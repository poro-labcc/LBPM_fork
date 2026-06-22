#include "models/SteadyPhaseModel.h"
#include "analysis/distance.h"
#include "common/ReadMicroCT.h"
#include <cmath>
#include <chrono>
#include <vector>

typedef Array<float> FloatArray;

extern "C" {
    void ScaLBL_D3Q19_AAeven_SteadyPhase(double *dist, int start, int finish,
                                         int Np, double Fx, double Fy, double Fz,
                                         double *ColorGrad, int *Phi, double tau_A, double tau_B, double *Velocity, double factor);

    void ScaLBL_D3Q19_AAodd_SteadyPhase(int *neighborList, double *dist,
                                        int start, int finish, int Np,
                                        double Fx, double Fy, double Fz,
                                        double *ColorGrad, int *Phi, double tau_A, double tau_B, double *Velocity, double factor);

    void SteadyComputeVelocity(double *dist, double *vel, double *pressure_out, int Np, double Fx, double Fy, double Fz, double *ColorGrad, double factor);
}

namespace {
    std::vector<double> get_gaussian_kernel(double sigma) {
        if (sigma <= 0.0) return {1.0};
        int radius = int(2.0 * sigma);
        int size = 2 * radius + 1;
        std::vector<double> k(size);
        double sum = 0.0;
        for (int i = -radius; i <= radius; i++) {
            double val = exp(-(i * i) / (2.0 * sigma * sigma));
            k[i + radius] = val;
            sum += val;
        }
        for (auto &v : k) v /= sum;
        return k;
    }
}

ScaLBL_SteadyPhaseModel::ScaLBL_SteadyPhaseModel(int RANK, int NP, const Utilities::MPI &COMM)
    : rank(RANK), nprocs(NP), Restart(false), timestep(0), timestepMax(0),
      ANALYSIS_INTERVAL(1000), BoundaryCondition(0), tau_A(1.0), tau_B(1.0),
      mu_A(0.0), mu_B(0.0), sigma(2.0), tolerance(1e-8), Fx(0.0), Fy(0.0), Fz(0.0), flux(0.0),
      din(1.0), dout(1.0), Nx(0), Ny(0), Nz(0), N(0), Np(0), nprocx(0), nprocy(0), nprocz(0),
      Lx(0), Ly(0), Lz(0), NeighborList(nullptr), fq(nullptr), Velocity(nullptr), Pressure(nullptr),
      Phi(nullptr), ColorGrad(nullptr), comm(COMM) {}

ScaLBL_SteadyPhaseModel::~ScaLBL_SteadyPhaseModel() {}

void ScaLBL_SteadyPhaseModel::ReadParams(string filename) {
    db = std::make_shared<Database>(filename);
    domain_db = db->getDatabase("Domain");
    steady_db = db->getDatabase("Steady");
    vis_db = db->getDatabase("Visualization");

    timestepMax = steady_db->getWithDefault<int>("timestepMax", 100000);
    ANALYSIS_INTERVAL = steady_db->getWithDefault<int>("analysis_interval", 1000);
    tolerance = steady_db->getWithDefault<double>("tolerance", 1.0e-8);
    tau_A = steady_db->getWithDefault<double>("tauA", 1.0);
    tau_B = steady_db->getWithDefault<double>("tauB", 1.0);
    sigma = steady_db->getWithDefault<double>("sigma", 2.0);
    factor = steady_db->getWithDefault<double>("factor", -2.0);

    if (steady_db->keyExists("din"))
        din = steady_db->getScalar<double>("din");
    if (steady_db->keyExists("dout"))
        dout = steady_db->getScalar<double>("dout");
    if (steady_db->keyExists("flux"))
        flux = steady_db->getScalar<double>("flux");

    if (steady_db->keyExists("F")) {
        Fx = steady_db->getVector<double>("F")[0];
        Fy = steady_db->getVector<double>("F")[1];
        Fz = steady_db->getVector<double>("F")[2];
    }

    if (steady_db->keyExists("BoundaryCondition")) {
        BoundaryCondition = steady_db->getScalar<int>("BC");
    } else if (domain_db->keyExists("BC")) {
        BoundaryCondition = domain_db->getScalar<int>("BC");
    }

    mu_A = (tau_A - 0.5) / 3.0;
    mu_B = (tau_B - 0.5) / 3.0;
}

void ScaLBL_SteadyPhaseModel::SetDomain() {
    Dm = std::shared_ptr<Domain>(new Domain(domain_db, comm));
    Mask = std::shared_ptr<Domain>(new Domain(domain_db, comm));
    Nx = Dm->Nx; Ny = Dm->Ny; Nz = Dm->Nz;
    Lx = Dm->Lx; Ly = Dm->Ly; Lz = Dm->Lz;
    N = Nx * Ny * Nz;
    SignDistance.resize(Nx, Ny, Nz);

    for (int i = 0; i < Nx * Ny * Nz; i++)
        Dm->id[i] = 1;
    
    comm.barrier();
    Dm->CommInit();
    comm.barrier();

    rank = Dm->rank();
    nprocx = Dm->nprocx();
    nprocy = Dm->nprocy();
    nprocz = Dm->nprocz();
}

void ScaLBL_SteadyPhaseModel::ReadInput() {

    sprintf(LocalRankString, "%05d", Dm->rank());
    sprintf(LocalRankFilename, "%s%s", "ID.", LocalRankString);
    sprintf(LocalRestartFile, "%s%s", "Restart.", LocalRankString);

    if (domain_db->keyExists("Filename")) {
        auto Filename = domain_db->getScalar<std::string>("Filename");
        Mask->Decomp(Filename);
    } else if (domain_db->keyExists("GridFile")) {
        auto input_id = readMicroCT(*domain_db, comm);
        array<int, 3> size0 = {(int)input_id.size(0), (int)input_id.size(1),
                               (int)input_id.size(2)};
        ArraySize size1 = {(size_t)Mask->Nx, (size_t)Mask->Ny,
                           (size_t)Mask->Nz};
        ASSERT((int)size1[0] == size0[0] + 2 && (int)size1[1] == size0[1] + 2 &&
               (int)size1[2] == size0[2] + 2);
        fillHalo<signed char> fill(comm, Mask->rank_info, size0, {1, 1, 1}, 0,
                                   1);
        Array<signed char> id_view;
        id_view.viewRaw(size1, Mask->id.data());
        fill.copy(input_id, id_view);
        fill.fill(id_view);
    } else {
        Mask->ReadIDs();
    }

    Array<char> id_solid(Nx, Ny, Nz);
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int n = k * Nx * Ny + j * Nx + i;
                if (Mask->id[n] > 0)
                    id_solid(i, j, k) = 1;
                else
                    id_solid(i, j, k) = 0;
            }
        }
    }
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                SignDistance(i, j, k) = 2.0 * double(id_solid(i, j, k)) - 1.0;
            }
        }
    }
    if (rank == 0)
        printf("Initialized solid phase -- Converting to Signed Distance "
               "function \n");
    CalcDist(SignDistance, id_solid, *Dm);
    if (rank == 0)
        cout << "Domain set." << endl;
}


void ScaLBL_SteadyPhaseModel::Create() {
    int rank = Mask->rank();
    for (int i = 0; i < Nx * Ny * Nz; i++)
        Dm->id[i] = Mask->id[i];
    Mask->CommInit();
    Np = Mask->PoreCount();
    if (rank == 0)
        printf("Create ScaLBL_Communicator \n");
    ScaLBL_Comm =
        std::shared_ptr<ScaLBL_Communicator>(new ScaLBL_Communicator(Mask));

    int Npad = (Np / 16 + 2) * 16;
    if (rank == 0)
        printf("Set up memory efficient layout \n");
    Map.resize(Nx, Ny, Nz);
    Map.fill(-2);
    auto neighborList = new int[18 * Npad];
    Np = ScaLBL_Comm->MemoryOptimizedLayoutAA(Map, neighborList,
                                              Mask->id.data(), Np, 1);
    comm.barrier();

    if (rank == 0)
        printf("Allocating distributions \n");
    dist_mem_size = Np * sizeof(double);
    neighborSize = 18 * (Np * sizeof(int));

    ScaLBL_AllocateDeviceMemory((void **)&Phi, Np * sizeof(int));
    ScaLBL_AllocateDeviceMemory((void **)&ColorGrad, 3 * dist_mem_size);
    ComputeNormals();

    ScaLBL_AllocateDeviceMemory((void **)&NeighborList, neighborSize);
    ScaLBL_AllocateDeviceMemory((void **)&fq, 19 * dist_mem_size);
    ScaLBL_AllocateDeviceMemory((void **)&Pressure, sizeof(double) * Np);
    ScaLBL_AllocateDeviceMemory((void **)&Velocity, 3 * sizeof(double) * Np);


    if (rank == 0)
        printf("Setting up device map and neighbor list \n");
    
    ScaLBL_CopyToDevice(NeighborList, neighborList, neighborSize);
    comm.barrier();
    delete [] neighborList;
    double MLUPS = ScaLBL_Comm->GetPerformance(NeighborList, fq, Np);
    printf("  MLPUS=%f from rank %i\n", MLUPS, rank);
}


void ScaLBL_SteadyPhaseModel::ComputeNormals() {
    IntArray HostPhi(Nx, Ny, Nz);
    DoubleArray NormX(Nx, Ny, Nz), NormY(Nx, Ny, Nz), NormZ(Nx, Ny, Nz);

    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int n = k * Nx * Ny + j * Nx + i;
                if (Mask->id[n] == 0) HostPhi(i, j, k) = 0.0;
                else HostPhi(i, j, k) = (Mask->id[n] == 1) ? 1.0 : -1.0;
            }
        }
    }

    DoubleArray SmoothPhi(Nx, Ny, Nz);

    if (sigma > 0.0) {
        auto k_gauss = get_gaussian_kernel(sigma);
        int r = k_gauss.size() / 2;
        
        int R = r + 1; 
        int exNx = (Nx - 2) + 2 * R;
        int exNy = (Ny - 2) + 2 * R;
        int exNz = (Nz - 2) + 2 * R;

        DoubleArray ExtendedPhi(exNx, exNy, exNz);
        ExtendedPhi.fill(0.0);

        for (int k = 1; k < Nz - 1; k++) {
            for (int j = 1; j < Ny - 1; j++) {
                for (int i = 1; i < Nx - 1; i++) {
                    ExtendedPhi(i - 1 + R, j - 1 + R, k - 1 + R) = HostPhi(i, j, k);
                }
            }
        }

        int pFlag = (BoundaryCondition == 0) ? 1 : 0;
        std::array<int, 3> subSize = {Nx - 2, Ny - 2, Nz - 2};
        std::array<int, 3> haloSize = {R, R, R};
        fillHalo<double> haloSync(Dm->Comm, Dm->rank_info, subSize, haloSize, pFlag, 1);
        haloSync.fill(ExtendedPhi);

        DoubleArray Tmp1(exNx, exNy, exNz);
        DoubleArray Tmp2(exNx, exNy, exNz);
        DoubleArray SmoothPhiExt(exNx, exNy, exNz);

        for (int k = 0; k < exNz; k++) {
            for (int j = 0; j < exNy; j++) {
                for (int i = R - 1; i <= exNx - R; i++) {
                    double sum = 0;
                    for (int d = -r; d <= r; d++) {
                        sum += k_gauss[d + r] * ExtendedPhi(i + d, j, k);
                    }
                    Tmp1(i, j, k) = sum;
                }
            }
        }

        for (int k = 0; k < exNz; k++) {
            for (int j = R - 1; j <= exNy - R; j++) {
                for (int i = R - 1; i <= exNx - R; i++) {
                    double sum = 0;
                    for (int d = -r; d <= r; d++) {
                        sum += k_gauss[d + r] * Tmp1(i, j + d, k);
                    }
                    Tmp2(i, j, k) = sum;
                }
            }
        }

        for (int k = R - 1; k <= exNz - R; k++) {
            for (int j = R - 1; j <= exNy - R; j++) {
                for (int i = R - 1; i <= exNx - R; i++) {
                    double sum = 0;
                    for (int d = -r; d <= r; d++) {
                        sum += k_gauss[d + r] * Tmp2(i, j, k + d);
                    }
                    SmoothPhiExt(i, j, k) = sum;
                }
            }
        }

        for (int k = 0; k < Nz; k++) {
            for (int j = 0; j < Ny; j++) {
                for (int i = 0; i < Nx; i++) {
                    SmoothPhi(i, j, k) = SmoothPhiExt(i - 1 + R, j - 1 + R, k - 1 + R);
                }
            }
        }

    } else {
        for (int k = 0; k < Nz; k++)
            for (int j = 0; j < Ny; j++)
                for (int i = 0; i < Nx; i++)
                    SmoothPhi(i, j, k) = HostPhi(i, j, k);
    }

    int cx[19] = {0, 1, -1, 0,  0, 0,  0, 1, -1,  1, -1, 1, -1,  1, -1, 0, 0,   0,  0};
    int cy[19] = {0, 0,  0, 1, -1, 0,  0, 1, -1, -1,  1, 0,  0,  0,  0, 1, -1,  1, -1};
    int cz[19] = {0, 0,  0, 0,  0, 1, -1, 0,  0,  0,  0, 1, -1, -1,  1, 1, -1, -1,  1};
    double w[19] = {1./3., 1./18.,1./18.,1./18.,1./18.,1./18.,1./18., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36.};

    for (int k = 1; k < Nz - 1; k++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                int n = k * Nx * Ny + j * Nx + i;
                if (Mask->id[n] == 0) continue;
                
                double nx = 0, ny = 0, nz = 0;
                for(int d = 1; d < 19; d++) {
                    double p = SmoothPhi(i - cx[d], j - cy[d], k - cz[d]);
                    nx += cx[d] * w[d] * 18.0 * p; 
                    ny += cy[d] * w[d] * 18.0 * p; 
                    nz += cz[d] * w[d] * 18.0 * p;
                }
                
                double mag = sqrt(nx*nx + ny*ny + nz*nz);
                if (mag > 1e-12) {
                    NormX(i,j,k) = nx / mag; 
                    NormY(i,j,k) = ny / mag; 
                    NormZ(i,j,k) = nz / mag;
                } else {
                    NormX(i,j,k) = 0.0;
                    NormY(i,j,k) = 0.0;
                    NormZ(i,j,k) = 0.0;
                }
            }
        }
    }

    for (int k = 1; k < Nz - 1; k++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                int n = k * Nx * Ny + j * Nx + i;
                if (Mask->id[n] == 0) continue;
                
                bool near_f_opp = false;
                double p_val = HostPhi(i, j, k);
                
                for(int d = 1; d < 19; d++) {
                    double neighbor_p = HostPhi(i - cx[d], j - cy[d], k - cz[d]);
                    if ((p_val == 1.0 && neighbor_p == -1.0) || (p_val == -1.0 && neighbor_p == 1.0)) {
                        near_f_opp = true;
                        break;
                    }
                }

                if (!near_f_opp) {
                    NormX(i,j,k) = 0.0;
                    NormY(i,j,k) = 0.0;
                    NormZ(i,j,k) = 0.0;
                }
            }
        }
    }

    std::vector<int> h_Phi(Np, 0);
    std::vector<double> h_Grad(3*Np, 0.0);

    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int id = Map(i,j,k);
                if (id >= 0) {
                    h_Phi[id] = (int)HostPhi(i,j,k);
                    h_Grad[id] = NormX(i,j,k);
                    h_Grad[Np+id] = NormY(i,j,k);
                    h_Grad[2*Np+id] = NormZ(i,j,k);
                }
            }
        }
    }
    ScaLBL_CopyToDevice(Phi, h_Phi.data(), Np * sizeof(int));
    ScaLBL_CopyToDevice(ColorGrad, h_Grad.data(), 3*dist_mem_size);
}

void ScaLBL_SteadyPhaseModel::Initialize() {
    if (rank == 0)
        printf("Initializing distributions \n");
    ScaLBL_D3Q19_Init(fq, Np);
}

void ScaLBL_SteadyPhaseModel::Run() {

    Minkowski Morphology(Mask);

    if (rank == 0) {
        bool WriteHeader = false;
        FILE *log_file = fopen("Permeability.csv", "r");
        if (log_file != NULL)
            fclose(log_file);
        else
            WriteHeader = true;

        if (WriteHeader) {
            log_file = fopen("Permeability.csv", "a+");
            fprintf(log_file, "time satA Fx Fy Fz muA muB keffA(Da) keffA*(Da) keffB(Da) keffB*(Da)\n");
            fclose(log_file);
        }
    }



    if (rank == 0) {
        printf("Beginning AA timesteps, timestepMax = %i \n", timestepMax);
        printf("********************************************************\n");
    }

    timestep = 0;
    double error = 1.0;
    double flow_rate_previous_A = 0.0, flow_rate_previous_B = 0.0;

    auto t1 = std::chrono::system_clock::now();
    while (timestep < timestepMax && error > tolerance) {
        timestep++;
        ScaLBL_Comm->SendD3Q19AA(fq);
        ScaLBL_D3Q19_AAodd_SteadyPhase(NeighborList, fq, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np, Fx, Fy, Fz, ColorGrad, Phi, tau_A, tau_B, Velocity, factor);
        ScaLBL_Comm->RecvD3Q19AA(fq);
        if (BoundaryCondition == 3) {
            ScaLBL_Comm->D3Q19_Pressure_BC_z(NeighborList, fq, din, timestep);
            ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
        } else if (BoundaryCondition == 4) {
            din = ScaLBL_Comm->D3Q19_Flux_BC_z(NeighborList, fq, flux, timestep);
            ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
        } else if (BoundaryCondition == 5) {
            ScaLBL_Comm->D3Q19_Reflection_BC_z(fq);
            ScaLBL_Comm->D3Q19_Reflection_BC_Z(fq);
        }
        ScaLBL_D3Q19_AAodd_SteadyPhase(NeighborList, fq, 0, ScaLBL_Comm->LastExterior(), Np, Fx, Fy, Fz, ColorGrad, Phi, tau_A, tau_B, Velocity, factor);
        
        ScaLBL_DeviceBarrier();
        comm.barrier();

        timestep++;
        ScaLBL_Comm->SendD3Q19AA(fq);
        ScaLBL_D3Q19_AAeven_SteadyPhase(fq, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np, Fx, Fy, Fz, ColorGrad, Phi, tau_A, tau_B, Velocity, factor);
        ScaLBL_Comm->RecvD3Q19AA(fq);
        if (BoundaryCondition == 3) {
            ScaLBL_Comm->D3Q19_Pressure_BC_z(NeighborList, fq, din, timestep);
            ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
        } else if (BoundaryCondition == 4) {
            din =
                ScaLBL_Comm->D3Q19_Flux_BC_z(NeighborList, fq, flux, timestep);
            ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
        } else if (BoundaryCondition == 5) {
            ScaLBL_Comm->D3Q19_Reflection_BC_z(fq);
            ScaLBL_Comm->D3Q19_Reflection_BC_Z(fq);
        }
        ScaLBL_D3Q19_AAeven_SteadyPhase(fq, 0, ScaLBL_Comm->LastExterior(), Np, Fx, Fy, Fz, ColorGrad, Phi, tau_A, tau_B, Velocity, factor);
        
        ScaLBL_DeviceBarrier();
        comm.barrier();

        if (timestep % ANALYSIS_INTERVAL == 0) {
            SteadyComputeVelocity(fq, Velocity, Pressure, Np, Fx, Fy, Fz, ColorGrad, factor);
            ScaLBL_DeviceBarrier();
            comm.barrier();

            FloatArray Velocity_x(Nx, Ny, Nz);
            FloatArray Velocity_y(Nx, Ny, Nz);
            FloatArray Velocity_z(Nx, Ny, Nz);
            
            Velocity_x.fill(0.0f);
            Velocity_y.fill(0.0f);
            Velocity_z.fill(0.0f);

            for (int k = 0; k < Nz; k++) {
                for (int j = 0; j < Ny; j++) {
                    for (int i = 0; i < Nx; i++) {
                        int id = Map(i, j, k);
                        if (id >= 0) {
                            Velocity_x(i, j, k) = (float)Velocity[id];
                            Velocity_y(i, j, k) = (float)Velocity[Np + id];
                            Velocity_z(i, j, k) = (float)Velocity[2 * Np + id];
                        }
                    }
                }
            }

            double vax_loc_A = 0.0, vay_loc_A = 0.0, vaz_loc_A = 0.0;
            double count_loc_A = 0.0;
            
            double vax_loc_B = 0.0, vay_loc_B = 0.0, vaz_loc_B = 0.0;
            double count_loc_B = 0.0;

            for (int k = 1; k < Nz - 1; k++) {
                for (int j = 1; j < Ny - 1; j++) {
                    for (int i = 1; i < Nx - 1; i++) {
                        int n = k * Nx * Ny + j * Nx + i;
                        if (Mask->id[n] == 1) {
                            vax_loc_A += Velocity_x(i, j, k);
                            vay_loc_A += Velocity_y(i, j, k);
                            vaz_loc_A += Velocity_z(i, j, k);
                            count_loc_A += 1.0;
                        } else if (Mask->id[n] == 2) {
                            vax_loc_B += Velocity_x(i, j, k);
                            vay_loc_B += Velocity_y(i, j, k);
                            vaz_loc_B += Velocity_z(i, j, k);
                            count_loc_B += 1.0;
                        }
                    }
                }
            }

            double vax_A = Dm->Comm.sumReduce(vax_loc_A);
            double vay_A = Dm->Comm.sumReduce(vay_loc_A);
            double vaz_A = Dm->Comm.sumReduce(vaz_loc_A);
            double count_A  = Dm->Comm.sumReduce(count_loc_A);
            
            double vax_B = Dm->Comm.sumReduce(vax_loc_B);
            double vay_B = Dm->Comm.sumReduce(vay_loc_B);
            double vaz_B = Dm->Comm.sumReduce(vaz_loc_B);
            double count_B  = Dm->Comm.sumReduce(count_loc_B);

            vax_A  = (count_A > 0) ? vax_A / count_A : 0.0;
            vay_A  = (count_A > 0) ? vay_A / count_A : 0.0;
            vaz_A  = (count_A > 0) ? vaz_A / count_A : 0.0;

            vax_B  = (count_B > 0) ? vax_B / count_B : 0.0;
            vay_B  = (count_B > 0) ? vay_B / count_B : 0.0;
            vaz_B  = (count_B > 0) ? vaz_B / count_B : 0.0;

            double force_mag = sqrt(Fx * Fx + Fy * Fy + Fz * Fz);
            double dir_x = Fx / force_mag;
            double dir_y = Fy / force_mag;
            double dir_z = Fz / force_mag;
            if (force_mag == 0.0) {
                dir_x = 0.0;
                dir_y = 0.0;
                dir_z = 1.0;
                force_mag = 1.0;
            }
            double flow_rate_A = (vax_A * dir_x + vay_A * dir_y + vaz_A * dir_z);
            double flow_rate_B = (vax_B * dir_x + vay_B * dir_y + vaz_B * dir_z);

            double error_A = fabs(flow_rate_A - flow_rate_previous_A) / fabs(flow_rate_A);
            double error_B = fabs(flow_rate_B - flow_rate_previous_B) / fabs(flow_rate_B);

            if (std::isnan(error_A)) {
                error = error_B;
            }
            else if (std::isnan(error_B)) {
                error = error_A;
            }
            else {
                error = error_A > error_B ? error_A : error_B;
            }

            flow_rate_previous_A = flow_rate_A;
            flow_rate_previous_B = flow_rate_B;


            double total_nodes_global = Dm->Comm.sumReduce((double)(Nx-2)*(Ny-2)*(Nz-2));

            double frac_A  = count_A / total_nodes_global;
            double frac_B = count_B / total_nodes_global;
            double sat_A = count_A / (count_A + count_B);

            double h = Dm->voxel_length;

            double absperm_A = h * h * mu_A * frac_A * flow_rate_A / force_mag;
            // absperm_A *= 1013.0; 

            double absperm_B = h * h * mu_B * frac_B * flow_rate_B / force_mag;
            // absperm_B *= 1013.0;

            if (rank == 0) {
                printf("Step %d | kA = %.5f | kB = %.5f | error = %.5f\n", timestep, absperm_A, absperm_B, error);
                FILE *log_file = fopen("Permeability.csv", "a");
                if (log_file != NULL) {
                    fprintf(log_file,
                        "%i %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g\n",
                        timestep, sat_A, Fx, Fy, Fz, mu_A, mu_B, absperm_A, absperm_A * Mask->Porosity(), absperm_B, absperm_B * Mask->Porosity());
                    fclose(log_file);
                }
            }
        }
    }

    if (rank == 0) printf("--------------------------------------------------------\n");
    auto t2 = std::chrono::system_clock::now();
    double cputime = std::chrono::duration<double>(t2 - t1).count() / timestep;
    double MLUPS = double(Np) / cputime / 1000000.0;

    if (rank == 0) {
        printf("********************************************************\n");
        printf("CPU time = %f \n", cputime);
        printf("Lattice update rate (per core)= %f MLUPS \n", MLUPS);
        MLUPS *= nprocs;
        printf("Lattice update rate (total)= %f MLUPS \n", MLUPS);
        printf("********************************************************\n");
    }
}

void ScaLBL_SteadyPhaseModel::VelocityField() {
    auto format = vis_db->getWithDefault<std::string>("format", "silo");
    if (!vis_db->getWithDefault<bool>("write_silo", false)) return;

    bool save_phase    = vis_db->getWithDefault<bool>("save_phase", true);
    bool save_velocity = vis_db->getWithDefault<bool>("save_velocity", true);
    bool save_normal   = vis_db->getWithDefault<bool>("save_normal", false);
    bool save_signdist = vis_db->getWithDefault<bool>("save_signdistance", false);
    bool save_pressure = vis_db->getWithDefault<bool>("save_pressure", true);

    IO::initialize("", format, false);
    std::vector<IO::MeshDataStruct> visData(1);
    visData[0].meshName = "domain";
    visData[0].mesh = std::make_shared<IO::DomainMesh>(Dm->rank_info, Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2, Dm->Lx, Dm->Ly, Dm->Lz);
    
    int var_index = 0;
    auto add_variable = [&](const std::string& name) {
        auto var = std::make_shared<IO::Variable>();
        var->name = name;
        var->type = IO::VariableType::VolumeVariable;
        var->dim = 1;
        var->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
        visData[0].vars.push_back(var);
        return var_index++;
    };

    int idx_phase = -1, idx_velx = -1, idx_vely = -1, idx_velz = -1;
    int idx_normx = -1, idx_normy = -1, idx_normz = -1;
    int idx_signdist = -1, idx_pressure = -1;

    if (save_phase)     idx_phase = add_variable("Phase");
    if (save_velocity) {
        idx_velx = add_variable("Velocity_x");
        idx_vely = add_variable("Velocity_y");
        idx_velz = add_variable("Velocity_z");
    }
    if (save_normal) {
        idx_normx = add_variable("Normal_x");
        idx_normy = add_variable("Normal_y");
        idx_normz = add_variable("Normal_z");
    }
    if (save_signdist)  idx_signdist = add_variable("SignDistance");
    if (save_pressure)  idx_pressure = add_variable("Pressure");

    fillHalo<double> fillData(Dm->Comm, Dm->rank_info, {Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2}, {1, 1, 1}, 0, 1);
    
    if (save_phase && idx_phase >= 0) {
        Array<unsigned char> PhaseHost(Nx, Ny, Nz);
        std::vector<int> h_Phi_1D(Np);
        ScaLBL_CopyToHost(h_Phi_1D.data(), Phi, Np * sizeof(int));
        
        PhaseHost.fill(0);
        for (int k = 0; k < Nz; k++) {
            for (int j = 0; j < Ny; j++) {
                for (int i = 0; i < Nx; i++) {
                    int id = Map(i, j, k);
                    if (id >= 0) {
                        int raw_phi = h_Phi_1D[id];
                        if (raw_phi == -1) PhaseHost(i, j, k) = 2;
                        else if (raw_phi == 1) PhaseHost(i, j, k) = 1;
                        else PhaseHost(i, j, k) = 0;
                    }
                }
            }
        }
        fillData.copy(PhaseHost, visData[0].vars[idx_phase]->data);
    } 

    if (save_velocity) {
        {
            DoubleArray Vel_x(Nx, Ny, Nz);
            ScaLBL_Comm->RegularLayout(Map, &Velocity[0], Vel_x);
            fillData.copy(Vel_x, visData[0].vars[idx_velx]->data);
        }
        {
            DoubleArray Vel_y(Nx, Ny, Nz);
            ScaLBL_Comm->RegularLayout(Map, &Velocity[Np], Vel_y);
            fillData.copy(Vel_y, visData[0].vars[idx_vely]->data);
        }
        {
            DoubleArray Vel_z(Nx, Ny, Nz);
            ScaLBL_Comm->RegularLayout(Map, &Velocity[2 * Np], Vel_z);
            fillData.copy(Vel_z, visData[0].vars[idx_velz]->data);
        }
    }

    if (save_normal) {
        {
            DoubleArray Norm_x(Nx, Ny, Nz);
            ScaLBL_Comm->RegularLayout(Map, &ColorGrad[0], Norm_x);
            fillData.copy(Norm_x, visData[0].vars[idx_normx]->data);
        }
        {
            DoubleArray Norm_y(Nx, Ny, Nz);
            ScaLBL_Comm->RegularLayout(Map, &ColorGrad[Np], Norm_y);
            fillData.copy(Norm_y, visData[0].vars[idx_normy]->data);
        }
        {
            DoubleArray Norm_z(Nx, Ny, Nz);
            ScaLBL_Comm->RegularLayout(Map, &ColorGrad[2 * Np], Norm_z);
            fillData.copy(Norm_z, visData[0].vars[idx_normz]->data);
        }
    }

    if (save_signdist && idx_signdist >= 0) {
        fillData.copy(SignDistance, visData[0].vars[idx_signdist]->data);
    }

    if (save_pressure && idx_pressure >= 0) {
        {
            DoubleArray Pressure_field(Nx, Ny, Nz);
            ScaLBL_Comm->RegularLayout(Map, &Pressure[0], Pressure_field);
            fillData.copy(Pressure_field, visData[0].vars[idx_pressure]->data);
        }
    }
    
    IO::writeData(timestep, visData, Dm->Comm);
}