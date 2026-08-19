/*
  Copyright 2013--2018 James E. McClure, Virginia Polytechnic & State University
  Copyright Equnior ASA

  This file is part of the Open Porous Media project (OPM).
  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
/*
 * Multi-relaxation time LBM Model
 */
#include "models/MRTModel.h"
#include "analysis/distance.h"
#include "common/ReadMicroCT.h"


ScaLBL_MRTModel::ScaLBL_MRTModel(int RANK, int NP, const Utilities::MPI &COMM)
    : rank(RANK), nprocs(NP), Restart(0), timestep(0), timestepMax(0), tau(0),
      Fx(0), Fy(0), Fz(0), flux(0), din(0), dout(0), dp(0), mu(0), Nx(0), Ny(0),
      Nz(0), N(0), Np(0), nprocx(0), nprocy(0), nprocz(0), BoundaryCondition(0), Lx(0),
      Ly(0), Lz(0), comm(COMM) {}


ScaLBL_MRTModel::~ScaLBL_MRTModel() {}


void ScaLBL_MRTModel::ReadParams(string filename) {
    // read the input database
    db = std::make_shared<Database>(filename);
    domain_db = db->getDatabase("Domain");
    mrt_db = db->getDatabase("MRT");
    vis_db = db->getDatabase("Visualization");
    ana_db = db->getDatabase("Analysis");

    tau                 = 1.0;
    timestepMax         = 100000;
    ANALYSIS_INTERVAL   = 1000;
    VISUAL_INTERVAL     = 100001;
    save_velocity       = true;
    save_pressure       = false;
    tolerance           = 1.0e-8;
    Fx = Fy = 0.0;
    Fz = 1.0e-5;
    dout    = 1.0;
    din     = 1.0;
    dp      = 0.0;

    // Color Model parameters
    if (mrt_db->keyExists("timestepMax")) {
        timestepMax = mrt_db->getScalar<int>("timestepMax");
    }
    if (ana_db->keyExists("analysis_interval")) {
        ANALYSIS_INTERVAL = ana_db->getScalar<int>("analysis_interval");
    }
    if (ana_db->keyExists("visualization_interval")) {
        VISUAL_INTERVAL = ana_db->getScalar<int>("visualization_interval");
    }
    else{
        VISUAL_INTERVAL = timestepMax+1; // Saving only in the end by default
    }
    if (vis_db->keyExists("save_pressure")) {
        save_pressure = vis_db->getScalar<bool>("save_pressure");
    }
    if (vis_db->keyExists("save_velocity")) {
        save_velocity = vis_db->getScalar<bool>("save_velocity");
    }
    if (mrt_db->keyExists("tolerance")) {
        tolerance = mrt_db->getScalar<double>("tolerance");
    }
    if (mrt_db->keyExists("tau")) {
        tau = mrt_db->getScalar<double>("tau");
    }
    if (mrt_db->keyExists("F")) {
        Fx = mrt_db->getVector<double>("F")[0];
        Fy = mrt_db->getVector<double>("F")[1];
        Fz = mrt_db->getVector<double>("F")[2];
    }
    if (mrt_db->keyExists("Restart")) {
        Restart = mrt_db->getScalar<bool>("Restart");
    }
    if (mrt_db->keyExists("dp")) {
        dp = mrt_db->getScalar<double>("dp");
    }
    if (mrt_db->keyExists("din")) {
        din = mrt_db->getScalar<double>("din");
    }
    if (mrt_db->keyExists("dout")) {
        dout = mrt_db->getScalar<double>("dout");
    }
    if (mrt_db->keyExists("flux")) {
        flux = mrt_db->getScalar<double>("flux");
    }

    // Read domain parameters
    if (mrt_db->keyExists("BoundaryCondition")) {
        BoundaryCondition = mrt_db->getScalar<int>("BC");
    } else if (domain_db->keyExists("BC")) {
        BoundaryCondition = domain_db->getScalar<int>("BC");
    }

    mu       = (tau - 0.5) / 3.0;
    rlx_setA = 1.0 / tau;
    rlx_setB = 8.f * (2.f - rlx_setA) / (8.f - rlx_setA);
}
void ScaLBL_MRTModel::SetDomain() {
    Dm = std::shared_ptr<Domain>(
        new Domain(domain_db, comm)); // full domain for analysis
    Mask = std::shared_ptr<Domain>(
        new Domain(domain_db, comm)); // mask domain removes immobile phases

    // domain parameters
    Nx = Dm->Nx;
    Ny = Dm->Ny;
    Nz = Dm->Nz;
    Lx = Dm->Lx;
    Ly = Dm->Ly;
    Lz = Dm->Lz;

    N = Nx * Ny * Nz;
    Distance.resize(Nx, Ny, Nz);
    Velocity_x.resize(Nx, Ny, Nz);
    Velocity_y.resize(Nx, Ny, Nz);
    Velocity_z.resize(Nx, Ny, Nz);
    Pressure_f.resize(Nx, Ny, Nz);

    for (int i = 0; i < Nx * Ny * Nz; i++)
        Dm->id[i] = 1; // initialize this way
    //Averages = std::shared_ptr<TwoPhase> ( new TwoPhase(Dm) ); // TwoPhase analysis object
    comm.barrier();
    Dm->CommInit();
    comm.barrier();

    rank = Dm->rank();
    nprocx = Dm->nprocx();
    nprocy = Dm->nprocy();
    nprocz = Dm->nprocz();
}

void ScaLBL_MRTModel::ReadInput() {

    sprintf(LocalRankString, "%05d", Dm->rank());
    sprintf(LocalRankFilename, "%s%s", "ID.", LocalRankString);
    sprintf(LocalRestartFile, "%s%s", "Restart.", LocalRankString);

    if (domain_db->keyExists("Filename")) {
        auto Filename = domain_db->getScalar<std::string>("Filename");
        Mask->Decomp(Filename);
    } else if (domain_db->keyExists("GridFile")) {
        // Read the local domain data
        auto input_id = readMicroCT(*domain_db, comm);
        // Fill the halo (assuming GCW of 1)
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

    // Generate the signed distance map
    // Initialize the domain and communication
    Array<char> id_solid(Nx, Ny, Nz);
    // Solve for the position of the solid phase
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int n = k * Nx * Ny + j * Nx + i;
                // Initialize the solid phase
                if (Mask->id[n] > 0)
                    id_solid(i, j, k) = 1;
                else
                    id_solid(i, j, k) = 0;
            }
        }
    }
    // Initialize the signed distance function
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                // Initialize distance to +/- 1
                Distance(i, j, k) = 2.0 * double(id_solid(i, j, k)) - 1.0;
            }
        }
    }
    //	MeanFilter(Averages->SDs);
    if (rank == 0)
        printf("Initialized solid phase -- Converting to Signed Distance "
               "function \n");
    CalcDist(Distance, id_solid, *Dm);
    if (rank == 0)
        cout << "Domain set." << endl;
}

void ScaLBL_MRTModel::Create() {
    /*
	 *  This function creates the variables needed to run a LBM 
	 */
    int rank = Mask->rank();
    //.........................................................
    // Initialize communication structures in averaging domain
    for (int i = 0; i < Nx * Ny * Nz; i++)
        Dm->id[i] = Mask->id[i];
    Mask->CommInit();
    Np = Mask->PoreCount();
    //...........................................................................
    if (rank == 0)
        printf("Create ScaLBL_Communicator \n");
    // Create a communicator for the device (will use optimized layout)
    // ScaLBL_Communicator ScaLBL_Comm(Mask); // original
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

    //...........................................................................
    //                MAIN  VARIABLES ALLOCATED HERE
    //...........................................................................
    // LBM variables
    if (rank == 0)
        printf("Allocating distributions \n");
    //......................device distributions.................................
    size_t dist_mem_size = Np * sizeof(double);
    size_t neighborSize = 18 * (Np * sizeof(int));
    //...........................................................................
    ScaLBL_AllocateDeviceMemory((void **)&NeighborList, neighborSize);
    ScaLBL_AllocateDeviceMemory((void **)&fq, 19 * dist_mem_size);
    ScaLBL_AllocateDeviceMemory((void **)&Pressure, sizeof(double) * Np);
    ScaLBL_AllocateDeviceMemory((void **)&Velocity, 3 * sizeof(double) * Np);
    //...........................................................................
    // Update GPU data structures
    if (rank == 0)
        printf("Setting up device map and neighbor list \n");
    // copy the neighbor list
    ScaLBL_CopyToDevice(NeighborList, neighborList, neighborSize);
    comm.barrier();
    double MLUPS = ScaLBL_Comm->GetPerformance(NeighborList, fq, Np);
    printf("  MLPUS=%f from rank %i\n", MLUPS, rank);
}


void ScaLBL_MRTModel::Initialize() {
    /*
	 * This function initializes model with equilibium distribution
     * for a null velocity field
	 */
    if (rank == 0)
        printf("Initializing distributions \n");
    ScaLBL_D3Q19_Init(fq, Np);
}

void ScaLBL_MRTModel::Initialize_Dist() {
    /*
	 * This function initializes model with custom distributions
	 */
    
     // Initialize distributions as no velocity Equilibrium
    ScaLBL_D3Q19_Init(fq, Np);


    // If there is a .raw file
    char raw_filename[256];
    sprintf(raw_filename, "StartF.%05d.raw", rank);
    std::ifstream binaryFile(raw_filename, std::ios::binary);

    if (binaryFile.good() && Start) {
        // Remove halo extra voxels (added in Domain class)
        unsigned int nx = Nx - 2;
        unsigned int ny = Ny - 2;
        unsigned int nz = Nz - 2;
        unsigned int n_items = 19;

        // Announces the start of the process
        if (rank == 0) printf("Reading start file (19 items float64): %s (%dx%dx%d)\n", raw_filename, nx, ny, nz);

        // Allocate buffer for the domain
        // n_items doubles per voxel (Ux, Uy, Uz, Pr)
        size_t total_voxels = (size_t)nx * ny * nz;
        std::vector<double> file_data(total_voxels * n_items);
        binaryFile.read(reinterpret_cast<char*>(file_data.data()), file_data.size() * sizeof(double));
        binaryFile.close();

        // Allocate auxiliary distributions
        double* temp_fq = new double[19 * Np];
        memset(temp_fq, 0, 19 * Np * sizeof(double)); // Initialize it with zeros


        for (unsigned int k = 1; k < nz+1; k++) {
            for (unsigned int j = 1; j < ny+1; j++) {
                for (unsigned int i = 1; i < nx+1; i++) {


                    // Get index in flatten array (solid + fluid + extra cells)
                    int cell_offset = Map(i, j, k);

                    // If is a fluid cell
                    if (cell_offset >= 0) {
                        // Calculate flatten index for file_data (discounting offset of 1)
                        size_t flat_idx = ((size_t)(k - 1) * nx * ny + (size_t)(j - 1) * nx + (i - 1)) * n_items;
                        // SET EQUILIBRIUM
                        for (int q = 0; q < 19; q++) {
                            double fq_data = file_data[flat_idx + q];
                            // Save in flatten array
                            temp_fq[q*Np + cell_offset] = fq_data;
                        }
                    }
                }
            }
        }

        ScaLBL_CopyToDevice(fq, temp_fq, 19 * Np * sizeof(double));
        delete[] temp_fq;
    }
    else {
        if (rank == 0) printf("No start file. Initializing with null velocity case.\n");
    }

    // Update Velocity state from fq
    ScaLBL_D3Q19_Momentum(fq,Velocity,Np);
    ScaLBL_DeviceBarrier();
    comm.barrier();
    ScaLBL_Comm->RegularLayout(Map, &Velocity[0   ], Velocity_x);
    ScaLBL_Comm->RegularLayout(Map, &Velocity[Np  ], Velocity_y);
    ScaLBL_Comm->RegularLayout(Map, &Velocity[2*Np], Velocity_z);

    // Update Pressure State from fq
    ScaLBL_D3Q19_Pressure(fq, Pressure, Np);    // Calculate Pressure Field
    ScaLBL_DeviceBarrier();                     // Sync
    comm.barrier();                             // Sync
    ScaLBL_Comm->RegularLayout(Map, &Pressure[0   ], Pressure_f);  // Transform Pressure Field in 3D domain

}

void ScaLBL_MRTModel::Initialize_fEq() {
    //
	// This function initializes model with equilibrium distributions
    // given by custom velocity and pressure fields
	//
    // Initialize distributions as no velocity Equilibrium
    ScaLBL_D3Q19_Init(fq, Np);


    // If there is a .raw file
    char raw_filename[256];
    sprintf(raw_filename, "Start.%05d.raw", rank);
    std::ifstream binaryFile(raw_filename, std::ios::binary);

    if (binaryFile.good() && Start) {
        // Remove halo extra voxels (added in Domain class)
        unsigned int nx = Nx - 2;
        unsigned int ny = Ny - 2;
        unsigned int nz = Nz - 2;
        unsigned int n_items = 4;

        // Announces the start of the process
        if (rank == 0) printf("Reading start file: %s (%dx%dx%d)\n", raw_filename, nx, ny, nz);

        // --- D3Q19 CONSTANTS ---

        // Allocate buffer for the domain
        // n_items doubles per voxel (Ux, Uy, Uz, Pr)
        size_t total_voxels = (size_t)nx * ny * nz;
        std::vector<double> file_data(total_voxels * n_items);
        binaryFile.read(reinterpret_cast<char*>(file_data.data()), file_data.size() * sizeof(double));
        binaryFile.close();

        // Allocate auxiliary distributions
        double* temp_fq = new double[19 * Np];
        memset(temp_fq, 0, 19 * Np * sizeof(double)); // Initialize it with zeros

        // MRT constants
        constexpr double mrt_V1 = 0.05263157894736842;
        constexpr double mrt_V2 = 0.012531328320802;
        constexpr double mrt_V3 = 0.04761904761904762;
        constexpr double mrt_V4 = 0.004594820384294068;
        constexpr double mrt_V5 = 0.01587301587301587;
        constexpr double mrt_V6 = 0.0555555555555555555555555;
        constexpr double mrt_V7 = 0.02777777777777778;
        constexpr double mrt_V8 = 0.08333333333333333;
        constexpr double mrt_V9 = 0.003341687552213868;
        constexpr double mrt_V10 = 0.003968253968253968;
        constexpr double mrt_V11 = 0.01388888888888889;
        constexpr double mrt_V12 = 0.04166666666666666;

        // For each cell
        for (unsigned int k = 1; k < nz+1; k++) {
            for (unsigned int j = 1; j < ny+1; j++) {
                for (unsigned int i = 1; i < nx+1; i++) {


                    // Get index in flatten array (solid + fluid + extra cells)
                    int cell_offset = Map(i, j, k);


                    // If is a fluid cell
                    if (cell_offset >= 0) {
                        // Calculate flatten index for file_data (discounting offset of 1)
                        size_t flat_idx = ((size_t)(k - 1) * nx * ny + (size_t)(j - 1) * nx + (i - 1)) * n_items;
                        // Get velocity data from flatten
                        double ux   = file_data[flat_idx + 0];
                        double uy   = file_data[flat_idx + 1];
                        double uz   = file_data[flat_idx + 2];
                        double rho  = file_data[flat_idx + 3]*3.0;

                        
                        // Macroscopic momentums from read file
                        double jx   = rho*ux;
                        double jy   = rho*uy;
                        double jz   = rho*uz;

                        // MRT Equilibrium momentums
                        double m1 = (19 * (jx * jx + jy * jy + jz * jz) / rho - 11 * rho);
                        double m2 = (3 * rho - 5.5 * (jx * jx + jy * jy + jz * jz) / rho);
                        //double m_eq3 = //?
                        double m4 = (-0.6666666666666666 * jx);
                        //double m_eq5 = //?
                        double m6 = (-0.6666666666666666 * jy);
                        //double m_eq7 = //?
                        double m8 = (-0.6666666666666666 * jz);
                        double m9 = ((2 * jx * jx - jy * jy - jz * jz) / rho);
                        double m10 = -0.5 * ((2 * jx * jx - jy * jy - jz * jz) / rho);
                        double m11 = ((jy * jy - jz * jz) / rho);
                        double m12 = -0.5 * ((jy * jy - jz * jz) / rho);
                        double m13 = (jx * jy / rho);
                        double m14 = (jy * jz / rho);
                        double m15 = (jx * jz / rho);
                        double m16 = 0.0;
                        double m17 = 0.0;
                        double m18 = 0.0;

                        // MRT Inverse: converting initialized momemtum as distributions
                        // q=0
                        double f_value = 0.0;
                        f_value = mrt_V1 * rho - mrt_V2 * m1 + mrt_V3 * m2;
                        temp_fq[cell_offset] = f_value;

                        // q = 1
                        f_value = mrt_V1 * rho - mrt_V4 * m1 - mrt_V5 * m2 + 0.1 * (jx - m4) +
                            mrt_V6 * (m9 - m10) + 0.16666666 * Fx;
                        temp_fq[1 * Np + cell_offset] = f_value;

                        // q=2
                        f_value = mrt_V1 * rho - mrt_V4 * m1 - mrt_V5 * m2 + 0.1 * (m4 - jx) +
                            mrt_V6 * (m9 - m10) - 0.16666666 * Fx;
                        temp_fq[2 * Np + cell_offset] = f_value;

                        // q = 3
                        f_value = mrt_V1 * rho - mrt_V4 * m1 - mrt_V5 * m2 + 0.1 * (jy - m6) +
                            mrt_V7 * (m10 - m9) + mrt_V8 * (m11 - m12) + 0.16666666 * Fy;
                        temp_fq[3 * Np + cell_offset] = f_value;

                        // q = 4
                        f_value = mrt_V1 * rho - mrt_V4 * m1 - mrt_V5 * m2 + 0.1 * (m6 - jy) +
                            mrt_V7 * (m10 - m9) + mrt_V8 * (m11 - m12) - 0.16666666 * Fy;
                        temp_fq[4 * Np + cell_offset] = f_value;

                        // q = 5
                        f_value = mrt_V1 * rho - mrt_V4 * m1 - mrt_V5 * m2 + 0.1 * (jz - m8) +
                            mrt_V7 * (m10 - m9) + mrt_V8 * (m12 - m11) + 0.16666666 * Fz;
                        temp_fq[5 * Np + cell_offset] = f_value;

                        // q = 6
                        f_value = mrt_V1 * rho - mrt_V4 * m1 - mrt_V5 * m2 + 0.1 * (m8 - jz) +
                            mrt_V7 * (m10 - m9) + mrt_V8 * (m12 - m11) - 0.16666666 * Fz;
                        temp_fq[6 * Np + cell_offset] = f_value;

                        // q = 7
                        f_value = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jx + jy) +
                            0.025 * (m4 + m6) + mrt_V7 * m9 + mrt_V11 * m10 + mrt_V8 * m11 +
                            mrt_V12 * m12 + 0.25 * m13 + 0.125 * (m16 - m17) +
                            0.08333333333 * (Fx + Fy);
                        temp_fq[7 * Np + cell_offset] = f_value;

                        // q = 8
                        f_value = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 - 0.1 * (jx + jy) -
                            0.025 * (m4 + m6) + mrt_V7 * m9 + mrt_V11 * m10 + mrt_V8 * m11 +
                            mrt_V12 * m12 + 0.25 * m13 + 0.125 * (m17 - m16) -
                            0.08333333333 * (Fx + Fy);
                        temp_fq[8 * Np + cell_offset] = f_value;

                        // q = 9
                        f_value = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jx - jy) +
                            0.025 * (m4 - m6) + mrt_V7 * m9 + mrt_V11 * m10 + mrt_V8 * m11 +
                            mrt_V12 * m12 - 0.25 * m13 + 0.125 * (m16 + m17) +
                            0.08333333333 * (Fx - Fy);
                        temp_fq[9 * Np + cell_offset] = f_value;

                        // q = 10
                        f_value = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jy - jx) +
                            0.025 * (m6 - m4) + mrt_V7 * m9 + mrt_V11 * m10 + mrt_V8 * m11 +
                            mrt_V12 * m12 - 0.25 * m13 - 0.125 * (m16 + m17) -
                            0.08333333333 * (Fx - Fy);
                        temp_fq[10 * Np + cell_offset] = f_value;

                        // q = 11
                        f_value = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jx + jz) +
                            0.025 * (m4 + m8) + mrt_V7 * m9 + mrt_V11 * m10 - mrt_V8 * m11 -
                            mrt_V12 * m12 + 0.25 * m15 + 0.125 * (m18 - m16) +
                            0.08333333333 * (Fx + Fz);
                        temp_fq[11 * Np + cell_offset] = f_value;

                        // q = 12
                        f_value = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 - 0.1 * (jx + jz) -
                            0.025 * (m4 + m8) + mrt_V7 * m9 + mrt_V11 * m10 - mrt_V8 * m11 -
                            mrt_V12 * m12 + 0.25 * m15 + 0.125 * (m16 - m18) -
                            0.08333333333 * (Fx + Fz);
                        temp_fq[12 * Np + cell_offset] = f_value;

                        // q = 13
                        f_value = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jx - jz) +
                            0.025 * (m4 - m8) + mrt_V7 * m9 + mrt_V11 * m10 - mrt_V8 * m11 -
                            mrt_V12 * m12 - 0.25 * m15 - 0.125 * (m16 + m18) +
                            0.08333333333 * (Fx - Fz);
                        temp_fq[13 * Np + cell_offset] = f_value;

                        // q= 14
                        f_value = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jz - jx) +
                            0.025 * (m8 - m4) + mrt_V7 * m9 + mrt_V11 * m10 - mrt_V8 * m11 -
                            mrt_V12 * m12 - 0.25 * m15 + 0.125 * (m16 + m18) -
                            0.08333333333 * (Fx - Fz);

                        temp_fq[14 * Np + cell_offset] = f_value;

                        // q = 15
                        f_value = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jy + jz) +
                            0.025 * (m6 + m8) - mrt_V6 * m9 - mrt_V7 * m10 + 0.25 * m14 +
                            0.125 * (m17 - m18) + 0.08333333333 * (Fy + Fz);
                        temp_fq[15 * Np + cell_offset] = f_value;

                        // q = 16
                        f_value = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 - 0.1 * (jy + jz) -
                            0.025 * (m6 + m8) - mrt_V6 * m9 - mrt_V7 * m10 + 0.25 * m14 +
                            0.125 * (m18 - m17) - 0.08333333333 * (Fy + Fz);
                        temp_fq[16 * Np + cell_offset] = f_value;

                        // q = 17
                        f_value = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jy - jz) +
                            0.025 * (m6 - m8) - mrt_V6 * m9 - mrt_V7 * m10 - 0.25 * m14 +
                            0.125 * (m17 + m18) + 0.08333333333 * (Fy - Fz);
                        temp_fq[17 * Np + cell_offset] = f_value;

                        // q = 18
                        f_value = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jz - jy) +
                            0.025 * (m8 - m6) - mrt_V6 * m9 - mrt_V7 * m10 - 0.25 * m14 -
                            0.125 * (m17 + m18) - 0.08333333333 * (Fy - Fz);
                        temp_fq[18 * Np + cell_offset] = f_value;



                    }
                }
            }
        }

        ScaLBL_CopyToDevice(fq, temp_fq, 19 * Np * sizeof(double));
        delete[] temp_fq;
    }
    else {
        if (rank == 0) printf("No start file. Initializing with null velocity case.\n");
    }

    // Update Velocity state from fq
    ScaLBL_D3Q19_Momentum(fq,Velocity,Np);
    ScaLBL_DeviceBarrier();
    comm.barrier();
    ScaLBL_Comm->RegularLayout(Map, &Velocity[0   ], Velocity_x);
    ScaLBL_Comm->RegularLayout(Map, &Velocity[Np  ], Velocity_y);
    ScaLBL_Comm->RegularLayout(Map, &Velocity[2*Np], Velocity_z);

    // Update Pressure State from fq
    ScaLBL_D3Q19_Pressure(fq, Pressure, Np);    // Calculate Pressure Field
    ScaLBL_DeviceBarrier();                     // Sync
    comm.barrier();                             // Sync
    ScaLBL_Comm->RegularLayout(Map, &Pressure[0   ], Pressure_f);  // Transform Pressure Field in 3D domain

}

void ScaLBL_MRTModel::Initialize_fEqNeq() {
    //
	// This function initializes model with equilibrium distributions
    // given by custom velocity and pressure fields.
    // The non-equilibrium moments of: 'e', 'pxx', 'pww', 'pxy', 'pxz' and 'pyz' are initialized
	//
    // Initialize distributions as no velocity Equilibrium
    ScaLBL_D3Q19_Init(fq, Np);


    // If there is a .raw file
    char raw_filename[256];
    sprintf(raw_filename, "Start.%05d.raw", rank);
    std::ifstream binaryFile(raw_filename, std::ios::binary);

    if (binaryFile.good() && Start) {
        // Remove halo extra voxels (added in Domain class)
        unsigned int nx = Nx - 2;
        unsigned int ny = Ny - 2;
        unsigned int nz = Nz - 2;
        unsigned int n_items = 4;

        // Announces the start of the process
        if (rank == 0) printf("Reading start file: %s (%dx%dx%d)\n", raw_filename, nx, ny, nz);


        // Allocate buffer for the domain
        // n_items doubles per voxel (Ux, Uy, Uz, Pr)
        size_t total_voxels = (size_t)nx * ny * nz;
        std::vector<double> file_data(total_voxels * n_items);
        binaryFile.read(reinterpret_cast<char*>(file_data.data()), file_data.size() * sizeof(double));
        binaryFile.close();

        // Allocate auxiliary distributions
        double* temp_fq = new double[19 * Np];
        memset(temp_fq, 0, 19 * Np * sizeof(double)); // Initialize it with zeros

        // MRT constants
        constexpr double mrt_V1 = 0.05263157894736842;
        constexpr double mrt_V2 = 0.012531328320802;
        constexpr double mrt_V3 = 0.04761904761904762;
        constexpr double mrt_V4 = 0.004594820384294068;
        constexpr double mrt_V5 = 0.01587301587301587;
        constexpr double mrt_V6 = 0.0555555555555555555555555;
        constexpr double mrt_V7 = 0.02777777777777778;
        constexpr double mrt_V8 = 0.08333333333333333;
        constexpr double mrt_V9 = 0.003341687552213868;
        constexpr double mrt_V10 = 0.003968253968253968;
        constexpr double mrt_V11 = 0.01388888888888889;
        constexpr double mrt_V12 = 0.04166666666666666;

        // For each cell
        for (unsigned int k = 1; k < nz+1; k++) {
            for (unsigned int j = 1; j < ny+1; j++) {
                for (unsigned int i = 1; i < nx+1; i++) {


                    // Get index in flatten array (solid + fluid + extra cells)
                    int cell_offset = Map(i, j, k);


                    // If is a fluid cell
                    if (cell_offset >= 0) {
                        // Calculate flatten index for file_data (discounting offset of 1)
                        size_t flat_idx = ((size_t)(k - 1) * nx * ny + (size_t)(j - 1) * nx + (i - 1)) * n_items;
                        // Get velocity data from flatten
                        double ux   = file_data[flat_idx + 0];
                        double uy   = file_data[flat_idx + 1];
                        double uz   = file_data[flat_idx + 2];
                        double rho  = file_data[flat_idx + 3]*3.0;

                        // Get derivatives
                        // ADD NON-EQUILIBRIUM
                        // Gradients (default values)
                        double dux_x    = 0.0;
                        double duy_x    = 0.0;
                        double duz_x    = 0.0;
                        double dux_y    = 0.0;
                        double duy_y    = 0.0;
                        double duz_y    = 0.0;
                        double dux_z    = 0.0;
                        double duy_z    = 0.0;
                        double duz_z    = 0.0;
                        // If x direction is free to flow
                        if (Map(i+1, j, k) >= 0 && Map(i-1, j, k) >= 0){
                            size_t id_pox   = ((size_t)(k - 1 +0) * nx * ny + (size_t)(j - 1 + 0) * nx + (i - 1 + 1)) * n_items; // Index in Start.Raw
                            size_t id_prx   = ((size_t)(k - 1 -0) * nx * ny + (size_t)(j - 1 - 0) * nx + (i - 1 - 1)) * n_items; // Index in Start.Raws
                            double ux_pox   = file_data[id_pox + 0]; // Indexes of velocities from next cell
                            double uy_pox   = file_data[id_pox + 1];
                            double uz_pox   = file_data[id_pox + 2];
                            double ux_prx   = file_data[id_prx + 0]; // Indexes of velocities from previous cell
                            double uy_prx   = file_data[id_prx + 1];
                            double uz_prx   = file_data[id_prx + 2];
                            dux_x           = (ux_pox - ux_prx)/(2.0);   // Ux derivate along x
                            duy_x           = (uy_pox - uy_prx)/(2.0);   // Uy derivate along x
                            duz_x           = (uz_pox - uz_prx)/(2.0);   // Uz derivate along x
                        }
                        // If only next x cell is free to flow
                        else if(Map(i+1, j, k) >= 0){
                            size_t id_pox   = ((size_t)(k - 1 +0) * nx * ny + (size_t)(j - 1 + 0) * nx + (i - 1 + 1)) * n_items; // Index in Start.Raw
                            double ux_pox   = file_data[id_pox + 0]; // Indexes of velocities from next cell
                            double uy_pox   = file_data[id_pox + 1];
                            double uz_pox   = file_data[id_pox + 2];
                            dux_x           = (ux_pox - ux);   // Ux derivate along x
                            duy_x           = (uy_pox - uy);   // Uy derivate along x
                            duz_x           = (uz_pox - uz);   // Uz derivate along x
                        }
                        // If only prev. x cell is free to flow
                        else if(Map(i-1, j, k) >= 0){
                            size_t id_prx   = ((size_t)(k - 1 -0) * nx * ny + (size_t)(j - 1 - 0) * nx + (i - 1 - 1)) * n_items; // Index in Start.Raw
                            double ux_prx   = file_data[id_prx + 0]; // Indexes of velocities from previous cell
                            double uy_prx   = file_data[id_prx + 1];
                            double uz_prx   = file_data[id_prx + 2];
                            dux_x           = (ux - ux_prx);   // Ux derivate along x
                            duy_x           = (uy - uy_prx);   // Uy derivate along x
                            duz_x           = (uz - uz_prx);   // Uz derivate along x

                        }

                        // If y direction is free to flow
                        if (Map(i, j+1, k) >= 0 && Map(i, j-1, k) >= 0){
                            size_t id_poy   = ((size_t)(k - 1 + 0) * nx * ny + (size_t)(j - 1 + 1) * nx + (i - 1 + 0)) * n_items; // Index in Start.Raw
                            size_t id_pry   = ((size_t)(k - 1 - 0) * nx * ny + (size_t)(j - 1 - 1) * nx + (i - 1 - 0)) * n_items; // Index in Start.Raw
                            double ux_poy   = file_data[id_poy + 0]; // Indexes of velocities from next cell
                            double uy_poy   = file_data[id_poy + 1];
                            double uz_poy   = file_data[id_poy + 2];
                            double ux_pry   = file_data[id_pry + 0]; // Indexes of velocities from previous cell
                            double uy_pry   = file_data[id_pry + 1];
                            double uz_pry   = file_data[id_pry + 2];
                            dux_y           = (ux_poy - ux_pry)/2.0;   // Ux derivate along y
                            duy_y           = (uy_poy - uy_pry)/2.0;   // Uy derivate along y
                            duz_y           = (uz_poy - uz_pry)/2.0;   // Uz derivate along y
                        }else if(Map(i, j+1, k) >= 0){
                            size_t id_poy   = ((size_t)(k - 1 + 0) * nx * ny + (size_t)(j - 1 + 1) * nx + (i - 1 + 0)) * n_items; // Index in Start.Raw
                            double ux_poy   = file_data[id_poy + 0]; // Indexes of velocities from next cell
                            double uy_poy   = file_data[id_poy + 1];
                            double uz_poy   = file_data[id_poy + 2];
                            dux_y           = (ux_poy - ux);   // Ux derivate along y
                            duy_y           = (uy_poy - uy);   // Uy derivate along y
                            duz_y           = (uz_poy - uz);   // Uz derivate along y
                        }else if(Map(i, j-1, k) >= 0){
                            size_t id_pry   = ((size_t)(k - 1 - 0) * nx * ny + (size_t)(j - 1 - 1) * nx + (i - 1 - 0)) * n_items; // Index in Start.Raw
                            double ux_pry   = file_data[id_pry + 0]; // Indexes of velocities from previous cell
                            double uy_pry   = file_data[id_pry + 1];
                            double uz_pry   = file_data[id_pry + 2];
                            dux_y           = (ux - ux_pry);   // Ux derivate along y
                            duy_y           = (uy - uy_pry);   // Uy derivate along y
                            duz_y           = (uz - uz_pry);   // Uz derivate along y
                        }

                        // If z direction is free to flow
                        if (Map(i, j, k+1) >= 0 && Map(i, j, k-1) >= 0){
                            size_t id_poz   = ((size_t)(k - 1 + 1) * nx * ny + (size_t)(j - 1 + 0) * nx + (i - 1 + 0)) * n_items; // Index in Start.Raw
                            size_t id_prz   = ((size_t)(k - 1 - 1) * nx * ny + (size_t)(j - 1 - 0) * nx + (i - 1 - 0)) * n_items; // Index in Start.Raw
                            double ux_poz   = file_data[id_poz + 0]; // Indexes of velocities from next cell
                            double uy_poz   = file_data[id_poz + 1];
                            double uz_poz   = file_data[id_poz + 2];
                            double ux_prz   = file_data[id_prz + 0]; // Indexes of velocities from previous cell
                            double uy_prz   = file_data[id_prz + 1];
                            double uz_prz   = file_data[id_prz + 2];
                            dux_z           = (ux_poz - ux_prz)/(2.0);   // Ux derivate along z
                            duy_z           = (uy_poz - uy_prz)/(2.0);   // Uy derivate along z
                            duz_z           = (uz_poz - uz_prz)/(2.0);   // Uz derivate along z
                        } else if(Map(i, j, k+1) >= 0){
                            size_t id_poz   = ((size_t)(k - 1 + 1) * nx * ny + (size_t)(j - 1 + 0) * nx + (i - 1 + 0)) * n_items; // Index in Start.Raw
                            double ux_poz   = file_data[id_poz + 0]; // Indexes of velocities from next cell
                            double uy_poz   = file_data[id_poz + 1];
                            double uz_poz   = file_data[id_poz + 2];
                            dux_z           = (ux_poz - ux);   // Ux derivate along z
                            duy_z           = (uy_poz - uy);   // Uy derivate along z
                            duz_z           = (uz_poz - uz);   // Uz derivate along z

                        } else if(Map(i, j, k-1) >= 0){
                            size_t id_prz   = ((size_t)(k - 1 - 1) * nx * ny + (size_t)(j - 1 - 0) * nx + (i - 1 - 0)) * n_items; // Index in Start.Raw
                            double ux_prz   = file_data[id_prz + 0]; // Indexes of velocities from previous cell
                            double uy_prz   = file_data[id_prz + 1];
                            double uz_prz   = file_data[id_prz + 2];
                            dux_z           = (ux - ux_prz);   // Ux derivate along z
                            duy_z           = (uy - uy_prz);   // Uy derivate along z
                            duz_z           = (uz - uz_prz);   // Uz derivate along z
                        }
                        // Macroscopic momentums from read file
                        double jx   = rho*ux;
                        double jy   = rho*uy;
                        double jz   = rho*uz;
                        double divergent = (dux_x+duy_y+duz_z);

                        // MRT Equilibrium momentums
                        double m_eq1 = (19 * (jx * jx + jy * jy + jz * jz) / rho - 11 * rho);
                        double m_eq2 = (3 * rho - 5.5 * (jx * jx + jy * jy + jz * jz) / rho);
                        double m_eq4 = (-0.6666666666666666 * jx);
                        double m_eq6 = (-0.6666666666666666 * jy);
                        double m_eq8 = (-0.6666666666666666 * jz);
                        double m_eq9 = ((2 * jx * jx - jy * jy - jz * jz) / rho);
                        double m_eq10 = -0.5 * ((2 * jx * jx - jy * jy - jz * jz) / rho);
                        double m_eq11 = ((jy * jy - jz * jz) / rho);
                        double m_eq12 = -0.5 * ((jy * jy - jz * jz) / rho);
                        double m_eq13 = (jx * jy / rho);
                        double m_eq14 = (jy * jz / rho);
                        double m_eq15 = (jx * jz / rho);
                        double m_eq16 = 0.0;
                        double m_eq17 = 0.0;
                        double m_eq18 = 0.0;

                        // MRT Non-equilibrium momentums
                        // Relaxations according to collisor
                        double relax_e = rlx_setA; //m1
                        double relax_p = rlx_setA; //m9, m10, m13, m14, m15 
                        double relax_q = rlx_setB; // m16, m17, m18
                        // e
                        double m_neq1 = - 19* divergent /  relax_e; 
                        m_neq1 *= (1-relax_e); // Convert to post-collision

                        double m_neq2 = 0.0;   // Epsilon
                        double m_neq4 = 0.0;   // q_x
                        double m_neq6 = 0.0;   // q_y
                        double m_neq8 = 0.0;   // q_z
                        
                        // 3p_xx
                        double m_neq9 = - 2.0 * rho * (2*dux_x-duy_y-duz_z) / (3.0*relax_p);
                        m_neq9 *= (1-relax_p); // Convert to post-collision
                        double m_neq10 =  - 0.5 * m_neq9;  // pi_xx

                        // p_ww
                        double m_neq11 = - 2.0 * rho * (duy_y - duz_z)/ (3.0*relax_p);
                        m_neq11 *= (1-relax_p); // Convert to post-collision
                        double m_neq12 = -0.5 * m_neq11;  // pi_ww

                        // p_xy
                        double m_neq13 = - rho *(dux_y+duy_x)/ (3.0*relax_p);
                        m_neq13 *= (1-relax_p); // Convert to post-collision

                        // p_yz
                        double m_neq14 = - rho *(duy_z+duz_y)/ (3.0*relax_p);
                        m_neq14 *= (1-relax_p); // Convert to post-collision

                        // p_xz
                        double m_neq15 = - rho *(dux_z+duz_x)/ (3.0*relax_p);
                        m_neq15 *= (1-relax_p); // Convert to post-collision

                        double m_neq16 = 0.0*relax_q;  // m_x
                        m_neq16 *= (1-relax_q);
                        
                        double m_neq17 = 0.0*relax_q;  // m_y
                        m_neq17 *= (1-relax_q);

                        double m_neq18 = 0.0*relax_q;  // m_z
                        m_neq18 *= (1-relax_q);


                        // Define total of initialized moments 
                        double m1 = m_eq1 + m_neq1;
                        double m2 = m_eq2 + m_neq2;
                        double m4 = m_eq4 + m_neq4;
                        double m6 = m_eq6 + m_neq6;
                        double m8 = m_eq8 + m_neq8;
                        double m9 = m_eq9 + m_neq9;
                        double m10 = m_eq10 + m_neq10;
                        double m11 = m_eq11 + m_neq11;
                        double m12 = m_eq12 + m_neq12;
                        double m13 = m_eq13 + m_neq13;
                        double m14 = m_eq14 + m_neq14;
                        double m15 = m_eq15 + m_neq15;
                        double m16 = m_eq16 + m_neq16;
                        double m17 = m_eq17 + m_neq17;
                        double m18 = m_eq18 + m_neq18;


                        // MRT Inverse: converting initialized momemtum as distributions
                        // q=0
                        double f_value = 0.0;
                        f_value = mrt_V1 * rho - mrt_V2 * m1 + mrt_V3 * m2;
                        temp_fq[cell_offset] = f_value;

                        // q = 1
                        f_value = mrt_V1 * rho - mrt_V4 * m1 - mrt_V5 * m2 + 0.1 * (jx - m4) +
                            mrt_V6 * (m9 - m10) + 0.16666666 * Fx;
                        temp_fq[1 * Np + cell_offset] = f_value;

                        // q=2
                        f_value = mrt_V1 * rho - mrt_V4 * m1 - mrt_V5 * m2 + 0.1 * (m4 - jx) +
                            mrt_V6 * (m9 - m10) - 0.16666666 * Fx;
                        temp_fq[2 * Np + cell_offset] = f_value;

                        // q = 3
                        f_value = mrt_V1 * rho - mrt_V4 * m1 - mrt_V5 * m2 + 0.1 * (jy - m6) +
                            mrt_V7 * (m10 - m9) + mrt_V8 * (m11 - m12) + 0.16666666 * Fy;
                        temp_fq[3 * Np + cell_offset] = f_value;

                        // q = 4
                        f_value = mrt_V1 * rho - mrt_V4 * m1 - mrt_V5 * m2 + 0.1 * (m6 - jy) +
                            mrt_V7 * (m10 - m9) + mrt_V8 * (m11 - m12) - 0.16666666 * Fy;
                        temp_fq[4 * Np + cell_offset] = f_value;

                        // q = 5
                        f_value = mrt_V1 * rho - mrt_V4 * m1 - mrt_V5 * m2 + 0.1 * (jz - m8) +
                            mrt_V7 * (m10 - m9) + mrt_V8 * (m12 - m11) + 0.16666666 * Fz;
                        temp_fq[5 * Np + cell_offset] = f_value;

                        // q = 6
                        f_value = mrt_V1 * rho - mrt_V4 * m1 - mrt_V5 * m2 + 0.1 * (m8 - jz) +
                            mrt_V7 * (m10 - m9) + mrt_V8 * (m12 - m11) - 0.16666666 * Fz;
                        temp_fq[6 * Np + cell_offset] = f_value;

                        // q = 7
                        f_value = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jx + jy) +
                            0.025 * (m4 + m6) + mrt_V7 * m9 + mrt_V11 * m10 + mrt_V8 * m11 +
                            mrt_V12 * m12 + 0.25 * m13 + 0.125 * (m16 - m17) +
                            0.08333333333 * (Fx + Fy);
                        temp_fq[7 * Np + cell_offset] = f_value;

                        // q = 8
                        f_value = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 - 0.1 * (jx + jy) -
                            0.025 * (m4 + m6) + mrt_V7 * m9 + mrt_V11 * m10 + mrt_V8 * m11 +
                            mrt_V12 * m12 + 0.25 * m13 + 0.125 * (m17 - m16) -
                            0.08333333333 * (Fx + Fy);
                        temp_fq[8 * Np + cell_offset] = f_value;

                        // q = 9
                        f_value = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jx - jy) +
                            0.025 * (m4 - m6) + mrt_V7 * m9 + mrt_V11 * m10 + mrt_V8 * m11 +
                            mrt_V12 * m12 - 0.25 * m13 + 0.125 * (m16 + m17) +
                            0.08333333333 * (Fx - Fy);
                        temp_fq[9 * Np + cell_offset] = f_value;

                        // q = 10
                        f_value = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jy - jx) +
                            0.025 * (m6 - m4) + mrt_V7 * m9 + mrt_V11 * m10 + mrt_V8 * m11 +
                            mrt_V12 * m12 - 0.25 * m13 - 0.125 * (m16 + m17) -
                            0.08333333333 * (Fx - Fy);
                        temp_fq[10 * Np + cell_offset] = f_value;

                        // q = 11
                        f_value = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jx + jz) +
                            0.025 * (m4 + m8) + mrt_V7 * m9 + mrt_V11 * m10 - mrt_V8 * m11 -
                            mrt_V12 * m12 + 0.25 * m15 + 0.125 * (m18 - m16) +
                            0.08333333333 * (Fx + Fz);
                        temp_fq[11 * Np + cell_offset] = f_value;

                        // q = 12
                        f_value = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 - 0.1 * (jx + jz) -
                            0.025 * (m4 + m8) + mrt_V7 * m9 + mrt_V11 * m10 - mrt_V8 * m11 -
                            mrt_V12 * m12 + 0.25 * m15 + 0.125 * (m16 - m18) -
                            0.08333333333 * (Fx + Fz);
                        temp_fq[12 * Np + cell_offset] = f_value;

                        // q = 13
                        f_value = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jx - jz) +
                            0.025 * (m4 - m8) + mrt_V7 * m9 + mrt_V11 * m10 - mrt_V8 * m11 -
                            mrt_V12 * m12 - 0.25 * m15 - 0.125 * (m16 + m18) +
                            0.08333333333 * (Fx - Fz);
                        temp_fq[13 * Np + cell_offset] = f_value;

                        // q= 14
                        f_value = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jz - jx) +
                            0.025 * (m8 - m4) + mrt_V7 * m9 + mrt_V11 * m10 - mrt_V8 * m11 -
                            mrt_V12 * m12 - 0.25 * m15 + 0.125 * (m16 + m18) -
                            0.08333333333 * (Fx - Fz);

                        temp_fq[14 * Np + cell_offset] = f_value;

                        // q = 15
                        f_value = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jy + jz) +
                            0.025 * (m6 + m8) - mrt_V6 * m9 - mrt_V7 * m10 + 0.25 * m14 +
                            0.125 * (m17 - m18) + 0.08333333333 * (Fy + Fz);
                        temp_fq[15 * Np + cell_offset] = f_value;

                        // q = 16
                        f_value = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 - 0.1 * (jy + jz) -
                            0.025 * (m6 + m8) - mrt_V6 * m9 - mrt_V7 * m10 + 0.25 * m14 +
                            0.125 * (m18 - m17) - 0.08333333333 * (Fy + Fz);
                        temp_fq[16 * Np + cell_offset] = f_value;

                        // q = 17
                        f_value = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jy - jz) +
                            0.025 * (m6 - m8) - mrt_V6 * m9 - mrt_V7 * m10 - 0.25 * m14 +
                            0.125 * (m17 + m18) + 0.08333333333 * (Fy - Fz);
                        temp_fq[17 * Np + cell_offset] = f_value;

                        // q = 18
                        f_value = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jz - jy) +
                            0.025 * (m8 - m6) - mrt_V6 * m9 - mrt_V7 * m10 - 0.25 * m14 -
                            0.125 * (m17 + m18) - 0.08333333333 * (Fy - Fz);
                        temp_fq[18 * Np + cell_offset] = f_value;



                    }
                }
            }
        }

        ScaLBL_CopyToDevice(fq, temp_fq, 19 * Np * sizeof(double));
        delete[] temp_fq;
    }
    else {
        if (rank == 0) printf("No start file. Initializing with null velocity case.\n");
    }

    // Update Velocity state from fq
    ScaLBL_D3Q19_Momentum(fq,Velocity,Np);
    ScaLBL_DeviceBarrier();
    comm.barrier();
    ScaLBL_Comm->RegularLayout(Map, &Velocity[0   ], Velocity_x);
    ScaLBL_Comm->RegularLayout(Map, &Velocity[Np  ], Velocity_y);
    ScaLBL_Comm->RegularLayout(Map, &Velocity[2*Np], Velocity_z);

    // Update Pressure State from fq
    ScaLBL_D3Q19_Pressure(fq, Pressure, Np);    // Calculate Pressure Field
    ScaLBL_DeviceBarrier();                     // Sync
    comm.barrier();                             // Sync
    ScaLBL_Comm->RegularLayout(Map, &Pressure[0   ], Pressure_f);  // Transform Pressure Field in 3D domain

}



//This script was made for validation only, and inteded to be removed afterwards
void ScaLBL_MRTModel::Run_Timesteps(const std::vector<int>& coords) {

    // Warn about viability of positions
    for (size_t n = 0; n < coords.size(); n += 3){
        int i = coords[n]+1;
        int j = coords[n+1]+1;
        int k = coords[n+2]+1;
        if (Distance(i, j, k) <= 0){
            printf("Position (x=%d,y=%d,z=%d) is not valid fluid. It will not be considered.\n", coords[n], coords[n+1], coords[n+2]);
        }
    }

    double rlx_setA = 1.0 / tau;
    double rlx_setB = 8.f * (2.f - rlx_setA) / (8.f - rlx_setA);

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
            fprintf(log_file, "time Fx Fy Fz mu Vs As Js Xs vx vy vz absperm\n");
            fclose(log_file);
        }
    }

    //.......create and start timer............
    ScaLBL_DeviceBarrier();
    comm.barrier();
    if (rank == 0)
        printf("Beginning AA timesteps, timestepMax = %i \n", timestepMax);
    if (rank == 0)
        printf("********************************************************\n");
    timestep = 0;
    auto t1 = std::chrono::system_clock::now();


    double count_loc = 0;
    double kin_energy_init = 0;
    for (int k = 1; k < Nz - 1; k++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                if (Distance(i, j, k) > 0) {
                    kin_energy_init += 2*std::pow(Velocity_x(i, j, k),2) + std::pow(Velocity_y(i, j, k),2) + std::pow(Velocity_z(i, j, k),2);
                    count_loc  += 1.0;
                }
            }
        }
    }
    kin_energy_init /= count_loc;


    // --- 1. DYNAMIC HEADER GENERATION ---
    if (rank == 0) {
        // 10 (Time) + 20 (Global) = 30 spaces padding
        printf("\n%-10s %-20s", "", "");
        for (size_t n = 0; n < coords.size(); n += 3) {
            int i = coords[n]+1, j = coords[n+1]+1, k = coords[n+2]+1;
            if (Distance(i, j, k) > 0) {
                char label[32];
                sprintf(label, "Point (%d,%d,%d)", i, j, k);
                // MATH: '| ' (2) + label + remaining spaces to hit 33
                // The %-31s ensures the point block is exactly 33 chars wide (including the |)
                printf("| %-31s", label);
            }
        }

        // --- LEVEL 2: VARIABLE LABELS ---
        printf("\n%-10s %-20s", "Timestep", "Global Kin. Energy");
        for (size_t n = 0; n < coords.size(); n += 3) {
            if (Distance(coords[n]+1, coords[n+1]+1, coords[n+2]+1) > 0) {
                // MATH: '| ' (2) + 14 (Kin) + ' | ' (3) + 14 (Press) = 33 characters
                printf("| %-14s | %-14s", "Kin. Energy", "Pressure");
            }
        }
        printf("\n----------------------------------------------------------------------------------------------------------\n");
    }

    while (timestep < timestepMax) {
        //************************************************************************/
        timestep++;
        ScaLBL_Comm->SendD3Q19AA(fq); //READ FROM NORMAL
        ScaLBL_D3Q19_AAodd_MRT(NeighborList, fq, ScaLBL_Comm->FirstInterior(),
                               ScaLBL_Comm->LastInterior(), Np, rlx_setA,
                               rlx_setB, Fx, Fy, Fz);
        ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
        // Set boundary conditions
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

        ScaLBL_D3Q19_AAodd_MRT(NeighborList, fq, 0, ScaLBL_Comm->LastExterior(),
                               Np, rlx_setA, rlx_setB, Fx, Fy, Fz);
        ScaLBL_DeviceBarrier();
        comm.barrier();
        //************************************************************************/
        timestep++;
        ScaLBL_Comm->SendD3Q19AA(fq); //READ FORM NORMAL
        ScaLBL_D3Q19_AAeven_MRT(fq, ScaLBL_Comm->FirstInterior(),
                                ScaLBL_Comm->LastInterior(), Np, rlx_setA,
                                rlx_setB, Fx, Fy, Fz);
        ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
        // Set boundary conditions
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
        ScaLBL_D3Q19_AAeven_MRT(fq, 0, ScaLBL_Comm->LastExterior(), Np,
                                rlx_setA, rlx_setB, Fx, Fy, Fz);
        ScaLBL_DeviceBarrier();
        comm.barrier();
        //************************************************************************/

        if (timestep % ANALYSIS_INTERVAL == 0) {
            ScaLBL_D3Q19_Momentum(fq, Velocity, Np);
            ScaLBL_DeviceBarrier();
            comm.barrier();
            ScaLBL_Comm->RegularLayout(Map, &Velocity[0   ], Velocity_x);
            ScaLBL_Comm->RegularLayout(Map, &Velocity[Np  ], Velocity_y);
            ScaLBL_Comm->RegularLayout(Map, &Velocity[2*Np], Velocity_z);

            ScaLBL_D3Q19_Pressure(fq, Pressure, Np);
            ScaLBL_DeviceBarrier();                     // Sync
            comm.barrier();                             // Sync
            ScaLBL_Comm->RegularLayout(Map, &Pressure[0   ], Pressure_f);  // Transform Pressure Field in 3D domain

            // Calculate global kinectic energy
            double kin_energy   = 0;
            for (int k = 1; k < Nz - 1; k++) {
                for (int j = 1; j < Ny - 1; j++) {
                    for (int i = 1; i < Nx - 1; i++) {
                        if (Distance(i, j, k) > 0) {
                            kin_energy    += std::pow(Velocity_x(i, j, k),2) + std::pow(Velocity_y(i, j, k),2) + std::pow(Velocity_z(i, j, k),2);
                        }
                    }
                }
            }
            kin_energy /= count_loc;

            // Checking analyzed points
            if (rank == 0) {
                printf("%-10d %-20.6e", timestep, kin_energy);
                for (size_t n = 0; n < coords.size(); n += 3) {
                    int i = coords[n]+1, j = coords[n+1]+1, k = coords[n+2]+1;
                    if (Distance(i, j, k) > 0) {
                        double u_loc = std::pow(Velocity_x(i, j, k), 2) + std::pow(Velocity_y(i, j, k), 2) + std::pow(Velocity_z(i, j, k), 2);
                        double p_loc = Pressure_f(i, j, k);
                        // Matches the 33-char header block exactly
                        printf("| %-14.6e | %-14.6e", u_loc, p_loc);
                    }
                }
                printf("\n");
            }

            if (timestep % VISUAL_INTERVAL == 0) SaveFields();

        }
    }
        
    printf("------------------------------------------------------------\n");
    //************************************************************************/
    if (rank == 0)
        printf("--------------------------------------------------------\n");
    // Compute the walltime per timestep
    auto t2 = std::chrono::system_clock::now();
    double cputime = std::chrono::duration<double>(t2 - t1).count() / timestep;
    // Performance obtained from each node
    double MLUPS = double(Np) / cputime / 1000000;

    if (rank == 0)
        printf("********************************************************\n");
    if (rank == 0)
        printf("CPU time = %f \n", cputime);
    if (rank == 0)
        printf("Lattice update rate (per core)= %f MLUPS \n", MLUPS);
    MLUPS *= nprocs;
    if (rank == 0)
        printf("Lattice update rate (total)= %f MLUPS \n", MLUPS);
    if (rank == 0)
        printf("********************************************************\n");
}
//Stop removing here



void ScaLBL_MRTModel::Run() {
    

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
            fprintf(log_file, "time Fx Fy Fz mu Vs As Js Xs vx vy vz absperm(mDa) absperm*(mDa)\n");
            fclose(log_file);
        }
    }

    //.......create and start timer............
    ScaLBL_DeviceBarrier();
    comm.barrier();
    if (rank == 0)
        printf("Beginning AA timesteps, timestepMax = %i \n", timestepMax);
    if (rank == 0)
        printf("********************************************************\n");
    timestep = 0;
    double error = 1.0;
    double flow_rate_previous = 0.0;
    auto t1 = std::chrono::system_clock::now();
    while (timestep < timestepMax && error > tolerance) {
        //************************************************************************/
        timestep++;
        ScaLBL_Comm->SendD3Q19AA(fq); //READ FROM NORMAL
        ScaLBL_D3Q19_AAodd_MRT(NeighborList, fq, ScaLBL_Comm->FirstInterior(),
                               ScaLBL_Comm->LastInterior(), Np, rlx_setA,
                               rlx_setB, Fx, Fy, Fz);
        ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
        // Set boundary conditions
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
        } else if (BoundaryCondition == 6) {
            ScaLBL_Comm->D3Q19_PeriodicPressure_BC_z(NeighborList, fq, dp,
                                                     timestep);
            ScaLBL_Comm->D3Q19_PeriodicPressure_BC_Z(NeighborList, fq, dp,
                                                     timestep);
        } else if (BoundaryCondition == 7) {
            din =
                ScaLBL_Comm->D3Q19_FluxCalculate_BC_z(NeighborList, fq, flux, timestep);
            dout =
                ScaLBL_Comm->D3Q19_FluxCalculate_BC_Z(NeighborList, fq, flux, timestep);
            dp = (din - dout)/3.0f;
            ScaLBL_Comm->D3Q19_PeriodicPressure_BC_z(NeighborList, fq, dp,
                                                     timestep);
            ScaLBL_Comm->D3Q19_PeriodicPressure_BC_Z(NeighborList, fq, dp,
                                                     timestep);
        }
        ScaLBL_D3Q19_AAodd_MRT(NeighborList, fq, 0, ScaLBL_Comm->LastExterior(),
                               Np, rlx_setA, rlx_setB, Fx, Fy, Fz);
        ScaLBL_DeviceBarrier();
        comm.barrier();
        timestep++;
        ScaLBL_Comm->SendD3Q19AA(fq); //READ FORM NORMAL
        ScaLBL_D3Q19_AAeven_MRT(fq, ScaLBL_Comm->FirstInterior(),
                                ScaLBL_Comm->LastInterior(), Np, rlx_setA,
                                rlx_setB, Fx, Fy, Fz);
        ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
        // Set boundary conditions
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
        } else if (BoundaryCondition == 6) {
            ScaLBL_Comm->D3Q19_PeriodicPressure_BC_z(NeighborList, fq, dp,
                                                     timestep);
            ScaLBL_Comm->D3Q19_PeriodicPressure_BC_Z(NeighborList, fq, dp,
                                                     timestep);
        } else if (BoundaryCondition == 7) {
            din =
                ScaLBL_Comm->D3Q19_FluxCalculate_BC_z(NeighborList, fq, flux, timestep);
            dout =
                ScaLBL_Comm->D3Q19_FluxCalculate_BC_Z(NeighborList, fq, flux, timestep);
            dp = (din - dout)/3.0f;
            ScaLBL_Comm->D3Q19_PeriodicPressure_BC_z(NeighborList, fq, dp,
                                                     timestep);
            ScaLBL_Comm->D3Q19_PeriodicPressure_BC_Z(NeighborList, fq, dp,
                                                     timestep);
        }
        ScaLBL_D3Q19_AAeven_MRT(fq, 0, ScaLBL_Comm->LastExterior(), Np,
                                rlx_setA, rlx_setB, Fx, Fy, Fz);
        ScaLBL_DeviceBarrier();
        comm.barrier();
        //************************************************************************/

        if (timestep % ANALYSIS_INTERVAL == 0) {
            ScaLBL_D3Q19_Momentum_2nd_order(fq, Velocity, Np, Fx, Fy, Fz);
            ScaLBL_D3Q19_Pressure(fq, Pressure, Np);
            ScaLBL_DeviceBarrier();
            comm.barrier();
            ScaLBL_Comm->RegularLayout(Map, &Velocity[0], Velocity_x);
            ScaLBL_Comm->RegularLayout(Map, &Velocity[Np], Velocity_y);
            ScaLBL_Comm->RegularLayout(Map, &Velocity[2 * Np], Velocity_z);
            ScaLBL_Comm->RegularLayout(Map, &Pressure[0], Pressure_f);

            double count_loc = 0;
            double count;
            double vax, vay, vaz;
            double vax_loc, vay_loc, vaz_loc;
            vax_loc = vay_loc = vaz_loc = 0.f;
            for (int k = 1; k < Nz - 1; k++) {
                for (int j = 1; j < Ny - 1; j++) {
                    for (int i = 1; i < Nx - 1; i++) {
                        if (Distance(i, j, k) > 0) {
                            vax_loc += Velocity_x(i, j, k);
                            vay_loc += Velocity_y(i, j, k);
                            vaz_loc += Velocity_z(i, j, k);
                            count_loc += 1.0;
                        }
                    }
                }
            }
            vax = Dm->Comm.sumReduce(vax_loc);
            vay = Dm->Comm.sumReduce(vay_loc);
            vaz = Dm->Comm.sumReduce(vaz_loc);
            count = Dm->Comm.sumReduce(count_loc);

            vax /= count;
            vay /= count;
            vaz /= count;

            double force_mag = sqrt(Fx * Fx + Fy * Fy + Fz * Fz);
            double dir_x = Fx / force_mag;
            double dir_y = Fy / force_mag;
            double dir_z = Fz / force_mag;
            if (force_mag == 0.0) {
                // default to z direction
                dir_x = 0.0;
                dir_y = 0.0;
                dir_z = 1.0;
                force_mag = 1.0;
            }
            double flow_rate = (vax * dir_x + vay * dir_y + vaz * dir_z);

            error = fabs(flow_rate - flow_rate_previous) / fabs(flow_rate);
            flow_rate_previous = flow_rate;

            //if (rank==0) printf("Computing Minkowski functionals \n");
            Morphology.ComputeScalar(Distance, 0.f);
            //Morphology.PrintAll();
            double mu = (tau - 0.5) / 3.f;
            double Vs = Morphology.V();
            double As = Morphology.A();
            double Hs = Morphology.H();
            double Xs = Morphology.X();
            Vs = Dm->Comm.sumReduce(Vs);
            As = Dm->Comm.sumReduce(As);
            Hs = Dm->Comm.sumReduce(Hs);
            Xs = Dm->Comm.sumReduce(Xs);

            double h = Dm->voxel_length;


            // Calculate Permeability
            double absperm = 0.0;
            if (BoundaryCondition == 3){
                absperm = h * h * mu * Mask->Porosity() * flow_rate / (Fz*(dout+din)/2   - (dout-din)/((Nz-2)*nprocz*3.0));
            }
            else if (BoundaryCondition == 6 || BoundaryCondition == 7){
                absperm = h * h * mu * Mask->Porosity() * flow_rate / (dp/((Nz-2)*nprocz));
            }
            else{
                absperm = h * h * mu * Mask->Porosity() * flow_rate / (force_mag);
            }
            absperm *= 1013.0; // Convert to mDarcy

            if (timestep % VISUAL_INTERVAL == 0) SaveFields();

            if (rank == 0) {
                printf("     %f\n", absperm);
                FILE *log_file = fopen("Permeability.csv", "a");
                fprintf(log_file,
                        "%i %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g "
                        "%.8g %.8g %.8g\n",
                        timestep, Fx, Fy, Fz, mu, h * h * h * Vs, h * h * As,
                        h * Hs, Xs, vax, vay, vaz, absperm, absperm * Mask->Porosity());
                fclose(log_file);
            }
        }
    }
    //************************************************************************/
    if (rank == 0)
        printf("--------------------------------------------------------\n");
    // Compute the walltime per timestep
    auto t2 = std::chrono::system_clock::now();
    double cputime = std::chrono::duration<double>(t2 - t1).count() / timestep;
    // Performance obtained from each node
    double MLUPS = double(Np) / cputime / 1000000;

    if (rank == 0)
        printf("********************************************************\n");
    if (rank == 0)
        printf("CPU time = %f \n", cputime);
    if (rank == 0)
        printf("Lattice update rate (per core)= %f MLUPS \n", MLUPS);
    MLUPS *= nprocs;
    if (rank == 0)
        printf("Lattice update rate (total)= %f MLUPS \n", MLUPS);
    if (rank == 0)
        printf("********************************************************\n");
}

void ScaLBL_MRTModel::SaveFields() {
    // This function saves the current content of Velocity and Pressure
    // The content must be previously computed by ScaLBL_D3Q19_Momentum() and ScaLBL_D3Q19_Pressure()


    // Define output format
    auto format = vis_db->getWithDefault<string>("format", "silo");

    if (vis_db->getWithDefault<bool>("write_silo", false)) {
        // Create Mesh
        std::vector<IO::MeshDataStruct> visData;
        fillHalo<double> fillData(  Dm->Comm, Dm->rank_info,
                                    {Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2},
                                    {1, 1, 1}, 0, 1);
        auto SignDistVar    = std::make_shared<IO::Variable>();
        auto VxVar          = std::make_shared<IO::Variable>();
        auto VyVar          = std::make_shared<IO::Variable>();
        auto VzVar          = std::make_shared<IO::Variable>();
        auto Press          = std::make_shared<IO::Variable>();

        IO::initialize("", format, false);

        // Create the MeshDataStruct
        visData.resize(1);
        visData[0].meshName = "domain";
        visData[0].mesh = std::make_shared<IO::DomainMesh>(
            Dm->rank_info, Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2, Dm->Lx, Dm->Ly,
            Dm->Lz);

        // SAVE VARIABLES
        unsigned int c_vars = 0;

        // Save Distance Transform
        SignDistVar->name = "SignDist";
        SignDistVar->type = IO::VariableType::VolumeVariable;
        SignDistVar->dim = 1;
        SignDistVar->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
        visData[0].vars.push_back(SignDistVar);
        Array<double> &SignData = visData[0].vars[c_vars]->data;
        ASSERT(visData[0].vars[c_vars]->name == "SignDist");
        fillData.copy(Distance, SignData);
        c_vars++;

        if(save_velocity){
            VxVar->name = "Velocity_x";
            VxVar->type = IO::VariableType::VolumeVariable;
            VxVar->dim = 1;
            VxVar->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
            visData[0].vars.push_back(VxVar);
            Array<double> &VelxData = visData[0].vars[c_vars]->data;
            ASSERT(visData[0].vars[c_vars]->name == "Velocity_x");
            fillData.copy(Velocity_x, VelxData);
            c_vars++;

            VyVar->name = "Velocity_y";
            VyVar->type = IO::VariableType::VolumeVariable;
            VyVar->dim = 1;
            VyVar->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
            visData[0].vars.push_back(VyVar);
            Array<double> &VelyData = visData[0].vars[c_vars]->data;
            ASSERT(visData[0].vars[c_vars]->name == "Velocity_y");
            fillData.copy(Velocity_y, VelyData);
            c_vars++;

            VzVar->name = "Velocity_z";
            VzVar->type = IO::VariableType::VolumeVariable;
            VzVar->dim = 1;
            VzVar->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
            visData[0].vars.push_back(VzVar);
            Array<double> &VelzData = visData[0].vars[c_vars]->data;
            ASSERT(visData[0].vars[c_vars]->name == "Velocity_z");
            fillData.copy(Velocity_z, VelzData);
            c_vars++;
        }

        if (save_pressure){
            Press->name = "Pressure";
            Press->type = IO::VariableType::VolumeVariable;
            Press->dim = 1;
            Press->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
            visData[0].vars.push_back(Press);
            Array<double> &PressData = visData[0].vars[c_vars]->data;
            ASSERT(visData[0].vars[c_vars]->name == "Pressure");
            fillData.copy(Pressure_f,  PressData);
            c_vars++;
        }

        IO::writeData(timestep, visData, Dm->Comm);
    }

}

void ScaLBL_MRTModel::VelocityField() {

    auto format = vis_db->getWithDefault<string>("format", "silo");

    /*	memcpy(Morphology.SDn.data(), Distance.data(), Nx*Ny*Nz*sizeof(double));
	Morphology.Initialize();
	Morphology.UpdateMeshValues();
	Morphology.ComputeLocal();
	Morphology.Reduce();
	
	double count_loc=0;
	double count;
	double vax,vay,vaz;
	double vax_loc,vay_loc,vaz_loc;
	vax_loc = vay_loc = vaz_loc = 0.f;
	for (int n=0; n<ScaLBL_Comm->LastExterior(); n++){
		vax_loc += VELOCITY[n];
		vay_loc += VELOCITY[Np+n];
		vaz_loc += VELOCITY[2*Np+n];
		count_loc+=1.0;
	}
	
	for (int n=ScaLBL_Comm->FirstInterior(); n<ScaLBL_Comm->LastInterior(); n++){
		vax_loc += VELOCITY[n];
		vay_loc += VELOCITY[Np+n];
		vaz_loc += VELOCITY[2*Np+n];
		count_loc+=1.0;
	}
	MPI_Allreduce(&vax_loc,&vax,1,MPI_DOUBLE,MPI_SUM,Mask->Comm);
	MPI_Allreduce(&vay_loc,&vay,1,MPI_DOUBLE,MPI_SUM,Mask->Comm);
	MPI_Allreduce(&vaz_loc,&vaz,1,MPI_DOUBLE,MPI_SUM,Mask->Comm);
	MPI_Allreduce(&count_loc,&count,1,MPI_DOUBLE,MPI_SUM,Mask->Comm);
	
	vax /= count;
	vay /= count;
	vaz /= count;
	
	double mu = (tau-0.5)/3.f;
	if (rank==0) printf("Fx Fy Fz mu Vs As Js Xs vx vy vz\n");
	if (rank==0) printf("%.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g\n",Fx, Fy, Fz, mu, 
						Morphology.V(),Morphology.A(),Morphology.J(),Morphology.X(),vax,vay,vaz);
						*/
    vis_db = db->getDatabase("Visualization");
    if (vis_db->getWithDefault<bool>("write_silo", false)) {

        std::vector<IO::MeshDataStruct> visData;
        fillHalo<double> fillData(Dm->Comm, Dm->rank_info,
                                  {Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2},
                                  {1, 1, 1}, 0, 1);

        auto VxVar = std::make_shared<IO::Variable>();
        auto VyVar = std::make_shared<IO::Variable>();
        auto VzVar = std::make_shared<IO::Variable>();
        auto SignDistVar = std::make_shared<IO::Variable>();

        IO::initialize("", format, false);
        // Create the MeshDataStruct
        visData.resize(1);
        visData[0].meshName = "domain";
        visData[0].mesh = std::make_shared<IO::DomainMesh>(
            Dm->rank_info, Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2, Dm->Lx, Dm->Ly,
            Dm->Lz);
        SignDistVar->name = "SignDist";
        SignDistVar->type = IO::VariableType::VolumeVariable;
        SignDistVar->dim = 1;
        SignDistVar->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
        visData[0].vars.push_back(SignDistVar);

        VxVar->name = "Velocity_x";
        VxVar->type = IO::VariableType::VolumeVariable;
        VxVar->dim = 1;
        VxVar->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
        visData[0].vars.push_back(VxVar);
        VyVar->name = "Velocity_y";
        VyVar->type = IO::VariableType::VolumeVariable;
        VyVar->dim = 1;
        VyVar->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
        visData[0].vars.push_back(VyVar);
        VzVar->name = "Velocity_z";
        VzVar->type = IO::VariableType::VolumeVariable;
        VzVar->dim = 1;
        VzVar->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
        visData[0].vars.push_back(VzVar);

        Array<double> &SignData = visData[0].vars[0]->data;
        Array<double> &VelxData = visData[0].vars[1]->data;
        Array<double> &VelyData = visData[0].vars[2]->data;
        Array<double> &VelzData = visData[0].vars[3]->data;

        ASSERT(visData[0].vars[0]->name == "SignDist");
        ASSERT(visData[0].vars[1]->name == "Velocity_x");
        ASSERT(visData[0].vars[2]->name == "Velocity_y");
        ASSERT(visData[0].vars[3]->name == "Velocity_z");

        fillData.copy(Distance, SignData);
        fillData.copy(Velocity_x, VelxData);
        fillData.copy(Velocity_y, VelyData);
        fillData.copy(Velocity_z, VelzData);

        IO::writeData(timestep, visData, Dm->Comm);
    }
}
