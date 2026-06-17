#ifndef ScaLBL_SteadyPhaseModel_H
#define ScaLBL_SteadyPhaseModel_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "common/ScaLBL.h"
#include "common/Communication.h"
#include "common/MPI.h"
#include "analysis/TwoPhase.h"
#include "analysis/Minkowski.h"
#include "ProfilerApp.h"

class ScaLBL_SteadyPhaseModel {
public:
    ScaLBL_SteadyPhaseModel(int RANK, int NP, const Utilities::MPI &COMM);
    ~ScaLBL_SteadyPhaseModel();
    
    void ReadParams(std::string filename);
    void ReadParams(std::shared_ptr<Database> db0);
    void SetDomain();
    void ReadInput();
    void Create();
    void Initialize();
    void ComputeNormals();
    void Run();
    void VelocityField();

    bool Restart, pBC;
    int timestep, timestepMax;
    int ANALYSIS_INTERVAL;
    int BoundaryCondition;
    
    double tau_A, tau_B;
    double mu_A, mu_B;
    double sigma;
    double tolerance;
    double Fx, Fy, Fz, flux;
    double din, dout;
    double factor;

    int Nx, Ny, Nz, N, Np;
    int rank, nprocx, nprocy, nprocz, nprocs;
    double Lx, Ly, Lz;

    std::shared_ptr<Domain> Dm;   
    std::shared_ptr<Domain> Mask; 
    std::shared_ptr<ScaLBL_Communicator> ScaLBL_Comm;
    
    std::shared_ptr<Database> db;
    std::shared_ptr<Database> domain_db;
    std::shared_ptr<Database> steady_db;
    std::shared_ptr<Database> vis_db;

    IntArray Map;
    DoubleArray SignDistance;
    
    int *NeighborList;
    double *fq;
    double *Velocity;
    double *Pressure;
    int *Phi;
    double *ColorGrad;

private:
    Utilities::MPI comm;
    size_t dist_mem_size;
    size_t neighborSize;
    char LocalRankString[8];
    char LocalRankFilename[40];
    char LocalRestartFile[40];

    void LoadParams(std::shared_ptr<Database> db0);
};

#endif