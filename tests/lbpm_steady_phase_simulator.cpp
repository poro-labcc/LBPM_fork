#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "common/ScaLBL.h"
#include "common/Communication.h"
#include "analysis/TwoPhase.h"
#include "common/MPI.h"
#include "models/SteadyPhaseModel.h"

using namespace std;

int main(int argc, char **argv)
{
    Utilities::startup( argc, argv, false );
    Utilities::MPI comm( MPI_COMM_WORLD );
    int rank = comm.getRank();
    int nprocs = comm.getSize();
    {
        if (rank == 0){
            printf("********************************************************\n");
            printf("Running Steady Phase Permeability Calculation \n");
            printf("********************************************************\n");
        }
        int device = ScaLBL_SetDevice(rank);
        NULL_USE( device );
        ScaLBL_DeviceBarrier();
        comm.barrier();
        
        ScaLBL_SteadyPhaseModel Steady(rank, nprocs, comm);
        auto filename = argv[1];
        
        Steady.ReadParams(filename);
        Steady.SetDomain();    
        Steady.ReadInput();
        Steady.Create();       
        Steady.Initialize();   
        Steady.Run();   
        Steady.VelocityField();
        cout << flush;
    }
    Utilities::shutdown();
    return 0;
}