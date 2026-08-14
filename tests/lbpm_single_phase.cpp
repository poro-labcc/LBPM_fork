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
#include "models/MRTModel.h"
//#define WRITE_SURFACES

/*
 * Veryfication script to validate the single phase initialization routines. Intended to be removed afterwards:
 * Remove:
 *  This script
 *  The function ScaLBL_MRTModel::Run_Timesteps
 *  This script's definition in MRTModel.h
 *  The cmake item in LBPM_source/tests/CMakeList.txt
 */

using namespace std;


int main(int argc, char **argv)
{
	// Initialize MPI
    Utilities::startup( argc, argv, false );
    Utilities::MPI comm( MPI_COMM_WORLD );
    int rank = comm.getRank();
    int nprocs = comm.getSize();
	{
		if (rank == 0){
			printf("*************************************************************\n");
			printf("Running Single Phase with Custom Initialization (MRT Model) \n");
			printf("*************************************************************\n");
		}
		// Initialize compute device
		int device=ScaLBL_SetDevice(rank);
        NULL_USE( device );
		ScaLBL_DeviceBarrier();
		comm.barrier();

		auto filename = argv[1];
		std::string coord_input = (argc >= 3) ? argv[2] : "";

		std::vector<int> flat_coords;
		if (!coord_input.empty()) {
			std::string input = coord_input;
			for (char &c : input) if (c == '(' || c == ')' || c == ',') c = ' ';
			std::stringstream ss(input);
			int val;
			while (ss >> val) flat_coords.push_back(val);
		}

		// 2. Parse the --init flag
		int init = 0; // Default case
		for (int idx = 1; idx < argc; idx++) {
			if (std::string(argv[idx]) == "--init" && idx + 1 < argc) {
				init = std::stoi(argv[idx + 1]);
				break;
			}
		}

		// 3. Logic branches for initialization
		ScaLBL_MRTModel MRT(rank, nprocs, comm);
		MRT.ReadParams(filename);
		MRT.SetDomain();
		MRT.ReadInput();
		MRT.Create();
		
		if (init == 1) {
			if (rank == 0) printf("--> Init. Mode 1: Collecting Ux,Uy,Uz and P from file and calculating Feq.\n");
			MRT.Initialize_fEq();
		}
		else if (init == 2) {
			if (rank == 0) printf("--> Init. Mode 2: Collecting Ux,Uy,Uz and P from file and calculating Feq + Fneq\n");
			MRT.Initialize_fEqNeq();
		}
		else if (init == 3) {
			if (rank == 0) printf("--> Init. Mode 3: Collecting F from file\n");
			MRT.Initialize_Dist();
		}
		else {
			if (rank == 0) printf("--> Default Init: Calculating Feq for null velocity field\n", init);
			MRT.Initialize();
		}

		MRT.SaveFields();
		MRT.Run_Timesteps(flat_coords);
		MRT.SaveFields();
		cout << flush;
	}
    Utilities::shutdown();
}
