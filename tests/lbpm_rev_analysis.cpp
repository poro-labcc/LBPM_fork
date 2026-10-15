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
#include "../analysis/rev_analysis.h"

//#define WRITE_SURFACES

/*
 * Standalone representative elementary volume analysis
 * James E. McClure 2013-2014
 */

using namespace std;


int main(int argc, char **argv)
{
	if (argc < 2) {
		std::cerr << "Usage: lbpm_rev_analysis params.db\n";
		return 1;
	}
	// Initialize MPI
    Utilities::startup( argc, argv, false );
    Utilities::MPI comm( MPI_COMM_WORLD );
    int rank = comm.getRank();
    int nprocs = comm.getSize();
	{
		if (rank == 0){
			printf("********************************************************\n");
			printf("Running REV Analysis \n");
			printf("********************************************************\n");
		}
		comm.barrier();
		
		ScaLBL_REV REV(rank,nprocs,comm);

		auto filename = argv[1];
		REV.ReadParams(filename);
		REV.SetDomain();    // this reads in the domain 
		REV.DetRevAnalysis();
		REV.StatRevAnalysis();
		cout << flush;
	}
    Utilities::shutdown();
}
