/** NEMDMain.cpp --
 **
 ** Copyright (C) 2003
 ** Centre for Molecular Simulation (CMS)
 ** School of Information Technology
 ** Swinburne University of Technology
 ** PO Box 21 Hawthorn, Vic 3122, Australia
 **
 ** Author: Zhongwu Zhou
 ** Email: zzhou@it.swin.edu.au
 **
 ** Edited by
 ** Author: Edoardo Savoia
 ** Email: esavoia@swin.edu.au
 **/

#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <mpi.h>          // MPI
#include "utils/Mpi_V.h"  // MPI mpi_size_g, mpi_rank_g

#include "Simulation.h"
#include "Errors.h" // NO_ERR
#include "Defaults.h" // DEFAULT_INFILE

int mpi_size_g;
int mpi_rank_g;


int main(int argc, char **argv) {

    // Timing
    time_t elapsed_time;
    elapsed_time = time(NULL);

    // Initialize MPI.
    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_g);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_g);

    // Before everything check the cluster size
    if(mpi_size_g < 1){
		MPI_Finalize();
		FATALMSG("At least 2 processors must be initialized!")
    }
	// Log something nice:
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len); // Get the name of the processor
    // printf("Hello from processor %s, rank %d out of %d processors\n", processor_name, mpi_rank_g, mpi_size_g);
    LOG("Hello from processor " << processor_name <<", rank " << mpi_rank_g << " out of "<< mpi_size_g << " processors")


    // Initialise Simulation Object
    Simulation *sim = NULL;
	/* THIS IS BAD!!! */
    // With no arguments use default input file name
    //sim = new Simulation( (argc <=1) ? DEFAULT_INFILE : argv[1] ); // One liner
    /* Explicit Version
     if(argc <= 1) {
        sim = new Simulation(DEFAULT_INFILE);
     } else {  // argc > 1
        sim = new Simulation(argv[1]);
     }
     */
	 /* BAD ENDS HERE */

	/* OPTION 1: Using a Builder */
	// SimBuilder *builder = new SimBuilder( (argc <=1) ? DEFAULT_INFILE : argv[1] );
 	// Simulation *sim = builder->buildSimulation();
 	/* OPTION 2: Using a static build() method */
    /* OPTION 3: Forget about checking the result and let the build() static method fail and exit */
    
    /* Today's menu: OPTION 2 */
    sim = Simulation::build((argc <=1) ? DEFAULT_INFILE : argv[1]);
	if (sim == NULL) {
        MPI_Finalize();
        FATALMSG("Can't create Simulation object")
    }
    
    
    // sim->setup_simulation(); // Absorbed by Simulation::build()
    // Let's go!
    //sim->run();

    // Done with MPI
    MPI_Finalize();

    // Run timing
    elapsed_time = time(NULL) - elapsed_time;

    if(mpi_rank_g == 0){
        LOG("Computation time is...... " << elapsed_time/3600<<"hr: "\
        <<(elapsed_time%3600)/60<<"min"\
        )
    }

    return NO_ERR;
}
