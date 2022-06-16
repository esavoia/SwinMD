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


#define DEFAULT_INFILE "config.txt"

#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <mpi.h>          // MPI
#include "utils/Mpi_V.h"  // MPI ssize, rrank

#include "Simulation.h"
#include "Errors.h"


int main(int argc, char **argv) {

    // Timing
    time_t elapsed_time;
    elapsed_time = time(NULL);
    
    // Initialize MPI.
    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rrank);
    MPI_Comm_size(MPI_COMM_WORLD, &ssize);
    
    // Before everything check the cluster size
    if(ssize < 1){
        cout <<"[FATAL ERROR] !!!--- at least 2 processors must be initialized ---!!!"<<endl;
        MPI_Finalize();
        exit(FATAL_ERROR);
    }
    
    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    // Where am I?
    printf("Hello from processor %s, rank %d out of %d processors\n", processor_name, rrank, ssize);
    
    
    // Initialise Simulation Object
    Simulation *sim = NULL;
    
    // With no arguments use default input file name
    sim = new Simulation( (argc <=1) ? DEFAULT_INFILE : argv[1] ); // One liner
    /* Explicit Version
     if(argc <= 1) {
        sim = new Simulation(DEFAULT_INFILE);
     } else {  // argc > 1
        sim = new Simulation(argv[1]);
     }
     */
    if (sim == NULL) {
        cerr << "[FATAL ERROR] !!!--- Can't to create Simulation object ---!!!" << endl;
        exit(FATAL_ERROR);
    }
    
    // Let's go!
    sim->setup_simulation();
    sim->run();
    
    // Done with MPI
    MPI_Finalize();
    
    // Run timing
    elapsed_time = time(NULL) - elapsed_time;
    
    if(rrank == 0){
        cout<<"Computation time is...... " << elapsed_time/3600<<"hr: "\
        <<(elapsed_time%3600)/60<<"min"\
        <<endl;
    }
    
    return 0;
}
