/** MPI variables definition
 **
 ** principal variables:rank size are defined in this head file
 **
 **
 **/
#include <mpi.h>

#ifndef MPI_V
#define MPI_V

// FIX_THIS: These are global variables are used only by main
// Now are used also by Simulation
extern int mpi_size_g;
extern int mpi_rank_g;

#endif
