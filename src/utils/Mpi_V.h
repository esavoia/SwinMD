/** MPI variables definition
 **
 ** principal variables:rank size are defined in this head file
 **
 **
 **/
#include <mpi.h>

#ifndef MPI_V
#define MPI_V

// These are global variables used only by main
int mpi_size_g;
int mpi_rank_g;

#endif
