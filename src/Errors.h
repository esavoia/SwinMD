//
//  Errors.h
//  SwinMD
//
//  Created by Edoardo Savoia on 15/6/2022.
//

#ifndef Errors_h
#define Errors_h

#include "utils/Mpi_V.h"

// Error Codes
#define NO_ERR 0
#define JUST_ERR 1
#define FATAL_ERROR 53
#define RUNTIME_ERROR 42
#define UNEXPECTED_ERROR 99

// Colours
#ifdef USE_COLOURS
#define COLOUR_RESET "\x1B[0m"
#define COLOUR_RED   "\x1B[31m"
#define COLOUR_YELLOW "\x1b[33m"
#define COLOUR_BLUE "\x1b[34m"
#define COLOUR_GREEN "\x1b[32m"
#define COLOUR_MAGENTA "\x1b[35m"
#define COLOUR_CYAN  "\x1b[36m"
#define COLOUR_WHITE "\x1b[37m"
#else

#define COLOUR_RESET
#define COLOUR_RED
#define COLOUR_YELLOW
#define COLOUR_BLUE
#define COLOUR_GREEN
#define COLOUR_MAGENTA
#define COLOUR_CYAN
#define COLOUR_WHITE

#endif

// Message Intro
#define FATAL_INTRO   "[" COLOUR_RED "FATAL" COLOUR_RESET "] "
#define ERROR_INTRO   "[" COLOUR_YELLOW "ERROR" COLOUR_RESET "] "
#define RUNTIME_INTRO "[" COLOUR_MAGENTA "RUNTIME" COLOUR_RESET "] "
#define WARN_INTRO    "[" COLOUR_GREEN "WARNING" COLOUR_RESET "] "
#define DEBUG_INTRO   "[" COLOUR_BLUE "DEBUG" COLOUR_RESET "] "

#ifndef USE_FILE_AND_LINE

// Always returns FATAL_ERROR code [No LOG_FATAL]
#define FATALMSG(msg)   {  cerr << FATAL_INTRO << __FILE__ << '(' << __LINE__ << "): " << msg << endl; exit(FATAL_ERROR); }
// #define LOG_FATAL(msg)   cerr << FATAL_INTRO << __FILE__ << '(' << __LINE__ << "): " << msg << endl;

// Always returns UNEXPECTED_ERROR code [No LOG_ERROR]
#define ERRORMSG(msg)   {  cerr << ERROR_INTRO << __FILE__ << '(' << __LINE__ << "): " << msg << endl; exit(UNEXPECTED_ERROR); }
// #define LOG_ERROR(msg)  cerr << ERROR_INTRO << __FILE__ << '(' << __LINE__ << "): " << msg << endl;

// Can lead to returning RUNTIME_ERROR code [Has LOG_ERR]
#define ERRMSG(msg)   {  cerr << RUNTIME_INTRO << __FILE__ << '(' << __LINE__ << "): " << msg << endl; exit(RUNTIME_ERROR); }
#define LOG_ERR(msg)  cerr << RUNTIME_INTRO << __FILE__ << '(' << __LINE__ << "): " << msg << endl;

// No error code is returned after this
#define WARNMSG(msg)   { cerr << WARN_INTRO << __FILE__ << " (" << __LINE__  << ") " << msg << endl; }

// #define DEBUGMSG(msg)   { cerr << DEBUG_INTRO << __FILE__ << " (" << __LINE__  << ") " << msg << endl; }
#ifdef DEBUG
#define DEBUGMSG(msg)   { cerr << DEBUG_INTRO << __FILE__  << " ("  <<  __LINE__  << ") " << msg << endl; }
#else
#define DEBUGMSG(msg) {;}
#endif

#else // USE_MPI_RANK_AND_SIZE

// Always returns FATAL_ERROR code [No LOG_FATAL]
#define FATALMSG(msg)   {  cerr << FATAL_INTRO << mpi_rank_g << '(' << mpi_size_g << "): " << msg << endl; exit(FATAL_ERROR); }
// #define LOG_FATAL(msg)   cerr << FATAL_INTRO << mpi_rank_g << '(' << mpi_size_g << "): " << msg << endl;

// Always returns UNEXPECTED_ERROR code [No LOG_ERROR]
#define ERRORMSG(msg)   {  cerr << ERROR_INTRO << mpi_rank_g << '(' << mpi_size_g << "): " << msg << endl; exit(UNEXPECTED_ERROR); }
// #define LOG_ERROR(msg)  cerr << ERROR_INTRO << mpi_rank_g << '(' << mpi_size_g << "): " << msg << endl;

// Can lead to returning RUNTIME_ERROR code [Has LOG_ERR]
#define ERRMSG(msg)   {  cerr << RUNTIME_INTRO << mpi_rank_g << '(' << mpi_size_g << "): " << msg << endl; exit(RUNTIME_ERROR); }
#define LOG_ERR(msg)  cerr << RUNTIME_INTRO << mpi_rank_g << '(' << mpi_size_g << "): " << msg << endl;

// No error code is returned after this
#define WARNMSG(msg)   { cerr << WARN_INTRO << mpi_rank_g << " (" << mpi_size_g  << ") " << msg << endl; }

// #define DEBUGMSG(msg)   { cerr << DEBUG_INTRO << mpi_rank_g << " (" << mpi_size_g  << ") " << msg << endl; }
#ifdef DEBUG
#define DEBUGMSG(msg)   { cerr << DEBUG_INTRO << mpi_rank_g  << " ("  <<  mpi_size_g  << ") " << msg << endl; }
#else
#define DEBUGMSG(msg) {;}
#endif


#endif


// LOG function just for sake of completeness

#define LOG(msg) cout << "[ LOG ] " msg << endl;


#endif /* Errors_h */
