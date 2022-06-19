/* NEMD_defs.h -- definitions and macros for NEMD codes
   Copyright (C) 2002

   Zhongwu Zhou
   Centre for Molecular Simulation (CMS)
   School of Information Technology
   Swinburne University of Technology
   PO Box 21 Hawthorn, Vic 3122, Australia
   Email: zzhou@it.swin.edu.au

*/


#ifndef NEMD_DEFS_H
#define NEMD_DEFS_H

using namespace std;

#ifdef  MAX
#undef  MAX
#endif
#define MAX(A,B)		(((A)>(B))?(A):(B))
#define MAX3(A,B,C)	 MAX((A), MAX((B),(C)))

#ifdef  MIN
#undef  MIN
#endif
#define MIN(A,B)		(((A)<(B))?(A):(B))
#define MIN3(A,B,C)	 MIN((A), MIN((B),(C)))

#ifdef ANINT
#undef ANINT
#endif
#define ANINT(A) ((int)((A)+((A)>=0?.5:-.5)))
//#define ANINT(A) ((int)((A)+((A)<(-0.5) ? -.5: ((A) > 0.5 ? 0.5 : 0))))

#define SQR(A)		((A)*(A))
#define CUBE(A)		((A)*(A)*(A))
#define RIJSQR(A, B)	(SQR(A.x-B.x)+SQR(A.y-B.y)+SQR(A.z-B.z))
#define RIJ(A, B)		sqrt(RIJSQR((A), (B)))
#define CBRT(x)		(pow((x),1./3.))
#define SIGN(x, y)	((y) >= 0 ? fabs(x) : -fabs(x))

/* index for virial array */
#define XX      0
#define XY      1
#define XZ      2
#define YX      3
#define YY      4
#define YZ      5
#define ZX      6
#define ZY      7
#define ZZ      8

/* definition and typedef for md codes */
#define SHORT_STR_LEN	8
#define MD_STR_LEN	64
#define LONG_STR_LEN	1024
#define FALSE		0
#define TRUE		1

typedef double		Double;		/* assume 8 bytes double and 4 bytes int */
typedef int			Int;

typedef unsigned	long int 	time_mt;    /* Possible largest integer */
typedef unsigned	long		size_mt;    /* Wide type for passing sizeof	      */
typedef char	SHORT_STR[SHORT_STR_LEN];
typedef char	MD_STR[MD_STR_LEN];
typedef char	LONG_STR[LONG_STR_LEN];
typedef void	gptr;                   /* generic ptr */

/**  Constants **/
//--------------------
// NUMERICAL CONSTANTS - from PROTOMOL/base/protomol_constants.h
//--------------------
static const Double ALPHAMIN            = 1e-7;
static const Double PI                  = 3.14159265358979323846;
static const Double SQRT_PI             = 1.7724538509055160273;
// static const Double DTOR                = (PI/180.0);

static const Double SI_C                = 299792458.0;    // [m/s]
static const Double SI_COULOMB_FACTOR   = SI_C*SI_C*1e-7; // [Vm/C]
static const Double SI_ELECTRON_CHARGE  = 1.6021892e-19;  // [C]
static const Double SI_LENGTH           = 1e+9;          //  [nm]
static const Double SI_AVOGADRO         = 6.022045e+23;   // [1/mol]
static const Double SI_AMU              = 1.6605655e-27;  // [kg]
static const Double SI_KCAL             = 1/4184.0;         // [J]
static const Double SI_TIME             = 1e+12;          // [ps]
static const Double SI_BOLTZMANN        = 1.380662e-23;   // [J/K]

static const Double FORCE_TO_SI         = SI_LENGTH/(SI_AVOGADRO*SI_KCAL);
static const Double ENERGY_TO_SI        = 1.0/(SI_AVOGADRO*SI_KCAL);
static const Double LENGTH_TO_SI        = 1.0/SI_LENGTH;
static const Double TIME_TO_SI          = 1.0/SI_TIME;


//program unit & constants
static const Double MUNIT               = SI_AMU;         // atomic mass
static const Double LUNIT               = 1.0e-9;         // nm

// sqrt(1.0/SI_LENGTH*SI_C*SI_C*1e-7*SI_ELECTRON_CHARGE*
//      SI_ELECTRON_CHARGE*SI_AVOGADRO/1000) (sqrt[kJ.mol-1.nm.e-2])
static const Double SQRTCOULOMBCONSTANT = 11.78708976;    // 1/sqrt(4PI*EPS0)
static const Double COUlOMBS     =  138.93548501;           // 1/(4PI*EPS0) 86.71955;
// static const Double SQRTCOULOMBCONSTANT = 1.33207;     // 1/sqrt(4PI*EPS0)

// SI_AVOGADRO*SI_BOLTZMANN/1000;
static const Double BOLTZMAN            = 0.0083145112;


// sqrt(1.0/SI_LENGTH*SI_C*SI_C*1e-7*SI_ELECTRON_CHARGE*
//      SI_ELECTRON_CHARGE*SI_AVOGADRO*SI_KCAL)
// static const Double SQRTCOULOMBCONSTANT = 18.2226123264;    // 1/sqrt(4PI*EPS0)
// SI_AVOGADRO*SI_BOLTZMANN*SI_KCAL;
// = 0.0019871913704086990
// static const Double BOLTZMAN            = 0.00198719137;

// SI_TIME*sqrt(1e-3*1.0/SI_LENGTH*1.0/SI_LENGTH*SI_KCAL)
// = 48.888213315262199
static const Double TIMEFACTOR          = 48.88821290839616;   // TIMEUNIT is 1/sqrt(4.184e-4)


static const Int  FASTDELTAMAX          = 32;
static const Double PDBVELSCALINGFACTOR = 20.45482706;

/*

// #define DEBUGMSG(ofp, args) { ofp << "DBG: " << __FILE__ << " (" << __LINE__ << ") " << args << endl;}

#define DEBUGMSG(args)   { cerr << "DBG: " << __FILE__ << " (" << __LINE__ << ") " << args << endl;}
        // cerr << "DBG: " << __FILE__ << '(' << __LINE__ << "): ";
        //cerr << "at timestep [" << timestepnr <<"] "<< args << endl;}

#define ERRORMSG(args)   {  \
        cerr << "ERROR: " << __FILE__ << '(' << __LINE__ << "): "; \
        cerr <<  args << endl;    \
        exit (1);   }

#define ERRORMSG(ofp, args)   {  \
        ofp << "ERROR: " << __FILE__ << '(' << __LINE__ << "): "; \
        ofp <<  args << endl;    \
        exit (1);   }
*/


//--------------------
// UNITS
//--------------------
//
//  r : 1 [AAm]          = 1e-10 [m]
//  t : 1 [s']           = 1e+15 * sqrt((1.660e-27 * 1e-20 * 6.022035e+23)/4184) [fs]
//                       = 1e+15 * sqrt((1e-3 * 1e-20)/4184) [fs]
//  v : 1 [AAm/s']       = 1e+15 * sqrt((1.660e-27 * 1e-20 * 6.022035e+23)/4184) [m/fs]
//                       = 1e+15 * sqrt((1e-3 * 1e-20)/4184) [m/fs]
//  F : 1 [kcal/mol AAm] = 4184 / (6.022045e+23 * 1e-10) [N]
//  E : 1 [kcal/mol]     = 4184 / 6.022045e+23 [J]
//  m : 1 [amu]          = 1.6605e-27 [kg]
//  e : 1 [e]            = 1.602e-19 [C]


#endif
