/** Induction.h -- 
 **
 ** Copyright (C) 2003
 ** Centre for Molecular Simulation (CMS)
 ** School of Information Technology
 ** Swinburne University of Technology
 ** PO Box 21 Hawthorn, Vic 3122, Australia
 **
 ** Author: eoyarzua
 ** Email: eoyarzua@swin.edu.au
 **/

#ifndef Induction_H
#define Induction_H

#include <vector>
#include <cmath>
///#include <math.h>
#include "../Force.h" 
#include "../utils/Local2Global.h"
#include "../utils/Vector3.h"
#include "../utils/CellManager.h"
#include <mpi.h>

using namespace std;

// #define USE_LOOK_UP_TABLE
#define TAB_SIZE  15000
#define DX  1.0e-12
//E_Inter

#define R_OP 0.505783

//
#define Aph  14.3170660 //59.942665 // short distance corrective coefficient 1[kj/mol]
#define Bph  -2.7567793 // short distance corrective coefficient 2
#define Apo  51.1496400 //214.15322 // short distance corrective coefficient 3[kj/mol]
#define Bpo  -2.1794969 // short distance corrective coefficient 4
#define Bohr_Radius 5.29189379e-11 //[Meter] atomic length unit
#define Hartree 4.3597438134e-18   //[Joule] atomic energy unit
#define Charge 1.60217646263e-11   //[C] atomic electrical charge unit
#define E_Mass 9.10e-31            //[Kg] atomic mass unit
#define Avogadro 6.02e+23          //Avogadro Constant, the number of atoms comtaining in one Mole 
#define Amu 1.66053e-27            //[Kg] Atmic Mass Unit
#define F_Length (1.0e-9/Bohr_Radius) // Length Scaled Factor
#define Aa0   0.05291772                // Bohr radius [nm]
#define Aa0met   5.29177210903E-11      // Bohr radius [m]


/**
 ** Forward declarations
 **/

class Induction : public Force {
    /**
     ** Data Member - 
     **/
    private:
    vector<Vector3>  kVector;     // reciprocal lattice vector(2PI*nx/Lx, 2PI*ny/Ly, 2PI*nz/Lz)
    vector<double>   kModulu;     // modulus of kVector
    Double *expK;                 // pre-computed exp( ) value for kVector
    // Int nxMax, nyMax, nzMax; // maximum values of (nx, ny, nz)

    Int numAtoms, counter,numMols, nAtomsx3;
    Molecule* myMols;
    Local2Global* myL2G;
    Double realCut, realCut2;     // cutoff and cut off square for real term
    Double cutoffInd, cutoffInd2; // cutoff for Induction computation
    Double kCut;                  // cutoff of reciprocal term
    Double alpha;                 // real/reciprocal space particion parameter, or Ewald convergence parameter
    Double alpha2, alpha2PI, alphaR4; 
    Double stepr;
    Double a;                     // a damping parameter

    Double boxLx, boxLy, boxLz;
    Double volume, volume2, volumer;

    double *Pi;  
    double *Ei0;  
    double **A;  
    double Tij[3][3];		 // The historically defined 3x3 Dipole-dipole tensor

    double *qQsin_i;  
    double *qQcos_i;  
    double *mksin_i;  
    double *mkcos_i;  
    double *pksin_i;  
    double *pkcos_i;  
    Double* forceTable;           // look up table for real term
    Double* potentialTable;       
    Int tableSize;

    Double realEnergy, longEnergy, selfEnergy, chargedEnergy, intraMolEnergy, surfDipoleEnergy;
    double intraVircorrect, selfdUdV, longderiv, realderiv, long2ndDeriv, real2ndDeriv;
    double surfderiv, surf2ndderiv, molderiv, mol2ndDeriv;
    Double longVirial[9];
    int prank, psize;
    int statEnsem;

    bool doUpdate, correctSurfDipole, useConstraint, doMinimization, dampTij, computeTorques;                 
    bool dampAMOEBA;
    double *Dummy;
    Vector3 *trq, *kXmsina, *kXmcosa, *KXQcosa, *KXQsina; 
    double totTrqI;

    /**
     ** constructor and destructor
     **/
    public:
    Induction(Ensemble* ensemble);
    ~Induction();

    /**
     ** methods
     **/
    public:
    virtual void compute();     // from class Force
// JC    virtual string get_force_id(); // in order to comply with the cluster compiler

    virtual void DummyPosition();
    virtual void InducedDipoles();

    virtual void write_force_info(ofstream& of);
    virtual void write_energy(ofstream& of);
    void compute_long();        // new to this class
    Double get_long_energy() { return longEnergy; }
    void get_long_virial(Double *vir)
    {
        for(Int i = XX; i <= ZZ; i++)
            vir[i] = longVirial[i];
    }

    /** new methods of this class **/
    private: 
    void set_up();
    void init_table(Int size, Double** table);
    void build_look_up_table();
    void self_correct(double &selfEnergy, double &selfderiv);
    void real_term(Double &realEnergy, double &realderiv, double &real2ndDeriv);
    void reciprocal_term(Double &reciproEnergy, double &longderiv, double &long2ndDeriv);
    void intraMolCorrect(Double &intraMolEnergy, double &intraVircorrect, double &molderiv, double &mol2ndDeriv);
    void surfDipoleCorrect(Double &surfaceDipoleEnergy, double &surfderiv, double &surf2ndderiv);
    void ConjugateGradient();
    void ParallelCG();
    void torques2forces(double &dUdVol);

    
    inline  Double lookup(Double r2, Double* table) const 
    {
        Double n, val;
        int nt;
    
        n = r2*stepr;
        nt=(int)n;
        
        val=table[nt];
        
        return val + (table[nt+1] - val)*(n - (Double)nt);
    }
};

#endif
