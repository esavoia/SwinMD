/** NccForce.h -- 
 **
 ** Copyright (C) 2004
 ** Centre for Molecular Simulation (CMS)
 ** School of Information Technology
 ** Swinburne University of Technology
 ** PO Box 21 Hawthorn, Vic 3122, Australia
 **
 ** Author: Jianhui Li
 ** Email: jli@it.swin.edu.au
 **/
//JC this is a child class of Force and on the same class level as the LFForce and Ewarld class
//JC it is called in LFIntegrator, which is actually the one currently used in the calculation!
//JC accordingly the LFIntegrator needs certain change,when the NCCForce is implemented in future 
//JC calculation, 
//JC NccForce is defined as a child class of Force.
#ifndef NccForce_H
#define NccForce_H

#include "../Force.h" 
#include <iomanip>
#include <time.h>
#include "mpi.h"
#include "../NEMD_defs.h"

#include <iostream>
#include <math.h>

using namespace std;


#define TAB_SIZE  15000
#define DX   1.0e-12

//E_Inter
#define Q    0.400    //0.717484  // Sqrt(Q2)************************************************ WAS 0.4 ************************************************
#define Q2   0.514783  //0.565117   //(0.5525372*0.5525372) //(0.4238*0.4238) /SQRTCOULOMBCONSTANT) position positive point scaled charge [Coulomb]
#define Aoo  1734.196   //1864.271482   //2.8682660 //12.008851 // well depth of O-O intermolecular interaction [Hartree]
#define Boo  -2.726696  //-2.753110  //-1.2390449 // exponential coefficient of Pairwise potential of O-O interaction[1/Bohr] 
#define Ahh  1.061887  //0.662712   //0.1531258 //0.6411068 // well depth of H-H intermolecular interaction[Hartree]
#define Bhh  -1.460975 // -1.299982 //-1.0133393 // exponential coefficient of Pairwise potential of H-H interaction
#define Aoh  2.319395  //2.684452   //41.5624830 //174.01372 // well depth of O-H intermolecular interaction
#define Boh  -1.567367 //-1.439787  //-1.2227832 // exponential coefficient of Pairwise potential of O-H interaction 
#define A_oh 0.436006  //0.675342   //40.0772130 //167.77952 // corrective well depth coefficient of O-H intermolecular interaction[kj/mol]
#define B_oh -1.181792 //-1.141494  //-1.2111490 // corrective exponential coefficient of O-H intermolecular interaction
#define R_OP 0.505783   

//E_CI
/*
#define Q    0.751743  // Sqrt(Q2)
#define Q2   0.565117     //(0.5525372*0.5525372) //(0.4238*0.4238) /SQRTCOULOMBCONSTANT) position positive point scaled charge [Coulomb]
#define Aoo  1864.271482  //2.8682660 //12.008851 // well depth of O-O intermolecular interaction [Hartree]
#define Boo  -2.753110    //-1.2390449 // exponential coefficient of Pairwise potential of O-O interaction[1/Bohr] 
#define Ahh  0.662712     //0.1531258 //0.6411068 // well depth of H-H intermolecular interaction[Hartree]
#define Bhh  -1.299982    //-1.0133393 // exponential coefficient of Pairwise potential of H-H interaction
#define Aoh  2.684452     //41.5624830 //174.01372 // well depth of O-H intermolecular interaction
#define Boh  -1.439787    //-1.2227832 // exponential coefficient of Pairwise potential of O-H interaction 
#define A_oh 0.675342     //40.0772130 //167.77952 // corrective well depth coefficient of O-H intermolecular interaction[kj/mol]
#define B_oh -1.141494    //-1.2111490 // corrective exponential coefficient of O-H intermolecular interaction
#define R_OP 0.487741
*/

#define Aph  14.3170660 //59.942665 // short distance corrective coefficient 1[kj/mol]
#define Bph  -2.7567793 // short distance corrective coefficient 2
#define Apo  51.1496400 //214.15322 // short distance corrective coefficient 3[kj/mol]
#define Bpo  -2.1794969 // short distance corrective coefficient 4
#define Bohr_Radius 5.29189379e-11 //[Meter] atomic length unit
#define Hartree 4.3597438134e-18   //[Joule] atomic energy unit
#define E_Charge 1.60217646263e-19 //[C] atomic electrical charge unit
#define E_Mass 9.10e-31            //[Kg] atomic mass unit
#define Avogadro 6.02e+23          //Avogadro Constant, the number of atoms comtaining in one Mole 
#define Amu 1.66053e-27            //[Kg] Atmic Mass Unit
#define F_Length (1.0e-9/Bohr_Radius) // Length Scaled Factor
#define F_Energy 2625.49992        // Energy Scaled Unit(from Hartree to Kilo Joule/Mole)
#define F_Force  (Hartree/Bohr_Radius)*(1.0/Amu)*(1.0/1.0e-9)/(1.0/1.0e-24) // Force scaled factor

// Jc: following parameters are defined for the three new 3-body interactions
#define eps3B   0.64852                                             // from ZW SPC/E potential in [KJ/Mol]
#define sigma3B 0.3166*0.3166*0.3166*(F_Length*F_Length*F_Length)   // from ZW SPC/E potential in [Nm3]
#define alpha3B 1.45e-3*(F_Length*F_Length*F_Length)                // from ZW SPC/E potential in [Nm3]
//#define Vddd 0.75*(alpha3B)*(sigma3B*sigma3B)*4*(eps3B/F_Energy)    // from Y.singh, J. Phys.B: Atom. Molec.Phys, 1971. Vol.4 
#define Vddd 2.459E-06  // kJnm9/mol // 287.95 //518.3  //[a.u.] from  P.J. Leonard, Theoretical Chemistry: Advances &PrespectivesVol.1,1975  
#define Coef3F Vddd*F_Force        //scaling coeffcient for force in 3-body interaction
#define Coef3E Vddd*F_Energy       //coefficient of energy
// Jc: define polarization coefficients
#define polar  0.00144000  //polarizability of H2O[nm3], taken from Literature 20050421001    
#define Pf_Energy 400.18400          // 1 Kcal/Mol = 4.18400 KJ/Mol
#define Pf_Force  4100.855           // (Kcal*4.184/(Avogadro*Angstrom))/(Amu*1.0e-9/1.0e-24) cf, calcualtion on 07Mar2006

class NccForce : public Force {
    /**
     ** Data Member - 
     **/
    private:
    vector<Vector3>  kVector;     // reciprocal lattice vector(2PI*nx/Lx, 2PI*ny/Ly, 2PI*nz/Lz)
    vector<double>   kModulu;     // modulus of kVector
    Vector3  DipoleSum;           //Jc: added by Jianhui
    double A[1500][1500];
    double KSwitch[500][4500];    // storing [1-K(r)] for all the Rij
    double Pi[1500];
    double Tij[3][3];             // Jc: T matrix 
    double Ei0[1500],b[1500];     // Jc: b vector 
    Double *expK;                 // pre-computed exp( ) value for kVector
    // Int nxMax, nyMax, nzMax;   // maximum values of (nx, ny, nz)
    Parameters * params;
    Double cutOff,cutOff2,switchDist;
    Double c1,c2,c3,c4;
    Double eLrc, vLrc;
    int    *p;
    bool   switchOn;
    bool   computeCoulomb;
    double Dummy[500*3];
//    double Dummy[5322];
//    double Dummy[1500];

// JC other variables in the NccForce calculation
    Int numAtoms, counter;        // number of atoms in the system = nAtoms in the ensamble
    Int numMols;

    Molecule *mols;
    Double realCut, realCut2;     // cutoff and cut off square for real term
    Double kCut;                  // cutoff of reciprocal term
    Double alpha;                 // real/reciprocal space particion parameter, or Ewald convergence parameter
    Double alpha2, alpha2PI, alphaR4; 
    Double stepr;

    Double boxLx, boxLy, boxLz;
    Double volume, volume2, volumer;

    Double* sinCosA;              // look up table for reciprocal term
    Double* forceTable;           // look up table for real term
    Double* potentialTable;       
    Int tableSize;

    Double realEnergy, longEnergy, selfEnergy, chargedEnergy, intraMolEnergy, surfDipoleEnergy;
    Double longVirial[9];
    Double WPotential,CPotential; // Wpotential = van der walss, Cpotential + columbic Potential

    bool doUpdate, correctSurfDipole;
    double R13,R14,R23,R24,R78,R81,R82,R73,R74,R56,R53,R54,R61,R62,R87,R76,R85;
    double P2Sum,PSum;      // Jc: the time average of dipole moments
    ofstream *ofDummy,*ofDummy1;
    ofstream *ofCoordinate;
    ofstream *ofpairlist;  // JC test the pair generation of the system
    ofstream *ofEnergy;    // JC save the Potential of the system
    ofstream *ofMolEng;    // JC save potential energy of the certain Molecule
    ofstream *ofminDist;   // JC save the minimum distance of firt O atom to other O atoms in the system
    ofstream *outTimeFile; // ("procsTime.txt",ios::out); // Jc:save computational time
    int CNum;              // JC compute calcualtion times
    int rank,size;         // Jc: MPI variables
    int lBond,uBond,lJBreak,lKBreak,uJBreak,uKBreak; // Jc: seperators in the precise method 
    bool computeEwald; 
    double halfLx,halfLy, halfLz;               
    int particleIndexXL;
    int particleIndexXU;
    int particleIndexYL;
    int particleIndexYU;
    int particleIndexZL;
    int particleIndexZU;

    public:                       // Constructor and Destructor
    NccForce(Ensemble *ensemble);
    ~NccForce() { ; };

    /**
     ** methods
     **/
    public:
    virtual void compute();     // from class Force
    virtual void DummyPosition();
    virtual void seperator();   // calculate the seperators for precise load balance
    virtual void write_force_info(ofstream& of);
    virtual void write_energy(ofstream& of);
    virtual void long_range_correct();
    void compute_long();        

    private:
    void GaussElimination();
    void ConjugateGradient();
// Jc: these two Switch functions L(r) and K(r) are taken from the by F. H. Stillinger and Carl W. David
// Jc: length unit is Angstrom and the Energy unit is Kcal/Mol cf. reference 20060221001PE 
    inline double switchFuncK (double r) // JC designed by Jianhui Li to switch the induced dipole -charge interaction 
    {
      double kValue;
      r = r*10; // Jc: change the length unit(Nm) to Angstrom 
      double r2 = r*r;
      double r4 = r2*r2;
//      double cutoffA = cutOff*10;
       
//      double cutOffA2 = cufOff*cutOff*100;
//      double cutOffA4 = cufOff*cutOff*100*cufOff*cutOff*100;
//      kValue = (1 - 1.666667*r/cutOffA + 1.666667*r4/cutOffA4 - r*r4/(cutOffA4*cutOffA)); //Jc: from ZW'Code
/*      double r3 = r*r*r;
      double r_cutoff2 = (r - 0.9584)*(r - 0.9584); 
      kValue = r3/(r3 + 1.855785*r_cutoff2*exp(-8.0*r_cutoff2) + 16.95145727*exp(-2.702563425*r));
*/
      kValue = 1.0;
      return kValue;
    }

    inline double switchFuncL (double r) // JC designed by Jianhui Li to switch the induced dipole -charge energy 
    {
      double lValue;
      r = r*10; // Jc: change the length unit(Nm) to Angstrom 
      double r2 = r*r;
//Jc:      lValue = 1.0 - exp(-3.169888166*r)*(1 + 3.169888166*r + 5.024095492*r2 - 17.99599078*r*r2 + 23.92285*r2*r2);
      lValue = 1.0;
      return lValue;
    }

    inline void apply_pbc1(Vector3 &r)
    {
      while (r.x >= halfLx)  r.x -= boxLx;
      while (r.x < -halfLx)  r.x += boxLx;
      while (r.y >= halfLy)  r.y -= boxLy;
      while (r.y < -halfLy)  r.y += boxLy;
      while (r.z >= halfLz)  r.z -= boxLz;
      while (r.z < -halfLz)  r.z += boxLz;  
    }


    /** new methods of this class **/
    private: 
//    void set_up();
//    void init_table(Int size, Double** table);
//    void build_look_up_table();
//    void real_term(Double &realEnergy);
//    void reciprocal_term(Double &reciproEnergy);
//    void intraMolCorrect(Double &intraMolEnergy);
//    void surfDipoleCorrect(Double &surfaceDipoleEnergy);
    
};

#endif
