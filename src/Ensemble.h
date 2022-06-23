/** Ensemble.h -- 
 **
 ** Copyright (C) 2003
 ** Centre for Molecular Simulation (CMS)
 ** School of Information Technology
 ** Swinburne University of Technology
 ** PO Box 21 Hawthorn, Vic 3122, Australia
 **
 ** Author: Zhongwu Zhou
 ** Email: zzhou@it.swin.edu.au
 **/
// JC commented by Jianhui Li
#ifndef Ensemble_H
#define Ensemble_H

#include <iostream>
#include <fstream>
#include <time.h>
#include "NEMD_defs.h"
#include "SimConfiguration.h"
#include "Parameters.h"
#include "Molecule.h"
#include "Atom.h"
#include "utils/Vector3.h" 
#include "utils/CellManager.h"
#include <iomanip>
#include <math.h>

// STL vector 
#include <vector>
#include <iomanip>

using namespace std;

using std::vector;

// STL set 
#include <set>
using std::set;

// define the apply_pbc function
/*#define apply_pbc(Vector3 &r)
{
   while (r.x >= halfLx)  r.x -= boxLx;
   while (r.x < -halfLx)  r.x += boxLx;
   while (r.y >= halfLy)  r.y -= boxLy;
   while (r.y < -halfLy)  r.y += boxLy;
   while (r.z >= halfLz)  r.z -= boxLz;
   while (r.z < -halfLz)  r.z += boxLz;    
} 
*/
typedef struct  
{
    Double ljEnergy;
    Double realEnergy;
    Double longEnergy;
    Double correctEnergy;
    Double molCorrectEnergy;
    Double surfCorrectEnergy;
    Double angleEnergy;
    Double bondEnergy;
    Double atomKinEnergy;
    Double molKinEnergy;
    Double miscEnergy;
    Double totEnergy;
    Double inducEnergy;
    Double TBEnergy;
    Double realPolEnergy;
    Double longPolEnergy;
    Double correctPolEnergy;
    Double molCorrectPolEnergy;
    Double surfacePolEnergy;
    Double ubEnergy;
    Double TotTorq;
} Energy;

typedef struct
{
    Double ljUlus;
    Double ljdUlus;
    Double ljd2Ulus;
    Double virialcorrect;
    Double intravirCorrect;
    Double longDeriv;
    Double realDeriv;
    Double realSecDeriv;
    Double longSecDeriv;
    Double surfDeriv;
    Double surfSecDeriv;
    Double inducEn;
    Double inducDeriv;
    Double molDeriv;
    Double molSecDeriv;
    Double BondDeriv;
    Double BondSecDeriv;
    Double intravirPolCorrect;
    Double longPolDeriv;
    Double realPolDeriv;
    Double realSecPolDeriv;
    Double longSecPolDeriv;
    Double surfPolDeriv;
    Double surfSecPolDeriv;
    Double molPolDeriv;
    Double molSecPolDeriv;
    Double selfPolDeriv;
    Double ubDeriv;
    Double ubSecDeriv;
} Lustig;                    //add here Lustig averages for future potentials 

typedef struct
{
    Double bUF;
    Double bUF2;
    Double bUFt;
    Double bUF2t;
    Double dbUdV;
    Double Kinav;
    Double invKin;
    Double Uen;
    Double Uen2;
    Double Kin2;
    Double Entot;
    Double Entot2;
    Double Kinavat;
    Double invKinat;
    Double Kin2at;
    Double bdUdv;
    Double bdUdvsq;
    Double bd2Udv2;
    Double UdUdvF;
    Double Volnm;
    Double Volsq;
    Double Hconf;
    Double Hconf2;
    Double HV;
    Double Hinter;
    Double Hinter2;
    Double HintV;
    Double P;

} AveragLus; 
    
typedef struct
{
    Double pressure;
    Double shearStress;
    Double antiShearStress;
    Double firstNormStress;
    Double secondNormStress;
} Stress;

/**
 ** Forward declarations
 **/

class Ensemble {
    /**
     ** Data Member - 
     **/
    // molecule and atom components for this simulation system - the topology  
    public:  
    Int     nAtoms;             // total number of atoms (or sites)
    Int     nA3;            
    Int     nMols;              // total number of molecules            
    Int     nBonds;             // total number of bonds            
    Int     nAngles;            // total number of angles
    Int     nDihedrals;         // total number of dihedrals
    Int     nImpropers;         // total number of Impropers

    Int     nAtomTypes;         // number of atom types
    Int     nMolTypes;          // number of molecule types
    Int     nBondTypes;         // number of bond types
    Int     nAngleTypes;        // number of angle types
    Int     nDihedralTypes;     // number of dihedral types
    Int     nImproperTypes;     // number of improper types
    Int     rank, size;         // Jc: for Paralle computation

    AtomStruct* atomStructs;
    MolStruct*  molStructs;
    Bond*       bonds;          // array of Bonds
    Angle*      angles;
    Dihedral*   dihedrals;
    Improper*   impropers;
    Int**       exclusionTable;

    // atom & molecule objects
    // array of atom objects, should use Atom**, array of ptr to Atom objs??
    Atom*       atoms;          
    Molecule*   molecules;      // array of molecule objects

    // Simulation box, boundary, or cell manager data members 
    Double    boxLx, boxLy, boxLz, halfLx, halfLy, halfLz,  cellLen;
    Double    volume, density;
    Double    rBoxLx, rBoxLy, rBoxLz;           // reverse box lenth
    Double    pairCutOff;            // square of cutoff plus buffer size for determining pairs.
    Double    mySeed;

    double totsteps;
    Parameters*         myParams;
    SimConfiguration*   myConfig;
    CellManager*        myCellManager;

    ofstream            *sysdataFile;     // used for debuging
    ofstream            *resultFile;
    ofstream*   ofo;   // JC: added by Jianhui

    Energy myEnergy;
    Lustig myLustig;
    AveragLus myAvLus;
    Double virial[9];
    Double longVirial[9];
    Double molVirial[9];
    Double atomPressure[9];
    Double molPressure[9];
    Double atomTemp;
    Double molTemp;

    // for accumulation
    Double energy;
    double ulus;
    double dudv;
    double d2udv;
    Stress atomStress;
    Stress molStress;

    // for MSD and/or VACF 
    Int         counter;        // number of production timesteps completed
    int		luscount;     
    Int         nSamples;       // total number of samplings or measurements
    Int         interval;       // interval number of timesteps for sampling
    Int         sampleCounter;  // number of sampling completed
    Double      deltaT;         // timestep
    Double      *time;
    Vector3     **trajectory;
    Vector3     **velocity;
    Int         *numPerType;    // number of atoms or mols for each type of species     
    bool        computeMSD, computeVACF;
    bool        useCellPairList;
    double      temperature;
    /**
     ** constructor and destructor
     **/
    public:
    Ensemble();
    ~Ensemble();

    /**
     ** methods
     **/
    public:
    void set_mol_system(Parameters* params);
    void set_configuration(SimConfiguration* config);
    void set_pairlist();
    void compute_bound_box();
    void update_pairlist();
    void read_coordinate(char* fname);
    void read_velocity(char* fname);
    void read_coordinate(FILE* fptr);
    void read_realPos(FILE* fptr);
    void read_velocity(FILE* fptr);
    void random_velocity(unsigned seed);    
    void init_velocity();      // Translational velocitys from Maxwell-Boltzmann distribution
    Double gauss();                         // generate uniform random number from standard normal Gaussian distribution
    Double drand();                         // generate random number between 0 - 1;
    void build_exclusion_table();
    void build_pairlist();          // build pairlist using brute force method
    void set_sysdata_file(ofstream *ofp);
    void set_result_file(ofstream *ofp);

    void zero_sum();
    void accumulate();
    void sampling();            // store sampling data (trajectory and/or velocity) for bunch output

    void write_ensemble_info(ofstream &of);
    void write_sys_data(ofstream &of);
    void write_pairlist(ofstream &of);
    void write_trajectory(ofstream &of, Int numSteps);
    void write_sampling_data();
    void write_result(ofstream &of, Int numSteps);
    void write_resultInduction(ofstream &of, Int numSteps);
    void write_bond(Molecule *mol, ofstream &of);
    void write_Potential(ofstream &of);      //Jc: added by Jianhui
    void write_systemInfo(ofstream &of);     //Jc: created by Jianhui
    void write_pressureTensor(ofstream &of,int counter); //Jc:added by Jianhui
    void write_press_res(ofstream &ofpr, int counter);
    void write_pressac(ofstream &of, int counter);
    void write_lustig(ofstream &of, int counter);
    void write_TPbylustig(ofstream &of, int counter);
    void write_averaglus(ofstream &of, int counter);
    void write_velbehav(ofstream &of, int counter);
    void write_position(int timeStep); // Jc: added byu Jianhui 
    // check_exclusion(): checks if a pair of atoms (atom1 and atom2) should be
    // excluded from non-bonded interactions based on exclusion table
    inline bool exclusion_check(Int atom1, Int atom2)
    {
        Int i = 0;
        Int val; 

        while ((val = exclusionTable[atom1][i++]) != -1)
        {
            if (atom2 == val)
                return true;             // excluded 
        }
        return false;                    
    }

    inline bool exclusion_check(Atom* atom1, Atom* atom2)
    {
        return exclusion_check(atom1->atomID, atom2->atomID);
    }

    inline void box_dimenssion(Double* lx, Double* ly, Double* lz)
    {
        *lx = boxLx;
        *ly = boxLy;
        *lz = boxLz;
    }

    inline void apply_pbc(Vector3 &r)
    {
        /* if (r.x >= halfLx)  r.x -= boxLx;
        if (r.x < -halfLx)  r.x += boxLx;
        if (r.y >= halfLy)  r.y -= boxLy;
        if (r.y < -halfLy)  r.y += boxLy;
        if (r.z >= halfLz)  r.z -= boxLz;
        if (r.z < -halfLz)  r.z += boxLz; */

        while (r.x >= halfLx)  r.x -= boxLx;
        while (r.x < -halfLx)  r.x += boxLx;
        while (r.y >= halfLy)  r.y -= boxLy;
        while (r.y < -halfLy)  r.y += boxLy;
        while (r.z >= halfLz)  r.z -= boxLz;
        while (r.z < -halfLz)  r.z += boxLz;  
    }

    private:
        
};

#endif


/******************************************************************************
    // check_exclusion() checks if a pair of atoms (atom1 and atom2) should be
    // excluded from non-bonded interactions based on exclusion policy specified
    // by exclusionDelta. The function returns exclusion type as following:
    // 0 - none
    // 1 - 1,2 exclusion
    // 2 - 1,3 exclusion
    // 3 - 1,4 exclusion
    inline Int check_exclusion(Atom* atom1, Atom* atom2, Int exclusionDelta)
    {
       Int delta;

       if (atom1->molID != atom2->molID)
          return 0;                     
       if (atom2->atomID > atom1->atomID)
          delta = atom2->atomID - atom1->atomID;
       else 
          delta = atom1->atomID - atom2->atomID;
       if (delta > exclusionDelta)
          return 0;      
       return delta;
    } 

    inline Int check_exclusion(Int atom1, Int atom2, Int exclusionDelta)
    {
        Int delta;

        delta = atom2 > atom1?(atom2 - atom1):(atom1 - atom2);
        if (delta > exclusionDelta)
            return 0;
        if (atoms[atom1].molID != atoms[atom2].molID)
            return 0;
        return delta;
    }
**********************************************************************************/

