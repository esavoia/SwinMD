/** Parameters.h -- This head file actually defines data structures for molecule parameters 
 **     and topology of the whole molecular and atom systems. Structs are used for possible 
 **     improvement in both performance and memory usage.
 **
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

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "NEMD_defs.h"
#include <mpi.h>

/** 
 ** infomation and parameters of the molecular systems
 **/
typedef struct atomParam 
{
    Int atomType;                   // ID of this atom type 
    char typeName[4];
    char atomName[4];
    Double mass;                    
    Double charge;
    Double sigma;                   // Sigma value
    Double eps;                     // Epsilon value
    Double dipx;		    // permanent dipole x
    Double dipy;		    // permanent dipole y
    Double dipz;		    // permanent dipole z
    Double quadxx;		    // permanent quadrupole xx
    Double quadxy;		    // permanent quadrupole xy
    Double quadxz;		    // permanent quadrupole xz
    Double quadyy;		    // permanent quadrupole yy
    Double quadyz;		    // permanent quadrupole yz
    Double quadzz;		    // permanent quadrupole zz
    Double polar;		    // polarizability value
    // Double sigma14;              // Sigma value for 1-4 interactions for later use
    // Double eps14;                // Epsilon value for 1-4 interactions
} AtomParam;

typedef struct bondParam
{
    Int bondType;                   // ID of this bond type
    Int atomType1;                  // first type atom forming the bond 
    Int atomType2;                  // second type atom forming the bond       
	Double k;		                // Force constant for the bond
	Double r0;	                    // Rest distance for the bond
} BondParam;

typedef struct angleParam
{
    Int angleType;                  // ID of this angle type 
    Int atomType1;                  // 1st type atom at one end in this angle
    Int atomType2;                  // 2nd type atom at the vertex in this angle
    Int atomType3;                  // 3rd type atom at the other end in this angle
	Double k;               		// Force constant for this angle
	Double theta0;	                // Rest angle for this angle
	// Double k_ub;	                // Urey-Bradley force constant
	// Double r_ub;	                // Urey-Bradley rest length
} AngleParam;

typedef struct torsionParam
{
    Int torsionType;                // ID of this torsion
    Int atomType1;
    Int atomType2;
    Int atomType3;
    Int atomType4;
    Double k;                       // force constant for this torsion
    // Double phaseShift;           // phase shift for this torsion
    Int d;
    Int n;                          // periodicity for this interaction
} TorsionParam;

// ljPairParam is used to search L-J force parameters for a pair of atoms
// from atom type index.
typedef struct ljPairParam
{
   	Int atomType1;		    //  Index for first atom type
   	Int atomType2;		    //  Index for second atom type
   	Double sigma;	    	//  Parameter A (Sigma*Sigma) for this pair
	Double eps;		        //  Parameter B (4*Epsilon) for this pair
    // double A14;          //  Parameter A for 1-4 interactions
	// double B14;		    //  Parameter B for 1-4 interactions
} LJPairParam;

typedef struct molType
{
    Int molType;
    Int numOfAtoms;
    char name[16];
} MolType;


/**
 ** topology of the molecular system
 **/
typedef struct atomstruct
{
        Int atomID;
        Int atomType;
        Int molID;
} AtomStruct;

typedef struct molstruct
{
        Int molID;
        Int molType;
} MolStruct;

typedef struct bond
{
        Int atom1;
        Int atom2;
        Int bondType;
} Bond;

typedef struct angle
{
        Int atom1;
        Int atom2;
        Int atom3;
        Int angleType;
} Angle;

typedef struct dihedral
{
        Int atom1;
        Int atom2;
        Int atom3;
        Int atom4;
        Int dihedralType;
} Dihedral;

typedef struct improper
{
        Int atom1;
        Int atom2;
        Int atom3;
        Int atom4;
        Int improperType;
} Improper;

class Parameters
{
    // parameters
    public:
    Int     nAtomTypes;         // number of atom types
    Int     nMolTypes;          // number of molecule types
    Int     nBondTypes;         // number of bond types
    Int     nAngleTypes;        // number of angle types
    Int     nDihedralTypes;     // number of dihedral types
    Int     nImproperTypes;     // number of improper types
    Int     ljTableSize;          

    AtomParam*      atomParams;            // array of AtomParam
    BondParam*      bondParams;            // array of BondParam
    AngleParam*     angleParams;
    TorsionParam*   dihedralParams;
    TorsionParam*   improperParams; 
    MolType*        molTypes; 
   
    // ljParamTable - linearised triangular array contains L-J params for all atom type pairs.
    // lj-pair-parameters are duplicated to give a symmetric array for easy data retrieve.
    // ljTableSize is actually number of atomTypes invoved in lj pair parameters therefore 
    // the total array size is ljTableSize*ljTableSize.
    LJPairParam*   ljParamTable;  
 
    // molecular system topology
    Int     nAtoms;             // total number of atoms
    Int     nMolecules;         // total number of molecules
    Int     nBonds;             // total number of bonds            
    Int     nAngles;            // total number of angles
    Int     nDihedrals;         // total number of dihedrals
    Int     nImpropers;         // total number of Impropers
    
    AtomStruct* atomArray;      // atom and mol arrays are used for storing read in data 
    MolStruct*  molArray;       // atom and mol objs to be created in Ensemble with more data
    Bond*       bondArray;
    Angle*      angleArray;
    Dihedral*   dihedralArray;
    Improper*   improperArray;
    int rank,size;              // Jc: MPI features
    // constructor & destructor
public:
    Parameters(const char* sysDataFile);
    ~Parameters();

    // read parameter data and structure data
    void read_sys_data(const char* sysDataFile); 

    // for retrieving parameters - given the index of a bond or angle, etc.
    // first retrive its type-index as the entry to the params array
    // and then retrieve relevant parameter values from the params array.
    inline void get_bond_params(Int bondID, Double* k, Double* r0)
    {
        Int bondType = bondArray[bondID].bondType;
        *k  = bondParams[bondType].k;
        *r0 = bondParams[bondType].r0;
    }

    // given ids of two bonded atoms, this method retrieves the bond length
    Double get_bond_length(Int a1, Int a2)
    {
        int atom1, atom2;
        for (Int i = 0; i < nBonds; i++)
        {
            atom1 = bondArray[i].atom1;
            atom2 = bondArray[i].atom2;

            if(((a1 == atom1) && (a2 == atom2)) || ((a1 == atom2) && (a2 == atom1)))
            {
                Int type = bondArray[i].bondType;
                return bondParams[type].r0;
            }
        }
        return -1;             // atoms a1, a2 are not bonded
    }

    // given ids of three angled atoms, this method retrieves the length between atom1 & atom3
    inline Double get_bond_length(Int a1, Int a2, Int a3)
    {
        Double len1, len2, angle;
        len1 = get_bond_length(a1, a2);
        len2 = get_bond_length(a2, a3);

        angle = get_angle(a1, a2, a3);
        return sqrt(len1*len1 + len2*len2 - 2*len1*len2*cos(angle));

    }

    inline Double get_angle(Int a1, Int a2, Int a3)
    {
        Int atom1, atom2, atom3, type;
        for (Int i = 0; i < nAngles; i++)
        {
            atom1 = angleArray[i].atom1;
            atom2 = angleArray[i].atom2;
            atom3 = angleArray[i].atom3;
        
            if(((a1==atom1) && (a2==atom2) && (a3==atom3)) || ((a1==atom3) && (a2==atom2) && (a3==atom1)))
            {
                type = angleArray[i].angleType;
                return angleParams[type].theta0;
            }
        }
        return -1;            // atoms a1, a2 and a3 are not angle
    }

    inline void get_angle_params(Int angleID, Double* k, Double* theta0)
    {
        Int angleType = angleArray[angleID].angleType;
        *k = angleParams[angleType].k;
        *theta0 = angleParams[angleType].theta0;
    }

    inline void get_improper_params(Int improperID, Double* k, Int* d, Int* n)
    {
        Int improperType = improperArray[improperID].improperType;
        *k = improperParams[improperType].k;
        *d = improperParams[improperType].d;
        *n = improperParams[improperType].n;
    }

    inline void get_dihedral_params(Int dihedralID, Double* k, Int* d, Int* n)
    {
        Int dihedralType = dihedralArray[dihedralID].dihedralType;
        *k = dihedralParams[dihedralType].k;
        *d = dihedralParams[dihedralType].d;
        *n = dihedralParams[dihedralType].n;
    }

    // methods for handling LJ-table
    void build_lj_param_table();
    void set_lj_param(LJPairParam *ljTable, const AtomParam &ap1, const AtomParam &ap2); 
    const LJPairParam& get_lj_parameter(Int atomType1, Int atomType2)
    {
        return ljParamTable[atomType1*ljTableSize + atomType2];
    }

private:
    int check_heading(const char* heading, FILE* fptr);

};

#endif
