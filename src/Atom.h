/** Atom.h -- 
 **
 ** Copyright (C) 2002
 ** Centre for Molecular Simulation (CMS)
 ** School of Information Technology
 ** Swinburne University of Technology
 ** PO Box 21 Hawthorn, Vic 3122, Australia
 **
 ** Author: Zhongwu Zhou
 ** Email: zzhou@it.swin.edu.au
 **/

#ifndef Atom_H
#define Atom_H

#include "NEMD_defs.h"
#include "Parameters.h"
#include "utils/Vector3.h"

/**
 ** Forward declarations
 **/


class Atom {
    /**
     ** Data Member - 
     **/
    // atom info & parameters
    public:
    Int     atomID;
    Int     atomType;
    Int     molID;              // parent molecule
    Double  scaledCharge;       // scaled by the square root of Coulomb constant
    Double  mass;
    char*   typeName;
    double dipx;		// Dipole in the x direction of local frame
    double dipy;		// Dipole in the y direction of local frame
    double dipz;		// Dipole in the z direction of local frame
    double quadxx;		// Quadrupole xx in the local frame
    double quadxy;		// Quadrupole xy in the local frame
    double quadxz;		// Quadrupole xz in the local frame
    double quadyy;		// Quadrupole yy in the local frame
    double quadyz;		// Quadrupole yz in the local frame
    double quadzz;		// Quadrupole zz in the local frame
    double polar;		// Polarizability of atom


    Vector3 position;           // current position
    Vector3 displacement;       // displacement since making of pairlist
    Vector3 momentum;
    Vector3 velocity;
    Vector3 longForce;          // long range electric forces
    Vector3 force;              // total force on this atom

    Vector3 realPos;            // real position of the atom for computing MSD
    // Vector3 acceleration;
    // Vector3 shortForce;         // short range forces

    Atom**  myPairList;         // array of points to Atoms

    private:
    Int     listSize;           // size of the pairList
    Int     currentSize;        // current number of pairs

    /**
     ** constructor and destructor
     **/
    public:
    Atom();  
    ~Atom();


    /**
     ** methods
     **/
    public:    
    Int get_list_size() {return currentSize; }
    void clear_pairlist();
    void set_pair(Atom* a);
    void resize();
};

#endif
