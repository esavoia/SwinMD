/** Molecule.h -- 
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

#ifndef MOLECULE_H
#define MOLECULE_H

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "Atom.h"
#include "NEMD_defs.h"
#include "Parameters.h"
#include "utils/Vector3.h"

class Molecule {

    /**
     ** Data Member - 
     **/
    // molecular info
    public:
    Int     molID;
    Int     molType;
    char*   molName;    
    Double  mass;           // Mass of whole molecule
    Vector3 massCenter;     // Mass center positions
    // Vector3 molCenter;      // geometry center
    Vector3 momenta;        // Molecular momenta
    Vector3 force;

    // component info
    Int   numAtoms;         // number of atoms
    Int   currentSize;      // current number of atoms in the array of myAtoms
    Atom**  myAtoms;          // array of atoms(ID) owned by this molecule

    /**
     ** constructor and destructor
     **/
    public:
    Molecule();
    ~Molecule(); 

    /**
     ** methods
     **/
    public:
    void set_molecule(Int id, Int type, Int nAtoms, char* name);
    void add_atom(Atom* atom);

    // inline Vector3& get_mol_center() { return molCenter; }

    inline void mass_center()
    {
        massCenter = 0.0;
        for (Int i = 0; i < numAtoms; i++)
            massCenter += myAtoms[i]->mass*myAtoms[i]->position;
        massCenter /= mass;
    }

    /* inline void mol_center()
    {
        molCenter = 0.0;
        for (Int i = 0; i < numAtoms; i++)
            molCenter += myAtoms[i]->position;
        molCenter /= numAtoms;
    } */
            
    inline void mol_momenta()
    {
        momenta.x = momenta.y = momenta.z = 0;
        for (Int i = 0; i < numAtoms; i++)
            momenta += myAtoms[i]->momentum;
    }

    inline void mol_force()
    {
        force.x = force.y = force.z = 0;
        for (Int i = 0; i < numAtoms; i++)
            force += myAtoms[i]->force;
    }
            
};

#endif
