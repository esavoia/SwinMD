/** Molecule.cpp -- 
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

#include "Molecule.h"

Molecule::Molecule()  
{ 
    molID = -1;
    molType = -1;
    numAtoms = 0;
    mass = 0.0;
    molName = NULL;
    myAtoms = NULL;
    currentSize = 0;
}

Molecule::~Molecule()  
{   if (myAtoms != NULL)    delete [] myAtoms;   } 

    /**
     ** methods
     **/

void Molecule::set_molecule(Int id, Int type, Int nAtoms, char* name) 
{
    molID = id;
    molType = type;
    numAtoms = nAtoms;
    molName = name;

    if(numAtoms)
    {
       myAtoms = new Atom* [numAtoms+1];
       for (int i=0; i <= nAtoms; i++)
           myAtoms[i] = NULL;
    }
} 

void Molecule::add_atom(Atom* atom)
{
    if ((myAtoms == NULL) || (currentSize > numAtoms))
    {
       printf("%s", "no memory allocated for adding new atom");
       exit(1);
    }
    myAtoms[currentSize] = atom;
    currentSize++;
}


    


