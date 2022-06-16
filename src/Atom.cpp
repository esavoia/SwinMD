/** Atom.cpp -- 
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

#include "Atom.h"

Atom::Atom()
{
   Int i;

   atomID        = -1;
   atomType      = -1;
   molID         = -1;
   scaledCharge  = 0.0;
   dipx          = 0.0;              // Dipole in the x direction of local frame
   dipy          = 0.0;              // Dipole in the y direction of local frame
   dipz          = 0.0;              // Dipole in the z direction of local frame
   quadxx        = 0.0;              // Quadrupole xx in the local frame
   quadxy        = 0.0;              // Quadrupole xy in the local frame
   quadxz        = 0.0;              // Quadrupole xz in the local frame
   quadyy        = 0.0;              // Quadrupole yy in the local frame
   quadyz        = 0.0;              // Quadrupole yz in the local frame
   quadzz        = 0.0; 
   polar         = 0.0; 	     // Polarizability of atom
   listSize      = 20;    // take it & size increament when resize() called as input?? 
   currentSize   = 0;
   mass = 0.0;

   myPairList    = new Atom*[listSize];
   for (i = 0; i < listSize; i++)
      myPairList[i] = NULL;
}

Atom::~Atom()
{
    if (myPairList != NULL)   delete [] myPairList;
}

void Atom::set_pair(Atom* a)
{
   if(currentSize >= listSize)        
      resize();

   myPairList[currentSize] = a;
   currentSize++;
}

void Atom::resize()
{
   Int i;
   Atom** tempList;

   listSize += 10;
   tempList = new Atom*[listSize];
   if(tempList == NULL)
   {
      printf("%s", "resize neighbour list error");
      exit(1);
   }

   for ( i = 0; i < currentSize; i++)
      tempList[i] = myPairList[i];
   for (i = currentSize; i < listSize; i++)
      tempList[i] = NULL;
   
   delete [] myPairList;
   myPairList = tempList;
}

// clear & reset pair list
void Atom::clear_pairlist()
{
    if (myPairList != NULL)
        delete [] myPairList;

    myPairList = new Atom*[20];
    listSize = 20;
    currentSize = 0;
    for (int i = 0; i < listSize; i++)
        myPairList[i] = NULL;
}
