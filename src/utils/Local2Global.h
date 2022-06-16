/** Local2Global.h -- 
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

#ifndef Local2Global_H
#define Local2Global_H

#include <vector>
#include <cmath>
#include "../Ensemble.h"
#include "../Atom.h"
#include "../Molecule.h"
#include "Vector3.h"
#include "CellManager.h"

using namespace std;

class Local2Global {
    /**
    ** Data Member
    **/
    private:

    int numAtoms, counter,numMols, nAtomsx3;
    Molecule* myMols;
    Ensemble* myEnsemble;
    Atom        *atoms;

    Vector3 **pp1, **pp2;
    Vector3 **uu, **vv, **ww;

    public:
    Vector3 *dip;
    double *Qxx, *Qxy, *Qxz, *Qyy, *Qyz, *Qzz;

    /**
     ** constructor and destructor
     **/
    public:
    Local2Global(Ensemble* ensemble);
    ~Local2Global();

    /**
     ** methods
     **/
    public:
    virtual void FrameofReference();


};

#endif
