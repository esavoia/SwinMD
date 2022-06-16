/** Force.cpp -- 
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

#include "Force.h"

// string Force::forceIdentifier("Force");

Force::Force(Ensemble *ensemble)
{
    myEnsemble = ensemble;
    // myParams = myEnsemble->myParams;
    atoms = myEnsemble->atoms;   
// JC    forceIdentifier = "Force";
    sysdataFile = ensemble->sysdataFile;
};

Force::~Force()
{
    ;
}

//JC string Force::get_force_id() { // in order to comply the cluster compiler
//JC  return Force::forceIdentifier;// in order to comply the cluster compiler
//JC }                              // in order to comply the cluster compiler