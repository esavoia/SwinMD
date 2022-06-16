/** Force.h -- 
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

#ifndef Force_H
#define Force_H

#include <iostream>
#include <fstream>
#include <string.h>
#include "NEMD_defs.h"
#include "Ensemble.h"
#include "Parameters.h"
#include "Atom.h"
#include "Molecule.h"

/**
 ** Forward declarations
 **/


class Force {
    /**
     ** Data Member - 
     **/
    protected:
    Ensemble    *myEnsemble;
    // Parameters  *myParams;
    Atom        *atoms;
    Double      energy;
    double ulus;
    double dudv;
    double d2udv;
    Double      virial[9];
    ofstream    *sysdataFile;       //  used for out put force info

    public:
//JC    string forceIdentifier; // in order to comply the cluster compiler

    /**
     ** constructor and destructor
     **/
    public:
    Force(Ensemble *ensemble);
    virtual ~Force();

    /**
     ** methods
     **/
    public:
    virtual void compute()=0;
//JC    virtual string get_force_id() =0; // in order to comply the cluster compiler
    // for debug use
    virtual void write_force_info(ofstream &of) =0;
    virtual void write_energy(ofstream &of) =0;

    Double get_energy() { return energy; }
    void get_virial(Double* vir)
    {
        for(Int i = XX; i <= ZZ; i++)
            vir[i] = virial[i];
    }
    
};

#endif
