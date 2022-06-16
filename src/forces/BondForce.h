/** BondForce.h -- 
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

#ifndef BondForce_H
#define BondForce_H

#include "../Force.h" 

/**
 ** Forward declarations
 **/


class BondForce : public Force {
    /**
     ** New Data Member - 
     **/
    Int numBonds;
    Bond* bonds;
    BondParam* params;
    bool use_harmonic, use_amoeba;

    /**
     ** constructor and destructor
     **/
    public:
    BondForce(Ensemble* ensemble);
    ~BondForce();

    /**
     ** methods
     **/
    public:
    virtual void compute();
//JC    virtual string get_force_id();// in order to comply the cluster compiler
    virtual void write_force_info(ofstream& of);
    virtual void write_energy(ofstream& of);
};

#endif
