/** UreyBradleyForce.h -- 
 **
 ** Copyright (C) 2021
 ** Centre for Molecular Simulation (CMS)
 ** School of Information Technology
 ** Swinburne University of Technology
 ** PO Box 21 Hawthorn, Vic 3122, Australia
 **
 ** Author: Elton Oyarzua
 ** Email: eoyarzua@swin.edu.au
 **/

#ifndef UreyBradleyForce_H
#define UreyBradleyForce_H

#include "../Force.h" 

/**
 ** Forward declarations
 **/


class UreyBradleyForce : public Force {
    /**
     ** New Data Member - 
     **/
    Int numBonds;
    Bond* bonds;
    BondParam* params;
    bool use_harmonic, use_amoeba;

    private:
    int numAtoms, numMols;
    Molecule* myMols;
    double volume, kl, l0;

    /**
     ** constructor and destructor
     **/
    public:
    UreyBradleyForce(Ensemble* ensemble);
    ~UreyBradleyForce();

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
