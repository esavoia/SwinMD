/** MieForce.h -- 
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

#ifndef MieForce_H
#define MieForce_H

#include "../Force.h" 
#include <cmath>
// #define SWITCHDIST  0.085

/**
 ** Forward declarations
 **/


class MieForce : public Force {
    /**
     ** My Data Member - 
     **/
    Int numAtoms;
    Int numMols;
    // LJPairParam* ljTable;
    Parameters* params;
    Atom*     myAtoms;
    Double cutOff, cutOff2, switchDist;
    double nMie, mMie, Cn;
    Double c1, c2, c3, c4;              // parameters for switch function 
    Double eLrc, vLrc, dUlrc, d2Ulrc;                  // energy & virial long range correction
    bool    switchOn;           
    bool   computeCoulomb;

    /**
     ** constructor and destructor
     **/
    public:
    MieForce(Ensemble* ensemble);
    ~MieForce()  { ;  }

    /**
     ** methods
     **/
    public:
    virtual void compute();     // from class Force
// JC    virtual string get_force_id();
    virtual void write_force_info(ofstream& of);
    virtual void write_energy(ofstream& of);

    // new method
    void long_range_correct();

};

#endif
