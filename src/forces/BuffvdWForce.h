/** BuffvdWForce.h -- 
 **
 ** Copyright (C) 2003
 ** Centre for Molecular Simulation (CMS)
 ** School of Information Technology
 ** Swinburne University of Technology
 ** PO Box 21 Hawthorn, Vic 3122, Australia
 **
 ** Author: Elton Oyarzua
 ** Email: eoyarzua@swin.edu.au
 ** Based on Force.h by Zhongwu Zhou
 **/

#ifndef BuffvdWForce_H
#define BuffvdWForce_H

#include "../Force.h" 
// #define SWITCHDIST  0.085

/**
 ** Forward declarations
 **/

typedef struct bufftable
{
        int atomType1;
        int atomType2;
        double sigma;
        double eps;
} BuffTable;


class BuffvdWForce : public Force {
    /**
     ** My Data Member - 
     **/
    Int numAtoms;
    Int numMols;
    // LJPairParam* ljTable;
    Parameters* params;
    Atom*     myAtoms;
    Double cutOff, cutOff2, switchDist, switchd1;
    double c0, c1, c2, c3, c4, c5, denom;              // parameters for switch function 
    Double eLrc, vLrc, dUlrc, d2Ulrc;                  // energy & virial long range correction
    bool    switchOn, useConstraint;           
    bool   computeCoulomb;
    int statEnsem;

    private:
    Molecule* myMols;
    int prank, psize;
    int numTypes;
    double factor, n, m, delta, gamma;
    Vector3 *virt;
    BuffTable* buffpairs;

    /**
     ** constructor and destructor
     **/
    public:
    BuffvdWForce(Ensemble* ensemble);
    ~BuffvdWForce()  { ;  }

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
    void VirtualPositions();
    void set_buff_param();
    void potFunc(double &r, double &sigma, double &eps, double &tmpE, double &tmpf);

};

#endif
