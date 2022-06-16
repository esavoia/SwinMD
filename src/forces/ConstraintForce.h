/** ConstraintForce.h -- 
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

#ifndef ConstraintForce_H
#define ConstraintForce_H

#include "../Force.h"
// #include 

#define TINY  1.0e-20
/**
 ** Forward declarations
 **/


class ConstraintForce : public Force {
    /**
     ** Data Member - 
     **/
    Int nAtoms;
    Int nMolecules;
    Int ***matrixM;             // selector matrix
    Double ***matrixL;          // different matrix from the selector matrix
    Double **bondLength;
    Molecule *molecules;
    Parameters *params;

    Vector3 **rm;               // distance difference matrix
    Vector3 **vm;               // veclocity difference matrix
    Vector3 **fm;               // force difference matrix
    Vector3 *velfb;             // velocity feed back for bonded atoms
    Vector3 *momfb;             // momentum feed back
    Double ***matrixC;          // coefficient matrix for linear equation set
    Double **vecB;              // right hand vector for linear equation set
    Int    *index;
    Double d;                   // index & d are used if the complete LU solver is used

    Double  CFB, DFB;           // feed back factors for bond velocty and momentum
    int counter;

    /**
     ** constructor and destructor
     **/
    public:
    ConstraintForce(Ensemble* ensemble, Double cfb, Double dfb);
    ~ConstraintForce();

    /**
     ** methods
     **/
    public:
    virtual void compute();
//JC    virtual string get_force_id(); // in order to comply with the cluster compiler
    virtual void write_force_info(ofstream& of);
    virtual void write_energy(ofstream& of);

    // new methods
    void init();
    void set_matrix();
    Vector3* get_velocity_feedback();
    Vector3* get_momentum_feedback();
    void lu_decomposite(Double** mat, Int size);
    void lu_decomposite(Double** mat, Int* indx, Double* d, Int size);
    void lu_solver(Double** mat, Double* vec, Int size);
    void lu_solver(Double** mat, Double* vec, Int* indx, Int size);

};


#endif
