/** MCYForce.h -- 
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

#ifndef MCYForce_H
#define MCYForce_H

#include "../Force.h" 
// #define SWITCHDIST  0.085

// Parameters from J. Chem. Phys., Vol. 64, No.4, 15 February 1976 in a.u.
// or J. Chem. Phys., Vol. 64, No.6. 15 March 1976 in kcal/nm and A
#define a1mcy  4553131.00302996 //  kJ/mol // 1734.196
#define b1mcy  51.5270866511259 //  1/nm   // 2.726696
#define a2mcy  2787.98395418654 //  kJ/mol // 1.061887
#define b2mcy  27.6084262492513 //  1/nm   // 1.460975
#define a3mcy  6089.57077675919 //  kJ/mol // 2.319395
#define b3mcy  29.6189436677632 //  1/nm   // 1.567367
#define a4mcy  1144.73360341454 //  kJ/mol // 0.436006
#define b4mcy  22.3326321627374 //  1/nm   // 1.181792  

typedef struct mcypotParam
{
        Int atomType1;          //  Index for first atom type
        Int atomType2;          //  Index for second atom type
        Double AA1;              //  Parameter a1 for this pair
        Double BB1;              //  Parameter b1 for this pair
        Double AA2;              //  Parameter a2 for this pair
        Double BB2;              //  Parameter b2 for this pair
} MCYParam;


/**
 ** Forward declarations
 **/


class MCYForce : public Force {
    /**
     ** My Data Member - 
     **/
    Int numAtoms;
    Int numAtomTypes;
    Int numMols;
    // LJPairParam* ljTable;
    Parameters* params;
    Atom*     myAtoms;
    Molecule* myMols;
    MCYParam*  mcyParams;
    Double cutOff, cutOff2, switchDist;
    Double c1, c2, c3, c4;              // parameters for switch function 
    Double eLrc, vLrc, dUlrc, d2Ulrc;                  // energy & virial long range correction
    bool    switchOn;           
    bool   computeCoulomb;

    /**
     ** constructor and destructor
     **/
    public:
    MCYForce(Ensemble* ensemble);
    ~MCYForce()  { ;  }

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
    void set_mcy_param();

};

#endif
