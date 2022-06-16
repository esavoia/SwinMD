/** BuckinghamForce.h -- 
 **
 ** Copyright (C) 2020
 ** Centre for Molecular Simulation (CMS)
 ** School of Information Technology
 ** Swinburne University of Technology
 ** PO Box 21 Hawthorn, Vic 3122, Australia
 **
 ** Author: 
 ** Email: 
 **/

#ifndef BuckinghamForce_H
#define BuckinghamForce_H

#include "../Force.h" 
// #define SWITCHDIST  0.085

//  Parameter from model 3 in  J. Phys. Chem., Vol. 100, No. 11, 1996
#define AAcc   1.489504E-03 // 356.0 kcal A^6/mol    		//
#define AAco   1.556448E-03 // 372.0 kcal A^6/mol  
#define AAoo   1.627576E-03 // 389.0 kcal A^6/mol    		//
#define BBcc   907.928	    // 217.0 kcal/mol
#define BBco   36651.84     // 8760.0 kcal/mol   		//
#define BBoo   1481136      // 354000.0 kcal/mol
#define CCcc   22.7 	    // 2.27 A^-1 
#define CCco   31.6 	    // 3.16 A^-1    		//
#define CCoo   44.0 	    // 4.40 A^-1

typedef struct buckinghamParam
{
        Int atomType1;          //  Index for first atom type
        Int atomType2;          //  Index for second atom type
        Double AA;           	//  Parameter A for this pair
        Double BB;              //  Parameter B for this pair
	Double CC;		//  Parameter C for this pair
} BuckinghamParam;
         

/**
 ** Forward declarations
 **/


class BuckinghamForce : public Force {
    /**
     ** My Data Member - 
     **/
    Int numAtoms;
    Int numAtomTypes;
    Int numMols;
    Vector3 **virt1;
    // LJPairParam* ljTable;
    Parameters* params;
    Molecule* myMols;
    BuckinghamParam*  buckParams;
    Double cutOff, cutOff2, switchDist;
    Double c1, c2, c3, c4;              // parameters for switch function 
    Double eLrc, vLrc, dUlrc, d2Ulrc, vgromLRC;                  // energy & virial long range correction
    double co2bond;
    bool    switchOn;           
    bool   computeCoulomb;

    /**
     ** constructor and destructor
     **/
    public:
    BuckinghamForce(Ensemble* ensemble);
    ~BuckinghamForce()  { ;  }

    /**
     ** methods
     **/
    public:
    virtual void compute();     // from class Force
// JC    virtual string get_force_id();
    virtual void write_force_info(ofstream& of);
    virtual void write_energy(ofstream& of);
    virtual void VirtualPositions();

    // new method
    void long_range_correct();
    void set_buck_param();

};

#endif
