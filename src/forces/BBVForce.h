/** BBVForce.h -- 
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

#ifndef BBVForce_H
#define BBVForce_H

#include "../Force.h" 
// #define SWITCHDIST  0.085

//  Parameters from S. Bock et al. / Chemical Physics 257 (2000) 147-156
#define a1cc   730230.214728524	 // 278.13    		//
#define a1co   36158.3810349162  // 13.772 		//  
#define a1oo   60730.4321608079  //23.131    		//
#define a1bc   0.0    		//
#define a1bo   0.0   		//
#define a1bb   0.0    		//

#define a2cc   41.4246863779215   // 2.1921    		//
#define a2co   33.8676716055431   // 1.7922 		//  
#define a2oo   31.8759002701875   // 1.6868    		//
#define a2bc   0.0    		//
#define a2bo   0.0   		//
#define a2bb   0.0    		//

#define a3cc   -4.56460396927518E-06 // -3.6056E6   		//
#define a3co   -1.54385803764452E-06 // -1.2195E6	//  
#define a3oo   3.17620781520815E-07  // 2.5089E5		//
#define a3bc   -2.06569344815462E-06  // -1.6317E6		//
#define a3bo   -4.86970424567725E-07 // -3.8466E5		//
#define a3bb   4.06315077640024E-06  // 3.2095E6		//

#define a4cc   0.000591828827475  // 1.3091E6    		//
#define a4co   0.000112660410822  // 2.4920E5		//  
#define a4oo   -2.59575554099094E-05 // -5.7417E4  		//
#define a4bc   -2.29579493673654E-05 // -5.0782E4		//
#define a4bo   4.14031530652999E-05 // 9.1582E4		//
#define a4bb   -0.000225112861013  // -4.9794E5		//

#define a5cc   -0.021759378098886  // -1.3478E5   		//
#define a5co   -0.003041921665968  // -1.8842E4	//  
#define a5oo   0.000891265935099   // 5.5206E3		//
#define a5bc   0.004967783228081   // 3.0771E4		//
#define a5bo   -0.000981658246997  // -6.0805E3		//
#define a5bb   0.003566290712304   // 2.2090E4		//

#define a6cc   0.330141744489845   // 5.7264E3    		//
#define a6co   0.020094230864852   // 3.4854E2		//  
#define a6oo   -0.017044409286983  // -2.9564E2  		//
#define a6bc   -0.118447229333331  // -2.0545E3		//
#define a6bo   0.014663934183954  // 2.5435E2		//
#define a6bb   0.003839777253862  // 6.6602E1		//

#define Aa0   0.05291772		// Bohr radius [nm]
#define Aa0met   5.29177210903E-11	// Bohr radius [m]


typedef struct bbvParam
{
        Int atomType1;          //  Index for first atom type
        Int atomType2;          //  Index for second atom type
        Double Aa1;           	//  Parameter a1 for this pair
        Double Aa2;             //  Parameter a2 for this pair
	Double Aa3;		//  Parameter a3 for this pair
        Double Aa4;           	//  Parameter a4 for this pair
        Double Aa5;             //  Parameter a5 for this pair
	Double Aa6;		//  Parameter a6 for this pair
} BBVParam;
         

/**
 ** Forward declarations
 **/


class BBVForce : public Force {
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
    BBVParam*  bbvParams;
    Double cutOff, cutOff2, switchDist;
    Double c1, c2, c3, c4;              // parameters for switch function 
    Double eLrc, vLrc, dUlrc, d2Ulrc, vgromLRC;                  // energy & virial long range correction
    double co2bond, dudvq;
    bool    switchOn;           
    bool   computeCoulomb;

    /**
     ** constructor and destructor
     **/
    public:
    BBVForce(Ensemble* ensemble);
    ~BBVForce()  { ;  }

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
    void set_bbv_param();

};

#endif
