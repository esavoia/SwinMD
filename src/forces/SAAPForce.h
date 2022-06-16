/** SAAPForce.h -- 
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

#ifndef SAAPForce_H
#define SAAPForce_H

#include "../Force.h" 
// #define SWITCHDIST  0.085

//  Parameters from Re-Parametrization by Ulrich Deiters   /
#define sa0cc   2647.3278943132   		//
#define sa0co   2955.40298345344  // 13.772 		//  
#define sa0oo   434.785828619992  //23.131    		//
#define sa0bc   330.977327095528   		//
#define sa0bo   -156.959014370464   		//
#define sa0bb   210.145778485296    		//

#define sa1cc   -12.87744	 // 278.13    		//
#define sa1co   -26.86706  // 13.772 		//  
#define sa1oo   -1.083306  //23.131    		//
#define sa1bc   -0.1054953	//
#define sa1bo   -2.51832		//
#define sa1bb   -10.6224		//

#define sa2cc   -38.0847515837792   // 2.1921    		//
#define sa2co   -1239.54146436504  // 1.7922 		//  
#define sa2oo   -1104.83723696272  // 1.6868    		//
#define sa2bc   -1962.70522194736 //
#define sa2bo   -1089.13611401264		//
#define sa2bb   116.465598779112	//

#define sa3cc   -2.088866  // -3.6056E6   		//
#define sa3co   -3.633782   // -1.2195E6	//  
#define sa3oo   -2.444001  // 2.5089E5		//
#define sa3bc   -0.05936309  // -1.6317E6		//
#define sa3bo   -2.656448  // -3.8466E5		//
#define sa3bb   -2.701114  // 3.2095E6		//

#define sa4cc   92.170181327152  // 1.3091E6    		//
#define sa4co   397.636343143056  // 2.4920E5		//  
#define sa4oo   -880.71543031112 // -5.7417E4  		//
#define sa4bc   -539.110661289976 // -5.0782E4		//
#define sa4bo   1777.13281502344 // 9.1582E4		//
#define sa4bb   14.6959983201344  // -4.9794E5		//

#define sa5cc   282.1627  // -1.3478E5   		//
#define sa5co   20914.05  // -1.8842E4	//  
#define sa5oo   51233.58   // 5.5206E3		//
#define sa5bc   18413.18   // 3.0771E4		//
#define sa5bo   123054.9  // -6.0805E3		//
#define sa5bb   3138.018   // 2.2090E4		//

#define sa6cc   0.0   // 5.7264E3    		//
#define sa6co   0.0   // 3.4854E2		//  
#define sa6oo   0.0  // -2.9564E2  		//
#define sa6bc   0.0  // -2.0545E3		//
#define sa6bo   0.0  // 2.5435E2		//
#define sa6bb   0.0  // 6.6602E1		//

typedef struct saapParam
{
        Int atomType1;          //  Index for first atom type
        Int atomType2;          //  Index for second atom type
        Double Aa0;           	//  Parameter a0 for this pair
        Double Aa1;           	//  Parameter a1 for this pair
        Double Aa2;             //  Parameter a2 for this pair
	Double Aa3;		//  Parameter a3 for this pair
        Double Aa4;           	//  Parameter a4 for this pair
        Double Aa5;             //  Parameter a5 for this pair
	Double Aa6;		//  Parameter a6 for this pair
} SAAPParam;
         

/**
 ** Forward declarations
 **/


class SAAPForce : public Force {
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
    SAAPParam*  saapParams;
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
    SAAPForce(Ensemble* ensemble);
    ~SAAPForce()  { ;  }

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
    void set_saap_param();

};

#endif
