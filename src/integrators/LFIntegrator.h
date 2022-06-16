/** LFIntegrator.h -- 
 **
 ** Copyright (C) 2002
 ** Centre for Molecular Simulation (CMS)
 ** School of Information Technology
 ** Swinburne University of Technology
 ** PO Box 21 Hawthorn, Vic 3122, Australia
 **
 ** Author: Zhongwu Zhou
 ** Email: zzhou@it.swin.edu.au
 **/
// JC modified by Jianhui Li in 2004 according to the 
#ifndef LFIntegrator_H
#define LFIntegrator_H

#include <iostream>
#include <random>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../Integrator.h"
#include "../utils/Vector3.h"
#include "../Atom.h"
#include <time.h>

/**
 ** Forward declarations
 **/


class LFIntegrator : public Integrator{
    /**
     **  New Data Member - 
     **/
    Double *mass_r;         // reverse atom mass
    Vector3 *vtPlus, *vtMinus;
    Vector3 **rOld, **rNew;
    Vector3 **virt1, **virt2, **virtfor;
    bool    **moving, **moved;
    Double  **bondLength;  
    Int     maxCount;          // only used when SHAKE is applied
    Double  tol;               // tolerance for SHAKE constraint
    double totsteps;
    ofstream *ofPotential;     // Jc: added by jinahui
    ofstream *ofNH;
    ofstream *outTimeFile;     // Jc: added by jinahui
    ofstream *ofHeatCapacity;  // Jc: added by jinahui store the E2Sum and ESum 
    ofstream *ofo;             // Jc: added by jinahui output vtMinus.out
    ofstream *ifo;             // Jc: added by jinahui input vtMinus.out

    int prank,psize,a,a3;      // Jc: added by jinahui
    double E2Sum, ESum, aTSum, mTSum, aT2Sum, mT2Sum;        // Jc: added by Jianhui Li E2Sum =<E2>, ESum = <E>
    double co2bond, llra, Lvir;
   // Double  shakeVirial[9];

    /**
     ** constructor and destructor
     **/
    public:
    LFIntegrator(Ensemble* ensemble, SimConfiguration* simConfig);
    ~LFIntegrator();

    /**
     ** methods
     **/
    public:
    virtual void initialise();
    virtual void run(Int numTimeSteps);
    virtual void VirtualPositions();

    // new function members
    private:
    void thermostat();
    void start();
    // void computePressure();
    void shake(Double tol, Int maxCount, Double *virial, int currentSteps);
    void shakeCFR(Double tol, Int maxCount, Double *virial, int currentSteps);
    void shakeBBV(Double tol, Int maxCount, Double *virial, int currentSteps);
    double urandf();
    double gaussVel(double sigma);
    void write_vtMinus();
    void read_vtMinus();


};

#endif
