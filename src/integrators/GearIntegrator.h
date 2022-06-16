/** GearIntegrator.h -- 
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

#ifndef GearIntegrator_H
#define GearIntegrator_H

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../Integrator.h"
#include "../utils/Vector3.h"
#include "../Atom.h"

/**
 ** Forward declarations
 **/

// Gear corrector cofficients for a second-order equation
static const Double     K20 = 19.0/120.0;      // or 19.0/90.0
static const Double     K21 = 3.0/4.0;
static const Double     K22 = 1.0;
static const Double     K23 = 0.5;
static const Double     K24 = 1.0/12.0;
// static const Double     TFB = 0.1;              // feed back for thermostat multiplier

// Gear corrector cofficients for a first-order equation
static const Double     K10 = 0.3486111111;   // 251.0/720.0;
static const Double     K11 = 1.0;
static const Double     K12 = 0.9166666667;   // 11.0/12.0;
static const Double     K13 = 0.3333333333;   // 1.0/3.0;
static const Double     K14 = 0.0416666667;   // 1.0/24.0;

class GearIntegrator : public Integrator{
    /**
     **  New Data Member - 
     **/
    Vector3* dr1;        // velocity, the first derivatives
    Vector3* dr2;        // accelerations, the second derivatives - dr2/dt2
    Vector3* dr3;        // the third derivatives - dr3/dt3
    Vector3* dr4;        // the forth derivatives - dr4/dt4
    // Vector3* position;   // previous positions of all atoms

    Vector3* pa1;       // first momenta derivatives  
    Vector3* pa2;       // second momenta derivatives
    Vector3* pa3;       // third momenta derivatives
    Vector3* pa4;       // forth momenta derivatives

    Vector3* pThmo;     // thermostat contribution to pa1
    Vector3* velfb;     // velocity feed back
    Vector3* momfb;     // momentum feed back

    bool init_derivatives;

    ConstraintForce* constraintForce;       // particularly Gaussan constraint is used with Gear Integrator

    /* Double temperature, TFB;
    Double drMax;       // maximum displacement for updating pair list
    Int numUpdates;     // accumulated number of pair list updates
    */

    /**
     ** constructor and destructor
     **/
    public:
    GearIntegrator(Ensemble* ensemble, SimConfiguration* simConfig);
    ~GearIntegrator();

    /**
     ** methods
     **/
    public:
    virtual void initialise();
    virtual void run(Int numTimeSteps);
    virtual void write_state();
    // new function members
    private:
    void predict();
    void correct();
    // void computePressure();
    void gaussThermostat(Vector3 *pThmo, bool atomThermo);
    void start();       // for new start, we initialise first derivatives for predictor-corrector
    void read_derivative(FILE* fptr);

};

#endif
