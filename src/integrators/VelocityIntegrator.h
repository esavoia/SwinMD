/** VelocityIntegrator.h -- 
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

#ifndef VelocityIntegrator_H
#define VelocityIntegrator_H

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


class VelocityIntegrator : public Integrator{
    /**
     **  New Data Member - 
     **/
    Double *mass_r;         // reverse atom mass

    /**
     ** constructor and destructor
     **/
    public:
    VelocityIntegrator(Ensemble* ensemble, SimConfiguration* simConfig);
    virtual ~VelocityIntegrator();

    /**
     ** methods
     **/
    public:
    virtual void initialise();
    virtual void run(Int numTimeSteps);
    // virtual     void write_state();
    // new function members
    private:
    void thermostat();
    void start();
    // void computePressure();


};

#endif
