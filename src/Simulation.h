/** Simulation.h -- handle input,  
 **     and build the simulation system and invoke the simulation
 **
 ** Copyright (C) 2003
 ** Centre for Molecular Simulation (CMS)
 ** School of Information Technology
 ** Swinburne University of Technology
 ** PO Box 21 Hawthorn, Vic 3122, Australia
 **
 ** Author: Zhongwu Zhou, modified by Jianhui LI for parallel computation
 ** Email: zzhou@it.swin.edu.au, jli@it.swin.edu.au
 **
 ** Edited by
 ** Author: Edoardo Savoia
 ** Email: esavoia@swin.edu.au
 **/

#ifndef SIMULATION_H
#define SIMULATION_H

#include <iostream>
#include <fstream>
#include <stdio.h> //??
#include <vector>

#include "NEMD_defs.h"
#include "Parameters.h"
#include "SimConfiguration.h"
#include "Ensemble.h"

#include "Integrator.h"
#include "integrators/GearIntegrator.h"
#include "integrators/LFIntegrator.h"
#include "integrators/NHIntegrator.h"
#include "integrators/VelocityIntegrator.h"
#include "observables/RDF.h"

#include <mpi.h>   // added by Jianhui // ?? Needed?
#include "utils/Vector3.h"
#include "utils/CellManager.h"


class Simulation {
    private:
    // I/O files
    const char*   sysConfigFile;    // system config data file
    const char*   coordinateFile;         // coordinate data file
    const char*   sysDataFile;            // molecular system parameters and structures file
    const char*   velocityFile;
    const char*   resultFile;
    const char*   trajectoryFile;
    
    const char*   dumpFile;               // have to consider what data will be dumped into it (?)
    const char*   restartFile;            // also used as a backup file ?
    
    FILE*   dumpFilePtr;
    ofstream *ofp;
    ofstream *ofpr;
    ofstream *ofac;
    ofstream *oflu;
    ofstream *oftp;
    ofstream *ofvel;
    ofstream *ofavl;
    ofstream *ofind;
    ofstream *ofpTrajectory;
    ofstream *ofPressureTensor;  //JC added by Jianhui Li
    ofstream *ofColPotential;    //JC added by Jianhui Li

    SimConfiguration  *myConfig;
    Ensemble    *myEnsemble;
    Parameters  *myParams;
    Integrator  *myIntegrator;
    RDF         *rdf;

// paralle setup parameters
    int size,rank;        

    /**
     ** constructor and destructor
     **/
    public:
        Simulation(const char configFile[]);
        ~Simulation();

    /**
     ** methods
     **/
    public:
        void setup_simulation(void);
        void run(void);
        void finish(void);       
};


#endif
