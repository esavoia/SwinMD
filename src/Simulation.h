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
    // I/O file paths/ file names
    const char*   sysConfigFile;    // system config data file
    const char*   coordinateFile;         // coordinate data file
    const char*   sysDataFile;            // molecular system parameters and structures file
    const char*   velocityFile;
    const char*   resultFile;
    const char*   trajectoryFile;

    const char*   dumpFile;               // have to consider what data will be dumped into it (?)
    const char*   restartFile;            // also used as a backup file ?
	// ofp -> result.out
	// ofvel -> velbehav.out
	// ofpr -> pressureResult.out
	// ofac -> pressureAcum.out
	// oflu -> Lustig.out
	// oftp -> ThermoProp.out
	// ofavl -> LustigAverages.out
	// ofind -> resultInduction.out
	// Output streams
    FILE*   dumpFilePtr;
    ofstream *ofp; // -> result.out
    ofstream *ofpr; // -> pressureResult.out
    ofstream *ofac; // -> pressureAcum.out
    ofstream *oflu; // -> Lustig.out
    ofstream *oftp; // -> ThermoProp.out
    ofstream *ofvel; // -> velbehav.out
    ofstream *ofavl; // -> LustigAverages.out
    ofstream *ofind; // -> resultInduction.out
    ofstream *ofpTrajectory;
    ofstream *ofPressureTensor;  // -> pTensor.out // JC added by Jianhui Li
	/* NEVER USED */
	// ofstream *ofColPotential;    //JC added by Jianhui Li

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
		// static unique_ptr<Simulation> build(const char configFile[]) {
		// 	unique_ptr<Simulation> a = make_unique(Simulation(configFile));
		// 	// return std::make_unique<Simulation>(configFile);
		// 	return a;
		// }
		static Simulation* build(const char[]);
};


#endif
