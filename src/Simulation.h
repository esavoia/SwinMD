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

#include "Parameters.h"
#include "SimConfiguration.h"
#include "Ensemble.h"
#include "Integrator.h"
#include "observables/RDF.h"

class Simulation {
private:
    // I/O file paths/ file names
    const char*   sysConfigFile;    	// System configuration file
    const char*   sysDataFile;          // Molecular system parameters and structures file
	const char*   trajectoryFile;		// Trajectory file
	
    // Output streams
	// [of] -> sysData.out
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
    
	// Simulation Configuration
    SimConfiguration  *myConfig;
	// Observables (???)
    Ensemble    *myEnsemble;
    // Molecular System Params
	Parameters  *myParams;
	Integrator  *myIntegrator;
    RDF         *rdf;
    
    
public:
    ~Simulation();
        
private:
	// Private Constructor: in this way to have a Simulation Object you have to use the static build function.
	Simulation(const char configFile[]);
    // This function setup all output files and initialise their headers when needed
	// Returns FALSE in case of FAILURE
	bool setup_outputs();
	// This function setup the simulation
	// Returns FALSE in case of FAILURE
	bool setup_simulation(void);
    
public:
	// run the simulation
	// NOTE: there is no way to escalate an error inside this function to the caller (main), it should return a bool or an err/int flag?
    void run(void);
	// Static build function, it build Simulation object.
	// Returns NULL in case of FAILURE
    static Simulation* build(const char[]);
};

#endif
