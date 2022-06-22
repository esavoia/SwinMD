/** Simulation.c -- implementation to handle input, maintain the global data structure
 **     and build the simulation system
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
// JC modified & commented by Jianhui Li, adding output file
#include "Simulation.h"
// #include <mpi.h> // added by Jianhui
#include "utils/Mpi_V.h" // added by E.S.
#include "Errors.h"
#include "Defaults.h"

#include "integrators/GearIntegrator.h"
#include "integrators/LFIntegrator.h"
#include "integrators/NHIntegrator.h"
#include "integrators/VelocityIntegrator.h"



// constructor - initialisation
Simulation::Simulation(const char configFile[]){
	// MPI params (size & rank) initilization
	// can't we use mpi_size_g and mpi_rank_g? [and rename them as mpi_size,mpi_rank]
	// MPI_Comm_size(MPI_COMM_WORLD, &size);
	// MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	sysConfigFile = configFile;
	// coordinateFile = NULL;
	sysDataFile = NULL;         // name of file of parameters and structures for molecular system
	// velocityFile = NULL;
	// resultFile = NULL;
	trajectoryFile = NULL;
	// dumpFile = NULL;            // have to consider what data will be dumped into it
	// restartFile = NULL;         // also used as to backup file
	
	// Moved to setup stage
	// TODO: error checking
	//    ofp = new ofstream(DEFAULT_result_FNAME, ios::out);
	//    ofvel = new ofstream(DEFAULT_velbehav_FNAME, ios::out);
	//    ofpr = new ofstream(DEFAULT_pressureResult_FNAME, ios::out);
	//    ofac = new ofstream(DEFAULT_pressureAcum_FNAME, ios::out);
	//    oflu = new ofstream(DEFAULT_Lustig_FNAME, ios::out);
	//    oftp = new ofstream(DEFAULT_ThermoProp_FNAME, ios::out);
	//    ofavl = new ofstream(DEFAULT_LustigAverages_FNAME, ios::out);
	//    ofind = new ofstream(DEFAULT_resultInduction_FNAME, ios::out);
	//    ofpTrajectory = NULL;
	//    ofPressureTensor = new ofstream(DEFAULT_pTensor_FNAME,ios::out);      // JC initialise the ofp pointer;
	//
	//	myConfig = new SimConfiguration(sysConfigFile);
	//	myEnsemble = new Ensemble();
	//    myParams = NULL;
	//    myIntegrator = NULL;
	//    rdf = NULL;
	
	ofp = NULL;
	ofvel = NULL;
	ofpr = NULL;
	ofac = NULL;
	oflu = NULL;
	oftp = NULL;
	ofavl = NULL;
	ofind = NULL;
	ofpTrajectory = NULL;
	ofPressureTensor = NULL;      // JC initialise the ofp pointer;
	
	myConfig = NULL;
	myEnsemble = NULL;
	myParams = NULL;
	myIntegrator = NULL;
	rdf = NULL;
}

Simulation::~Simulation(){
	// TODO: proper destructor
	;
}

Simulation* Simulation::build(const char configFile[]){
	Simulation *s = new Simulation(configFile);
	return (s->setup_simulation()) ? s : NULL;
}

// TODO: Proper headers with '\t' separators
bool Simulation::setup_outputs(){
	
	ofp = new ofstream(DEFAULT_result_FNAME, ios::out);
	if(!ofp){
		LOG_ERR("Can't open "<< DEFAULT_result_FNAME)
		return false;
	}
	// Header
	if(mpi_rank_g == 0){
		*ofp << "Energy" << endl;
		*ofp << "Steps" << '\t' << " L-J" << '\t' << " Angle" << '\t' << " Bond" << '\t' << " UB-hydrogens "<< "Real" << '\t';
		*ofp << "  Long" << '\t' << "  Correct" << '\t' << "  molCorrect" << '\t' << "  surfCorrect" << '\t';
		*ofp << "  AtomTemp" << '\t' << "  MolTemp" << '\t' << "  Total" << endl;
	}
	
	ofvel = new ofstream(DEFAULT_velbehav_FNAME, ios::out);
	if(!ofvel){
		// Error!
		LOG_ERR("Can't open "<< DEFAULT_velbehav_FNAME)
		return false;
	}
	// Header
	if(mpi_rank_g == 0){
		*ofvel << "# Step  " <<" velcomX  " << "  velcomY  "<< "  velcomZ  " <<"   Diff_vel  " <<"  Carbon_vel  "<< "   Density  ";
		*ofvel <<  " C-O0_dist  " <<  "  stdC-O0   "<<  "   C-O1_dist  " << "   stdC-O1  " << "   Angle  " << "    anglSTDdev  " <<endl;
	}
	
	ofpr = new ofstream(DEFAULT_pressureResult_FNAME, ios::out);
	if(!ofpr){
		// Error!
		LOG_ERR("Can't open "<< DEFAULT_pressureResult_FNAME)
		return false;
	}
	// Header
	if(mpi_rank_g == 0){
		*ofpr << "## Step  " <<"Atom_Pdiag  " << "shear_stress  "<< "antisymm_SS  " <<"first_normal_Sh  " <<"second_normal_Sh  "<< "Mol_Pdiag  ";
		*ofpr << "shear_stress  "<< "antisymm_SS  " <<"first_normal_Sh  " <<"second_normal_Sh  " << endl;
	}
	
	ofac = new ofstream(DEFAULT_pressureAcum_FNAME, ios::out);
	if(!ofac){
		// Error!
		LOG_ERR("Can't open "<< DEFAULT_pressureAcum_FNAME)
		return false;
	}
	// Header
	if(mpi_rank_g == 0){
		*ofac << "## Step  " <<"Atom_Pdiag  " << "shear_stress  "<< "antisymm_SS  " <<"first_normal_Sh  " <<"second_normal_Sh  "<< "Mol_Pdiag  ";
		*ofac << "shear_stress  "<< "antisymm_SS  " <<"first_normal_Sh  " <<"second_normal_Sh  " << endl;
	}
	
	oflu = new ofstream(DEFAULT_Lustig_FNAME, ios::out);
	if(!oflu){
		// Error!
		LOG_ERR("Can't open "<< DEFAULT_Lustig_FNAME)
		return false;
	}
	// Header
	if(mpi_rank_g == 0){
		*oflu << "# Step  " << "   lusStep    " << " Atomtemp  " << "  MolTemp  " << "  U_LJ    " << "  dUdV_LJ  " << "  d2UdV2_LJ   ";
		*oflu << "      -bdUdV    " << "    (-bdUdV)2  "<< "      -bd2UdV2   " << "     UdUdVF      "  <<"  real2ndDev  " <<"  long2ndDev  "<< endl;
	}
	
	oftp = new ofstream(DEFAULT_ThermoProp_FNAME, ios::out);
	if(!oftp){
		// Error!
		LOG_ERR("Can't open "<< DEFAULT_ThermoProp_FNAME)
		return false;
	}
	// Header
	if(mpi_rank_g == 0){
		*oftp << "## Step     " << " U (kJ/mol) "  << " Cv (kJ/molK)  " << "   P (Mpa)" << endl;
	}
	
	ofavl = new ofstream(DEFAULT_LustigAverages_FNAME, ios::out);
	if(!ofavl){
		// Error!
		LOG_ERR("Can't open "<< DEFAULT_LustigAverages_FNAME)
		return false;
	}
	// Header
	if(mpi_rank_g == 0){
		*ofavl << "# Step  " << "   count   "<<"     bUconf/F  " << "    (bUconf/F)2  " << "    bavU/F  " << "     (bavU/F)2  " << "   Kinav  " << "  1/Kin   ";
		*ofavl << "      Uen  " <<"      Uen2   " << "     Kin2  "<< "      Etot    " << "      Etot2   " << "      Kinatom   "<< "       1/Kinatom   " << "      Kinatom2  " <<endl;
	}
	
	ofind = new ofstream(DEFAULT_resultInduction_FNAME, ios::out);
	if(!ofind){
		// Error!
		LOG_ERR("Can't open "<< DEFAULT_resultInduction_FNAME)
		return false;
	}
	// Header
	if(mpi_rank_g == 0){
		*ofind << "## Step     " << "   Real "  << "   Long  " << "   SelfCorrect  " << "   MolCorrect " << "  SurfacePol   " << "    TotInduction  "<< endl;
	}
	// ofpTrajectory = NULL;
	ofPressureTensor = new ofstream(DEFAULT_pTensor_FNAME,ios::out);      // JC initialise the ofp pointer;
	if(!ofPressureTensor){
		// Error!
		LOG_ERR("Can't open "<< DEFAULT_pTensor_FNAME)
		return false;
	}
	return true;
}

bool Simulation::setup_simulation() {
	if(!setup_outputs()){
		return false;
	}
	// Initialise sysData output
	ofstream *of = new ofstream(DEFAULT_sysData_FNAME , ios::out);
	if(!of){
		// Error!
		LOG_ERR("Can't open "<< DEFAULT_sysData_FNAME)
		return false;
	}
	
	// Process "config.txt"
	myConfig = new SimConfiguration(sysConfigFile);
	// read basic data from config.txt by process 0 and distribute them to other processes
	// the distribution is carried out in the SimConfiguration.read_config_file() [???]
	// if(mpi_rank_g ==0 ){
	myConfig->read_config_file();
	// }
	
	// Handle trajectory storage and filename
	if(myConfig->write_trajectory()) {
		trajectoryFile = myConfig->get_trajectory_file();
		if (trajectoryFile == NULL){
			// ERRORMSG("no trajectory file name provided");
			LOG_ERR("No trajectory file name provided")
			return false;
		}
		// Open trajectory file
		ofpTrajectory = new ofstream (trajectoryFile, ios::out);
		if (!ofpTrajectory) {
			LOG_ERR("Can't open "<< trajectoryFile)
			return false;
		}
	}
	
	// Handle sysDataFile
	
	// sysDataFile contains parameters and topology data for building molecule & atom systems
	// these data will be read in and filled into system data structures
	sysDataFile = myConfig->get_sysData_file();
	myParams = new Parameters(sysDataFile);
	
	// Setup Ensemble
	myEnsemble = new Ensemble();
	myEnsemble->set_mol_system(myParams);
	myEnsemble->set_configuration(myConfig);
	// set up box, boundary, cell manager, cell list and pairlists.
	myEnsemble->compute_bound_box();
	myEnsemble->set_pairlist();
	
	myEnsemble->set_sysdata_file(of);
	myEnsemble->set_result_file(ofp);
	
	// Setup Integrator
	int integratorType = myConfig->get_integrator_type();
	switch (integratorType) {
		case GearIntegrator_TYPE:
			myIntegrator = new GearIntegrator(myEnsemble, myConfig);
			break;
		case VelocityIntegrator_TYPE:
			myIntegrator = new VelocityIntegrator(myEnsemble, myConfig);
			break;
		case LFIntegrator_TYPE:
			myIntegrator = new LFIntegrator(myEnsemble, myConfig);
			break;
		case NHIntegrator_TYPE:
			myIntegrator = new NHIntegrator(myEnsemble, myConfig);
			break;
		default:
			LOG_ERR("Unknown Integrator Type: "<< integratorType)
			return false;
	}
	
	// Last check in case Integrators construction via new failed (???)
	// it's too much?
	if (myIntegrator == NULL)
		ERRORMSG("Null integrator");
	
	myIntegrator->set_result_file(ofp);
	
	
	// Generate measuremental (??) objects
	Int nSamples = myConfig->get_n_steps() - myConfig->get_start_sampling();
	if (myConfig->get_sampling_ts_freq() != 0){
		nSamples /= myConfig->get_sampling_ts_freq(); // in the simConfiguration.cpp
	}
	if (myConfig->compute_rdf()){
		rdf = new RDF(myEnsemble, nSamples);
	}
	
	// Write configuration and ensemble info in sysData.out
	myConfig->write_config(*of);
	myEnsemble->write_ensemble_info(*of);
	
	return true;
}

//void Simulation::run(void) // tray transfer the rank

void Simulation::run(void){
	Int counter = 0;
	Int numSteps = myConfig->get_n_steps();
	Int printSteps =  myConfig->get_print_ts_freq();
	Int samplingSteps = myConfig->get_sampling_ts_freq();
	Int startSteps = 0;
	
	// Jc: MPI initialization
	/* WARNING */
	// int rank, size; // these definitions shadows the private fields and mpi_rank/mpi_size globals
	// MPI_Comm_size(MPI_COMM_WORLD, &size);
	// MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	// Jc:
	
	startSteps = myConfig->get_start_steps();
	
	// this check is unnecessary, ofp has been already opened.
	// What about all other of[something]?
	if (ofp == NULL)
		ERRORMSG ("null output stream");
	
	/* Write Headers Moved to setup_outputs()
	 if(mpi_rank_g == 0){
	 *ofp << "Energy" << endl;
	 *ofp << "Steps" << '\t' << " L-J" << '\t' << " Angle" << '\t' << " Bond" << '\t' << " UB-hydrogens "<< "Real" << '\t';
	 *ofp << "  Long" << '\t' << "  Correct" << '\t' << "  molCorrect" << '\t' << "  surfCorrect" << '\t';
	 *ofp << "  AtomTemp" << '\t' << "  MolTemp" << '\t' << "  Total" << endl;
	 
	 *ofpr << "## Step  " <<"Atom_Pdiag  " << "shear_stress  "<< "antisymm_SS  " <<"first_normal_Sh  " <<"second_normal_Sh  "<< "Mol_Pdiag  ";
	 *ofpr << "shear_stress  "<< "antisymm_SS  " <<"first_normal_Sh  " <<"second_normal_Sh  " << endl;
	 
	 *ofac << "## Step  " <<"Atom_Pdiag  " << "shear_stress  "<< "antisymm_SS  " <<"first_normal_Sh  " <<"second_normal_Sh  "<< "Mol_Pdiag  ";
	 *ofac << "shear_stress  "<< "antisymm_SS  " <<"first_normal_Sh  " <<"second_normal_Sh  " << endl;
	 
	 *ofvel << "# Step  " <<" velcomX  " << "  velcomY  "<< "  velcomZ  " <<"   Diff_vel  " <<"  Carbon_vel  "<< "   Density  ";
	 *ofvel <<  " C-O0_dist  " <<  "  stdC-O0   "<<  "   C-O1_dist  " << "   stdC-O1  " << "   Angle  " << "    anglSTDdev  " <<endl;
	 
	 *oflu << "# Step  " << "   lusStep    " << " Atomtemp  " << "  MolTemp  " << "  U_LJ    " << "  dUdV_LJ  " << "  d2UdV2_LJ   ";
	 *oflu << "      -bdUdV    " << "    (-bdUdV)2  "<< "      -bd2UdV2   " << "     UdUdVF      "  <<"  real2ndDev  " <<"  long2ndDev  "<< endl;
	 
	 *ofavl << "# Step  " << "   count   "<<"     bUconf/F  " << "    (bUconf/F)2  " << "    bavU/F  " << "     (bavU/F)2  " << "   Kinav  " << "  1/Kin   ";
	 *ofavl << "      Uen  " <<"      Uen2   " << "     Kin2  "<< "      Etot    " << "      Etot2   " << "      Kinatom   "<< "       1/Kinatom   " << "      Kinatom2  " <<endl;
	 
	 *oftp << "## Step     " << " U (kJ/mol) "  << " Cv (kJ/molK)  " << "   P (Mpa)" << endl;
	 
	 *ofind << "## Step     " << "   Real "  << "   Long  " << "   SelfCorrect  " << "   MolCorrect " << "  SurfacePol   " << "    TotInduction  "<< endl;
	 }
	 */
	
	if(startSteps > 0){
		myIntegrator->run(startSteps);
		counter += startSteps;
	}
	
	myEnsemble->zero_sum();
	
	if(mpi_rank_g == 0){
		//Jc: by Jianhui to output the system information to pTensor.dat
		myEnsemble->write_systemInfo(*ofPressureTensor);
	}
	
	while (counter < numSteps){
		for (Int i = 0; i < printSteps; i += samplingSteps){
			myIntegrator->run(samplingSteps);    // JC sampling steps too large in the system (???)
			myEnsemble->sampling();              // JC what is the function of sampling
			// RDF sampling
			if (rdf){
				rdf->sampling();                 // JC what is the function of sampling
			}
			if (mpi_rank_g == 0){
				myEnsemble->write_pressureTensor(*ofPressureTensor,counter + i); // to pTensor.out
			}
		}
		counter += printSteps;
		if(mpi_rank_g == 0){
			
			myEnsemble->write_result(*ofp, counter); // JC: pressure tensor to result.out
			
			myIntegrator->write_state();            // JC: data output is carried out in the Integrator print data every sampling time steps
			//myIntegrator->write_statexyz();
			
			myEnsemble->write_press_res(*ofpr, counter); // to pressureResult.out
			myEnsemble->write_pressac(*ofac, counter); // to pressureAcum.out
			if (myConfig->compute_induction()){
				myEnsemble->write_resultInduction(*ofind, counter); // to resultInduction.out
			}
			if (myConfig->compute_lustig()){
				myEnsemble->write_lustig(*oflu, counter); // to Lustig.out
			}
			if (myConfig->compute_lustig()){
				myEnsemble->write_TPbylustig(*oftp, counter); // to ThermoProp.out
			}
			if (myConfig->compute_lustig()){
				myEnsemble->write_averaglus(*ofavl, counter); // to LustigAverages.out
			}
			myEnsemble->write_velbehav(*ofvel, counter); // to velbehav.out
		}
	}
	
	
	if(mpi_rank_g == 0){
		myEnsemble->write_sampling_data();     // Jc: write Mol Trajectory and Velocity
	}
	
}


/* NEVER USED */
// void Simulation::finish(void) { ; }
