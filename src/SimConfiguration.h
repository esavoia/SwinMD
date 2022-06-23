/** SimConfiguration.h -- handle NEMD configuration input, maintain the global configuration 
 **    data structure for the simulation system.
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
// JC modified by Jianhui Li to cpmply with the new compiler
#ifndef SIMCONFIGURATION_H
#define SIMCONFIGURATION_H

#include <iostream>
#include <fstream>

// #include <stdio.h> //???
// #include <stdlib.h>
// #include "NEMD_defs.h"
#include "utils/Vector3.h"
// #include <mpi.h>
using namespace std;
/**
 ** Forward declarations
 **/
enum EnsembleStatus { NVT, NPT, NVE, NPE };
enum RunType {EQUIL, COMPR, SHEAR, ELONG};

class SimConfiguration {
	/**
	 ** Data Member -- most data are global, setup at the begining of the simulation
	 **      and should be not modified on the run-tume by other objects.
	 **/
private:
	// I/O file names
	const char*   sysConfigFile;          // name of system config data file
	char   sysDataFile[80];         // input file
	char   coordinateFile[80];      // input file
	char   velocityFile[80];        // input file, required if readVelocity flag is on
	char   resultFile[80];          // output file
	char   restartFile[80];         // used as backup file as well, optional
	char   trajectoryFile[80];      // atomic trajectory output, optional
	char   dumpFile[80];            // have to consider what data will be dumped into it
	char   molTrajectoryFile[80];   // molecular trajectory output for MSD analysis, optional
	char   molVelocityFile[80];     // molecular velocity output for VACF analysis, optional
	char   thermotype[80];          // Thermostat type name
	
	/**
	 ** timestep and controls
	 **/
	Double  timeStep;               // time length of one step, default 1 fs
	Int     nSteps,                 // total number of steps to loop, default 1
	startSteps;             // from which steps to start as "equlibration", default 0
	// TS for multi-time step integration
	Int     shortTsFreq,            // time steps btw two short force evaluations, default 1
	longTsFreq;             // time steps btw long force evaluations, default 1
	// TS for analysis and output control
	Int     averageTsFreq,          // TS Freq of average calculation, default nSteps
	printTsFreq,            // default nSteps
	backupTsFreq;           // default nSteps
	// TS control for sampling or measurement such as rdf or msd
	Int     startSampling,          // from which steps to start sampling, such as RDF, MSD...
	samplingTsFreq;         // Freq for sampling
	double  couplthermo;	    // copuling constant for NH or Berendsen thermostat
	double  couplBaros;		    // copuling constant for NH Barostat
	
	// TS control for TTCF
	Int     startTTCF,              // from which steps to start TTCF
	ttcfTsFreq;             // Freq to conduct TTCF and output data
	Int     nBins;                  // number of bins for computing RDF
	Int     integratorType;         // 0 - Gear; 1 - Velocity Verlet
	Double  elapsedTime;            // time since startSteps
	// this may be defined in other classes or method
	unsigned seed;                  // seed for random velocity, default 1
	int rank,size;                  // Jc: MPI features
	/**
	 ** System info
	 **/
	EnsembleStatus ensembleStatus;
	Double  temperature,
	pressure,
	volume,
	energy,
	density,
	cutoff,                 // cutoff for vdW forces
	cutoffEw,               // cutoff for Ewald real space if existing electric force
	cutoffInd,              // cutoff for Induction computation/minimization
	kCutoff,                // cutoff in k space for ewald sum
	alpha,                  // ewald convergent factor
	accuracy,               // accuracy for converging ewald sums
	cutBuff,                // buffer distance for building neighbour list
	limitRDF,               // limiting distance for RDF calculation
	nMie,                   // repulsive exponent of Mie potential
	mMie,                   // attractive exponent of Mie potential
	nBuff,                  // "n" exponent of Buff potential
	mBuff,                  // "m" exponent of Buff potential
	deltaBuff,              // "delta" parameter of Buff potential
	gammaBuff,              // "gamma" parameter of Buff potential
	Hfactor,                // reduction factor for Hydrogens in vdW AMOEBA
	klUB,                   // harmonic constant for UB potential for Hydrogens in AMOEBA
	l0UB,                   // Urey-Bradley ideal distance between Hydrogens in AMOEBA
	adamp,                  // controls the strength of damping
	switchDist;             // Distance at which shifting function is active
	Double  CFB, DFB, TFB;          // feed back constants if using Gaussain constraint methods
	Double  tolerance;              // tolerance for SHAKE constraint method
	Int     maxCount;               // max number of iterative loops for SHAKE
	Int     bendingType;            // Bending type: 1:Harmonic, 2:AMOEBA
	// Int     nOriginsMSD;            // number of origins for MSD statistic analysis
	// Int     nOriginsVACF;           // number of origins for VACF statistic analysis
	
	// box and cell parameters,  we may not need these as input
	Double  lxBox,                  // length of box in x-diraction
	lyBox,
	lzBox,
	cellSideLen;            // side length of a cell (using cubic cell)
	
	// run type - equilibrium, compress, shear or elongation
	RunType runType;
	Double  rate;                   // default 0 for equilibrium. It can be shear, compress
	bool    doEquilibrium,          // or elongation rate depending on testType.
	doCompression,
	doShear,
	doElongation;
	
	// flags
	bool    printSysConfig,         // print out system config data
	writeTrajectory,        // print out trajectory data
	readVelocity,           // read velocity otherwise random velocity
	isNewStart,             // read new sysConfig instead of reStart
	zeroAverage,            // set averages to zero
	useCellPairList,        // use cell method to build pair list
	// useGreenKubo,           // use Green-Kubo for transport properties, default false
	useTTCF,                // use TTCF for transport properties, default false
	useAtomThermo,          // use atom thermostat, default: true
	constraintOn,           // bond constraints, default: false.
	harmonicbend,           // Harmonic potential for both angle and bond, default: false.
	amoebabend,             // AMOEBA anharmonic potential for both angle and bond, default: false.
	binaryOutput,           // if writing output in binary formate, default: false
	switchOn,               // long-range L-J force correct with shifting fun,default: false
	lCorrectOn,             // long-range correction for L-J force, default: false
	computeCoulomb,         // use Coulomb term for electricstatic force
	computeEwald,           // use Ewald method for electricstatic force
	computeWolf,            // use Wolf term for electricstatic force
	computeMultiPol,        // use Multipolar electrostatic interctions
	computeInduction,       // compute Induction Energy for Multipole systems with Ewald method
	computeTorq,            // compute Torques and add/translate them into atomic forces
	dampReal,               // DEPRECATED!! apply damping function to the real-space part of Ewald
	dampInduction,          // apply damping function to Induction calculation
	dampAmoeba,             // Damping in the Induced part according to Tinker code
	doInductionMin,         // do Minimization of Induced Dipoles in the Induction calculation
	computeMie,             // use Mie potential for vdW force
	computeBuffvdW,         // use buffered potential for vdW force (AMOEBA)
	computeharmoUB,         // use a harmonic Urey-Bradley term for Hydrogens (AMOEBA)
	computeSurfCorrect,     // compute surface dipole correction for Ewald force
	hasExtForces,           // are there external forces, default: false
	hasFixedAtom,           // are there fixed atoms, mostly for confined sys, default false
	computeMSD,
	computeVACF,
	computeLustig,          // compute Lustig averages
	computeRDF;
	
	// flags for ensemble status - assume constant N particles
	bool    constantVolume,
	constantPressure,
	constantTemperature,
	constantEnergy;
	
	
	/**
	 ** constructor and destructor
	 **/
public:
	// SimConfiguration();
	SimConfiguration(const char configFile[]);
	~SimConfiguration() { ; }

private:
	void set_run_type(Int tType);
	void set_ensemble_status(Int e, Double v1, Double v2);
	/**
	 ** methods
	 **/
public:
	// initialise simulation parameters to default value, most of them will be
	// re-setup after reading imput files.
	void init();
	// read input files and set up simulation system
	void read_config_file();
	void read_config_file_old();
	// void write_config(FILE* fptr);
	void write_config(ofstream &ofp);
	void set_box(Double lx, Double ly, Double lz);
	void set_cutoff(Int nAtoms, Double volume);
	void set_use_cell_flag(bool useCellList);
	void set_cell_side(Double len);
	void set_density(Double d);
	void set_temperature(Double t);
	void set_density(Int nAtoms);
	void set_temperature(Int nAtoms, Double kinEnergy);
	void set_elapsedTime(Double currentSteps);
	
	
	// void set_freedom_degree();
	
	// inline Int get_freedom_degree() const { return freedomDegree; }
	Double get_timestep() const { return timeStep; }
	Int get_n_steps() const { return nSteps; }
	Int get_start_steps() const { return startSteps; }
	Int get_short_ts_freq() const { return shortTsFreq; }
	Int get_long_ts_freq() const { return longTsFreq; }
	Int get_avg_ts_freq() const { return averageTsFreq; }
	Int get_print_ts_freq() const { return printTsFreq; }
	Int get_backup_ts_freq() const { return backupTsFreq; }
	// inline Int get_dump_ts_freq() const { return dumpTsFreq; }
	Int get_start_sampling() const { return startSampling; }
	Int get_sampling_ts_freq() const { return samplingTsFreq; }
	Int get_start_TTCF() const { return startTTCF; }
	Int get_ttcf_ts_freq() const { return ttcfTsFreq; }
	Int get_integrator_type() { return integratorType; }
	int get_ensemble_status() {return ensembleStatus; }
	// Int get_msd_nOrigins() { return nOriginsMSD; }
	// Int get_vacf_nOrigins() { return nOriginsVACF; }
	Int get_rdf_nbins() { return nBins; }
	unsigned get_seed() { return seed; }
	
	// inline Double get_max_displacement() const { return maxDisplacement; }
	Double get_cutoff() const { return cutoff; }
	Double get_cutoffEw() const { return cutoffEw; }
	Double get_cutoffInd() const { return cutoffInd; }
	Double get_adamping() const { return adamp; }
	Double get_cutBuff() const { return cutBuff; }
	Double get_k_cutoff() const { return kCutoff; }
	Double get_alpha() const { return alpha; }
	Double get_accuracy() const { return accuracy; }
	Double get_limit_RDF() const { return limitRDF; }
	Double get_elapsed_time() const { return elapsedTime; }
	Double get_cell_side_len() const { return cellSideLen; }
	Double get_box_volume() const { return (lxBox*lyBox*lzBox);}
	Vector3 get_box_dimension() const { return Vector3(lxBox, lyBox, lzBox); }
	Double  get_density() const { return density; }
	Double get_density(Int nAtoms) const {
		double vol = lxBox*lyBox*lzBox;
		if (vol > 0.0)  return nAtoms/vol;
		return 0.0;
	}
	Double get_temperature() const { return temperature; }
	Double get_coupl_const() const { return couplthermo; }
	Double get_pressure() const { return pressure; }
	Double get_baros_const() const { return couplBaros; }
	Double get_CFB() const { return CFB; }
	Double get_DFB() const { return DFB; }
	Double get_TFB() const { return TFB; }
	Double get_tolerance() const { return tolerance; }
	Double get_max_count() const { return maxCount; }
	Double get_switch_distance() const { return switchDist; }
	Double get_kl_ub() const { return klUB; }
	Double get_l0_ub() const { return l0UB; }
	Double get_nMie_param() const { return nMie; }
	Double get_mMie_param() const { return mMie; }
	Double get_nBuff_param() const { return nBuff; }
	Double get_mBuff_param() const { return mBuff; }
	Double get_deltaBuff_param() const { return deltaBuff; }
	Double get_gammaBuff_param() const { return gammaBuff; }
	Double get_hfactor_reduction() const { return Hfactor; }
	
	
	char* get_coordinate_file() { return coordinateFile; }
	char* get_velocity_file() {return velocityFile; }
	char* get_sysData_file() { return sysDataFile; }
	char* get_result_file() { return resultFile; }
	char* get_restart_file() { return restartFile; }
	char* get_trajectory_file() { return trajectoryFile; }
	char* get_dump_file() { return dumpFile; }
	char* get_mol_traject_file() { return molTrajectoryFile; }
	char* get_mol_velocity_file() { return molVelocityFile;  }
	char* get_thermo_type() { return thermotype; }
	
	
	bool is_print_config_on() const { return printSysConfig; }
	bool is_read_velocity() const { return readVelocity; }
	bool is_new_start() const { return isNewStart; }
	bool is_zero_avg_on() { return zeroAverage; }
	bool is_constraint_on() { return constraintOn; }
	bool is_switch_on() { return switchOn; }
	bool write_trajectory() { return writeTrajectory; }
	bool use_long_correct() { return lCorrectOn; }
	bool use_cell_list() { return useCellPairList; }
	bool use_TTCF() { return useTTCF; }
	bool use_harmonic_potential() { return harmonicbend; }
	bool use_amoeba_bending() { return amoebabend; }
	// bool use_Green_Kubo() { return useGreenKubo; }
	bool use_bin_output() { return binaryOutput; }
	bool use_atom_thermo() { return useAtomThermo; }
	bool compute_coulomb() { return computeCoulomb; }
	bool compute_mie() { return computeMie; }
	bool compute_buffamoeba() { return computeBuffvdW; }
	bool compute_harmonic_ub() { return computeharmoUB; }
	bool compute_ewald_force() { return computeEwald; }
	bool compute_surf_correct() { return computeSurfCorrect; }
	bool compute_wolf() { return computeWolf; }
	bool compute_multipol() { return computeMultiPol; }
	bool compute_induction() { return computeInduction; }
	bool compute_torques() { return computeTorq; }
	bool use_damping_real() { return dampReal; }
	bool use_damping_induction() { return dampInduction; }
	bool do_induction_minimization() { return doInductionMin; }
	bool use_damping_indAmoeba() { return dampAmoeba; }
	// bool has_elec_forces() { return hasElecForces; }
	bool has_ext_forces() { return hasExtForces; }
	bool has_fixed_atom() { return hasFixedAtom; }
	bool compute_msd() { return computeMSD; }
	bool compute_rdf() { return computeRDF; }
	bool compute_lustig() { return computeLustig; }
	bool compute_vacf() { return computeVACF; }
};

#endif
