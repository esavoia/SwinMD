/** Integrator.h -- 
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

#ifndef Integrator_H
#define Integrator_H

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include "NEMD_defs.h"
#include "utils/Vector3.h"
#include "Atom.h"
#include "Molecule.h"
#include "Ensemble.h"
#include "SimConfiguration.h"
#include "forces/AngleForce.h"
#include "forces/BondForce.h"
#include "forces/UreyBradleyForce.h"
#include "forces/ConstraintForce.h"
// JC #include "DihedralForce.h" // in order to comply with the cluster compiler
// JC #include "ImproperForce.h" // in order to comply with the cluster compiler
#include "forces/LJForce.h"
#include "forces/BuffvdWForce.h"
#include "forces/MCYForce.h"
#include "forces/MieForce.h"
#include "forces/SAAPForce.h"
#include "forces/BuckinghamForce.h"
#include "forces/BBVForce.h"
#include "forces/EwaldForce.h"
#include "forces/Multipole.h"
#include "forces/Induction.h"
#include "forces/WolfForce.h"
#include "forces/NccForce.h"

#include <vector>
using namespace std;
 
/**
 ** Forward declarations
 **/


class Integrator{
    /**
     ** Data Member - 
     **/
    protected:
    Double timeStep;
    Int shortSteps, longSteps;              // steps btw two nonbonded evaluations / long electric evaluations
    Int numTimeSteps;
    Int numAtoms;
    Int numMols;
    Int numBonds;
    Int atomNDF;
    Int molNDF;
    Int currentSteps;
    Int nht;
    int statEnsem;
    Double fnh, nnh;
    Double fnhpl, nnhpl, Vpl;
    Double fnhmi, nnhmi, Vmin;
    Double Mtot;
    Vector3 Rm;

    Ensemble* myEnsemble;
    SimConfiguration* mySimConfig;
    Atom*     myAtoms;
    Molecule* myMols;
    vector<Force*> forceGroup;               // group of short and bonded forces
    AngleForce* angleForce;
    BondForce*  bondForce;
    UreyBradleyForce*  ubForce;
    // ConstraintForce* constraintForce;
//JC    DihedralForce* dihedralForce;
//JC    ImproperForce* improperForce;
    LJForce*    ljForce;
    BuffvdWForce*    buffForce;
    MCYForce*   mcyForce;
    MieForce*   mieForce;
    SAAPForce*    saapForce;
    BuckinghamForce*    buckinghamForce;
    BBVForce*   bbvForce;
    EwaldForce* ewaldForce;
    Multipole*  multiElec;
    Induction*  inductionForce;
    WolfForce*  wolfForce;
    NccForce* nccForce; // JC created by Jianhui Li, pointer to NCCForce
    
    bool has_angle_force;
    bool has_bond_force;
    bool has_constr_force;
    bool computeEwald;
    bool computeMultipole;
    bool computeInduction;
    bool computeWolf;
    bool computeBuff;
    bool has_elec_force;
    bool useAtomThermo;
    bool NpT;
    bool use_harmonic;
    bool use_amoeba;
    bool use_ub_harmonic;

    Double temperature, TFB, pressure;
    Double couplconst, BarosConst;
    char* thermotype;
    Double atomKinEnergy, molKinEnergy;
    Double volume;
    Double molTemp, atomTemp, molPt, molVirt;
    Double drMax;       // maximum displacement for updating pair list
    Int numUpdates;     // accumulated number of pair list updates

    char* dumpFile;               // for out put system data
    ofstream*   resultFile;
    ofstream*   ofo;   // JC
    /**
     ** constructor and destructor
     **/
    public:
    Integrator(Ensemble* ensemble, SimConfiguration* simConfig);
    virtual ~Integrator() { ; }

    /**
     ** methods
     **/
    public:
    virtual void initialise();
    virtual void run(Int numTimeSteps) = 0;
    virtual void write_state();
    virtual void write_statexyz();
    void computePressure();
    void write_mol_position(Molecule *mol, ofstream &of);
    void write_bond(ofstream &of);
    void set_dump_file(char* fname);
    void set_result_file(ofstream *ofp);
    // void set_force_group();
};

#endif
