/** Integrator.cpp --
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
// JC modified by Jianhui LI to comply with the new compiler
#include "Integrator.h"
#include "Errors.h"

Integrator::Integrator(Ensemble* ensemble, SimConfiguration* simConfig)
{
    myEnsemble = ensemble;
    mySimConfig = simConfig;
    myAtoms = myEnsemble->atoms;
    numAtoms = myEnsemble->nAtoms;
    myMols = myEnsemble->molecules;
    numMols = myEnsemble->nMols;
    numBonds = myEnsemble->nBonds;
    molNDF = 3*numMols - 3;                                          // Degrees of Freedom to Compute T
    if (mySimConfig->is_constraint_on())
        atomNDF = 6*numMols - 3;   // Degrees of Freedom for CO2 ??
    else
        atomNDF = 3*numAtoms - 3;                                    // Degrees of Freedom to Compute T
    useAtomThermo = mySimConfig->use_atom_thermo();
    timeStep = mySimConfig->get_timestep();
    shortSteps = mySimConfig->get_short_ts_freq();
    longSteps = mySimConfig->get_long_ts_freq();
    numTimeSteps = mySimConfig->get_n_steps();
    dumpFile = mySimConfig->get_dump_file();
    drMax = mySimConfig->get_cutBuff()*mySimConfig->get_cutBuff()/4.0;
    temperature = mySimConfig->get_temperature();
    thermotype = mySimConfig->get_thermo_type();
	molTemp = temperature;
    atomTemp = temperature;
    couplconst = mySimConfig->get_coupl_const();
    TFB = mySimConfig->get_TFB();
    statEnsem = mySimConfig->get_ensemble_status();
    pressure  = mySimConfig->get_pressure() / 1.66053904042716;
    molPt     = pressure;
    molVirt   = - pressure;

    NpT       = false;
    if (statEnsem == 1){
        NpT        = true;
        BarosConst = mySimConfig->get_baros_const();
	DEBUGMSG(" NpT status ON in Integrator.cpp");

        if (strcmp(thermotype, "NoseHoover")){
	    ERRORMSG(" NpT Ensamble must use NH thermostat");
        }
    }
    else {
	DEBUGMSG("   NVT ensemble  ");
    }

    use_harmonic = mySimConfig->use_harmonic_potential();
    use_amoeba   = mySimConfig->use_amoeba_bending();

    atomKinEnergy = 0.0;
    molKinEnergy = 0.0;
    currentSteps = 0;
    nht = 0.0;
    fnh = 0.0;
    fnhpl = 0.0;
    fnhmi = 0.0;
    nnh   = 0.0;
    nnhpl = 0.0;
    nnhmi = 0.0;
    Mtot = 0.0;
    volume = myEnsemble->boxLx*myEnsemble->boxLy*myEnsemble->boxLz;
    Vpl   = volume;
    Vmin  = volume;

    angleForce = NULL;
    bondForce  = NULL;
    ubForce    = NULL;
//JC    dihedralForce = NULL;
//JC    improperForce = NULL;
    ljForce  = NULL;
    buffForce = NULL;
    mcyForce  = NULL;
    mieForce  = NULL;
    saapForce  = NULL;
    buckinghamForce  = NULL;
    bbvForce = NULL;
    ewaldForce = NULL;
    multiElec = NULL;
    inductionForce = NULL;
    wolfForce = NULL;
    // constraintForce = NULL;
    nccForce = NULL; //JC
    has_angle_force = false;
    has_bond_force = false;
    has_constr_force = false;
    use_ub_harmonic  = false;
    computeEwald = false;
    computeMultipole = false;
    computeInduction = false;
    computeWolf  = false;
    computeBuff  = false;
    // has_elec_force = false;
    // ofo  = new ofstream("state.out", ios::out); // JC
}

void Integrator::initialise()
{
    // initialise forces

    //nccForce = new NccForce(myEnsemble);

    if (mySimConfig->compute_buffamoeba())
    {
	computeBuff = true;
        buffForce = new BuffvdWForce(myEnsemble);
        forceGroup.push_back(buffForce);
    }
    else {
        ljForce = new LJForce(myEnsemble);
        forceGroup.push_back(ljForce);
    }
    //mcyForce = new MCYForce(myEnsemble);
    //forceGroup.push_back(mcyForce);
    //mieForce = new MieForce(myEnsemble);
    //forceGroup.push_back(mieForce);
    //saapForce = new SAAPForce(myEnsemble);
    //forceGroup.push_back(saapForce);
    //bbvForce = new BBVForce(myEnsemble);
    //forceGroup.push_back(bbvForce);
    //buckinghamForce = new BuckinghamForce(myEnsemble);
    //forceGroup.push_back(buckinghamForce);
     // forceGroup.push_back(nccForce);
    // if (mySimConfig->has_elec_forces())
    if (mySimConfig->compute_ewald_force())
    {
        computeEwald = true;
        ewaldForce = new EwaldForce(myEnsemble);
        // forceGroup.push_back(ewaldForce);
    }
    if (mySimConfig->compute_multipol())
    {
        computeMultipole = true;
        multiElec = new Multipole(myEnsemble);

        if (mySimConfig->compute_induction())
        {
            computeInduction = true;
            inductionForce = new Induction(myEnsemble);
	}
    }
    if (mySimConfig->compute_wolf())
    {
        computeWolf = true;
        wolfForce = new WolfForce(myEnsemble);
        // forceGroup.push_back(ewaldForce);
    }
    if (mySimConfig->is_constraint_on()&&(numBonds > 0))
    {
        has_constr_force = true;
        // constraintForce = new ConstraintForce(myEnsemble, mySimConfig->get_CFB(), mySimConfig->get_DFB());
    }
    if (!mySimConfig->is_constraint_on()&&(numBonds > 0))
    {
        has_bond_force = true;
        bondForce = new BondForce(myEnsemble);
        // forceGroup.push_back(bondForce);
    }
    if (!mySimConfig->is_constraint_on() && (myEnsemble->nAngles > 0))
    {
        has_angle_force = true;
        angleForce = new AngleForce(myEnsemble);
        // forceGroup.push_back(angleForce);
    }
    if (mySimConfig->compute_harmonic_ub())
    {
	if (mySimConfig->is_constraint_on())
    	{
	    DEBUGMSG(" WARNING!! ACHTUNG!! PELIGRO!!");
	    DEBUGMSG(" WARNING!! ACHTUNG!! PELIGRO!!");
	    DEBUGMSG(" Are you using a flexible Urey-Bradley term between Hydrogens but with constraints ??");
	    ERRORMSG(" Check your model and/or turn off the constraint");
	}

        use_ub_harmonic = true;
        ubForce = new UreyBradleyForce(myEnsemble);
        // forceGroup.push_back(ubForce);
    }
    if (myEnsemble->nDihedrals > 0)
    {
        // has_dihedral_force = true;
        // dihedralForce = new DihedralForce(myEnsemble);
    }
    if (myEnsemble->nImpropers > 0)
    {
        // has_improper_force = true;
        // improperForce = new ImproperForce(myEnsemble);
    }
}

void Integrator :: computePressure()
{
    Double molCorrectVirial[9];
    Molecule* mol;
    Atom* atom;

    if (statEnsem == 1){
        volume = myEnsemble->boxLx*myEnsemble->boxLy*myEnsemble->boxLz;
    }
    for (Int i = XX; i <= ZZ; i++)
        molCorrectVirial[i] = 0.0;
    for (Int i = 0; i < numMols; i++)
    {
        mol = &myMols[i];
        Vector3 Ri = mol->massCenter;
        for (Int j = 0; j < mol->numAtoms; j++)
        {
            Vector3 fj = mol->myAtoms[j]->force;
            Vector3 rj = mol->myAtoms[j]->position;
            molCorrectVirial[XX] += fj.x*(rj.x - Ri.x);
            molCorrectVirial[XY] += fj.y*(rj.x - Ri.x);
            molCorrectVirial[XZ] += fj.z*(rj.x - Ri.x);
            molCorrectVirial[YX] += fj.x*(rj.y - Ri.y);
            molCorrectVirial[YY] += fj.y*(rj.y - Ri.y);
            molCorrectVirial[YZ] += fj.z*(rj.y - Ri.y);
            molCorrectVirial[ZX] += fj.x*(rj.z - Ri.z);
            molCorrectVirial[ZY] += fj.y*(rj.z - Ri.z);
            molCorrectVirial[ZZ] += fj.z*(rj.z - Ri.z);
        }
    }
    for (Int i = XX; i <= ZZ; i++)
        myEnsemble->molVirial[i] = myEnsemble->virial[i] - molCorrectVirial[i];

    for (Int i = 0; i < numAtoms; i++)
    {
        Atom* atom = &myAtoms[i];
        myEnsemble->atomPressure[XX] += atom->momentum.x*atom->momentum.x/atom->mass;
        myEnsemble->atomPressure[XY] += atom->momentum.x*atom->momentum.y/atom->mass;
        myEnsemble->atomPressure[XZ] += atom->momentum.x*atom->momentum.z/atom->mass;
        myEnsemble->atomPressure[YY] += atom->momentum.y*atom->momentum.y/atom->mass;
        myEnsemble->atomPressure[YZ] += atom->momentum.y*atom->momentum.z/atom->mass;
        myEnsemble->atomPressure[ZZ] += atom->momentum.z*atom->momentum.z/atom->mass;
    }
    myEnsemble->atomPressure[YX] = myEnsemble->atomPressure[XY];
    myEnsemble->atomPressure[ZX] = myEnsemble->atomPressure[XZ];
    myEnsemble->atomPressure[ZY] = myEnsemble->atomPressure[YZ];
    for (Int i = 0; i < numMols; i++)
    {
        mol = &myMols[i];
        myEnsemble->molPressure[XX] += mol->momenta.x*mol->momenta.x/mol->mass;
        myEnsemble->molPressure[XY] += mol->momenta.x*mol->momenta.y/mol->mass;
        myEnsemble->molPressure[XZ] += mol->momenta.x*mol->momenta.z/mol->mass;
        myEnsemble->molPressure[YY] += mol->momenta.y*mol->momenta.y/mol->mass;
        myEnsemble->molPressure[YZ] += mol->momenta.y*mol->momenta.z/mol->mass;
        myEnsemble->molPressure[ZZ] += mol->momenta.z*mol->momenta.z/mol->mass;
    }
    myEnsemble->molPressure[YX] = myEnsemble->molPressure[XY];
    myEnsemble->molPressure[ZX] = myEnsemble->molPressure[XZ];
    myEnsemble->molPressure[ZY] = myEnsemble->molPressure[YZ];

    for (Int i = XX; i <= ZZ; i++)
    {
        myEnsemble->atomPressure[i] = (myEnsemble->atomPressure[i] + myEnsemble->virial[i])/volume;
        myEnsemble->molPressure[i] = (myEnsemble->molPressure[i] + myEnsemble->molVirial[i])/volume;
    }
}

//Jc: both the position and real positions are output every sampling time steps
//Jc: real position are used in the next simulation for calcualtion of the diffusion
//Jc:
void Integrator :: write_state()
{
    // write configurations to default file
//    ofstream ofo = ofstream("state.out", ios::out);
    ofo = new ofstream("state.out", ios::out);
    *ofo << "Coordinate at step: " << currentSteps << endl;
    *ofo << numAtoms << endl;
    for(int j = 0; j < numAtoms; j++)
        *ofo<<myAtoms[j].position.x<<' '<<myAtoms[j].position.y<<' '<<myAtoms[j].position.z<<endl;

    *ofo << "Real positions at step: " << currentSteps << endl; // Jc: modified by Jianhui Li
    *ofo << numAtoms << endl;
    for(int j = 0; j < numAtoms; j++)
        *ofo<<myAtoms[j].realPos.x<<' '<<myAtoms[j].realPos.y<<' '<<myAtoms[j].realPos.z<<endl;

    *ofo << "velocity at steps: " << currentSteps << endl;
    *ofo << numAtoms << endl;
    for(int j = 0; j < numAtoms; j++)
    {
        *ofo << myAtoms[j].velocity.x;
        *ofo << ' ' << myAtoms[j].velocity.y;
        *ofo << ' ' << myAtoms[j].velocity.z;
        *ofo << endl;
    }
    *ofo << endl;
    ofo->close();
    // write_bond(ofo); //to comply with the new compiler
}

void Integrator :: write_statexyz()
{
    string base(".xyz");
    string namef("stateout");
    string posit(9,'0');
    string stepxyz = to_string(currentSteps);
    posit.replace(posit.end()-stepxyz.size(),posit.end(),stepxyz);
    ofo = new ofstream(posit+namef+base, ios::out);
    *ofo << numAtoms << endl;
    *ofo << "output xyz by SwinMD" << endl;
    for(int j = 0; j < numAtoms; j++)
        *ofo<<myAtoms[j].typeName<<' '<<myAtoms[j].position.x*10<<' '<<myAtoms[j].position.y*10<<' '<<myAtoms[j].position.z*10<<endl;
    *ofo << endl;
    ofo->close();
}


void Integrator :: write_mol_position(Molecule *mol, ofstream &of)
{
     of<<"mol["<< mol->molID <<"] ";
     for (Int j = 0; j < mol->numAtoms; j++)
     {
         of <<"  atom[" << mol->myAtoms[j]->atomID << "] " << mol->myAtoms[j]->position.x;
         of <<' '<<mol->myAtoms[j]->position.y<<' '<<mol->myAtoms[j]->position.z<<"   ";
     }
     of << endl;
}


void Integrator :: write_bond(ofstream &of)
{
    double l1, l2, l3, angle;
    Vector3 r12, r32, r13;
    Molecule* mol;

    for (Int i = 0; i < numMols; i++)
    {
        mol = &myMols[i];
        int ns = mol->numAtoms;
        if (ns <= 1)
            return;
        of << "mol[" << mol->molID << "]:  ";
        if (ns == 2)
        {
            r12 = mol->myAtoms[1]->position - mol->myAtoms[0]->position;
            of << r12.length2() ;
        }
        else
        for(int site = 0; site < (ns-2); site++)
        {
            r12 = mol->myAtoms[site + 1]->position - mol->myAtoms[site]->position;
            l1 = r12.length2();
            r32 = mol->myAtoms[site + 2]->position - mol->myAtoms[site + 1]->position;
            l2 = r32.length2();
            r13 = mol->myAtoms[site + 2]->position - mol->myAtoms[site]->position;
            l3 = r13.length2();
            of << l1 << ' ' << l2 << ' ' << l3;
        }
        of << endl;
    }
}

void Integrator::set_dump_file(char* fname)
{
    dumpFile = fname;
}

void Integrator::set_result_file(ofstream* ofp)
{
    resultFile = ofp;
}
