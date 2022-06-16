/** UreyBradleyForce.cpp -- 
 **
 ** Copyright (C) 2o21
 ** Centre for Molecular Simulation (CMS)
 ** School of Information Technology
 ** Swinburne University of Technology
 ** PO Box 21 Hawthorn, Vic 3122, Australia
 **
 ** Author: Elton Oyarzua
 ** Email: eoyarzua@swin.edu.au
 **/

#include "UreyBradleyForce.h" 

UreyBradleyForce::UreyBradleyForce(Ensemble* ensemble) : Force(ensemble)
{
    numAtoms = myEnsemble->nAtoms;
    numMols = myEnsemble->nMols;
    myMols = myEnsemble->molecules;
    volume = myEnsemble->volume;
    params = myEnsemble->myParams->bondParams;
    kl = myEnsemble->myConfig->get_kl_ub();
    l0 = myEnsemble->myConfig->get_l0_ub();

//JC    forceIdentifier = "UreyBradleyForce"; // in order to comply the cluster compiler
    write_force_info(*sysdataFile);
    #ifdef DEBUG
        DEBUGMSG("Creating Urey-Bradley Force ");
    #endif
    use_harmonic = myEnsemble->myConfig->use_harmonic_potential();
    use_amoeba   = myEnsemble->myConfig->use_amoeba_bending();
}

UreyBradleyForce::~UreyBradleyForce()
{
    ;
}

void UreyBradleyForce::compute()
{
    Molecule *mol;
    int id1, id2, bondType;
    Vector3 r12, force;
    double r, diff, diff2;
    double val, Ubond, dubond, d2ubond;    

    energy = 0.0;
    dubond = 0.0;
    d2ubond = 0.0;

    for (Int i = XX; i <= ZZ; i++)
        virial[i] = 0.0;

    for (int nm = 0; nm < numMols; nm++)
    {
        int ns = myMols[nm].numAtoms;
        mol = &myMols[nm];

        Atom *atomH1 = mol->myAtoms[0];
        Atom *atomO  = mol->myAtoms[1];
        Atom *atomH2 = mol->myAtoms[2];

	id1 = atomH1->atomID;
	id2 = atomH2->atomID;

        r12 = atomH1->position - atomH2->position;
        myEnsemble->apply_pbc(r12);
        r = r12.length();
        diff = r - l0;

        Ubond = kl*diff*diff;
        val = 2.0*kl*diff;

        force    = -val*r12/r;
        energy  += Ubond;
	dubond  += val*r;
        d2ubond += -2.0*r*val + r*r*kl;

        virial[XX] += -force.x * r12.x;
        virial[XY] += -force.x * r12.y;
        virial[XZ] += -force.x * r12.z;
        virial[YX] += -force.y * r12.x;
        virial[YY] += -force.y * r12.y;
        virial[YZ] += -force.y * r12.z;
        virial[ZX] += -force.z * r12.x;
        virial[ZY] += -force.z * r12.y;
        virial[ZZ] += -force.z * r12.z;
        
        // add to total force
        atoms[id1].force += force;
        atoms[id2].force -= force;
    }
    myEnsemble->myEnergy.ubEnergy   = energy;
    myEnsemble->myLustig.ubDeriv    = dubond;
    myEnsemble->myLustig.ubSecDeriv = d2ubond;

    for (Int k = XX; k <= ZZ; k++)
        myEnsemble->virial[k] -= virial[k];
}

void UreyBradleyForce::write_force_info(ofstream& of)
{
    of << "Urey-Bradley force" << endl;
    of << "force constant: " << kl << endl;
    of << "ideal distance: " << l0 << endl;
    of << endl;
}

void UreyBradleyForce::write_energy(ofstream& of)
{
    of << " UB Energy: " << energy << endl;
}

// JC string UreyBradleyForce::get_force_id()  // in order to comply with the cluster compiler
// JC {
// JC  return UreyBradleyForce::forceIdentifier;
// JC }
