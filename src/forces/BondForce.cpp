/** BondForce.cpp -- 
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

#include "BondForce.h" 

BondForce::BondForce(Ensemble* ensemble) : Force(ensemble)
{
    numBonds = myEnsemble->nBonds;
    bonds = myEnsemble->bonds;
    params = myEnsemble->myParams->bondParams;
//JC    forceIdentifier = "BondForce"; // in order to comply the cluster compiler
    write_force_info(*sysdataFile);
    #ifdef DEBUG
        DEBUGMSG("Create bond force ");
    #endif
    use_harmonic = myEnsemble->myConfig->use_harmonic_potential();
    use_amoeba   = myEnsemble->myConfig->use_amoeba_bending();

    if (use_amoeba) { DEBUGMSG(" BondForce: Use AMOEBA anhamornic potential.");}
        else { DEBUGMSG(" BondForce: harmonic potential"); }
}

BondForce::~BondForce()
{
    ;
}

void BondForce::compute()
{
    Int atom1, atom2, bondType;
    Vector3 r12;
    Double forceConst, r0, r, diff, diff2;
    Vector3 force;
    double val, Ubond, dubond, d2ubond;    
    double anharmo, c1, c2;
    if (use_amoeba) {
        c1 = 25.5;
	c2 = 379.3125; 
    } else {
	c1 = 0.0;
	c2 = 0.0;
    }


    energy = 0.0;
    dubond = 0.0;
    d2ubond = 0.0;

    for (Int i = XX; i <= ZZ; i++)
        virial[i] = 0.0;

    for(Int i = 0; i < numBonds; i++)
    {
        atom1 = bonds[i].atom1;
        atom2 = bonds[i].atom2;
        bondType = bonds[i].bondType;
        forceConst = params[bondType].k;
        r0 = params[bondType].r0;

        r12 = atoms[atom1].position - atoms[atom2].position;
        myEnsemble->apply_pbc(r12);
        r = r12.length();
        diff = r - r0;

	if (use_harmonic) {
            Ubond = 0.5*forceConst*diff*diff;
            val = forceConst*diff;
	} else if (use_amoeba) {
	    diff2   = diff*diff;
	    anharmo = 1.0 - c1*diff + c2*diff2;
	    Ubond   = forceConst*diff2*anharmo;
	    val     = 2.0*forceConst*diff*anharmo + forceConst*diff2*(-c1 + c2*2.0*diff);
	}

        force    = -val*r12/r;
        energy  += Ubond;
	dubond  += val*r;
        d2ubond += -2.0*r*val + r*r*forceConst;

        virial[XX] += force.x * r12.x;
        virial[XY] += force.x * r12.y;
        virial[XZ] += force.x * r12.z;
        virial[YX] += force.y * r12.x;
        virial[YY] += force.y * r12.y;
        virial[YZ] += force.y * r12.z;
        virial[ZX] += force.z * r12.x;
        virial[ZY] += force.z * r12.y;
        virial[ZZ] += force.z * r12.z;
        
        // add to total force
        atoms[atom1].force += force;
        atoms[atom2].force -= force;
    }
    myEnsemble->myEnergy.bondEnergy = energy;
    myEnsemble->myLustig.BondDeriv    = dubond;
    myEnsemble->myLustig.BondSecDeriv = d2ubond;

    for (Int k = XX; k <= ZZ; k++)
        myEnsemble->virial[k] -= virial[k];
}

void BondForce::write_force_info(ofstream& of)
{
    of << "Bond force" << endl;
    of << "Bond number: " << numBonds << endl;
    of << endl;
}

void BondForce::write_energy(ofstream& of)
{
    of << "Bond Energy: " << energy << endl;
}

// JC string BondForce::get_force_id()  // in order to comply with the cluster compiler
// JC {
// JC  return BondForce::forceIdentifier;
// JC }
