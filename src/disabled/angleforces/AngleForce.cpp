/** AngleForce.cpp -- 
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
//JC modified by Jianhui Li to comply with the linux compiler
#include "AngleForce.h" 
#include "../utils/Vector3.h"


AngleForce::AngleForce(Ensemble* ensemble) : Force(ensemble)
{
    numAngles = myEnsemble->nAngles;
    angles = myEnsemble->angles;
    params = myEnsemble->myParams->angleParams;
//JC    forceIdentifier= "AngleForce"; // in order to comply the cluster compiler
    write_force_info(*sysdataFile);
    #ifdef DEBUG
        DEBUGMSG("create angle force ");
    #endif
    use_harmonic = myEnsemble->myConfig->use_harmonic_potential();
    use_amoeba   = myEnsemble->myConfig->use_amoeba_bending();

    if (use_amoeba) { DEBUGMSG(" AngleForce: Use AMOEBA anhamornic potential.");}
        else { DEBUGMSG(" AngleForce: harmonic potential"); }
}

void AngleForce::compute()
{
    Int atom1, atom2, atom3;
    Int angleType;
    Double theta, sinTheta, cosTheta, diff;
    Double theta0, forceConst;
    double Uangle, v, invSin;
    double diff2, diff3, diff4, anharm, deg;
    double c1, c2, c3 ,c4;
    Vector3 r12, r32;        // vector from atom i to atomj
    Double d12r, d32r;         // distance between atom i and atom j
    Vector3 force1, force2, force3;
    if (use_amoeba){
	c1 = 0.014;
	c2 = 5.6E-5;
	c3 = 7.0E-7;
	c4 = 2.2E-8;
    } else {
	c1 = 0.0;
	c2 = 0.0;
	c3 = 0.0;
	c4 = 0.0;
    }
    

    energy = 0.0;
    for (Int i = XX; i <= ZZ; i++)
        virial[i] = 0.0;

    for (Int i = 0; i < numAngles; i++)
    {
        atom1 = angles[i].atom1;
        atom2 = angles[i].atom2;
        atom3 = angles[i].atom3;
        angleType = angles[i].angleType;

        theta0 = params[angleType].theta0;
        forceConst = params[angleType].k;
        
        r12 = atoms[atom2].position - atoms[atom1].position; 
        myEnsemble->apply_pbc(r12);
        d12r = 1.0/r12.length();                // distance inverse between atom1 & atom2 
        r32 = atoms[atom2].position - atoms[atom3].position;
        myEnsemble->apply_pbc(r32);
        d32r = 1.0/r32.length();

        cosTheta = (r12*r32)*d12r*d32r;
	cosTheta = MAX(cosTheta, -0.9999999);
//cout << cosTheta << endl;
        sinTheta = sqrt(1.0 - cosTheta*cosTheta);
	invSin   = 1.0/sinTheta;
        theta = acos(cosTheta); 
        diff = theta - theta0;

        // 2k(theta-thea0)/sin(theta)
        Vector3 r12_d12 = r12*d12r;
        Vector3 r32_d32 = r32*d32r;
	if (use_harmonic) {
            Uangle = 0.5*forceConst*diff*diff;
            v      = forceConst*diff;
	} else if (use_amoeba) {
	    deg   = 180/PI;
	    diff *= deg;
	    diff2 = diff*diff;
	    diff3 = diff*diff2;
	    diff4 = diff2*diff2;

	    anharm = 1.0 - c1*diff + c2*diff2 - c3*diff3 + c4*diff4;
	    Uangle = forceConst*diff2*anharm/(deg*deg);
	    v      = 2.0*forceConst*diff*anharm/deg + forceConst*diff2*(-c1 + c2*2.0*diff - c3*3.0*diff2 + c4*4.0*diff3)/deg;
	}        

        energy += Uangle;        
        force1 = v*d12r*invSin*(r12_d12*cosTheta - r32_d32);
        force3 = v*d32r*invSin*(r32_d32*cosTheta - r12_d12);
        force2 = -(force1 + force3);
        virial[XX] += force1.x * r12.x + force3.x * r32.x;
        virial[XY] += force1.x * r12.y + force3.x * r32.y;
        virial[XZ] += force1.x * r12.z + force3.x * r32.z;
        virial[YX] += force1.y * r12.x + force3.y * r32.x;
        virial[YY] += force1.y * r12.y + force3.y * r32.y;
        virial[YZ] += force1.y * r12.z + force3.y * r32.z;
        virial[ZX] += force1.z * r12.x + force3.z * r32.x;
        virial[ZY] += force1.z * r12.y + force3.z * r32.y;
        virial[ZZ] += force1.z * r12.z + force3.z * r32.z;

        // add to the total force;
        atoms[atom1].force += force1;
        atoms[atom2].force += force2;
        atoms[atom3].force += force3;       
    }

    myEnsemble->myEnergy.angleEnergy = energy;
    for (Int k = XX; k <= ZZ; k++)
        myEnsemble->virial[k] -= virial[k];
}

void AngleForce::write_force_info(ofstream& of)
{
    of << "Angle Force" << endl;
    of << "number of angles: " << numAngles << endl;
    of << endl;
}

void AngleForce::write_energy(ofstream& of)
{
    of << "Angle Energy: " << energy << endl;
}

//JC string AngleForce::get_force_id() // in order to comply the cluster compiler
//Jc{
//JC  return AngleForce::forceIdentifier;
//Jc}
