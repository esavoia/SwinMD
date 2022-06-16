/** MCYForce.cpp -- 
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

#include "MCYForce.h" 


MCYForce::MCYForce(Ensemble* ensemble) : Force(ensemble)
{
    #ifdef DEBUG
        DEBUGMSG("Creating MCYForce ");
    #endif
    
    mcyParams = NULL;
//JC    forceIdentifier = "LJForce";// in order to comply with the cluster compiler
    numAtoms = myEnsemble->nAtoms;
    numAtomTypes = myEnsemble->nAtomTypes;
    numMols = myEnsemble->nMols;
    params = myEnsemble->myParams;
    myMols = myEnsemble->molecules;
    cutOff = myEnsemble->myConfig->get_cutoff();
    cutOff2 = cutOff*cutOff;
    switchDist = cutOff2;           
    switchOn = myEnsemble->myConfig->is_switch_on();
    
    if (switchOn)
        switchDist = myEnsemble->myConfig->get_switch_distance()*myEnsemble->myConfig->get_switch_distance();
    if (switchOn && (cutOff2 > switchDist))
    {
        c1 = 1.0/(cutOff2 - switchDist);
        c1 = c1*c1*c1;
        c2 = cutOff2 - 3.0*switchDist;
        c3 = 4*c1;
    }
    else
    {
        c1 = 0.0;
        c2 = 0.0;
        c3 = 0.0;
    }

    computeCoulomb = myEnsemble->myConfig->compute_coulomb();

    mcyParams = new MCYParam[numAtomTypes*numAtomTypes];
          if (mcyParams  == NULL)
              ERRORMSG( "memory allocation of MCY Parameters error");

    set_mcy_param();


    if (myEnsemble->myConfig->use_long_correct())
        long_range_correct();           // calculate once at the beginning
    write_force_info(*sysdataFile);
}

// this implementation of compute() might be not efficient in parallelisation
// in case of OpenMP, data dependency may occur in the summation of force
// in case of MPI, many communication may be involved 
void MCYForce::compute()
{

    Atom *atomi, *atomj;
    Vector3 rij, fij;
    Double r, r_1, rij2, rij_2, sr2, sr6, sr12, tmpE, tmpf;
    double dulj, d2ulj, voll, vol2;
    Double qi, qij;
    Double rcut_2 = 1.0/cutOff2;
    Double elecE = 0.0;

    Vector3 rim, rjm, rijm;
    Double  molVirial[9];
    
    energy = 0.0;
    ulus   = 0.0;
    dudv   = 0.0;
    d2udv  = 0.0;

    voll = myEnsemble->boxLx*myEnsemble->boxLy*myEnsemble->boxLz;
    vol2 = voll*voll;
//cout << "volume lustig = " << voll << "  vol2  = " << vol2 << endl;
    
    for (Int i = XX; i <= ZZ; i++)
    {
        virial[i] = 0.0;
        molVirial[i] = 0.0;
    }

    // go through each atom
    for (Int i = 0; i < numAtoms; i++)
    {
        atomi = &atoms[i]; 
        Int id = atomi->molID;
        rim = myEnsemble->molecules[id].massCenter; //
        Int size = atomi->get_list_size();      // size of pairlist
        qi = atomi->scaledCharge;

// For CO2 the COM has to be the same as the Carbon position:
//if (atomi->atomType == 1 ){
//cout << atomi->position << " <- atom      Rij ->   " << rim  << endl;}

        // go through the pairlist of atomi
        for (Int j = 0; j < size; j++)   /// Int j = i+1; j<numAtoms; j++
        {
            atomj = atomi->myPairList[j];
//	if (atoms[i]->molID != atoms[j]->molID ){
	if (atomi->molID != atomj->molID){
            id = atomj->molID;    //
            rij = atomj->position - atomi->position;
            rijm = myEnsemble->molecules[id].massCenter - rim;  //


            myEnsemble->apply_pbc(rij);
            myEnsemble->apply_pbc(rijm);
            r = rij.length();
            rij2 = r*r;

            if (rij2 <= cutOff2)
            {
                int mcyid =  atomi->atomType * numAtomTypes + atomj->atomType;    // probably in the future we are goin to need numAtomtypes per molecule, as a vector.
                double aa1 = mcyParams[mcyid].AA1;
                double bb1 = mcyParams[mcyid].BB1;
	        double aa2 = mcyParams[mcyid].AA2;
	        double bb2 = mcyParams[mcyid].BB2;
                

//if (atomi->atomType == 1 && atomj->atomType == 1){
//cout << rij << " <- atom      Rij ->   " << rijm << endl;}

                r_1 = 1.0/r;
                rij_2 = r_1*r_1;
                tmpE = aa1*exp(-bb1*r) - aa2*exp(-bb2*r);
                tmpf = aa1*bb1*r_1*exp(-bb1*r) - aa2*bb2*r_1*exp(-bb2*r);
		double rijRij = rij*rijm;

		if (myEnsemble->myConfig->compute_lustig())
		{
		dulj  = tmpf*rijRij;                  // Lustig 2011 Eq. 15b
		d2ulj = ( aa1*bb1*bb1*exp(-bb1*r) - aa2*bb2*bb2*exp(-bb2*r) )*rijRij*rijRij*rij_2 + 2.0*dulj;	        // Lustig 2011 Eq. 15c
		}

                if(switchOn /*&& (rij2 > switchDist)*/)
                {/*
                    Double c = cutOff2 - rij2;
                    Double c4 = c*(2*rij2 + c2);
                                        
                    Double switchVal = c*c4*c1;
                    Double dSwitchVal = c3*(c*c-c4);
                    tmpf = tmpf*switchVal - tmpE*dSwitchVal;
                    tmpE *= switchVal;*/
//cout << "cut2 = " << rcut_2 << endl;
		    
		    double Eljcut   = aa1*exp(-bb1*cutOff) - aa2*exp(-bb2*cutOff);
		    tmpE = tmpE - Eljcut;
                }

                energy += tmpE;   // intra + inter
		if (myEnsemble->myConfig->compute_lustig())
		{
		ulus   += tmpE; 
    		dudv   += dulj;
    		d2udv  += d2ulj;
		}
                
                // computing coulomb term with NAMD shift function 
                /* if (computeCoulomb)
                {
                    qij = qi*atomj->scaledCharge;
                    Double efac = 1 - rij2*rcut_2;
                    Double prefac = qij*r_1*efac;
                    elecE += prefac*efac;
                    tmpf += prefac*rij_2*(1 + 3.0*rij2*rcut_2);
                } */

                if (computeCoulomb)
                {
                    qij = qi*atomj->scaledCharge;
                    Double e = qij*r_1;
                    Double f = e*rij_2;
                    // use GROMACS shift function
                    Double r_rcut = r/cutOff;
                    Double r4_rcut4 = rij2*rij2*rcut_2*rcut_2;
                    elecE += e*(1 - 1.666667*r_rcut + 1.666667*r4_rcut4 - r_rcut*r4_rcut4);
                    tmpf += f*(1 - 5*r4_rcut4 + 4*r_rcut*r4_rcut4);
                    // CHARMM's shift function
                    // elecE += e*(1 - rij2*rcut_2)*(1-rij2*rcut_2);
                    // tmpf += f*(1 - r/cutOff)*(1 - r/cutOff);
                } 

                fij = -tmpf*rij;   
                atomi->force += fij;
                atomj->force -= fij;
                virial[XX] += fij.x * rij.x;
                virial[XY] += fij.y * rij.x;
                virial[XZ] += fij.z * rij.x;
                virial[YX] += fij.x * rij.y;
                virial[YY] += fij.y * rij.y;
                virial[YZ] += fij.z * rij.y;
                virial[ZX] += fij.x * rij.z;
                virial[ZY] += fij.y * rij.z;
                virial[ZZ] += fij.z * rij.z;

                /* molVirial[XX] += fij.x * rijm.x;   // 
                molVirial[XY] += fij.y * rijm.x;
                molVirial[XZ] += fij.z * rijm.x;
                molVirial[YX] += fij.x * rijm.y;
                molVirial[YY] += fij.y * rijm.y;
                molVirial[YZ] += fij.z * rijm.y;
                molVirial[ZX] += fij.x * rijm.z;
                molVirial[ZY] += fij.y * rijm.z;
                molVirial[ZZ] += fij.z * rijm.z;   //  */

            }
          }
	}        
    }

    //   long range correction  
    energy += eLrc;
    virial[XX] += vLrc;
    virial[YY] += vLrc;
    virial[ZZ] += vLrc;
    myEnsemble->myEnergy.ljEnergy = energy + elecE;

    if (myEnsemble->myConfig->compute_lustig())
    {
    	dudv   = -dudv/(3.0*voll);
    	d2udv  = d2udv/(9.0*vol2);

//cout << "  dudv  =  " << dudv << "  d2udv =  " <<d2udv << endl;

    	myEnsemble->myLustig.ljUlus   = ulus + eLrc;
    	myEnsemble->myLustig.ljdUlus  = dudv + dUlrc;
    	myEnsemble->myLustig.ljd2Ulus = d2udv + d2Ulrc;
    }
    // myEnsemble->myEnergy.realEnergy = elecE;
    for (Int k = XX; k <= ZZ; k++)
    {
        myEnsemble->virial[k] -= virial[k];
        // myEnsemble->molVirial[k] -= molVirial[k];
    }
} 

// Beyond the cutoff, usually the repulsive term can be neglected.
// long_range_correct() will make correction to dispersion term 
// algorithm refering to GROMACS Appendix C
void MCYForce::long_range_correct()
{
    Int numTypes = myEnsemble->nAtomTypes;
    myAtoms = myEnsemble->atoms;
    Int *atomNum = new Int[numTypes];
    Double val = 0.0;
    double valU = 0.0, valdU = 0.0, vald2U = 0.0;
    double ljcutU = 0.0, ljcutdU = 0.0, ljcutd2U = 0.0;
    double sigma6, cutoff3, volrc, NNval, csix,c6excl,invvol;
    double a1,b1,a2,b2;
    int    npair = 0, nexcl = 0;

    for (Int i = 0; i < numTypes; i++)
        atomNum[i] = 0;

    // count atom number for each type
    for (Int i = 0; i < numAtoms; i++)
    {   
        int type = atoms[i].atomType;
        atomNum[type] += 1;
    }

    cutoff3 = cutOff*cutOff*cutOff;                                        //For Lustig 2011 
    volrc   = myEnsemble->boxLx*myEnsemble->boxLy*myEnsemble->boxLz;
    invvol  = 1.0/volrc; 

    for (Int i = 0; i < numTypes; i++)
    {
        for (Int j = i; j < numTypes; j++)
        {
            int mcyid = i*numTypes + j;
            double Aa1 = mcyParams[mcyid].AA1;
            double Bb1 = mcyParams[mcyid].BB1;
            double Aa2 = mcyParams[mcyid].AA2;
            double Bb2 = mcyParams[mcyid].BB2;
	    int npair_ij;

	   if (i != j) {
		npair_ij = atomNum[i]*atomNum[j];
	   } else {
		npair_ij = atomNum[i]*(atomNum[i] - 1)*0.5;
	   }

	   a1  += npair_ij*Aa1;
	   b1  += npair_ij*Bb1;
	   a2  += npair_ij*Aa2;
	   b2  += npair_ij*Bb2;
	   npair += npair_ij;
	}
    }

  cout << "  b2  =  "<< setprecision(6) << setw(14) << b2 << "	npair = "<< npair << endl;

   for (int i = 0; i < numAtoms; i++){
	for (int j = i; j < numAtoms; j++){
   	    if (myEnsemble->exclusion_check(i, j)){
                int mcyid = i*numTypes + j;
                double Aa1 = mcyParams[mcyid].AA1;
                double Bb1 = mcyParams[mcyid].BB1;
                double Aa2 = mcyParams[mcyid].AA2;
                double Bb2 = mcyParams[mcyid].BB2;

		a1 -= Aa1;
		b1 -= Bb1;
		a2 -= Aa2;
		b2 -= Bb2;
		nexcl += 1;
	    }
	}
   }

  cout << "  b2 =  "<< setprecision(6) << setw(14) << b2 << "	npair = "<< npair << "	nexcl =" << nexcl << endl;

   if (npair - nexcl <= 0) {
	a1 = 0;
	b1 = 0;
	a2 = 0;
	b2 = 0;
   }
   else {
	a1 /= npair - nexcl;
	b1 /= npair - nexcl;
	a2 /= npair - nexcl;
	b2 /= npair - nexcl;
   }

  cout << "  Csix_final =" << b2 << endl;


   ljcutU   =  2.0*PI*numAtoms*numAtoms*invvol*(a1*exp(-b1*cutOff)*(b1*b1*cutOff*cutOff + 2.0*b1*cutOff + 2.0)/(b1*b1*b1) - a2*exp(-b2*cutOff)*(b2*b2*cutOff*cutOff + 2.0*b2*cutOff + 2.0)/(b2*b2*b2)); 
   ljcutdU  = -2.0*PI*numAtoms*numAtoms*invvol*invvol*(a1*cutoff3*exp(-b1*cutOff) - a2*cutoff3*exp(-b2*cutOff))/(3.0) - ljcutU*invvol;
   ljcutd2U = 0.0; //-4.0*PI*numAtoms*numAtoms*invvol*invvol*invvol*csix/cutoff3;

   eLrc   = ljcutU;
   dUlrc  = ljcutdU;                   // Lustig 2011 LRC
   d2Ulrc = ljcutd2U;
 	
/*
            val += atomNum[i]*atomNum[j]*eps4*sigma2*sigma2*sigma2;    //Gromacs value
	    if (myEnsemble->myConfig->compute_lustig())
	    {
            	sigma6   = sigma2*sigma2*sigma2;                                        //For Lustig 2011 
            	ljcutU   = npair_ij*eps4*sigma6;
            	ljcutdU  = eps4*(4*sigma6*sigma6/(cutoff3*cutoff3*cutoff3) - 6*sigma6/cutoff3);
            	ljcutd2U = eps4*(20*sigma6*sigma6/(cutoff3*cutoff3*cutoff3) - 18*sigma6/cutoff3);

	    	NNval  = (atomNum[i] + atomNum[j])*(atomNum[i] + atomNum[j])*0.250 - (atomNum[i] + atomNum[j])*0.5;

cout << "  atomNum[i]  =  " << atomNum[i] <<  "  atomNum[j]  =  "<< atomNum[j] << " NNval  =  " << NNval << endl; 

            	valU   += NNval*ljcutU;
            	valdU  += NNval*ljcutdU;
            	vald2U += NNval*ljcutd2U;
            }
        }
    }
    eLrc   = -2*PI*val/(3*myEnsemble->boxLx*myEnsemble->boxLy*myEnsemble->boxLz*cutOff*cutOff*cutOff);     //Gromacs LRC
    if (myEnsemble->myConfig->compute_lustig())
    {
	eLrc   =  2*PI*valU/(9*volrc);                          // Lustig 2011 LRC
    	dUlrc  = -2*PI*valdU/(9*volrc*volrc);                   // Lustig 2011 LRC
    	d2Ulrc =  2*PI*vald2U/(9*volrc*volrc*volrc);            // Lustig 2011 LRC
    }
*/
cout << "  eLrc  =  " << eLrc <<  " dUlrc  =  "<< dUlrc << " d2Ulrc  =  " << d2Ulrc << endl;

    vLrc = ljcutdU*volrc;          
}

//-----------------------------------------------------------------------------------------------
// Pair table for MCY potential obtained as
// atomType1 * numtypes + atomType2
// Original Parameters from J. Chem. Phys., Vol. 64, No.4, 15 February 1976
// MC simulations J. Chem. Phys., Vol. 64, No.6, 15 March 1976
//-----------------------------------------------------------------------------------------------

void MCYForce::set_mcy_param()
{
        mcyParams[0].atomType1 = 0;
        mcyParams[0].atomType2 = 0;
        mcyParams[0].AA1 = a2mcy;
        mcyParams[0].BB1 = b2mcy;
        mcyParams[0].AA2 = 0.0;
        mcyParams[0].BB2 = 0.0;

        mcyParams[1].atomType1 = 0;
        mcyParams[1].atomType2 = 1;
        mcyParams[1].AA1 = a3mcy;
        mcyParams[1].BB1 = b3mcy;
        mcyParams[1].AA2 = a4mcy;
        mcyParams[1].BB2 = b4mcy;

        mcyParams[2].atomType1 = 1;
        mcyParams[2].atomType2 = 0;
        mcyParams[2].AA1 = a3mcy;
        mcyParams[2].BB1 = b3mcy;
        mcyParams[2].AA2 = a4mcy;
        mcyParams[2].BB2 = b4mcy;

        mcyParams[3].atomType1 = 1;
        mcyParams[3].atomType2 = 1;
        mcyParams[3].AA1 = a1mcy;
        mcyParams[3].BB1 = b1mcy;
        mcyParams[3].AA2 = 0.0;
        mcyParams[3].BB2 = 0.0;

}

void MCYForce::write_force_info(ofstream& of)
{
    of << "MCY Force" << endl;
    of << "cutoff2: " << cutOff2 << endl;
    if (switchOn)
        of << "compute MCY force with switch function" << endl;
    if (computeCoulomb)
        of << "Non-bonded force include coulomb" << endl;
    of << "MCY parameters:" << endl;
    of << "type1" << '\t' << "type2" << '\t' << "a1" << '\t' << " b1" << '\t' << "  a2" << '\t' << "  b2" << endl;
    for (Int i = 0; i <myEnsemble->nAtomTypes; i++)
    {
        for (Int j = 0; j < myEnsemble->nAtomTypes; j++)
        {
	    int mcyid = i*numAtomTypes + j;
            of << i << '\t' << j << '\t' << mcyParams[mcyid].AA1 << '\t' << mcyParams[mcyid].BB1 << '\t' <<  mcyParams[mcyid].AA2 << '\t' << mcyParams[mcyid].BB2  <<endl;
        }
    }
    of << "Long range correct energy: " << eLrc << endl;        
    of << endl;
}

void MCYForce::write_energy(ofstream& of)
{
    of << "MCY Energy: " << energy << endl;
}

// JC string LJForce::get_force_id() 
// JC {
// JC  return LJForce::forceIdentifier;
// JC }
