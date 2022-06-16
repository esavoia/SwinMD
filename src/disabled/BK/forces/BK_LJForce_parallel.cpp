/** BuffvdWForce.cpp -- 
 **
 ** Copyright (C) 2003
 ** Centre for Molecular Simulation (CMS)
 ** School of Information Technology
 ** Swinburne University of Technology
 ** PO Box 21 Hawthorn, Vic 3122, Australia
 **
 ** Author: Elton Oyarzua
 ** Email: eoyarzua@swin.edu.au
 ** Based on Force.h by Zhongwu Zhou
 **/

#include "BuffvdWForce.h" 


BuffvdWForce::BuffvdWForce(Ensemble* ensemble) : Force(ensemble)
{
    #ifdef DEBUG
        DEBUGMSG("Creating BuffvdWForce ");
    #endif

//JC    forceIdentifier = "LJForce";// in order to comply with the cluster compiler
    numAtoms = myEnsemble->nAtoms;
    numMols = myEnsemble->nMols;
    params = myEnsemble->myParams;
    cutOff = myEnsemble->myConfig->get_cutoff();
    cutOff2 = cutOff*cutOff;
    switchDist = cutOff2;           
    switchOn = myEnsemble->myConfig->is_switch_on();
    useConstraint = myEnsemble->myConfig->is_constraint_on();
    statEnsem = myEnsemble->myConfig->get_ensemble_status();

    virt = NULL;

    virt = new Vector3 [numAtoms];

    if ((virt==NULL))
        ERRORMSG("memory allocation error for BuffvdWForce ");


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
    if (myEnsemble->myConfig->use_long_correct())
        long_range_correct();           // calculate once at the beginning at constant Volume
    write_force_info(*sysdataFile);

    MPI_Comm_rank(MPI_COMM_WORLD, &prank);
    MPI_Comm_size(MPI_COMM_WORLD, &psize);
}

// this implementation of compute() might be not efficient in parallelisation
// in case of OpenMP, data dependency may occur in the summation of force
// in case of MPI, many communication may be involved 
void BuffvdWForce::compute()
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

    if (statEnsem == 1){
        if (myEnsemble->myConfig->use_long_correct()){
            long_range_correct();       }
    }
    
    for (Int i = XX; i <= ZZ; i++)
    {
        virial[i] = 0.0;
        molVirial[i] = 0.0;
    }

    // go through each atom
   int numDecompose = (int) numAtoms/psize;
   if ((numAtoms % psize) == 0) numDecompose--;

   for(int i = 0; i <= numDecompose;i++)
   {
     int atomIndex = i*psize + prank;
     if(atomIndex <numAtoms )
     {
        atomi = &atoms[atomIndex]; 
        Int id = atomi->molID;
        rim = myEnsemble->molecules[id].massCenter; //
        Int size = atomi->get_list_size();      // size of pairlist
        qi = atomi->scaledCharge;

        // go through the pairlist of atomi
        for (Int j = 0; j < size; j++)   /// Int j = i+1; j<numAtoms; j++
        {
            atomj = atomi->myPairList[j];
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
                LJPairParam ljParam = params->get_lj_parameter(atomi->atomType, atomj->atomType);
                Double sigma2 = ljParam.sigma;
                Double eps4 = ljParam.eps;                

//if (atomi->atomType == 1 && atomj->atomType == 1){
//cout << rij << " <- atom      Rij ->   " << rijm << endl;}

                r_1 = 1.0/r;
                rij_2 = r_1*r_1;
                sr2 = sigma2*rij_2;
                sr6 = sr2*sr2*sr2;
                sr12 = sr6*sr6;
                tmpE = eps4*(sr12 - sr6);
                tmpf = (sr12 + sr12 - sr6)*6*eps4*rij_2;
		double rijRij;
		if (useConstraint){
		    rijRij = rij*rijm;
		    //cout << "Constraint!!!!"  << endl;
		}
		else {
		    rijRij = rij2;
		    //cout << "APAGADO!!" << endl;
		}

		if (myEnsemble->myConfig->compute_lustig())
		{
		dulj  = tmpf*rijRij;                  // Lustig 2011 Eq. 15b
		d2ulj = eps4*(156.0*sr12 - 42.0*sr6)*rijRij*rijRij*rij_2*rij_2 + 2.0*dulj;	        // Lustig 2011 Eq. 15c
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
		    
		    double sr2cut   = sigma2*rcut_2;
		    double sr6cut   = sr2cut*sr2cut*sr2cut;
		    double sr12cut  = sr6cut*sr6cut;
		    double Eljcut   = eps4*(sr12cut - sr6cut);
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

            }  // if-sentence cut-off
          }  // if-different molIDs
	} // for j-loop       
      } // if sentence i-loop
    } // parallel decomposition 

   double temp;
   MPI_Reduce(&energy,&temp,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   energy = temp;

   temp = 0.0;
   MPI_Allreduce(&dudv,&temp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   dudv = temp;

   temp = 0.0;
   MPI_Reduce(&d2udv,&temp,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   d2udv = temp;


    //------------------------------------------------------------------
    // Add the long-range correction to the virial only at the root
    // processor. Remember that AllReduce of the virial is 
    // performed in LFIntegrator.cpp: WE DON'T want to add to
    // the virial psize-times corrections!!! 
    //------------------------------------------------------------------
    if (prank == 0) {
        virial[XX] += vLrc;
        virial[YY] += vLrc;
        virial[ZZ] += vLrc;
    }
    myEnsemble->myEnergy.ljEnergy = energy + eLrc;

    if (myEnsemble->myConfig->compute_lustig())
    {
    	dudv   = -dudv/(3.0*voll);
    	d2udv  = d2udv/(9.0*vol2);

//cout << "  dudv  =  " << dudv << "  d2udv =  " <<d2udv << endl;

    	myEnsemble->myLustig.ljUlus   = energy + eLrc;
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

//---------------------------------------------------------------------------------------------------
// Virtual position of Hydrogens following the "reduction factor" of
// the vdW interactions of the AMOEBA force field.
//
// Permutation of distances based from the Lectures of 
// Berendsen, H. J. C., van Gunsteren, W. F.  Molecular dynamics simulations: Techniques and
// approaches. In: Molecular Liquids-Dynamics and Interactions. et al., A. J. B. ed. NATO
// ASI C 135. Reidel Dordrecht, The Netherlands 1984 475–500. Figure 8.
//
// This method is computed in serial: All the processors must know 
// each virtual position and this method it is not time consuming.
//---------------------------------------------------------------------------------------------------

void BuffvdWForce::VirtualPositions()
{
    double factor;
 ;
}
// Beyond the cutoff, usually the repulsive term can be neglected.
// long_range_correct() will make correction to dispersion term 
// algorithm refering to GROMACS Appendix C
void BuffvdWForce::long_range_correct()
{
    Int numTypes = myEnsemble->nAtomTypes;
    myAtoms = myEnsemble->atoms;
    Int *atomNum = new Int[numTypes];
    Double val = 0.0;
    double valU = 0.0, valdU = 0.0, vald2U = 0.0;
    double ljcutU = 0.0, ljcutdU = 0.0, ljcutd2U = 0.0;
    double sigma6, cutoff3, volrc, NNval, csix,c6excl,invvol;
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
            LJPairParam ljParam = params->get_lj_parameter(i, j);
            Double sigma2 = ljParam.sigma;
            Double eps4 = ljParam.eps;
            sigma6  = sigma2*sigma2*sigma2;                                        //
	    int npair_ij;

	   if (i != j) {
		npair_ij = atomNum[i]*atomNum[j];
	   } else {
		npair_ij = atomNum[i]*(atomNum[i] - 1)*0.5;
	   }

	   csix  += npair_ij*eps4*sigma6;
	   npair += npair_ij;
	}
    }


   for (int i = 0; i < numAtoms; i++){
	for (int j = i; j < numAtoms; j++){
   	    if (myEnsemble->exclusion_check(i, j)){
		LJPairParam ljParam = params->get_lj_parameter(atoms[i].atomType, atoms[j].atomType);
                Double sigma2 = ljParam.sigma;
                Double eps4 = ljParam.eps;
                sigma6  = sigma2*sigma2*sigma2;

		csix -= eps4*sigma6;
		nexcl += 1;
	    }
	}
   }


   if (npair - nexcl <= 0) {
	csix = 0;
   }
   else {
	csix /= npair - nexcl;
   }


   ljcutU   = -2.0*PI*numAtoms*numAtoms*invvol*csix/(3.0 * cutoff3); 
   ljcutdU  =  4.0*PI*numAtoms*numAtoms*invvol*invvol*csix/(3.0 * cutoff3);
   ljcutd2U = -4.0*PI*numAtoms*numAtoms*invvol*invvol*invvol*csix/cutoff3;

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
    if (statEnsem == 0){
	cout << "  eLrc  =  " << eLrc <<  " dUlrc  =  "<< dUlrc << " d2Ulrc  =  " << d2Ulrc << endl;
    }

    vLrc = -2.0*eLrc;          
}

void BuffvdWForce::write_force_info(ofstream& of)
{
    of << "LJ Force" << endl;
    of << "cutoff2: " << cutOff2 << endl;
    if (switchOn)
        of << "compute LJ force with switch function" << endl;
    if (computeCoulomb)
        of << "Non-bonded force include coulomb" << endl;
    of << "LJ parameters:" << endl;
    of << "type1" << '\t' << "type2" << '\t' << "eps4" << '\t' << "sigma2" << endl;
    for (Int i = 0; i <myEnsemble->nAtomTypes; i++)
    {
        for (Int j = 0; j < myEnsemble->nAtomTypes; j++)
        {
            LJPairParam ljParam = params->get_lj_parameter(i, j);
            of << i << '\t' << j << '\t' << ljParam.eps << '\t' << ljParam.sigma << endl;
        }
    }
    of << "Long range correct energy: " << eLrc << endl;        
    of << endl;
}

void BuffvdWForce::write_energy(ofstream& of)
{
    of << "LJ Energy: " << energy << endl;
}

// JC string LJForce::get_force_id() 
// JC {
// JC  return LJForce::forceIdentifier;
// JC }
