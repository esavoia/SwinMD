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
    myMols = myEnsemble->molecules;
    params = myEnsemble->myParams;
    numTypes = myEnsemble->nAtomTypes;

    //-------------------------------------------------------------------------
    // Get parameters from SimConfiguration.h
    //-------------------------------------------------------------------------
    
    cutOff        = myEnsemble->myConfig->get_cutoff();
    cutOff2       = cutOff*cutOff;
    switchDist    = cutOff2;           
    switchOn      = myEnsemble->myConfig->is_switch_on();
    useConstraint = myEnsemble->myConfig->is_constraint_on();
    statEnsem     = myEnsemble->myConfig->get_ensemble_status();
    factor        = myEnsemble->myConfig->get_hfactor_reduction();
    n             = myEnsemble->myConfig->get_nBuff_param();
    m             = myEnsemble->myConfig->get_mBuff_param();
    delta         = myEnsemble->myConfig->get_deltaBuff_param();
    gamma         = myEnsemble->myConfig->get_gammaBuff_param();

    virt = NULL;
    buffpairs = NULL;

    virt = new Vector3 [numAtoms];
    buffpairs = new BuffTable[numTypes*numTypes];

    if ((virt==NULL)||(buffpairs==NULL))
        ERRORMSG("memory allocation error for BuffvdWForce ");


    if (switchOn)
        switchDist = myEnsemble->myConfig->get_switch_distance()*myEnsemble->myConfig->get_switch_distance();
	switchd1   = myEnsemble->myConfig->get_switch_distance();
    if (switchOn && (cutOff2 > switchDist))
    {
    /*    c1 = 1.0/(cutOff2 - switchDist);
        c1 = c1*c1*c1;
        c2 = cutOff2 - 3.0*switchDist;
        c3 = 4*c1;*/
	double offcut = cutOff - switchd1;
	denom = offcut*offcut*offcut*offcut*offcut;
	c0    = cutOff*cutOff2*( cutOff2 - 5.0*cutOff*switchd1 + 10.0*switchDist) / denom;
	c1    = -30.0*cutOff2*switchDist / denom;
	c2    = 30.0*(cutOff2*switchd1 + cutOff*switchDist) / denom;
	c3    = -10.0*(cutOff2 + 4.0*cutOff*switchd1 + switchDist) / denom;
	c4    = 15.0*(cutOff + switchd1) / denom;
	c5    = -6.0 / denom;
    }
    else
    {
        c1 = 0.0;
        c2 = 0.0;
        c3 = 0.0;
    }

    set_buff_param();

    computeCoulomb = myEnsemble->myConfig->compute_coulomb();
    if (myEnsemble->myConfig->use_long_correct())
        long_range_correct();           // calculate once at the beginning at constant Volume

    write_force_info(*sysdataFile);

    MPI_Comm_rank(MPI_COMM_WORLD, &prank);
    MPI_Comm_size(MPI_COMM_WORLD, &psize);
}

void BuffvdWForce::compute()
{

    Atom *atomi, *atomj;
    Molecule *moli, *molj;
    Vector3 ri, rj, rij, fij;
    Double r, r_1, rij2, tmpE, tmpf;
    double dulj, d2ulj, voll, vol2;
    int id, molid, index;
    double eps, sigma, rho;
    double par1, par2, rho_m, rho_m1;
    double expnm, par1_nm, invsigma, m_1;
    double rij3, rij4, rij5, switchVal, dSwitchVal;
    Double rcut_2 = 1.0/cutOff2;

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

    //-------------------------------------------------------------
    // First compute the vdW displacement of the Hydrogens
    // in serial at each processor. 
    //-------------------------------------------------------------
    
    VirtualPositions();

    // go through each atom
   int numDecompose = (int) numAtoms/psize;
   if ((numAtoms % psize) == 0) numDecompose--;

   for(int i = 0; i <= numDecompose;i++)
   {
     int atomIndex = i*psize + prank;
     if(atomIndex <numAtoms )
     {
        atomi = &atoms[atomIndex]; 
	id    = atomi->atomID;
        molid = atomi->molID;
        rim = myEnsemble->molecules[molid].massCenter; //
        Int size = atomi->get_list_size();      // size of pairlist
	ri = atomi->position; 

	if (atomi->atomType == 0){
	    ri =  virt[id];
	}

        // go through the pairlist of atomi
        for (Int j = 0; j < size; j++)   /// Int j = i+1; j<numAtoms; j++
        {
            atomj = atomi->myPairList[j];
	    if (atomi->molID == atomj->molID){
	        ERRORMSG(" BuffvdWForce: vdW interaction within the same molecule!! Check your PairList.");
            }  
            molid = atomj->molID;    
	    id    = atomj->atomID;
	    rj    = atomj->position;

	    if (atomj->atomType == 0){
		rj =  virt[id];
	    }

            rij = ri - rj;
            rijm = rim - myEnsemble->molecules[molid].massCenter;  
            myEnsemble->apply_pbc(rij);
            myEnsemble->apply_pbc(rijm);
            r = rij.length();
            rij2 = r*r;

            if (rij2 <= cutOff2)
            {
		index = atomi->atomType*numTypes + atomj->atomType;
		
		if ( (atomi->atomType != buffpairs[index].atomType1) || (atomj->atomType != buffpairs[index].atomType2) ){
	            ERRORMSG(" BuffvdWForce: The atomType count is wrong. Check pairTable for BuffvdW parameters or general atomTypes IDs");
		}

	        sigma    = buffpairs[index].sigma;
	        eps      = buffpairs[index].eps;

		if (eps == 0) {
		    tmpE = 0;
		    tmpf = 0;
		} else {
		    invsigma = 1.0/sigma;
		    m_1      = m - 1.0;
		    rho      = r*invsigma;
		    rho_m    = pow(rho, m);
		    rho_m1   = pow(rho, m_1);
		    expnm    = n - m;

                    r_1  = 1.0/r;
		    par1 = (1.0 + delta)/(rho + delta);
		    par2 = (1.0 + gamma)/(rho_m + gamma);
		    par1_nm = pow(par1, expnm);

                    tmpE = eps*par1_nm*(par2 - 2.0);
                    tmpf = -eps*expnm*par1_nm*(par2 - 2.0)*invsigma/(rho + delta) - eps*par1_nm*par2*m*rho_m1*invsigma/(rho_m + gamma);
		}

		double rijRij;
		if (useConstraint){
		    rijRij = rij*rijm;
		}
		else {
		    rijRij = rij2;
		}

		if (myEnsemble->myConfig->compute_lustig())
		{
		    d2ulj = 0.0;	        
		}

                if(switchOn && (rij2 > switchDist))
                {/*
                    Double c = cutOff2 - rij2;
                    Double c4 = c*(2*rij2 + c2);
                                        
                    Double switchVal = c*c4*c1;
                    Double dSwitchVal = c3*(c*c-c4);*/

		    rij3 = rij2*r;
                    rij4 = rij2*rij2;
                    rij5 = rij2*rij3;

                    switchVal  = c5*rij5 + c4*rij4 + c3*rij3 + c2*rij2 + c1*r + c0;
                    dSwitchVal = 5.0*c5*rij4 + 4.0*c4*rij3 + 3.0*c3*rij2 + 2.0*c2*r + c1;

                    tmpf = tmpf*switchVal + tmpE*dSwitchVal;
                    tmpE *= switchVal;
                }

		dulj    = tmpf*r;             
                energy += tmpE;   
    		dudv   += dulj;

		if (myEnsemble->myConfig->compute_lustig())
		{
    		    d2udv  += d2ulj;
		}
                
                fij = -tmpf*rij*r_1;   

		//----------------------------------------------------------
		//  Force redistribution due to relocation of Hydrogens
		//----------------------------------------------------------

		moli = &myMols[atomi->molID];
		molj = &myMols[atomj->molID];
		Atom *atomOi  = moli->myAtoms[1];
		Atom *atomOj  = molj->myAtoms[1];

		if (atomi->atomType == 0){
                    atomi->force  += factor*fij;
                    atomOi->force += (1.0 - factor)*fij;
		} else {
                    atomi->force += fij;
		}

		if (atomj->atomType == 0){
		    atomj->force  -= factor*fij;
		    atomOj->force -= (1.0 - factor)*fij;
		} else {
		    atomj->force -= fij;
		}

                virial[XX] += -fij.x * rij.x;
                virial[XY] += -fij.y * rij.x;
                virial[XZ] += -fij.z * rij.x;
                virial[YX] += -fij.x * rij.y;
                virial[YY] += -fij.y * rij.y;
                virial[YZ] += -fij.z * rij.y;
                virial[ZX] += -fij.x * rij.z;
                virial[ZY] += -fij.y * rij.z;
                virial[ZZ] += -fij.z * rij.z;

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

    dudv = dudv/(3.0*voll);

    myEnsemble->myEnergy.ljEnergy = energy + eLrc;
    myEnsemble->myLustig.ljUlus   = energy + eLrc;
    myEnsemble->myLustig.ljdUlus  = dudv + dUlrc;

    if (myEnsemble->myConfig->compute_lustig())
    {
    	d2udv  = d2udv/(9.0*vol2);
    	myEnsemble->myLustig.ljd2Ulus = d2udv + d2Ulrc;
    }

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
    int ns, molID, id;
    Molecule *mol;

    for (int i = 0; i<numMols; i++)
    {
        ns    = myEnsemble->molecules[i].numAtoms;
        molID = myEnsemble->molecules[i].molID;

        mol = &myMols[molID];
        Atom *atomO  = mol->myAtoms[1];
	
	if (atomO->atomType != 1){
	    ERRORMSG(" BuffvdWForce: atomType==1 is not the central Oxygen!!. I can't compute the virtual displacement of Hydrogens. ");
	}

        for (int i = 0; i < ns; i++)
        {
            Atom *atomi = mol->myAtoms[i];
	    id = atomi->atomID;

	    if (atomi->atomType == 0) {
	        virt[id] = factor*atomi->position + (1.0 - factor)*atomO->position;
	    }
	}
    }
}

//-------------------------------------------------------------------------------------------
// Pair table for Buff-vdW Parameters accessed by 
// atomType1 * numtypes + atomType2
//-------------------------------------------------------------------------------------------

void BuffvdWForce::set_buff_param() 
{
    int index;
    double sigmai, epsi;
    double sigmaj, epsj;
    double Rij, epsij, sigmai2, sigmaj2, epssqrt;
    DEBUGMSG("Creating Buff-vdW pairwise table");

    for (int i = 0; i < numTypes; i++)
    {
	epsi   = myEnsemble->myParams->atomParams[i].eps;
	sigmai = myEnsemble->myParams->atomParams[i].sigma;

	if (epsi == 0 && sigmai == 0){
    	    DEBUGMSG("WARNING!! BuffvdWForce: eps = 0 and Rij = 0, check your energies (just in case)");
	}

        for(int j = 0; j < numTypes; j++)
        {
	    epsj   = myEnsemble->myParams->atomParams[j].eps;
	    sigmaj = myEnsemble->myParams->atomParams[j].sigma;
	    index  = i*numTypes + j;

	    sigmai2 = sigmai*sigmai;
	    sigmaj2 = sigmaj*sigmaj;
	    epssqrt = sqrt(epsi) + sqrt(epsj);

	    if (epsi == 0 || epsj == 0) {
		Rij   = 0.0000001;
		epsij = 0;
	    } else {
	        Rij   = (sigmai2*sigmai + sigmaj2*sigmaj)/(sigmai2 + sigmaj2); 
	        epsij = (4.0*epsi*epsj)/(epssqrt*epssqrt);
	    }

	    buffpairs[index].atomType1 = i;
	    buffpairs[index].atomType2 = j;
	    buffpairs[index].sigma     = Rij;
	    buffpairs[index].eps       = epsij;
	}
    }
}

// Beyond the cutoff, usually the repulsive term can be neglected.
// long_range_correct() will make correction to dispersion term 
// algorithm refering to GROMACS Appendix C
//

void BuffvdWForce::potFunc(double &r, double &sigma, double &eps, double &tmpE, double &tmpf)
{

    double par1, par2, rho_m, rho_m1;
    double expnm, par1_nm, invsigma, m_1;
    double rho;

    if (eps == 0) {
        tmpE = 0;
        tmpf = 0;
    } else {
        invsigma = 1.0/sigma;
        m_1      = m - 1.0;
        rho      = r*invsigma;
        rho_m    = pow(rho, m);
        rho_m1   = pow(rho, m_1);
        expnm    = n - m;

        par1 = (1.0 + delta)/(rho + delta);
        par2 = (1.0 + gamma)/(rho_m + gamma);
        par1_nm = pow(par1, expnm);

        tmpE = eps*par1_nm*(par2 - 2.0);
        tmpf = -eps*expnm*par1_nm*(par2 - 2.0)*invsigma/(rho + delta) - eps*par1_nm*par2*m*rho_m1*invsigma/(rho_m + gamma);
    }

}

void BuffvdWForce::long_range_correct()
{
    Int numTypes = myEnsemble->nAtomTypes;
    myAtoms = myEnsemble->atoms;
    Int *atomNum = new Int[numTypes];
    double ljcutU = 0.0, ljcutdU = 0.0, ljcutd2U = 0.0;
    double sigma6, cutoff3, volrc, csix,invvol;
    int    index, npair = 0, nexcl = 0;
    double eps, sigma, r, rij2;
    int level = 20;
    int nsteps = (int) pow(2,level) - 1;
    double tmpE, tmpf, h, maxright = 5000;
    double Ulrc = 0.0, sumU = 0.0;
    double Plrc = 0.0, sumP = 0.0;

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

    h = (maxright - cutOff)/pow(2,level);

    for (Int i = 0; i < numTypes; i++)
    {
        for (Int j = i; j < numTypes; j++)
        {
	   int npair_ij;
	   sumU = 0.0;
	   sumP = 0.0;

	   if (i != j) {
		npair_ij = atomNum[i]*atomNum[j];
	   } else {
		npair_ij = atomNum[i]*(atomNum[i] - 1)*0.5;
	   }

           npair += npair_ij;
	   index = i*numTypes + j;
	
	   if ( (i != buffpairs[index].atomType1) || (j != buffpairs[index].atomType2) ){
	       ERRORMSG(" BuffvdWForce: The atomType count is wrong. Check pairTable for BuffvdW parameters or general atomTypes IDs");
	   }

	   sigma    = buffpairs[index].sigma;
	   eps      = buffpairs[index].eps;

	   r    = cutOff;
    	   rij2 = r*r;
	   tmpE = 0.0;
	   tmpf = 0.0;

           potFunc(r, sigma, eps, tmpE, tmpf);

	   sumU += 0.5*h*rij2*tmpE;
	   sumP += 0.5*h*rij2*r*tmpf;

           for (int k=1; k <= nsteps; k++)
            {
                r = cutOff + k*h;
    		rij2 = r*r;
	        tmpE = 0.0;
	        tmpf = 0.0;

		potFunc(r, sigma, eps, tmpE, tmpf);

	        sumU += h*rij2*tmpE;
	        sumP += h*rij2*r*tmpf;
	    }

	   r = maxright;
    	   rij2 = r*r;
	   tmpE = 0.0;
	   tmpf = 0.0;

	   potFunc(r, sigma, eps, tmpE, tmpf);

	   sumU += 0.5*h*rij2*tmpE;
	   sumP += 0.5*h*rij2*r*tmpf;

	   Ulrc += npair_ij*sumU;
	   Plrc += npair_ij*sumP;
	   
	}
    }

/*
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
*/

   if (npair <= 0) {
	Ulrc = 0;
	Plrc = 0;
   }
   else {
	Ulrc /= npair;
	Plrc /= npair;
   }


   ljcutU   = 2.0*PI*numAtoms*numAtoms*invvol*Ulrc; 
   ljcutdU  = 2.0*PI*numAtoms*numAtoms*invvol*invvol*Plrc/3.0;
   ljcutd2U = 0.0;

   eLrc   = ljcutU;
   dUlrc  = ljcutdU;                   
   d2Ulrc = ljcutd2U;
 	
    if (statEnsem == 0){
	cout << "  eLrc  =  " << eLrc <<  " dUlrc  =  "<< dUlrc << " d2Ulrc  =  " << d2Ulrc << endl;
    }

    vLrc = ljcutdU*volrc;  	// It's divided by V after in Integrator.cpp          
}

void BuffvdWForce::write_force_info(ofstream& of)
{
    of << "BuffvdW Force" << endl;
    of << "cutoff2: " << cutOff2 << endl;
    if (switchOn)
        of << "compute BuffvdW force with switch function" << endl;
    if (computeCoulomb)
        of << "Non-bonded force include coulomb" << endl;
    of << "BuffvdW parameters:" << endl;
    of << "type1" << '\t' << "type2" << '\t' << " sigma(Rij)" << '\t' << " eps " << endl;
    for (Int i = 0; i <myEnsemble->nAtomTypes; i++)
    {
        for (Int j = 0; j < myEnsemble->nAtomTypes; j++)
        {
	    int index =  i*numTypes + j; 
            of << buffpairs[index].atomType1 << '\t' << buffpairs[index].atomType2 << '\t' << buffpairs[index].sigma << '\t' << buffpairs[index].eps << endl;
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
