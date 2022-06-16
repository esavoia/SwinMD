/** BBVForce.cpp -- 
 **
 ** Copyright (C) 2020
 ** Centre for Molecular Simulation (CMS)
 ** School of Information Technology
 ** Swinburne University of Technology
 ** PO Box 21 Hawthorn, Vic 3122, Australia
 **
 ** Author: eltoyar
 ** Email: eoyarzua@.swin.edu.au
 **/

#include "BBVForce.h" 


BBVForce::BBVForce(Ensemble* ensemble) : Force(ensemble)
{
    #ifdef DEBUG
        DEBUGMSG("Creating BBV Force ");
    #endif

    bbvParams = NULL;
    virt1 = NULL;
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


    bbvParams = new BBVParam[numAtomTypes*numAtomTypes];
          if (bbvParams  == NULL)
              ERRORMSG( "memory allocation of BBV Parameters error");

    set_bbv_param();

    if (myEnsemble->myConfig->use_long_correct())
        long_range_correct();           // calculate once at the beginning

    write_force_info(*sysdataFile);
}

// this implementation of compute() might be not efficient in parallelisation
// in case of OpenMP, data dependency may occur in the summation of force
// in case of MPI, many communication may be involved 
void BBVForce::compute()
{

    Atom *atomi, *atomj;
    Atom* atom;
    Molecule *moli, *molj;
    Vector3 rij, fij, ri, rj;
    Double r, r_1, rij2, rij_2, sr2, sr6, sr12, tmpE, tmpf, rij6;
    double dulj, d2ulj, voll, vol2, polbbv;
    double expbbv, invexpbbv, expbbv2, expbbv8, expbbv15;
    Double qi, qij;
    Double rcut_2 = 1.0/cutOff2;
    Double elecE = 0.0;

    Vector3 rim, rjm, rijm;
    Double  molVirial[9];
    
    energy = 0.0;
    ulus   = 0.0;
    dudv   = 0.0;
    dudvq  = 0.0;
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
        Int size = atomi->get_list_size();          // size of pairlist
        qi = atomi->scaledCharge;
  	ri = atomi->position;

        // go through the pairlist of atomi
        for (Int j = 0; j < size; j++)   /// Int j = i+1; j<numAtoms; j++
        {
            atomj = atomi->myPairList[j];
//	if (atoms[i]->molID != atoms[j]->molID ){
	if (atomi->molID != atomj->molID){
            id = atomj->molID;    //
	    rj = atomj->position;
            rij = ri - rj;
            rijm = myEnsemble->molecules[id].massCenter - rim;  //
            rijm -= rij;   //

            myEnsemble->apply_pbc(rij);
            rijm += rij;
            r = rij.length();
            rij2 = r*r;

            if (rij2 <= cutOff2)
            {
		int bbvid = atomi->atomType*numAtomTypes + atomj->atomType;
		double Aa1 = bbvParams[bbvid].Aa1;
		double Aa2 = bbvParams[bbvid].Aa2;
		double Aa3 = bbvParams[bbvid].Aa3;
		double Aa4 = bbvParams[bbvid].Aa4;
		double Aa5 = bbvParams[bbvid].Aa5;
		double Aa6 = bbvParams[bbvid].Aa6;

                r_1 = 1.0/r;
                rij_2 = r_1*r_1;
		rij6 = rij_2*rij_2*rij_2;
		expbbv = 1.0 + exp(-2.0*(r/Aa0 - 2.0));
		invexpbbv = 1.0/expbbv;
		expbbv2 = invexpbbv*invexpbbv ;
		expbbv8 = expbbv2*expbbv2*expbbv2*expbbv2;
		expbbv15 = expbbv8*expbbv2*expbbv2*expbbv2*invexpbbv;
		polbbv  = (Aa3*rij6*rij6 + Aa4*rij6*rij_2*rij_2 + Aa5*rij6*rij_2 + Aa6*rij6);
		tmpE = Aa1*exp(-Aa2*r) + polbbv*expbbv15;
		//tmpf = -Aa1*Aa2*exp(-Aa2*r) - (12.0*Aa3*rij6*rij6 + 10.0*Aa4*rij6*rij_2*rij_2 + 8.0*Aa5*rij6*rij_2 + 6.0*Aa6*rij6)*r_1;
		tmpf = -Aa1*Aa2*r_1*exp(-Aa2*r) + (30/Aa0)*r_1*polbbv*expbbv8*expbbv8*exp(4.0 - 2.0*r/Aa0) - (12.0*Aa3*rij6*rij6 + 10.0*Aa4*rij6*rij_2*rij_2 + 8.0*Aa5*rij6*rij_2 + 6.0*Aa6*rij6)*rij_2*expbbv15;

//cout << atomi->atomType << "    " << atomj->atomType << "    " << Aa5 << endl;

		if (myEnsemble->myConfig->compute_lustig())
		{
		dulj  = rij2*tmpf;                  		// Lustig 2011 Eq. 14b
		d2ulj = -2.0*rij2*tmpf;	        // Lustig 2011 Eq. 14c
		}

                if(switchOn && (rij2 > switchDist))
                {
                    Double c = cutOff2 - rij2;
                    Double c4 = c*(2*rij2 + c2);
                                        
                    Double switchVal = c*c4*c1;
                    Double dSwitchVal = c3*(c*c-c4);
                    tmpf = tmpf*switchVal - tmpE*dSwitchVal;
                    tmpE *= switchVal;
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
                    Double e = qij*r_1*expbbv15;
                    Double f = (30/Aa0)*e*invexpbbv*exp(4.0 - 2.0*r/Aa0)*r_1 - e*rij_2;
                    // use GROMACS shift function
                    Double r_rcut = r/cutOff;
                    Double r4_rcut4 = rij2*rij2*rcut_2*rcut_2;
                    elecE += e*(1 - 1.666667*r_rcut + 1.666667*r4_rcut4 - r_rcut*r4_rcut4);
                    tmpf += f*(1 - 5*r4_rcut4 + 4*r_rcut*r4_rcut4);
		    dudvq += f*(1 - 5*r4_rcut4 + 4*r_rcut*r4_rcut4)*rij2;
                    // CHARMM's shift function
                    // elecE += e*(1 - rij2*rcut_2)*(1-rij2*rcut_2);
                    // tmpf += f*(1 - r/cutOff)*(1 - r/cutOff);
                } 

                fij = -tmpf*rij;  
 
		atomi->force += fij;
               	atomj->force -= fij;

                virial[XX] += - fij.x * rij.x;
                virial[XY] += - fij.y * rij.x;
                virial[XZ] += - fij.z * rij.x;
                virial[YX] += - fij.x * rij.y;
                virial[YY] += - fij.y * rij.y;
                virial[YZ] += - fij.z * rij.y;
                virial[ZX] += - fij.x * rij.z;
                virial[ZY] += - fij.y * rij.z;
                virial[ZZ] += - fij.z * rij.z;

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
    	dudv   = dudv/(3*voll);
    	d2udv  = d2udv/(9*vol2);

    	myEnsemble->myLustig.ljUlus   = ulus + eLrc + elecE;
    	myEnsemble->myLustig.ljdUlus  = dudv + dUlrc;
	if (computeCoulomb) {
	    myEnsemble->myLustig.realDeriv = dudvq; }
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
void BBVForce::long_range_correct()
{
    Int numTypes = myEnsemble->nAtomTypes;
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
	    int bbvid = i*numTypes + j;
	    double Aa = bbvParams[bbvid].Aa6;
            int npair_ij;

	   cout << " i = " << i << "    j = " << j << "    a6_ij = "<< Aa << endl;

           if (i != j) {
                npair_ij = atomNum[i]*atomNum[j];
           } else {
                npair_ij = atomNum[i]*(atomNum[i] - 1)*0.5;
           }

           csix  += npair_ij*Aa;
           npair += npair_ij;
        }
    }

  cout << "  Csix =  "<< setprecision(6) << setw(14) << csix << "       npair = "<< npair << endl;

   for (int i = 0; i < numAtoms; i++){
        for (int j = i; j < numAtoms; j++){
            if (myEnsemble->exclusion_check(i, j)){
		int bbvid = atoms[i].atomType*numTypes + atoms[j].atomType;
                double Aa = bbvParams[bbvid].Aa6;

                csix -= Aa;
                nexcl += 1;
            }
        }
   }

  cout << "  Csix =  "<< setprecision(6) << setw(14) << csix << "       npair = "<< npair << "  nexcl =" << nexcl << endl;

   if (npair - nexcl <= 0) {
        csix = 0;
   }
   else {
        csix /= npair - nexcl;
   }

  cout << "  Csix_final =" << csix << endl;


   ljcutU   =  2.0*PI*numAtoms*numAtoms*invvol*csix/(3.0 * cutoff3);
   ljcutdU  = -4.0*PI*numAtoms*numAtoms*invvol*invvol*csix/(3.0 * cutoff3);
   ljcutd2U = -4.0*PI*numAtoms*numAtoms*invvol*invvol*invvol*csix/cutoff3;

   eLrc   = ljcutU;
   dUlrc  = ljcutdU;                   // Lustig 2011 LRC
   d2Ulrc = -ljcutd2U;

   vLrc = -2.0*eLrc;		// It's divided by V after in Integrator.cpp

   cout << "  eLrc  =  " << eLrc <<  " dUlrc  =  "<< dUlrc << " d2Ulrc  =  " << d2Ulrc << endl;          
}

//-------------------------------------------------------------------------------------------
// Pair table for Buckingham Parameters obtained as 
// atomType1 * numtypes + atomType2
// Constants from Tsuzuki et al. J. Phys. Chem., Vol. 100, No. 11, 1996 
// model 3 of Table 2 in the article
// -----------------------------------------------------------------------------------------

void BBVForce::set_bbv_param()                              
{
	bbvParams[0].atomType1 = 0;
	bbvParams[0].atomType2 = 0;
	bbvParams[0].Aa1 = a1oo;
	bbvParams[0].Aa2 = a2oo;
	bbvParams[0].Aa3 = a3oo;
	bbvParams[0].Aa4 = a4oo;
	bbvParams[0].Aa5 = a5oo;
	bbvParams[0].Aa6 = a6oo;

	bbvParams[1].atomType1 = 0;
	bbvParams[1].atomType2 = 1;
	bbvParams[1].Aa1 = a1co;
	bbvParams[1].Aa2 = a2co;
	bbvParams[1].Aa3 = a3co;
	bbvParams[1].Aa4 = a4co;
	bbvParams[1].Aa5 = a5co;
	bbvParams[1].Aa6 = a6co;

	bbvParams[2].atomType1 = 0;
	bbvParams[2].atomType2 = 2;
	bbvParams[2].Aa1 = a1bo;
	bbvParams[2].Aa2 = a2bo;
	bbvParams[2].Aa3 = a3bo;
	bbvParams[2].Aa4 = a4bo;
	bbvParams[2].Aa5 = a5bo;
	bbvParams[2].Aa6 = a6bo;

	bbvParams[3].atomType1 = 1;
	bbvParams[3].atomType2 = 0;
	bbvParams[3].Aa1 = a1co;
	bbvParams[3].Aa2 = a2co;
	bbvParams[3].Aa3 = a3co;
	bbvParams[3].Aa4 = a4co;
	bbvParams[3].Aa5 = a5co;
	bbvParams[3].Aa6 = a6co;

	bbvParams[4].atomType1 = 1;
	bbvParams[4].atomType2 = 1;
	bbvParams[4].Aa1 = a1cc;
	bbvParams[4].Aa2 = a2cc;
	bbvParams[4].Aa3 = a3cc;
	bbvParams[4].Aa4 = a4cc;
	bbvParams[4].Aa5 = a5cc;
	bbvParams[4].Aa6 = a6cc;

	bbvParams[5].atomType1 = 1;
	bbvParams[5].atomType2 = 2;
	bbvParams[5].Aa1 = a1bc;
	bbvParams[5].Aa2 = a2bc;
	bbvParams[5].Aa3 = a3bc;
	bbvParams[5].Aa4 = a4bc;
	bbvParams[5].Aa5 = a5bc;
	bbvParams[5].Aa6 = a6bc;

	bbvParams[6].atomType1 = 2;
	bbvParams[6].atomType2 = 0;
	bbvParams[6].Aa1 = a1bo;
	bbvParams[6].Aa2 = a2bo;
	bbvParams[6].Aa3 = a3bo;
	bbvParams[6].Aa4 = a4bo;
	bbvParams[6].Aa5 = a5bo;
	bbvParams[6].Aa6 = a6bo;

	bbvParams[7].atomType1 = 2;
	bbvParams[7].atomType2 = 1;
	bbvParams[7].Aa1 = a1bc;
	bbvParams[7].Aa2 = a2bc;
	bbvParams[7].Aa3 = a3bc;
	bbvParams[7].Aa4 = a4bc;
	bbvParams[7].Aa5 = a5bc;
	bbvParams[7].Aa6 = a6bc;

	bbvParams[8].atomType1 = 2;
	bbvParams[8].atomType2 = 2;
	bbvParams[8].Aa1 = a1bb;
	bbvParams[8].Aa2 = a2bb;
	bbvParams[8].Aa3 = a3bb;
	bbvParams[8].Aa4 = a4bb;
	bbvParams[8].Aa5 = a5bb;
	bbvParams[8].Aa6 = a6bb;

}

void BBVForce::write_force_info(ofstream& of)
{
    of << "BBV Force" << endl;
    of << "cutoff2: " << cutOff2 << endl;
    if (switchOn)
        of << "compute BBV force with switch function" << endl;
    if (computeCoulomb)
        of << "Non-bonded force include coulomb" << endl;
    of << "BBV parameters:" << endl;
    of << "type1" << '\t' << "type2" << '\t' << "a1" << '\t' << "a2" << '\t' << "    a3" <<  '\t' << "    a4" << '\t' << "    a5" <<  '\t' << "    a6" << endl;
    for (Int i = 0; i <myEnsemble->nAtomTypes; i++)
    {
        for (Int j = 0; j < myEnsemble->nAtomTypes; j++)
        {
	    int buckid = i*numAtomTypes + j;

            of << i << '\t' << j << '\t' << bbvParams[buckid].Aa1 << '\t' << bbvParams[buckid].Aa2 <<  '\t' << bbvParams[buckid].Aa3 <<  '\t' << bbvParams[buckid].Aa4;
	    of  <<  '\t' << bbvParams[buckid].Aa5 <<  '\t' << bbvParams[buckid].Aa6 << endl;
        }
    }
    of << "Long range correct energy: " << eLrc << endl;        
    of << endl;
}

void BBVForce::write_energy(ofstream& of)
{
    of << "BBV Energy: " << energy << endl;
}

// JC string LJForce::get_force_id() 
// JC {
// JC  return LJForce::forceIdentifier;
// JC }
