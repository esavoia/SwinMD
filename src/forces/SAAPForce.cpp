/** SAAPForce.cpp -- 
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

#include "SAAPForce.h" 
#include <cmath>

SAAPForce::SAAPForce(Ensemble* ensemble) : Force(ensemble)
{
    #ifdef DEBUG
        DEBUGMSG("Creating SAAP Force ");
    #endif

    saapParams = NULL;
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


    saapParams = new SAAPParam[numAtomTypes*numAtomTypes];
          if (saapParams  == NULL)
              ERRORMSG( "memory allocation of SAAP Parameters error");

    set_saap_param();

    if (myEnsemble->myConfig->use_long_correct())
        long_range_correct();           // calculate once at the beginning

    write_force_info(*sysdataFile);
}

// this implementation of compute() might be not efficient in parallelisation
// in case of OpenMP, data dependency may occur in the summation of force
// in case of MPI, many communication may be involved 
void SAAPForce::compute()
{

    Atom *atomi, *atomj;
    Atom* atom;
    Molecule *moli, *molj;
    Vector3 rij, fij, ri, rj;
    Double r, r_1, rij2, rij_2, sr2, sr6, sr12, tmpE, tmpf, rij6;
    double dulj, d2ulj, voll, vol2;
    double expSAAP, a5r6; 
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
		int saapid = atomi->atomType*numAtomTypes + atomj->atomType;
		double Aa0 = saapParams[saapid].Aa0;
		double Aa1 = saapParams[saapid].Aa1;
		double Aa2 = saapParams[saapid].Aa2;
		double Aa3 = saapParams[saapid].Aa3;
		double Aa4 = saapParams[saapid].Aa4;
		double Aa5 = saapParams[saapid].Aa5;
		double Aa6 = saapParams[saapid].Aa6;

                r_1 = 1.0/r;
                rij_2 = r_1*r_1;
		expSAAP = Aa0*r_1*exp(Aa1*r) + Aa2*exp(Aa3*r) + Aa4;
		a5r6    = 1.0 + Aa5*rij2*rij2*rij2;
		tmpE = expSAAP/a5r6;
		tmpf = (-Aa0*rij_2*exp(Aa1*r) + Aa0*Aa1*r_1*exp(Aa1*r) + Aa2*Aa3*exp(Aa3*r))*r_1/a5r6 - 6.0*Aa5*rij2*rij2*expSAAP/(a5r6*a5r6);

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
                    Double e = qij*r_1;
                    Double f = -e*rij_2;
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
void SAAPForce::long_range_correct()
{
    Int numTypes = myEnsemble->nAtomTypes;
    Int *atomNum = new Int[numTypes];
    Double val = 0.0;
    double valU = 0.0, valdU = 0.0, vald2U = 0.0;
    double ljcutU = 0.0, ljcutdU = 0.0, ljcutd2U = 0.0;
    double sigma6, cutoff3, volrc, NNval, c4, c5, csix,invvol;
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
	    int saapid = i*numTypes + j;
	    double Aa4 = saapParams[saapid].Aa4;
	    double Aa5 = saapParams[saapid].Aa5;
            int npair_ij;

	   cout << " i = " << i << "    j = " << j << "    a4_ij = "<< Aa4 << "    a5_ij = "<< Aa5 << endl;

           if (i != j) {
                npair_ij = atomNum[i]*atomNum[j];
           } else {
                npair_ij = atomNum[i]*(atomNum[i] - 1)*0.5;
           }

           c4    += npair_ij*Aa4;
           c5    += npair_ij*Aa5;
           csix  += npair_ij*Aa4/Aa5;
           npair += npair_ij;
        }
    }

  cout << "  < a4 > =  "<< setprecision(6) << setw(14) << c4 << "       npair = "<< npair << endl;

   for (int i = 0; i < numAtoms; i++){
        for (int j = i; j < numAtoms; j++){
            if (myEnsemble->exclusion_check(i, j)){
		int saapid = atoms[i].atomType*numTypes + atoms[j].atomType;
	        double Aa4 = saapParams[saapid].Aa4;
	        double Aa5 = saapParams[saapid].Aa5;

                c4 -= Aa4;
                c5 -= Aa5;
                csix -= Aa4/Aa5;
                nexcl += 1;
            }
        }
   }

  cout << "  < a4 > =  "<< setprecision(6) << setw(14) << c4 << "       npair = "<< npair << "  nexcl =" << nexcl << endl;

   if (npair - nexcl <= 0) {
        c4 = 0;
        c5 = 0;
        csix = 0;
   }
   else {
        c4 /= npair - nexcl;
        c5 /= npair - nexcl;
        csix /= npair - nexcl;
   }

    if(c5 <= 0){
        ERRORMSG("a5 constants should be greater than 0"); }

  cout << "  < a4 >_final =" << c4 << endl;

   double  roota5   = sqrt(c5);
   double  atanSaap = atan(roota5 * cutoff3 );

  cout << "  arctan(ra5 * rc3) = " << atanSaap << endl;

   ljcutU   = 2.0*PI*numAtoms*numAtoms*invvol*c4/(3.0 * roota5) * (0.5*PI - atanSaap);
   ljcutdU  = 2.0*PI*numAtoms*numAtoms*invvol*invvol*c4/(3.0) * (atanSaap/roota5 - cutoff3/(c5*cutoff3*cutoff3 + 1.0) - 0.5*PI/roota5);
   ljcutd2U =  2.0*PI*numAtoms*numAtoms*invvol*(c4/c5)/(3.0 * cutoff3);
   //ljcutd2U = -4.0*PI*numAtoms*numAtoms*invvol*invvol*invvol*csix/cutoff3;

   eLrc   = ljcutU;
   dUlrc  = ljcutdU;                   // 
   //d2Ulrc = -ljcutd2U;

   vLrc = volrc*ljcutdU;	// It's divided by V after in Integrator.cpp 

   cout << "  eLrc  =  " << eLrc <<  " dUlrc  =  "<< dUlrc << " elrcTest  =  " << ljcutd2U << endl;          
}

//-------------------------------------------------------------------------------------------
// Pair table for SAAP Parameters obtained as 
// atomType1 * numtypes + atomType2
// -----------------------------------------------------------------------------------------

void SAAPForce::set_saap_param()                              
{
	saapParams[0].atomType1 = 0;
	saapParams[0].atomType2 = 0;
	saapParams[0].Aa0 = sa0oo;
	saapParams[0].Aa1 = sa1oo;
	saapParams[0].Aa2 = sa2oo;
	saapParams[0].Aa3 = sa3oo;
	saapParams[0].Aa4 = sa4oo;
	saapParams[0].Aa5 = sa5oo;
	saapParams[0].Aa6 = sa6oo;

	saapParams[1].atomType1 = 0;
	saapParams[1].atomType2 = 1;
	saapParams[1].Aa0 = sa0co;
	saapParams[1].Aa1 = sa1co;
	saapParams[1].Aa2 = sa2co;
	saapParams[1].Aa3 = sa3co;
	saapParams[1].Aa4 = sa4co;
	saapParams[1].Aa5 = sa5co;
	saapParams[1].Aa6 = sa6co;

	saapParams[2].atomType1 = 0;
	saapParams[2].atomType2 = 2;
	saapParams[2].Aa0 = sa0bo;
	saapParams[2].Aa1 = sa1bo;
	saapParams[2].Aa2 = sa2bo;
	saapParams[2].Aa3 = sa3bo;
	saapParams[2].Aa4 = sa4bo;
	saapParams[2].Aa5 = sa5bo;
	saapParams[2].Aa6 = sa6bo;

	saapParams[3].atomType1 = 1;
	saapParams[3].atomType2 = 0;
	saapParams[3].Aa0 = sa0co;
	saapParams[3].Aa1 = sa1co;
	saapParams[3].Aa2 = sa2co;
	saapParams[3].Aa3 = sa3co;
	saapParams[3].Aa4 = sa4co;
	saapParams[3].Aa5 = sa5co;
	saapParams[3].Aa6 = sa6co;

	saapParams[4].atomType1 = 1;
	saapParams[4].atomType2 = 1;
	saapParams[4].Aa0 = sa0cc;
	saapParams[4].Aa1 = sa1cc;
	saapParams[4].Aa2 = sa2cc;
	saapParams[4].Aa3 = sa3cc;
	saapParams[4].Aa4 = sa4cc;
	saapParams[4].Aa5 = sa5cc;
	saapParams[4].Aa6 = sa6cc;

	saapParams[5].atomType1 = 1;
	saapParams[5].atomType2 = 2;
	saapParams[5].Aa0 = sa0bc;
	saapParams[5].Aa1 = sa1bc;
	saapParams[5].Aa2 = sa2bc;
	saapParams[5].Aa3 = sa3bc;
	saapParams[5].Aa4 = sa4bc;
	saapParams[5].Aa5 = sa5bc;
	saapParams[5].Aa6 = sa6bc;

	saapParams[6].atomType1 = 2;
	saapParams[6].atomType2 = 0;
	saapParams[6].Aa0 = sa0bo;
	saapParams[6].Aa1 = sa1bo;
	saapParams[6].Aa2 = sa2bo;
	saapParams[6].Aa3 = sa3bo;
	saapParams[6].Aa4 = sa4bo;
	saapParams[6].Aa5 = sa5bo;
	saapParams[6].Aa6 = sa6bo;

	saapParams[7].atomType1 = 2;
	saapParams[7].atomType2 = 1;
	saapParams[7].Aa0 = sa0bc;
	saapParams[7].Aa1 = sa1bc;
	saapParams[7].Aa2 = sa2bc;
	saapParams[7].Aa3 = sa3bc;
	saapParams[7].Aa4 = sa4bc;
	saapParams[7].Aa5 = sa5bc;
	saapParams[7].Aa6 = sa6bc;

	saapParams[8].atomType1 = 2;
	saapParams[8].atomType2 = 2;
	saapParams[8].Aa0 = sa0bb;
	saapParams[8].Aa1 = sa1bb;
	saapParams[8].Aa2 = sa2bb;
	saapParams[8].Aa3 = sa3bb;
	saapParams[8].Aa4 = sa4bb;
	saapParams[8].Aa5 = sa5bb;
	saapParams[8].Aa6 = sa6bb;

}

void SAAPForce::write_force_info(ofstream& of)
{
    of << "SAAP Force" << endl;
    of << "cutoff2: " << cutOff2 << endl;
    if (switchOn)
        of << "compute SAAP force with switch function" << endl;
    if (computeCoulomb)
        of << "Non-bonded force include coulomb" << endl;
    of << "SAAP parameters:" << endl;
    of << "type1" << '\t' << "type2" << '\t' << "    a0" << '\t' << "   a1" << '\t' << "  a2" << '\t' << "      a3" <<  '\t' << "      a4" << '\t' << "      a5" <<  '\t' << "      a6" << endl;
    for (Int i = 0; i <myEnsemble->nAtomTypes; i++)
    {
        for (Int j = 0; j < myEnsemble->nAtomTypes; j++)
        {
	    int saapid = i*numAtomTypes + j;

            of << i << '\t' << j << '\t' << saapParams[saapid].Aa0 << '\t' << saapParams[saapid].Aa1 << '\t' << saapParams[saapid].Aa2 <<  '\t' << saapParams[saapid].Aa3 <<  '\t' << saapParams[saapid].Aa4;
	    of  <<  '\t' << saapParams[saapid].Aa5 <<  '\t' << saapParams[saapid].Aa6 << endl;
        }
    }
    of << "Long range correct energy: " << eLrc << endl;        
    of << endl;
}

void SAAPForce::write_energy(ofstream& of)
{
    of << "SAAP Energy: " << energy << endl;
}

// JC string LJForce::get_force_id() 
// JC {
// JC  return LJForce::forceIdentifier;
// JC }
