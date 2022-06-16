/** BuckinghamForce.cpp -- 
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

#include "BuckinghamForce.h" 


BuckinghamForce::BuckinghamForce(Ensemble* ensemble) : Force(ensemble)
{
    #ifdef DEBUG
        DEBUGMSG("Creating BuckinghamForce ");
    #endif

    buckParams = NULL;
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

    virt1 = new Vector3* [numMols];
    for (int nm = 0; nm < numMols; nm++)
        {
            int ns = myMols[nm].numAtoms;
            if (ns <= 1)
                continue;
            virt1[nm] = new Vector3 [ns];
	    if((virt1[nm]==NULL))
                ERRORMSG("memory allocation error for Buckingham Virtual positions");
        }

    buckParams = new BuckinghamParam[numAtomTypes*numAtomTypes];
          if (buckParams  == NULL)
              ERRORMSG( "memory allocation of Buckingham Parameters error");

    set_buck_param();

    if (myEnsemble->myConfig->use_long_correct())
        long_range_correct();           // calculate once at the beginning

    write_force_info(*sysdataFile);
}

// this implementation of compute() might be not efficient in parallelisation
// in case of OpenMP, data dependency may occur in the summation of force
// in case of MPI, many communication may be involved 
void BuckinghamForce::compute()
{

    Atom *atomi, *atomj;
    Atom* atom;
    Molecule *moli, *molj;
    Vector3 rij, fij, ri, rj;
    Double r, r_1, rij2, rij_2, sr2, sr6, sr12, tmpE, tmpf, rij6;
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

    VirtualPositions();		// Oxygen atoms of CO2  moved slightly to the center

    // go through each atom
    for (Int i = 0; i < numAtoms; i++)
    {
        atomi = &atoms[i]; 
        Int id = atomi->molID;
        rim = myEnsemble->molecules[id].massCenter; //
        Int size = atomi->get_list_size();          // size of pairlist
        qi = atomi->scaledCharge;
  	ri = atomi->position;
	if (!strcmp(myMols[id].molName, "CO2") && atomi->atomType == 0){ 
		ri = virt1[id][0];
	//cout << " ----- II ----- " << endl;
	}
	if (!strcmp(myMols[id].molName, "CO2") && atomi->atomType == 1){ 
		ri = virt1[id][1];
	}

        // go through the pairlist of atomi
        for (Int j = 0; j < size; j++)   /// Int j = i+1; j<numAtoms; j++
        {
            atomj = atomi->myPairList[j];
//	if (atoms[i]->molID != atoms[j]->molID ){
	if (atomi->molID != atomj->molID){
            id = atomj->molID;    //
	    rj = atomj->position;
	    if (!strcmp(myMols[id].molName, "CO2") && atomj->atomType == 0){ 
		rj = virt1[id][0];
	    }
	    if (!strcmp(myMols[id].molName, "CO2") && atomj->atomType == 1){ 
		rj = virt1[id][1];
	    }
            rij = ri - rj;
            rijm = myEnsemble->molecules[id].massCenter - rim;  //
            rijm -= rij;   //

            myEnsemble->apply_pbc(rij);
            rijm += rij;
            r = rij.length();
            rij2 = r*r;

            if (rij2 <= cutOff2)
            {
		int buckid = atomi->atomType*myMols[id].numAtoms + atomj->atomType;
		double Aa = buckParams[buckid].AA;
		double Bb = buckParams[buckid].BB;
		double Cc = buckParams[buckid].CC;

                r_1 = 1.0/r;
                rij_2 = r_1*r_1;
		rij6 = rij_2*rij_2*rij_2;
		tmpE = Bb*exp(-Cc*r) - Aa*rij6;
		tmpf = -Bb*Cc*r_1*exp(-Cc*r) + 6.0*Aa*rij6*rij_2;

//cout << atomi->atomType << "    " << atomj->atomType << "    " << Bb << endl;

		if (myEnsemble->myConfig->compute_lustig())
		{
		dulj  = -Bb*Cc*r*exp(-Cc*r) + 6.0*Aa*rij6;                  		// Lustig 2011 Eq. 14b
		d2ulj = (2*r + Cc*rij2)*Bb*Cc*exp(-Cc*r) - 54.0*Aa*rij6;	        // Lustig 2011 Eq. 14c
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

		moli = &myMols[atomi->molID];
		molj = &myMols[atomj->molID];
        	Atom *atomCi = moli->myAtoms[2];
        	Atom *atomCj = molj->myAtoms[2];
 
		if ( (atomi->atomType == 0) || (atomi->atomType == 1) ){ 
                	atomi->force += 0.91*fij;
			atomCi->force += 0.09*fij;
		}
		else if (atomi->atomType == 2){
			atomi->force += fij;
		}
		
	        if ( (atomj->atomType == 0) || (atomj->atomType == 1) ){ 
                	atomj->force -= 0.91*fij;
			atomCj->force -= 0.09*fij;
		}
	        else if (atomj->atomType == 2){ 
                	atomj->force -= fij;
		}

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
    	dudv   = dudv/(3*voll);
    	d2udv  = d2udv/(9*vol2);

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
void BuckinghamForce::long_range_correct()
{
    Int numTypes = myEnsemble->nAtomTypes;
    Int *atomNum = new Int[numTypes];
    Double val = 0.0;
    double valU = 0.0, valdU = 0.0, vald2U = 0.0;
    double ljcutU = 0.0, ljcutdU = 0.0, ljcutd2U = 0.0;
    double sigma6, cutoff3, volrc, NNval;

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

    for (Int i = 0; i < numTypes; i++)
    {
        for (Int j = 0; j < numTypes; j++)
        {
	    //if (j < i) {  continue; }

	    int buckid = i*numTypes + j;
	    double Aa = buckParams[buckid].AA;
	    cout << " i = " << i << "    j = " << j << "    Aij = "<< Aa << endl;

            val += atomNum[i]*atomNum[j]*Aa;    //Gromacs value

	    /*if (myEnsemble->myConfig->compute_lustig())
	    {
            	sigma6   = sigma2*sigma2*sigma2;                                        //For Lustig 2011 
            	ljcutU   = eps4*(sigma6*sigma6/(cutoff3*cutoff3*cutoff3) - 3*sigma6/cutoff3);
            	ljcutdU  = eps4*(4*sigma6*sigma6/(cutoff3*cutoff3*cutoff3) - 6*sigma6/cutoff3);
            	ljcutd2U = eps4*(20*sigma6*sigma6/(cutoff3*cutoff3*cutoff3) - 18*sigma6/cutoff3);

	    	NNval  = (atomNum[i] + atomNum[j])*(atomNum[i] + atomNum[j])*0.250 - (atomNum[i] + atomNum[j])*0.5;

		cout << "  atomNum[i]  =  " << atomNum[i] <<  "  atomNum[j]  =  "<< atomNum[j] << " NNval  =  " << NNval << endl; 

            	valU   += NNval*ljcutU;
            	valdU  += NNval*ljcutdU;
            	vald2U += NNval*ljcutd2U;
            }*/
        }
    }

    eLrc   = -2*PI*val/(3*volrc*cutoff3);     //Gromacs LRC

 /*   if (myEnsemble->myConfig->compute_lustig())
    {
	eLrc   =  2*PI*valU/(9*volrc);                          // Lustig 2011 LRC
    	dUlrc  = -2*PI*valdU/(9*volrc*volrc);                   // Lustig 2011 LRC
    	d2Ulrc = -2*PI*vald2U/(9*volrc*volrc*volrc);            // Lustig 2011 LRC
    }*/


    vLrc = -eLrc;
    vgromLRC = -4*PI*val/(3*volrc*volrc*cutoff3);	// Pressure, not virial
          
    cout << "  eLrc  =  " << eLrc <<  " vLrc  =  "<< vLrc << " vgromLRC  =  " << vgromLRC << endl;
}

//----------------------------------------------------------------------------------------------
// Virtual positions of Buckingham based on
// Tsuzuki et al. / Chemical Physics Letters 255 (1996) 347-349
//
// Permutation of distances based from the Lectures of 
// Berendsen, H. J. C., van Gunsteren, W. F.  Molecular dynamics simulations: Techniques and
// approaches. In: Molecular Liquids-Dynamics and Interactions. et al., A. J. B. ed. NATO
// ASI C 135. Reidel Dordrecht, The Netherlands 1984 475–500. Figure 8.
// Also available in Gromacs Manual.
// ----------------------------------------------------------------------------------------------

void BuckinghamForce :: VirtualPositions()
{
   Atom* atom;
   Molecule *mol;
   co2bond = myEnsemble->myParams->bondParams[0].r0;

   const double ll = co2bond*0.5;	// Only if the bond is defined as the distance between Oxygens
   const double cc = 0.91*ll;
   const double bb = ll - cc;
   for (int i = 0; i<numMols; i++)
   {
     if (!strcmp(myMols[i].molName, "CO2")) {

        mol = &myMols[i];
        Atom *atomO0 = mol->myAtoms[0];
        Atom *atomO1 = mol->myAtoms[1];
        Atom *atomC2 = mol->myAtoms[2];

     virt1[i][0]  = (bb/ll)*atomC2->position + (cc/ll)*atomO0->position;
     virt1[i][1]  = (bb/ll)*atomC2->position + (cc/ll)*atomO1->position;

             /*Vector3 co1dif =  virt1[i][1] - virt1[i][0];   // New oxygen separation 
             double diffco1 = co1dif.length();
             Vector3 co2dif =  virt1[i][0] - atomO0->position;      // how was moved 
             double diffco2 = co2dif.length();
             cout<<"  virtualSep   =  "<< diffco1 << " oxysep  =  "<< diffco2 << "   bb = " << bb << endl;
		*/
	 }
   }

}

//-------------------------------------------------------------------------------------------
// Pair table for Buckingham Parameters obtained as 
// atomType1 * numtypes + atomType2
// Constants from Tsuzuki et al. J. Phys. Chem., Vol. 100, No. 11, 1996 
// model 3 of Table 2 in the article
// -----------------------------------------------------------------------------------------

void BuckinghamForce::set_buck_param()                              
{
	buckParams[0].atomType1 = 0;
	buckParams[0].atomType2 = 0;
	buckParams[0].AA = AAoo;
	buckParams[0].BB = BBoo;
	buckParams[0].CC = CCoo;

	buckParams[1].atomType1 = 0;
	buckParams[1].atomType2 = 1;
	buckParams[1].AA = AAoo;
	buckParams[1].BB = BBoo;
	buckParams[1].CC = CCoo;

	buckParams[2].atomType1 = 0;
	buckParams[2].atomType2 = 2;
	buckParams[2].AA = AAco;
	buckParams[2].BB = BBco;
	buckParams[2].CC = CCco;

	buckParams[3].atomType1 = 1;
	buckParams[3].atomType2 = 0;
	buckParams[3].AA = AAoo;
	buckParams[3].BB = BBoo;
	buckParams[3].CC = CCoo;

	buckParams[4].atomType1 = 1;
	buckParams[4].atomType2 = 1;
	buckParams[4].AA = AAoo;
	buckParams[4].BB = BBoo;
	buckParams[4].CC = CCoo;

	buckParams[5].atomType1 = 1;
	buckParams[5].atomType2 = 2;
	buckParams[5].AA = AAco;
	buckParams[5].BB = BBco;
	buckParams[5].CC = CCco;

	buckParams[6].atomType1 = 2;
	buckParams[6].atomType2 = 0;
	buckParams[6].AA = AAco;
	buckParams[6].BB = BBco;
	buckParams[6].CC = CCco;

	buckParams[7].atomType1 = 2;
	buckParams[7].atomType2 = 1;
	buckParams[7].AA = AAco;
	buckParams[7].BB = BBco;
	buckParams[7].CC = CCco;

	buckParams[8].atomType1 = 2;
	buckParams[8].atomType2 = 2;
	buckParams[8].AA = AAcc;
	buckParams[8].BB = BBcc;
	buckParams[8].CC = CCcc;
}

void BuckinghamForce::write_force_info(ofstream& of)
{
    of << "Buckingham Force" << endl;
    of << "cutoff2: " << cutOff2 << endl;
    if (switchOn)
        of << "compute Buckingham force with switch function" << endl;
    if (computeCoulomb)
        of << "Non-bonded force include coulomb" << endl;
    of << "Buckingham parameters:" << endl;
    of << "type1" << '\t' << "type2" << '\t' << "A" << '\t' << "B" << '\t' << "C" << endl;
    for (Int i = 0; i <myEnsemble->nAtomTypes; i++)
    {
        for (Int j = 0; j < myEnsemble->nAtomTypes; j++)
        {
	    int buckid = i*numAtomTypes + j;

            of << i << '\t' << j << '\t' << buckParams[buckid].AA << '\t' << buckParams[buckid].BB <<  '\t' << buckParams[buckid].CC << endl;
        }
    }
    of << "Long range correct energy: " << eLrc << endl;        
    of << endl;
}

void BuckinghamForce::write_energy(ofstream& of)
{
    of << "Buckingham Energy: " << energy << endl;
}

// JC string LJForce::get_force_id() 
// JC {
// JC  return LJForce::forceIdentifier;
// JC }
