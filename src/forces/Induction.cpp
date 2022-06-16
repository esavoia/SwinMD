/** Induction.cpp -- 
 **
 ** Copyright (C) 2003
 ** Centre for Molecular Simulation (CMS)
 ** School of Information Technology
 ** Swinburne University of Technology
 ** PO Box 21 Hawthorn, Vic 3122, Australia
 **
 ** Ewald Method for Induced Multipole systems by
 ** Author: eoyarzua
 ** Email: eoyarzua@swin.edu.au
 **/

#include "Induction.h" 

///#include "erf.h"
#include <mpi.h>

Induction::Induction(Ensemble* ensemble) : Force(ensemble) 
{
    DEBUGMSG("Creating Induction");

//JC    forceIdentifier = "EwaldForce"; // in order to comply with the cluster compiler
    numAtoms = myEnsemble->nAtoms;
    numMols = myEnsemble->nMols;
    myMols = myEnsemble->molecules;
    volume = myEnsemble->volume;
    volume2 = volume*volume;
    volumer = 1.0/volume;
    nAtomsx3 = 3*numAtoms;

    //------------------------------------------------------------------------------------------
    // Get parameters and booleans from SimConfiguration.h
    //------------------------------------------------------------------------------------------

    alpha             = myEnsemble->myConfig->get_alpha();
    alpha2            = alpha*alpha;
    alpha2PI          = 2.0*alpha/SQRT_PI;
    alphaR4           = 1.0/(4.0*alpha2);                 
    realCut           = myEnsemble->myConfig->get_cutoffEw();
    realCut2          = realCut*realCut;
    cutoffInd         = myEnsemble->myConfig->get_cutoffInd();
    cutoffInd2        = cutoffInd*cutoffInd;
    a                 = myEnsemble->myConfig->get_adamping();
    kCut              = myEnsemble->myConfig->get_k_cutoff(); 
    correctSurfDipole = myEnsemble->myConfig->compute_surf_correct(); 
    useConstraint     = myEnsemble->myConfig->is_constraint_on();
    doMinimization    = myEnsemble->myConfig->do_induction_minimization();
    dampTij           = myEnsemble->myConfig->use_damping_induction();
    statEnsem         = myEnsemble->myConfig->get_ensemble_status();  
    computeTorques    = myEnsemble->myConfig->compute_torques(); 
    dampAMOEBA        = myEnsemble->myConfig->use_damping_indAmoeba(); 

   if (doMinimization) { DEBUGMSG(" Performs minimization of induced dipoles by means of CG method");}
	else { DEBUGMSG(" Doesn't Perform minimization of induced dipoles: iAMOEBA"); }

    if (computeTorques) { DEBUGMSG(" Induction: Compute the torques due to dipoles and translate them into atomic forces");}
        else { DEBUGMSG(" Torques are NOT included in the Induction computation"); }

    if (dampTij) { DEBUGMSG(" Induction: Tij dipole tensor is being damped.");}
        else { DEBUGMSG(" Induction: Tij dipole tensor is not damped, check the divergency of CG method"); }

    if (dampAMOEBA) { DEBUGMSG(" Induction: Real-space damped according to Tinker/AMOEBA code.");}
        else { DEBUGMSG(" Induction: Real-space damped according to Bl*lamb"); }

    stepr = 0.0;
    chargedEnergy = 0.0;
    longEnergy = 0.0;
    longderiv  = 0.0;
    long2ndDeriv  = 0.0;
    tableSize = TAB_SIZE + 2;

    myL2G  = NULL;

    trq  = NULL;
    expK = NULL;
    qQsin_i = NULL;
    qQcos_i = NULL;
    mksin_i = NULL;
    mkcos_i = NULL;
    pksin_i = NULL;
    pkcos_i = NULL;
    kXmsina = NULL;
    kXmcosa = NULL;
    KXQcosa = NULL;
    KXQsina = NULL;

    A   = NULL;
    Pi  = NULL;
    Ei0 = NULL;

    forceTable = NULL;
    potentialTable = NULL;
    doUpdate = false;
    if (statEnsem == 1){
        doUpdate = true;
        DEBUGMSG(" Volume Update for NPT ensemble on Induction Energy");
    }

    myL2G  = new Local2Global(myEnsemble);

    qQsin_i = new double[numAtoms];
    qQcos_i = new double[numAtoms];
    mksin_i = new double[numAtoms];
    mkcos_i = new double[numAtoms];
    pksin_i = new double[numAtoms];
    pkcos_i = new double[numAtoms];
    trq     = new Vector3 [numAtoms];
    kXmsina = new Vector3 [numAtoms];
    kXmcosa = new Vector3 [numAtoms];
    KXQcosa = new Vector3 [numAtoms];
    KXQsina = new Vector3 [numAtoms];

    Pi  = new double[3*numAtoms];
    Ei0 = new double[3*numAtoms];
    A   = new double* [3*numAtoms];

    Dummy = new Double[numMols*3];
    if ((qQsin_i==NULL)||(qQcos_i==NULL)||(mksin_i==NULL)||(mkcos_i==NULL) ||(pksin_i==NULL)||(pkcos_i==NULL))
        ERRORMSG("memory allocation error for Induction with Ewald");

    if ((trq==NULL)||(kXmsina==NULL)||(kXmcosa==NULL)||(KXQcosa==NULL)||(KXQsina==NULL))
        ERRORMSG("memory allocation error for Induction with Ewald (2)");

    if ((Pi==NULL)||(Ei0==NULL)||( A==NULL))
        ERRORMSG("memory allocation error for Induction with Ewald (3)");

    for(int i = 0; i < 3*numAtoms; i++)
        {
	    A[i] = new double [3*numAtoms]; 
            if((A[i]==NULL))
                ERRORMSG("memory allocation error for Induction with Ewald");
        }

    set_up();
    DEBUGMSG("Induction set up");

    // pre-compute self energy correction - Gaussian term
    Double qSum = 0.0;
    Double q2Sum = 0.0;
    for(Int i = 0; i < numAtoms; i++)
    {
        qSum += atoms[i].scaledCharge;
        q2Sum += atoms[i].scaledCharge*atoms[i].scaledCharge;
    }
    if(fabs(qSum) > 1.0e-5)
       chargedEnergy = -PI*qSum*qSum/(2.0*volume*alpha2);  // Volume needs to be updated in NpT cases with highly charged systems (Ions probably?)
    
    write_force_info(*sysdataFile);
 
   // prank = MPI::COMM_WORLD.Get_rank();  // Jc: get the MPI parameters // Updated to MPI-3 by elton
    //psize = MPI::COMM_WORLD.Get_size();

    MPI_Comm_rank(MPI_COMM_WORLD, &prank);
    MPI_Comm_size(MPI_COMM_WORLD, &psize);

}


Induction::~Induction()
{
    if(A   != NULL) delete A; 
    if(Pi  != NULL) delete Pi; 
    if(Ei0 != NULL) delete Ei0; 
    if(trq != NULL) delete [] trq;
    if(kXmsina != NULL) delete [] kXmsina;
    if(kXmcosa != NULL) delete [] kXmcosa;
    if(KXQcosa != NULL) delete [] KXQcosa;
    if(KXQsina != NULL) delete [] KXQsina;
    if(expK != NULL) delete [] expK;
    if(qQsin_i != NULL) delete [] qQsin_i;
    if(qQcos_i != NULL) delete [] qQcos_i;
    if(mksin_i != NULL) delete [] mksin_i;
    if(mkcos_i != NULL) delete [] mkcos_i;
    if(pksin_i != NULL) delete [] pksin_i;
    if(pkcos_i != NULL) delete [] pkcos_i;
    if(forceTable != NULL) delete [] forceTable;
    if(potentialTable != NULL) delete [] potentialTable;    
}

void Induction::InducedDipoles()
{
    int i3, molIDi, j3, molIDj, idj;
    double qi, dipxi, dipyi, dipzi;
    double qj, dipxj, dipyj, dipzj;
    double lamb3, lamb5, lamb7, u3, exp_au3;
    double r, rij2, rij_1, rij_3, rij_5, rij_7;
    double Qxxi, Qxyi, Qxzi, Qyyi, Qyzi, Qzzi; 
    double Qyxi, Qzxi, Qzyi; 
    double Qxxj, Qxyj, Qxzj, Qyyj, Qyzj, Qzzj; 
    double Qyxj, Qzxj, Qzyj; 
    double polarity, polarj;
    double polmin = 0.00000001;
    double QiI, QjI, QjRij, QiRji;
    Vector3 ri, rj, rij;
    Vector3 dipi, dipj;
    Vector3 bi, bj;
    
    //------------------------------------------------------------------------------
    // First, clean this up
    //------------------------------------------------------------------------------

    for(int i = 0; i < 3*numAtoms;i++)
    {
        Ei0[i] = 0.0;
        Pi[i]  = 0.0;
        for(int j = 0; j < 3*numAtoms;j++)
        {
            A[i][j] = 0.0;
        } 
    }

    //------------------------------------------------------------------------------
    // Now start with the Natoms loop. The 3*nAtoms is the solution
    // to implement the A matrix as a linear coupled system (3 components at each site).
    //------------------------------------------------------------------------------

    for (int i = 0; i < numAtoms; i++)
        {
	    i3 = 3 * i;
	    molIDi = atoms[i].molID;
            ri = atoms[i].position;
            qi = atoms[i].scaledCharge/SQRTCOULOMBCONSTANT;
	    polarity = MAX(polmin,atoms[i].polar);

            dipi = myL2G->dip[i];

	    //-----------------------------------------------------------------------
	    // Call Traceless Quadrupoles from Local2Global
	    //-----------------------------------------------------------------------

            Qxxi  = myL2G->Qxx[i];
            Qxyi  = myL2G->Qxy[i];
            Qxzi  = myL2G->Qxz[i];
            Qyyi  = myL2G->Qyy[i];
            Qyzi  = myL2G->Qyz[i];
            Qzzi  = myL2G->Qzz[i];

            Qyxi  = Qxyi;
            Qzxi  = Qxzi;
            Qzyi  = Qyzi;

            //----------------------------------------------------------------------------
            //  Before jumping to the j-index, Let's build the
            //  diagonal elements of the A matrix.
	    //----------------------------------------------------------------------------
	    
	    A[i3][i3]      = 1.0/polarity; 
	    A[i3][i3+1]    = 0.0; 
	    A[i3][i3+2]    = 0.0; 
	    A[i3+1][i3]    = 0.0; 
	    A[i3+1][i3+1]  = 1.0/polarity; 
	    A[i3+1][i3+2]  = 0.0; 
	    A[i3+2][i3]    = 0.0; 
	    A[i3+2][i3+1]  = 0.0; 
	    A[i3+2][i3+2]  = 1.0/polarity; 

            for(int j = i+1; j < numAtoms; j++)
        	{
		    j3     = 3 * j;
	    	    molIDj = atoms[j].molID;
		    idj    = atoms[j].atomID;

	    	    if (j != idj) { 
			ERRORMSG("Induction Minimization: ID counting is wrong. atomIDj =/= j or Check your water molecule IDs (3 sites per molecule). ");}

            	    rj = atoms[j].position;
                    qj = atoms[j].scaledCharge/SQRTCOULOMBCONSTANT;

            	    dipj = myL2G->dip[j];

	    	    //-----------------------------------------------------------------------
	    	    // Call Traceless Quadrupoles from Local2Global
	    	    //-----------------------------------------------------------------------

            	    Qxxj  = myL2G->Qxx[j];
            	    Qxyj  = myL2G->Qxy[j];
            	    Qxzj  = myL2G->Qxz[j];
            	    Qyyj  = myL2G->Qyy[j];
            	    Qyzj  = myL2G->Qyz[j];
            	    Qzzj  = myL2G->Qzz[j];

            	    Qyxj  = Qxyj;
            	    Qzxj  = Qxzj;
            	    Qzyj  = Qyzj;

		    rij = ri - rj;
		    myEnsemble->apply_pbc(rij);
		    r = rij.length();
            	    rij2 = r*r;

		    
		    if (rij2 <= cutoffInd2)
	                {
			    //----------------------------------------------------------------------------
			    //  inverse of the length of rij
			    //----------------------------------------------------------------------------
			    
			    rij_1 = 1.0/r;
			    rij_3 = rij_1*rij_1*rij_1;
			    rij_5 = rij_3*rij_1*rij_1;
			    rij_7 = rij_5*rij_1*rij_1;

			    //----------------------------------------------------------------------------
			    //  Dipole-Dipole Tensor
			    //----------------------------------------------------------------------------

	    		    polarj  = MAX(polmin,atoms[idj].polar);
			    u3      = rij2*r/(sqrt(polarity * polarj));
			    exp_au3 = exp(-a*u3);
			    lamb3   = 1.0 - exp_au3;
			    lamb5   = 1.0 - (1.0 + a*u3)*exp_au3;
			    lamb7   = 1.0 - (1.0 + a*u3 + 3.0*a*a*u3*u3/5.0)*exp_au3;

			    if (dampTij) 
			    {
		                Tij[0][0] = lamb5*3.0*rij_5*rij.x*rij.x - lamb3*rij_3;  
		                Tij[0][1] = lamb5*3.0*rij_5*rij.x*rij.y;        
		                Tij[0][2] = lamb5*3.0*rij_5*rij.x*rij.z;        
		                Tij[1][0] = Tij[0][1];
		                Tij[1][1] = lamb5*3.0*rij_5*rij.y*rij.y - lamb3*rij_3;
		                Tij[1][2] = lamb5*3.0*rij_5*rij.y*rij.z;
		                Tij[2][0] = Tij[0][2];
	   	                Tij[2][1] = Tij[1][2];
		                Tij[2][2] = lamb5*3.0*rij_5*rij.z*rij.z - lamb3*rij_3;

			    } else {

		                Tij[0][0] = 3.0*rij_5*rij.x*rij.x - rij_3;  
		                Tij[0][1] = 3.0*rij_5*rij.x*rij.y;        
		                Tij[0][2] = 3.0*rij_5*rij.x*rij.z;        
		                Tij[1][0] = Tij[0][1];
		                Tij[1][1] = 3.0*rij_5*rij.y*rij.y - rij_3;
		                Tij[1][2] = 3.0*rij_5*rij.y*rij.z;
		                Tij[2][0] = Tij[0][2];
	   	                Tij[2][1] = Tij[1][2];
		                Tij[2][2] = 3.0*rij_5*rij.z*rij.z - rij_3;

			    }

			    //----------------------------------------------------------------------------
			    //  Now let's build the A matrix.
			    //----------------------------------------------------------------------------

			    A[i3][j3]      = -Tij[0][0];
			    A[i3][j3+1]    = -Tij[0][1];
			    A[i3][j3+2]    = -Tij[0][2];
			    A[i3+1][j3]    = -Tij[1][0];
			    A[i3+1][j3+1]  = -Tij[1][1];
			    A[i3+1][j3+2]  = -Tij[1][2];
			    A[i3+2][j3]    = -Tij[2][0];
			    A[i3+2][j3+1]  = -Tij[2][1];
			    A[i3+2][j3+2]  = -Tij[2][2];
			    
			    A[j3][i3]      = -Tij[0][0];
			    A[j3][i3+1]    = -Tij[0][1];
			    A[j3][i3+2]    = -Tij[0][2];
			    A[j3+1][i3]    = -Tij[1][0];
			    A[j3+1][i3+1]  = -Tij[1][1];
			    A[j3+1][i3+2]  = -Tij[1][2];
			    A[j3+2][i3]    = -Tij[2][0];
			    A[j3+2][i3+1]  = -Tij[2][1];
			    A[j3+2][i3+2]  = -Tij[2][2];


		             if (!myEnsemble->exclusion_check(i, idj))
                		{
                    		    if (molIDi == molIDj) {
                        		ERRORMSG("Induction computation: atoms in the same molecule are interacting. Check exclusion_check() pairs.");}

				    //----------------------------------------------------------------------------
				    //  Vectors needed for bi/Ei0 terms
				    //----------------------------------------------------------------------------
		
				    Vector3 rji, QiXrji, QjXrij, RijXmuj, RjiXmui;

				    rji = -rij;

				    RijXmuj.x = rij.x*rij.x*dipj.x + rij.x*rij.y*dipj.y + rij.x*rij.z*dipj.z; 
				    RijXmuj.y = rij.y*rij.x*dipj.x + rij.y*rij.y*dipj.y + rij.y*rij.z*dipj.z; 
				    RijXmuj.z = rij.z*rij.x*dipj.x + rij.z*rij.y*dipj.y + rij.z*rij.z*dipj.z; 

				    RjiXmui.x = rji.x*rji.x*dipi.x + rji.x*rji.y*dipi.y + rji.x*rji.z*dipi.z;
				    RjiXmui.y = rji.y*rji.x*dipi.x + rji.y*rji.y*dipi.y + rji.y*rji.z*dipi.z;
				    RjiXmui.z = rji.z*rji.x*dipi.x + rji.z*rji.y*dipi.y + rji.z*rji.z*dipi.z;

				    QiXrji.x = Qxxi*rji.x + Qxyi*rji.y + Qxzi*rji.z; 
				    QiXrji.y = Qyxi*rji.x + Qyyi*rji.y + Qyzi*rji.z; 
				    QiXrji.z = Qzxi*rji.x + Qzyi*rji.y + Qzzi*rji.z; 

				    QjXrij.x = Qxxj*rij.x + Qxyj*rij.y + Qxzj*rij.z;
				    QjXrij.y = Qyxj*rij.x + Qyyj*rij.y + Qyzj*rij.z;
				    QjXrij.z = Qzxj*rij.x + Qzyj*rij.y + Qzzj*rij.z;

				    //----------------------------------------------------------------------------
				    //  Constants needed for bi/Ei0 terms
				    //----------------------------------------------------------------------------

				    QiI   = Qxxi + Qyyi + Qzzi;
				    QjI   = Qxxj + Qyyj + Qzzj;
				    QjRij = Qxxj*rij.x*rij.x + Qxyj*rij.x*rij.y + Qxzj*rij.x*rij.z \
				          + Qyxj*rij.y*rij.x + Qyyj*rij.y*rij.y + Qyzj*rij.y*rij.z \
				          + Qzxj*rij.z*rij.x + Qzyj*rij.z*rij.y + Qzzj*rij.z*rij.z ;

				    QiRji = Qxxi*rij.x*rij.x + Qxyi*rij.x*rij.y + Qxzi*rij.x*rij.z \
					  + Qyxi*rij.y*rij.x + Qyyi*rij.y*rij.y + Qyzi*rij.y*rij.z \ 
					  + Qzxi*rij.z*rij.x + Qzyi*rij.z*rij.y + Qzzi*rij.z*rij.z ; 

				    //----------------------------------------------------------------------------
				    //  bi/Ei0 computation
				    //----------------------------------------------------------------------------
				    
        			    bi = qj*lamb3*rij_3*rij + 3.0*lamb5*rij_5*RijXmuj - lamb3*rij_3*dipj + 15.0*lamb7*rij_7*QjRij*rij - 6.0*lamb5*rij_5*QjXrij - 3.0*lamb5*rij_5*QjI*rij;

        			    Ei0[i3]   += bi.x;
        			    Ei0[i3+1] += bi.y;
        			    Ei0[i3+2] += bi.z;

				    bj = qi*lamb3*rij_3*rji + 3.0*lamb5*rij_5*RjiXmui - lamb3*rij_3*dipi + 15.0*lamb7*rij_7*QiRji*rji - 6.0*lamb5*rij_5*QiXrji - 3.0*lamb5*rij_5*QiI*rji;

        			    Ei0[j3]   += bj.x;
        			    Ei0[j3+1] += bj.y;
        			    Ei0[j3+2] += bj.z;

			    	} // end molecular pairs check
		        } //end cutoffInd loop

		} // end j loop
	} // end i loop

    //------------------------------------------------------------------------------
    // Finally, solve the Ax=b system by means of the 
    // Conjugate Gradient method for a symmetric matrix or not.
    // The solution without minimizing the Induced Dipoles is based on the
    // iAMOEBA model (J. Phys. Chem. B 2013, 117, 9956−9972). 
    //------------------------------------------------------------------------------
    
    if (doMinimization) {
	if (psize > 1) {
            ParallelCG();
	}
	else {
            ConjugateGradient();
	}
    }
    else {
    	for(int i = 0; i < numAtoms;i++)
    	{
	    int i3 = 3 * i;
       	    Pi[i3]   = atoms[i].polar * Ei0[i3];
	    Pi[i3+1] = atoms[i].polar * Ei0[i3+1];
	    Pi[i3+2] = atoms[i].polar * Ei0[i3+2];
	}
    }
}

void Induction::DummyPosition()

{
   Vector3 Roh1, Roh2, Rhh, Rpo, Dummyi;
   Atom *atom1, *atom2, *atom3;
   Double  RRpo;

   for (int i = 0;i<numMols;i++)
   {
     if (!strcmp(myMols[i].molName, "MCY")) {
     int j = 3 * i;
     atom1 = &atoms[j];
     atom2 = &atoms[j + 1];
     atom3 = &atoms[j + 2];

     Rpo.x = 0.5*(atom3->position.x + atom1->position.x); 
     Rpo.y = 0.5*(atom3->position.y + atom1->position.y);
     Rpo.z = 0.5*(atom3->position.z + atom1->position.z);
     Rpo   = Rpo - atom2->position;
     RRpo  = Rpo.length(); 
     Dummyi = atom2->position + (R_OP/(RRpo*F_Length))*Rpo; // Jc: R_OP needs F_Length to nm  !
     Dummy[j]     = Dummyi.x;
     Dummy[j + 1] = Dummyi.y;
     Dummy[j + 2] = Dummyi.z;
     }
   }

}

void Induction::compute()
{
    energy = 0.0;
    selfEnergy = 0.0;
    selfdUdV   = 0.0;
    realEnergy = 0.0;
    realderiv = 0.0;
    real2ndDeriv = 0.0;
    molderiv = 0.0;
    mol2ndDeriv = 0.0;
    intraMolEnergy = 0.0;
    intraVircorrect = 0.0;
    surfDipoleEnergy = 0.0;
    surfderiv = 0.0;
    surf2ndderiv = 0.0;
    totTrqI = 0.0;

    DummyPosition();

    myL2G->FrameofReference();

    InducedDipoles();

    if (doMinimization) {
	double r[3*numAtoms];
	int div = 0;
	for(int i=0;i<3*numAtoms;i++)
	{
	    r[i] =  Ei0[i];
	    for(int j = 0;j<3*numAtoms;j++)
	    {
		r[i] -= A[i][j] * Pi[j];
	    }
	    if (r[i] > 1.0) {
		div += 1;
	    }
	}
	if (div > 10 && prank == 0) {
	    DEBUGMSG(" ACHTUNG/WARNING: CG is starting to diverge!");
	    cout << " ACHTUNG/WARNING: CG is starting to diverge: " << div << " out of: " << 3*numAtoms << endl;
	}
    }

    for (Int i = XX; i <= ZZ; i++)
        virial[i] = 0.0;

    if(doUpdate)
        set_up();

    self_correct(selfEnergy, selfdUdV);

    real_term(realEnergy, realderiv, real2ndDeriv);                     // Jc: the first term in Ewald

    intraMolCorrect(intraMolEnergy, intraVircorrect, molderiv, mol2ndDeriv);           // Molec correction in Ewald

    if(correctSurfDipole){
        surfDipoleCorrect(surfDipoleEnergy, surfderiv, surf2ndderiv);   // Jc: the first term in Ewald
    }

    myEnsemble->myEnergy.realPolEnergy = realEnergy;
    myEnsemble->myEnergy.correctPolEnergy = selfEnergy; //Jc: selfEnergy and charged Energy need no reduce
    myEnsemble->myEnergy.molCorrectPolEnergy = intraMolEnergy;
    myEnsemble->myEnergy.surfacePolEnergy = surfDipoleEnergy;

    myEnsemble->myLustig.intravirPolCorrect = intraVircorrect;
    myEnsemble->myLustig.realPolDeriv = realderiv;
    myEnsemble->myLustig.realSecPolDeriv = real2ndDeriv;

    myEnsemble->myLustig.molPolDeriv = molderiv;
    myEnsemble->myLustig.molSecPolDeriv = mol2ndDeriv;

    myEnsemble->myLustig.surfPolDeriv = surfderiv;
    myEnsemble->myLustig.surfSecPolDeriv = surf2ndderiv;

    myEnsemble->myLustig.selfPolDeriv = selfdUdV;

 
    for (Int k = XX; k <= ZZ; k++)
        myEnsemble->virial[k] -= virial[k];
}     
 
void Induction::compute_long()
{
    DummyPosition();    

    myL2G->FrameofReference();

    for (Int i = XX; i <= ZZ; i++)
        longVirial[i] = 0.0;

    for (Int i = 0; i < numAtoms; i++)
        atoms[i].longForce = 0.0;

    reciprocal_term(longEnergy, longderiv, long2ndDeriv);
    myEnsemble->myEnergy.longPolEnergy = longEnergy;
    
    myEnsemble->myLustig.longPolDeriv = longderiv;
    myEnsemble->myLustig.longSecPolDeriv = long2ndDeriv;


    myEnsemble->myEnergy.TotTorq += totTrqI;  // Assuming that trqperm is computed first.

    for(Int i = 0; i < ZZ; i++)
        myEnsemble->longVirial[i] += longVirial[i];
}

void Induction::set_up()
{
    Double kx, ky, kz;
    Double kCut2 = kCut*kCut; 
    Int nxMax, nyMax, nzMax;


    // build reciprocal vector - kmax=kcut=2*PI*nxMax/Lx
    myEnsemble->box_dimenssion(&boxLx, &boxLy, &boxLz);
    nxMax = (int) (kCut*boxLx/(2*PI));
    nyMax = (int) (kCut*boxLy/(2*PI));
    nzMax = (int) (kCut*boxLz/(2*PI));
    kVector.clear();
    kModulu.clear();

    Double vx = 2.0*PI/boxLx;
    Double vy = 2.0*PI/boxLy;
    Double vz = 2.0*PI/boxLz;

    // should we check if nx(ny,nz)Max <= 1
    for(int nx = 0; nx <= nxMax; nx++)
    {
        kx = nx*vx;
        for(int ny = (nx==0 ? 0:-nyMax); ny <= nyMax; ny++)
        {
            ky = ny*vy;
            for(int nz = ((nx==0 && ny==0)? 1:-nzMax); nz <= nzMax; nz++)
            {
                kz = nz*vz;
                Double k2 = kx*kx + ky*ky + kz*kz;
                if(k2 < kCut2)
                {
                    kVector.push_back(Vector3(kx, ky, kz));
                    kModulu.push_back(k2);
                }
            }
        }
    }   // done reciprocal vector

    // pre-compute expK[] for reciprocal term
    #ifdef DEBUG
    //    DEBUGMSG("init expK table");
    #endif
    init_table(kVector.size(), &expK);

    for(Int i = 0; i < kVector.size(); i++)
    {
        if (kModulu[i] == 0)
            ERRORMSG("Divided by zero");
        if (expK == NULL)
             ERRORMSG("expK - NULL pointer");
        expK[i] = exp(-kModulu[i]*alphaR4)/kModulu[i];        
    }
  
    // build force & potential look-up tables for real term
    #ifdef USE_LOOK_UP_TABLE                      
    build_look_up_table();
    #endif

}   // end set_up()

// build lookup table for real term
void Induction::build_look_up_table()
{
    Double r, r2;
    Double step = realCut2/(tableSize - 2);

    #ifdef DEBUG
        DEBUGMSG("build look_up table");
    #endif

    stepr = 1.0/step;
    
    init_table(tableSize, &forceTable);
    init_table(tableSize, &potentialTable);

    for(Int i = 1; i < tableSize; i++)
    {
        r2 = i*step;
        r = sqrt(r2);
        potentialTable[i] = erfc(alpha*r)/r;
        forceTable[i] = (potentialTable[i]+alpha2PI*exp(-alpha2*r2))/r2;
    }
    potentialTable[0] = 2.0*potentialTable[1];
    forceTable[0] = 2.0*forceTable[1];    
}

void Induction::init_table(Int size, Double** table)
{
    if (*table != NULL)
        delete [] *table;
    *table = new Double[size];
    if(*table == NULL) 
        ERRORMSG("Induction::init_table(), not enough memory");
}

void Induction::self_correct(double &selfEnergy, double &selfderiv)
{

    #ifdef DEBUG2
        DEBUGMSG("Computing Induction Self-Interaction Correction term");
    #endif

    double qi, dipxi, dipyi, dipzi;
    Vector3 dipi, indPi;
    int molIDi, i3;
    double selfterm, term;
    
    selfEnergy = 0.0;
    selfderiv  = 0.0;

    //--------------------------------------------------------------------------
    //  First clean up the torques for all the processors
    //--------------------------------------------------------------------------

    for (int i = 0; i < numAtoms; i++)
    {
        trq[i].x = 0;
        trq[i].y = 0;
        trq[i].z = 0;
    }

    for (Int i = 0; i < numAtoms; i++)
        {
	    molIDi = atoms[i].molID;
            qi = atoms[i].scaledCharge/SQRTCOULOMBCONSTANT;
            dipi = myL2G->dip[i];

	   //-------------------------------------------------------------------------------
	   // i-index Induced Dipole
	   //-------------------------------------------------------------------------------
	
            i3  = 3 * i;
	    indPi.x =  Pi[i3]; 
	    indPi.y =  Pi[i3+1]; 
	    indPi.z =  Pi[i3+2]; 

	   //-------------------------------------------------------------------------------
	   // Self-Correction Energy
	   //-------------------------------------------------------------------------------

            selfterm = (dipi*indPi)*alpha2PI*alpha2/3;
            selfEnergy += -selfterm*COUlOMBS;

	   //-------------------------------------------------------------------------------
	   // Self-Correction to the Torque
	   //-------------------------------------------------------------------------------
	   
	    if (computeTorques)
            {
		term = COUlOMBS*(2.0/3.0)*alpha2PI*alpha2;
		trq[i] = term*( cross(dipi, indPi) );

	        totTrqI += trq[i].length();
	    }
	}

	//----------------------------------------------------------------
	// Torques into Forces
	//----------------------------------------------------------------

	if (computeTorques){
	    double dUdV = 0.0;
            torques2forces(dUdV);

	    selfderiv += dUdV;
	}
}

//Jc:  compute real term in Ewald summation,it is necessary as well as the reciprocal term 
void Induction::real_term(Double &realEnergy, double &realderiv, double &real2ndDeriv)
{

    #ifdef DEBUG2
        DEBUGMSG("Computing Induction real space");
    #endif

    Atom *atomi, *atomj;
    Molecule *moli, *molj;
    Vector3 ri, rj, fi, fj, fij, rij, rji, rijm, rim, rjm;
    int size;
    int molIDi, molIDj, i3, j3;
    double qi, qj, qij;
    double dipxi, dipyi, dipzi;
    double dipxj, dipyj, dipzj;
    double Qxxi, Qxyi, Qxzi, Qyyi, Qyzi, Qzzi; 
    double Qyxi, Qzxi, Qzyi; 
    double Qxxj, Qxyj, Qxzj, Qyyj, Qyzj, Qzzj; 
    double Qyxj, Qzxj, Qzyj; 
    double G0_i, G1_i, G2_i, G3_i, G4_i;
    double G0_j, G1_j, G2_j, G3_j, G4_j;
    double G1_ip, G2_ip;
    double B0, B1, B2, B3, B4, B5, Uen;
    double lambda3, lambda5, lambda7, lambda9;
    double gradlam3, gradlam5, gradlam7, gradlam9;
    double gradB0, gradB1, gradB2, gradB3, scale;
    double gradB1u, gradB2u, gradB3u;
    double invr, invr2, invr3, invr5, invr7, invr9;
    double B0d, B1d, B2d, B3d, B4d;
    double B1du, B2du, B3du;
    double r2, rijRij, r, u3;
    double QiI, QjI, QiRij, QjRji, QipjXrij, QjpiXrji;
    double expau3, au3, a2u6, a3u9, a4u12; 
    double TorqueArray[numAtoms][3], tempArray[numAtoms][3];
    Vector3 dipi, dipj;
    Vector3 indPi, indPj;
    Vector3 gradG1_i, gradG2_i, gradG3_i;
    Vector3 gradG1_j, gradG2_j, gradG3_j;
    Vector3 gradG2_ip, gradG2_jp;
    Vector3 Gt1i, Gt2i, Gt3i;
    Vector3 Gt1j, Gt2j, Gt3j;
    Vector3 QiXrij, QjXrji, QiXpj, QjXpi;
    
    realEnergy = 0.0;
    realderiv  = 0.0;
    real2ndDeriv  = 0.0;

   //--------------------------------------------------------------------------
   //  First clean up the torques for all the processors
   //--------------------------------------------------------------------------

    for (int i = 0; i < numAtoms; i++)
    {
        trq[i].x = 0;
        trq[i].y = 0;
        trq[i].z = 0;
    }

// Jc: palallel code Created by Jinahui LI

   int numDecompose = (int) numAtoms/psize;
   if ((numAtoms % psize) == 0) numDecompose--;

   for(int i = 0; i <= numDecompose;i++)     
   { 
     int atomIndex = i*psize + prank;  
     if(atomIndex <numAtoms ) 
     {

//    for (Int i = 0; i < numAtoms; i++)    // Jc: ZW code
//    {                                     // Jc: ZW code
        Vector3 fi;
        atomi = &atoms[atomIndex];
        
        molIDi = atomi->molID;

        size = atomi->get_list_size();

        if(!strcmp(myMols[molIDi].molName, "MCY") && atomi->atomType == 1)            // Jc: O atom change to Dummy parameters
        { 
          ri.x = Dummy[molIDi*3];
          ri.y = Dummy[molIDi*3+1];
          ri.z = Dummy[molIDi*3+2];
        }
        else
        { 
          ri = atomi->position;
        }

  	rim = myEnsemble->molecules[molIDi].massCenter;
        qi = atomi->scaledCharge/SQRTCOULOMBCONSTANT;

	//--------------------------------------------------------------------------
	// Next we will transform the multipole from the local frame of reference 
	// of each site in a water molecule to the Global frame of reference
	// -------------------------------------------------------------------------
	
        dipi = myL2G->dip[atomIndex];

	//-----------------------------------------------------------------------
	// Call Traceless Quadrupoles from Local2Global
	//-----------------------------------------------------------------------

        Qxxi  = myL2G->Qxx[atomIndex];
        Qxyi  = myL2G->Qxy[atomIndex];
        Qxzi  = myL2G->Qxz[atomIndex];
        Qyyi  = myL2G->Qyy[atomIndex];
        Qyzi  = myL2G->Qyz[atomIndex];
        Qzzi  = myL2G->Qzz[atomIndex];

        Qyxi  = Qxyi;
        Qzxi  = Qxzi;
        Qzyi  = Qyzi;

	//-----------------------------------------------------------------------
	// i-index Induced Dipole	
	//-----------------------------------------------------------------------
	
        i3  = 3 * atomIndex;
	indPi.x =  Pi[i3]; 
	indPi.y =  Pi[i3+1]; 
	indPi.z =  Pi[i3+2]; 

        /*for (int j = 0; j < size; j++)
        {

            atomj = atomi->myPairList[j];*/
        for (int j = atomIndex +1; j < numAtoms; j++)
        {

            atomj = &atoms[j];
            molIDj = atomj->molID;

            if (!strcmp(myMols[molIDj].molName, "MCY") && atomj->atomType == 1)        // Jc: O atom change to Dummy parameters
            { 
              rj.x = Dummy[molIDj*3];
              rj.y = Dummy[molIDj*3+1];
              rj.z = Dummy[molIDj*3+2];
            }
            else
            { 
              rj = atomj->position;
            }

	    rijm = myEnsemble->molecules[molIDj].massCenter - rim;
            qj  = atomj->scaledCharge/SQRTCOULOMBCONSTANT;
            rij = ri - rj;

            dipj = myL2G->dip[atomj->atomID];

	    //-----------------------------------------------------------------------
	    // Call Traceless Quadrupoles from Local2Global
	    //-----------------------------------------------------------------------

            Qxxj  = myL2G->Qxx[atomj->atomID];
            Qxyj  = myL2G->Qxy[atomj->atomID];
            Qxzj  = myL2G->Qxz[atomj->atomID];
            Qyyj  = myL2G->Qyy[atomj->atomID];
            Qyzj  = myL2G->Qyz[atomj->atomID];
            Qzzj  = myL2G->Qzz[atomj->atomID];

            Qyxj  = Qxyj;
            Qzxj  = Qxzj;
            Qzyj  = Qyzj;

   	    //-----------------------------------------------------------------------
	    // j-index Induced Dipole	
	    //-----------------------------------------------------------------------
	
            j3  = 3 * atomj->atomID;
	    indPj.x =  Pi[j3]; 
	    indPj.y =  Pi[j3+1]; 
	    indPj.z =  Pi[j3+2]; 


            myEnsemble->apply_pbc(rij);
            myEnsemble->apply_pbc(rijm);

            r2 = rij.length2();


            if (r2 < realCut2)
            {
                rji = -rij;

		//----------------------------------------------------------------------------
		//  Constant needed for the G terms
		//----------------------------------------------------------------------------

		QiI   = Qxxi + Qyyi + Qzzi;
		QjI   = Qxxj + Qyyj + Qzzj;

		QiRij = Qxxi*rij.x*rij.x + Qxyi*rij.x*rij.y + Qxzi*rij.x*rij.z \
		      + Qyxi*rij.y*rij.x + Qyyi*rij.y*rij.y + Qyzi*rij.y*rij.z \ 
		      + Qzxi*rij.z*rij.x + Qzyi*rij.z*rij.y + Qzzi*rij.z*rij.z ; 

		QjRji = Qxxj*rij.x*rij.x + Qxyj*rij.x*rij.y + Qxzj*rij.x*rij.z \
		      + Qyxj*rij.y*rij.x + Qyyj*rij.y*rij.y + Qyzj*rij.y*rij.z \
		      + Qzxj*rij.z*rij.x + Qzyj*rij.z*rij.y + Qzzj*rij.z*rij.z ;

		QipjXrij = indPj.x*(Qxxi*rij.x + Qxyi*rij.y + Qxzi*rij.z) \ 
			 + indPj.y*(Qyxi*rij.x + Qyyi*rij.y + Qyzi*rij.z) \ 
			 + indPj.z*(Qzxi*rij.x + Qzyi*rij.y + Qzzi*rij.z) ;

		QjpiXrji = indPi.x*(Qxxj*rji.x + Qxyj*rji.y + Qxzj*rji.z) \ 
			 + indPi.y*(Qyxj*rji.x + Qyyj*rji.y + Qyzj*rji.z) \ 
			 + indPi.z*(Qzxj*rji.x + Qzyj*rji.y + Qzzj*rji.z) ;

		G0_i = 0.0;
		G1_i = - qi*indPj*rij - dipi*indPj;
		G2_i = - (dipi*rij)*(indPj*rij) - 2.0*QipjXrij - indPj*rij*QiI; 
		G3_i = - indPj*rij*QiRij;
		G4_i = 0.0;

		G0_j = 0.0;
		G1_j = - qj*indPi*rji - dipj*indPi;
		G2_j = - (dipj*rji)*(indPi*rji) - 2.0*QjpiXrji - indPi*rji*QjI; 
		G3_j = - indPi*rji*QjRji;
		G4_j = 0.0;

		//----------------------------------------------------------------------------
		//  Vectors needed for the gradG terms
		//----------------------------------------------------------------------------
		

		QiXrij.x = Qxxi*rij.x + Qxyi*rij.y + Qxzi*rij.z; 
		QiXrij.y = Qyxi*rij.x + Qyyi*rij.y + Qyzi*rij.z; 
		QiXrij.z = Qzxi*rij.x + Qzyi*rij.y + Qzzi*rij.z; 

		QjXrji.x = Qxxj*rji.x + Qxyj*rji.y + Qxzj*rji.z;
		QjXrji.y = Qyxj*rji.x + Qyyj*rji.y + Qyzj*rji.z;
		QjXrji.z = Qzxj*rji.x + Qzyj*rji.y + Qzzj*rji.z;

		QiXpj.x  = Qxxi*indPj.x + Qxyi*indPj.y + Qxzi*indPj.z; 
		QiXpj.y  = Qyxi*indPj.x + Qyyi*indPj.y + Qyzi*indPj.z; 
		QiXpj.z  = Qzxi*indPj.x + Qzyi*indPj.y + Qzzi*indPj.z; 

		QjXpi.x  = Qxxj*indPi.x + Qxyj*indPi.y + Qxzj*indPi.z;
		QjXpi.y  = Qyxj*indPi.x + Qyyj*indPi.y + Qyzj*indPi.z;
		QjXpi.z  = Qzxj*indPi.x + Qzyj*indPi.y + Qzzj*indPi.z;


		gradG1_i = - qi*indPj;
		gradG2_i = - (indPj*rij)*dipi - (dipi*rij)*indPj - 2*QiXpj - QiI*indPj; 
		gradG3_i = - QiRij*indPj - (2*indPj*rij)*QiXrij;

		gradG1_j = - qj*indPi;
		gradG2_j = - (indPi*rji)*dipj - (dipj*rji)*indPi - 2*QjXpi - QjI*indPi; 
		gradG3_j = - QjRji*indPi - (2*indPi*rji)*QjXrji;

		if (useConstraint){
		    rijRij = rij*rijm;
		}
		else {
		    rijRij = r2;
		}
                r  = sqrt(r2);
                u3 = r*r*r/(sqrt(atomi->polar * atomj->polar));

		B0 = erfc(alpha*r)/r;
		B1 = (-B0 - alpha2PI*exp(-alpha2*r2))/r2;
		B2 = (-3*B1 + 2*alpha2PI*alpha2*exp(-alpha2*r2) )/r2;
		B3 = (-5*B2 - 4*alpha2PI*alpha2*alpha2*exp(-alpha2*r2) )/r2;
		B4 = (-7*B3 + 8*alpha2PI*alpha2*alpha2*alpha2*exp(-alpha2*r2) )/r2;
		B5 = (-9*B4 - 16*alpha2PI*alpha2*alpha2*alpha2*alpha2*exp(-alpha2*r2) )/r2;

    		expau3  = exp(-a*u3);
    		au3     = a*u3;
    		a2u6    = au3*au3;
    		a3u9    = a2u6*au3;
    		a4u12   = a3u9*au3;

    		lambda3 = 1.0 - expau3;
		lambda5 = 1.0 - (1.0 + au3)*expau3;
		lambda7 = 1.0 - (1.0 + au3 + 3.0*a2u6/5.0)*expau3;
		lambda9 = 1.0 - (1.0 + au3 + (18.0*a2u6 + 9.0*a3u9)/35.0)*expau3;

    		gradlam3 = 3.0*au3*expau3/r2 ;
		gradlam5 = 3.0*a2u6*expau3/r2;
		gradlam7 = (-3.0*a2u6 + 9.0*a3u9)*expau3/(5.0*r2);
		gradlam9 = (-3.0*a2u6 - 27.0*a3u9 + 27.0*a4u12)*expau3/(35.0*r2);

                if (myEnsemble->exclusion_check(atomIndex, j))
                {
                    scale = 0.0;
                } else {
                    scale = 1.0;
                }

		if (dampAMOEBA) {
		    invr  = 1.0/r;
                    invr2 = invr*invr;
                    invr3 = invr2*invr;
                    invr5 = 3.0*invr3*invr2;
                    invr7 = 5.0*invr5*invr2;
                    invr9 = 7.0*invr7*invr2;

                    B1d = B1 + (1.0 - scale*lambda3)*invr3;
                    B2d = B2 - (1.0 - scale*lambda5)*invr5;
                    B3d = B3 + (1.0 - scale*lambda7)*invr7;

                    gradB1 = B2 - scale*gradlam3*invr3 - (1.0 - scale*lambda3)*invr5;
                    gradB2 = B3 + scale*gradlam5*invr5 + (1.0 - scale*lambda5)*invr7;
                    gradB3 = B4 - scale*gradlam7*invr7 - (1.0 - scale*lambda7)*invr9;

                    B1du = B1 + (1.0 - lambda3)*invr3;
                    B2du = B2 - (1.0 - lambda5)*invr5;
                    B3du = B3 + (1.0 - lambda7)*invr7;

                    gradB1u = B2 - gradlam3*invr3 - (1.0 - lambda3)*invr5;
                    gradB2u = B3 + gradlam5*invr5 + (1.0 - lambda5)*invr7;
                    gradB3u = B4 - gradlam7*invr7 - (1.0 - lambda7)*invr9;
		} else {
		    B1d = B1*lambda3; 
		    B2d = B2*lambda3 + B1*gradlam3;
		    B3d = B3*lambda3 + B2*2.0*gradlam3 + B1*(gradlam3 - 3.0*gradlam5)/r2;
		    B4d = B4*lambda3 + B3*3.0*gradlam3 + B2*(3.0*gradlam3 - 9.0*gradlam5)/r2 + B1*(-gradlam3 - 6.0*gradlam5 + 15.0*gradlam7)/(r2*r2);
		}

		//------------------------------------------------------------------------------------------------------
		//  Total Real-Space Energy as 0.5* sum_i sum_j=/=i sum_l G_lij * Blij
		//------------------------------------------------------------------------------------------------------
		if (dampAMOEBA) {
                    Uen = 0.5*((G1_i + G1_j)*B1d + (G2_i + G2_j)*B2d + (G3_i + G3_j)*B3d)*COUlOMBS;
		} else{	
		    Uen = 0.5*((G1_i + G1_j)*B1d + (G2_i + G2_j)*B2d + (G3_i + G3_j)*B3d)*scale*COUlOMBS;
		}

                realEnergy  += Uen;
		real2ndDeriv += 0.0;

		//---------------------------------------------------------------------------------------------------
		//  Total Real-Space Force as sum_j=/=i sum_l gradG_lij * Blij + G_lij * Bl+1 * rij
		//---------------------------------------------------------------------------------------------------
		
		G1_ip = - indPi*indPj;
		G2_ip = - (indPi*rij)*(indPj*rij); 
		gradG2_ip = - (indPj*rij)*indPi; 
		gradG2_jp = - (indPi*rji)*indPj; 

		if (dampAMOEBA) {
		    if (doMinimization) {
                        fi = - (B1d*gradG1_i + B2d*gradG2_i + B3d*gradG3_i + B2du*gradG2_ip) \
                             - (G1_i*gradB1 + G2_i*gradB2 + G3_i*gradB3 + G1_ip*gradB1u + G2_ip*gradB2u)*rij;

                        fj = - (B1d*gradG1_j + B2d*gradG2_j + B3d*gradG3_j + B2du*gradG2_jp) \
                             - (G1_j*gradB1 + G2_j*gradB2 + G3_j*gradB3)*rji;
                    } else {
                        fi = - (B1d*gradG1_i + B2d*gradG2_i + B3d*gradG3_i) \
                             - (G1_i*gradB1 + G2_i*gradB2 + G3_i*gradB3)*rij;

                        fj = - (B1d*gradG1_j + B2d*gradG2_j + B3d*gradG3_j) \
                             - (G1_j*gradB1 + G2_j*gradB2 + G3_j*gradB3)*rji;
                    }
		} else {
                    if (doMinimization) {
                        fi = - (B1d*gradG1_i*scale + B2d*gradG2_i*scale + B3d*gradG3_i*scale + B2d*gradG2_ip) \
                             - (G1_i*B2d*scale + G2_i*B3d*scale + G3_i*B4d*scale + G1_ip*B2d + G2_ip*B3d)*rij;

                        fj = - (B1d*gradG1_j*scale + B2d*gradG2_j*scale + B3d*gradG3_j*scale + B2d*gradG2_jp) \
                             - (G1_j*B2d*scale + G2_j*B3d*scale + G3_j*B4d*scale)*rji;
		    } else {
                        fi = - (B1d*gradG1_i + B2d*gradG2_i + B3d*gradG3_i)*scale \
                             - (G1_i*B2d + G2_i*B3d + G3_i*B4d)*rij*scale     ;

                        fj = - (B1d*gradG1_j + B2d*gradG2_j + B3d*gradG3_j)*scale \
                             - (G1_j*B2d + G2_j*B3d + G3_j*B4d)*rji*scale     ;
		    }
		}

		fij = (fi - fj); 
		realderiv += -fij*rij*COUlOMBS; 

// Jc: distribute force among molecule atoms, when dummy charge is calulated
// Jc: the calculated results are not improving. so there is no force distribution 
// Jc: in this calculation
		moli =  &myMols[molIDi]; 
		molj =  &myMols[molIDj]; 
                Atom *atomh1i = moli->myAtoms[0];
                Atom *atomh2i = moli->myAtoms[2];

                Atom *atomh1j = molj->myAtoms[0];
                Atom *atomh2j = molj->myAtoms[2];

                if ( !strcmp(myMols[molIDi].molName, "MCY") &&  (atomi->atomType == 1) ){
                        atomi->force += 0.5431676*fij*COUlOMBS;
                        atomh1i->force += 0.5*0.4568324*fij*COUlOMBS;
                        atomh2i->force += 0.5*0.4568324*fij*COUlOMBS;
                }
                else {
			atomi->force += fij*COUlOMBS;
                }

                if ( !strcmp(myMols[molIDj].molName, "MCY") && (atomj->atomType == 1) ){
                        atomj->force -= 0.5431676*fij*COUlOMBS;
                        atomh1j->force -= 0.5*0.4568324*fij*COUlOMBS;
                        atomh2j->force -= 0.5*0.4568324*fij*COUlOMBS;
                }
                else {
                        atomj->force -= fij*COUlOMBS;
                }

//                if(atomj->atomType == 0)  //Jc: H atoms
//                {
                  //atomj->force -= fij;
/*                }
                else
                {
                  double b    =  R_OP*Bohr_Radius/1.0e-9;
                  double gama =  -b/(0.09572 * cos(0.912109067)); // gama= 0.456840908

                  Vector3 p, rod, ro, rd;
                  ro = atomj->position; 
                  rd.x = Dummy[molIDj*3];
                  rd.y = Dummy[molIDj*3+1];
                  rd.z = Dummy[molIDj*3+2];
                  rod  = rd - ro;
                  p    = (rod*fij)/(b*b)*rod;
                  p    = gama*(fij - p);
                  atomj->force             -= fij -  p;
                  atoms[molIDj*3].force    -= 0.5* p  ;
                  atoms[molIDj*3 +2].force -= 0.5* p  ;
*/
/*                  atoms[molIDj*3].force     -= fij*0.068455;
                  atomj->force              -= fij*0.86309;
                  atoms[molIDj*3 + 2].force -= fij*0.068455;   
*/
//                } 
                virial[XX] += fij.x * rij.x;
                virial[XY] += fij.y * rij.x;
                virial[XZ] += fij.z * rij.x;
                virial[YX] += fij.x * rij.y;
                virial[YY] += fij.y * rij.y;
                virial[YZ] += fij.z * rij.y;
                virial[ZX] += fij.x * rij.z;
                virial[ZY] += fij.y * rij.z;
                virial[ZZ] += fij.z * rij.z;

	        //-------------------------------------------------------------------------------------------
	        // Computation of Torques
	        //-------------------------------------------------------------------------------------------
 
	        if (computeTorques)
                {
		    double Qi[3][3], Rij[3][3],pjrij[3][3], rijpj[3][3];
		    double Qj[3][3], pirji[3][3], rjipi[3][3];

		    Vector3 QiXpjrij, QiXrijpj, QiXRij;
		    Vector3 QjXpirji, QjXrjipi, QjXRji;
		    Vector3 ones(1,1,1);

		    Qi[0][0] = Qxxi;
		    Qi[0][1] = Qxyi;
		    Qi[0][2] = Qxzi;
		    Qi[1][0] = Qyxi;
		    Qi[1][1] = Qyyi;
		    Qi[1][2] = Qyzi;
		    Qi[2][0] = Qzxi;
		    Qi[2][1] = Qzyi;
		    Qi[2][2] = Qzzi;

		    Qj[0][0] = Qxxj;
		    Qj[0][1] = Qxyj;
		    Qj[0][2] = Qxzj;
		    Qj[1][0] = Qyxj;
		    Qj[1][1] = Qyyj;
		    Qj[1][2] = Qyzj;
		    Qj[2][0] = Qzxj;
		    Qj[2][1] = Qzyj;
		    Qj[2][2] = Qzzj;

		    Rij[0][0] = rij.x*rij.x;
		    Rij[0][1] = rij.x*rij.y;
		    Rij[0][2] = rij.x*rij.z;
		    Rij[1][0] = rij.y*rij.x;
		    Rij[1][1] = rij.y*rij.y;
		    Rij[1][2] = rij.y*rij.z;
		    Rij[2][0] = rij.z*rij.x;
		    Rij[2][1] = rij.z*rij.y;
		    Rij[2][2] = rij.z*rij.z;

		    pjrij[0][0] = indPj.x*rij.x;
		    pjrij[0][1] = indPj.x*rij.y;
		    pjrij[0][2] = indPj.x*rij.z;
		    pjrij[1][0] = indPj.y*rij.x;
		    pjrij[1][1] = indPj.y*rij.y;
		    pjrij[1][2] = indPj.y*rij.z;
		    pjrij[2][0] = indPj.z*rij.x;
		    pjrij[2][1] = indPj.z*rij.y;
		    pjrij[2][2] = indPj.z*rij.z;

		    rijpj[0][0] = pjrij[0][0];
		    rijpj[0][1] = pjrij[1][0];
		    rijpj[0][2] = pjrij[2][0];
		    rijpj[1][0] = pjrij[0][1];
		    rijpj[1][1] = pjrij[1][1];
		    rijpj[1][2] = pjrij[2][1];
		    rijpj[2][0] = pjrij[0][2];
		    rijpj[2][1] = pjrij[1][2];
		    rijpj[2][2] = pjrij[2][2];

		    pirji[0][0] = indPi.x*rji.x;
		    pirji[0][1] = indPi.x*rji.y;
		    pirji[0][2] = indPi.x*rji.z;
		    pirji[1][0] = indPi.y*rji.x;
		    pirji[1][1] = indPi.y*rji.y;
		    pirji[1][2] = indPi.y*rji.z;
		    pirji[2][0] = indPi.z*rji.x;
		    pirji[2][1] = indPi.z*rji.y;
		    pirji[2][2] = indPi.z*rji.z;

		    rjipi[0][0] = pirji[0][0];
		    rjipi[0][1] = pirji[1][0];
		    rjipi[0][2] = pirji[2][0];
		    rjipi[1][0] = pirji[0][1];
		    rjipi[1][1] = pirji[1][1];
		    rjipi[1][2] = pirji[2][1];
		    rjipi[2][0] = pirji[0][2];
		    rjipi[2][1] = pirji[1][2];
		    rjipi[2][2] = pirji[2][2];

		    QiXpjrij = cross(Qi, pjrij, ones);
		    QiXrijpj = cross(Qi, rijpj, ones);
		    QiXRij   = cross(Qi, Rij, ones);		    


		    QjXpirji = cross(Qj, pirji, ones);
		    QjXrjipi = cross(Qj, rjipi, ones);
		    QjXRji   = cross(Qj, Rij, ones);		// because Rji = Rij

    		    Gt1i = cross(dipi, indPj);
		    Gt2i = (indPj*rij)*(cross(dipi, rij)) + 2.0*QiXpjrij + 2.0*QiXrijpj;
		    Gt3i = 2.0*(indPj*rij)*QiXRij;

    		    Gt1j = cross(dipj, indPi);
		    Gt2j = (indPi*rji)*(cross(dipj, rji)) + 2.0*QjXpirji + 2.0*QjXrjipi;
		    Gt3j = 2.0*(indPi*rji)*QjXRji;

		    if (dampAMOEBA) {
			trq[atomi->atomID] += COUlOMBS*(Gt1i*B1d + Gt2i*B2d + Gt3i*B3d);
                        trq[atomj->atomID] += COUlOMBS*(Gt1j*B1d + Gt2j*B2d + Gt3j*B3d);
		    } else {
        	        trq[atomi->atomID] += COUlOMBS*(Gt1i*B1d + Gt2i*B2d + Gt3i*B3d)*scale;
        	        trq[atomj->atomID] += COUlOMBS*(Gt1j*B1d + Gt2j*B2d + Gt3j*B3d)*scale;
		    }

		}  // If torques

            } //If cut-off
        }  // j-loop

//        if(atomi->atomType == 0)
//        {
         // atomi->force += fi;
/*        }
        else
        {
          double b    = R_OP*Bohr_Radius/1.0e-9;         // 0.026765499
          double gama = -b/(0.09572 * cos(0.912109067)); // gama= 0.456840908

          Vector3 p, rod, ro, rd;
          ro = atomi->position; 
          rd.x = Dummy[molIDi*3];
          rd.y = Dummy[molIDi*3+1];
          rd.z = Dummy[molIDi*3+2];
          rod  = rd - ro;
          p    = (rod*fi)/(b*b)*rod;
          p    = gama*(fi - p);

          atomi->force             += fi -  p;
          atoms[molIDi*3].force    += 0.5* p;
          atoms[molIDi*3 +2].force += 0.5* p;
*/
/*          atoms[molIDi*3].force     += fi*0.068455;
          atomi->force              += fi*0.86309;
          atoms[molIDi*3 + 2].force += fi*0.068455;   
*/
//        } 
     }
   }

   double temp;
   //MPI::COMM_WORLD.Reduce(&realEnergy,&temp,1,MPI::DOUBLE,MPI::SUM,0); Updated to MPI-3 by elton
   MPI_Reduce(&realEnergy,&temp,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   realEnergy = temp;

   temp = 0.0;
   MPI_Allreduce(&realderiv,&temp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   realderiv = temp;

   temp = 0.0;
   MPI_Reduce(&real2ndDeriv,&temp,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   real2ndDeriv = temp;


    if (computeTorques){
        for (int i = 0; i < numAtoms; i++)
        {
            TorqueArray[i][0] = trq[i].x;
            TorqueArray[i][1] = trq[i].y;
            TorqueArray[i][2] = trq[i].z;

            tempArray[i][0] = 0.0;
            tempArray[i][1] = 0.0;
            tempArray[i][2] = 0.0;
        }

        MPI_Allreduce(&TorqueArray,&tempArray,nAtomsx3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);

        for (int i = 0; i < numAtoms; i++)
        {
            trq[i].x = tempArray[i][0];
            trq[i].y = tempArray[i][1];
            trq[i].z = tempArray[i][2];

	    totTrqI += trq[i].length();
        }

	//----------------------------------------------------------------
	// Torques into Forces
	//----------------------------------------------------------------

	double dUdV = 0.0;
        torques2forces(dUdV);

	realderiv += dUdV;

    }

}

void Induction::reciprocal_term(Double &pe, double &longderiv, double &long2ndDeriv)
{    
    if(doUpdate) {
        volume = myEnsemble->volume;
        volume2 = volume*volume;
        volumer = 1.0/volume;
    }

    Double pi4Volumer = 4.0*PI*volumer;
    Double pi8Volumer = 2.0*pi4Volumer;
    Double tmpenergy = 0.0;
    Molecule *moli;
    clock_t start,end;
    Vector3 ri,rj,rim;
    pe = 0.0;
    longderiv = 0.0;
    long2ndDeriv = 0.0;
    Vector3 dipi, indPi;
    Vector3 kj, cosrR, sinrR;
    double Qxxi, Qxyi, Qxzi, Qyyi, Qyzi, Qzzi; 
    double Qyxi, Qzxi, Qzyi; 
    double UmqQ, Ump, dUdV3, dUdVpQ;
    double Upp, Uup;
    double qi;
    double qiQisinSum, qiQicosSum, miksinSum, mikcosSum;
    double piksinSum, pikcosSum, Q2sinSum, Q2cosSum;
    double FpiqQ, FqQpj, Fpmj, Fupj, val;
    double Sk, dUdV1, dUdV2;
    double TorqueArray[numAtoms][3], tempArray[numAtoms][3];
    double rk, QiK, cosAqiQi, sinAqiQi;
    double mikcos, miksin, pikcos, piksin, Q2sin, Q2cos;
    double b, Ff, Ffk2, d2UV1, d2UrR;
    double d2u1, d2u2;

    int molIDi, i3;

   //--------------------------------------------------------------------------
   //  First clean up the torques for all the processors
   //--------------------------------------------------------------------------

    for (int i = 0; i < numAtoms; i++)
    {
        trq[i].x = 0;
        trq[i].y = 0;
        trq[i].z = 0;
    }

// Jc: parallel code  created by Jianhui
    int numDecompose = (int) kVector.size()/psize;

    if ((kVector.size() % psize) == 0) numDecompose--;
    for(int k = 0; k <= numDecompose;k++)     
    { 
      int j = k*psize + prank; 
      if(j < kVector.size()) 
      {
    
//    for (unsigned int j = 0; j < kVector.size(); j++)   // Jc: ZW code
//    {                                                   // Jc: ZW code  
        qiQisinSum = 0.0; 
        qiQicosSum = 0.0;
        miksinSum  = 0.0; 
        mikcosSum  = 0.0;
        piksinSum  = 0.0; 
        pikcosSum  = 0.0;
        Q2sinSum   = 0.0;
        Q2cosSum   = 0.0;
        kj = kVector[j];
        Double coeff = 2.0*(alphaR4 + 1.0/kModulu[j]);         // for virial compute

	cosrR.x = 0.0, cosrR.y = 0.0, cosrR.z = 0.0;
	sinrR.x = 0.0, sinrR.y = 0.0, sinrR.z = 0.0;
 
        for (Int i = 0; i < numAtoms; i++)
        {
	    molIDi = atoms[i].molID;
	    rim = myEnsemble->molecules[molIDi].massCenter;
            qi = atoms[i].scaledCharge/SQRTCOULOMBCONSTANT;

            dipi = myL2G->dip[i];

	    //-----------------------------------------------------------------------
	    // Call Traceless Quadrupoles from Local2Global
	    //-----------------------------------------------------------------------

            Qxxi  = myL2G->Qxx[i];
            Qxyi  = myL2G->Qxy[i];
            Qxzi  = myL2G->Qxz[i];
            Qyyi  = myL2G->Qyy[i];
            Qyzi  = myL2G->Qyz[i];
            Qzzi  = myL2G->Qzz[i];

            Qyxi  = Qxyi;
            Qzxi  = Qxzi;
            Qzyi  = Qyzi;

	    //-----------------------------------------------------------------------
	    // i-index Induced Dipole	
	    //-----------------------------------------------------------------------
	
            i3  = 3 * i;
	    indPi.x =  Pi[i3]; 
	    indPi.y =  Pi[i3+1]; 
	    indPi.z =  Pi[i3+2]; 


            if(!strcmp(myMols[atoms[i].molID].molName, "MCY") && atoms[i].atomType == 1) // Jc: oxygen atom
            {
              ri.x = Dummy[3*(atoms[i].molID)];
              ri.y = Dummy[3*(atoms[i].molID)+1];
              ri.z = Dummy[3*(atoms[i].molID)+2];
            }
            else
            { 
              ri = atoms[i].position;
            }
            // na q er myEnsemble->apply_pbc(ri);
            // na q er myEnsemble->apply_pbc(rim);
            rk   = ri*kj;
	    QiK = Qxxi*kj.x*kj.x + Qxyi*kj.x*kj.y + Qxzi*kj.x*kj.z \
                + Qyxi*kj.y*kj.x + Qyyi*kj.y*kj.y + Qyzi*kj.y*kj.z \
                + Qzxi*kj.z*kj.x + Qzyi*kj.z*kj.y + Qzzi*kj.z*kj.z ; 

            cosAqiQi = (qi - QiK)*cos(rk);
            sinAqiQi = (qi - QiK)*sin(rk);
	    mikcos   = (dipi*kj)*cos(rk);
	    miksin   = (dipi*kj)*sin(rk);
	    pikcos   = (indPi*kj)*cos(rk);
	    piksin   = (indPi*kj)*sin(rk);
            Q2sin    = 2.0*QiK*sin(rk);
            Q2cos    = 2.0*QiK*cos(rk);

    	    qQcos_i[i] = cosAqiQi;
    	    qQsin_i[i] = sinAqiQi;
    	    mkcos_i[i] = mikcos;
    	    mksin_i[i] = miksin;
    	    pkcos_i[i] = pikcos;
    	    pksin_i[i] = piksin;

            qiQicosSum += cosAqiQi;
            qiQisinSum += sinAqiQi; 
            mikcosSum  += mikcos;
            miksinSum  += miksin; 
            pikcosSum  += pikcos;
            piksinSum  += piksin; 
	    Q2sinSum   += Q2sin;
            Q2cosSum   += Q2cos;


	    cosrR  += cosAqiQi*(ri - rim); // deprecated
	    sinrR  += sinAqiQi*(ri - rim); // deprecated

	    //------------------------------------------------------------------------------------------
	    // cross products for torque computations 
	    //------------------------------------------------------------------------------------------

	    if (computeTorques)
            {
                double Kj[3][3], Qi[3][3];
                Vector3 ones(1,1,1);

                Kj[0][0] = kj.x*kj.x;
                Kj[0][1] = kj.x*kj.y;
                Kj[0][2] = kj.x*kj.z;
                Kj[1][0] = kj.y*kj.x;
                Kj[1][1] = kj.y*kj.y;
                Kj[1][2] = kj.y*kj.z;
                Kj[2][0] = kj.z*kj.x;
                Kj[2][1] = kj.z*kj.y;
                Kj[2][2] = kj.z*kj.z;

                Qi[0][0] = Qxxi;
                Qi[0][1] = Qxyi;
                Qi[0][2] = Qxzi;
                Qi[1][0] = Qyxi;
                Qi[1][1] = Qyyi;
                Qi[1][2] = Qyzi;
                Qi[2][0] = Qzxi;
                Qi[2][1] = Qzyi;
                Qi[2][2] = Qzzi;

		Vector3 kXmi = cross(kj, dipi);
		Vector3 KXQi = cross(Kj, Qi, ones);

    		kXmsina[i] = kXmi*sin(rk);
    		kXmcosa[i] = kXmi*cos(rk);
    		KXQcosa[i] = KXQi*cos(rk);
    		KXQsina[i] = KXQi*sin(rk);
	    }
        }

	//-------------------------------------------------------------------------------------------------
	// Reciprocal-space Energy computation between all multipoles
	//-------------------------------------------------------------------------------------------------
	
	UmqQ  = pi4Volumer*expK[j]*(qiQisinSum*pikcosSum - qiQicosSum*piksinSum);
	Ump   = pi4Volumer*expK[j]*(pikcosSum*mikcosSum + piksinSum*miksinSum);

        tmpenergy = (UmqQ + Ump)*COUlOMBS;
        pe += tmpenergy;

	//-------------------------------------------------------------------------------------------------
	// Reciprocal-space First Volume derivative 
	//-------------------------------------------------------------------------------------------------

	dUdVpQ = pi4Volumer*expK[j]*(Q2sinSum*pikcosSum - Q2cosSum*piksinSum);
	if (doMinimization) {
	    Upp  =  pi4Volumer*expK[j]*(pikcosSum*pikcosSum + piksinSum*piksinSum);
	} else {
	    Upp  = 0.0;
	}
	
	Sk = 2.0*UmqQ + 2.0*Ump + Upp;
	dUdV1 = -Sk;
	dUdV2 = Sk*kj*kj*2.0*alphaR4;
	dUdV3  = -2.0*UmqQ + 2.0*dUdVpQ - 4.0*Ump - 2.0*Upp;

	longderiv += (dUdV1 + dUdV2 + dUdV3)*COUlOMBS;
	long2ndDeriv += 4.0*tmpenergy - 7.0*tmpenergy*kj*kj*2.0*alphaR4 + tmpenergy*kj*kj*2.0*alphaR4*kj*kj*2.0*alphaR4;

        longVirial[XX] += tmpenergy - tmpenergy*coeff*kj.x*kj.x;
        longVirial[YY] += tmpenergy - tmpenergy*coeff*kj.y*kj.y;
        longVirial[ZZ] += tmpenergy - tmpenergy*coeff*kj.z*kj.z;
        longVirial[XY] -= tmpenergy*coeff*kj.x*kj.y;
        longVirial[XZ] -= tmpenergy*coeff*kj.x*kj.z;
        longVirial[YZ] -= tmpenergy*coeff*kj.y*kj.z;
        

	//-----------------------------------------------------------------------------------------------
	// Force Computation of the Reciprocal-space contribution
	//-----------------------------------------------------------------------------------------------

        b = pi8Volumer*expK[j];
	Ff = 0.0; Ffk2 = 0.0; d2UV1 = 0.0; d2UrR = 0.0;

        for (Int i = 0; i < numAtoms; i++)
        {
	  FpiqQ =  b*(pkcos_i[i]*qiQicosSum + pksin_i[i]*qiQisinSum); 
	  FqQpj = -b*(qQcos_i[i]*pikcosSum + qQsin_i[i]*piksinSum) ;
	  Fupj  =  b*(mksin_i[i]*pikcosSum - mkcos_i[i]*piksinSum) ;

	  if (doMinimization) {
	      Fpmj  =  b*(pksin_i[i]*(mikcosSum + pikcosSum) - pkcos_i[i]*(miksinSum + piksinSum) ) ;
	  } else {
	      Fpmj  =  b*(pksin_i[i]*mikcosSum - pkcos_i[i]*miksinSum) ;
	  }

          val = COUlOMBS*(FpiqQ + FqQpj + Fpmj + Fupj);  // Total Force

          d2u1 = b*kj*kj*(qQcos_i[i]*qiQicosSum + qQsin_i[i]*qiQisinSum);  // Equation for monopole only
	  molIDi = atoms[i].molID;
	  rim = myEnsemble->molecules[molIDi].massCenter;
	  ri = atoms[i].position;

          d2u2 = b*kj*kj*(qQcos_i[i]*cosrR + qQsin_i[i]*sinrR)*(ri - rim);  // Equation for monopole only

	      //------------------------------------------------------------------------------------------
	      // if the MCY model is tested, we need to redistribute the forces in order
	      // to keep thermal equilibrium
	      //------------------------------------------------------------------------------------------

	        moli =  &myMols[molIDi];
                Atom *atomh1i = moli->myAtoms[0];
                Atom *atomh2i = moli->myAtoms[2];

                if (!strcmp(myMols[atoms[i].molID].molName, "MCY") && (atoms[i].atomType == 1) ){
                        atoms[i].force += 0.5431676*kVector[j]*val;
                        atomh1i->force += 0.5*0.4568324*kVector[j]*val;
                        atomh2i->force += 0.5*0.4568324*kVector[j]*val;
                }
                else {
                        atoms[i].force += val*kVector[j];
                }

	     Ff +=  val * ( kVector[j]*(ri - rim) );
	     Ffk2 += 2.0*val*kj*kj*2.0*alphaR4 * ( kVector[j]*(ri - rim) ); 
	     d2UV1 += d2u1*(ri - rim)*(ri - rim);
	     d2UrR += d2u2;

	    //-------------------------------------------------------------------------------------------
	    // Computation of Torques
	    //-------------------------------------------------------------------------------------------
 
	    if (computeTorques)
            {
		Vector3 taum = b*( kXmcosa[i]*pikcosSum + kXmsina[i]*piksinSum);
		Vector3 tauQ = 2.0*b*( KXQcosa[i]*piksinSum - KXQsina[i]*pikcosSum);

		trq[i] += COUlOMBS*(taum + tauQ); 
	    }
        }

	if (useConstraint){
	    longderiv += Ff;
	    long2ndDeriv += -6.0*Ff + Ffk2 - d2UV1 + d2UrR;
	}

      }                                                                  
    }

    longVirial[YX] = longVirial[XY];
    longVirial[ZX] = longVirial[XZ];
    longVirial[ZY] = longVirial[YZ];

    double temp;
//    MPI::COMM_WORLD.Reduce(&pe,&temp,1,MPI::DOUBLE,MPI::SUM,0); MPI-3 by elton
    MPI_Reduce(&pe,&temp,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    pe = temp;
   
    temp = 0.0; 
    MPI_Allreduce(&longderiv,&temp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    longderiv = temp;

    temp = 0.0; 
    MPI_Reduce(&long2ndDeriv,&temp,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    long2ndDeriv = temp;


    if (computeTorques){
        for (int i = 0; i < numAtoms; i++)
        {
            TorqueArray[i][0] = trq[i].x;
            TorqueArray[i][1] = trq[i].y;
            TorqueArray[i][2] = trq[i].z;

            tempArray[i][0] = 0.0;
            tempArray[i][1] = 0.0;
            tempArray[i][2] = 0.0;
        }

        MPI_Allreduce(&TorqueArray,&tempArray,nAtomsx3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);

        for (int i = 0; i < numAtoms; i++)
        {
            trq[i].x = tempArray[i][0];
            trq[i].y = tempArray[i][1];
            trq[i].z = tempArray[i][2];

	    totTrqI += trq[i].length();
        }

	//----------------------------------------------------------------
	// Torques into Forces
	//----------------------------------------------------------------

	double dUdV = 0.0;
        torques2forces(dUdV);

	longderiv += dUdV;

    }

}

// intra-molecular self energy correction 
// Jc: it is one of the biggest energy contribution but some times, it is balanced off by 
// Jc: the point self-energy, the self point self-energy is calcualted only once in the EwaldForce.cpp
// Jc: by the compute() method
// Jc: the self-energy is not related to the any positions in the system, but only the chrages 
// Jc: however the intraMolecularEnergy, also called nolCorrection are charge-position depedent
// Jc: Since the positions of charges in MCY model is different from SPC/E, Dummy positions need to be 
// Jc: introduced to account for the molCorrection calculation
 
void Induction::intraMolCorrect(Double &intraMolEnergy, double &intraVircorrect,  double &molderiv, double &mol2ndDeriv)
{
    Double intraMolecularEnergy = 0.0;
    Double r, r2, qi, qj;
    Int idi, idj, ns;
    int molIDi, molIDj; 
    int i3, j3;
    double dipxi, dipyi, dipzi;
    double dipxj, dipyj, dipzj;
    double Qxxi, Qxyi, Qxzi, Qyyi, Qyzi, Qzzi; 
    double Qyxi, Qzxi, Qzyi; 
    double Qxxj, Qxyj, Qxzj, Qyyj, Qyzj, Qzzj; 
    double Qyxj, Qzxj, Qzyj; 
    double G0_i, G1_i, G2_i, G3_i, G4_i;
    double G0_j, G1_j, G2_j, G3_j, G4_j;
    double D0, D1, D2, D3, D4, D5;
    double QiI, QjI, QiRij, QjRji, QipjXrij, QjpiXrji;
    double TorqueArray[numAtoms][3], tempArray[numAtoms][3];
    Vector3 dipi, dipj;
    Vector3 indPi, indPj;
    Vector3 ri, rj, rij, rji;
    Vector3 gradG1_i, gradG2_i, gradG3_i;
    Vector3 gradG1_j, gradG2_j, gradG3_j;
    Vector3 Gt1i, Gt2i, Gt3i;
    Vector3 Gt1j, Gt2j, Gt3j;

    intraMolEnergy = 0;
    intraVircorrect = 0.0;
    molderiv = 0.0;
    mol2ndDeriv = 0.0;
    if(doUpdate) {
        volume = myEnsemble->volume;
        volume2 = volume*volume;
        volumer = 1.0/volume;
    }

   //--------------------------------------------------------------------------
   //  First clean up the torques for all the processors
   //--------------------------------------------------------------------------

    for (int i = 0; i < numAtoms; i++)
    {
        trq[i].x = 0;
        trq[i].y = 0;
        trq[i].z = 0;
    }


    // intra-molecular term. 
    // Molecular virial correction should be computed after all forces have been computed.  

    Int nmol = myEnsemble->nMols;

// Jc: parallel code
    int numDecompose = (int) nmol/psize;
    if ((nmol % psize) == 0) numDecompose--;
    for(int i = 0; i <= numDecompose;i++)     
    { 
      int molIndex = i*psize + prank; // Jc:i*psize + psize - 1- prank;
      if(molIndex <nmol ) 
      {

//    for(Int nm = 0; nm < nmol; nm++)        //Jc: ZW code
//    {                                       //Jc: ZW code
        ns = myEnsemble->molecules[molIndex].numAtoms;
        for (Int i = 0; i < ns; i++)
        {
            idi = myEnsemble->molecules[molIndex].myAtoms[i]->atomID;
            molIDi = myEnsemble->molecules[molIndex].myAtoms[i]->molID;
            qi = atoms[idi].scaledCharge/SQRTCOULOMBCONSTANT;

	    //--------------------------------------------------------------------------
	    // Next we will transform the multipole from the local frame of reference 
	    // of each site in a water molecule to the Global frame of reference
	    // -------------------------------------------------------------------------

            dipi = myL2G->dip[idi];

	    //-----------------------------------------------------------------------
	    // Call Traceless Quadrupoles from Local2Global
	    //-----------------------------------------------------------------------

            Qxxi  = myL2G->Qxx[idi];
            Qxyi  = myL2G->Qxy[idi];
            Qxzi  = myL2G->Qxz[idi];
            Qyyi  = myL2G->Qyy[idi];
            Qyzi  = myL2G->Qyz[idi];
            Qzzi  = myL2G->Qzz[idi];

            Qyxi  = Qxyi;
            Qzxi  = Qxzi;
            Qzyi  = Qyzi;

	    //-----------------------------------------------------------------------
	    // i-index Induced Dipole	
	    //-----------------------------------------------------------------------
	
            i3  = 3 * idi;
	    indPi.x =  Pi[i3]; 
	    indPi.y =  Pi[i3+1]; 
	    indPi.z =  Pi[i3+2]; 


// Jc: substitute the Oxygen position with the Dummdy Charge position
            ri = atoms[idi].position;
            if(!strcmp(myMols[molIDi].molName, "MCY") && atoms[idi].atomType ==1)
            { 
              ri.x = Dummy[molIDi*3];
              ri.y = Dummy[molIDi*3 +1];
              ri.z = Dummy[molIDi*3 +2]; 
            }
            for (Int j = i+1; j < ns; j++)
            {
                idj = myEnsemble->molecules[molIndex].myAtoms[j]->atomID;
                molIDj = myEnsemble->molecules[molIndex].myAtoms[j]->molID;
                if (myEnsemble->exclusion_check(idi, idj))
                {
	            if (molIDi != molIDj) {
                        ERRORMSG("Molecular Correction term: Different molecules are interacting. Check exclusion_check() pairs.");}

                    rj = atoms[idj].position;
                    if(!strcmp(myMols[molIDj].molName, "MCY") && atoms[idj].atomType ==1)
                    { 
                     rj.x = Dummy[3*molIDj];
                     rj.y = Dummy[3*molIDj +1];
                     rj.z = Dummy[3*molIDj +2]; 
                    }

                    rij = ri - rj;
                    myEnsemble->apply_pbc(rij);
                    r2 = rij.length2();
                    r = sqrt(r2);

                    qj = atoms[idj].scaledCharge/SQRTCOULOMBCONSTANT;
            	    dipj = myL2G->dip[idj];

	    	    //-----------------------------------------------------------------------
	    	    // Call Traceless Quadrupoles from Local2Global
	    	    //-----------------------------------------------------------------------

            	    Qxxj  = myL2G->Qxx[idj];
            	    Qxyj  = myL2G->Qxy[idj];
            	    Qxzj  = myL2G->Qxz[idj];
            	    Qyyj  = myL2G->Qyy[idj];
            	    Qyzj  = myL2G->Qyz[idj];
            	    Qzzj  = myL2G->Qzz[idj];

            	    Qyxj  = Qxyj;
            	    Qzxj  = Qxzj;
            	    Qzyj  = Qyzj;

   	    	    //-----------------------------------------------------------------------
	    	    // j-index Induced Dipole	
	    	    //-----------------------------------------------------------------------
	
            	    j3  = 3 * idj;
	    	    indPj.x =  Pi[j3]; 
	    	    indPj.y =  Pi[j3+1]; 
	    	    indPj.z =  Pi[j3+2]; 

		    //----------------------------------------------------------------------------
		    //  Constant needed for the G terms
		    //----------------------------------------------------------------------------

		    rji = -rij;

		    QiI   = Qxxi + Qyyi + Qzzi;
		    QjI   = Qxxj + Qyyj + Qzzj;

		    QiRij = Qxxi*rij.x*rij.x + Qxyi*rij.x*rij.y + Qxzi*rij.x*rij.z \
			  + Qyxi*rij.y*rij.x + Qyyi*rij.y*rij.y + Qyzi*rij.y*rij.z \ 
			  + Qzxi*rij.z*rij.x + Qzyi*rij.z*rij.y + Qzzi*rij.z*rij.z ; 

		    QjRji = Qxxj*rij.x*rij.x + Qxyj*rij.x*rij.y + Qxzj*rij.x*rij.z \
		          + Qyxj*rij.y*rij.x + Qyyj*rij.y*rij.y + Qyzj*rij.y*rij.z \
		          + Qzxj*rij.z*rij.x + Qzyj*rij.z*rij.y + Qzzj*rij.z*rij.z ;

		    QipjXrij = indPj.x*(Qxxi*rij.x + Qxyi*rij.y + Qxzi*rij.z) \ 
			     + indPj.y*(Qyxi*rij.x + Qyyi*rij.y + Qyzi*rij.z) \ 
			     + indPj.z*(Qzxi*rij.x + Qzyi*rij.y + Qzzi*rij.z) ;

		    QjpiXrji = indPi.x*(Qxxj*rji.x + Qxyj*rji.y + Qxzj*rji.z) \ 
			     + indPi.y*(Qyxj*rji.x + Qyyj*rji.y + Qyzj*rji.z) \ 
			     + indPi.z*(Qzxj*rji.x + Qzyj*rji.y + Qzzj*rji.z) ;

		    G0_i = 0.0;
		    G1_i = - qi*indPj*rij - dipi*indPj;
		    G2_i = - (dipi*rij)*(indPj*rij) - 2.0*QipjXrij - indPj*rij*QiI ; 
		    G3_i = - indPj*rij*QiRij;
		    G4_i = 0.0;

		    G0_j = 0.0;
		    G1_j = - qj*indPi*rji - dipj*indPi;
		    G2_j = - (dipj*rji)*(indPi*rji) - 2.0*QjpiXrji - indPi*rji*QjI ; 
		    G3_j = - indPi*rji*QjRji;
		    G4_j = 0.0;

		    //---------------------------------------------------------------------------------------------------
		    //  D constants for Induction Multipole Correction
		    //---------------------------------------------------------------------------------------------------

		    D0 = erf(alpha*r)/r;
		    D1 = (-D0 + alpha2PI*exp(-alpha2*r2))/r2;
		    D2 = (-3*D1 - 2*alpha2PI*alpha2*exp(-alpha2*r2) )/r2;
		    D3 = (-5*D2 + 4*alpha2PI*alpha2*alpha2*exp(-alpha2*r2) )/r2;
		    D4 = (-7*D3 - 8*alpha2PI*alpha2*alpha2*alpha2*exp(-alpha2*r2) )/r2;
		    D5 = (-9*D4 + 16*alpha2PI*alpha2*alpha2*alpha2*alpha2*exp(-alpha2*r2) )/r2;

                    intraMolecularEnergy += (G1_i + G1_j)*D1 + (G2_i + G2_j)*D2 + (G3_i + G3_j)*D3;

                    // intra-molecular self force, may be no need when using constraint method
		    if (!useConstraint){

			Vector3 QiXrij, QjXrji, QiXpj, QjXpi;

			QiXrij.x = Qxxi*rij.x + Qxyi*rij.y + Qxzi*rij.z; 
			QiXrij.y = Qyxi*rij.x + Qyyi*rij.y + Qyzi*rij.z; 
			QiXrij.z = Qzxi*rij.x + Qzyi*rij.y + Qzzi*rij.z; 

			QjXrji.x = Qxxj*rji.x + Qxyj*rji.y + Qxzj*rji.z;
			QjXrji.y = Qyxj*rji.x + Qyyj*rji.y + Qyzj*rji.z;
			QjXrji.z = Qzxj*rji.x + Qzyj*rji.y + Qzzj*rji.z;

			QiXpj.x  = Qxxi*indPj.x + Qxyi*indPj.y + Qxzi*indPj.z; 
			QiXpj.y  = Qyxi*indPj.x + Qyyi*indPj.y + Qyzi*indPj.z; 
			QiXpj.z  = Qzxi*indPj.x + Qzyi*indPj.y + Qzzi*indPj.z; 

			QjXpi.x  = Qxxj*indPi.x + Qxyj*indPi.y + Qxzj*indPi.z;
			QjXpi.y  = Qyxj*indPi.x + Qyyj*indPi.y + Qyzj*indPi.z;
			QjXpi.z  = Qzxj*indPi.x + Qzyj*indPi.y + Qzzj*indPi.z;


			gradG1_i = - qi*indPj;
			gradG2_i = - (indPj*rij)*(dipi /*+ indPi*/) - (dipi*rij)*indPj - 2*QiXpj - QiI*indPj; 
			gradG3_i = - QiRij*indPj - (2*indPj*rij)*QiXrij;

			gradG1_j = - qj*indPi;
			gradG2_j = - (indPi*rji)*(dipj /*+ indPj*/) - (dipj*rji)*indPi - 2*QjXpi - QjI*indPi; 
			gradG3_j = - QjRji*indPi - (2*indPi*rji)*QjXrji;

		        //-------------------------------------------------------------------------------------------
		        // Computation of Forces and dU/dV for pressure
		        // Note: The E_ind field contribution within the molecules is computed
		        // in the reciprocal part, and do counts as the mutual force, 
		        // therefore it should not be substracted.
		        //-------------------------------------------------------------------------------------------
		        /*G1_i += - indPi*indPj;
		        G2_i += - (indPi*rij)*(indPj*rij); */

                        Vector3 fi   = D1*gradG1_i + D2*gradG2_i + D3*gradG3_i + (G1_i*D2 + G2_i*D3 + G3_i*D4)*rij;
                        Vector3 fj   = D1*gradG1_j + D2*gradG2_j + D3*gradG3_j + (G1_j*D2 + G2_j*D3 + G3_j*D4)*rji;
                        Vector3 fij  = fi - fj;

    			molderiv    += -fij*rij*COUlOMBS/(3.0*volume);
    			mol2ndDeriv += 0.0;

                        atoms[idi].force += fij*COUlOMBS;
                        atoms[idj].force -= fij*COUlOMBS;

		        //-------------------------------------------------------------------------------------------
		        // Computation of Torques
		        //-------------------------------------------------------------------------------------------
	 
		        if (computeTorques)
		        {
			    double Qi[3][3], Qj[3][3], pjrij[3][3], rijpj[3][3], Rij[3][3];
			    double pirji[3][3], rjipi[3][3] ;
			    Vector3 QiXpjrij, QiXrijpj, QiXRij;
			    Vector3 QjXpirji, QjXrjipi, QjXRji;
			    Vector3 ones(1,1,1);

			    Qi[0][0] = Qxxi;
			    Qi[0][1] = Qxyi;
			    Qi[0][2] = Qxzi;
			    Qi[1][0] = Qyxi;
			    Qi[1][1] = Qyyi;
			    Qi[1][2] = Qyzi;
			    Qi[2][0] = Qzxi;
			    Qi[2][1] = Qzyi;
			    Qi[2][2] = Qzzi;

			    Qj[0][0] = Qxxj;
			    Qj[0][1] = Qxyj;
			    Qj[0][2] = Qxzj;
			    Qj[1][0] = Qyxj;
			    Qj[1][1] = Qyyj;
			    Qj[1][2] = Qyzj;
			    Qj[2][0] = Qzxj;
			    Qj[2][1] = Qzyj;
			    Qj[2][2] = Qzzj;

			    Rij[0][0] = rij.x*rij.x;
			    Rij[0][1] = rij.x*rij.y;
			    Rij[0][2] = rij.x*rij.z;
			    Rij[1][0] = rij.y*rij.x;
			    Rij[1][1] = rij.y*rij.y;
			    Rij[1][2] = rij.y*rij.z;
			    Rij[2][0] = rij.z*rij.x;
			    Rij[2][1] = rij.z*rij.y;
			    Rij[2][2] = rij.z*rij.z;

			    pjrij[0][0] = indPj.x*rij.x;
			    pjrij[0][1] = indPj.x*rij.y;
			    pjrij[0][2] = indPj.x*rij.z;
			    pjrij[1][0] = indPj.y*rij.x;
			    pjrij[1][1] = indPj.y*rij.y;
			    pjrij[1][2] = indPj.y*rij.z;
			    pjrij[2][0] = indPj.z*rij.x;
			    pjrij[2][1] = indPj.z*rij.y;
			    pjrij[2][2] = indPj.z*rij.z;

			    rijpj[0][0] = pjrij[0][0];
			    rijpj[0][1] = pjrij[1][0];
			    rijpj[0][2] = pjrij[2][0];
			    rijpj[1][0] = pjrij[0][1];
			    rijpj[1][1] = pjrij[1][1];
			    rijpj[1][2] = pjrij[2][1];
			    rijpj[2][0] = pjrij[0][2];
			    rijpj[2][1] = pjrij[1][2];
			    rijpj[2][2] = pjrij[2][2];

			    pirji[0][0] = indPi.x*rji.x;
			    pirji[0][1] = indPi.x*rji.y;
			    pirji[0][2] = indPi.x*rji.z;
			    pirji[1][0] = indPi.y*rji.x;
			    pirji[1][1] = indPi.y*rji.y;
			    pirji[1][2] = indPi.y*rji.z;
			    pirji[2][0] = indPi.z*rji.x;
			    pirji[2][1] = indPi.z*rji.y;
			    pirji[2][2] = indPi.z*rji.z;
		
			    rjipi[0][0] = pirji[0][0];
			    rjipi[0][1] = pirji[1][0];
			    rjipi[0][2] = pirji[2][0];
			    rjipi[1][0] = pirji[0][1];
			    rjipi[1][1] = pirji[1][1];
			    rjipi[1][2] = pirji[2][1];
			    rjipi[2][0] = pirji[0][2];
			    rjipi[2][1] = pirji[1][2];
			    rjipi[2][2] = pirji[2][2];

			    QiXpjrij = cross(Qi, pjrij, ones);
			    QiXrijpj = cross(Qi, rijpj, ones);
			    QiXRij   = cross(Qi, Rij, ones);

			    QjXpirji = cross(Qj, pirji, ones);
			    QjXrjipi = cross(Qj, rjipi, ones);
			    QjXRji   = cross(Qj, Rij, ones);    // Because Rji = Rij.

    			    Gt1i = cross(dipi, indPj);
			    Gt2i = (indPj*rij)*(cross(dipi, rij)) + 2.0*QiXpjrij + 2.0*QiXrijpj;
			    Gt3i = 2.0*(indPj*rij)*QiXRij;
			    
    			    Gt1j = cross(dipj, indPi);
			    Gt2j = (indPi*rji)*(cross(dipj, rji)) + 2.0*QjXpirji + 2.0*QjXrjipi;
			    Gt3j = 2.0*(indPi*rji)*QjXRji;

			    trq[idi] += -COUlOMBS*(Gt1i*D1 + Gt2i*D2 + Gt3i*D3); 
			    trq[idj] += -COUlOMBS*(Gt1j*D1 + Gt2j*D2 + Gt3j*D3); 

		        } // If-sentence Torques

                    }  // if-sentence not constraints
                }
            }          
        }
      }
    }
    intraMolEnergy = -0.5*intraMolecularEnergy*COUlOMBS;
    intraVircorrect = 0.0;
    //cout << U_i << endl;
    double temp;
//    MPI::COMM_WORLD.Reduce(&intraMolEnergy,&temp,1,MPI::DOUBLE,MPI::SUM,0);
    MPI_Reduce(&intraMolEnergy,&temp,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    intraMolEnergy = temp;

    temp = 0.0;
    MPI_Reduce(&intraVircorrect,&temp,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    intraVircorrect = temp;

    temp = 0.0;
    MPI_Allreduce(&molderiv,&temp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    molderiv = temp;

    temp = 0.0;
    MPI_Reduce(&mol2ndDeriv,&temp,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    mol2ndDeriv = temp;


    if (computeTorques){
        for (int i = 0; i < numAtoms; i++)
        {
            TorqueArray[i][0] = trq[i].x;
            TorqueArray[i][1] = trq[i].y;
            TorqueArray[i][2] = trq[i].z;

            tempArray[i][0] = 0.0;
            tempArray[i][1] = 0.0;
            tempArray[i][2] = 0.0;
        }

        MPI_Allreduce(&TorqueArray,&tempArray,nAtomsx3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);

        for (int i = 0; i < numAtoms; i++)
        {
            trq[i].x = tempArray[i][0];
            trq[i].y = tempArray[i][1];
            trq[i].z = tempArray[i][2];

	    totTrqI += trq[i].length();
        }

	//----------------------------------------------------------------
	// Torques into Forces
	//----------------------------------------------------------------

	double dUdV = 0.0;
        torques2forces(dUdV);

	molderiv += dUdV/(3.0*volume);

    }

}

//  surface dipole correction
void Induction::surfDipoleCorrect(Double &surfaceDipoleEnergy, double &surfderiv, double &surf2ndderiv)
{
    if(doUpdate) {
        volume = myEnsemble->volume;
        volume2 = volume*volume;
        volumer = 1.0/volume;
    }
    int molIDi, i3; 
    int idi, idj, ns;
    int atomID;
    double r, r2, qi;
    double dipxi, dipyi, dipzi;
    double  selfSum = 0.0;
    double  selfpi = 0.0;
    double dUdV1, dUdV2;
    double selfterm;
    double TorqueArray[numAtoms][3], tempArray[numAtoms][3];
    const double coeff = 2.0*PI*volumer/3.0;
    Vector3 dipi, indPi;
    Vector3 qrSum, qRsum, rim;
    Vector3 miSum, piSum;
    Vector3 ri, fi, temp;
    Molecule *moli;

    surfaceDipoleEnergy = 0.0;
    surfderiv           = 0.0;
    surf2ndderiv        = 0.0;

    //--------------------------------------------------------------------------
    //  First clean up the torques for all the processors
    //--------------------------------------------------------------------------

    for (int i = 0; i < numAtoms; i++)
    {
        trq[i].x = 0;
        trq[i].y = 0;
        trq[i].z = 0;
    }

    for (Int i = 0; i < numAtoms; i++)
    {
        Vector3 r = atoms[i].position;

	molIDi = atoms[i].molID;
        rim = myEnsemble->molecules[molIDi].massCenter;

        if(!strcmp(myMols[molIDi].molName, "MCY") && atoms[i].atomType ==1)   // Jc: added by Jianhui
        { 
            r.x = Dummy[molIDi*3];
            r.y = Dummy[molIDi*3 +1];
            r.z = Dummy[molIDi*3 +2]; 
        }
        qi = atoms[i].scaledCharge/SQRTCOULOMBCONSTANT;
        dipi = myL2G->dip[i];

	//-----------------------------------------------------------------------
	// i-index Induced Dipole	
	//-----------------------------------------------------------------------
	
        i3  = 3 * i;
	indPi.x =  Pi[i3]; 
	indPi.y =  Pi[i3+1]; 
	indPi.z =  Pi[i3+2]; 

        // myEnsemble->apply_pbc(r);  //Jc: why?
        qrSum += qi*r;
        qRsum += qi*rim;
        miSum += dipi;
        piSum += indPi;
        selfSum += dipi*indPi;
	selfpi  += indPi*indPi;

    }

    selfterm = qrSum*piSum + miSum*piSum - selfSum;

    surfaceDipoleEnergy = coeff*selfterm*COUlOMBS;

    if (useConstraint){
	if (doMinimization) {
            dUdV1 = -coeff*volumer*(2.0*qrSum*piSum + 2.0*miSum*piSum + piSum*piSum - 2.0*selfSum - selfpi) ;
	} else {
            dUdV1 = -coeff*volumer*(2.0*qrSum*piSum + 2.0*miSum*piSum - 2.0*selfSum) ;
	}
	dUdV2 = 2.0*coeff*volumer*qrSum*piSum/3.0 ;
        surfderiv      = COUlOMBS*(dUdV1 + dUdV2);
        surf2ndderiv   = 2.0*coeff*volumer*volumer*qrSum*qrSum - 16.0*coeff*volumer*volumer*qrSum*qRsum/9.0 + 2.0*coeff*volumer*volumer*qRsum*qRsum/9.0;
    }
    else {
	if (doMinimization) {
            dUdV1 = -coeff*volumer*(2.0*qrSum*piSum + 2.0*miSum*piSum + piSum*piSum - 2.0*selfSum - selfpi) ;
	} else {
            dUdV1 = -coeff*volumer*(2.0*qrSum*piSum + 2.0*miSum*piSum - 2.0*selfSum) ;
	}
	dUdV2 = 2.0*coeff*volumer*qrSum*piSum/3.0 ;
        surfderiv      = COUlOMBS*(dUdV1 + dUdV2);
        surf2ndderiv   = 4.0*coeff*volumer*volumer*qrSum*qrSum/9.0;
    }

// Jc: parallel code;
    int numDecompose = (int) numAtoms/psize;
    if ((numAtoms % psize) == 0) numDecompose--;
    for(int i = 0; i <= numDecompose;i++)     
    { 
      int atomIndex = i*psize + psize - 1- prank;
      if(atomIndex <numAtoms ) 
      {

	atomID = atoms[atomIndex].atomID ;
	molIDi = atoms[atomIndex].molID;

        if (atomIndex != atomID ) { 
     	        ERRORMSG("Surface term Correction: ID counting is probably wrong. Check your atom index in sysData.in. ");}

        qi = atoms[atomID].scaledCharge/SQRTCOULOMBCONSTANT;
        dipi = myL2G->dip[atomID];

	//-----------------------------------------------------------------------
	// i-index Induced Dipole	
	//-----------------------------------------------------------------------
	
        i3  = 3 * atomID;
	indPi.x =  Pi[i3]; 
	indPi.y =  Pi[i3+1]; 
	indPi.z =  Pi[i3+2]; 

        ri = atoms[atomID].position;
        if(!strcmp(myMols[molIDi].molName, "MCY") && atoms[atomID].atomType ==1)   // Oxygen location MCYna
        {
            ri.x = Dummy[molIDi*3];
            ri.y = Dummy[molIDi*3 +1];
            ri.z = Dummy[molIDi*3 +2];
        }

        fi = -2.0*coeff*qi*piSum;
        temp = -atoms[atomIndex].scaledCharge*qrSum;

	moli =  &myMols[atoms[atomID].molID];
        Atom *atomh1i = moli->myAtoms[0];
        Atom *atomh2i = moli->myAtoms[2];

        if (!strcmp(myMols[atoms[atomID].molID].molName, "MCY") && (atoms[atomID].atomType == 1) ){
                atoms[atomID].force += 0.5431676*fi*COUlOMBS;
                atomh1i->force += 0.5*0.4568324*fi*COUlOMBS;
                atomh2i->force += 0.5*0.4568324*fi*COUlOMBS;
        }
        else {
                atoms[atomID].force += fi*COUlOMBS;
        }
          //atoms[atomIndex].force -= 2.0*atoms[atomIndex].scaledCharge*qrSum;        // Gromacs. This may be right
                virial[XX] += temp.x * ri.x;
                virial[XY] += temp.y * ri.x;
                virial[XZ] += temp.z * ri.x;
                virial[YX] += temp.x * ri.y;
                virial[YY] += temp.y * ri.y;
                virial[YZ] += temp.z * ri.y;
                virial[ZX] += temp.x * ri.z;
                virial[ZY] += temp.y * ri.z;
                virial[ZZ] += temp.z * ri.z;

        //-------------------------------------------------------------------------------------------
        // Computation of Torques
        //-------------------------------------------------------------------------------------------

        if (computeTorques)
        {
	    Vector3 SumpjXmi = cross(piSum, dipi);
	    Vector3 piXmi    = cross(indPi, dipi);

	    trq[atomID] += COUlOMBS*2.0*coeff*(SumpjXmi - piXmi); 

        }
      }
    }
 

    if (computeTorques){
        for (int i = 0; i < numAtoms; i++)
        {
            TorqueArray[i][0] = trq[i].x;
            TorqueArray[i][1] = trq[i].y;
            TorqueArray[i][2] = trq[i].z;

            tempArray[i][0] = 0.0;
            tempArray[i][1] = 0.0;
            tempArray[i][2] = 0.0;
        }

        MPI_Allreduce(&TorqueArray,&tempArray,nAtomsx3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);

        for (int i = 0; i < numAtoms; i++)
        {
            trq[i].x = tempArray[i][0];
            trq[i].y = tempArray[i][1];
            trq[i].z = tempArray[i][2];

	    totTrqI += trq[i].length();
        }

	//----------------------------------------------------------------
	// Torques into Forces
	//----------------------------------------------------------------

	double dUdV = 0.0;
        torques2forces(dUdV);

	surfderiv += dUdV*volumer/3.0;
    }

}

void Induction:: ConjugateGradient(){

  int dim = numAtoms*3;
  const double tolerance = 1.0e-14;
  double c,t,d,dotV;
  double r[dim], v[dim], z[dim];

  c = 0.0; d = 0.0;t = 0.0,dotV=0.0;

  for(int i =0;i<dim;i++)
  {
     Pi[i]  = 1.0;
  }

  for(int i=0;i<dim;i++)
  {
    r[i] =  Ei0[i];
    for(int j = 0;j<dim;j++)
    {
      r[i] -= A[i][j] * Pi[j];
    }
    v[i] = r[i];
    c    += r[i]*r[i];
    dotV += v[i]*v[i];
  }

  //-----------------------------------------------------------------------
  // CG Iteration
  //-----------------------------------------------------------------------
  
  for(int i=0;i<dim;i++){
    if(sqrt(dotV)<tolerance){
      cerr << "An error has occurred in ConjugateGradient: execution of function terminated" << endl;
      break;
    }
    double dotVZ = 0.0;
    for(int j =0;j<dim;j++)
    {
       z[j] = 0.0;
       for(int k = 0; k<dim;k++)
       {
         z[j] += A[j][k]*v[k];
       }
       dotVZ += v[j]*z[j];
    }
    t = c/dotVZ;

    d = 0.0;
    for(int j =0;j<dim;j++)
    {
      Pi[j] += t*v[j];
      r[j]  -= t*z[j];
      d     += r[j]*r[j];
    }

    if(sqrt(d) < tolerance)
      break;

    dotV = 0.0;
    for(int j =0; j<dim;j++)
    {
      v[j] = r[j] + (d/c)*v[j];
      dotV += v[j]*v[j];
    }
    c = d;
  }
}

void Induction:: ParallelCG(){

  int dim = numAtoms*3;
  const double tolerance = 1.0e-14;
  double c,t,d,dotV, sum, local_sum;
  int offset;
  double v[dim];

  c = 0.0; d = 0.0;t = 0.0,dotV=0.0;

  for(int i =0;i<dim;i++)
  {
     Pi[i]  = 1.0;
     v[i]   = 1.0;
  }

  int Nrows_per_proc = (int) dim/psize;
  double r[Nrows_per_proc], x[Nrows_per_proc], z[Nrows_per_proc];

  offset = prank * Nrows_per_proc;

  for(int i=0;i<Nrows_per_proc;i++)
  {
    r[i] =  Ei0[offset + i];
    for(int j = 0;j<dim;j++)
    {
      r[i] -= A[offset + i][j] * Pi[j];
    }
    local_sum += r[i]*r[i];
    x[i] = Pi[offset + i];
  }

   sum = 0.0;
   MPI_Allreduce(&local_sum,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   c = sum;

   MPI_Allgather(r,Nrows_per_proc,MPI_DOUBLE,v,Nrows_per_proc,MPI_DOUBLE,MPI_COMM_WORLD);

  //-----------------------------------------------------------------------
  // CG Iteration
  //-----------------------------------------------------------------------
  
  for(int i=0;i<dim;i++){

    double dotVZ = 0.0;
    for(int j =0;j<Nrows_per_proc;j++)
    {
       z[j] = 0.0;
       for(int k = 0; k<dim;k++)
       {
         z[j] += A[offset + j][k]*v[k];
       }
       dotVZ += v[offset + j]*z[j];
    }

    sum = 0.0;
    MPI_Allreduce(&dotVZ,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    t = c/sum;

    d = 0.0, local_sum = 0.0;
    for(int j =0;j<Nrows_per_proc;j++)
    {
      x[j]      += t*v[offset + j];
      r[j]      -= t*z[j];
      local_sum += r[j]*r[j];
    }

    sum = 0.0;
    MPI_Allreduce(&local_sum,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    d = sum;

    if(sqrt(d) < tolerance)
      break;

    for(int j =0; j<Nrows_per_proc;j++)
    {
      z[j] = r[j] + (d/c)*v[offset + j];
    }

    MPI_Allgather(z,Nrows_per_proc,MPI_DOUBLE,v,Nrows_per_proc,MPI_DOUBLE,MPI_COMM_WORLD);

    c = d;
  }

  MPI_Allgather(x,Nrows_per_proc,MPI_DOUBLE,Pi,Nrows_per_proc,MPI_DOUBLE,MPI_COMM_WORLD);

}

void Induction::torques2forces(double &dUdVol)
{
    //----------------------------------------------------------------------
    // Variables for this method
    //----------------------------------------------------------------------

    int ns, molID, id;
    double val;
    Molecule *mol;
    Vector3 dUdu, dUdv;
    Vector3 vxu, wxu, wxv;
    Vector3 u, v, w, tau;
    double ul, vl, wl;
    double vxuL, wxuL, wxvL;
    double wrot_f = 0.5;

    dUdVol = 0.0;

    //----------------------------------------------------------------------
    // Parallel decomposition as number of molecules
    //----------------------------------------------------------------------

    int numDecompose = (int) numMols/psize;
    if ((numMols % psize) == 0) numDecompose--;
    for(int i = 0; i <= numDecompose;i++)
    {
      int molIndex = i*psize + prank; // Jc:i*psize + psize - 1- prank;
      if(molIndex < numMols)
      {
        ns    = myEnsemble->molecules[molIndex].numAtoms;
        molID = myEnsemble->molecules[molIndex].molID;

        mol = &myMols[molID];

        Atom *atomH1 = mol->myAtoms[0];
        Atom *atomO  = mol->myAtoms[1];
        Atom *atomH2 = mol->myAtoms[2];

        //-------------------------------------------------------------------
        //  Transformation for H1. The distances are ri - rj 
        //  in order to mantain consistency of the code.
        //-------------------------------------------------------------------

        id  = atomH1->atomID;
        tau = trq[id];

        //:::::::::::::::: local vectors ::::::::::::::::::::::::::::::::::
        
	u   = atomH1->position - atomO->position;
        ul  = u.length();

        v.x = 1.0;
        v.y = 0.0;
        v.z = 0.0;
        val = u.x/ul;
        if (fabs(val) > 0.8666){
            v.x = 0.0;
            v.y = 1.0;
        }
        vl = 1.0;

        w  = cross(u, v);
        wl = w.length();

        //::::::::::::::::: unit vectors ::::::::::::::::::::::::::::::::::

        u = u/ul;
        v = v/vl;
        w = w/wl;

        //::: infinitesimal displacement of u due to v and w ::::::::::::::
        
	vxu  = cross(v, u);
        vxuL = vxu.length();

        wxu  = cross(w, u);
        wxuL = wxu.length();

        dUdu = (-v*tau)*vxu/(vxuL*vxuL*ul) + (-w*tau)*wxu/(wxuL*wxuL*ul);

        atomH1->force += -dUdu;
        atomO->force  +=  dUdu;

        dUdVol += dUdu*u*ul;       // 1/3V factor is applied outside.

        //-------------------------------------------------------------------
        //  Transformation for O. The distances are ri - rj 
        //  in order to mantain consistency of the code.
        //-------------------------------------------------------------------

	id  = atomO->atomID;
        tau = trq[id];

        //:::::::::::::::: local vectors ::::::::::::::::::::::::::::::::::

        u  = atomO->position - atomH1->position;
        v  = atomO->position - atomH2->position;
        w  = cross(u, v);

        ul = u.length();
        vl = v.length();
        wl = w.length();

        //::::::::::::::::: unit vectors ::::::::::::::::::::::::::::::::::
        
	u = u/ul;
        v = v/vl;
        w = w/wl;

        //::: infinitesimal displacement of u due to v and w ::::::::::::::

	vxu  = cross(v, u);
        vxuL = vxu.length();

        wxu  = cross(w, u);
        wxuL = wxu.length();

        dUdu = (-v*tau)*vxu/(vxuL*vxuL*ul) + wrot_f*(-w*tau)*wxu/(wxuL*wxuL*ul);

        //::: infinitesimal displacement of v due to u and w ::::::::::::::

	wxv  = cross(w, v);
        wxvL = wxv.length();

        dUdv = (-u*tau)*(-vxu)/(vxuL*vxuL*vl) + wrot_f*(-w*tau)*wxv/(wxvL*wxvL*vl);

        atomO->force  += - dUdu - dUdv;
        atomH1->force += dUdu;
        atomH2->force += dUdv;

        dUdVol += dUdu*u*ul + dUdv*v*vl;       // 1/3V factor is applied outside.

        //-------------------------------------------------------------------
        //  Transformation for H2. The distances are ri - rj 
        //  in order to mantain consistency of the code.
        //-------------------------------------------------------------------

        id  = atomH2->atomID;
        tau = trq[id];

        //:::::::::::::::: local vectors ::::::::::::::::::::::::::::::::::
        
	u   = atomH2->position - atomO->position;
        ul  = u.length();

        v.x = 1.0;
        v.y = 0.0;
        v.z = 0.0;
        val = u.x/ul;
        if (fabs(val) > 0.8666){
            v.x = 0.0;
            v.y = 1.0;
        }
        vl = 1.0;

        w  = cross(u, v);
        wl = w.length();

        //::::::::::::::::: unit vectors ::::::::::::::::::::::::::::::::::

        u = u/ul;
        v = v/vl;
        w = w/wl;

        //::: infinitesimal displacement of u due to v and w ::::::::::::::
        
	vxu  = cross(v, u);
        vxuL = vxu.length();

        wxu  = cross(w, u);
        wxuL = wxu.length();

        dUdu = (-v*tau)*vxu/(vxuL*vxuL*ul) + (-w*tau)*wxu/(wxuL*wxuL*ul);

        atomH2->force += -dUdu;
        atomO->force  +=  dUdu;

        dUdVol += dUdu*u*ul;       // 1/3V factor is applied outside.

      }
    }

   double temp = 0.0;
   MPI_Allreduce(&dUdVol,&temp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   dUdVol = temp;

}


void Induction::write_force_info(ofstream& of)
{   
    of << "Induction Force created" << endl;
    of<<"alpha: "<<alpha<<'\t'<<"real cut off: "<<realCut<<'\t'<<"K space cut off: "<<kCut<<endl;
    of << endl;
    of << 'I' << '\t' << "kModulu[i]" << '\t' << "expK[i]" << endl;
    for (int i = 0; i < kModulu.size(); i++)
        of << i << '\t' << kModulu[i] << '\t' << expK[i] << endl;
    of << endl;

    #ifdef USE_LOOK_UP_TABLE 
    /**
    of << "r2" << '\t' << "potentialTable[i]" << '\t' << "forceTable[i]" << endl; 
    for (int i = 0; i < tableSize; i++)
        of << i/stepr << '\t' << potentialTable[i] << '\t' << forceTable[i] << endl;  
    of << endl;
    **/
    #endif
    
    of << "Self correction energy: " << selfEnergy << '\t' << "Charged energy: " << chargedEnergy << endl;
    of << endl;
    
}

void Induction::write_energy(ofstream& of)
{
    of << "Induction Energy: " << energy+longEnergy << endl;
    of << "realEnergy: " << realEnergy << '\t' << "longEnergy: " << longEnergy << endl ;
    // of << "introMolEnergy: " << introMolEnergy << '\t' << "surfaceEnergy" << surfDipoleEnergy << endl;
}


