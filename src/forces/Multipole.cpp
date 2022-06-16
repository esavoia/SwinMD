/** Multipole.cpp -- 
 **
 ** Copyright (C) 2003
 ** Centre for Molecular Simulation (CMS)
 ** School of Information Technology
 ** Swinburne University of Technology
 ** PO Box 21 Hawthorn, Vic 3122, Australia
 **
 ** Author: Zhongwu Zhou
 ** Email: zzhou@it.swin.edu.au
 ** Ewald Method for Multipole systems by
 ** Author: eoyarzua
 ** Email: eoyarzua@swin.edu.au
 **/

#include "Multipole.h" 

///#include "erf.h"
#include <mpi.h>

Multipole::Multipole(Ensemble* ensemble) : Force(ensemble) 
{
    DEBUGMSG("Creating Multipole");

//JC    forceIdentifier = "EwaldForce"; // in order to comply with the cluster compiler
    numAtoms = myEnsemble->nAtoms;
    numMols = myEnsemble->nMols;
    myMols = myEnsemble->molecules;
    volume = myEnsemble->volume;
    volume2 = volume*volume;
    volumer = 1.0/volume;
    alpha = myEnsemble->myConfig->get_alpha();
    alpha2 = alpha*alpha;
    alpha2PI = 2.0*alpha/SQRT_PI;
    alphaR4 = 1.0/(4.0*alpha2);                 
    realCut = myEnsemble->myConfig->get_cutoffEw();
    realCut2 = realCut*realCut;
    kCut = myEnsemble->myConfig->get_k_cutoff(); 
    correctSurfDipole = myEnsemble->myConfig->compute_surf_correct(); 
    useConstraint = myEnsemble->myConfig->is_constraint_on();
    statEnsem = myEnsemble->myConfig->get_ensemble_status();  
    computeTorques = myEnsemble->myConfig->compute_torques();

    if (computeTorques) { DEBUGMSG("Compute the torques due to dipoles and translate them into atomic forces");}
        else { DEBUGMSG("Torques are NOT included in the force computation"); }

    nAtomsx3 = 3*numAtoms;
    stepr = 0.0;
    chargedEnergy = 0.0;
    longEnergy = 0.0;
    longderiv  = 0.0;
    long2ndDeriv  = 0.0;
    tableSize = TAB_SIZE + 2;

    myL2G  = NULL;

    trq  = NULL;
    kXmsina = NULL; 
    kXmcosa = NULL;
    KXQcosa = NULL;
    KXQsina = NULL;
    expK = NULL;
    qQsin_i = NULL;
    qQcos_i = NULL;
    mksin_i = NULL;
    mkcos_i = NULL;
    forceTable = NULL;
    potentialTable = NULL;
    doUpdate = false;
    if (statEnsem == 1){
        doUpdate = true;
        DEBUGMSG(" Volume Update for NPT ensemble on Multipole Electrostatic");
    }

    myL2G  = new Local2Global(myEnsemble);

    qQsin_i = new double[numAtoms];
    qQcos_i = new double[numAtoms];
    mksin_i = new double[numAtoms];
    mkcos_i = new double[numAtoms];
    trq     = new Vector3 [numAtoms];
    kXmsina = new Vector3 [numAtoms];
    kXmcosa = new Vector3 [numAtoms];
    KXQcosa = new Vector3 [numAtoms];
    KXQsina = new Vector3 [numAtoms];
    Dummy = new Double[numMols*3];
    if ((qQsin_i==NULL)||(qQcos_i==NULL)||(mksin_i==NULL)||(mkcos_i==NULL))
        ERRORMSG("memory allocation error for Multipole Ewald");
    if ((trq==NULL)||(kXmsina==NULL)||(kXmcosa==NULL)||(KXQcosa==NULL)||(KXQsina==NULL))
        ERRORMSG("memory allocation error for Multipole Ewald (2)");

    set_up();
    DEBUGMSG("Multipole set up");

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


Multipole::~Multipole()
{
    if(trq != NULL) delete [] trq;
    if(expK != NULL) delete [] expK;
    if(KXQcosa != NULL) delete [] KXQcosa;
    if(KXQsina != NULL) delete [] KXQsina;
    if(kXmsina != NULL) delete [] kXmsina;
    if(kXmcosa != NULL) delete [] kXmcosa;
    if(qQsin_i != NULL) delete [] qQsin_i;
    if(qQcos_i != NULL) delete [] qQcos_i;
    if(mksin_i != NULL) delete [] mksin_i;
    if(mkcos_i != NULL) delete [] mkcos_i;
    if(forceTable != NULL) delete [] forceTable;
    if(potentialTable != NULL) delete [] potentialTable;    
}


void Multipole::DummyPosition()

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

void Multipole::compute()
{
    energy = 0.0;
    selfEnergy = 0.0;
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
    totTrq = 0.0;

    DummyPosition();

    myL2G->FrameofReference();

    for (Int i = XX; i <= ZZ; i++)
        virial[i] = 0.0;

    if(doUpdate)
        set_up();

    self_correct(selfEnergy);

    real_term(realEnergy, realderiv, real2ndDeriv);                     // Jc: the first term in Ewald

    intraMolCorrect(intraMolEnergy, intraVircorrect, molderiv, mol2ndDeriv);           // Molec correction in Ewald

    if(correctSurfDipole){
        surfDipoleCorrect(surfDipoleEnergy, surfderiv, surf2ndderiv);   // Jc: the first term in Ewald
    }

    myEnsemble->myEnergy.realEnergy = realEnergy;
    myEnsemble->myEnergy.correctEnergy = selfEnergy+chargedEnergy; //Jc: selfEnergy and charged Energy need no reduce
    myEnsemble->myEnergy.molCorrectEnergy = intraMolEnergy;
    myEnsemble->myEnergy.surfCorrectEnergy = surfDipoleEnergy;

    myEnsemble->myLustig.intravirCorrect = intraVircorrect;
    myEnsemble->myLustig.realDeriv = realderiv;
    myEnsemble->myLustig.realSecDeriv = real2ndDeriv;

    myEnsemble->myLustig.molDeriv = molderiv;
    myEnsemble->myLustig.molSecDeriv = mol2ndDeriv;

    myEnsemble->myLustig.surfDeriv = surfderiv;
    myEnsemble->myLustig.surfSecDeriv = surf2ndderiv;

    for (Int k = XX; k <= ZZ; k++)
        myEnsemble->virial[k] -= virial[k];
}     
 
void Multipole::compute_long()
{
    DummyPosition();    

    myL2G->FrameofReference();

    for (Int i = XX; i <= ZZ; i++)
        longVirial[i] = 0.0;

    for (Int i = 0; i < numAtoms; i++)
        atoms[i].longForce = 0.0;

    reciprocal_term(longEnergy, longderiv, long2ndDeriv);
    myEnsemble->myEnergy.longEnergy = longEnergy;
    
    myEnsemble->myLustig.longDeriv = longderiv;
    myEnsemble->myLustig.longSecDeriv = long2ndDeriv;

    myEnsemble->myEnergy.TotTorq = totTrq;
    
    for(Int i = 0; i < ZZ; i++)
        myEnsemble->longVirial[i] = longVirial[i];
}

void Multipole::set_up()
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
void Multipole::build_look_up_table()
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

void Multipole::init_table(Int size, Double** table)
{
    if (*table != NULL)
        delete [] *table;
    *table = new Double[size];
    if(*table == NULL) 
        ERRORMSG("Multipole::init_table(), not enough memory");
}

void Multipole::self_correct(double &selfEnergy)
{

    #ifdef DEBUG2
        DEBUGMSG("Computing Multipole Self-Interaction Correction term");
    #endif

    int molIDi;
    double qi, dipxi, dipyi, dipzi;
    double Qxxi, Qxyi, Qxzi, Qyyi, Qyzi, Qzzi; 
    double Qyxi, Qzxi, Qzyi; 
    double QiI, QiQi, selfterm; 
    Vector3 dipi;
    
    selfEnergy = 0.0;

    for (Int i = 0; i < numAtoms; i++)
        {
	    molIDi = atoms[i].molID;
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

	//-------------------------------------------------------------------------------
	// Let's compute the scalar product between matrices.
	//-------------------------------------------------------------------------------
	    QiI   = Qxxi + Qyyi + Qzzi;

	    QiQi  = Qxxi*Qxxi + Qxyi*Qxyi + Qxzi*Qxzi \
		  + Qyxi*Qyxi + Qyyi*Qyyi + Qyzi*Qyzi \
		  + Qzxi*Qzxi + Qzyi*Qzyi + Qzzi*Qzzi ;


            selfterm = qi*qi*alpha2PI + (dipi*dipi - 2*qi*QiI)*2*alpha2PI*alpha2/3 + (QiI*QiI + 2*QiQi)*4*alpha2PI*alpha2*alpha2/5;
            selfEnergy += -0.5*selfterm*COUlOMBS;

	}

    //---------------------------------------------------------------------------
    // NOTE: there is no self-correction to forces or Torques in
    // the permanent multipole contribution.
    //---------------------------------------------------------------------------

}

//Jc:  compute real term in Ewald summation,it is necessary as well as the reciprocal term 
void Multipole::real_term(Double &realEnergy, double &realderiv, double &real2ndDeriv)
{

    #ifdef DEBUG2
        DEBUGMSG("Computing Multipole real space");
    #endif

    Atom *atomi, *atomj;
    Molecule *moli, *molj;
    Vector3 ri, rj, fij, rij, rijm, rim, rjm;
    int size;
    int molIDi, molIDj;
    double qi, qj;
    double dipxi, dipyi, dipzi;
    double dipxj, dipyj, dipzj;
    double Qxxi, Qxyi, Qxzi, Qyyi, Qyzi, Qzzi; 
    double Qyxi, Qzxi, Qzyi; 
    double Qxxj, Qxyj, Qxzj, Qyyj, Qyzj, Qzzj; 
    double Qyxj, Qzxj, Qzyj; 
    double G0, G1, G2, G3, G4, G5;
    double B0, B1, B2, B3, B4, B5;
    double r, r2, rijRij;
    double QiI, QjI, QiRij, QjRij, QiQj;
    double QimjXrij, QjmiXrij, QiXrijXQj;
    double TorqueArray[numAtoms][3], tempArray[numAtoms][3];
    Vector3 dipi, dipj;
    Vector3 gradG1, gradG2, gradG3, gradG4;
    Vector3 Gt1i, Gt2i, Gt3i, Gt4i;
    Vector3 Gt1j, Gt2j, Gt3j, Gt4j;
    Vector3 QiXrij, QjXrij, QiXmj, QjXmi;
    Vector3 QiXQjrij, QjXQirij;
    
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


        for (int j = 0; j < size; j++)
        {

            atomj = atomi->myPairList[j];
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


            myEnsemble->apply_pbc(rij);
            myEnsemble->apply_pbc(rijm);

            r2 = rij.length2();


            if (r2 < realCut2)
            {
		//----------------------------------------------------------------------------
		//  Constant needed for the G terms
		//----------------------------------------------------------------------------

		QiI   = Qxxi + Qyyi + Qzzi;
		QjI   = Qxxj + Qyyj + Qzzj;
		QiRij = Qxxi*rij.x*rij.x + Qxyi*rij.x*rij.y + Qxzi*rij.x*rij.z \
		      + Qyxi*rij.y*rij.x + Qyyi*rij.y*rij.y + Qyzi*rij.y*rij.z \ 
		      + Qzxi*rij.z*rij.x + Qzyi*rij.z*rij.y + Qzzi*rij.z*rij.z ; 

		QjRij = Qxxj*rij.x*rij.x + Qxyj*rij.x*rij.y + Qxzj*rij.x*rij.z \
		      + Qyxj*rij.y*rij.x + Qyyj*rij.y*rij.y + Qyzj*rij.y*rij.z \
		      + Qzxj*rij.z*rij.x + Qzyj*rij.z*rij.y + Qzzj*rij.z*rij.z ;

		QiQj  = Qxxi*Qxxj + Qxyi*Qxyj + Qxzi*Qxzj\
		      + Qyxi*Qyxj + Qyyi*Qyyj + Qyzi*Qyzj\
		      + Qzxi*Qzxj + Qzyi*Qzyj + Qzzi*Qzzj;

		QimjXrij = dipj.x*(Qxxi*rij.x + Qxyi*rij.y + Qxzi*rij.z) \ 
			 + dipj.y*(Qyxi*rij.x + Qyyi*rij.y + Qyzi*rij.z) \ 
			 + dipj.z*(Qzxi*rij.x + Qzyi*rij.y + Qzzi*rij.z) ;

		QjmiXrij = dipi.x*(Qxxj*rij.x + Qxyj*rij.y + Qxzj*rij.z) \ 
			 + dipi.y*(Qyxj*rij.x + Qyyj*rij.y + Qyzj*rij.z) \ 
			 + dipi.z*(Qzxj*rij.x + Qzyj*rij.y + Qzzj*rij.z) ;

		QiXrijXQj = rij.x*Qxxj*(Qxxi*rij.x + Qxyi*rij.y + Qxzi*rij.z) + rij.y*Qxyj*(Qxxi*rij.x + Qxyi*rij.y + Qxzi*rij.z) + rij.z*Qxzj*(Qxxi*rij.x + Qxyi*rij.y + Qxzi*rij.z) \ 
			  + rij.x*Qyxj*(Qyxi*rij.x + Qyyi*rij.y + Qyzi*rij.z) + rij.y*Qyyj*(Qyxi*rij.x + Qyyi*rij.y + Qyzi*rij.z) + rij.z*Qyzj*(Qyxi*rij.x + Qyyi*rij.y + Qyzi*rij.z) \ 
			  + rij.x*Qzxj*(Qzxi*rij.x + Qzyi*rij.y + Qzzi*rij.z) + rij.y*Qzyj*(Qzxi*rij.x + Qzyi*rij.y + Qzzi*rij.z) + rij.z*Qzzj*(Qzxi*rij.x + Qzyi*rij.y + Qzzi*rij.z) ;
		G0 = qi*qj;
		G1 = qj*dipi*rij - qi*dipj*rij - dipi*dipj + qj*QiI + qi*QjI;
		G2 = qj*QiRij + qi*QjRij - (dipi*rij)*(dipj*rij) - dipj*rij*QiI + dipi*rij*QjI -2*QimjXrij + 2*QjmiXrij + QiI*QjI + 2*QiQj; 
		G3 = dipi*rij*QjRij - dipj*rij*QiRij + QiI*QjRij + QjI*QiRij + 4*QiXrijXQj;
		G4 = QiRij*QjRij;

		//----------------------------------------------------------------------------
		//  Vectors needed for the gradG terms
		//----------------------------------------------------------------------------
		

		QiXrij.x = Qxxi*rij.x + Qxyi*rij.y + Qxzi*rij.z; 
		QiXrij.y = Qyxi*rij.x + Qyyi*rij.y + Qyzi*rij.z; 
		QiXrij.z = Qzxi*rij.x + Qzyi*rij.y + Qzzi*rij.z; 

		QjXrij.x = Qxxj*rij.x + Qxyj*rij.y + Qxzj*rij.z;
		QjXrij.y = Qyxj*rij.x + Qyyj*rij.y + Qyzj*rij.z;
		QjXrij.z = Qzxj*rij.x + Qzyj*rij.y + Qzzj*rij.z;

		QiXmj.x  = Qxxi*dipj.x + Qxyi*dipj.y + Qxzi*dipj.z; 
		QiXmj.y  = Qyxi*dipj.x + Qyyi*dipj.y + Qyzi*dipj.z; 
		QiXmj.z  = Qzxi*dipj.x + Qzyi*dipj.y + Qzzi*dipj.z; 

		QjXmi.x  = Qxxj*dipi.x + Qxyj*dipi.y + Qxzj*dipi.z;
		QjXmi.y  = Qyxj*dipi.x + Qyyj*dipi.y + Qyzj*dipi.z;
		QjXmi.z  = Qzxj*dipi.x + Qzyj*dipi.y + Qzzj*dipi.z;

		QiXQjrij.x  = (Qxxi*Qxxj + Qxyi*Qyxj + Qxzi*Qzxj)*rij.x + (Qxxi*Qxyj + Qxyi*Qyyj + Qxzi*Qzyj)*rij.y + (Qxxi*Qxzj + Qxyi*Qyzj + Qxzi*Qzzj)*rij.z;
		QiXQjrij.y  = (Qyxi*Qxxj + Qyyi*Qyxj + Qyzi*Qzxj)*rij.x + (Qyxi*Qxyj + Qyyi*Qyyj + Qyzi*Qzyj)*rij.y + (Qyxi*Qxzj + Qyyi*Qyzj + Qyzi*Qzzj)*rij.z;
		QiXQjrij.z  = (Qzxi*Qxxj + Qzyi*Qyxj + Qzzi*Qzxj)*rij.x + (Qzxi*Qxyj + Qzyi*Qyyj + Qzzi*Qzyj)*rij.y + (Qzxi*Qxzj + Qzyi*Qyzj + Qzzi*Qzzj)*rij.z;

		QjXQirij.x  = (Qxxj*Qxxi + Qxyj*Qyxi + Qxzj*Qzxi)*rij.x + (Qxxj*Qxyi + Qxyj*Qyyi + Qxzj*Qzyi)*rij.y + (Qxxj*Qxzi + Qxyj*Qyzi + Qxzj*Qzzi)*rij.z;
		QjXQirij.y  = (Qyxj*Qxxi + Qyyj*Qyxi + Qyzj*Qzxi)*rij.x + (Qyxj*Qxyi + Qyyj*Qyyi + Qyzj*Qzyi)*rij.y + (Qyxj*Qxzi + Qyyj*Qyzi + Qyzj*Qzzi)*rij.z;
		QjXQirij.z  = (Qzxj*Qxxi + Qzyj*Qyxi + Qzzj*Qzxi)*rij.x + (Qzxj*Qxyi + Qzyj*Qyyi + Qzzj*Qzyi)*rij.y + (Qzxj*Qxzi + Qzyj*Qyzi + Qzzj*Qzzi)*rij.z;


		gradG1 = qj*dipi - qi*dipj;
		gradG2 = 2*qj*QiXrij - (dipj*rij)*dipi - (dipi*rij)*dipj - 2*QiXmj - QiI*dipj + 2*qi*QjXrij + 2*QjXmi + QjI*dipi; 
		gradG3 = -QiRij*dipj - (2*dipj*rij)*QiXrij + QjRij*dipi + (2*dipi*rij)*QjXrij + (2*QiI)*QjXrij + (2*QjI)*QiXrij + 4*QiXQjrij + 4*QjXQirij;
		gradG4 = (2*QjRij)*QiXrij + (2*QiRij)*QjXrij;

		    if (useConstraint){
			rijRij = rij*rijm;
			//cout << "ENCENDIDO" << endl;
		    }
		    else {
			rijRij = r2;
			//cout << "APAGADO" << endl;
		    }
                    r = sqrt(r2);

		    B0 = erfc(alpha*r)/r;
		    B1 = (-B0 - alpha2PI*exp(-alpha2*r2))/r2;
		    B2 = (-3*B1 + 2*alpha2PI*alpha2*exp(-alpha2*r2) )/r2;
		    B3 = (-5*B2 - 4*alpha2PI*alpha2*alpha2*exp(-alpha2*r2) )/r2;
		    B4 = (-7*B3 + 8*alpha2PI*alpha2*alpha2*alpha2*exp(-alpha2*r2) )/r2;
		    B5 = (-9*B4 - 16*alpha2PI*alpha2*alpha2*alpha2*alpha2*exp(-alpha2*r2) )/r2;


                    /* double e = erfc(alpha*r)/r;
                       double f = (e + alpha2PI*exp(-alpha2*r2))/r2;
		    if (!strcmp(myMols[molIDi].molName, "BBV")  && !strcmp(myMols[molIDj].molName, "BBV")) {
			double expBBV = 1.0 + exp(-2.0*(r/Aa0 - 2.0));
			double invexp = 1.0/expBBV;
			double Dij2   = invexp*invexp;
			double Dij8   = Dij2*Dij2*Dij2*Dij2;
			double Dij15  = Dij8*Dij2*Dij2*Dij2*invexp;
			   e = e*Dij15;
			   f = f*Dij15 - (30.0/Aa0)*exp(-2.0*(r/Aa0 - 2.0))*e*invexp/r;
		    }*/ 

		    //---------------------------------------------------------------------------
		    //  Total Real-Space Energy as sum_i sum_j>i sum_l G_lij * Blij
		    //---------------------------------------------------------------------------
		    
                    realEnergy  += (G0*B0 + G1*B1 + G2*B2 + G3*B3 + G4*B4)*COUlOMBS;
		    real2ndDeriv += 0.0;

		    //---------------------------------------------------------------------------------------------------
		    //  Total Real-Space Force as sum_j=/=i sum_l gradG_lij * Blij + G_lij * Bl+1 * rij
		    //---------------------------------------------------------------------------------------------------

		    fij = -(B1*gradG1 + B2*gradG2 + B3*gradG3 + B4*gradG4) - (G0*B1 + G1*B2 + G2*B3 + G3*B4 + G4*B5)*rij;
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

		//----------------------------------------------------------------------------
		//  Computation of Torques
		//----------------------------------------------------------------------------
		
		if (computeTorques){

		    double Qi[3][3], Qj[3][3], Rij[3][3];
		    double mjrij[3][3], rijmj[3][3], rijQjrij[3][3], Qjrijrij[3][3];
		    double mirji[3][3], rjimi[3][3], rjiQirji[3][3], Qirjirji[3][3];

    		    Vector3 ones(1,1,1);
		    Vector3 rijXmi, mjXmi, RijXQi, QjXQi, mjrijXQi, rijmjXQi;
		    Vector3 QjrijXmi, rijQjrijXQi, QjrijrijXQi;
		    Vector3 rji, rjiXmj, RjiXQj, mirjiXQj, rjimiXQj;
		    Vector3 QirijXmj, rjiQirjiXQj, QirjirjiXQj;

		    rji = -rij;
		
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

		    mjrij[0][0] = dipj.x*rij.x; 
		    mjrij[0][1] = dipj.x*rij.y; 
		    mjrij[0][2] = dipj.x*rij.z; 
		    mjrij[1][0] = dipj.y*rij.x; 
		    mjrij[1][1] = dipj.y*rij.y; 
		    mjrij[1][2] = dipj.y*rij.z; 
		    mjrij[2][0] = dipj.z*rij.x; 
		    mjrij[2][1] = dipj.z*rij.y; 
		    mjrij[2][2] = dipj.z*rij.z; 

		    rijmj[0][0] = mjrij[0][0];
		    rijmj[0][1] = mjrij[1][0];
		    rijmj[0][2] = mjrij[2][0];
		    rijmj[1][0] = mjrij[0][1];
		    rijmj[1][1] = mjrij[1][1];
		    rijmj[1][2] = mjrij[2][1];
		    rijmj[2][0] = mjrij[0][2];
		    rijmj[2][1] = mjrij[1][2];
		    rijmj[2][2] = mjrij[2][2];

		    Qjrijrij[0][0] = QjXrij.x*rij.x;
		    Qjrijrij[0][1] = QjXrij.x*rij.y;
		    Qjrijrij[0][2] = QjXrij.x*rij.z;
		    Qjrijrij[1][0] = QjXrij.y*rij.x;
		    Qjrijrij[1][1] = QjXrij.y*rij.y;
		    Qjrijrij[1][2] = QjXrij.y*rij.z;
		    Qjrijrij[2][0] = QjXrij.z*rij.x;
		    Qjrijrij[2][1] = QjXrij.z*rij.y;
		    Qjrijrij[2][2] = QjXrij.z*rij.z;

		    rijQjrij[0][0] = Qjrijrij[0][0];
		    rijQjrij[0][1] = Qjrijrij[1][0];
		    rijQjrij[0][2] = Qjrijrij[2][0];
		    rijQjrij[1][0] = Qjrijrij[0][1];
		    rijQjrij[1][1] = Qjrijrij[1][1];
		    rijQjrij[1][2] = Qjrijrij[2][1];
		    rijQjrij[2][0] = Qjrijrij[0][2];
		    rijQjrij[2][1] = Qjrijrij[1][2];
		    rijQjrij[2][2] = Qjrijrij[2][2];


		    mirji[0][0] = dipi.x*rji.x; 
		    mirji[0][1] = dipi.x*rji.y; 
		    mirji[0][2] = dipi.x*rji.z; 
		    mirji[1][0] = dipi.y*rji.x; 
		    mirji[1][1] = dipi.y*rji.y; 
		    mirji[1][2] = dipi.y*rji.z; 
		    mirji[2][0] = dipi.z*rji.x; 
		    mirji[2][1] = dipi.z*rji.y; 
		    mirji[2][2] = dipi.z*rji.z; 

		    rjimi[0][0] = mirji[0][0]; 
		    rjimi[0][1] = mirji[1][0]; 
		    rjimi[0][2] = mirji[2][0]; 
		    rjimi[1][0] = mirji[0][1]; 
		    rjimi[1][1] = mirji[1][1]; 
		    rjimi[1][2] = mirji[2][1]; 
		    rjimi[2][0] = mirji[0][2]; 
		    rjimi[2][1] = mirji[1][2]; 
		    rjimi[2][2] = mirji[2][2]; 

		    Qirjirji[0][0] = QiXrij.x*rij.x;
		    Qirjirji[0][1] = QiXrij.x*rij.y;
		    Qirjirji[0][2] = QiXrij.x*rij.z;
		    Qirjirji[1][0] = QiXrij.y*rij.x;
		    Qirjirji[1][1] = QiXrij.y*rij.y;
		    Qirjirji[1][2] = QiXrij.y*rij.z;
		    Qirjirji[2][0] = QiXrij.z*rij.x;
		    Qirjirji[2][1] = QiXrij.z*rij.y;
		    Qirjirji[2][2] = QiXrij.z*rij.z;

		    rjiQirji[0][0] = Qirjirji[0][0];
		    rjiQirji[0][1] = Qirjirji[1][0];
		    rjiQirji[0][2] = Qirjirji[2][0];
		    rjiQirji[1][0] = Qirjirji[0][1];
		    rjiQirji[1][1] = Qirjirji[1][1];
		    rjiQirji[1][2] = Qirjirji[2][1];
		    rjiQirji[2][0] = Qirjirji[0][2];
		    rjiQirji[2][1] = Qirjirji[1][2];
		    rjiQirji[2][2] = Qirjirji[2][2];

		    rijXmi      = cross(rij,dipi);
		    mjXmi       = cross(dipj,dipi);
		    QjrijXmi    = cross(QjXrij, dipi);
		    RijXQi      = cross(Rij, Qi, ones);
		    mjrijXQi    = cross(mjrij, Qi, ones);
		    rijmjXQi    = cross(rijmj, Qi, ones); 
		    QjXQi       = cross(Qj, Qi, ones);
		    rijQjrijXQi = cross(rijQjrij, Qi, ones);
		    QjrijrijXQi = cross(Qjrijrij, Qi, ones);

		    rjiXmj      = cross(rji, dipj);
		    QirijXmj    = cross(QiXrij, dipj);       // should be rji, just don't forget the minus! 
		    RjiXQj      = cross(Rij, Qj, ones);      // where: Rij = Rji
		    mirjiXQj    = cross(mirji, Qj, ones);
		    rjimiXQj    = cross(rjimi, Qj, ones);
		    rjiQirjiXQj = cross(rjiQirji, Qj, ones);
		    QirjirjiXQj = cross(Qirjirji, Qj, ones);

		    Gt1i = qj*rijXmi - mjXmi;
		    Gt2i = -(dipj*rij)*rijXmi + 2.0*QjrijXmi + QjI*rijXmi + 2.0*qj*RijXQi - 2.0*mjrijXQi - 2.0*rijmjXQi + 4.0*QjXQi;  
		    Gt3i = QjRij*rijXmi - 2.0*(dipj*rij)*RijXQi + 2.0*QjI*RijXQi + 4.0*rijQjrijXQi + 4.0*QjrijrijXQi;
		    Gt4i = 2.0*QjRij*RijXQi;
		
		    Gt1j = qi*rjiXmj + mjXmi;
		    Gt2j = -(dipi*rji)*rjiXmj - 2.0*QirijXmj + QiI*rjiXmj + 2.0*qi*RjiXQj - 2.0*mirjiXQj - 2.0*rjimiXQj - 4.0*QjXQi;
		    Gt3j = QiRij*rjiXmj - 2.0*(dipi*rji)*RjiXQj + 2.0*QiI*RjiXQj + 4.0*rjiQirjiXQj + 4.0*QirjirjiXQj;
		    Gt4j = 2.0*QiRij*RjiXQj;

		    trq[atomi->atomID] += COUlOMBS*( B1*Gt1i + B2*Gt2i + B3*Gt3i + B4*Gt4i ); 
		    trq[atomj->atomID] += COUlOMBS*( B1*Gt1j + B2*Gt2j + B3*Gt3j + B4*Gt4j ); 


		}	// If Torques
            }  // cut-off If sentence
        }   // j-loop

     }  // i-loop
   }  // parallel decomposition

   //-----------------------------------------------------------------------
   // NOTE: 1/3V factor of dUdV (realderiv) is applied in Ensemble.cpp
   // and in molVirt in LFIntegrator.cpp.
   //-----------------------------------------------------------------------
   
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

	    totTrq += trq[i].length();
	}

        //-----------------------------------------------------
        // Torques into forces.
        //-----------------------------------------------------

	double dUdV = 0.0;
        torques2forces(dUdV);

	realderiv = realderiv + dUdV;

    }
}

void Multipole::reciprocal_term(Double &pe, double &longderiv, double &long2ndDeriv)
{    
    if(doUpdate) {
        volume = myEnsemble->volume;
        volume2 = volume*volume;
        volumer = 1.0/volume;
    }
    int molIDi;
    clock_t start,end;
    Double pi4Volumer = 4.0*PI*volumer;
    Double pi8Volumer = 2.0*pi4Volumer;
    Double tmpenergy = 0.0;
    double dipxi, dipyi, dipzi;
    double qi, a;
    double Qxxi, Qxyi, Qxzi, Qyyi, Qyzi, Qzzi; 
    double Qyxi, Qzxi, Qzyi; 
    double UqQ, UmqQ, Umm;
    double dUdVqQ, dUdVqQm, dUdV3rd;
    double qiQisinSum, qiQicosSum, miksinSum, mikcosSum, qi2QicosSumV, qi2QisinSumV;
    double QiK, cosAqiQi, sinAqiQi, mikcos, miksin, qi2Qicos, qi2Qisin; 
    double b, Ff, Ffk2, d2UV1, d2UrR;
    double FqQ, FmiqQ, FqQmi, Fmm;
    double val, d2u1, d2u2;
    double TorqueArray[numAtoms][3], tempArray[numAtoms][3];
    Vector3 ri, rim;
    Vector3 kj, cosrR, sinrR;
    Vector3 dipi;
    Molecule *moli;

    pe = 0.0;
    longderiv = 0.0;
    long2ndDeriv = 0.0;

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
        qi2QicosSumV = 0.0;
        qi2QisinSumV = 0.0;
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
            a   = ri*kj;
	    QiK = Qxxi*kj.x*kj.x + Qxyi*kj.x*kj.y + Qxzi*kj.x*kj.z \
                + Qyxi*kj.y*kj.x + Qyyi*kj.y*kj.y + Qyzi*kj.y*kj.z \
                + Qzxi*kj.z*kj.x + Qzyi*kj.z*kj.y + Qzzi*kj.z*kj.z ; 

            cosAqiQi = (qi - QiK)*cos(a);
            sinAqiQi = (qi - QiK)*sin(a);
	    mikcos   = (dipi*kj)*cos(a);
	    miksin   = (dipi*kj)*sin(a);
            qi2Qicos = 2.0*QiK*cos(a);
            qi2Qisin = 2.0*QiK*sin(a);

    	    qQcos_i[i] = cosAqiQi;
    	    qQsin_i[i] = sinAqiQi;
    	    mkcos_i[i] = mikcos;
    	    mksin_i[i] = miksin;

            qiQicosSum += cosAqiQi;
            qiQisinSum += sinAqiQi; 
            mikcosSum  += mikcos;
            miksinSum  += miksin; 
            qi2QicosSumV += qi2Qicos;
            qi2QisinSumV += qi2Qisin;

	    cosrR  += cosAqiQi*(ri - rim); // deprecated
	    sinrR  += sinAqiQi*(ri - rim); // deprecated

	    //-----------------------------------------------------------------------------
	    // cross products for torque computation
	    //-----------------------------------------------------------------------------

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

		Vector3 kjXmi   = cross(kj, dipi);
		Vector3 kjkjXQi = cross(Kj, Qi, ones);

		kXmsina[i] = kjXmi*sin(a); 
		kXmcosa[i] = kjXmi*cos(a);
		KXQcosa[i] = kjkjXQi*cos(a);
		KXQsina[i] = kjkjXQi*sin(a);
	    }
	    

        }

	//-------------------------------------------------------------------------------------------------
	// Reciprocal-space Energy computation between all multipoles
	//-------------------------------------------------------------------------------------------------
	
    	UqQ   = pi4Volumer*expK[j]*(qiQicosSum*qiQicosSum + qiQisinSum*qiQisinSum); 
	UmqQ  = pi4Volumer*expK[j]*(qiQisinSum*mikcosSum - qiQicosSum*miksinSum);
	Umm   = pi4Volumer*expK[j]*(mikcosSum*mikcosSum + miksinSum*miksinSum);

        tmpenergy = (UqQ + 2*UmqQ + Umm)*COUlOMBS;
        pe += tmpenergy;

        dUdVqQ  = pi8Volumer*expK[j]*(qi2QicosSumV*qiQicosSum + qi2QisinSumV*qiQisinSum);
	dUdVqQm = pi8Volumer*expK[j]*(qi2QisinSumV*mikcosSum - qi2QicosSumV*miksinSum);
	dUdV3rd = (dUdVqQ - 2.0*UmqQ + dUdVqQm - 2.0*Umm)*COUlOMBS;

	longderiv += -tmpenergy + tmpenergy*kj*kj*2.0*alphaR4 + dUdV3rd;
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
	Ff = 0.0;
	Ffk2 = 0.0;
	d2UV1 = 0.0;
	d2UrR = 0.0;

        for (Int i = 0; i < numAtoms; i++)
        {
	  FqQ   =  b*(qQsin_i[i]*qiQicosSum - qQcos_i[i]*qiQisinSum);
	  FmiqQ =  b*(mkcos_i[i]*qiQicosSum + mksin_i[i]*qiQisinSum);
	  FqQmi = -b*(qQcos_i[i]*mikcosSum  + qQsin_i[i]*miksinSum);
	  Fmm   =  b*(mksin_i[i]*mikcosSum  - mkcos_i[i]*miksinSum);

          val = COUlOMBS*(FqQ + FmiqQ + FqQmi + Fmm);  // Total Force

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
 
		//-------------------------------------------------------------------------------------
		//  Computation of Torques
		//-------------------------------------------------------------------------------------

		if (computeTorques)
		{
		    Vector3 taum_qQ =  b*(qiQisinSum*kXmcosa[i] - qiQicosSum*kXmsina[i]) ;
		    Vector3 taum_m  =  b*(mikcosSum*kXmcosa[i]  + miksinSum*kXmsina[i]) ;
		    Vector3 tauQ_qQ = -2.0*b*(qiQicosSum*KXQcosa[i] + qiQisinSum*KXQsina[i]);
		    Vector3 tauQ_m  =  2.0*b*(miksinSum*KXQcosa[i]  - mikcosSum*KXQsina[i]) ;

		    trq[i] += COUlOMBS*( taum_qQ + taum_m + tauQ_qQ + tauQ_m );

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

   //-----------------------------------------------------------------------
   // NOTE: 1/3V factor of dUdV (longderiv) is applied in Ensemble.cpp
   // and in molVirt in LFIntegrator.cpp.
   //-----------------------------------------------------------------------
   
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

	    totTrq += trq[i].length();
	}

        //-----------------------------------------------------
        // Torques into forces.
        //-----------------------------------------------------

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
 
void Multipole::intraMolCorrect(Double &intraMolEnergy, double &intraVircorrect,  double &molderiv, double &mol2ndDeriv)
{
    Double intraMolecularEnergy = 0.0;
    Double r, r2, qi, qj;
    Int idi, idj, ns;
    int molIDi, molIDj; 
    double dipxi, dipyi, dipzi;
    double dipxj, dipyj, dipzj;
    double Qxxi, Qxyi, Qxzi, Qyyi, Qyzi, Qzzi; 
    double Qyxi, Qzxi, Qzyi; 
    double Qxxj, Qxyj, Qxzj, Qyyj, Qyzj, Qzzj; 
    double Qyxj, Qzxj, Qzyj; 
    double G0, G1, G2, G3, G4, G5;
    double D0, D1, D2, D3, D4, D5;
    double QiI, QjI, QiRij, QjRij;
    double QiQj, QimjXrij, QjmiXrij, QiXrijXQj; 
    double TorqueArray[numAtoms][3], tempArray[numAtoms][3];
    Vector3 dipi, dipj;
    Vector3 gradG1, gradG2, gradG3, gradG4;
    Vector3 Gt1i, Gt2i, Gt3i, Gt4i;
    Vector3 Gt1j, Gt2j, Gt3j, Gt4j;
    Vector3 QiXrij, QjXrij, QiXmj, QjXmi;
    Vector3 QiXQjrij, QjXQirij;
    Vector3 ri, rj, rij;

    intraMolEnergy = 0;
    intraVircorrect = 0.0;
    molderiv = 0.0;
    mol2ndDeriv = 0.0;

    if(doUpdate) {
        volume = myEnsemble->volume;
        volume2 = volume*volume;
        volumer = 1.0/volume;
    }

    // intra-molecular term. 
    // Molecular virial correction should be computed after all forces have been computed.  

    Int nmol = myEnsemble->nMols;

   //--------------------------------------------------------------------------
   //  First clean up the torques for all the processors
   //--------------------------------------------------------------------------

    for (int i = 0; i < numAtoms; i++)
    {
	trq[i].x = 0;
	trq[i].y = 0;
	trq[i].z = 0;
    }

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

		    //----------------------------------------------------------------------------
		    //  Constant needed for the G terms
		    //----------------------------------------------------------------------------

		    QiI   = Qxxi + Qyyi + Qzzi;
		    QjI   = Qxxj + Qyyj + Qzzj;

		    QiRij = Qxxi*rij.x*rij.x + Qxyi*rij.x*rij.y + Qxzi*rij.x*rij.z \
			  + Qyxi*rij.y*rij.x + Qyyi*rij.y*rij.y + Qyzi*rij.y*rij.z \ 
			  + Qzxi*rij.z*rij.x + Qzyi*rij.z*rij.y + Qzzi*rij.z*rij.z ; 

		    QjRij = Qxxj*rij.x*rij.x + Qxyj*rij.x*rij.y + Qxzj*rij.x*rij.z \
		          + Qyxj*rij.y*rij.x + Qyyj*rij.y*rij.y + Qyzj*rij.y*rij.z \
		          + Qzxj*rij.z*rij.x + Qzyj*rij.z*rij.y + Qzzj*rij.z*rij.z ;

		    QiQj  = Qxxi*Qxxj + Qxyi*Qxyj + Qxzi*Qxzj\
			  + Qyxi*Qyxj + Qyyi*Qyyj + Qyzi*Qyzj\
			  + Qzxi*Qzxj + Qzyi*Qzyj + Qzzi*Qzzj;

		    QimjXrij = dipj.x*(Qxxi*rij.x + Qxyi*rij.y + Qxzi*rij.z) \ 
			     + dipj.y*(Qyxi*rij.x + Qyyi*rij.y + Qyzi*rij.z) \ 
			     + dipj.z*(Qzxi*rij.x + Qzyi*rij.y + Qzzi*rij.z) ;

		    QjmiXrij = dipi.x*(Qxxj*rij.x + Qxyj*rij.y + Qxzj*rij.z) \ 
			     + dipi.y*(Qyxj*rij.x + Qyyj*rij.y + Qyzj*rij.z) \ 
			     + dipi.z*(Qzxj*rij.x + Qzyj*rij.y + Qzzj*rij.z) ;

		    QiXrijXQj = rij.x*Qxxj*(Qxxi*rij.x + Qxyi*rij.y + Qxzi*rij.z) + rij.y*Qxyj*(Qxxi*rij.x + Qxyi*rij.y + Qxzi*rij.z) + rij.z*Qxzj*(Qxxi*rij.x + Qxyi*rij.y + Qxzi*rij.z) \ 
			      + rij.x*Qyxj*(Qyxi*rij.x + Qyyi*rij.y + Qyzi*rij.z) + rij.y*Qyyj*(Qyxi*rij.x + Qyyi*rij.y + Qyzi*rij.z) + rij.z*Qyzj*(Qyxi*rij.x + Qyyi*rij.y + Qyzi*rij.z) \ 
			      + rij.x*Qzxj*(Qzxi*rij.x + Qzyi*rij.y + Qzzi*rij.z) + rij.y*Qzyj*(Qzxi*rij.x + Qzyi*rij.y + Qzzi*rij.z) + rij.z*Qzzj*(Qzxi*rij.x + Qzyi*rij.y + Qzzi*rij.z) ;

		    G0 = qi*qj;
		    G1 = qj*dipi*rij - qi*dipj*rij - dipi*dipj + qj*QiI + qi*QjI;
		    G2 = qj*QiRij + qi*QjRij - (dipi*rij)*(dipj*rij) - dipj*rij*QiI + dipi*rij*QjI -2*QimjXrij + 2*QjmiXrij + QiI*QjI + 2*QiQj; 
		    G3 = dipi*rij*QjRij - dipj*rij*QiRij + QiI*QjRij + QjI*QiRij + 4*QiXrijXQj;
		    G4 = QiRij*QjRij;

		    //---------------------------------------------------------------------------------------------------
		    //  D constants for Multipole Correction
		    //---------------------------------------------------------------------------------------------------

		    D0 = erf(alpha*r)/r;
		    D1 = (-D0 + alpha2PI*exp(-alpha2*r2))/r2;
		    D2 = (-3*D1 - 2*alpha2PI*alpha2*exp(-alpha2*r2) )/r2;
		    D3 = (-5*D2 + 4*alpha2PI*alpha2*alpha2*exp(-alpha2*r2) )/r2;
		    D4 = (-7*D3 - 8*alpha2PI*alpha2*alpha2*alpha2*exp(-alpha2*r2) )/r2;
		    D5 = (-9*D4 + 16*alpha2PI*alpha2*alpha2*alpha2*alpha2*exp(-alpha2*r2) )/r2;

                    intraMolecularEnergy += G0*D0 + G1*D1 + G2*D2 + G3*D3 + G4*D4;

                    // intra-molecular self force, may be no need when using constraint method
		    if (!useConstraint){

			QiXrij.x = Qxxi*rij.x + Qxyi*rij.y + Qxzi*rij.z; 
			QiXrij.y = Qyxi*rij.x + Qyyi*rij.y + Qyzi*rij.z; 
			QiXrij.z = Qzxi*rij.x + Qzyi*rij.y + Qzzi*rij.z; 

			QjXrij.x = Qxxj*rij.x + Qxyj*rij.y + Qxzj*rij.z;
			QjXrij.y = Qyxj*rij.x + Qyyj*rij.y + Qyzj*rij.z;
			QjXrij.z = Qzxj*rij.x + Qzyj*rij.y + Qzzj*rij.z;

			QiXmj.x  = Qxxi*dipj.x + Qxyi*dipj.y + Qxzi*dipj.z; 
			QiXmj.y  = Qyxi*dipj.x + Qyyi*dipj.y + Qyzi*dipj.z; 
			QiXmj.z  = Qzxi*dipj.x + Qzyi*dipj.y + Qzzi*dipj.z; 

			QjXmi.x  = Qxxj*dipi.x + Qxyj*dipi.y + Qxzj*dipi.z;
			QjXmi.y  = Qyxj*dipi.x + Qyyj*dipi.y + Qyzj*dipi.z;
			QjXmi.z  = Qzxj*dipi.x + Qzyj*dipi.y + Qzzj*dipi.z;

			QiXQjrij.x  = (Qxxi*Qxxj + Qxyi*Qyxj + Qxzi*Qzxj)*rij.x + (Qxxi*Qxyj + Qxyi*Qyyj + Qxzi*Qzyj)*rij.y + (Qxxi*Qxzj + Qxyi*Qyzj + Qxzi*Qzzj)*rij.z;
			QiXQjrij.y  = (Qyxi*Qxxj + Qyyi*Qyxj + Qyzi*Qzxj)*rij.x + (Qyxi*Qxyj + Qyyi*Qyyj + Qyzi*Qzyj)*rij.y + (Qyxi*Qxzj + Qyyi*Qyzj + Qyzi*Qzzj)*rij.z;
			QiXQjrij.z  = (Qzxi*Qxxj + Qzyi*Qyxj + Qzzi*Qzxj)*rij.x + (Qzxi*Qxyj + Qzyi*Qyyj + Qzzi*Qzyj)*rij.y + (Qzxi*Qxzj + Qzyi*Qyzj + Qzzi*Qzzj)*rij.z;

			QjXQirij.x  = (Qxxj*Qxxi + Qxyj*Qyxi + Qxzj*Qzxi)*rij.x + (Qxxj*Qxyi + Qxyj*Qyyi + Qxzj*Qzyi)*rij.y + (Qxxj*Qxzi + Qxyj*Qyzi + Qxzj*Qzzi)*rij.z;
			QjXQirij.y  = (Qyxj*Qxxi + Qyyj*Qyxi + Qyzj*Qzxi)*rij.x + (Qyxj*Qxyi + Qyyj*Qyyi + Qyzj*Qzyi)*rij.y + (Qyxj*Qxzi + Qyyj*Qyzi + Qyzj*Qzzi)*rij.z;
			QjXQirij.z  = (Qzxj*Qxxi + Qzyj*Qyxi + Qzzj*Qzxi)*rij.x + (Qzxj*Qxyi + Qzyj*Qyyi + Qzzj*Qzyi)*rij.y + (Qzxj*Qxzi + Qzyj*Qyzi + Qzzj*Qzzi)*rij.z;


			gradG1 = qj*dipi - qi*dipj;
			gradG2 = 2*qj*QiXrij - (dipj*rij)*dipi - (dipi*rij)*dipj - 2*QiXmj - QiI*dipj + 2*qi*QjXrij + 2*QjXmi + QjI*dipi; 
			gradG3 = -QiRij*dipj - (2*dipj*rij)*QiXrij + QjRij*dipi + (2*dipi*rij)*QjXrij + (2*QiI)*QjXrij + (2*QjI)*QiXrij + 4*QiXQjrij + 4*QjXQirij;
			gradG4 = (2*QjRij)*QiXrij + (2*QiRij)*QjXrij;

                        Vector3 fij  = D1*gradG1 + D2*gradG2 + D3*gradG3 + D4*gradG4 + (G0*D1 + G1*D2 + G2*D3 + G3*D4 + G4*D5)*rij;
    			molderiv    += -(fij*rij)*COUlOMBS/(3.0*volume);
    			mol2ndDeriv += 0.0;

                        atoms[idi].force += fij*COUlOMBS;
                        atoms[idj].force -= fij*COUlOMBS;

			if (computeTorques){

			    double Qi[3][3], Qj[3][3], Rij[3][3];
			    double mjrij[3][3], rijmj[3][3], rijQjrij[3][3], Qjrijrij[3][3];
			    double mirji[3][3], rjimi[3][3], rjiQirji[3][3], Qirjirji[3][3];

			    Vector3 ones(1,1,1);
			    Vector3 rijXmi, mjXmi, RijXQi, QjXQi, mjrijXQi, rijmjXQi;
			    Vector3 QjrijXmi, rijQjrijXQi, QjrijrijXQi;
			    Vector3 rji, rjiXmj, RjiXQj, mirjiXQj, rjimiXQj;
			    Vector3 QirijXmj, rjiQirjiXQj, QirjirjiXQj;

			    rji = -rij;

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

			    mjrij[0][0] = dipj.x*rij.x; 
			    mjrij[0][1] = dipj.x*rij.y; 
			    mjrij[0][2] = dipj.x*rij.z; 
			    mjrij[1][0] = dipj.y*rij.x; 
			    mjrij[1][1] = dipj.y*rij.y; 
			    mjrij[1][2] = dipj.y*rij.z; 
			    mjrij[2][0] = dipj.z*rij.x; 
			    mjrij[2][1] = dipj.z*rij.y; 
			    mjrij[2][2] = dipj.z*rij.z; 

			    rijmj[0][0] = mjrij[0][0];
			    rijmj[0][1] = mjrij[1][0];
			    rijmj[0][2] = mjrij[2][0];
			    rijmj[1][0] = mjrij[0][1];
			    rijmj[1][1] = mjrij[1][1];
			    rijmj[1][2] = mjrij[2][1];
			    rijmj[2][0] = mjrij[0][2];
			    rijmj[2][1] = mjrij[1][2];
			    rijmj[2][2] = mjrij[2][2];

			    Qjrijrij[0][0] = QjXrij.x*rij.x;
			    Qjrijrij[0][1] = QjXrij.x*rij.y;
			    Qjrijrij[0][2] = QjXrij.x*rij.z;
			    Qjrijrij[1][0] = QjXrij.y*rij.x;
			    Qjrijrij[1][1] = QjXrij.y*rij.y;
			    Qjrijrij[1][2] = QjXrij.y*rij.z;
			    Qjrijrij[2][0] = QjXrij.z*rij.x;
			    Qjrijrij[2][1] = QjXrij.z*rij.y;
			    Qjrijrij[2][2] = QjXrij.z*rij.z;

			    rijQjrij[0][0] = Qjrijrij[0][0];
			    rijQjrij[0][1] = Qjrijrij[1][0];
			    rijQjrij[0][2] = Qjrijrij[2][0];
			    rijQjrij[1][0] = Qjrijrij[0][1];
			    rijQjrij[1][1] = Qjrijrij[1][1];
			    rijQjrij[1][2] = Qjrijrij[2][1];
			    rijQjrij[2][0] = Qjrijrij[0][2];
			    rijQjrij[2][1] = Qjrijrij[1][2];
			    rijQjrij[2][2] = Qjrijrij[2][2];


			    mirji[0][0] = dipi.x*rji.x; 
			    mirji[0][1] = dipi.x*rji.y; 
			    mirji[0][2] = dipi.x*rji.z; 
			    mirji[1][0] = dipi.y*rji.x; 
			    mirji[1][1] = dipi.y*rji.y; 
			    mirji[1][2] = dipi.y*rji.z; 
			    mirji[2][0] = dipi.z*rji.x; 
			    mirji[2][1] = dipi.z*rji.y; 
			    mirji[2][2] = dipi.z*rji.z; 

			    rjimi[0][0] = mirji[0][0]; 
			    rjimi[0][1] = mirji[1][0]; 
			    rjimi[0][2] = mirji[2][0]; 
			    rjimi[1][0] = mirji[0][1]; 
			    rjimi[1][1] = mirji[1][1]; 
			    rjimi[1][2] = mirji[2][1]; 
			    rjimi[2][0] = mirji[0][2]; 
			    rjimi[2][1] = mirji[1][2]; 
			    rjimi[2][2] = mirji[2][2]; 

			    Qirjirji[0][0] = QiXrij.x*rij.x;
			    Qirjirji[0][1] = QiXrij.x*rij.y;
			    Qirjirji[0][2] = QiXrij.x*rij.z;
			    Qirjirji[1][0] = QiXrij.y*rij.x;
			    Qirjirji[1][1] = QiXrij.y*rij.y;
			    Qirjirji[1][2] = QiXrij.y*rij.z;
			    Qirjirji[2][0] = QiXrij.z*rij.x;
			    Qirjirji[2][1] = QiXrij.z*rij.y;
			    Qirjirji[2][2] = QiXrij.z*rij.z;

			    rjiQirji[0][0] = Qirjirji[0][0];
			    rjiQirji[0][1] = Qirjirji[1][0];
			    rjiQirji[0][2] = Qirjirji[2][0];
			    rjiQirji[1][0] = Qirjirji[0][1];
			    rjiQirji[1][1] = Qirjirji[1][1];
			    rjiQirji[1][2] = Qirjirji[2][1];
			    rjiQirji[2][0] = Qirjirji[0][2];
			    rjiQirji[2][1] = Qirjirji[1][2];
			    rjiQirji[2][2] = Qirjirji[2][2];


			    rijXmi      = cross(rij,dipi);
			    mjXmi       = cross(dipj,dipi);
			    QjrijXmi    = cross(QjXrij, dipi);
			    RijXQi      = cross(Rij, Qi, ones);
			    mjrijXQi    = cross(mjrij, Qi, ones);
			    rijmjXQi    = cross(rijmj, Qi, ones); 
			    QjXQi       = cross(Qj, Qi, ones);
			    rijQjrijXQi = cross(rijQjrij, Qi, ones);
			    QjrijrijXQi = cross(Qjrijrij, Qi, ones);

			    rjiXmj      = cross(rji, dipj);
			    QirijXmj    = cross(QiXrij, dipj);       // should be rji, just don't forget the minus! 
			    RjiXQj      = cross(Rij, Qj, ones);      // where: Rij = Rji
			    mirjiXQj    = cross(mirji, Qj, ones);
			    rjimiXQj    = cross(rjimi, Qj, ones);
			    rjiQirjiXQj = cross(rjiQirji, Qj, ones);
			    QirjirjiXQj = cross(Qirjirji, Qj, ones);

			    Gt1i = qj*rijXmi - mjXmi;
			    Gt2i = -(dipj*rij)*rijXmi + 2.0*QjrijXmi + QjI*rijXmi + 2.0*qj*RijXQi - 2.0*mjrijXQi - 2.0*rijmjXQi + 4.0*QjXQi;  
			    Gt3i = QjRij*rijXmi - 2.0*(dipj*rij)*RijXQi + 2.0*QjI*RijXQi + 4.0*rijQjrijXQi + 4.0*QjrijrijXQi;
			    Gt4i = 2.0*QjRij*RijXQi;
			
			    Gt1j = qi*rjiXmj + mjXmi;
			    Gt2j = -(dipi*rji)*rjiXmj - 2.0*QirijXmj + QiI*rjiXmj + 2.0*qi*RjiXQj - 2.0*mirjiXQj - 2.0*rjimiXQj - 4.0*QjXQi;
			    Gt3j = QiRij*rjiXmj - 2.0*(dipi*rji)*RjiXQj + 2.0*QiI*RjiXQj + 4.0*rjiQirjiXQj + 4.0*QirjirjiXQj;
			    Gt4j = 2.0*QiRij*RjiXQj;


			    //--------------------------------------------------------------------------------------------------------------------
			    // Computation of Torques
			    //--------------------------------------------------------------------------------------------------------------------

		    	    trq[idi] += -COUlOMBS*( D1*Gt1i + D2*Gt2i + D3*Gt3i + D4*Gt4i ); 
		    	    trq[idj] += -COUlOMBS*( D1*Gt1j + D2*Gt2j + D3*Gt3j + D4*Gt4j ); 

			}  // if Torques
                    }  // If constraints (forces and torques)
                }
            }          
        }
      }
    }
    intraMolEnergy = -intraMolecularEnergy*COUlOMBS;
    intraVircorrect = 0.0;

    //----------------------------------------------------------------
    // Note: the 1/3V factor for dUdV (molderiv) 
    // is applied here!!
    //----------------------------------------------------------------
    
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

	    totTrq += trq[i].length();
	}

        //-----------------------------------------------------
        // Torques into forces.
        //-----------------------------------------------------

	double dUdV = 0.0;
        torques2forces(dUdV);

	molderiv += dUdV/(3.0*volume);
    }

}

//  surface dipole correction
void Multipole::surfDipoleCorrect(Double &surfaceDipoleEnergy, double &surfderiv, double &surf2ndderiv)
{
    if(doUpdate) {
        volume = myEnsemble->volume;
        volume2 = volume*volume;
        volumer = 1.0/volume;
    }
    Double coeff = 2.0*PI*volumer/3.0;
    Double r, r2, qi;
    int idi, idj, ns;
    int molIDi, atomID; 
    double dipxi, dipyi, dipzi;
    double Qxxi, Qyyi, Qzzi; 
    double  selfSum = 0.0;
    double QiI, selfterm;
    double TorqueArray[numAtoms][3], tempArray[numAtoms][3];
    Vector3 dipi;
    Vector3 qrSum, qRsum, rim;
    Vector3 miSum, qirimi;
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
	// Call Traceless Quadrupoles from Local2Global 	
	//-----------------------------------------------------------------------
	
	Qxxi  = myL2G->Qxx[i];
	Qyyi  = myL2G->Qyy[i];
	Qzzi  = myL2G->Qzz[i];

	QiI = Qxxi + Qyyi + Qzzi;
        
        // myEnsemble->apply_pbc(r);  //Jc: why?
        qrSum += qi*r;
        qRsum += qi*rim;
        miSum += dipi;
	qirimi += qi*r + dipi;
        selfSum += 2*qi*QiI - dipi*dipi;

    }

    selfterm = qirimi*qirimi + selfSum;

    surfaceDipoleEnergy = coeff*selfterm*COUlOMBS;

    if (useConstraint){
        surfderiv           = -surfaceDipoleEnergy*volumer + 2.0*coeff*volumer*qirimi*qrSum*COUlOMBS/3.0;
        surf2ndderiv        = 2.0*coeff*volumer*volumer*qrSum*qrSum - 16.0*coeff*volumer*volumer*qrSum*qRsum/9.0 + 2.0*coeff*volumer*volumer*qRsum*qRsum/9.0;
    }
    else {
        surfderiv           = -surfaceDipoleEnergy*volumer + 2.0*coeff*volumer*qirimi*qrSum*COUlOMBS/3.0;
        surf2ndderiv        = 4.0*coeff*volumer*volumer*qrSum*qrSum/9.0;
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

        ri = atoms[atomID].position;
        if(!strcmp(myMols[molIDi].molName, "MCY") && atoms[atomID].atomType ==1)   // Oxygen location MCY
        { 
            ri.x = Dummy[molIDi*3];
            ri.y = Dummy[molIDi*3 +1];
            ri.z = Dummy[molIDi*3 +2]; 
        }

        fi = -2.0*coeff*qi*qirimi;
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

	//-------------------------------------------------------------------------------------
	//  Computation of Torques
	//-------------------------------------------------------------------------------------

	if (computeTorques)
	{
	    trq[atomID] += 2.0*coeff*COUlOMBS*( cross(qirimi, dipi) ); 
	}
      }
    }
 
    //--------------------------------------------------------------------------
    // Note: the 1/3V factor of dU/dV (surfderiv) is applied here!
    //--------------------------------------------------------------------------

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

	    totTrq += trq[i].length();
	}

        //-----------------------------------------------------
        // Torques into forces.
        //-----------------------------------------------------

	double dUdV = 0.0;
        torques2forces(dUdV);

	surfderiv += dUdV*volumer/3.0;
    }

}

// JC string EwaldForce::get_force_id() // in order to comply with the cluster compiler
// JC {
// JC  return EwaldForce::forceIdentifier;
// JC}

void Multipole::torques2forces(double &dUdVol)
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

	//--------------------------------------------------------------------
	//  Transformation for H1. The distances are ri - rj 
	//  in order to mantain consistency of the code.
	//  Note: this method can be written using the .atomType object and
	//  using a "for-loop" with ns, if you are reading this, plz go ahead.
	//--------------------------------------------------------------------

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


	//-----------------------------------------------------------------
	//  Transformation for O. The distances are ri - rj 
	//  in order to mantain consistency of the code.
	//-----------------------------------------------------------------

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

	//-----------------------------------------------------------------
	//  Transformation for H2. The distances are ri - rj 
	//  in order to mantain consistency of the code.
	//-----------------------------------------------------------------

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


void Multipole::write_force_info(ofstream& of)
{   
    of << "Multipole Force created" << endl;
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

void Multipole::write_energy(ofstream& of)
{
    of << "Multipole Energy: " << energy+longEnergy << endl;
    of << "realEnergy: " << realEnergy << '\t' << "longEnergy: " << longEnergy << endl ;
    // of << "introMolEnergy: " << introMolEnergy << '\t' << "surfaceEnergy" << surfDipoleEnergy << endl;
}


