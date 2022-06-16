/** EwaldForce.cpp -- 
 **
 ** Copyright (C) 2003
 ** Centre for Molecular Simulation (CMS)
 ** School of Information Technology
 ** Swinburne University of Technology
 ** PO Box 21 Hawthorn, Vic 3122, Australia
 **
 ** Author: Zhongwu Zhou
 ** Email: zzhou@it.swin.edu.au
 ** Ewald Volume Derivatives and NpT Update by
 ** Author: eoyarzua
 ** Email: eoyarzua@swin.edu.au
 **/

#include "EwaldForce.h" 

#include "../utils/erf.h"
#include <mpi.h>

EwaldForce::EwaldForce(Ensemble* ensemble) : Force(ensemble) 
{
    DEBUGMSG("Creating EwaldForce");

//JC    forceIdentifier = "EwaldForce"; // in order to comply with the cluster compiler
    numAtoms = myEnsemble->nAtoms;
    numMols = myEnsemble->nMols;
    myMols = myEnsemble->molecules;
    volume = myEnsemble->volume;
cout << "VOLUME = " << volume << endl;
    volume2 = volume*volume;
    volumer = 1.0/volume;
    alpha = myEnsemble->myConfig->get_alpha();
    alpha2 = alpha*alpha;
    alpha2PI = 2.0*alpha/SQRT_PI;
    alphaR4 = 1.0/(4.0*alpha2);                 
    realCut = myEnsemble->myConfig->get_cutoff();
    realCut2 = realCut*realCut;
    kCut = myEnsemble->myConfig->get_k_cutoff(); 
    correctSurfDipole = myEnsemble->myConfig->compute_surf_correct(); 
    useConstraint = myEnsemble->myConfig->is_constraint_on();
    statEnsem = myEnsemble->myConfig->get_ensemble_status();  

    stepr = 0.0;
    selfEnergy = 0.0;
    chargedEnergy = 0.0;
    longEnergy = 0.0;
    longderiv  = 0.0;
    long2ndDeriv  = 0.0;
    tableSize = TAB_SIZE + 2;

    expK = NULL;
    sinCosA = NULL;
    forceTable = NULL;
    potentialTable = NULL;
    doUpdate = false;
    if (statEnsem == 1){
        doUpdate = true;
        DEBUGMSG(" Ewald Volume Update for NpT ensemble");
    }


    Dummy = new Double[numMols*3];

    set_up();
    DEBUGMSG("EwaldForce set up");
    init_table(2*numAtoms, &sinCosA);

    // pre-compute self energy correction - Gaussian term
    Double qSum = 0.0;
    Double q2Sum = 0.0;
    for(Int i = 0; i < numAtoms; i++)
    {
        qSum += atoms[i].scaledCharge;
        q2Sum += atoms[i].scaledCharge*atoms[i].scaledCharge;
    }
    selfEnergy = -q2Sum*alpha/SQRT_PI;
    if(fabs(qSum) > 1.0e-5)
       chargedEnergy = -PI*qSum*qSum/(2.0*volume*alpha2);  // Volume needs to be updated in NpT cases with highly charged systems (Ions probably?)
    
    write_force_info(*sysdataFile);
 
   // prank = MPI::COMM_WORLD.Get_rank();  // Jc: get the MPI parameters // Updated to MPI-3 by elton
    //psize = MPI::COMM_WORLD.Get_size();

    MPI_Comm_rank(MPI_COMM_WORLD, &prank);
    MPI_Comm_size(MPI_COMM_WORLD, &psize);

}

EwaldForce::~EwaldForce()
{
    if(expK != NULL) delete [] expK;
    if(sinCosA != NULL) delete [] sinCosA;
    if(forceTable != NULL) delete [] forceTable;
    if(potentialTable != NULL) delete [] potentialTable;    
}

void EwaldForce::DummyPosition()

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

void EwaldForce::compute()
{
    energy = 0.0;
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

    DummyPosition();

    for (Int i = XX; i <= ZZ; i++)
        virial[i] = 0.0;

    if(doUpdate)
        set_up();

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
 
void EwaldForce::compute_long()
{
    DummyPosition();    

    for (Int i = XX; i <= ZZ; i++)
        longVirial[i] = 0.0;

    for (Int i = 0; i < numAtoms; i++)
        atoms[i].longForce = 0.0;

    reciprocal_term(longEnergy, longderiv, long2ndDeriv);
    myEnsemble->myEnergy.longEnergy = longEnergy;
    
    myEnsemble->myLustig.longDeriv = longderiv;
    myEnsemble->myLustig.longSecDeriv = long2ndDeriv;

    for(Int i = 0; i < ZZ; i++)
        myEnsemble->longVirial[i] = longVirial[i];
}

void EwaldForce::set_up()
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
void EwaldForce::build_look_up_table()
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

void EwaldForce::init_table(Int size, Double** table)
{
    if (*table != NULL)
        delete [] *table;
    *table = new Double[size];
    if(*table == NULL) 
        ERRORMSG("EwaldForce::init_table(), not enough memory");
}


//Jc:  compute real term in Ewald summation,it is necessary as well as the reciprocal term 
void EwaldForce::real_term(Double &realEnergy, double &realderiv, double &real2ndDeriv)
{

    #ifdef DEBUG2
        DEBUGMSG("Computing Ewald real space");
    #endif

    Atom *atomi, *atomj;
    Molecule *moli, *molj;
    Vector3 ri, rj, fij, rij, rijm, rim, rjm;
    Int size;
    Double qi, qij;
    
    realEnergy = 0.0;
    realderiv  = 0.0;
    real2ndDeriv  = 0.0;

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
        
        int molIDi = atomi->molID;

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
        qi = atomi->scaledCharge;

        for (int j = 0; j < size; j++)
        {

            atomj = atomi->myPairList[j];
            int molIDj = atomj->molID;

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
            qij = qi*atomj->scaledCharge;
            rij = rj - ri;


            myEnsemble->apply_pbc(rij);
            myEnsemble->apply_pbc(rijm);

            double r2 = rij.length2();

//if (atomi->atomType == 1 && atomj->atomType == 1){
//cout << rij << " <- atom      Rij ->   " << endl;
//cout << rijm  << endl;}


            if (r2 < realCut2)
            {
                #ifdef USE_LOOK_UP_TABLE
                {
                    realEnergy += qij*lookup(r2, potentialTable);
                    fij = -qij*lookup(r2, forceTable)*rij;
                }
                #else
                {
		    double rijRij;
		    if (useConstraint){
			rijRij = rij*rijm;
			//cout << "ENCENDIDO" << endl;
		    }
		    else {
			rijRij = r2;
			//cout << "APAGADO" << endl;
		    }
                    double r = sqrt(r2);
                    double e = erfc(alpha*r)/r;
                    double f = (e + alpha2PI*exp(-alpha2*r2))/r2;
		   /* if (!strcmp(myMols[molIDi].molName, "BBV")  && !strcmp(myMols[molIDj].molName, "BBV")) {
			double expBBV = 1.0 + exp(-2.0*(r/Aa0 - 2.0));
			double invexp = 1.0/expBBV;
			double Dij2   = invexp*invexp;
			double Dij8   = Dij2*Dij2*Dij2*Dij2;
			double Dij15  = Dij8*Dij2*Dij2*Dij2*invexp;
			   e = e*Dij15;
			   f = f*Dij15 - (30.0/Aa0)*exp(-2.0*(r/Aa0 - 2.0))*e*invexp/r;
		    }*/ 
                    realEnergy  += qij*e;
		    realderiv   += -qij*f*rijRij; 
		    real2ndDeriv += qij*(2.0*f + 2.0*alpha2PI*alpha2*exp(-alpha2*r2))*rijRij*rijRij/r2 + 2.0*qij*f*rijRij;

                    fij = -qij*f*rij;
                }
                #endif

// Jc: distribute force among molecule atoms, when dummy charge is calulated
// Jc: the calculated results are not improving. so there is no force distribution 
// Jc: in this calculation
                fi += fij;
		moli =  &myMols[molIDi]; 
		molj =  &myMols[molIDj];; 
                Atom *atomh1i = moli->myAtoms[0];
                Atom *atomh2i = moli->myAtoms[2];

                Atom *atomh1j = molj->myAtoms[0];
                Atom *atomh2j = molj->myAtoms[2];

                if ( !strcmp(myMols[molIDi].molName, "MCY") &&  (atomi->atomType == 1) ){
                        atomi->force += 0.5431676*fij;
                        atomh1i->force += 0.5*0.4568324*fij;
                        atomh2i->force += 0.5*0.4568324*fij;
                }
                else {
			atomi->force += fij;
                }

                if ( !strcmp(myMols[molIDj].molName, "MCY") && (atomj->atomType == 1) ){
                        atomj->force -= 0.5431676*fij;
                        atomh1j->force -= 0.5*0.4568324*fij;
                        atomh2j->force -= 0.5*0.4568324*fij;
                }
                else {
                        atomj->force -= fij;
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
            }
        }

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
}

void EwaldForce::reciprocal_term(Double &pe, double &longderiv, double &long2ndDeriv)
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
        Double sinSum = 0.0; 
        Double cosSum = 0.0;
        Vector3 kj = kVector[j];
        Vector3 cosrR, sinrR;
        Double tmpVir;
        Double coeff = 2.0*(alphaR4 + 1.0/kModulu[j]);         // for virial compute

	cosrR.x = 0.0, cosrR.y = 0.0, cosrR.z = 0.0;
	sinrR.x = 0.0, sinrR.y = 0.0, sinrR.z = 0.0;
 
        for (Int i = 0; i < numAtoms; i++)
        {
	    int molIDi = atoms[i].molID;
	    rim = myEnsemble->molecules[molIDi].massCenter;
            Double qi = atoms[i].scaledCharge;
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
            Double a  = ri*kj;
            Double sinAqi = qi*sin(a);
            Double cosAqi = qi*cos(a);
            sinCosA[i*2  ] = sinAqi;
            sinCosA[i*2+1] = cosAqi;
            sinSum += sinAqi;
            cosSum += cosAqi;
	    cosrR  += cosAqi*(ri - rim);
	    sinrR  += sinAqi*(ri - rim);
        }

        tmpenergy = pi4Volumer*expK[j]*(sinSum*sinSum + cosSum*cosSum);
        pe += tmpenergy;

	longderiv += -tmpenergy + tmpenergy*kj*kj*2.0*alphaR4;
	long2ndDeriv += 4.0*tmpenergy - 7.0*tmpenergy*kj*kj*2.0*alphaR4 + tmpenergy*kj*kj*2.0*alphaR4*kj*kj*2.0*alphaR4;

        longVirial[XX] += tmpenergy - tmpenergy*coeff*kj.x*kj.x;
        longVirial[YY] += tmpenergy - tmpenergy*coeff*kj.y*kj.y;
        longVirial[ZZ] += tmpenergy - tmpenergy*coeff*kj.z*kj.z;
        longVirial[XY] -= tmpenergy*coeff*kj.x*kj.y;
        longVirial[XZ] -= tmpenergy*coeff*kj.x*kj.z;
        longVirial[YZ] -= tmpenergy*coeff*kj.y*kj.z;
        
        Double b = pi8Volumer*expK[j];
	double Ff = 0.0, Ffk2 = 0.0, d2UV1 = 0.0, d2UrR = 0.0;
        for (Int i = 0; i < numAtoms; i++)
        {
          Double val = b*(sinCosA[i*2]*cosSum - sinCosA[i*2+1]*sinSum);  // sinCosA[]=qi*sin(k*ri)
          double d2u1 = b*kj*kj*(sinCosA[i*2+1]*cosSum + sinCosA[i*2]*sinSum);  // sinCosA[]=qi*sin(k*ri)
	  int molIDi = atoms[i].molID;
	  rim = myEnsemble->molecules[molIDi].massCenter;
	  ri = atoms[i].position;

          double d2u2 = b*kj*kj*(sinCosA[i*2+1]*cosrR + sinCosA[i*2]*sinrR)*(ri - rim);  // sinCosA[]=qi*sin(k*ri)
//Jc: if Dummy, then distribute force among molecular memebers

/*          if(atoms[i].atomType == 1)
          {
             double bi    = R_OP*Bohr_Radius/1.0e-9;        // 0.026765499 ;//
             double gama = -bi/(0.09572 * cos(0.912109067)); // gama= 0.456840908//
             Vector3 p, rod, ro, rd;
             ro = atoms[i].position; 
             rd.x = Dummy[atoms[i].molID*3];
             rd.y = Dummy[atoms[i].molID*3+1];
             rd.z = Dummy[atoms[i].molID*3+2];
             rod  = rd - ro;
             p    = (rod*kVector[j]*val)/(bi*bi)*rod;
             p    = gama*(kVector[j]*val - p);

             atoms[i].force      += (kVector[j]*val -  p);
             atoms[i-1].force    += 0.5 * p;
             atoms[i+1].force    += 0.5 * p;
*/
/*
            atoms[i-1].force     += kVector[j]*val*0.068455;
            atoms[i].force       += kVector[j]*val*0.86309;
            atoms[i+1].force     += kVector[j]*val*0.068455;   
*/
//          }
//          else
//          { 

	        moli =  &myMols[molIDi];
                Atom *atomh1i = moli->myAtoms[0];
                Atom *atomh2i = moli->myAtoms[2];

                if (!strcmp(myMols[atoms[i].molID].molName, "MCY") && (atoms[i].atomType == 1) ){
                        atoms[i].force += 0.5431676*kVector[j]*val;
                        atomh1i->force += 0.5*0.4568324*kVector[j]*val;
                        atomh2i->force += 0.5*0.4568324*kVector[j]*val;
                }
                else {
                        atoms[i].force += kVector[j]*val;
                }

             //atoms[i].force += kVector[j]*val;  // The q is included in SinCosA
	     Ff +=  val * ( kVector[j]*(ri - rim) );
	     Ffk2 += 2.0*val*kj*kj*2.0*alphaR4 * ( kVector[j]*(ri - rim) ); 
	     d2UV1 += d2u1*(ri - rim)*(ri - rim);
	     d2UrR += d2u2;
 
//          } 
//Jc: ZW code            atoms[i].longForce += kVector[j]*val;
        }

	if (useConstraint){
	    longderiv += Ff;
	//long2ndDeriv += -4.0*(-tmpenergy + tmpenergy*kj*kj*2.0*alphaR4 + Ff) - 3.0*tmpenergy*kj*kj*2.0*alphaR4 + tmpenergy*kj*kj*2.0*alphaR4*kj*kj*2.0*alphaR4 - 2.0*Ff + Ffk2 - d2UV1 + d2UrR;
	    long2ndDeriv += -6.0*Ff + Ffk2 - d2UV1 + d2UrR;
	   //cout << "  !!!!!!  ENCENDIDO" << endl;
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

}

// intra-molecular self energy correction 
// Jc: it is one of the biggest energy contribution but some times, it is balanced off by 
// Jc: the point self-energy, the self point self-energy is calcualted only once in the EwaldForce.cpp
// Jc: by the compute() method
// Jc: the self-energy is not related to the any positions in the system, but only the chrages 
// Jc: however the intraMolecularEnergy, also called nolCorrection are charge-position depedent
// Jc: Since the positions of charges in MCY model is different from SPC/E, Dummy positions need to be 
// Jc: introduced to account for the molCorrection calculation
 
void EwaldForce::intraMolCorrect(Double &intraMolEnergy, double &intraVircorrect,  double &molderiv, double &mol2ndDeriv)
{
    Double intraMolecularEnergy = 0.0;
    double U_i = 0.0;
    Double r, r2, qi, qij;
    Int idi, idj, ns;
    int midi, midj; 
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
            midi = myEnsemble->molecules[molIndex].myAtoms[i]->molID;
            qi = atoms[idi].scaledCharge;
// Jc: substitute the Oxygen position with the Dummdy Charge position
            Vector3 ri = atoms[idi].position;
            if(!strcmp(myMols[midi].molName, "MCY") && atoms[idi].atomType ==1)
            { 
              ri.x = Dummy[midi*3];
              ri.y = Dummy[midi*3 +1];
              ri.z = Dummy[midi*3 +2]; 
            }
            for (Int j = i+1; j < ns; j++)
            {
                idj = myEnsemble->molecules[molIndex].myAtoms[j]->atomID;
                midj = myEnsemble->molecules[molIndex].myAtoms[j]->molID;
                if (myEnsemble->exclusion_check(idi, idj))
                {
                    Vector3 rj = atoms[idj].position;
                    if(!strcmp(myMols[midj].molName, "MCY") && atoms[idj].atomType ==1)
                    { 
                     rj.x = Dummy[3*midj];
                     rj.y = Dummy[3*midj +1];
                     rj.z = Dummy[3*midj +2]; 
                    }

                    Vector3 rij = rj - ri;
                    myEnsemble->apply_pbc(rij);
                    r2 = rij.length2();
                    r = sqrt(r2);
                    qij = qi*atoms[idj].scaledCharge;
                    Double e = erf(alpha*r)/r;
                    intraMolecularEnergy += qij*e;
		    U_i += alpha2PI*qij*exp(-alpha2*r2);
		    //cout << "vol = "  << volume2  << endl;
                    // intra-molecular self force, may be no need when using constraint method
		    if (!useConstraint){
			double dudr  = e - alpha2PI*exp(-alpha2*r2);
    			molderiv    += qij*dudr/(3.0*volume);
    			mol2ndDeriv += ( -2.0*qij*dudr + qij*(-2.0*dudr + 2.0*alpha2PI*alpha2*exp(-alpha2*r2)*r2) ) / (9.0*volume2);
                        Vector3 fij  = qij*dudr/r2*rij;
                        atoms[idi].force += fij;
                        atoms[idj].force -= fij;
		        //cout << " HOLAA CTM!!!!"  << volume2  << endl;
                    }
                }
            }          
        }
      }
    }
    intraMolEnergy = -intraMolecularEnergy;
    intraVircorrect = -U_i;
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
}

//  surface dipole correction
void EwaldForce::surfDipoleCorrect(Double &surfaceDipoleEnergy, double &surfderiv, double &surf2ndderiv)
{
    if(doUpdate) {
        volume = myEnsemble->volume;
        volume2 = volume*volume;
        volumer = 1.0/volume;
    }
    Double coeff = 2.0*PI*volumer/3.0;
    Double r, r2, qi, qij;
    Int idi, idj, ns;
    Vector3 qrSum, qRsum, rim;
    int midi,midj;
    Molecule *moli;

    surfaceDipoleEnergy = 0.0;
    surfderiv           = 0.0;
    surf2ndderiv        = 0.0;

    // surface dipole term, no virial apply to this term

    for (Int i = 0; i < numAtoms; i++)
    {
        Vector3 r = atoms[i].position;

        midi = atoms[i].molID;
        rim = myEnsemble->molecules[midi].massCenter;
        if(!strcmp(myMols[midi].molName, "MCY") && atoms[i].atomType ==1)   // Jc: added by Jianhui
        { 
          r.x = Dummy[midi*3];
          r.y = Dummy[midi*3 +1];
          r.z = Dummy[midi*3 +2]; 
        }
        
        // myEnsemble->apply_pbc(r);  //Jc: why?
        qrSum += atoms[i].scaledCharge*r;
        qRsum += atoms[i].scaledCharge*rim;
    }

    surfaceDipoleEnergy = coeff*qrSum.length2();
    if (useConstraint){
        surfderiv           = -coeff*volumer*qrSum*qrSum + 2.0*coeff*volumer*qrSum*qRsum/3.0;
        surf2ndderiv        = 2.0*coeff*volumer*volumer*qrSum*qrSum - 16.0*coeff*volumer*volumer*qrSum*qRsum/9.0 + 2.0*coeff*volumer*volumer*qRsum*qRsum/9.0;
        //cout << "ENCENDIDO!!"<< endl;
    }
    else {
        surfderiv           = -coeff*volumer*qrSum*qrSum/3.0;
        surf2ndderiv        = 4.0*coeff*volumer*volumer*qrSum*qrSum/9.0;
        //cout << "APAGADO!!"<< endl;
    }
    qrSum = coeff*qrSum;                         // need to check force calculation 

// Jc: parallel code;
    Vector3 temp1;
    int numDecompose = (int) numAtoms/psize;
    if ((numAtoms % psize) == 0) numDecompose--;
    for(int i = 0; i <= numDecompose;i++)     
    { 
      int atomIndex = i*psize + psize - 1- prank;
      if(atomIndex <numAtoms ) 
      {
        Vector3 temp = -atoms[atomIndex].scaledCharge*qrSum;
        Vector3 ri = atoms[atomIndex].position;

//    for (Int i = 0; i < numAtoms; i++)                           // Jc: ZW code
        // atoms[i].force += atoms[i].scaledCharge*qrSum;          // Moldy
/*        if(atoms[atomIndex].atomType ==1)
        {  
           double b    =  R_OP*Bohr_Radius/1.0e-9;  //
           double gama = -b/(0.09572 * cos(0.912109067)); // gama= 0.456840908

           Vector3 p, rod, ro, rd;
           ro = atoms[atomIndex].position; 
           rd.x = Dummy[atoms[atomIndex].molID*3];
           rd.y = Dummy[atoms[atomIndex].molID*3+1];
           rd.z = Dummy[atoms[atomIndex].molID*3+2];
           rod  = rd - ro;
           p    = (rod*temp)/(b*b)*rod;
           p    = gama*(temp - p);
           atoms[atomIndex].force        -= temp -  p;
           atoms[atomIndex - 1].force    -= 0.5* p;
           atoms[atomIndex + 1].force    -= 0.5* p;
*/
/*
          atoms[atomIndex -1].force     -= temp*0.068455;
          atoms[atomIndex].force        -= temp*0.86309;
          atoms[atomIndex +1].force     -= temp*0.068455;   
*/
//        }  
//        else
//        {
	        moli =  &myMols[atoms[atomIndex].molID];
                Atom *atomh1i = moli->myAtoms[0];
                Atom *atomh2i = moli->myAtoms[2];

                if (!strcmp(myMols[atoms[atomIndex].molID].molName, "MCY") && (atoms[atomIndex].atomType == 1) ){
                        atoms[atomIndex].force -= 0.5431676*2.0*atoms[atomIndex].scaledCharge*qrSum;
                        atomh1i->force -= 0.5*0.4568324*2.0*atoms[atomIndex].scaledCharge*qrSum;
                        atomh2i->force -= 0.5*0.4568324*2.0*atoms[atomIndex].scaledCharge*qrSum;
                }
                else {
                        atoms[atomIndex].force -= 2.0*atoms[atomIndex].scaledCharge*qrSum;
                }
          //atoms[atomIndex].force -= 2.0*atoms[atomIndex].scaledCharge*qrSum;        // Gromacs. This may be right
//        } 
//        temp1 -= 2*atoms[atomIndex].scaledCharge*qrSum; // Jc: ZW code
                virial[XX] += temp.x * ri.x;
                virial[XY] += temp.y * ri.x;
                virial[XZ] += temp.z * ri.x;
                virial[YX] += temp.x * ri.y;
                virial[YY] += temp.y * ri.y;
                virial[YZ] += temp.z * ri.y;
                virial[ZX] += temp.x * ri.z;
                virial[ZY] += temp.y * ri.z;
                virial[ZZ] += temp.z * ri.z;
      }
    }
 
//    double temp,temp2;
//    temp2 = temp1.x;
//    MPI::COMM_WORLD.Reduce(&temp2,&temp,1,MPI::DOUBLE,MPI::SUM,0);

}

// JC string EwaldForce::get_force_id() // in order to comply with the cluster compiler
// JC {
// JC  return EwaldForce::forceIdentifier;
// JC}

void EwaldForce::write_force_info(ofstream& of)
{   
    of << "Ewald Force created" << endl;
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

void EwaldForce::write_energy(ofstream& of)
{
    of << "Ewald Energy: " << energy+longEnergy << endl;
    of << "realEnergy: " << realEnergy << '\t' << "longEnergy: " << longEnergy << endl ;
    // of << "introMolEnergy: " << introMolEnergy << '\t' << "surfaceEnergy" << surfDipoleEnergy << endl;
}


