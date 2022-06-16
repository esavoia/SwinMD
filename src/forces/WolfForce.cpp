/** WolfForce.cpp -- 
 **
 ** Copyright (C) 2020
 ** Centre for Molecular Simulation (CMS)
 ** School of Information Technology
 ** Swinburne University of Technology
 ** PO Box 21 Hawthorn, Vic 3122, Australia
 **
 ** Author: Elton
 ** Email: eoyarzua@
 **/

#include "WolfForce.h" 

//#include "erfW.h"
#include <mpi.h>

WolfForce::WolfForce(Ensemble* ensemble) : Force(ensemble) 
{
    DEBUGMSG("Creating WolfForce for electrostatics");

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
    realCut = myEnsemble->myConfig->get_cutoff();
    realCut2 = realCut*realCut;

    stepr = 0.0;
    selfEnergy = 0.0;
    chargedEnergy = 0.0;


    Dummy = new Double[numMols*3];


    // pre-compute self energy correction - Gaussian term
    Double qSum = 0.0;
    Double q2Sum = 0.0;
    for(Int i = 0; i < numAtoms; i++)
    {
        qSum += atoms[i].scaledCharge;
        q2Sum += atoms[i].scaledCharge*atoms[i].scaledCharge;
    }
    selfEnergy = -q2Sum*(alpha/SQRT_PI + 0.5*erfc(alpha * realCut)/realCut);
    if(fabs(qSum) > 1.0e-5)
       chargedEnergy = -PI*qSum*qSum/(2.0*volume*alpha2);
    
    write_force_info(*sysdataFile);
 
   // prank = MPI::COMM_WORLD.Get_rank();  // Jc: get the MPI parameters // Updated to MPI-3 by elton
    //psize = MPI::COMM_WORLD.Get_size();

    MPI_Comm_rank(MPI_COMM_WORLD, &prank);
    MPI_Comm_size(MPI_COMM_WORLD, &psize);

}

WolfForce::~WolfForce()
{
 ;
}

void WolfForce::DummyPosition()

{
   Vector3 Roh1, Roh2, Rhh, Rpo, Dummyi;
   Atom *atom1, *atom2, *atom3;
   Double  RRpo;

   for (int i = 0;i<numMols;i++)
   {
     if (!strcmp(myMols[i].molName, "WAT")) {
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

void WolfForce::compute()
{
    energy = 0.0;
    realEnergy = 0.0;
    realderiv = 0.0;
    real2ndDeriv = 0.0;

    DummyPosition();

    for (Int i = XX; i <= ZZ; i++)
        virial[i] = 0.0;


    real_term(realEnergy, realderiv, real2ndDeriv);                     // Jc: the first term in Ewald



    myEnsemble->myEnergy.realEnergy = realEnergy;
    //myEnsemble->myEnergy.correctEnergy = selfEnergy+chargedEnergy; //Jc: selfEnergy and charged Energy need no reduce

    myEnsemble->myLustig.realDeriv = realderiv;
    myEnsemble->myLustig.realSecDeriv = real2ndDeriv;
 
    for (Int k = XX; k <= ZZ; k++)
        myEnsemble->virial[k] -= virial[k];
}     
 

//Jc:  compute real term in Ewald summation,it is necessary as well as the reciprocal term 
void WolfForce::real_term(Double &realEnergy, double &realderiv, double &real2ndDeriv)
{

    #ifdef DEBUG2
        DEBUGMSG("Computing Wolf Electrostatics with damping erf()");
    #endif

    Atom *atomi, *atomj;
    Vector3 ri, rj, fij, rij, rijm, rim, rjm;
    Int size;
    Double qi, qij;
    
    realEnergy = 0.0;
    realderiv  = 0.0;
    real2ndDeriv = 0.0;

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

        if(!strcmp(myMols[molIDi].molName, "WAT") && atomi->atomType == 1)            // Jc: O atom change to Dummy parameters
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

// For CO2 the COM has to be the same as the Carbon position:
//if (atomi->atomType == 1 ){
//cout << ri << " <- atom      Rij ->   " << rim  << endl;}


        for (int j = 0; j < size; j++)
        {

            atomj = atomi->myPairList[j];
            int molIDj = atomj->molID;

            if (!strcmp(myMols[molIDj].molName, "WAT") && atomj->atomType == 1)        // Jc: O atom change to Dummy parameters
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
//if (atomi->atomType == 1 && atomj->atomType == 1){
//cout << rij << " <- atom      Rij ->   " << endl;
//cout << rijm  << endl;}


            myEnsemble->apply_pbc(rij);
            myEnsemble->apply_pbc(rijm);

            double r2 = rij.length2();

            if (r2 < realCut2)
            {
		    double rijRij = rij*rijm;
                    double r = sqrt(r2);
		    double ea = erfc(alpha*r)/r;
		    double dv = (ea + alpha2PI*exp(-alpha2*r2))/r2;
		    double eR = erfc(alpha*realCut)/realCut; 
		    double Rc = eR/realCut + alpha2PI*exp(-alpha2*realCut2)/realCut;
                    double e = ea - eR;
                    double f = dv - Rc/realCut ;
		    double dv2 = ( 2.0*dv + 2.0*alpha2PI*alpha2*exp(-alpha2*r2) )/r2 - (2.0*Rc/realCut + 2.0*alpha2PI*alpha2*exp(-alpha2*realCut2) )/realCut2;
		    if (!strcmp(myMols[molIDi].molName, "BBV")  && !strcmp(myMols[molIDj].molName, "BBV")) {
			double expBBV = 1.0 + exp(-2.0*(r/Aa0 - 2.0));
                        double invexp = 1.0/expBBV;
                        double Dij2   = invexp*invexp;
                        double Dij8   = Dij2*Dij2*Dij2*Dij2;
                        double Dij15  = Dij8*Dij2*Dij2*Dij2*invexp;
//cout << "Holaaaa" << endl;
			double expRc  = 1.0 + exp(-2.0*(realCut/Aa0 - 2.0));
                        double invRc  = 1.0/expRc;
                        double Dij2Rc = invRc*invRc;
                        double Dij8Rc = Dij2Rc*Dij2Rc*Dij2Rc*Dij2Rc;
                        double Dij15Rc = Dij8Rc*Dij2Rc*Dij2Rc*Dij2Rc*invRc;

                        double dRc    = Rc*Dij15Rc - (30.0/Aa0)*exp(-2.0*(realCut/Aa0 - 2.0))*eR*Dij8Rc*Dij8Rc;
                           e = ea*Dij15 - eR*Dij15Rc ;
                           f = (ea + alpha2PI*exp(-alpha2*r2))/r2*Dij15 - (30.0/Aa0)*exp(-2.0*(r/Aa0 - 2.0))*ea*Dij8*Dij8/r - dRc/realCut;
		    } 
                    realEnergy  += qij*e;
		    realderiv   += -qij*f*rijRij; 
		    real2ndDeriv  += qij*dv2*rijRij*rijRij + 2.0*qij*f*rijRij ; 

                    fij = -qij*f*rij;

// Jc: distribute force among molecule atoms, when dummy charge is calulated
// Jc: the calculated results are not improving. so there is no force distribution 
// Jc: in this calculation
                  atomi->force += fij;
//                if(atomj->atomType == 0)  //Jc: H atoms
//                {
                  atomj->force -= fij;
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
//          atomi->force += fi;
//        } 
     }
   }

   double temp;
   //MPI::COMM_WORLD.Reduce(&realEnergy,&temp,1,MPI::DOUBLE,MPI::SUM,0); Updated to MPI-3 by elton
   MPI_Reduce(&realEnergy,&temp,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   realEnergy = temp;

   temp = 0.0;
   MPI_Reduce(&realderiv,&temp,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   realderiv = temp;

   temp = 0.0;
   MPI_Reduce(&real2ndDeriv,&temp,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   real2ndDeriv = temp;
}

void WolfForce::write_force_info(ofstream& of)
{   
    of << "Wolf Force created" << endl;
    of<<"alpha: "<<alpha<<'\t'<<"real cut off: "<<realCut<<'\t'<<"K space cut off: "<<kCut<<endl;
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

void WolfForce::write_energy(ofstream& of)
{
    of << "Wolf Energy: " << energy+longEnergy << endl;
    of << "realEnergy: " << realEnergy << '\t' << "longEnergy: " << longEnergy << endl ;
    // of << "introMolEnergy: " << introMolEnergy << '\t' << "surfaceEnergy" << surfDipoleEnergy << endl;
}


