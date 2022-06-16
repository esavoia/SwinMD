/** NccForce.cpp -- 
 **
 ** Copyright (C) 2004
 ** Centre for Molecular Simulation (CMS)
 ** School of Information Technology
 ** Swinburne University of Technology
 ** PO Box 21 Hawthorn, Vic 3122, Australia
 **
 ** Author: Jianhui Li
 ** Email: jli@it.swin.edu.au
 **/
//JC this is a child class of Force, at the similar level of the LFForce and Ewarld class
//JC nowthat the unit system is unclear till now, the reasonable method to verify is to assume
//JC all the parameters in the potential in atomic unit, 

//JC Energy unit: Hartree (Eh) is the physical constant = 4.3597438134e-18 J =27.2144eV
//JC Eh=H*H/(Me*A0), where H is the dirac constant; Me is the electron rest mass;A0 is the Bohr radius

//JC Length: Bohr radius, having a value of 5.29189379e-11m, about 53 picometers(half of Angstrom) is used
//JC         in a perturbative expansions of wave function
//JC Electrical Charge; unit is elementary charge, the fundamental physical constant, and the 
//JC having a value of 1.60217646263e-19 C,
//JC 

//JC Mass: taking the rest mass of the electron, 9.10e-31 Kg 

//JC should be at the same level as Zhongwu's calculation and the historical results.29/10/2004 
//JC ******************************************************************************************
//JC ** the current version was changed to the parameters for MCY                            **
//JC ** there is no the vand der Waals uinteraction between the dummy and other positions    **
//JC ******************************************************************************************

#include "NccForce.h" 
#include <mpicxx.h>

NccForce::NccForce(Ensemble* ensemble) : Force(ensemble) //Constructor

{
    #ifdef DEBUG
        DEBUGMSG("Creating NccForce  ");  
    #endif

//Jc: forceIdentifier = "NccForce";

    mols = myEnsemble->molecules;
    numAtoms = myEnsemble->nAtoms;  
    numMols = myEnsemble->nMols;
    params = myEnsemble->myParams; 
    cutOff = myEnsemble->myConfig->get_cutoff();
    cutOff2 = cutOff*cutOff;
//    cout<<" the system number is ,,,,,,  "<<numAtoms<<endl;
    switchDist = cutOff2; // JC:          
    switchOn = myEnsemble->myConfig->is_switch_on();
    CNum = 0;

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

    ofDummy = new ofstream("dummyPos.dat",ios::out);
    ofEnergy = new ofstream("Energy.dat",ios::out);
    ofMolEng = new ofstream("MolEng.dat",ios::out);
    *ofEnergy<<left;
    *ofEnergy<<"3BD           "<<"CHH           "<<"CHD           "<<"CDD           "<<"WHH           "<<"WOO           "\
             <<"WHO           "<<"Total         "<<endl;

    *ofMolEng<<" i   j  "<<"CHH     "<<"CHD       "<<"CDD     "<<"WHH         "<<"WOO         "<<"WHD         "<<\
             "WOD         "<<"WHO          "<<"Total     "<<endl;
    ofpairlist = new ofstream("pairlist.dat",ios::out);
    *ofpairlist<<left; 
    ofminDist  = new ofstream("minDist.dat",ios::out);
    outTimeFile = new ofstream("procsTime.txt",ios::out);
// Jc:    computeCoulomb = myEnsemble->myConfig->compute_coulomb();
// Jc:    if (myEnsemble->myConfig->use_long_correct())
// Jc:        long_range_correct();           // calculate once at the beginning
// Jc:    write_force_info(*sysdataFile);

}

void NccForce::DummyPosition()

{
   Vector3 Roh1, Roh2, Rhh, Rpo, Dummyi;
   Atom *atom1, *atom2, *atom3;
   Double  RRpo;

   for (int i = 0;i<numMols;i++)
   {
    int j = 3 * i;
    atom1 = &atoms[j];
    atom2 = &atoms[j + 1];
    atom3 = &atoms[j + 2];

    Rpo.x = 0.5*(atom3->position.x + atom1->position.x); 
    Rpo.y = 0.5*(atom3->position.y + atom1->position.y);
    Rpo.z = 0.5*(atom3->position.z + atom1->position.z);
    Rpo   = Rpo - atom2->position;
    RRpo  = Rpo.length(); 
    Dummyi = atom2->position + (0.4481389/(RRpo*F_Length))*Rpo; // Jc: unit system change!
    Dummy[j]     = Dummyi.x;
    Dummy[j + 1] = Dummyi.y;
    Dummy[j + 2] = Dummyi.z;

    *ofDummy<<scientific<<setw(15)<<Dummy[i*3] <<setw(15)<< Dummy[i*3 + 1]<<setw(15)<< Dummy[i*3 + 2]<<endl;
   }

   ofDummy->close();      // Jc: closing the stream  
}

// Jc: force calcualtion, but the integration process is done 
void NccForce::seperator()
{
  ;
}

// Jc: in another class named Integrator

void NccForce::compute()  
{

    Atom *atomi,*atomj;
    Atom *atomih1,*atomih2,*atomio, atomid;
    Atom *atomjh1,*atomjh2,*atomjo, atomjd;
    Atom *atomt;
    Atom *atomA,  *atomB,  *atomC;

    Vector3 rij, fij, distance,rtemp,atomPosition;   // JC: instantiate 2 vectors for whatever reasons??
    Vector3 rA, rB, rC;                              // JC: vectors for 3 body calcualtion
    
    double dA1,dA2,dA3,dA4,dA5,dA6;
    double dB1,dB2,dB3,dB4,dB5,dB6;
    double dC1,dC2,dC3,dC4,dC5,dC6;
    double dVdRa,dVdRb,dVdRc;
    double u3 = 0.0;
    double Alpha = 1.0;

    double r, rij2;

    double r_1,rij_2,sr2,sr6,sr12,tmpE,tmpf;


    double B_cutoff  = sqrt(cutOff2) * F_Length;     // Jc: change to the larger cutoff
    double B_cutoff2 = B_cutoff * B_cutoff;
    double T_cutoff  = 0.5*B_cutoff;

// Jc: end of parameterization
// Jc: timing variables
    clock_t start1,end1,start2,end2,start3,start4,end4;

    start1 = clock();

// Jc: save computation time in each processor
    double p_time[100], pp_time[100];

    for(int numProcs=0; numProcs<101; numProcs++)
    {
      p_time[numProcs] = 0.0, pp_time[numProcs] = 0.0;      
    }     
    
// Jc: MPI variables
    int rank,size;
    MPI::Intracomm Intracomm1;

    Vector3 rim, rjm, rijm;
    double  molVirial[9];
    double  temf= 0.0, teme = 0.0;
    double  Ctemp, Wtemp;

// Jc: test MPI
    rank = MPI::COMM_WORLD.Get_rank();
    size = MPI::COMM_WORLD.Get_size();
// Jc: end of test MPI

    Ctemp = 0.0, Wtemp = 0.0;    
    energy = 0.0;
    WPotential = 0.0,CPotential = 0.0;

// Jc: generate the file name each interation!!

    char fileNum[] = "0123456789";
    char fileName[] = "traj0000.dat";
    int firstDigit  = int(CNum / 1000);
    int secondDigit = int((CNum%1000)/100);
    int thirdDigit  = int((CNum%100)/10);
    int fourthDigit = int(CNum%10);
    fileName[4] = fileNum[firstDigit];
    fileName[5] = fileNum[secondDigit];
    fileName[6] = fileNum[thirdDigit];
    fileName[7] = fileNum[fourthDigit];
//    cout<<"test ---------NCCForce ----------" << endl; // JC test the 
    DummyPosition();               // Jc: update every time step to keep the dummy position 

    for (int i = XX; i <= ZZ; i++) // Jc: XX and ZZ Wirial components cf. NEMD_defs.h
    {
        virial[i] = 0.0;           // Jc: reset the virial components
        molVirial[i] = 0.0;       
    }

// JC calculation of assumed negative charges
   double CHH = 0.0; // the coulombic energy between H-H pairs
   double CHD = 0.0; // the coulombic energy between H-Negative pairs
   double CDD = 0.0; // the coulombic energy between Dummy-Dummy pairs
   double WHH = 0.0; // the Van der Waals interaction between H-H pairs
   double WOO = 0.0; // the Van dre Waals interaction between O-O pairs
   double WHD = 0.0; // the Van der Waals interaction between H-D pairs
   double WOD = 0.0; // the Van der Waals interaction between O-D pairs
   double WHO = 0.0; // the Van der Waals interaction between O-D pairs

// temp Potential of Certain Molecule
   double MCHH = 0.0;// the coulombic energy between H-H pairs
   double MCHD = 0.0;// the coulombic energy between H-Negative pairs
   double MCDD = 0.0;// the coulombic energy between Dummy-Dummy pairs
   double MWHH = 0.0;// the Van der Waals interaction between H-H pairs
   double MWOO = 0.0;// the Van dre Waals interaction between O-O pairs
   double MWHD = 0.0;// the Van der Waals interaction between H-D pairs
   double MWOD = 0.0;// the Van der Waals interaction between O-D pairs
   double MWHO = 0.0;// the Van der Waals interaction between O-D pairs
   int flag = 1;     //flag

//   cout<<"the cutoff is ,,,,,,,,,,,____  "<< B_cutoff2<<endl;
   cout<<"the number of atoms in the system is ,,,,,, "<<numAtoms<<" "<<size<<endl;
   int numDecompose = (int) numAtoms/size;
   cout<<" the numDecompose is ...... "<<numDecompose<<endl;
   if ((numAtoms % size) == 0) numDecompose--;

   for(int i = 0; i <= numDecompose;i++)                         // Jc: molID = i

   {
//     int atomIndex = i *size - rank + size -1;    
//     if((i % 2)==0)                                            // Jc: define the first atom: atomA
//     {
       int atomIndex = i *size + rank;
//     }

     if(atomIndex <= (numAtoms -2)) 
     {

       int i3 =  atomIndex; //(i*size + rank);

       atomi = &atoms[i3];

       for(int j = (atomIndex + 1); j<numAtoms;j++) // Jc: molId = j
       { 

         MCHH = 0.0; // Jc: the coulombic energy between H-H pairs
         MCHD = 0.0; // Jc: the coulombic energy between H-Negative pairs
         MCDD = 0.0; // Jc: the coulombic energy between Dummy-Dummy pairs
         MWHH = 0.0; // Jc: the Van der Waals interaction between H-H pairs
         MWOO = 0.0; // Jc: the Van dre Waals interaction between O-O pairs
         MWHD = 0.0; // Jc: the Van der Waals interaction between H-D pairs
         MWOD = 0.0; // Jc: the Van der Waals interaction between O-D pairs
         MWHO = 0.0; // Jc: the Van der Waals interaction between O-D pairs
//         if((atomIndex==0)&&(j==1))
//        {
//         cout<<"the cutoff is ,,,,,,,,,,,____  "<< B_cutoff2<<endl;
//         cout<<"the R(0,1)   is ,,,,,,,,,,,____  "<< rij2<<endl;
//         }

         flag = 1;

//         int delay = 1;
//         start4 = clock();
//         do{
//           end4 = clock();
//           if (end4-start4 >= 10)
//           { delay = 0;}
//         } while(delay == 1);

         int j3 = j;
         atomj = &atoms[j3];
     
         rij = (atomj->position - atomi->position); 
         myEnsemble->apply_pbc(rij);
         rij = F_Length*rij;
         r   = rij.length();
         rij2= r*r;
         double r_cutoff = r/B_cutoff;
         double r4_cutoff4 = rij2*rij2/(B_cutoff2 * B_cutoff2);
         if (rij2 <= (B_cutoff2))
         {
       
           Double sigma2 = 0.36;                    // Jc: for copper
           Double eps4 = 5.479447e-23 * 6.022e+23;  // Jc: 5.479447e-20               
//           cout<<"calcualting the force ........."<<endl;
           r_1 = 1.0/r;
           rij_2 = r_1*r_1;
           sr2 = sigma2*rij_2;
           sr6 = sr2*sr2*sr2;
           sr12 = sr6*sr6;
           tmpE = eps4*(sr12 - sr6);
           tmpf = (sr12 + sr12 - sr6)*6*eps4*rij_2;
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
       }
     }
   }                                             // Jc: end of i-loop in 2-body calculation  


   double energyTerm[6],tempEnergyTerm[6];
   
   energyTerm[0] = CHH;
   energyTerm[1] = CHD;
   energyTerm[2] = CDD;
   energyTerm[3] = WHH;
   energyTerm[4] = WOO;
   energyTerm[5] = WHO;

   MPI::COMM_WORLD.Allreduce(energyTerm,tempEnergyTerm,6,MPI::DOUBLE,MPI::SUM);

   CHH = tempEnergyTerm[0];
   CHD = tempEnergyTerm[1];
   CDD = tempEnergyTerm[2];
   WHH = tempEnergyTerm[3];
   WOO = tempEnergyTerm[4];
   WHO = tempEnergyTerm[5];   
   
   
   end1 = clock();

   p_time[rank] = ((double)(end1-start1))/CLOCKS_PER_SEC;
//   cout<<"the processing time in processor "<<rank<<" is "<<p_time[rank]<<" at time step "<<CNum<<endl;

   MPI::COMM_WORLD.Allreduce(&p_time,&pp_time,size,MPI::DOUBLE,MPI::SUM);
   
   for(int numProcs = 0; numProcs<size;numProcs++)
   {
     p_time[numProcs] = pp_time[numProcs];
   }   

   if(rank == 0)
   { 
     for(int numProcs = 0; numProcs<size;numProcs++){     
     *outTimeFile<<"  "<<numProcs<<"   "<<p_time[numProcs]<<endl;}
   }

// Jc: every processor should wait at end of each time steps , while processor 0 reduce data from others 
// Jc: force, energy components are to be reduced
// Jc: and output on processor 0

// Jc: create force Array 

   start2 = clock();    
//   cout<<"test ---------NCCForce %%%----------" << endl; // JC test the    
/*   double forceArray[5500][3],tempForceArray[5500][3];
   
   for (int numIndex = 0;numIndex < numAtoms; numIndex++)
   {
     forceArray[numIndex][0]     = atoms[numIndex].force.x;
     forceArray[numIndex][1]     = atoms[numIndex].force.y;
     forceArray[numIndex][2]     = atoms[numIndex].force.z;

     tempForceArray[numIndex][0] = atoms[numIndex].force.x; 
     tempForceArray[numIndex][1] = atoms[numIndex].force.y;   
     tempForceArray[numIndex][2] = atoms[numIndex].force.z;     
   }
     
   MPI::COMM_WORLD.Barrier();                                   // Jc: process synchronization 

// Jc: reduce and distribute 3-body energy component

   MPI::COMM_WORLD.Allreduce(forceArray,tempForceArray,4500,MPI::DOUBLE,MPI::SUM);
   
   double temp = 0.0;

   MPI::COMM_WORLD.Allreduce(&u3,&temp,1,MPI::DOUBLE,MPI::SUM);

   u3 = temp;

   for (int numIndex = 0;numIndex < numAtoms; numIndex++)
   {
     forceArray[numIndex]     = tempForceArray[numIndex];

     atoms[numIndex].force.x  = forceArray[numIndex][0];
     atoms[numIndex].force.y  = forceArray[numIndex][1];
     atoms[numIndex].force.z  = forceArray[numIndex][2];
   }
     
   end2 = clock();

   cout<<"the communication time in each iteration is "<<((double)(end2-start2))/CLOCKS_PER_SEC<<" seconds in processor "<<rank<<endl; 
*/
// Jc: output in processor 0
        
   if(rank ==0) {
 
//1      if(CNum == 0)
//2      {
//3       *ofMolEng<<"-------------------total energy of the system --------------------------"<<endl;
//4       *ofMolEng<<"  SUM   "<<"CHH     "<<"CHD       "<<"CDD     "<<"WHH         "<<"WOO         "<<"WHD         "<<\
//5        "WOD         "<<"WHO          "<<"Total     "<<endl;      
//6       *ofMolEng<<"        "<<setw(13)<<CHH/numMols<<setw(13)<<CHD/numMols<<setw(13)<<CDD/numMols<<setw(13)<<\
//7        setw(12)<<WOD<<setw(13)<<WHO<<setw(12)<<CHH+CHD+CDD+WHH+WOO+WHO+WOD+WHD<<endl;             
//8      }

      if((CNum % 2) == 0)
      {
        CHH /= numMols; CHD /= numMols; CDD /= numMols; WHH /= numMols;WOO /=numMols;
        WHD /= numMols; WOD /= numMols; WHO /= numMols; u3  /= numMols;

        *ofEnergy<<scientific<<setw(14)<<u3<<setw(14)<<CHH<<setw(14)<<CHD<<setw(14)<<CDD<<setw(14)<<WHH<< \
          setw(14)<<WOO<<setw(14)<<WHO<<setw(14)<<CHH+CHD+CDD+WHH+WOO+WHO<<endl;

// Jc: test the minimum disrance of first O atom to other O atoms
// Jc: this is to test the structural properties

        atomi = &atoms[1];double temDis = 5.0;
        for (int hj = 0;hj< numAtoms;hj++)
        {
          atomt = &atoms[hj]; 
          rtemp = atomt->position;       
          if((hj >2)&&((atomt->atomType)==1))
          {
             rij = rtemp - atomi->position;
             r   = rij.length();
             if (r<temDis) temDis = r;             
          }
        }

        *ofminDist<<setw(5)<<CNum<<scientific<<setw(15)<<temDis<<endl;
      }

      int printFren = 50;

      if((CNum % printFren) ==0)
      {
        ofCoordinate = new ofstream(fileName,ios::out);

        for (int hj =0;hj< numAtoms;hj++) // Jc: save the cooridnates,forces and velocities of particels 
        {
          atomt = &atoms[hj]; 
          rtemp = atomt->position; rij = atomt->force;       
          *ofCoordinate<<fixed<<setw(5)<<hj<<setw(10)<<rtemp.x<<setw(10)<<rtemp.y<<setw(10)<<rtemp.z<<\
           setw(14)<<rij.x<<setw(14)<<rij.y<<setw(14)<<rij.z<<endl;
        }

      }
      
   }                                       // Jc: end of data output in process 0
 

   CNum++;

   for (Int k = XX; k <= ZZ; k++)
   {
      myEnsemble->virial[k] -= virial[k];
        // myEnsemble->molVirial[k] -= molVirial[k];
   }

}

// Beyond the cutoff, usually the repulsive term can be neglected.
// long_range_correct() will make correction to dispersion term 
// algorithm refering to GROMACS Appendix C

void NccForce::long_range_correct()
{
    Int numTypes = myEnsemble->nAtomTypes;
    Int *atomNum = new Int[numTypes];
    Double val = 0.0;

    for (Int i = 0; i < numTypes; i++)
        atomNum[i] = 0;

    // count atom number for each type
    for (Int i = 0; i < numAtoms; i++)
    {   
        int type = atoms[i].atomType;
        atomNum[type] += 1;
    }

    for (Int i = 0; i < numTypes; i++)
    {
        for (Int j = 0; j < numTypes; j++)
        {
            LJPairParam ljParam = params->get_lj_parameter(i, j);
            Double sigma2 = ljParam.sigma;
            Double eps4 = ljParam.eps;

            val += atomNum[i]*atomNum[j]*eps4*sigma2*sigma2*sigma2;
        }
    }

    eLrc = -2*PI*val/(3*myEnsemble->boxLx*myEnsemble->boxLy*myEnsemble->boxLz*cutOff*cutOff*cutOff);
    vLrc = -eLrc;          
}

void NccForce::write_force_info(ofstream& of)
{
    of << "Ncc Force" << endl;
    of << "cutoff2: " << cutOff2 << endl;
    if (switchOn)
        of << "compute LJ force with switch function" << endl;
    if (computeCoulomb)
        of << "Non-bonded force include coulomb" << endl;
    of << "Nc parameters:" << endl;
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

void NccForce::write_energy(ofstream& of) // JC designed by Jianhui Li to store the 
{
      of << WPotential <<'\t' << CPotential <<endl; 
      WPotential = 0.0;
      CPotential = 0.0;
}

//JC string NccForce::get_force_id() 
//JC{
//Jc  return NccForce::forceIdentifier;
//JC}

