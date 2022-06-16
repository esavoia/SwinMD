/** LFIntegrator.cpp -- 
 **
 ** Copyright (C) 2003
 ** Centre for Molecular Simulation (CMS)
 ** School of Information Technology
 ** Swinburne University of Technology
 ** PO Box 21 Hawthorn, Vic 3122, Australia
 **
 ** Author: Zhongwu Zhou
 ** Email: zzhou@it.swin.edu.au
 ** Commented and modified by Jianhui LI in 2004/jli@it.swin.edu.au
 **/
// JC accordingly change have been implemented in the new Class LJIntegrator1
// JC in the primary step, there is no LJ force calculation and Ewald, long distance correction
// JC long distance correction were done in the LJforce class, see Zhongwu's code
#include "LFIntegrator.h"
#include <time.h>
#include "mpi.h"


// JC: constructor but inherit from the class Integrator, note that they have the same argu
LFIntegrator :: LFIntegrator(Ensemble* ensemble, SimConfiguration* simConfig) : Integrator(ensemble, simConfig)
{
    #ifdef DEBUG
        DEBUGMSG("Creating LFIntegrator");
    #endif
    a = numAtoms;
    a3 = 3*a; 
    E2Sum = 0.0; ESum = 0.0; // Jc: give the initial values E2Av = <E2>, EAv = <E>
    aTSum = 0.0; mTSum = 0.0; aT2Sum = 0.0; mT2Sum = 0.0;
    numUpdates = 0;
//    prank = MPI::COMM_WORLD.Get_rank();  // Jc: get the MPI parameters //  Updated to MPI-3 elton
//    psize = MPI::COMM_WORLD.Get_size();
    MPI_Comm_rank(MPI_COMM_WORLD, &prank);
    MPI_Comm_size(MPI_COMM_WORLD, &psize);

        
    if(prank ==0)   outTimeFile = new ofstream("proForce.out",ios::out);
    if(prank ==0) 
    {
        ofPotential = NULL;
        ofPotential = new ofstream("potential.out",ios::out);
        *ofPotential<<"##             atomTemp "<<"      molTemp  "<<"       Pressure "<<"   molCor"<<\
         " sCorr"<<"    corEn"<<"   induc"<<"   3BDEn"<<"   totEn"<<"        Cv"<<endl;
	ofNH = new ofstream("frictionNH.out",ios::out);
	*ofNH << "## step" << "      " << "friction" << "    CouplingConst "<< "    Q_nh" << endl;
    } 
    //Jc: reduce the output frequency
    initialise();
}

LFIntegrator :: ~LFIntegrator()
{
    if (mass_r != NULL)     delete mass_r;
    if (vtPlus != NULL)     delete vtPlus;
    if (vtMinus != NULL)    delete vtMinus;
    if (rOld != NULL)       delete rOld;
    if (virtfor != NULL)    delete virtfor;
    if (virt1 != NULL)      delete virt1;
    if (virt2 != NULL)      delete virt2;
    if (rNew != NULL)       delete rNew;
    if (moving != NULL)     delete moving;
    if (moved  != NULL)     delete moved;
    if (bondLength != NULL) delete bondLength;

    rOld = NULL;
    virtfor = NULL;
    virt1 = NULL;
    virt2 = NULL;
    rNew = NULL;
    moving = NULL;
    moved  = NULL;
    bondLength = NULL;
    mass_r = NULL;
    vtPlus = NULL;
    vtMinus = NULL;   
}

void LFIntegrator :: initialise()
{
    #ifdef DEBUG
        DEBUGMSG("initialise LFIntegrator");// JC 54 message from Debug
    #endif
    
    Integrator :: initialise(); 

    maxCount = (int)mySimConfig->get_max_count();      //10000;
    tol = mySimConfig->get_tolerance();                //0.00001;
    totsteps = mySimConfig->get_n_steps();

    mass_r = new Double[numAtoms];
    vtPlus = new Vector3[numAtoms];
    vtMinus = new Vector3[numAtoms];
    if ((mass_r == NULL) || (vtPlus == NULL) || (vtMinus == NULL))
        ERRORMSG("memory allocation error for LFIntegrator");
    for (int i = 0; i < numAtoms; i++)
        mass_r[i] = 1.0/myAtoms[i].mass;

    rOld = NULL;
    virtfor  = NULL;
    virt1 = NULL;
    virt2 = NULL;
    rNew = NULL;
    moving = NULL;
    moved  = NULL;
    bondLength = NULL;
    if (has_constr_force) 
    {
        #ifdef DEBUG
            DEBUGMSG("Constraint with SHAKE");
        #endif 

        rOld = new Vector3* [numMols];
	virtfor = new Vector3* [numMols];
        virt1 = new Vector3* [numMols];
        virt2 = new Vector3* [numMols];
        rNew = new Vector3* [numMols];
        moving = new bool* [numMols];
        moved  = new bool* [numMols];
        bondLength = new Double* [numMols];

        if ((rOld==NULL)||(virt1==NULL)||(virt2==NULL)||(virtfor==NULL)||(rNew==NULL)||(moving==NULL)||(moved==NULL)||(bondLength==NULL))
            ERRORMSG("memory allocation error for LFIntegrator");
        for (int nm = 0; nm < numMols; nm++)
        {
            int ns = myMols[nm].numAtoms;
            if (ns <= 1)
                continue;
            rOld[nm] = new Vector3 [ns];
	    virtfor[nm]  = new Vector3 [ns];
            virt1[nm] = new Vector3 [ns];
            virt2[nm] = new Vector3 [ns];
            rNew[nm] = new Vector3 [ns];
            moving[nm] = new bool [ns];
            moved[nm]  = new bool [ns];
            bondLength[nm] = new Double [ns];    //Somehow this is no necessary because the lenght can be obtained from params once, not for every pair 

            if((rOld[nm]==NULL)||(virt1[nm]==NULL)||(virt2[nm]==NULL)||(virtfor[nm]==NULL)||(rNew[nm]==NULL)||(moving[nm]==NULL)||(moved[nm]==NULL)||(bondLength[nm]==NULL))
                ERRORMSG("memory allocation error for LFIntegrator");           
        }
        // set bond length, we only consider cyclic bond but not include angle

        for (int nm = 0; nm < numMols; nm++)
        {
            int ns = myMols[nm].numAtoms;
            if (ns <= 1)
                continue;
            //for(int site = 0; site < (ns-1); site++)
            //{
            //    int a1 = myMols[nm].myAtoms[site]->atomID;
            //    int a2 = myMols[nm].myAtoms[site+1]->atomID;
            //    Double l1 = myEnsemble->myParams->get_bond_length(a1, a2);
            //    bondLength[nm][site] = l1*l1;                
            //}
            // the cyclic bond, only consider Water here
            if ( ( !strcmp(myMols[nm].molName, "WAT") ) || ( !strcmp(myMols[nm].molName, "MCY") ) )
            { 
                int a1 = myMols[nm].myAtoms[0]->atomID;
                int a2 = myMols[nm].myAtoms[1]->atomID;
                int a3 = myMols[nm].myAtoms[2]->atomID;
		Double l1 = myEnsemble->myParams->get_bond_length(a1, a2);
                bondLength[nm][0] = l1*l1;
	        Double l2 = myEnsemble->myParams->get_bond_length(a2, a3);
                bondLength[nm][1] = l2*l2;
                Double l3 = myEnsemble->myParams->get_bond_length(a1, a2, a3);
                bondLength[nm][2] = l3*l3;
            }
	    if (!strcmp(myMols[nm].molName, "CO2"))
            {
                int a1 = myMols[nm].myAtoms[0]->atomID;
                int a2 = myMols[nm].myAtoms[1]->atomID;
                Double l1 = myEnsemble->myParams->get_bond_length(a1, a2);
                bondLength[nm][0] = l1*l1;				// Check this matrix when water is included  
            }
	    if (!strcmp(myMols[nm].molName, "DAM"))
            {
                int a1 = myMols[nm].myAtoms[0]->atomID;
                int a2 = myMols[nm].myAtoms[1]->atomID;
                Double l1 = myEnsemble->myParams->get_bond_length(a1, a2);
                bondLength[nm][0] = l1*l1;				// Check this matrix when water is included  
            }

   	   //cout << " nmol  = " << nm  << "  bondLength = " << bondLength[nm][2] << endl;
        }
    }
    
    if (mySimConfig->is_new_start())
       start();
    else
    {
       read_vtMinus();
          
    }


    if (!strcmp(thermotype, "NoseHoover")){
            DEBUGMSG("Apply Nose-Hoover Thermostat");
    }
    if (!strcmp(thermotype, "Berendsen")){
            DEBUGMSG("Apply Berendsen Thermostat");
    }
    if (!strcmp(thermotype, "Andersen")){
            DEBUGMSG("Apply Andersen Thermostat");
    }
        char*  strm = myEnsemble->myParams->molTypes[0].name;
        if (!strcmp(strm, "BBV"))
        {
                double Imol;
                const double Omass = myEnsemble->myParams->atomParams[0].mass;
                const double Cmass = myEnsemble->myParams->atomParams[1].mass;
                const double Mbbv = 2*Omass + Cmass;
                co2bond = myEnsemble->myParams->bondParams[0].r0;

                Imol = 2.0*Omass*(co2bond)*(co2bond);           // Change here if the masses are changed at sysData.in

                Lvir = sqrt(4.0*Imol/Mbbv);
                llra = 0.5*(2.0*co2bond - Lvir);
                Lvir = 0.1979944;
                llra = 0.0171028;
                cout << " CO2_bondlength = "<< co2bond << endl;
                cout << " Increment =  " << llra << endl;
        }

}

void LFIntegrator :: VirtualPositions()
{
   Atom* atom;
   Molecule *mol;

   const double ll = 2.0*co2bond;
   const double cc = llra;
//   const double ll = 0.2298;
//   const double cc = 0.016926032;
   const double bb = ll - cc;
//cout << myEnsemble->myParams->bondParams[0].r0 <<endl;
   for (int i = 0; i<numMols; i++)
   {
     if (!strcmp(myMols[i].molName, "BBV")) {

        mol = &myMols[i];
        Atom *atomO0 = mol->myAtoms[0];
        Atom *atomC2 = mol->myAtoms[2];
        Atom *atomO1 = mol->myAtoms[4];

     virt1[i][0]  = (cc/ll)*atomO1->position + (bb/ll)*atomO0->position;
     virt1[i][1]  = (bb/ll)*atomO1->position + (cc/ll)*atomO0->position;
     }
   }
}

void LFIntegrator :: run(Int numTimeSteps)// this is called by the Simulation.run()
{
    Double shakeVirial[9];            
    bool doUpdate;
    Vector3 position, diff;
    Molecule* mol;
    Atom* atom;
    double vomx, vomy, vomz, masstott;
    int ns, NHcount;
    double qq, invqq, ww, invww, boxLen, density;
    double forceArray[a][3],tempForceArray[a][3];
 
    clock_t start1,end1,start2,end2,start3,startlong,endlong;  //Jc: timing

    double p_time, pp_time[100];

    for(int numProcs=0; numProcs<101; numProcs++)
    {
      pp_time[numProcs] = 0.0;      
    }     

    #ifdef DEBUGXX 
        if (dumpFile == NULL)   dumpFile = "dump.out";
        ofstream of = ofstream(dumpFile, ios::out);
    #endif 

    for (Int i = 0; i < numTimeSteps; i++)
    {
        atomKinEnergy = 0.0;
        molKinEnergy = 0.0;
        Rm.x = 0.0;
        Rm.y = 0.0;
        Rm.z = 0.0;
        Mtot = 0.0;

        currentSteps++;

/*        if(prank==0)   //Jc: for testing the continue simualtion
        {
          if(mySimConfig->is_new_start())
          { 
             if(currentSteps ==101)
             myEnsemble->write_position(currentSteps); 
          }
          else
          {
           if(currentSteps == 1)
             myEnsemble->write_position(currentSteps);
          } 
        } 
*/
        if((prank == 0)&&((currentSteps%20)==0)) 
            cerr << "integration steps: " << currentSteps << endl;
                               
        // check if need to update pairlist
        doUpdate = false;
        for (Int j = 0; j < numAtoms; j++)
        {
            if(myAtoms[j].displacement.length2() > drMax)
            {
                doUpdate = true;
                break;
            }
        }

        if(doUpdate)
        {
            numUpdates++;
            for (Int j = 0; j < numAtoms; j++)
                myAtoms[j].displacement = 0.0; 
            
            #ifdef DEBUGXX //JC blocked by Jianhui Li
            // trace mol & atom position before apply_pbc
                of << "integration steps: " << currentSteps << "  updating list for number: " << numUpdates<< endl;
                write_bond(of);    
            #endif

            // apply PBC to molecules
            for (Int nm = 0; nm < numMols; nm++)
            {  
                mol = &myMols[nm];

                position = mol->massCenter; 
                myEnsemble->apply_pbc(mol->massCenter);        
                diff = mol->massCenter - position;
                for (Int j = 0; j < mol->numAtoms; j++)
                    mol->myAtoms[j]->position += diff;
            }
            myEnsemble->update_pairlist();                                   
        }

        // update box if const. pressure
            
        // zero force & virial befor force computation
        for (Int j = XX; j <= ZZ; j++)
        {
            myEnsemble->virial[j] = 0.0;
            myEnsemble->molVirial[j] = 0.0;
            myEnsemble->atomPressure[j] = 0.0;
            myEnsemble->molPressure[j] = 0.0;
        }
	
	vomx = 0.0;
	vomy = 0.0;
	vomz = 0.0;
	masstott = 0.0;
// Jc: reset the Atom force before each iteration
        for (Int j = 0; j < numAtoms; j++){
            myAtoms[j].force = 0.0;
	    vomx += vtMinus[j].x*myAtoms[j].mass;
	    vomy += vtMinus[j].y*myAtoms[j].mass;
	    vomz += vtMinus[j].z*myAtoms[j].mass;
	    masstott += myAtoms[j].mass; 
	}
	vomx = vomx/masstott;
	vomy = vomy/masstott;
	vomz = vomz/masstott;

	if (fabs(vomx) > 1e-11 || fabs(vomy) > 1e-11 || fabs(vomz) > 1e-11){
                for (Int j = 0; j < numAtoms; j++){
                                vtMinus[j].x -=  vomx;
                                vtMinus[j].y -=  vomy;
                                vtMinus[j].z -=  vomz;
				myAtoms[j].velocity.x -= vomx;
				myAtoms[j].velocity.y -= vomy;
				myAtoms[j].velocity.z -= vomz;
                }
         }


        double tensor[9],tensorTemp[9];

        // compute forces at time t+dt
        startlong = clock();
        start1 = clock();

        if (computeBuff) {
            buffForce->compute();
        }
        else {
            if(prank ==0)
            {
                ljForce->compute();
                // mcyForce->compute();
                //mieForce->compute();
                // saapForce->compute();
	        //buckinghamForce->compute();
	        //bbvForce->compute();
            }
            //nccForce->compute();    //Jc: to integrate the MCY force into the simulation
	    double tmp = 0.0;
	    tmp = myEnsemble->myLustig.ljdUlus;
	    MPI_Bcast(&tmp,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	    myEnsemble->myLustig.ljdUlus = tmp;
	}
        endlong   = clock();

//cout <<" r = " << prank << "  ljDudv = " << myEnsemble->myLustig.ljdUlus << endl; 

        if (computeEwald)                             
        {
            startlong = clock();              
            ewaldForce->compute(); 
            endlong   = clock();
                  
            if ((i%longSteps) == 0)
            {   
              startlong = clock();      
              ewaldForce->compute_long();
              endlong   = clock();
            }

// Jc: commented for 2-body properties calcualtion
// Jc: reduce virial calculated EwaldForce respectively

            for (int tensorIndex = XX;tensorIndex <= ZZ; tensorIndex++)
            {
              tensor[tensorIndex]     = myEnsemble->longVirial[tensorIndex];
              tensorTemp[tensorIndex] = myEnsemble->longVirial[tensorIndex];
            }
     
//            MPI::COMM_WORLD.Allreduce(tensor,tensorTemp,9,MPI::DOUBLE,MPI::SUM);
            MPI_Allreduce(tensor,tensorTemp,9,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

            for (int tensorIndex =XX;tensorIndex <= ZZ; tensorIndex++)
            {
              tensor[tensorIndex]                  = tensorTemp[tensorIndex];
              myEnsemble->longVirial[tensorIndex]  = tensor[tensorIndex];
            } 
        }
        else if (computeMultipole)                             
        {
            startlong = clock();              
            multiElec->compute(); 
            endlong   = clock();
                  
            startlong = clock();      
            multiElec->compute_long();
            endlong   = clock();

            if (computeInduction)                             
            {
                inductionForce->compute(); 
                  
                inductionForce->compute_long();
	    }
// Jc: reduce virial calculated EwaldForce respectively

            for (int tensorIndex = XX;tensorIndex <= ZZ; tensorIndex++)
            {
              tensor[tensorIndex]     = myEnsemble->longVirial[tensorIndex];
              tensorTemp[tensorIndex] = myEnsemble->longVirial[tensorIndex];
            }
     
            MPI_Allreduce(tensor,tensorTemp,9,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

            for (int tensorIndex =XX;tensorIndex <= ZZ; tensorIndex++)
            {
              tensor[tensorIndex]                  = tensorTemp[tensorIndex];
              myEnsemble->longVirial[tensorIndex]  = tensor[tensorIndex];
            } 
        }
        else if (computeWolf)                             
        {
	    wolfForce->compute();
	}
// Jc: force reduce
   
         for (int numIndex = 0;numIndex < numAtoms; numIndex++)
         {
            forceArray[numIndex][0]     = myAtoms[numIndex].force.x;
            forceArray[numIndex][1]     = myAtoms[numIndex].force.y;
            forceArray[numIndex][2]     = myAtoms[numIndex].force.z;
  
         }
     
//         MPI::COMM_WORLD.Allreduce(forceArray,tempForceArray,a3,MPI::DOUBLE,MPI::SUM);
         MPI_Allreduce(forceArray,tempForceArray,a3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

//         MPI::COMM_WORLD.Barrier();
         MPI_Barrier(MPI_COMM_WORLD);

         for (int numIndex = 0;numIndex < numAtoms; numIndex++)
         {
            myAtoms[numIndex].force.x  = tempForceArray[numIndex][0];
            myAtoms[numIndex].force.y  = tempForceArray[numIndex][1];
            myAtoms[numIndex].force.z  = tempForceArray[numIndex][2];
         } 

// Jc: reduce virial calculated in NccForce
         for (int tensorIndex = XX;tensorIndex <= ZZ; tensorIndex++)
         {
            tensor[tensorIndex]     = myEnsemble->virial[tensorIndex];
            tensorTemp[tensorIndex] = myEnsemble->virial[tensorIndex];
         }

//         MPI::COMM_WORLD.Allreduce(tensor,tensorTemp,9,MPI::DOUBLE,MPI::SUM);
         MPI_Allreduce(tensor,tensorTemp,9,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

         for (int tensorIndex =XX;tensorIndex <= ZZ; tensorIndex++)
         {
            tensor[tensorIndex]              = tensorTemp[tensorIndex];
            myEnsemble->virial[tensorIndex]  = tensor[tensorIndex];
         }
 
/*           end1 = clock();
            p_time = ((double)(end1-start1))/CLOCKS_PER_SEC;

            MPI::COMM_WORLD.Gather(&p_time,1,MPI::DOUBLE,&pp_time,1,MPI::DOUBLE,0);   

            if(prank == 0)
            { 
              for(int numProcs = 0; numProcs<psize;numProcs++){     
              *outTimeFile<<"  "<<numProcs<<"   "<<pp_time[numProcs]<<endl;}
            }
*/

            // Jc: reduce the output frequency
            // sum long force and virial to the total.  // sum directly in the Ewald reciprocal
//Jc:            for (Int j = 0; j < numAtoms; j++)
//Jc:                myAtoms[j].force += myAtoms[j].longForce;
            for (Int j = XX; j <= ZZ; j++)
                myEnsemble->virial[j] += myEnsemble->longVirial[j];

        if(has_bond_force)
        {
          //cout<<".......calcualting bond force ............."<<endl;         
          bondForce->compute();
        }
        if(has_angle_force)
        {
          //cout<<".......calcualting angle force ............."<<endl; 
          angleForce->compute();
        }
        if (use_ub_harmonic)
        {
          ubForce->compute();
        }

	myEnsemble->myEnergy.inducEnergy = myEnsemble->myEnergy.realPolEnergy + myEnsemble->myEnergy.longPolEnergy \
	  +  myEnsemble->myEnergy.correctPolEnergy + myEnsemble->myEnergy.molCorrectPolEnergy + myEnsemble->myEnergy.surfacePolEnergy;
        myEnsemble->myEnergy.totEnergy =  myEnsemble->myEnergy.ljEnergy + myEnsemble->myEnergy.realEnergy\
          +  myEnsemble->myEnergy.longEnergy + myEnsemble->myEnergy.molCorrectEnergy     \
          +  myEnsemble->myEnergy.surfCorrectEnergy + myEnsemble->myEnergy.correctEnergy \
	  +  myEnsemble->myEnergy.angleEnergy + myEnsemble->myEnergy.bondEnergy +  myEnsemble->myEnergy.ubEnergy   \
          +  myEnsemble->myEnergy.inducEnergy + myEnsemble->myEnergy.TBEnergy;

// Jc: compute new position and velocity
//cout << " Mol types  =   " << myEnsemble->nMolTypes << endl;
        if (has_constr_force)
        {
	  // Compute here Rm for molecular systems
	  //
	    for (int i = 0; i < myEnsemble->nMolTypes; i++)
 	     {
		char*  str = myEnsemble->myParams->molTypes[i].name;
//cout << " moltype =  "  << str <<endl;
	        if ( (!strcmp(str, "WAT") )  ||  ( !strcmp(str, "MCY") )  )
                 {
                	shake(tol, maxCount, shakeVirial,currentSteps);
	         }
	        if (!strcmp(str, "DAM"))
                 {
                	shake(tol, maxCount, shakeVirial,currentSteps);
	         }
         	if (!strcmp(str, "CO2"))
           	 {
                	shakeCFR(tol, maxCount, shakeVirial,currentSteps);
           	 }
		if (!strcmp(str, "BBV"))
                 {
                        shakeBBV(tol, maxCount, shakeVirial,currentSteps);
                 }
	     }
        }               
        else {
	 for (Int i = 0; i < numAtoms; i++)
         {
            Atom* atom = &myAtoms[i];
	    Vector3 velvect;
	    double velcheck = 1.0;
	    NHcount = 0;
	    if (!strcmp(thermotype, "NoseHoover")/* && (currentSteps <= 0.25*totsteps)*/ ){
		while ( (velcheck > 1E-10) && (NHcount < 11) ){
            	    vtPlus[i] = vtMinus[i] + timeStep*(mass_r[i]*atom->force - fnh*atom->velocity); // NH eqns   
		    if (NpT){
            	        vtPlus[i] +=  - timeStep*nnh*atom->velocity;   
		    }  
	            velvect  = (vtPlus[i] + vtMinus[i])*0.5  - atom->velocity;
		    velcheck = velvect.length();
		    atom->velocity = (vtPlus[i] + vtMinus[i])*0.5;
	        //if (prank==0){ cout << " NHcount =  "  << NHcount  << "      velcheck =  "<<  velcheck <<  "   velvect.x =  "<< velvect.x << endl;}
		    NHcount++;
		}
		if (NHcount == 10){
		    ERRORMSG("NH iteration Error");
		}
		
	    } else {
	    	vtPlus[i] = vtMinus[i] + timeStep*mass_r[i]*atom->force;
	    }

            Vector3 dr = timeStep*vtPlus[i]; // JC displacement in one Time step

	    velcheck = 1.0;
	    NHcount = 0;

            if (NpT){
                Vector3 rplus = atom->position;
                Vector3 rhalfp = atom->position;
		while ( (velcheck > 1E-10) && (NHcount < 11) ){
		    rplus = atom->position + timeStep*(vtPlus[i] + nnhpl*(rhalfp - Rm));
	            velvect  = 0.5*(atom->position + rplus) - rhalfp;
		    velcheck = velvect.length();
		    rhalfp = 0.5*(atom->position + rplus);
		//if(prank==0){
		//cout << " NHcount =  "  << NHcount  << "      velcheck =  "<<  velcheck <<  "   velvect.x =  "<< velvect.x << endl;}
		    NHcount++;
		}
		if (NHcount == 10){
		    ERRORMSG("NpT barostat integration Error");
		}

		dr = rplus - atom->position;
	    }
            atom->position += dr;
            atom->realPos += dr;             // JC i didn't find any member in the Atom class named realPos
            atom->displacement += dr;
	    atom->velocity = (vtPlus[i] + vtMinus[i])*0.5;
            vtMinus[i] = vtPlus[i];
         }
	}

	// -------------------------------------------------------------------------
        //  Now compute the molecular centre, momentum, & kinetic energy
        //  with the updated positions and velocities at the corresponding 
        //  timestep t.
        // -------------------------------------------------------------------------
        for (Int nm = 0; nm < numMols; nm++)
        {  
            mol = &myMols[nm];

            mol->momenta = 0.0;
            mol->massCenter = 0.0;
	    if (!strcmp(mol->molName, "BBV")) {
                Atom *atom0  = mol->myAtoms[0];
                Atom *atomB1 = mol->myAtoms[1];
                Atom *atomc  = mol->myAtoms[2];
                Atom *atomB2 = mol->myAtoms[3];
                Atom *atom1  = mol->myAtoms[4];

                atom0->momentum = 0.5*(atom0->velocity)*mol->mass;
                atom1->momentum = 0.5*(atom1->velocity)*mol->mass;
                atomc->momentum = 0.0;
                atomB1->momentum = 0.0;
                atomB2->momentum = 0.0;

                mol->momenta = atom0->momentum + atom1->momentum;
                mol->massCenter = atom0->mass*atom0->position + atom1->mass*atom1->position + atomc->mass*atomc->position + atomB1->mass*atomB1->position +  atomB2->mass*atomB2->position;
                atomKinEnergy += 0.5*atom0->velocity*atom0->velocity*mol->mass + 0.5*atom1->velocity*atom1->velocity*mol->mass; // The 0.5 is because half the tot mass
//cout << "   posC.x = " << atomc->position.x << "    posC.y = " << atomc->position.y << endl;
            }
            else{          
            for (Int j = 0; j < mol->numAtoms; j++)
            {
                atom = mol->myAtoms[j];
                atom->momentum = atom->velocity*atom->mass;
                mol->momenta += atom->momentum;
                mol->massCenter += atom->mass*atom->position;
                atomKinEnergy += atom->velocity*atom->velocity*atom->mass; // JC compute the Kinetic energy
	    }  //if (j==1){cout << "   posC.x = " << atom->position.x << "    posC.y = " << atom->position.y << endl;}
            }            
            mol->massCenter /= mol->mass;
            molKinEnergy += mol->momenta*mol->momenta/mol->mass;
        }

	if (!strcmp(thermotype, "NoseHoover")){
		if(useAtomThermo){
			qq = atomNDF*BOLTZMAN*temperature*couplconst*couplconst;
		}
		else{
			qq = molNDF*BOLTZMAN*temperature*couplconst*couplconst;
		}

		invqq = 1/qq;

		//if (currentSteps <= 0.5*totsteps)
		{
			if(useAtomThermo){
			nht = nht + 1;
			fnhpl = fnhmi + timeStep*invqq*(atomKinEnergy - atomNDF*BOLTZMAN*temperature);
		}
			else{
			fnhpl = fnhmi + timeStep*invqq*(molKinEnergy - molNDF*BOLTZMAN*temperature);
			nht = nht + 1;
			}
		    fnh = 0.5*(fnhpl + fnhmi);
		    fnhmi = fnhpl;	
		}

		if(((currentSteps%100)==0)&&(prank==0))
        	{
        	*ofNH << fixed << setw(6) << setprecision(4) << currentSteps; 
		*ofNH << fixed << setprecision(6) << setw(14) << fnh << fixed << setprecision(6) << setw(14) << couplconst;
		*ofNH << fixed << setprecision(6) << setw(14) << qq;
		if (!NpT){
		    *ofNH << endl;
		}
		}
	}
	else { fnh = 0.0;}

        myEnsemble->myEnergy.atomKinEnergy = 0.5*atomKinEnergy;
        myEnsemble->myEnergy.molKinEnergy = 0.5*molKinEnergy;
        molTemp = molKinEnergy/(molNDF*BOLTZMAN);
        atomTemp = atomKinEnergy/(atomNDF*BOLTZMAN); //JC compute the atomic temperature
        myEnsemble->atomTemp = atomTemp;
        myEnsemble->molTemp = molTemp;

	molVirt   =  myEnsemble->myLustig.ljdUlus + myEnsemble->myLustig.inducDeriv - 3.0*myEnsemble->myEnergy.TBEnergy/volume \ 
		   + myEnsemble->myLustig.surfDeriv + myEnsemble->myLustig.molDeriv + myEnsemble->myLustig.surfPolDeriv +  myEnsemble->myLustig.molPolDeriv \
		   + (myEnsemble->myLustig.longDeriv + myEnsemble->myLustig.realDeriv + myEnsemble->myLustig.ubDeriv + myEnsemble->myLustig.BondDeriv \
		   + myEnsemble->myLustig.longPolDeriv + myEnsemble->myLustig.realPolDeriv + myEnsemble->myLustig.selfPolDeriv)/(3.0*volume) ;
 
        if (has_constr_force) {
	     molPt     = molKinEnergy/(3.0*volume) - molVirt;
	}
	else {
	     molPt     = atomKinEnergy/(3.0*volume) - molVirt;
	}

	myEnsemble->myAvLus.P = molPt;

vomx = 0.0;
vomy = 0.0;
vomz = 0.0;
masstott = 0.0;


	if (!strcmp(thermotype, "Andersen")){
 	 	//if ((currentSteps%40)==0) {
 	 	double prob = urandf();
	        MPI_Bcast(&prob,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    		MPI_Barrier(MPI_COMM_WORLD);
		
 	 	if ( prob < couplconst*timeStep ) {
		 nht = nht + 1;
		if (prank == 0) {
		    cout << "Andersen Massive collisions:  " << nht << "   P " << couplconst*timeStep << "      step = " << currentSteps << endl; }
		for (Int i = 0; i < numAtoms; i++)
        	{
            		Atom* atom = &myAtoms[i];
			qq = sqrt(BOLTZMAN*temperature/atom->mass); 

			double vtemp[3];
			vtemp[0] = gaussVel(qq); 
			vtemp[1] = gaussVel(qq); 
			vtemp[2] = gaussVel(qq); 
	        
			MPI_Bcast(&vtemp,3,MPI_DOUBLE,0,MPI_COMM_WORLD);
    			MPI_Barrier(MPI_COMM_WORLD);

	    		vtPlus[i].x = vtemp[0];
	    		vtPlus[i].y = vtemp[1];
	    		vtPlus[i].z = vtemp[2];

            		//atom->velocity = vtPlus[i];
            		vtMinus[i] = vtPlus[i];
       	 	}
		for (Int j = 0; j < numAtoms; j++){
                        Atom* atom = &myAtoms[j];
                        vomx += vtMinus[j].x*atom->mass;
                        vomy += vtMinus[j].y*atom->mass;
                        vomz += vtMinus[j].z*atom->mass;
                        masstott += atom->mass;
                }

                vomx = vomx/masstott;
                vomy = vomy/masstott;
                vomz = vomz/masstott;

                for (Int j = 0; j < numAtoms; j++){
                        vtMinus[j].x -=  vomx;
                        vtMinus[j].y -=  vomy;
                        vtMinus[j].z -=  vomz;
                       // atom->velocity.x -= vomx;
                       // atom->velocity.y -= vomy;
                      //  atom->velocity.z -= vomz;
                }
		}
	}

        
        computePressure();                           // JC where's the computePressure(Integrator.cpp)

        if (has_constr_force)                        // correct atom virial and pressure due to shake virial
        {
           for (Int i = XX; i <= ZZ; i++)
            {
                myEnsemble->virial[i] += shakeVirial[i];
                myEnsemble->atomPressure[i] += (shakeVirial[i]/volume);
            }
        }
 
	// virial for Lustig:
	//

	Double molCorrectLustig = 0.0;

    for (Int i = 0; i < numMols; i++)
    {
        mol = &myMols[i];
        Vector3 Ri = mol->massCenter;
        for (Int j = 0; j < mol->numAtoms; j++)
        {
            Vector3 fj = mol->myAtoms[j]->force;
            Vector3 rj = mol->myAtoms[j]->position;
	    molCorrectLustig += fj*(rj - Ri);
            //molCorrectVirial[XX] += fj.x*(rj.x - Ri.x);
            //molCorrectVirial[YY] += fj.y*(rj.y - Ri.y);
            //molCorrectVirial[ZZ] += fj.z*(rj.z - Ri.z);
        }
    }

	myEnsemble->myLustig.virialcorrect = molCorrectLustig / (3.0 * volume);


// Jc : calculation in EwaldSummation in the shakeVirial[i],     

        myEnsemble->accumulate();    // Jc: calculate the pressure tensor in the 
	
	if (!strcmp(thermotype, "Berendsen")){
		if ((currentSteps%10) == 0){    //Jc: tempearture coupling every 10 steps
     		//if (currentSteps <= 0.25*totsteps){
        		 thermostat();    //Berendsen THERMOSTAT
      		//}
	    }
	}

	// -----------------------------------------------------------------------------------
	// Perform the integration of the Barostat viariable and
	// update the volume for the next time step
	// -----------------------------------------------------------------------------------
	
	if (NpT){
            if (has_constr_force) {
	        ww     = molNDF*BOLTZMAN*temperature*BarosConst*BarosConst;
	    }
	    else {
	        ww     = atomNDF*BOLTZMAN*temperature*BarosConst*BarosConst;
	    }
	
	    invww = 1.0/ww;
    	    nnhpl = nnhmi + 3*timeStep*invww*volume*(molPt - pressure);
            nnh   = (nnhpl + nnhmi)*0.5;
	    nnhmi = nnhpl;

	   Vpl    = Vmin*exp(3*timeStep*nnhpl);
	   boxLen = pow(Vpl, 1.0/3.0);
           myEnsemble->boxLx = boxLen;
           myEnsemble->boxLy = boxLen;
           myEnsemble->boxLz = boxLen;
	   myEnsemble->halfLx = 0.5*boxLen;
	   myEnsemble->halfLy = 0.5*boxLen;
	   myEnsemble->halfLz = 0.5*boxLen;

	   volume  = myEnsemble->boxLx*myEnsemble->boxLy*myEnsemble->boxLz;
	   density = numAtoms/volume;
           myEnsemble->volume = volume;
	   mySimConfig->set_box(boxLen, boxLen, boxLen);
	   mySimConfig->set_density(density);
	   Vmin    = Vpl;

           if(((currentSteps%100)==0)&&(prank==0)){
  	   	*ofNH << fixed << setprecision(6) << setw(14) << nnh << fixed << setprecision(6) << setw(14) << BarosConst;
	   	*ofNH << fixed << setprecision(6) << setw(14) << ww << endl;
	    }
	} 
	// -------------------------------------------------------------------------------
        // Apply PBC to mass center & then shift position to new center
        // and compute the center of mass of the systems usefull for 
        // NH barostat with Melchionna/Ciccotti correction
        // MOLECULAR PHYSICS, 1993, VOL. 78, No. 3, 533-544
        // -------------------------------------------------------------------------------

        for (Int nm = 0; nm < numMols; nm++)
        {  
            mol = &myMols[nm];

            position = mol->massCenter; 
            myEnsemble->apply_pbc(mol->massCenter);        
            diff = mol->massCenter - position;
            for (Int j = 0; j < mol->numAtoms; j++) {
                mol->myAtoms[j]->position += diff;
	    }
            Rm   += mol->mass * mol->massCenter;
            Mtot += mol->mass;

        }

	Rm  = Rm / Mtot ;

// Jc: calculate the Heat Capacity designed by Jianhui LI on 31 Mar 2006
        if((currentSteps>0 )&&(prank ==0)){

          E2Sum += (myEnsemble->myEnergy.totEnergy)*(myEnsemble->myEnergy.totEnergy);
          ESum +=  myEnsemble->myEnergy.totEnergy;
          aTSum += atomTemp;
          mTSum += molTemp;
          aT2Sum += atomTemp*atomTemp;
          mT2Sum += molTemp*molTemp;
// Jc: Heat capacity unit is [KJ/(Kg*K)], details see Jianhui Notebook on 1st April 2006
        if(((currentSteps%100)==0)&&(prank==0))
        {
            *ofPotential<<fixed<<setw(6)<<setprecision(2)<<currentSteps;
            *ofPotential<< setprecision(8) << setw(17)<< atomTemp;
            *ofPotential<< setprecision(8) << setw(17)<< molTemp;
            *ofPotential<< setprecision(8) << setw(17)<< molPt;
            *ofPotential<< setprecision(6) << setw(17)<< volume;
            *ofPotential<< setprecision(6) << setw(17)<< myEnsemble->boxLx;
            *ofPotential<< setprecision(3) << setw(9)<<myEnsemble->myEnergy.correctEnergy/numMols;
            *ofPotential<<setw(7)<<myEnsemble->myEnergy.inducEnergy/numMols;
            *ofPotential<<setw(6)<<myEnsemble->myEnergy.TBEnergy/numMols;                              
            *ofPotential<<setw(8)<<myEnsemble->myEnergy.totEnergy/numMols << endl;
        }
        }
//Jc: save the the vtMinus for continue simulation
        write_vtMinus();

    }
}

void LFIntegrator :: thermostat()
{
    Molecule *mol;
    int ns;
    Double factor;        // scale factor
    double vomx, vomy, vomz, masstott;
    double dttau = couplconst*timeStep;

    vomx = 0.0;
    vomy = 0.0;
    vomz = 0.0;
    masstott = 0.0;

    for (Int nm = 0; nm < numMols; nm++)     // 
    {                                        // 
        mol = &myMols[nm];
        ns = mol->numAtoms;
        for (Int i = 0; i < ns; i++)
        {
            Atom *atom = mol->myAtoms[i];
            vomx += atom->velocity.x*atom->mass;
            vomy += atom->velocity.y*atom->mass;
            vomz += atom->velocity.z*atom->mass;
            masstott += atom->mass;
        }

     }

vomx = vomx/masstott;
vomy = vomy/masstott;
vomz = vomz/masstott;

    if(useAtomThermo)
        factor = sqrt(1.0 + (timeStep/dttau)*(temperature/atomTemp - 1.0));
    else
        factor = sqrt(1.0 + (timeStep/dttau)*(temperature/molTemp - 1.0));

    for (Int i = 0; i < numAtoms; i++)
    {
        myAtoms[i].velocity *= factor;
        vtPlus[i] *= factor;
        vtMinus[i] *= factor;
    }


}

void LFIntegrator :: start()
{
    double forceArray[a][3],tempForceArray[a][3];

    clock_t startlong,endlong;  

    #ifdef DEBUG
        DEBUGMSG("start()");
    #endif
    
    for (Int i = 0; i < numAtoms; i++)
        myAtoms[i].force = 0.0;

    double tensor[9],tensorTemp[9];

    startlong = clock();

    if (computeBuff) {
        buffForce->compute();
    }
    else {
        if(prank ==0)
        {
    	    ljForce->compute();
    	    // mcyForce->compute();
    	    //mieForce->compute();
    	    //  saapForce->compute();
    	    //buckinghamForce->compute();
    	    //bbvForce->compute();
        }
        //nccForce->compute();            // JC changed by Jianhui LI since there is no LJ force in NCC
    }
    endlong   = clock();

    if (computeEwald )         
    {
	if (computeMultipole || computeWolf){
	    ERRORMSG("More than one Electrostatic method is ON, WHY BROOO???");
	}
        DEBUGMSG("compute Ewald");
        startlong = clock();
        ewaldForce->compute();      // Jc: force changed in real_term and intraMolecluar Energy, surfcorrection
        endlong = clock();

        startlong = clock();      
        ewaldForce->compute_long();
        endlong   = clock();

// Jc: reduce virial calculated in NccForce and EwaldForce respectively
        for (int tensorIndex = XX;tensorIndex <= ZZ; tensorIndex++)
        {
          tensor[tensorIndex]     = myEnsemble->longVirial[tensorIndex];
          tensorTemp[tensorIndex] = myEnsemble->longVirial[tensorIndex];
        }
     
//        MPI::COMM_WORLD.Allreduce(tensor,tensorTemp,9,MPI::DOUBLE,MPI::SUM);
        MPI_Allreduce(tensor,tensorTemp,9,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

        for (int tensorIndex =XX;tensorIndex <= ZZ; tensorIndex++)
        {
          tensor[tensorIndex]                  = tensorTemp[tensorIndex];
          myEnsemble->longVirial[tensorIndex]  = tensor[tensorIndex];
        } 

    }
    else if (computeMultipole)         
    {
	if (computeEwald || computeWolf){
	    ERRORMSG("More than one Electrostatic method is ON, WHY BROOO???");
	}
        DEBUGMSG("compute Multipole Electrostatic with Ewald");
        startlong = clock();
        multiElec->compute();      // Jc: force changed in real_term and intraMolecluar Energy, surfcorrection
        endlong = clock();

        startlong = clock();      
        multiElec->compute_long();
        endlong   = clock();

        if (computeInduction)                             
        {
            inductionForce->compute(); 
                
            inductionForce->compute_long();
            DEBUGMSG("Include to the Multipole Electrostatic the Induction computation with Ewald");
	}

// Jc: reduce virial calculated in Multipole 
        for (int tensorIndex = XX;tensorIndex <= ZZ; tensorIndex++)
        {
          tensor[tensorIndex]     = myEnsemble->longVirial[tensorIndex];
          tensorTemp[tensorIndex] = myEnsemble->longVirial[tensorIndex];
        }
     
//        MPI::COMM_WORLD.Allreduce(tensor,tensorTemp,9,MPI::DOUBLE,MPI::SUM);
        MPI_Allreduce(tensor,tensorTemp,9,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

        for (int tensorIndex =XX;tensorIndex <= ZZ; tensorIndex++)
        {
          tensor[tensorIndex]                  = tensorTemp[tensorIndex];
          myEnsemble->longVirial[tensorIndex]  = tensor[tensorIndex];
        } 

    }
    else if (computeWolf)                             
    {
        DEBUGMSG("compute Wolf Electrostatic");
	wolfForce->compute();
    }
    else                              
    {
        DEBUGMSG("NO ELECTROSTATIC METHOD HAS BEEN IMPLEMENTED!!!!!!!!");
        DEBUGMSG("NO ELECTROSTATIC METHOD HAS BEEN IMPLEMENTED!!!  IS THIS OK??");
    }
//Jc: commented for 2-body calcualtion
 
    for (int numIndex = 0;numIndex < numAtoms; numIndex++)
    {
      forceArray[numIndex][0]     = myAtoms[numIndex].force.x;
      forceArray[numIndex][1]     = myAtoms[numIndex].force.y;
      forceArray[numIndex][2]     = myAtoms[numIndex].force.z;
    }
   
//    MPI::COMM_WORLD.Allreduce(forceArray,tempForceArray,a3,MPI::DOUBLE,MPI::SUM);
    MPI_Allreduce(forceArray,tempForceArray,a3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    for (int numIndex = 0;numIndex < numAtoms; numIndex++)
    {

      myAtoms[numIndex].force.x  = tempForceArray[numIndex][0];
      myAtoms[numIndex].force.y  = tempForceArray[numIndex][1];
      myAtoms[numIndex].force.z  = tempForceArray[numIndex][2];
    } 

// Jc: reduce virial calculated in NccForce and EwaldForce

    for (int tensorIndex = XX;tensorIndex <= ZZ; tensorIndex++)
    {
      tensor[tensorIndex]     = myEnsemble->virial[tensorIndex];
      tensorTemp[tensorIndex] = myEnsemble->virial[tensorIndex];
    }
     
//    MPI::COMM_WORLD.Allreduce(tensor,tensorTemp,9,MPI::DOUBLE,MPI::SUM);
    MPI_Allreduce(tensor,tensorTemp,9,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    for (int tensorIndex =XX;tensorIndex <= ZZ; tensorIndex++)
    {
      tensor[tensorIndex]              = tensorTemp[tensorIndex];
      myEnsemble->virial[tensorIndex]  = tensor[tensorIndex];
    }
 
// Jc: end of partial result reduce

//    for (Int i = 0; i < numAtoms; i++)
//   { 
//      myAtoms[i].force += myAtoms[i].longForce;
//      if(myAtoms[i].longForce.x != 0) cout<<" long force X of atom "<<i<<" is"<<myAtoms[i].longForce.x<<endl;
//      if(myAtoms[i].longForce.y != 0) cout<<" long force Y of atom "<<i<<" is"<<myAtoms[i].longForce.y<<endl;
//      if(myAtoms[i].longForce.z != 0) cout<<" long force Z of atom "<<i<<" is"<<myAtoms[i].longForce.z<<endl;
//    }
    if (has_angle_force)
    {
        DEBUGMSG("compute angle" );
        cout<<".......calcualting angle force ............."<<endl;
        angleForce->compute();// JC is there any angle force in the NCC calculation
    }

    if (has_bond_force)
    {
        DEBUGMSG("compute bond");
        cout<<".......calcualting bond force ............."<<endl;
        bondForce->compute();
    }
    if (use_ub_harmonic)
    {
        DEBUGMSG("compute Urey-Bradley Force");
        cout<<".......calcualting Urey-Bradley Force ............."<<endl;
        ubForce->compute();
    }
//    MPI::COMM_WORLD.Barrier();   
    MPI_Barrier(MPI_COMM_WORLD);
    for (Int i = 0; i < numAtoms; i++)
        vtMinus[i] = myAtoms[i].velocity - 0.5*timeStep*mass_r[i]*myAtoms[i].force ;
}

void LFIntegrator :: shakeCFR(Double tol, Int maxCount, Double *virial, int currentSteps)
{
    Atom* atom;
    Molecule *mol;
    Int ns, nBonds;
    
    double vomx, vomy, vomz, masstott;
    Double tol_2 = 2*tol;

    for (Int i = XX; i <= ZZ; i++)
    {    
      virial[i] = 0.0;
    }
    Int count = 0;
    Int next;
    Vector3 rijNew, rijOld;
    Double rijNew2, rijNewOld, diff, factor; // rijOld2,commented by Jianhui
    bool done;


//    int numDecompose = (int) numMols/psize;
//    if ((numMols % psize) == 0) numDecompose--;

//    for(int i = 0; i <= numDecompose;i++)   
//    { 
//     int nm = i*psize + prank;
//     if(nm <= (numMols -1)) 
//     {

    for (Int nm = 0; nm < numMols; nm++)     // 
    {                                        // 
        mol = &myMols[nm];
        ns = mol->numAtoms;
        for (Int i = 0; i < ns; i++) 
        {
            Atom *atom = mol->myAtoms[i];
            vomx += atom->velocity.x*atom->mass;
            vomy += atom->velocity.y*atom->mass;
            vomz += atom->velocity.z*atom->mass;
            masstott += atom->mass;
	}
     }	

vomx = vomx/masstott;
vomy = vomy/masstott;
vomz = vomz/masstott;

//cout << "vomx  = " << vomx << " vomy = " << vomy << "  vomz  =  " << vomz << endl;

	if (fabs(vomx) > 1e-11 || fabs(vomy) > 1e-11 || fabs(vomz) > 1e-11)
		for (Int nm = 0; nm < numMols; nm++)     // 
    		{                                        // 
        		mol = &myMols[nm];
        		ns = mol->numAtoms;
        		for (Int i = 0; i < ns; i++)
        		{
            			Atom *atom = mol->myAtoms[i];
				Int id = atom->atomID;
            			vtMinus[id].x -=  vomx;
            			vtMinus[id].y -=  vomy;
            			vtMinus[id].z -=  vomz;
        		}
    		}


    for (Int nm = 0; nm < numMols; nm++)     // apply constraint to each molecule By ZhongWu
    {                                        // created by ZhongWu
        mol = &myMols[nm]; 
        ns = mol->numAtoms;

        // we only treat water as cyclic molecule here and water must have a name of "CO2"
        if (!strcmp(mol->molName, "CO2"))     // Only CO2 
           {
            nBonds = ns - 2; 

        // update atomic positions
            double MM = mol->mass; 
//cout<<"  MolMass = "<< MM  << "  MolType = "<<mol->molType <<"  MolID = "<<mol->molID<<"  ns = "<<ns<<"  nBonds = "<<nBonds<<endl;
        //for (Int i = 0; i < (ns-1); i++) // calculate the new assumed positions of atoms on nmth molecule
       // {
            Atom *atom0 = mol->myAtoms[0];
            Atom *atomc = mol->myAtoms[1];
            Atom *atom2 = mol->myAtoms[2];
            Int id0 = atom0->atomID;
            Int id1 = atomc->atomID;
            Int id2 = atom2->atomID;
	    double massC = atomc->mass;
	    double massO = atom0->mass;
       //cout<<"  id0 =  "<< id0 << "   id1 =  "<< id1 << "  mass0 = "<< massO << "  massC = "<< massC <<endl;
            rOld[nm][0] = atom0->position;
            if (!strcmp(thermotype, "NoseHoover")){
               // if (currentSteps <= 0.25*totsteps){
	            vtPlus[id0] = vtMinus[id0] + timeStep*mass_r[id0]*(atom0->force - (massC/MM)*0.5*atom0->force + (massO/MM)*atomc->force - (massC/MM)*0.5*atom2->force - fnh*mol->momenta/ns);
	/*	}
		else {
	            vtPlus[id0] = vtMinus[id0] + timeStep*mass_r[id0]*(atom0->force - (massC/MM)*0.5*atom0->force + (massO/MM)*atomc->force - (massC/MM)*0.5*atom2->force);
		}*/
	    }
	    else {
	            vtPlus[id0] = vtMinus[id0] + timeStep*mass_r[id0]*(atom0->force - (massC/MM)*0.5*atom0->force + (massO/MM)*atomc->force - (massC/MM)*0.5*atom2->force);
	    }
            rNew[nm][0] = rOld[nm][0] + timeStep*vtPlus[id0];
            //moving[nm][0] = false;
            //moved[nm][0]  = true;

            rOld[nm][2] = atom2->position;
            if (!strcmp(thermotype, "NoseHoover")){
              //  if (currentSteps <= 0.25*totsteps){
            	     vtPlus[id2] = vtMinus[id2] + timeStep*mass_r[id2]*(atom2->force - (massC/MM)*0.5*atom2->force + (massO/MM)*atomc->force - (massC/MM)*0.5*atom0->force - fnh*mol->momenta/ns );
	/*	}
	    	else {
            	     vtPlus[id2] = vtMinus[id2] + timeStep*mass_r[id2]*(atom2->force - (massC/MM)*0.5*atom2->force + (massO/MM)*atomc->force - (massC/MM)*0.5*atom0->force);
	        }*/
	    }
	    else {
            	     vtPlus[id2] = vtMinus[id2] + timeStep*mass_r[id2]*(atom2->force - (massC/MM)*0.5*atom2->force + (massO/MM)*atomc->force - (massC/MM)*0.5*atom0->force);
	    }
            rNew[nm][2] = rOld[nm][2] + timeStep*vtPlus[id2];
            //moving[nm][1] = false;
            //moved[nm][1]  = true;

            rOld[nm][1] = atomc->position;
            //vtPlus[id2] = vtMinus[id2] + timeStep*mass_r[id2]*(atomc->force + (massC/MM)*(atom->force + atom1->force) - 2*(massO/MM)*atomc->force);
            //vtPlus[id2] = (vtPlus[id0] + vtPlus[id1])*0.5 ;
            rNew[nm][1] = (rNew[nm][2] + rNew[nm][0])*0.5;
            //moving[nm][2] = false;
            //moved[nm][2]  = true;
       // }

             Vector3 co1dif = rOld[nm][0] - rOld[nm][1];      // Jc the actual bond length 
             double diffco1 = co1dif.length();
             Vector3 co2dif =  rOld[nm][1] -  rOld[nm][2];      // Jc the actual bond length 
             double diffco2 = co2dif.length();
//cout<<"  diffco1   =  "<< diffco1 << "  diffco2 =  "<< diffco2 << endl;

//cout << bondLength[nm][0] << endl;

        // iterative application of constraints
        count = 0;
        done = false;
        while ((done == false) && (count < maxCount))           // JC after MaxCount times and the program breaks down!
        {   
            done = true;
                    rijNew = rNew[nm][0] - rNew[nm][2];      // Jc the actual bond length
	 	    myEnsemble->apply_pbc(rijNew); 
                    rijNew2 = rijNew.length2(); 
                    diff = 4.0*bondLength[nm][0] - rijNew2;         // JC the difference between 

                    if (fabs(diff) > (bondLength[nm][0]*tol_2)) // JC start Shake algorithm
                    {                        
                        rijOld = rOld[nm][0] - rOld[nm][2];  // 
	 	        myEnsemble->apply_pbc(rijOld); 
                        rijNewOld = rijOld*rijNew;              // inner dot product the criterion
                        if (rijNewOld < 4.0*bondLength[nm][0]*tol)
                        {
                          cout<<" Mol not converged " <<nm<<" at pos "<<mol->myAtoms[0]->position<<endl;
                          cout<<" Mol not converged " <<nm<<" at pos "<<mol->myAtoms[1]->position<<endl;
                          cout<<" Mol not converged " <<nm<<" at pos "<<mol->myAtoms[2]->position<<endl;
 
                          ERRORMSG("Shake constraint failure");
                        }
 
                        factor = diff/(2*rijNewOld*(mass_r[id0] + mass_r[id2]));                         
                        //factor = diff/(2*rijNewOld*(2/MM + 2/MM));                         
                        // negative virial or stress contribution

                        virial[XX] += factor*rijNew.x*rijOld.x;  // Jc: these are the Shake virial from the Shake algorithm
                        virial[XY] += factor*rijNew.x*rijOld.y;
                        virial[XZ] += factor*rijNew.x*rijOld.z;
                        virial[YX] += factor*rijNew.y*rijOld.x;
                        virial[YY] += factor*rijNew.y*rijOld.y;
                        virial[YZ] += factor*rijNew.y*rijOld.z;
                        virial[ZX] += factor*rijNew.z*rijOld.x;
                        virial[ZY] += factor*rijNew.z*rijOld.y;
                        virial[ZZ] += factor*rijNew.z*rijOld.z; 

                        rNew[nm][0] += factor*rijOld*mass_r[id0];
                        rNew[nm][2] -= factor*rijOld*mass_r[id2];
		        rNew[nm][1] = (rNew[nm][0] + rNew[nm][2])*0.5;
                        //moving[nm][0] = true;
                        //moving[nm][1] = true;
                        //moving[nm][2] = true;
                        done = false;
                    }                        
            //for (Int i = 0; i < ns; i++)
            //{
            //    moved[nm][i] = moving[nm][i];
            //    moving[nm][i] = false;
            //} 
            count++;
        }   // complete iterative loop for one molecule
       

        if(!done)
        {
            for (Int i = 0; i < ns; i++) // calculate the new assumed positions of atoms on nmth molecule
            {
              Atom *atom = mol->myAtoms[i];
              Int id = atom->atomID;
              rOld[nm][i] = atom->position;
              cout<<"the number of atoms "<<id<<endl;
              cout<<"the position of the atom are x=  "<<rOld[nm][i].x<<endl;
              cout<<"the position of the atom are y=  "<<rOld[nm][i].y<<endl;
              cout<<"the position of the atom are z=  "<<rOld[nm][i].z<<endl;
             }
            cerr << "SHAKE not converge after " << maxCount << " iterative loops on Molecule: " << nm << endl;
            ERRORMSG("try large iterative loops");
           
        }
        // update position and velocity for this molecule
        for(Int i = 0; i < ns; i++)
        {
            atom = mol->myAtoms[i];
            Int id = atom->atomID;
            Vector3 dr = rNew[nm][i] - atom->position;
            atom->position += dr;
            atom->realPos += dr;
            atom->displacement += dr;
            vtPlus[id] = dr/timeStep;
            atom->velocity = (vtPlus[id] + vtMinus[id])*0.5;
            vtMinus[id] = vtPlus[id];
        }
      }
    }
    
    for (Int i = XX; i <= ZZ; i++)
        virial[i] /= (timeStep*timeStep);
//  }      //Jc: parallel decomposition end
}

void LFIntegrator :: shakeBBV(Double tol, Int maxCount, Double *virial, int currentSteps)
{
    Atom* atom;
    Molecule *mol;
    Int ns, nBonds;

    double vomx, vomy, vomz, masstott;
    Double tol_2 = 2*tol;
    const double ll = 2*co2bond;
    const double CC = llra;
    const double BB = ll - 2*CC;                // same as Lvir
    const double adist = 0.5*(Lvir - co2bond);
    const double adiff = Lvir - adist;
 //cout << " CC =" << CC << "     BB = " << BB << "    Lvir =" << Lvir << "   adist = " << adist << "  1-adist " << adiff << endl;

    for (Int i = XX; i <= ZZ; i++)
    {
      virial[i] = 0.0;
    }
    Int count = 0;
    Int next;
    Vector3 rijNew, rijOld;
    Double rijNew2, rijNewOld, diff, factor; // rijOld2,commented by Jianhui
    bool done;
    
//    int numDecompose = (int) numMols/psize;
//    if ((numMols % psize) == 0) numDecompose--;
    
//    for(int i = 0; i <= numDecompose;i++)   
//    { 
//     int nm = i*psize + prank;
//     if(nm <= (numMols -1)) 
//     {

    for (Int nm = 0; nm < numMols; nm++)     // 
    {                                        // 
        mol = &myMols[nm];
        ns = mol->numAtoms;
        //for (Int i = 0; i < ns; i++)
        {
            Atom *atom0 = mol->myAtoms[0];
            Atom *atom1 = mol->myAtoms[4];
            vomx += 0.5*atom1->velocity.x*mol->mass + 0.5*atom0->velocity.x*mol->mass;
            vomy += 0.5*atom1->velocity.y*mol->mass + 0.5*atom0->velocity.y*mol->mass;
            vomz += 0.5*atom1->velocity.z*mol->mass + 0.5*atom0->velocity.z*mol->mass;
            masstott += mol->mass;

        }
     }

vomx = vomx/masstott;
vomy = vomy/masstott;
vomz = vomz/masstott;

//cout << "vomx  = " << vomx << " vomy = " << vomy << "  vomz  =  " << vomz << endl;

	if (fabs(vomx) > 1e-11 || fabs(vomy) > 1e-11 || fabs(vomz) > 1e-11)
                for (Int nm = 0; nm < numMols; nm++)     // 
                {                                        // 
                        mol = &myMols[nm];
                        ns = mol->numAtoms;
                        for (Int i = 0; i < ns; i++)
                        {
                                Atom *atom = mol->myAtoms[i];
                                Int id = atom->atomID;
                                vtMinus[id].x -=  vomx;
                                vtMinus[id].y -=  vomy;
                                vtMinus[id].z -=  vomz;
                        }
                }

     VirtualPositions();

    for (Int nm = 0; nm < numMols; nm++)     // apply constraint to each molecule By ZhongWu
    {                                        // created by ZhongWu
        mol = &myMols[nm];
        ns = mol->numAtoms;

        // we only treat water as cyclic molecule here and water must have a name of "CO2"
        if (!strcmp(mol->molName, "BBV"))     // Only CO2 5 particle molecule
           {
            nBonds = ns - 2;

        // update atomic positions
            double MM = mol->mass;
//cout<<"  MolMass = "<< MM  << "  MolType = "<<mol->molType <<"  MolID = "<<mol->molID<<"  ns = "<<ns<<"  nBonds = "<<nBonds<<endl;
	//for (Int i = 0; i < (ns-1); i++) // calculate the new assumed positions of atoms on nmth molecule
	// {
	    Atom *atom0  = mol->myAtoms[0];
            Atom *atomB1 = mol->myAtoms[1];
            Atom *atomc  = mol->myAtoms[2];
            Atom *atomB2 = mol->myAtoms[3];
            Atom *atom1  = mol->myAtoms[4];

            Int id0 = atom0->atomID;
            Int id1 = atom1->atomID;
            Int id2 = atomc->atomID;
            Int id3 = atomB1->atomID;
            Int id4 = atomB2->atomID;
            double massC = atomc->mass;
            double massO = atom0->mass;
//       cout<<"  id0 =  "<< id0 << "   id1 =  "<< id1 << "  mass0 = "<< massO << "  massC = "<< massC <<endl;

	    rOld[nm][0] = virt1[nm][0];
            virt2[nm][0] = atom0->force;// - (massC/MM)*0.5*atom->force + (massO/MM)*atomc->force - (massC/MM)*0.5*atom1->force;

            //vtPlus[id0] = vtMinus[id0] + timeStep*mass_r[id0]*(atom->force - (massC/MM)*0.5*atom->force + (massO/MM)*atomc->force - (massC/MM)*0.5*atom1->force);
            //rNew[nm][0] = rOld[nm][0] + timeStep*vtPlus[id0];

	    rOld[nm][1] = virt1[nm][1];
            virt2[nm][1] = atom1->force;// - (massC/MM)*0.5*atom1->force + (massO/MM)*atomc->force - (massC/MM)*0.5*atom->force;

            //rOld[nm][2] = virt1[nm][0];
	    virt2[nm][2] = atomc->force;        // do I need this??

            virtfor[nm][0] = -CC/BB*virt2[nm][1] + (1+CC/BB)*virt2[nm][0] + 0.5*atomc->force + adiff/Lvir*atomB1->force + adist/Lvir*atomB2->force;
            virtfor[nm][1] =  (1+CC/BB)*virt2[nm][1] - CC/BB*virt2[nm][0] + 0.5*atomc->force + adist/Lvir*atomB1->force + adiff/Lvir*atomB2->force;

            vtPlus[id0] = vtMinus[id0] + timeStep*2.0/MM*virtfor[nm][0];    // Half the total mass
            rNew[nm][0] = rOld[nm][0] + timeStep*vtPlus[id0];

            vtPlus[id1] = vtMinus[id1] + timeStep*2.0/MM*virtfor[nm][1];
            rNew[nm][1] = rOld[nm][1] + timeStep*vtPlus[id1];

            rNew[nm][2] = (rNew[nm][1] + rNew[nm][0])*0.5;
       // }
       //cout << "adiff/Lvir = " << adiff/Lvir << "     adist/Lvir = " << adist/Lvir << endl;
       //             Vector3 co1dif = rOld[nm][0] - rOld[nm][2];      // Jc the actual bond length 
       //             double diffco1 = co1dif.length();
       //             Vector3 co2dif =  rOld[nm][1] -  rOld[nm][2];      // Jc the actual bond length 
       //             double diffco2 = co2dif.length();
       //cout<<"  diffco1   =  "<< diffco1 << "  diffco2 =  "<< diffco2 << endl;
      
       
       // iterative application of constraints
        count = 0;
        done = false;
        while ((done == false) && (count < maxCount))           // JC after MaxCount times and the program breaks down!
        {
            done = true;
                    rijNew = rNew[nm][0] - rNew[nm][1];      // Jc the actual bond length
                    myEnsemble->apply_pbc(rijNew);
                    rijNew2 = rijNew.length2();
                    diff = BB*BB - rijNew2;         // JC the difference between 

                    if (fabs(diff) > (BB*BB*tol_2)) // JC start Shake algorithm
                    {
                        rijOld = rOld[nm][0] - rOld[nm][1];  // 
                        myEnsemble->apply_pbc(rijOld);
                        rijNewOld = rijOld*rijNew;              // inner dot product the criterion
                        if (rijNewOld < BB*BB*tol)
                        {
                          cout<<" Mol not converged " <<nm<<" at pos "<<mol->myAtoms[0]->position<<endl;
                          cout<<" Mol not converged " <<nm<<" at pos "<<mol->myAtoms[1]->position<<endl;
                          cout<<" Mol not converged " <<nm<<" at pos "<<mol->myAtoms[2]->position<<endl;

                          ERRORMSG("Shake constraint failure");
                        }

                        //factor = diff/(2*rijNewOld*(mass_r[id0] + mass_r[id1]));                  
                        factor = diff/(2*rijNewOld*(2/MM + 2/MM));
                        // negative virial or stress contribution

			virial[XX] += factor*rijNew.x*rijOld.x;  // Jc: these are the Shake virial from the Shake algorithm
                        virial[XY] += factor*rijNew.x*rijOld.y;
                        virial[XZ] += factor*rijNew.x*rijOld.z;
                        virial[YX] += factor*rijNew.y*rijOld.x;
                        virial[YY] += factor*rijNew.y*rijOld.y;
                        virial[YZ] += factor*rijNew.y*rijOld.z;
                        virial[ZX] += factor*rijNew.z*rijOld.x;
                        virial[ZY] += factor*rijNew.z*rijOld.y;
                        virial[ZZ] += factor*rijNew.z*rijOld.z;

                        rNew[nm][0] += factor*rijOld*2/MM;//mass_r[id0];
                        rNew[nm][1] -= factor*rijOld*2/MM;//mass_r[id1];
                        rNew[nm][2] = (rNew[nm][0] + rNew[nm][1])*0.5;
                        //moving[nm][0] = true;
                        //moving[nm][1] = true;
                        //moving[nm][2] = true;
                        done = false;
                    }
            //for (Int i = 0; i < ns; i++)
            //{
            //    moved[nm][i] = moving[nm][i];
            //    moving[nm][i] = false;
            //} 
            count++;
        }   // complete iterative loop for one molecule


        if(!done)
        {
            for (Int i = 0; i < ns; i++) // calculate the new assumed positions of atoms on nmth molecule
            {
              Atom *atom = mol->myAtoms[i];
              Int id = atom->atomID;
              rOld[nm][i] = atom->position;
              cout<<"the number of atoms "<<id<<endl;
              cout<<"the position of the atom are x=  "<<rOld[nm][i].x<<endl;
              cout<<"the position of the atom are y=  "<<rOld[nm][i].y<<endl;
              cout<<"the position of the atom are z=  "<<rOld[nm][i].z<<endl;
             }
            cerr << "SHAKE not converge after " << maxCount << " iterative loops on Molecule: " << nm << endl;
            ERRORMSG("try large iterative loops");

        }
        // update position and velocity for this molecule
        //for(Int i = 0; i < ns; i++)
        //{
        //

	    Vector3 newpos0  = -CC/BB*rNew[nm][1] + (1+CC/BB)*rNew[nm][0];
            Vector3 newpos1  = -CC/BB*rNew[nm][0] + (1+CC/BB)*rNew[nm][1];
            //Vector3 newposB1 = adiff/Lvir*rNew[nm][0] + adist/Lvir*rNew[nm][1];
            //Vector3 newposB2 = adist/Lvir*rNew[nm][0] + adiff/Lvir*rNew[nm][1];
            Vector3 newposB1 = (newpos0 + rNew[nm][2])*0.5;
            Vector3 newposB2 = (newpos1 + rNew[nm][2])*0.5;

            Vector3 drr0  = newpos0 - atom0->position;
            Vector3 drr1  = newpos1 - atom1->position;
            Vector3 drrB1 = newposB1 - atomB1->position;
            Vector3 drrB2 = newposB2 - atomB2->position;
            Vector3 drrc  = rNew[nm][2] - atomc->position;

            atom0->position += drr0;
            atom0->realPos += drr0;
            atom0->displacement += drr0;

            atom1->position += drr1;
            atom1->realPos += drr1;
            atom1->displacement += drr1;

            atomB1->position += drrB1;
            atomB1->realPos += drrB1;
            atomB1->displacement += drrB1;

            atomB2->position += drrB2;
            atomB2->realPos += drrB2;
            atomB2->displacement += drrB2;

	    atomc->position += drrc;
            atomc->realPos += drrc;
            atomc->displacement += drrc;

            vtPlus[id0] = (rNew[nm][0] - rOld[nm][0])/timeStep;
            vtPlus[id1] = (rNew[nm][1] - rOld[nm][1])/timeStep;
            vtPlus[id2] = 0.0;
            vtPlus[id3] = 0.0;
            vtPlus[id4] = 0.0;

            atom0->velocity = (vtPlus[id0] + vtMinus[id0])*0.5;
            atom1->velocity = (vtPlus[id1] + vtMinus[id1])*0.5;
            atomc->velocity = 0.0;
            atomB1->velocity = 0.0;
            atomB2->velocity = 0.0;

            vtMinus[id0] = vtPlus[id0];
            vtMinus[id1] = vtPlus[id1];
            vtMinus[id2] = vtPlus[id2];
            vtMinus[id3] = vtPlus[id3];
            vtMinus[id4] = vtPlus[id4];

        //}
      }
    }

    for (Int i = XX; i <= ZZ; i++)
        virial[i] /= (timeStep*timeStep);
//  }      //Jc: parallel decomposition end
}


void LFIntegrator :: shake(Double tol, Int maxCount, Double *virial, int currentSteps)
{
    Atom* atom;
    Molecule *mol;
    Int ns, nBonds;

    Double tol_2 = 2*tol;

    for (Int i = XX; i <= ZZ; i++)
    {    
      virial[i] = 0.0;
    }
    Int count = 0;
    Int next;
    Vector3 rijNew, rijOld;
    Double rijNew2, rijNewOld, diff, factor; // rijOld2,commented by Jianhui
    bool done;


    for (Int nm = 0; nm < numMols; nm++)     // apply constraint to each molecule By ZhongWu
    {                                        // created by ZhongWu
        mol = &myMols[nm]; 
        ns = mol->numAtoms;

        // we only treat water as cyclic molecule here and water must have a name of "WAT" or "MCY"
        if ( (!strcmp(mol->molName, "WAT"))  || (!strcmp(mol->molName, "MCY"))  )     //  water 
           {
            nBonds = ns; 
           }
        else
           {    
            nBonds = ns - 1; 
           }

        // update atomic positions
      if (!NpT){
        for (Int i = 0; i < ns; i++) // calculate the new assumed positions of atoms on nmth molecule
        {
            Atom *atom = mol->myAtoms[i];
            Int id = atom->atomID;
            rOld[nm][i] = atom->position;
	 if (!strcmp(thermotype, "NoseHoover")){
     	    //if (currentSteps <= 0.25*totsteps){
            	vtPlus[id] = vtMinus[id] + timeStep*mass_r[id]*( atom->force - fnh*mol->momenta/ns );
	    /*}
	    else{
	    	vtPlus[id] = vtMinus[id] + timeStep*mass_r[id]*atom->force;
	    }*/
	 }
	 else{
	    vtPlus[id] = vtMinus[id] + timeStep*mass_r[id]*atom->force;
	 }
            rNew[nm][i] = rOld[nm][i] + timeStep*vtPlus[id];
            moving[nm][i] = false;
            moved[nm][i]  = true;
        }

        // iterative application of constraints
        count = 0;
        done = false;
        while ((done == false) && (count < maxCount))           // JC after MaxCount times and the program breaks down!
        {   
            done = true;
            for (Int i = 0; i < nBonds; i++)
            {   
                next = i + 1;
                if (next == ns)
                    next = 0;
                // if (moved[nm][i] || moved[nm][next])
                {
                    rijNew = rNew[nm][i] - rNew[nm][next];      // Jc the actual bond length 
		    myEnsemble->apply_pbc(rijNew);
                    rijNew2 = rijNew.length2(); 
                    diff = bondLength[nm][i] - rijNew2;         // JC the difference between 

                    if (fabs(diff) > (bondLength[nm][i]*tol_2)) // JC start Shake algorithm
                    {                        
                        rijOld = rOld[nm][i] - rOld[nm][next];  // 
		        myEnsemble->apply_pbc(rijOld);
                        rijNewOld = rijOld*rijNew;              // inner dot product the criterion
                        if (rijNewOld < bondLength[nm][i]*tol)
                        {
                          cout<<" Mol not converge are " <<nm<<" "<<mol->myAtoms[0]->position<<endl;
                          cout<<" Mol not converge are " <<nm<<" "<<mol->myAtoms[1]->position<<endl;
                          cout<<" Mol not converge are " <<nm<<" "<<mol->myAtoms[2]->position<<endl;
 
                          ERRORMSG("Shake constraint failure");
                        } 
                        Int id1 = mol->myAtoms[i]->atomID;
                        Int id2 = mol->myAtoms[next]->atomID;
                        factor = diff/(2*rijNewOld*(mass_r[id1] + mass_r[id2]));                         
                        // negative virial or stress contribution

                        virial[XX] += factor*rijNew.x*rijOld.x;  // Jc: these are the Shake virial from the Shake algorithm
                        virial[XY] += factor*rijNew.x*rijOld.y;
                        virial[XZ] += factor*rijNew.x*rijOld.z;
                        virial[YX] += factor*rijNew.y*rijOld.x;
                        virial[YY] += factor*rijNew.y*rijOld.y;
                        virial[YZ] += factor*rijNew.y*rijOld.z;
                        virial[ZX] += factor*rijNew.z*rijOld.x;
                        virial[ZY] += factor*rijNew.z*rijOld.y;
                        virial[ZZ] += factor*rijNew.z*rijOld.z; 

                        rNew[nm][i] += factor*rijOld*mass_r[id1];
                        rNew[nm][next] -= factor*rijOld*mass_r[id2];
                        moving[nm][i] = true;
                        moving[nm][next] = true;
                        done = false;
                    }                        
                }
            }
            for (Int i = 0; i < nBonds; i++)
            {
                moved[nm][i] = moving[nm][i];
                moving[nm][i] = false;
            } 
            count++;
        }   // complete iterative loop for one molecule
       

        if(!done)
        {
            for (Int i = 0; i < ns; i++) // calculate the new assumed positions of atoms on nmth molecule
            {
              Atom *atom = mol->myAtoms[i];
              Int id = atom->atomID;
              rOld[nm][i] = atom->position;
              cout<<"the number of atoms "<<id<<endl;
              cout<<"the position of the atom are x=  "<<rOld[nm][i].x<<endl;
              cout<<"the position of the atom are y=  "<<rOld[nm][i].y<<endl;
              cout<<"the position of the atom are z=  "<<rOld[nm][i].z<<endl;
             }
            cerr << "SHAKE not converge after " << maxCount << " iterative loops on Molecule: " << nm << endl;
            ERRORMSG("try large iterative loops");
           
        }
        // update position and velocity for this molecule
        for(Int i = 0; i < ns; i++)
        {
            atom = mol->myAtoms[i];
            Int id = atom->atomID;
            Vector3 dr = rNew[nm][i] - atom->position;
            atom->position += dr;
            atom->realPos += dr;
            atom->displacement += dr;
            vtPlus[id] = dr/timeStep;
            atom->velocity = (vtPlus[id] + vtMinus[id])*0.5;
            vtMinus[id] = vtPlus[id];
        }

    }  // end of NVE/NVT ensemble loop
    else {

        Vector3 velvect;
        Vector3 posvect;
	Vector3 rhalfp;
	Vector3 momentime;
        double  velcheck = 1.0;
        double  poscheck = 1.0;
        double  NHcount = 0;

	momentime = mol->momenta;
	rhalfp = mol->massCenter;

	while ( !( (poscheck < 1E-05) && (velcheck < 1E-05) ) &&  (NHcount < 11) ){
            for (Int i = 0; i < ns; i++) // calculate the new assumed positions of atoms on nmth molecule
            {
                Atom *atom = mol->myAtoms[i];
                Int id = atom->atomID;
                rOld[nm][i] = atom->position;

		if (useAtomThermo) {
                    vtPlus[id] = vtMinus[id] + timeStep*mass_r[id]*( atom->force - fnh*atom->velocity*atom->mass - atom->mass*nnh*momentime/mol->mass); 
		}
		else {
                    vtPlus[id] = vtMinus[id] + timeStep*mass_r[id]*( atom->force - atom->mass*fnh*momentime/mol->mass - atom->mass*nnh*momentime/mol->mass);
		}

                rNew[nm][i] = rOld[nm][i] + timeStep*(vtPlus[id] + nnhpl*(rhalfp - Rm) );
                moving[nm][i] = false;
                moved[nm][i]  = true;
            }

		// iterative application of constraints
		count = 0;
		done = false;
		while ((done == false) && (count < maxCount))           // JC after MaxCount times and the program breaks down!
		{   
		    done = true;
		    for (Int i = 0; i < nBonds; i++)
		    {   
			next = i + 1;
			if (next == ns)
			    next = 0;
			// if (moved[nm][i] || moved[nm][next])
			{
			    rijNew = rNew[nm][i] - rNew[nm][next];      // Jc the actual bond length 
			    myEnsemble->apply_pbc(rijNew);
			    rijNew2 = rijNew.length2(); 
			    diff = bondLength[nm][i] - rijNew2;         // JC the difference between 

			    if (fabs(diff) > (bondLength[nm][i]*tol_2)) // JC start Shake algorithm
			    {                        
				rijOld = rOld[nm][i] - rOld[nm][next];  // 
				myEnsemble->apply_pbc(rijOld);
				rijNewOld = rijOld*rijNew;              // inner dot product the criterion
				if (rijNewOld < bondLength[nm][i]*tol)
				{
				  cout<<" Mol not converge are " <<nm<<" "<<mol->myAtoms[0]->position<<endl;
				  cout<<" Mol not converge are " <<nm<<" "<<mol->myAtoms[1]->position<<endl;
				  cout<<" Mol not converge are " <<nm<<" "<<mol->myAtoms[2]->position<<endl;
	 
				  ERRORMSG("Shake constraint failure");
				} 
				Int id1 = mol->myAtoms[i]->atomID;
				Int id2 = mol->myAtoms[next]->atomID;
				factor = diff/(2*rijNewOld*(mass_r[id1] + mass_r[id2]));                         
				// negative virial or stress contribution

				virial[XX] += factor*rijNew.x*rijOld.x;  // Jc: these are the Shake virial from the Shake algorithm
				virial[XY] += factor*rijNew.x*rijOld.y;
				virial[XZ] += factor*rijNew.x*rijOld.z;
				virial[YX] += factor*rijNew.y*rijOld.x;
				virial[YY] += factor*rijNew.y*rijOld.y;
				virial[YZ] += factor*rijNew.y*rijOld.z;
				virial[ZX] += factor*rijNew.z*rijOld.x;
				virial[ZY] += factor*rijNew.z*rijOld.y;
				virial[ZZ] += factor*rijNew.z*rijOld.z; 

				rNew[nm][i] += factor*rijOld*mass_r[id1];
				rNew[nm][next] -= factor*rijOld*mass_r[id2];
                		vtPlus[id1] += factor*rijOld*mass_r[id1]/timeStep;
                		vtPlus[id2] -= factor*rijOld*mass_r[id2]/timeStep;
				moving[nm][i] = true;
				moving[nm][next] = true;
				done = false;
			    }                        
			}
		    }
		    for (Int i = 0; i < nBonds; i++)
		    {
			moved[nm][i] = moving[nm][i];
			moving[nm][i] = false;
		    } 
		    count++;
		}   // complete iterative loop for one molecule
	       

		if(!done)
		{
		    for (Int i = 0; i < ns; i++) // calculate the new assumed positions of atoms on nmth molecule
		    {
		      Atom *atom = mol->myAtoms[i];
		      Int id = atom->atomID;
		      rOld[nm][i] = atom->position;
		      cout<<"the number of atoms "<<id<<endl;
		      cout<<"the position of the atom are x=  "<<rOld[nm][i].x<<endl;
		      cout<<"the position of the atom are y=  "<<rOld[nm][i].y<<endl;
		      cout<<"the position of the atom are z=  "<<rOld[nm][i].z<<endl;
		     }
		    cerr << "SHAKE not converge after " << maxCount << " iterative loops on Molecule: " << nm << endl;
		    ERRORMSG("try large iterative loops");
		   
		}

	    Vector3 rNewCM;
	    Vector3 vNewCM;
	    Vector3 halfMomen;
            for (Int i = 0; i < ns; i++) // calculate the new mol momenta and center of mass on nmth molecule
            {
                Atom *atom = mol->myAtoms[i];
                Int id = atom->atomID;
		rNewCM    += rNew[nm][i] * atom->mass;
		vNewCM    += vtPlus[id]  * atom->mass;
		halfMomen += vtMinus[id] * atom->mass;

	        if(useAtomThermo){
                    atom->velocity = (vtPlus[id] + vtMinus[id])*0.5;}
            }

	    rNewCM = rNewCM/mol->mass;

	    velvect  = (vNewCM + halfMomen)*0.5  - momentime;
	    velcheck = velvect.length();
	    momentime = (vNewCM + halfMomen)*0.5;

	    posvect  = (mol->massCenter + rNewCM)*0.5 - rhalfp;
	    poscheck = posvect.length();
	    rhalfp   = (mol->massCenter + rNewCM)*0.5;


            // if(prank==0){
            //    cout << " NHcount =  "  << NHcount  << "      velcheck =  "<<  velcheck <<  "     poscheck =  "<< poscheck << endl;}
	    NHcount++;

         } // NpT while

         if (NHcount == 10 || velcheck > 1E-05  || poscheck > 1E-05){
            ERRORMSG("NpT barostat integration Error with SHAKE: Check tolerance of both v(t) and r(t+dt/2) or Use another initial Conf. with an approximated Pressure/Volume");
         }


        // update position and velocity for this molecule
        for(Int i = 0; i < ns; i++)
        {
            atom = mol->myAtoms[i];
            Int id = atom->atomID;
            Vector3 dr = rNew[nm][i] - atom->position;
            atom->position += dr;
            atom->realPos += dr;
            atom->displacement += dr;
            atom->velocity = (vtPlus[id] + vtMinus[id])*0.5;
            vtMinus[id] = vtPlus[id];
        }
    }  // NpT ensemble

    }
    
    for (Int i = XX; i <= ZZ; i++)
        virial[i] /= (timeStep*timeStep);
//  }      //Jc: parallel decomposition end
}


double LFIntegrator :: gaussVel(double sigma)
{
	double r, invr, v1, v2, l;
	
	r = 2.0;
	while (r >= 1.0) {
		v1 = 2.0*urandf() - 1.0;
		v2 = 2.0*urandf() - 1.0;
		r = v1*v1 + v2*v2;
	}
	
	invr = 1.0/r;
	l = v1*sqrt(-2.0*invr*log(r));
	l = 0.0 + sigma*l;
	return l;
}

double LFIntegrator :: urandf()
{
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dis(0.0, 1.0);
	// Use dis to transform the random unsigned int generated by gen into a 
        // double in [0, 1). Each call to dis(gen) generates a new random double
        return dis(gen);
}

// Jc: ouput VtMinus for the continue simualtion
void LFIntegrator :: write_vtMinus()
{
//Jc:    if((prank ==0)&&(currentSteps ==100))
    if((prank ==0)&&(currentSteps%100 == 0))
    {
      ofo = new ofstream("vtMinus.out", ios::out);

      *ofo << "vtMinus" << currentSteps << endl;
      *ofo << numAtoms << endl;
      for(int j = 0; j < numAtoms; j++)
      {
        *ofo << vtMinus[j].x;
        *ofo << ' ' << vtMinus[j].y;
        *ofo << ' ' << vtMinus[j].z;
        *ofo << endl;
      }
      *ofo << endl;
      ofo->close();
    }
}

void LFIntegrator :: read_vtMinus()
{       
    int numCoordinates = 0;     
    int i, v1, ret;
    double d1, d2, d3;
    char buf[128];
    double posTemp[numAtoms][3];    

    if(prank == 0)
    {
      FILE* fptr = fopen("vtMinus.in", "r");

      if (fptr == NULL)
         ERRORMSG("fail to open vtMinus.in file" );

      fgets(buf, 512, fptr);             // skip one line

      DEBUGMSG("Reading vtMinus in LF integrator");

      fgets(buf, 512, fptr);

      if(((ret=sscanf(buf,"%d", &numCoordinates))!= 1)||(numCoordinates!= numAtoms))
          ERRORMSG("read the number of vtMinus error \n");
        
      for (i = 0; i < numCoordinates; i++)
      {
        fgets(buf, 512, fptr);
        ret = sscanf(buf, "%lf%lf%lf",  &d1, &d2, &d3);
        posTemp[i][0] = d1;          
        posTemp[i][1] = d2;          
        posTemp[i][2] = d3;          
      }
      fclose(fptr);
    }

//    MPI::COMM_WORLD.Bcast(&numCoordinates,1,MPI::INT,0);
//    MPI::COMM_WORLD.Bcast(posTemp,4500,MPI::DOUBLE,0);
    MPI_Bcast(&numCoordinates,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(posTemp,a3,MPI_DOUBLE,0,MPI_COMM_WORLD);

    for (i = 0; i < numCoordinates; i++)
    {
       vtMinus[i].x = posTemp[i][0];
       vtMinus[i].y = posTemp[i][1];
       vtMinus[i].z = posTemp[i][2];
    }  
}
