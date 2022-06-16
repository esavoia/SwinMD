/** VelocityIntegrator.cpp -- 
 **
 ** Copyright (C) 2002
 ** Centre for Molecular Simulation (CMS)
 ** School of Information Technology
 ** Swinburne University of Technology
 ** PO Box 21 Hawthorn, Vic 3122, Australia
 **
 ** Author: Zhongwu Zhou
 ** Email: zzhou@it.swin.edu.au
 **/

#include "VelocityIntegrator.h"



VelocityIntegrator :: VelocityIntegrator(Ensemble* ensemble, SimConfiguration* simConfig) : Integrator(ensemble, simConfig)
{
    #ifdef DEBUG
        DEBUGMSG("Creating VelocityIntegrator");
    #endif

    numUpdates = 0;
    initialise();
}

VelocityIntegrator :: ~VelocityIntegrator() 
{
    if (mass_r != NULL) delete mass_r;
    mass_r = NULL;
}

void VelocityIntegrator :: initialise()
{
    #ifdef DEBUG
        DEBUGMSG("initialise VelocityIntegrator");
    #endif

    Integrator :: initialise();

    mass_r = new Double[numAtoms];
    if (mass_r == NULL)
        ERRORMSG("memory allocation error for VelocityIntegrator");
    for (int i = 0; i < numAtoms; i++)
        mass_r[i] = 1.0/myAtoms[i].mass;

    // if (mySimConfig->is_new_start())
        start();
}

void VelocityIntegrator :: run(Int numTimeSteps)
{
//    Double molCorrectVirial[9]; // never used commented by Jianhui Li            
    bool doUpdate;
    Vector3 position, diff;
    Molecule* mol;
    Atom* atom;
    Double timeStep2 = timeStep*timeStep;

    #ifdef DEBUGXX
        if (dumpFile == NULL)   dumpFile = "dump.out";
        ofstream of = ofstream(dumpFile, ios::out);// JC no 'of' any more
    #endif 

    for (Int i = 0; i < numTimeSteps; i++)
    {
        atomKinEnergy = 0.0;
        molKinEnergy = 0.0;

        currentSteps++;
        #ifdef DEBUG
            cerr << "integration steps: " << currentSteps << endl;
        #endif

        // compute new position
        for (Int j = 0; j< numAtoms; j++)
        {
            Atom* atom = &myAtoms[i];        
            Vector3 dr = atom->velocity*timeStep + 0.5*timeStep2*atom->force*mass_r[i];
            atom->position += dr;
            atom->realPos += dr;
            atom->displacement += dr;
            atom->velocity += 0.5*timeStep*atom->force*mass_r[i];
        }
                               
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
            
            #ifdef DEBUGXX //JC
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
        for (Int j = 0; j < numAtoms; j++)
            myAtoms[j].force = 0.0;
 
        // compute forces at time t+dt
        ljForce->compute();         // assume always have ljForce
        if(has_bond_force)
            bondForce->compute();
        if(has_angle_force)
            angleForce->compute();
        if (computeEwald)        
        {
            ewaldForce->compute();
            if ((i%longSteps) == 0)
                ewaldForce->compute_long();

            // sum long force and virial to the total.
            for (Int j = 0; j < numAtoms; j++)
                myAtoms[j].force += myAtoms[j].longForce;
            for (Int j = XX; j <= ZZ; j++)
                myEnsemble->virial[j] += myEnsemble->longVirial[j];
        }            
        /* if (has_constr_force)             // must after all other forces are computed
        {
            constraintForce->compute();
            velfb = constraintForce->get_velocity_feedback();
            momfb = constraintForce->get_momentum_feedback();         
        } */

        // cmplete velocity computing
        for (int i = 0; i < numAtoms; i++)
            myAtoms[i].velocity += 0.5*timeStep*myAtoms[i].force*mass_r[i];

        // compute molecular centre, momentum, & kinetic energy
        for (Int nm = 0; nm < numMols; nm++)
        {  
            mol = &myMols[nm];

            mol->momenta = 0.0;
            mol->massCenter = 0.0;          
            for (Int j = 0; j < mol->numAtoms; j++)
            {
                atom = mol->myAtoms[j];
                atom->momentum = atom->velocity*atom->mass;
                mol->momenta += atom->momentum;
                mol->massCenter += atom->mass*atom->position;
                atomKinEnergy += atom->velocity*atom->velocity*atom->mass;
            }            
            mol->massCenter /= mol->mass;
            molKinEnergy += mol->momenta*mol->momenta/mol->mass;
        }
        myEnsemble->myEnergy.atomKinEnergy = 0.5*atomKinEnergy;
        myEnsemble->myEnergy.molKinEnergy = 0.5*molKinEnergy;
        molTemp = molKinEnergy/(molNDF*BOLTZMAN);
        atomTemp = atomKinEnergy/(atomNDF*BOLTZMAN);
        
        computePressure();
        thermostat();


        // apply PBC to mass center & then shift position to new center
        for (Int nm = 0; nm < numMols; nm++)
        {  
            mol = &myMols[nm];

            position = mol->massCenter; 
            myEnsemble->apply_pbc(mol->massCenter);        
            diff = mol->massCenter - position;
            for (Int j = 0; j < mol->numAtoms; j++)
                mol->myAtoms[j]->position += diff;
        }
    }
}


void VelocityIntegrator :: thermostat()
{
    Double factor;        // scale factor
      
    if(useAtomThermo)
        factor = sqrt(temperature/atomTemp);
    else
        factor = sqrt(temperature/molTemp);

    for (Int i = 0; i < numAtoms; i++)
        myAtoms[i].velocity *= factor;
}

void VelocityIntegrator :: start()
{
    #ifdef DEBUG
        cerr << "New Start() " << endl;
    #endif
    
    for (Int i = 0; i < numAtoms; i++)
        myAtoms[i].force = 0.0;
    ljForce->compute();
    if (has_angle_force)
        angleForce->compute();
    if (has_bond_force)
        bondForce->compute();
    if (computeEwald )         
    {
        ewaldForce->compute();
        ewaldForce->compute_long();
        for (Int i = 0; i < numAtoms; i++)
            myAtoms[i].force += myAtoms[i].longForce;
    }

    // what constraint should be used?
}


