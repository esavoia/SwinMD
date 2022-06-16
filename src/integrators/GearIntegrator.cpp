/** GearIntegrator.cpp -- 
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

#include "GearIntegrator.h"



GearIntegrator :: GearIntegrator(Ensemble* ensemble, SimConfiguration* simConfig) : Integrator(ensemble, simConfig)
{
    #ifdef DEBUG
        DEBUGMSG("Creating GearIntegrator");
    #endif

    dr1 = NULL;
    dr2 = NULL;
    dr3 = NULL;
    dr4 = NULL;    
    pa1 = NULL;
    pa2 = NULL;
    pa3 = NULL;
    pa4 = NULL;
    pThmo = NULL;
    velfb = NULL;
    momfb = NULL;
    constraintForce = NULL;

    numUpdates = 0;
    init_derivatives = false;
    initialise();
}

GearIntegrator :: ~GearIntegrator()
{
    if(dr1 != NULL)     delete [] dr1;
    if(dr2 != NULL)     delete [] dr2;
    if(dr3 != NULL)     delete [] dr3;
    if(dr4 != NULL)     delete [] dr4;    
    if(pa1 != NULL)     delete [] pa1;
    if(pa2 != NULL)     delete [] pa2;
    if(pa3 != NULL)     delete [] pa3;
    if(pa4 != NULL)     delete [] pa4;
    if(pThmo != NULL) delete [] pThmo;
}

void GearIntegrator :: initialise()
{
    #ifdef DEBUG
        DEBUGMSG("initialise GearIntegrator");
    #endif

    Integrator :: initialise();
    if (has_constr_force)
        constraintForce = new ConstraintForce(myEnsemble, mySimConfig->get_CFB(), mySimConfig->get_DFB());

    dr1 = new Vector3[numAtoms];
    dr2 = new Vector3[numAtoms];
    dr3 = new Vector3[numAtoms];
    dr4 = new Vector3[numAtoms];
    pa1 = new Vector3[numAtoms];
    pa2 = new Vector3[numAtoms];
    pa3 = new Vector3[numAtoms];
    pa4 = new Vector3[numAtoms];
    pThmo = new Vector3[numAtoms];

    if((dr1==NULL)||(dr2==NULL)||(dr3==NULL)||(dr4==NULL)  \
       ||(pa1==NULL)||(pa2==NULL)||(pa3==NULL)||(pa4==NULL)||(pThmo==NULL))
        ERRORMSG("memory allocation error for GearIntegrator");

    if (mySimConfig->is_new_start())
        init_derivatives = true;
    else
    {
        FILE* fptr;
        // if fail to read derivatives for Gear predictor then need initialse the first derivatives. 
        if ((fptr = fopen("derivative.in", "r")) == NULL)
            init_derivatives = true;  
        else
        {
            read_derivative(fptr);
/* JC           ofstream ofp = ofstream("Dout", ios::out);
            ofp << numAtoms << endl;
            for(int j = 0; j < numAtoms; j++)
            {
                ofp << dr1[j].x << ' ' << dr1[j].y << ' ' << dr1[j].z << endl;
                ofp << dr2[j].x << ' ' << dr2[j].y << ' ' << dr2[j].z << endl;
                ofp << dr3[j].x << ' ' << dr3[j].y << ' ' << dr3[j].z << endl;
                ofp << dr4[j].x << ' ' << dr4[j].y << ' ' << dr4[j].z << endl;
            }

            ofp << numAtoms << endl;
            for(int j = 0; j < numAtoms; j++)
            {
                ofp << pa1[j].x << ' ' << pa1[j].y << ' ' << pa1[j].z << endl;
                ofp << pa2[j].x << ' ' << pa2[j].y << ' ' << pa2[j].z << endl;
                ofp << pa3[j].x << ' ' << pa3[j].y << ' ' << pa3[j].z << endl;
                ofp << pa4[j].x << ' ' << pa4[j].y << ' ' << pa4[j].z << endl;
            } */
        }     
            
    }
    if (init_derivatives)
        start();
}

void GearIntegrator :: run(Int numTimeSteps)
{
    Double molCorrectVirial[9];            
    bool doUpdate;
    Vector3 position, diff;
    Molecule* mol;
    Atom* atom;

    #ifdef DEBUGXX // JC change to the DEBUG to DEBUGXX ??
        if (dumpFile == NULL)   dumpFile = "dump.out";
        ofstream of = ofstream(dumpFile, ios::out);
    #endif 

    for (Int i = 0; i < numTimeSteps; i++)
    {
        atomKinEnergy = 0.0;
        molKinEnergy = 0.0;

        currentSteps++;
        #ifdef DEBUG
            cerr << "integration steps: " << currentSteps << endl;
        #endif

        predict();
                        
        // compute molecular centre, momentum, & kinetic energy
        for (Int nm = 0; nm < numMols; nm++)
        {  
            mol = &myMols[nm];

            mol->momenta = 0.0;
            mol->massCenter = 0.0;          
            for (Int j = 0; j < mol->numAtoms; j++)
            {
                atom = mol->myAtoms[j];
                mol->momenta += atom->momentum;
                mol->massCenter += atom->mass*atom->position;
                atomKinEnergy += atom->momentum*atom->momentum/atom->mass;
            }            
            mol->massCenter /= mol->mass;
            molKinEnergy += mol->momenta*mol->momenta/mol->mass;
        }
        myEnsemble->myEnergy.atomKinEnergy = 0.5*atomKinEnergy;
        myEnsemble->myEnergy.molKinEnergy = 0.5*molKinEnergy;
        myEnsemble->atomTemp = atomKinEnergy/(atomNDF*BOLTZMAN);
        myEnsemble->molTemp = molKinEnergy/(molNDF*BOLTZMAN);

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
            
            #ifdef DEBUGXX // JC change from debug to debugxx??
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
 
        // compute forces
        ljForce->compute();         // assume always have ljForce
        if(has_bond_force)
            bondForce->compute();
        if(has_angle_force)
            angleForce->compute();
        //if (has_elec_force ) 
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
        if (has_constr_force)             // must after all other forces are computed
        {
            constraintForce->compute();
            velfb = constraintForce->get_velocity_feedback();
            momfb = constraintForce->get_momentum_feedback();         
        }

        computePressure();
        myEnsemble->accumulate();
        gaussThermostat(pThmo, useAtomThermo);
        correct();

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

// we use first order motion equation & scaled time derivatives for Gear method
// dr1[i] = deltat*(dr[i]/dt), dr2[i] = 1/2 * delta2 * (d2r[i]/dt2), ...
void GearIntegrator :: predict()
{
    for (Int i = 0; i < numAtoms; i++)
    {
        Atom* atom = &myAtoms[i];
        // predict new position
        Vector3 dr = dr1[i] + dr2[i] + dr3[i] + dr4[i];
        atom->position += dr;
        atom->realPos += dr;
        atom->displacement += dr;
        // predict new velocity
        dr1[i] += 2.0*dr2[i] + 3.0*dr3[i] + 4.0*dr4[i];
        // predict new acceleration
        dr2[i] += 3.0*dr3[i] + 6.0*dr4[i];
        // predict the third derivative
        dr3[i] += 4.0*dr4[i];
        // the fourth derivative d4[i] = d4[i]

        // momenta
        atom->momentum += pa1[i] + pa2[i] + pa3[i] + pa4[i];
        pa1[i] += 2.0*pa2[i] + 3.0*pa3[i] + 4.0*pa4[i];
        pa2[i] += 3.0*pa3[i] + 6.0*pa4[i];
        pa3[i] += 4.0*pa4[i];
        // predicted velocity will be used for computing constraint force & correct()
        atom->velocity = atom->momentum/atom->mass;
    }
}

void GearIntegrator :: correct()
{
    for (Int i = 0; i < numAtoms; i++)
    {
        Atom* atom = &myAtoms[i];

        Vector3 pa1Cor = (atom->force - pThmo[i])*timeStep - pa1[i];
        if (has_constr_force)       // feedback from constraint
            pa1Cor -= momfb[i]*timeStep;
        atom->momentum += K10*pa1Cor;
        atom->velocity = atom->momentum/atom->mass;
        pa1[i] += pa1Cor;
        pa2[i] += K12*pa1Cor;
        pa3[i] += K13*pa1Cor;
        pa4[i] += K14*pa1Cor;

        Vector3 dr1Cor = atom->velocity*timeStep - dr1[i];
        if (has_constr_force)       // feedback from constraint
            dr1Cor -= velfb[i]*timeStep;
        atom->position += K10*dr1Cor;
        atom->realPos += K10*dr1Cor;
        atom->displacement += K10*dr1Cor;
        dr1[i] += dr1Cor;
        dr2[i] += K12*dr1Cor;
        dr3[i] += K13*dr1Cor;
        dr4[i] += K14*dr1Cor;
    }
}


/* void GearIntegrator :: correct()
{
    for (Int i = 0; i < numAtoms; i++)
    {
        Atom* atom = &myAtoms[i];
        Vector3 dr1Cor = atom->velocity*timeStep - dr1[i];
        Vector3 pa1Cor = (atom->force - pThmo[i])*timeStep - pa1[i];
        
        if (has_constr_force)       // feedback for constraint
        {
            dr1Cor -= velfb[i]*timeStep;
            pa1Cor -= momfb[i]*timeStep;
        }

        atom->position += K10*dr1Cor;
        atom->realPos += K10*dr1Cor;
        atom->displacement += K10*dr1Cor;
        dr1[i] += dr1Cor;
        dr2[i] += K12*dr1Cor;
        dr3[i] += K13*dr1Cor;
        dr4[i] += K14*dr1Cor;

        atom->momentum += K10*pa1Cor;
        pa1[i] += pa1Cor;
        pa2[i] += K12*pa1Cor;
        pa3[i] += K13*pa1Cor;
        pa4[i] += K14*pa1Cor;
    }
}
*/

 


void GearIntegrator :: gaussThermostat(Vector3* pThmo, bool atomThermo)
{
    Double zeta;        // thermostat multiplier
    Double sumPF = 0.0;
    Double molTemp;
  
    if(atomThermo)
    {
        Atom* atom;
        for (Int i = 0; i < numAtoms; i++)
        {
            atom = &myAtoms[i];
            sumPF += atom->momentum*atom->force/atom->mass;
        }
        molTemp = atomKinEnergy/(atomNDF*BOLTZMAN);
        if (atomKinEnergy == 0.0)
            ERRORMSG("Devided by zero");
        zeta = sumPF/atomKinEnergy + TFB*(molTemp - temperature);    // include feedback for temperature
        #ifdef DEBUG
            cerr << "Zeta, T, targetT " << zeta << "; " << molTemp << "; " << temperature << endl;
        #endif
        for (Int i = 0; i < numAtoms; i++)
            pThmo[i] = zeta*myAtoms[i].momentum;
    }

    else                // using molecular thermostat
    {
        Int id;
        Molecule *mol;
        for (Int i = 0; i < numMols; i++)
        {  
            mol = &myMols[i];
            mol->force = 0.0;  
            for (Int j = 0; j < mol->numAtoms; j++)
                mol->force += mol->myAtoms[j]->force;
            sumPF += mol->momenta*mol->force/mol->mass;
        }
        molTemp = molKinEnergy/(molNDF*BOLTZMAN);

        if (molKinEnergy == 0.0)
            ERRORMSG("Devided by zero");
        zeta = sumPF/molKinEnergy + TFB*(molTemp - temperature);    // include feedback for temperature
        #ifdef DEBUG
            cerr << "Zeta, molT, targetT " << zeta << "; " << molTemp << "; " << temperature << endl;
        #endif

        for (Int i = 0; i < numAtoms; i++)
        {
            id = myAtoms[i].molID;
            mol = &myMols[id];
            pThmo[i] = zeta*myAtoms[i].mass/mol->mass*mol->momenta;
        }
    }
    
}

void GearIntegrator :: start()
{
    #ifdef DEBUG
        cerr << "New Start() " << endl;
    #endif
    
    for (Int i = 0; i < numAtoms; i++)
        myAtoms[i].force = 0.0;
    ljForce->compute();
    if (has_angle_force)
    {
        cerr << "compute angle" << endl;
        angleForce->compute();
    }
    if (has_bond_force)
    {
        cerr << "compute bond" << endl;
        bondForce->compute();
    }
    if (computeEwald )         
    {
        ewaldForce->compute();
        ewaldForce->compute_long();
        for (Int i = 0; i < numAtoms; i++)
            myAtoms[i].force += myAtoms[i].longForce;
    }
    for (Int i = 0; i < numAtoms; i++)
    {
        dr1[i] = myAtoms[i].velocity*timeStep;
        pa1[i] = myAtoms[i].force*timeStep;
    }        
    if (has_constr_force)            
    {
        constraintForce->compute();
        velfb = constraintForce->get_velocity_feedback();
        momfb = constraintForce->get_momentum_feedback();
        for (Int i = 0; i < numAtoms; i++)
        {
            dr1[i] -= velfb[i]*timeStep;
            pa1[i] -= momfb[i]*timeStep;
        }
    }
}

void GearIntegrator :: read_derivative(FILE* fptr)
{
    Int nAtoms, ret;
    Double vx, vy, vz;
    char buf[256];

    fgets(buf, 512, fptr);        // skip one line
    fgets(buf, 512, fptr); 
    if((ret = sscanf(buf, "%d", &nAtoms))!=1 || (nAtoms != numAtoms))
    {
        DEBUGMSG("reading atom number error" );
        return;
    }
    for(int j = 0; j < nAtoms; j++)
    {
        fgets(buf, 512, fptr);
        sscanf(buf, "%lf%lf%lf", &vx, &vy, &vz);
        dr1[j].x = vx;      dr1[j].y = vy;       dr1[j].z = vz;
        fgets(buf, 512, fptr);
        sscanf(buf, "%lf%lf%lf", &vx, &vy, &vz);
        dr2[j].x = vx;      dr2[j].y = vy;       dr2[j].z = vz;
        fgets(buf, 512, fptr);
        sscanf(buf, "%lf%lf%lf", &vx, &vy, &vz);
        dr3[j].x = vx;      dr3[j].y = vy;       dr3[j].z = vz;
        fgets(buf, 512, fptr);
        sscanf(buf, "%lf%lf%lf", &vx, &vy, &vz);
        dr4[j].x = vx;      dr4[j].y = vy;       dr4[j].z = vz;
    }

    fgets(buf, 512, fptr);        // skip one line
    fgets(buf, 512, fptr); 
    if((ret = sscanf(buf, "%d", &nAtoms))!=1 || (nAtoms != numAtoms))
    {
        DEBUGMSG ( "reading atom number error" );
        return;
    }
    for(int j = 0; j < nAtoms; j++)
    {
        fgets(buf, 512, fptr);
        sscanf(buf, "%lf%lf%lf", &vx, &vy, &vz);
        pa1[j].x = vx;      pa1[j].y = vy;       pa1[j].z = vz;
        fgets(buf, 512, fptr);
        sscanf(buf, "%lf%lf%lf", &vx, &vy, &vz);
        pa2[j].x = vx;      pa2[j].y = vy;       pa2[j].z = vz;
        fgets(buf, 512, fptr);
        sscanf(buf, "%lf%lf%lf", &vx, &vy, &vz);
        pa3[j].x = vx;      pa3[j].y = vy;       pa3[j].z = vz;
        fgets(buf, 512, fptr);
        sscanf(buf, "%lf%lf%lf", &vx, &vy, &vz);
        pa4[j].x = vx;      pa4[j].y = vy;       pa4[j].z = vz;
    }
    
}


void GearIntegrator :: write_state()
{
    Integrator :: write_state();

    // write Gear predictor derivatives to file
    /* JC ofstream ofp = ofstream("derivative.out", ios::out);
    ofp << "velocity derivatives at steps: " << currentSteps << endl;
    ofp << numAtoms << endl;
    for(int j = 0; j < numAtoms; j++)
    {
        ofp << dr1[j].x << ' ' << dr1[j].y << ' ' << dr1[j].z << endl;
        ofp << dr2[j].x << ' ' << dr2[j].y << ' ' << dr2[j].z << endl;
        ofp << dr3[j].x << ' ' << dr3[j].y << ' ' << dr3[j].z << endl;
        ofp << dr4[j].x << ' ' << dr4[j].y << ' ' << dr4[j].z << endl;
    }
    ofp << "momentum derivatives at steps: " << currentSteps << endl;
    ofp << numAtoms << endl;
    for(int j = 0; j < numAtoms; j++)
    {
        ofp << pa1[j].x << ' ' << pa1[j].y << ' ' << pa1[j].z << endl;
        ofp << pa2[j].x << ' ' << pa2[j].y << ' ' << pa2[j].z << endl;
        ofp << pa3[j].x << ' ' << pa3[j].y << ' ' << pa3[j].z << endl;
        ofp << pa4[j].x << ' ' << pa4[j].y << ' ' << pa4[j].z << endl;
    } */

}

