/** Simulation.c -- implementation to handle input, maintain the global data structure 
 **     and build the simulation system 
 **
 ** Copyright (C) 2003
 ** Centre for Molecular Simulation (CMS)
 ** School of Information Technology
 ** Swinburne University of Technology
 ** PO Box 21 Hawthorn, Vic 3122, Australia
 **
 ** Author: Zhongwu Zhou
 ** Email: zzhou@it.swin.edu.au
 **/
// JC modified & commented by Jianhui Li, adding output file 
#include "Simulation.h"
#include <mpi.h> // added by Jianhui 

// constructor - initialisation
Simulation::Simulation(const char configFile[])
{   

//  MPI initilization
//    size = MPI::COMM_WORLD.Get_size();
//    rank = MPI::COMM_WORLD.Get_rank();
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    coordinateFile = NULL;     
    sysDataFile = NULL;         // name of file of parameters and structures for molecular system
    velocityFile = NULL;
    resultFile = NULL;
    restartFile = NULL;         // also used as to backup file
    trajectoryFile = NULL;
    dumpFile = NULL;            // have to consider what data will be dumped into it 

    sysConfigFile = configFile;
    myConfig = new SimConfiguration(sysConfigFile); 
    myEnsemble = new Ensemble(); 
    myParams = NULL;
    myIntegrator = NULL;  
    rdf = NULL;

    ofp = NULL; 
    ofpTrajectory = NULL;  
    ofp = new ofstream("result.out", ios::out);
    ofvel = new ofstream("velbehav.out", ios::out);
    ofpr = new ofstream("pressureResult.out", ios::out);
    ofac = new ofstream("pressureAcum.out", ios::out);
    oflu = new ofstream("Lustig.out", ios::out);
    oftp = new ofstream("ThermoProp.out", ios::out);
    ofavl = new ofstream("LustigAverages.out", ios::out);
    ofind = new ofstream("resultInduction.out", ios::out);
//    ofPotential =NULL;                                          // JC initialise the ofp pointer; 
    ofPressureTensor = new ofstream("pTensor.out",ios::out);      // JC initialise the ofp pointer;
 
}

Simulation::~Simulation()
{
    ;
}


void Simulation::setup_simulation()
{
// read basic data from config.txt by process 0 and distribute them to other porcess
// the distribution is carried out in the SimConfiguration.read_config_file();
//    if(rank ==0 ){ 
      myConfig->read_config_file();
//    }
      sysDataFile = myConfig->get_sysData_file();
//      #ifdef DEBUG
//          DEBUGMSG("setup --  simulation -- system  from process 0  ");
//      #endif    
      // sysDataFile contains parameters and topology data for building molecule & atom systems
      // these data will be read in and filled into system data structures
      myParams = new Parameters(sysDataFile);
      myEnsemble->set_mol_system(myParams);
      myEnsemble->set_configuration(myConfig);
      // set up box, boundary, cell manager, cell list and pairlists.
      myEnsemble->compute_bound_box();
      myEnsemble->set_pairlist();
  
      // initialise output
      ofstream of("sysData.out" , ios::out);
      myConfig->write_config(of);
      myEnsemble->set_sysdata_file(&of);
      myEnsemble->set_result_file(ofp);
      myEnsemble->write_ensemble_info(of);
      if(myConfig->write_trajectory())
      {
          trajectoryFile = myConfig->get_trajactory_file();
          if (trajectoryFile == NULL)
              ERRORMSG("no trajectory file name provided");
          ofpTrajectory = new ofstream (trajectoryFile, ios::out);
      }

    // generate Integrator
      if (myConfig->get_integrator_type() == 0)
          myIntegrator = new GearIntegrator(myEnsemble, myConfig);
      else if (myConfig->get_integrator_type() == 1)
          myIntegrator = new VelocityIntegrator(myEnsemble, myConfig);
      else if (myConfig->get_integrator_type() == 2)
          myIntegrator = new LFIntegrator(myEnsemble, myConfig);    
      else if (myConfig->get_integrator_type() == 3)
          myIntegrator = new NHIntegrator(myEnsemble, myConfig);    
 
      if (myIntegrator == NULL)
          ERRORMSG("Null integrator");
      myIntegrator->set_result_file(ofp);

      // generate measuremental objects

      Int nSamples = myConfig->get_n_steps() - myConfig->get_start_sampling();
      if (myConfig->get_sampling_ts_freq() != 0)        
          nSamples /= myConfig->get_sampling_ts_freq(); // in the simConfiguration.cpp
      if (myConfig->compute_rdf())    rdf = new RDF(myEnsemble, nSamples);
     
} 

//void Simulation::run(void) // tray transfer the rank 

  void Simulation::run(void)
{ 
    Int counter = 0;
    Int numSteps = myConfig->get_n_steps();
    Int printSteps =  myConfig->get_print_ts_freq();
    Int samplingSteps = myConfig->get_sampling_ts_freq();
    Int startSteps = 0;

// Jc: MPI initialization 
    int rank, size;
//    rank = MPI::COMM_WORLD.Get_rank();
//    size = MPI::COMM_WORLD.Get_size();
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

// Jc: 

    startSteps = myConfig->get_start_steps();

    if (ofp == NULL)
        ERRORMSG ("null output stream");
    if(rank == 0){  
      *ofp << "Energy" << endl;
      *ofp << "Steps" << '\t' << " L-J" << '\t' << " Angle" << '\t' << " Bond" << '\t' << " UB-hydrogens "<< "Real" << '\t';
      *ofp << "  Long" << '\t' << "  Correct" << '\t' << "  molCorrect" << '\t' << "  surfCorrect" << '\t';
      *ofp << "  AtomTemp" << '\t' << "  MolTemp" << '\t' << "  Total" << endl;
      *ofpr << "## Step  " <<"Atom_Pdiag  " << "shear_stress  "<< "antisymm_SS  " <<"first_normal_Sh  " <<"second_normal_Sh  "<< "Mol_Pdiag  ";
      *ofpr << "shear_stress  "<< "antisymm_SS  " <<"first_normal_Sh  " <<"second_normal_Sh  " << endl;
      *ofac << "## Step  " <<"Atom_Pdiag  " << "shear_stress  "<< "antisymm_SS  " <<"first_normal_Sh  " <<"second_normal_Sh  "<< "Mol_Pdiag  ";
      *ofac << "shear_stress  "<< "antisymm_SS  " <<"first_normal_Sh  " <<"second_normal_Sh  " << endl;
      *ofvel << "# Step  " <<" velcomX  " << "  velcomY  "<< "  velcomZ  " <<"   Diff_vel  " <<"  Carbon_vel  "<< "   Density  ";
      *ofvel <<  " C-O0_dist  " <<  "  stdC-O0   "<<  "   C-O1_dist  " << "   stdC-O1  " << "   Angle  " << "    anglSTDdev  " <<endl;
      *oflu << "# Step  " << "   lusStep    " << " Atomtemp  " << "  MolTemp  " << "  U_LJ    " << "  dUdV_LJ  " << "  d2UdV2_LJ   ";
      *oflu << "      -bdUdV    " << "    (-bdUdV)2  "<< "      -bd2UdV2   " << "     UdUdVF      "  <<"  real2ndDev  " <<"  long2ndDev  "<< endl;
      *ofavl << "# Step  " << "   count   "<<"     bUconf/F  " << "    (bUconf/F)2  " << "    bavU/F  " << "     (bavU/F)2  " << "   Kinav  " << "  1/Kin   ";
      *ofavl << "      Uen  " <<"      Uen2   " << "     Kin2  "<< "      Etot    " << "      Etot2   " << "      Kinatom   "<< "       1/Kinatom   " << "      Kinatom2  " <<endl;
      *oftp << "## Step     " << " U (kJ/mol) "  << " Cv (kJ/molK)  " << "   P (Mpa)" << endl;
      *ofind << "## Step     " << "   Real "  << "   Long  " << "   SelfCorrect  " << "   MolCorrect " << "  SurfacePol   " << "    TotInduction  "<< endl;
    }
    if (startSteps > 0)
    {
       myIntegrator->run(startSteps);
       counter += startSteps;
    }

    myEnsemble->zero_sum();

    if(rank == 0)                                //Jc: by Jianhui to output the system information to pTensor.dat
      myEnsemble->write_systemInfo(*ofPressureTensor);

    while (counter < numSteps)
    {   

        for (Int i = 0; i < printSteps; i += samplingSteps)
        {
            myIntegrator->run(samplingSteps);    // JC sampling steps too large in the system
            myEnsemble->sampling();              // JC what is the function of sampling 
            // RDF sampling
            if (rdf)    rdf->sampling();         // JC what is the function of sampling
            if (rank == 0 )
            {        
              myEnsemble->write_pressureTensor(*ofPressureTensor,counter + i);
            }
        }        

        counter += printSteps;
        if(rank ==0)
        { 
          myEnsemble->write_result(*ofp, counter); // JC: pressure tensor to result.out
//          if(counter ==100)
           myIntegrator->write_state();            // JC: data output is carried out in the Integrator print data every sampling time steps
           //myIntegrator->write_statexyz();
           myEnsemble->write_press_res(*ofpr, counter);
           myEnsemble->write_pressac(*ofac, counter);
           if (myConfig->compute_induction()) {myEnsemble->write_resultInduction(*ofind, counter);}
           if (myConfig->compute_lustig()) {myEnsemble->write_lustig(*oflu, counter);}
           if (myConfig->compute_lustig()) {myEnsemble->write_TPbylustig(*oftp, counter);}
           if (myConfig->compute_lustig()) {myEnsemble->write_averaglus(*ofavl, counter);}
           myEnsemble->write_velbehav(*ofvel, counter);
        }
    }


    if(rank == 0){                        
      myEnsemble->write_sampling_data();     // Jc: write Mol Trajectory and Velocity
    }

}

void Simulation::finish(void) { ; }  


