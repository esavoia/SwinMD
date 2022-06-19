/** Ensemble.c --
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

#include "Ensemble.h"
#include "Errors.h" 

Ensemble::Ensemble()
{
    DEBUGMSG("Creating Ensemble object");

    nAtoms         = 0;
    nAtomTypes     = 0;
    nMols          = 0;
    nMolTypes      = 0;
    nBonds         = 0;
    nBondTypes     = 0;
    nAngles        = 0;
    nAngleTypes    = 0;
    nDihedrals     = 0;
    nDihedralTypes = 0;
    nImpropers     = 0;
    nImproperTypes = 0;

    atomStructs    = NULL;
    molStructs     = NULL;
    bonds          = NULL;
    angles         = NULL;
    dihedrals      = NULL;
    impropers      = NULL;
    exclusionTable = NULL;

    // atom and molecule objects
    atoms          = NULL;
    molecules      = NULL;

    boxLx  = 0.0;      boxLy  = 0.0;       boxLz  = 0.0;
    halfLx = 0.0;      halfLy = 0.0;       halfLz = 0.0;
    rBoxLx = 0.0;      rBoxLy = 0.0;       rBoxLz = 0.0;
    volume = 0.0;      cellLen = 0.0;
    mySeed = 0.0;

    myConfig = NULL;
    myParams = NULL;
    myCellManager = NULL;
    sysdataFile = NULL;

    time = NULL;
    numPerType = NULL;
    trajectory = NULL;
    velocity   = NULL;

    myEnergy.ljEnergy = 0.0;
    myEnergy.realEnergy = 0.0;
    myEnergy.longEnergy = 0.0;
    myEnergy.correctEnergy = 0.0;
    myEnergy.molCorrectEnergy = 0.0;
    myEnergy.surfCorrectEnergy = 0.0;
    myEnergy.angleEnergy = 0.0;
    myEnergy.bondEnergy = 0.0;
    myEnergy.atomKinEnergy = 0.0;
    myEnergy.molKinEnergy = 0.0;
    myEnergy.miscEnergy = 0.0;
    myEnergy.totEnergy = 0.0;
    myEnergy.inducEnergy = 0.0;
    myEnergy.TBEnergy = 0.0;
    myEnergy.ubEnergy = 0.0;
    myEnergy.TotTorq = 0.0;

    myEnergy.realPolEnergy = 0.0;
    myEnergy.longPolEnergy = 0.0;
    myEnergy.correctPolEnergy = 0.0;
    myEnergy.molCorrectPolEnergy = 0.0;
    myEnergy.surfacePolEnergy = 0.0;
    myLustig.intravirPolCorrect = 0.0;
    myLustig.longPolDeriv = 0.0;
    myLustig.longSecPolDeriv = 0.0;
    myLustig.realPolDeriv = 0.0;
    myLustig.realSecPolDeriv = 0.0;
    myLustig.molPolDeriv = 0.0;
    myLustig.molSecPolDeriv = 0.0;
    myLustig.surfPolDeriv = 0.0;
    myLustig.surfSecPolDeriv = 0.0;
    myLustig.ubDeriv = 0.0;
    myLustig.ubSecDeriv = 0.0;


    myLustig.ljUlus   = 0.0;
    myLustig.ljdUlus  = 0.0;
    myLustig.ljd2Ulus = 0.0;
    myLustig.virialcorrect = 0.0;
    myLustig.intravirCorrect = 0.0;
    myLustig.longDeriv = 0.0;
    myLustig.longSecDeriv = 0.0;
    myLustig.realDeriv = 0.0;
    myLustig.realSecDeriv = 0.0;
    myLustig.molDeriv = 0.0;
    myLustig.molSecDeriv = 0.0;
    myLustig.surfDeriv = 0.0;
    myLustig.surfSecDeriv = 0.0;
    myLustig.inducEn = 0.0;
    myLustig.inducDeriv = 0.0;
    myLustig.BondDeriv = 0.0;
    myLustig.BondSecDeriv  = 0.0;
    myLustig.selfPolDeriv  = 0.0;

    atomTemp = 0.0;
    molTemp = 0.0;

    counter = 0;
    luscount = 0.0;
    sampleCounter = 0;

    zero_sum();

   /// rank = MPI::COMM_WORLD.Get_rank();  Updated to MPI-3 by elton
   /// size = MPI::COMM_WORLD.Get_size();

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


}

Ensemble::~Ensemble()
{
    if (exclusionTable != NULL) delete [] exclusionTable;
    if (atoms != NULL)  delete [] atoms;
    if (molecules != NULL)  delete [] molecules;
}

void Ensemble::set_mol_system(Parameters* params)
{
    Int i, type, id, num;
    Double charge, mass;
    char* str, typeName;

    #ifdef DEBUG
        DEBUGMSG("setting up molecular system");
    #endif

    myParams = params;
    if(myParams == NULL)
        ERRORMSG("fail to access parameter object");

    nAtoms         = myParams->nAtoms;
    nAtomTypes     = myParams->nAtomTypes;
    nMols          = myParams->nMolecules;
    nMolTypes      = myParams->nMolTypes;
    nBonds         = myParams->nBonds;
    nBondTypes     = myParams->nBondTypes;
    nAngles        = myParams->nAngles;
    nAngleTypes    = myParams->nAngleTypes;
    nDihedrals     = myParams->nDihedrals;
    nDihedralTypes = myParams->nDihedralTypes;
    nImpropers     = myParams->nImpropers;
    nImproperTypes = myParams->nImproperTypes;

    atomStructs    = myParams->atomArray;
    molStructs     = myParams->molArray;
    bonds          = myParams->bondArray;
    angles         = myParams->angleArray;
    dihedrals      = myParams->dihedralArray;
    impropers      = myParams->improperArray;

    nA3   = 3*nAtoms;
cout << " 3*nAtom = " << nA3 << endl;
    atoms = new Atom[nAtoms];
    molecules = new Molecule[nMols];

    if (atoms == NULL || molecules == NULL)
        ERRORMSG("memory allocation of atoms or moleclues error");

    // fill atom info to atom object array
    for (i = 0; i < nAtoms; i++)
    {
        type = atomStructs[i].atomType;
        charge = myParams->atomParams[type].charge;

        atoms[i].atomID = atomStructs[i].atomID;    // should be same as i
        atoms[i].atomType = type;
        atoms[i].molID = atomStructs[i].molID;
        atoms[i].scaledCharge = charge*SQRTCOULOMBCONSTANT;
        atoms[i].dipx = myParams->atomParams[type].dipx;
        atoms[i].dipy = myParams->atomParams[type].dipy;
        atoms[i].dipz = myParams->atomParams[type].dipz;
        atoms[i].quadxx = myParams->atomParams[type].quadxx;
        atoms[i].quadxy = myParams->atomParams[type].quadxy;
        atoms[i].quadxz = myParams->atomParams[type].quadxz;
        atoms[i].quadyy = myParams->atomParams[type].quadyy;
        atoms[i].quadyz = myParams->atomParams[type].quadyz;
        atoms[i].quadzz = myParams->atomParams[type].quadzz;
        atoms[i].polar  = myParams->atomParams[type].polar;
        atoms[i].mass = myParams->atomParams[type].mass;
        atoms[i].typeName = myParams->atomParams[type].typeName;
    }

    // fill molecular info to molecule objects
    for (i=0; i < nMols; i++)
    {
        type = molStructs[i].molType;
        num = myParams->molTypes[type].numOfAtoms;
        str = myParams->molTypes[type].name;

        molecules[i].set_molecule(i, type, num, str);
    }

    // set atoms and mass to molecules
    for (i = 0; i < nAtoms; i++)
    {
        type = atoms[i].atomType;
        id = atoms[i].molID;

        molecules[id].mass += atoms[i].mass;
        molecules[id].add_atom(&atoms[i]);
    }

    build_exclusion_table();
}

void Ensemble::set_configuration(SimConfiguration* config)
{
    char dump[256];

    myConfig = config;
    mySeed = myConfig->get_seed();
    nSamples = myConfig->get_n_steps() - myConfig->get_start_sampling();
    interval = myConfig->get_sampling_ts_freq();
    totsteps = myConfig->get_n_steps();
    if(interval != 0)
        nSamples /= interval;
    deltaT = myConfig->get_timestep();
    computeMSD = myConfig->compute_msd();
    computeVACF = myConfig->compute_vacf();
    temperature = myConfig->get_temperature(); // Jc: for output the temperature data;
    if (computeMSD || computeVACF)
    {
        time = new Double[nSamples + 1];
        numPerType = new Int[nMolTypes];
        if (time == NULL || numPerType == NULL)
            ERRORMSG("Error in memory allocation");

        for (Int i = 0; i < nMols; i++)
        {
            Int typeID = molecules[i].molType;
            numPerType[typeID]++;
        }
    }

    if (computeMSD)
    {
        trajectory = new Vector3* [nSamples + 1];
        if (trajectory == NULL)
            ERRORMSG("Error in memory allocation");
        for (Int j = 0; j <= nSamples; j++)
        {
            trajectory[j] = new Vector3[nAtoms];//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            if (trajectory[j] == NULL)
                ERRORMSG("Fail to allocate memory for MSD");
        }
    }
    if (computeVACF)
    {
        velocity = new Vector3* [nSamples + 1];
        if (velocity == NULL)
            ERRORMSG("Error in memory allocation");
        for (Int j = 0; j <= nSamples; j++)
        {
            velocity[j] = new Vector3[nMols];
            if (velocity[j] == NULL)
                ERRORMSG("Fail to allocate memory for MSD");
        }
    }

    if (myConfig ->is_new_start())
    {
        char*  coordinateFile = myConfig->get_coordinate_file();
        FILE*  fptrc = fopen(coordinateFile, "r");
        if (fptrc == NULL)
            ERRORMSG("fail to open coordinate file");
        read_coordinate(fptrc);
        if (myConfig->is_read_velocity())
        {
            char* velocityFile = myConfig->get_velocity_file();
            FILE* fptrv = fopen(velocityFile, "r");
            if (fptrv == NULL)
                ERRORMSG("fail to open velocity file");
            read_velocity(fptrv);
            fclose(fptrv);
        }
        // else random_velocity(mySeed);
        else init_velocity(); // Jc: in this calculation, velocity are assigned according to Gaussian
        fclose(fptrc);
    }
    else
    {
        // char* fname = myConfig->get_restart_file();
        FILE* fptr = fopen("state.in", "r");
        if (fptr == NULL)
            ERRORMSG("fail to open restart file" );

        fgets(dump, 512, fptr);             // skip one line
        read_coordinate(fptr);
        fgets(dump, 512, fptr);             // skip one line
        read_realPos(fptr);
        fgets(dump, 512, fptr);             // skip one line
        read_velocity(fptr);
        fclose(fptr);
    }
}

// Currently, exclusion table is built basing on bonded atom pairs without consideration of 1-4 exclusion.
// For atom i, exclusionTable[i] is an array containing all other atoms directly bonded with it and angle as well.
// Therefor, an exclusion pair for atomI and atomJ is duplicated in the table.
void Ensemble::build_exclusion_table()
{
    Int a1, a2;

    DEBUGMSG("building exclusion table");

    // reforming bond list in the order of atom ID
    vector < set<int> > bondlist(nAtoms);
    for (Int i = 0; i < nBonds; i++)
    {
        a1=bonds[i].atom1;
        a2=bonds[i].atom2;
        bondlist[a1].insert(a2);
        bondlist[a2].insert(a1);
    }

    // angle between i-j-k. i-j and j-k are excluded as bonded pairs as above.
    // So for excluding angle we also need to exclude i-k.
    vector < set<int> > anglelist(nAtoms);
    for (Int i = 0; i < nAngles; i++)
    {
        a1 = angles[i].atom1;
        a2 = angles[i].atom3;
        bondlist[a1].insert(a2);
        bondlist[a2].insert(a1);
    }

    // allocate memory for exclusionTable
    exclusionTable = new Int* [nAtoms];
    for (int i = 0; i < nAtoms; i++)
    {
        unsigned int size = bondlist[i].size();   // + anglelist[i].size();
        exclusionTable[i] = new Int[size+1];
        exclusionTable[i][size] = -1;
    }

    // build the table
    for(Int i = 0; i < nAtoms; i++)
    {
        Int j = 0;
        for (set<int>::iterator atom=bondlist[i].begin();atom!=bondlist[i].end();atom++)
            exclusionTable[i][j++] = *atom;
    }
    // need free bondlist??
}

// set dimensions of the simulation box from atom's positions and scale the box to a density if specified
void Ensemble::compute_bound_box()
{
    Int i;
    Double boxLen = 0.0;
    Vector3   boxMin, boxMax;

    if (nAtoms <= 0)
        return;

    boxMin = atoms[0].position;
    boxMax = atoms[0].position;

    for (i = 1; i < nAtoms; i++)
    {
        boxMin.x = MIN(boxMin.x, atoms[i].position.x);
        boxMin.y = MIN(boxMin.y, atoms[i].position.y);
        boxMin.z = MIN(boxMin.z, atoms[i].position.z);

        boxMax.x = MAX(boxMax.x, atoms[i].position.x);
        boxMax.y = MAX(boxMax.y, atoms[i].position.y);
        boxMax.z = MAX(boxMax.z, atoms[i].position.z);
    }

    /*
    // shift atoms to middle of the box
    Double xshift = (boxMax.x + boxMin.x)/2;
    Double yshift = (boxMax.y + boxMin.y)/2;
    Double zshift = (boxMax.z + boxMin.z)/2;
    for (i = 0; i < nAtoms; i++)
    {
        atoms[i].position.x -= xshift;
        atoms[i].position.y -= yshift;
        atoms[i].position.z -= zshift;
    }
    */

    boxLx = boxMax.x - boxMin.x;
    boxLy = boxMax.y - boxMin.y;
    boxLz = boxMax.z - boxMin.z;

    // volume = boxLx*boxLy*boxLz;

    if ((density = myConfig->get_density()) > 0)
    {
        Double vol = nAtoms/density;
        // default cubic box
        boxLen = pow(vol, 1.0/3.0);
        // boxLen = MAX(MAX(boxLx, boxLy), MAX(boxLz, boxLen));
        boxLx = boxLen;
        boxLy = boxLen;
        boxLz = boxLen;

        // using rectangle box that matchs with configration
        /* if (vol > volume)             // scaling the box to match density
        {
            Double factor = pow(vol/volume, 1.0/3.0);
            boxLx *= factor;
            boxLy *= factor;
            boxLz *= factor;
        } */
    }
    volume = boxLx*boxLy*boxLz;
    density = nAtoms/volume;

    halfLx = 0.5*boxLx;
    halfLy = 0.5*boxLy;
    halfLz = 0.5*boxLz;
    if(boxLx > 0.0) rBoxLx = 1.0/boxLx;
    if(boxLy > 0.0) rBoxLy = 1.0/boxLy;
    if(boxLz > 0.0) rBoxLz = 1.0/boxLz;
    myConfig->set_box(boxLx, boxLy, boxLz);
    myConfig->set_density(density);
//    cout<<"the boxLx is : "<<boxLx<<endl;
//    cout<<"the boxLy is : "<<boxLy<<endl;
//    cout<<"the boxLz is : "<<boxLz<<endl;
}


// set_pairlist() -- a wraped function to set cell number and dimmensions,
//                 -- generate cell list, set index, origin and x, y, z ranks,
//                 -- build neighbor cell list for each cell
//                 -- and build pair list for each atom
// This function usually called at the starting of simulation
// but can also be called at running of the simulation if the dimension of sim-box is changed
void Ensemble::set_pairlist()
{
    pairCutOff =  myConfig->get_cutoff()+myConfig->get_cutBuff();
    pairCutOff *= pairCutOff;
    useCellPairList = myConfig->use_cell_list();
    if (useCellPairList)
    {
        if (myCellManager == NULL)
            myCellManager = new CellManager();
        cellLen = myConfig->get_cutoff() + myConfig->get_cutBuff();
        // should reassign 'useCellPairList' flag in case the box is too small for using cell method
        useCellPairList = myCellManager->set_cellManager(cellLen, boxLx, boxLy, boxLz);
    }
    if(useCellPairList)         // continue to set cell list & cell neighbor list
    {
        myCellManager->build_cell_list();
        myCellManager->build_cell_neighbor();
    }
    update_pairlist();
}

void Ensemble::update_pairlist()
{
    // clear & reinitialise pairlist for each atom
    for(Int i = 0; i < nAtoms; i++)
        atoms[i].clear_pairlist();

    if (useCellPairList)
    {
        myCellManager->clear_atomlist();
        for(Int i = 0; i < nAtoms; i++)
            myCellManager->set_atom(&atoms[i]);
        myCellManager->update_pair_list(exclusionTable, pairCutOff);
    }
    else        // build pair list using brute force method
        build_pairlist();
}

void Ensemble::build_pairlist()
{
    Atom* atomi;
    Atom* atomj;
    Vector3 rij;

    // pair lists built with following loop are uneven lists
    // when parallelise the code have to pay attention on the load balance problem
    for (Int i = 0; i < (nAtoms - 1); i++)
    {
        atomi = &atoms[i];
        #ifdef DEBUG
            if (atomi == NULL)
                ERRORMSG("NULL atom pointer");
        #endif
        for(Int j = i+1; j < nAtoms; j++)
        {
            atomj = &atoms[j];
            #ifdef DEBUG
            if (atomj == NULL)
                ERRORMSG("NULL atom pointer");
            #endif
            if (!exclusion_check(atomi->atomID, atomj->atomID))
            {
                rij = atomj->position - atomi->position;
                apply_pbc(rij);
                if(rij.length2() < pairCutOff)
                {
                    atomi->set_pair(atomj);
                }
            }
        }
    }
}

void Ensemble::read_coordinate(char* fname)  // Jc: never used for the coordinate reading
{
    Int numCoordinates = 0;
    Int i, v1, ret;
    Double d1, d2, d3;
    char buf[128];
    FILE *fptr;
    double posTemp[nAtoms][3];   //Jc: temporary storage needed for reading coordinates in Paralle computating

    if(rank ==0)
    {
      DEBUGMSG("Reading coordinate using the first method");

      if ((fptr = fopen(fname, "r")) == NULL)
          ERRORMSG("open coordinate file error");

      fgets(buf, 512, fptr);
      if(((ret=sscanf(buf,"%d", &numCoordinates))!= 1)||(numCoordinates!= nAtoms))
          ERRORMSG("read the number of coordiantes error \n");
      for (i = 0; i < numCoordinates; i++)
      {
          fgets(buf, 512, fptr);
          // ret = sscanf(buf, "%d%lf%lf%lf", &v1, &d1, &d2, &d3);
          ret = sscanf(buf, "%lf%lf%lf",  &d1, &d2, &d3);
          posTemp[i][0] = d1;
          posTemp[i][1] = d2;
          posTemp[i][2] = d3;

//ZW          atoms[i].position.x = d1;
//ZW          atoms[i].position.y = d2;
//ZW          atoms[i].position.z = d3;
//ZW          atoms[i].realPos = atoms[i].position;
      }
      fclose(fptr);
    }
// broadcast coordinates from rank 0
    //MPI::COMM_WORLD.Bcast(&numCoordinates,1,MPI::INT,0);   Updated to MPI-3 by elton
    //MPI::COMM_WORLD.Bcast(posTemp,4500,MPI::DOUBLE,0);
    MPI_Bcast(&numCoordinates,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(posTemp,nA3,MPI_DOUBLE,0,MPI_COMM_WORLD);
    for (i = 0; i < numCoordinates; i++)
    {

      atoms[i].position.x = posTemp[i][0];
      atoms[i].position.y = posTemp[i][1];
      atoms[i].position.z = posTemp[i][2];
      atoms[i].realPos = atoms[i].position;
    }

}

void Ensemble::read_velocity(char* fname)  // Jc: never used for the velocity reading
{
    Int numVelocities = 0;
    Int i, v1, ret;
    Double d1, d2, d3;
    char buf[128];
    FILE *fptr;
    double posTemp[nAtoms][3];   //Jc: temporary storage needed in reading velocity

    if(rank ==0)
    {
      DEBUGMSG("Reading velocity from the first method");

      if ((fptr = fopen(fname, "r")) == NULL)
          ERRORMSG("open velocity file error");

      fgets(buf, 512, fptr);
      if(((ret=sscanf(buf,"%d", &numVelocities))!= 1)||(numVelocities!= nAtoms))
          ERRORMSG("read the number of velocities error \n");

      for (i = 0; i < numVelocities; i++)
      {
        fgets(buf, 512, fptr);
//Jc;        ret = sscanf(buf, "%d%lf%lf%lf", &v1, &d1, &d2, &d3);
        ret = sscanf(buf, "%lf%lf%lf",  &d1, &d2, &d3);

        posTemp[i][0] = d1;
        posTemp[i][1] = d2;
        posTemp[i][2] = d3;

//ZW        atoms[i].velocity.x = d1;
//ZW        atoms[i].velocity.y = d2;
//ZW        atoms[i].velocity.z = d3;
//ZW        atoms[i].momentum = atoms[i].mass*atoms[i].velocity;
      }
    }
// broadcast volecities from rank 0
    //MPI::COMM_WORLD.Bcast(&numVelocities,1,MPI::INT,0); Updated to MPI-3 by elton
    //MPI::COMM_WORLD.Bcast(posTemp,4500,MPI::DOUBLE,0);
    MPI_Bcast(&numVelocities,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(posTemp,nA3,MPI_DOUBLE,0,MPI_COMM_WORLD);
    for (i = 0; i < numVelocities; i++)
    {
      atoms[i].velocity.x = posTemp[i][0];
      atoms[i].velocity.y = posTemp[i][1];
      atoms[i].velocity.z = posTemp[i][2];
      atoms[i].momentum = atoms[i].mass*atoms[i].velocity;
    }
}

//JC:  below methods used in the data reading

void Ensemble::read_coordinate(FILE* fptr)  // Jc: used for the coordinate reading
{
    Int numCoordinates = 0;
    Int i, v1, ret;
    Double d1, d2, d3;
    char buf[128];
    double posTemp[nAtoms][3];

    if(rank == 0)
    {
      DEBUGMSG("Reading coordinate from the 2nd method");

      fgets(buf, 512, fptr);
      if(((ret=sscanf(buf,"%d", &numCoordinates))!= 1)||(numCoordinates!= nAtoms))
          ERRORMSG("read the number of coordiantes error \n");

      for (i = 0; i < numCoordinates; i++)
      {
        fgets(buf, 512, fptr);
//ZW        ret = sscanf(buf, "%d%lf%lf%lf", &v1, &d1, &d2, &d3);
        ret = sscanf(buf, "%lf%lf%lf",  &d1, &d2, &d3);
        posTemp[i][0] = d1;
        posTemp[i][1] = d2;
        posTemp[i][2] = d3;
//ZW        atoms[i].position.x = d1;
//ZW        atoms[i].position.y = d2;
//ZW        atoms[i].position.z = d3;
        // atoms[i].position = atoms[i].position0;

      }
    }
// broadcast coordinates from rank 0
    //MPI::COMM_WORLD.Bcast(&numCoordinates,1,MPI::INT,0);    Updated to MPI-3 by elton
    //MPI::COMM_WORLD.Bcast(posTemp,4500,MPI::DOUBLE,0);
    MPI_Bcast(&numCoordinates,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(posTemp,nA3,MPI_DOUBLE,0,MPI_COMM_WORLD);
    for (i = 0; i < numCoordinates; i++)
    {
      atoms[i].position.x = posTemp[i][0];
      atoms[i].position.y = posTemp[i][1];
      atoms[i].position.z = posTemp[i][2];
      atoms[i].realPos = atoms[i].position; // Jc: only effective for new start simulation
    }
}

// Jc: created by Jianhui Li to read the real position from state.in
void Ensemble::read_realPos(FILE* fptr)  // Jc: used for the reaPos reading
{
    Int numCoordinates = 0;
    Int i, v1, ret;
    Double d1, d2, d3;
    char buf[128];
    double posTemp[nAtoms][3];

    if(rank == 0)
    {
      DEBUGMSG("Reading realPos from the 2nd method");

      fgets(buf, 512, fptr);
      if(((ret=sscanf(buf,"%d", &numCoordinates))!= 1)||(numCoordinates!= nAtoms))
          ERRORMSG("read the number of coordiantes error \n");

      for (i = 0; i < numCoordinates; i++)
      {
        fgets(buf, 512, fptr);
        ret = sscanf(buf, "%lf%lf%lf",  &d1, &d2, &d3);
        posTemp[i][0] = d1;
        posTemp[i][1] = d2;
        posTemp[i][2] = d3;
      }
    }
// broadcast coordinates from rank 0
    //MPI::COMM_WORLD.Bcast(&numCoordinates,1,MPI::INT,0);   Updated to MPI-3 by elton
    //MPI::COMM_WORLD.Bcast(posTemp,4500,MPI::DOUBLE,0);
    MPI_Bcast(&numCoordinates,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(posTemp,nA3,MPI_DOUBLE,0,MPI_COMM_WORLD);
    for (i = 0; i < numCoordinates; i++)
    {
      atoms[i].realPos.x = posTemp[i][0];
      atoms[i].realPos.y = posTemp[i][1];
      atoms[i].realPos.z = posTemp[i][2];
    }
}

void Ensemble::read_velocity(FILE* fptr)  // Jc: used in the velocity reading in continue simulation
{
    Int numVelocities = 0;
    Int i, v1, ret;
    Double d1, d2, d3;
    char buf[128];
    double posTemp[nAtoms][3];              // Jc: temporary storage needed in reading velocity

    if(rank ==0)
    {

      DEBUGMSG("Reading velocity from the 2nd method");

      fgets(buf, 512, fptr);
      if(((ret=sscanf(buf,"%d", &numVelocities))!= 1)||(numVelocities!= nAtoms))
        ERRORMSG("read the numbert of velocities error \n");

      for (i = 0; i < numVelocities; i++)
      {
        fgets(buf, 512, fptr);
//ZW       ret = sscanf(buf, "%d%lf%lf%lf", &v1, &d1, &d2, &d3);
        ret = sscanf(buf, "%lf%lf%lf", &d1, &d2, &d3);
        posTemp[i][0] = d1;
        posTemp[i][1] = d2;
        posTemp[i][2] = d3;

//ZW        atoms[i].velocity.x = d1;
//ZW        atoms[i].velocity.y = d2;
//ZW        atoms[i].velocity.z = d3;
//ZW        atoms[i].momentum = atoms[i].mass*atoms[i].velocity;
      }
    }

   // MPI::COMM_WORLD.Bcast(&numVelocities,1,MPI::INT,0);  Updated to MPI-3 by elton
   // MPI::COMM_WORLD.Bcast(posTemp,4500,MPI::DOUBLE,0);
   MPI_Bcast(&numVelocities,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(posTemp,nA3,MPI_DOUBLE,0,MPI_COMM_WORLD);
    for (i = 0; i < numVelocities; i++)
    {
      atoms[i].velocity.x = posTemp[i][0];
      atoms[i].velocity.y = posTemp[i][1];
      atoms[i].velocity.z = posTemp[i][2];
      atoms[i].momentum = atoms[i].mass*atoms[i].velocity;
    }
}


void Ensemble::set_sysdata_file(ofstream *ofp)
{
    sysdataFile = ofp;
}

void Ensemble::set_result_file(ofstream *ofp)
{
    resultFile = ofp;
}

void Ensemble :: zero_sum()
{
    counter = 0;
    energy = 0.0;
    atomStress.pressure = 0.0;
    atomStress.shearStress = 0.0;
    atomStress.antiShearStress = 0.0;
    atomStress.firstNormStress = 0.0;
    atomStress.secondNormStress = 0.0;
    molStress.pressure = 0.0;
    molStress.shearStress = 0.0;
    molStress.antiShearStress = 0.0;
    molStress.firstNormStress = 0.0;
    molStress.secondNormStress = 0.0;
    myAvLus.bUF = 0.0;
    myAvLus.bUF2 = 0.0;
    myAvLus.bUFt = 0.0;
    myAvLus.bUF2t = 0.0;
    myAvLus.Kinav = 0.0;
    myAvLus.invKin = 0.0;
    myAvLus.Uen = 0.0;
    myAvLus.Uen2 = 0.0;
    myAvLus.Kin2 = 0.0;
    myAvLus.Entot = 0.0;
    myAvLus.Entot2 = 0.0;
    myAvLus.Kinavat = 0.0;
    myAvLus.invKinat = 0.0;
    myAvLus.Kin2at = 0.0;
    myAvLus.bdUdv = 0.0;
    myAvLus.bdUdvsq = 0.0;
    myAvLus.bd2Udv2 = 0.0;
    myAvLus.UdUdvF = 0.0;
    myAvLus.Volnm = 0.0;
    myAvLus.Volsq = 0.0;
    myAvLus.Hconf = 0.0;
    myAvLus.Hconf2 = 0.0;
    myAvLus.HV = 0.0;
    myAvLus.Hinter = 0.0;
    myAvLus.Hinter2 = 0.0;
    myAvLus.HintV = 0.0;
    myAvLus.P = 0.0;
}

void Ensemble :: accumulate() // Jc: calculate the Average stress tensor
{
    counter++;
    temperature = myConfig->get_temperature();
    double Tpress = myConfig->get_pressure() / 1.66053904042716;
    //double Tpress = myAvLus.P;
    double voll = boxLx*boxLy*boxLz;
    double vol2 = voll*voll;
    double Upot, UelecPerm, UelecPol, Etot;
    double Hconf, Uinter, Hinter;

    if (counter > 0.15*totsteps)
    {

    int molNDF;
        if(myConfig->is_constraint_on()){
    	    molNDF = 6*nMols - 3;
	}
	else {
	    molNDF = 3*nAtoms - 3;
	}
    luscount++;

    UelecPerm = myEnergy.realEnergy + myEnergy.longEnergy + myEnergy.correctEnergy + myEnergy.molCorrectEnergy + myEnergy.surfCorrectEnergy;
    UelecPol  = myEnergy.realPolEnergy + myEnergy.longPolEnergy + myEnergy.correctPolEnergy + myEnergy.molCorrectPolEnergy + myEnergy.surfacePolEnergy;
    Upot      = UelecPerm + UelecPol + myEnergy.ljEnergy + myEnergy.angleEnergy + myEnergy.bondEnergy + myEnergy.ubEnergy;
    Uinter    = UelecPerm + UelecPol + myEnergy.ljEnergy;
    Hconf     = Upot + Tpress*voll;
    Hinter    = Upot + myEnergy.atomKinEnergy + Tpress*voll;

    energy += Upot;

    myAvLus.Hconf  += Hconf;
    myAvLus.Hconf2 += Hconf*Hconf;
    myAvLus.HV     += Hconf*voll;

    myAvLus.Hinter  += Hinter;
    myAvLus.Hinter2 += Hinter*Hinter;
    myAvLus.HintV   += Hinter*voll;

    myAvLus.Uen  += (myEnergy.ljEnergy + myEnergy.TBEnergy + myEnergy.inducEnergy + myEnergy.angleEnergy + myEnergy.bondEnergy + myEnergy.ubEnergy + myEnergy.realEnergy + myEnergy.longEnergy + myEnergy.correctEnergy + myEnergy.molCorrectEnergy + myEnergy.surfCorrectEnergy);
    myAvLus.Uen2 += (myEnergy.ljEnergy + myEnergy.TBEnergy + myEnergy.inducEnergy + myEnergy.angleEnergy + myEnergy.bondEnergy + myEnergy.ubEnergy + myEnergy.realEnergy + myEnergy.longEnergy + myEnergy.correctEnergy + myEnergy.molCorrectEnergy + myEnergy.surfCorrectEnergy)*(myEnergy.ljEnergy + myEnergy.TBEnergy + myEnergy.inducEnergy + myEnergy.angleEnergy + myEnergy.bondEnergy + myEnergy.ubEnergy + myEnergy.realEnergy + myEnergy.longEnergy + myEnergy.correctEnergy + myEnergy.molCorrectEnergy + myEnergy.surfCorrectEnergy);

    myAvLus.Entot  += (myEnergy.ljEnergy + myEnergy.TBEnergy + myEnergy.inducEnergy + myEnergy.angleEnergy + myEnergy.bondEnergy + myEnergy.ubEnergy + myEnergy.realEnergy + myEnergy.longEnergy + myEnergy.correctEnergy + myEnergy.molCorrectEnergy + myEnergy.surfCorrectEnergy + myEnergy.atomKinEnergy);
    myAvLus.Entot2 += (myEnergy.ljEnergy + myEnergy.TBEnergy + myEnergy.inducEnergy + myEnergy.angleEnergy + myEnergy.bondEnergy + myEnergy.ubEnergy + myEnergy.realEnergy + myEnergy.longEnergy + myEnergy.correctEnergy + myEnergy.molCorrectEnergy + myEnergy.surfCorrectEnergy + myEnergy.atomKinEnergy)*(myEnergy.ljEnergy + myEnergy.TBEnergy + myEnergy.inducEnergy + myEnergy.angleEnergy + myEnergy.bondEnergy + myEnergy.ubEnergy + myEnergy.realEnergy + myEnergy.longEnergy + myEnergy.correctEnergy + myEnergy.molCorrectEnergy + myEnergy.surfCorrectEnergy + myEnergy.atomKinEnergy);

    myAvLus.bUF += Upot/(temperature*BOLTZMAN*molNDF*0.5);
    myAvLus.bUF2 += Upot*Upot/(temperature*BOLTZMAN*molNDF*0.5)/(temperature*BOLTZMAN*molNDF*0.5);

    myAvLus.bUFt  += (myLustig.ljUlus + myEnergy.TBEnergy + myEnergy.inducEnergy + myEnergy.angleEnergy + myEnergy.bondEnergy + myEnergy.ubEnergy + myEnergy.realEnergy + myEnergy.longEnergy + myEnergy.correctEnergy + myEnergy.molCorrectEnergy + myEnergy.surfCorrectEnergy)/(temperature*BOLTZMAN*molNDF*0.5);
    myAvLus.bUF2t += (myLustig.ljUlus + myEnergy.TBEnergy + myEnergy.inducEnergy + myEnergy.angleEnergy + myEnergy.bondEnergy + myEnergy.ubEnergy + myEnergy.realEnergy + myEnergy.longEnergy + myEnergy.correctEnergy + myEnergy.molCorrectEnergy + myEnergy.surfCorrectEnergy)/(temperature*BOLTZMAN*molNDF*0.5)*(myLustig.ljUlus + myEnergy.TBEnergy + myEnergy.inducEnergy + myEnergy.angleEnergy + myEnergy.bondEnergy + myEnergy.ubEnergy + myEnergy.realEnergy + myEnergy.longEnergy + myEnergy.correctEnergy + myEnergy.molCorrectEnergy + myEnergy.surfCorrectEnergy)/(temperature*BOLTZMAN*molNDF*0.5);

    myAvLus.Kinav += myEnergy.molKinEnergy;
    myAvLus.invKin += 1/myEnergy.molKinEnergy;
    myAvLus.Kin2 += (myEnergy.molKinEnergy)*(myEnergy.molKinEnergy);
    myAvLus.Kinavat += myEnergy.atomKinEnergy;
    myAvLus.invKinat += 1/myEnergy.atomKinEnergy;
    myAvLus.Kin2at += (myEnergy.atomKinEnergy)*(myEnergy.atomKinEnergy);
    myAvLus.Volnm += voll;
    myAvLus.Volsq += vol2;

    myAvLus.bdUdv   +=   -(myLustig.ljdUlus + myLustig.inducDeriv - 3.0*myEnergy.TBEnergy/voll + myLustig.surfDeriv + myLustig.molDeriv + myLustig.surfPolDeriv + myLustig.molPolDeriv + (myLustig.longDeriv + myLustig.realDeriv + myLustig.longPolDeriv + myLustig.realPolDeriv + myLustig.selfPolDeriv + myLustig.ubDeriv + myLustig.BondDeriv)/(3.0*voll))/(temperature*BOLTZMAN) ;
    myAvLus.bdUdvsq += ( -(myLustig.ljdUlus + myLustig.inducDeriv - 3.0*myEnergy.TBEnergy/voll + myLustig.surfDeriv + myLustig.molDeriv + myLustig.surfPolDeriv + myLustig.molPolDeriv + (myLustig.longDeriv + myLustig.realDeriv + myLustig.longPolDeriv + myLustig.realPolDeriv + myLustig.selfPolDeriv + myLustig.ubDeriv + myLustig.BondDeriv)/(3.0*voll))/(temperature*BOLTZMAN) )*( -(myLustig.ljdUlus + myLustig.inducDeriv - 3.0*myEnergy.TBEnergy/voll + myLustig.surfDeriv + myLustig.molDeriv + myLustig.surfPolDeriv + myLustig.molPolDeriv + (myLustig.longDeriv + myLustig.realDeriv + myLustig.longPolDeriv + myLustig.realPolDeriv + myLustig.selfPolDeriv + myLustig.ubDeriv + myLustig.BondDeriv)/(3.0*voll))/(temperature*BOLTZMAN) );

    myAvLus.bd2Udv2 += -( myLustig.ljd2Ulus + 12.0*myEnergy.TBEnergy/vol2 + myLustig.surfSecDeriv + myLustig.molSecDeriv + (myLustig.longSecDeriv + myLustig.realSecDeriv + myLustig.BondSecDeriv)/(9.0*vol2) )/(temperature*BOLTZMAN) ;

    myAvLus.UdUdvF  += (myLustig.ljUlus + myEnergy.TBEnergy + myEnergy.inducEnergy + myEnergy.angleEnergy + myEnergy.bondEnergy + myEnergy.ubEnergy + myEnergy.realEnergy + myEnergy.longEnergy + myEnergy.correctEnergy + myEnergy.molCorrectEnergy + myEnergy.surfCorrectEnergy)/(temperature*BOLTZMAN*molNDF*0.5) * ( -(myLustig.ljdUlus + myLustig.inducDeriv - 3.0*myEnergy.TBEnergy/voll + myLustig.surfDeriv + myLustig.molDeriv + (myLustig.longDeriv + myLustig.realDeriv + myLustig.selfPolDeriv + myLustig.ubDeriv + myLustig.BondDeriv)/(3.0*voll))/(temperature*BOLTZMAN) );

    atomStress.pressure += myAvLus.P;
    atomStress.shearStress += (atomPressure[XY] + atomPressure[YX])*0.5;
    atomStress.antiShearStress += (atomPressure[YX] - atomPressure[XY])*0.5;
    atomStress.firstNormStress += (atomPressure[XX] - atomPressure[YY]);
    atomStress.secondNormStress += (atomPressure[YY] - atomPressure[ZZ]);
    molStress.pressure += (molPressure[XX]+molPressure[YY]+molPressure[ZZ])/3.0;
    molStress.shearStress += (molPressure[XY] + molPressure[YX])*0.5;
    molStress.antiShearStress += (molPressure[YX] - molPressure[XY])*0.5;
    molStress.firstNormStress += (molPressure[XX] - molPressure[YY]);
    molStress.secondNormStress += (molPressure[YY] - molPressure[ZZ]);
//    cout<<" the calculated Atom pressure ZZ is from Esemble.cpp *** "<<atomPressure[ZZ]<<endl;
    }
}

// the sampling() function simply stores molecular trajectory and/or velocity data for
// MSD/VACF which will be dumped into files at the end of simulation for later analysis
void Ensemble :: sampling()
{
    if (sampleCounter > nSamples)
        return;
    if (computeMSD == NULL && computeVACF == NULL)
        return;

/*    for (Int i = 0; i < nMols; i++)
    {
        Vector3 Rmt, Vmt;
        for (Int j = 0; j < molecules[i].numAtoms; j++)
        {
            Rmt += molecules[i].myAtoms[j]->realPos*molecules[i].myAtoms[j]->mass;
            Vmt += molecules[i].myAtoms[j]->momentum;
        }
        if (computeMSD)
            trajectory[sampleCounter][i] = Rmt/molecules[i].mass;
        if (computeVACF)
            velocity[sampleCounter][i] = Vmt/molecules[i].mass;
    }
*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(Int i = 0; i < nAtoms; i++)
    {
        Vector3 Rmt, Vmt;

        if (computeMSD)
            trajectory[sampleCounter][i] = atoms[i].realPos;
        if (computeVACF)
            velocity[sampleCounter][i] = Vmt/molecules[i].mass;
    }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    time[sampleCounter] = sampleCounter*interval*deltaT + time[0];
    sampleCounter++;
}

void Ensemble :: write_ensemble_info(ofstream &of)
{

    DEBUGMSG("write ensemble info");

    char* tab = "       ";

    write_sys_data(of);
    of << endl;

    of << "***Box dimension***" << endl;
    of << boxLx << tab << boxLy << tab << boxLz << tab << endl;
    of <<"half box length: "<< halfLx << tab << halfLy << tab << halfLz  << endl;
    of << "density: " << density << endl;
    of << endl;

    #ifdef DEBUG
    of << "*** Initial coordinates of atoms ***"<<endl;
    of << nAtoms << endl;
    for (int i = 0; i < nAtoms; i++)
        of << atoms[i].position.x << tab << atoms[i].position.y << tab << atoms[i].position.z << endl;
    of << endl;
    of << " *** Initial atom velocity ***" << endl;
    of << nAtoms << endl;
    for (Int i = 0; i < nAtoms; i++)
        of << atoms[i].velocity.x << tab << atoms[i].velocity.y << tab << atoms[i].velocity.z << endl;
    of << endl;
    of << "*** Initial bond lenght *** " << endl;
    for (Int i = 0; i < nMols; i++)
        write_bond(&molecules[i], of);
    write_pairlist(of);
    #endif

    DEBUGMSG("end of write ensemble info");
}

void Ensemble::write_pairlist(ofstream &of)
{
    of<< "***atom pair list***" << endl;
    for (Int i = 0; i < nAtoms; i++)
    {
        of << "atom[" << i << "]: ";
        for(Int j = 0; j < atoms[i].get_list_size(); j++)
        {
            of << atoms[i].myPairList[j]->atomID << ' ' ;
        }
        of << endl;
    }
    of << endl;
}

void Ensemble::write_sys_data(ofstream& of)
{
    char* tab = "       ";

    of << "!sysData file " << endl;
    of << "!Size_of_system" << endl;
    of << nMols << tab << nMolTypes << "       ! number of molecules and number of molecule types" << endl;
    of << nAtoms << tab << nAtomTypes << "       ! number of atoms and number of atom types" << endl;
    of << nBonds << tab << nBondTypes << "         ! number of bonds and number of bond types" << endl;
    of << nAngles << tab << nAngleTypes << "        ! number of angles and number of angle types" << endl;
    of << nDihedrals << tab << nDihedralTypes << "        ! number of dihedrals and number of dihedrals types" << endl;
    of << nImpropers << tab << nImproperTypes << "        ! number of impropers and number of improper types" << endl;
    of << endl;

    of << "!AtomParams          ! section header " << endl;
    of << "!typeID  Mass   Charge   Epsilon  Sigma  typeName  atomName" << endl;
    for (int i = 0; i < nAtomTypes; i++)
    {
        of<<myParams->atomParams[i].atomType<<tab<<myParams->atomParams[i].mass<<tab<<myParams->atomParams[i].charge;
        of<<tab<<myParams->atomParams[i].eps<<tab<<myParams->atomParams[i].sigma<<tab<<myParams->atomParams[i].typeName<<endl;
    }
    of << endl;

    of << "!BondParams" << endl;
    of << "!bondTypeID  atomType1  atomType2  bondLength  forceConstant" << endl;
    for (int i = 0; i < nBondTypes; i++)
    {
        of<<myParams->bondParams[i].bondType<<tab<<myParams->bondParams[i].atomType1<<tab<<myParams->bondParams[i].atomType2;
        of<<tab<<myParams->bondParams[i].r0<<tab<<myParams->bondParams[i].k<<endl;
    }
    of << endl;

    of << "!AngleParams" << endl;
    of << "!angTypeID    aType1    aType2    aType3    degree    K" << endl;
    for (int i = 0; i < nAngleTypes; i++)
    {
        double val = myParams->angleParams[i].theta0*180/PI;
        of<<myParams->angleParams[i].angleType<<tab<<myParams->angleParams[i].atomType1;
        of<<tab<<myParams->angleParams[i].atomType2<<tab<<myParams->angleParams[i].atomType3;
        of<<tab<<myParams->angleParams[i].theta0*180/PI<<tab<<myParams->angleParams[i].k<<endl;
    }
    of << endl;

    of << "!DihedralParams" << endl;
    of << "!TypeID       aType1    aType2    aType3    aType4   K   d   n" << endl;
    of << endl;
    of << "!ImproperParams " << endl;
    of << "!TypeID       aType1    aType2    aType3    aType4   K   d   n " << endl;
    of << endl;

    of << "!MolTypes" << endl;
    of << "!molTypeID       numOfAtoms      molNames" << endl;
    for (int i = 0; i < nMolTypes; i++)
    {
        of<<myParams->molTypes[i].molType<<tab<<myParams->molTypes[i].numOfAtoms;
        of<<tab<<myParams->molTypes[i].name<<endl;
    }
    of << endl;

    of << "!Atoms                    following are topology data" << endl;
    of << "!atomID    atomType      molID       molType        " << endl;
    for (int i = 0; i < nAtoms; i++)
        of<<atoms[i].atomID<<tab<<atoms[i].atomType<<tab<<atoms[i].molID<<endl;
    of << endl;

    of << "!Bonds" << endl;
    of << "!bondID  atom1          atom2        bondType" << endl;
    for (int i = 0; i < nBonds; i++)
        of<<i<<tab<<bonds[i].atom1<<tab<<bonds[i].atom2<<tab<<bonds[i].bondType<< endl;
    of << endl;

    of << "!Angles" << endl;
    of << "!angleID  atom1   atom2    atom3    angleType " << endl;
    for (int i = 0; i < nAngles; i++)
        of<<i<<tab<<angles[i].atom1<<tab<<angles[i].atom2<<tab<<angles[i].atom3<<tab<<angles[i].angleType<<endl;
    of << endl;

    of << "!Dihedrals" << endl;
    of << "!dihedralID atom1      atom2     atom3      atom4      dihedralType" << endl;
    of << endl;
    of << "!Impropers" << endl;
    of << "!ImproperID atom1      atom2     atom3      atom4      dihedralType" << endl;
    of << endl;

    of << "LJ parameter table" << endl;
    of << "atomType1    atomType2    Epsilon  Sigma" << endl;
    for (int i = 0; i < myParams->ljTableSize*myParams->ljTableSize; i++)
    {
        of<<myParams->ljParamTable[i].atomType1<<tab<<myParams->ljParamTable[i].atomType2<<tab;
        of<<myParams->ljParamTable[i].eps<<tab<<myParams->ljParamTable[i].sigma<<endl;
    }
    of << endl;

    of << "!Dipoles" << endl;
    for (int i = 0; i < nAtomTypes; i++)
    {
        of<<myParams->atomParams[i].atomType<<tab<< myParams->atomParams[i].typeName <<tab<< myParams->atomParams[i].dipx;
        of<<tab<<myParams->atomParams[i].dipy<<tab<<myParams->atomParams[i].dipz<<endl;
    }
    of << endl;

    of << "!Quadrupoles" << endl;
    for (int i = 0; i < nAtomTypes; i++)
    {
        of<<myParams->atomParams[i].atomType<<tab<< myParams->atomParams[i].typeName <<tab<< myParams->atomParams[i].quadxx;
        of<<tab<<myParams->atomParams[i].quadxy<<tab<<myParams->atomParams[i].quadxz<<endl;
        of<<tab<< tab<< tab<< myParams->atomParams[i].quadyy<<tab<<myParams->atomParams[i].quadyz<< tab<< myParams->atomParams[i].quadzz <<endl;
    }
    of << endl;

    of << "exclusion table" << endl;
    for (int i = 0; i < nAtoms; i++)
    {
        of << "atom[" << i << "]: ";
        int j = 0;
        while(exclusionTable[i][j] >= 0)
        {
            of << exclusionTable[i][j] << ' ' ;
                j++;
        }
        of << endl;
    }

}

void Ensemble::write_result(ofstream& of, Int numSteps) // Jc: output pressure to the <result.out>
{
    // instant values
    of << numSteps;
    of << fixed << setprecision(6) << setw(14) << myEnergy.ljEnergy;
    of << fixed << setprecision(6) << setw(14) << myEnergy.angleEnergy;
    of << fixed << setprecision(6) << setw(14) << myEnergy.bondEnergy;
    of << fixed << setprecision(6) << setw(14) << myEnergy.ubEnergy;
    of << fixed << setprecision(6) << setw(16) << myEnergy.realEnergy;
    of << fixed << setprecision(6) << setw(14) << myEnergy.longEnergy;
    of << fixed << setprecision(5) << setw(16) << myEnergy.correctEnergy;
    of << fixed << setprecision(5) << setw(14) << myEnergy.molCorrectEnergy;
    of << scientific << setprecision(6) << setw(14) << myEnergy.surfCorrectEnergy;
    of << fixed << setprecision(6) << setw(14) << atomTemp;
    of << fixed << setprecision(6) << setw(14) << molTemp;
//    myEnergy.totEnergy = myEnergy.ljEnergy + myEnergy.realEnergy + myEnergy.longEnergy + myEnergy.correctEnergy;
//    myEnergy.totEnergy += myEnergy.molCorrectEnergy + myEnergy.surfCorrectEnergy;
    of << fixed << setprecision(6) << setw(16) << myEnergy.totEnergy;
    if (counter > 0)            // average energy
        of << fixed << setprecision(6) << setw(16) << energy/luscount << endl;

}

void Ensemble::write_resultInduction(ofstream& of, Int numSteps) // Output for resulInduction.out
{
    // instant values
    of << numSteps;
    of << fixed << setprecision(6) << setw(16) << myEnergy.realPolEnergy;
    of << fixed << setprecision(6) << setw(14) << myEnergy.longPolEnergy;
    of << fixed << setprecision(5) << setw(16) << myEnergy.correctPolEnergy;
    of << fixed << setprecision(5) << setw(14) << myEnergy.molCorrectPolEnergy;
    of << scientific << setprecision(6) << setw(14) << myEnergy.surfacePolEnergy;
    of << fixed << setprecision(6) << setw(16) <<   myEnergy.realPolEnergy + myEnergy.longPolEnergy + myEnergy.correctPolEnergy + myEnergy.molCorrectPolEnergy + myEnergy.surfacePolEnergy;
    of << fixed << setprecision(6) << setw(16) << myEnergy.inducEnergy;
    of << fixed << setprecision(6) << setw(16) << myEnergy.TotTorq;
    of << fixed << setprecision(6) << setw(16) << myLustig.surfPolDeriv;
    of << fixed << setprecision(6) << setw(16) << myLustig.molPolDeriv;
    of << fixed << setprecision(6) << setw(16) << myLustig.selfPolDeriv;
    of << fixed << setprecision(6) << setw(16) << myLustig.inducDeriv;
    if (counter > 0)            // average energy
        of << fixed << setprecision(6) << setw(16) << energy/luscount << endl;
}

void Ensemble :: write_press_res(ofstream& ofpr, Int numSteps)
{
    ofpr<< numSteps << setw(14)<<(atomPressure[XX]+atomPressure[YY]+atomPressure[ZZ])/3.0<<setw(14)<<(atomPressure[XY] + atomPressure[YX])*0.5;
    ofpr<<setw(14)<<(atomPressure[YX] - atomPressure[XY])*0.5<<setw(14)<<(atomPressure[XX] - atomPressure[YY]) <<setw(14)<<(atomPressure[YY] - atomPressure[ZZ]);
    ofpr<<setw(14)<<(molPressure[XX]+molPressure[YY]+molPressure[ZZ])/3.0<<setw(14)<<(molPressure[XY] + molPressure[YX])*0.5;
    ofpr<<setw(14)<<(molPressure[YX] - molPressure[XY])*0.5<<setw(14)<<(molPressure[XX] - molPressure[YY])<<setw(14)<<(molPressure[YY] - molPressure[ZZ]);
    ofpr << endl;
}

void Ensemble :: write_pressac(ofstream& of, Int numSteps)
{
    if (luscount > 0 )
{
    of<< numSteps << setw(14)<<atomStress.pressure/luscount<<setw(14)<<atomStress.shearStress/luscount;
    of<<setw(14)<<atomStress.antiShearStress/luscount <<setw(14)<<atomStress.firstNormStress/luscount <<setw(14)<<atomStress.secondNormStress/luscount;
    of<<setw(14)<<molStress.pressure/luscount <<setw(14)<<molStress.shearStress/luscount;
    of<<setw(14)<<molStress.antiShearStress/luscount  <<setw(14)<<molStress.firstNormStress/luscount<<setw(14)<<molStress.secondNormStress/luscount;
    of << endl;
}
}

void Ensemble :: write_lustig(ofstream& of, Int numSteps)
{
    of<< numSteps << "   " <<luscount;
    of<< scientific << setprecision(8) << setw(17) << myLustig.ljdUlus;
    of<< scientific << setprecision(8) << setw(17) << myLustig.inducDeriv;
    of<< scientific << setprecision(8) << setw(17) << myLustig.surfDeriv;       // please sum here future lustig averages for other potentials
    of<< scientific << setprecision(8) << setw(17) << myLustig.molDeriv;      // please sum here future lustig first derivate averages for other potentials
    of<< scientific << setprecision(8) << setw(17) << myLustig.surfPolDeriv;     // please sum here future lustig second derivate  averages for other potentials
    of<< scientific << setprecision(8) << setw(17) << myLustig.molPolDeriv;
    of<< scientific << setprecision(8) << setw(17) << myLustig.longDeriv/(3.0*volume);
    of<< scientific << setprecision(8) << setw(17) << myLustig.realDeriv/(3.0*volume);
    of<< scientific << setprecision(8) << setw(17) << myLustig.ubDeriv/(3.0*volume);
    of<< scientific << setprecision(8) << setw(17) << myLustig.BondDeriv/(3.0*volume);
    of<< scientific << setprecision(8) << setw(17) << myLustig.longPolDeriv/(3.0*volume);
    of<< scientific << setprecision(8) << setw(17) << myLustig.realPolDeriv/(3.0*volume);
    of<< scientific << setprecision(8) << setw(17) << myLustig.selfPolDeriv/(3.0*volume);
    of<< scientific << setprecision(8) << setw(17) << 2.0*myEnergy.atomKinEnergy/(3.0*volume);
    of<< endl;
}

void Ensemble :: write_TPbylustig(ofstream& of, Int numSteps)
{
    double M10, M20, M01, M02, M11;
    double bb, Uu, Cv, Pp, Pvir, Cvens;
    double Cp, ap, Bt, Bs, w2, w;
    double Cpint, Cvint, Bsint, w2int, wint;
    double Hconf, Hconf2, HV, avV, avV2, Z01_bt;
    double Hint, Hint2, HintV, Ti;
    double voll, vol2, HV_HV, HV_HVconf, volFluc;
    temperature = myConfig->get_temperature();
    int statEnsem = myConfig->get_ensemble_status();
    bb = 1.0/(temperature*BOLTZMAN);
    int NDF, Nn;
    if(myConfig->is_constraint_on()){
        NDF = 6*nMols - 3;
	Nn  = nMols;	// with nmols Lustig will use the TranslationalKinetic energy
    }
    else {
        NDF = 3*nAtoms - 3;
	Nn  = nAtoms;
    }

    if (luscount > 0 ) {
	voll = boxLx*boxLy*boxLz;
	vol2 = voll*voll;

	//---------------------------------------------------------------------------
	//  Let's compute the derivatives of the partition
	//  function accodring to Eq.(12) of the NVT Lustig method
	//---------------------------------------------------------------------------
	M10 = -0.5*NDF/bb*(1.0 + myAvLus.bUFt/luscount);
	M20 = (-0.5*NDF/bb)*(-0.5*NDF/bb)*( 1.0 + 2.0/NDF + 2.0*myAvLus.bUFt/luscount + myAvLus.bUF2t/luscount );
	M01 = Nn/voll + myAvLus.bdUdv/luscount;

        if(myConfig->is_constraint_on()){
            Pvir = 2.0*myAvLus.Kinav/luscount/(3.0 * voll) + myAvLus.bdUdv/luscount/bb;
        }
        else {
            Pvir = 2.0*myAvLus.Kinavat/luscount/(3.0 * voll) + myAvLus.bdUdv/luscount/bb;
        }

	Cvens = myAvLus.Entot2/luscount - ( myAvLus.Entot/luscount * myAvLus.Entot/luscount );
	Cvint = myAvLus.Uen2/luscount - ( myAvLus.Uen*myAvLus.Uen/luscount/luscount );

	//---------------------------------------------------------------------------
	//  Now compute the thermodynamic properties
	//---------------------------------------------------------------------------
	Ti = 2.0*myAvLus.Kinavat/NDF/BOLTZMAN/luscount;
	Uu = myAvLus.Uen/luscount/nMols;		     // in kJ/mol
	Cv = BOLTZMAN*bb*bb*(M20 - M10*M10)/nMols*1000;      // in J/mol.K
	Pp = 1.66053904042716*atomStress.pressure/luscount;  // in MPa
	Cvens = Cvens/(temperature*temperature*BOLTZMAN*nMols)*1000;
	Cvint = Cvint/(temperature*temperature*BOLTZMAN*nMols)*1000;   // configurational Cv

        //---------------------------------------------------------------------------
        //  Compute the thermodynamic properties in the NPT ensamble
        //---------------------------------------------------------------------------
	if (statEnsem == 1){

            Hconf  = myAvLus.Hconf/luscount;
            Hconf2 = myAvLus.Hconf2/luscount;
            HV     = myAvLus.HV/luscount;
            avV    = myAvLus.Volnm/luscount;
            avV2   = myAvLus.Volsq/luscount;

    	    Hint  = myAvLus.Hinter/luscount;
	    Hint2 = myAvLus.Hinter2/luscount;
	    HintV = myAvLus.HintV/luscount;

            M10 = -0.5*NDF/bb - Hconf;
            M20 = 0.5*NDF*(0.5*NDF + 1.0)/(bb*bb) + NDF*Hconf/bb + Hconf2;
            M01 = -bb*avV;
            M02 = bb*bb*avV2;
            M11 = (0.5*NDF - 1.0)*avV + bb*HV;

	    volFluc = avV2 - avV*avV;
            Z01_bt = M01 - bb*(M11 - M10*M01);
	    HV_HVconf = HV - Hconf*avV;
            Cp  = (Hconf2 - Hconf*Hconf)/(BOLTZMAN*Ti*Ti);

	    if (volFluc == 0) {
                Cv  = Cp;
	    } else {
                Cv  = Cp - HV_HVconf*HV_HVconf/volFluc/(BOLTZMAN*Ti*Ti);
	    }
	    ap  = HV_HVconf/avV/(BOLTZMAN*Ti*Ti);
            Bt  = volFluc/(avV*BOLTZMAN*Ti);
            Bs  = Bt*Cv/Cp;
            w2  = avV/(nMols*18.015*Bs);        // for Water only!!!
            w   = sqrt(w2);

	    // Properties using the interatomic energy
	    HV_HV = HintV - Hint*avV;
	    Cpint = (Hint2 - Hint*Hint)/(BOLTZMAN*Ti*Ti);
	    Cvint = Cpint - HV_HV*HV_HV/(avV2 - avV*avV)/(BOLTZMAN*Ti*Ti);
	    Bsint = Bt*Cvint/Cpint;
            w2int = avV/(nMols*18.015*Bsint);        // for Water only!!!
            wint  = sqrt(w2int);

	    // Now in MPa and J/molK units
	    Cp  *= 1000/nMols;
            Cv  *= 1000/nMols;
	    ap  *= 1E04;	// in  K-1 10-4
            Bt  *= 1000/(1.66053904042716);  // in GPa
            Bs  *= 1000/(1.66053904042716);  // in GPa
            w   *= 1000;    // in m/s

	    // Hint in MPa and J/molK units
	    Cpint *= 1000/nMols;
            Cvint *= 1000/nMols;
            Bsint *= 1000/(1.66053904042716);  // in GPa
            wint  *= 1000;    // in m/s
	}

	of<< luscount;
        of<< scientific << setprecision(8) << setw(17) << Uu;
        of<< scientific << setprecision(8) << setw(17) << Ti;
        of<< scientific << setprecision(8) << setw(17) << Pp;
        of<< scientific << setprecision(8) << setw(17) << Cvens;
        of<< scientific << setprecision(8) << setw(17) << 1.66053904042716*Pvir;
        of<< scientific << setprecision(8) << setw(17) << myAvLus.Volnm/luscount;
        of<< scientific << setprecision(8) << setw(17) << myAvLus.Volsq/luscount;
	of<< scientific << setprecision(8) << setw(17) << Cp;
        of<< scientific << setprecision(8) << setw(17) << Cv;
        of<< scientific << setprecision(8) << setw(17) << ap;
        of<< scientific << setprecision(8) << setw(17) << Bt;
        of<< scientific << setprecision(8) << setw(17) << Bs;
        of<< scientific << setprecision(8) << setw(17) << w;
	of<< scientific << setprecision(8) << setw(17) << Cpint;
        of<< scientific << setprecision(8) << setw(17) << Cvint;
        of<< scientific << setprecision(8) << setw(17) << Bsint;
        of<< scientific << setprecision(8) << setw(17) << wint;
	of<< endl;
    }
}

void Ensemble :: write_averaglus(ofstream& of, Int numSteps)
{

    of<< numSteps << "   " <<luscount;
    of<< scientific << setprecision(8) << setw(17) << myAvLus.bUF/luscount;
    of<< scientific << setprecision(8) << setw(17) << myAvLus.bUF2/luscount;
    of<< scientific << setprecision(8) << setw(17) << myAvLus.bUFt/luscount;
    of<< scientific << setprecision(8) << setw(17) << myAvLus.bUF2t/luscount;
    of<< scientific << setprecision(8) << setw(17) << myAvLus.Kinav/luscount;
    of<< scientific << setprecision(8) << setw(17) << myAvLus.invKin/luscount;
    of<< scientific << setprecision(8) << setw(17) << myAvLus.Uen/luscount;
    of<< scientific << setprecision(8) << setw(17) << myAvLus.Uen2/luscount;
    of<< scientific << setprecision(8) << setw(17) << myAvLus.Kin2/luscount;
    of<< scientific << setprecision(8) << setw(17) << myAvLus.Entot/luscount;
    of<< scientific << setprecision(8) << setw(17) << myAvLus.Entot2/luscount;
    of<< scientific << setprecision(8) << setw(17) << myAvLus.Kinavat/luscount;
    of<< scientific << setprecision(8) << setw(17) << myAvLus.invKinat/luscount;
    of<< scientific << setprecision(8) << setw(17) << myAvLus.Kin2at/luscount;
    of<< endl;
}

void Ensemble :: write_velbehav(ofstream& of, Int numSteps)
{

int ns, ccountt;
double vomx, vomy, vomz, masstott, mmden;
Vector3 veldiff, ccvelv, co1diff, co2diff;
double diffv, ccvel, difco1, difco2, stdco1, stdco2;
double difval1, difval2, dr12, dr32, cosTheta, anglt;
double stdan, stdangle;
//Atom* atom;

vomx = 0.0;
vomy = 0.0;
vomz = 0.0;
masstott = 0.0;
for (int i = 0; i < nAtoms; i++){
        vomx += atoms[i].velocity.x*atoms[i].mass;
        vomy += atoms[i].velocity.y*atoms[i].mass;
        vomz += atoms[i].velocity.z*atoms[i].mass;
        masstott += atoms[i].mass;
}
ccountt = 0.0;
diffv = 0.0;
ccvel = 0.0;
difco1 = 0.0;
difco2 = 0.0;
stdco1 = 0.0;
stdco2 = 0.0;

of<< numSteps << setw(14) << vomx/masstott <<  setw(14) << vomy/masstott <<  setw(14) << vomz/masstott;
for (int nm = 0; nm < nMols; nm++){
        ns = molecules[nm].numAtoms;
        if (!strcmp(molecules[nm].molName, "CO2")) {
                double MM = molecules[nm].mass;
                ccountt = ccountt + 1;
                Atom *atom0 = molecules[nm].myAtoms[0];
                Atom *atom1 = molecules[nm].myAtoms[1];
                Atom *atom2 = molecules[nm].myAtoms[2];

                Int id0 = atom0->atomID;
                Int id1 = atom1->atomID;
                Int id2 = atom2->atomID;

                veldiff = atom0->velocity - atom2->velocity;
                diffv  += veldiff.length();
                ccvelv  = atom1->velocity;
                ccvel  += ccvelv.length();

                co1diff = atom0->position - atom1->position;
                co2diff = atom2->position - atom1->position;
                difco1 += co1diff.length();
                difco2 += co2diff.length();
                mmden  += MM;

		dr12 = 1.0/co1diff.length();
		dr32 = 1.0/co2diff.length();
		cosTheta = (co1diff*co2diff)*dr12*dr32;
		anglt  += acos(cosTheta)*180.00/PI;

		//vomx += 0.5*(atom0->velocity.x)*MM + 0.5*(atom1->velocity.x)*MM;
		//vomy += 0.5*(atom0->velocity.y)*MM + 0.5*(atom1->velocity.y)*MM;
		//vomz += 0.5*(atom0->velocity.z)*MM + 0.5*(atom1->velocity.z)*MM;

       }
        if (!strcmp(molecules[nm].molName, "DAM")) {
                double MM = molecules[nm].mass;
                ccountt = ccountt + 1;
                Atom *atom0 = molecules[nm].myAtoms[0];
                Atom *atom1 = molecules[nm].myAtoms[1];

                Int id0 = atom0->atomID;
                Int id1 = atom1->atomID;

                veldiff = atom0->velocity - atom1->velocity;
                diffv  += veldiff.length();
                ccvelv  = atom1->velocity;
                ccvel  += ccvelv.length();

                co1diff = atom0->position - atom1->position;
                difco1 += co1diff.length();
                mmden  += MM;
       }
       if (!strcmp(molecules[nm].molName, "BBV")) {
                double MM = molecules[nm].mass;
                ccountt = ccountt + 1;
                Atom *atom0 = molecules[nm].myAtoms[0];
                Atom *atom2 = molecules[nm].myAtoms[2];
                Atom *atom1 = molecules[nm].myAtoms[4];

                Int id0 = atom0->atomID;
                Int id1 = atom1->atomID;
                Int id2 = atom2->atomID;

                veldiff = atom0->velocity - atom1->velocity;
                diffv  += veldiff.length();
                ccvelv  = atom2->velocity;
                ccvel  += ccvelv.length();

                co1diff = atom0->position - atom2->position;
                co2diff = atom2->position - atom1->position;
                difco1 += co1diff.length();
                difco2 += co2diff.length();
                mmden  += MM;

                dr12 = 1.0/co1diff.length();
                dr32 = 1.0/co2diff.length();
                cosTheta = (co1diff*co2diff)*dr12*dr32;
                anglt  += acos(cosTheta)*180.00/PI;

                //vomx += 0.5*(atom0->velocity.x)*MM + 0.5*(atom1->velocity.x)*MM;
                //vomy += 0.5*(atom0->velocity.y)*MM + 0.5*(atom1->velocity.y)*MM;
                //vomz += 0.5*(atom0->velocity.z)*MM + 0.5*(atom1->velocity.z)*MM;

        }
	if ( (!strcmp(molecules[nm].molName, "WAT")) ||  (!strcmp(molecules[nm].molName, "MCY"))) {
		double MM = molecules[nm].mass;
                ccountt = ccountt + 1;
                Atom *atom0 = molecules[nm].myAtoms[0];
                Atom *atom1 = molecules[nm].myAtoms[1];
                Atom *atom2 = molecules[nm].myAtoms[2];

                Int id0 = atom0->atomID;
                Int id1 = atom1->atomID;
                Int id2 = atom2->atomID;

                veldiff = atom0->velocity - atom2->velocity;
                diffv  += veldiff.length();
                ccvelv  = atom1->velocity;
                ccvel  += ccvelv.length();

                co1diff = atom0->position - atom1->position;
                co2diff = atom2->position - atom1->position;
                difco1 += co1diff.length();
                difco2 += co2diff.length();
                mmden  += MM;

		dr12 = 1.0/co1diff.length();
		dr32 = 1.0/co2diff.length();
		cosTheta = (co1diff*co2diff)*dr12*dr32;
		anglt  += acos(cosTheta)*180.00/PI;
	}

}

stdco1   = 0.0;
stdco2   = 0.0;
stdangle = 0.0;

for (int nm = 0; nm < nMols; nm++){
        ns = molecules[nm].numAtoms;
        if (!strcmp(molecules[nm].molName, "CO2") ||  !strcmp(molecules[nm].molName, "WAT") || !strcmp(molecules[nm].molName, "MCY")) {
                Atom *atom0 = molecules[nm].myAtoms[0];
                Atom *atom1 = molecules[nm].myAtoms[1];
                Atom *atom2 = molecules[nm].myAtoms[2];

                co1diff = atom0->position - atom1->position;
                co2diff = atom2->position - atom1->position;

                difval1 = co1diff.length() - difco1/ccountt;
                difval2 = co2diff.length() - difco2/ccountt;
                stdco1 += difval1*difval1;
                stdco2 += difval2*difval2;

		dr12 = 1.0/co1diff.length();
                dr32 = 1.0/co2diff.length();
                cosTheta = (co1diff*co2diff)*dr12*dr32;
                stdan = acos(cosTheta)*180.00/PI - anglt/ccountt;
		stdangle += stdan*stdan;

        }
        if (!strcmp(molecules[nm].molName, "DAM")) {
                Atom *atom0 = molecules[nm].myAtoms[0];
                Atom *atom1 = molecules[nm].myAtoms[1];

                co1diff = atom0->position - atom1->position;

                difval1 = co1diff.length() - difco1/ccountt;
                stdco1 += difval1*difval1;
        }
	if (!strcmp(molecules[nm].molName, "BBV")) {
                Atom *atom0 = molecules[nm].myAtoms[0];
                Atom *atom2 = molecules[nm].myAtoms[2];
                Atom *atom1 = molecules[nm].myAtoms[4];

                co1diff = atom0->position - atom2->position;
                co2diff = atom2->position - atom1->position;

                difval1 = co1diff.length() - difco1/ccountt;
                difval2 = co2diff.length() - difco2/ccountt;
                stdco1 += difval1*difval1;
                stdco2 += difval2*difval2;

                dr12 = 1.0/co1diff.length();
                dr32 = 1.0/co2diff.length();
                cosTheta = (co1diff*co2diff)*dr12*dr32;
                stdan = acos(cosTheta)*180.00/PI - anglt/ccountt;
                stdangle += stdan*stdan;

        }
}

   of << setw(10) << diffv/ccountt << setw(10) << ccvel/ccountt << setw(10) << mmden/(boxLx*boxLy*boxLz)/6.0221408e+23/1000/1e-27;
   of << setw(12) << difco1/ccountt << setw(12) << sqrt(fabs(stdco1)/(ccountt-1)) << setw(12) << difco2/ccountt << setw(12) << sqrt(fabs(stdco2)/(ccountt-1));
   of << setw(10) << anglt/ccountt << setw(12) << sqrt(fabs(stdangle)/(ccountt-1));
   of << endl;

}

// atom trajectory data usually for visiualisation purpose
void Ensemble::write_trajectory(ofstream &of, Int numSteps)
{
    of << numSteps << endl;
    for (Int i = 0; i < nAtoms; i++)
        of << atoms[i].position.x << '\t'<< atoms[i].position.y << '\t'<< atoms[i].position.z << endl;
}

void Ensemble :: write_bond(Molecule *mol, ofstream &of)
{
    double l1, l2, l3, angle;
    Vector3 r12, r32, r13;

    int ns = mol->numAtoms;
    if (ns <= 1)
        return;
    of << "mol[" << mol->molID << "]:  ";
    if (ns == 2)
    {
        r12 = mol->myAtoms[1]->position - mol->myAtoms[0]->position;
        of << l1 ;
    }
    else
    for(int site = 0; site < (ns-2); site++)
    {
        r12 = mol->myAtoms[site + 1]->position - mol->myAtoms[site]->position;
        l1 = r12.length2();
        r32 = mol->myAtoms[site + 2]->position - mol->myAtoms[site + 1]->position;
        l2 = r32.length2();
        r13 = mol->myAtoms[site + 2]->position - mol->myAtoms[site]->position;
        l3 = r13.length2();
        of << l1 << ' ' << l2 << ' ' << l3;
    }
    of << endl;
}

// store MSD/VACF sampling data to files for analysing after simulation
// currently the results are writen into default files and the data include
// necessary parameters about each species and the simulation
void Ensemble::write_sampling_data()
{
    ofstream *msd = NULL;
    ofstream *vacf = NULL;

    if (computeMSD)     msd = new ofstream(myConfig->get_mol_traject_file(), ios::out);
    if (computeVACF)    vacf = new ofstream(myConfig->get_mol_velocity_file(), ios::out);

    if(msd)
    {
//        *msd << nMols << '\t' << "// number of molecules" << endl;
//        *msd << nMolTypes << '\t' << "// number of molecular types" << endl;
//        *msd << nSamples  << '\t' << "// total number of samplings" << endl;
//        *msd << interval  << '\t' << "// interval number of timesteps for sampling" << endl;
//        *msd << deltaT    << '\t' << "// timeStep of simulation" << endl;
        // molecular info
//        for (Int i = 0; i < nMols; i++)
//            *msd << i << '\t' << molecules[i].molType << endl;
        // molecular positions of all samples
        for (Int i = 0; i < nSamples; i++)
        {
            *msd << nAtoms << endl;/////////////////////////////////////////////////////////////////////////////////////
            *msd << i << endl;
            for (Int nm = 0; nm < nMols; nm++)
                *msd<<"1"/*myParams->atomParams[i].atomType*/<<'\t'<<trajectory[i][nm].x<<'\t'<<trajectory[i][nm].y<<'\t'<<trajectory[i][nm].z<<endl;
        }
    }

    if(vacf)
    {
        *vacf << nMols << '\t' << "// number of molecules" << endl;
        *vacf << nMolTypes << '\t' << "// number of molecular types" << endl;
        *vacf << nSamples  << '\t' << "// total number of samplings" << endl;
        *vacf << interval  << '\t' << "// interval number of timesteps for sampling" << endl;
        *vacf << deltaT    << '\t' << "// timeStep of simulation" << endl;
        // molecular info
        for (Int i = 0; i < nMols; i++)
            *vacf << i << '\t' << molecules[i].molType << endl;
        // molecular velocities of all samples
        for (Int i = 0; i < nSamples; i++)
        {
            *vacf << i << endl;
            for (Int nm = 0; nm < nMols; nm++)
                *vacf<<velocity[i][nm].x<<'\t'<<velocity[i][nm].y<<'\t'<<velocity[i][nm].z<<endl;
        }
    }
}

void Ensemble::random_velocity(unsigned seed)
{
    Int i, type, n, ndf;
    Double mass, pxavg, pyavg, pzavg, tavg, scale;
    Vector3* momentum;

    DEBUGMSG("random velocity");

    if(myConfig->use_atom_thermo())
    {
        n = nAtoms;
        if(myConfig->is_constraint_on())
            ndf = 3*nAtoms - nBonds - nAngles - 3;
        else
            ndf = 3*nAtoms - 3;
    }
    else
    {
        n = nMols;
        ndf = 3*nMols - 3;
    }

    momentum = new Vector3[n];
    pxavg = 0.0;
    pyavg = 0.0;
    pzavg = 0.0;
    tavg = 0.0;

    srand(seed);
    for(i = 0; i < n; i++)
    {
        momentum[i].x = (Double)rand()/RAND_MAX;         // random 0 - 1
        momentum[i].x = momentum[i].x - 0.5;             // random -0.5 - 0.5
        pxavg += momentum[i].x;
        momentum[i].y = (Double)rand()/RAND_MAX;
        momentum[i].y = momentum[i].y - 0.5;
        pyavg += momentum[i].y;
        momentum[i].z = (Double)rand()/RAND_MAX;
        momentum[i].z = momentum[i].z - 0.5;
        pzavg += momentum[i].z;
    }
    pxavg = pxavg/n;
    pyavg = pyavg/n;
    pzavg = pzavg/n;


    // zero total momentum
    if(myConfig->use_atom_thermo())
    {
        for(i = 0; i < n; i++)
        {
            atoms[i].momentum.x = momentum[i].x - pxavg;
            atoms[i].momentum.y = momentum[i].y - pyavg;
            atoms[i].momentum.z = momentum[i].z - pzavg;

            // sum of momentum*velocity
            tavg += atoms[i].momentum.length2()/atoms[i].mass;
        }
        tavg = tavg/(ndf* BOLTZMAN);
        scale = sqrt(myConfig->get_temperature()/tavg);
        for(i = 0; i < n; i++)
        {
            atoms[i].momentum *= scale;
            atoms[i].velocity = atoms[i].momentum/atoms[i].mass;
        }
    }
    else
    {
        for(i = 0; i < n; i++)
        {
            molecules[i].momenta.x = momentum[i].x - pxavg;
            molecules[i].momenta.y = momentum[i].y - pyavg;
            molecules[i].momenta.z = momentum[i].z - pzavg;

            // sum of momentum*velocity
            tavg += molecules[i].momenta.length2()/molecules[i].mass;
        }
        // temperature  - assuming reduced units used
        // tavg = tavg/(3.0 * nMols - 3);

        // temperature - Using Kelvin temperature
        tavg = tavg/(ndf * BOLTZMAN);

        // scale velocities to initial temperature
        scale = sqrt(myConfig->get_temperature()/tavg);
        for(i = 0; i < n; i++)
        {
            molecules[i].momenta *= scale;
            Vector3 velocity = molecules[i].momenta/molecules[i].mass;
            for (int j = 0; j < molecules[i].numAtoms; j++)
            {
                molecules[i].myAtoms[j]->velocity = velocity;
                molecules[i].myAtoms[j]->momentum = molecules[i].myAtoms[j]->mass*velocity;
            }
        }
    }
    delete [] momentum;
}


/** init_velocity() - initialise translational velocity from Maxwell-Boltzmann Distribution
 ** The distribution is determined by temperature and (real) mass.
 **
 ** Reference:
 **     Allen & Tildesley, Computer Simulation of Liquids, 1989
 **/
void Ensemble::init_velocity()
{
    Int i, n, ndf;
    Double pxavg, pyavg, pzavg, tavg, scale;
    Vector3* momentum;

    DEBUGMSG("init velocity");

    if(myConfig->use_atom_thermo())
    {
        n = nAtoms;
        if(myConfig->is_constraint_on())
            ndf = 3*nAtoms - nBonds - nAngles - 3;
        else
            ndf = 3*nAtoms - 3;
    }
    else
    {
        n = nMols;
        ndf = 3*nMols - 3;
    }

    momentum = new Vector3[n];
    pxavg = 0.0;
    pyavg = 0.0;
    pzavg = 0.0;
    tavg = 0.0;

    for(i = 0; i < n; i++)
    {
        momentum[i].x = gauss();
        momentum[i].y = gauss();
        momentum[i].z = gauss();

        pxavg += momentum[i].x;
        pyavg += momentum[i].y;
        pzavg += momentum[i].z;
    }
    pxavg = pxavg/n;
    pyavg = pyavg/n;
    pzavg = pzavg/n;

    for(i = 0; i < n; i++)
    {
        momentum[i].x -= pxavg;
        momentum[i].y -= pyavg;
        momentum[i].z -= pzavg;

        if(myConfig->use_atom_thermo())
            tavg += momentum[i].length2()/atoms[i].mass;
        else
            tavg += momentum[i].length2()/molecules[i].mass;
    }
    tavg /= (ndf* BOLTZMAN);
    scale = sqrt(myConfig->get_temperature()/tavg);

    for (i = 0; i < n; i++)
    {
        if(myConfig->use_atom_thermo())
        {
            atoms[i].momentum = scale*momentum[i];
            atoms[i].velocity = atoms[i].momentum/atoms[i].mass;
        }
        else
        {
            molecules[i].momenta = scale*momentum[i];
            Vector3 velocity = molecules[i].momenta/molecules[i].mass;
            for (int j = 0; j < molecules[i].numAtoms; j++)
            {
                molecules[i].myAtoms[j]->velocity = velocity;
                molecules[i].myAtoms[j]->momentum = molecules[i].myAtoms[j]->mass*velocity;
            }
        }
    }
    delete [] momentum;
}

// generate uniform random number from standard normal Gaussian distribution
// with zero mean and unit variance
Double Ensemble::gauss()
{
    Double A1, A3, A5, A7, A9;
    Double sum, r, r2, rst;

    // constants
    A1 = 3.949846138;
    A3 = 0.252408784;
    A5 = 0.076542912;
    A7 = 0.008355968;
    A9 = 0.029899776;

    sum = 0.0;
    for (int i = 0; i < 12; i++)
        sum += drand();
    r = (sum - 6.0)/4.0;
    r2 = r*r;
    rst = (((A9*r2 + A7)*r2 + A5)*r2 + A3)*r2 + A1;

    return rst*r;
}

Double Ensemble::drand()
{
    Int l, c, m;
    static unsigned int seed = mySeed;

    l = 1029;
    c = 221591;
    m = 1048576;

    seed = (seed * l + c)%m;
    return (Double) (seed)/m;
}

//JC function member created by Jianhui LI in accordence with the
void Ensemble:: write_Potential(ofstream &of)
{
//    of << WPotential <<'\t' << CPotential <<endl;
}


// Jc: function member created by Jianhui LI to output pressure tensor
// Jc: for shear viscosity calcualtion
void Ensemble:: write_systemInfo(ofstream &of)
{
  of<<temperature<<endl;
  of<<volume<<endl;
  of<<nSamples<<endl;
  of<<interval<<endl;
  of<<deltaT<<endl;
}

// Jc: function member created by Jianhui LI to output pressure tensor
// Jc: for shear viscosity calcualtion

void Ensemble:: write_pressureTensor(ofstream &of, int counter1)
{

  of<<counter1<<'\t'<<endl;
  of<<scientific<<"  XX = "<<atomPressure[XX]<<"  YY = "<<atomPressure[YY]<<"  ZZ = "<<atomPressure[ZZ];////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  of<<scientific<<"  XY = "<<atomPressure[XY]<<"  YX = "<<atomPressure[YX]<<"  XZ = "<<atomPressure[XZ]<<"  ZX = "<<atomPressure[ZX];///////////////////////////////////////////////////////////////////////////////////
  of<<scientific<<"  YZ = "<<atomPressure[YZ]<<"  ZY = "<<atomPressure[ZY];/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  of<<scientific<<setw(14)<<(atomPressure[XY] + atomPressure[YX])*0.5;
  of<<scientific<<setw(14)<<(atomPressure[ZY] + atomPressure[YZ])*0.5;
  of<<scientific<<setw(14)<<(atomPressure[XZ] + atomPressure[ZX])*0.5;
  of<<scientific<<setw(14)<<(atomPressure[XZ]+atomPressure[ZX] + \
    atomPressure[ZY] + atomPressure[YZ]+atomPressure[XY] + atomPressure[YX])/6.0<<endl;

  of<<scientific<<"  molXX = "<<molPressure[XX]<<"  molYY = "<<molPressure[YY]<<"  molZZ = "<<molPressure[ZZ];////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  of<<scientific<<"  molXY = "<<molPressure[XY]<<"  molYX = "<<molPressure[YX]<<"  molXZ = "<<molPressure[XZ]<<"  molZX = "<<molPressure[ZX];/////////////////////////////////////////////////////////////////////////////////
  of<<scientific<<"  molYZ = "<<molPressure[YZ]<<"  molZY = "<<molPressure[ZY];///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  of<<scientific<<setw(14)<<(molPressure[XY]  + molPressure[YX])*0.5;
  of<<scientific<<setw(14)<<(molPressure[ZY]  + molPressure[YZ])*0.5;
  of<<scientific<<setw(14)<<(molPressure[XZ]  + molPressure[ZX])*0.5;
  of<<scientific<<setw(14)<<(molPressure[XZ]+molPressure[ZX] + \
    molPressure[ZY] + molPressure[YZ] + molPressure[XY] + molPressure[YX])/6.0<<endl;

  of<<scientific<<setw(14)<<atomStress.shearStress/counter;
  of<<scientific<<setw(14)<<molStress.shearStress/counter<<endl;

}

void Ensemble:: write_position(int timeStep)
{
    char* tab = "       ";
    ofo = new ofstream("position.out", ios::out);
    *ofo << "Coordinate at step: " << timeStep << endl;
    for (int i = 0; i < nAtoms; i++)
        *ofo << i<< tab <<atoms[i].position.x << tab << atoms[i].position.y << tab << atoms[i].position.z << endl;
    *ofo << endl;
}
