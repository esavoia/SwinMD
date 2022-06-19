/** Parameters.cpp --
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

#include "Parameters.h"
#include <iostream>
#include "Errors.h"

Parameters :: Parameters(const char* sysDataFile)
{
    #ifdef DEBUG
        DEBUGMSG("Creating Parameter object");
    #endif

//    rank = MPI::COMM_WORLD.Get_rank();  Updated to MPI-3 elton
//    size = MPI::COMM_WORLD.Get_size();
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    atomParams = NULL;
    bondParams = NULL;
    angleParams = NULL;
    dihedralParams = NULL;
    improperParams = NULL;
    molTypes = NULL;
    ljParamTable = NULL;

    atomArray = NULL;
    molArray = NULL;
    bondArray = NULL;
    angleArray = NULL;
    dihedralArray = NULL;
    improperArray = NULL;

    nAtomTypes = 0;
    nMolTypes = 0;
    nBondTypes = 0;
    nAngleTypes = 0;
    nDihedralTypes = 0;
    nImproperTypes = 0;
    ljTableSize = 0;

    nAtoms = 0;
    nMolecules = 0;
    nBonds = 0;
    nAngles = 0;
    nDihedrals = 0;
    nImpropers = 0;

    read_sys_data(sysDataFile);
    // done_reading_sysData();
}

Parameters :: ~Parameters()
{
    if (atomParams != NULL)     delete [] atomParams;
    if (bondParams != NULL)     delete [] bondParams;
    if (angleParams != NULL)    delete [] angleParams;
    if (dihedralParams != NULL) delete [] dihedralParams;
    if (improperParams != NULL) delete [] improperParams;
    if (molTypes != NULL)       delete [] molTypes;
    if (ljParamTable != NULL)   delete [] ljParamTable;

    if (atomArray != NULL)      delete [] atomArray;
    if (molArray != NULL)       delete [] molArray;
    if (bondArray != NULL)      delete [] bondArray;
    if (angleArray != NULL)     delete [] angleArray;
    if (dihedralArray != NULL)  delete [] dihedralArray;
    if (improperArray != NULL)  delete [] improperArray;
}

void Parameters::read_sys_data(const char* sysDataFile)
{
    int  i, ret, v1, v2, v3, v4, v5, v6;
    double d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14;
    char buf[256];
    char str1[32];
    char str2[32];

//    #ifdef DEBUG
//        DEBUGMSG("Reading system data");
//    #endif

    FILE *fptr;

    if(rank ==0)
    {
      if ((fptr = fopen(sysDataFile, "r")) == NULL)
          ERRORMSG("open sysDataFile error");

    // skip first two lines
      fgets(buf, 512, fptr);
      fgets(buf, 512, fptr);

    // next 6 lines should contain the number about molecular system sizes
    // and each line contains two numbers
      fgets(buf, 512, fptr);
      if((ret = sscanf(buf, "%d%d", &nMolecules, &nMolTypes)) != 2)
          ERRORMSG("read nMolecules & nMolTypes error ");

      fgets(buf, 512, fptr);
      if((ret = sscanf(buf, "%d%d", &nAtoms, &nAtomTypes)) != 2)
          ERRORMSG("read nAtom & nAtomTypes error \n");

      fgets(buf, 512, fptr);
      if((ret = sscanf(buf, "%d%d", &nBonds, &nBondTypes)) != 2)
          ERRORMSG("read nBonds & nBondTypes error \n");

      fgets(buf, 512, fptr);
      if((ret = sscanf(buf, "%d%d", &nAngles, &nAngleTypes)) != 2)
          ERRORMSG("read nAngles & nAngleTypes error \n");

      fgets(buf, 512, fptr);
      if((ret = sscanf(buf, "%d%d", &nDihedrals, &nDihedralTypes)) != 2)
          ERRORMSG("read nDihedrals & nDihedralTypes error");

      fgets(buf, 512, fptr);
      if((ret = sscanf(buf, "%d%d", &nImpropers, &nImproperTypes)) != 2)
          ERRORMSG("read nImpropers & nImproperTypes error");
    }   // Jc: end of read sysData in processor 0

// Jc: Broadcast to other processors

//    MPI::COMM_WORLD.Bcast(&nMolecules,1,MPI::DOUBLE,0);
    MPI_Bcast(&nMolecules,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
//    MPI::COMM_WORLD.Bcast(&nMolTypes,1,MPI::DOUBLE,0);
    MPI_Bcast(&nMolTypes,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
//    MPI::COMM_WORLD.Bcast(&nAtoms,1,MPI::DOUBLE,0);
    MPI_Bcast(&nAtoms,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
//    MPI::COMM_WORLD.Bcast(&nAtomTypes,1,MPI::DOUBLE,0);
    MPI_Bcast(&nAtomTypes,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
//    MPI::COMM_WORLD.Bcast(&nBonds,1,MPI::DOUBLE,0);
    MPI_Bcast(&nBonds,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
//    MPI::COMM_WORLD.Bcast(&nBondTypes,1,MPI::DOUBLE,0);
    MPI_Bcast(&nBondTypes,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
//    MPI::COMM_WORLD.Bcast(&nAngles,1,MPI::DOUBLE,0);
    MPI_Bcast(&nAngles,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
//    MPI::COMM_WORLD.Bcast(&nAngleTypes,1,MPI::DOUBLE,0);
    MPI_Bcast(&nAngleTypes,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
//    MPI::COMM_WORLD.Bcast(&nDihedrals,1,MPI::DOUBLE,0);
    MPI_Bcast(&nDihedrals,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
//    MPI::COMM_WORLD.Bcast(&nDihedralTypes,1,MPI::DOUBLE,0);
    MPI_Bcast(&nDihedralTypes,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
//    MPI::COMM_WORLD.Bcast(&nImpropers,1,MPI::DOUBLE,0);
    MPI_Bcast(&nImpropers,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
//    MPI::COMM_WORLD.Bcast(&nImproperTypes,1,MPI::DOUBLE,0);
    MPI_Bcast(&nImproperTypes,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

//    cout<< "the received data are.................... "<<rank<<" "<<nAngleTypes<<endl;


    // allocate memory to data structures
      if(nAtomTypes)
      {
          atomParams = new AtomParam[nAtomTypes];
          if (atomParams == NULL)
              ERRORMSG("memory allocation of atomParams error");
      }
      if(nMolTypes)
      {
          molTypes = new MolType[nMolTypes];
          if (molTypes == NULL)
              ERRORMSG("memory allocation of nMolTypes error");
      }
      if(nBondTypes)
      {
          bondParams = new BondParam[nBondTypes];
          if (bondParams == NULL)
              ERRORMSG("memory allocation of bondParams error");
      }

      if(nAngleTypes)
      {
          angleParams = new AngleParam[nAngleTypes];
          if (angleParams == NULL)
              ERRORMSG("memory allocation of angleParams error");
      }
      if(nDihedralTypes)
      {
          dihedralParams = new TorsionParam[nDihedralTypes];
          if (dihedralParams == NULL)
              ERRORMSG("memory allocation of dihedralParams error");
      }
      if(nImproperTypes)
      {
          improperParams = new TorsionParam[nImproperTypes];
          if (improperParams == NULL)
              ERRORMSG("memory allocation of improperParams error");
      }
      if(nAtoms)
      {
          atomArray = new AtomStruct[nAtoms];
          if (atomArray == NULL)
              ERRORMSG("memory allocation of atomArray error");
      }
      if(nMolecules)
      {
          molArray = new MolStruct[nMolecules];
          if (molArray == NULL)
              ERRORMSG("memory allocation of molArray error");
      }
      if(nBonds)
      {
          bondArray = new Bond[nBonds];
          if (bondArray == NULL)
              ERRORMSG("memory allocation of bondArray error");
      }
      if(nAngles)
      {
          angleArray = new Angle[nAngles];
          if (angleArray == NULL)
              ERRORMSG("memory allocation of angleArray error");
      }
      if(nDihedrals)
      {
          dihedralArray = new Dihedral[nDihedrals];
          if (dihedralArray == NULL)
              ERRORMSG("memory allocation of dihedralArray error");
      }
      if(nImpropers)
      {
          improperArray = new Improper[nImpropers];
          if (improperArray == NULL)
              ERRORMSG("memory allocation of improperArray error");
      }

      ljTableSize = nAtomTypes;
      if(ljTableSize)
      {
          ljParamTable = new LJPairParam[ljTableSize*ljTableSize];
          if (ljParamTable  == NULL)
              ERRORMSG( "memory allocation of ljParamTable error");
      }


/*Jc:    #ifdef DEBUG
        cerr << "DBG: nMols: "<< nMolecules <<", nMolType: "<< nMolTypes<< endl;
        cerr << "DBG: nAtom: " << nAtoms << ", nAtomType: " << nAtomTypes<< endl;
        cerr << "DBG: nBond: " << nBonds << ", nBondType: " << nBondTypes << endl;
        cerr << "DBG: nAngle: " << nAngles << ", nAngleType: " << nAngleTypes << endl;
        cerr << "DBG: nDihedral: " << nDihedrals << ", nDihedralType: " << nDihedralTypes << endl;
        cerr << "DBG: nImproper: " << nImpropers << ", nImproperType: " << nImproperTypes << endl;
    #endif
*/
    // read atom parameters
    if(rank == 0)
    {
      if (!check_heading("!AtomParams", fptr))
         ERRORMSG("head not match with AtomParams");
      if(nAtomTypes)
      {
       #ifdef DEBUG
          DEBUGMSG("read atom parameters");
       #endif
         for (i = 0; i < nAtomTypes; i++)
         {
            fgets(buf, 512, fptr);
            for(int desRank = 1;desRank<size;desRank++)
            {
//              MPI::COMM_WORLD.Send(buf,256,MPI_CHAR,desRank,0);
              MPI_Send(buf,256,MPI_CHAR,desRank,0,MPI_COMM_WORLD);
            }
            sscanf(buf, "%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%s%s", &v1, &d1, &d2, &d3, &d4, &d5, &d6, &d7, &d8, &d9, &d10, &d11, &d12, &d13, &d14, str1, str2);
            atomParams[i].atomType = v1;
            atomParams[i].mass = d1;
            atomParams[i].charge = d2;
            atomParams[i].eps = d3;
            atomParams[i].sigma = d4;
            atomParams[i].dipx = d5;
            atomParams[i].dipy = d6;
            atomParams[i].dipz = d7;
            atomParams[i].quadxx = d8;
            atomParams[i].quadxy = d9;
            atomParams[i].quadxz = d10;
            atomParams[i].quadyy = d11;
            atomParams[i].quadyz = d12;
            atomParams[i].quadzz = d13;
            atomParams[i].polar  = d14;
            strcpy(atomParams[i].typeName, str1);
            strcpy(atomParams[i].atomName, str2);
         }
      }
    }
    else
    {
      for (i = 0; i < nAtomTypes; i++)
      {
//        MPI::COMM_WORLD.Recv(buf,256,MPI::CHAR,0,0);
        MPI_Recv(buf,256,MPI_CHAR,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);   // not sure about Ignore
        sscanf(buf, "%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%s%s", &v1, &d1, &d2, &d3, &d4, &d5, &d6, &d7, &d8, &d9, &d10, &d11, &d12, &d13, &d14, str1, str2);
        atomParams[i].atomType = v1;
        atomParams[i].mass = d1;
        atomParams[i].charge = d2;
        atomParams[i].eps = d3;
        atomParams[i].sigma = d4;
        atomParams[i].dipx = d5;
        atomParams[i].dipy = d6;
        atomParams[i].dipz = d7;
        atomParams[i].quadxx = d8;
        atomParams[i].quadxy = d9;
        atomParams[i].quadxz = d10;
        atomParams[i].quadyy = d11;
        atomParams[i].quadyz = d12;
        atomParams[i].quadzz = d13;
        atomParams[i].polar  = d14;
        strcpy(atomParams[i].typeName, str1);
        strcpy(atomParams[i].atomName, str2);
      }
    }

    // read bond parameters
    if(rank == 0)
    {
      if (!check_heading("!BondParams", fptr))
          ERRORMSG("head not match with BondParams");
      if(nBondTypes)
      {
        #ifdef DEBUG
           DEBUGMSG("read bond parameters");
        #endif
        for (i = 0; i < nBondTypes; i++)
        {
          fgets(buf, 512, fptr);
          for(int desRank = 1;desRank<size;desRank++)
          {
//            MPI::COMM_WORLD.Send(buf,256,MPI::CHAR,desRank,0);
            MPI_Send(buf,256,MPI_CHAR,desRank,0,MPI_COMM_WORLD);
          }
          ret = sscanf(buf, "%d%d%d%lf%lf", &v1, &v2, &v3, &d1, &d2);
//         cout <<"sent data nBOnde  "<<rank<<" "<<v1<<" "<<v3 <<" "<<" "<<d1<<" "<<d2<<endl;

          bondParams[i].bondType = v1;
          bondParams[i].atomType1 = v2;
          bondParams[i].atomType2 = v3;
          bondParams[i].r0 = d1;
          bondParams[i].k = d2;
        }
      }
    }
    else
    {
      for (i = 0; i < nBondTypes; i++)
      {

//        MPI::COMM_WORLD.Recv(buf,256,MPI::CHAR,0,0);
        MPI_Recv(buf,256,MPI_CHAR,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        ret = sscanf(buf, "%d%d%d%lf%lf", &v1, &v2, &v3, &d1, &d2);
//        cout <<"received data nBOnde  "<<rank<<" "<<v1<<" "<<v3 <<" "<<" "<<d1<<" "<<d2<<endl;

        bondParams[i].bondType = v1;
        bondParams[i].atomType1 = v2;
        bondParams[i].atomType2 = v3;
        bondParams[i].r0 = d1;
        bondParams[i].k = d2;
      }
    }

    // read angle parameters

    if(rank ==0)
    {
      if(!check_heading("!AngleParams", fptr))
         ERRORMSG("head not match with AngleParams");
      if(nAngleTypes)
      {
        #ifdef DEBUG
          DEBUGMSG("read angle parameters");
        #endif

        for (i = 0; i < nAngleTypes; i++)
        {
          fgets(buf, 512, fptr);
//          cout <<"sent data nAangle at "<<rank<<" are "<<buf<<endl;

          for(int desRank = 1;desRank<size;desRank++)
          {
//            MPI::COMM_WORLD.Send(buf,256,MPI::CHAR,desRank,0);
            MPI_Send(buf,256,MPI_CHAR,desRank,0,MPI_COMM_WORLD);
          }
          ret= sscanf(buf, "%d%d%d%d%lf%lf", &v1, &v2, &v3, &v4, &d1,&d2);
//          cout <<"sent data nAangle at "<<rank<<" "<<v1<<" "<<v2<<" "<<v3 <<" "<<v4<<" "<<d1<<" "<<d2<<endl;

          angleParams[i].angleType = v1;
          angleParams[i].atomType1 = v2;
          angleParams[i].atomType2 = v3;
          angleParams[i].atomType3 = v4;
          angleParams[i].theta0 = d1*2*PI/360.0;           // assum the unit of d1 is 'degree'
          angleParams[i].k = d2;
        }
      }
    }
    else
    {
      for (i = 0; i < nAngleTypes; i++)
      {
//        MPI::COMM_WORLD.Recv(buf,256,MPI::CHAR,0,0);
        MPI_Recv(buf,256,MPI_CHAR,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);   // Not sure about Ignore
        ret = sscanf(buf, "%d%d%d%d%lf%lf", &v1, &v2, &v3, &v4, &d1,&d2);
//        cout <<"received data nAngle at "<<rank<<" "<<v1<<" "<<v2 <<" "<<v3<<" "<<v4<<" "<<d1<<" "<<d2<<endl;

        angleParams[i].angleType = v1;
        angleParams[i].atomType1 = v2;
        angleParams[i].atomType2 = v3;
        angleParams[i].atomType3 = v4;
        angleParams[i].theta0 = d1*2*PI/360.0;           // assum the unit of d1 is 'degree'
        angleParams[i].k = d2;
      }
    }

    // read dihedral parameters
    if(rank ==0)
    {
      if(!check_heading("!DihedralParams", fptr))
         ERRORMSG("head not match with DihedralParams");

      if(nDihedralTypes)
      {
        #ifdef DEBUG
          DEBUGMSG("read dihedral parameters");
        #endif

        for (i = 0; i < nDihedralTypes; i++)
        {
          fgets(buf, 512, fptr);
          for(int desRank = 1;desRank<size;desRank++)
          {
//            MPI::COMM_WORLD.Send(buf,256,MPI::CHAR,desRank,0);
            MPI_Send(buf,256,MPI_CHAR,desRank,0,MPI_COMM_WORLD);
          }

          sscanf(buf, "%d%d%d%d%d%lf%lf%lf", &v1,&v2,&v3,&v4,&v5,&d1,&d2,&d3);
          dihedralParams[i].torsionType = v1;
          dihedralParams[i].atomType1 = v2;
          dihedralParams[i].atomType2 = v3;
          dihedralParams[i].atomType3 = v4;
          dihedralParams[i].atomType4 = v5;
          dihedralParams[i].k = d1;
          dihedralParams[i].d = (int) d2;
          dihedralParams[i].n = (int) d3;
        }
      }
    }
    else
    {
      for (i = 0; i < nDihedralTypes; i++)
      {
//        MPI::COMM_WORLD.Recv(buf,256,MPI::CHAR,0,0);
        MPI_Recv(buf,256,MPI_CHAR,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);   // Not sure about Ignore, elton

        sscanf(buf, "%d%d%d%d%d%lf%lf%lf", &v1,&v2,&v3,&v4,&v5,&d1,&d2,&d3);
//       cout <<"received data nDihedral  "<<rank<<" "<<v1<<" "<<d1 <<" "<<d2<<endl;

        dihedralParams[i].torsionType = v1;
        dihedralParams[i].atomType1 = v2;
        dihedralParams[i].atomType2 = v3;
        dihedralParams[i].atomType3 = v4;
        dihedralParams[i].atomType4 = v5;
        dihedralParams[i].k = d1;
        dihedralParams[i].d = (int) d2;
        dihedralParams[i].n = (int) d3;
      }
    }

    // read improper parameters
    if(rank == 0)
    {
      if(!check_heading("!ImproperParams", fptr))
        ERRORMSG("head not match with ImproperParams");
      if(nImproperTypes)
      {
        #ifdef DEBUG
          DEBUGMSG("read improper parameters");
        #endif

        for (i = 0; i < nImproperTypes; i++)
        {
          fgets(buf, 512, fptr);
          for(int desRank = 1;desRank<size;desRank++)
          {
//            MPI::COMM_WORLD.Send(buf,256,MPI::CHAR,desRank,0);
            MPI_Send(buf,256,MPI_CHAR,desRank,0,MPI_COMM_WORLD);
          }
          sscanf(buf, "%d%d%d%d%d%lf%lf%lf", &v1,&v2,&v3,&v4,&v5,&d1,&d2,&d3);
          improperParams[i].torsionType = v1;
          improperParams[i].atomType1 = v2;
          improperParams[i].atomType2 = v3;
          improperParams[i].atomType3 = v4;
          improperParams[i].atomType4 = v5;
          improperParams[i].k = d1;
          improperParams[i].d = (int) d2;
          improperParams[i].n = (int) d3;
        }
      }
    }
    else
    {
      for (i = 0; i < nImproperTypes; i++)
      {
//        MPI::COMM_WORLD.Recv(buf,256,MPI::CHAR,0,0);
        MPI_Recv(buf,256,MPI_CHAR,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);   // Not sure about Ignore, elton

        sscanf(buf, "%d%d%d%d%d%lf%lf%lf", &v1,&v2,&v3,&v4,&v5,&d1,&d2,&d3);
//        cout <<"received data nImpro  "<<rank<<" "<<v1<<" "<<d1 <<" "<<d3<<endl;

        improperParams[i].torsionType = v1;
        improperParams[i].atomType1 = v2;
        improperParams[i].atomType2 = v3;
        improperParams[i].atomType3 = v4;
        improperParams[i].atomType4 = v5;
        improperParams[i].k = d1;
        improperParams[i].d = (int) d2;
        improperParams[i].n = (int) d3;
      }
    }

    // read molecule types
    if(rank ==0)
    {
      if(!check_heading("!MolTypes", fptr))
        ERRORMSG("head not match with MolTypes");

      if(nMolTypes)
      {
        #ifdef DEBUG
          DEBUGMSG("read molecular types");
        #endif

        for (i = 0; i < nMolTypes; i++)
        {
          fgets(buf, 512, fptr);
          for(int desRank = 1;desRank<size;desRank++)
          {
//            MPI::COMM_WORLD.Send(buf,256,MPI::CHAR,desRank,0);
            MPI_Send(buf,256,MPI_CHAR,desRank,0,MPI_COMM_WORLD);
          }

          sscanf(buf, "%d%d%s", &v1, &v2, str1);
//          cout <<"sent data molType  "<<rank<<" "<<v1<<" "<<v2 <<" "<<str1<<endl;

          molTypes[i].molType = v1;
          molTypes[i].numOfAtoms = v2;
          strcpy(molTypes[i].name, str1);
        }
      }
    }
    else
    {
      for (i = 0; i < nMolTypes; i++)
      {
//        MPI::COMM_WORLD.Recv(buf,256,MPI::CHAR,0,0);
        MPI_Recv(buf,256,MPI_CHAR,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        sscanf(buf, "%d%d%s", &v1, &v2, str1);
//        cout <<"received data molType  "<<rank<<" "<<v1<<" "<<v2 <<" "<<str1<<endl;

        molTypes[i].molType = v1;
        molTypes[i].numOfAtoms = v2;
        strcpy(molTypes[i].name, str1);
      }
    }

    // read topology
    // read atoms and molecules
    if(rank == 0)
    {
      if(!check_heading("!Atoms", fptr))
        ERRORMSG("head not match with Atoms");
      if (nAtoms)
      {
        #ifdef DEBUG
          DEBUGMSG("read atoms");
        #endif

        for (i = 0; i < nAtoms; i++)
        {
          fgets(buf, 512, fptr);

          for(int desRank = 1;desRank<size;desRank++)
          {
//            MPI::COMM_WORLD.Send(buf,256,MPI::CHAR,desRank,0);
            MPI_Send(buf,256,MPI_CHAR,desRank,0,MPI_COMM_WORLD);
          }

          sscanf(buf, "%d%d%d%d", &v1, &v2, &v3, &v4);
//          if(i==1499)
//          cout <<"sent out data are  "<<rank<<" "<<buf<<endl;

          atomArray[v1].atomID = v1;
          atomArray[v1].atomType = v2;
          atomArray[v1].molID = v3;
          molArray[v3].molID = v3;
          molArray[v3].molType = v4;
        }
      }
    }
    else
    {
      for (i = 0; i < nAtoms; i++)
      {
//        MPI::COMM_WORLD.Recv(buf,256,MPI::CHAR,0,0);
        MPI_Recv(buf,256,MPI_CHAR,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        sscanf(buf, "%d%d%d%d", &v1, &v2, &v3, &v4);
//        if(i==1499)
//        cout <<"receive out data are  "<<rank<<" "<<buf<<endl;
        atomArray[v1].atomID = v1;
        atomArray[v1].atomType = v2;
        atomArray[v1].molID = v3;
        molArray[v3].molID = v3;
        molArray[v3].molType = v4;
      }
    }
//    MPI::COMM_WORLD.Barrier();
    // read bonds
    if(rank ==0)
      {
      if(!check_heading("!Bonds", fptr))
        ERRORMSG("head not match with Bonds");
      if(nBonds)
      {
        #ifdef DEBUG
          DEBUGMSG("read bonds");
        #endif

        for (i = 0; i < nBonds; i++)
        {
          fgets(buf, 512, fptr);
          for(int desRank = 1;desRank<size;desRank++)
          {
//            MPI::COMM_WORLD.Send(buf,256,MPI::CHAR,desRank,0);
            MPI_Send(buf,256,MPI_CHAR,desRank,0,MPI_COMM_WORLD);
          }

          sscanf(buf, "%d%d%d%d", &v1, &v2, &v3, &v4);
//          if(i==nBonds-1)
//          cout <<"sent out data are  "<<rank<<" "<<v1<<" "<<v2<<" "<<v3<<" "<<v4<<endl;

          bondArray[v1].atom1 = v2;
          bondArray[v1].atom2 = v3;
          bondArray[v1].bondType = v4;
        }
      }
    }
    else
    {
        for (i = 0; i < nBonds; i++)
        {
//          MPI::COMM_WORLD.Recv(buf,256,MPI::CHAR,0,0);
          MPI_Recv(buf,256,MPI_CHAR,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
          sscanf(buf, "%d%d%d%d", &v1, &v2, &v3, &v4);
//          if(i==nBonds-1)
//          cout <<"received data are  "<<rank<<" "<<v1<<" "<<v2<<" "<<v3<<" "<<v4<<endl;

          bondArray[v1].atom1 = v2;
          bondArray[v1].atom2 = v3;
          bondArray[v1].bondType = v4;
        }
    }

    // read angles
    if(rank == 0)
    {
      if(!check_heading("!Angles", fptr))
        ERRORMSG("head not match with Angles");
      if(nAngles)
      {
        #ifdef DEBUG
          DEBUGMSG("read angles");
        #endif

        for (i = 0; i < nAngles; i++)
        {
          fgets(buf, 512, fptr);
          for(int desRank = 1;desRank<size;desRank++)
          {
//            MPI::COMM_WORLD.Send(buf,256,MPI::CHAR,desRank,0);
            MPI_Send(buf,256,MPI_CHAR,desRank,0,MPI_COMM_WORLD);
          }
          sscanf(buf, "%d%d%d%d%d", &v1, &v2, &v3, &v4, &v5);
//          if(i==nAngles-1)
//          cout <<"sent out data are  "<<rank<<" "<<v1<<" "<<v2<<" "<<v3<<" "<<v4<<"  "<<v5<<endl;

          angleArray[v1].atom1 = v2;
          angleArray[v1].atom2 = v3;
          angleArray[v1].atom3 = v4;
          angleArray[v1].angleType = v5;
        }
      }
    }
    else
    {
        for (i = 0; i < nAngles; i++)
        {
//          MPI::COMM_WORLD.Recv(buf,256,MPI::CHAR,0,0);
          MPI_Recv(buf,256,MPI_CHAR,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
          sscanf(buf, "%d%d%d%d%d", &v1, &v2, &v3, &v4, &v5);
//          if(i==nAngles-1)
//          cout <<"received data are  "<<rank<<" "<<v1<<" "<<v2<<" "<<v3<<" "<<v4<<"  "<<v5<<endl;

          angleArray[v1].atom1 = v2;
          angleArray[v1].atom2 = v3;
          angleArray[v1].atom3 = v4;
          angleArray[v1].angleType = v5;
        }
    }

//    cout<<" the system finished the reading process 999999........."<<rank<<endl;

    // read dihedrals
    if(rank ==0)
    {
      if(!check_heading("!Dihedrals", fptr))
        ERRORMSG("head not match with Dihedrals");
      if(nDihedrals)
      {
        #ifdef DEBUG
          DEBUGMSG("read dihedrals");
        #endif

        for (i = 0; i < nDihedrals; i++)
        {
          fgets(buf, 512, fptr);
          for(int desRank = 1;desRank<size;desRank++)
          {
//            MPI::COMM_WORLD.Send(buf,256,MPI::CHAR,desRank,0);
            MPI_Send(buf,256,MPI_CHAR,desRank,0,MPI_COMM_WORLD);
          }
          sscanf(buf, "%d%d%d%d%d%d", &v1, &v2, &v3, &v4, &v5, &v6);
          dihedralArray[v1].atom1 = v2;
          dihedralArray[v1].atom2 = v3;
          dihedralArray[v1].atom3 = v4;
          dihedralArray[v1].atom4 = v5;
          dihedralArray[v1].dihedralType = v6;
        }
      }
    }
    else
    {
        for (i = 0; i < nDihedrals; i++)
        {
//          MPI::COMM_WORLD.Recv(buf,256,MPI::CHAR,0,0);
          MPI_Recv(buf,256,MPI_CHAR,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
          sscanf(buf, "%d%d%d%d%d%d", &v1, &v2, &v3, &v4, &v5, &v6);
          dihedralArray[v1].atom1 = v2;
          dihedralArray[v1].atom2 = v3;
          dihedralArray[v1].atom3 = v4;
          dihedralArray[v1].atom4 = v5;
          dihedralArray[v1].dihedralType = v6;
        }
    }
//    cout<<" the system finished the reading process 888888......... "<<rank<<endl;

    // read impropers

    if(rank ==0)
    {
      if(!check_heading("!Impropers", fptr))
        ERRORMSG("head not match with Impropers");
      if(nImpropers)
      {
        #ifdef DEBUG
          DEBUGMSG("read impropers");
        #endif

        for (i = 0; i < nImpropers; i++)
        {
          fgets(buf, 512, fptr);
          for(int desRank = 1;desRank<size;desRank++)
          {
//            MPI::COMM_WORLD.Send(buf,256,MPI::CHAR,desRank,0);
            MPI_Send(buf,256,MPI_CHAR,desRank,0,MPI_COMM_WORLD);
          }
          sscanf(buf, "%d%d%d%d%d%d", &v1, &v2, &v3, &v4, &v5, &v6);
          improperArray[v1].atom1 = v2;
          improperArray[v1].atom2 = v3;
          improperArray[v1].atom3 = v4;
          improperArray[v1].atom4 = v5;
          improperArray[v1].improperType = v6;
        }
      }
    }
    else
    {
        for (i = 0; i < nImpropers; i++)
        {
//          MPI::COMM_WORLD.Recv(buf,256,MPI::CHAR,0,0);
          MPI_Recv(buf,256,MPI_CHAR,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
          sscanf(buf, "%d%d%d%d%d%d", &v1, &v2, &v3, &v4, &v5, &v6);
          improperArray[v1].atom1 = v2;
          improperArray[v1].atom2 = v3;
          improperArray[v1].atom3 = v4;
          improperArray[v1].atom4 = v5;
          improperArray[v1].improperType = v6;
        }
    }
//    cout<<" the system finished the reading process 7777......... "<<rank<<endl;
    if(rank ==0 )
      fclose(fptr);

//    cout<<" the system finished the reading process 6666........ ."<<rank<<endl;

//    MPI::COMM_WORLD.Barrier();
    build_lj_param_table();

//    #ifdef DEBUG
//        DEBUGMSG("Done reading system data");
//    #endif

}   // end read_sys_data()

// check_heading checks section head. if it is right return 1 else return 0
// error will be handled by the caller.
int Parameters::check_heading(const char* heading, FILE* fptr)
{
    char buf[256];
    char str[32];

    // skip a blank line
    fgets(buf, 512, fptr);
    // check section head
    fgets(buf, 512, fptr);
    sscanf(buf, "%s", str);
    if(strcmp(str, heading))     // may need case insensitive comparision
       return 0;
    // skip following notation line
    fgets(buf, 512, fptr);
    return 1;
}

void Parameters::build_lj_param_table()
{
    Int i, j;

    #ifdef DEBUG
        DEBUGMSG("build LJ table ");
    #endif

    if((ljTableSize == 0) || (ljParamTable == NULL))
        ERRORMSG("build_lj_param_table error");

    // lj_table is a triangle table. But duplicat values are set for the table.
    // This simplifies and speeds up the data retrieve in simulation
    for (i = 0; i < ljTableSize; i++)
    {
       for(j = 0; j < ljTableSize; j++)
       {
          set_lj_param(ljParamTable, atomParams[i], atomParams[j]);
       }
    }
//    cout<<" the program is coming here in the calcualation   .....................  "<<rank<<endl;
}

void Parameters::set_lj_param(LJPairParam *ljTable, const AtomParam &ap1,const AtomParam &ap2)
{
    Double sigma, eps;
    Int index;

    sigma = (ap1.sigma + ap2.sigma)/2;
    eps = sqrt(ap1.eps * ap2.eps);
    index = ap1.atomType*ljTableSize + ap2.atomType;

    ljTable[index].atomType1 = ap1.atomType;
    ljTable[index].atomType2 = ap2.atomType;
    ljTable[index].sigma = sigma*sigma;
    ljTable[index].eps = 4.0*eps;
}
