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
#include <fstream>
#include "Errors.h"

Parameters :: Parameters(const char* sysDataFile){

	DEBUGMSG("Creating Parameter object");

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

Parameters :: ~Parameters(){
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

#define SKIP_REST ifsys->ignore(std::numeric_limits<std::streamsize>::max(), '\n');


void Parameters::read_sys_data(const char* sysDataFile){
	DEBUGMSG("Reading system data '" << sysDataFile << "'");
	ifstream *ifsys = new ifstream(sysDataFile,ios::in);
	if (!ifsys || !ifsys->is_open()){
		ERRMSG("Can't open " << sysDataFile)
	}
	// skip first two lines
	SKIP_REST
	SKIP_REST
	// next 6 lines should contain the number about molecular system sizes
	// and each line contains two numbers
	*ifsys >> nMolecules >> nMolTypes;
	SKIP_REST
	*ifsys >> nAtoms >> nAtomTypes;
	SKIP_REST
	*ifsys >> nBonds >> nBondTypes;
	SKIP_REST
	*ifsys >> nAngles >> nAngleTypes;
	SKIP_REST
	*ifsys >> nDihedrals >> nDihedralTypes;
	SKIP_REST
	*ifsys >> nImpropers >> nImproperTypes;
	SKIP_REST
	// Allocate memory for data structures and populate them
	
	if(nAtomTypes){
		atomParams = new AtomParam[nAtomTypes];
		if (atomParams == NULL)
			ERRORMSG("memory allocation of atomParams error");
		
		// read atom parameters
		DEBUGMSG("read atom parameters");
		// check section header
		if (!check_heading("!AtomParams", ifsys))
			ERRORMSG("head not match with AtomParams");
		for (int i = 0; i < nAtomTypes; i++){
			AtomParam ap = atomParams[i];
			// oneliner
			// *ifsys >> ap.atomType >> ap.mass >> ap.charge >> ap.eps >> ap.sigma >> ap.dipx >> ap.dipy >> ap.dipz >> ap.quadxx >> ap.quadxy >> ap.quadxz >> ap.quadyy >> ap.quadyz >> ap.quadzz >> ap.polar  >> ap.typeName >> ap.atomName;
			// SKIP_REST
			*ifsys >> ap.atomType;
			*ifsys >> ap.mass;
			*ifsys >> ap.charge;
			*ifsys >> ap.eps;
			*ifsys >> ap.sigma;
			*ifsys >> ap.dipx;
			*ifsys >> ap.dipy;
			*ifsys >> ap.dipz;
			*ifsys >> ap.quadxx;
			*ifsys >> ap.quadxy;
			*ifsys >> ap.quadxz;
			*ifsys >> ap.quadyy;
			*ifsys >> ap.quadyz;
			*ifsys >> ap.quadzz;
			*ifsys >> ap.polar ;
			*ifsys >> ap.typeName;
			*ifsys >> ap.atomName;
			// DEBUGMSG(ap.atomType);
			// DEBUGMSG(ap.mass);
			SKIP_REST
		}
		
	} else {
		// check section header
		DEBUGMSG("check atom parameters");
		check_heading("!AtomParams", ifsys);
	}
	
	if(nBondTypes){
		bondParams = new BondParam[nBondTypes];
		if (bondParams == NULL)
			ERRORMSG("memory allocation of bondParams error");
		
		// read bond parameters
		DEBUGMSG("read bond parameters");
		if (!check_heading("!BondParams", ifsys))
			ERRORMSG("head not match with BondParams");
		
		
		for (int i = 0; i < nBondTypes; i++){
			*ifsys >> bondParams[i].bondType;
			*ifsys >> bondParams[i].atomType1;
			*ifsys >> bondParams[i].atomType2;
			*ifsys >> bondParams[i].r0;
			*ifsys >> bondParams[i].k;
			SKIP_REST
			// DEBUGMSG(bondParams[i].bondType);
			// DEBUGMSG(bondParams[i].r0);
		}
	} else { DEBUGMSG("check bond parameters"); check_heading("!BondParams", ifsys); }
	
	if(nAngleTypes){
		angleParams = new AngleParam[nAngleTypes];
		if (angleParams == NULL)
			ERRORMSG("memory allocation of angleParams error");
		// read angle parameters
		DEBUGMSG("read angle parameters");
		if(!check_heading("!AngleParams", ifsys))
			ERRORMSG("head not match with AngleParams");
		
		for (int i = 0; i < nAngleTypes; i++){
			double tmp = 0.0;
			*ifsys >> angleParams[i].angleType;
			*ifsys >> angleParams[i].atomType1;
			*ifsys >> angleParams[i].atomType2;
			*ifsys >> angleParams[i].atomType3;
			*ifsys >> tmp;
			angleParams[i].theta0 = tmp*2*PI/360.0;           // assum the unit of d1 is 'degree'
			*ifsys >> angleParams[i].k;
			SKIP_REST
			// DEBUGMSG(tmp);
		}
	} else { DEBUGMSG("check angle parameters"); check_heading("!AngleParams", ifsys); }
	
	if(nDihedralTypes){
		dihedralParams = new TorsionParam[nDihedralTypes];
		if (dihedralParams == NULL)
			ERRORMSG("memory allocation of dihedralParams error");
		// read dihedral parameters
		DEBUGMSG("read dihedral parameters");
		if(!check_heading("!DihedralParams", ifsys))
			ERRORMSG("head not match with DihedralParams");

		for (int i = 0; i < nDihedralTypes; i++){
			*ifsys >> dihedralParams[i].torsionType;
			*ifsys >> dihedralParams[i].atomType1;
			*ifsys >> dihedralParams[i].atomType2;
			*ifsys >> dihedralParams[i].atomType3;
			*ifsys >> dihedralParams[i].atomType4;
			*ifsys >> dihedralParams[i].k;
			*ifsys >> dihedralParams[i].d;
			*ifsys >> dihedralParams[i].n;
			SKIP_REST
			// DEBUGMSG(dihedralParams[i].k);
		}
	} else {
		DEBUGMSG("check dihedral parameters");
		check_heading("!DihedralParams", ifsys);
	}
	
	if(nImproperTypes){
		improperParams = new TorsionParam[nImproperTypes];
		if (improperParams == NULL)
			ERRORMSG("memory allocation of improperParams error");
		// read improper parameters
		DEBUGMSG("read improper parameters");
		if(!check_heading("!ImproperParams", ifsys))
			ERRORMSG("head not match with ImproperParams");
		
		for (int i = 0; i < nImproperTypes; i++){
			*ifsys >> improperParams[i].torsionType;
			*ifsys >> improperParams[i].atomType1;
			*ifsys >> improperParams[i].atomType2;
			*ifsys >> improperParams[i].atomType3;
			*ifsys >> improperParams[i].atomType4;
			*ifsys >> improperParams[i].k;
			*ifsys >> improperParams[i].d;
			*ifsys >> improperParams[i].n;
			SKIP_REST
		}
	} else {
		DEBUGMSG("check improper parameters");
		check_heading("!ImproperParams", ifsys);
	}
	
	if(nMolTypes){
		molTypes = new MolType[nMolTypes];
		if (molTypes == NULL)
			ERRORMSG("memory allocation of nMolTypes error");
		// read molecule types
		DEBUGMSG("read molecular types");
		if(!check_heading("!MolTypes", ifsys))
			ERRORMSG("head not match with MolTypes");

		for (int i = 0; i < nMolTypes; i++){
			*ifsys >> molTypes[i].molType;
			*ifsys >> molTypes[i].numOfAtoms;
			*ifsys >> molTypes[i].name;
			SKIP_REST
		}
	}

	if(nMolecules){
		molArray = new MolStruct[nMolecules];
		if (molArray == NULL)
			ERRORMSG("memory allocation of molArray error");
	}
	
	if(nAtoms){
		atomArray = new AtomStruct[nAtoms];
		if (atomArray == NULL)
			ERRORMSG("memory allocation of atomArray error");
		// read topology
		// read atoms and molecules
		DEBUGMSG("read atoms");
		if(!check_heading("!Atoms", ifsys))
			ERRORMSG("head not match with Atoms");
		
		for (int i = 0; i < nAtoms; i++){
			unsigned int atom_id, mol_id;
			*ifsys >> atom_id;
			atomArray[atom_id].atomID = atom_id;
			*ifsys >> atomArray[atom_id].atomType;
			*ifsys >> mol_id;
			atomArray[atom_id].molID = mol_id;
			// Here we waste a bit of time resetting molID and molTypes each time a molecule has more than one atom.
			molArray[mol_id].molID = mol_id;
			*ifsys >> molArray[mol_id].molType;
			SKIP_REST
		}
		
	}
	
	if(nBonds){
		bondArray = new Bond[nBonds];
		if (bondArray == NULL)
			ERRORMSG("memory allocation of bondArray error");
		// read bonds
		DEBUGMSG("read bonds");
		if(!check_heading("!Bonds", ifsys))
			ERRORMSG("head not match with Bonds");
		
		for (int i = 0; i < nBonds; i++){
			unsigned int bond_id;
			*ifsys >> bond_id;
			*ifsys >> bondArray[bond_id].atom1;
			*ifsys >> bondArray[bond_id].atom2;
			*ifsys >> bondArray[bond_id].bondType;
			SKIP_REST
		}
	}
	if(nAngles){
		angleArray = new Angle[nAngles];
		if (angleArray == NULL)
			ERRORMSG("memory allocation of angleArray error");
		// read angles
		DEBUGMSG("read angles");
		if(!check_heading("!Angles", ifsys))
			ERRORMSG("head not match with Angles");

		for (int i = 0; i < nAngles; i++){
			unsigned int angle_id;
			*ifsys >> angle_id;
			*ifsys >> angleArray[angle_id].atom1;
			*ifsys >> angleArray[angle_id].atom2;
			*ifsys >> angleArray[angle_id].atom3;
			*ifsys >> angleArray[angle_id].angleType;
			SKIP_REST
		}
	}

	if(nDihedrals){
		dihedralArray = new Dihedral[nDihedrals];
		if (dihedralArray == NULL)
			ERRORMSG("memory allocation of dihedralArray error");
		// read dihedrals
		DEBUGMSG("read dihedrals");
		if(!check_heading("!Dihedrals", ifsys))
			ERRORMSG("head not match with Dihedrals");
		
		for (int i = 0; i < nDihedrals; i++){
			unsigned int dihedral_id;
			*ifsys >> dihedral_id;
			*ifsys >> dihedralArray[dihedral_id].atom1;
			*ifsys >> dihedralArray[dihedral_id].atom2;
			*ifsys >> dihedralArray[dihedral_id].atom3;
			*ifsys >> dihedralArray[dihedral_id].atom4;
			*ifsys >> dihedralArray[dihedral_id].dihedralType;
			SKIP_REST
		}
	}
	if(nImpropers){
		improperArray = new Improper[nImpropers];
		if (improperArray == NULL)
			ERRORMSG("memory allocation of improperArray error");
		// read impropers
		DEBUGMSG("read impropers");
		if(!check_heading("!Impropers", ifsys))
			ERRORMSG("head not match with Impropers");
		
		for (int i = 0; i < nImpropers; i++){
			unsigned int improper_id;
			*ifsys >> improper_id;
			*ifsys >> improperArray[improper_id].atom1;
			*ifsys >> improperArray[improper_id].atom2;
			*ifsys >> improperArray[improper_id].atom3;
			*ifsys >> improperArray[improper_id].atom4;
			*ifsys >> improperArray[improper_id].improperType;
			SKIP_REST
		}
	}
	
	ljTableSize = nAtomTypes;
	if(ljTableSize){
		ljParamTable = new LJPairParam[ljTableSize*ljTableSize];
		if (ljParamTable  == NULL)
			ERRORMSG( "memory allocation of ljParamTable error");
		
		// lj_table is a triangle table. But duplicat values are set for the table.
		// This simplifies and speeds up the data retrieve in simulation
		for (int i = 0; i < ljTableSize; i++){
			for(int j = 0; j < ljTableSize; j++){
				set_lj_param(ljParamTable, atomParams[i], atomParams[j]);
			}
		}
	}
}   // end read_sys_data() 298

void Parameters::read_sys_data_test(const char* sysDataFile){
	LOG(sysDataFile)
	DEBUGMSG("Reading system data '" << sysDataFile << "'");

	ifstream *ifsys = new ifstream(sysDataFile,ios::in);
	if (!ifsys || !ifsys->is_open()){
		ERRMSG("Can't open " << sysDataFile)
	}

	// skip first two lines
	SKIP_REST
	SKIP_REST
	// next 6 lines should contain the number about molecular system sizes
	// and each line contains two numbers
	*ifsys >> nMolecules >> nMolTypes;
	SKIP_REST
	*ifsys >> nAtoms >> nAtomTypes;
	SKIP_REST
	*ifsys >> nBonds >> nBondTypes;
	SKIP_REST
	*ifsys >> nAngles >> nAngleTypes;
	SKIP_REST
	*ifsys >> nDihedrals >> nDihedralTypes;
	SKIP_REST
	*ifsys >> nImpropers >> nImproperTypes;
	SKIP_REST
	// Jc: Broadcast to other processors
	// unsigned int buf[] = {nMolecules, nMolTypes, nAtoms, nAtomTypes, nBonds, nBondTypes, nAngles, nAngleTypes, nDihedrals, nDihedralTypes, nImpropers, nImproperTypes};
	// MPI_Bcast(&buf,12,MPI_UNSIGNED,0,MPI_COMM_WORLD);

	// allocate memory to data structures
	if(nAtomTypes){
		atomParams = new AtomParam[nAtomTypes];
		if (atomParams == NULL)
			ERRORMSG("memory allocation of atomParams error");
	}
	if(nMolTypes){
		molTypes = new MolType[nMolTypes];
		if (molTypes == NULL)
			ERRORMSG("memory allocation of nMolTypes error");
	}
	if(nBondTypes){
		bondParams = new BondParam[nBondTypes];
		if (bondParams == NULL)
			ERRORMSG("memory allocation of bondParams error");
	}
	if(nAngleTypes){
		angleParams = new AngleParam[nAngleTypes];
		if (angleParams == NULL)
			ERRORMSG("memory allocation of angleParams error");
	}
	if(nDihedralTypes){
		dihedralParams = new TorsionParam[nDihedralTypes];
		if (dihedralParams == NULL)
			ERRORMSG("memory allocation of dihedralParams error");
	}
	if(nImproperTypes){
		improperParams = new TorsionParam[nImproperTypes];
		if (improperParams == NULL)
			ERRORMSG("memory allocation of improperParams error");
	}
	if(nAtoms){
		atomArray = new AtomStruct[nAtoms];
		if (atomArray == NULL)
			ERRORMSG("memory allocation of atomArray error");
	}
	if(nMolecules){
		molArray = new MolStruct[nMolecules];
		if (molArray == NULL)
			ERRORMSG("memory allocation of molArray error");
	}
	if(nBonds){
		bondArray = new Bond[nBonds];
		if (bondArray == NULL)
			ERRORMSG("memory allocation of bondArray error");
	}
	if(nAngles){
		angleArray = new Angle[nAngles];
		if (angleArray == NULL)
			ERRORMSG("memory allocation of angleArray error");
	}
	if(nDihedrals){
		dihedralArray = new Dihedral[nDihedrals];
		if (dihedralArray == NULL)
			ERRORMSG("memory allocation of dihedralArray error");
	}
	if(nImpropers){
		improperArray = new Improper[nImpropers];
		if (improperArray == NULL)
			ERRORMSG("memory allocation of improperArray error");
	}
	ljTableSize = nAtomTypes;
	if(ljTableSize){
		ljParamTable = new LJPairParam[ljTableSize*ljTableSize];
		if (ljParamTable  == NULL)
			ERRORMSG( "memory allocation of ljParamTable error");
	}

	// read atom parameters
	// check section header
	if (!check_heading("!AtomParams", ifsys))
		ERRORMSG("head not match with AtomParams");

	if(nAtomTypes){
		DEBUGMSG("read atom parameters");
		for (int i = 0; i < nAtomTypes; i++){
			AtomParam ap = atomParams[i];
			// oneliner
			// *ifsys >> ap.atomType >> ap.mass >> ap.charge >> ap.eps >> ap.sigma >> ap.dipx >> ap.dipy >> ap.dipz >> ap.quadxx >> ap.quadxy >> ap.quadxz >> ap.quadyy >> ap.quadyz >> ap.quadzz >> ap.polar  >> ap.typeName >> ap.atomName;
			// SKIP_REST
			*ifsys >> ap.atomType;
			*ifsys >> ap.mass;
			*ifsys >> ap.charge;
			*ifsys >> ap.eps;
			*ifsys >> ap.sigma;
			*ifsys >> ap.dipx;
			*ifsys >> ap.dipy;
			*ifsys >> ap.dipz;
			*ifsys >> ap.quadxx;
			*ifsys >> ap.quadxy;
			*ifsys >> ap.quadxz;
			*ifsys >> ap.quadyy;
			*ifsys >> ap.quadyz;
			*ifsys >> ap.quadzz;
			*ifsys >> ap.polar ;
			*ifsys >> ap.typeName;
			*ifsys >> ap.atomName;
			// DEBUGMSG(ap.atomType);
			// DEBUGMSG(ap.mass);
			SKIP_REST
		}
	}

	// read bond parameters
	if (!check_heading("!BondParams", ifsys))
		ERRORMSG("head not match with BondParams");

	if(nBondTypes){
		DEBUGMSG("read bond parameters");
		for (int i = 0; i < nBondTypes; i++){
			*ifsys >> bondParams[i].bondType;
			*ifsys >> bondParams[i].atomType1;
			*ifsys >> bondParams[i].atomType2;
			*ifsys >> bondParams[i].r0;
			*ifsys >> bondParams[i].k;
			SKIP_REST
			// DEBUGMSG(bondParams[i].bondType);
			// DEBUGMSG(bondParams[i].r0);
		}
	}

	// read angle parameters
	if(!check_heading("!AngleParams", ifsys))
		ERRORMSG("head not match with AngleParams");

	if(nAngleTypes){
		DEBUGMSG("read angle parameters");
		for (int i = 0; i < nAngleTypes; i++){
			double tmp = 0.0;
			*ifsys >> angleParams[i].angleType;
			*ifsys >> angleParams[i].atomType1;
			*ifsys >> angleParams[i].atomType2;
			*ifsys >> angleParams[i].atomType3;
			*ifsys >> tmp;
			angleParams[i].theta0 = tmp*2*PI/360.0;           // assum the unit of d1 is 'degree'
			*ifsys >> angleParams[i].k;
			SKIP_REST
			// DEBUGMSG(tmp);
		}
	}

	// read dihedral parameters
	if(!check_heading("!DihedralParams", ifsys))
		ERRORMSG("head not match with DihedralParams");

	if(nDihedralTypes){
		DEBUGMSG("read dihedral parameters");
		for (int i = 0; i < nDihedralTypes; i++){
			*ifsys >> dihedralParams[i].torsionType;
			*ifsys >> dihedralParams[i].atomType1;
			*ifsys >> dihedralParams[i].atomType2;
			*ifsys >> dihedralParams[i].atomType3;
			*ifsys >> dihedralParams[i].atomType4;
			*ifsys >> dihedralParams[i].k;
			*ifsys >> dihedralParams[i].d;
			*ifsys >> dihedralParams[i].n;
			SKIP_REST
			// DEBUGMSG(dihedralParams[i].k);
		}
	}

	// read improper parameters
	if(!check_heading("!ImproperParams", ifsys))
		ERRORMSG("head not match with ImproperParams");
	if(nImproperTypes){
		DEBUGMSG("read improper parameters");
		for (int i = 0; i < nImproperTypes; i++){
			*ifsys >> improperParams[i].torsionType;
			*ifsys >> improperParams[i].atomType1;
			*ifsys >> improperParams[i].atomType2;
			*ifsys >> improperParams[i].atomType3;
			*ifsys >> improperParams[i].atomType4;
			*ifsys >> improperParams[i].k;
			*ifsys >> improperParams[i].d;
			*ifsys >> improperParams[i].n;
			SKIP_REST
		}
	}

	// read molecule types
	if(!check_heading("!MolTypes", ifsys))
		ERRORMSG("head not match with MolTypes");

	if(nMolTypes){
		DEBUGMSG("read molecular types");
		for (int i = 0; i < nMolTypes; i++){
			*ifsys >> molTypes[i].molType;
			*ifsys >> molTypes[i].numOfAtoms;
			*ifsys >> molTypes[i].name;
			SKIP_REST
		}
	}

	// read topology
	// read atoms and molecules
	if(!check_heading("!Atoms", ifsys))
		ERRORMSG("head not match with Atoms");
	if (nAtoms){
		DEBUGMSG("read atoms");
		for (int i = 0; i < nAtoms; i++){
			unsigned int atom_id, mol_id;
			*ifsys >> atom_id;
			atomArray[atom_id].atomID = atom_id;
			*ifsys >> atomArray[atom_id].atomType;
			*ifsys >> mol_id;
			atomArray[atom_id].molID = mol_id;
			// Here we waste a bit of time resetting molID and molTypes each time a molecule has more than one atom.
			molArray[mol_id].molID = mol_id;
			*ifsys >> molArray[mol_id].molType;
			SKIP_REST
		}
	}

	// read bonds
	if(!check_heading("!Bonds", ifsys))
		ERRORMSG("head not match with Bonds");
	
	if(nBonds){
		DEBUGMSG("read bonds");
		for (int i = 0; i < nBonds; i++){
			unsigned int bond_id;
			*ifsys >> bond_id;
			*ifsys >> bondArray[bond_id].atom1;
			*ifsys >> bondArray[bond_id].atom2;
			*ifsys >> bondArray[bond_id].bondType;
			SKIP_REST
		}
	}

	// read angles
	if(!check_heading("!Angles", ifsys))
		ERRORMSG("head not match with Angles");
	
	if(nAngles){
		DEBUGMSG("read angles");
		for (int i = 0; i < nAngles; i++){
			unsigned int angle_id;
			*ifsys >> angle_id;
			*ifsys >> angleArray[angle_id].atom1;
			*ifsys >> angleArray[angle_id].atom2;
			*ifsys >> angleArray[angle_id].atom3;
			*ifsys >> angleArray[angle_id].angleType;
			SKIP_REST
		}
	}

	// read dihedrals
	if(!check_heading("!Dihedrals", ifsys))
		ERRORMSG("head not match with Dihedrals");

	if(nDihedrals){
		DEBUGMSG("read dihedrals");
		for (int i = 0; i < nDihedrals; i++){
			unsigned int dihedral_id;
			*ifsys >> dihedral_id;
			*ifsys >> dihedralArray[dihedral_id].atom1;
			*ifsys >> dihedralArray[dihedral_id].atom2;
			*ifsys >> dihedralArray[dihedral_id].atom3;
			*ifsys >> dihedralArray[dihedral_id].atom4;
			*ifsys >> dihedralArray[dihedral_id].dihedralType;
			SKIP_REST
		}
	}

	// read impropers
	if(!check_heading("!Impropers", ifsys))
		ERRORMSG("head not match with Impropers");

	if(nImpropers){
		DEBUGMSG("read impropers");
		for (int i = 0; i < nImpropers; i++){
			unsigned int improper_id;
			*ifsys >> improper_id;
			*ifsys >> improperArray[improper_id].atom1;
			*ifsys >> improperArray[improper_id].atom2;
			*ifsys >> improperArray[improper_id].atom3;
			*ifsys >> improperArray[improper_id].atom4;
			*ifsys >> improperArray[improper_id].improperType;
			SKIP_REST
		}
	}
	build_lj_param_table();
}   // end read_sys_data() 298


void Parameters::read_sys_data_orig(const char* sysDataFile){
	int  i, ret, v1, v2, v3, v4, v5, v6;
	double d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14;
	char buf[256];
	char str1[32];
	char str2[32];
	FILE *fptr;
	DEBUGMSG("Reading system data");

	if(mpi_rank_g == 0){
		if ((fptr = fopen(sysDataFile, "r")) == NULL)
			ERRORMSG("open sysDataFile error");

		// skip first two lines
		fgets(buf, 512, fptr);
		fgets(buf, 512, fptr);
		// next 6 lines should contain the number about molecular system sizes
		// and each line contains two numbers
		fgets(buf, 512, fptr);
		if((sscanf(buf, "%d%d", &nMolecules, &nMolTypes)) != 2)
			ERRORMSG("read nMolecules & nMolTypes error ");

		fgets(buf, 512, fptr);
		if((sscanf(buf, "%d%d", &nAtoms, &nAtomTypes)) != 2)
			ERRORMSG("read nAtom & nAtomTypes error \n");

		fgets(buf, 512, fptr);
		if((sscanf(buf, "%d%d", &nBonds, &nBondTypes)) != 2)
			ERRORMSG("read nBonds & nBondTypes error \n");

		fgets(buf, 512, fptr);
		if((sscanf(buf, "%d%d", &nAngles, &nAngleTypes)) != 2)
			ERRORMSG("read nAngles & nAngleTypes error \n");

		fgets(buf, 512, fptr);
		if((sscanf(buf, "%d%d", &nDihedrals, &nDihedralTypes)) != 2)
			ERRORMSG("read nDihedrals & nDihedralTypes error");

		fgets(buf, 512, fptr);
		if((sscanf(buf, "%d%d", &nImpropers, &nImproperTypes)) != 2)
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

	// allocate memory to data structures
	if(nAtomTypes){
		atomParams = new AtomParam[nAtomTypes];
		if (atomParams == NULL)
			ERRORMSG("memory allocation of atomParams error");
	}
	if(nMolTypes){
		molTypes = new MolType[nMolTypes];
		if (molTypes == NULL)
			ERRORMSG("memory allocation of nMolTypes error");
	}
	if(nBondTypes){
		bondParams = new BondParam[nBondTypes];
		if (bondParams == NULL)
			ERRORMSG("memory allocation of bondParams error");
	}

	if(nAngleTypes){
		angleParams = new AngleParam[nAngleTypes];
		if (angleParams == NULL)
			ERRORMSG("memory allocation of angleParams error");
	}
	if(nDihedralTypes){
		dihedralParams = new TorsionParam[nDihedralTypes];
		if (dihedralParams == NULL)
			ERRORMSG("memory allocation of dihedralParams error");
	}
	if(nImproperTypes){
		improperParams = new TorsionParam[nImproperTypes];
		if (improperParams == NULL)
			ERRORMSG("memory allocation of improperParams error");
	}
	if(nAtoms){
		atomArray = new AtomStruct[nAtoms];
		if (atomArray == NULL)
			ERRORMSG("memory allocation of atomArray error");
	}
	if(nMolecules){
		molArray = new MolStruct[nMolecules];
		if (molArray == NULL)
			ERRORMSG("memory allocation of molArray error");
	}
	if(nBonds){
		bondArray = new Bond[nBonds];
		if (bondArray == NULL)
			ERRORMSG("memory allocation of bondArray error");
	}
	if(nAngles){
		angleArray = new Angle[nAngles];
		if (angleArray == NULL)
			ERRORMSG("memory allocation of angleArray error");
	}
	if(nDihedrals){
		dihedralArray = new Dihedral[nDihedrals];
		if (dihedralArray == NULL)
			ERRORMSG("memory allocation of dihedralArray error");
	}
	if(nImpropers){
		improperArray = new Improper[nImpropers];
		if (improperArray == NULL)
			ERRORMSG("memory allocation of improperArray error");
	}

	ljTableSize = nAtomTypes;
	if(ljTableSize){
		ljParamTable = new LJPairParam[ljTableSize*ljTableSize];
		if (ljParamTable  == NULL)
			ERRORMSG( "memory allocation of ljParamTable error");
	}

	// read atom parameters
	if(mpi_rank_g == 0){
		if (!check_heading_orig("!AtomParams", fptr))
			ERRORMSG("head not match with AtomParams");
		if(nAtomTypes){
			DEBUGMSG("read atom parameters");
			for (int i = 0; i < nAtomTypes; i++){
				fgets(buf, 512, fptr);
				for(int desRank = 1;desRank<mpi_size_g;desRank++){
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
	else {
		for (int i = 0; i < nAtomTypes; i++){
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
	if(mpi_rank_g == 0){
		if (!check_heading_orig("!BondParams", fptr))
			ERRORMSG("head not match with BondParams");
		if(nBondTypes){
			DEBUGMSG("read bond parameters");
			for (int i = 0; i < nBondTypes; i++){
				fgets(buf, 512, fptr);
				for(int desRank = 1;desRank<mpi_size_g;desRank++){
					//            MPI::COMM_WORLD.Send(buf,256,MPI::CHAR,desRank,0);
					MPI_Send(buf,256,MPI_CHAR,desRank,0,MPI_COMM_WORLD);
				}
				ret = sscanf(buf, "%d%d%d%lf%lf", &v1, &v2, &v3, &d1, &d2);

				bondParams[i].bondType = v1;
				bondParams[i].atomType1 = v2;
				bondParams[i].atomType2 = v3;
				bondParams[i].r0 = d1;
				bondParams[i].k = d2;
			}
		}
	}
	else {
		for (int i = 0; i < nBondTypes; i++){
			//        MPI::COMM_WORLD.Recv(buf,256,MPI::CHAR,0,0);
			MPI_Recv(buf,256,MPI_CHAR,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

			ret = sscanf(buf, "%d%d%d%lf%lf", &v1, &v2, &v3, &d1, &d2);

			bondParams[i].bondType = v1;
			bondParams[i].atomType1 = v2;
			bondParams[i].atomType2 = v3;
			bondParams[i].r0 = d1;
			bondParams[i].k = d2;
		}
	}

	// read angle parameters

	if(mpi_rank_g == 0){
		if(!check_heading_orig("!AngleParams", fptr))
			ERRORMSG("head not match with AngleParams");

		if(nAngleTypes){
			DEBUGMSG("read angle parameters");

			for (int i = 0; i < nAngleTypes; i++){
				fgets(buf, 512, fptr);
				for(int desRank = 1;desRank<mpi_size_g;desRank++){
					//            MPI::COMM_WORLD.Send(buf,256,MPI::CHAR,desRank,0);
					MPI_Send(buf,256,MPI_CHAR,desRank,0,MPI_COMM_WORLD);
				}
				ret= sscanf(buf, "%d%d%d%d%lf%lf", &v1, &v2, &v3, &v4, &d1,&d2);

				angleParams[i].angleType = v1;
				angleParams[i].atomType1 = v2;
				angleParams[i].atomType2 = v3;
				angleParams[i].atomType3 = v4;
				angleParams[i].theta0 = d1*2*PI/360.0;           // assum the unit of d1 is 'degree'
				angleParams[i].k = d2;
			}
		}
	}
	else {
		for (int i = 0; i < nAngleTypes; i++){
			//        MPI::COMM_WORLD.Recv(buf,256,MPI::CHAR,0,0);
			MPI_Recv(buf,256,MPI_CHAR,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);   // Not sure about Ignore
			ret = sscanf(buf, "%d%d%d%d%lf%lf", &v1, &v2, &v3, &v4, &d1,&d2);

			angleParams[i].angleType = v1;
			angleParams[i].atomType1 = v2;
			angleParams[i].atomType2 = v3;
			angleParams[i].atomType3 = v4;
			angleParams[i].theta0 = d1*2*PI/360.0;           // assum the unit of d1 is 'degree'
			angleParams[i].k = d2;
		}
	}

	// read dihedral parameters
	if(mpi_rank_g == 0){
		if(!check_heading_orig("!DihedralParams", fptr))
			ERRORMSG("head not match with DihedralParams");

		if(nDihedralTypes){

			DEBUGMSG("read dihedral parameters");


			for (int i = 0; i < nDihedralTypes; i++){
				fgets(buf, 512, fptr);
				for(int desRank = 1;desRank<mpi_size_g;desRank++){
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
	else {
		for (int i = 0; i < nDihedralTypes; i++){
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
	if(mpi_rank_g == 0){
		if(!check_heading_orig("!ImproperParams", fptr))
			ERRORMSG("head not match with ImproperParams");
		if(nImproperTypes){
			DEBUGMSG("read improper parameters");

			for (int i = 0; i < nImproperTypes; i++){
				fgets(buf, 512, fptr);
				for(int desRank = 1;desRank<mpi_size_g;desRank++){
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
	else {
		for (int i = 0; i < nImproperTypes; i++){
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
	if(mpi_rank_g == 0){
		if(!check_heading_orig("!MolTypes", fptr))
			ERRORMSG("head not match with MolTypes");

		if(nMolTypes){

			DEBUGMSG("read molecular types");


			for (int i = 0; i < nMolTypes; i++){
				fgets(buf, 512, fptr);
				for(int desRank = 1;desRank<mpi_size_g;desRank++){
					//            MPI::COMM_WORLD.Send(buf,256,MPI::CHAR,desRank,0);
					MPI_Send(buf,256,MPI_CHAR,desRank,0,MPI_COMM_WORLD);
				}

				sscanf(buf, "%d%d%s", &v1, &v2, str1);


				molTypes[i].molType = v1;
				molTypes[i].numOfAtoms = v2;
				strcpy(molTypes[i].name, str1);
			}
		}
	}
	else {
		for (int i = 0; i < nMolTypes; i++){
			//        MPI::COMM_WORLD.Recv(buf,256,MPI::CHAR,0,0);
			MPI_Recv(buf,256,MPI_CHAR,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

			sscanf(buf, "%d%d%s", &v1, &v2, str1);


			molTypes[i].molType = v1;
			molTypes[i].numOfAtoms = v2;
			strcpy(molTypes[i].name, str1);
		}
	}

	// read topology
	// read atoms and molecules
	if(mpi_rank_g == 0){
		if(!check_heading_orig("!Atoms", fptr))
			ERRORMSG("head not match with Atoms");
		if (nAtoms){

			DEBUGMSG("read atoms");


			for (int i = 0; i < nAtoms; i++){
				fgets(buf, 512, fptr);

				for(int desRank = 1;desRank<mpi_size_g;desRank++){
					//            MPI::COMM_WORLD.Send(buf,256,MPI::CHAR,desRank,0);
					MPI_Send(buf,256,MPI_CHAR,desRank,0,MPI_COMM_WORLD);
				}

				sscanf(buf, "%d%d%d%d", &v1, &v2, &v3, &v4);

				atomArray[v1].atomID = v1;
				atomArray[v1].atomType = v2;
				atomArray[v1].molID = v3;
				molArray[v3].molID = v3;
				molArray[v3].molType = v4;
			}
		}
	}
	else {
		for (int i = 0; i < nAtoms; i++){
			//        MPI::COMM_WORLD.Recv(buf,256,MPI::CHAR,0,0);
			MPI_Recv(buf,256,MPI_CHAR,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			sscanf(buf, "%d%d%d%d", &v1, &v2, &v3, &v4);
			atomArray[v1].atomID = v1;
			atomArray[v1].atomType = v2;
			atomArray[v1].molID = v3;
			molArray[v3].molID = v3;
			molArray[v3].molType = v4;
		}
	}
	//    MPI::COMM_WORLD.Barrier();
	// read bonds
	if(mpi_rank_g == 0){
		if(!check_heading_orig("!Bonds", fptr))
			ERRORMSG("head not match with Bonds");
		if(nBonds){

			DEBUGMSG("read bonds");


			for (int i = 0; i < nBonds; i++){
				fgets(buf, 512, fptr);
				for(int desRank = 1;desRank<mpi_size_g;desRank++){
					//            MPI::COMM_WORLD.Send(buf,256,MPI::CHAR,desRank,0);
					MPI_Send(buf,256,MPI_CHAR,desRank,0,MPI_COMM_WORLD);
				}

				sscanf(buf, "%d%d%d%d", &v1, &v2, &v3, &v4);

				bondArray[v1].atom1 = v2;
				bondArray[v1].atom2 = v3;
				bondArray[v1].bondType = v4;
			}
		}
	}
	else {
		for (int i = 0; i < nBonds; i++){
			//          MPI::COMM_WORLD.Recv(buf,256,MPI::CHAR,0,0);
			MPI_Recv(buf,256,MPI_CHAR,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			sscanf(buf, "%d%d%d%d", &v1, &v2, &v3, &v4);

			bondArray[v1].atom1 = v2;
			bondArray[v1].atom2 = v3;
			bondArray[v1].bondType = v4;
		}
	}

	// read angles
	if(mpi_rank_g == 0){
		if(!check_heading_orig("!Angles", fptr))
			ERRORMSG("head not match with Angles");
		if(nAngles){

			DEBUGMSG("read angles");


			for (int i = 0; i < nAngles; i++){
				fgets(buf, 512, fptr);
				for(int desRank = 1;desRank<mpi_size_g;desRank++){
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
	else {
		for (int i = 0; i < nAngles; i++){
			//          MPI::COMM_WORLD.Recv(buf,256,MPI::CHAR,0,0);
			MPI_Recv(buf,256,MPI_CHAR,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			sscanf(buf, "%d%d%d%d%d", &v1, &v2, &v3, &v4, &v5);

			angleArray[v1].atom1 = v2;
			angleArray[v1].atom2 = v3;
			angleArray[v1].atom3 = v4;
			angleArray[v1].angleType = v5;
		}
	}



	// read dihedrals
	if(mpi_rank_g == 0){
		if(!check_heading_orig("!Dihedrals", fptr))
			ERRORMSG("head not match with Dihedrals");
		if(nDihedrals){

			DEBUGMSG("read dihedrals");

			for (int i = 0; i < nDihedrals; i++){
				fgets(buf, 512, fptr);
				for(int desRank = 1;desRank<mpi_size_g;desRank++){
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
	else {
		for (int i = 0; i < nDihedrals; i++){
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


	// read impropers

	if(mpi_rank_g == 0)    {
		if(!check_heading_orig("!Impropers", fptr))
			ERRORMSG("head not match with Impropers");
		if(nImpropers){

			DEBUGMSG("read impropers");


			for (int i = 0; i < nImpropers; i++){
				fgets(buf, 512, fptr);
				for(int desRank = 1;desRank<mpi_size_g;desRank++){
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
	else {
		for (int i = 0; i < nImpropers; i++){
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

	if(mpi_rank_g == 0)
		fclose(fptr);

	build_lj_param_table();

}   // end read_sys_data()


// check_heading checks section head. if it is right return 1 else return 0
// error will be handled by the caller.
bool Parameters::check_heading(const char* heading, ifstream* ifsys) {
	// skip a blank line
	ifsys->ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	// check section head
	string s;
	*ifsys >> s;
	// skip the rest of the line
	ifsys->ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	LOG(s);
	if (s != heading) {
		return false;
	}
	// skip following notation line
	ifsys->ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	return true;
}



// check_heading checks section head. if it is right return 1 else return 0
// error will be handled by the caller.
int Parameters::check_heading_orig(const char* heading, FILE* fptr) {
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

void Parameters::build_lj_param_table(){
	DEBUGMSG("build LJ table ");
	if((ljTableSize == 0) || (ljParamTable == NULL))
		ERRORMSG("build_lj_param_table error");

	// lj_table is a triangle table. But duplicat values are set for the table.
	// This simplifies and speeds up the data retrieve in simulation
	for (int i = 0; i < ljTableSize; i++){
		for(int j = 0; j < ljTableSize; j++){
			set_lj_param(ljParamTable, atomParams[i], atomParams[j]);
		}
	}
}

void Parameters::set_lj_param(LJPairParam *ljTable, const AtomParam &ap1,const AtomParam &ap2){
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
