/** SimConfiguration.cpp -- implementation to handle input, maintain the global data structure
 **     and build the simulation system
 **
 ** Copyright (C) 2003
 ** Centre for Molecular SimConfiguration (CMS)
 ** School of Information Technology
 ** Swinburne University of Technology
 ** PO Box 21 Hawthorn, Vic 3122, Australia
 **
 ** Author: Zhongwu Zhou
 ** Email: zzhou@it.swin.edu.au
 **/

#include "SimConfiguration.h"
// #include <mpi.h>
#include "Errors.h"
#include "Defaults.h"

SimConfiguration::SimConfiguration(const char configFile[])
{

	if (configFile == NULL)
		ERRORMSG("missing simulation configuration file");
	sysConfigFile = configFile;
	// Initialise data members [Good luck!]
	init();
}

void SimConfiguration::init()
{
	nSteps = 1;   startSteps = 0;   timeStep = 1.0;
	shortTsFreq = 1;                longTsFreq = 4;
	averageTsFreq = nSteps;         printTsFreq = nSteps;
	backupTsFreq = nSteps;
	startSampling = startSteps;     samplingTsFreq = nSteps;
	startTTCF = startSteps;         ttcfTsFreq = nSteps;
	integratorType = 0;
	nBins = 100;
	couplthermo = 1.0;              couplBaros = 1.0;

	lxBox = 1.0;         lyBox = 1.0;         lzBox = 1.0;
	cutoff = 1.2;        kCutoff = 5;         cutBuff = 0.0;
	cutoffEw = 1.0;
	rate = 0.0;          cutoffInd = 1.0;     adamp = 200.0;
	density = 1.0;       limitRDF = 1.0;      switchDist = 1.0;
	nMie = 12.0;	 mMie = 6.0;
	nBuff = 14.0;	 mBuff = 7.0;         Hfactor = 0.91;
	deltaBuff = 0.07;	 gammaBuff = 0.12;
	cellSideLen  =1.0;   elapsedTime = 0.0;   temperature = 1.0;
	CFB = 0.0;           DFB = 0.0;        //   pressure  = 1.0;
	tolerance = 0.00001;    maxCount = 10000;   seed = 1;
	klUB = 0.0;          l0UB = 0.0;
	bendingType = 1;

	printSysConfig = false;   isNewStart = false;     zeroAverage = true;
	doEquilibrium = true;     doCompression = false;  doShear = false;
	doElongation = false;     computeRDF = false;     computeMSD = false;
	computeVACF = false;      useTTCF = false;        useCellPairList = false;
	binaryOutput = false;     constraintOn = false;   useAtomThermo = true;
	switchOn = false;         lCorrectOn = false;
	computeCoulomb = false;   computeEwald = false;   computeSurfCorrect = false;
	computeWolf = false;      computeMie = false;     computeMultiPol = false;
	computeInduction = false; dampReal = false;       dampInduction = false;
	dampAmoeba = false;
	doInductionMin = false;   computeTorq = false;    computeBuffvdW = false;
	hasExtForces = false;     hasFixedAtom = false;   computeharmoUB = false;
	constantVolume = false;   constantPressure = false;
	constantEnergy = false;   constantTemperature = false;
	readVelocity = false;     writeTrajectory = false;   computeLustig = false;
	harmonicbend = false;     amoebabend = false;

	coordinateFile[0] = '\0';
	velocityFile[0] = '\0';
	sysDataFile[0] = '\0';;
	resultFile[0] = '\0';
	restartFile[0] = '\0';
	trajectoryFile[0] = '\0';
	dumpFile[0] = '\0';
	molTrajectoryFile[0] = '\0';
	molVelocityFile[0] = '\0';
	thermotype[0] = '\0';

	DEBUGMSG("SimConfiguration initialised");
}

#define SKIP_REST config_s->ignore(std::numeric_limits<std::streamsize>::max(), '\n');

void SimConfiguration::read_config_file(){
	LOG("Opening config file: " << sysConfigFile);
	ifstream *config_s = new ifstream(sysConfigFile,ios::in);
	if (!config_s){
		ERRMSG("Can't open "<< sysConfigFile)
	}

	*config_s >> sysDataFile >> coordinateFile;
	SKIP_REST

	*config_s >> resultFile >> restartFile >> dumpFile;
	SKIP_REST

	string boolFlag;
	*config_s >> boolFlag; // must output trajectory file
	if (boolFlag == "t") {
		writeTrajectory = true;
		*config_s >> trajectoryFile;
	}
	SKIP_REST

	*config_s >> boolFlag; // isNewStart
	if (boolFlag == "t") {
		isNewStart = true;
	}

	*config_s >> seed;

	*config_s >> boolFlag; // readVelocity
	if (boolFlag == "t") {
		readVelocity = true;
	}
	*config_s >> velocityFile;
	SKIP_REST

	*config_s >> timeStep >> nSteps >> startSteps;
	SKIP_REST

	*config_s >> shortTsFreq >> longTsFreq;
	SKIP_REST

	*config_s >> averageTsFreq >> printTsFreq >> backupTsFreq;
	SKIP_REST

	*config_s >> integratorType;
	SKIP_REST

	int ensStatus,atomFlag;
	double ensParam1,ensParam2;

	*config_s >> ensStatus >> ensParam1 >> ensParam2 >> atomFlag;
	set_ensemble_status(ensStatus, ensParam1, ensParam2);
	if (atomFlag != 1){
		useAtomThermo = false;
	}
	*config_s >> couplthermo >> thermotype >> couplBaros;
	SKIP_REST

	int runType_t;
	*config_s >> runType_t >> rate;
	set_run_type(runType_t);
	SKIP_REST

	*config_s >> cutoff >> cutoffEw >> kCutoff >> cutBuff;
	SKIP_REST

	*config_s >> density;
	SKIP_REST

	*config_s >> lxBox >> lyBox >> lzBox;
	SKIP_REST

	*config_s >> boolFlag >> cellSideLen;
	if (boolFlag == "t") {
		useCellPairList = true;
	}
	SKIP_REST

	int constraintMode;
	*config_s >> constraintMode;

	if (constraintMode == GAUSSIAN_CONSTRAINT){ // Gaussan constraint
		constraintOn = true;
		*config_s >> CFB >> DFB >> TFB;
		SKIP_REST
	} else if (constraintMode == SHAKE_CONSTRAINT){
		constraintOn = true;
		*config_s >> tolerance >> maxCount;
		SKIP_REST
	} else if (constraintMode == NO_CONSTRAINT){
		constraintOn = false;
		*config_s >> tolerance >> bendingType;
		if (bendingType == 1){
			harmonicbend = true;
		}
		else if (bendingType == 2){
			amoebabend = true;
		}
		SKIP_REST
	}
	*config_s >> boolFlag >> klUB >> l0UB;

	if (boolFlag == "t"){
		computeharmoUB = true;
	}
	SKIP_REST

	*config_s >> boolFlag;
	if (boolFlag == "t"){
		lCorrectOn = true;
	}
	*config_s >> boolFlag;
	if (boolFlag == "t"){
		switchOn = true;
	}
	*config_s >> switchDist;
	SKIP_REST

	*config_s >> boolFlag;
	if (boolFlag == "t"){
		computeCoulomb = true;
	}
	*config_s >> boolFlag;
	if (boolFlag == "t"){
		computeEwald = true;
	}
	*config_s >> boolFlag;
	if (boolFlag == "t"){
		computeSurfCorrect = true;
	}
	*config_s >> boolFlag;
	if (boolFlag == "t"){
		computeWolf = true;
	}
	*config_s >> alpha >> kCutoff;
	SKIP_REST

	*config_s >> boolFlag;
	if (boolFlag == "t"){
		computeMultiPol  = true;
	}
	*config_s >> boolFlag;
	if (boolFlag == "t"){
		computeInduction = true;
	}
	*config_s >> boolFlag;
	if (boolFlag == "t"){
		computeTorq      = true;
	}
	*config_s >> boolFlag;
	if (boolFlag == "t"){
		dampInduction    = true;
	}
	*config_s >> boolFlag;
	if (boolFlag == "t"){
		doInductionMin   = true;
	}
	*config_s >> boolFlag;
	if (boolFlag == "t"){
		dampAmoeba       = true;
	}
	*config_s >> cutoffInd >> adamp;
	SKIP_REST

	*config_s >> boolFlag;
	if (boolFlag == "t"){
		computeBuffvdW = true;
	}
	*config_s >> nBuff >> mBuff >> deltaBuff >> gammaBuff >> Hfactor;
	SKIP_REST

	*config_s >> boolFlag;
	if (boolFlag == "t"){
		computeMie = true;
	}
	*config_s >> nMie >> mMie;
	SKIP_REST

	*config_s >> boolFlag;
	if (boolFlag == "t"){
		hasExtForces = true;
	}
	*config_s >> boolFlag;
	if (boolFlag == "t"){
		hasFixedAtom = true;
	}
	SKIP_REST
	
	
	*config_s >> startSampling >> samplingTsFreq;
	SKIP_REST

	*config_s >> boolFlag;
	if (boolFlag == "t"){
		useTTCF = true;
	}
	*config_s >> startTTCF >> ttcfTsFreq;
	SKIP_REST

	*config_s >> boolFlag;
	if (boolFlag == "t"){
		computeRDF = true;
	}
	*config_s >> nBins;
	SKIP_REST

	*config_s >> boolFlag;
	if (boolFlag == "t"){
		computeMSD = true;
	}
	*config_s >> molTrajectoryFile;
	SKIP_REST
	
	*config_s >> boolFlag;
	if (boolFlag == "t"){
		computeVACF = true;
	}
	*config_s >> molVelocityFile;
	SKIP_REST
	
	*config_s >> boolFlag;
	if (boolFlag == "t"){
		computeLustig = true;
	}
	SKIP_REST
	
	DEBUGMSG("close config file" );
}


// config.txt data is read by every process
void SimConfiguration::read_config_file_old(){
	int  v1, v2;
	double d1, d2;
	char buf[256];
	char str1[32];
	char str2[32];
	char str3[32];
	char str4[32];
	char str5[32];
	char str6[32];
	// int rank,size; // rank and size of the MPI

	FILE *fptr;

	//  setup rank, size for MPI
	//    rank = MPI::COMM_WORLD.Get_rank();
	//    size = MPI::COMM_WORLD.Get_size();

	//  process 0 reading data from config.txt
	//    if(rank == 0) {
	// cerr << "Opening config file: " << sysConfigFile << endl;
	LOG("Opening config file: " << sysConfigFile);
	ifstream *config_s = new ifstream(sysConfigFile,ios::in);
	if (!config_s){
		ERRMSG("Can't open "<< sysConfigFile)
	}
	if ((fptr = fopen(sysConfigFile, "r")) == NULL)
		ERRORMSG("fail to open configuration file");

	char   sysDataFile_t[80];
	char   coordinateFile_t[80];

	*config_s >> sysDataFile_t >> coordinateFile_t;
	config_s->ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	fgets(buf, 512, fptr);
	sscanf(buf, "%s%s", sysDataFile, coordinateFile);
	DEBUGMSG("sscan: '" << sysDataFile << "' _t: '" << sysDataFile_t << "'");
	DEBUGMSG("sscan: '" << coordinateFile << "' _t: '" << coordinateFile_t << "'");

	if(strcmp(sysDataFile, sysDataFile_t) != 0 || strcmp(coordinateFile, coordinateFile_t) != 0){
		ERRMSG("Difference between sscan and >>!")
	}

	fgets(buf, 512, fptr);
	sscanf(buf, "%s%s%s", resultFile, restartFile, dumpFile);

	char   resultFile_t[80];
	char   restartFile_t[80];
	char   dumpFile_t[80];
	*config_s >> resultFile_t >> restartFile_t >> dumpFile_t;
	config_s->ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	DEBUGMSG("sscan: '" << resultFile << "' _t: '" << resultFile_t << "'");
	DEBUGMSG("sscan: '" << restartFile << "' _t: '" << restartFile_t << "'");
	DEBUGMSG("sscan: '" << dumpFile << "' _t: '" << dumpFile_t << "'");
	if(strcmp(resultFile, resultFile_t) != 0 || strcmp(restartFile, restartFile_t) != 0 || strcmp(dumpFile,dumpFile_t) != 0){
		ERRMSG("Difference between sscan and >>!")
	}

	char   trajectoryFile_t[80];
	string s;
	*config_s >> s;
	if (s == "t") {
		*config_s >> trajectoryFile_t;
	}
	config_s->ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	fgets(buf, 512, fptr);
	sscanf(buf, "%s", str1);
	DEBUGMSG("sscan: '" << str1 << "' _t: '" << s << "'");
	if (str1[0] == 't'){
		writeTrajectory = true;
		sscanf(buf, "%s", trajectoryFile);
		DEBUGMSG("sscan: '" << trajectoryFile << "' _t: '" << trajectoryFile_t << "'");
		if(strcmp(trajectoryFile, trajectoryFile_t) != 0){
			ERRMSG("Difference between sscan and >>!")
		}
	}
	if( s != (std::string) str1 ){
		ERRMSG("Difference between sscan and >>!")
	}


	fgets(buf, 512, fptr);
	sscanf(buf, "%s%u%s%s", str1, &seed, str2, velocityFile);
	if (str1[0] == 't') isNewStart = true;
	if (str2[0] == 't') readVelocity = true;

	*config_s >> s;
	DEBUGMSG("sscan: '" << str1 << "' _t: '" << s << "'");
	if( s != (std::string) str1 ){
		ERRMSG("Difference between sscan and >>!")
	}
	unsigned seed_t;
	*config_s >> seed_t;
	DEBUGMSG("sscan: '" << seed << "' _t: '" << seed_t << "'");
	if (seed !=seed_t) {
		ERRMSG("Difference between sscan and >>!")
	}

	*config_s >> s;
	DEBUGMSG("sscan: '" << str2 << "' _t: '" << s << "'");
	if( s != (std::string) str2 ){
		ERRMSG("Difference between sscan and >>!")
	}
	char   velocityFile_t[80];
	*config_s >> velocityFile_t;
	config_s->ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	DEBUGMSG("sscan: '" << velocityFile << "' _t: '" << velocityFile_t << "'");
	if(strcmp(velocityFile, velocityFile_t) != 0){
		ERRMSG("Difference between sscan and >>!")
	}

	fgets(buf, 512, fptr);
	sscanf(buf, "%lf%d%d", &timeStep, &nSteps, &startSteps);
	Double timeStep_t;
	Int nSteps_t,	startSteps_t;
	*config_s >> timeStep_t >> nSteps_t >> startSteps_t;
	config_s->ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	DEBUGMSG("sscan: '" << timeStep << "' _t: '" << timeStep_t << "'");
	DEBUGMSG("sscan: '" << nSteps << "' _t: '" << nSteps_t << "'");
	DEBUGMSG("sscan: '" << startSteps << "' _t: '" << startSteps_t << "'");
	if (timeStep != timeStep_t || nSteps != nSteps_t ||startSteps != startSteps_t) {
		ERRMSG("Difference between sscan and >>!")
	}


	fgets(buf, 512, fptr);
	sscanf(buf, "%d%d", &shortTsFreq, &longTsFreq);
	Int     shortTsFreq_t, longTsFreq_t;
	*config_s >> shortTsFreq_t >> longTsFreq_t;
	config_s->ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	DEBUGMSG("sscan: '" << shortTsFreq << "' _t: '" << shortTsFreq_t << "'");
	DEBUGMSG("sscan: '" << longTsFreq << "' _t: '" << longTsFreq_t << "'");
	if (shortTsFreq != shortTsFreq_t || longTsFreq != longTsFreq_t) {
		ERRMSG("Difference between sscan and >>!")
	}


	fgets(buf, 512, fptr);
	sscanf(buf, "%d%d%d", &averageTsFreq, &printTsFreq, &backupTsFreq);
	Int averageTsFreq_t, printTsFreq_t, backupTsFreq_t;
	*config_s >> averageTsFreq_t >> printTsFreq_t >> backupTsFreq_t;
	config_s->ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	DEBUGMSG("sscan: '" << averageTsFreq << "' _t: '" << averageTsFreq_t << "'");
	DEBUGMSG("sscan: '" << printTsFreq << "' _t: '" << printTsFreq_t << "'");
	DEBUGMSG("sscan: '" << backupTsFreq << "' _t: '" << backupTsFreq_t << "'");
	if (averageTsFreq != averageTsFreq_t || printTsFreq != printTsFreq_t ||backupTsFreq != backupTsFreq_t) {
		ERRMSG("Difference between sscan and >>!")
	}

	fgets(buf, 512, fptr);
	sscanf(buf, "%d", &integratorType);
	Int     integratorType_t;
	*config_s >> integratorType_t;
	config_s->ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	if (integratorType != integratorType_t) {
		ERRMSG("Difference between sscan and >>!")
	}

	fgets(buf, 512, fptr);
	sscanf(buf, "%d%lf%lf%d%lf%s%lf", &v1, &d1, &d2, &v2, &couplthermo, thermotype, &couplBaros);
	set_ensemble_status(v1, d1, d2);
	if (v2 != 1)    useAtomThermo = false;
	int v;
	double d;
	*config_s >> v >> d;
	if (v1 != v || d1 != d ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> d >> v;
	if (v2 != v || d2 != d ) {
		ERRMSG("Difference between sscan and >>!")
	}
	double  couplthermo_t, couplBaros_t;
	string thermotype_t;
	*config_s >> couplthermo_t >> thermotype_t >> couplBaros_t;
	if (couplthermo != couplthermo_t ||  thermotype_t != (std::string) thermotype || couplBaros != couplBaros_t) {
		ERRMSG("Difference between sscan and >>!")
	}
	config_s->ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	fgets(buf, 512, fptr);
	sscanf(buf, "%d%lf", &v1, &rate);
	set_run_type(v1);
	*config_s >> v >> d;
	if (v1 != v || rate != d ) {
		ERRMSG("Difference between sscan and >>!")
	}
	config_s->ignore(std::numeric_limits<std::streamsize>::max(), '\n');


	fgets(buf, 512, fptr);
	sscanf(buf, "%lf%lf%lf%lf", &cutoff, &cutoffEw, &kCutoff, &cutBuff);
	*config_s >> d;
	if (d != cutoff ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> d;
	if (d != cutoffEw ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> d;
	if (d != kCutoff ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> d;
	if (d != cutBuff ) {
		ERRMSG("Difference between sscan and >>!")
	}
	config_s->ignore(std::numeric_limits<std::streamsize>::max(), '\n');


	fgets(buf, 512, fptr);
	sscanf(buf, "%lf", &density);

	*config_s >> d;
	if (d != density ) {
		ERRMSG("Difference between sscan and >>!")
	}
	config_s->ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	fgets(buf, 512, fptr);
	sscanf(buf, "%lf%lf%lf", &lxBox, &lyBox, &lzBox);
	*config_s >> d;
	if (d != lxBox ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> d;
	if (d != lyBox ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> d;
	if (d != lzBox ) {
		ERRMSG("Difference between sscan and >>!")
	}
	config_s->ignore(std::numeric_limits<std::streamsize>::max(), '\n');


	fgets(buf, 512, fptr);
	sscanf(buf, "%s%lf", str1, &cellSideLen);
	if (str1[0] == 't')  useCellPairList = true;

	*config_s >> s >> d;
	if (d != cellSideLen || s != (std::string) str1 ) {
		ERRMSG("Difference between sscan and >>!")
	}
	config_s->ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	fgets(buf, 512, fptr);
	sscanf(buf, "%d", &v1);
	*config_s >> v;
	if (v != v1 ) {
		ERRMSG("Difference between sscan and >>!")
	}
	// config_s->ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	if (v1 == 1)        // Gaussan constraint
	{
		constraintOn = true;
		sscanf(buf, "%d%lf%lf%lf", &v1, &CFB, &DFB, &TFB);
		*config_s >> v;
		if (v != v1 ) {
			ERRMSG("Difference between sscan and >>!")
		}
		*config_s >> d;
		if (d != CFB ) {
			ERRMSG("Difference between sscan and >>!")
		}
		*config_s >> d;
		if (d != DFB ) {
			ERRMSG("Difference between sscan and >>!")
		}
		*config_s >> d;
		if (d != TFB ) {
			ERRMSG("Difference between sscan and >>!")
		}
		config_s->ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
	else if (v1 == 2)
	{
		constraintOn = true;
		sscanf(buf, "%d%lf%d", &v1, &tolerance, &maxCount);
		*config_s >> v;
		if (v != v1 ) {
			ERRMSG("Difference between sscan and >>!")
		}
		*config_s >> d;
		if (d != tolerance ) {
			ERRMSG("Difference between sscan and >>!")
		}
		*config_s >> v;
		if (v != maxCount ) {
			ERRMSG("Difference between sscan and >>!")
		}
		config_s->ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
	else if (v1 == 0)
	{
		constraintOn = false;
		sscanf(buf, "%d%lf%d", &v1, &tolerance, &bendingType);
		if (bendingType == 1){
			harmonicbend = true;
		}
		else if (bendingType == 2){
			amoebabend = true;
		}
		*config_s >> v;
		if (v != v1 ) {
			ERRMSG("Difference between sscan and >>!")
		}

		*config_s >> d;
		DEBUGMSG("sscan: '" << tolerance << "' _t: '" << d << "'");
		if (d != tolerance ) {
			ERRMSG("Difference between sscan and >>!")
		}
		*config_s >> v;
		DEBUGMSG("sscan: '" << bendingType << "' _t: '" << v << "'");
		if (v != bendingType ) {
			ERRMSG("Difference between sscan and >>!")
		}
		config_s->ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}

	fgets(buf, 512, fptr);
	sscanf(buf, "%s%lf%lf", str1, &klUB, &l0UB);
	if (str1[0] == 't') computeharmoUB = true;

	*config_s >> s >> d;
	if (d != klUB || s != (std::string) str1 ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> d;
	if (d != l0UB ) {
		ERRMSG("Difference between sscan and >>!")
	}
	config_s->ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	fgets(buf, 512, fptr);
	sscanf(buf, "%s%s%lf", str1, str2, &switchDist);
	if (str1[0] == 't') lCorrectOn = true;
	if (str2[0] == 't') switchOn = true;
	*config_s >> s;
	if (s != (std::string) str1 ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> s;
	if (s != (std::string) str2 ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> d;
	if (d != switchDist ) {
		ERRMSG("Difference between sscan and >>!")
	}
	config_s->ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	fgets(buf, 512, fptr);
	sscanf(buf, "%s%s%s%s%lf%lf", str1, str2, str3, str4, &alpha, &kCutoff);
	if (str1[0] == 't'){
		computeCoulomb = true;
	}
	if (str2[0] == 't'){
		computeEwald = true;
	}
	if (str3[0] == 't'){
		computeSurfCorrect = true;
	}
	if (str4[0] == 't'){
		computeWolf = true;
	}
	*config_s >> s;
	if (s != (std::string) str1 ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> s;
	if (s != (std::string) str2 ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> s;
	if (s != (std::string) str3 ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> s;
	if (s != (std::string) str4 ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> d;
	if (d != alpha ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> d;
	if (d != kCutoff ) {
		ERRMSG("Difference between sscan and >>!")
	}
	config_s->ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	fgets(buf, 512, fptr);
	sscanf(buf, "%s%s%s%s%s%s%lf%lf", str1, str2, str3, str4, str5, str6, &cutoffInd, &adamp);
	if (str1[0] == 't') computeMultiPol  = true;
	if (str2[0] == 't') computeInduction = true;
	if (str3[0] == 't') computeTorq      = true;
	if (str4[0] == 't') dampInduction    = true;
	if (str5[0] == 't') doInductionMin   = true;
	if (str6[0] == 't') dampAmoeba       = true;

	*config_s >> s;
	if (s != (std::string) str1 ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> s;
	if (s != (std::string) str2 ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> s;
	if (s != (std::string) str3 ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> s;
	if (s != (std::string) str4 ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> s;
	if (s != (std::string) str5 ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> s;
	if (s != (std::string) str6 ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> d;
	if (d != cutoffInd ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> d;
	if (d != adamp ) {
		ERRMSG("Difference between sscan and >>!")
	}
	config_s->ignore(std::numeric_limits<std::streamsize>::max(), '\n');


	fgets(buf, 512, fptr);
	sscanf(buf, "%s%lf%lf%lf%lf%lf", str1, &nBuff, &mBuff, &deltaBuff, &gammaBuff, &Hfactor);
	if (str1[0] == 't') computeBuffvdW = true;
	*config_s >> s;
	if (s != (std::string) str1 ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> d;
	if (d != nBuff ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> d;
	if (d != mBuff ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> d;
	if (d != deltaBuff ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> d;
	if (d != gammaBuff ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> d;
	if (d != Hfactor ) {
		ERRMSG("Difference between sscan and >>!")
	}
	config_s->ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	fgets(buf, 512, fptr);
	sscanf(buf, "%s%lf%lf", str1, &nMie, &mMie);
	if (str1[0] == 't') computeMie = true;

	*config_s >> s;
	if (s != (std::string) str1 ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> d;
	if (d != nMie ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> d;
	if (d != mMie ) {
		ERRMSG("Difference between sscan and >>!")
	}
	config_s->ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	fgets(buf, 512, fptr);
	sscanf(buf, "%s%s", str1, str2);
	if (str1[0] == 't') hasExtForces = true;
	if (str2[0] == 't') hasFixedAtom = true;

	*config_s >> s;
	if (s != (std::string) str1 ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> s;
	if (s != (std::string) str2 ) {
		ERRMSG("Difference between sscan and >>!")
	}
	config_s->ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	fgets(buf, 512, fptr);
	sscanf(buf, "%d%d", &startSampling, &samplingTsFreq);
	*config_s >> v;
	if (v != startSampling ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> v;
	if (v != samplingTsFreq ) {
		ERRMSG("Difference between sscan and >>!")
	}
	config_s->ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	fgets(buf, 512, fptr);
	sscanf(buf, "%s%d%d", str1, &startTTCF, &ttcfTsFreq);
	if (str1[0] == 't') useTTCF = true;
	*config_s >> s;
	if (s != (std::string) str1 ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> v;
	if (v != startTTCF ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> v;
	if (v != ttcfTsFreq ) {
		ERRMSG("Difference between sscan and >>!")
	}
	config_s->ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	fgets(buf, 512, fptr);
	sscanf(buf, "%s%d", str1, &nBins);
	if (str1[0] == 't') computeRDF = true;
	*config_s >> s;
	if (s != (std::string) str1 ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> v;
	if (v != nBins ) {
		ERRMSG("Difference between sscan and >>!")
	}
	config_s->ignore(std::numeric_limits<std::streamsize>::max(), '\n');


	fgets(buf, 512, fptr);
	sscanf(buf, "%s%s", str1, molTrajectoryFile);
	if (str1[0] == 't') computeMSD = true;
	*config_s >> s;
	if (s != (std::string) str1 ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> s;
	if (s != (std::string) molTrajectoryFile ) {
		ERRMSG("Difference between sscan and >>!")
	}
	config_s->ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	fgets(buf, 512, fptr);
	sscanf(buf, "%s%s", str1, molVelocityFile);
	if (str1[0] == 't') computeVACF = true;
	*config_s >> s;
	if (s != (std::string) str1 ) {
		ERRMSG("Difference between sscan and >>!")
	}
	*config_s >> s;
	if (s != (std::string) molVelocityFile ) {
		ERRMSG("Difference between sscan and >>!")
	}
	config_s->ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	fgets(buf, 512, fptr);
	sscanf(buf, "%s", str1);
	if (str1[0] == 't') computeLustig = true;

	*config_s >> s;
	if (s != (std::string) str1 ) {
		ERRMSG("Difference between sscan and >>!")
	}
	config_s->ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	DEBUGMSG("close config file" );
	//    }                               // end of reading process on process 0

	//    MPI::COMM_WORLD.Bcast(sysDataFile,80,MPI::CHAR,0);
	//    MPI::COMM_WORLD.Bcast(coordinateFile,80,MPI::CHAR,0);
	//    MPI::COMM_WORLD.Bcast(resultFile,80,MPI::CHAR,0);
	//    MPI::COMM_WORLD.Bcast(restartFile,80,MPI::CHAR,0);
	//    MPI::COMM_WORLD.Bcast(dumpFile,80,MPI::CHAR,0);
	//    MPI::COMM_WORLD.Bcast(coordinateFile,80,MPI::CHAR,0);
	//    MPI::COMM_WORLD.Bcast(coordinateFile,80,MPI::CHAR,0);
	//    MPI::COMM_WORLD.Bcast(coordinateFile,80,MPI::CHAR,0);
	//    MPI::COMM_WORLD.Barrier();                         // good
	//    if(rank != 0) {cout<<"the data stored in the file ...... "<<sysDataFile<<endl;}
	fclose(fptr);
}

void SimConfiguration::set_box(Double lx, Double ly, Double lz)
{
	lxBox=lx;   lyBox=ly;   lzBox=lz;
}

void SimConfiguration::set_cell_side(Double len)
{
	cellSideLen=len;
}

// set_cutoff() will be called when need to evaluate ewald force
// first we estimate a cutoff for ewald real space and then the larger cutoff is chosen
// for both LJ force and ewald real space
void SimConfiguration::set_cutoff(Int nAtoms, Double volume)
{
	Double tRate = 5.5;         // default value as established by Moldy
	Double pVal = -log(accuracy);
	pVal = sqrt(pVal);
	alpha = SQRT_PI*pow(tRate*nAtoms/(volume*volume), 1.0/6.0);
	Double realCut = pVal/alpha;

	cutoff = MAX(cutoff, realCut);        // cutOff should be adequate for both LJ force & ewald real space
	alpha = pVal/cutoff;                  // re-correct alpha based on chosed cutOff
	kCutoff = 2.0*alpha*pVal;
}

void SimConfiguration::set_density(Double d)
{
	density = d;
}

void SimConfiguration::set_temperature(Double t)
{
	temperature = t;
}

void SimConfiguration::set_density(Int nAtoms)
{
	Double vol;
	vol = lxBox*lyBox*lzBox;
	if (vol > 0.0) density = nAtoms/vol;
	else density = 0.0;
}

void SimConfiguration::set_temperature(Int nAtoms, Double kinEnergy)
{
	if (nAtoms > 0) temperature = kinEnergy/(1.5*nAtoms);
	else temperature = 0.0;
}

void SimConfiguration::set_elapsedTime(Double currentSteps)
{
	elapsedTime = currentSteps*timeStep;
}

void SimConfiguration::set_ensemble_status(Int e, Double v1, Double v2)
{
	ensembleStatus = (EnsembleStatus)e;
	switch(ensembleStatus)
	{
		case NVT:
			constantTemperature = true;
			constantVolume = true;
			temperature = v1;
			volume = v2;
			break;
		case NPT:
			constantTemperature = true;
			constantPressure = true;
			temperature = v1;
			pressure = v2;
			break;
		case NVE:
			constantVolume = true;
			constantEnergy = true;
			energy = v1;
			volume = v2;
			break;
		case NPE:
			constantPressure = true;
			constantEnergy = true;
			energy = v1;
			pressure = v2;
			break;
	}
}

void SimConfiguration::set_run_type(Int t)
{
	runType =(RunType) t;
	switch(runType)
	{
		case EQUIL:
			doEquilibrium = true;
			break;
		case COMPR:
			doCompression = true;
			break;
		case SHEAR:
			doShear = true;
			break;
		case ELONG:
			doElongation = true;
			break;
	}
}

void SimConfiguration::write_config(ofstream &ofp)
{
#ifdef DEBUG
	DEBUGMSG("write config to file");
#endif

	char* tab = "     ";
	ofp << sysDataFile << tab << coordinateFile << tab <<  "// input file names" << endl;
	ofp << resultFile << tab << restartFile << tab << trajectoryFile << tab << dumpFile << endl;
	ofp << timeStep << tab << nSteps << tab << startSteps << tab << "//timeStep  nSteps  startSteps"<< endl;
	ofp << shortTsFreq << tab << longTsFreq << tab << "//nonbondTsFreq, longElectTsFreq" << endl;
	ofp << averageTsFreq << tab << printTsFreq << tab << backupTsFreq << tab;
	ofp << "// computeAvgTsFreq  printTsFreq  backupTsFreq" << endl;
	ofp << integratorType << tab << "//integrator type: 0 - Gear; 1 - Velocity Verlet" << endl;
	switch(ensembleStatus)
	{
		case NVT:
			ofp << ensembleStatus << tab << temperature << tab << volume << tab;
			break;
		case NPT:
			ofp << ensembleStatus << tab << temperature << tab << pressure << tab;
			break;
		case NVE:
			break;
		case NPE:
			break;
	}
	ofp << couplthermo << tab << thermotype << "-thermostat" << tab << couplBaros << tab;
	ofp << "//ensembleStatus(0-NVT, 1-NPT)  Temperature volume/pressure" << endl;
	if (useAtomThermo)  ofp << "use Atom thermostat" << endl;
	else ofp << "use molecular thermostat" << endl;

	ofp << runType << tab << rate << tab << "//runType(0-Equil, 1-Comps, 2-Shear, 3-Elong)  rate" << endl;
	ofp << cutoff << tab << cutoffEw << tab << kCutoff << tab << cutBuff<< tab;
	ofp << "//cutoff, cutoffEw,  kCutoff, cutBuff" << endl;
	ofp <<  density << tab << "// density" << endl;
	ofp << lxBox << tab << lyBox << tab << lzBox << tab << "// lxBox  lyBox  lzBox" << endl;
	if (useCellPairList)
		ofp << tab << tab  << cellSideLen <<  tab <<"// useCellPairList   cellSideLen" << endl;
	else
		ofp <<  "f" << tab << cellSideLen << tab << "// useCellPairList   cellSideLen" << endl;
	if (constraintOn) { ofp << "t" << tab;
		ofp << CFB << tab << DFB << tab << TFB << tab << "//constraintOn, CFB, DFB, TFB" << endl;}
	else { ofp << "f" << tab << tolerance << tab << bendingType << tab;
		if (harmonicbend) {
			ofp << " Harmonic_pot" << tab ;
		}
		else if (amoebabend) {
			ofp << " AMOEBA " << tab ;
		}
		else {
			DEBUGMSG("ACHTUNG: are you ussing a flexible molecule??,");
			ERRORMSG("if so, please assing a bending potential in config.txt: 1:Harmonic, 2:Amoeba.");
		}
		ofp << "//constraintOn, tol, 1:Harmonic 2:amoeba, bool" << endl;}

	if (computeharmoUB) ofp << "t" << tab;
	else ofp << "f" << tab;
	ofp <<  klUB << tab << l0UB <<  "// computeharmoUB;  klUB;  l0UB" << endl;

	if (lCorrectOn)  ofp << "t" << tab;
	else  ofp << "f" << tab;
	if(switchOn)
		ofp <<  "t" << tab << switchDist << tab ;
	else
		ofp << "f" << tab << switchDist << tab;
	ofp << "// longRangeCorrect; shiftOn; shiftingDist" << endl;
	if (computeCoulomb) ofp << "t" << tab;
	else ofp << "f" << tab;
	if (computeEwald) ofp << "t" << tab << accuracy << tab;
	else  ofp << "f" << tab << accuracy << tab;
	if (computeWolf) ofp << "t" << tab;
	else ofp << "f" << tab;
	ofp << "// computeCoulomb;  computeEwald; accuracy; computeWolf" << endl;
	if (computeMultiPol) ofp << "t" << tab;
	else ofp << "f" << tab;
	if (computeInduction) ofp << "t" << tab;
	else ofp << "f" << tab;
	if (computeTorq) ofp << "t" << tab;
	else ofp << "f" << tab;
	if (dampInduction) ofp << "t" << tab;
	else ofp << "f" << tab;
	if (doInductionMin) ofp << "t" << tab;
	else ofp << "f" << tab;
	if (dampAmoeba) ofp << "t" << tab;
	else ofp << "f" << tab;
	ofp << cutoffInd << tab;
	ofp << adamp << tab;
	ofp << "// computeMultipoles;  computeInduction;  computeTorques;  dampInduction; doInductionMin; dampAmoeba; cutoffInd; a_damping " << endl;
	if (computeBuffvdW) ofp << "t" << tab;
	else ofp << "f" << tab;
	ofp << nBuff << tab << mBuff << tab << deltaBuff << tab << gammaBuff << tab << Hfactor << tab << "// computeBuffvdW;  nBuff;  mBuff;  deltaBuff;  gammaBuff;  Hfactor;  " << endl;
	if (hasExtForces) ofp << "t" << tab;
	else  ofp << "f" << tab;
	if (hasFixedAtom) ofp << "t" << tab;
	else  ofp << "f" << tab;
	ofp << "// hasExtForces   hasFixedAtom" << endl;
	ofp << startSampling << tab << samplingTsFreq << tab;
	ofp << "// start sampling steps;  samplingTsFreq" << endl;
	if (useTTCF)
		ofp << "t" << tab << startTTCF << tab << ttcfTsFreq << tab << "// useTTCF  startTTCF   ttcfTsFreq"<< endl;
	else
		ofp << "f" << tab << startTTCF << tab << ttcfTsFreq << tab << "// useTTCF  startTTCF   ttcfTsFreq"<< endl;
	if (computeVACF) ofp << "t" << endl;
	else  ofp << "f" << endl;
	if (computeRDF)
		ofp << "t" << tab << nBins << tab;
	else
		ofp << "f" << tab << nBins << tab;
	ofp << "// computeRDF  numOfBins" << endl;
	if (computeMSD) ofp << "t          // computeMSD" << endl;
	else            ofp << "f          // computeMSD" << endl;

	ofp << endl;
}
