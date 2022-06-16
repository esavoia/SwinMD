/** RDF.cpp -- 
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

#include "RDF.h"

RDF :: RDF(Ensemble *ensemble, Int nSampling)
{
    myEnsemble = ensemble;
    nAtoms = myEnsemble->nAtoms;
    nTypes = myEnsemble->nAtomTypes;
    atoms  = myEnsemble->atoms;
    halfLx = myEnsemble->halfLx;
    halfLy = myEnsemble->halfLy;
    halfLz = myEnsemble->halfLz;
    nBins = myEnsemble->myConfig->get_rdf_nbins();
    counter = 0;
    interval = 1;
    outInterval = 1;    
    nSamples = nSampling;
    initialised = false;
    rdfTable = NULL;
    numPerType = NULL;

    init();
}

RDF :: ~RDF()
{
    clear();
}

void RDF::init()
{
    clear();
    maxR = halfLx;
    maxR = MIN(maxR, halfLy);
    maxR = MIN(maxR, halfLz);
    maxR2 = maxR*maxR;
    deltaR = maxR/nBins;
    deltaRr = nBins/maxR;

    if(maxR <= 0)
        ERRORMSG("rdf found boxdimension <= 0");

    rdfTable = new Double** [nTypes];
    if (rdfTable == NULL)
        ERRORMSG("could not allocate rdfTable");
    for (Int i = 0; i < nTypes; i++)
    {
        rdfTable[i] = new Double* [nTypes];
        if (rdfTable[i] == NULL)
            ERRORMSG("could not allocate rdfTable");
    }
    for (Int i = 0; i < nTypes; i++)
    {
        for (Int j = 0; j < nTypes; j++)
        {
            rdfTable[i][j] = new Double[nBins];
            if (rdfTable[i][j] == NULL)
                ERRORMSG("could not allocate rdfTable");
        }
    }
    numPerType = new Int[nTypes];
    if (numPerType == NULL)
        ERRORMSG("Fail to allocate numPerType");

    for (Int i = 0; i < nTypes; i++)
        for (Int j = 0; j < nTypes; j++)
            for (Int k = 0; k < nBins; k++)
                rdfTable[i][j][k] = 0;
    for (Int i = 0; i < nTypes; i++)
        numPerType[i] = 0;

    for (Int i = 0; i < nAtoms; i++)
    {
        Int type = atoms[i].atomType;
        if (type >= nTypes)
            ERRORMSG("undefined atom type encountered in RDF");
        numPerType[type] ++;
    }

    //check number
    Int num = 0;
    for (Int i = 0; i < nTypes; i++)
        num += numPerType[i];
    if (num != nAtoms)
        ERRORMSG("wrong number of atoms");

    initialised = true;
}

void RDF::clear()
{
    if (rdfTable != NULL)       delete [] rdfTable;
    if (numPerType != NULL)     delete [] numPerType;
    rdfTable = NULL;
    numPerType = NULL;
}

void RDF :: sampling()
{
    Atom *atomi; //, *atomj; // never used and commented by Jianhui Li
    Int type1, type2;

    if(!initialised)    init();
    
    for (Int i = 0; i < nAtoms; i++)
    {
        atomi = &myEnsemble->atoms[i];
        type1 = atomi->atomType;
        for (Int j = i+1; j < nAtoms; j++)
        {
            if (myEnsemble->exclusion_check(atomi->atomID, myEnsemble->atoms[j].atomID))
                continue;
            type2 = myEnsemble->atoms[j].atomType;
            Vector3 diff = atomi->position - myEnsemble->atoms[j].position;
            myEnsemble->apply_pbc(diff);
            Double r = diff.length();
            Double r2 = r*r;
            if (r2 < maxR2)
            {
                Int binid = (Int) (deltaRr*r);
                rdfTable[type1][type2][binid] ++;
                rdfTable[type2][type1][binid] ++;
            }
        }                
    } 
    counter++;
//    cout<<"counter = llllllllllllllllllllllllllllllllllllllllllllll "<<counter<<"  nSamples =  "<<nSamples <<endl;
    if (counter >= nSamples)
    {
//        cout<<"counter = "<<counter<<"   nSamples 3273072390572390523759023573290124712907=  "<<nSamples <<endl;   
        write_result();
    }
}

void RDF :: write_result()
{
    Double pi4dr;
    ofstream of("RDF.out", ios::out);
    if(myEnsemble->volume <= 0)
        ERRORMSG("zero or negtive volume");

    pi4dr = 4.0*PI*deltaR;
    pi4dr /= myEnsemble->volume;

    of << "#RDF output after " << counter*interval << " timesteps" << endl;
    of << "Types: " << '\t';
    for (Int t1 = 0; t1 < nTypes; t1++)
    {
        for (Int t2 = t1; t2 < nTypes; t2++)
            of << t1 << '-' << t2 << '\t';
    }
    of << endl;
//    cout<<"nBins ()()()()()()()()()()()()() =  "<<nBins<<endl;
    for (Int i = 0; i < nBins; i++)
    {
        Double r = i*deltaR + 0.5*deltaR;
        Double r2 = r*r;
        of <<fixed<<setw(14)<< r;
        for (Int t1 = 0; t1 < nTypes; t1++)
        {
            for (Int t2 = t1; t2 < nTypes; t2++)
            {
                Double rdf = (rdfTable[t1][t2][i] + rdfTable[t2][t1][i])/pi4dr;
                rdf = 0.5*rdf/(r2*numPerType[t1]*numPerType[t2]*counter);
                of <<setw(14)<< rdf;
            }
        }
        of << endl;
    }
}
