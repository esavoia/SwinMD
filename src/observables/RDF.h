/** RDF.h -- 
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

#ifndef RDF_H
#define RDF_H

#include "../NEMD_defs.h"
#include "../Ensemble.h"
#include "../utils/Vector3.h"

/**
 ** Forward declarations
 **/


class RDF {
    /**
     ** Data Member - 
     **/
    Double  ***rdfTable;           
    Int     nBins;              // number of histogram bins
    Int     interval;           // interval number of timesteps for sampling 
    Int     outInterval;        // interval number of timesteps for output 
    bool    initialised;

    Int     nAtoms, nTypes;    
    Double  halfLx, halfLy, halfLz;
    Double  maxR, maxR2, deltaR, deltaRr;               // 

    Int     nSamples;           // number of sampling or measurements
    Int     counter;            // number of sampling completed
    Int     *numPerType;         // number of atoms or mols for each type of species
    Atom    *atoms;

    Ensemble *myEnsemble;

    /**
     ** constructor and destructor
     **/
    public:
    RDF(Ensemble *ensemble, Int nSampling);
    ~RDF();


    /**
     ** methods
     **/
    public:
    void init();
    void sampling();
    void clear();
    void write_result();
    void set_bin_num(Int n);

};

#endif
