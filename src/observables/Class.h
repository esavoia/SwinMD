/** MSD.h -- 
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

#ifndef MSD_H
#define MSD_H

#include "Ensemble.h"
#include "utils/Vector3.h"

/**
 ** Forward declarations
 **/


class MSD {
    /**
     ** Data Member - 
     **/
    Double *msdSum;         // acumulated msd for every sampling time
    Int    counter;
    Int    interval;        // interval number of timesteps for sampling 
    Int    outInterval;     // interval number of timesteps for output 
    Double halfLx, halfLy, halfLz;

    Vector3 *R0;
    Vector3 *Rt;

    /**
     ** constructor and destructor
     **/
    public:
    MSD();
    ~MSD();


    /**
     ** methods
     **/
    public:
    void set_R0();
    void sampling();

    inline void undo_pbc(Vector3 rt, Vector3 ref)
    {
        Vector3 diff = rt - ref;

        if (diff.x > halfLx)        rt.x -= halfLx;
        else if (diff.x < halfLx)   rt.x += halfLx;
        if (diff.y > halfLy)        rt.y -= halfLy;
        else if (diff.y < halfLy)   rt.y += halfLy;
        if (diff.z > halfLz)        rt.z -= halfLz;
        else if (diff.z < halfLz)   rt.z += halfLz;
    }


};

#endif
