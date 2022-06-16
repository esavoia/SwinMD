/** Local2Global.cpp -- 
 **
 ** Copyright (C) 2003
 ** Centre for Molecular Simulation (CMS)
 ** School of Information Technology
 ** Swinburne University of Technology
 ** PO Box 21 Hawthorn, Vic 3122, Australia
 **
 ** Local to Global coordinante transformation
 ** Author: eoyarzua
 ** Email: eoyarzua@swin.edu.au
 **/

#include "Local2Global.h"
Local2Global::Local2Global(Ensemble* ensemble)
{
    DEBUGMSG("Creating Local2Global");

    myEnsemble = ensemble;
    // myParams = myEnsemble->myParams;
    atoms = myEnsemble->atoms;
    numAtoms = myEnsemble->nAtoms;
    numMols = myEnsemble->nMols;
    myMols = myEnsemble->molecules;

    pp1  = NULL;
    pp2  = NULL;
    uu   = NULL;
    vv   = NULL;
    ww   = NULL;
    dip  = NULL;
    Qxx  = NULL;
    Qxy  = NULL;
    Qxz  = NULL;
    Qyy  = NULL;
    Qyz  = NULL;
    Qzz  = NULL;

    pp1 = new Vector3* [numMols];
    pp2 = new Vector3* [numMols];
    uu  = new Vector3* [numMols];
    vv  = new Vector3* [numMols];
    ww  = new Vector3* [numMols];
    dip  = new Vector3[numAtoms];
    Qxx  = new double[numAtoms];
    Qxy  = new double[numAtoms];
    Qxz  = new double[numAtoms];
    Qyy  = new double[numAtoms];
    Qyz  = new double[numAtoms];
    Qzz  = new double[numAtoms];

    if ((pp1==NULL)||(pp2==NULL)||(uu==NULL)||(vv==NULL)||(ww==NULL)||(dip==NULL))
        ERRORMSG("memory allocation error in Local2Global");

    if ((Qxx==NULL)||(Qxy==NULL)||(Qxz==NULL)||(Qyy==NULL)||(Qyz==NULL)||(Qzz==NULL))
        ERRORMSG("memory allocation error 2 in Local2Global");

    for (int nm = 0; nm < numMols; nm++)
        {
            int ns = myMols[nm].numAtoms;
            if (ns <= 1)
                continue;
            pp1[nm] = new Vector3 [ns];
            pp2[nm] = new Vector3 [ns];
            uu[nm]  = new Vector3 [ns];
            vv[nm]  = new Vector3 [ns];
            ww[nm]  = new Vector3 [ns];
            if((pp1[nm]==NULL)||(pp2[nm]==NULL)||(uu[nm]==NULL)||(vv[nm]==NULL)||(ww[nm]==NULL))
                ERRORMSG("memory allocation error for Multipole Ewald");
        }
}

Local2Global::~Local2Global()
{
    if(Qxx != NULL) delete Qxx;
    if(Qxy != NULL) delete Qxy;
    if(Qxz != NULL) delete Qxz;
    if(Qyy != NULL) delete Qyy;
    if(Qyz != NULL) delete Qyz;
    if(Qzz != NULL) delete Qzz;
    if(uu != NULL) delete uu;
    if(vv != NULL) delete vv;
    if(ww != NULL) delete ww;
    if(pp1 != NULL) delete pp1;
    if(pp2 != NULL) delete pp2;
    if(dip != NULL) delete dip;
}


void Local2Global::FrameofReference()
{
    Molecule *mol;
    Vector3 oneX(1,0,0);
    Vector3 oneY(0,1,0);
    double val;
    double dipxi, dipyi, dipzi;
    double Qxxi, Qxyi, Qxzi, Qyyi, Qyzi, Qzzi;
    int idi;

    //-------------------------------------------------------------------------------
    //  First Let's build the Reference axes p1 and p2 for each site 
    //  at each water molecule
    //-------------------------------------------------------------------------------

    for (int nm = 0; nm < numMols; nm++)
        {
            int ns = myMols[nm].numAtoms;
            mol = &myMols[nm];

            Atom *atomH1 = mol->myAtoms[0];
            Atom *atomO  = mol->myAtoms[1];
            Atom *atomH2 = mol->myAtoms[2];

            // Local Axes of Left Hydrogen1
            pp1[nm][0] = atomO->position - atomH1->position;
            val        = pp1[nm][0].x / pp1[nm][0].length();
            pp2[nm][0] = oneX;
            if ( fabs(val) > 0.8660) {
                pp2[nm][0] = oneY;
            }

            // Lolal Axes of Central Oxygen
	    pp1[nm][1] = atomH1->position - atomO->position;
            pp2[nm][1] = atomH2->position - atomO->position;

            pp1[nm][1] /= pp1[nm][1].length();
            pp2[nm][1] /= pp2[nm][1].length();

            pp1[nm][1] = pp1[nm][1] + pp2[nm][1];

            // Local Axes of Right Hydrogen2
	    pp1[nm][2] = atomO->position - atomH2->position;
            val        = pp1[nm][2].x / pp1[nm][2].length();
            pp2[nm][2] = oneX;
            if ( fabs(val) > 0.8660) {
                pp2[nm][2] = oneY;
            }
	}

    //-------------------------------------------------------------------------------
    // Next Let's build the Local coordinate system u,v,w
    // for each water molecule at each site
    //-------------------------------------------------------------------------------
    
    Vector3 ss;

    for (int nm = 0; nm < numMols; nm++)
        {
            int ns = myMols[nm].numAtoms;
            mol = &myMols[nm];
            for (int i = 0; i < ns; i++)
                {
                    pp1[nm][i] /= pp1[nm][i].length();

                    ww[nm][i] = pp1[nm][i];
                    ss        = pp2[nm][i] - ( ( pp2[nm][i] * ww[nm][i] ) * ww[nm][i] );
                    uu[nm][i] = ss / ss.length();
                    vv[nm][i] = cross( ww[nm][i], uu[nm][i] );
                }
        }

    //-------------------------------------------------------------------------------
    // Now we transform the Multipoles in the Local Frame
    // of reference to the Global Frame of reference.
    //-------------------------------------------------------------------------------
    for (int nm = 0; nm < numMols; nm++)
        {
            int ns = myMols[nm].numAtoms;
            mol = &myMols[nm];
            for (int i = 0; i < ns; i++)
                {
		    idi = myEnsemble->molecules[nm].myAtoms[i]->atomID;

		    dipxi = uu[nm][i].x*atoms[idi].dipx + vv[nm][i].x*atoms[idi].dipy + ww[nm][i].x*atoms[idi].dipz;
            	    dipyi = uu[nm][i].y*atoms[idi].dipx + vv[nm][i].y*atoms[idi].dipy + ww[nm][i].y*atoms[idi].dipz;
            	    dipzi = uu[nm][i].z*atoms[idi].dipx + vv[nm][i].z*atoms[idi].dipy + ww[nm][i].z*atoms[idi].dipz;

		    dip[idi].x = dipxi;
		    dip[idi].y = dipyi;
		    dip[idi].z = dipzi;

		    Qxxi  = uu[nm][i].x*( uu[nm][i].x*atoms[idi].quadxx + vv[nm][i].x*atoms[idi].quadxy + ww[nm][i].x*atoms[idi].quadxz )  \
                          + vv[nm][i].x*( uu[nm][i].x*atoms[idi].quadxy + vv[nm][i].x*atoms[idi].quadyy + ww[nm][i].x*atoms[idi].quadyz )  \
                          + ww[nm][i].x*( uu[nm][i].x*atoms[idi].quadxz + vv[nm][i].x*atoms[idi].quadyz + ww[nm][i].x*atoms[idi].quadzz )  ;

            	    Qxyi  = uu[nm][i].x*( uu[nm][i].y*atoms[idi].quadxx + vv[nm][i].y*atoms[idi].quadxy + ww[nm][i].y*atoms[idi].quadxz )  \
                          + vv[nm][i].x*( uu[nm][i].y*atoms[idi].quadxy + vv[nm][i].y*atoms[idi].quadyy + ww[nm][i].y*atoms[idi].quadyz )  \
                          + ww[nm][i].x*( uu[nm][i].y*atoms[idi].quadxz + vv[nm][i].y*atoms[idi].quadyz + ww[nm][i].y*atoms[idi].quadzz )  ;

            	    Qxzi  = uu[nm][i].x*( uu[nm][i].z*atoms[idi].quadxx + vv[nm][i].z*atoms[idi].quadxy + ww[nm][i].z*atoms[idi].quadxz )  \
                          + vv[nm][i].x*( uu[nm][i].z*atoms[idi].quadxy + vv[nm][i].z*atoms[idi].quadyy + ww[nm][i].z*atoms[idi].quadyz )  \
                          + ww[nm][i].x*( uu[nm][i].z*atoms[idi].quadxz + vv[nm][i].z*atoms[idi].quadyz + ww[nm][i].z*atoms[idi].quadzz )  ;

            	    Qyyi  = uu[nm][i].y*( uu[nm][i].y*atoms[idi].quadxx + vv[nm][i].y*atoms[idi].quadxy + ww[nm][i].y*atoms[idi].quadxz )  \
                          + vv[nm][i].y*( uu[nm][i].y*atoms[idi].quadxy + vv[nm][i].y*atoms[idi].quadyy + ww[nm][i].y*atoms[idi].quadyz )  \
                          + ww[nm][i].y*( uu[nm][i].y*atoms[idi].quadxz + vv[nm][i].y*atoms[idi].quadyz + ww[nm][i].y*atoms[idi].quadzz )  ;

                    Qyzi  = uu[nm][i].y*( uu[nm][i].z*atoms[idi].quadxx + vv[nm][i].z*atoms[idi].quadxy + ww[nm][i].z*atoms[idi].quadxz )  \
                          + vv[nm][i].y*( uu[nm][i].z*atoms[idi].quadxy + vv[nm][i].z*atoms[idi].quadyy + ww[nm][i].z*atoms[idi].quadyz )  \
                          + ww[nm][i].y*( uu[nm][i].z*atoms[idi].quadxz + vv[nm][i].z*atoms[idi].quadyz + ww[nm][i].z*atoms[idi].quadzz )  ;

            	    Qzzi  = uu[nm][i].z*( uu[nm][i].z*atoms[idi].quadxx + vv[nm][i].z*atoms[idi].quadxy + ww[nm][i].z*atoms[idi].quadxz )  \
                          + vv[nm][i].z*( uu[nm][i].z*atoms[idi].quadxy + vv[nm][i].z*atoms[idi].quadyy + ww[nm][i].z*atoms[idi].quadyz )  \
                          + ww[nm][i].z*( uu[nm][i].z*atoms[idi].quadxz + vv[nm][i].z*atoms[idi].quadyz + ww[nm][i].z*atoms[idi].quadzz )  ;

            	    //----------------------------------------------------------------------
            	    //  Traceless Quadrupoles into traced formulae 
            	    //----------------------------------------------------------------------

		    Qxxi  = Qxxi/3;
                    Qxyi  = Qxyi/3;
                    Qxzi  = Qxzi/3;
                    Qyyi  = Qyyi/3;
                    Qyzi  = Qyzi/3;
                    Qzzi  = Qzzi/3;

		    Qxx[idi] = Qxxi;
		    Qxy[idi] = Qxyi;
		    Qxz[idi] = Qxzi;
		    Qyy[idi] = Qyyi;
		    Qyz[idi] = Qyzi;
		    Qzz[idi] = Qzzi;
		}
	}
}
