/** ConstraintForce.cpp -- 
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

#include "ConstraintForce.h"
#include <iostream>
#include <fstream>

ConstraintForce :: ConstraintForce(Ensemble *ensemble, Double cfb, Double dfb) : Force(ensemble)
{
    nAtoms = myEnsemble->nAtoms;
    molecules = myEnsemble->molecules;
    nMolecules = myEnsemble->nMols;
    params = myEnsemble->myParams;

    matrixM = NULL;
    matrixL = NULL;
    matrixC = NULL;
    bondLength = NULL;
    rm = NULL;
    vm = NULL;
    fm = NULL;
    velfb = NULL;
    momfb = NULL;
    index = NULL;
    
    CFB = cfb;
    DFB = dfb;

    init();
    set_matrix();

    write_force_info(*sysdataFile);
}

ConstraintForce :: ~ConstraintForce()
{
    if (matrixM != NULL) delete [] matrixM;
    if (matrixL != NULL) delete [] matrixL;
    if (matrixC != NULL) delete [] matrixC;
    if (bondLength != NULL) delete [] bondLength;
    if (vecB != NULL) delete [] vecB;
    if (rm != NULL) delete [] rm;
    if (vm != NULL) delete [] vm;
    if (fm != NULL) delete [] fm;
    if (velfb != NULL) delete [] velfb;
    if (momfb != NULL) delete [] momfb;
    if (index != NULL) delete [] index;
}

// allocate memory and zero matrix
void ConstraintForce :: init()
{
    Int ns;

    #ifdef DEBUG
        DEBUGMSG("Initialise ConstraintForce");
    #endif

    matrixM = new Int**[nMolecules];
    matrixL = new Double**[nMolecules];
    matrixC = new Double**[nMolecules];
    bondLength = new Double*[nMolecules];
    vecB = new Double*[nMolecules];
    rm = new Vector3* [nMolecules];
    vm = new Vector3* [nMolecules];
    fm = new Vector3* [nMolecules];
    index = new Int [nMolecules];
    velfb = new Vector3 [nAtoms];
    momfb = new Vector3 [nAtoms];

    if((matrixM==NULL)||(matrixL==NULL)||(matrixC==NULL)||(bondLength==NULL)    \
        ||(vecB==NULL)||(rm==NULL)||(vm==NULL)||(fm==NULL)||(velfb==NULL)||(momfb==NULL))
        ERRORMSG("fail to allocate memory for ConstraintForce");
       

    for (int nm = 0; nm < nMolecules; nm++)
    {
        ns = molecules[nm].numAtoms;
        if (ns <= 1)
            continue;

        matrixM[nm] = new Int*[ns];
        matrixL[nm] = new Double*[2*ns - 3];
        matrixC[nm] = new Double*[2*ns - 3];
        bondLength[nm] = new Double[2*ns - 3];
        vecB[nm] = new Double[2*ns - 3];
        rm[nm] = new Vector3 [2*ns - 3];
        vm[nm] = new Vector3 [2*ns - 3];
        fm[nm] = new Vector3 [2*ns - 3];

        if((matrixM[nm]==NULL)||(matrixL[nm]==NULL)||(matrixC[nm]==NULL)||(bondLength[nm]==NULL)    \
           ||(vecB[nm]==NULL)||(rm[nm]==NULL)||(vm[nm]==NULL)||(fm[nm]==NULL))
            ERRORMSG("fail to allocate memory for ConstraintForce");

        for (int j = 0; j < ns; j++)
        {
            matrixM[nm][j] = new Int[2*ns - 3];
            if(matrixM[nm][j]==NULL)
                ERRORMSG("fail to allocate memory for ConstraintForce");
            for(int k = 0; k < (2*ns - 3); k++)
                matrixM[nm][j][k] = 0;
        }
        for (int j = 0; j < (2*ns - 3); j++)
        {
            matrixL[nm][j] = new Double[2*ns - 3];
            if(matrixL[nm][j]==NULL)
                ERRORMSG("fail to allocate memory for ConstraintForce");
            for(int k = 0; k < (2*ns - 3); k++)            
                matrixL[nm][j][k] = 0.0;
        }
        for (int j = 0; j < (2*ns - 3); j++)
        {
            matrixC[nm][j] = new Double[2*ns - 3];
            if(matrixM[nm][j]==NULL)
                ERRORMSG("fail to allocate memory for ConstraintForce");
            for(int k = 0; k < (2*ns - 3); k++)            
                matrixC[nm][j][k] = 0.0;
        }
    }
}
    
void ConstraintForce :: set_matrix()
{
    Int ns;
    Int a0, a1, a2;
    Double m0, m1, m2;
    Double l1, l2, angle;

    #ifdef DEBUG
        DEBUGMSG("Setting Constraint matrix");
    #endif

    // set the selector matrix
    for (int nm = 0; nm < nMolecules; nm++)
    {
        ns = molecules[nm].numAtoms;
        if (ns <= 1)
            continue;
        int site = 0;
        for (int col = 0; col < (2*ns - 3); col += 2)
        {
            matrixM[nm][site++][col] = -1;
            matrixM[nm][site][col] = 1;
        }
        site = 0;
        for (int col = 1; col <(2*ns - 3); col += 2)
        {
            matrixM[nm][site++][col] = -1;
            matrixM[nm][site+1][col] = 1;
        }
    }

    // set the  matrixL
    for (int nm = 0; nm < nMolecules; nm++)
    {
        ns = molecules[nm].numAtoms;
        if (ns <= 1)
            continue;
        for (int col = 0; col < (2*ns - 3); col++)
        {
            for(int site = 0; site < (ns-2); site++)
            {   
                m0 = molecules[nm].myAtoms[site]->mass; 
                m1 = molecules[nm].myAtoms[site+1]->mass; 
                m2 = molecules[nm].myAtoms[site+2]->mass; 

                // bond for i molecule between 'site' and 'site+1'
                matrixL[nm][2*site][col] = matrixM[nm][site+1][col]/m1 - matrixM[nm][site][col]/m0;
                // angle for i molecule among 'site', 'site+1' and 'site+2' constraint as a bond
                // between 'site' and 'site+2'
                matrixL[nm][2*site+1][col] = matrixM[nm][site+2][col]/m2 - matrixM[nm][site][col]/m0;
            }
            // the last bond
            m1 = molecules[nm].myAtoms[ns - 2]->mass;
            m2 = molecules[nm].myAtoms[ns - 1]->mass;
            matrixL[nm][2*ns - 4][col] = matrixM[nm][ns-1][col]/m2 - matrixM[nm][ns-2][col]/m1;
        }
    }

    // set bond length
    for (int nm = 0; nm < nMolecules; nm++)
    {
        ns = molecules[nm].numAtoms;
        if (ns <= 1)
            continue;
        for(int site = 0; site < (ns-2); site++)
        {
            a0 = molecules[nm].myAtoms[site]->atomID;
            a1 = molecules[nm].myAtoms[site+1]->atomID;
            a2 = molecules[nm].myAtoms[site+2]->atomID;
            l1 = params->get_bond_length(a0, a1);
            l2 = params->get_bond_length(a1, a2);
            angle = params->get_angle(a0, a1, a2);
            bondLength[nm][2*site] = l1*l1;
            bondLength[nm][2*site+1] = l1*l1 + l2*l2 - 2*l1*l2*cos(angle);
        }
        // the last bond
        a0 = molecules[nm].myAtoms[ns-2]->atomID;
        a1 = molecules[nm].myAtoms[ns-1]->atomID;
        l1 = params->get_bond_length(a0, a1);
        bondLength[nm][2*ns-4] = l1*l1;
    }


}

void ConstraintForce :: compute()
{
    Double val, vcfb, mcfb;

    energy = 0.0;
    for (Int i = XX; i <= ZZ; i++)
        virial[i] = 0.0;

    // compute rm, vm, fm & vecB
    for (int nm = 0; nm < nMolecules; nm++)
    {
        int ns = molecules[nm].numAtoms;
        if (ns <= 1)
            continue;
        Atom** atoms = molecules[nm].myAtoms;

        for (int i = 0; i < (ns - 2); i++)
        {
            rm[nm][2*i] = atoms[i + 1]->position - atoms[i]->position;
            rm[nm][2*i + 1] = atoms[i + 2]->position - atoms[i]->position;
            vm[nm][2*i] = atoms[i + 1]->velocity - atoms[i]->velocity;
            vm[nm][2*i + 1] = atoms[i + 2]->velocity - atoms[i]->velocity;
            fm[nm][2*i] = atoms[i + 1]->force/atoms[i + 1]->mass - atoms[i]->force/atoms[i]->mass;
            fm[nm][2*i + 1] = atoms[i + 2]->force/atoms[i + 2]->mass - atoms[i]->force/atoms[i]->mass;
        }
        rm[nm][2*ns - 4] = atoms[ns - 1]->position - atoms[ns - 2]->position;
        vm[nm][2*ns - 4] = atoms[ns - 1]->velocity - atoms[ns - 2]->velocity;
        fm[nm][2*ns - 4] = atoms[ns - 1]->force/atoms[ns - 1]->mass - atoms[ns - 2]->force/atoms[ns - 2]->mass;

        for (int i = 0; i < (2*ns - 3); i++)
        {
            vecB[nm][i] = -(fm[nm][i]*rm[nm][i] + vm[nm][i].length2());
        }
    }
    // compute the coefficient matrix
    for (int nm = 0; nm < nMolecules; nm++)
    {
        int ns = molecules[nm].numAtoms;
        if (ns <= 1)
            continue;
        for (int i = 0; i< (2*ns - 3); i++)
            for (int j = 0; j < (2*ns - 3); j++)
                matrixC[nm][i][j] = (rm[nm][i]*rm[nm][j])*matrixL[nm][i][j];
    }

    // solve linear equations for multipliers
    for (int nm = 0; nm < nMolecules; nm++)
    {
        if (molecules[nm].numAtoms <= 1)
            continue;
        int num = 2*molecules[nm].numAtoms - 3;

        // lu_decomposite(matrixC[nm], num);
        // lu_solver(matrixC[nm], vecB[nm], num);
        lu_decomposite(matrixC[nm], index, &d, num);
        lu_solver(matrixC[nm], vecB[nm], index, num);
    }

    for (int i = 0; i < nAtoms; i++)
    {
        velfb[i] = 0.0;
        momfb[i] = 0.0;
    }

    // compute constraint force & feedback  
    for (Int nm = 0; nm < nMolecules; nm++)
    {
        int ns = molecules[nm].numAtoms;
        if (ns <= 1)
            continue;

        for (int i = 0; i < ns; i++)
        {
            Vector3 fc, vc, mc;
            Atom* atom = molecules[nm].myAtoms[i];
            Double Rij2, Rij;

            for (int j = 0; j < (2*ns - 3); j++)
            {
                // constraint force part
                val = matrixM[nm][i][j]*vecB[nm][j];
                // velocity feedback part
                Rij2 = rm[nm][j]*rm[nm][j];
                vcfb = matrixM[nm][i][j]*(1.0 - bondLength[nm][j]/Rij2);
                // momentum feedback part
                mcfb = matrixM[nm][i][j]*(rm[nm][j]*vm[nm][j])/Rij2;
                fc += val*rm[nm][j];
                vc += vcfb*rm[nm][j];
                mc += mcfb*rm[nm][j];
            }
            atom->force += fc;
            velfb[atom->atomID] = CFB*vc;
            momfb[atom->atomID] = DFB*mc;
            virial[XX] += fc.x * atom->position.x; 
            virial[XY] += fc.y * atom->position.x;
            virial[XZ] += fc.z * atom->position.x;
            virial[YX] += fc.x * atom->position.y;
            virial[YY] += fc.y * atom->position.y;
            virial[YZ] += fc.z * atom->position.y;
            virial[ZX] += fc.x * atom->position.z;
            virial[ZY] += fc.y * atom->position.z;
            virial[ZZ] += fc.z * atom->position.z;
        }
    }

    for (Int k = XX; k <= ZZ; k++)
        myEnsemble->virial[k] += virial[k];
}

// given a matrix - mat[size][size], this method replaces it using LU decomposition method.
// The modified matrix 'mat' will be as an input to lu_solver() to solve linear equations.
// This is a simplified implementation by omiting the search for the largest pivot element
// therefore particularly efficient for the matrix having the largest elements in the diagonal
// positions, such as used here, but it may not stable if some diagonal elements become too small.
void ConstraintForce :: lu_decomposite(Double** mat, Int size)
{
    #ifdef DEBUG
        if (mat==NULL)
            ERRORMSG("'mat' - null pointer");
    #endif
    for (int col = 0; col < size; col++)
    {
        for (int row = 0; row <= col; row++)
        {
            Double sum = mat[row][col];
            for (int k = 0; k < row; k++)
                sum -= mat[row][k]*mat[k][col];
            mat[row][col] = sum;
        }
        
        Double pivot = 1.0/mat[col][col];
        for (int row = col + 1; row < size; row++)
        {
            Double sum = mat[row][col];
            for (int k = 0; k < col; k++)
                sum -= mat[row][k]*mat[k][col];
            mat[row][col] = pivot* sum;
        }
    }
}

// given a matrix - mat[size][size], this method replaces it using LU decomposition method.
// The modified matrix 'mat' will be as an input to lu_solver() to solve linear equations.
// This is a complete version copied from Numerical Recipes in C. Corrections are made on
// loop boundaries.
// Param: mat[size][size] - input as coefficient matrix of linear equations, output as 
//        LU decomposed matrix.
// Param: indx[size] - output recording the row index interchanged.
// Param: d - output as 1 or -1 indicating even or odd number of raws being interchanged.
void ConstraintForce :: lu_decomposite(Double** mat, Int* indx, Double* d, Int size)
{
    Int row, col, k, imax;
    Double big, dum, sum, temp;
    Double *vv;                 // store the implicit scaling of each row

    #ifdef DEBUG2
        DEBUGMSG("Performing LU decomposition");
    #endif

    vv = new Double [size];
    *d = 1;                     // no row interchanged yet

    // get implicit scaling values for each row
    for (row = 0; row < size; row++)
    {
        big = 0.0;
        for (col = 0; col < size; col++)
            big = MAX(fabs(mat[row][col]), big);
        if (big == 0)
            ERRORMSG("devided by zero error");
        vv[row] = 1.0/big;       // save the scaling
    }

    // for each column, apply Crout's method
    for (col = 0; col < size; col++)
    {
        for (row = 0; row < col; row++)
        {
            sum = mat[row][col];
            for (k = 0; k < row; k++)
                sum -= mat[row][k]*mat[k][col];
            mat[row][col] = sum;
        }
        big = 0.0;              // initialize for the search for the largest pivot element.
        for (row = col; row < size; row++)
        {
            sum = mat[row][col];
            for (k = 0; k < col; k++)
                sum -= mat[row][k]*mat[k][col];
            mat[row][col] = sum;
            if ((dum = vv[row]*fabs(sum)) >= big)
            {
                big = dum;
                imax = row;           // index of a row to be interchanged
            }
        }
        if (col != imax)            // need to interchange rows
        {
            for (k = 0; k < size; k++)
            {
                temp = mat[imax][k];
                mat[imax][k] = mat[col][k];
                mat[col][k] = temp;
            }
            *d = -(*d);
            vv[imax] = vv[col];     // change the scale factor
        }
        indx[col] = imax;             // record the index of the row interchanged
        if (mat[col][col] == 0.0)
            mat[col][col] = TINY;   // substitute zero with TINY for the diagonal element
        // finally, divide by the pivot element
        if (col != (size - 1))
        {
            dum = 1.0/mat[col][col];
            for (row = col+1; row < size; row++)
                mat[row][col] *= dum;
        }
    }
    delete [] vv;
}


// this method solves a set of linear equations, A*X = B, in combination of lu_decomposite().
// Input mat[size][size] is LU decomposition matrix determined by lu_decomposite().
// indx[size] is input as a vector indicating row interchanged by lu_decomposite().
// vec[size], as input, is the right hand vector B and as output contains the solutions of X.  
void ConstraintForce :: lu_solver(Double** mat, Double* vec, Int* indx, Int size)
{    
    Double temp, sum;
    Int ip;

    #ifdef DEBUG2
        DEBUGMSG("lu_solver");
    #endif

    // interchange 'vec[]'
    for (int row = 0; row < size; row++)
    {
        ip= indx[row];
        temp = vec[ip];
        vec[ip] = vec[row];
        vec[row] = temp;
    }
    // forward substitution
    for (int row = 1; row < size; row++)
    {
        sum = vec[row];
        for (int j = 0; j < row - 1; j++)
            sum -= mat[row][j]*vec[j];
        vec[row] = sum;
    }
    // backward substitution
    for (int row = size -1; row >= 0; row--)
    {
        sum = vec[row];
        for (int j = row+1; j < size; j++)
            sum -= mat[row][j]*vec[j];
        vec[row] = sum/mat[row][row];
    }
}

// this simplified method solves a set of linear equations, A*X = B, in combination of 
// the simplified function of lu_decomposite().
void ConstraintForce :: lu_solver(Double** mat, Double* vec, Int size)
{
    Double sum;

    #ifdef DEBUG2
        DEBUGMSG("lu_solver");
    #endif

    // forward substitution
    for (int row = 1; row < size; row++)
    {
        sum = vec[row];
        for (int j = 0; j < row - 1; j++)
            sum -= mat[row][j]*vec[j];
        vec[row] = sum;
    }
    // backward substitution
    for (int row = size -1; row >= 0; row--)
    {
        sum = vec[row];
        for (int j = row+1; j < size; j++)
            sum -= mat[row][j]*vec[j];
        vec[row] = sum/mat[row][row];
    }    
}

Vector3* ConstraintForce::get_velocity_feedback()
{
    return velfb;
}

Vector3* ConstraintForce::get_momentum_feedback()
{
    return momfb;
}

void ConstraintForce::write_force_info(ofstream& of)
{
    of << "Constraint Force created" << endl;
    of << "CFB , DFB:  " << CFB << ",  " << DFB << endl;
    of << "selector matrixM" << endl;
    for (int i = 0; i < nMolecules; i++)
    {
        int ns = molecules[i].numAtoms;
        of << "mol[" << i << "]:" << endl;
        for (int row = 0; row < ns;  row++)
        {
            for (int col = 0; col < (2*ns-3); col++)
                of << matrixM[i][row][col] << ' ';
            of << endl;
        }
    }

    of<< "difference matrixL" << endl;
    for (int i = 0; i < nMolecules; i++)
    {
        int ns = molecules[i].numAtoms;
        of << "mol[" << i << "]:" << endl;
        for (int row = 0; row < (2*ns-3);  row++)
        {
            for (int col = 0; col < (2*ns-3); col++)
                of << matrixL[i][row][col] << ' ';
            of << endl;
        }
    }
        
    of << "bond length vector" << endl;
    for (int i = 0; i < nMolecules; i++)
    {
        int ns = molecules[i].numAtoms;
        of << "mol[" << i << "]:" << endl;
        for (int row = 0; row < (2*ns-3);  row++)
            of << bondLength[i][row] << ' ';
        of << endl;
    }
}

void ConstraintForce::write_energy(ofstream& of)
{
    ;
}

// JC string ConstraintForce::get_force_id()// in order to comply with the cluster compiler
// JC {
// JC    return ConstraintForce::forceIdentifier;
// JC}
