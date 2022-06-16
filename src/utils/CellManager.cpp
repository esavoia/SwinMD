/** CellManager.cpp -- generate & maintain cell and cell list
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

#include "CellManager.h"

CellManager::CellManager()
{
    cellNum = 0;    numOnX = 0;    numOnY = 0;    numOnZ = 0;
    cellLen = 0.0;  lenOnX = 0.0;  lenOnY = 0.0;  lenOnZ = 0.0;
    rLenOnX = 0.0;  rLenOnY = 0.0; rLenOnZ = 0.0;
    nbrCellNum = 13;         // default value
    cellList = NULL;
}

CellManager::~CellManager()
{
    if (cellList != NULL)
    {
        for (Int i = 0; i < cellNum; i++)
        {
            if (cellList[i].nbrCellList != NULL)
                delete [] cellList[i].nbrCellList;
            if (cellList[i].atomList != NULL)
                delete [] cellList[i].atomList;
        }
        delete [] cellList;
    }
}

// given a cell length value and the dimensions of the simulation box, we first 
// compute the number of cells in each directions and then the accurate size of the cell.
// The cubic cell is only true when the simulation box is cubic.
// In future, we need consider shear & elongation flows, ref. FJC5 init_neighbour & cell
bool CellManager::set_cellManager(Double len, Double lxBox, Double lyBox, Double lzBox)
{
    #ifdef DEBUG
        DEBUGMSG("set_cellManager");
    #endif

    // check if cell length has a proper value
    if ((len <= 0)||(len > lxBox/3.0)||(len > lyBox/3.0)||(len > lzBox/3.0))
        return false;

    boxLx = lxBox;
    boxLy = lyBox;
    boxLz = lzBox;
    halfLx = boxLx/2.0;
    halfLy = boxLy/2.0;
    halfLz = boxLz/2.0;

    cellLen = len;
    numOnX  = (Int)(lxBox/len);
    numOnY  = (Int)(lyBox/len);
    numOnZ  = (Int)(lzBox/len);
    numOnXY = numOnX*numOnY;
    cellNum = numOnXY*numOnZ;

    lenOnX  = lxBox/numOnX;
    lenOnY  = lyBox/numOnY;
    lenOnZ  = lzBox/numOnZ;
    rLenOnX = 1.0/lenOnX;
    rLenOnY = 1.0/lenOnY;
    rLenOnZ = 1.0/lenOnZ;

    origin0 = Vector3(-lxBox/2.0, -lyBox/2.0, -lzBox/2.0);

    #ifdef DEBUG
        cerr << "cell length: " << cellLen << endl;
        cerr << "box length: " << boxLx << ' ' << boxLy << ' ' << boxLz << endl;
        cerr << "cell number: " << numOnX << ' ' << numOnY << ' ' << numOnZ << ' ' << cellNum<<endl;
        cerr << "cell origin: " << origin0.x << ' ' << origin0.y << ' ' << origin0.z << endl;
    #endif

    return true;
}

// nbrCellNum is set to 13 by default, but also allowed to set other values when required. 
void CellManager::set_nbr_cell_size(Int numNbrCells)
{
    nbrCellNum = numNbrCells;
}

// build_cell_list() generates cell list and set index, origin and ranks in x, y and z 
// directions for each cell. 
void CellManager::build_cell_list()
{
    Int xi, yi, zi, id;
    Vector3 tmp;

    #ifdef DEBUG
        DEBUGMSG("building cell list");
    #endif

    if (cellList != NULL)
        delete [] cellList;

    cellList = new Cell[cellNum];
    if (cellList == NULL)
        ERRORMSG("memory allocation error for cell list");
         
    id = 0;
    for (zi = 0; zi < numOnZ; zi++)
    {
       for (yi = 0; yi < numOnY; yi++)
       {
          for (xi = 0; xi < numOnX; xi++)
          {
             // id = cell_index(i, j, k);
             cellList[id].x = xi;
             cellList[id].y = yi;
             cellList[id].z = zi;
             cellList[id].nbrCellList = new Int[nbrCellNum+1];
             if(cellList[id].nbrCellList == NULL)
                ERRORMSG("memory allocation error for neighbor cell list");

             tmp = Vector3(xi*lenOnX, yi*lenOnY, zi*lenOnZ);
             // cellList[id].origin = tmp + origin0;             
             cellList[id].listSize = 20;        
             cellList[id].currentSize = 0;
             cellList[id].atomList = new Atom*[20];

             for (Int i = 0; i < 20; i++)
                 cellList[id].atomList[i] = NULL;

             // id = xi + yi*numOnX + zi*numOnX*numOnY
             id++;
          }
       }
    }
    for (Int i = 0; i < cellNum; i++)
        for (Int j = 0; j <= nbrCellNum; j++)
            cellList[i].nbrCellList[j] = -1;
}

// set_atom() sets an atom into a cell's atomList
void CellManager::set_atom(Atom* atom)
{
    #ifdef DEBUG
        if (atom == NULL)
            ERRORMSG("NULL atom in set_atom()");
    #endif

    // find which cell this atom belongs to
    Int id = find_cell(atom->position);
    #ifdef DEBUG
        if (id >= cellNum)
        {
            cerr << "ERO: wrong cell ID: " << id << endl;
            cerr << "atom position: "<<atom->position.x<<' '<<atom->position.y<<' '<<atom->position.z<<endl;
            ERRORMSG("Cell ID outoff range");
        }
    #endif

    Int size = cellList[id].currentSize;

    // add the atom into the atomList of this cell
    if(size >= cellList[id].listSize)
       resize_atom_list(cellList[id]);
    cellList[id].atomList[size] = atom;
    cellList[id].currentSize++;
}

// clear & reset atom list for each cell
void CellManager::clear_atomlist()
{
    Int i;
    Int num = numOnXY*numOnZ;

    for ( i = 0; i < num; i++)
    {
        if (cellList[i].atomList != NULL)
           delete [] cellList[i].atomList;
        // reset stom list
        cellList[i].atomList = new Atom*[20];
        cellList[i].listSize = 20;
        cellList[i].currentSize = 0;

        for(Int j = 0; j < 20; j++)
            cellList[i].atomList[j] = NULL;
    }
}

// resize cell's atom list
void CellManager::resize_atom_list(Cell& cell)
{
   Int i;
   Atom** tempList;

   cell.listSize += 10;
   tempList = new Atom*[cell.listSize];
   if(tempList == NULL)
   {
      printf("%s", "resize atom list error");
      exit(1);
   }

   for ( i = 0; i < cell.currentSize; i++)
      tempList[i] = cell.atomList[i];
   for (i = cell.currentSize; i < cell.listSize; i++)
      tempList[i] = NULL;
   
   delete [] cell.atomList;
   cell.atomList = tempList;
}


void CellManager::update_pair_list(Int** exclusionTable, Double delta)
{
    Int i, j, k, l, selfSize, nbrSize, nbrCellID;
    Int atomID1, atomID2, atomID;
    Double r2;
    Vector3 dif;
    Cell *selfCell, *nbrCell;
    Atom *atom1, *atom2;                  // pair atoms
    bool isExcluded;    

    #ifdef DEBUG2
        cerr << "DBG: cell lists" << endl;
        for (i = 0; i < cellNum; i++)         // for each cell
        {
            cerr<< "cell[" << i << "]: ";
            selfCell = &cellList[i];      
            if (selfCell == NULL)
            {
                cerr << "ERO: on cell list [" << i <<']' << endl;
                ERRORMSG("NULL point");
            }
            for (j = 0; j < selfCell->currentSize; j++)
                cerr << selfCell->atomList[j]->atomID << ' ';
            cerr<< endl;
        }
        cerr << endl; 
        cerr << endl;
    #endif

    for (i = 0; i < cellNum; i++)         // for each cell
    {
       selfCell = &cellList[i];
       selfSize = selfCell->currentSize;      // number of atoms in this cell's atomlist
       for (j = 0; j < selfSize; j++)         // for each atom in this cell
       {
           atom1 = selfCell->atomList[j];
           atomID1 = atom1->atomID;

           // pairs in the same cell, no need apply PBC on distance
           // for (k = 0; k < selfSize; k++)
           // atomList should be in increament order if use Newton's third law, "for" loop should be
           for (k = j+1 ; k < selfSize; k++)   
           {
               atom2 = selfCell->atomList[k];
               atomID2 = atom2->atomID;
               isExcluded = false;

               Int m = 0;
               while ((atomID = exclusionTable[atomID1][m++]) != -1)
               {
                   if (atomID2 == atomID)
                   {
                       isExcluded = true;           // excluded 
                       break;
                   }
               }
               if (!isExcluded)
               {
                  dif = atom2->position - atom1->position;
                  r2 = dif.length2();                  // square distance
                  if (r2 < delta)
                     atom1->set_pair(atom2);
               }
           }

           // pairs in neighbour cells 
           k = 0;
           while ((nbrCellID = selfCell->nbrCellList[k]) > -1)
           {
               nbrCell = &cellList[nbrCellID];
               nbrSize = nbrCell->currentSize;
               for (l = 0; l < nbrSize; l++)              // for atoms in this neighbor cell
               {
                   atom2 = nbrCell->atomList[l];
                   atomID2 = atom2->atomID;
                   isExcluded = false;
                   Int m = 0;
                   while ((atomID = exclusionTable[atomID1][m++]) != -1)
                   {
                       if (atomID2 == atomID)
                       {
                           isExcluded = true;           // excluded
                           break;
                       } 
                   }
                   if (!isExcluded)                     // should apply PBC for pairs in neighbour cells
                   {
                      dif = atom2->position - atom1->position;
                      apply_pbc(dif);
                      r2 = dif.length2();              
                      if (r2 < delta)
                         atom1->set_pair(atom2);
                   }
               }
               k++;                                     // next neighbor cell
           }
        }
    }    
}


// update pairList for each atom in every cell
// @params: exclusionDelta - value specifying exclusion policy
//          delta - square of cutoff plus buffer size for determining pairs.
/*  void CellManager::update_pair_list(Int exclusionDelta, Double delta)
{
    Int i, j, k, l, selfSize, nbrSize, nbrCellID;
    Double r2;
    Vector3 dif;
    Cell *selfCell, *nbrCell;
    Atom *atom1, *atom2;                  // pair atoms

    #ifdef DEBUG2
        DEBUGMSG("CellManager updates pair list");
    #endif

    #ifdef DEBUG2
        cerr << "DBG: cell lists" << endl;
        for (i = 0; i < cellNum; i++)         // for each cell
        {
            cerr << "cell[" << i << "]: ";
            selfCell = &cellList[i];      
            if (selfCell == NULL)
            {
                cerr << "ERO: on cell list [" << i <<']' << endl;
                ERRORMSG("NULL point");
            }
            for (j = 0; j < selfCell->currentSize; j++)
                cerr << selfCell->atomList[j]->atomID << ' ';
            cerr << endl;
        }
    #endif

    for (i = 0; i < cellNum; i++)         // for each cell
    {
       selfCell = &cellList[i];

       selfSize = selfCell->currentSize;      // number of atoms in this cell's atomlist
       for (j = 0; j < selfSize; j++)         // for each atom in this cell
       {
           atom1 = selfCell->atomList[j];

           // pairs in the same cell
           // for (k = 0; k < selfSize; k++)
           // atomList should be in increament order if use Newton's third law, "for" loop should be
           for (k = j+1; k < selfSize; k++)   
           {
               atom2 = selfCell->atomList[k];
               if (!exclusion(atom1, atom2, exclusionDelta))
               {
                  dif = atom2->position - atom1->position;
                  r2 = dif.length2();                  // square distance
                  if (r2 < delta)
                     atom1->set_pair(atom2);
               }
           }

           // pairs in neighbour cells 
           k = 0;
           while ((nbrCellID = selfCell->nbrCellList[k]) > -1)
           {
               nbrCell = &cellList[nbrCellID];
               nbrSize = nbrCell->currentSize;
               for (l = 0; l < nbrSize; l++)              // for atoms in this neighbor cell
               {
                   atom2 = nbrCell->atomList[l];
                   if (!exclusion(atom1, atom2, exclusionDelta))
                   {
                      dif = atom2->position - atom1->position;
                      r2 = dif.length2();              
                      if (r2 < delta)
                         atom1->set_pair(atom2);
                   }
               }
               k++;                                     // next neighbor cell
           }
        }
    }       
}
**/

