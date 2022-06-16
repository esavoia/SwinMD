/** CellManager.h -- generate & maintain cell and cell list,
 **     update pair list and neighbor cell list
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

#ifndef CellManager_H
#define CellManager_H

#include <stdlib.h>
#include "../NEMD_defs.h"
#include "Vector3.h"
#include "../Atom.h"
using namespace std;
/**
 ** Forward declarations
 **/

// struct Cell specifies a cell data structure in the simulation box
typedef struct cell {
    // location of theis cell
    // Vector3 origin;
    // rank in x, y & z directions
    Int     x;            
    Int     y;
    Int     z;
    // neighbour cell list - store cell index of neighbour cells
    Int*    nbrCellList;
    // atom list belongs to this cell
    Atom**  atomList;
    Int     listSize;
    Int     currentSize;
} Cell; 

class CellManager {
    /**
     ** Data Member - 
     **/
    public:
    Int     numOnX, numOnY, numOnZ, numOnXY;      // numOnXY = numOnX*numOnY, cellNum = numOnXY*numOnZ
    Int     cellNum,nbrCellNum;                   // number of neighbour cells, default 13
    Double  cellLen, lenOnX, lenOnY, lenOnZ;      // length of the cell
    Double  rLenOnX, rLenOnY, rLenOnZ;            // reverse of length
    Vector3 origin0;                              // origin of the Cell(0, 0, 0)
    Cell*   cellList;                             // list of all cells consisting of the sim box
    Double boxLx, boxLy, boxLz;
    Double halfLx, halfLy, halfLz;

    /**
     ** constructor and destructor
     **/
    public:
    CellManager();
    ~CellManager();

    /**
     ** methods
     **/
    public:   
    // set_cellManager sets the number of cells and cell dimmensions, return 'true' when success
    // if cell length is negative or larger than 1/3 of box length return 'false' & cell will be not set
    bool set_cellManager(Double len, Double lxBox, Double lyBox, Double lzBox);
    void set_nbr_cell_size(Int numNbrCells);
    void build_cell_list();
    void set_atom(Atom* atom);
    void clear_atomlist();
    void resize_atom_list(Cell& cell);
    void update_pair_list(Int** exclusionTable, Double delta);
    void update_pair_list(Int exclusionDelta, Double delta);

    void get_cell_size(Vector3& v) { v = Vector3(lenOnX, lenOnY, lenOnZ); }
    Double get_cell_volume() { return (lenOnX*lenOnY*lenOnZ); }

    // given ranks of a cell in x, y & z directions, get the cell's index as in cell list
    inline Int cell_index(Int x, Int y, Int z)
    {
        return (z*numOnXY + y*numOnX + x);
    }

    // given an atom's position, find the cell where the atom sits 
    // and return the index of the cell 
    // we let Ensemble to apply PBC on the atom before calling this function
    // this guarantees the atom is within the sim-box when this function is called
    inline Int find_cell(const Vector3 &position)
    {
        Int x, y, z;           // cell ranks in each directions         
        Vector3 delta;                  
                                        
        delta = position - origin0;     
        x = (Int)(delta.x*rLenOnX);
        y = (Int)(delta.y*rLenOnY);
        z = (Int)(delta.z*rLenOnZ);

        #ifdef DEBUG
            if ((x>= numOnX)||(x<0)||(y>= numOnY)||(y<0)||(z >= numOnZ)||(z<0))
            {
                DEBUGMSG(" ");
                cerr << "improper cell ranks: " << x << ' ' << y << ' ' << z << ' ' << endl;
                cerr << "position: " << position.x << ' ' << position.y << ' ' << position.z << endl;
            }
        #endif

        // due to the computing precision, after applying PBC on atom it's still possible that 
        // an atom is little bit out of boundary, this will cause program crash when setting atoms
        // to cells, here we force the cell ranks of any atom within the range
        
        if(x >= numOnX)     x = numOnX - 1;
        else if(x < 0)      x = 0;
        if(y >= numOnY)     y = numOnY - 1;
        else if(y < 0)      y = 0;
        if(z >= numOnZ)     z = numOnZ - 1;
        else if(z < 0)      z = 0; 
        return (z*numOnXY + y*numOnX + x);
    }  

    // this function build neighbor cell list for each central cell
    // minimum image and Newton's third law are used and each cell will
    // have 13 neighbor cells.
    inline void build_cell_neighbor()
    {
        #ifdef DEBUG
            DEBUGMSG("building cell neighbor");
        #endif

        Int xPlus, yPlus, zPlus;
        Int xMinus, yMinus, zMinus;
        Int myID = 0;
        Cell *myCell;
        Int maxX = numOnX - 1;
        Int maxY = numOnY - 1;
        Int maxZ = numOnZ - 1;
        
        for(int zi = 0; zi < numOnZ; zi++)       // for all cells
        {
           for(int yi = 0; yi < numOnY; yi++)
           {
              for(int xi = 0; xi < numOnX; xi++)
              {
                  myCell = &cellList[myID];      // the central cell

                  if(xi < maxX)  xPlus = xi + 1;
                  else xPlus = 0;
                  if(yi < maxY)  yPlus = yi + 1;
                  else yPlus = 0;
                  if(zi < maxZ)  zPlus = zi + 1;
                  else zPlus = 0;
                  if(xi > 0)  xMinus = xi - 1;
                  else xMinus = maxX; 
                  if(yi > 0)  yMinus = yi - 1;
                  else yMinus = maxY;
                  if(zi > 0)  zMinus = zi - 1;
                  else zMinus = maxZ;

                  // shift (0 0 -1)
                  myCell->nbrCellList[0]  = cell_index(xi, yi, zMinus);
                  // shift (1 0 -1)
                  myCell->nbrCellList[1]  = cell_index(xPlus, yi, zMinus);
                  // shift (1 0 0)
                  myCell->nbrCellList[2]  = cell_index(xPlus, yi, zi);
                  // shift (1 0 1)
                  myCell->nbrCellList[3]  = cell_index(xPlus, yi, zPlus);
                  // shift (-1 1 -1)
                  myCell->nbrCellList[4]  = cell_index(xMinus, yPlus, zMinus);
                  // shift (0 1 -1)
                  myCell->nbrCellList[5]  = cell_index(xi, yPlus, zMinus);
                  // shift (1 1 -1)
                  myCell->nbrCellList[6]  = cell_index(xPlus, yPlus, zMinus);
                  // shift (-1 1 0)
                  myCell->nbrCellList[7]  = cell_index(xMinus, yPlus, zi);
                  // shift (0 1 0)
                  myCell->nbrCellList[8]  = cell_index(xi, yPlus, zi);
                  // shift (1 1 0)
                  myCell->nbrCellList[9]  = cell_index(xPlus, yPlus, zi);
                  // shift (-1 1 1)
                  myCell->nbrCellList[10] = cell_index(xMinus, yPlus, zPlus);
                  // shift (0 1 1)
                  myCell->nbrCellList[11] = cell_index(xi, yPlus, zPlus);
                  // shift (1 1 1)               
                  myCell->nbrCellList[12] = cell_index(xPlus, yPlus, zPlus);

                  myCell->nbrCellList[13] = -1;     // Sentinel for end of neighbors
                  myID++;                           // next central cell
              }
           }
        }
    }

    inline void apply_pbc(Vector3 &r)
    {
        while (r.x >= halfLx)  r.x -= boxLx;
        while (r.x < -halfLx)  r.x += boxLx;
        while (r.y >= halfLy)  r.y -= boxLy;
        while (r.y < -halfLy)  r.y += boxLy;
        while (r.z >= halfLz)  r.z -= boxLz;
        while (r.z < -halfLz)  r.z += boxLz; 
    }
};



    /**
     ** below function build neighbor cell lists without forwllowing Minimum image convention
     ** therefore the cell at the  boundary or cooner of the simulation box does not necessary
     ** to have 13 neighbours.
    inline void build_cell_neighbor()
    {
        Int myID = 0;
        Cell *myCell;
        
        for(int zi = 0; zi < numOnZ; zi++)
        {
           for(int yi = 0; yi < numOnY; yi++)
           {
              for(int xi = 0; xi < numOnX; xi++)
              {
                  Int n = 0;                        // number of neighbor cells found so far
                  myCell = &cellList[myID];      // the central cell

                  if(xi<numOnX-1)  myCell->nbrCellList[n++] = myID + 1;
                  if(yi<numOnY-1)  myCell->nbrCellList[n++] = myID + numOnX;
                  if(zi<numOnZ-1)  myCell->nbrCellList[n++] = myID + numOnXY;
                  if(xi<(numOnX-1)&&yi<(numOnY-1))  myCell->nbrCellList[n++] = myID + numOnX + 1;
                  if(xi<(numOnX-1)&&zi<(numOnZ-1))  myCell->nbrCellList[n++] = myID + numOnXY + 1;
                  if(yi<(numOnY-1)&&zi<(numOnZ-1))  myCell->nbrCellList[n++] = myID + numOnXY + numOnX;
                  if(xi<(numOnX-1)&&yi>0)  myCell->nbrCellList[n++] = myID - numOnX + 1;
                  if(xi>0&&zi<(numOnZ-1))  myCell->nbrCellList[n++] = myID + numOnXY - 1;
                  if(yi>0&&zi<(numOnZ-1))  myCell->nbrCellList[n++] = myID + numOnXY - numOnX;
                  if(xi<(numOnX-1)&&yi<(numOnY-1)&&zi<(numOnZ-1))  
                      myCell->nbrCellList[n++] = myID + numOnX + numOnXY + 1;
                  if(xi>0&&yi<(numOnY-1)&&zi<(numOnZ-1))  
                      myCell->nbrCellList[n++] = myID + numOnX + numOnXY - 1;
                  if(xi<(numOnX-1)&&(yi>0)&&zi<(numOnZ-1))  
                      myCell->nbrCellList[n++] = myID - numOnX + numOnXY + 1;
                  if(xi>0&&(yi>0)&&zi<(numOnZ-1))  
                      myCell->nbrCellList[n++] = myID - numOnX + numOnXY - 1;

                  myCell->nbrCellList[n] = -1;     // Sentinel for end of neighbors
                  myID++;                          // next central cell
              }
           }
        }
    }
    **/

#endif
