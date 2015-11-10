/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2000 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): P.ferraro (pascal.ferraro@cirad.fr)
 *
 *       $Source$
 *       $Id$
 *
 *       Forum for AMAPmod developers    : amldevlp@cirad.fr
 *
 *  ----------------------------------------------------------------------------
 *
 *                      GNU General Public Licence
 *
 *       This program is free software; you can redistribute it and/or
 *       modify it under the terms of the GNU General Public License as
 *       published by the Free Software Foundation; either version 2 of
 *       the License, or (at your option) any later version.
 *
 *       This program is distributed in the hope that it will be useful,
 *       but WITHOUT ANY WARRANTY; without even the implied warranty of
 *       MERCHANTABILITY or FITNESS For A PARTICULAR PURPOSE. See the
 *       GNU General Public License for more details.
 *
 *       You should have received a copy of the GNU General Public
 *       License along with this program; see the file COPYING. If not,
 *       write to the Free Software Foundation, Inc., 59
 *       Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 *  ----------------------------------------------------------------------------
 */


/**
 *\class IntTable
 *\brief Defines an array of integer
 *\par Requirements
 * - r number of rows;
 * - c number of column.
 *\author Pascal ferraro
 *\date 04/2000
 */
// -------------------------------
// Classe: IntTable
// Cette classe gere un tableau
// d entier de r ligne et c
// colonnes
// -------------------------------

#ifndef SB_INT_TABLE_HEADER
#define SB_INT_TABLE_HEADER

#include <vector>
#include "indexmanager.h"

typedef std::vector<int> IntVector;
typedef std::vector<IntVector> IntVectorTable;


class IntTable
{

 public:
  IntTable(){};
  IntTable(int ,int ,int );
  void resize(int ,int ,int );
  ~IntTable();
  int getInt(int , int ) const ;
  void putInt(int , int , int );
  void openIntVector(int );
  void closeIntVector(int );
  int getColumnSize() const { return(_columnSize); };
  int getRowSize() const { return(_rowSize); };
  int getSimulatedSize() const { return(_simulatedSize); };

 private:
  IndexManager _indices;
  IntVectorTable _intTable;
  int _columnSize;
  int _rowSize;
  int _simulatedSize;
};

#endif



