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


#ifndef SB_TREE_MATCH_EXTRACT_HEADER
#define SB_TREE_MATCH_EXTRACT_HEADER

#include "definitions.h"
#include "treematch.h"

typedef float Rate;
typedef RWTValVector<Rate>   RateVector;
typedef RWTValDlist<int>    NumberList;

enum FileType { DAT=0, PLOT, PRINT };

enum EditOperation { DEL=0, INS , MAT , SUB };

class TreeMatchExtract
{
  public :
    TreeMatchExtract(MTG& mtg, TreeMatch& treematch );
  void statistics() ;
  void plot_statistics(const char* prefix_file_name) ;
  void plot_statistic(const char* prefix_file_name,int imput_tree,int reference_tree) ;
  void text_statistics(const char* prefix_file_name) ;
  void text_statistic(const char* prefix_file_name,int imput_tree,int reference_tree) ;
  void mapping() ;
  void editOpStat();
  private :
    // data members
    MTG* _mtg;
  TreeMatch* _treematch;
  RateVector _delVertexRate;
  RateVector _insVertexRate;
  RateVector _subVertexRate;
  RateVector _matVertexRate;
  NumberList _newOrder;
  // private functions
  EditOperation editOp(Sequence*, VId, int);
};

#endif	

