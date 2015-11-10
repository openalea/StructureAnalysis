/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       TreeMatching : Comparison of Tree Structures
 *
 *       Copyright 1995-2009 UMR LaBRI
 *
 *       File author(s): P.ferraro (pascal.ferraro@labri.fr)
 *
 *       $Source$
 *       $Id: mdtable.h 3258 2007-06-06 13:18:26Z dufourko $
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


#ifndef SB_MATCHING_DISTANCE_TABLE
#define SB_MATCHING_DISTANCE_TABLE

#include "definitions.h"
#include "nodecost.h"
#include "treegraph.h"
#include "distancetable.h"
/**
 *\class MatchingDistanceTable
 *\brief Retains the distance between all the subtrees of two tre-graphs.
 * The distance is computing from a given NodeCost.
 *\par Requirements
 *  - Two TreeGraph;
 *  - A NodeCost
 *\author Pascal ferraro
 *\date 2009
 */
class MatchingDistanceTable;
/** Type of  MatchingDistanceTable */
#define STD 0
#define COMPACT 1
//enum MDTableType { STD, COMPACT};

typedef int MDTableType;

typedef std::vector<double> DistanceVector;
typedef std::vector<DistanceVector> DistanceVectorTable;




class TREEMATCH_API MatchingDistanceTable 
{
//   friend class Matching;
//   friend class MatchPath;

  public :
  MatchingDistanceTable(){};
  //MatchingDistanceTable(const TreeGraphPtr& ,const TreeGraphPtr& ,const NodeCostPtr);
  virtual ~MatchingDistanceTable(){};

  virtual DistanceType getDBT(int ,int ) const = 0;
  virtual DistanceType getDBF(int ,int ) const = 0;
  virtual void putDBT(int ,int ,DistanceType)   = 0;
  virtual void putDBF(int ,int ,DistanceType )  = 0;
  //Initialisation
  virtual void putInputTreeToEmpty(int,DistanceType){};
  virtual DistanceType getInputTreeToEmpty(int){return 0.;};
  virtual void putReferenceTreeFromEmpty(int,DistanceType){};
  virtual DistanceType getReferenceTreeFromEmpty(int){return 0.;};

  virtual void openDistancesVector(int input_vertex){};
  virtual void closeDistancesVector(int input_vertex){};

  
  virtual DistanceType inputTreeToEmpty(int ){return 0.;};
  virtual DistanceType inputForestToEmpty(int ){return 0.;};
  virtual DistanceType referenceTreeFromEmpty(int ){return 0.;};
  virtual DistanceType referenceForestFromEmpty(int ){return 0.;};
  virtual DistanceType getICost(int& ) const {return 0.;};
  virtual DistanceType getDCost(int& ) const {return 0.;};
  virtual DistanceType getCCost(int ,int ) const {return 0.;};
  virtual DistanceType getMCost(vector<int> ,int ) const {return 0.;};
  virtual DistanceType getSCost(int ,vector<int> ) const {return 0.;};
  virtual DistanceTable &getTreeDistanceTable()   = 0;
  virtual DistanceVectorTable &getDistanceTable() = 0; 
  MDTableType getType() const {return _type;};
  
 protected :
  MDTableType _type;

  // public :
//    int _n1;
//    int _n2;
//    TreeGraphPtr T1;
//    TreeGraphPtr T2;
//    NodeCostPtr ND;


};

class TREEMATCH_API StdMatchingDistanceTable : public MatchingDistanceTable
{
//   friend class Matching;
//   friend class MatchPath;

  public :
  StdMatchingDistanceTable(){};
  StdMatchingDistanceTable(const TreeGraphPtr& ,const TreeGraphPtr& ,const NodeCostPtr);
  ~StdMatchingDistanceTable(){};
  DistanceType getDBT(int ,int ) const ;
  DistanceType getDBF(int ,int ) const ;
  void putDBT(int ,int ,DistanceType ) ;
  void putDBF(int ,int ,DistanceType ) ;
  //Initialisation
  void putInputTreeToEmpty(int,DistanceType);
  DistanceType getInputTreeToEmpty(int);
  void putReferenceTreeFromEmpty(int,DistanceType);
  DistanceType getReferenceTreeFromEmpty(int);
  DistanceType getICost(int& ) const ;
  DistanceType getDCost(int& ) const ;
  DistanceType getCCost(int ,int ) const;
  DistanceType getMCost(vector<int> ,int ) const;
  DistanceType getSCost(int ,vector<int> ) const;
 
  DistanceType inputTreeToEmpty(int );
  DistanceType inputForestToEmpty(int );
  DistanceType referenceTreeFromEmpty(int );
  DistanceType referenceForestFromEmpty(int );

  virtual DistanceVectorTable& getDistanceTable() { return(_treeDistTable); }
  virtual DistanceTable& getTreeDistanceTable() { assert(1);
    // just a compilation trick. unuseful.
	static DistanceTable DUMMYTABLE;
	return DUMMYTABLE; 
  };
 
    private :
    int _n1;
    int _n2;
    TreeGraphPtr T1;
    TreeGraphPtr T2;
    NodeCostPtr ND;


  protected :
  //DistanceVector _inputTreeToEmpty;
  //DistanceVector _referenceTreeFromEmpty;
  DistanceVectorTable _treeDistTable;
  DistanceVectorTable _forestDistTable;
};


class TREEMATCH_API CompactMatchingDistanceTable : public MatchingDistanceTable
{
//   friend class Matching;
//   friend class MatchPath;

  public :
  CompactMatchingDistanceTable(){};
  CompactMatchingDistanceTable(const TreeGraphPtr& ,const TreeGraphPtr& ,const NodeCostPtr);
  ~CompactMatchingDistanceTable(){};
  DistanceType getDBT(int ,int ) const ;
  DistanceType getDBF(int ,int ) const ;
  void putDBT(int ,int ,DistanceType ) ;
  void putDBF(int ,int ,DistanceType ) ;
  //Initialisation
  void putInputTreeToEmpty(int,DistanceType);
  DistanceType getInputTreeToEmpty(int);
  void putReferenceTreeFromEmpty(int,DistanceType);
  DistanceType getReferenceTreeFromEmpty(int);
  DistanceType getICost(int& ) const ;
  DistanceType getDCost(int& ) const ;
  DistanceType getCCost(int ,int ) const;
  DistanceType getMCost(vector<int> ,int ) const;
  DistanceType getSCost(int ,vector<int> ) const;

  void openDistancesVector(int input_vertex);
  void closeDistancesVector(int input_vertex);
 
  DistanceType inputTreeToEmpty(int );
  DistanceType inputForestToEmpty(int );
  DistanceType referenceTreeFromEmpty(int );
  DistanceType referenceForestFromEmpty(int );
  virtual DistanceTable &getTreeDistanceTable() {
    return(_treeDistTable);};
  virtual DistanceVectorTable &getDistanceTable() { exit(1); }

    private :
    int _n1;
    int _n2;
    TreeGraphPtr T1;
    TreeGraphPtr T2;
    NodeCostPtr ND;

  protected :
  DistanceVector _inputTreeToEmpty;
  DistanceVector _referenceTreeFromEmpty;
  DistanceTable _treeDistTable;
  DistanceTable _forestDistTable;
};




#endif 

