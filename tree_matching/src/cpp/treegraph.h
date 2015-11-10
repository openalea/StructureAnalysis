/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2000 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): P.ferraro (pascal.ferraro@cirad.fr)
 *
 *       $Source: /usr/cvsmaster/AMAPmod/src/TreeMatching/treegraph.h,v $
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


#ifndef SB_TREE_GRAPH_HEADER
#define SB_TREE_GRAPH_HEADER

#include "definitions.h"
#include "mtg/mtg.h"
#include "treenode.h"
#include "distancetable.h"
#include <list>
#include <vector>
#include <iterator>

typedef std::vector<int> NodeList;
typedef std::vector<int> ClassList;
typedef std::vector<NodeList*> NodeTable;
typedef  int node;
typedef  int edge;


/**
 * List of TreeNode
 */
// typedef RWTPtrSlist<TreeNode> TreeNodeList;
typedef std::vector<TreeNode*> TreeNodeList;

enum TreeType  { TOPO , COMPO };
using namespace std;


/**
 *\class TreeGraph
 *\brief Definition of a tree graph
 *\author Pascal ferraro
 *\date 1999
 */

class TreeGraph
{

  public :
  typedef const TreeNode* const_iterator;
    /** Constructor */
    TreeGraph();
    TreeGraph(MTG& ,VId ,TreeType,EType,AmlBoolean valued=FALSE,NodeFunctionList* funcions=NULL);
    TreeGraph(MTG& ,VIdList* ,TreeType,EType,AmlBoolean valued=FALSE,NodeFunctionList* funcions=NULL);
    
    TreeGraph(MTG& ,VId ,TreeType,AmlBoolean valued=FALSE,NodeFunctionList* funcions=NULL);
    TreeGraph(MTG& ,VId ,int,AmlBoolean valued=FALSE,NodeFunctionList* functions=NULL);
    TreeGraph(char* fich_name, MTG&, VId root, AmlBoolean valued=FALSE,NodeFunctionList* functions=NULL);
  
  TreeGraph(TreeGraph*,int);

    /** Destructor */
    ~TreeGraph();

  int getRealNumber(int postfix) const;

    /** Return the father of a node */
    int father(int ) const;

    /** Give a child of the ChildrenList */
    int child(int ,int) const ;

    /** Give a the left brother of child */
    int leftBrother(int) const ;

    /** Give a the right brother of child */
    int rightBrother(int) const ;

    /** return 1 if a given child of a vertex is
     *  child is in the same complex of vertex */
    int childIsInComplex(int ,int ) const;
	int childIsInAxis(int ,int ) const;

    /**father isin the same complex of vertex */
    int fatherIsInComplex(int) const;

    /** return TreeNode corresponding to vertex */
    TreeNode* getTreeNode(int vertex) const;

//   /** remove A subTree corresponding to a given vertex **/
     void delSubTree(int vertex);

    /** Give the SonsList */
    const NodeList* sons(int ) const ;

  /** Return the root number */
  int getRoot() const;
  
  /** Return the depth of the tree */
  int getDepth() const;
  
  /** Return the degree of the tree */
    int getDegree() const;

    /** Return the order of the tree */
    int getOrder() const;

    /** Return the number of child of a node  */
    int getNbChild(int ) const;
    int getNbDesc(int ) const;

    /** Return the number of node in the tree */
    int getNbVertex() const;
    int getNbClass() const;
    int getNbAxisVertex() const;

    /** Return wether the node is a leaf or not          */
  int isLeaf(int ) const;
  int getChildLess( int ) const;
  void printNodedClass() const;
    /** Return a node */
    TreeNode* getNode(int ) const;
	int getNumber(int vertex) const;

    /** Return the edge type */
    EType getEdgeType() const { return(_edge); }

    /** Return a pointer to the mtg the tree is linked with  */
    MTG* getMTG() const { return(_mtg); }

    /** Return wheter the tree is null or not  */
    int isNull();

    /** Method for printing */
    void print();
    ostream& mtg_write(ostream &os) const;
  bool mtg_write( const char *path ) const;

    void upDate(MTG& );

    int outdeg(const node);
    std::list<int>::const_iterator getOutNodes(const node);
   
    void removeSon(int vertex,int son,int decalage);
  void shiftSonNumber(int vertex,int son,int decalage);
  void addSubTree(TreeGraph* T,int root, int insertion_point) ;

  int minimalClass(int node) const;
  int rightClasse(int clas);
  NodeList* getRootsInForestClass(int clas);
  
  NodeList* getRootsInClass(int clas);
  int getNbRootsInClass(int clas);

  private :
  MTG* _mtg;
  VId _root;
  int _nbVertex;
  int _nbClass;
  int _nbAxisVertex;
  int _degree;
  int _depth;
  int _number;
  int _numPostfix;
  int _numClass;
  std::vector<int> _postfix;
  std::vector<int> _class;
  int _order;
  EType _edge;
  AmlBoolean _valued;
  NodeFunctionList* _functions;
  int depth(int );
  NodeTable _outNodeTable;
  NodeTable _classNodeTable;
  TreeNodeList _treenodes;
  void putNode(int,int);
  void putNodeinClass(int,int);
  void putTreeNode(TreeNode*);
  void putNodeList(NodeList* );
  void putClassNodeList(NodeList* );
  void makeNode(VId ,int ,VId,TreeType );
  void makeNode(VId ,int ,VId ,TreeType ,VIdList*);	
  void makeComplexNode(VId ,int ,VId,TreeType ); 
  void makeAxeNode(VId ,int ,int ,int );	
  
};


#endif


