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


#ifndef SB_TREE_NODE_HEADER
#define SB_TREE_NODE_HEADER

#include <list>
#include <iostream>
#include "definitions.h"
#include "mtg/mtg.h"
#include "amlnodefunction.h"


typedef DistanceType ValueType;

/**
 * A type for a list of function that can be applied to a TreeNode.
 */
typedef std::list<NodeFunction> NodeFunctionList;

/**
 * A type for attributes of TreeNodes.
 */
typedef std::vector<DistanceType> ValueVector;


/**
 *\class TreeNode
 *\brief Definition of a vertex in a tree graph
 *\author Pascal ferraro
 *\date 1999
 */


class TreeNode
{

  public :

    /** Default constructor
     *\par Remarks
     * Never used
     */
    TreeNode() : _number(0),
		 _numPostfix(0),
		 _depth(0),
		 _father(0),
		 _vertex(0),
		 _complex(0),
		 _value(0),
		 _order(0),
		 _position(0),
		 _nb_value(0),
		 _mtg(NULL),
		 _values(NULL) {};

    /** Copy constructor. */
  TreeNode(const TreeNode& inode):  
    _number(inode.getNumber()),
    _numPostfix(inode.getNumPostfix()),
    _depth(inode.depth()),
    _father(inode.father()),
    _vertex(inode.getVertex()),
    _complex(inode.getComplex()),
    _value(inode.getValue()),
    _order(inode.getOrder()),
    _position(inode.getPosition()),
    _nb_value(inode.getValueSize()),
    _mtg(inode.getMTG()),
    _values(NULL){
    //    cout <<"_nbvalue =" <<_nbvalue<<endl;
    resize(_nb_value);
    for(int i=0;i<_nb_value;i++)
      putValue(i,inode.getValue(i));
  }

    /** Constructs by default a TreeNode from a MTG /e mtg. */

  TreeNode(MTG& mtg) : _mtg(&mtg),_values(0),_nb_value(0){}

    /** constructs a TreeNode with /e mtg , /e number, /e depth, /e father, /e vertex, /e complex. */
    TreeNode(MTG& ,int ,int ,int ,int ,int );

    /** Destructor. */
    ~TreeNode();

    /** \par Reading Functions. */


    /** Returns the /e i value of TreeNode. */
    ValueType  value(int ) const ;

    /** Returns the depth of the TreeNode in the TreeGraph. */
    int depth() const { return(_depth);}

    /** Returns the father of TreeNode. */
    int father() const { return(_father);}

    /** Returns the the Vertex of the TreeNode in MTG.
     *  vertex is the reference of TreeNode in a MTG */
    int getVertex() const { return(_vertex);}

    /** Returns the number of the TreeNode in MTG.
     *  number is the reference of TreeNode in a TreeGraph */
  int getNumber() const { return(_number);}
  int getId() const { return(_id);}
 

  int getNumPostfix() const { return(_numPostfix);}
  
  /** Returns the complex of the TreeNode in MTG. */
    int getComplex() const { return(_complex);}

    /** Returns the value of the TreeNode in MTG. */
    int getValue() const { return(_value);}

    /** Returns the order of the TreeNode in MTG. */
    int getOrder() const { return(_order);}

    /** Returns the position of the TreeNode in MTG. */
    int getPosition() const { return(_position);}

    /** \par Writing Functions. */

    /** Puts the position of the TreeNode in MTG. */
    void putPosition(int position) { _position=position;}

    /** Puts the order of the TreeNode in MTG. */
    void putOrder(int order) { _order=order;}

    /** Puts the value of the TreeNode in MTG. */
    void putValue(int value) { _value=value;}

    /** Puts the complex of the TreeNode in MTG. */
    void putComplex(int complex) { _complex=complex;}

    /** Puts the number of the TreeNode in a TreeGraph.
     *  number is the reference of TreeNode in a TreeGraph */
  void putNumber(int number) { _number=number;}
void putId(int id) { _id=id;}
  void putNumPostfix(int number) { _numPostfix=number;}
  
  
  /** Puts the vertex of the TreeNode in MTG.
   *  vertex is the reference of TreeNode in a MTG */
  void putVertex(int vertex) { _vertex=vertex;}
  
  /** Puts the father of the TreeNode in MTG. */
    void putFather(int father ){ _father=father;}
  
  /** Puts the depth of the TreeNode in MTG. */
  void putDepth(int depth){ _depth=depth;}
    void put(std::ostream& os);
  
  /** Print a TreeNode*/
  void print();
  
    /** Changes the number of value affected to a TreeNode.*/
    void resize(int new_size);

    /** Returns the number of value affected to a TreeNode.*/
    int getValueSize() const;
    DistanceType getValue(int index) const;

    /** Puts a value /e new_value at teh index /e index to a TreeNode.*/
    void putValue(int index, DistanceType new_value);
    MTG* getMTG() const { return(_mtg);};
  void upDate(MTG& );
  
  int getAt1() const {return at1;};
  int getAt2() const {return at2;};
  int getAt3() const {return at3;};
  char getAt4() const {return at4;};
  
  
  private :
  int _number;
  int _numPostfix;
  int _depth;
  int _father;
  int _vertex;
  int _id;
  int _complex;
  int _value;
  int _order;
  int _position;
  int _nb_value;
  ValueVector* _values;
  MTG* _mtg;
  int at1;
  int at2;
  int at3;
  char at4;
};

#endif

