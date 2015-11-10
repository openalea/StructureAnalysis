/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       TreeMatching : Comparison of Tree Structures
 *
 *       Copyright 1995-2009 UMR LaBRI
 *
 *       File author(s): P.ferraro (pascal.ferraro@labri.fr)
 *
 *       $Source: /usr/cvsmaster/AMAPmod/src/TreeMatching/treegraph.h,v $
 *       $Id: treegraph.h 3264 2007-06-06 14:22:22Z dufourko $
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
#include "treenode.h"
#include <list>
#include <vector>
#include <stack>
#include <iterator>

/**
 * List of TreeNode
 */
typedef std::vector<TreeNodePtr> TreeNodeList;

/**
 * List of Children
 */

typedef std::vector<int> NodeList;


using namespace std;


/**
 *\class TreeGraph
 *\brief Definition of a tree graph
 *\author Pascal ferraro
 *\date 2009
 */

class TREEMATCH_API TreeGraph
{

  public :
  typedef TreeNodeList::const_iterator const_iterator;

  /** Constructor */
  TreeGraph();
  
  /** Destructor */
  ~TreeGraph();

  void addNode(int, int);

  void addNode(TreeNodePtr);

  /** return TreeNode corresponding to vertex */
  inline TreeNodePtr getTreeNode(int vertex) { return getNode(vertex); }
  TreeNodePtr getNode(int vertex);

  int getDepth() const { return _depth;};

   /** Return the father of a node */
  int father(int ) const;

  /** Give a child of the ChildrenList */
  int child(const int ,const int) const ;

  void addValue(int, DistanceType);
  DistanceType getValue(int, int) const;

  /** Give the SonsList */
  const NodeList& childList(const int ) const ;

  /** Return the number of child of a node  */
  int getNbChild(int ) const;
  int getNbDesc(int ) const;

  vector<int> getPath(int, int) const;

  /** Return the number of node in the tree */
  inline int getNbVertex() const { return _treenodes.size(); }

  /** Return the root id */
  inline int getRoot() const { return _rootId; }

  /** Return the root id */
  inline int getDegree() const{ return _degree; }
  
  /** Return whether the node is a leaf or not          */
  inline bool isLeaf(int node) const { return (getNbChild(node)==0); }

  /** Return wheter the tree is null or not  */
  inline bool isNull() {  return _treenodes.empty(); }

  /** Method for printing */
  void print() const;
  ostream& mtg_write(ostream &os) const;
  bool mtg_write(  char *path ) ;   



  private :
  // int _nbVertex;
  TreeNodeList _treenodes;
  int _depth;
  int _degree;
  int _rootId;

public:

  /*
    Definition of an iterator of the node if of the graph in pre order mode
  */
  struct TREEMATCH_API const_pre_order_iterator {
        friend class TreeGraph;
    protected:
        // pointer to the graph
        const TreeGraph * __treegraph;

        typedef NodeList::const_iterator TIter;

        typedef int value_t ;
        typedef const int * pointer_t;

        // pointer to actual node id in the graph
        TIter __current;
        // pointer to end of list of actual node id in the graph
        TIter __current_end;

        // type definition of a pair of pointer and a stack of pairs
        typedef std::pair<TIter,TIter> TIterPair;
        typedef std::stack<TIterPair> IteratorStack;
        // stack of higher value of pointers in the graph hierarchy
        IteratorStack __stack;

        // Pop pointers value from the stack
        inline void pop() 
        {  
           TIterPair p = __stack.top(); 
           __current = p.first;  __current_end = p.second; 
           __stack.pop();
        }

        // Push pointers value into the stack
        inline void push() { __stack.push( TIterPair(__current, __current_end)); }

        const_pre_order_iterator(TIter current, TIter current_end, const TreeGraph * treegraph):
            __current(current), __current_end(current_end), __treegraph(treegraph) {}

        // advance in the traversal
        void increment();

    public:
        bool  atEnd() const ;

        // pre increment operator for iterator
        inline const_pre_order_iterator& operator++() 
        { increment(); return *this; }

        // post incremental operator for iterator
        inline const_pre_order_iterator operator++(int i) 
        { const_pre_order_iterator original = *this; for (int j = 0 ; j < i; ++j) increment(); return original; }

        // value access operator for iterator
        inline value_t operator*() const { return *__current; }

        // return pointer to class object
        inline pointer_t operator->() const
		{	return __current.operator->(); }

        const_pre_order_iterator& operator+=(int i)
		{	// increment by integer
            for (int j = 0 ; j < i; ++j) increment();
            return *this;
        }

        const_pre_order_iterator operator+(int i)
		{	// increment by integer
            const_pre_order_iterator newiter = *this;
            for (int j = 0 ; j < i; ++j) newiter.increment();
            return newiter;
        }

        // compare two iterators
        inline bool operator==(const_pre_order_iterator other) {
            return __treegraph == other.__treegraph && __current == other.__current;
        }

        // compare two iterators
        inline bool operator!=(const_pre_order_iterator other) {
            return __treegraph != other.__treegraph || __current != other.__current;
        }
  };

  /// Get a begin iterator for the traversal of the subtree rooted in nodeid (nodeid is not given)
  inline const_pre_order_iterator subtree_iterator_begin(int nodeid) { 
      TreeNodePtr node = getNode(nodeid); 
      return const_pre_order_iterator(node->getChildList().begin(), node->getChildList().end(),this);
  }

  /// Get an end iterator for the traversal of the subtree rooted in nodeid
  inline const_pre_order_iterator subtree_iterator_end(int nodeid) { 
      TreeNodePtr node = getNode(nodeid); 
      return const_pre_order_iterator(node->getChildList().end(),  node->getChildList().end(),this);
  }

  /*
    To access to all the descendants of a node, one can use the previous iterator in the following way:
    TreeGraph t;
    ...
    int root = 0;
    for(TreeGraph::const_pre_order_iterator it = t.subtree_iterator_begin(root); it != t.subtree_iterator_begin(root); ++it)
    {
        int next_ascendant_it = *it;
        ...
     }
  */
};

// A smart pointer on a TreeGraph
typedef boost::shared_ptr<TreeGraph> TreeGraphPtr;
// A weak reference to avoid loop.
typedef boost::weak_ptr<TreeGraph> TreeGraphWeakPtr;

#endif


