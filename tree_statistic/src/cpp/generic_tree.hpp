/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2010 CIRAD/INRIA Virtual Plants
 *
 *       File author(s): J.-B. Durand (jean-baptiste.durand@imag.fr)
 *
 *       $Source$
 *       $Id: generic_tree.hpp 3186 2007-05-25 15:10:30Z dufourko $
 *
 *       Forum for V-Plants developers:
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
 *       MERCHANTABILITY or FITNESS for A PARTICULAR PURPOSE. See the
 *       GNU General Public License for more details.
 *
 *       You should have received a copy of the GNU General Public
 *       License along with this program; see the file COPYING. If not,
 *       write to the Free Software Foundation, Inc., 59
 *       Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 *  ----------------------------------------------------------------------------
 */

#ifndef GENERIC_TREE_TCC
#define GENERIC_TREE_TCC

using Tree_tie::tie;

/*--------------------------------------------------------------*
 *
 *  Generic algorithm for subtree selection, using a tree,
 *  the subtree root and a flag on keeping or pruning
 *  the selected subtree
 *
 *--------------------------------------------------------------*/

template <class Tree>
Tree* select_subtree(Tree& t,
                     typename tree_traits<Tree>::vertex_descriptor v,
                     bool keep_flag)
{
   typedef typename tree_traits<Tree>::vertex_descriptor vid;
   typename tree_traits<Tree>::children_iterator it, end;
   std::deque<vid> traversal_nodes,
                   subtree_nodes;
    // list of nodes used for tree traversal and for building result tree
   vid csnode, // current source node id
       *vid_array;
   int nvertices= 0; // number of subtree vertices
   unsigned int nvertices_init= 0; // number of initial vertices as set by constructor
   Tree* res= NULL;

   if (keep_flag)
      traversal_nodes.push_back(v);
   else
      traversal_nodes.push_back(t.root());
   vid_array= new vid[t.get_size()];

   while (!traversal_nodes.empty())
   {
      csnode= traversal_nodes.front();
      traversal_nodes.pop_front();
      subtree_nodes.push_back(csnode);
      nvertices++;
      if ((keep_flag) || (csnode != v))
      {
         // if !keep_flag, the vertices belonging to subtree
         // rooted at node v have to be ignored
         Tree_tie::tie(it, end)= t.children(csnode);
         while (it < end)
         {
            if ((keep_flag) || (*it != v))
               traversal_nodes.push_back(*it++);
            else
               it++;
         }
      }
   }

   res= new Tree(0, 0); // (root, number of vertices)
   nvertices_init=res->get_size();
   while (!subtree_nodes.empty())
   {
      csnode= subtree_nodes.front();
      subtree_nodes.pop_front();
      vid_array[csnode]= res->add_vertex(t.get(csnode));
   }

   if (keep_flag)
      traversal_nodes.push_back(v);
   else
      traversal_nodes.push_back(t.root());

   while (!traversal_nodes.empty())
   {
      csnode= traversal_nodes.front();
      traversal_nodes.pop_front();
      if ((keep_flag) || (csnode != v))
      {
         // if !keep_flag, the vertex belonging to subtree
         // rooted at node v have to be ignored
         Tree_tie::tie(it, end)= t.children(csnode);
         while (it < end)
         {
            if ((keep_flag) || (*it != v))
            {
               traversal_nodes.push_back(*it);
               res->add_edge(vid_array[csnode], vid_array[*it++]);
            }
            else
               it++;
         }
      }
   }

   assert(nvertices == res->get_size()-nvertices_init);

   delete [] vid_array;

   return res;
}

/*--------------------------------------------------------------*
 *
 *  Default constructor of Generic_tree class
 *
 *--------------------------------------------------------------*/

template <typename Type>
Generic_tree<Type>::Generic_tree(key root,
                                 int n,
                                 const value& default_value)
 : tree_simple<Type>(root, n, default_value)
{}


/*--------------------------------------------------------------*
 *
 *  Copy constructor of Generic_tree class
 *
 *--------------------------------------------------------------*/

template <typename Type>
Generic_tree<Type>::Generic_tree(const Generic_tree& tree)
  :tree_simple<Type>(tree) {}

/*--------------------------------------------------------------*
 *
 *  Constructor of Generic_tree class using a structure and
 *  a default value
 *
 *--------------------------------------------------------------*/

template <typename Type>
Generic_tree<Type>::Generic_tree(Unlabelled_tree& utree, const value& default_value)
{
  typedef typename tree_simple<char>::vertex_iterator char_vertex_iterator;
  char_vertex_iterator it, end;

  this->_root= utree._root;

  Tree_tie::tie(it,end)= utree.vertices();
  while( it != end )
    {
      this->add_vertex(default_value);
      it++;
    }

  this->_children= utree._children;
  this->_parent= utree._parent;
}

/*--------------------------------------------------------------*
 *
 *  Destructor of Generic_tree class
 *
 *--------------------------------------------------------------*/

template <typename Type>
Generic_tree<Type>::~Generic_tree()
 {}


/*--------------------------------------------------------------*
 *
 *  Assignement operator of Generic_tree class
 *
 *--------------------------------------------------------------*/

template <typename Type>
Generic_tree<Type>& Generic_tree<Type>::operator=(const Generic_tree<Type> &tree)
{
   copy(tree);
   return *this;
}

/*--------------------------------------------------------------*
 *
 *  Copy of the tree structure
 *
 *--------------------------------------------------------------*/

template <typename Type>
Unlabelled_tree* Generic_tree<Type>::get_structure() // const
{
  Unlabelled_tree *utree= new Unlabelled_tree(*this);
  return utree;
}

/*--------------------------------------------------------------*
 *
 *  Setting of the tree structure
 *
 *--------------------------------------------------------------*/

template <typename Type>
void Generic_tree<Type>::set_structure(Unlabelled_tree& utree,
                                       const value& default_value)
{
  typedef typename tree_simple<char>::vertex_iterator char_vertex_iterator;
  // typedef typename tree_simple<Type>::vertex_iterator vertex_iterator;
  char_vertex_iterator ch_it, ch_end;
  vertex_iterator it, end;

  this->_size= utree.get_size();
  this->_vertices.resize(this->_size);
  this->_root= utree._root;
  this->_children= utree._children;
  this->_parent= utree._parent;

  Tree_tie::tie(ch_it, ch_end)= utree.vertices();
  Tree_tie::tie(it, end)= this->vertices();
  while (ch_it != ch_end)
  {
      this->put(*it++, default_value);
      ch_it++;
  }
}

/*--------------------------------------------------------------*
 *
 *  Number of vertices
 *
 *--------------------------------------------------------------*/

template <typename Type>
unsigned int Generic_tree<Type>::get_size() const
{ return this->_size; }

/*--------------------------------------------------------------*
 *
 *  Tree depth
 *
 *--------------------------------------------------------------*/

template <typename Type>
unsigned int Generic_tree<Type>::get_depth()
{
    typedef typename tree_traits<Generic_tree<Type> >::vertex_descriptor vid;
    typedef typename generic_visitor<Generic_tree<Type> >::vertex_array vertex_array;

    generic_visitor<Generic_tree<Type> > visitor;
    // std::deque<vid> _nodes;
    // vid v;
    vertex_array va= visitor.get_breadthorder(*this);

    return get_depth(va.back());
}

/*--------------------------------------------------------------*
 *
 *  Depth of a vertex
 *
 *--------------------------------------------------------------*/

template <typename Type>
unsigned int Generic_tree<Type>::get_depth(key v)
{
    if(v == this->root())
         return 0;
    else
         return(get_depth(this->parent(v))+1);
}

/*--------------------------------------------------------------*
 *
 *  Number of children of a vertex
 *
 *--------------------------------------------------------------*/

template <typename Type>
unsigned int Generic_tree<Type>::get_nb_children(key v)
{ return(this->_children[v].size()); }

/*--------------------------------------------------------------*
 *
 *  Indicate whether the given couple defines an existing edge
 *
 *--------------------------------------------------------------*/

template <typename Type>
bool Generic_tree<Type>::is_edge(key parent, key child)
{
   // typename tree_traits<Tree>::
   children_iterator it, end;
   bool res= false;

   Tree_tie::tie(it, end)= this->children(parent);
   while (it < end)
      if (*it++ == child)
         res= true;

   return res;
}

/*--------------------------------------------------------------*
 *
 *  Display
 *
 *--------------------------------------------------------------*/

template <typename Type>
void Generic_tree<Type>::display(ostream& os,
                                 key v)
{
  bool lineskip= 0;

  display_skip(os, v, "", lineskip);
}

template <typename Type>
void Generic_tree<Type>::display_skip(ostream& os,
                                      key v,
                                      std::string tabulation,
                                      bool& lineskip)
{
   std::string tab= "|-";

   os << this->get(v) << endl;
   lineskip= 0;
   children_iterator i, end;
   Tree_tie::tie(i, end) = this->children(v);
   if (i != end)
   {
      int s= tabulation.size();
      if (s > 0)
      {
         tabulation.erase(tabulation.end()-1,tabulation.end());
         tabulation+= " ";
      }
      tabulation+= tab;
      os << tabulation;

      if (i == end-1)
      // last child of node v is being written
      {
         tabulation.erase(tabulation.end()-2,tabulation.end());
         tabulation+= "  ";
      }

      display_skip(os, *i++, tabulation, lineskip);

      while(i != end)
      {
         os << tabulation;
         if (i == end-1)
         // last child of node v is being written
         {
            tabulation.erase(tabulation.end()-2,tabulation.end());
            tabulation+= "  ";
         }
         display_skip(os, *i++, tabulation, lineskip);
      }
      if ((tabulation.size() >= 3) && !(lineskip))
      {
         lineskip= 1;
         tabulation.erase(tabulation.end()-3,tabulation.end());
         tabulation+= "   ";
         os << tabulation << endl;
      }
   }
   else // i == end
   {
      if( tabulation.size() >= 2 )
      {
       // lineskip= 1;
         std::string::iterator j=tabulation.end(),i=j-2;
         tabulation.erase(i,j);
      }
   }
}


/*--------------------------------------------------------------*
 *
 *  Copy operator of Generic_tree class
 *
 *--------------------------------------------------------------*/

template <typename Type>
void Generic_tree<Type>::copy(const Generic_tree<Type> &tree)
{
  this->_root= tree._root;
  this->_children= tree._children;
  this->_parent= tree._parent;
  this->_vertices= tree._vertices;

  this->_size= tree._size;
  this->_capacity= tree._capacity;
}

/*--------------------------------------------------------------*
 *
 *  Left (bit) shift operator of Generic_tree class
 *
 *--------------------------------------------------------------*/

template <typename Type>
std::ostream& operator<<(std::ostream& os,
                         Generic_tree<Type>& tree)
{
   tree.display(os, tree.root());
   return os;
}

/*--------------------------------------------------------------*
 *
 *  Constructor of Unlabelled_tree class from Generic_tree
 *
 *--------------------------------------------------------------*/

template<typename Type>
Unlabelled_tree::Unlabelled_tree(Generic_tree<Type>& gtree)
{
   typedef typename tree_simple<Type>::vertex_iterator generic_vertex_iterator;
   generic_vertex_iterator it, end;

   _root= gtree._root;

   Tree_tie::tie(it,end)= gtree.vertices();
   while( it < end )
   {
      add_vertex(C_DEFAULT_CHAR);
      it++;
   }

   _children= gtree._children;
   _parent= gtree._parent;
}

/*--------------------------------------------------------------*
 *
 *  Default constructor of Random_unlabelled_tree_generator class
 *
 *--------------------------------------------------------------*/

template <class Tree>
Random_unlabelled_tree_generator<Tree>::Random_unlabelled_tree_generator(const Distribution& distrib,
                                                                         int max_size,
                                                                         int max_depth) :
    _curr_width(0),
    _nb_children(0),
    _max_size(max_size),
    _max_depth(max_depth),
    _leaves()
{ _distribution= new Distribution(distrib); }


/*--------------------------------------------------------------*
 *
 *  Copy constructors of Random_unlabelled_tree_generator class
 *
 *--------------------------------------------------------------*/

template <class Tree>
Random_unlabelled_tree_generator<Tree>::Random_unlabelled_tree_generator(const Random_unlabelled_tree_generator& rutg)
{
  _curr_width= rutg._curr_width;
  _nb_children= rutg._nb_children;
  _max_size= rutg._nb_children;
  _max_depth= rutg._max_depth;
  _distribution= new Distribution(*(rutg._distribution));
}

/*--------------------------------------------------------------*
 *
 *  Destructor of Random_unlabelled_tree_generator class
 *
 *--------------------------------------------------------------*/

template <class Tree>
Random_unlabelled_tree_generator<Tree>::~Random_unlabelled_tree_generator()
{
   if(_distribution)
      delete _distribution;
}

/*--------------------------------------------------------------*
 *
 * Generation of the random structure
 *
 *--------------------------------------------------------------*/


template <class Tree>
Unlabelled_tree* Random_unlabelled_tree_generator<Tree>::run()
{
    Unlabelled_tree *_res_tree= new Unlabelled_tree;

    _curr_width= 1;
    // vid root= _res_tree->add_vertex(C_DEFAULT_CHAR);
    // _res_tree->root()= root;
    _leaves.push_front(_res_tree->root());
    while(( _max_size > 0 ) && (_max_depth > 0) && (!_leaves.empty()))
        random_children(_leaves.back(), _res_tree);
    // processing of next leaf

    return _res_tree;
}


/*--------------------------------------------------------------*
 *
 * Generation of the children of a given parent
 *
 *--------------------------------------------------------------*/

template <class Tree>
void Random_unlabelled_tree_generator<Tree>::random_children(vid parent,
                                                             Unlabelled_tree* const _tree)
{

    if (!_curr_width)
      {
        _max_depth-= 1;
        _curr_width= _nb_children;
        _nb_children= 0;
      }

    if (_max_depth)
    {
         _leaves.pop_back();
         _curr_width-= 1;
         // current vertex need not be processed any more
         int n= min(_distribution->simulation(), _max_size-_tree->size());
         for( int i= 0; i < n; i++ )
           {
               vid v= _tree->add_vertex(C_DEFAULT_CHAR);
               _tree->add_edge(parent, v);
               _leaves.push_front(v);
           }
         _max_size-= n;
         _nb_children+= n;

    }
}


#endif
