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
 *       $Id: generic_typed_edge_tree.hpp 3186 2007-05-25 15:10:30Z dufourko $
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

#ifndef GENERIC_TYPED_EDGE_TREE_TCC
#define GENERIC_TYPED_EDGE_TREE_TCC

/*--------------------------------------------------------------*
 *
 *  Selection of a subtree of Typed_edge_tree
 *  using the subtree root and a flag on keeping
 *  or pruning the selected subtree
 *
 *--------------------------------------------------------------*/

template <class Tree>
Tree* select_typed_edge_subtree(Tree& t,
                                typename tree_traits<Tree>::vertex_descriptor v,
                                bool keep_flag)
{
   typedef typename tree_traits<Tree>::vertex_descriptor key;
   typename tree_traits<Tree>::children_iterator it, end;
   // typedef typename Typed_edge_tree<Type>::key key;
   // typename Typed_edge_tree<Type>::children_iterator it, end;
   std::deque<key> traversal_nodes,
                   subtree_nodes;
    // list of nodes used for tree traversal and for building result tree
   bool type;
   key csnode, // current source node id
       *vid_array;
   int nvertices= 0; // number of subtree vertices
   int nvertices_init= 0; // number of initial vertices as set by constructor
   Tree* res= NULL;

   if (keep_flag)
      traversal_nodes.push_back(v);
   else
      traversal_nodes.push_back(t.root());
   vid_array= new key[t.size()];

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
   nvertices_init=res->size();
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
               // type= t.edge_type(vid_array[csnode], vid_array[*it]);
               type= t.edge_type(csnode, *it);
               res->add_edge(vid_array[csnode], vid_array[*it++], type);
            }
            else
               it++;
         }
      }
   }

   assert(nvertices == res->size()-nvertices_init);

   delete [] vid_array;

   return res;
}

/*--------------------------------------------------------------*
 *
 *  Default constructor of Typed_edge_tree class
 *
 *--------------------------------------------------------------*/

template <typename Type>
Typed_edge_tree<Type>::Typed_edge_tree(key root,
                                       int n,
                                       const value& default_value)
 : Generic_tree<Type>(root, n, default_value)
 , edge_types()
{
   if (n > 0)
      edge_types.reserve(n-1);
}


/*--------------------------------------------------------------*
 *
 *  Copy constructor of Typed_edge_tree class
 *
 *--------------------------------------------------------------*/

template <typename Type>
Typed_edge_tree<Type>::Typed_edge_tree(const Typed_edge_tree& tree)
 : Generic_tree<Type>(tree)

{ edge_types= tree.edge_types; }

/*--------------------------------------------------------------*
 *
 *  Constructor of Typed_edge_tree class using a structure and
 *  a default value
 *
 *--------------------------------------------------------------*/

template <typename Type>
Typed_edge_tree<Type>::Typed_edge_tree(Unlabelled_typed_edge_tree& utree,
                                       const value& default_value)
 : Generic_tree<Type>(utree, default_value)
 , edge_types()
{
   // const unsigned int size= utree.get_size();

   edge_types= utree.edge_types;
}


/*--------------------------------------------------------------*
 *
 *  Destructor of Typed_edge_tree class
 *
 *--------------------------------------------------------------*/

template <typename Type>
Typed_edge_tree<Type>::~Typed_edge_tree()
{}

/*--------------------------------------------------------------*
 *
 *  Assignement operator of Typed_edge_tree class
 *
 *--------------------------------------------------------------*/

template <typename Type>
Typed_edge_tree<Type>& Typed_edge_tree<Type>::operator=(const Typed_edge_tree<Type> &tree)
{
   copy(tree);
   edge_types= tree.edge_types;
   return *this;
}

/*--------------------------------------------------------------*
 *
 *  Return the tree structure from a typed-edge tree
 *
 *--------------------------------------------------------------*/
template <typename Type>
Unlabelled_typed_edge_tree* Typed_edge_tree<Type>::get_structure()
{
   Unlabelled_typed_edge_tree *res= new Unlabelled_typed_edge_tree(*this);
   // res->edge_types= edge_types;

   return res;
}

/*--------------------------------------------------------------*
 *
 *  Set the tree structure from a typed-edge tree
 *
 *--------------------------------------------------------------*/

template <typename Type>
void Typed_edge_tree<Type>::set_structure(Unlabelled_typed_edge_tree& utree,
                                          const value& default_value)
{
   Generic_tree<Type>::set_structure(utree, default_value);
   edge_types= utree.edge_types;
}

/*--------------------------------------------------------------*
 *
 * Add a vertex and return its identifier
 *
 *--------------------------------------------------------------*/

template <typename Type>
typename Typed_edge_tree<Type>::key Typed_edge_tree<Type>::add_vertex(const value& data)
{
   edge_types.push_back(false);
   return Generic_tree<Type>::add_vertex(data);
}

/*--------------------------------------------------------------*
 *
 *  Add an edge to a Typed_edge_tree
 *
 *--------------------------------------------------------------*/

template <typename Type>
bool Typed_edge_tree<Type>::add_edge(key parent,
                                     key child,
                                     bool type)
{
   bool res= true;
   if (!this->is_edge(parent, child))
      res= Generic_tree<Type>::add_edge(parent, child);
   edge_types[child]= type;
   return res;
}


/*--------------------------------------------------------------*
 *
 *  Return the type of a given edge
 *
 *--------------------------------------------------------------*/

template <typename Type>
bool Typed_edge_tree<Type>::edge_type(key parent, key child)
{
   assert(this->is_edge(parent, child));
   return edge_types[child];
}

/*--------------------------------------------------------------*
 *
 *  Set the type of a given edge and return true
 *  iif the edge exists
 *
 *--------------------------------------------------------------*/

template <typename Type>
bool Typed_edge_tree<Type>::set_edge_type(key parent, key child, bool type)
{
   const bool res = this->is_edge(parent, child);

   if (res)
      edge_types[child] = type;

   return res;
}

/*--------------------------------------------------------------*
 *
 *  Return the branching order (depth) of a typed-edge tree
 *
 *--------------------------------------------------------------*/

template <typename Type>
unsigned int Typed_edge_tree<Type>::get_branching_order()
{
   typedef typename tree_traits<Generic_tree<Type> >::vertex_descriptor vid;
   typedef typename generic_visitor<tree_type>::vertex_array vertex_array;

   unsigned int res= 1; //, vorder;
   vid node;
   generic_visitor<Generic_tree<Type> > visitor;
   vertex_array va= visitor.get_breadthorder(*this);

   for(node= 0; node < (vid)va.size(); node++)
      if (this->get_nb_children(va[node]) == 0)
         res= max(res, get_branching_order(va[node]));


   return res;
}

/*--------------------------------------------------------------*
 *
 * Return the branching order (depth) of a given vertex
 *
 *--------------------------------------------------------------*/

template <typename Type>
unsigned int Typed_edge_tree<Type>::get_branching_order(key v)
{
   unsigned int res= 1; //, pdepth;
   key p;

   if(v != this->root())
   {
      p= this->parent(v);
      if (!this->edge_type(p, v))
         // branching
         res= get_branching_order(this->parent(v))+1;
      else
         res= get_branching_order(this->parent(v));
   }
   return res;
}

/*--------------------------------------------------------------*
 *
 *  Display
 *
 *--------------------------------------------------------------*/
/*
template <typename Type>
void Typed_edge_tree<Type>::display(ostream& os,
                                 key v)
{
  bool lineskip= 0;

  display_skip(os, v, "", lineskip);
}

template <typename Type>
void Typed_edge_tree<Type>::display_skip(ostream& os,
                                      key v,
                                      std::string tabulation,
                                      bool& lineskip)
{
   std::string tab= "|-";

   os << get(v) << endl;
   lineskip= 0;
   children_iterator i, end;
   Tree_tie::tie(i, end) = children(v);
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
*/

/*--------------------------------------------------------------*
 *
 *  Copy operator of Typed_edge_tree class
 *
 *--------------------------------------------------------------*/

template <typename Type>
void Typed_edge_tree<Type>::copy(const Typed_edge_tree<Type> &tree)
{
   Generic_tree<Type>::copy(tree);
   edge_types= tree.edge_types;
}

/*--------------------------------------------------------------*
 *
 *  Left (bit) shift operator of Typed_edge_tree class
 *
 *--------------------------------------------------------------*/

template <typename Type>
std::ostream& operator<<(std::ostream& os,
                         Typed_edge_tree<Type>& tree)
{
   tree.display(os, tree.root());
   return os;
}

/*--------------------------------------------------------------*
 *
 *  Constructor of Unlabelled_typed_edge_tree class from Generic_tree
 *
 *--------------------------------------------------------------*/

template<typename Type>
Unlabelled_typed_edge_tree::Unlabelled_typed_edge_tree(Generic_tree<Type>& gtree)
 : Typed_edge_tree<char>(gtree.Generic_tree<Type>::_root, 1)
 , Unlabelled_tree()
{
   typedef typename Generic_tree<Type>::tree_type tree_type;
   typedef typename tree_type::vertex_descriptor vid;
   typedef generic_visitor<tree_type> visitor;
   typedef typename visitor::vertex_array vertex_array;
   int u;
   vid cnode, *dest_vids= NULL;
   visitor *v= NULL;
   vertex_array va;

   Typed_edge_tree<char>::_size= Unlabelled_tree::_size;
   Typed_edge_tree<char>::_root= Unlabelled_tree::_root;
   assert((unsigned int)this->Typed_edge_tree<char>::size() == this->get_size());
   dest_vids= new vid[gtree.get_size()];
   // add root vertex
   dest_vids[gtree.root()]= Unlabelled_tree::root();

   v= new generic_visitor<tree_type>;
   traverse_tree(gtree.root(), gtree, *v);
   va= v->get_breadthorder(gtree);
   delete v;
   v= NULL;

   for(u= 0; u < (int)gtree.get_size(); u++)
   {
      cnode= va[u];
      if (cnode != this->root())
      {
         dest_vids[cnode]= Unlabelled_tree::add_vertex(C_DEFAULT_CHAR);
         Unlabelled_tree::add_edge(gtree.parent(dest_vids[cnode]),
                                                dest_vids[cnode]);
      }
   }

   edge_types.resize(gtree.size(), false);
   Typed_edge_tree<char>::_size= Unlabelled_tree::_size;
   Typed_edge_tree<char>::_root= Unlabelled_tree::_root;
   Typed_edge_tree<char>::_children= Unlabelled_tree::_children;
   Typed_edge_tree<char>::_parent= Unlabelled_tree::_parent;
   Typed_edge_tree<char>::_vertices= Unlabelled_tree::_vertices;

   delete [] dest_vids;
   dest_vids= NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructor of Unlabelled_typed_edge_tree class from Typed_edge_tree
 *
 *--------------------------------------------------------------*/

template<typename Type>
Unlabelled_typed_edge_tree::Unlabelled_typed_edge_tree(Typed_edge_tree<Type>& gtree)
 : Typed_edge_tree<char>(gtree.Generic_tree<Type>::_root, 1)
 , Unlabelled_tree()
{
   typedef typename Generic_tree<Type>::tree_type tree_type;
   typedef typename tree_type::vertex_descriptor vid;
   typedef generic_visitor<tree_type> visitor;
   typedef typename visitor::vertex_array vertex_array;
   int u;
   vid cnode, *dest_vids= NULL;
   visitor *v= NULL;
   vertex_array va;

   Typed_edge_tree<char>::_size= Unlabelled_tree::_size;
   Typed_edge_tree<char>::_root= Unlabelled_tree::_root;
   assert((unsigned int)this->Typed_edge_tree<char>::size() == this->get_size());
   dest_vids= new vid[gtree.get_size()];
   // add root vertex
   dest_vids[gtree.root()]= Unlabelled_tree::root();

   v= new generic_visitor<tree_type>;
   traverse_tree(gtree.root(), gtree, *v);
   va= v->get_breadthorder(gtree);
   delete v;
   v= NULL;

   for(u= 0; u < (int)gtree.get_size(); u++)
   {
      cnode= va[u];
      if (cnode != this->root())
      {
         dest_vids[cnode]= Unlabelled_tree::add_vertex(C_DEFAULT_CHAR);
         Unlabelled_tree::add_edge(gtree.parent(dest_vids[cnode]),
                                                dest_vids[cnode]);
      }
   }

   edge_types= gtree.Typed_edge_tree<Type>::edge_types;
   Typed_edge_tree<char>::_size= Unlabelled_tree::_size;
   Typed_edge_tree<char>::_root= Unlabelled_tree::_root;
   Typed_edge_tree<char>::_children= Unlabelled_tree::_children;
   Typed_edge_tree<char>::_parent= Unlabelled_tree::_parent;
   Typed_edge_tree<char>::_vertices= Unlabelled_tree::_vertices;

   delete [] dest_vids;
   dest_vids= NULL;
}

#endif
