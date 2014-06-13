// Includes ====================================================================
#include "tree/basic_visitors.h"
#include "tree/tree_simple.h"
#include "tree/tree_traits.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/distribution.h"
#include "stat_tool/vectors.h"

#include "sequence_analysis/sequences.h"
#include "sequence_analysis/sequence_label.h"

#include "tree_statistic/int_fl_containers.h"
#include "tree_statistic/tree_labels.h"
#include "tree_statistic/generic_typed_edge_tree.h"
#include "tree_statistic/typed_edge_trees.h"

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
// definition of boost::python::len
#include <boost/python/make_constructor.hpp>
// definition of boost::python::make_constructor
#include "../errors.h"

// Using =======================================================================
using namespace boost::python;
using namespace Stat_trees;

// Declarations ================================================================
namespace  {
//
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Default_tree_add_vertex_overloads_0_1, add_vertex, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Unlabelled_typed_edge_tree_add_vertex_overloads_0_1, add_vertex, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Unlabelled_typed_edge_tree_simulation_overloads_1_3, simulation, 1, 3)
//

template<int num> struct UniqueInt { int v; enum { value=num };
UniqueInt(int _v) : v(_v) { } operator int() const { return v; } };

/*************************************************************
 *
 *  Exporting constants as constant functions
 */

int NB_TREES_() { return NB_TREES; }
int I_DEFAULT_TREE_SIZE_() { return I_DEFAULT_TREE_SIZE; }
int I_DEFAULT_TREE_DEPTH_() { return I_DEFAULT_TREE_DEPTH; }

/*************************************************************
 *
 *  Generic wrappers:
 */

template <class Tree>
typename tree_traits<Tree>::vertex_descriptor
Tree_wrapper_Parent(Tree& tree, typename tree_traits<Tree>::vertex_descriptor v)
{
   // int parent;
   ostringstream error_message;
   bool status= true;

   if ((v < 0) || (v > (int)tree.get_size()-1))
   {
      status= false;
      error_message << STAT_TREES_error[TREESTATR_VERTEX_ID] << ": "
                    << v << endl;
   }
   if (v == tree.root())
   {
      status= false;
      error_message << "Parent of root node is undefined" << endl;
   }
   if (!status)
   {
      PyErr_SetString(PyExc_IndexError, (error_message.str()).c_str());
      throw_error_already_set();
   }
   // else
   return tree.parent(v);
}

template <class Tree> class Tree_iterator
{
private :

   typedef typename tree_traits<Tree>::vertex_iterator vertex_iterator;
   typedef typename tree_traits<Tree>::vertex_descriptor key;
   vertex_iterator it, end;

public :

   Tree_iterator(Tree& tree) { Tree_tie::tie(it, end)= tree.vertices(); }
   ~Tree_iterator() {}
   Tree_iterator iter() { return *this; }
   key next()
   {
      if (it == end)
      {
         PyErr_SetString(PyExc_StopIteration, "No more vertices.");
         boost::python::throw_error_already_set();
      }
      return *it++;
   }
};

template <class Tree>
Tree_iterator<Tree> Tree_wrapper_vertex_iterator(Tree& tree)
{
   Tree_iterator<Tree> it(tree);
   return it;
}

template <class Tree> class Tree_iterator_preorder
{
private :

   typedef typename generic_visitor<Tree>::vertex_array vertex_array;
   typedef typename tree_traits<Tree>::vertex_descriptor key;
   generic_visitor<Tree> visitor;
   vertex_array va;
   int current_vertex, end;

public :

   Tree_iterator_preorder(Tree& tree)
   {
      // visitor= new generic_visitor<Tree>;
      traverse_tree(tree.root(), tree, visitor);
      va= visitor.get_preorder(tree);
      current_vertex= 0;
      end= tree.size();
   }
   ~Tree_iterator_preorder() {}
   Tree_iterator_preorder iter() { return *this; }
   key next()
   {
      if (current_vertex == end)
      {
         PyErr_SetString(PyExc_StopIteration, "No more vertices.");
         boost::python::throw_error_already_set();
      }
      return va[current_vertex++];
   }
};

template <class Tree>
Tree_iterator_preorder<Tree> Tree_wrapper_vertex_iterator_preorder(Tree& tree)
{
   Tree_iterator_preorder<Tree> it(tree);
   return it;
}

template <class Tree> class Tree_iterator_inorder
{
private :

   typedef typename generic_visitor<Tree>::vertex_array vertex_array;
   typedef typename tree_traits<Tree>::vertex_descriptor key;
   generic_visitor<Tree> visitor;
   vertex_array va;
   int current_vertex, end;

public :

   Tree_iterator_inorder(Tree& tree)
   {
      traverse_tree(tree.root(), tree, visitor);
      va= visitor.get_inorder(tree);
      current_vertex= 0;
      end= tree.size();
   }
   ~Tree_iterator_inorder() {}
   Tree_iterator_inorder iter() { return *this; }
   key next()
   {
      if (current_vertex == end)
      {
         PyErr_SetString(PyExc_StopIteration, "No more vertices.");
         boost::python::throw_error_already_set();
      }
      return va[current_vertex++];
   }
};

template <class Tree>
Tree_iterator_inorder<Tree> Tree_wrapper_vertex_iterator_inorder(Tree& tree)
{
   Tree_iterator_inorder<Tree> it(tree);
   return it;
}

template <class Tree> class Tree_iterator_postorder
{
private :

   typedef typename generic_visitor<Tree>::vertex_array vertex_array;
   typedef typename tree_traits<Tree>::vertex_descriptor key;
   generic_visitor<Tree> visitor;
   vertex_array va;
   int current_vertex, end;

public :

   Tree_iterator_postorder(Tree& tree)
   {
      traverse_tree(tree.root(), tree, visitor);
      va= visitor.get_postorder(tree);
      current_vertex= 0;
      end= tree.size();
   }
   ~Tree_iterator_postorder() {}
   Tree_iterator_postorder iter() { return *this; }
   key next()
   {
      if (current_vertex == end)
      {
         PyErr_SetString(PyExc_StopIteration, "No more vertices.");
         boost::python::throw_error_already_set();
      }
      return va[current_vertex++];
   }
};

template <class Tree>
Tree_iterator_postorder<Tree> Tree_wrapper_vertex_iterator_postorder(Tree& tree)
{
   Tree_iterator_postorder<Tree> it(tree);
   return it;
}

template <class Tree> class Tree_iterator_breadthorder
{
private :

   typedef typename generic_visitor<Tree>::vertex_array vertex_array;
   typedef typename tree_traits<Tree>::vertex_descriptor key;
   generic_visitor<Tree> visitor;
   vertex_array va;
   int current_vertex, end;

public :

   Tree_iterator_breadthorder(Tree& tree)
   {
      va= visitor.get_breadthorder(tree);
      current_vertex= 0;
      end= tree.size();
   }
   ~Tree_iterator_breadthorder() {}
   Tree_iterator_breadthorder iter() { return *this; }
   key next()
   {
      if (current_vertex == end)
      {
         PyErr_SetString(PyExc_StopIteration, "No more vertices.");
         boost::python::throw_error_already_set();
      }
      return va[current_vertex++];
   }
};

template <class Tree>
Tree_iterator_breadthorder<Tree> Tree_wrapper_vertex_iterator_breadthorder(Tree& tree)
{
   Tree_iterator_breadthorder<Tree> it(tree);
   return it;
}

template <class Tree> class Tree_iterator_leavesfirstorder
{
private :

   typedef typename generic_visitor<Tree>::vertex_array vertex_array;
   typedef typename tree_traits<Tree>::vertex_descriptor key;
   generic_visitor<Tree> visitor;
   vertex_array va;
   int current_vertex, end;

public :

   Tree_iterator_leavesfirstorder(Tree& tree)
   {
      std::vector<int> depth;

      va= visitor.get_leavesfirstorder(tree, depth);
      current_vertex= 0;
      end= tree.size();
   }
   ~Tree_iterator_leavesfirstorder() {}
   Tree_iterator_leavesfirstorder iter() { return *this; }
   key next()
   {
      if (current_vertex == end)
      {
         PyErr_SetString(PyExc_StopIteration, "No more vertices.");
         boost::python::throw_error_already_set();
      }
      return va[current_vertex++];
   }
};

template <class Tree>
Tree_iterator_leavesfirstorder<Tree> Tree_wrapper_vertex_iterator_leavesfirstorder(Tree& tree)
{
   Tree_iterator_leavesfirstorder<Tree> it(tree);
   return it;
}

template <class Tree> class Children_iterator
{
private :

   typedef typename tree_traits<Tree>::children_iterator children_iterator;
   typedef typename tree_traits<Tree>::vertex_descriptor key;
   children_iterator it, end;

public :

   Children_iterator(Tree& tree, key v)
   {
      ostringstream error_message;
      bool status= true;

      if ((v < 0) || (v > (int)tree.get_size()-1))
      {
         status= false;
         error_message << STAT_TREES_error[TREESTATR_VERTEX_ID] << ": "
                       << v << endl;
      }
      if (!status)
      {
         PyErr_SetString(PyExc_IndexError, (error_message.str()).c_str());
         throw_error_already_set();
      }
      else
          Tree_tie::tie(it, end)= tree.children(v);
   }
   ~Children_iterator() {}
   Children_iterator iter() { return *this; }
   key next()
   {
      if (it == end)
      {
         PyErr_SetString(PyExc_StopIteration, "No more children.");
         boost::python::throw_error_already_set();
      }
      return *it++;
   }
};

template <class Tree>
Children_iterator<Tree> Tree_wrapper_children_iterator(Tree& tree,
                                                       typename tree_traits<Tree>::vertex_descriptor v)
{
   Children_iterator<Tree> it(tree, v);
   return it;
}


template <class Tree> int Generic_tree_wrapper_Depth0(Tree& tree)
{ return tree.get_depth(); }

template <class Tree> int Generic_tree_wrapper_Depth1(Tree& tree,
                                                      typename Tree::key v)
{ return tree.get_depth(v); }

template <class Tree> int Generic_tree_wrapper_Order0(Tree& tree)
{ return tree.get_branching_order(); }

template <class Tree> int Generic_tree_wrapper_Order1(Tree& tree,
                                                      typename Tree::key v)
{ return tree.get_branching_order(v); }

template <class Tree>
std::string Generic_tree_wrapper_Display0(Tree& tree)
{
   std::stringstream s;
   std::string res;
   if (tree.get_size() > 0)
   {
      tree.display(s, tree.root());
      res= s.str();
   }
   else
      res= "(empty tree)";
   return res;
}

template <class Tree>
std::string Generic_tree_wrapper_Display1(Tree& tree,
                                          typename Tree::key v)
{
   std::stringstream s;
   std::string res;
   ostringstream error_message;
   bool status= true;

   if ((v < 0) || (v > (int)tree.get_size()-1))
   {
      status= false;
      error_message << STAT_TREES_error[TREESTATR_VERTEX_ID] << ": "
                    << v << endl;
   }
   if (!status)
   {
      PyErr_SetString(PyExc_IndexError, (error_message.str()).c_str());
      throw_error_already_set();
   }

   if (tree.get_size() > 0)
   {
      tree.display(s, v);
      res= s.str();
   }
   else
      res= "(empty tree)";
   return res;
}

/*************************************************************
 *
 *  Wrappers for Python class Tree:
 */


Default_tree* Tree_wrapper_init1(Default_tree::key root,
                                 int n, object default_value)
{
   Int_fl_container cvalue;
   boost::python::list values, types;
   int nb_integral, nb_float, index, type,
       int_index, float_index; // length,
   object current_object;
   ostringstream error_message;
   // bool status= true;
   Default_tree* tree;

   values= extract<boost::python::list>(default_value.attr("Values")());
   types= extract<boost::python::list>(default_value.attr("Types")());
   nb_integral= extract<int>(default_value.attr("NbInt")());
   nb_float= extract<int>(default_value.attr("NbFloat")());
   int_index= 0;
   float_index= 0;

   cvalue.reset(nb_integral, nb_float);
   for(index=0; index < nb_integral+nb_float; index++)
   {
      type= extract<int>(types[index]);
      if (type==REAL_VALUE)
         cvalue.Double(float_index++)= extract<double>(values[index]);
      else
      {
         if ((type==INT_VALUE) || (type==STATE) || (type==NB_INTERNODE))
            cvalue.Int(int_index++)= extract<int>(values[index]);
         else
         {
            PyErr_SetString(PyExc_IndexError, "unknown variable type");
            throw_error_already_set();
         }
      }
   }
   tree= new Typed_edge_int_fl_tree<Int_fl_container>(root, n, cvalue);
   return tree;
}

tree_traits<Default_tree>::vertex_descriptor
Tree_wrapper_AddVertex(Default_tree& tree, object o)
{
   Int_fl_container cvalue(0, 0);
   boost::python::list values, types;
   int nb_integral, nb_float, index, type, //length,
       int_index, float_index, tree_nb_integral, tree_nb_float;
   object current_object;
   ostringstream error_message;
   bool status= true;

   tree_nb_integral= tree.get_nb_int();
   tree_nb_float= tree.get_nb_float();

   values= extract<boost::python::list>(o.attr("Values")());
   types= extract<boost::python::list>(o.attr("Types")());
   nb_integral= extract<int>(o.attr("NbInt")());
   nb_float= extract<int>(o.attr("NbFloat")());
   int_index= 0;
   float_index= 0;

   if (nb_integral != tree_nb_integral)
   {
      status= false;
      error_message << "incompatible number of integral variables: "
                    << tree_nb_integral << " for tree and "
                    << nb_integral << " for TreeValue object." << endl;
   }

   if (nb_float != tree_nb_float)
   {
      status= false;
      error_message << "incompatible number of floating variables: "
                    << tree_nb_float << " for tree and "
                    << nb_float << " for TreeValue object." << endl;
   }

   if (!status)
   {
      PyErr_SetString(PyExc_IndexError, (error_message.str()).c_str());
      throw_error_already_set();
   }
   else
   {
      cvalue.reset(nb_integral, nb_float);
      for(index=0; index < nb_integral+nb_float; index++)
      {
         type= extract<int>(types[index]);
         if (type==REAL_VALUE)
            cvalue.Double(float_index++)= extract<double>(values[index]);
         else
         {
            if (type==INT_VALUE)
               cvalue.Int(int_index++)= extract<int>(values[index]);
            else
            {
               PyErr_SetString(PyExc_IndexError, "unknown variable type");
               throw_error_already_set();
            }
         }
      }
   }
   return tree.add_vertex(cvalue);
}

void Tree_wrapper_IidSimulation(Default_tree& tree, const boost::python::list& distribution_list)
{
   int length, index, nb_int= tree.get_nb_int();
   Default_tree::pt_Distribution_array distributions;
   ostringstream error_message;
   bool status= true;
   object current_object;

   length= boost::python::len(distribution_list);
   distributions= new Distribution*[nb_int];
   if (length != nb_int)
   {
      status= false;
      error_message << "bad number of distributions; should be " << nb_int << endl;
   }

   for(index= 0; index < length; index++)
   {
      current_object= distribution_list[index];
      distributions[index]= NULL;

      extract<Distribution> x(current_object);
      if (x.check())
         distributions[index]= new Distribution(x());
      else
      {
         status= false;
         error_message << "bad type for object number " << index
                       << " - should be a Distribution object" << endl;
      }
   }

   if (!status)
   {
      PyErr_SetString(PyExc_IndexError, (error_message.str()).c_str());
      throw_error_already_set();
   }
   else
      tree.iid_simulation(distributions);

   for(index= 0; index < length; index++)
      if (distributions[index] == NULL)
         delete distributions[index];

   delete [] distributions;
}

void Tree_wrapper_Put(Default_tree& tree, Default_tree::key v, object o)
{
   Int_fl_container cvalue(0, 0);
   boost::python::list values, types;
   int nb_integral, nb_float, index, type, //length,
       int_index, float_index, tree_nb_integral, tree_nb_float;
   object current_object;
   ostringstream error_message;
   bool status= true;

   tree_nb_integral= tree.get_nb_int();
   tree_nb_float= tree.get_nb_float();

   values= extract<boost::python::list>(o.attr("Values")());
   types= extract<boost::python::list>(o.attr("Types")());
   nb_integral= extract<int>(o.attr("NbInt")());
   nb_float= extract<int>(o.attr("NbFloat")());
   int_index= 0;
   float_index= 0;

   if (nb_integral != tree_nb_integral)
   {
      status= false;
      error_message << "incompatible number of integral variables: "
                    << tree_nb_integral << " for tree and "
                    << nb_integral << " for TreeValue object." << endl;
   }

   if (nb_float != tree_nb_float)
   {
      status= false;
      error_message << "incompatible number of floating variables: "
                    << tree_nb_float << " for tree and "
                    << nb_float << " for TreeValue object." << endl;
   }

   if ((v < 0) || (v > (int)tree.get_size()-1))
   {
      status= false;
      error_message << STAT_TREES_error[TREESTATR_VERTEX_ID] << ": "
                    << v << endl;
   }


   if (!status)
   {
      PyErr_SetString(PyExc_IndexError, (error_message.str()).c_str());
      throw_error_already_set();
   }
   else
   {
      cvalue.reset(nb_integral, nb_float);
      for(index=0; index < nb_integral+nb_float; index++)
      {
         type= extract<int>(types[index]);
         if (type==REAL_VALUE)
            cvalue.Double(float_index++)= extract<double>(values[index]);
         else
         {
            if (type==INT_VALUE)
               cvalue.Int(int_index++)= extract<int>(values[index]);
            else
            {
               PyErr_SetString(PyExc_IndexError, "unknown variable type");
               throw_error_already_set();
            }
         }
      }
      // cout << cvalue << endl;
      tree.put(v, cvalue);
   }
}

/*************************************************************
 *
 *  Wrappers for Python class CharTree:
 */

int CharTree_wrapper_Depth0(Unlabelled_tree& tree)
{ return tree.get_depth(); }

int CharTree_wrapper_Depth1(Unlabelled_tree& tree, Unlabelled_tree::key v)
{
   ostringstream error_message;
   bool status= true;

   if ((v < 0) || (v > (int)tree.get_size()-1))
   {
      status= false;
      error_message << STAT_TREES_error[TREESTATR_VERTEX_ID] << ": "
                    << v << endl;
   }
   if (!status)
   {
      PyErr_SetString(PyExc_IndexError, (error_message.str()).c_str());
      throw_error_already_set();
   }

   return tree.get_depth(v);
}

std::string CharTree_wrapper_Display0(Unlabelled_tree& tree)
{
   std::stringstream s;
   std::string res;
   if (tree.get_size() > 0)
   {
      tree.display(s, tree.root());
      res= s.str();
   }
   else
      res= "(empty tree)";
   return res;
}

std::string CharTree_wrapper_Display1(Unlabelled_tree& tree, Default_tree::key v)
{
   std::stringstream s;
   std::string res;
   ostringstream error_message;
   bool status= true;

   if ((v < 0) || (v > (int)tree.get_size()-1))
   {
      status= false;
      error_message << STAT_TREES_error[TREESTATR_VERTEX_ID] << ": "
                    << v << endl;
   }
   if (!status)
   {
      PyErr_SetString(PyExc_IndexError, (error_message.str()).c_str());
      throw_error_already_set();
   }

   if (tree.get_size() > 0)
   {
      tree.display(s, v);
      res= s.str();
   }
   else
      res= "(empty tree)";
   return res;
}

Unlabelled_typed_edge_tree::vertex_descriptor
CharTree_wrapper_AddVertex(Unlabelled_typed_edge_tree& tree)
{ return tree.Typed_edge_tree<char>::add_vertex(C_DEFAULT_CHAR); }

class MyClass : public PyObject
{
public :

   int i;
   MyClass() {}
   ~MyClass() {}
   std::string f() const { return "Hello, world !"; }

};
// std::string MyClass_wrapper_f(const MyClass& m)
// { return m.f(); }

class TreeValue
{
   friend std::ostream& operator<<(std::ostream& os, const TreeValue& v);

protected :

   int _nb_integral, _nb_float;
   boost::python::list values;
   std::vector<int> types;

public :

   TreeValue(boost::python::list lvalues_or_types)
   : _nb_integral(0)
   , _nb_float(0)
   , values(lvalues_or_types)
   , types()
   {
      int length, index;
      object current_object;
      ostringstream error_message;

      length= boost::python::len(lvalues_or_types);

      // equivalent to:
      // object check_length( values.attr("__len__")());
      // extract<int> x(check_length);

      for(index= 0; index < length; index++)
      {
         current_object= values[index];

         extract<int> x(current_object);
         if (x.check())
         {
            types.push_back(INT_VALUE);
            _nb_integral++;
         }
         else
         {
            extract<double> x(current_object);
            if (x.check())
            {
               types.push_back(REAL_VALUE);
               _nb_float++;
            }
            else
            {
               error_message << STAT_TREES_error[TREESTATR_VARIABLE_TYPE];
               PyErr_SetString(PyExc_IndexError, (error_message.str()).c_str());
               throw_error_already_set();
            }
         }
      }
   }
   ~TreeValue() {};
   TreeValue& operator=(const TreeValue& tree)
   {
      values= tree.values;
      types= tree.types;
      _nb_integral= tree._nb_integral;
      _nb_float= tree._nb_float;
      return *this;
   }
   const object& operator[](int index) const
   {
      const object& res= values[index];
      return res;
   }

   int get_type(int index) const
   { return types[index]; }

   boost::python::list get_values() const { return values; }
   boost::python::list get_types() const
   {
      boost::python::list res;
      unsigned int index;

      for(index= 0; index < types.size(); index++)
         res.append(types[index]);

      return res;
   }

   void setitem(int index, const object& value)
   {
      object o;
      StatError error;
      ostringstream error_message;

      error.init();
      if (index >= 0 && index < _nb_integral+_nb_float)
      {
         extract<int> x(value);
         if (x.check())
         {
            if ((types[index] == INT_VALUE) || (types[index] == STATE)
                || (types[index] == POSITION) || (types[index] == NB_INTERNODE))
               values[index]= x();
            else
            {
               if (types[index] == REAL_VALUE)
                  // we accept an int where a double is expected
                  values[index]= (double)x();
               else
               {
                  error_message << STAT_label[STATL_VARIABLE] << " " << index
                                << " " << STAT_TREES_error[TREESTATR_VARIABLE_TYPE]
                                << ": should be " << STAT_variable_word[INT_VALUE]
                                << " or " << STAT_word[STATE]
                                << " or " << STAT_word[POSITION]
                                << " or " << STAT_word[NB_INTERNODE];

                  PyErr_SetString(PyExc_IndexError, (error_message.str()).c_str());
                  throw_error_already_set();
                  // should use format_error instead
               }
            }
         }
         else
         {
            extract<double> x(value);
            if (x.check())
            {
               if (types[index] == REAL_VALUE)
                  values[index]= x();
            }
            else
            {
               ostringstream error_message;
               error_message << STAT_label[STATL_VARIABLE] << " " << index
                             << ": " << STAT_TREES_error[TREESTATR_VARIABLE_TYPE]
                             << ": should be " << STAT_variable_word[REAL_VALUE];
               error.update((error_message.str()).c_str());
               // PyErr_SetString(PyExc_IndexError, s.str());
               PyErr_SetString(PyExc_IndexError, (error_message.str()).c_str());
               throw_error_already_set();
            }
         }
      }
      else
      {
         PyErr_SetString(PyExc_IndexError, "index out of range");
         throw_error_already_set();
      }
   }

   object getitem(int index) const
   {
      if ((index < 0) || (index >= _nb_integral+_nb_float))
      {
         PyErr_SetString(PyExc_IndexError, "index out of range");
         throw_error_already_set();
      }
      return values[index];
   }

   int get_nb_integral() const { return _nb_integral; }
   int get_nb_float() const { return _nb_float; }
};

std::ostream& operator<<(std::ostream& os, const TreeValue& v)
{
   str s( v.values.attr("__str__")());
   std::string res;
   extract<std::string> x(s);
   if (x.check())
   {
      res= x();
      os << res;
   }

   return os;
}

}// namespace


// Module ======================================================================
BOOST_PYTHON_MODULE(ctree)
{
    /* class_< MyClass >("MyClass", init< const MyClass& >())
        .def(init< >())
        .def_readwrite("i", &MyClass::i)
        .def("f", &MyClass::f)
    ;*/

    class_< TreeValue >("TreeValue", init< const TreeValue& >())
        .def(init<boost::python::list>())
        .def("NbInt", &TreeValue::get_nb_integral)
        .def("NbFloat", &TreeValue::get_nb_float)
        .def("Type", &TreeValue::get_type)
        .def("Types", &TreeValue::get_types)
        .def("Values", &TreeValue::get_values)
        .def("NbFloat", &TreeValue::get_nb_float)
        .def("__getitem__", &TreeValue::getitem)
        .def("__setitem__", &TreeValue::setitem)
        .def(self_ns::str(self))
    ;

    class_< Default_tree >("CTree", init< const Default_tree& >())
        .def(init< optional< int, int, Default_tree::value > > ())
        .def(init< int, int, Default_tree::key, int >())
        .def(init< Unlabelled_typed_edge_tree&, optional< const Default_tree::value& > >())
        .def(init< const Typed_edge_tree<Int_fl_container>& >())
        .def(init< int, int, Unlabelled_typed_edge_tree& >())
        .def("__init__", make_constructor(Tree_wrapper_init1))
        // .def("Int", (const int& (Int_fl_container::*)(int) const)&Int_fl_container::Int, return_value_policy< copy_const_reference >())
        .def("Root", (Default_tree::vertex_descriptor (Default_tree::*)() const)&Default_tree::root,
                     "Root(self) -> vid. \n\n"
                     "Return the tree root vertex identifier (vid).")
        .def("IsRoot", &Default_tree::is_root,
                     "IsRoot(self, vid) -> Boolean. \n\n"
                     "Return true if and only if argument is the tree root "
                     "vertex identifier (vid).")
        .def("Size", &Default_tree::get_size,
                     "Size(self) -> integer. \n\n"
                     "Return the number of (allocated) vertices.")
        .def("Depth", &Generic_tree_wrapper_Depth0<Default_tree>,
                     "Depth(self) -> integer. \n\n"
                     "Return the tree depth.")
        .def("Depth", &Generic_tree_wrapper_Depth1<Default_tree>,
                     "Depth(self, vid) -> integer. \n\n"
                     "Return the depth of a given vertex, referred to "
                     "by its identifier (vid).")
        .def("NbInt", &Default_tree::get_nb_int,
                      "NbInt(self) -> integer. \n\n"
                      "Return the number of variables with integer type.")
        .def("NbFloat", &Default_tree::get_nb_float,
                        "NbFloat(self) -> integer. \n\n"
                        "Return the number of variables with floating type.")
        .def("NbChildren", &Default_tree::get_nb_children,
                           "NbChildren(self, vid) -> integer. \n\n"
                           "Return the number of children of a given vertex,\n"
                           "referred to by its identifier (vid).")
        .def("Order", &Generic_tree_wrapper_Order0<Default_tree>,
                     "Order(self) -> integer. \n\n"
                     "Return the tree branching order.")
        .def("Order", &Generic_tree_wrapper_Order1<Default_tree>,
                     "Order(self, vid) -> integer. \n\n"
                     "Return the branching order of a given vertex, "
                     "referred to by its identifier (vid).")
        .def("AddVertex", &Tree_wrapper_AddVertex,
                          "AddVertex(self, value) -> vid. \n\n"
                          "Add a vertex with given value to the tree and "
                          "return its identifier.")
        .def("AddEdge", &Default_tree::add_edge,
                        "AddEdge(self, vid, vid, type) -> True. \n\n"
                        "Add an edge of a given type between two vertices referred to "
                        "by their descriptors (vids).")
        .def("EdgeType", &Default_tree::edge_type,
                         "EdgeType(self, vid, vid) -> bool. \n\n"
                         "Return the type of an edge defined by two vertices, "
                         "referred to by their descriptors (vids).")
        .def("IsEdge", &Default_tree::is_edge,
                       "IsEdge(self, vid, vid) -> bool. \n\n"
                       "Return whether a given edge (defined by two vertices, "
                       "referred to by their vids) exists.")
        .def("Get", &Default_tree::get, return_value_policy< copy_const_reference >())
        .def("Put", &Default_tree::put)
        .def("Put", &Tree_wrapper_Put)
        .def("Parent", &Tree_wrapper_Parent<Default_tree>,
                       "Parent(self, vid) -> vid. \n\n"
                       "Return the parent of a vertex referred to "
                       "by its descriptor (vid).")
        .def("GetStructure", &Default_tree::get_structure, return_value_policy< manage_new_object>(),
                             "GetStructure(self) \n\n"
                             "Return the tree structure, i.e. the tree without the labels.")
        .def("SetStructure", &Default_tree::set_structure,
                             "SetStructure(self, structure, default_value) \n\n"
                             "Modify the tree structure.")
        .def("IidSimulation", &Tree_wrapper_IidSimulation,
                              "IidSimulation(self, list) \n\n"
                              "Simulation of the integer tree attributes using\n"
                              "a list of distributions.")
        .def("Display", &Generic_tree_wrapper_Display0<Default_tree>,
                        "Display(self) -> str \n\n"
                        "Display the whole tree.")
        .def("Display", &Generic_tree_wrapper_Display1<Default_tree>,
                        "Display(self, vid) -> str \n\n"
                        "Display the subtree rooted at given vertex, \n "
                        "referred to by its identifier (vid).")
        .def("Vertices", &Tree_wrapper_vertex_iterator<Default_tree>,
                         "Vertices(self) -> iterator \n\n"
                         "Get an iterator on the tree vertices.")
        .def("Children", &Tree_wrapper_children_iterator<Default_tree>,
                         "Children(self, vid) -> iterator \n\n"
                         "Get an iterator on the children of a given vertex,\n"
                         "referred to by its identifier (vid).")
        .def("Preorder", &Tree_wrapper_vertex_iterator_preorder<Default_tree>,
                         "Preorder(self) -> iterator \n\n"
                         "Get an iterator on the tree vertices, using a preorder traversing.")
        .def("Inorder", &Tree_wrapper_vertex_iterator_inorder<Default_tree>,
                        "Preorder(self) -> iterator \n\n"
                        "Get an iterator on the tree vertices, using an inorder traversing.")
        .def("Postorder", &Tree_wrapper_vertex_iterator_postorder<Default_tree>,
                         "Postorder(self) -> iterator \n\n"
                         "Get an iterator on the tree vertices, using a postorder traversing.")
        .def("Breadthorder", &Tree_wrapper_vertex_iterator_breadthorder<Default_tree>,
                             "Breadthorder(self) -> iterator \n\n"
                             "Get an iterator on the tree vertices, using a "
                             "breadth-first tree \n traversing.")
        .def("LeavesFirstorder", &Tree_wrapper_vertex_iterator_leavesfirstorder<Default_tree>,
                             "LeavesFirstorder(self) -> iterator \n\n"
                             "Get an iterator on the tree vertices where the leaf nodes "
                             "come first, then their parents, etc.")
        .def("__iter__", &Tree_wrapper_vertex_iterator<Default_tree>,
                         "__iter__(self) -> iterator \n\n"
                         "Get an iterator on the tree vertices.")
        .def("__str__", &Generic_tree_wrapper_Display0<Default_tree>,
                        "__str__(self) -> str \n\n"
                        "Display the whole tree.")
    ;

    class_< Tree_iterator<Default_tree> >("TreeIterator", no_init)
        .def("next", &Tree_iterator<Default_tree>::next)
        .def("__iter__", &Tree_iterator<Default_tree>::iter)
    ;

    class_< Tree_iterator_preorder<Default_tree> >("TreeIteratorPreorder", no_init)
        .def("next", &Tree_iterator_preorder<Default_tree>::next)
        .def("__iter__", &Tree_iterator_preorder<Default_tree>::iter)
    ;

    class_< Tree_iterator_inorder<Default_tree> >("TreeIteratorInorder", no_init)
        .def("next", &Tree_iterator_inorder<Default_tree>::next)
        .def("__iter__", &Tree_iterator_inorder<Default_tree>::iter)
    ;

    class_< Tree_iterator_postorder<Default_tree> >("TreeIteratorPostorder", no_init)
        .def("next", &Tree_iterator_postorder<Default_tree>::next)
        .def("__iter__", &Tree_iterator_postorder<Default_tree>::iter)
    ;

    class_< Tree_iterator_breadthorder<Default_tree> >("TreeIteratorBreadthorder", no_init)
        .def("next", &Tree_iterator_breadthorder<Default_tree>::next)
        .def("__iter__", &Tree_iterator_breadthorder<Default_tree>::iter)
    ;

    class_< Tree_iterator_leavesfirstorder<Default_tree> >("TreeIteratorLeavesFirstorder", no_init)
        .def("next", &Tree_iterator_leavesfirstorder<Default_tree>::next)
        .def("__iter__", &Tree_iterator_leavesfirstorder<Default_tree>::iter)
    ;

    class_< Children_iterator<Default_tree> >("ChildrenIterator", no_init)
        .def("next", &Children_iterator<Default_tree>::next)
        .def("__iter__", &Children_iterator<Default_tree>::iter)
    ;


    class_< Unlabelled_typed_edge_tree>("TreeStructure",
                                        "An implementation of tree structures,"
                                        "i.e. trees with no attributes.",
                                        init< const Unlabelled_typed_edge_tree& >())
        .def(init< >())
        .def(init< Unlabelled_typed_edge_tree::key, int, optional< Unlabelled_typed_edge_tree::value > >())
        .def(init< Distribution, optional<int, int> >())
        .def(init< Generic_tree<Int_fl_container>& >())
        .def(init< Typed_edge_tree<Int_fl_container>& >())
        .def(init< Default_tree& >())
        .def("Root", (Unlabelled_typed_edge_tree::vertex_descriptor (Unlabelled_typed_edge_tree::*)() const)&Unlabelled_typed_edge_tree::root,
                     "Root(self) -> vid. \n\n"
                     "Return the tree root vertex identifier (vid).")
        .def("IsRoot", &Unlabelled_typed_edge_tree::is_root,
                     "IsRoot(self, vid) -> Boolean. \n\n"
                     "Return true if and only if argument is the tree root "
                     "vertex identifier (vid).")
        .def("Size", &Unlabelled_typed_edge_tree::get_size,
                     "Size(self) -> integer. \n\n"
                     "Return the number of (allocated) vertices.")
        .def("Depth", Generic_tree_wrapper_Depth0<Unlabelled_typed_edge_tree>,
                     "Depth(self) -> integer. \n\n"
                     "Return the tree depth.")
        .def("Depth", Generic_tree_wrapper_Depth1<Unlabelled_typed_edge_tree>,
                     "Depth(self, vid) -> integer. \n\n"
                     "Return the depth of a given vertex, referred to "
                     "by its identifier (vid).")
        // .def("GetSize", &Unlabelled_typed_edge_tree::get_size)
        .def("AddVertex", &CharTree_wrapper_AddVertex,
                          "AddVertex(self, value) -> vid. \n\n"
                          "Add a vertex with given value to the tree and "
                          "return its identifier.")
        .def("AddEdge", &Unlabelled_typed_edge_tree::add_edge,
                        // &Typed_edge_tree<char>::add_edge,
                        "AddEdge(self, vid, vid, type) -> True. \n\n"
                        "Add an edge of a given type between two vertices referred to "
                        "by their descriptors (vids).")
        .def("EdgeType", &Typed_edge_tree<char>::edge_type,
                         "EdgeType(self, vid, vid) -> bool. \n\n"
                         "Return the type of an edge defined by two vertices, "
                         "referred to by their descriptors (vids).")
        .def("IsEdge", &Unlabelled_tree::is_edge,
                       "IsEdge(self, vid, vid) -> bool. \n\n"
                       "Return whether a given edge (defined by two vertices, "
                       "referred to by their vids) exists.")
    // .def("Get", &Unlabelled_typed_edge_tree::get, return_value_policy< copy_const_reference >())
        .def("Parent", &Tree_wrapper_Parent<Unlabelled_typed_edge_tree>,
                       "Parent(self, vid) -> vid. \n\n"
                       "Return the parent of a vertex referred to "
                       "by its descriptor (vid).")
        .def("Simulate", (void (Unlabelled_typed_edge_tree::*)(const Distribution&, int, int) )&Unlabelled_typed_edge_tree::simulation, Unlabelled_typed_edge_tree_simulation_overloads_1_3())
        .def("Display", &Generic_tree_wrapper_Display0<Unlabelled_typed_edge_tree>,
                        "Display(self) -> str \n\n"
                        "Display the whole tree.")
        .def("Display", &Generic_tree_wrapper_Display1<Unlabelled_typed_edge_tree>,
                        "Display(self, vid) -> str \n\n"
                        "Display the subtree rooted at given vertex, \n "
                        "referred to by its identifier (vid).")
        .def("Vertices", &Tree_wrapper_vertex_iterator<Unlabelled_typed_edge_tree>,
                         "Vertices(self) -> iterator \n\n"
                         "Get an iterator on the tree vertices.")
        .def("Children", &Tree_wrapper_children_iterator<Unlabelled_typed_edge_tree>,
                         "Children(self, vid) -> iterator \n\n"
                         "Get an iterator on the children of a given vertex,\n"
                         "referred to by its identifier (vid).")
        .def("__iter__", &Tree_wrapper_vertex_iterator<Unlabelled_typed_edge_tree>,
                         "__iter__(self) -> iterator \n\n"
                         "Get an iterator on the tree vertices.")
        .def("__str__", &Generic_tree_wrapper_Display0<Unlabelled_typed_edge_tree>,
                        "__str__(self) -> str \n\n"
                        "Display the whole tree.")
    ;

    class_< Tree_iterator<Unlabelled_typed_edge_tree> >("TreeStructureIterator", no_init)
        .def("next", &Tree_iterator<Unlabelled_typed_edge_tree>::next)
        .def("__iter__", &Tree_iterator<Unlabelled_typed_edge_tree>::iter)
    ;

    class_< Children_iterator<Unlabelled_typed_edge_tree> >("TreeStructureChildrenIterator", no_init)
        .def("next", &Children_iterator<Unlabelled_typed_edge_tree>::next)
        .def("__iter__", &Children_iterator<Unlabelled_typed_edge_tree>::iter)
    ;


    def("NB_TREES", NB_TREES_);
    def("I_DEFAULT_TREE_SIZE", I_DEFAULT_TREE_SIZE_);
    def("I_DEFAULT_TREE_DEPTH", I_DEFAULT_TREE_DEPTH_);

    enum_<UniqueInt<6> >("Characteristic")
        .value("FIRST_OCCURRENCE_ROOT", Stat_trees::FIRST_OCCURRENCE_ROOT)
        .value("FIRST_OCCURRENCE_LEAVES", Stat_trees::FIRST_OCCURRENCE_LEAVES)
        .value("SOJOURN_SIZE", Stat_trees::SOJOURN_SIZE)
        .value("NB_ZONES", Stat_trees::NB_ZONES)
        .value("NB_OCCURRENCES", Stat_trees::NB_OCCURRENCES)
        .value("OBSERVATION", Stat_trees::OBSERVATION)
        .export_values()
    ;

    def("SelectSubTree", select_typed_edge_subtree<Default_tree>,
                         return_value_policy< manage_new_object >());

    def("SelectSubTreeStructure", select_subtree<Unlabelled_typed_edge_tree>,
                                  return_value_policy< manage_new_object >());

}
