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
using namespace tree_statistic;

using boost::python::list;

// Declarations ================================================================

namespace  {
//
// BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Trees_extract_sequences_overloads_2_3, Trees_wrapper_extract_sequences, 2, 3)
//

template<int num, int id> struct UniqueInt
{
   int v;
   enum { value=num };

   UniqueInt(int _v) : v(_v) { }
   operator int() const { return v; }
};

/*************************************************************
 *
 *  Exporting constants as constant functions
 */

// int NB_TREES_() { return NB_TREES; }

/*************************************************************
 *
 *  Wrappers for Python class Trees:
 */

Trees* Trees_wrapper_init1(boost::python::list tree_list)
{
   int nb_trees, t, nb_variables, var, nbint, findex, iindex;
   int *itype= NULL, *rawtype= NULL;
   Trees::pt_tree_type_array otrees;
   object o;
   boost::python::list ltype;
   ostringstream error_message;
   bool status= true, several_errors= false;
   Trees *trees= NULL;
   // string s;

   nb_trees= boost::python::len(tree_list);
   if (nb_trees > 0)
   {
      otrees= new Default_tree*[nb_trees];
      for (t= 0; t < nb_trees; t++)
      {
         o= tree_list[t];
         try
         {
            extract<Default_tree> x(o.attr("_ctree")());
            if (x.check())
               otrees[t]= new Default_tree(x());
            else
               status=false;
         }
         catch (...)
         {
            status= false;
         }
         if (!status)
         {
            if (several_errors)
               error_message << endl;
            else
               several_errors= true;
            error_message << "incorrect type for element " << t
                          << " of argument list: expecting a Tree object"; // << endl;
            otrees[t]= NULL;
         }
      }
      if (!status)
      {
         for (t= 0; t < nb_trees; t++)
            if (otrees[t] != NULL)
            {
               delete otrees[t];
               otrees[t]= NULL;
            }
         delete [] otrees;
         otrees= NULL;

         throw_python_error(PyExc_TypeError, error_message);
         // PyErr_SetString(PyExc_TypeError, (error_message.str()).c_str());
         // throw_error_already_set();
      }
      else
      {
         ltype= extract<boost::python::list>(o.attr("Types")());
         // s= extract<string>(ltype.attr("__str__")());
         nb_variables= boost::python::len(ltype);
         itype= new int[nb_variables];
         rawtype= new int[nb_variables];

         // extraction of the variable types
         // and computation of the number of integer values
         nbint= 0;
         for (var=0; var < nb_variables; var++)
         {
            rawtype[var]= extract<int>(ltype[var]);
            if (rawtype[var] != REAL_VALUE)
               nbint++;
         }

         // The variables in Trees must be packed as follows:
         // [integral_variables, real_variables]
         iindex= 0;
         findex= nbint;
         for (var=0; var < nb_variables; var++)
         {
            if (rawtype[var] == REAL_VALUE)
               itype[findex++]= REAL_VALUE;
            else
               itype[iindex++]= rawtype[var];
         }

         trees= new Trees(nb_trees, itype, otrees, true);
         for (t= 0; t < nb_trees; t++)
            if (otrees[t] != NULL)
            {
               delete otrees[t];
               otrees[t]= NULL;
            }
         delete [] otrees;
         delete [] itype;
         delete [] rawtype;
         otrees= NULL;
         itype= NULL;
         rawtype= NULL;
         if (trees == NULL)
         {
            status= false;
            error_message << "could not initialize a Trees object from argument"; // << endl;
            throw_stat_tree_error(error_message);
         }
      }
   }
   else
   {
      status= false;
      error_message << "at least one tree required to initialize a Trees object"; // << endl;
      throw_python_error(PyExc_IndexError, error_message);
   }
   return trees;
}

void Trees_wrapper_to_int_type(Trees& reftree, int variable)
{
   // requires that "variable" is valid
   ostringstream error_message;
   StatError error;

   reftree.to_int_type(error, variable);
   if (error.get_nb_error() > 0)
   {
      error_message << error;
      throw_stat_tree_error(error_message);
   }
}

Trees* Trees_wrapper_transcode(const Trees& reftree, int variable,
                               boost::python::list symbols)
{
   // requires that "variable" is valid
   int c, nb_classes, nb_values;
   object o;
   ostringstream error_message;
   bool status= true;
   bool several_errors= false;
   StatError error;
   Trees *trees= NULL;
   Trees::int_array isymbols= NULL;

   nb_classes= boost::python::len(symbols);
   if (nb_classes > 0)
      isymbols= new int[nb_classes];
   else
   {
      status= false;
      error_message << "bad number of classes: " << nb_classes
                    << "; must be positive"; // << endl;
      throw_python_error(PyExc_ValueError, error_message);
   }
   nb_values= reftree.get_max_int_value(variable-1)
      - reftree.get_min_int_value(variable-1) + 1;

   if ((status) && (nb_classes < nb_values))
   {
      status= false;
      error_message << "bad number of classes: " << nb_classes
                    << "; must be at least the number of values ("
                    << nb_values << ")"; // << endl;
      throw_python_error(PyExc_ValueError, error_message);
   }
   if (status)
      for (c= 0; c < nb_classes; c++)
      {
         o= symbols[c];
         extract<int> x(o);
         if (x.check())
            isymbols[c]= x();
         else
            status= false;
         if (!status)
         {
            if (several_errors)
               error_message << endl;
            else
               several_errors= true;
            error_message << "incorrect type for element " << c
                          << " of argument list: expecting an integer"; // << endl;
         }
      }
   if (!status)
   {
      if (isymbols != NULL)
      {
         delete [] isymbols;
         isymbols= NULL;
      }
      throw_python_error(PyExc_TypeError, error_message);
   }
   else
   {
      trees= reftree.transcode(error, variable, isymbols);
      delete [] isymbols;
      isymbols= NULL;
      if (trees == NULL)
      {
         status= false;
         error_message << error; // << endl;
         throw_stat_tree_error(error_message);
      }
   }
   return trees;
}

Trees* Trees_wrapper_cluster_step(const Trees& reftree, int variable,
                                  int step)
{
   StatError error;
   Trees *trees= NULL;
   ostringstream error_message;

   trees= reftree.cluster(error, variable, step);
   if (trees == NULL)
   {
      error_message << error; // << endl;
      throw_stat_tree_error(error_message);
      // PyErr_SetString(PyExc_RuntimeError, (error_message.str()).c_str());
      // throw_error_already_set();
   }
   return trees;
}

Trees* Trees_wrapper_cluster_limit(const Trees& reftree, int variable,
                                   boost::python::list limits)
{
   int c, nb_class;
   object o;
   ostringstream error_message;
   bool status= true;
   bool several_errors= false;
   StatError error;
   Trees *trees= NULL;
   Trees::int_array ilimits= NULL;

   nb_class= boost::python::len(limits);
   if (nb_class > 0)
      ilimits= new int[nb_class];
   for (c= 0; c < nb_class; c++)
   {
      o= limits[c];
      extract<int> x(o);
      if (x.check())
         ilimits[c]= x();
      else
         status=false;
      if (!status)
      {
         if (several_errors)
            error_message << endl;
         else
            several_errors= true;
         error_message << "incorrect type for element " << c
                       << " of argument list: expecting an integer"; // << endl;
      }
   }
   if (!status)
   {
      delete [] ilimits;
      ilimits= NULL;
      throw_python_error(PyExc_TypeError, error_message);
   }
   else
   {
      trees= reftree.cluster(error, variable, nb_class, ilimits);
      delete [] ilimits;
      ilimits= NULL;
      if (trees == NULL)
      {
         status= false;
         error_message << error; // << endl;
         throw_stat_tree_error(error_message);
      }
   }
   return trees;
}

Trees* Trees_wrapper_difference(const Trees& reftree, int variable)
{
   ostringstream error_message;
   bool status= true;
   StatError error;
   Trees *trees= NULL;

   trees= reftree.difference(error, variable);
   if (trees == NULL)
   {
      status= false;
      error_message << error;
      throw_stat_tree_error(error_message);
   }
   return trees;
}

std::string Trees_wrapper_ascii_write(const Trees& reftree,
                                      bool exhaustive= false)
{
   std::stringstream ss;
   std::string res;

   reftree.ascii_write(ss, exhaustive);
   res= ss.str();

   return res;
}

void Trees_wrapper_build_sequences(const Trees& reftree, const char* prefix,
                                   bool all_paths= true, bool auto_axis= false)
{
   ostringstream error_message;
   bool status= true;
   StatError error;
   Sequences* seq= NULL;


   seq= reftree.build_sequences(error, all_paths, auto_axis);
   if (seq == NULL)
      status= false;

   if (status)
   {
      status= seq->ascii_data_write(error, prefix);
      delete seq;
      seq= NULL;
   }
   if (!status)
   {
      error_message << error;
      throw_stat_tree_error(error_message);
   }
}

Sequences* Trees_wrapper_build_py_sequences(const Trees& reftree,
                                            bool all_paths = true,
                                            bool auto_axis = true)
{
   ostringstream error_message;
   bool status= true;
   StatError error;
   Sequences* seq= NULL;


   seq= reftree.build_sequences(error, all_paths, auto_axis);
   if (seq == NULL)
      status= false;
   if (!status)
   {
      error_message << error;
      throw_stat_tree_error(error_message);
   }
   else
       return seq;
}

void Trees_wrapper_build_vectors_path(const Trees& reftree, const char* prefix)
{
   ostringstream error_message;
   bool status= true;
   StatError error;
   Vectors* vec= NULL;


   vec= reftree.build_vectors(error);
   if (vec == NULL)
      status= false;

   if (status)
   {
      status= vec->ascii_data_write(error, prefix);
      delete vec;
      vec= NULL;
   }
   if (!status)
   {
      error_message << error;
      throw_stat_tree_error(error_message);
   }
}

Vectors* Trees_wrapper_build_vectors_object(const Trees& reftree)
{
   ostringstream error_message;
   bool status= true;
   StatError error;
   Vectors* vec= NULL;


   vec = reftree.build_vectors(error);
   if (vec == NULL)
      status= false;

   if (!status)
   {
      error_message << error;
      throw_stat_tree_error(error_message);
   }
   return vec;
}

DiscreteDistributionData* Trees_wrapper_extract_value(const Trees& reftree, int variable)
{
   ostringstream error_message;
   // bool status= true;
   StatError error;
   DiscreteDistributionData *histo= NULL;


   histo= reftree.extract(error, variable);
   if (histo == NULL)
   {
      error_message << error;
      throw_stat_tree_error(error_message);
   }
   return histo;
}

DiscreteDistributionData* Trees_wrapper_extract_feature(const Trees& reftree,
                                                        int type,
                                                        int variable,
                                                        int value)
{
   ostringstream error_message;
   StatError error;
   DiscreteDistributionData *histo= NULL;

   histo= reftree.extract(error, type, variable, value);
   if (histo == NULL)
   {
      error_message << error;
      throw_stat_tree_error(error_message);
   }
   return histo;
}

boost::python::list Trees_wrapper_max(const Trees& reftree)
{
   boost::python::list res;
   int i;
   Trees::value val;

   val= reftree.get_max_value();
   for(i= 0; i < reftree.get_nb_int(); i++)
      res.append(val.Int(i));
   for(i= 0; i < reftree.get_nb_float(); i++)
      res.append(val.Double(i));

   return res;
}

boost::python::list Trees_wrapper_min(const Trees& reftree)
{
   boost::python::list res;
   int i;
   Trees::value val;

   val= reftree.get_min_value();
   for(i= 0; i < reftree.get_nb_int(); i++)
      res.append(val.Int(i));
   for(i= 0; i < reftree.get_nb_float(); i++)
      res.append(val.Double(i));

   return res;
}

Trees* Trees_wrapper_select_variable(const Trees& reftree,
                                     boost::python::list variables, bool keep)
{
   int nb_variables, var;
   int *variable_list;
   bool several_errors= false;
   object o;
   StatError error;
   ostringstream error_message;
   bool status= true;
   Trees *trees= NULL;

   nb_variables= boost::python::len(variables);
   if (nb_variables > 0)
   {
      variable_list= new int[nb_variables];
      for (var= 0; var < nb_variables; var++)
      {
         o= variables[var];
         extract<int> x(o);
         if (x.check())
            variable_list[var]= x();
         else
            status=false;
         if (!status)
         {
            if (several_errors)
               error_message << endl;
            else
               several_errors= true;
            error_message << "incorrect type for element " << var
                          << " of argument list: expecting an integer";
            variable_list[var]= I_DEFAULT;
         }
      }
      if (!status)
      {
         delete [] variable_list;
         variable_list= NULL;
         throw_python_error(PyExc_TypeError, error_message);
      }
      else
      {
         trees= reftree.select_variable(error, nb_variables,
                                        variable_list, keep);
         delete [] variable_list;
         variable_list= NULL;
         if (trees == NULL)
         {
            status= false;
            error_message << error;
            throw_stat_tree_error(error_message);
         }
      }
   }
   else
   {
      status= false;
      trees= new Trees(reftree);
      error_message << "at least one variable required to select variables";
      throw_python_error(PyExc_UserWarning, error_message);
   }
   return trees;
}

Trees* Trees_wrapper_select_individual(const Trees& reftree,
                                       boost::python::list identifiers, bool keep)
{
   register int nb_ids;
   int t= 0;
   int *id_list;
   bool several_errors= false;
   object o;
   StatError error;
   ostringstream error_message;
   bool status= true;
   Trees *trees= NULL;

   nb_ids= boost::python::len(identifiers);
   if (nb_ids > 0)
   {
      id_list= new int[nb_ids];
      for (t= 0; t < nb_ids; t++)
      {
         o= identifiers[t];
         extract<int> x(o);
         if (x.check())
            id_list[t]= x();
         else
            status=false;
         if (!status)
         {
            if (several_errors)
               error_message << endl;
            else
               several_errors= true;
            error_message << "incorrect type for element " << t
                          << " of argument list: expecting an integer";
            id_list[t]= I_DEFAULT;
         }
      }
      if (!status)
      {
         delete [] id_list;
         id_list= NULL;
         throw_python_error(PyExc_TypeError, error_message);
      }
      else
      {
         trees= reftree.select_individual(error, nb_ids, id_list, keep);
         delete [] id_list;
         id_list= NULL;
         if (trees == NULL)
         {
            status= false;
            error_message << error;
            throw_stat_tree_error(error_message);
         }
      }
   }
   else
   {
      status= false;
      trees= new Trees(reftree);
      error_message << "at least one identifier required to select individuals";
      throw_python_error(PyExc_UserWarning, error_message);
   }
   return trees;
}

Trees* Trees_wrapper_segmentation_extract_value(const Trees& reftree,
                                                int variable,
                                                int value,
                                                bool keep)
{
   int *ivalue= new int[1];
   StatError error;
   ostringstream error_message;
   bool status= true;
   Trees *trees= NULL;

   ivalue[0]= value;

   trees= reftree.segmentation_extract(error, variable,
                                       1, ivalue, keep);
   delete [] ivalue;
   ivalue= NULL;
   if (trees == NULL)
   {
      status= false;
      error_message << error;
      throw_stat_tree_error(error_message);
   }
   return trees;
}

Trees* Trees_wrapper_segmentation_extract_values(const Trees& reftree,
                                                 int variable,
                                                 boost::python::list values,
                                                 bool keep)
{
   int nb_values, val;
   int *ivalue= NULL;
   StatError error;
   ostringstream error_message;
   bool status= true;
   object o;
   bool several_errors= false;
   Trees *trees= NULL;

   nb_values= boost::python::len(values);

   if (nb_values > 0)
   {
      ivalue= new int[nb_values];
      for(val= 0; val < nb_values; val++)
      {
         o= values[val];
         extract<int> x(o);
         if (x.check())
            ivalue[val]= x();
         else
            status= false;
         if (!status)
         {
            if (several_errors)
               error_message << endl;
            else
               several_errors= true;
            error_message << "incorrect type for element " << val
                          << " of argument list: expecting an integer";
            ivalue[val]= I_DEFAULT;
         }
      }
      if (!status)
      {
         delete [] ivalue;
         ivalue= NULL;
         throw_python_error(PyExc_TypeError, error_message);
      }
      else
      {
         trees= reftree.segmentation_extract(error, variable,
                                             nb_values, ivalue, keep);
         delete [] ivalue;
         ivalue= NULL;
         if (trees == NULL)
         {
            status= false;
            error_message << error;
            throw_stat_tree_error(error_message);
         }
      }
   }
   else
   {
      status= false;
      error_message << "at least one value required to extract segmentation";
      throw_stat_tree_error(error_message);
   }
   return trees;
}

Trees* Trees_wrapper_shift_int(const Trees& reftree, int variable, int shift)
{
   ostringstream error_message;
   StatError error;
   Trees *trees= NULL;

   trees= reftree.shift(error, variable, shift);
   if (trees == NULL)
   {
      error_message << error;
      throw_stat_tree_error(error_message);
   }
   return trees;
}

Trees* Trees_wrapper_shift_float(const Trees& reftree, int variable, double shift)
{
   ostringstream error_message;
   StatError error;
   Trees *trees= NULL;

   trees= reftree.shift(error, variable, shift);
   if (trees == NULL)
   {
      error_message << error;
      throw_stat_tree_error(error_message);
   }
   return trees;
}

Trees* Trees_wrapper_merge(const Trees& reftree, StatError& error,
                           boost::python::list tree_list)
{
   int nb_trees, t; //, nb_variables, var;
   Trees::pt_observed_trees_array otrees;
   object o;
   ostringstream error_message;
   bool several_errors= false;
   bool status= true;
   Trees *trees= NULL;

   nb_trees= boost::python::len(tree_list);
   if (nb_trees > 0)
   {
      otrees= new Trees*[nb_trees];
      for (t= 0; t < nb_trees; t++)
      {
         o= tree_list[t];
         extract<Trees> x(o);
         if (x.check())
            otrees[t]= new Trees(x());
         else
            status=false;
         if (!status)
         {
            if (several_errors)
               error_message << endl;
            else
               several_errors= true;
            error_message << "incorrect type for element " << t
                          << " of argument list: expecting a Trees object"; // << endl;
            otrees[t]= NULL;
         }
      }
      if (!status)
      {
         for (t= 0; t < nb_trees; t++)
            if (otrees[t] != NULL)
            {
               delete otrees[t];
               otrees[t]= NULL;
            }
         delete [] otrees;
         otrees= NULL;
         throw_python_error(PyExc_TypeError, error_message);
      }
      else
      {
         trees= reftree.merge(error, nb_trees, otrees);
         for (t= 0; t < nb_trees; t++)
            if (otrees[t] != NULL)
            {
               delete otrees[t];
               otrees[t]= NULL;
            }
         delete [] otrees;
         otrees= NULL;
         if (trees == NULL)
         {
            status= false;
            error_message << "could not merge Trees object from arguments"; // << endl;
            throw_stat_tree_error(error_message);
         }
      }
   }
   else
   {
      status= false;
      trees= new Trees(reftree);
      error_message << "at least one tree required to merge Trees objects"; // << endl;
      throw_python_error(PyExc_UserWarning, error_message);
      return trees;
   }
   return trees;
}

unsigned int Trees_wrapper_max_order(const Trees& reftree)
{
   typedef Trees::tree_type tree_type;
   register int t;
   unsigned int res= 1;
   tree_type *tree;

   for(t= 0; t < reftree.get_nb_trees(); t++)
   {
      tree= reftree.get_tree(t);
      res= max(res, tree->get_branching_order());
      delete tree;
      tree= NULL;
   }

   return res;
}


Trees* Trees_wrapper_merge2(const Trees& reftree, boost::python::list tree_list)
{
   int nb_trees, t; //, nb_variables, var;
   Trees::pt_observed_trees_array otrees;
   bool several_errors= false;
   object o;
   StatError error;
   ostringstream error_message;
   bool status= true, current_status;
   Trees *trees= NULL;

#ifdef   DEBUG
    cerr << "entering Trees_wrapper_merge2" << endl;
#endif
   nb_trees= boost::python::len(tree_list);
   if (nb_trees > 0)
   {
      otrees= new Trees*[nb_trees];
      for (t= 0; t < nb_trees; t++)
      {
         current_status= true;
         o= tree_list[t];
         extract<Trees> x(o);
         if (x.check())
            otrees[t]= new Trees(x());
         else
         {
            status= false;
            current_status= false;
         }
         if (!current_status)
         {
            if (several_errors)
               error_message << endl;
            else
               several_errors= true;
            error_message << "incorrect type for element " << t
                          << " of argument list: expecting a Trees object"; // << endl;
            otrees[t]= NULL;
         }
      }
      if (!status)
      {
         for (t= 0; t < nb_trees; t++)
            if (otrees[t] != NULL)
            {
               delete otrees[t];
               otrees[t]= NULL;
            }
         delete [] otrees;
         otrees= NULL;
         throw_python_error(PyExc_TypeError, error_message);
      }
      else
      {
#ifdef   DEBUG
         cerr << "Calling Trees::merge(" << error << ", " << nb_trees << ", "
              << otrees << ")" << endl;
#endif
         trees= reftree.merge(error, nb_trees, otrees);
#ifdef   DEBUG
         cerr << "Finished Trees::merge" << endl;
#endif
         for (t= 0; t < nb_trees; t++)
            if (otrees[t] != NULL)
            {
               delete otrees[t];
               otrees[t]= NULL;
            }
         delete [] otrees;
         otrees= NULL;
         if (trees == NULL)
         {
            status= false;
            error_message << error; // << endl;
            throw_stat_tree_error(error_message);
#ifdef      DEBUG
            cerr << "Sending exception." << endl;
#endif
            #ifdef      DEBUG
            cerr << "Exception sent." << endl;
#endif
         }
      }
   }
   else
   {
      status= false;
      trees= new Trees(reftree);
      error_message << "at least one tree required to merge Trees objects"; // << endl;
      throw_python_error(PyExc_UserWarning, error_message);
      return trees;
   }
   return trees;
}

Trees* Trees_wrapper_merge_variable(const Trees& reftree, boost::python::list tree_list)
{
   int nb_trees, t;
   Trees::pt_observed_trees_array otrees;
   bool several_errors= false;
   object o;
   StatError error;
   ostringstream error_message;
   bool status= true, current_status;
   Trees *trees= NULL;

   nb_trees= boost::python::len(tree_list);
   if (nb_trees > 0)
   {
      otrees= new Trees*[nb_trees];
      for (t= 0; t < nb_trees; t++)
      {
         current_status= true;
         o= tree_list[t];
         extract<Trees> x(o);
         if (x.check())
            otrees[t]= new Trees(x());
         else
         {
            status= false;
            current_status= false;
         }
         if (!current_status)
         {
            if (several_errors)
               error_message << endl;
            else
               several_errors= true;
            error_message << "incorrect type for element " << t
                          << " of argument list: expecting a Trees object"; // << endl;
            otrees[t]= NULL;
         }
      }
      if (!status)
      {
         for (t= 0; t < nb_trees; t++)
            if (otrees[t] != NULL)
            {
               delete otrees[t];
               otrees[t]= NULL;
            }
         delete [] otrees;
         otrees= NULL;
         throw_python_error(PyExc_TypeError, error_message);
      }
      else
      {
         trees= reftree.merge_variable(error, nb_trees, otrees);
         for (t= 0; t < nb_trees; t++)
            if (otrees[t] != NULL)
            {
               delete otrees[t];
               otrees[t]= NULL;
            }
         delete [] otrees;
         otrees= NULL;
         if (trees == NULL)
         {
            status= false;
            error_message << error; // << endl;
            throw_stat_tree_error(error_message);
         }
      }
   }
   else
   {
      status= false;
      trees= new Trees(reftree);
      error_message << "at least one tree required to merge variables"; // << endl;
      throw_python_error(PyExc_UserWarning, error_message);
      return trees;
   }
   return trees;
}

MultiPlotSet* Trees_wrapper_get_plotable(const Trees& reftree,
                                         int plot_type,
                     int variable)
{
   StatError error;
   ostringstream error_message;
   MultiPlotSet *plotset= NULL;

   plotset= reftree.get_plotable(error, plot_type, variable);

   if (plotset == NULL)
   {
      error_message << error;
      throw_stat_tree_error(error_message);
   }

   return plotset;

}
void Trees_wrapper_plot_write(const Trees& reftree,
                              const char* prefix,
                              const char* title)
{
   bool status= true;
   ostringstream error_message;
   StatError error;

   status= reftree.plot_write(error, prefix, title);
   if (not status)
   {
      error_message << error;
      throw_stat_tree_error(error_message);
   }
}

std::string Trees_wrapper_line_write(const Trees& reftree)
{
   std::stringstream ss;
   std::string res;

   reftree.line_write(ss);
   res= ss.str();

   return res;
}

}// namespace


// Module ======================================================================
BOOST_PYTHON_MODULE(ctrees)
{

    // Error initialisation
    object stat_tree_errors = import("openalea.tree_statistic._errors");
    // Import StatTreeError
    StatTreeError = stat_tree_errors.attr("StatTreeError");

    class_< Trees >("CTrees", init< const Trees& >())
        .def(init< optional< int, int, int > > ())
        .def(init< int, int, const FrequencyDistribution&, const FrequencyDistribution&,
                   optional< bool, bool > >())
        .def(init< const Trees&, int, Trees::int_array>())
        .def("__init__", make_constructor(Trees_wrapper_init1))
        .def("Cluster", &Trees_wrapper_cluster_step,
                        return_value_policy< manage_new_object >(),
                        "Cluster(self, variable, step) ->CTrees \n\n"
                        "Cluster the values of a variable, using a given step.")
        .def("Cluster", &Trees_wrapper_cluster_limit,
                        return_value_policy< manage_new_object >(),
                        "Cluster(self, variable, limits) ->CTrees \n\n"
                        "Cluster the values of a variable, using given limits.")
        .def("Difference", &Trees_wrapper_difference,
                           return_value_policy< manage_new_object >(),
                           "Difference(self, variable) -> CTrees \n\n"
                           "First-order differentiation of Trees.")
        .def("Display", &Trees_wrapper_ascii_write)
        .def("BuildSequences", &Trees_wrapper_build_sequences,
                               "BuildSequences(self, bool, bool) -> void \n\n"
                               "Build sequences from trees and print them"
                               "into a file.")
        .def("BuildPySequences", &Trees_wrapper_build_py_sequences,
                                 return_value_policy< manage_new_object >(),
                                 "BuildPySequences(self, bool, bool) -> void \n\n"
                                 "Build sequences from trees into a "
                                 "sequence_analysis.Sequence object.")
        .def("BuildVectors", &Trees_wrapper_build_vectors_path,
                             "BuildVectors(self, path) -> void \n\n"
                             "Build vectors from trees and print them"
                             "into a file.")
        .def("BuildVectors", &Trees_wrapper_build_vectors_object,
                          return_value_policy< manage_new_object >(),
                             "BuildVectors(self, path) -> void \n\n"
                             "Build vectors object from trees.")
        .def("ExtractSizeHistogram", &Trees::extract_size,
                                                 return_value_policy< manage_new_object >(),
                                                "ExtractSizeHistogram(self) -> FrequencyDistribution \n\n"
                                                "Extract the frequency distribution of tree size.")
        .def("ExtractNbChildrenHistogram", &Trees::extract_nb_children,
                                            return_value_policy< manage_new_object >())
        .def("ExtractValueHistogram", &Trees_wrapper_extract_value,
                                      return_value_policy< manage_new_object >())
        .def("ExtractFeatureHistogram", &Trees_wrapper_extract_feature,
                                        return_value_policy< manage_new_object >())
        .def("IsCharacteristic", &Trees::is_characteristic,
                             "IsCharacteristic(self, variable, charac) -> bool \n\n",
                             "Check whether a characteristic is present"
                             "for a given variable")
        .def("NbInt", &Trees::get_nb_int,
                      "NbInt(self) -> int \n\n"
                      "return the number of variables "
                      "with integer type\n")
        .def("NbFloat", &Trees::get_nb_float,
                      "NbInt(self) -> int \n\n"
                      "return the number of variables "
                      "with floating type\n")
        .def("NbValues", &Trees::get_nb_values,
                         "NbValues(self, variable) -> int \n\n"
                         "return the number of values of a variable "
                         "with finite values\n")
        .def("NbTrees", &Trees::get_nb_trees)
        .def("Max", &Trees_wrapper_max,
                    "Max(self) -> list. \n\n"
                    "Return the maximal value of each variable.")
        .def("Min", &Trees_wrapper_min,
                    "Max(self) -> list. \n\n"
                    "Return the minimal value of each variable.")
        .def("MaxDepth", &Trees::get_max_depth,
                         "MaxDepth(self) -> int. \n\n"
                         "Return the maximal depth of the trees.")
        .def("MaxOrder", &Trees_wrapper_max_order,
                         "MaxOrder(self) -> int. \n\n"
                         "Return the maximal branching order of the trees.")
        .def("Merge", &Trees_wrapper_merge2,
                      return_value_policy< manage_new_object >(),
                      "Merge(self, tree_list) -> CTrees. \n\n"
                      "Merge the trees of self and that in the list "
                      "given as an argument.")
        .def("MergeVariable", &Trees_wrapper_merge_variable,
                              return_value_policy< manage_new_object >(),
                              "MergeVariable(self, tree_list) -> CTrees. \n\n"
                              "Merge the variables of self and that of the trees "
                              "in the list given as an argument.")
        .def("get_plotable", &Trees_wrapper_get_plotable,
                              return_value_policy< manage_new_object >(),
                         "Fill MultiPlotSet structure.")
        .def("plot_write", &Trees_wrapper_plot_write,
                       "Write into a gnuplot file.")
        .def("SegmentationExtract", &Trees_wrapper_segmentation_extract_value,
                                    return_value_policy< manage_new_object >(),
                                    "SegmentationExtract(self, variable, value, mode) -> CTrees. \n\n"
                                    "Select homogeneous subtrees of self according "
                                    "to a certain variable and value of that variable.")
        .def("SegmentationExtract", &Trees_wrapper_segmentation_extract_values,
                                    return_value_policy< manage_new_object >(),
                                    "SegmentationExtract(self, variable, values, mode) -> CTrees. \n\n"
                                    "Select homogeneous subtrees of self according "
                                    "to a certain variable and range of values.")
        .def("SelectVariable", &Trees_wrapper_select_variable,
                               return_value_policy< manage_new_object >(),
                               "SelectVariable(self, variable, mode) -> CTrees. \n\n"
                               "Select the given variables of self.")
        .def("SelectIndividual", &Trees_wrapper_select_individual,
                                 return_value_policy< manage_new_object >(),
                                 "SelectIndividual(self, identifiers, mode) -> CTrees. \n\n"
                                 "Select the given trees of self.")
        .def("Shift", &Trees_wrapper_shift_float,
                      return_value_policy< manage_new_object >(),
                      "Shift(self, variable, shift) -> CTrees \n\n"
                      "Shift (i.e. translate) the values of a Trees object.")
        .def("Shift", &Trees_wrapper_shift_int,
                      return_value_policy< manage_new_object >(),
                      "Shift(self, variable, shift) -> CTrees \n\n"
                      "Shift (i.e. translate) the values of a Trees object.")
        // note that the order between both shift wrappers is relevant,
        // since an integer is a kind of float
        .def("Size", &Trees::get_total_size,
                     "Size(self) -> int \n\n"
                     "Return total number of vertices.")
        .def("Size", &Trees::get_size,
                     "Size(self, int) -> int \n\n"
                     "Return the number of vertices of a given tree.")
        .def("ToIntType", &Trees_wrapper_to_int_type,
                     "ToIntType(self, variable) -> void \n\n"
                     "Switch a variable type to INT_VALUE.")
        .def("Transcode", &Trees_wrapper_transcode,
                        return_value_policy< manage_new_object >(),
                        "Transcode(self, variable, new_values) ->CTrees \n\n"
                        "Transcode the values of a variable "
                        "using a given list of the new values.")
        .def("Tree", &Trees::get_tree,
                     return_value_policy< manage_new_object >(),
                     "Tree(self, id) -> CTree \n\n"
                     "Return a given tree from the Trees object.")
        .def("Type", &Trees::get_type)
        .def("__str__", &Trees_wrapper_line_write)
        // .def(self_ns::str(self))
    ;

}
