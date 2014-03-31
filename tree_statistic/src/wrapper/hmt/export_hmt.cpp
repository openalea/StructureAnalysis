/*------------------------------------------------------------------------------
 *
 *        VPlants.Tree-Statistic : VPlants Tree-Statistic module
 *        HiddenMarkovTrees
 *
 *        Copyright 2006-2007 INRIA - CIRAD - INRA
 *
 *        File author(s): Jean-Baptiste Durand <Jean-Baptiste.Durand@imag.fr>
 *
 *        Distributed under the GPL 2.0 License.
 *        See accompanying file LICENSE.txt or copy at
 *           http://www.gnu.org/licenses/gpl-2.0.txt
 *
 *        OpenAlea WebSite : http://openalea.gforge.inria.fr
 *
 *        $Id: export_hmt.cpp 9099 2010-06-08 09:03:00Z pradal $
 *
 *-----------------------------------------------------------------------------*/
// Includes ====================================================================
#include "tree/basic_visitors.h"
#include "tree/tree_traits.h"
#include "tree/tree_simple.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/distribution.h"
#include "stat_tool/vectors.h"
#include "stat_tool/discrete_mixture.h"

#include "tree_statistic/tree_labels.h"
#include "tree_statistic/int_fl_containers.h"
#include "tree_statistic/generic_typed_edge_tree.h"
#include "tree_statistic/typed_edge_trees.h"
#include "tree_statistic/hidden_markov_tree.h"
#include "tree_statistic/hidden_markov_ind_out_tree.h"

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

// Declarations ================================================================
template<int num> struct UniqueInt { int v; enum { value=num };
UniqueInt(int _v) : v(_v) { } operator int() const { return v; } };

/*************************************************************
 *
 *  Wrapper for Error class StatTreeError:
 */

#define WRAP StatTreeErrorWrap
void class_errors()
{
    // Error initialisation
    object stat_tree_errors = import("openalea.tree_statistic._errors");
    // Import StatError
    StatTreeError = stat_tree_errors.attr("StatTreeError");
};
#undef WRAP

/*************************************************************
 *
 *  Wrapper for enum. EntropyAlgorithm:
 */

#define WRAP EntropyAlgorithmWrap
void enum_entropy_algorithms()
{
    enum_<UniqueInt<2> >("EntropyAlgorithm")
        .value("UPWARD", Stat_trees::UPWARD)
        .value("DOWNWARD", Stat_trees::DOWNWARD)
        .export_values()
    ;
};
#undef WRAP

/*************************************************************
 *
 *  Wrappers for Python class CHmt:
 */

#define WRAP HiddenMarkovTreeWrap
class WRAP
{

public :

   static int CHmt_wrapper_nb_values(const HiddenMarkovTree& hmt, int variable)
   {
      int res= 0;
      ostringstream error_message;

      if ((variable >= 0) && (variable < hmt.get_nb_ioutput_process()))
         res= hmt.get_nb_values(variable);
      else
      {
         error_message << "Bad variable: " << variable << endl;
         throw_stat_tree_error(error_message);
      }
      return res;
   }

   static HiddenMarkovTreeData*
   CHmt_wrapper_extract_data(const HiddenMarkovTree& hmt)
   {
      HiddenMarkovTreeData *res;
      StatError error;
      ostringstream error_message;

      res= hmt.extract_data(error);
      if (res == NULL)
      {
         error_message << error;
         throw_stat_tree_error(error_message);
      }
      return res;
   }

   static DiscreteParametricModel*
   CHmt_wrapper_extract(const HiddenMarkovTree& hmt, int itype,
                        int ivariable, int ivalue)
   {
      DiscreteParametricModel *res;
      StatError error;
      ostringstream error_message;

      res= hmt.extract(error, itype, ivariable, ivalue);
      if (res == NULL)
      {
         error_message << error;
         throw_stat_tree_error(error_message);
      }
      return res;
   }

   static double CHmt_wrapper_likelihood(const HiddenMarkovTree& hmt,
                                         const Trees& trees)
   {
      double res;
      ostringstream error_message;

      res= hmt.likelihood_computation(trees, I_DEFAULT);
      if (res <= D_INF)
      {
         error_message << "Cannot compute the likelihood";
         throw_stat_tree_error(error_message);
      }
      return res;
   }

   static str CHmt_wrapper_ascii_write0(const HiddenMarkovTree& hmt)
   {
      std::stringstream s;
      str res;

      hmt.line_write(s);
      res= str(s.str());
      return res;
   }

}; // end class WRAP

void class_hmt()
{
    class_< HiddenMarkovTree >
    ("CHmt", init< const HiddenMarkovTree&, optional< bool, bool> >())
        .def("ExtractData", WRAP::CHmt_wrapper_extract_data,
                            return_value_policy< manage_new_object >())
        .def("Extract", WRAP::CHmt_wrapper_extract,
                        return_value_policy< manage_new_object >(),
                        "Extract(self, type, variable, value) -> _DiscreteParametricModel \n\n"
                        "Extract some characteristic or observation distribution")
        .def("IsParametric", &HiddenMarkovTree::is_parametric,
                            "IsParametric(self, variable) -> bool \n\n"
                            "Return True if process 'variable' "
                            "is parametric")
        .def("Likelihood", WRAP::CHmt_wrapper_likelihood,
                            "Likelihood(self, trees) -> float \n\n"
                            "Return the likelihood of the parameters"
                            "for the given trees")
        .def("NbInt", &HiddenMarkovTree::get_nb_ioutput_process)
        .def("NbFloat", &HiddenMarkovTree::get_nb_doutput_process)
        .def("NbStates", &HiddenMarkovTree::get_nb_state,
                         "NbStates(self) -> int \n\n"
                         "return the number of states "
                         "of the hidden Markov tree\n")
        .def("NbValues", WRAP::CHmt_wrapper_nb_values,
                         "NbValues(self) -> int \n\n"
                         "return the number of values for a given process "
                         "of the hidden Markov tree\n")
        .def("__str__", WRAP::CHmt_wrapper_ascii_write0,
                        "__str__(self) -> str \n\n"
                        "Display the hidden Markov tree.")
    ;
};
#undef WRAP

/*************************************************************
 *
 *  Wrappers for Python class CiHmot:
 */

#define WRAP HiddenMarkovIndOutTreeWrap
class WRAP
{

public :

   static HiddenMarkovTreeData*
   CiHmot_wrapper_extract_data(const HiddenMarkovIndOutTree& hmt)
   {
      HiddenMarkovTreeData *res;
      StatError error;
      ostringstream error_message;

      res= hmt.extract_data(error);
      if (res == NULL)
      {
         error_message << error;
         throw_stat_tree_error(error_message);
      }
      return res;
   }

   static str CiHmot_wrapper_display(const HiddenMarkovIndOutTree& hmt,
                                     bool exhaustive)
   {
      std::stringstream s;
      str res;

      hmt.ascii_write(s, exhaustive);
      res= str(s.str());
      return res;
   }

   static bool CiHmot_wrapper_file_ascii_write2(const HiddenMarkovIndOutTree& hmt,
                                                const char * path, bool exhaustive)
   {
      StatError error;
      bool res;
      ostringstream error_message;

      res= hmt.ascii_write(error, path, exhaustive);
      if (!res)
      {
         error_message << error;
         throw_stat_tree_error(error_message);
      }
      return res;
   }

   static bool CiHmot_wrapper_file_ascii_write1(const HiddenMarkovIndOutTree& hmt,
                                                const char * path)
   { return CiHmot_wrapper_file_ascii_write2(hmt, path, false); }

   static void CiHmot_wrapper_state_permutation(const HiddenMarkovIndOutTree& hmt,
                                                boost::python::list perm)
   {
      bool status= true, several_errors= false;
      int llength, i;
      ostringstream error_message;
      object o;
      StatError error;
      int *iperm= NULL;

      llength= boost::python::len(perm);
      if ((llength > 0) && (llength == hmt.get_nb_state()))
      {
         iperm= new int[llength];
         for (i= 0; i < llength; i++)
         {
            o= perm[i];
            try
            {
               extract<int> x(o);
               if (x.check())
                  iperm[i]= x();
               else
                  status= false;
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
               error_message << "incorrect type for element " << i
                             << " of argument list: expecting an int ";
            }
         }
         if (!status)
         {
            delete [] iperm;
            iperm= NULL;
            throw_python_error(PyExc_TypeError, error_message);
         }
      }
      else
      {
         status= false;
         error_message << "incorrect permutation" << endl;
         throw_stat_tree_error(error_message);
      }
      if (status)
      {
         hmt.state_permutation(error, iperm);
         if (error.get_nb_error() > 0)
         {
            error_message << error;
            throw_stat_tree_error(error_message);
         }
         delete [] iperm;
         iperm= NULL;
      }
   }

   static bool CiHmot_wrapper_spreadsheet_write1(const HiddenMarkovIndOutTree& hmt,
                                                 const char * path)
   {
      StatError error;
      bool res;
      ostringstream error_message;

      res= hmt.spreadsheet_write(error, path);
      if (!res)
      {
         error_message << error;
         throw_stat_tree_error(error_message);
      }
      return res;
   }

   static HiddenMarkovTreeData*
   CiHmot_wrapper_simulate_histo(const HiddenMarkovIndOutTree& hmt,
                                 const FrequencyDistribution& ihsize,
                                 const FrequencyDistribution& ihnb_children,
                                 bool counting_flag= true,
                                 bool divergence_flag= false)
   {
      StatError error;
      HiddenMarkovTreeData* markov_data= NULL;
      ostringstream error_message;

      markov_data= hmt.simulation(error, ihsize, ihnb_children,
                                  counting_flag, divergence_flag);
      if (markov_data == NULL)
      {
         error_message << error;
         throw_stat_tree_error(error_message);
      }
      return markov_data;
   }

   static HiddenMarkovTreeData*
   CiHmot_wrapper_simulate_size(const HiddenMarkovIndOutTree& hmt,
                                int inb_trees,
                                int isize, int inb_children_max,
                                bool counting_flag= true)
   {
      StatError error;
      HiddenMarkovTreeData* markov_data= NULL;
      ostringstream error_message;

      markov_data= hmt.simulation(error, inb_trees, isize, inb_children_max,
                                  counting_flag);
      if (markov_data == NULL)
      {
         error_message << error;
         throw_stat_tree_error(error_message);
      }
      return markov_data;
   }

   static HiddenMarkovTreeData*
   CiHmot_wrapper_simulate_trees(const HiddenMarkovIndOutTree& hmt,
                                 int inb_trees,
                                 const Trees& otrees,
                                 bool counting_flag= true)
   {
      StatError error;
      HiddenMarkovTreeData* markov_data= NULL;
      ostringstream error_message;

      markov_data= hmt.simulation(error, inb_trees, otrees, counting_flag);
      if (markov_data == NULL)
      {
         error_message << error;
         throw_stat_tree_error(error_message);
      }
      return markov_data;
   }

   static HiddenMarkovTreeData*
   CiHmot_wrapper_simulate_forest(const HiddenMarkovIndOutTree& hmt,
                                  const Trees& otrees,
                                  bool counting_flag= true)
   {
      HiddenMarkovTreeData* markov_data= NULL;

      markov_data= hmt.simulation(otrees, counting_flag);

      return markov_data;
   }

   static boost::python::list CiHmot_wrapper_state_profile5(const HiddenMarkovIndOutTree& hmt,
                                                            int viterbi_algorithm,
                                                            int nb_state_trees,
                                                            int index,
                                                            int entropy_algorithm,
                                                            int root)
   {
      typedef HiddenMarkovTreeData::tree_type tree_type;
      // typedef tree_type::vertex_descriptor vid;
      bool status= true;
      unsigned int cpt_msg= 0;
      StatError error;
      boost::python::list res;
      str msg;
      ostringstream error_message;
      tree_type *t= NULL;
      std::vector<ostringstream*> messages;
      HiddenMarkovTreeData *markov_data= NULL, *smoothed= NULL, *nstate_trees= NULL,
                           *viterbi_upward_downward= NULL,
                           *generalized_restoration= NULL;

      markov_data= hmt.extract_data(error);
      if (markov_data == NULL)
      {
         error_message << error;
         throw_stat_tree_error(error_message);
      }
      else
      {
         if (index == I_DEFAULT)
         {
            error_message << "Bad value for argument 'TreeId'" << endl;
            status= false;
         }
         if ((index < 0) || (index >= markov_data->get_nb_trees()))
         {
            error_message << "Bad value for argument 'TreeId': " << index << endl;
            status= false;
         }
         if ((status) && (root != I_DEFAULT) &&
             ((root < 0) || (root >= markov_data->get_size(index))))
         {
            error_message << "Bad value for argument 'RootVertex': " << root << endl;
            status= false;
         }
         if (status)
         {
            t= markov_data->get_tree(index);

            if (!(t->is_root(root)) && (root != I_DEFAULT)
                 && (t->get_nb_children(root) < 2))
            {
               error_message << "Bad number of children for subtree rooted at: "
                             << root << endl;
               status= false;
            }
            delete t;
            t= NULL;
         }
         if (!status)
         {
            delete markov_data;
            throw_python_error(PyExc_ValueError, error_message);
         }
         else
            status= hmt.state_profile(error, *markov_data, index, smoothed,
                                      nstate_trees, viterbi_upward_downward,
                                      generalized_restoration, messages, viterbi_algorithm,
                                      nb_state_trees, entropy_algorithm, root);

         if (!status)
         // if ((smoothed == NULL) || (viterbi_upward_downward == NULL)
         //     || (generalized_restoration == NULL)) // || (nstate_trees == NULL)
         {
            error_message << error;
            delete markov_data;
            if (smoothed != NULL)
            {
               delete smoothed;
               smoothed = NULL;
            }
            if (nstate_trees != NULL)
            {
               delete nstate_trees;
               nstate_trees= NULL;
            }
            if (viterbi_upward_downward != NULL)
            {
               delete viterbi_upward_downward;
               viterbi_upward_downward= NULL;
            }
            if (generalized_restoration != NULL)
            {
               delete generalized_restoration;
               generalized_restoration= NULL;
            }
            throw_stat_tree_error(error_message);
         }
         else
         {
            msg= (messages[cpt_msg]->str()).c_str();
            delete messages[cpt_msg];
            messages[cpt_msg++]= NULL;
            res.append(msg);
            res.append(smoothed);
            // messages for smoothed probabilities
            while (len(msg) > 0)
            {
               msg= (messages[cpt_msg]->str()).c_str();
               delete messages[cpt_msg];
               messages[cpt_msg++]= NULL;
               res.append(msg);
            }
            // res.append(nstate_trees);
            // head messages for viterbi upward-downward
            while (len(msg) > 0)
            {
               msg= (messages[cpt_msg]->str()).c_str();
               delete messages[cpt_msg];
               messages[cpt_msg++]= NULL;
               res.append(msg);
            }
            msg= (messages[cpt_msg]->str()).c_str();
            delete messages[cpt_msg];
            messages[cpt_msg++]= NULL;
            res.append(msg);

            res.append(viterbi_upward_downward);

            // tail message for viterbi upward-downward
            msg= (messages[cpt_msg]->str()).c_str();
            delete messages[cpt_msg];
            messages[cpt_msg++]= NULL;
            res.append(msg);

            res.append(generalized_restoration);

            // tail messages for generalized state restoration
            while(cpt_msg < messages.size())
            {
               msg= (messages[cpt_msg]->str()).c_str();
               delete messages[cpt_msg];
               messages[cpt_msg++]= NULL;
               res.append(msg);
            }
         }
         delete markov_data;
         markov_data= NULL;
      }
      return res;
   }

   static boost::python::list CiHmot_wrapper_state_profile4(const HiddenMarkovIndOutTree& hmt,
                                                            int viterbi_algorithm,
                                                            int nb_state_trees,
                                                            int index,
                                                            int entropy_algorithm)
   { return CiHmot_wrapper_state_profile5(hmt, viterbi_algorithm, nb_state_trees,
                                         index, entropy_algorithm, I_DEFAULT); }


   static boost::python::list CiHmot_wrapper_state_profile3(const HiddenMarkovIndOutTree& hmt,
                                                            int viterbi_algorithm,
                                                            int nb_state_trees,
                                                            int index)
   { return CiHmot_wrapper_state_profile5(hmt, viterbi_algorithm, nb_state_trees,
                                         index, Stat_trees::UPWARD, I_DEFAULT); }

   static HiddenMarkovIndOutTree*
   CiHmot_wrapper_ascii_read(const char * path, int size= I_DEFAULT_TREE_SIZE,
                             bool counting_flag= true, double cumul_threshold= OCCUPANCY_THRESHOLD)
   {
      StatError error;
      HiddenMarkovTree* hmt= NULL;
      HiddenMarkovIndOutTree* hmot= NULL;
      ostringstream error_message;

      hmt= Stat_trees::hidden_markov_ind_out_tree_ascii_read(error, path, size,
                                                             counting_flag, cumul_threshold);
      if (hmt == NULL)
      {
         error_message << error;
         throw_stat_tree_error(error_message);
      }
      else
         hmot= new HiddenMarkovIndOutTree(*hmt);

      delete hmt;
      hmt= NULL;
      return hmot;
   }

   static boost::python::list CiHmot_wrapper_entropy_computation1(const HiddenMarkovIndOutTree& hmt,
                                                                  int index)
   {
      typedef double** double_array_2d;
      typedef Stat_trees::HiddenMarkovTree::double_array_3d double_array_3d;
      typedef Stat_trees::HiddenMarkovTree::double_array_4d double_array_4d;

      bool status = true;
      int t, i, j, nb_trees, nb_states;
      double state_tree_entropy, state_tree_entropy_lik, upward_entropy,
             marginal_entropy_value, entropy1, max_marginal_entropy;
      ostringstream error_message;
      boost::python::list res;
      StatError error;
      HiddenMarkovTreeData *markov_data= NULL;
      double *marginal_entropy = NULL;
      double_array_3d state_marginal = NULL, output_cond = NULL,
                      upward_prob = NULL, upward_parent_prob = NULL,
                      state_entropy = NULL, downward_prob = NULL;
      double_array_4d downward_pair_prob = NULL;

      markov_data = hmt.extract_data(error);
      if (markov_data == NULL)
      {
         error_message << error;
         throw_stat_tree_error(error_message);
      }
      else
      {
         if ((index != I_DEFAULT) &&
             ((index < 0) || (index >= markov_data->get_nb_trees())))
         {
            error_message << "Bad value for argument 'TreeId': " << index << endl;
            status= false;
         }
         if (!status)
         {
            delete markov_data;
            throw_python_error(PyExc_ValueError, error_message);
         }
         else
         {
            nb_trees = markov_data->get_nb_trees();
            nb_states = markov_data->get_nb_states();
            hmt.get_state_marginal_distribution(*markov_data, state_marginal);
            hmt.get_output_conditional_distribution(*markov_data, output_cond);

            hmt.get_upward_step(*markov_data, upward_prob, upward_parent_prob,
                                state_entropy, state_marginal, output_cond,
                                state_tree_entropy, index);

            hmt.get_downward_step(*markov_data, downward_prob, downward_pair_prob,
                                  upward_prob, upward_parent_prob, state_marginal,
                                  output_cond, state_tree_entropy_lik, index);

            if (index != I_DEFAULT)
            {
               marginal_entropy_value = hmt.get_marginal_entropy(*markov_data, index,
                                                                 downward_prob, marginal_entropy,
                                                                 max_marginal_entropy);
            }
            else
            {
               marginal_entropy_value = 0.;
               for(t = 0; t < nb_trees; t++)
               {
                  marginal_entropy_value += hmt.get_marginal_entropy(*markov_data, t,
                                                                     downward_prob,
                                                                     marginal_entropy,
                                                                     max_marginal_entropy);
                  delete [] marginal_entropy;
                  marginal_entropy = NULL;
               }
            }

            for(t = 0; t < nb_trees; t++)
            {
               for(j = 0; j < nb_states; j++)
               {
                  delete [] state_marginal[t][j];
                  delete [] upward_prob[t][j];
                  delete [] upward_parent_prob[t][j];
                  delete [] downward_prob[t][j];
                  for(i= 0; i < nb_states; i++)
                     delete [] downward_pair_prob[t][j][i];
                  delete [] downward_pair_prob[t][j];
                  delete [] state_entropy[t][j];
                  delete [] output_cond[t][j];
               }
               delete [] state_marginal[t];
               delete [] output_cond[t];
               delete [] upward_prob[t];
               delete [] upward_parent_prob[t];
               delete [] downward_prob[t];
               delete [] downward_pair_prob[t];
               delete [] state_entropy[t];
            }

            delete [] state_marginal;
            delete [] output_cond;
            delete [] upward_prob;
            delete [] upward_parent_prob;
            delete [] downward_prob;
            delete [] downward_pair_prob;
            delete [] state_entropy;

            delete markov_data;
            markov_data = NULL;

            res.append(state_tree_entropy);
            res.append(marginal_entropy_value);
         }
      }
      return res;
   }

   static boost::python::list CiHmot_wrapper_entropy_computation0(const HiddenMarkovIndOutTree& hmt)
   { return CiHmot_wrapper_entropy_computation1(hmt, I_DEFAULT); }

   static boost::python::list CiHmot_wrapper_upward_entropy_computation1(const HiddenMarkovIndOutTree& hmt,
                                                                         int index)
   {
      typedef double** double_array_2d;
      typedef Stat_trees::HiddenMarkovTree::double_array_3d double_array_3d;
      typedef Stat_trees::HiddenMarkovTree::double_array_4d double_array_4d;

      bool status = true;
      int t, i, j, nb_trees, nb_states;
      double state_tree_entropy, state_tree_entropy2, upward_entropy,
             marginal_entropy_value, entropy1, max_marginal_entropy;
      ostringstream error_message;
      boost::python::list res;
      StatError error;
      HiddenMarkovTreeData *markov_data= NULL;
      double *marginal_entropy = NULL;
      double_array_2d conditional_upward_entropy = NULL;
      double_array_3d state_marginal = NULL, output_cond = NULL,
                      upward_prob = NULL, upward_parent_prob = NULL,
                      state_entropy = NULL, downward_prob = NULL;
      double_array_4d downward_pair_prob = NULL;

      markov_data = hmt.extract_data(error);
      if (markov_data == NULL)
      {
         error_message << error;
         throw_stat_tree_error(error_message);
      }
      else
      {
         if ((index != I_DEFAULT) &&
             ((index < 0) || (index >= markov_data->get_nb_trees())))
         {
            error_message << "Bad value for argument 'TreeId': " << index << endl;
            status= false;
         }
         if (!status)
         {
            delete markov_data;
            throw_python_error(PyExc_ValueError, error_message);
         }
         else
         {
            nb_trees = markov_data->get_nb_trees();
            nb_states = markov_data->get_nb_states();
            hmt.get_state_marginal_distribution(*markov_data, state_marginal);
            hmt.get_output_conditional_distribution(*markov_data, output_cond);

            hmt.get_upward_step(*markov_data, upward_prob, upward_parent_prob,
                                state_entropy, state_marginal, output_cond,
                                state_tree_entropy, index);
            hmt.get_downward_step(*markov_data, downward_prob, downward_pair_prob,
                                  upward_prob, upward_parent_prob, state_marginal,
                                  output_cond, state_tree_entropy2, index);

            upward_entropy = hmt.get_upward_conditional_entropy(*markov_data, state_marginal,
                                                                upward_prob, upward_parent_prob,
                                                                downward_prob,
                                                                conditional_upward_entropy,
                                                                index);
            if (index != I_DEFAULT)
            {
               marginal_entropy_value = hmt.get_marginal_entropy(*markov_data, index,
                                                                 downward_prob, marginal_entropy,
                                                                 max_marginal_entropy);
            }
            else
            {
               marginal_entropy_value = 0.;
               for(t = 0; t < nb_trees; t++)
               {
                  marginal_entropy_value += hmt.get_marginal_entropy(*markov_data, t,
                                                                     downward_prob,
                                                                     marginal_entropy,
                                                                     max_marginal_entropy);
                  delete [] marginal_entropy;
                  marginal_entropy = NULL;
               }
            }

            for(t = 0; t < nb_trees; t++)
            {
               for(j = 0; j < nb_states; j++)
               {
                  delete [] state_marginal[t][j];
                  delete [] upward_prob[t][j];
                  delete [] upward_parent_prob[t][j];
                  delete [] downward_prob[t][j];
                  for(i= 0; i < nb_states; i++)
                     delete [] downward_pair_prob[t][j][i];
                  delete [] downward_pair_prob[t][j];
                  delete [] state_entropy[t][j];
                  delete [] output_cond[t][j];
               }
               delete [] state_marginal[t];
               delete [] output_cond[t];
               delete [] upward_prob[t];
               delete [] upward_parent_prob[t];
               delete [] downward_prob[t];
               delete [] downward_pair_prob[t];
               delete [] conditional_upward_entropy[t];
               delete [] state_entropy[t];
            }

            delete [] state_marginal;
            delete [] output_cond;
            delete [] upward_prob;
            delete [] upward_parent_prob;
            delete [] downward_prob;
            delete [] downward_pair_prob;
            delete [] state_entropy;
            delete [] conditional_upward_entropy;

            delete markov_data;
            markov_data = NULL;

            res.append(state_tree_entropy);
            res.append(upward_entropy);
            res.append(marginal_entropy_value);
         }
      }
      return res;
   }

   static boost::python::list CiHmot_wrapper_upward_entropy_computation0(const HiddenMarkovIndOutTree& hmt)
   { return CiHmot_wrapper_upward_entropy_computation1(hmt, I_DEFAULT); }

   static HiddenMarkovTreeData*
   CiHmot_wrapper_compute_state_trees(const HiddenMarkovIndOutTree& hmt,
                                      const Trees& trees,
                                      int algorithm,
                                      bool characteristic_flag)
   {
      double gini_index, information;
      StatError error;
      HiddenMarkovTreeData* markov_data= NULL;
      ostringstream error_message;

      markov_data= hmt.state_tree_computation(error, trees, algorithm, characteristic_flag);
      if (markov_data == NULL)
      {
         error_message << error;
         throw_stat_tree_error(error_message);
      }
      return markov_data;
   }

   static void CiHmot_wrapper_plot_write(const HiddenMarkovIndOutTree& hmt,
                                         const char* prefix,
                                         const char* title)
   {
      bool status= true;
      ostringstream error_message;
      StatError error;

      status= hmt.plot_write(error, prefix, title);
      if (not status)
      {
         error_message << error;
         throw_stat_tree_error(error_message);
      }
   }

   static void CiHmot_wrapper_state_profile_plot_write(const HiddenMarkovIndOutTree& hmt,
                                                       const char* prefix,
                                                       const char* title,
                                                       int identifier, int vertex,
                                                       int entropy_algorithm)
   {
      bool status= true;
      ostringstream error_message;
      StatError error;

      status= hmt.state_profile_plot_write(error, prefix, identifier, vertex,
                                           title, entropy_algorithm);
      if (not status)
      {
         error_message << error;
         throw_stat_tree_error(error_message);
      }
   }

   static HiddenMarkovIndOutTree*
   CiHmot_wrapper_ascii_read1(const char * path)
   { return CiHmot_wrapper_ascii_read(path); }

}; // end WRAP

void class_hmiot()
{
    class_< HiddenMarkovIndOutTree, bases<HiddenMarkovTree> >
    ("CiHmot", init< const HiddenMarkovIndOutTree&, optional< bool, bool> >())
        // .def("AsciiWrite", &CiHmot_wrapper_ascii_write0)
        .def("ExtractData", WRAP::CiHmot_wrapper_extract_data,
                            return_value_policy< manage_new_object >())
        .def("Display", WRAP::CiHmot_wrapper_display)
        .def("ComputeStateTrees", WRAP::CiHmot_wrapper_compute_state_trees,
                                  return_value_policy< manage_new_object >())
        .def("EntropyComputation", WRAP::CiHmot_wrapper_entropy_computation1,
                                  "EntropyComputation(self, int) -> list \n\n"
                                  "Return state tree entropy "
                                  "and sum of marginal entropies "
                                  "for a given tree")
        .def("EntropyComputation", WRAP::CiHmot_wrapper_entropy_computation0,
                                  "EntropyComputation(self) -> list \n\n"
                                  "Return state tree entropy "
                                  "and sum of marginal entropies "
                                  "summed over every tree")

        .def("FileAsciiWrite", WRAP::CiHmot_wrapper_file_ascii_write1)
        .def("FileAsciiWrite", WRAP::CiHmot_wrapper_file_ascii_write2)
        .def("plot_write", WRAP::CiHmot_wrapper_plot_write,
                       "Write into a gnuplot file.")
        .def("Simulate", WRAP::CiHmot_wrapper_simulate_histo,
                         return_value_policy< manage_new_object >())
        .def("Simulate", WRAP::CiHmot_wrapper_simulate_size,
                         return_value_policy< manage_new_object >())
        .def("Simulate", WRAP::CiHmot_wrapper_simulate_trees,
                         return_value_policy< manage_new_object >())
        .def("Simulate", WRAP::CiHmot_wrapper_simulate_forest,
                         return_value_policy< manage_new_object >())
        .def("SpreadsheetWrite", WRAP::CiHmot_wrapper_spreadsheet_write1)
        .def("StatePermutation", WRAP::CiHmot_wrapper_state_permutation,
                                 "StatePermutation(self, list) \n\n"
                                 "permutation of the model states\n")
        .def("StateProfile", WRAP::CiHmot_wrapper_state_profile5,
                             "StateProfile(self, int, int, int, int, int) -> list \n\n"
                             "return trees object and strings"
                             "for the state tree analysis\n")
        .def("StateProfile", WRAP::CiHmot_wrapper_state_profile4,
                             "StateProfile(self, int, int, int, int) -> list \n\n"
                             "return trees object and strings"
                             "for the state tree analysis\n")
        .def("StateProfile", WRAP::CiHmot_wrapper_state_profile3,
                             "StateProfile(self, int, int, int) -> list \n\n"
                             "return trees object and strings"
                             "for the state tree analysis\n")
        .def("UpwardEntropyComputation", WRAP::CiHmot_wrapper_upward_entropy_computation1,
                                         "UpwardEntropyComputation(self, int) -> list \n\n"
                                         "Return state tree entropy, sum of upward "
                                         "entropies and sum of marginal entropies "
                                         "for a given tree")
        .def("UpwardEntropyComputation", WRAP::CiHmot_wrapper_upward_entropy_computation0,
                                         "UpwardEntropyComputation(self) -> list \n\n"
                                         "Return state tree entropy, sum of upward "
                                         "entropies and sum of marginal entropies "
                                         "summed over every tree")
        .def("state_profile_plot_write", WRAP::CiHmot_wrapper_state_profile_plot_write,
                                 "state_profile_plot_write(self, str, str, int, int, int) -> void \n\n"
                                 "write Gnuplot files for state and entropy"
                                 "profiles for the state tree analysis\n")
    ;

    def("HmtAsciiRead", WRAP::CiHmot_wrapper_ascii_read1,
                        return_value_policy< manage_new_object >());

};
#undef WRAP

/*************************************************************
 *
 *  Wrappers for Python class CHmt_data:
 */

#define WRAP HiddenMarkovTreeDataWrap
class WRAP
{

public :

   static bool Chmt_data_wrapper_is_parametric(const HiddenMarkovTreeData& reftree,
                                               int variable)
   {
      ostringstream error_message;
      HiddenMarkovTree* hmarkovt= reftree.get_markov();
      bool res= false;

      if (hmarkovt == NULL)
      {
         error_message << "No model associated with data";
         throw_stat_tree_error(error_message);
      }
      else
      {
         if ((variable <= 0) ||
             (variable > reftree.get_nb_int() + reftree.get_nb_float()))
         {
            error_message << STAT_TREES_error[TREESTATR_OUTPUT_PROCESS_INDEX];
            throw_python_error(PyExc_IndexError, error_message);
         }
         else
         {
            res= hmarkovt->is_parametric(variable-1);
            delete hmarkovt;
            hmarkovt= NULL;
            return res;
         }
      }
   }

   static DiscreteMixtureData* Chmt_data_wrapper_extract_marginal(const HiddenMarkovTreeData& reftree,
                                                                  int variable)
   {
      ostringstream error_message;
      StatError error;
      DiscreteMixtureData *histo = NULL;

      histo = reftree.extract_marginal(error, variable);
      if (histo == NULL)
      {
         error_message << error;
         throw_stat_tree_error(error_message);
      }
      return histo;
   }

   static HiddenMarkovTree*
   Chmt_data_wrapper_extract_model(const HiddenMarkovTreeData& tree)
   {
      HiddenMarkovTree *res;
      StatError error;
      ostringstream error_message;

      res = tree.extract_model(error);
      if (res == NULL)
      {
         error_message << error;
         throw_stat_tree_error(error_message);
      }
      return res;
   }

   static HiddenMarkovIndOutTree*
   Chmt_data_wrapper_hidden_markov_ind_out_tree_estimation_markov(const HiddenMarkovTreeData& hmtd,
                                                              const HiddenMarkovIndOutTree& ihmarkov,
                                                              bool counting_flag,
                                                              int state_trees,
                                                              int algorithm,
                                                              double saem_exponent,
                                                              int nb_iter,
                                                              bool force_param)
   {
      HiddenMarkovIndOutTree *hmt= NULL;
      StatError error;
      ostringstream error_message;

      hmt= hmtd.hidden_markov_ind_out_tree_estimation(error, cout, ihmarkov,
                                                  counting_flag, state_trees,
                                                  algorithm, saem_exponent,
                                                  nb_iter, force_param);
      if (hmt == NULL)
      {
         error_message << error;
         throw_stat_tree_error(error_message);
      }
      return hmt;
   }

   static HiddenMarkovIndOutTree*
   Chmt_data_wrapper_hidden_markov_ind_out_tree_estimation_self_transition(const HiddenMarkovTreeData& hmtd,
                                                                       int nb_state,
                                                                       bool left_right,
                                                                       bool counting_flag,
                                                                       int state_trees,
                                                                       int algorithm,
                                                                       double saem_exponent,
                                                                       double self_transition,
                                                                       int nb_iter,
                                                                       boost::python::list force_param)
   {
      bool status= true, several_errors= false;
      const int nb_variables= hmtd.get_nb_int() + hmtd.get_nb_float();
      int nb_fparam, p;
      char type= 'o'; // or may be not ...
      ostringstream error_message;
      object o;
      StatError error;
      bool *fparam= NULL;
      HiddenMarkovIndOutTree *hmt= NULL;

      nb_fparam= boost::python::len(force_param);
      if (nb_fparam > 0)
      {
         if (nb_fparam != nb_variables)
         {
            status= false;
            error_message << "bad size of argument list: " << nb_fparam
                          << ": should be the number of variables ("
                          << nb_variables << ")";
            throw_python_error(PyExc_ValueError, error_message);
         }
         else
         {
            fparam= new bool[nb_fparam];
            for (p= 0; p < nb_fparam; p++)
            {
               o= force_param[p];
               try
               {
                  extract<bool> x(o);
                  if (x.check())
                     fparam[p]= x();
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
                  error_message << "incorrect type for element " << p
                                << " of argument list: expecting a boolean"; // << endl;
               }
            }
            if (!status)
            {
               delete [] fparam;
               fparam= NULL;
               throw_python_error(PyExc_TypeError, error_message);
            }
         }
      }
      if (status)
      {
         hmt= hmtd.hidden_markov_ind_out_tree_estimation(error, cout, type, nb_state,
                                                     left_right, counting_flag,
                                                     state_trees, algorithm,
                                                     saem_exponent, self_transition,
                                                     nb_iter, fparam);

         if (hmt == NULL)
         {
            error_message << error;
            throw_stat_tree_error(error_message);
         }
         if (fparam != NULL)
         {
            delete [] fparam;
            fparam= NULL;
         }
      }
      return hmt;
   }

   static str CHmt_data_wrapper_ascii_write0(const HiddenMarkovTreeData& hmtd)
   {
      std::stringstream ss;
      str res;

      hmtd.line_write(ss);
      res= str(ss.str());

      return res;
   }

   static void CHmt_data_wrapper_plot_write(const HiddenMarkovTreeData& reftree,
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

   static DiscreteDistributionData* Chmt_data_wrapper_extract_value(const HiddenMarkovTreeData& reftree,
                                                                    int type,
                                                                    int variable,
                                                                    int value)
   {
      ostringstream error_message;
      // bool status= true;
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

}; // end WRAP


void class_hmt_data()
{
    class_< HiddenMarkovTreeData, bases<Trees> >
    ("CHmt_data", init< const HiddenMarkovTreeData&, optional< bool> >())
        .def(init< const Trees&>())
        // .def(init< const Trees&, int>())
        // .def("__init__", make_constructor(CHmt_data_wrapper_init1))
        .def("IsParametric", WRAP::Chmt_data_wrapper_is_parametric,
                            "IsParametric(self, variable) -> bool \n\n"
                            "Return True if process 'variable' "
                            "is parametric")
        .def("ExtractMarginal", WRAP::Chmt_data_wrapper_extract_marginal,
             return_value_policy< manage_new_object>(),
             "ExtractMarginal(self, variable) -> DiscreteMixtureData \n\n"
             "Return the mixture of observation distributions"
             "for a given variable")
        .def("ExtractMarkov", WRAP::Chmt_data_wrapper_extract_model,
             return_value_policy< manage_new_object>(),
             "ExtractMarginal(self) -> HiddenMarkovTree \n\n"
             "Return the model part of self \n")
        .def("StateTrees",
             &HiddenMarkovTreeData::get_state_hidden_markov_tree_data,
             return_value_policy< manage_new_object >(),
             "StateTrees(self) -> CHmt_data \n\n"
             "Return a CHmt_data containing the states as a variable \n")
        .def("SmoothedProbaTrees",
             &HiddenMarkovTreeData::get_state_smoothed_hidden_markov_tree_data,
             return_value_policy< manage_new_object >(),
             "SmoothedProbaTrees(self, index) -> CHmt_data \n\n"
             "For a given tree, return a CHmt_data containing the states "
             "and the smoothed probabilities as variables \n")
        .def("EstimationCiHmot",
             WRAP::Chmt_data_wrapper_hidden_markov_ind_out_tree_estimation_markov,
             return_value_policy< manage_new_object >())
        .def("EstimationCiHmot",
             WRAP::Chmt_data_wrapper_hidden_markov_ind_out_tree_estimation_self_transition,
             return_value_policy< manage_new_object >())
        .def("ExtractValueHistogram", WRAP::Chmt_data_wrapper_extract_value,
             return_value_policy< manage_new_object>())
        .def("plot_write", WRAP::CHmt_data_wrapper_plot_write)
        .def("__str__", WRAP::CHmt_data_wrapper_ascii_write0)
    ;
};



