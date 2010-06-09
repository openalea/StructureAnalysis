// Includes ====================================================================
#include "tree/basic_visitors.h"
#include "tree/tree_traits.h"
#include "tree/tree_simple.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/distribution.h"
#include "stat_tool/vectors.h"

#include "tree_statistic/tree_labels.h"
#include "tree_statistic/int_fl_containers.h"
#include "tree_statistic/generic_typed_edge_tree.h"
#include "tree_statistic/typed_edge_trees.h"
#include "tree_statistic/hidden_markov_tree.h"
#include "tree_statistic/hidden_markov_out_tree.h"

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
namespace  {
//
// BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(CiHmot_wrapper_ascii_write_overloads_1_2, CiHmot_wrapper_ascii_write, 1, 2)
//

template<int num> struct UniqueInt { int v; enum { value=num };
UniqueInt(int _v) : v(_v) { } operator int() const { return v; } };

/*************************************************************
 *
 *  Exporting constants as constant functions
 */

// int NB_TREES_() { return NB_TREES; }

/*************************************************************
 *
 *  Wrappers for Python class CHmt:
 */

int CHmt_wrapper_nb_values(const HiddenMarkovTree& hmt, int variable)
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

HiddenMarkovTreeData*
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

DiscreteParametricModel*
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

double CHmt_wrapper_likelihood(const HiddenMarkovTree& hmt,
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

str CHmt_wrapper_ascii_write0(const HiddenMarkovTree& hmt)
{
   std::stringstream s;
   str res;

   hmt.line_write(s);
   res= str(s.str());
   return res;
}

/*************************************************************
 *
 *  Wrappers for Python class CiHmot:
 */

HiddenMarkovTreeData*
CiHmot_wrapper_extract_data(const HiddenMarkovOutTree& hmt)
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

str CiHmot_wrapper_display(const HiddenMarkovOutTree& hmt,
                           bool exhaustive)
{
   std::stringstream s;
   str res;

   hmt.ascii_write(s, exhaustive);
   res= str(s.str());
   return res;
}

bool CiHmot_wrapper_file_ascii_write2(const HiddenMarkovOutTree& hmt,
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

bool CiHmot_wrapper_file_ascii_write1(const HiddenMarkovOutTree& hmt,
                                      const char * path)
{ return CiHmot_wrapper_file_ascii_write2(hmt, path, false); }

void CiHmot_wrapper_state_permutation(const HiddenMarkovOutTree& hmt,
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

bool CiHmot_wrapper_spreadsheet_write1(const HiddenMarkovOutTree& hmt,
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

HiddenMarkovTreeData*
CiHmot_wrapper_simulate_histo(const HiddenMarkovOutTree& hmt,
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

HiddenMarkovTreeData*
CiHmot_wrapper_simulate_size(const HiddenMarkovOutTree& hmt,
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

HiddenMarkovTreeData*
CiHmot_wrapper_simulate_trees(const HiddenMarkovOutTree& hmt,
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

boost::python::list CiHmot_wrapper_state_profile5(const HiddenMarkovOutTree& hmt,
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

         if (t->get_nb_children(root) < 2)
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
         throw_python_error(PyExc_ValueError, error_message);
      }
      else
      {
      status= hmt.state_profile(error, *markov_data, index, smoothed,
                                nstate_trees, viterbi_upward_downward,
                                generalized_restoration, messages, viterbi_algorithm,
                                nb_state_trees, entropy_algorithm, root);
      }

      if (!status)
      // if ((smoothed == NULL) || (viterbi_upward_downward == NULL)
      //     || (generalized_restoration == NULL)) // || (nstate_trees == NULL)
      {
         error_message << error;
         throw_stat_tree_error(error_message);
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

boost::python::list CiHmot_wrapper_state_profile4(const HiddenMarkovOutTree& hmt,
                                                 int viterbi_algorithm,
                                                 int nb_state_trees,
                                                 int index,
                                                 int entropy_algorithm)
{ return CiHmot_wrapper_state_profile5(hmt, viterbi_algorithm, nb_state_trees,
                                      index, entropy_algorithm, I_DEFAULT); }


boost::python::list CiHmot_wrapper_state_profile3(const HiddenMarkovOutTree& hmt,
                                                 int viterbi_algorithm,
                                                 int nb_state_trees,
                                                 int index)
{ return CiHmot_wrapper_state_profile5(hmt, viterbi_algorithm, nb_state_trees,
                                      index, Stat_trees::UPWARD, I_DEFAULT); }

HiddenMarkovOutTree*
Hmt_wrapper_ascii_read(const char * path, int size= I_DEFAULT_TREE_SIZE,
                       bool counting_flag= true, double cumul_threshold= OCCUPANCY_THRESHOLD)
{
   StatError error;
   HiddenMarkovTree* hmt= NULL;
   HiddenMarkovOutTree* hmot= NULL;
   ostringstream error_message;

   hmt= Stat_trees::hidden_markov_tree_ascii_read(error, path, size,
                                                  counting_flag, cumul_threshold);
   if (hmt == NULL)
   {
      error_message << error;
      throw_stat_tree_error(error_message);
   }
   else
      hmot= new HiddenMarkovOutTree(*hmt);

   delete hmt;
   hmt= NULL;
   return hmot;
}

HiddenMarkovTreeData*
CiHmot_wrapper_compute_state_trees(const HiddenMarkovOutTree& hmt,
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

void CiHmot_wrapper_plot_write(const HiddenMarkovOutTree& hmt,
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

void CiHmot_wrapper_state_profile_plot_write(const HiddenMarkovOutTree& hmt,
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

HiddenMarkovOutTree*
Hmt_wrapper_ascii_read1(const char * path)
{ return Hmt_wrapper_ascii_read(path); }


/*************************************************************
 *
 *  Wrappers for Python class CHmt_data:
 */

bool Chmt_data_wrapper_is_parametric(const HiddenMarkovTreeData& reftree,
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
         error_message << STAT_TREES_error[STATR_OUTPUT_PROCESS_INDEX];
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

HiddenMarkovOutTree*
Chmt_data_wrapper_hidden_markov_out_tree_estimation_markov(const HiddenMarkovTreeData& hmtd,
                                                           const HiddenMarkovOutTree& ihmarkov,
                                                           bool counting_flag,
                                                           int state_trees,
                                                           int algorithm,
                                                           double saem_exponent,
                                                           int nb_iter,
                                                           bool force_param)
{
   HiddenMarkovOutTree *hmt= NULL;
   StatError error;
   ostringstream error_message;

   hmt= hmtd.hidden_markov_out_tree_estimation(error, cout, ihmarkov,
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

HiddenMarkovOutTree*
Chmt_data_wrapper_hidden_markov_out_tree_estimation_self_transition(const HiddenMarkovTreeData& hmtd,
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
   HiddenMarkovOutTree *hmt= NULL;

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
      hmt= hmtd.hidden_markov_out_tree_estimation(error, cout, type, nb_state,
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

str CHmt_data_wrapper_ascii_write0(const HiddenMarkovTreeData& hmtd)
{
   std::stringstream ss;
   str res;

   hmtd.line_write(ss);
   res= str(ss.str());

   return res;
}

void CHmt_data_wrapper_plot_write(const HiddenMarkovTreeData& reftree,
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

DiscreteDistributionData* Chmt_data_wrapper_extract_value(const HiddenMarkovTreeData& reftree,
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

}// namespace


// Module ======================================================================
BOOST_PYTHON_MODULE(chmt)
{
    // Error initialisation
    object stat_tree_errors = import("openalea.tree_statistic._errors");
    // Import StatError
    StatTreeError = stat_tree_errors.attr("StatTreeError");

    class_< HiddenMarkovTree >
    ("CHmt", init< const HiddenMarkovTree&, optional< bool, bool> >())
        .def("ExtractData", &CHmt_wrapper_extract_data,
                            return_value_policy< manage_new_object >())
        .def("Extract", &CHmt_wrapper_extract,
                        return_value_policy< manage_new_object >(),
                        "Extract(self, type, variable, value) -> _DiscreteParametricModel \n\n"
                        "Extract some characteristic or observation distribution")
        .def("IsParametric", &HiddenMarkovTree::is_parametric,
                            "IsParametric(self, variable) -> bool \n\n"
                            "Return True if process 'variable' "
                            "is parametric")
        .def("Likelihood", &CHmt_wrapper_likelihood,
                            "Likelihood(self, trees) -> float \n\n"
                            "Return the likelihood of the parameters"
                            "for the given trees")
        .def("NbInt", &HiddenMarkovTree::get_nb_ioutput_process)
        .def("NbFloat", &HiddenMarkovTree::get_nb_doutput_process)
        .def("NbStates", &HiddenMarkovTree::get_nb_state,
                         "NbStates(self) -> int \n\n"
                         "return the number of states "
                         "of the hidden Markov tree\n")
        .def("NbValues", &CHmt_wrapper_nb_values,
                         "NbValues(self) -> int \n\n"
                         "return the number of values for a given process "
                         "of the hidden Markov tree\n")
        .def("__str__", &CHmt_wrapper_ascii_write0,
                        "__str__(self) -> str \n\n"
                        "Display the hidden Markov tree.")
    ;

    class_< HiddenMarkovOutTree, bases<HiddenMarkovTree> >
    ("CiHmot", init< const HiddenMarkovOutTree&, optional< bool, bool> >())
        // .def("AsciiWrite", &CiHmot_wrapper_ascii_write0)
        .def("ExtractData", &CiHmot_wrapper_extract_data,
                            return_value_policy< manage_new_object >())
        .def("Display", &CiHmot_wrapper_display)
        .def("ComputeStateTrees", &CiHmot_wrapper_compute_state_trees,
                                  return_value_policy< manage_new_object >())
        .def("FileAsciiWrite", &CiHmot_wrapper_file_ascii_write1)
        .def("FileAsciiWrite", &CiHmot_wrapper_file_ascii_write2)
        .def("plot_write", &CiHmot_wrapper_plot_write,
                       "Write into a gnuplot file.")
        .def("Simulate", &CiHmot_wrapper_simulate_histo,
                         return_value_policy< manage_new_object >())
        .def("Simulate", &CiHmot_wrapper_simulate_size,
                         return_value_policy< manage_new_object >())
        .def("Simulate", &CiHmot_wrapper_simulate_trees,
                         return_value_policy< manage_new_object >())
        .def("SpreadsheetWrite", &CiHmot_wrapper_spreadsheet_write1)
        .def("StatePermutation", &CiHmot_wrapper_state_permutation,
                                 "StatePermutation(self, list) \n\n"
                                 "permutation of the model states\n")
        .def("StateProfile", &CiHmot_wrapper_state_profile5,
                             "StateProfile(self, int, int, int, int, int) -> list \n\n"
                             "return trees object and strings"
                             "for the state tree analysis\n")
        .def("StateProfile", &CiHmot_wrapper_state_profile4,
                             "StateProfile(self, int, int, int, int) -> list \n\n"
                             "return trees object and strings"
                             "for the state tree analysis\n")
        .def("StateProfile", &CiHmot_wrapper_state_profile3,
                             "StateProfile(self, int, int, int) -> list \n\n"
                             "return trees object and strings"
                             "for the state tree analysis\n")
        .def("state_profile_plot_write", &CiHmot_wrapper_state_profile_plot_write,
                                 "state_profile_plot_write(self, str, str, int, int, int) -> void \n\n"
                                 "write Gnuplot files for state and entropy"
                                 "profiles for the state tree analysis\n")
    ;

    class_< HiddenMarkovTreeData, bases<Trees> >
    ("CHmt_data", init< const HiddenMarkovTreeData&, optional< bool> >())
        .def(init< const Trees&>())
        // .def(init< const Trees&, int>())
        // .def("__init__", make_constructor(CHmt_data_wrapper_init1))
        .def("IsParametric", &Chmt_data_wrapper_is_parametric,
                            "IsParametric(self, variable) -> bool \n\n"
                            "Return True if process 'variable' "
                            "is parametric")
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
             &Chmt_data_wrapper_hidden_markov_out_tree_estimation_markov,
             return_value_policy< manage_new_object >())
        .def("EstimationCiHmot",
             &Chmt_data_wrapper_hidden_markov_out_tree_estimation_self_transition,
             return_value_policy< manage_new_object >())
        .def("ExtractValueHistogram", &Chmt_data_wrapper_extract_value,
             return_value_policy< manage_new_object>())
        .def("plot_write", &CHmt_data_wrapper_plot_write)
        .def("__str__", &CHmt_data_wrapper_ascii_write0)
    ;

    def("HmtAsciiRead", Hmt_wrapper_ascii_read1,
                        return_value_policy< manage_new_object >());

    enum_<UniqueInt<2> >("EntropyAlgorithm")
        .value("UPWARD", Stat_trees::UPWARD)
        .value("DOWNWARD", Stat_trees::DOWNWARD)
        .export_values()
    ;
     // def("NB_TREES", NB_TREES_);

}
