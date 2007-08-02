// Includes ====================================================================
#include "tree/basic_visitors.h"
#include "tree/tree_traits.h"
#include "tree/tree_simple.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/distribution.h"

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

// Using =======================================================================
using namespace boost::python;
using namespace Stat_trees;

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

int CHmt_wrapper_nb_values(const Hidden_markov_tree& hmt, int variable)
{
   int res= 0;
   ostringstream error_message;

   if ((variable >= 0) && (variable < hmt.get_nb_ioutput_process()))
      res= hmt.get_nb_values(variable);
   else
   {
      error_message << "Bad variable: " << variable << endl;
      PyErr_SetString(PyExc_RuntimeError, (error_message.str()).c_str());
      throw_error_already_set();
   }
   return res;
}

Hidden_markov_tree_data*
CHmt_wrapper_extract_data(const Hidden_markov_tree& hmt)
{
   Hidden_markov_tree_data *res;
   Format_error error;
   ostringstream error_message;

   res= hmt.extract_data(error);
   if (res == NULL)
   {
      error_message << error;
      PyErr_SetString(PyExc_RuntimeError, (error_message.str()).c_str());
      throw_error_already_set();
   }
   return res;
}

double CHmt_wrapper_likelihood(const Hidden_markov_tree& hmt,
                               const Trees& trees)
{
   double res;
   ostringstream error_message;

   res= hmt.likelihood_computation(trees, I_DEFAULT);
   if (res <= D_INF)
   {
      error_message << "Cannot compute the likelihood";
      PyErr_SetString(PyExc_RuntimeError, (error_message.str()).c_str());
      throw_error_already_set();
   }
   return res;
}

str CHmt_wrapper_ascii_write0(const Hidden_markov_tree& hmt)
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

Hidden_markov_tree_data*
CiHmot_wrapper_extract_data(const Hidden_markov_out_tree& hmt)
{
   Hidden_markov_tree_data *res;
   Format_error error;
   ostringstream error_message;

   res= hmt.extract_data(error);
   if (res == NULL)
   {
      error_message << error;
      PyErr_SetString(PyExc_RuntimeError, (error_message.str()).c_str());
      throw_error_already_set();
   }
   return res;
}

str CiHmot_wrapper_display(const Hidden_markov_out_tree& hmt,
                           bool exhaustive)
{
   std::stringstream s;
   str res;

   hmt.ascii_write(s, exhaustive);
   res= str(s.str());
   return res;
}

bool CiHmot_wrapper_file_ascii_write2(const Hidden_markov_out_tree& hmt,
                                      const char * path, bool exhaustive)
{
   Format_error error;
   bool res;
   ostringstream error_message;

   res= hmt.ascii_write(error, path, exhaustive);
   if (!res)
   {
      error_message << error;
      PyErr_SetString(PyExc_RuntimeError, (error_message.str()).c_str());
      throw_error_already_set();
   }
   return res;
}

bool CiHmot_wrapper_file_ascii_write1(const Hidden_markov_out_tree& hmt,
                                      const char * path)
{ return CiHmot_wrapper_file_ascii_write2(hmt, path, false); }

bool CiHmot_wrapper_spreadsheet_write1(const Hidden_markov_out_tree& hmt,
                                       const char * path)
{
   Format_error error;
   bool res;
   ostringstream error_message;

   res= hmt.spreadsheet_write(error, path);
   if (!res)
   {
      error_message << error;
      PyErr_SetString(PyExc_RuntimeError, (error_message.str()).c_str());
      throw_error_already_set();
   }
   return res;
}

Hidden_markov_tree_data*
CiHmot_wrapper_simulate_histo(const Hidden_markov_out_tree& hmt,
                              const Histogram& ihsize,
                              const Histogram& ihnb_children,
                              bool counting_flag= true,
                              bool divergence_flag= false)
{
   Format_error error;
   Hidden_markov_tree_data* markov_data= NULL;
   ostringstream error_message;

   markov_data= hmt.simulation(error, ihsize, ihnb_children,
                               counting_flag, divergence_flag);
   if (markov_data == NULL)
   {
      error_message << error;
      PyErr_SetString(PyExc_RuntimeError, (error_message.str()).c_str());
      throw_error_already_set();
   }
   return markov_data;
}

Hidden_markov_tree_data*
CiHmot_wrapper_simulate_size(const Hidden_markov_out_tree& hmt,
                             int inb_trees,
                             int isize, int inb_children_max,
                             bool counting_flag= true)
{
   Format_error error;
   Hidden_markov_tree_data* markov_data= NULL;
   ostringstream error_message;

   markov_data= hmt.simulation(error, inb_trees, isize, inb_children_max,
                               counting_flag);
   if (markov_data == NULL)
   {
      error_message << error;
      PyErr_SetString(PyExc_RuntimeError, (error_message.str()).c_str());
      throw_error_already_set();
   }
   return markov_data;
}

Hidden_markov_tree_data*
CiHmot_wrapper_simulate_trees(const Hidden_markov_out_tree& hmt,
                              int inb_trees,
                              const Trees& otrees,
                              bool counting_flag= true)
{
   Format_error error;
   Hidden_markov_tree_data* markov_data= NULL;
   ostringstream error_message;

   markov_data= hmt.simulation(error, inb_trees, otrees, counting_flag);
   if (markov_data == NULL)
   {
      error_message << error;
      PyErr_SetString(PyExc_RuntimeError, (error_message.str()).c_str());
      throw_error_already_set();
   }
   return markov_data;
}

list CiHmot_wrapper_state_profile(const Hidden_markov_out_tree& hmt,
                                  int viterbi_algorithm,
                                  int nb_state_trees,
                                  int index,
                                  int entropy_algorithm= Stat_trees::UPWARD)
{
   bool status= false;
   unsigned int cpt_msg= 0;
   Format_error error;
   list res;
   str msg;
   ostringstream error_message;
   std::vector<ostringstream*> messages;
   Hidden_markov_tree_data *markov_data= NULL, *smoothed= NULL, *nstate_trees= NULL,
                           *viterbi_upward_downward= NULL,
                           *generalized_restoration= NULL;

   markov_data= hmt.extract_data(error);
   if (markov_data == NULL)
   {
      error_message << error;
      PyErr_SetString(PyExc_RuntimeError, (error_message.str()).c_str());
      throw_error_already_set();
   }
   else
   {
      status= hmt.state_profile(error, *markov_data, index, smoothed, nstate_trees,
                                viterbi_upward_downward, generalized_restoration,
                                messages, viterbi_algorithm, nb_state_trees,
                                entropy_algorithm);

      if (!status)
      // if ((smoothed == NULL) || (viterbi_upward_downward == NULL)
      //     || (generalized_restoration == NULL)) // || (nstate_trees == NULL)
      {
         error_message << error;
         PyErr_SetString(PyExc_RuntimeError, (error_message.str()).c_str());
         throw_error_already_set();
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

Hidden_markov_out_tree*
Hmt_wrapper_ascii_read(const char * path, int size= I_DEFAULT_TREE_SIZE,
                       bool counting_flag= true, double cumul_threshold= OCCUPANCY_THRESHOLD)
{
   Format_error error;
   Hidden_markov_tree* hmt= NULL;
   Hidden_markov_out_tree* hmot= NULL;
   ostringstream error_message;

   hmt= Stat_trees::hidden_markov_tree_ascii_read(error, path, size,
                                                  counting_flag, cumul_threshold);
   if (hmt == NULL)
   {
      error_message << error;
      PyErr_SetString(PyExc_RuntimeError, (error_message.str()).c_str());
      throw_error_already_set();
   }
   else
      hmot= new Hidden_markov_out_tree(*hmt);

   delete hmt;
   hmt= NULL;
   return hmot;
}

Hidden_markov_tree_data*
CiHmot_wrapper_compute_state_trees(const Hidden_markov_out_tree& hmt,
                                   const Trees& trees,
                                   int algorithm,
                                   bool characteristic_flag)
{
   double gini_index, information;
   Format_error error;
   Hidden_markov_tree_data* markov_data= NULL;
   ostringstream error_message;

   markov_data= hmt.state_tree_computation(error, trees, algorithm, characteristic_flag);
   if (markov_data == NULL)
   {
      error_message << error;
      PyErr_SetString(PyExc_RuntimeError, (error_message.str()).c_str());
      throw_error_already_set();
   }
   return markov_data;
}

void CiHmot_wrapper_plot_write(const Hidden_markov_out_tree& hmt,
                               const char* prefix,
                               const char* title)
{
   bool status= true;
   ostringstream error_message;
   Format_error error;

   status= hmt.plot_write(error, prefix, title);
   if (not status)
   {
      error_message << error;
      PyErr_SetString(PyExc_RuntimeError, (error_message.str()).c_str());
      throw_error_already_set();
   }
}

void CiHmot_wrapper_state_profile_plot_write(const Hidden_markov_out_tree& hmt,
                                             const char* prefix,
                                             const char* title,
                                             int identifier, int vertex,
                                             int entropy_algorithm)
{
   bool status= true;
   ostringstream error_message;
   Format_error error;

   status= hmt.state_profile_plot_write(error, prefix, identifier, vertex,
                                        title, entropy_algorithm);
   if (not status)
   {
      error_message << error;
      PyErr_SetString(PyExc_RuntimeError, (error_message.str()).c_str());
      throw_error_already_set();
   }
}

Hidden_markov_out_tree*
Hmt_wrapper_ascii_read1(const char * path)
{ return Hmt_wrapper_ascii_read(path); }


/*************************************************************
 *
 *  Wrappers for Python class CHmt_data:
 */


Hidden_markov_out_tree*
Chmt_data_wrapper_hidden_markov_out_tree_estimation_markov(const Hidden_markov_tree_data& hmtd,
                                                           const Hidden_markov_out_tree& ihmarkov,
                                                           bool counting_flag,
                                                           int state_trees,
                                                           int nb_iter,
                                                           bool characteristic_flag)
{
   Hidden_markov_out_tree *hmt= NULL;
   Format_error error;
   ostringstream error_message;

   hmt= hmtd.hidden_markov_out_tree_estimation(error, cout, ihmarkov,
                                               counting_flag, state_trees,
                                               nb_iter, characteristic_flag);
   if (hmt == NULL)
   {
      error_message << error;
      PyErr_SetString(PyExc_RuntimeError, (error_message.str()).c_str());
      throw_error_already_set();
   }
   return hmt;
}

Hidden_markov_out_tree*
Chmt_data_wrapper_hidden_markov_out_tree_estimation_self_transition(const Hidden_markov_tree_data& hmtd,
                                                                    int nb_state,
                                                                    bool left_right,
                                                                    bool counting_flag,
                                                                    int state_trees,
                                                                    double self_transition,
                                                                    int nb_iter,
                                                                    list force_param)
{
   bool status= true, several_errors= false;
   const int nb_variables= hmtd.get_nb_int() + hmtd.get_nb_float();
   int nb_fparam, p;
   char type= 'o'; // or may be not ...
   ostringstream error_message;
   object o;
   Format_error error;
   bool *fparam= NULL;
   Hidden_markov_out_tree *hmt= NULL;

   nb_fparam= boost::python::len(force_param);
   if (nb_fparam > 0)
   {
      if (nb_fparam != nb_variables)
      {
         status= false;
         error_message << "bad size of argument list: " << nb_fparam
                       << ": should be the number of variables ("
                       << nb_variables << ")";
         PyErr_SetString(PyExc_ValueError, (error_message.str()).c_str());
         throw_error_already_set();
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
            PyErr_SetString(PyExc_TypeError, (error_message.str()).c_str());
            throw_error_already_set();
         }
      }
   }
   if (status)
   {
      hmt= hmtd.hidden_markov_out_tree_estimation(error, cout, type, nb_state,
                                                  left_right, counting_flag,
                                                  state_trees, self_transition,
                                                  nb_iter, fparam);

      if (hmt == NULL)
      {
         error_message << error;
         PyErr_SetString(PyExc_RuntimeError, (error_message.str()).c_str());
         throw_error_already_set();
      }
   }
   return hmt;
}

str CHmt_data_wrapper_ascii_write0(const Hidden_markov_tree_data& hmtd)
{
   std::stringstream ss;
   str res;

   hmtd.line_write(ss);
   res= str(ss.str());

   return res;
}

void CHmt_data_wrapper_plot_write(const Hidden_markov_tree_data& reftree,
                                  const char* prefix,
                                  const char* title)
{
   bool status= true;
   ostringstream error_message;
   Format_error error;

   status= reftree.plot_write(error, prefix, title);
   if (not status)
   {
      error_message << error;
      PyErr_SetString(PyExc_RuntimeError, (error_message.str()).c_str());
      throw_error_already_set();
   }
}

Distribution_data* Chmt_data_wrapper_extract_value(const Hidden_markov_tree_data& reftree,
                                                   int type,
                                                   int variable,
                                                   int value)
{
   ostringstream error_message;
   // bool status= true;
   Format_error error;
   Distribution_data *histo= NULL;

   histo= reftree.extract(error, type, variable, value);
   if (histo == NULL)
   {
      error_message << error;
      PyErr_SetString(PyExc_RuntimeError, (error_message.str()).c_str());
      throw_error_already_set();
   }
   return histo;
}

}// namespace


// Module ======================================================================
BOOST_PYTHON_MODULE(chmt)
{
    class_< Hidden_markov_tree >
    ("CHmt", init< const Hidden_markov_tree&, optional< bool, bool> >())
        .def("ExtractData", &CHmt_wrapper_extract_data,
                            return_value_policy< manage_new_object >())
        .def("Likelihood", &CHmt_wrapper_likelihood,
                            "Likelihood(self, trees) -> float \n\n"
                            "Return the likelihood of the parameters"
                            "for the given trees")
        .def("NbInt", &Hidden_markov_tree::get_nb_ioutput_process)
        .def("NbFloat", &Hidden_markov_tree::get_nb_doutput_process)
        .def("NbStates", &Hidden_markov_tree::get_nb_state,
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

    class_< Hidden_markov_out_tree, bases<Hidden_markov_tree> >
    ("CiHmot", init< const Hidden_markov_out_tree&, optional< bool, bool> >())
        // .def("AsciiWrite", &CiHmot_wrapper_ascii_write0)
        .def("ExtractData", &CiHmot_wrapper_extract_data,
                            return_value_policy< manage_new_object >())
        .def("Display", &CiHmot_wrapper_display)
        .def("ComputeStateTrees", &CiHmot_wrapper_compute_state_trees,
                                  return_value_policy< manage_new_object >())
        .def("FileAsciiWrite", &CiHmot_wrapper_file_ascii_write1)
        .def("FileAsciiWrite", &CiHmot_wrapper_file_ascii_write2)
        .def("Plot", &CiHmot_wrapper_plot_write)
        .def("Simulate", &CiHmot_wrapper_simulate_histo,
                         return_value_policy< manage_new_object >())
        .def("Simulate", &CiHmot_wrapper_simulate_size,
                         return_value_policy< manage_new_object >())
        .def("Simulate", &CiHmot_wrapper_simulate_trees,
                         return_value_policy< manage_new_object >())
        .def("SpreadsheetWrite", &CiHmot_wrapper_spreadsheet_write1)
        .def("StateProfile", &CiHmot_wrapper_state_profile,
                             "StateProfile(self, int, int, int, int) -> list \n\n"
                             "return trees object and strings"
                             "for the state tree analysis\n")
        .def("StateProfilePlot", &CiHmot_wrapper_state_profile_plot_write,
                                 "StateProfilePlot(self, str, str, int, int) -> void \n\n"
                                 "write Gnuplot files for state and entropy"
                                 "profiles for the state tree analysis\n")
    ;

    class_< Hidden_markov_tree_data, bases<Trees> >
    ("CHmt_data", init< const Hidden_markov_tree_data&, optional< bool> >())
        .def(init< const Trees&>())
        // .def(init< const Trees&, int>())
        // .def("__init__", make_constructor(CHmt_data_wrapper_init1))
        .def("StateTrees",
             &Hidden_markov_tree_data::get_state_hidden_markov_tree_data,
             return_value_policy< manage_new_object >(),
             "StateTrees(self) -> CHmt_data \n\n"
             "Return a CHmt_data containing the states as a variable \n")
        .def("SmoothedProbaTrees",
             &Hidden_markov_tree_data::get_state_smoothed_hidden_markov_tree_data,
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
        .def("Plot", &CHmt_data_wrapper_plot_write)
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
