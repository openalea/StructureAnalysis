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
 *       $Id: hidden_markov_out_tree.cpp 3193 2007-05-29 10:03:19Z dufourko $
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

#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"

#include "tree/tree_simple.h"
#include "tree/tree_traits.h"
#include "tree/basic_visitors.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/distribution.h"   // definition of DiscreteParametricModel class
#include "stat_tool/vectors.h"
#include "sequence_analysis/sequences.h"
#include "int_fl_containers.h"
#include "generic_typed_edge_tree.h"
#include "typed_edge_trees.h"
#include "hidden_markov_tree.h"
#include "hidden_markov_out_tree.h"
#include "tree_labels.h"
// #include <sstream>

extern int cumul_method(int nb_value, const double *cumul, double scale= 1.);

using namespace Stat_trees;

/*****************************************************************
 *
 *  Default constructor of class HiddenMarkovOutTree
 *
 **/

HiddenMarkovOutTree::HiddenMarkovOutTree()
 : HiddenMarkovTree()
{}

/*****************************************************************
 *
 *  Constructor of HiddenMarkovOutTree class
 *  using the type of markov tree ('o'rdinary or 'e'quilibrium),
 *  the number of states, the number of observed integral
 *  and floating processes and the number of observed values
 *  for each variable (ending with the floating variables)
 *
 **/

HiddenMarkovOutTree::HiddenMarkovOutTree(char itype, int inb_state,
                                         int inb_ioutputprocess,
                                         int inb_doutputprocess,
                                         int* nb_value,
                                         bool* force_param)
 : HiddenMarkovTree(itype, inb_state, 1, inb_ioutputprocess,
                    inb_doutputprocess, nb_value, force_param)
{}

/*****************************************************************
 *
 *  Constructor of HiddenMarkovOutTree class
 *  using a Chain object, the number of observed integral processes,
 *  the processes themselves, the tree size
 *  and a flag on the counting distribution computation
 *
 **/

HiddenMarkovOutTree::HiddenMarkovOutTree(const Chain * pchain,
                                         int inb_ioutput_process,
                                         NonparametricProcess** pobservation,
                                         int size, bool counting_flag)
 : HiddenMarkovTree(pchain, 1, inb_ioutput_process, pobservation, size, counting_flag)
{}

/*****************************************************************
 *
 *  Constructor of HiddenMarkovOutTree class
 *  using a Chain object, the number of observed integral and
 *  floating processes,the processes themselves, the tree size
 *  and a flag on the counting distribution computation
 *
 **/

HiddenMarkovOutTree::HiddenMarkovOutTree(const Chain * pchain,
                                         int inb_ioutput_process,
                                         int inb_doutput_process,
                                         NonparametricProcess** np_observation,
                                         DiscreteParametricProcess** ip_observation,
                                         DiscreteParametricProcess** dp_observation,
                                         int size, bool counting_flag)
 : HiddenMarkovTree(pchain, 1, inb_ioutput_process, inb_doutput_process,
                    np_observation, ip_observation, dp_observation,
                    size, counting_flag)
 {}

/*****************************************************************
 *
 *  Copy constructor of HiddenMarkovOutTree class
 *  using a flag on the HiddenMarkovTreeData copy
 *  and one on the characteristic distribution copy
 *
 **/

HiddenMarkovOutTree::HiddenMarkovOutTree(const HiddenMarkovOutTree& markov,
                                         bool data_flag,
                                         bool characteristic_flag)
 : HiddenMarkovTree(markov, data_flag, characteristic_flag)
{}

/*****************************************************************
 *
 *  Constructor of HiddenMarkovOutTree class
 *  converting a HiddenMarkovTree to a HiddenMarkovOutTree,
 *  assuming that _ch_order is equal to one
 *
 **/

HiddenMarkovOutTree::HiddenMarkovOutTree(const HiddenMarkovTree& markov,
                                         bool data_flag,
                                         bool characteristic_flag)
 : HiddenMarkovTree()
{
   if (_ch_order == 1)
   {
       Chain::copy(markov);
       HiddenMarkovTree::copy(markov, data_flag, characteristic_flag);
   }
}

/*****************************************************************
 *
 *  Copy constructor of HiddenMarkovOutTree class
 *  using a default self transition probability
 *
 **/

HiddenMarkovOutTree::HiddenMarkovOutTree(const HiddenMarkovOutTree& markov,
                                         double self_transition)
 : HiddenMarkovTree(markov, false)
{
   register int i;

   // initialization of the hidden Markov tree parameters
   init(self_transition);

   // initialization of the observation distributions
   for(i= 0; i < _nb_ioutput_process; i++)
   {
      if (npprocess[i+1] != NULL)
         npprocess[i+1]->init();
      else
         piprocess[i+1]->init();
   }

   for(i= 0; i < _nb_doutput_process; i++)
      if (pdprocess[i] != NULL)
         pdprocess[i]->init();
}

/*****************************************************************
 *
 *  Copy constructor for HiddenMarkovOutTree class
 *  with a structural transformation
 *  using the type of transformation ('s' : state, 'o' : "children order"),
 *  and a parameter (reference state/ "children order")
 *
 **/

HiddenMarkovOutTree::HiddenMarkovOutTree(const HiddenMarkovOutTree& markov,
                                         char manip,
                                         int param)
 : HiddenMarkovTree(markov, manip, param)
{}

/*****************************************************************
 *
 *  Destructor for HiddenMarkovOutTree class
 *
 **/

HiddenMarkovOutTree::~HiddenMarkovOutTree()
{}
// calls HiddenMarkovTree::~HiddenMarkovTree()

/*****************************************************************
 *
 *  Return the data part of a HiddenMarkovOutTree
 *  using a StatError object, keeping a reference on self
 *
 **/

HiddenMarkovTreeData* HiddenMarkovOutTree::extract_data(StatError& error) const
{
   bool status= true;
   HiddenMarkovTreeData *tree= NULL;

   error.init();

   if (markov_data == NULL)
   {
      status= false;
      error.update(STAT_error[STATR_NO_DATA]);
   }

   if (status)
      if ((_nb_ioutput_process != markov_data->_nb_integral)
          || (_nb_doutput_process != markov_data->_nb_float))
      {
         status= false;
         error.update(STAT_TREES_error[STATR_STATE_TREES]);
      }

   if (status)
   {
      tree= new HiddenMarkovTreeData(*markov_data);
      tree->markov= new HiddenMarkovOutTree(*this, false);
   }

   return tree;
}

/*****************************************************************
 *
 *  Prints a HiddenMarkovOutTree on a single line
 *  using an output stream
 *
 **/

std::ostream& HiddenMarkovOutTree::line_write(std::ostream& os) const
{
   os << nb_state << " " << STAT_word[STATW_STATES];

   return os;
}

/*****************************************************************
 *
 *  Prints a HiddenMarkovOutTree
 *  using an output stream and a flag on the level of detail
 *
 **/

ostream& HiddenMarkovOutTree::ascii_write(ostream& os, bool exhaustive) const
{ return ascii_write(os, markov_data, exhaustive, false); }


/*****************************************************************
 *
 *  Prints a HiddenMarkovOutTree into a file
 *  using a StatError object, the path
 *  and a flag on the level of detail
 *
 **/

bool HiddenMarkovOutTree::ascii_write(StatError& error, const char * path,
                                      bool exhaustive) const
{
   bool status;
   ofstream out_file(path);

   error.init();

   if (!out_file)
   {
      status= false;
      error.update(STAT_error[STATR_FILE_NAME]);
   }

   else
   {
      status= true;
      ascii_write(out_file, markov_data, exhaustive, true);
   }

   return status;
}

/*****************************************************************
 *
 *  Prints a HiddenMarkovOutTree in a spreadsheet fashion
 *  using a StatError object and the path
 *
 **/

bool HiddenMarkovOutTree::spreadsheet_write(StatError& error,
                                            const char * path) const
{
   bool status;
   ofstream out_file(path);


   error.init();

   if (!out_file)
   {
      status= false;
      error.update(STAT_error[STATR_FILE_NAME]);
   }

   else
   {
      status= true;
      spreadsheet_write(out_file, markov_data);
   }

   return status;
}

/*****************************************************************
 *
 *  Gnuplot output for HiddenMarkovOutTree class
 *  using a StatError object, a prefix for the files
 *  and the title for figures
 *
 **/

bool HiddenMarkovOutTree::plot_write(StatError& error,
                                     const char * prefix,
                                     const char * title) const
{
   bool status= plot_write(prefix, title, markov_data);

   error.init();

   if (!status)
     error.update(STAT_error[STATR_FILE_PREFIX]);

   return status;
}

/*****************************************************************
 *
 *  Prints a HiddenMarkovOutTree and the corresponding data structure
 *  using an output stream, a HiddenMarkovTreeData object,
 *  a flag on the level of detail, a flag on the file use
 *  and a Test object
 *
 **/

ostream& HiddenMarkovOutTree::ascii_write(ostream& os,
                                          const HiddenMarkovTreeData * otrees,
                                          bool exhaustive,
                                          bool file_flag,
                                          const Test* test) const
{
   // register int i;
   // int variable, cumul_size, nb_output_process= _nb_ioutput_process+_nb_doutput_process;
   // FrequencyDistribution **observation= NULL;
   // TreeCharacteristics *characteristics= NULL;

   switch (type)
   {
      case 'o' :
         os << STAT_TREES_word[STATW_HIDDEN_MARKOV_OUT_TREE] << endl;
         break;
      case 'e' :
         os << STAT_TREES_word[STATW_EQUILIBRIUM_HIDDEN_MARKOV_OUT_TREE] << endl;
         break;
   }

   // printing of the hidden Markov tree parameters
   HiddenMarkovTree::ascii_write(os, otrees, exhaustive, file_flag, test, false);

   return os;
}

/*****************************************************************
 *
 *  Prints a HiddenMarkovOutTree and the corresponding data structure
 *  in a spreadsheet fashion
 *  using an output stream, a HiddenMarkovTreeData object
 *  and a Test object
 *
 **/

ostream& HiddenMarkovOutTree::spreadsheet_write(ostream& os, const HiddenMarkovTreeData * otrees,
                                                const Test * test) const
{
   register int i;
   int variable= 0, cumul_size, nb_output_process= _nb_ioutput_process+_nb_doutput_process;
   FrequencyDistribution **observation= NULL;
   TreeCharacteristics *characteristics= NULL;

   switch (type)
   {
      case 'o' :
         os << STAT_TREES_word[STATW_HIDDEN_MARKOV_OUT_TREE] << endl;
         break;
      case 'e' :
         os << STAT_TREES_word[STATW_EQUILIBRIUM_HIDDEN_MARKOV_OUT_TREE] << endl;
         break;
   }

   // printing of the Markov tree parameters

   spreadsheet_print(os);

   if ((otrees != NULL) && (otrees->_type[0] == STATE))
   {
      variable= 0;
      characteristics= otrees->characteristics[variable];
   }

   npprocess[0]->spreadsheet_print(os, 0, 0, characteristics);

   // printing of the (characteristic ?) distributions of each
   // observed process

   if (nb_output_process > 0)
   {
      os << "\n" << nb_output_process << "\t"
         << STAT_word[nb_output_process == 1 ? STATW_OUTPUT_PROCESS : STATW_OUTPUT_PROCESSES] << endl;

      for(i= 1; i <= _nb_ioutput_process; i++)
      {
         os << "\n" << STAT_word[STATW_OUTPUT_PROCESS];
         os << "\t" << i;

         if (npprocess[i] != NULL)
            os << "\t" << STAT_word[STATW_NONPARAMETRIC];
         else
            os << "\t" << STAT_word[STATW_PARAMETRIC];

         os << endl;

         if (otrees != NULL)
         {
            switch (otrees->_type[0])
            {
               case INT_VALUE :
                  variable= i - 1;
                  break;
               case STATE :
                  variable= i;
                  break;
            }

            if (otrees->observation != NULL)
               observation= otrees->observation[variable];
            if (otrees->characteristics[variable] != NULL)
               characteristics= otrees->characteristics[variable];
         }
         if (npprocess[i] != NULL)
            npprocess[i]->spreadsheet_print(os, i, observation, characteristics);
         else
            piprocess[i]->spreadsheet_print(os, observation);
      }

      for(i= 0; i < _nb_doutput_process; i++)
      {
         os << "\n" << STAT_word[STATW_OUTPUT_PROCESS];
         os << "\t" << i;

         os << "\t" << STAT_word[STATW_PARAMETRIC] << endl;

         if (otrees != NULL)
         {
            if (otrees->observation != NULL)
               observation= otrees->observation[i];
            // warning : akward !
         }
         pdprocess[i]->spreadsheet_print(os, observation);
      }
   }

   if (otrees != NULL)
   {
      int nb_parameter= nb_parameter_computation(MIN_PROBABILITY);
      double information, likelihood;


      // printing of the quantities for which the characteristic distributions
      // are invariant - if any

      os << "\n" << STAT_label[STATL_TREE_SIZE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
      otrees->hsize->spreadsheet_characteristic_print(os);

      os << "\n\t" << STAT_label[STATL_TREE_SIZE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
      otrees->hsize->spreadsheet_print(os);

      cumul_size= otrees->cumul_size_computation();
      os << "\n" << STAT_label[STATL_CUMULATIVE_SIZE] << "\t" << cumul_size << endl;

      // printing of the tree information quantity in the iid case

      information= otrees->iid_information_computation();

      os << "\n" << STAT_TREES_label[STATL_TREES_IID_INFORMATION] << "\t" << information << "\t"
         << information / cumul_size << endl;

      // printing of the likelihood

      if (otrees->likelihood != D_INF)
      {
         os << "\n" << STAT_label[STATL_STATE_TREES_LIKELIHOOD] << "\t" << otrees->likelihood << "\t"
            << STAT_label[STATL_NORMALIZED] << "\t" << otrees->likelihood / cumul_size << endl;
      }

      likelihood= otrees->hidden_likelihood;

      if (likelihood != D_INF)
      {
         os << "\n" << STAT_label[STATL_OBSERVED_TREES_LIKELIHOOD] << "\t" << likelihood << "\t"
            << STAT_label[STATL_NORMALIZED] << "\t" << likelihood / cumul_size << endl;
      }

      if ((likelihood != D_INF) && (nb_component == 1))
      {
         if (nb_parameter < cumul_size-1)
         {
            os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
               << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AICc] << ")\t"
               << 2 * (likelihood - (double)(nb_parameter * cumul_size) /
                  (double)(cumul_size-nb_parameter-1)) << endl;
         }

         os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
            << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BIC] << ")\t"
            << 2 * likelihood - nb_parameter * log((double)cumul_size) << endl;
      }
   }

   if (test)
   {
      os << "\n";
      test->spreadsheet_print(os);
   }

   return os;
}

/*****************************************************************
 *
 *  Gnuplot output for HiddenMarkovOutTree class
 *  using a prefix for the files, the title of figures
 *  and the observed trees
 *
 **/

bool HiddenMarkovOutTree::plot_write(const char * prefix, const char * title,
                                     const HiddenMarkovTreeData * otrees) const
{
   bool status;
   register int i;
   int variable= 0; //, cumul_size, nb_output_process= _nb_ioutput_process+_nb_doutput_process;
   FrequencyDistribution *hsize= NULL, **observation= NULL;
   TreeCharacteristics *characteristics= NULL;

   // print characteristic distributions for state process
   if ((otrees != NULL) && (otrees->_type[0] == STATE))
   {
      variable= 0;
      characteristics= otrees->characteristics[variable];
      hsize= otrees->hsize;
   }

   status= npprocess[0]->plot_print(prefix, title, 0, NULL, characteristics, hsize);

   // print characteristic distributions of each observed process

   if (status)
   {
      if (otrees != NULL)
         hsize= otrees->hsize;

      for(i= 1; i <= _nb_ioutput_process; i++)
      {
         if (otrees != NULL)
         {
            switch (otrees->_type[0])
            {
               case INT_VALUE :
                  variable= i - 1;
                  break;
               case STATE :
                  variable= i;
                  break;
            }

            if (otrees->observation != NULL)
               observation= otrees->observation[variable];

            if (otrees->characteristics[variable] != NULL)
               characteristics= otrees->characteristics[variable];
         }

         if (npprocess[i] != NULL)
            npprocess[i]->plot_print(prefix, title, i, observation,
                                     characteristics, hsize);
         else
            piprocess[i]->plot_print(prefix, title, i, observation);
      }

      for(i= 0; i < _nb_doutput_process; i++)
      {
         if (otrees->observation != NULL)
            observation= otrees->observation[i];

         piprocess[i]->plot_print(prefix, title, i, observation);
      }
   }
   return status;
}

void HiddenMarkovOutTree::get_state_marginal_distribution (const HiddenMarkovTreeData& trees,
                                                           double_array_3d& res) const
{ state_marginal_distribution(trees, res); }

void HiddenMarkovOutTree::get_output_conditional_distribution(const HiddenMarkovTreeData& trees,
                                                              double_array_3d& res) const
{ output_conditional_distribution(trees, res); }

double HiddenMarkovOutTree::get_upward_step(const HiddenMarkovTreeData& trees,
                                            double_array_3d& upward_prob,
                                            double_array_3d& upward_parent_prob,
                                            double_array_3d& state_entropy,
                                            double_array_3d marginal_prob,
                                            double_array_3d output_cond_prob,
                                            double& entropy1) const
{ return upward_step(trees, upward_prob, upward_parent_prob, state_entropy,
                     marginal_prob, output_cond_prob, entropy1); }

void HiddenMarkovOutTree::get_downward_step(const HiddenMarkovTreeData& trees,
                                            double_array_3d& downward_prob,
                                            double_array_3d*& downward_pair_prob,
                                            double_array_3d upward_prob,
                                            double_array_3d upward_parent_prob,
                                            double_array_3d marginal_prob,
                                            double_array_3d output_cond_prob,
                                            double& entropy2) const
{ return downward_step(trees, downward_prob, downward_pair_prob, upward_prob, upward_parent_prob,
                       marginal_prob, output_cond_prob, entropy2);}

double HiddenMarkovOutTree::get_viterbi(const HiddenMarkovTreeData& trees,
                                        int index)
{
   double res;

   create_cumul();
   log_computation();
   res= viterbi(trees, index);
   remove_cumul();
   return res;
}

double HiddenMarkovOutTree::get_upward_downward(const HiddenMarkovTreeData& trees,
                                                int index,
                                                ostream* os,
                                                char format) const
{
   double max_marginal_entropy, entropy1, likelihood, state_likelihood;
   std::deque<int> *path= NULL;
   state_likelihood= upward_downward(trees, max_marginal_entropy, entropy1,
                                     likelihood, path, index, os, format, 0);
   if (path != NULL)
   {
      delete path;
      path= NULL;
   }
   return state_likelihood;
}

/*****************************************************************
 *
 *  Computation of the parameter number of a HiddenMarkovOutTree
 *
 **/

int HiddenMarkovOutTree::nb_parameter_computation(double min_probability) const
{
   int nb_parameter= Chain::nb_parameter_computation(min_probability);
   register int var; //, j, val


   if (_nb_ioutput_process > 0)
      for(var= 0; var < _nb_ioutput_process; var++)
      {
         if (npprocess[var+1] != NULL)
            nb_parameter+= npprocess[var+1]->nb_parameter_computation(min_probability);
         else
            nb_parameter+= piprocess[var+1]->nb_parameter_computation();
      }
   // if (_nb_ioutput_process > 0)
   // ... idem
   return nb_parameter;
}

/*****************************************************************
 *
 *  Likelihood correction for a HiddenMarkovOutTree,
 *  using a given set of trees (HiddenMarkovTreeData)
 *
 **/

double HiddenMarkovOutTree::likelihood_correction(const HiddenMarkovTreeData& otrees) const
{
   register int j, var, t, s, val;
   int *cinitial;
   double correction, *pinitial;


   correction= 0.;

   cinitial= otrees.chain_data->initial;
   pinitial= initial;
   for(j= 0; j < nb_state; j++)
   {
     if (*cinitial > 0)
        correction+= *cinitial * log(*pinitial);
     cinitial++;
     pinitial++;
   }

   if (_nb_ioutput_process > 0)
      for(t= 0; t < otrees._nb_trees; t++)
         for(var= 0; var < _nb_ioutput_process; var++)
            if (npprocess[var+1]!=NULL)
            {
               s= (otrees.state_trees[t])->get((otrees.state_trees[t])->root()).Int();
               val= ((otrees.trees[t])->get(otrees.trees[t]->root())).Int(j);
               correction+= log(npprocess[var+1]->observation[s]->mass[val]);
            }

   return correction;
}


HiddenMarkovOutTree* Stat_trees::hidden_markov_out_tree_ascii_read(StatError& error,
                                                                  const char * path,
                                                                  int size,
                                                                  bool counting_flag,
                                                                  double cumul_threshold)
{
   RWLocaleSnapshot locale("en");
   RWCString buffer, token;
   size_t position;
   char type= 'v';
   bool status;
   register int i;
   int line;
   HiddenMarkovTree *markov= NULL;
   HiddenMarkovOutTree *res= NULL;
   ifstream in_file(path);

   error.init();

   if (!in_file)
      error.update(STAT_error[STATR_FILE_NAME]);
   else
   {
      status= true;
      line= 0;

      if (size < 2)
      {
         status= false;
         error.update(STAT_error[STATR_SMALL_TREE_SIZE]);
      }

      if (size > MAX_SIZE)
      {
         status= false;
         error.update(STAT_error[STATR_BIG_TREE_SIZE]);
      }

      while (buffer.readLine(in_file, false))
      {
         line++;

#        ifdef DEBUG
         cout << line << "  " << buffer << endl;
#        endif

         position= buffer.first('#');
         if (position != RW_NPOS)
            buffer.remove(position);

         i= 0;

         RWCTokenizer next(buffer);

         while (!((token = next()).isNull()))
         {
           // test keyword (EQUILIBRIUM) HIIDEN_MARKOV_OUT_TREE

            if (i == 0)
            {
               if (token == STAT_TREES_word[STATW_HIDDEN_MARKOV_OUT_TREE])
                  type= 'o';
               else
                  if (token == STAT_TREES_word[STATW_EQUILIBRIUM_HIDDEN_MARKOV_OUT_TREE])
                     type= 'e';
                  else
                  {
                     status= false;
                     ostringstream correction_message;
                     correction_message << STAT_TREES_word[STATW_HIDDEN_MARKOV_OUT_TREE] << " or "
                                        << STAT_TREES_word[STATW_EQUILIBRIUM_HIDDEN_MARKOV_OUT_TREE];
                     error.correction_update(STAT_parsing[STATP_KEY_WORD], (correction_message.str()).c_str(), line);
                  }
            }
           i++;
         }

         if (i > 0)
         {
            if (i != 1)
            {
               status= false;
               error.update(STAT_parsing[STATP_FORMAT], line);
            }
            break;
         }
      }

      if (status)
      {
         // analysis of the hidden Markov tree format and processing

         markov= Stat_trees::hidden_markov_tree_ascii_read(error, path, size,
                                                           counting_flag, cumul_threshold);
         if ((markov == NULL) || (error.get_nb_error() > 0))
            status= false;
         else
            if (markov->get_ch_order() != 1)
            {
               status= false;
               error.update(STAT_parsing[STATP_ORDER], line);
            }
      }

      if (status)
         res= new HiddenMarkovOutTree(*markov, false, false);

      if (markov != NULL)
         delete markov;

   }
   return res;
}

/*****************************************************************
 *
 *  Left (bit) shift operator of HiddenMarkovOutTree
 *
 **/

ostream& Stat_trees::operator<<(ostream &os, const HiddenMarkovOutTree& hmarkov)
{ return hmarkov.ascii_write(os); }
