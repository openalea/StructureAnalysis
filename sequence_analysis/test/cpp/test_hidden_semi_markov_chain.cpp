/****************************************************************
 *
 *  Test hidden semi-Markov tree
 */

#include "stat_tool/stat_tools.h"
#include "stat_tool/vectors.h"
#include "stat_tool/discrete_mixture.h"
#include "stat_tool/distribution.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "sequence_analysis/sequences.h"
#include "sequence_analysis/semi_markov.h"
#include "sequence_analysis/hidden_semi_markov.h"

using namespace stat_tool;
using namespace sequence_analysis;

int main(void)
{

   register int var, t, u, nb_integral, cptm= 0; // j, i, nb_states,

   bool status, geometric_poisson=false, common_dispersion=false, counting_flag=true, state_sequence=true;
   const int nb_sequence = 30 , length = 100;
   HiddenSemiMarkov *hsmc= NULL, *hsmc_ref= NULL, *hsmc_est_file= NULL;
   SemiMarkovData *hsmd= NULL;
   stat_tool::censoring_estimator estimator=stat_tool::COMPLETE_LIKELIHOOD;
   // Hidden_variable_order_markov *hmc= NULL, *hmc_init= NULL;
   MultiPlotSet *plotable=NULL;
   MarkovianSequences *seq_estim= NULL;
   StatError error;
   std::vector< int > select;
   const char * hsmcrefpath= "../../share/data/test_hidden_semi_markov_param.dat";
   // const char * hmcinitpath= "./hmc_init.hvom";

   // reading and printing of a hidden Markov out tree
   hsmc_ref = HiddenSemiMarkov::ascii_read(error, hsmcrefpath);
   cout << error;

   if (hsmc_ref != NULL)
   {
      cout << "Reference hidden semi Markov model : " << endl;
      hsmc_ref->ascii_write(cout, false);
      cout << endl;

      set_seed(1);

      plotable = hsmc_ref->get_plotable();
      delete plotable;
      // simulation of hidden Markov out trees
      hsmd= hsmc_ref->simulation(error, nb_sequence, length);
      cout << error;

      select.resize(1);
      select[0] = 2;
      // discard state variable
      seq_estim = hsmd->select_variable(error, 1, select, true);
    		  
#     ifdef DEBUG
      // printing of the simulated trees
      cout << "Simulated sequences: " << endl;
      hsmd->MarkovianSequences::ascii_data_write(cout);
      cout << "End simulated sequences" << endl;
#     endif

    if (seq_estim != NULL)
		 hsmc_est_file = seq_estim->hidden_semi_markov_estimation(error, &cout, *hsmc_ref, geometric_poisson , common_dispersion, estimator, counting_flag, state_sequence, 300);
	 else {
		 cout << error;
		 return 1;
	 }
     plotable = hsmc_est_file->get_plotable();
     delete plotable;
   }

   if (hsmc_ref != NULL)
		delete hsmc_ref;
	 hsmc_ref = NULL;
	 if (hsmd != NULL)
		delete hsmd;
	 if (seq_estim != NULL)
		 delete seq_estim;
	 seq_estim = NULL;
	 hsmd = NULL;
	 if (hsmc_est_file != NULL) {
		 cout << "Estimated model:" << endl;
		 hsmc_est_file->ascii_write(cout);
		 delete hsmc_est_file;
	 }
	 hsmc_est_file = NULL;


   return 0;
}
