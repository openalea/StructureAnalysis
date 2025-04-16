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
   const MarkovianSequences **seq_index_markov= new const MarkovianSequences*[1];
   const MarkovianSequences **seq_index_markov2= new const MarkovianSequences*[1];
   Sequences *seq_merge= NULL, **seq_index= NULL;
   int *ilength= NULL;
   int ***index = NULL;
   stat_tool::variable_nature itype= INT_VALUE;
   StatError error;
   std::vector< int > select;
   const char * hsmcrefpath= "../../share/data/switching_lmm_irred.hsc";
   DiscreteParametric *d= NULL;
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

#     ifdef DEBUG
      // printing of the simulated sequences
      cout << "Simulated sequences: " << endl;
      hsmd->MarkovianSequences::ascii_data_write(cout);
      cout << "End simulated sequences" << endl;
#     endif

      select.resize(1);
      select[0] = 2;
      // discard state variable
      seq_estim = hsmd->select_variable(error, 1, select, true);
      assert(seq_estim->get_type(0) == REAL_VALUE);
      delete hsmd;
      hsmd = NULL;
      ilength= new int[nb_sequence];
      index = new int**[nb_sequence];
      // build random invalid index
      d = new DiscreteParametric(NEGATIVE_BINOMIAL, 1, I_DEFAULT, 10., 0.5);
      for (u = 0; u < nb_sequence; u++) {
    	  index[u] = new int*[1];
    	  ilength[u] = length;
		  index[u][0] = new int[length];
    	  for (t = 0; t < length; t++) {
    		  index[u][0][t] = d->simulation();
    		  assert(index[u][0][t] < d->alloc_nb_value);
    	  }
      }
      seq_index = new Sequences*[1];
      seq_index[0] = new Sequences(nb_sequence , NULL , ilength , IMPLICIT_TYPE ,
    		  	     	  	  	  	1 , itype , index);

      delete [] ilength;
      ilength = NULL;
      // add some index
      for (u = 0; u < nb_sequence; u++) {
		  delete [] index[u][0];
		  delete [] index[u];
      }
      delete [] index;
	  delete [] ilength;
	  index = NULL;
	  ilength = NULL;

      seq_index_markov[0] = new MarkovianSequences(*seq_index[0]);
      delete seq_index[0];
      seq_index[0] = NULL;
      seq_merge = seq_estim->merge_variable(error, 1, seq_index_markov);
      if (seq_merge == NULL) {
  		 cout << error;
  		 return 1;
  	  }
      assert(seq_merge->get_type(0) == REAL_VALUE);
      delete seq_index_markov[0];
      delete [] seq_index_markov;
      seq_index[0] = seq_merge->set_variable_as_index_parameter(error, 2, TIME);

      delete seq_index[0];
      seq_index[0] = NULL;

      // Rebuild index without errors
      error.init();
      ilength= new int[nb_sequence];
      index = new int**[nb_sequence];
      // build random invalid index
      d = new DiscreteParametric(NEGATIVE_BINOMIAL, 1, I_DEFAULT, 10., 0.5);
      for (u = 0; u < nb_sequence; u++) {
    	  index[u] = new int*[1];
    	  ilength[u] = length;
		  index[u][0] = new int[length];
		  index[u][0][0] = d->simulation();
    	  for (t = 1; t < length; t++) {
    		  index[u][0][t] = index[u][0][t-1] + d->simulation();
    	  }
      }
      seq_index = new Sequences*[1];
      seq_index[0] = new Sequences(nb_sequence , NULL , ilength , IMPLICIT_TYPE ,
    		  	     	  	  	  	1 , itype , index);

      delete [] ilength;
      ilength = NULL;
      // add some index
      for (u = 0; u < nb_sequence; u++) {
		  delete [] index[u][0];
		  delete [] index[u];
      }
      delete [] index;
	  delete [] ilength;
	  index = NULL;
	  ilength = NULL;

      seq_index_markov2[0] = new MarkovianSequences(*seq_index[0]);
      delete seq_index[0];
      seq_index[0] = NULL;
      seq_merge = seq_estim->merge_variable(error, 1, seq_index_markov2);
      if (seq_merge == NULL) {
  		 cout << error;
  		 return 1;
  	  }
      assert(seq_merge->get_type(0) == REAL_VALUE);
      delete seq_index_markov2[0];
      delete [] seq_index_markov2;
      delete seq_estim;
      seq_index[0] = seq_merge->set_variable_as_index_parameter(error, 2, TIME);

      delete d;
      d = NULL;

      if (seq_index[0] != NULL) {
    	  cout << "Sequence with index: " << endl;
    	  seq_index[0]->ascii_write(cout);
    	  delete seq_merge;
    	  seq_estim = new MarkovianSequences(*seq_index[0]);
    	  delete seq_index[0];
    	  delete [] seq_index;
      } else {
   		 cout << error;
   		 return 1;
      }

      // resimulate with new index
     hsmd = hsmc_ref->semi_markov_switching_lm_simulation(error, 1, *seq_estim, I_DEFAULT);
     delete seq_estim;
     seq_estim = NULL;
     if (hsmd == NULL) {
		 cout << error;
		 return 1;
     }
     assert(hsmd->get_type(1) == REAL_VALUE);
     // discard state variable
     select[0] = 2;
     seq_estim = hsmd->select_variable(error, 1, select, true);
     delete hsmd;
	 hsmd == NULL;
     if (seq_estim == NULL) {
		 cout << error;
		 return 1;
     }
     assert(seq_estim->get_type(0) == REAL_VALUE);

	 hsmc_est_file = seq_estim->hidden_semi_markov_estimation(error, &cout, *hsmc_ref, geometric_poisson , common_dispersion, estimator, counting_flag, state_sequence, 300);
     if (hsmc_est_file != NULL) {
    	 plotable = hsmc_est_file->get_plotable();
    	 delete plotable;
     } else {
		 cout << error;
		 return 1;
	 }

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
