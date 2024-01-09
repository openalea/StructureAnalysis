/****************************************************************
 *
 *  Comparison of a hidden Markov chain and a hidden Markov line-tree
 */

#include "stat_tool/stat_tools.h"
#include "stat_tool/vectors.h"
#include "stat_tool/discrete_mixture.h"
#include "stat_tool/distribution.h"
#include "stat_tool/markovian.h"
#include "sequence_analysis/sequences.h"
#include "sequence_analysis/semi_markov.h"
#include "sequence_analysis/hidden_semi_markov.h"

using namespace stat_tool;
using namespace sequence_analysis;

int main(void)
{
   const int nb_sequence=10, length=300;
   bool status;
   const int tid= 0; // id of the sequence to be analyzed
   HiddenSemiMarkov *hmc= NULL, *hmc_init= NULL;
   SemiMarkovData *seq= NULL;
   StatError error;
   const char * hmcinitpath= "../../../test/cpp/hmc_init.hmc";

   hmc_init= hidden_semi_markov_ascii_read(error, hmcinitpath);
   cout << error;
   if (hmc_init != NULL)
   {

      cout << endl << "Initialization HMC : " << endl;
      hmc_init->ascii_write(cout, false);
      cout << endl;
      seq= hmc_init->simulation(error, nb_sequence, length, true);
      cout << endl << error;

      if (seq != NULL)
      {
         hmc= seq->hidden_semi_markov_estimation(error, cout, *hmc_init, true, COMPLETE_LIKELIHOOD,
                                                 true, VITERBI, 100);
         cout << error;
         if (hmc != NULL)
         {

            // check the likelihood computation
            cout << "Estimated HMC : " << endl;
            hmc->ascii_write(cout, false);
            cout << endl;

            // check the state profile computation
            status= hmc->state_profile_ascii_write(error, cout, tid, SSTATE, GENERALIZED_VITERBI,
                                                   5);
            if (!status)
               cout << error << endl;

            delete hmc;
            hmc= NULL;
         }
         delete seq;
         seq= NULL;
      }

      delete hmc_init;
      hmc_init= NULL;

   }
   return 0;
}
