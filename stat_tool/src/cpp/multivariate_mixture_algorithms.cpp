/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2015 CIRAD/INRA/Inria Virtual Plants
 *
 *       File author(s): J.-B. Durand (jean-baptiste.durand@imag.fr) and
 *                       Y. Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id: multivariate_mixture_algorithms.cpp 3256 2007-06-06 12:32:54Z dufourko $
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
 *       MERCHANTABILITY or FITNESS For A PARTICULAR PURPOSE. See the
 *       GNU General Public License for more details.
 *
 *       You should have received a copy of the GNU General Public
 *       License along with this program; see the file COPYING. If not,
 *       write to the Free Software Foundation, Inc., 59
 *       Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 *  ----------------------------------------------------------------------------
 */



#include <math.h>
#include <sstream>
#include <iomanip>
#include <assert.h>

#include <boost/scoped_array.hpp>

#include "stat_tools.h"
#include "distribution.h"
#include "discrete_mixture.h"
#include "markovian.h"
#include "vectors.h"
#include "multivariate_mixture.h"
#include "stat_label.h"

#include "stat_tool/distribution_reestimation.hpp"

using namespace std;


namespace stat_tool {



/*--------------------------------------------------------------*
 *
 *  Calcul de la quantite d'information d'un objet MultivariateMixtureData.
 *
 *--------------------------------------------------------------*/

double MultivariateMixtureData::information_computation() const

{
  register int i, var;
  double information , buff;

  information = weight->information_computation();

  if (information != D_INF) {
      for (i = 0;i < nb_component;i++) {
    if (weight->frequency[i] > 0) {
      buff = 0.;
      for (var = 0;var < nb_variable;var++)
        if (type[var] != STATE) {
          buff += component[var][i]->information_computation();
          if (buff != D_INF) {
        information += buff;
          }
          else {
        information = D_INF;
        break;
          }
      }
    }
      }
  }

  return information;
}

/*--------------------------------------------------------------*
 *
 *  Calcul de la densite conditionnelle d'un objet Vectors pour un melange de lois.
 *
 *  argument : reference sur un objet Vectors, sur les densites conditionnelles,
 *  et flag sur le calcul de la log densite.
 *
 *--------------------------------------------------------------*/

void MultivariateMixture::get_output_conditional_distribution(const Vectors &mixt_data,
                                                              double** &output_cond,
                                                              bool log_computation) const
{
  int n, i, var;
  const int nb_vector = mixt_data.get_nb_vector(),
    nb_variable = mixt_data.get_nb_variable();

  if (output_cond == NULL) {
    output_cond = new double*[nb_vector];
    for (n = 0; n < nb_vector; n++)
      output_cond[n] = NULL;
  }

  for (n = 0; n < nb_vector; n++)
    if (output_cond[n] == NULL)
      output_cond[n] = new double[nb_component];

  for (n = 0; n < nb_vector; n++) {
    for (i = 0; i < nb_component; i++) {
      if (log_computation)
	output_cond[n][i] = 0.;
      else
	output_cond[n][i] = 1.;
      for (var = 0; var < nb_variable; var++) {
	if (pcomponent[var] != NULL) {
	  if (log_computation) {
	    if (pcomponent[var]->observation[i]->mass[mixt_data.int_vector[n][var]] > 0)
	      output_cond[n][i] += log(pcomponent[var]->observation[i]->mass[mixt_data.int_vector[n][var]]);
	    else {
	      output_cond[n][i] = D_INF;
	      break;
	    }
	  }
	  else
	    output_cond[n][i]*= pcomponent[var]->observation[i]->mass[mixt_data.int_vector[n][var]];
	}
	else {
	  if (log_computation) {
	    if (npcomponent[var]->observation[i]->mass[mixt_data.int_vector[n][var]] > 0)
	      output_cond[n][i] += log(npcomponent[var]->observation[i]->mass[mixt_data.int_vector[n][var]]);
	    else {
	      output_cond[n][i] = D_INF;
	      break;
	    }
	  }
	  else
	    output_cond[n][i]*= npcomponent[var]->observation[i]->mass[mixt_data.int_vector[n][var]];
	}
      }
    }
  }
}

/*--------------------------------------------------------------*
 *
 *  Calcul de la loi des etats sachant un objet Vectors, pour un melange de lois.
 *
 *  argument : reference sur un objet Vectors, lois des observations sachant
 *  les etats et reference sur lois des etats
 *
 *--------------------------------------------------------------*/

void MultivariateMixture::get_posterior_distribution(const Vectors &mixt_data,
                                                     double** output_cond,
                                                     double** &posterior_dist) const {
  int n, i;
  const int nb_vector = mixt_data.get_nb_vector();
  double marginal;

  if (posterior_dist == NULL) {
    posterior_dist = new double*[nb_vector];
    for (n = 0; n < nb_vector; n++)
      posterior_dist[n] = NULL;
  }

  for (n = 0; n < nb_vector; n++)
    if (posterior_dist[n] == NULL)
      posterior_dist[n] = new double[nb_component];

  for (n = 0; n < nb_vector; n++) {
    marginal = 0.;
    for (i = 0; i < nb_component; i++) {
      posterior_dist[n][i] = output_cond[n][i] * weight->mass[i];
      marginal += posterior_dist[n][i];
    }
    for (i = 0; i < nb_component; i++) {
      if (marginal > 0)
    posterior_dist[n][i] /= marginal;
      else
    posterior_dist[n][i] = 0.;
    }
  }
}

/*--------------------------------------------------------------*
 *
 *  Restauration des etats d'un objet Vectors, pour un melange de lois.
 *
 *  arguments : reference sur un objet StatError, sur un object Vectors,
 *  l' algorithme de restauration, l'index (a partir de 0)
 *  et la loi a posteriori des etats (mise a jour au besoin)
 *
 *--------------------------------------------------------------*/

std::vector<int>* MultivariateMixture::state_computation(StatError &error, const Vectors &vec,
                                                         double** &posterior_dist,
                                                         int algorithm, int index) const {

  bool status = true, delete_cond = false;
  int n, i, s;
  double pmax;
  // output_cond n'est pas utile (pas alloue) si posterior_dist est donne en argument
  // Dans le cas contraire, posterior_dist est mis a jour
  double **output_cond = NULL, **cp_posterior_dist = NULL;
  int nb_vector = vec.get_nb_vector(),
    nb_variable = vec.get_nb_variable();
  std::vector<int> *res_vec = NULL;

  error.init();

  if (index != I_DEFAULT) {
    if ((index < 0) || (index >= nb_vector)) {
      status = false;
      error.update(STAT_error[STATR_SAMPLE_INDEX]); }
    nb_vector = 1;
  }
  if (nb_variable != nb_var) {
    status = false;
    error.update(STAT_error[STATR_NB_VARIABLE]);
  }
  if (status) {
    if (posterior_dist == NULL) {
      delete_cond = true;
      get_output_conditional_distribution(vec, output_cond, false);
      get_posterior_distribution(vec, output_cond, cp_posterior_dist);
      posterior_dist = cp_posterior_dist;
    }
    res_vec = new std::vector<int>(nb_vector);
    for (n = 0; n < nb_vector; n++) {
      if ((index == I_DEFAULT) || (index == n)) {
    s = 0;
    pmax =  posterior_dist[n][0];
    for (i = 1; i < nb_component; i++)
      if (posterior_dist[n][i] > pmax) {
        pmax = posterior_dist[n][i];
        s = i;
      }
    if (index == I_DEFAULT)
      (*res_vec)[n] = s;
    else
      (*res_vec)[0] = s;
      }
    }
    if (delete_cond) {
      for (n = 0; n < nb_vector; n++) {
	delete [] output_cond[n];
	output_cond[n] = NULL;
      }
      delete [] output_cond;
      output_cond = NULL;
    }
  }
  return res_vec;
}

/*--------------------------------------------------------------*
 *
 *  Calcul de la vraisemblance d'un objet MultivariateMixtureData
 *  pour un melange de lois.
 *
 *  argument : reference sur un objet MultivariateMixtureData.
 *
 *--------------------------------------------------------------*/

double MultivariateMixture::likelihood_computation(const Vectors &mixt_data,
                                                   bool log_computation) const

{
  register int n, i, var, nb_vector = mixt_data.get_nb_vector();
  double log_likelihood = 1. , buff , **output_cond = NULL;
  double *pmass = NULL;

  get_output_conditional_distribution(mixt_data, output_cond, false);

  for (n = 0; n < nb_vector; n++) {
    buff = 0.;
    for (i = 0;i < nb_component;i++)
      buff += weight->mass[i] * output_cond[n][i];

    if (buff > 0)
      log_likelihood += log(buff);
    else {
      if (log_computation)
    return D_INF;
      else
    return 0.;
    }
  }

  if (!log_computation) {
    if (log_likelihood > D_INF)
      log_likelihood = exp(log_likelihood);
    else
      log_likelihood = 0.;
  }

  for (n = 0; n < nb_vector; n++) {
    delete [] output_cond[n];
    output_cond[n] = NULL;
    }

  delete [] output_cond;
  output_cond = NULL;

  return log_likelihood;
}

/*--------------------------------------------------------------*
 *
 *  Initialisation d'un melange de lois
 *
 *--------------------------------------------------------------*/

void MultivariateMixture::init() {

  unsigned int j, var;

  if (weight == NULL)
    weight = new DiscreteParametric(nb_component);

  if (pcomponent == NULL)
    pcomponent = new DiscreteParametricProcess*[nb_var];

  if (npcomponent == NULL)
    npcomponent = new CategoricalProcess*[nb_var];


  for(j = 0; j < nb_component; j++)
    weight->mass[j]= 1. / nb_component;

  for(var = 0; var < nb_var; var++) {
    if (pcomponent[var] != NULL)
      pcomponent[var]->init();
    else
      npcomponent[var]->init();
  }
}

/*--------------------------------------------------------------*
 *
 *  Estimation du nombre de composantes et des parametres
 *  d'un melange de lois par l'algorithme EM.
 *
 *  arguments : reference sur un objet StatError, modele initial,
 *  nombre d'iteration maximal, flag sur l'utilisation forcee ou
 *  non de lois d'observation de meme nature que le modele initial
 *  (parametriques ou non)
 *
 *--------------------------------------------------------------*/

MultivariateMixture* Vectors::mixture_estimation(StatError &error, ostream& os,
                                                 const MultivariateMixture &imixture,
                                                 int nb_iter, bool *force_param) const {

  bool status, state_simulation, all_states_used;
  register int i , j , var;
  unsigned int n;
  int max_nb_value, iter, nb_likelihood_decrease, val;
  double likelihood= D_INF, previous_likelihood, observation_likelihood ,
    min_likelihood= 0,  *reestim= NULL, saem_coef= 0., state_likelihood= D_INF,
    best_likelihood= D_INF; // best likelihood for SEM
  double *state_array; // simulated states for SEM
  double **output_cond = NULL, **cond_prob = NULL;
  StatError error_v;
  Reestimation<double> ***observation_reestim = NULL;
  Reestimation<double> *weight_reestim = NULL;
  FrequencyDistribution *hobservation= NULL;
  MultivariateMixture *mixt = NULL,  *best_mixt= NULL; // best model for SEM
  MultivariateMixtureData *mixt_data = NULL,
    *state_restoration= NULL;

# ifdef DEBUG
  double test[NB_STATE][4];
# endif

  error.init();

  /*   if (algorithm == VITERBI)
       error_v.init();*/

  // test of the number of observed values per variable

  status= false;
  for(var = 0; var < nb_variable; var++)
    if ((int)get_max_value(var)-(int)get_min_value(var) > 0) {
      status= true;
      break;
    }

  if (!status)
    error.update(STAT_error[STATR_NB_VALUE]);
  else {
    if (imixture.nb_var != nb_variable) {
      // number of variables for observations and model
      // must match
      status = false;
      error.update(STAT_error[STATR_NB_VARIABLE]);
    }
    else {
      for(var = 0; var < nb_variable; var++) {
    if (((imixture.npcomponent[var] != NULL) &&
         (imixture.npcomponent[var]->nb_value != (int)get_max_value(var)+1)) ||
        ((imixture.pcomponent[var] != NULL) &&
         (imixture.pcomponent[var]->nb_value < (int)get_max_value(var)+1))) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_OUTPUT_PROCESS] << " " << var << ": "
            << STAT_error[STATR_NB_VALUE];
      error.update((error_message.str()).c_str());
    }
    else
      if ((imixture.npcomponent[var] != NULL) && (marginal_distribution != NULL))
        for(j = 0; j < (int)get_max_value(var); j++)
          if (marginal_distribution[var]->frequency[j] == 0) {
        status= false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << var << ": "
                  << STAT_error[STATR_MISSING_VALUE] << " " << j;
        error.update((error_message.str()).c_str());
          }
      }
    }
  }

  if ((nb_iter != I_DEFAULT) && (nb_iter < 1)) {
    status= false;
    error.update(STAT_error[STATR_NB_ITERATION]);
  }

  /*   if (!((algorithm == VITERBI) || (algorithm == FORWARD_BACKWARD_SAMPLING) ||
       (algorithm == GIBBS_SAMPLING) || (algorithm == FORWARD_BACKWARD)))
       {
       status= false;
       error.update(STAT_TREES_error[STATR_EM_ALGORITHM]);
       }

       if ((algorithm == VITERBI) || (algorithm == FORWARD_BACKWARD_SAMPLING) ||
       (algorithm == GIBBS_SAMPLING))
       {
       if ((saem_exponent < 0) || (saem_exponent >= 1))
       {
       status= false;
       ostringstream error_message;
       error_message << STAT_TREES_error[STATR_SAEM_EXP] << ": " << saem_exponent;
       error.update((error_message.str()).c_str());
       }
       state_simulation= true;
       }
       else */
  state_simulation= false;

  if (status) {

    // create mixture

    mixt = new MultivariateMixture(imixture, false);

    if (state_simulation) {
      best_mixt = new MultivariateMixture(*mixt, false);
      // initialize hidden states for state simulation

      /*         if (algorithm == GIBBS_SAMPLING)
         mixt_data= mixt->state_tree_computation(error_v, *this, VITERBI, false); */

      if (mixt_data == NULL)
    mixt_data = new MultivariateMixtureData(*this, false);

      /* if (mixt_data->state_trees == NULL)
     mixt_data->build_state_trees(); */
    }

#   ifdef DEBUG
    cout << *mixt;
#   endif

    // create data structures for estimation

    observation_reestim = new Reestimation<double>**[mixt->nb_var];
    weight_reestim = new Reestimation<double>(mixt->nb_component);
    for(var = 0; var < mixt->nb_var; var++) {
      observation_reestim[var] = new Reestimation<double>*[mixt->nb_component];
      for(j = 0; j < mixt->nb_component; j++)
    observation_reestim[var][j] = new Reestimation<double>((int)get_max_value(var)+1);
    }

    max_nb_value = 0;
    for(var = 0; var < mixt->nb_var; var++)
      if ((mixt->pcomponent[var] != NULL) && (max_nb_value < (int)get_max_value(var)+1))
    max_nb_value= (int)get_max_value(var) + 1;

    if (max_nb_value > 0)
      hobservation= new FrequencyDistribution(max_nb_value);
    else
      hobservation= NULL;

    iter = 0;
    nb_likelihood_decrease = 0;

    if (state_simulation)
      state_array = new double[this->nb_vector];

    do {
      iter++;
      previous_likelihood = likelihood;
      likelihood = 0.;

      // initialization of the reestimation quantities

      for(var = 0; var < mixt->nb_var; var++)
    for(j = 0; j < mixt->nb_component; j++) {
      reestim = observation_reestim[var][j]->frequency;
      for(val = 0; val < (int)get_max_value(var) + 1; val++)
        reestim[val] = 0.; // *reestim++ = 0.;
    }

#     ifdef DEBUG
      for(i = 0; i < mixt->nb_component; i++)
    for(j = 0; j < 4; j++)
      test[i][j] = 0.;
#     endif

      // if ((algorithm == FORWARD_BACKWARD) || (saem_exponent != .0))
      // {
      mixt->get_output_conditional_distribution(*this, output_cond);

      mixt->get_posterior_distribution(*this,output_cond, cond_prob);

      likelihood = mixt->likelihood_computation(*this, true);
      if (likelihood == D_INF)
    break;
      // }
      /* else
     likelihood= mixt->likelihood_computation(*this);*/

      /*         if (algorithm == VITERBI)
         {
         state_restoration= mixt->state_tree_computation(error_v, *mixt_data, VITERBI, false);
         if (error_v.get_nb_error() > 0)
         break;
         }

         if (algorithm == FORWARD_BACKWARD_SAMPLING)
         {
         state_restoration= mixt->sstate_simulation(*mixt_data,
         state_likelihood,
         false);
         if (state_likelihood <= D_INF)
         break;
         }

         if (algorithm == GIBBS_SAMPLING)
         {
         state_restoration= mixt->gibbs_state_simulation(*mixt_data,
         state_likelihood,
         false);
         if (state_likelihood <= D_INF)
         break;
         } */

      /* if (state_restoration != NULL) {
      // transition counts
      for(t= 0; t < state_restoration->_nb_trees; t++) {
      for(j= 0; j < mixt->nb_component; j++) {
      if (j == cstate)
      state_array[t][j][cnode]= 1.0;
      else
      state_array[t][j][cnode]= .0;
      for(i= 0; i < mixt->nb_component; i++)
      state_pair_array[t][i][j][cnode]= .0;
      // this quantity will be added to transition reestimation quantities
      }

      } // end for t
      delete state_restoration;
      state_restoration= NULL;
      }*/

      /* if (saem_exponent != .0)
     saem_coef= 1./(double)pow(iter+1, saem_exponent); */

      // accumulation of the reestimation quantities for initial distribution
      // and transition probabilities

      for(i = 0; i < mixt->nb_component; i++)
    for(n = 0; n < nb_vector; n++) {
      /* if ((algorithm != FORWARD_BACKWARD) && (saem_exponent != .0))
         (1-saem_coef)*cond_prob[n][i]
         + saem_coef*state_array[n][i];
         if ((algorithm != FORWARD_BACKWARD) && (saem_exponent == .0))
         state_array[n][i];
         }*/
      weight_reestim->frequency[i]+= cond_prob[n][i];

      // accumulation of the reestimation quantities for observation distributions
      for(var = 0; var < mixt->nb_var; var++) {
        val = int_vector[n][var];
        // if ((algorithm == FORWARD_BACKWARD))
        observation_reestim[var][i]->frequency[val] += cond_prob[n][i];
        /* if ((algorithm != FORWARD_BACKWARD) && (saem_exponent != .0))
           (1-saem_coef)* cond_prob[n][i] + saem_coef*state_array[n][i];
           if ((algorithm != FORWARD_BACKWARD) && (saem_exponent == .0))
           state_array[n][i];*/
      }
    }

      if (likelihood != D_INF) {
    if (likelihood < previous_likelihood)
      nb_likelihood_decrease++;
    else
      nb_likelihood_decrease = 0;
    // save best parameter for restoration algorithms
    if ((state_simulation) && (likelihood > best_likelihood)) {
      *best_mixt= *mixt;
      best_likelihood= likelihood;
    }
      }

      // reestimation of the weights
      /* for (i = 0; i < nb_component; i++) {
     mixt->weight->mass[i] = (double)mixt_histo->weight->frequency[k] /
     (double)mixt_histo->weight->nb_element;
     }*/

      reestimation(mixt->nb_component, weight_reestim->frequency ,
           mixt->weight->mass, MIN_PROBABILITY, false);


      // reestimation of the observation distributions

      for(var = 0; var < mixt->nb_var; var++) {
    if (mixt->npcomponent[var] != NULL)
      for(j = 0; j < mixt->nb_component; j++)
        reestimation((int)get_max_value(var)+1, observation_reestim[var][j]->frequency,
             mixt->npcomponent[var]->observation[j]->mass,
             MIN_PROBABILITY, false);

    else { // (mixt->pcomponent[var] != NULL)
      for(j = 0; j < mixt->nb_component; j++) {
        observation_reestim[var][j]->nb_value_computation();
        observation_reestim[var][j]->offset_computation();
        observation_reestim[var][j]->nb_element_computation();
        observation_reestim[var][j]->max_computation();
        observation_reestim[var][j]->mean_computation();
        observation_reestim[var][j]->variance_computation();

        hobservation->update(observation_reestim[var][j],
                 MAX((int)(observation_reestim[var][j]->nb_element
                       * MAX(sqrt(observation_reestim[var][j]->variance), 1.)
                       * MIXTURE_COEFF),
                     MIN_NB_ELEMENT));
        if ((force_param == NULL) || (!force_param[var]))
          observation_likelihood
        = hobservation->Reestimation<int>::type_parametric_estimation(mixt->pcomponent[var]->observation[j],
                                          0, true,
                                          OBSERVATION_THRESHOLD);
        // above instruction allows the type of the distribution to vary

        else
          observation_likelihood
        = hobservation->Reestimation<int>::parametric_estimation(mixt->pcomponent[var]->observation[j],
                                     0, true, OBSERVATION_THRESHOLD);
        // above instruction prevents the type of the distribution to vary
        // (not suitable for an automatic initialization of pcomponent of UNIFORM type

        if (observation_likelihood == D_INF)
          min_likelihood = D_INF;
        else {
          mixt->pcomponent[var]->observation[j]->computation((int)get_max_value(var)+1,
                                 OBSERVATION_THRESHOLD);

          if (mixt->pcomponent[var]->observation[j]->ident == BINOMIAL)
        for(i = mixt->pcomponent[var]->observation[j]->nb_value; i < (int)get_max_value(var)+1; i++)
          mixt->pcomponent[var]->observation[j]->mass[i]= 0.;
        }
      }
    }
      } // end for (var)

#     ifdef MESSAGE
      os << STAT_label[STATL_ITERATION] << " " << iter << "   "
     << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << endl;
#     endif

#     ifdef DEBUG
      if (iter % 5 == 0)
    cout << *mixt;
#     endif

    } // end do
    while ((likelihood != D_INF) && (((nb_iter == I_DEFAULT) && (iter < DISCRETE_MIXTURE_NB_ITER) &&
                      (((likelihood - previous_likelihood) / -likelihood > MVMIXTURE_LIKELIHOOD_DIFF) ||
                       (min_likelihood == D_INF) || (nb_likelihood_decrease == 1))) ||
                     ((nb_iter != I_DEFAULT) && (iter < nb_iter))));

    if (likelihood != D_INF) {

#     ifdef MESSAGE
      if (iter > 1)
    os << "\n" << iter << " " << STAT_label[STATL_ITERATIONS] << endl;
      else
    os << "\n" << iter << " " << STAT_label[STATL_ITERATION] << endl;
#     endif

      // reestimation of the weights
      /* for (i = 0; i < nb_component; i++) {
     mixt->weight->mass[i] = (double)mixt_histo->weight->frequency[k] /
     (double)mixt_histo->weight->nb_element;
     }*/

      reestimation(mixt->nb_component, weight_reestim->frequency ,
           mixt->weight->mass, MIN_PROBABILITY, false);


      // reestimation of the observation distributions

      for(var = 0; var < mixt->nb_var; var++) {
    if (mixt->npcomponent[var] != NULL)
      for(j = 0; j < mixt->nb_component; j++)
        reestimation((int)get_max_value(var)+1, observation_reestim[var][j]->frequency,
             mixt->npcomponent[var]->observation[j]->mass,
             MIN_PROBABILITY, false);

    else { // (mixt->pcomponent[var] != NULL)
      for(j = 0; j < mixt->nb_component; j++) {
        observation_reestim[var][j]->nb_value_computation();
        observation_reestim[var][j]->offset_computation();
        observation_reestim[var][j]->nb_element_computation();
        observation_reestim[var][j]->max_computation();
        observation_reestim[var][j]->mean_computation();
        observation_reestim[var][j]->variance_computation();

        hobservation->update(observation_reestim[var][j],
                 MAX((int)(observation_reestim[var][j]->nb_element
                       * MAX(sqrt(observation_reestim[var][j]->variance), 1.)
                       * MIXTURE_COEFF),
                     MIN_NB_ELEMENT));
        if ((force_param == NULL) || (!force_param[var]))
          observation_likelihood
        = hobservation->Reestimation<int>::type_parametric_estimation(mixt->pcomponent[var]->observation[j],
                                          0, true,
                                          OBSERVATION_THRESHOLD);
        // above instruction allows the type of the distribution to vary

        else
          observation_likelihood
        = hobservation->Reestimation<int>::parametric_estimation(mixt->pcomponent[var]->observation[j],

                                     0, true, OBSERVATION_THRESHOLD);
        // above instruction prevents the type of the distribution to vary
        // (not suitable for an automatic initialization of pcomponent of UNIFORM type

        if (observation_likelihood == D_INF)
          min_likelihood = D_INF;
        else {
          mixt->pcomponent[var]->observation[j]->computation((int)get_max_value(var)+1,
                                 OBSERVATION_THRESHOLD);

          if (mixt->pcomponent[var]->observation[j]->ident == BINOMIAL)
        for(i = mixt->pcomponent[var]->observation[j]->nb_value; i < (int)get_max_value(var)+1; i++)
          mixt->pcomponent[var]->observation[j]->mass[i]= 0.;
        }
      }
    }
      } // end for (var)

    } // end if (likelihood != D_INF)

      // deallocation of the arrays

    /* if ((algorithm == VITERBI) || (algorithm == FORWARD_BACKWARD_SAMPLING) ||
       (algorithm == GIBBS_SAMPLING)) {
       for(n = 0; n < nb_vector; n++) {
       for(j = 0; j < mixt->nb_component; j++) {
       delete [] state_array[n];
       state_array[n] = NULL;
       }
       delete [] state_array;
       state_array = NULL;
       }
       }*/

    // if ((algorithm == FORWARD_BACKWARD) || (saem_exponent != .0))
    // {
    for(n = 0; n < nb_vector; n++) {
      delete [] output_cond[n];
      output_cond[n] = NULL;
    }

    delete [] output_cond;
    output_cond = NULL;
    // }

    for(var = 0; var < mixt->nb_var; var++) {
      for(j = 0; j < mixt->nb_component; j++) {
    delete observation_reestim[var][j];
    observation_reestim[var][j] = NULL;
      }
      delete [] observation_reestim[var];
      observation_reestim[var] = NULL;
    }
    delete [] observation_reestim;
    observation_reestim = NULL;
    delete weight_reestim;
    weight_reestim = NULL;

    delete hobservation;
    hobservation= NULL;

    if (likelihood == D_INF) {
      delete mixt;
      mixt = NULL;
      error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
    }
    else {
      if ((state_simulation) && (likelihood > best_likelihood)) {
    *best_mixt = *mixt;
    best_likelihood = likelihood;
      }
      if (state_simulation) {
    *mixt= *best_mixt;
    delete mixt_data;
    mixt_data = NULL;
      }

      /* if ((state_trees == FORWARD_BACKWARD) || (state_trees == VITERBI)) {
     if (mixt->mixture_data != NULL)
     delete mixt->mixture_data;*/

      mixt->mixture_data = mixt->cluster(error,  *this, VITERBI);
      mixt_data = mixt->mixture_data;

      /* switch (state_trees) {

      case FORWARD_BACKWARD :
      {
      mixt_data->build_state_trees();
      mixt_data->hidden_likelihood=
      mixt->upward_downward(*mixt_data, max_marginal_entropy, entropy1,
      likelihood, path, I_DEFAULT, NULL, 'a', 0);
      if (path != NULL)
      {
      delete path;
      path= NULL;
      }
      break;
      };

      case VITERBI :
      {
      mixt->create_cumul();
      mixt->log_computation();

      mixt_data->build_state_trees();
      mixt_data->hidden_likelihood= mixt->viterbi(*mixt_data);

      mixt->remove_cumul();
      break;
      }
      }*/

      mixt_data->nb_component = mixt->nb_component;

      // computation of the parametric observation distributions
      // for mixt_data->component[0] is not to be used
      delete [] mixt_data->component[0];
      mixt_data->component[0] = NULL;
      for(var = 0; var < mixt->nb_var; var++)
    if (mixt->pcomponent[var] != NULL) {
      for(j = 0; j < mixt->nb_component; j++)
        mixt->pcomponent[var]->observation[j]->computation(mixt_data->component[var+1][j]->nb_value,
                                   OBSERVATION_THRESHOLD);
    }

#       ifdef MESSAGE
      cout << "\n" << STAT_label[STATL_LIKELIHOOD] << ": " << mixt->likelihood_computation(*this, true);
      cout << endl;
      // << " | " << mixt->Hidden_markov_out_tree::state_likelihood_computation(*mixt_data) << endl;
#       endif

    }
    /*
      else // state_trees != FORWARD_BACKWARD and VITERBI
      {
      if (mixt->mixture_data != NULL)
      delete mixt->mixture_data;

      mixt->mixture_data= new Hidden_markov_tree_data(*this);
      mixt_data= mixt->mixture_data;
      // mixt_data->state_variable_init(INT_VALUE);
      }*/
    all_states_used = true;
    for (j = 0; j < mixt->nb_component; j++) {
      if (mixt_data->component[1][j]->nb_value == 0) {
    all_states_used = false;
    break;
      }
    }


    if (!all_states_used) {

      ostringstream error_message;
      error_message << STAT_label[STATL_STATE] << " " << j << ": "
            << STAT_error[STATR_NOT_PRESENT];
      error.update((error_message.str()).c_str());
      delete mixt->mixture_data;
      mixt->mixture_data = NULL;
    }


    for(var = 0; var < mixt->nb_var; var++)
      if (mixt->npcomponent[var] != NULL)
    for(j = 0; j < mixt->nb_component; j++) {
      mixt->npcomponent[var]->observation[j]->cumul_computation();
      mixt->npcomponent[var]->observation[j]->max_computation();
    }
  }

  if (state_simulation) {
    if (best_mixt != NULL) {
      delete best_mixt;
      best_mixt = NULL;
    }
  }

  mixt_data= NULL;
  reestim= NULL;

  return mixt;
}

/*--------------------------------------------------------------*
 *
 *  Estimation du nombre de composantes et des parametres
 *  d'un melange de lois par l'algorithme EM.
 *
 *  arguments : reference sur un objet StatError, nombre de composantes,
 *  nombre d'iterations maximal, flag sur l'utilisation forcee ou
 *  non de lois d'observation parametriques (pour chaque variable)
 *
 *--------------------------------------------------------------*/

MultivariateMixture* Vectors::mixture_estimation(StatError &error, std::ostream& os ,
                                                 int nb_component, int nb_iter,
                                                 bool *force_param) const {

  // note: length of force_param must be checked before call
  bool status= true;
  register int var;
  boost::scoped_array<int> nb_value(new int[nb_variable]);
  MultivariateMixture *imixt = NULL, *mixt = NULL;

  error.init();

  if ((nb_component < 2) || (nb_component > DISCRETE_MIXTURE_NB_COMPONENT)) {
    status= false;
    error.update(STAT_error[STATR_NB_DISTRIBUTION]);
  }

  if (status) {
    for(var = 0; var < nb_variable; var++) {
      nb_value[var] = marginal_distribution[var]->nb_value;
    }

    // initial Mixture

    imixt = new MultivariateMixture(nb_component, nb_variable,
                                    nb_value.get(), force_param);

    imixt->init();

    mixt = mixture_estimation(error, os, *imixt, nb_iter);

    delete imixt;
    imixt = NULL;
  }
  return mixt;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par un melange de lois.
 *
 *  arguments : reference sur un objet StatError, effectif.
 *
 *--------------------------------------------------------------*/

MultivariateMixtureData* MultivariateMixture::simulation(StatError &error ,
                                                         int nb_element) const

{
  int k , n, var;
  int value;
  int *iidentifier, **iint_vector;
  MultivariateMixtureData *mixt_data = NULL;
  Vectors *vec = NULL;
  FrequencyDistribution *hweight = NULL;
  FrequencyDistribution ***hcomponent = NULL;


  error.init();

  if ((nb_element < 1) || (nb_element > DIST_NB_ELEMENT))
    error.update(STAT_error[STATR_SAMPLE_SIZE]);
  else {

    // creation d'un objet MultivariateMixtureData
    iidentifier = new int[nb_element];
    iint_vector = new int*[nb_element];
    hweight = new FrequencyDistribution(nb_component);
    hcomponent = new FrequencyDistribution**[nb_var];

    for (var = 0; var < nb_var; var++) {
      hcomponent[var] = new FrequencyDistribution*[nb_component];
      for (k = 0; k < nb_component; k++)
    hcomponent[var][k] = new FrequencyDistribution(nb_element);
    }

    for (n = 0; n < nb_element; n++) {

      iint_vector[n] = new int[nb_var];
      iidentifier[n] = n+1;

      // ponderation (simulation de l'etat)

      k = weight->simulation();
      (hweight->frequency[k])++;

      // composante
      for (var = 0; var < nb_var; var++) {
    if (pcomponent[var] != NULL)
      value = pcomponent[var]->observation[k]->simulation();
    else
      value = npcomponent[var]->observation[k]->simulation();
    (hcomponent[var][k]->frequency[value])++;
    iint_vector[n][var] = value;
      }
    } // end for (n)

    vec = new Vectors(nb_element, iidentifier, nb_var, iint_vector);

    for (n = 0; n < nb_element; n++) {
      delete [] iint_vector[n];
      iint_vector[n] = NULL;
    }

    // calcul des caracteristiques des lois empiriques
    for (var = 0; var < nb_var; var++) {
      for (k = 0; k < nb_component; k++) {
    hcomponent[var][k]->nb_value_computation();
    hcomponent[var][k]->offset_computation();
    hcomponent[var][k]->nb_element_computation();
    hcomponent[var][k]->max_computation();
    hcomponent[var][k]->mean_computation();
    hcomponent[var][k]->variance_computation();
      }
    }

    hweight->nb_value_computation();
    hweight->offset_computation();
    hweight->nb_element_computation();
    hweight->max_computation();
    hweight->mean_computation();
    hweight->variance_computation();

    delete [] iint_vector;
    iint_vector = NULL;
    delete [] iidentifier;
    iidentifier = NULL;

    // extraction des caracteristiques des lois empiriques
    mixt_data = new MultivariateMixtureData(*vec, nb_component);
    mixt_data->mixture = new MultivariateMixture(*this , false);

    mixt_data->weight = hweight;
    mixt_data->component = hcomponent;
  }
  return mixt_data;
}

/*--------------------------------------------------------------*
 *
 *  Ajout des etats restaures en tant que variable
 *
 *  arguments : reference sur un objet StatError, sur un objet Vectors,
 *  algorithme de restauration, et flag pour ajouter ou non l'entropie
 *
 *--------------------------------------------------------------*/

MultivariateMixtureData* MultivariateMixture::cluster(StatError &error,  const Vectors &vec,
                                                      int algorithm, bool add_state_entropy) const {

  int n, var, s, k;
  int nb_vector = vec.get_nb_vector(),
    nb_variable = vec.get_nb_variable(),
    nb_res_variable = nb_variable + 1, // nombre de variables du resultat
    // nombre de variables entieres et reelles du resultat
    nb_int_variable = 0, nb_real_variable = 0; 
  int **iint_vector = NULL;
  double **ireal_vector = NULL;  
  double **posterior_dist = NULL; // probabilites a posteriori des etats
  int *itypes = NULL; // type des variables
  std::vector<int> *states = NULL;
  MultivariateMixtureData* clusters_vec = NULL;
  Vectors* state_vec = NULL;
  FrequencyDistribution *hweight = NULL, ***hcomponent = NULL;

  if (add_state_entropy) {
    nb_res_variable++;
    nb_real_variable++;
  }

  itypes = new int[nb_res_variable];
  itypes[0] = STATE; // states
  for (var = 0; var < nb_variable; var++) {
    itypes[var+1] = vec.get_type(var);
    if ((itypes[var+1] == INT_VALUE) || (itypes[var+1] == STATE))
      nb_int_variable++;
    else
      nb_real_variable++;
  }

  nb_int_variable++; // state variable

  if (add_state_entropy)  
    itypes[nb_res_variable-1] = REAL_VALUE;
  
  // calcul des etats et leur loi a posteriori
  states = state_computation(error, vec, posterior_dist, algorithm, I_DEFAULT);

  if (states != NULL) {
    iint_vector = new int*[nb_vector];
    for (n = 0; n < nb_vector; n++) {
      iint_vector[n] = new int[nb_res_variable];
      iint_vector[n][0] = (*states)[n];
    }
    if (nb_real_variable > 0) {
      ireal_vector = new double*[nb_vector];
      for (n = 0; n < nb_vector; n++) {
	ireal_vector[n] = new double[nb_res_variable];
      }
    }

    for (var = 1; var < nb_res_variable; var++) {
      for (n = 0; n < nb_vector; n++) {
	if ((itypes[var] == INT_VALUE) || (itypes[var] == STATE))
	  iint_vector[n][var] = vec.int_vector[n][var-1];
	else {
	  if (add_state_entropy && (var == nb_res_variable-1)) {
	    // compute entropy
	    ireal_vector[n][nb_real_variable-1] = 0;
	    for (s = 0; s < nb_component; s++) {
	      if (log(posterior_dist[n][s]) > D_INF) 
		ireal_vector[n][nb_real_variable-1] -= 
		  posterior_dist[n][s] * log(posterior_dist[n][s]); 
	      }
	  }
	  else
	    ireal_vector[n][var] = vec.real_vector[n][var-1];
	}
      }
    }
    
    state_vec = new Vectors(nb_vector, vec.identifier, nb_res_variable, itypes, 
			    iint_vector, ireal_vector);

    for (n = 0; n < nb_vector; n++) {
      delete [] iint_vector[n];
      iint_vector[n] = NULL;
    }
    delete [] iint_vector;
    iint_vector = NULL;
    delete [] itypes;
    itypes = NULL;
    if (ireal_vector != NULL) {
      for (n = 0; n < nb_vector; n++) {
	delete [] ireal_vector[n];
	ireal_vector[n] = NULL;
      }
      delete [] ireal_vector;
      ireal_vector = NULL; 
    }
    if (add_state_entropy) {
      for (n = 0; n < nb_vector; n++) {
	delete [] posterior_dist[n]; 
	posterior_dist[n] = NULL; }
      delete [] posterior_dist;
      posterior_dist = NULL;
    }

    state_vec->type[0] = STATE;
    clusters_vec = new MultivariateMixtureData(*state_vec, nb_component);
    clusters_vec->mixture = new MultivariateMixture(*this, false);

    assert(clusters_vec->type[0] = STATE);

    // calcul des lois empiriques d'observation

    hweight = clusters_vec->weight;
    hcomponent = clusters_vec->component;

    for (n = 0; n < nb_vector; n++) {
      s = (*states)[n];
      (hweight->frequency[s])++;
      for (var = 1; var < nb_var+1; var++) {
	(hcomponent[var][s]->frequency[vec.int_vector[n][var-1]])++;
      }
    } // end for (n)

    delete [] hcomponent[0];
    hcomponent[0] = NULL;

    for (var = 1; var <= nb_var; var++) {
      for (k = 0; k < nb_component; k++) {
    hcomponent[var][k]->nb_value_computation();
    hcomponent[var][k]->offset_computation();
    hcomponent[var][k]->nb_element_computation();
    hcomponent[var][k]->max_computation();
    hcomponent[var][k]->mean_computation();
    hcomponent[var][k]->variance_computation();
      }
    }

    hweight->nb_value_computation();
    hweight->offset_computation();
    hweight->nb_element_computation();
    hweight->max_computation();
    hweight->mean_computation();
    hweight->variance_computation();

    delete states;
    states = NULL;
    delete state_vec;
    state_vec = NULL;
  }
  return clusters_vec;

}


};  // namespace stat_tool
