/*  ----------------------------------------------------------------------------
 *
 *
 *       Copyright 2005-2009 UMR DAP
 *
 *       File author(s): F. Chaubert-Pereira (chaubert@cirad.fr)
 *
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


/*---------------------------------------------------------------------
 *
 * Correspond aux fonctions propres aux semi-Markov switching 
 * linear mixed model
 * (modèle avec effets aléatoires temporels)
 *  
 *--------------------------------------------------------------------*/

#include <math.h>
#include <iostream>
#include <iomanip>
#include <limits.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_eigen.h>
#include "stat_tool/stat_label.h"
#include "tool/util_math.h"
#include "stat_tool/stat_tools.h"
#include "tool/config.h"
#include "stat_tool/curves.h"
#include "sequence_analysis/renewal.h"
#include "stat_tool/markovian.h"
#include "sequence_analysis/sequences.h"
#include "sequence_analysis/sequence_label.h"
#include "stat_tool/vectors.h"
#include "sequence_analysis/tops.h"
#include "continuous_histo.h"
#include "switching_sequence.h"
#include "switching_process.h"
#include "markov_switching.h"

#include "stat_tool/distribution_reestimation.h"

using namespace std;

extern void convert_matrix_int (int **matrix, int nb_row, int nb_col, gsl_matrix* convert_matrix);
extern void convert_matrix_diag_int (int *matrix, int nb_row, gsl_matrix* convert_matrix); 
extern void convert_vector_int (int *vect, int length, gsl_vector* convert_vect);
extern void convert_array2D_int (gsl_matrix *tab, int nb_row, int nb_col, int **convert_tzb);
extern void convert_array1D_diag_int (gsl_matrix *tab, int length, int *convert_tab); 
extern void convert_array1D_int (gsl_vector *tab, int length, int *convert_tab);

extern void convert_matrix_double (double **matrix, int nb_row, int nb_col, gsl_matrix* convert_matrix);
extern void convert_matrix_diag_double (gsl_vector *matrix, int nb_row, gsl_matrix* convert_matrix); 
extern void convert_vector_double (double *vect, int length, gsl_vector* convert_vect);
extern void convert_array2D_double (gsl_matrix *tab, int nb_row, int nb_col, double **convert_tab);
extern void convert_array1D_diag_double (gsl_matrix *tab, int length, double *convert_tab);
extern void convert_array1D_double (gsl_vector *tab, int length, double *convert_tab);

extern void cumul_computation(int nb_value, const double *pmass, double *pcumul);
extern int cumul_method(int nb_value, const double *cumul, double scale=1.);
extern void log_computation(int nb_value, const double *pmass, double *plog);
extern char* label(const char *file_name);

extern double random_normal(double imean, double ivariance);


extern double random_unif();

extern double interval_bisection(Reestimation<double> *distribution_reestim, 
 				 Reestimation<double> *length_bias_reestim);


/*--------------------------------------------------------------*
 *
 *  Simulation par un semi-Markov switching linear mixed model.
 *
 *  arguments : reference sur un objet Format_error, histogramme 
 *              des longueurs des sequences, nombre de covariables,
 *              nombre d'effets aléatoires, présence/absence de 
 *              constante, variance résiduelle, paramètres de 
 *              regression, variance aléatoire, covariables,
 *              effets aléatoires temporels, paramètres d'index,
 *              flag sur le calcul des lois de comptage, flag 
 *              calcul d'une divergence de Kullback-Leibler.
 *
 *--------------------------------------------------------------*/

Markov_switching_data* Markov_switching::simulation_year_effect(Format_error &error, const Histogram &hlength, int inb_covariable, 
								int inb_random, int iconstant, double *iresidual_variance, 
								double **iregression, double **irandom_variance, double ***icovar,  
								double *iyear_effect, int **iindex, int T, bool counting_flag, 
								bool divergence_flag) const
{
  bool status = true;
  register int i , j , k, m;
  int cumul_length ,occupancy, *pstate ;
  double **poutput;
  double *pcovar, *ccovar, *pyear_effect, *cyear_effect;
  int *pindex, *cindex;
  int nb_covar;
  Markov_switching *sw_markov;
  Markov_switching_data *seq;
  
  seq = 0;
  error.init();

  if ((hlength.nb_element < 1) || (hlength.nb_element > NB_SEQUENCE)) {
    status = false;
    error.update(SEQ_error[SEQR_NB_SEQUENCE]);
  }
  if (hlength.offset < 2) {
    status = false;
    error.update(SEQ_error[SEQR_SHORT_SEQUENCE_LENGTH]);
  }
  if (hlength.nb_value - 1 > MAX_LENGTH) {
    status = false;
    error.update(SEQ_error[SEQR_LONG_SEQUENCE_LENGTH]);
  }

  if (status) {
    cumul_length = 0;
    for (i = hlength.offset;i < hlength.nb_value;i++) {
      cumul_length += i * hlength.frequency[i];
    }

    if (cumul_length > CUMUL_LENGTH) {
      status = false;
      error.update(SEQ_error[SEQR_CUMUL_SEQUENCE_LENGTH]);
    }
  }

  if (status) {

    // initialisations

    seq = new Markov_switching_data(2, hlength , inb_covariable, inb_random, 0, T, iconstant, false);
    seq->type[0] = STATE;
    seq->sw_markov = new Markov_switching(*this , false);

    sw_markov = seq->sw_markov;

    sw_markov->create_cumul();
    sw_markov->cumul_computation();

    poutput = new double*[1];
    
    for (i = 0; i < sw_markov->nb_state; i++) {
      sw_markov->sw_process[1]->observation[i]->Set_random_variance(irandom_variance[i]);
      sw_markov->sw_process[1]->observation[i]->Set_residual_variance(iresidual_variance[i]);
      if (!iconstant){
	sw_markov->sw_process[1]->observation[i]->nb_covariable = inb_covariable +1;
	nb_covar = inb_covariable + 1;
      }
      else {
	sw_markov->sw_process[1]->observation[i]->nb_covariable = inb_covariable;
	nb_covar = inb_covariable;
      }
      sw_markov->sw_process[1]->observation[i]->Set_regression(iregression[i]);
    }

    for (i = 0;i <seq->nb_sequence;i++) {
      for (j = 0;j < seq->length[i];j++) {
	pcovar = seq->covar[i][j];
	ccovar = icovar[i][j];
	if (!iconstant) {
	  *pcovar++ = 1;
	}

	if(inb_covariable > 0){
	  for (k = 0;k < nb_covar;k++) {
	    *pcovar++ = *ccovar++;
	  }
	}
      }
    }
    

    for (i = 0;i < seq->nb_sequence;i++) {
      pindex = seq->index[i];
      cindex = iindex[i];
      for (j = 0;j < seq->length[i];j++) {
	*pindex++ = *cindex++;
      }
    }

    pyear_effect = seq->year_effect;
    cyear_effect = iyear_effect; 

    for (i = 0;i < T;i++) {
      *pyear_effect++ = *cyear_effect++;  
    }
    
    for (i = 0;i < seq->nb_sequence;i++) {
      pstate = seq->int_sequence[i][0];
      *pstate = cumul_method(sw_markov->nb_state , sw_markov->cumul_initial);
      poutput[0] = seq->real_sequence[i][1];
            
      j = 0;
      do {
	if (j > 0) {
	  pstate++;
	  *pstate = cumul_method(sw_markov->nb_state , sw_markov->cumul_transition[*(pstate - 1)]);
	}
	
	switch (sw_markov->state_subtype[*pstate]) {
	  
	case SEMI_MARKOVIAN : {
	  if ((sw_markov->type == 'e') && (j == 0)) {
	    occupancy = sw_markov->forward[*pstate]->simulation();
	  }
	  else {
	    occupancy = sw_markov->nonparametric_process[0]->sojourn_time[*pstate]->simulation();
	  }
	  
	  if (j + occupancy > seq->length[i]) {
	    occupancy = seq->length[i] - j;
	  }
	  break;
	}
	  
	case MARKOVIAN : {
	  if (sw_markov->transition[*pstate][*pstate] < 1.) {
	    occupancy = 1;
	  }
	  else {
	    occupancy = seq->length[i] - j;
	  }
	  break;
	}
	}
	
	j += occupancy;
	
	for (k = 1;k < occupancy;k++) {
	  pstate++;
	  *pstate = *(pstate - 1);
	}
	
	for(k = 0; k < occupancy; k++){
	  if (sw_markov->sw_process[1]) {
	    *poutput[0]++ = sw_markov->sw_process[1]->observation[*pstate]->simulation(seq->covar[i][j])
	      + sqrt(irandom_variance[*pstate][inb_random-1])*seq->year_effect[seq->index[i][j]-1];
	  }
	}
      }
      while (j < seq->length[i]);
    }
	  
    sw_markov->remove_cumul();
    
    delete [] poutput;
      
    // extraction des caracteristiques des sequences simulees

    for (i = 0;i < seq->nb_variable;i++) {
      seq->max_value_computation(i);
      seq->build_marginal_histogram(i);
    }

    seq->build_transition_count(*sw_markov);

    if (!divergence_flag) {
      sw_markov->characteristic_computation(*seq , counting_flag);

      // calcul de la vraisemblance

      //      seq->likelihood = sw_markov->likelihood_computation(*seq); 

    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par un semi-Markov switching linear mixed model.
 *
 *  arguments : reference sur un objet Format_error, nombre et  
 *              longueurs des sequences, nombre de covariables,
 *              nombre d'effets aléatoires, présence/absence de 
 *              constante, variance résiduelle, paramètres de 
 *              regression, variance aléatoire, covariables,
 *              effets aléatoires temporels, paramètres d'index, 
 *              flag sur le calcul des lois de comptage.
 *
 *--------------------------------------------------------------*/

Markov_switching_data* Markov_switching::simulation_year_effect(Format_error &error , int nb_sequence ,
								int length , int inb_covariable, int inb_random,
								int iconstant,  double *iresidual_variance,
								double **iregression, double **irandom_variance, 
								double ***icovar, double *iyear_effect, int **iindex, 
								int T, bool counting_flag) const
{
  
  register int i, j;
  bool status = true;
  Markov_switching_data *seq;
  
  seq = 0;
  error.init();
  
  if ((nb_sequence < 1) || (nb_sequence > NB_SEQUENCE)) {
    status = false;
    error.update(SEQ_error[SEQR_NB_SEQUENCE]);
  }
  if (length < 2) {
    status = false;
    error.update(SEQ_error[SEQR_SHORT_SEQUENCE_LENGTH]);
  }
  if (length > MAX_LENGTH) {
    status = false;
    error.update(SEQ_error[SEQR_LONG_SEQUENCE_LENGTH]);
  }

  if (status) {
    Histogram hlength(length + 1);

    hlength.nb_element = nb_sequence;
    hlength.offset = length;
    hlength.max = nb_sequence;
    hlength.mean = length;
    hlength.variance = 0.;
    hlength.frequency[length] = nb_sequence;

    seq = simulation_year_effect(error , hlength , inb_covariable, inb_random, iconstant, iresidual_variance, iregression, 
				 irandom_variance,  icovar, iyear_effect, iindex, T, counting_flag, false);
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par un semi-Markov switching linear mixed model.
 *
 *  arguments : reference sur un objet Format_error, nombre de 
 *              séquences, référence sur un objet Switching_sequence,
 *              nombre de covariables, nombre d'effets aléatoires, 
 *              présence/absence de constante, variance résiduelle, 
 *              paramètres de regression, variance aléatoire, covariables,
 *              effets aléatoires temporels, paramètres d'index, nombre 
 *              de temps distincts, flag sur le calcul des lois de comptage.
 *
 *--------------------------------------------------------------*/

Markov_switching_data* Markov_switching::simulation_year_effect(Format_error &error , int nb_sequence ,
								const Switching_sequence &iseq ,
								int inb_covariable, int inb_random, int iconstant,
								double *iresidual_variance, double **iregression,
								double **irandom_variance, double ***icovar,
								double *iyear_effect, int **iindex, int T, bool counting_flag) const
{
  Histogram *hlength;
  Markov_switching_data *seq;

  error.init();

  if ((nb_sequence < 1) || (nb_sequence > NB_SEQUENCE)) {
    seq = 0;
    error.update(SEQ_error[SEQR_NB_SEQUENCE]);
  }

  else {
    hlength = iseq.hlength->frequency_scale(nb_sequence);

    seq = simulation_year_effect(error , *hlength , inb_covariable, inb_random, iconstant, iresidual_variance, iregression, 
				 irandom_variance, icovar, iyear_effect, iindex, T, counting_flag, false);
    delete hlength;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par un semi-Markov switching linear mixed model.
 *
 *  arguments : reference sur un objet Format_error, histogramme
 *              des longueurs des séquences, nombre de covariables, 
 *              nombre d'effets aléatoires, présence/absence de 
 *              constante, variance résiduelle, paramètres de 
 *              regression, variance aléatoire, covariables,
 *              effets aléatoires temporels, paramètres d'index, 
 *              nombre de temps distincts, flag sur le calcul 
 *              des lois de comptage, flag calcul d'une divergence 
 *              de Kullback-Leibler.
 *
 *--------------------------------------------------------------*/

Markov_switching_data* Markov_switching::simulation_year_effect(Format_error &error, const Histogram &hlength, int inb_covariable, 
								int inb_random, int iconstant, double *iresidual_variance, 
								double **iregression, double **irandom_variance, double ***icovar,
								double *iyear_effect, int **iindex, int T, bool counting_flag ,
								bool divergence_flag, bool hidden) const
{
  Switching_sequence *observ_seq;
  Markov_switching_data *seq;

  seq = Markov_switching::simulation_year_effect(error, hlength, inb_covariable, inb_random, iconstant, iresidual_variance, 
						 iregression, irandom_variance, icovar, iyear_effect, iindex, T, counting_flag, 
						 divergence_flag);

  if ((seq) && (!divergence_flag)) {
    observ_seq = seq->remove_variable_1();
    seq->hidden_likelihood = likelihood_computation(*observ_seq, 0,I_DEFAULT);
    delete observ_seq;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par un semi-Markov switching linear mixed model.
 *
 *  arguments : reference sur un objet Format_error, nombre et  
 *              longueurs des sequences, nombre de covariables,
 *              nombre d'effets aléatoires, présence/absence de 
 *              constante, variance résiduelle, paramètres de 
 *              regression, variance aléatoire, covariables,
 *              effets aléatoires temporels, paramètres d'index, 
 *              nombre de temps distincts, flag sur le calcul 
 *              des lois de comptage.
 *
 *--------------------------------------------------------------*/

Markov_switching_data* Markov_switching::simulation_year_effect(Format_error &error, int nb_sequence, int length, int inb_covariable, 
								int inb_random, int iconstant, double *iresidual_variance, 
								double **iregression, double **irandom_variance, double ***icovar,
								double *iyear_effect, int **iindex, int T, bool counting_flag, 
								bool hidden) const
{
  Switching_sequence *observ_seq;
  Markov_switching_data *seq;
    
  seq = Markov_switching::simulation_year_effect(error, nb_sequence, length, inb_covariable, inb_random, iconstant, iresidual_variance, 
						 iregression, irandom_variance, icovar, iyear_effect, iindex, T, counting_flag);
  
  if (seq) {
    observ_seq = seq->remove_variable_1();
    seq->hidden_likelihood = likelihood_computation(*observ_seq, 0,  I_DEFAULT);
    delete observ_seq;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par un semi-Markov switching linear mixed model.
 *
 *  arguments : reference sur un objet Format_error, nombre de 
 *              séquences, référence sur un objet Switching_sequence,
 *              nombre de covariables, nombre d'effets aléatoires, 
 *              présence/absence de constante, variance résiduelle, 
 *              paramètres de regression, variance aléatoire, covariables,
 *              effets aléatoires temporels, paramètres d'index, nombre 
 *              de temps distincts, flag sur le calcul des lois de comptage.
 *
 *--------------------------------------------------------------*/

Markov_switching_data* Markov_switching::simulation_year_effect(Format_error &error, int nb_sequence, const Switching_sequence &iseq ,
								int inb_covariable, int inb_random, int iconstant, 
								double *iresidual_variance, double **iregression,
								double **irandom_variance, double ***icovar, double *iyear_effect, 
								int **iindex, int T, bool counting_flag, bool hidden) const

{
  Switching_sequence *observ_seq;
  Markov_switching_data *seq;
  
  seq = Markov_switching::simulation_year_effect(error , nb_sequence , iseq , inb_covariable, inb_random,iconstant, iresidual_variance, 
						 iregression, irandom_variance, icovar, iyear_effect, iindex,T,  counting_flag);
  
  if (seq) {
    observ_seq = seq->remove_variable_1();
    seq->hidden_likelihood = likelihood_computation(*observ_seq, 0, I_DEFAULT);
    delete observ_seq;
  }
  
  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'un semi-Markov switching linear
 *  mixed model a partir d'un echantillon de sequences par 
 *  l'algorithme SEM/MCEM.
 *
 *  State sequences : "forward-backward" 
 *  Random effects : sampling by Metropolis-Hastings algorithm 
 *
 *  arguments : reference sur un objet Format_error, stream, 
 *              semi-Markov switching linear mixed model initial,
 *              nombre d'itérations pour convergence, parametres 
 *              pour le nombre d'effets aléatoires  simules, type 
 *              d'estimateur pour la reestimation des lois d'occupation 
 *              des etats, variance commune à tous les états, sortie 
 *              de type R, flags sur le calcul des lois de comptage 
 *              et sur le calcul des sequences d'etats optimales, 
 *              nombre d'iterations.
 *
 *--------------------------------------------------------------*/

Markov_switching* Switching_sequence::markov_switching_stochastic_estimation_MH_year(Format_error &error, ostream &os,
										     const Markov_switching &ihsw_markov,
										     int nb_iter_conv, int min_nb_random,
										     int max_nb_random, double VarMH,
										     int estimator, bool VarCommune, 
										     bool output_R, double parameter, bool counting_flag,
										     bool state_sequence, int nb_iter,
										     int mean_computation) const
{
  bool status;
  register int i , j , k , m, n, r;
  int max_nb_value , nb_likelihood_decrease, offset, nb_value, *occupancy_nb_value,  *censored_occupancy_nb_value;
  int nb_random_simul, iter ;
  double **poutput;
  int tmp = 0;
  double likelihood = D_INF , previous_likelihood , occupancy_likelihood, observation_likelihood , min_likelihood ,
    obs_product, buff, sum, occupancy_mean, **observation, *norm, *state_norm, **forward, **state_in, ***backward, **backward1,
    *auxiliary, *reestim, *ofrequency, *lfrequency, *occupancy_survivor, *censored_occupancy_survivor;
  
  double *temp;

  Chain_reestimation<double> *chain_reestim;
  Parametric *occupancy;
  Reestimation<double> **occupancy_reestim, **length_bias_reestim, **censored_occupancy_reestim;
  Histogram *hoccupancy;
  Continuous_histo *hobservation;
  Markov_switching *hsw_markov;
  Markov_switching_data *hsw_markov_data;

  double *back;
  double **predicted;

  hsw_markov = 0;
  error.init();
 
# ifdef DEBUG
  double test[NB_STATE][4];
# endif
 
  // test nombre de valeurs observees par variable
  
  status = false;
  for (i = 0;i < nb_variable;i++) {
    if (max_value[i] > min_value[i]) {
      status = true;
      break;
    }
  }

  if (!status) {
    error.update(SEQ_error[SEQR_VARIABLE_NB_VALUE]);
  }


  if ((nb_iter != I_DEFAULT) && (nb_iter < 1)) {
    status = false;
    error.update(STAT_error[STATR_NB_ITERATION]);
  }


  if((min_nb_random < 0) || (min_nb_random>max_nb_random)){
    status = false;
    error.update(SEQ_error[SEQR_MIN_NB_STATE_SEQUENCE]);
  }


  if (status) {

    // creation de la chaine de Markov cachee

    hsw_markov = new Markov_switching(ihsw_markov , false, (int)(max_length * SAMPLE_NB_VALUE_COEFF));

#   ifdef DEBUG
    cout << *hsw_markov;
#   endif

    // initialisations

    // creation des structures de donnees de l'algorithme

    temp = new double[hsw_markov->nb_state];

    observation = new double*[max_length];
    for (i = 0;i < max_length;i++) {
      observation[i] = new double[hsw_markov->nb_state];
    }

    norm = new double[max_length];
    state_norm = new double[hsw_markov->nb_state];

    forward = new double*[max_length];
    for (i = 0;i < max_length;i++) {
      forward[i] = new double[hsw_markov->nb_state];
    }

    state_in = new double*[max_length - 1];
    for (i = 0;i < max_length - 1;i++) {
      state_in[i] = new double[hsw_markov->nb_state];
    }

    back = new double[hsw_markov->nb_state];

    backward1 = new double*[max_length];
    for (i = 0;i < max_length;i++) {
      backward1[i] = new double[hsw_markov->nb_state];
    }

    auxiliary = new double[hsw_markov->nb_state];

    predicted = new double*[max_length];
    for (i = 0;i < max_length;i++) {
      predicted[i] = new double[hsw_markov->nb_state];
    }
    
    chain_reestim = new Chain_reestimation<double>('o' , hsw_markov->nb_state , hsw_markov->nb_state);

    occupancy_nb_value = new int[hsw_markov->nb_state];
    occupancy_reestim = new Reestimation<double>*[hsw_markov->nb_state];
    
    for (i = 0;i < hsw_markov->nb_state;i++) {
      switch (hsw_markov->state_subtype[i]) {
	
      case SEMI_MARKOVIAN : {
	occupancy_nb_value[i] = hsw_markov->nonparametric_process[0]->sojourn_time[i]->alloc_nb_value;
        occupancy_reestim[i] = new Reestimation<double>(occupancy_nb_value[i]);
        break;
      }

      case MARKOVIAN : {
        occupancy_reestim[i] = 0;
        break;
      }
      }
    }

    max_nb_value = 0;
    for (i = 0;i < hsw_markov->nb_state;i++) {
      if ((hsw_markov->state_subtype[i] == SEMI_MARKOVIAN) && (occupancy_nb_value[i] > max_nb_value)) {
        max_nb_value = occupancy_nb_value[i];
      }
    }

    hoccupancy = new Histogram(max_nb_value);

    max_nb_value = 0;
    if ((hsw_markov->sw_process[1]) && (max_nb_value < marginal[0]->nb_value)) {
      max_nb_value = marginal[0]->nb_value;
    }

    poutput = new double*[nb_variable];
    

    int T = 0;//renvoie le nombre de temps distincts
    double **year_effect_predict;
    
    for (i = 0; i < nb_sequence; i++){
      if( T < index[i][length[i]-1]){
	T = index[i][length[i]-1];
      }
    }	
    
    year_effect_predict = new double *[max_nb_random+1];
    for ( i = 0; i <= max_nb_random; i++){
      year_effect_predict[i] = new double [T];
    }

    //HMM-LM ou HSMM-LM
    // initialisation des paramètres de régression, variances résiduelles, probabilités initiales et probabilités de transition,

    hsw_markov = markov_switching_estimation(error, cout, ihsw_markov, VarCommune, estimator, false, false, false);

    iter = 0;
    nb_likelihood_decrease = 0;

    do {
      iter++;
      previous_likelihood = likelihood;
      likelihood = 0.;

      // calcul du nombre d'effets aléatoires simulés

      if(min_nb_random + (int)::round(parameter * iter) < max_nb_random){
	nb_random_simul = min_nb_random + (int)::round(parameter * iter);
      }
      else {
	nb_random_simul = max_nb_random;
      } 

      backward = new double**[nb_sequence];
      for (i = 0; i < nb_sequence; i++) {
	backward[i] = new double*[max_length];
	for (j = 0;j < max_length;j++) {
	  backward[i][j] = new double[hsw_markov->nb_state];
	}
      }

      // initialisation des quantites de reestimation
      
      chain_reestim->init();
      
      for (i = 0;i < hsw_markov->nb_state;i++) {
        if (hsw_markov->state_subtype[i] == SEMI_MARKOVIAN) {
          reestim = occupancy_reestim[i]->frequency;
          for (j = 0;j < occupancy_nb_value[i];j++) {
            *reestim++ = 0.;
          }
	}
      }

#     ifdef DEBUG
      for (i = 0;i < hsw_markov->nb_state;i++) {
        for (j = 0;j < 4;j++) {
          test[i][j] = 0.;
        }
      }
#     endif

      for (i = 0;i < nb_sequence;i++) {
	
	for (j = 0;j < nb_variable;j++) {
	  poutput[j] = real_sequence[i][j];
	  
	}

	// recurrence "forward"
	
	for (j = 0;j < length[i];j++) {
          norm[j] = 0.;

          for (k = 0;k < hsw_markov->nb_state;k++) {

            // calcul des probabilites d'observation

            observation[j][k] = 1.;
	    observation[j][k] *= hsw_markov->sw_process[1]->observation[k]->density_computation(*poutput[0], covar[i][j], 0, 
												year_effect[index[i][j]-1]);
	   
	    switch (hsw_markov->state_subtype[k]) {

            // cas etat semi-markovien

            case SEMI_MARKOVIAN : {
              if (j == 0) {
                state_norm[k] = hsw_markov->initial[k];
              }
              else {
                state_norm[k] += state_in[j - 1][k] - forward[j - 1][k];
              }

	      state_norm[k] *= observation[j][k];
	      norm[j] += state_norm[k];
              break;
            }

            // cas etat markovien

            case MARKOVIAN : {
              if (j == 0) {
                forward[j][k] = hsw_markov->initial[k];
              }
              else {
                forward[j][k] = state_in[j - 1][k];
              }

              forward[j][k] *= observation[j][k];

              norm[j] += forward[j][k];

              break;
            }
            }
          }
  
	  if (norm[j] > 0.) {

	    for (k = 0;k < hsw_markov->nb_state;k++) {
              switch (hsw_markov->state_subtype[k]) {
              case SEMI_MARKOVIAN :
                state_norm[k] /= norm[j];
                break;
              case MARKOVIAN :
                forward[j][k] /= norm[j];
                break;
              }
            }

            likelihood += log(norm[j]);
          }

          else {
            likelihood = D_INF;
            break;
          }

          for (k = 0;k < hsw_markov->nb_state;k++) {

            // cas etat semi-markovien

            if (hsw_markov->state_subtype[k] == SEMI_MARKOVIAN) {
              
	      occupancy = hsw_markov->nonparametric_process[0]->sojourn_time[k];
              obs_product = 1.;
              forward[j][k] = 0.;

              if (j < length[i] - 1) {
                for (m = 1;m <= MIN(j + 1 , occupancy->nb_value   - 1);m++) {
                  obs_product *= observation[j - m + 1][k] / norm[j - m + 1];
                  if (obs_product == 0.) {
                    break;
                  }
                  if (m < j + 1) {
                    forward[j][k] += obs_product * occupancy->mass[m] * state_in[j - m][k];
                  }
                  else {
		    forward[j][k] += obs_product * occupancy->mass[m] * hsw_markov->initial[k];
                  }
                }
              }

              else {
                for (m = 1;m <= MIN(j + 1 , occupancy->nb_value - 1);m++) {
                  obs_product *= observation[j - m + 1][k] / norm[j - m + 1];
                  if (obs_product == 0.) {
                    break;
                  }

                  if (m < j + 1) {
                    forward[j][k] += obs_product * (1. - occupancy->cumul[m - 1]) * state_in[j - m][k];
                  }
                  else {
		    forward[j][k] += obs_product * (1. - occupancy->cumul[m - 1]) *
		      hsw_markov->initial[k];
                  }
                }
              }
            }
          }

          if (j < length[i] - 1) {
            for (k = 0;k < hsw_markov->nb_state;k++) {
              state_in[j][k] = 0.;
              for (m = 0;m < hsw_markov->nb_state;m++) {
                state_in[j][k] += hsw_markov->transition[m][k] * forward[j][m];
  	      }
            }
          }

          for (k = 0;k < nb_variable;k++) {
            poutput[k]++;
          }
        }

        if (likelihood == D_INF) {
          break;
        }

#       ifdef DEBUG
	if(i==60){
	  for (j = 0;j < length[i];j++) {
	    double susu = 0;
	    cout << j << " : ";
	    for (k = 0;k < hsw_markov->nb_state;k++) {
	      susu = susu + forward[j][k];
	      cout << forward[j][k] << " ";
	    }
	    cout<<": "<<susu;
	    cout << endl;
	  }
	  cout << endl;
	}
#       endif

	// recurrence "backward"

        for (j = 0;j < nb_variable;j++) {
          poutput[j]--;
        }

        j = length[i] - 1;
        for (k = 0;k < hsw_markov->nb_state;k++) {
          backward[i][j][k] = forward[j][k];

	  back[k] = forward[j][k];
          backward1[j][k] = back[k];
        }

        for (j = length[i] - 2;j >= 0;j--) {
          for (k = 0;k < nb_variable;k++) {
            poutput[k]--;
          }

          for (k = 0;k < hsw_markov->nb_state;k++) {
            auxiliary[k] = 0.;

            switch (hsw_markov->state_subtype[k]) {

            // cas etat semi-markovien

            case SEMI_MARKOVIAN : {
              occupancy = hsw_markov->nonparametric_process[0]->sojourn_time[k];
              obs_product = 1.;

              for (m = 1;m < MIN(length[i] - j , occupancy->nb_value);m++) {
                obs_product *= observation[j + m][k] / norm[j + m];
                if (obs_product == 0.) {
                  break;
                }

                if (backward1[j + m][k] > 0.) {
                  if (m < length[i] - j - 1) {
                    buff = backward1[j + m][k] * obs_product * occupancy->mass[m] /
                           forward[j + m][k];

                    // accumulation des quantites de reestimation des lois d'occupation des etats

                    occupancy_reestim[k]->frequency[m] += buff * state_in[j][k];
                  }

                  else {
                    buff = obs_product * (1. - occupancy->cumul[m - 1]);

                    // accumulation des quantites de reestimation des lois d'occupation des etats

		    for (n = m;n < occupancy->nb_value;n++) {
		      occupancy_reestim[k]->frequency[n] += obs_product * occupancy->mass[n] *
			state_in[j][k];
		    }
                  }
		  
                  auxiliary[k] += buff;
                }
              }
              break;
            }

            // cas etat markovien

            case MARKOVIAN : {
              if (backward1[j + 1][k] > 0.) {
                auxiliary[k] = backward1[j + 1][k] / state_in[j][k];
              }
              break;
            }
            }
          }

          for (k = 0;k < hsw_markov->nb_state;k++) {
            backward1[j][k] = 0.;

            for (m = 0;m < hsw_markov->nb_state;m++) {
              buff = auxiliary[m] * hsw_markov->transition[k][m] * forward[j][k];
              backward1[j][k] += buff;

              // accumulation des quantites de reestimation des probabilites de transition

              chain_reestim->transition[k][m] += buff;
            }

            switch (hsw_markov->state_subtype[k]) {

            // cas etat semi-markovien

            case SEMI_MARKOVIAN : {
              backward[i][j][k] = backward[i][j+1][k] + backward1[j][k] - auxiliary[k] * state_in[j][k];
              if (backward[i][j][k] < 0.) {
                backward[i][j][k] = 0.;
              }
              if (backward[i][j][k] > 1.) {
                backward[i][j][k] = 1.;
              }

	      back[k] = back[k] + backward1[j][k] - auxiliary[k] * state_in[j][k];
              if (back[k] < 0.) {
                back[k] = 0.;
              }
              if (back[k] > 1.) {
                back[k] = 1.;
              }
              break;
            }

            // cas etat markovien

            case MARKOVIAN : {
              backward[i][j][k] = backward1[j][k];
	      back[k] = backward1[j][k];
              break;
            }
            }
          }
        }

        // accumulation des quantites de reestimation des probabilites initiales

	for (j = 0;j < hsw_markov->nb_state;j++) {
	  chain_reestim->initial[j] += back[j];
	}
        
        // accumulation des quantites de reestimation des lois d'occupation des etats initiaux

	for (j = 0;j < hsw_markov->nb_state;j++) {
	  if ((hsw_markov->state_subtype[j] == SEMI_MARKOVIAN) && (hsw_markov->initial[j] > 0.)) {
	    occupancy = hsw_markov->nonparametric_process[0]->sojourn_time[j];
	    obs_product = 1.;
	    
	    for (k = 1;k < MIN(length[i] + 1 , occupancy->nb_value);k++) {
	      obs_product *= observation[k - 1][j] / norm[k - 1];
	      if (obs_product == 0.) {
		break;
	      }
	      
	      if (backward1[k - 1][j] > 0.) {
		if (k < length[i]) {
		  occupancy_reestim[j]->frequency[k] += backward1[k - 1][j] * obs_product *
		    occupancy->mass[k] * hsw_markov->initial[j] / forward[k - 1][j];
		}

		else {
		  for (m = k;m < occupancy->nb_value;m++) {
		    occupancy_reestim[j]->frequency[m] += obs_product * occupancy->mass[m] *
		      hsw_markov->initial[j];
		  }
		}
	      }
	    }
	  }
	}
	

#       ifdef DEBUG
	for (j = length[i] - 1;j >= 0;j--) {
	  cout << j << " : ";
	  double sum = 0.;
	  for (k = 0;k < hsw_markov->nb_state;k++) {
	    sum += back[k];
	    cout << back[k];
	    if ((hsw_markov->state_subtype[k] == SEMI_MARKOVIAN) && (j < length[i] - 1)){
	      cout << " (" << backward1[j][k] << ") ";
	    }
	  }
	  cout << "| " << sum << endl;
        }
#       endif
      }

# ifdef DEBUG
      cout<<*chain_reestim<<endl;
# endif

      //  simulation des réalisations de l'effet aléatoire
            
      for(i = 0; i < T; i++){
	year_effect_predict[0][i] = year_effect[i];
      }
      

      double U;
      for (r = 1; r <= nb_random_simul; r++){
	
       	U = random_unif();

	for(j = 0; j < T; j++){
	  year_effect_predict[r][j] = year_effect_predict[r-1][j]+ random_normal(0.,VarMH);
	}

	double prob;
	double likelihoodnum = 0.;
	double likelihoodden = likelihood;

	for (i = 0;i < nb_sequence;i++) {

	  // recurrence "forward"
	  for (j = 0;j < nb_variable;j++) {
	    poutput[j] = real_sequence[i][j];
	  }

	  for (j = 0;j < length[i];j++) {
	    norm[j] = 0.;

	    for (k = 0;k < hsw_markov->nb_state;k++) {

	      // calcul des probabilites d'observation

	      observation[j][k] = 1.;
	      observation[j][k] *= hsw_markov->sw_process[1]->observation[k]->
		density_computation(*poutput[0], covar[i][j], 0, year_effect_predict[r][index[i][j]-1]);


	      switch (hsw_markov->state_subtype[k]) {

		// cas etat semi-markovien

	      case SEMI_MARKOVIAN : {
		if (j == 0) {
		  state_norm[k] = hsw_markov->initial[k];
		}
		else {
		  state_norm[k] += state_in[j - 1][k] - forward[j - 1][k];
		}

		state_norm[k] *= observation[j][k];
		norm[j] += state_norm[k];
		break;
	      }

		// cas etat markovien

	      case MARKOVIAN : {
		if (j == 0) {
		  forward[j][k] = hsw_markov->initial[k];
		}
		else {
		  forward[j][k] = state_in[j - 1][k];
		}

		forward[j][k] *= observation[j][k];

		norm[j] += forward[j][k];

		break;
	      }
	      }
	    }
  
	    if (norm[j] > 0.) {

	      for (k = 0;k < hsw_markov->nb_state;k++) {
		switch (hsw_markov->state_subtype[k]) {
		case SEMI_MARKOVIAN :
		  state_norm[k] /= norm[j];
		  break;
		case MARKOVIAN :
		  forward[j][k] /= norm[j];
		  break;
		}
	      }

	      likelihoodnum += log(norm[j]);
	    }

	    else {
	      likelihoodnum = D_INF;
	      break;
	    }

	    for (k = 0;k < hsw_markov->nb_state;k++) {

	      // cas etat semi-markovien

	      if (hsw_markov->state_subtype[k] == SEMI_MARKOVIAN) {
              
		occupancy = hsw_markov->nonparametric_process[0]->sojourn_time[k];
		obs_product = 1.;
		forward[j][k] = 0.;

		if (j < length[i] - 1) {
		  for (m = 1;m <= MIN(j + 1 , occupancy->nb_value   - 1);m++) {
		    obs_product *= observation[j - m + 1][k] / norm[j - m + 1];
		    if (obs_product == 0.) {
		      break;
		    }
		    if (m < j + 1) {
		      forward[j][k] += obs_product * occupancy->mass[m] * state_in[j - m][k];
		    }

		    else {
		      forward[j][k] += obs_product * occupancy->mass[m] * hsw_markov->initial[k];
		    }
		  }
		}

		else {
		  for (m = 1;m <= MIN(j + 1 , occupancy->nb_value - 1);m++) {
		    obs_product *= observation[j - m + 1][k] / norm[j - m + 1];
		    if (obs_product == 0.) {
		      break;
		    }

		    if (m < j + 1) {
		      forward[j][k] += obs_product * (1. - occupancy->cumul[m - 1]) * state_in[j - m][k];
		    }

		    else {
		      forward[j][k] += obs_product * (1. - occupancy->cumul[m - 1]) *
			hsw_markov->initial[k];
		    }
		  }
		}
	      }
	    }

	    if (j < length[i] - 1) {
	      for (k = 0;k < hsw_markov->nb_state;k++) {
		state_in[j][k] = 0.;
		for (m = 0;m < hsw_markov->nb_state;m++) {
		  state_in[j][k] += hsw_markov->transition[m][k] * forward[j][m];
              
		}
	      }
	    }

	    for (k = 0;k < nb_variable;k++) {
	      poutput[k]++;
	    }
	  }

	  if (likelihood == D_INF) {
	    break;
	  }

	}

	for (j = 0; j < T; j++){
	  likelihoodnum = likelihoodnum + log(hsw_markov->sw_process[1]->observation[0]->dnorm(year_effect_predict[r][j],0.,1.));
	}

	for (j = 0; j < T; j++){
	  likelihoodden = likelihoodden + log(hsw_markov->sw_process[1]->observation[0]->dnorm(year_effect_predict[r-1][j],0.,1.));
	}


	prob = fmin(0, likelihoodnum - likelihoodden);

	if(U > exp(prob)){
	  for(j = 0; j < T; j++){
	    year_effect_predict[r][j] = year_effect_predict[r-1][j];
	  }


	  for (j = 0; j < T; j++){
	    likelihoodden = likelihoodden - log(hsw_markov->sw_process[1]->observation[0]->dnorm(year_effect_predict[r-1][j],0.,1.));
	  }
	  
	  
	  likelihood = likelihoodden;
	}
	else{
	  for (j = 0; j < T; j++){
	    likelihoodnum = likelihoodnum - log(hsw_markov->sw_process[1]->observation[0]->dnorm(year_effect_predict[r][j],0.,1.));
	  }


	  likelihood = likelihoodnum;
	}

      }



      for(i = 0; i < T; i++){
	year_effect[i] = 0.;
      }

      for ( r = 1; r <= nb_random_simul; r++){
	for(i = 0; i < T; i++){
	  year_effect[i] = year_effect[i] + year_effect_predict[r][i];
	}
      }
      
      for(i = 0; i < T; i++){
	year_effect[i] = year_effect[i] / (double) nb_random_simul;
      }



      // estimation des parametres des modèles linéaires pour chacun des états 

      // Cas avec Variance différente pour chaque état	
      if (!VarCommune){
	for (k = 0; k < hsw_markov->nb_state; k++) {
	  
	  double *beta_k, *tau_k;
	  double sigma_k = 0.;
	  
	  beta_k = new double[nb_covariable];
	  tau_k = new double[nb_random];
	  
	  // estimation des parametres de regression pour chacun des états
	  
	  beta_k = regression_parameter_temporal_stochastic(hsw_markov, nb_sequence, length, real_sequence, nb_covariable, covar, 
							    year_effect, backward, index, nb_random, k); 

	  // Calcul des variances résiduelles associées à chaque état

	  sigma_k = residual_variance_temporal_stochastic(hsw_markov, nb_sequence, length, real_sequence, nb_covariable, covar, 
							  year_effect, backward, index, beta_k, nb_random, k); 

	  // Calcul des variances associées aux effets aléatoires 
	  
	  tau_k = random_variance_temporal_stochastic(hsw_markov, nb_sequence, length, real_sequence, nb_covariable, covar, 
						      year_effect, backward, index, beta_k, nb_random, k);

#   ifdef DEBUG
	  cout<<"beta pour l'état "<<k<<" :  "<<endl;
	  for(i = 0; i< nb_covariable; i++){
	    cout<<beta_k[i]<<"   ";
	  }
	  cout<<endl;
#endif

#   ifdef DEBUG
	  cout<<"variance pour l'état "<<k<<" :  "<<sigma_k<<endl;
#   endif

	  // mise à jour des lois d'observations
	  
	  hsw_markov->sw_process[1]->observation[k]->Set_regression(beta_k);
	  hsw_markov->sw_process[1]->observation[k]->Set_residual_variance(sigma_k);
	  hsw_markov->sw_process[1]->observation[k]->Set_random_variance(tau_k);
	  
	  delete [] tau_k;
	  delete [] beta_k;
	}
      }

      else{
	// Cas avec Variance commune pour chaque état	
	
	double **beta;
	double *tau_k;
 	double sigma = 0.;
	
 	beta = new double*[hsw_markov->nb_state];
	tau_k = new double[nb_random];	

	// estimation des parametres de regression pour chacun des états

 	for (k = 0; k < hsw_markov->nb_state; k++){
	  beta[k] = new double[nb_covariable];
	  beta[k] = regression_parameter_temporal_stochastic(hsw_markov, nb_sequence, length, real_sequence, nb_covariable, covar, 
							     year_effect, backward, index, nb_random, k); 
	  
#   ifdef DEBUG
	  cout<<"beta pour l'état "<<k<<" :  "<<endl;
	  for(i = 0; i< nb_covariable; i++){
	    cout<<beta_k[i]<<"   ";
	  }
	  cout<<endl;
#endif
	}

	// Calcul de la variance résiduelle commune à tous les états

	sigma = residual_variance_temporal_VarCommune_stochastic(hsw_markov, nb_sequence, length, real_sequence, nb_covariable, 
								 covar, year_effect, backward, index, beta, nb_random); 


	// Calcul de la variance aléatoire commune à tous les états
	
	tau_k = random_variance_temporal_VarCommune_stochastic(hsw_markov, nb_sequence, length, real_sequence, nb_covariable, covar, 
							       year_effect, backward, index, beta, nb_random);


#   ifdef DEBUG
	cout<<"variance résiduelle  :  "<< sigma << endl;
#   endif
	
	// mise à jour des lois d'observations
	
	for (k = 0; k < hsw_markov->nb_state; k++){
	  hsw_markov->sw_process[1]->observation[k]->Set_regression(beta[k]);
	  hsw_markov->sw_process[1]->observation[k]->Set_residual_variance(sigma);
	  hsw_markov->sw_process[1]->observation[k]->Set_random_variance(tau_k);
	}

	delete [] tau_k;
	delete [] beta;
      }
 
      if (likelihood != D_INF) {
	
	// reestimation des probabilites initiales
	
	if (hsw_markov->type == 'o') {
	  reestimation(hsw_markov->nb_state, chain_reestim->initial, hsw_markov->initial , MIN_PROBABILITY , false);
	}
	
	// reestimation des probabilites de transition
	
	for (i = 0;i < hsw_markov->nb_state; i++) {
	  reestimation(hsw_markov->nb_state, chain_reestim->transition[i], hsw_markov->transition[i] , MIN_PROBABILITY , false);
	}

	// reestimation des lois d'occupation des etats

        min_likelihood = 0.;

        for (i = 0;i < hsw_markov->nb_state;i++) {
          if (hsw_markov->state_subtype[i] == SEMI_MARKOVIAN) {
            occupancy = hsw_markov->nonparametric_process[0]->sojourn_time[i];

#           ifdef DEBUG
            cout << STAT_label[STATL_STATE] << " " << i << " (";
#           endif

            if ((hsw_markov->type == 'o') || (estimator == PARTIAL_LIKELIHOOD)) {
              occupancy_reestim[i]->nb_value_computation();
              occupancy_reestim[i]->offset_computation();
              occupancy_reestim[i]->nb_element_computation();
              occupancy_reestim[i]->max_computation();
              occupancy_reestim[i]->mean_computation();
              occupancy_reestim[i]->variance_computation();

#             ifdef DEBUG
              if (hsw_markov->type == 'o') {
                switch (estimator) {
                case COMPLETE_LIKELIHOOD :
                  cout << test[i][0] + test[i][2] << " | " << test[i][1] + test[i][3];
                  break;
                case PARTIAL_LIKELIHOOD :
                  cout << test[i][0];
                  break;
                }
                cout << " | " << occupancy_reestim[i]->nb_element << ") - ";
                occupancy_reestim[i]->ascii_characteristic_print(cout);
              }
#             endif

            }

            else {
              offset = 1;
              nb_value = occupancy_nb_value[i];

              ofrequency = occupancy_reestim[i]->frequency + occupancy_nb_value[i];
              lfrequency = length_bias_reestim[i]->frequency + occupancy_nb_value[i];
              while ((*--ofrequency == 0) && (*--lfrequency == 0) && (nb_value > 2)) {
                nb_value--;
              }
              occupancy_reestim[i]->nb_value = nb_value;
              length_bias_reestim[i]->nb_value = nb_value;

              ofrequency = occupancy_reestim[i]->frequency + offset;
              lfrequency = length_bias_reestim[i]->frequency + offset;
              while ((*ofrequency++ == 0) && (*lfrequency++ == 0) && (offset < nb_value - 1)) {
                offset++;
              }
              occupancy_reestim[i]->offset = offset;
              length_bias_reestim[i]->offset = offset;

              occupancy_reestim[i]->nb_element_computation();
              length_bias_reestim[i]->nb_element_computation();

#             ifdef DEBUG
              occupancy_reestim[i]->max_computation();
              occupancy_reestim[i]->mean_computation();
              occupancy_reestim[i]->variance_computation();

              cout << test[i][1] << " | " << occupancy_reestim[i]->nb_element << ") - ";
              occupancy_reestim[i]->ascii_characteristic_print(cout);

              length_bias_reestim[i]->max_computation();
              length_bias_reestim[i]->mean_computation();
              length_bias_reestim[i]->variance_computation();

              cout << STAT_label[STATL_STATE] << " " << i << " (" << test[i][3] << " | "
                   << length_bias_reestim[i]->nb_element << ") - ";
              length_bias_reestim[i]->ascii_characteristic_print(cout);
#             endif

              switch (mean_computation) {
              case COMPUTED :
                occupancy_mean = interval_bisection(occupancy_reestim[i] , length_bias_reestim[i]);
                break;
              case ONE_STEP_LATE :
                occupancy_mean = occupancy->mean;
                break;
              }

#             ifdef DEBUG
              cout << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_MEAN] << ": "
                   << occupancy_mean << endl;
#             endif

              occupancy_reestim[i]->equilibrium_process_combination(length_bias_reestim[i] , occupancy_mean);

#             ifdef DEBUG
              cout << test[i][0] + test[i][2] << " | " << test[i][1] + test[i][3] << " | "
                   << occupancy_reestim[i]->nb_element << ") - ";
              occupancy_reestim[i]->ascii_characteristic_print(cout);
#             endif
            }

            hoccupancy->update(occupancy_reestim[i] ,
                               MAX((int)(occupancy_reestim[i]->nb_element *
                                         MAX(sqrt(occupancy_reestim[i]->variance) , 1.) * OCCUPANCY_COEFF) , MIN_NB_ELEMENT));

            if (iter <= EXPLORATION_NB_ITER) {
              occupancy_likelihood = hoccupancy->Reestimation<int>::parametric_estimation(occupancy, 1, true, OCCUPANCY_THRESHOLD);
            }
            else {
              occupancy_likelihood = hoccupancy->Reestimation<int>::type_parametric_estimation(occupancy, 1, true, OCCUPANCY_THRESHOLD);
            }

            if (occupancy_likelihood == D_INF) {
              min_likelihood = D_INF;
            }
            else {
              occupancy->computation(hoccupancy->nb_value , OCCUPANCY_THRESHOLD);
            }

#           ifdef DEBUG
            cout << STAT_word[STATW_STATE] << " " << i << " " << STAT_word[STATW_OCCUPANCY_DISTRIBUTION] << endl;
            occupancy->ascii_print(cout);
#           endif

          }
        }
      }
      
      likelihood = hsw_markov->likelihood_computation(*this, 0,I_DEFAULT);
      
      //#     ifdef MESSAGE
      os << STAT_label[STATL_ITERATION] << " " << iter << "   "
	 << SEQ_label[SEQL_OBSERVED_SEQUENCES_LIKELIHOOD] << ": " << likelihood << " ( "<< nb_random_simul << " )" << endl;
      //#     endif
      
    
      for (i = 0;i < nb_sequence;i++) {
	for (j = 0; j < max_length; j++) {
	  delete [] backward[i][j];
	}
	delete [] backward[i];
      }
      delete [] backward;
    
      if(likelihood - previous_likelihood <= 0.01){tmp = tmp+1;}
      
    }
    while  (iter < nb_iter);
    //((likelihood != D_INF) && ((likelihood - previous_likelihood) / -likelihood > MARKOV_SWITCHING_LIKELIHOOD_DIFF) &&
    //	   ((nb_iter == I_DEFAULT) && (iter < MARKOV_SWITCHING_NB_ITER) || ((nb_iter != I_DEFAULT) && (iter < nb_iter)))
    //	   || (tmp < nb_iter_conv));
    

  
    if (likelihood != D_INF) {
      
#     ifdef MESSAGE
      os << "\n" << iter << " " << STAT_label[STATL_ITERATIONS] << endl;
#     endif
      
      // reestimation des probabilites initiales
      
      if (hsw_markov->type == 'o') {
	reestimation(hsw_markov->nb_state , chain_reestim->initial ,
		     hsw_markov->initial , MIN_PROBABILITY , true);
      }
      
      // reestimation des probabilites de transition
      
      for (i = 0;i < hsw_markov->nb_state;i++) {
	reestimation(hsw_markov->nb_state , chain_reestim->transition[i] ,
		     hsw_markov->transition[i] , MIN_PROBABILITY , true);
      }
      
      for (i = 0;i < hsw_markov->nb_state;i++) {
        if ((hsw_markov->state_subtype[i] == SEMI_MARKOVIAN) &&
            (hsw_markov->nonparametric_process[0]->sojourn_time[i]->mean == 1.)) {
          hsw_markov->state_subtype[i] = MARKOVIAN;
          delete hsw_markov->nonparametric_process[0]->sojourn_time[i];
          hsw_markov->nonparametric_process[0]->sojourn_time[i] = 0;
          delete hsw_markov->forward[i];
          hsw_markov->forward[i] = 0;
        }
      }
   
      // reestimation des lois d'observation
      
      for (j = 0;j < hsw_markov->nb_state;j++) {
	hsw_markov->sw_process[1]->observation[j]->parametric_variance_computation();
	hsw_markov->sw_process[1]->observation[j]->min_value_computation();
	hsw_markov->sw_process[1]->observation[j]->max_value_computation();
	hsw_markov->sw_process[1]->observation[j]->nb_value_computation();
      }
      
    }
    

#   ifdef DEBUG 
    cout<<"AFFICHAGE DES VARIANCES RESIDUELLES:  "<<endl;
    for (j = 0; j< hsw_markov->nb_state; j++) {
      cout<< hsw_markov->sw_process[1]->observation[j]->residual_variance<<"   ";
    }
    cout<<endl;

    cout<<"AFFICHAGE DES PARAMETRES DE REGRESSION: "<<endl;
    for (j = 0; j< hsw_markov->nb_state; j++) {
      for(k=0; k<nb_covariable; k++) {
	cout << hsw_markov->sw_process[1]->observation[j]->regression[k]<< "   ";
      }
      cout<<endl;
    }
    cout<<endl ;
#   endif

    // destruction des structures de donnees de l'algorithme

    for (i = 0;i < max_length;i++) {
      delete [] observation[i];
    }
    delete [] observation;

    delete [] norm;
    delete [] state_norm;

    for (i = 0;i < max_length;i++) {
      delete [] forward[i];
    }
    delete [] forward;

    for (i = 0;i < max_length - 1;i++) {
      delete [] state_in[i];
    }
    delete [] state_in;

    delete [] back;

    for (i = 0;i < max_length;i++) {
      delete [] backward1[i];
    }
    delete [] backward1;

    delete [] auxiliary;

    delete chain_reestim;

    delete [] occupancy_nb_value;

    for (i = 0;i < hsw_markov->nb_state;i++) {
      delete occupancy_reestim[i];
    }
    delete [] occupancy_reestim;

    if (estimator == KAPLAN_MEIER) {
      delete [] censored_occupancy_nb_value;

      for (i = 0;i < hsw_markov->nb_state;i++) {
        delete censored_occupancy_reestim[i];
      }
      delete [] censored_occupancy_reestim;

      delete [] occupancy_survivor;
      delete [] censored_occupancy_survivor;
    }

    delete hoccupancy;


    
    for (i = 0;i < max_length;i++) {
      delete [] predicted[i];
    }
    delete [] predicted;
    
    delete [] temp;

    delete [] poutput;

    
    if (likelihood == D_INF) {
      delete hsw_markov;
      hsw_markov = 0;
      error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
    }
    
    else {
      if (state_sequence) {

	hsw_markov->sw_markov_data = new Markov_switching_data(*this ,0);

	hsw_markov_data = hsw_markov->sw_markov_data;
	hsw_markov_data->type[0] = STATE;
	
	if ((hsw_markov->sw_process[1]) && (hsw_markov_data->characteristics[1])) {
	  delete hsw_markov_data->characteristics[1];
	  hsw_markov_data->characteristics[1] = 0;
	}

	hsw_markov->create_cumul();	
	
#   ifdef DEBUG
	for (i = 0;i < hsw_markov->nb_state;i++) {
	  for (j = 0; j < hsw_markov->nb_state; j++){
	    cout<<hsw_markov->initial[i]<<"   ";
	  }
	  cout<<endl;
	}
#   endif

	hsw_markov_data->posterior_probability = new double[hsw_markov_data->nb_sequence];
	
	if (!output_R){
	  hsw_markov_data->likelihood = hsw_markov->viterbi(*hsw_markov_data, hsw_markov_data->posterior_probability, 
							    true, false);
	}
	else{
	  hsw_markov_data->likelihood = hsw_markov->viterbi(*hsw_markov_data, hsw_markov_data->posterior_probability, 
							    false, true);
	}

	hsw_markov->remove_cumul();
        hsw_markov_data->max_value[0] = hsw_markov->nb_state - 1;
	hsw_markov_data->build_marginal_histogram(0);
	hsw_markov_data->build_characteristic(0);
	hsw_markov_data->build_transition_count(*hsw_markov);

        // calcul des lois d'occupation des etats

        for (i = 0;i < hsw_markov->nb_state;i++) {
          if (hsw_markov->state_subtype[i] == SEMI_MARKOVIAN) {
            hsw_markov->nonparametric_process[0]->sojourn_time[i]->computation((hsw_markov_data->characteristics[0] ? 
										hsw_markov_data->characteristics[0]->
										sojourn_time[i]->nb_value : 1) , OCCUPANCY_THRESHOLD);
            if (hsw_markov->state_type[i] == 'r') {
	      hsw_markov->forward[i]->copy(*(hsw_markov->nonparametric_process[0]->sojourn_time[i]));
	      hsw_markov->forward[i]->computation(*(hsw_markov->nonparametric_process[0]->sojourn_time[i]));
            }
          }
        }

#       ifdef MESSAGE
	cout << "\n" << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << ": " << hsw_markov_data->likelihood      <<endl;
#       endif

      }
      
      else {
	for (i = 0;i < hsw_markov->nb_state;i++) {
	  if ((hsw_markov->state_subtype[i] == SEMI_MARKOVIAN) && (hsw_markov->state_type[i] == 'r')) {
	    hsw_markov->forward[i]->copy(*(hsw_markov->nonparametric_process[0]->sojourn_time[i]));
	    hsw_markov->forward[i]->computation(*(hsw_markov->nonparametric_process[0]->sojourn_time[i]));
	  }
	}
     
	hsw_markov->sw_markov_data = new Markov_switching_data(*this , false);
	hsw_markov_data = hsw_markov->sw_markov_data;
	hsw_markov_data->state_variable_init(REAL_VALUE); 

	if ((hsw_markov->sw_process[1]) && (hsw_markov_data->characteristics[0])) {
	  delete hsw_markov_data->characteristics[0];
	  hsw_markov_data->characteristics[0] = 0;
	}
      }
 
      // calcul de la vraisemblance et des lois caracteristiques du modele

      hsw_markov_data->hidden_likelihood = hsw_markov->likelihood_computation(*this, hsw_markov_data->posterior_probability, I_DEFAULT);
      
#     ifdef DEBUG
      cout << *hsw_markov;
      cout << "iteration " << iter << "  " << endl
	   << SEQ_label[SEQL_OBSERVED_SEQUENCES_LIKELIHOOD] << ": " << hsw_markov_data->hidden_likelihood << endl;
#     endif
      
#     ifdef MESSAGE
      if  ((state_sequence) && (seq->nb_sequence <= POSTERIOR_PROBABILITY_NB_SEQUENCE)) {
        int *pstate;

        os << "\n" << SEQ_label[SEQL_POSTERIOR_STATE_SEQUENCE_PROBABILITY] << endl;
        for (i = 0;i < seq->nb_sequence;i++) {
          os << SEQ_label[SEQL_SEQUENCE] << " " << seq->identifier[i] << ": "
             << seq->posterior_probability[i];

          if (hsw_markov->nb_component == hsw_markov->nb_state) {
            os << " | " << SEQ_label[SEQL_STATE_BEGIN] << ": ";

            pstate = seq->int_sequence[i][0] + 1;
            if (seq->index_parameter) {
              for (j = 1;j < seq->length[i];j++) {
                if (*pstate != *(pstate - 1)) {
                  os << seq->index_parameter[i][j] << ", ";
                }
                pstate++;
              }
            }

            else {
              for (j = 1;j < seq->length[i];j++) {
                if (*pstate != *(pstate - 1)) {
                  os << j << ", ";
                }
                pstate++;
              }
            }
          }

          os << endl;
        }

	delete [] pstate;
      }
#     endif
 
      hsw_markov->component_computation();
      hsw_markov->characteristic_computation(*hsw_markov_data , counting_flag , I_DEFAULT , false);
    }
    
    cout<<"Realisation de l'effet aléatoire environnement commun"<<endl;
    cout<<"c( ";
    for ( i = 0; i < (T-1); i++){
      cout<<year_effect[i]<<" ,";
    }
    cout<<year_effect[T-1]<<")"<<endl<<endl;

  }
  
  return hsw_markov;
  
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'un semi-Markov switching linear 
 *  mixed model a partir d'un echantillon de sequences par 
 *  l'algorithme SEM/MCEM.
 *
 *  State sequences : "forward-backward" 
 *  Random effects : Sampling by Metropolis-Hastings algorithm 
 *
 *  arguments : reference sur un objet Format_error, stream, type 
 *              de processus ('o' : ordinaire, 'e' : en equilibre), 
 *              nombre d'etats de la chaine de Markov, flag sur la 
 *              nature de la chaine de Markov, variance résiduelle, 
 *              variance aléatoire, paramètres de regression, nombre 
 *              de covariables, nombre d'effets aléatoires, presence/
 *              absence de constante, ordre, nombre d'iterations pour 
 *              convergence, parametres pour nombre de sequences d'etats 
 *              simulees, type Markovien, variance commune à tous les 
 *              états, sortie de type R, flags sur le calcul des lois 
 *              de comptage et sur le calcul des sequences d'etats 
 *              optimales, probabilite de rester dans un etat initiale, 
 *              nombre d'iterations.
 *
 *--------------------------------------------------------------*/

Markov_switching* Switching_sequence::markov_switching_stochastic_estimation_MH_year(Format_error &error, ostream &os,
										     char type, int nb_state, bool left_right,
										     double *residual_variance, double **random_variance,
										     double **regression, int nb_covariable,
										     int nb_random, int constant, int order,
										     int nb_iter_conv, int min_nb_random,
										     int max_nb_random, double VarMH,
										     int estimator, bool markov,
										     bool VarCommune, bool output_R, double parameter, 
										     bool counting_flag, bool state_sequence, 
										     double occupancy_mean, double self_transition, 
										     int nb_iter, int mean_computation) const
{
  bool status = true;
  register int i, j;
  int nb_value[SEQUENCE_NB_VARIABLE];
  Chain *pchain;
  Continuous_parametric **cont_param;
  Switching_process *ihsw_process;
  Markov_switching *ihsw_markov , *hsw_markov;
  double proba;
  Parametric **occupancy;

  hsw_markov = 0;
  error.init();

  if ((nb_state < 2) || (nb_state > NB_STATE)) {
    status = false;
    error.update(SEQ_error[SEQR_NB_STATE]);
  }

  if ((occupancy_mean != D_DEFAULT) && (occupancy_mean <= 1.)) {
    status = false;
    error.update(SEQ_error[SEQR_OCCUPANCY]);
  }

  if (status) {
    for (i = 0;i < nb_variable;i++) {
      nb_value[i] = marginal[i]->nb_value;
    }

    if (!constant){
      nb_covariable = nb_covariable + 1;
    }

    pchain = new Chain('o', nb_state, nb_state, true);
    if (self_transition == D_DEFAULT) {
      if (hlength->mean != 0) {
      self_transition = MAX(1. - 1. / hlength->mean , SELF_TRANSITION);
      }
    }

    if (!markov){
      pchain->init(left_right,self_transition);
    }
    else {
      for(i = 0; i < nb_state; i++){
	pchain->transition[i][i] = 0.7;
	for (j = 0; j < nb_state; j++){
	  pchain->transition[i][j] = 0.3/((double)nb_state - 1);
	}
      }
    }

    for (i = 0; i < nb_state; i++){
      pchain->initial[i] = 1./((double)nb_state);
    }

    cont_param = new Continuous_parametric*[nb_state];
    for (i = 0; i < nb_state; i++) {
      cont_param[i] = new Continuous_parametric(NB_VALUE, nb_covariable, nb_random, residual_variance[i], random_variance[i], 
						regression[i], HETERO_RANDOM);
    }
    ihsw_process = new Switching_process(nb_state, cont_param);

 // initialisation des lois d'occupations des etats

    occupancy = new Parametric*[nb_state];
    if (occupancy_mean == D_DEFAULT) {
      occupancy_mean = MAX(hlength->mean , OCCUPANCY_MEAN);
    }

    proba = 1./occupancy_mean;

    if (!markov){
      for(i = 0; i < (nb_state-1); i++){
	occupancy[i] = new Parametric(NEGATIVE_BINOMIAL, 1, I_DEFAULT, 1., proba, OCCUPANCY_THRESHOLD);
      }
      occupancy[nb_state-1] = 0;
    }
    else {
      for (i = 0; i < nb_state; i++){
	occupancy[i] = 0;
      }
    }

    Nonparametric_sequence_process *poccupancy;
    poccupancy = new Nonparametric_sequence_process(nb_state, occupancy);
    ihsw_markov = new Markov_switching(pchain, poccupancy, ihsw_process, 20);

    // initialisation des lois d'observation

    ihsw_markov->sw_process[1]->init();
    hsw_markov = markov_switching_stochastic_estimation_MH_year(error, os, *ihsw_markov, nb_iter_conv, min_nb_random, max_nb_random, 
								VarMH, estimator, VarCommune, output_R, parameter, counting_flag, 
								state_sequence, nb_iter, mean_computation);
    delete ihsw_markov;
  }

  return hsw_markov;
}

