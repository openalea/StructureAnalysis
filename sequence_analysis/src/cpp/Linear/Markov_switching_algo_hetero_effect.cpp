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
 * (modèle avec effets aléatoires individuels)
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
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
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
 *              effets aléatoires individuels, flag sur le calcul 
 *              des lois de comptage, flag calcul d'une divergence de Kullback-Leibler.
 *
 *--------------------------------------------------------------*/

Markov_switching_data* Markov_switching::simulation_hetero_effect(Format_error &error , const Histogram &hlength ,
								  int inb_covariable, int inb_random, 
								  int iconstant, double *iresidual_variance, 
								  double **iregression, double **irandom_variance, 
								  double ***icovar, double **ieffect,
								  bool counting_flag , bool divergence_flag) const
{
  bool status = true;
  register int i , j , k, m;
  int cumul_length , occupancy, *pstate ;
  double **poutput;
  double *pcovar, *ccovar, *peffect, *ceffect;
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

    seq = new Markov_switching_data(2 , hlength , inb_covariable, inb_random, nb_state, 0, iconstant, false);

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
	for (k = 0;k < nb_covar;k++) {
	  *pcovar++ = *ccovar++;
	}
      }
    }

    for (i = 0;i < seq->nb_sequence;i++) {
      peffect = seq->effect[i];
      ceffect = ieffect[i];
      for (j = 0;j < sw_markov->nb_state;j++) {
	*peffect++ = *ceffect++;
      }
    }

    for (i = 0;i < seq->nb_sequence;i++) {
      pstate = seq->int_sequence[i][0];
      *pstate = cumul_method(sw_markov->nb_state , sw_markov->cumul_initial);

      poutput[j] = seq->real_sequence[i][1];
      
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
	      + sqrt(irandom_variance[*pstate][0])*seq->effect[i][*pstate];
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

      // seq->likelihood = sw_markov->likelihood_computation(*seq); 

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
 *              effets aléatoires individuels, flag sur le calcul 
 *              des lois de comptage.
 *
 *--------------------------------------------------------------*/

Markov_switching_data* Markov_switching::simulation_hetero_effect(Format_error &error , int nb_sequence ,
								  int length , int inb_covariable, int inb_random,
								  int iconstant,  double *iresidual_variance,
								  double **iregression, double **irandom_variance, 
								  double ***icovar, double **ieffect, bool counting_flag) const
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

    seq = simulation_hetero_effect(error , hlength , inb_covariable, inb_random, iconstant, iresidual_variance, iregression, 
				   irandom_variance,  icovar, ieffect, counting_flag, false);
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
 *              effets aléatoires individuels, flag sur le calcul 
 *              des lois de comptage.
 *
 *--------------------------------------------------------------*/

Markov_switching_data* Markov_switching::simulation_hetero_effect(Format_error &error , int nb_sequence ,
								  const Switching_sequence &iseq ,
								  int inb_covariable, int inb_random, int iconstant,
								  double *iresidual_variance, double **iregression,
								  double **irandom_variance, double ***icovar,
								  double **ieffect, bool counting_flag) const
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

    seq = simulation_hetero_effect(error , *hlength , inb_covariable, inb_random, iconstant, iresidual_variance,iregression, 
				   irandom_variance, icovar, ieffect, counting_flag, false);
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
 *              effets aléatoires individuels, flag sur le calcul 
 *              des lois de comptage, flag calcul d'une divergence 
 *              de Kullback-Leibler.
 *
 *--------------------------------------------------------------*/

Markov_switching_data* Markov_switching::simulation_hetero_effect(Format_error &error ,
								  const Histogram &hlength ,
								  int inb_covariable, int inb_random, int iconstant, 
								  double *iresidual_variance, double **iregression,
								  double **irandom_variance, double ***icovar,
								  double **ieffect, bool counting_flag ,
								  bool divergence_flag, bool hidden) const
{
  Switching_sequence *observ_seq;
  Markov_switching_data *seq;


  seq = Markov_switching::simulation_hetero_effect(error , hlength , inb_covariable, inb_random, iconstant,  iresidual_variance, 
				     iregression, irandom_variance, icovar, ieffect, counting_flag , divergence_flag);

  if ((seq) && (!divergence_flag)) {
    observ_seq = seq->remove_variable_1();
    seq->hidden_likelihood = likelihood_computation(*observ_seq,0, I_DEFAULT);
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
 *              effets aléatoires individuels, flag sur le calcul 
 *              des lois de comptage.
 *
 *--------------------------------------------------------------*/

Markov_switching_data* Markov_switching::simulation_hetero_effect(Format_error &error , int nb_sequence , int length ,
								  int inb_covariable, int inb_random, int iconstant, 
								  double *iresidual_variance, double **iregression,
								  double **irandom_variance, double ***icovar,
								  double **ieffect,  bool counting_flag, bool hidden) const
{
  Switching_sequence *observ_seq;
  Markov_switching_data *seq;


  seq = Markov_switching::simulation_hetero_effect(error , nb_sequence , length , inb_covariable, inb_random, iconstant,
						   iresidual_variance, iregression, irandom_variance, icovar, ieffect, counting_flag);

  if (seq) {
    observ_seq = seq->remove_variable_1();
    seq->hidden_likelihood = likelihood_computation(*observ_seq,0, I_DEFAULT);
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
 *              effets aléatoires individuels, flag sur le calcul 
 *              des lois de comptage.
 *
 *--------------------------------------------------------------*/

Markov_switching_data* Markov_switching::simulation_hetero_effect(Format_error &error ,
								  int nb_sequence ,
								  const Switching_sequence &iseq ,
								  int inb_covariable, int inb_random, int iconstant, 
								  double *iresidual_variance, double **iregression,
								  double **irandom_variance, double ***icovar,
								  double **ieffect, 
								  bool counting_flag, bool hidden) const

{
  Switching_sequence *observ_seq;
  Markov_switching_data *seq;
  
  seq = Markov_switching::simulation_hetero_effect(error , nb_sequence , iseq , inb_covariable, inb_random,iconstant,
						   iresidual_variance, iregression, irandom_variance, icovar, ieffect, counting_flag);
  
  if (seq) {
    observ_seq = seq->remove_variable_1();
    seq->hidden_likelihood = likelihood_computation(*observ_seq,0, I_DEFAULT);
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
 *  State sequences : "forward-backward" sampling
 *  Random effects : Prediction by conditionnal expectation 
 *
 *
 *  arguments : reference sur un objet Format_error, stream, 
 *              semi-Markov switching linear mixed model initial,
 *              type d'effet aléatoire, nombre d'itérations pour convergence
 *              parametres pour le nombre de sequences d'etats simulees, 
 *              variance résiduelle commune à tous les états, sortie de 
 *              type R, type d'estimateur pour la reestimation 
 *              des lois d'occupation des etats, flags sur le calcul
 *              des lois de comptage et sur le calcul des sequences 
 *              d'etats optimales, nombre d'iterations.
 *
 *--------------------------------------------------------------*/

Markov_switching* Switching_sequence::markov_switching_stochastic_estimation_hetero(Format_error &error , ostream &os ,
										    const Markov_switching &ihsw_markov ,
										    int type_random, int nb_iter_conv,
										    int min_nb_state_sequence, int max_nb_state_sequence, 
										    bool VarCommune, bool output_R, double parameter , 
										    int estimator, bool counting_flag ,
										    bool state_sequence , int nb_iter) const
{
  bool status;
  register int i, j, k, l, m, n, r, kk, tmp;
  int max_nb_value , iter , nb_state_sequence , state_occupancy , nb_likelihood_decrease ,
    *occupancy_nb_value , *state_seq , *pstate ;
  double **poutput;
  double likelihood = D_INF , previous_likelihood , occupancy_likelihood , observation_likelihood ,
    min_likelihood , obs_product , **observation , *norm , *state_norm , **forward ,
    **state_in , *backward , *cumul_backward , *reestim , *occupancy_survivor ,
    *censored_occupancy_survivor;
  double ****indic;
  double ***random_predict;

  Chain_reestimation<double> *chain_reestim;
  Parametric *occupancy;
  Reestimation<double> *bcomplete_run , *censored_run , **complete_run , **final_run , **initial_run ,
    **single_run;
  Markov_switching *hsw_markov;
  Markov_switching_data *hsw_markov_data;
  const Reestimation<double> *prun[3];


# ifdef DEBUG
  double sum;
# endif

  hsw_markov = 0;
  error.init();

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

  if ((min_nb_state_sequence < 1) || (min_nb_state_sequence > max_nb_state_sequence)) {
    status = false;
    error.update(SEQ_error[SEQR_MIN_NB_STATE_SEQUENCE]);
  }

  tmp = 0 ;


  if (status) {
    if (max_length > COUNTING_MAX_LENGTH) {
      counting_flag = false;
    }

    // creation de la chaine de Markov cachee

    hsw_markov = new Markov_switching(ihsw_markov , false, (int)(max_length*SAMPLE_NB_VALUE_COEFF));

    Parametric **sojourn;
    sojourn = new Parametric*[hsw_markov->nb_state];


#   ifdef DEBUG
    cout << *hsw_markov;
#   endif

    // initialisations
    
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
    
    backward = new double[max_length + 1];
    cumul_backward = new double[max_length + 1];
    
    state_seq = new int[max_length];

    chain_reestim = new Chain_reestimation<double>('o', hsw_markov->nb_state , hsw_markov->nb_state);
    
    occupancy_nb_value = new int[hsw_markov->nb_state];
    complete_run = new Reestimation<double>*[hsw_markov->nb_state];
    final_run = new Reestimation<double>*[hsw_markov->nb_state];
    
    for (i = 0;i < hsw_markov->nb_state;i++) {
      switch (hsw_markov->state_subtype[i]) {
	
      case SEMI_MARKOVIAN : {
        occupancy_nb_value[i] = MIN(hsw_markov->nonparametric_process[0]->sojourn_time[i]->alloc_nb_value , max_length + 1);
	complete_run[i] = new Reestimation<double>(occupancy_nb_value[i]);
        final_run[i] = new Reestimation<double>(occupancy_nb_value[i]);
        break;
      }
	
      case MARKOVIAN : {
        complete_run[i] = 0;
        final_run[i] = 0;
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
    
    if (estimator != PARTIAL_LIKELIHOOD) {
      occupancy_survivor = new double[max_nb_value];
      censored_occupancy_survivor = new double[max_nb_value + 1];
    }

    poutput = new double*[nb_variable];
    
    random_predict = new double **[nb_sequence];
    for(i = 0; i < nb_sequence; i++){
      random_predict[i] = new double*[max_nb_state_sequence];
      for(j = 0; j < max_nb_state_sequence; j++){
	random_predict[i][j] = new double[hsw_markov->nb_state];
      }
    }
    

    // HMM-LM ou HSMM-LM 
    // initialisation des paramètres de régression, variances résiduelles, probabilités initiales et probabilités de transition,

    hsw_markov = markov_switching_estimation(error, cout, ihsw_markov, VarCommune, estimator, false, false, false);

    tmp = 0;
    iter = 0;
    nb_likelihood_decrease = 0;

    do {
      previous_likelihood = likelihood;
      likelihood = 0.;

      // calcul du nombre de sequences d'etats simulees
      int nb_state_sequence_previous = 1;

      if(iter>0){
	nb_state_sequence_previous = nb_state_sequence;
      }

      if (min_nb_state_sequence + (int)::round(parameter * iter) < max_nb_state_sequence) {
	nb_state_sequence = min_nb_state_sequence + (int)::round(parameter * iter);
      }
      else {
	nb_state_sequence = max_nb_state_sequence;
      }

      indic = new double ***[nb_sequence];
      for(i = 0; i < nb_sequence; i++){
	indic[i] = new double**[nb_state_sequence];
	for(j = 0; j < nb_state_sequence; j++){
	  indic[i][j] = new double*[max_length];
	  for (k = 0; k < max_length; k++) {
	    indic[i][j][k] = new double[hsw_markov->nb_state];
	  }
	}
      }

      if( iter == 0) {
	for(i = 0; i < nb_sequence; i++){
	  for(j = 0; j < nb_state_sequence; j++){
	    for (k = 0; k < hsw_markov->nb_state; k++) {
	      random_predict[i][j][k] = effect[i][k];
	    }
	  }
	}
      }
      
      iter++;
      
      // initialisation des quantites de reestimation
      
      chain_reestim->init();

      for (i = 0;i < hsw_markov->nb_state;i++) {
        if (hsw_markov->state_subtype[i] == SEMI_MARKOVIAN) {
          reestim = complete_run[i]->frequency;
          for (j = 0;j < occupancy_nb_value[i];j++) {
            *reestim++ = 0.;
          }
	  
          reestim = final_run[i]->frequency;
          for (j = 0;j < occupancy_nb_value[i];j++) {
            *reestim++ = 0.;
          }
	  
	}
      }
      
      for (i = 0;i < nb_sequence;i++) {
	for (j = 0;j < nb_state_sequence;j++) {

	  // tirage des M^{(k+1)} effets aléatoires parmi les M^{(k)} effets aléatoires
 
	  int choixeffet = j+1;
	  if(nb_state_sequence != nb_state_sequence_previous){
	    choixeffet = rand()%nb_state_sequence_previous+1;
	  }
  
	  for (k = 0;k < nb_variable;k++) {
	    poutput[k] = real_sequence[i][k];
	  }
	  	  
	  // recurrence "forward"
	  
	  for (k = 0;k < length[i];k++) {
	    norm[k] = 0.;
	    
	    for (l = 0;l < hsw_markov->nb_state;l++) {
	      
	      // calcul des probabilites d'observation
	      
	      observation[k][l] = 1.;
	      if (iter == 1){
		observation[k][l] *= hsw_markov->sw_process[1]->observation[l]->density_computation(*poutput[0], covar[i][k], 
												    effect[i][l]);
	      }
	      else {
		observation[k][l] *= hsw_markov->sw_process[1]->observation[l]->density_computation(*poutput[0], covar[i][k], 
												    random_predict[i][choixeffet-1][l]);
	      }
	      
	      switch (hsw_markov->state_subtype[l]) {
		
		// cas etat semi-markovien
		
	      case SEMI_MARKOVIAN : {
		if (k == 0) {
		  state_norm[l] = hsw_markov->initial[l];
		}
		else {
		  state_norm[l] += state_in[k - 1][l] - forward[k - 1][l];
		}
		state_norm[l] *= observation[k][l];
		
		norm[k] += state_norm[l];
		break;
	      }
		
		// cas etat markovien
		
	      case MARKOVIAN : {
		if (k == 0) {
		  forward[k][l] = hsw_markov->initial[l];
		}
		else {
		  forward[k][l] = state_in[k - 1][l];
		}
		forward[k][l] *= observation[k][l];
		
		norm[k] += forward[k][l];
		break;
	      }
	      }
	    }
	    
	    if (norm[k] > 0.) {
	      for (l = 0;l < hsw_markov->nb_state;l++) {
		switch (hsw_markov->state_subtype[l]) {
		case SEMI_MARKOVIAN :
		  state_norm[l] /= norm[k];
		  break;
		case MARKOVIAN :
		  forward[k][l] /= norm[k];
		  break;
		}
	      }
	      
	      likelihood += log(norm[k]);
	    }
	    
	    else {
	      likelihood = D_INF;
	      break;
	    }
	    
  
	    for (l = 0;l < hsw_markov->nb_state;l++) {
	      
	      // cas etat semi-markovien
	      
	      if (hsw_markov->state_subtype[l] == SEMI_MARKOVIAN) {
		occupancy = hsw_markov->nonparametric_process[0]->sojourn_time[l];
		
		obs_product = 1.;
		forward[k][l] = 0.;
		
		if (k < length[i] - 1) {
		  for (m = 1;m <= MIN(k + 1 , occupancy->nb_value - 1);m++) {
		    obs_product *= observation[k - m + 1][l] / norm[k - m + 1];
		    if (obs_product == 0.) {
		      break;
		    }
		    
		    if (m < k + 1) {
		      forward[k][l] += obs_product * occupancy->mass[m] * state_in[k - m][l];
		    }
		    
		    else {
		      forward[k][l] += obs_product * occupancy->mass[m] * hsw_markov->initial[l];
		    }
		  }
		}
		
		else {
		  for (m = 1;m <= MIN(k + 1 , occupancy->nb_value - 1);m++) {
		    obs_product *= observation[k - m + 1][l] / norm[k - m + 1];
		    if (obs_product == 0.) {
		      break;
		    }
		    
		    if (m < k + 1) {
		      forward[k][l] += obs_product * (1. - occupancy->cumul[m - 1]) * state_in[k - m][l];
		    }
		    
		    else {
		      forward[k][l] += obs_product * (1. - occupancy->cumul[m - 1]) * hsw_markov->initial[l];
		    }
		  }
		}
	      }
	    }

	    if (k < length[i] - 1) {
	      for (l = 0;l < hsw_markov->nb_state;l++) {
		state_in[k][l] = 0.;
		for (m = 0;m < hsw_markov->nb_state;m++) {
		  state_in[k][l] += hsw_markov->transition[m][l] * forward[k][m];
		}
	      }
	    }
	    
	    for (l = 0;l < nb_variable;l++) {
	      poutput[l]++;
	    }
	  }
	  
	  if (likelihood == D_INF) {
	    break;
	  }
	  
#       ifdef DEBUG
	  for (k = 0;k < length[i];k++) {
	    cout << k << " : ";
	    for (l = 0;l < hsw_markov->nb_state;l++) {
	      cout << forward[k][l] << " ";
	    }
	    cout << endl;
	  }
	  cout << endl;
#       endif


	  // passes "backward"

	  k = length[i] - 1;
	  pstate = state_seq + k;
	  
	  for (m = 0;m < nb_variable;m++) {
	    poutput[m] = real_sequence[i][m] + k;
	  }

	  cumul_computation(hsw_markov->nb_state , forward[k] , cumul_backward);
	  *pstate = cumul_method(hsw_markov->nb_state , cumul_backward);
	  
	  do {
	    
            // cas etat semi-markovien
	    
            if (hsw_markov->state_subtype[*pstate] == SEMI_MARKOVIAN) {
              occupancy = hsw_markov->nonparametric_process[0]->sojourn_time[*pstate];
	      obs_product = 1.;
	      
              if (k < length[i] - 1) {
                for (m = 1;m <= MIN(k + 1 , occupancy->nb_value - 1);m++) {
                  obs_product *= observation[k - m + 1][*pstate] / norm[k - m + 1];
                  if (obs_product == 0.) {
                    break;
                  }
		  
                  if (m < k + 1) {
                    backward[m] = obs_product * occupancy->mass[m] * state_in[k - m][*pstate] / forward[k][*pstate];
                  }

                  else {
		    backward[m] = obs_product * occupancy->mass[m] * hsw_markov->initial[*pstate] / forward[k][*pstate];
		  }
                }
              }

              else {
                for (m = 1;m <= MIN(k + 1 , occupancy->nb_value - 1);m++) {
                  obs_product *= observation[k - m + 1][*pstate] / norm[k - m + 1];
                  if (obs_product == 0.) {
                    break;
                  }

                  if (m < k + 1) {
                    backward[m] = obs_product * (1. - occupancy->cumul[m - 1]) * state_in[k - m][*pstate] / forward[k][*pstate];
                  }

                  else {
		    backward[m] = obs_product * (1. - occupancy->cumul[m - 1]) * hsw_markov->initial[*pstate] / forward[k][*pstate];
                  }
                }
              }

              cumul_computation(m - 1 , backward + 1 , cumul_backward);
              state_occupancy = 1 + cumul_method(m - 1 , cumul_backward);

#             ifdef DEBUG
              double sum = 0.;
              for (n = 1;n < m;n++) {
                sum += backward[n];
              }
              if ((sum < 1. - DOUBLE_ERROR) || (sum > 1. + DOUBLE_ERROR)) {
                cout << "\nERROR: " << k << " " << sum << endl;
              }
#             endif

	      // accumulation des quantites de reestimation des lois d'occupation des etats

              if (k < length[i] - 1) {
                if (state_occupancy < k + 1) {
                  (complete_run[*pstate]->frequency[state_occupancy])++;
                }

                else {
		  (complete_run[*pstate]->frequency[state_occupancy])++;
                }
              }

              else {
                if (state_occupancy < k + 1) {
                  (final_run[*pstate]->frequency[state_occupancy])++;
                }

                else {
		  (final_run[*pstate]->frequency[state_occupancy])++;
                }
              }

              for (m = 1;m < state_occupancy;m++) {
                pstate--;
                *pstate = *(pstate + 1);
	      }

	      kk = k;
	      k -= (state_occupancy - 1);

	      for(m = (k+1); m <= kk; m++){
		for (n = 0; n < hsw_markov->nb_state; n++){
		  if ( *(pstate+1) == n) {
		    indic[i][j][m][n] = 1.;
		  }
		  else {
		    indic[i][j][m][n] = 0.;
		  }
		}
	      }

              if (k == 0) {
                break;
              }
	    }

            k--;

            for (m = 0;m < hsw_markov->nb_state;m++) {
              backward[m] = hsw_markov->transition[m][*pstate] * forward[k][m] / state_in[k][*pstate];
            }
            cumul_computation(hsw_markov->nb_state , backward , cumul_backward);
            *--pstate = cumul_method(hsw_markov->nb_state , cumul_backward);

	    for(m = 0; m < hsw_markov->nb_state; m++) {
	      if ( *(pstate+1) == m) {
		indic[i][j][k+1][m] = 1.;
	      }
	      else {
		indic[i][j][k+1][m] = 0.;
	      }
	    }

            // accumulation des quantites de reestimation des probabilites de transition

            (chain_reestim->transition[*pstate][*(pstate + 1)])++;

#           ifdef DEBUG
            double sum = 0.;
            for (m = 0;m < hsw_markov->nb_state;m++) {
              sum += backward[m];
	      cout<<backward[m]<<"  ";
            }
	    cout<<sum<<endl;
            if ((sum < 1. - DOUBLE_ERROR) || (sum > 1. + DOUBLE_ERROR)) {
              cout << "\nERROR: " << k << " " << sum << endl;
            }
#           endif

          }
          while (k > 0);


	  for(m = 0; m < hsw_markov->nb_state; m++) {
	    if ( *pstate == m) {
	      indic[i][j][0][m] = 1.;
	    }
	    else {
	      indic[i][j][0][m] = 0.;
	    }
	  }

          // accumulation des quantites de reestimation des probabilites initiales

	  (chain_reestim->initial[*pstate])++;



	}// fin de la boucle du nb_state_sequence
	
      }// fin de la boucle des individus


      // Prédiction des effets aléatoires


      double **mediane = new double*[hsw_markov->nb_state];
      for(i = 0; i < hsw_markov->nb_state; i++){
	mediane[i] = new double[nb_state_sequence];
      }  
      
      for (i = 0; i < nb_sequence; i++) {
	
	if(type_random == 0){
	  random_predict[i] =  conditional_expectation_heterogeneity_wise(hsw_markov, length, nb_covariable, covar, effect,
									  real_sequence, indic, nb_state_sequence, i);
	}
	else{
	  random_predict[i] = conditional_expectation_heterogeneity_state_wise(hsw_markov, length, nb_covariable, covar, effect,
									       real_sequence, indic, nb_state_sequence, i);
	}
	
	// il faut prendre la médiane de chaque état pour chaque inidvidu
	for (k = 0; k < hsw_markov->nb_state; k++) {
	  effect[i][k] = 0.;
	  for (j = 0; j < nb_state_sequence; j++){
	    mediane[k][j] = random_predict[i][j][k];
	  }
	}
	
	for (k = 0; k < hsw_markov->nb_state; k++){
	  gsl_sort(mediane[k],1,nb_state_sequence);
	  effect[i][k] = gsl_stats_median_from_sorted_data(mediane[k], 1, nb_state_sequence);
	} 
      }
      
      for (i = 0;i < hsw_markov->nb_state;i++) {
	delete [] mediane[i];
      }
      delete [] mediane;


      // estimation des parametres des modèles linéaires mixtes pour chacun des états 
 
      for (k = 0; k < hsw_markov->nb_state; k++) {
	double *beta_k,*tau_k, *sebeta_k, *setau_k;
	double  sigma_k, sesigma_k;
	
	beta_k = new double[nb_covariable];
	tau_k = new double[nb_random];
	sebeta_k = new double[nb_covariable];
	setau_k = new double[nb_random];

	// estimation des parametres de regression pour chacun des états 

	beta_k = regression_parameter_heterogeneity_stochastic(hsw_markov, nb_sequence, length, real_sequence, nb_covariable, 
							       covar, random_predict, indic, k, nb_state_sequence); 
	 
	// Calcul des variances résiduelles associées à chaque état
     
 	sigma_k = residual_variance_heterogeneity_stochastic(hsw_markov, nb_sequence, length, real_sequence, nb_covariable, covar, 
							     random_predict, indic, beta_k, k, nb_state_sequence);

	// Calcul des variances associées aux effets aléatoires 

	tau_k = random_variance_heterogeneity_stochastic(hsw_markov, nb_sequence, length, real_sequence, nb_covariable, covar, 
							 random_predict, indic, beta_k, k, nb_state_sequence);

#   ifdef DEBUG
	cout<<"beta pour l'état "<<k<<" :  "<<endl;
	for(i = 0; i < nb_covariable; i++){
	  cout<<beta_k[i]<<"   ";
	}
	cout<<endl;
#   endif
	 
	  
#   ifdef DEBUG
	cout<<"variance pour l'état "<<k<<" :  "<<sigma_k<<endl;
#   endif

 	// mise à jour des lois d'observations
	  
	hsw_markov->sw_process[1]->observation[k]->Set_regression(beta_k);
	hsw_markov->sw_process[1]->observation[k]->Set_residual_variance(sigma_k);
	hsw_markov->sw_process[1]->observation[k]->Set_random_variance(tau_k);	

	// calcul des standards errors des parametres de regression 
	sebeta_k = standard_error_regression_parameter_heterogeneity(hsw_markov, nb_sequence, length, real_sequence, nb_covariable, 
								   covar, nb_random, random_predict, indic, k, nb_state_sequence);

	// Calcul des standard errors des variances résiduelles
	sesigma_k = standard_error_residual_variance_heterogeneity(hsw_markov, nb_sequence, length, real_sequence, nb_covariable, 
								   covar, nb_random, random_predict, indic, k, nb_state_sequence);

	// Calcul des standard errors des ecarts types associées aux effets aléatoires 
	setau_k = standard_error_random_variance_heterogeneity(hsw_markov, nb_sequence, length, real_sequence, nb_covariable, 
							       covar, nb_random, random_predict, indic, k, nb_state_sequence);
	
	if (iter == nb_iter){
	  for (i = 0; i < nb_covariable; i++){
	    cout<<"standard deviation du paramètre de régression " << i << " pour l'état "<< k <<" :  "	<< sebeta_k[i] <<endl;
	  }
	  cout<<"standard deviation de la variance résiduelle pour l'état "<<k<<" :  "<< sesigma_k <<endl;
 	  cout<<"standard deviation du paramètre aléatoire pour l'état "<<k<<" :  "<< setau_k[0] <<endl;
	  cout<<endl;
	}

	delete [] setau_k;
	delete [] sebeta_k;
	delete [] tau_k;
	delete [] beta_k;
      }
      
  
    if (likelihood != D_INF) {
	  
	// reestimation des probabilites initiales

	reestimation(hsw_markov->nb_state , chain_reestim->initial ,
		     hsw_markov->initial , MIN_PROBABILITY , false);
  
	
	// reestimation des probabilites de transition
	  
	for (i = 0;i < hsw_markov->nb_state;i++) {
	  reestimation(hsw_markov->nb_state , chain_reestim->transition[i] ,
		       hsw_markov->transition[i] , MIN_PROBABILITY , false);
	}
	  
	// reestimation des lois d'occupation des etats
	  
	min_likelihood = 0.;
	  
	for (i = 0;i < hsw_markov->nb_state;i++) {
	  if (hsw_markov->state_subtype[i] == SEMI_MARKOVIAN) {
	    occupancy = hsw_markov->nonparametric_process[0]->sojourn_time[i];
	    
	    complete_run[i]->nb_value_computation();
	    complete_run[i]->offset_computation();
	    complete_run[i]->nb_element_computation();
	    
#           ifdef DEBUG
	    cout << "\n" << STAT_label[STATL_STATE] << " " << i << " ";
	    complete_run[i]->print(cout);
#           endif
	    
	    final_run[i]->nb_value_computation();
	    final_run[i]->offset_computation();
	    final_run[i]->nb_element_computation();
	    
	    if (final_run[i]->nb_element > 0.) {
	      complete_run[i]->state_occupancy_estimation(final_run[i] , complete_run[i] , occupancy_survivor ,
							  censored_occupancy_survivor , false);
	    }
	    
            complete_run[i]->max_computation();
            complete_run[i]->mean_computation();
            complete_run[i]->variance_computation();
	    
	    occupancy_likelihood = complete_run[i]->type_parametric_estimation(occupancy , 1 , true , OCCUPANCY_THRESHOLD);
 
#           ifdef DEBUG
	    cout << "\n" << STAT_label[STATL_STATE] << " " << i << " ";
   	    complete_run[i]->print(cout);
#           endif

#           ifdef DEBUG
	    cout << "\n" << STAT_label[STATL_STATE] << " " << i << " ";
	    final_run[i]->print(cout);
#           endif

#           ifdef DEBUG
            if (i == 1) {
              occupancy->print(cout);
            }
#           endif

            if (occupancy_likelihood == D_INF) {
              min_likelihood = D_INF;
            }
            else {
              occupancy->computation(complete_run[i]->nb_value , OCCUPANCY_THRESHOLD);
	    }

#           ifdef DEBUG
	    cout << endl<<STAT_word[STATW_STATE] << " " << i << endl;
            occupancy->ascii_print(cout);
#           endif
	    
          }
        }
      }

      likelihood=hsw_markov->likelihood_computation(*this, 0, I_DEFAULT);  
   
      //#     ifdef MESSAGE
      os << STAT_label[STATL_ITERATION] << " " << iter << "   "
	 << SEQ_label[SEQL_OBSERVED_SEQUENCES_LIKELIHOOD] << ": " << likelihood
	 << "   (" << nb_state_sequence << ")" << endl;
      //#     endif
   

      for(i = 0; i < nb_sequence; i++){
	for (j = 0; j < nb_state_sequence; j++){
	  for (k = 0; k < max_length; k++) {
	    delete [] indic[i][j][k];
	  }
	  delete [] indic[i][j];
	}
	delete [] indic[i];
      }
      delete [] indic;

      
#     ifdef DEBUG
      if (iter % 5 == 0) {
	cout << *hsw_markov;
      }
#     endif

    }
    while (iter<nb_iter);
    /*((likelihood != D_INF) && (((nb_iter == I_DEFAULT) && (iter < MARKOV_SWITCHING_NB_ITER) &&
				     ((likelihood - previous_likelihood) / -likelihood > MARKOV_SWITCHING_LIKELIHOOD_DIFF)
				        || (tmp < nb_iter_conv))
					|| ((nb_iter != I_DEFAULT) && (iter < nb_iter))));*/

    if (likelihood != D_INF) {

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

      // reestimation des lois d'observation non-parametriques

      for (j = 0;j < hsw_markov->nb_state;j++) {
	hsw_markov->sw_process[1]->observation[j]->parametric_variance_computation();
	hsw_markov->sw_process[1]->observation[j]->min_value_computation();
	hsw_markov->sw_process[1]->observation[j]->max_value_computation();
	hsw_markov->sw_process[1]->observation[j]->nb_value_computation();
      }
    }

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

    delete [] backward;
    delete [] cumul_backward;

    delete [] state_seq;

    delete chain_reestim;

    for (i = 0;i < hsw_markov->nb_state;i++) {
      delete complete_run[i];
    }
    delete [] complete_run;

    for (i = 0;i < hsw_markov->nb_state;i++) {
      delete final_run[i];
    }
    delete [] final_run;

    delete [] occupancy_nb_value;

    if (estimator != PARTIAL_LIKELIHOOD) {
      delete [] occupancy_survivor;
      delete [] censored_occupancy_survivor;
    }

    delete [] poutput;


    for(i = 0; i < nb_sequence; i++){
      for (j = 0; j < nb_state_sequence; j++){
	delete [] random_predict[i][j];
      }
      delete [] random_predict[i];
    }
    delete [] random_predict;


    if (likelihood == D_INF) {
      delete hsw_markov;
      hsw_markov = 0;
      error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
    }
    else {
      if  (state_sequence) {
        hsw_markov->sw_markov_data = new Markov_switching_data(*this , 0);

        hsw_markov_data = hsw_markov->sw_markov_data;
        hsw_markov_data->type[0] = STATE;

	if ((hsw_markov->sw_process[1]) && (hsw_markov_data->characteristics[1])) {
	  delete hsw_markov_data->characteristics[1];
	  hsw_markov_data->characteristics[1] = 0;
	}
        
	hsw_markov->create_cumul();
	hsw_markov_data->posterior_probability = new double[hsw_markov_data->nb_sequence];

	if(!output_R){
	  hsw_markov_data->likelihood = hsw_markov->viterbi(*hsw_markov_data, hsw_markov_data->posterior_probability, true, false);
	}
	else {
	  hsw_markov_data->likelihood = hsw_markov->viterbi(*hsw_markov_data, hsw_markov_data->posterior_probability, false, true);
	}

	hsw_markov->remove_cumul();
        hsw_markov_data->max_value[0] = hsw_markov->nb_state - 1;
        hsw_markov_data->build_marginal_histogram(0);
        hsw_markov_data->build_characteristic(0);
        hsw_markov_data->build_transition_count(*hsw_markov);
        
#       ifdef MESSAGE
	cout << "\n" << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << ": " << hsw_markov_data->likelihood <<endl;
#       endif

	// calcul des lois d'occupation des etats

        for (i = 0;i < hsw_markov->nb_state;i++) {
          if (hsw_markov->state_subtype[i] == SEMI_MARKOVIAN) {
            hsw_markov->nonparametric_process[0]->sojourn_time[i]->computation((hsw_markov_data->characteristics[0] ? 
										hsw_markov_data->characteristics[0]->sojourn_time[i]
										->nb_value : 1) ,
									       OCCUPANCY_THRESHOLD);
            if (hsw_markov->state_type[i] == 'r') {
	      hsw_markov->forward[i]->copy(*(hsw_markov->nonparametric_process[0]->sojourn_time[i]));
              hsw_markov->forward[i]->computation(*(hsw_markov->nonparametric_process[0]->sojourn_time[i]));
            }
          }
        }
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

      hsw_markov_data->hidden_likelihood = hsw_markov->likelihood_computation(*this , hsw_markov_data->posterior_probability, 
									      I_DEFAULT);

      hsw_markov->component_computation();
      hsw_markov->characteristic_computation(*hsw_markov_data , counting_flag , I_DEFAULT , false);

#     ifdef MESSAGE
      if  ((state_sequence) && (hsw_markov_data->nb_sequence <= POSTERIOR_PROBABILITY_NB_SEQUENCE)) {
        os << "\n" << SEQ_label[SEQL_POSTERIOR_STATE_SEQUENCE_PROBABILITY] << endl;
        for (i = 0;i < hsw_markov_data->nb_sequence;i++) {
          os << SEQ_label[SEQL_SEQUENCE] << " " << hsw_markov_data->identifier[i] << ": "
             << hsw_markov_data->posterior_probability[i];

          if (hsw_markov->nb_component == hsw_markov->nb_state) {
            os << " | " << SEQ_label[SEQL_STATE_BEGIN] << ": ";

            pstate = hsw_markov_data->int_sequence[i][0] + 1;
            if (hsw_markov_data->index_parameter) {
              for (j = 1;j < hsw_markov_data->length[i];j++) {
                if (*pstate != *(pstate - 1)) {
                  os << hsw_markov_data->index_parameter[i][j] << ", ";
                }
                pstate++;
              }
            }

            else {
              for (j = 1;j < hsw_markov_data->length[i];j++) {
                if (*pstate != *(pstate - 1)) {
                  os << j << ", ";
                }
                pstate++;
              }
            }
          }

          os << endl;
        }
      }
#     endif

    }
  }
  return hsw_markov;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'un semi-Markov switching linear 
 *  mixed model a partir d'un echantillon de sequences par 
 *  l'algorithme SEM/MCEM.
 *
 *  State sequences : "forward-backward" sampling
 *  Random effects : Prediction by conditionnal expectation 
 *
 *  arguments : reference sur un objet Format_error, stream, type de processus
 *              ('o' : ordinaire, 'e' : en equilibre), nombre d'etats de la chaine de Markov,
 *              flag sur la nature de la chaine de Markov, variance résiduelle, variance
 *              aléatoire, paramètres de regression, nombre de covariables, nombre
 *              d'effets aléatoires, presence/absence de constante, ordre, type d'effet
 *              aléatoire, nombre d'iterations pour convergence, parametres pour nombre de 
 *              sequences d'etats simulees, type Markovien, variance résiduelle commune à tous les états,
 *              sortie de type R, flags sur le calcul des lois de comptage et sur le calcul
 *              des sequences d'etats optimales, probabilite de rester dans un etat initiale,
 *              nombre d'iterations.
 *
 *--------------------------------------------------------------*/

Markov_switching* Switching_sequence::markov_switching_stochastic_estimation_hetero(Format_error &error, ostream &os, char type, 
										    int nb_state, bool left_right,
										    double *residual_variance, double **random_variance,
										    double **regression, int nb_covariable, int nb_random,
										    int constant, int order, int type_random, 
										    int nb_iter_conv, int min_nb_state_sequence,
										    int max_nb_state_sequence, bool markov, 
										    bool VarCommune, bool output_R, double parameter, 
										    bool counting_flag, bool state_sequence, 
										    double occupancy_mean, double self_transition, 
										    int nb_iter) const
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
    hsw_markov = markov_switching_stochastic_estimation_hetero(error, os, *ihsw_markov, type_random, nb_iter_conv, min_nb_state_sequence, 
							       max_nb_state_sequence, VarCommune, output_R,parameter, COMPLETE_LIKELIHOOD,
							       counting_flag, state_sequence, nb_iter);
    delete ihsw_markov;
  }

  return hsw_markov;
}






















