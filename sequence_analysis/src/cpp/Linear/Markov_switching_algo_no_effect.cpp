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
 * Correspond aux fonctions propres aux semi-Markov switching linear model
 * (modèle sans effets aléatoires)
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


/*--------------------------------------------------------------*
 *  Déjà dans renewal_algorithms.cpp
 *
 *  Calcul de la moyenne d'une loi par bissection d'intervalle.
 *
 *  arguments : pointeurs sur les quantites de reestimation de la 
 *              loi et de la loi biaisee par la longueur.
 *
 *--------------------------------------------------------------*/

double interval_bisection(Reestimation<double> *distribution_reestim ,
                          Reestimation<double> *length_bias_reestim)

{
  register int i;
  double ratio , inf_ratio , sup_ratio , mean , inf_mean , sup_mean , *dfrequency , *lfrequency;

# ifdef DEBUG
  int iter = 0;
# endif

  // initialisations : calculs des 2 premieres valeurs
  dfrequency = distribution_reestim->frequency + distribution_reestim->offset;
  lfrequency = length_bias_reestim->frequency + distribution_reestim->offset;
  inf_ratio = 0.;
  sup_ratio = 0.;
  inf_mean = distribution_reestim->offset;
  sup_mean = distribution_reestim->nb_value - 1;

  for (i = distribution_reestim->offset;i < distribution_reestim->nb_value;i++) {
    inf_ratio += (*dfrequency + *lfrequency) * i /
                 (distribution_reestim->nb_element * sup_mean + length_bias_reestim->nb_element * i);
    sup_ratio += (*dfrequency++ + *lfrequency++) * i /
                 (distribution_reestim->nb_element * inf_mean + length_bias_reestim->nb_element * i);
  }

  do {
    dfrequency = distribution_reestim->frequency + distribution_reestim->offset;
    lfrequency = length_bias_reestim->frequency + distribution_reestim->offset;
    ratio = 0.;
    mean = (inf_mean + sup_mean) / 2.;

    for (i = distribution_reestim->offset;i < distribution_reestim->nb_value;i++) {
      ratio += (*dfrequency++ + *lfrequency++) * i /
               (distribution_reestim->nb_element * mean + length_bias_reestim->nb_element * i);
    }

#   ifdef DEBUG
    cout << STAT_label[STATL_ITERATION] << " " << iter++ << ": " << mean << " " << ratio << endl;
#   endif

    if (ratio < 1.) {
      inf_ratio = ratio;
      sup_mean = mean;
    }
    else {
      sup_ratio = ratio;
      inf_mean = mean;
    }
  }
  while (sup_ratio - inf_ratio > BISECTION_RATIO_THRESHOLD);

  mean = (inf_mean + sup_mean) / 2.;

# ifdef DEBUG
  cout << STAT_label[STATL_MEAN] << ": " << mean << " " << ratio << endl;
# endif

  return mean;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par un semi-Markov switching linear model.
 *
 *  arguments : reference sur un objet Format_error, histogramme 
 *              des longueurs des sequences, nombre de covariables, 
 *              présence/absence de constante, variance résiduelle,
 *              paramètres de regression, covariables,
 *              flag sur le calcul des lois de comptage,
 *              flag calcul d'une divergence de Kullback-Leibler.
 *
 *--------------------------------------------------------------*/

Markov_switching_data* Markov_switching::simulation(Format_error &error , const Histogram &hlength ,
						    int inb_covariable, int iconstant, 
						    double *iresidual_variance, double **iregression,
						    double ***icovar, bool counting_flag , bool divergence_flag) const
{
  bool status = true;
  register int i , j , k, m;
  int cumul_length , occupancy, *pstate;
  int nb_covar;
  double **poutput;
  double *pcovar, *ccovar;
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

    seq = new Markov_switching_data( 2 , hlength , inb_covariable, 0, 0, 0, iconstant, false);
    seq->type[0] = STATE;

    seq->sw_markov = new Markov_switching(*this , false);
    sw_markov = seq->sw_markov;


    sw_markov->create_cumul();
    sw_markov->cumul_computation();

    poutput = new double*[1];
    
    for (i = 0; i < sw_markov->nb_state; i++) {
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
  
    for (i = 0;i < seq->nb_sequence;i++) {
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

        for (k = 0;k < occupancy;k++) {
	  if (sw_markov->sw_process[1]) {
	    *poutput[0]++ = sw_markov->sw_process[1]->observation[*pstate]->simulation(seq->covar[i][j]);
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
      //  seq->likelihood = sw_markov->likelihood_computation(*seq); 

    }

  }

  return seq;
}



/*--------------------------------------------------------------*
 *
 *  Simulation par un semi-Markov switching linear model.
 *
 *  arguments : reference sur un objet Format_error, nombre et 
 *              longueur des sequences, nombre de covariables, 
 *              présence/absence de constante, variance résiduelle,
 *              paramètres de régression, covariables,
 *              flag sur le calcul des lois de comptage.
 *
 *--------------------------------------------------------------*/

Markov_switching_data* Markov_switching::simulation(Format_error &error , int nb_sequence ,
						    int length , int inb_covariable, int iconstant, 
						    double *iresidual_variance, double **iregression,
						    double ***icovar, bool counting_flag) const
{
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

    seq = simulation(error , hlength , inb_covariable, iconstant, iresidual_variance, iregression, icovar, counting_flag, false);
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par un semi-Markov switching linear model.
 *
 *  arguments : reference sur un objet Format_error, nombre de sequences,
 *              reference sur un objet Switching_sequences,
 *              nombre de covariable,s présence/absence de constante,
 *              variance résiduelle, paramètres de régression,
 *              covariables, flag sur le calcul des lois de comptage.
 *
 *--------------------------------------------------------------*/

Markov_switching_data* Markov_switching::simulation(Format_error &error , int nb_sequence ,
						    const Switching_sequence &iseq ,
						    int inb_covariable, int iconstant, 
						    double *iresidual_variance, double **iregression,
						    double ***icovar, bool counting_flag) const
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

    seq = simulation(error , *hlength , inb_covariable, iconstant, iresidual_variance,
		     iregression, icovar, counting_flag, false);
    delete hlength;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par un semi-Markov switching linear model.
 *
 *  arguments : reference sur un objet Format_error, histogramme 
 *              des longueurs des sequences, nombre de covariables,
 *              présence/absence de constante, variance résiduelle,
 *              paramètres de régression, covariables,
 *              flag sur le calcul des lois de comptage,
 *              flag calcul d'une divergence de Kullback-Leibler.
 *
 *-------------------------------------------------------------*/

Markov_switching_data* Markov_switching::simulation(Format_error &error ,
                                                    const Histogram &hlength ,
						    int inb_covariable, int iconstant, 
						    double *iresidual_variance, double **iregression,
						    double ***icovar, bool counting_flag ,
                                                    bool divergence_flag, bool hidden) const
{
  Switching_sequence *observ_seq;
  Markov_switching_data *seq;


  seq = Markov_switching::simulation(error , hlength , inb_covariable, iconstant,
				     iresidual_variance, iregression, icovar, counting_flag , divergence_flag);


  if ((seq) && (!divergence_flag)) {
    observ_seq = seq->remove_variable_1();
    seq->hidden_likelihood = likelihood_computation(*observ_seq, 0, I_DEFAULT);
    delete observ_seq;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par un semi-Markov switching linear model.
 *
 *  arguments : reference sur un objet Format_error, nombre et 
 *              longueur des sequences, nombre de covariables, 
 *              présence/absence de constante, variance résiduelle,
 *              paramètres de régression, covariables, flag sur 
 *              le calcul des lois de comptage.
 *
 *--------------------------------------------------------------*/

Markov_switching_data* Markov_switching::simulation(Format_error &error ,
                                                    int nb_sequence , int length ,
						    int inb_covariable, int iconstant, 
						    double *iresidual_variance, double **iregression,
						    double ***icovar, bool counting_flag, bool hidden) const
{
  Switching_sequence *observ_seq;
  Markov_switching_data *seq;

  seq = Markov_switching::simulation(error , nb_sequence , length , inb_covariable, iconstant,
				     iresidual_variance, iregression, icovar, counting_flag);

  if (seq) {
    observ_seq = seq->remove_variable_1();
    seq->hidden_likelihood = likelihood_computation(*observ_seq, 0, I_DEFAULT);
    delete observ_seq;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par un semi-Markov switching linear model.
 *
 *  arguments : reference sur un objet Format_error, nombre de sequences,
 *              reference sur un objet Switching_sequences, nombre
 *              de covariables, présence/absence de constante, 
 *              variance résiduelle, paramètres de régression, 
 *              covariables, flag sur le calcul des lois de comptage.
 *
 *--------------------------------------------------------------*/

Markov_switching_data* Markov_switching::simulation(Format_error &error ,
                                                    int nb_sequence ,
                                                    const Switching_sequence &iseq ,
						    int inb_covariable, int iconstant, 
						    double *iresidual_variance, double **iregression,
						    double ***icovar, bool counting_flag, bool hidden) const
{
  Switching_sequence *observ_seq;
  Markov_switching_data *seq;

  seq = Markov_switching::simulation(error , nb_sequence , iseq , inb_covariable, iconstant,
				     iresidual_variance, iregression, icovar,  counting_flag);

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
 *  a partir d'un echantillon de sequences par l'algorithme EM.
 *
 *  arguments : reference sur un objet Format_error, stream, semi-Markov 
 *              linear switching initial, variance résiduelle commune
 *              à tous les états, type d'estimateur pour la 
 *              reestimation des lois d'occupation des etats, sortie
 *              de type R, flags sur le calcul des lois de comptage 
 *              et sur le calcul des sequences d'etats optimales, 
 *              nombre d'iterations.
 *
 *--------------------------------------------------------------*/

Markov_switching* Switching_sequence::markov_switching_estimation(Format_error &error , ostream &os ,
								  const Markov_switching &isw_markov , bool VarCommune, 
								  int estimator, bool output_R,
								  bool counting_flag , bool state_sequence ,
								  int nb_iter, int mean_computation) const
{
  bool status;
  register int i , j , k , m, n;
  int max_nb_value , iter, nb_likelihood_decrease, offset, nb_value, *occupancy_nb_value,
    *censored_occupancy_nb_value;
  double **poutput;
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


  if (status) {
    if (max_length > COUNTING_MAX_LENGTH) {
      counting_flag = false;
    }

    // creation de la chaine de Markov cachee

    hsw_markov = new Markov_switching(isw_markov , false, (int)(max_length * SAMPLE_NB_VALUE_COEFF));


#   ifdef DEBUG
    cout << *hsw_markov;
#   endif


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
    
    iter = 0;
    nb_likelihood_decrease = 0;

    do {
      iter++;

      previous_likelihood = likelihood;
      likelihood = 0.;
      
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

	    observation[j][k] *= hsw_markov->sw_process[1]->observation[k]->density_computation(*poutput[0], covar[i][j]);

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
	    double sum1 = 0;
	    cout << j << " : ";
	    for (k = 0;k < hsw_markov->nb_state;k++) {
	      sum1 = sum1 + forward[j][k];
	      cout << forward[j][k] << " ";
	    }
	    cout<<": "<<sum1;
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
	  
	  for (k = 0;k < hsw_markov->nb_state;k++) {
	    if (hsw_markov->state_subtype[k] == SEMI_MARKOVIAN) {
              if (j < length[i] - 1) {
                test[k][0] += backward1[j][k];
                test[k][1] += auxiliary[k] * state_in[j][k];
              }
              else {
                test[k][2] += backward[i][j][k];
              }
              if (j == 0) {
                test[k][3] += backward[i][j][k];
              }
            }
          }
        }
#       endif

      }

# ifdef DEBUG
      cout<<*chain_reestim<<endl;
# endif


      // estimation des parametres des modèles linéaires pour chacun des états 

      // Cas avec Variance différente pour chaque état	
      if (!VarCommune){
	for (k = 0; k < hsw_markov->nb_state; k++) {

	  double *beta_k, *sebeta_k;
	  double  sigma_k, sesigma_k;
	
	  beta_k = new double[nb_covariable];
	  sebeta_k = new double[nb_covariable];
	  
	  beta_k = regression_parameter_no_effect(hsw_markov, nb_sequence, length, real_sequence, nb_covariable, covar, backward, k); 
	  
#   ifdef DEBUG
	  cout<<"beta pour l'état "<<k<<" :  "<<endl;
	  for(i = 0; i< nb_covariable; i++){
	    cout<<beta_k[i]<<"   ";
	  }
	  cout<<endl;
#endif
	  
	  sigma_k = residual_variance_no_effect(hsw_markov, nb_sequence, length, real_sequence, nb_covariable, 
						covar, backward, beta_k, k); 

#   ifdef DEBUG
	  cout<<"variance pour l'état "<<k<<" :  "<<sigma_k<<endl;
#   endif
	  
	  // mise à jour des lois d'observations
	  
	  hsw_markov->sw_process[1]->observation[k]->Set_regression(beta_k);
	  hsw_markov->sw_process[1]->observation[k]->Set_residual_variance(sigma_k);
	  
	  // calcul des standards errors des parametres de regression 
	  sebeta_k = standard_error_regression_parameter(hsw_markov, nb_sequence, length, real_sequence, nb_covariable,covar, backward, k);

	  // Calcul des standard errors des variances résiduelles
	  sesigma_k = standard_error_residual_variance(hsw_markov, nb_sequence, length, real_sequence, nb_covariable, covar, backward, k);

	
	  if (iter == nb_iter){
	    for (i = 0; i < nb_covariable; i++){
	      cout<<"standard deviation du paramètre de régression " << i << " pour l'état "<< k <<" :  "	<< sebeta_k[i] <<endl;
	    }
	    cout<<"standard deviation de la variance résiduelle pour l'état "<<k<<" :  "<< sesigma_k <<endl;
	    cout<<endl;
	  }

	  delete [] sebeta_k;
	  delete [] beta_k;
	}
      }

      else{
	// Cas avec Variance commune pour chaque état	
	
	double **beta;
	double sigma = 0.;
	
	beta = new double*[hsw_markov->nb_state];
	
	for (k = 0; k < hsw_markov->nb_state; k++){
	  beta[k] = new double[nb_covariable];
	  beta[k] = regression_parameter_no_effect(hsw_markov, nb_sequence, length, real_sequence, nb_covariable, covar, backward, k); 
	  
#   ifdef DEBUG
	  cout<<"beta pour l'état "<<k<<" :  "<<endl;
	  for(i = 0; i< nb_covariable; i++){
	    cout<<beta_k[i]<<"   ";
	  }
	  cout<<endl;
#endif
	  
	}
	
	sigma = residual_variance_no_effect_VarCommune(hsw_markov, nb_sequence, length, real_sequence, nb_covariable, 
						       covar, backward, beta); 
	
#   ifdef DEBUG
	cout<<"variance résiduelle  :  "<< sigma << endl;
#   endif
	
	// mise à jour des lois d'observations
	
	for (k = 0; k < hsw_markov->nb_state; k++){
	  hsw_markov->sw_process[1]->observation[k]->Set_regression(beta[k]);
	  hsw_markov->sw_process[1]->observation[k]->Set_residual_variance(sigma);
	}
	


	for (k = 0; k < hsw_markov->nb_state; k++){
	  double *sebeta_k;
	  sebeta_k = new double[nb_covariable];
	  
	  // calcul des standards errors des parametres de regression 
	  
	  sebeta_k = standard_error_regression_parameter(hsw_markov, nb_sequence, length, real_sequence, nb_covariable,covar, backward, k);

	  if (iter == nb_iter){
	    for (i = 0; i < nb_covariable; i++){
	      cout<<"standard deviation du paramètre de régression " << i << " pour l'état "<< k <<" :  "	<< sebeta_k[i] <<endl;
	    }
	  }
	  
	  delete [] sebeta_k;
	}
	  

	// Calcul des standard errors des variances résiduelles
	double sesigma;
	sesigma = standard_error_residual_variance_VarCommune(hsw_markov, nb_sequence, length, real_sequence, nb_covariable, 
							      covar, backward);
	  
	  
	if (iter == nb_iter){
	  cout<<"standard deviation de la variance résiduelle :  "<< sesigma <<endl;
	  cout<<endl;
	}
	 

	for (k = 0; k < hsw_markov->nb_state; k++){
	  delete [] beta[k];
	}
	delete [] beta;
      }
 



      if (likelihood != D_INF) {
	if (likelihood < previous_likelihood) {
          nb_likelihood_decrease++;
        }
        else {
          nb_likelihood_decrease = 0;
        }
	
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
              occupancy_likelihood = hoccupancy->Reestimation<int>::parametric_estimation(occupancy , 1 , true ,
                                                                                          OCCUPANCY_THRESHOLD);
            }
            else {
              occupancy_likelihood = hoccupancy->Reestimation<int>::type_parametric_estimation(occupancy , 1 , true ,
                                                                                               OCCUPANCY_THRESHOLD);
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
      
      //#     ifdef MESSAGE
      os << STAT_label[STATL_ITERATION] << " " << iter << "   "
	 << SEQ_label[SEQL_OBSERVED_SEQUENCES_LIKELIHOOD] << ": " << likelihood << endl;
      //#     endif



      
#     ifdef DEBUG
      if (iter % 5 == 0) {
	cout << *hsw_markov;
      }
#     endif

    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < max_length;j++) {
        delete backward[i][j];
      }
      delete [] backward[i];
    }
    delete [] backward;

    }
    while ((likelihood != D_INF) && ((likelihood - previous_likelihood) / -likelihood > MARKOV_SWITCHING_LIKELIHOOD_DIFF) &&
	   ((nb_iter == I_DEFAULT) && (iter < MARKOV_SWITCHING_NB_ITER) || ((nb_iter != I_DEFAULT) && (iter < nb_iter))));


 
    if (likelihood != D_INF) {
      
#     ifdef MESSAGE
      os << "\n" << iter << " " << STAT_label[STATL_ITERATIONS] << endl;
#     endif
      
      // reestimation des probabilites initiales
      
      reestimation(hsw_markov->nb_state , chain_reestim->initial ,
		   hsw_markov->initial , MIN_PROBABILITY , true);
      
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
    for(i = 0; i< hsw_markov->nb_output_process; i++){
      for (j = 0; j< hsw_markov->nb_state; j++) {
	cout<< hsw_markov->sw_process[i + 1]->observation[j]->residual_variance<<"   ";
      }
      cout<<endl;
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
										sojourn_time[i]->nb_value : 1) ,
									       OCCUPANCY_THRESHOLD);
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

      hsw_markov_data->hidden_likelihood = hsw_markov->likelihood_computation(*this, hsw_markov_data->posterior_probability,
									      I_DEFAULT);
      
      
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
      }
#     endif

 
      hsw_markov->component_computation();
      hsw_markov->characteristic_computation(*hsw_markov_data , counting_flag , I_DEFAULT , false);
    
    }
  }
  
  return hsw_markov;
  
}    


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'un semi-Markov switching linear a partir
 *  d'un echantillon de sequences.
 *
 *  arguments : reference sur un objet Format_error, stream, nombre d'etats
 *              de la semi-chaine de Markov, flag sur la nature de la 
 *              semi-chaine de Markov, variance résiduelle, paramètres 
 *              de regression, nombre de covariables, présence/absence
 *              de constante, flag sur le type Markovien, variance commune 
 *              à tous les états, sortie de type R, ordre,flags sur le calcul
 *              des lois de comptage et sur le calcul des sequences d'etats optimales,
 *              probabilite de rester dans un etat initiale, nombre d'iterations.
 *
 *--------------------------------------------------------------*/

Markov_switching* Switching_sequence::markov_switching_estimation(Format_error &error , ostream &os , 
								  int nb_state , bool left_right ,  
								  double *residual_variance, double **regression,
								  int nb_covariable, int constant, bool markov, 
								  bool VarCommune, bool output_R, int order, bool counting_flag , 
								  bool state_sequence , double occupancy_mean, double self_transition , 
								  int nb_iter, int mean_computation) const
{
  bool status = true;
  register int i,j;
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

    if(!constant) {
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
      cont_param[i] = new Continuous_parametric(NB_VALUE, nb_covariable, 0, 
						residual_variance[i], 0, regression[i], NO_RANDOM); 
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

    hsw_markov = new Markov_switching(*ihsw_markov, true);
    hsw_markov = markov_switching_estimation(error , os , *ihsw_markov , VarCommune, COMPLETE_LIKELIHOOD, output_R, counting_flag ,
					     state_sequence , nb_iter, mean_computation);
    delete ihsw_markov;
  }

  return hsw_markov;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'un semi-Markov switching linear model
 *  a partir d'un echantillon de sequences par l'algorithme SEM/MCEM.
 *
 *  arguments : reference sur un objet Format_error, stream, 
 *              semi-Markov switching linear model initial,
 *              parametres pour le nombre de sequences d'etats simulees, 
 *              variance résiduelle commune à tous les états, sortie de 
 *              type R, type d'estimateur pour la reestimation 
 *              des lois d'occupation des etats, flags sur le calcul
 *              des lois de comptage et sur le calcul des sequences 
 *              d'etats optimales, nombre d'iterations.
 *
 *--------------------------------------------------------------*/

Markov_switching* Switching_sequence::markov_switching_stochastic_estimation(Format_error &error , ostream &os ,
									     const Markov_switching &ihsw_markov ,
									     int min_nb_state_sequence , int max_nb_state_sequence ,
									     bool VarCommune, bool output_R, double parameter , 
									     int estimator, bool counting_flag , bool state_sequence , 
									     int nb_iter) const
{
  bool status;
  register int i , j , k , m , n, r;
  int max_nb_value , iter , nb_state_sequence , state_occupancy , nb_likelihood_decrease ,
    *occupancy_nb_value , *state_seq , *pstate ;
  double **poutput;
  double likelihood = D_INF , previous_likelihood , occupancy_likelihood , observation_likelihood ,
         min_likelihood , obs_product , **observation , *norm , *state_norm , **forward ,
         **state_in , *backward , *cumul_backward , *reestim , *occupancy_survivor ,
         *censored_occupancy_survivor;
  double ****indic;
  Chain_reestimation<double> *chain_reestim;
  Parametric *occupancy;
  Reestimation<double> *bcomplete_run , *censored_run , **complete_run , **final_run , **initial_run ,
                       **single_run;
  Markov_switching *hsmarkov;
  Markov_switching_data *seq;
  const Reestimation<double> *prun[3];


# ifdef DEBUG
  double sum;
# endif

  hsmarkov = 0;
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

  if (status) {
    if (max_length > COUNTING_MAX_LENGTH) {
      counting_flag = false;
    }

    // creation de la chaine de Markov cachee

    hsmarkov = new Markov_switching(ihsw_markov , false, (int)(max_length*SAMPLE_NB_VALUE_COEFF));

#   ifdef DEBUG
    cout << *hsmarkov;
#   endif

    // initialisations

    observation = new double*[max_length];
    for (i = 0;i < max_length;i++) {
      observation[i] = new double[hsmarkov->nb_state];
    }

    norm = new double[max_length];
    state_norm = new double[hsmarkov->nb_state];

    forward = new double*[max_length];
    for (i = 0;i < max_length;i++) {
      forward[i] = new double[hsmarkov->nb_state];
    }

    state_in = new double*[max_length - 1];
    for (i = 0;i < max_length - 1;i++) {
      state_in[i] = new double[hsmarkov->nb_state];
    }

    backward = new double[max_length + 1];
    cumul_backward = new double[max_length + 1];

    state_seq = new int[max_length];

    chain_reestim = new Chain_reestimation<double>('o', hsmarkov->nb_state , hsmarkov->nb_state);

    occupancy_nb_value = new int[hsmarkov->nb_state];
    complete_run = new Reestimation<double>*[hsmarkov->nb_state];
    final_run = new Reestimation<double>*[hsmarkov->nb_state];

    for (i = 0;i < hsmarkov->nb_state;i++) {
      switch (hsmarkov->state_subtype[i]) {

      case SEMI_MARKOVIAN : {
        occupancy_nb_value[i] = MIN(hsmarkov->nonparametric_process[0]->sojourn_time[i]->alloc_nb_value ,
                                    max_length + 1);

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
    for (i = 0;i < hsmarkov->nb_state;i++) {
      if ((hsmarkov->state_subtype[i] == SEMI_MARKOVIAN) && (occupancy_nb_value[i] > max_nb_value)) {
        max_nb_value = occupancy_nb_value[i];
      }
    }

    if (estimator != PARTIAL_LIKELIHOOD) {
      occupancy_survivor = new double[max_nb_value];
      censored_occupancy_survivor = new double[max_nb_value + 1];
    }


    poutput = new double*[nb_variable];


    iter = 0;
    nb_likelihood_decrease = 0;

    do {
      previous_likelihood = likelihood;
      likelihood = 0.;

      // calcul du nombre de sequences d'etats simulees

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
	    indic[i][j][k] = new double[hsmarkov->nb_state];
	  }
	}
      }

      iter++;

      // initialisation des quantites de reestimation

      chain_reestim->init();

      for (i = 0;i < hsmarkov->nb_state;i++) {
        if (hsmarkov->state_subtype[i] == SEMI_MARKOVIAN) {
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
        for (j = 0;j < nb_variable;j++) {
          poutput[j] = real_sequence[i][j];
        }


        // recurrence "forward"

	for (j = 0;j < length[i];j++) {
          norm[j] = 0.;

          for (k = 0;k < hsmarkov->nb_state;k++) {

            // calcul des probabilites d'observation

            observation[j][k] = 1.;
	    observation[j][k] *= hsmarkov->sw_process[1]->observation[k]->density_computation(*poutput[0], covar[i][j]);
        

            switch (hsmarkov->state_subtype[k]) {

            // cas etat semi-markovien

            case SEMI_MARKOVIAN : {
              if (j == 0) {
                state_norm[k] = hsmarkov->initial[k];
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
                forward[j][k] = hsmarkov->initial[k];
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
            for (k = 0;k < hsmarkov->nb_state;k++) {
              switch (hsmarkov->state_subtype[k]) {
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

          for (k = 0;k < hsmarkov->nb_state;k++) {

            // cas etat semi-markovien

            if (hsmarkov->state_subtype[k] == SEMI_MARKOVIAN) {
              occupancy = hsmarkov->nonparametric_process[0]->sojourn_time[k];
              obs_product = 1.;
              forward[j][k] = 0.;

              if (j < length[i] - 1) {
                for (m = 1;m <= MIN(j + 1 , occupancy->nb_value - 1);m++) {
                  obs_product *= observation[j - m + 1][k] / norm[j - m + 1];
                  if (obs_product == 0.) {
                    break;
                  }

                  if (m < j + 1) {
                    forward[j][k] += obs_product * occupancy->mass[m] * state_in[j - m][k];
                  }

                  else {
                      forward[j][k] += obs_product * occupancy->mass[m] * hsmarkov->initial[k];
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
		    forward[j][k] += obs_product * (1. - occupancy->cumul[m - 1]) * hsmarkov->initial[k];
		  }
                }
              }
            }
          }

          if (j < length[i] - 1) {
            for (k = 0;k < hsmarkov->nb_state;k++) {
              state_in[j][k] = 0.;
              for (m = 0;m < hsmarkov->nb_state;m++) {
                state_in[j][k] += hsmarkov->transition[m][k] * forward[j][m];
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
        for (j = 0;j < length[i];j++) {
          cout << j << " : ";
          for (k = 0;k < hsmarkov->nb_state;k++) {
            cout << forward[j][k] << " ";
          }
          cout << endl;
        }
        cout << endl;
#       endif

        // passes "backward"

        for (j = 0;j < nb_state_sequence;j++) {
          k = length[i] - 1;
          pstate = state_seq + k;

          for (m = 0;m < nb_variable;m++) {
            poutput[m] = real_sequence[i][m] + k;
          }

	  cumul_computation(hsmarkov->nb_state , forward[k] , cumul_backward);
	  *pstate = cumul_method(hsmarkov->nb_state , cumul_backward);


	  do {

            // cas etat semi-markovien

            if (hsmarkov->state_subtype[*pstate] == SEMI_MARKOVIAN) {
              occupancy = hsmarkov->nonparametric_process[0]->sojourn_time[*pstate];
              obs_product = 1.;

              if (k < length[i] - 1) {
                for (m = 1;m <= MIN(k + 1 , occupancy->nb_value - 1);m++) {
                  obs_product *= observation[k - m + 1][*pstate] / norm[k - m + 1];
                  if (obs_product == 0.) {
                    break;
                  }

                  if (m < k + 1) {
                    backward[m] = obs_product * occupancy->mass[m] * state_in[k - m][*pstate] /
                                  forward[k][*pstate];
                  }

                  else {
		    backward[m] = obs_product * occupancy->mass[m] * hsmarkov->initial[*pstate] /
		      forward[k][*pstate];
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
                    backward[m] = obs_product * (1. - occupancy->cumul[m - 1]) * state_in[k - m][*pstate] /
                                  forward[k][*pstate];
                  }

                  else {
		    backward[m] = obs_product * (1. - occupancy->cumul[m - 1]) *
		      hsmarkov->initial[*pstate] / forward[k][*pstate];
                  }
                }
              }

              cumul_computation(m - 1 , backward + 1 , cumul_backward);
              state_occupancy = 1 + cumul_method(m - 1 , cumul_backward);

#             ifdef DEBUG
              sum = 0.;
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

	      r = k;
	      k -= (state_occupancy - 1);

	      for(m = k; m <= r; m++){
		for (n = 0; n < hsmarkov->nb_state; n++){
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
            for (m = 0;m < hsmarkov->nb_state;m++) {
              backward[m] = hsmarkov->transition[m][*pstate] * forward[k][m] / state_in[k][*pstate];
            }
            cumul_computation(hsmarkov->nb_state , backward , cumul_backward);
            *--pstate = cumul_method(hsmarkov->nb_state , cumul_backward);

	    for(m = 0; m < hsmarkov->nb_state; m++) {
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
            sum = 0.;
            for (m = 0;m < hsmarkov->nb_state;m++) {
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


	  for(m = 0; m < hsmarkov->nb_state; m++) {
	    if ( *pstate == m) {
	      indic[i][j][0][m] = 1.;
	    }
	    else {
	      indic[i][j][0][m] = 0.;
	    }
	  }

          // accumulation des quantites de reestimation des probabilites initiales

	  (chain_reestim->initial[*pstate])++;
        }
      }

      // estimation des parametres de regression pour chacun des états 
      
      // Cas avec Variance différente pour chaque état	
      if (!VarCommune){
	for (k = 0; k < hsmarkov->nb_state; k++) {
	  
	  double *beta_k;
	  double sigma_k = 0.;
	  
	  beta_k = new double[nb_covariable];
	  
	  beta_k = regression_parameter_no_effect_stochastic(hsmarkov, nb_sequence, length, real_sequence, nb_covariable, 
							     covar, indic, k, nb_state_sequence); 
	  
#   ifdef DEBUG
	  cout<<"beta pour l'état "<<k<<" :  "<<endl;
	  for(i = 0; i< nb_covariable; i++){
	    cout<<beta_k[i]<<"   ";
	  }
	  cout<<endl;
#   endif
	  
	  sigma_k = residual_variance_no_effect_stochastic(hsmarkov, nb_sequence, length, real_sequence, nb_covariable, 
							   covar, indic, beta_k, k, nb_state_sequence); 

#   ifdef DEBUG
	  cout<<"variance pour l'état "<<k<<" :  "<<sigma_k<<endl;
#   endif
	  
	  // mise à jour des lois d'observations
	  
	  hsmarkov->sw_process[1]->observation[k]->Set_regression(beta_k);
	  hsmarkov->sw_process[1]->observation[k]->Set_residual_variance(sigma_k);
	  
	  delete [] beta_k;
	}
      }

      else{
	// Cas avec Variance commune pour chaque état	
	
	double **beta;
	double sigma = 0.;
	
	beta = new double*[hsmarkov->nb_state];
	
	for (k = 0; k < hsmarkov->nb_state; k++){
	  beta[k] = new double[nb_covariable];
	  beta[k] = regression_parameter_no_effect_stochastic(hsmarkov, nb_sequence, length, real_sequence, nb_covariable, 
							     covar, indic, k, nb_state_sequence); 
	  
#   ifdef DEBUG
	  cout<<"beta pour l'état "<<k<<" :  "<<endl;
	  for(i = 0; i< nb_covariable; i++){
	    cout<<beta_k[i]<<"   ";
	  }
	  cout<<endl;
#endif
	  
	}
	
	sigma = residual_variance_no_effect_VarCommune_stochastic(hsmarkov, nb_sequence, length, real_sequence, nb_covariable, 
								  covar, indic, beta, nb_state_sequence); 
	
#   ifdef DEBUG
	cout<<"variance résiduelle  :  "<< sigma << endl;
#   endif
	
	// mise à jour des lois d'observations
	
	for (k = 0; k < hsmarkov->nb_state; k++){
	  hsmarkov->sw_process[1]->observation[k]->Set_regression(beta[k]);
	  hsmarkov->sw_process[1]->observation[k]->Set_residual_variance(sigma);
	}
	
	delete [] beta;
      }
 
	      


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


      if (likelihood != D_INF) {
        if (likelihood < previous_likelihood) {
          nb_likelihood_decrease++;
        }
        else {
          nb_likelihood_decrease = 0;
        }

        // reestimation des probabilites initiales

	reestimation(hsmarkov->nb_state , chain_reestim->initial ,
		     hsmarkov->initial , MIN_PROBABILITY , false);

        // reestimation des probabilites de transition

        for (i = 0;i < hsmarkov->nb_state;i++) {
          reestimation(hsmarkov->nb_state , chain_reestim->transition[i] ,
                       hsmarkov->transition[i] , MIN_PROBABILITY , false);
        }

        // reestimation des lois d'occupation des etats

        min_likelihood = 0.;

        for (i = 0;i < hsmarkov->nb_state;i++) {
          if (hsmarkov->state_subtype[i] == SEMI_MARKOVIAN) {
            occupancy = hsmarkov->nonparametric_process[0]->sojourn_time[i];

            complete_run[i]->nb_value_computation();
            complete_run[i]->offset_computation();
            complete_run[i]->nb_element_computation();

#           ifdef DEBUG
            cout << "\n" << STAT_label[STATL_STATE] << " " << i << " ";

            complete_run[i]->max_computation();
            complete_run[i]->mean_computation();
            complete_run[i]->variance_computation();

            complete_run[i]->print(cout);
#           endif

	    if ((iter > STOCHASTIC_EXPLORATION_NB_ITER) && (estimator == COMPLETE_LIKELIHOOD)) {
              final_run[i]->nb_value_computation();
              final_run[i]->offset_computation();
              final_run[i]->nb_element_computation();

	      if (final_run[i]->nb_element > 0.) {
		complete_run[i]->state_occupancy_estimation(final_run[i] , complete_run[i] , occupancy_survivor ,
							    censored_occupancy_survivor , false);
	      }

	    }

            complete_run[i]->max_computation();
            complete_run[i]->mean_computation();
            complete_run[i]->variance_computation();

            if (iter <= EXPLORATION_NB_ITER) {
              occupancy_likelihood = complete_run[i]->parametric_estimation(occupancy , 1 , true , OCCUPANCY_THRESHOLD);
            }
            else {
              occupancy_likelihood = complete_run[i]->type_parametric_estimation(occupancy , 1 , true , OCCUPANCY_THRESHOLD);
	    }

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
            cout << STAT_word[STATW_STATE] << " " << i << " " << STAT_word[STATW_OCCUPANCY_DISTRIBUTION] << endl;
            occupancy->ascii_print(cout);
#           endif

          }
        }
      }


      //#     ifdef MESSAGE
      os << STAT_label[STATL_ITERATION] << " " << iter << "   "
         << SEQ_label[SEQL_OBSERVED_SEQUENCES_LIKELIHOOD] << ": " << likelihood
         << "   (" << nb_state_sequence << ")" << endl;
      //#     endif

#     ifdef DEBUG
      if (iter % 5 == 0) {
        cout << *hsmarkov;
      }
#     endif

    }
    while (iter<100);
//     ((likelihood != D_INF) && ((iter < STOCHASTIC_EXPLORATION_NB_ITER + 2) ||
// 			       ((nb_iter == I_DEFAULT) && (iter < MARKOV_SWITCHING_NB_ITER) &&
// 				(((likelihood - previous_likelihood) / -likelihood > MARKOV_SWITCHING_LIKELIHOOD_DIFF) ||
// 				 (min_likelihood == D_INF) || (nb_likelihood_decrease == 1))) ||
// 			       ((nb_iter != I_DEFAULT) && (iter < nb_iter))));


    if (likelihood != D_INF) {

#     ifdef MESSAGE
      os << "\n" << iter << " " << STAT_label[STATL_ITERATIONS] << endl;

#     endif

      // reestimation des probabilites initiales

      if (hsmarkov->type == 'o') {
        reestimation(hsmarkov->nb_state , chain_reestim->initial ,
                     hsmarkov->initial , MIN_PROBABILITY , true);
      }

      // reestimation des probabilites de transition

      for (i = 0;i < hsmarkov->nb_state;i++) {
        reestimation(hsmarkov->nb_state , chain_reestim->transition[i] ,
                     hsmarkov->transition[i] , MIN_PROBABILITY , true);
      }


      for (i = 0;i < hsmarkov->nb_state;i++) {
        if ((hsmarkov->state_subtype[i] == SEMI_MARKOVIAN) &&
            (hsmarkov->nonparametric_process[0]->sojourn_time[i]->mean == 1.)) {
          hsmarkov->state_subtype[i] = MARKOVIAN;
          delete hsmarkov->nonparametric_process[0]->sojourn_time[i];
          hsmarkov->nonparametric_process[0]->sojourn_time[i] = 0;
          delete hsmarkov->forward[i];
          hsmarkov->forward[i] = 0;
        }
      }

      // reestimation des lois d'observation non-parametriques

      for (j = 0;j < hsmarkov->nb_state;j++) {
	
	hsmarkov->sw_process[1]->observation[j]->parametric_variance_computation();
	hsmarkov->sw_process[1]->observation[j]->min_value_computation();
	hsmarkov->sw_process[1]->observation[j]->max_value_computation();
	hsmarkov->sw_process[1]->observation[j]->nb_value_computation();
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

   for (i = 0;i < hsmarkov->nb_state;i++) {
     delete complete_run[i];
   }
   delete [] complete_run;

   for (i = 0;i < hsmarkov->nb_state;i++) {
     delete final_run[i];
   }
   delete [] final_run;

   if (hsmarkov->type == 'e') {
     for (i = 0;i < hsmarkov->nb_state;i++) {
       delete initial_run[i];
     }
     delete [] initial_run;

     for (i = 0;i < hsmarkov->nb_state;i++) {
       delete single_run[i];
     }
     delete [] single_run;
   }

   delete [] occupancy_nb_value;

   if (estimator != PARTIAL_LIKELIHOOD) {
     delete [] occupancy_survivor;
     delete [] censored_occupancy_survivor;
   }

   delete [] poutput;


    if (likelihood == D_INF) {
      delete hsmarkov;
      hsmarkov = 0;
      error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
    }

    else {
      if  (state_sequence) {
        hsmarkov->sw_markov_data = new Markov_switching_data(*this , 0);
        seq = hsmarkov->sw_markov_data;
        seq->type[0] = STATE;

	if ((hsmarkov->sw_process[1]) && (seq->characteristics[1])) {
	  delete seq->characteristics[1];
	  seq->characteristics[1] = 0;
	}

	hsmarkov->create_cumul();

	seq->posterior_probability = new double[seq->nb_sequence];
	if(!output_R){
	  seq->likelihood = hsmarkov->viterbi(*seq, seq->posterior_probability,true, false);
	}
	else {
	  seq->likelihood = hsmarkov->viterbi(*seq, seq->posterior_probability, false, true);
	}

 	hsmarkov->remove_cumul();

        seq->max_value[0] = hsmarkov->nb_state - 1;

        seq->build_marginal_histogram(0);
        seq->build_characteristic(0);

        seq->build_transition_count(*hsmarkov);
        
#       ifdef MESSAGE
        cout << "\n" << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << ": " << seq->likelihood <<endl;
#       endif

	// calcul des lois d'occupation des etats

        for (i = 0;i < hsmarkov->nb_state;i++) {
          if (hsmarkov->state_subtype[i] == SEMI_MARKOVIAN) {
            hsmarkov->nonparametric_process[0]->sojourn_time[i]->computation((seq->characteristics[0] ? 
									      seq->characteristics[0]->sojourn_time[i]->nb_value : 1) , 
									     OCCUPANCY_THRESHOLD);
            if (hsmarkov->state_type[i] == 'r') {
	      hsmarkov->forward[i]->copy(*(hsmarkov->nonparametric_process[0]->sojourn_time[i]));
              hsmarkov->forward[i]->computation(*(hsmarkov->nonparametric_process[0]->sojourn_time[i]));
            }
          }
        }
      }
      else {
	for (i = 0;i < hsmarkov->nb_state;i++) {
	  if ((hsmarkov->state_subtype[i] == SEMI_MARKOVIAN) && (hsmarkov->state_type[i] == 'r')) {
	    hsmarkov->forward[i]->copy(*(hsmarkov->nonparametric_process[0]->sojourn_time[i]));
	    hsmarkov->forward[i]->computation(*(hsmarkov->nonparametric_process[0]->sojourn_time[i]));
	  }
	}
	
        hsmarkov->sw_markov_data = new Markov_switching_data(*this , false);
        seq = hsmarkov->sw_markov_data;
        seq->state_variable_init(REAL_VALUE);

	if ((hsmarkov->sw_process[1]) && (seq->characteristics[0])) {
	  delete seq->characteristics[0];
	  seq->characteristics[0] = 0;
	}
      }

      // calcul de la vraisemblance et des lois caracteristiques du modele

      seq->hidden_likelihood = hsmarkov->likelihood_computation(*this , seq->posterior_probability, I_DEFAULT);

      hsmarkov->component_computation();
      hsmarkov->characteristic_computation(*seq , counting_flag , I_DEFAULT , false);

#     ifdef MESSAGE
      if  ((state_sequence) && (seq->nb_sequence <= POSTERIOR_PROBABILITY_NB_SEQUENCE)) {
        os << "\n" << SEQ_label[SEQL_POSTERIOR_STATE_SEQUENCE_PROBABILITY] << endl;
        for (i = 0;i < seq->nb_sequence;i++) {
          os << SEQ_label[SEQL_SEQUENCE] << " " << seq->identifier[i] << ": "
             << seq->posterior_probability[i];

          if (hsmarkov->nb_component == hsmarkov->nb_state) {
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
      }
#     endif

    }
  }

  return hsmarkov;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'un semi-Markov switching linear 
 *  model a partir d'un echantillon de sequences par l'algorithme 
 *  SEM/MCEM.
 *
 *  arguments : reference sur un objet Format_error, stream, type de processus
 *              ('o' : ordinaire, 'e' : en equilibre), nombre d'etats de la chaine de Markov,
 *              flag sur la nature de la chaine de Markov, variance résiduelle, paramètres
 *              de regression, nombre de covariables, presence/absence de constante,
 *              ordre, parametres pour nombre de sequences d'etats simulees, 
 *              type Markovien, variance résiduelle commune à tous les états,
 *              sortie de type R, flags sur le calcul des lois de comptage et sur le calcul
 *              des sequences d'etats optimales, probabilite de rester dans un etat initiale,
 *              nombre d'iterations.
 *
 *--------------------------------------------------------------*/

Markov_switching* Switching_sequence::markov_switching_stochastic_estimation(Format_error &error , ostream &os ,
									     char type , int nb_state , bool left_right ,
									     double *residual_variance, double **regression, 
									     int nb_covariable, int constant, int order,
									     int min_nb_state_sequence , int max_nb_state_sequence ,
									     bool markov, bool VarCommune, bool output_R,
									     double parameter , bool counting_flag ,
									     bool state_sequence , double occupancy_mean, 
									     double self_transition ,
									     int nb_iter) const
{
  bool status = true;
  register int i,j;
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

    if(!constant) {
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
      cont_param[i] = new Continuous_parametric(NB_VALUE, nb_covariable, 0, 
						residual_variance[i], 0, regression[i], NO_RANDOM); 
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

    hsw_markov = new Markov_switching(*ihsw_markov, true);
    hsw_markov = markov_switching_stochastic_estimation(error , os , *ihsw_markov , min_nb_state_sequence ,
							max_nb_state_sequence , VarCommune, false, parameter ,
							COMPLETE_LIKELIHOOD, counting_flag , state_sequence , nb_iter);
    delete ihsw_markov;
  }

  return hsw_markov;
}
