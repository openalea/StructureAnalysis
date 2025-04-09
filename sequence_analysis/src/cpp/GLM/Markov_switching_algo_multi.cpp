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
#include <gsl/gsl_rng.h>
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
#include "switching_sequence.h"
#include "switching_process.h"
#include "markov_switching.h"

#include "stat_tool/distribution_reestimation.h"

using namespace std;

extern void convert_matrix_int (int **essai, int nb_row, int nb_col, gsl_matrix* convert_essai);
extern void convert_matrix_diag_int (int *essai, int nb_row, gsl_matrix* convert_essai); 
extern void convert_vector_int (int *essai, int length, gsl_vector* essai_convert);
extern void convert_array2D_int (gsl_matrix *essai, int nb_row, int nb_col, int **convert_essai);
extern void convert_array1D_diag_int (gsl_matrix *essai, int length, int *essai_convert); 
extern void convert_array1D_int (gsl_vector *essai, int length, int *essai_convert);

extern void convert_matrix_double (double **essai, int nb_row, int nb_col, gsl_matrix* convert_essai);
extern void convert_matrix_diag_double (gsl_vector *essai, int nb_row, gsl_matrix* convert_essai); 
extern void convert_vector_double (double *essai, int length, gsl_vector* essai_convert);
extern void convert_array2D_double (gsl_matrix *essai, int nb_row, int nb_col, double **convert_essai);
extern void convert_array1D_diag_double (gsl_matrix *essai, int length, double *essai_convert);
extern void convert_array1D_double (gsl_vector *essai, int length, double *essai_convert);

extern void cumul_computation(int nb_value, const double *pmass, double *pcumul);
extern int cumul_method(int nb_value, const double *cumul, double scale=1.);
extern void log_computation(int nb_value, const double *pmass, double *plog);
extern char* label(const char *file_name);


/*--------------------------------------------------------------*
 *  Dans renewal_algorithms.cpp
 *
 *  Calcul de la moyenne d'une loi par bissection d'intervalle.
 *
 *  arguments : pointeurs sur les quantites de reestimation de la loi et
 *              de la loi biaisee par la longueur.
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



//--------------------------------------------------------------
// Calcul de la vraisemblance de sequences pour une semi-chaine de Markov.
//--------------------------------------------------------------
double Markov_switching::likelihood_computation(const Switching_sequence &isw_seq , int index) const
{
  register int i , j , k, m, l;
  int nb_value , occupancy, *pstate , **poutput;
  double likelihood = 0. , proba;


  // verification de la compatibilite entre le modele et les donnees

  if (nb_output_process +1 == isw_seq.nb_variable) {
    for (i = 0;i <= nb_output_process;i++) {
      nb_value = sw_process[i]->nb_value;
      if ((isw_seq.marginal[i]) && (nb_value < isw_seq.marginal[i]->nb_value)) {
	likelihood = D_INF;
	break;
      }
    }
  }
  
  else {
    likelihood = D_INF;
  }
  
  if (likelihood != D_INF) {
    if (nb_output_process > 0) {
      poutput = new int*[nb_output_process];
    }

    for (i = 0;i < isw_seq.nb_sequence;i++) {
      if ((index == I_DEFAULT) || (index == i)) {
        pstate = isw_seq.int_sequence[i][0];

        proba = initial[*pstate];
        if (proba > 0.) {
          likelihood += log(proba);
        }
        else {
          likelihood = D_INF;
          break;
        }

        for (j = 0;j < nb_output_process;j++) {
          poutput[j] = isw_seq.int_sequence[i][j + 1];
	}

	j = 0;
	do {
	  if (j > 0){
	    pstate++;

	    proba = transition[*(pstate-1)][*pstate];
	    if (proba > 0.){
	      likelihood += log(proba);
	    }
	    else {
	      likelihood = D_INF;
	      break;
	    }
	  }
	   
	  if (transition[*pstate][*pstate] < 1.) {
	    j++;
	    occupancy = 1;
	    
	    if (state_subtype[*pstate] == SEMI_MARKOVIAN) {
              while ((j < isw_seq.length[i]) && (*(pstate + 1) == *pstate)) {
                j++;
                occupancy++;
                pstate++;
              }
	      
              proba = 0.;
              if ((type == 'e') && (j == occupancy)) {
                if (occupancy < forward[*pstate]->nb_value) {
                  if (j < isw_seq.length[i]) {
                    proba = forward[*pstate]->mass[occupancy];
                  }
                  else {
                    proba = (1. - forward[*pstate]->cumul[occupancy - 1]);
                  }
                }
              }

              else {
                if (occupancy < nonparametric_process[0]->sojourn_time[*pstate]->nb_value) {
                  if (j < isw_seq.length[i]) {
                    proba = nonparametric_process[0]->sojourn_time[*pstate]->mass[occupancy];
                  }
                  else {
                    proba = (1. - nonparametric_process[0]->sojourn_time[*pstate]->cumul[occupancy - 1]);
                  }
                }
              }

              if (proba > 0.) {
                likelihood += log(proba);
              }
              else {
                likelihood = D_INF;
                break;
              }
            }
          }

          else {
            occupancy = isw_seq.length[i] - j;
            j += occupancy;
          }

	  for (m = 0; m < occupancy; m++){
	    for (k = 0;k < nb_output_process;k++) {
	      proba = sw_process[k + 1]->observation[*pstate]->mass_computation(*++poutput[k], isw_seq.covar[i][j]);
	      
	      if (proba > 0.) {
		likelihood += log(proba);
	      }
	      else {
		likelihood = D_INF;
		break;
	      }
	    }
	  
	    if (likelihood == D_INF) {
	      break;
	    }
	  }

	  if (likelihood == D_INF) {
	    break;
	  }
	}
	while (j < isw_seq.length[i]);


        if (likelihood == D_INF) {
          break;
        }
      }
    }

    if (nb_output_process > 0) {
      delete [] poutput;
    }
  }
  
  return likelihood;
}

//--------------------------------------------------------------
// Calcul de la vraisemblance de sequences pour une chaine de Markov.
//
// A REVOIR , probleme au niveau du calcul de la vraisemblance 
// de l'histogramme des observations
//--------------------------------------------------------------
double Markov_switching::likelihood_computation(const Markov_switching_data &isw_markov_data) const
{
  register int i , j;
  int nb_value;
  double buff , likelihood = 0.;
  Histogram **initial_run, **final_run, **single_run;



  // verification de la compatibilite entre le modele et les donnees

  if (nb_output_process +1  == isw_markov_data.nb_variable) { 
    for (i = 0;i <= nb_output_process;i++) {
      nb_value = sw_process[i]->nb_value;
      
      if (nb_value < isw_markov_data.marginal[i]->nb_value) {
        likelihood = D_INF;
        break;
      }
    }
  }
  
  else {
    likelihood = D_INF;
  }
  
  if (likelihood != D_INF) {
    likelihood = Chain::likelihood_computation(*(isw_markov_data.chain_data));

    for (i = 0;i < nb_state;i++) {
      if (state_subtype[i] == SEMI_MARKOVIAN) {
	buff = nonparametric_process[0]->sojourn_time[i]->likelihood_computation(*(isw_markov_data.characteristics[0]->sojourn_time[i]));
	
	if (buff != D_INF) {
	  likelihood += buff;
	}
	else {
	  likelihood = D_INF;
	  break;
	}
	
	buff = nonparametric_process[0]->sojourn_time[i]
	  ->survivor_likelihood_computation(*(isw_markov_data.characteristics[0]->final_run[i]));
	
	if (buff != D_INF) {
	  likelihood += buff;
	}
	else {
	  likelihood = D_INF;
	}

	if (likelihood == D_INF) {
	  break;
	}
      }
    }

    
    if (likelihood != D_INF) {
      for (i = 1;i <= nb_output_process;i++) {
	for (j = 0;j < nb_state;j++) {
	  buff = sw_process[i]->observation[j]->likelihood_computation(*(isw_markov_data.observation[i][j]));

	  if (buff != D_INF) {
	    likelihood += buff;
	  }
	  else { 
	    likelihood = D_INF;
	    break;
	  }
	}
	
	if (likelihood == D_INF) {
	  break;
	}
      }
    }
  }
  
  
  return likelihood;
}

//--------------------------------------------------------------
// Calcul de la vraisemblance de sequences pour une chaine de Markov
// cachee par l'algorithme forward.
//--------------------------------------------------------------
double Markov_switching::likelihood_computation(const Switching_sequence &isw_seq ,
						double *posterior_probability,
                                                int index, bool hidden) const
{
  register int i , j , k , l, m, p;
  int nb_value , length,**poutput;
  double likelihood = 0. , seq_likelihood, obs_product, **observation, *state_norm, *forward1 , **state_in , *norm;
  Parametric *occupancy;
  
  if (!hidden){
    likelihood = (*this).likelihood_computation(isw_seq, index);
  }
  else {
    if (likelihood != D_INF) {
      
      // initialisations
      
      length = (index == I_DEFAULT ? isw_seq.max_length : isw_seq.length[index]);

      observation = new double*[length];
      for (i = 0;i < length;i++) {
	observation[i] = new double[nb_state];
      }
      
      norm = new double[length];
      state_norm = new double[nb_state];
      forward1 = new double[nb_state];
      
      state_in = new double*[length - 1];
      for (i = 0;i < length - 1;i++) {
	state_in[i] = new double[nb_state];
      }

      poutput = new int*[isw_seq.nb_variable];
      
      for (i = 0;i < isw_seq.nb_sequence;i++) {
	if ((index == I_DEFAULT) || (index == i)) {
	  for (j = 0;j < isw_seq.nb_variable;j++) {
	    poutput[j] = isw_seq.int_sequence[i][j];
	  }
	  seq_likelihood = 0.;
	  
	  for (j = 0; j < isw_seq.length[i];j++){
	    norm[j] = 0.;

	    for (k = 0;k < nb_state;k++) {

	      // calcul des probabilites d'observation
	      observation[j][k] = 1.; 
	    
	      for (m = 0;m < nb_output_process;m++) {
		observation[j][k] *= sw_process[m + 1]->observation[k]->mass_computation(*poutput[m], isw_seq.covar[i][j]);
	      }
	

	      switch (state_subtype[k]) {
		
		// cas etat semi-markovien
		
	      case SEMI_MARKOVIAN : {
		if (j == 0) {
		  state_norm[k] = initial[k];
		}
		else {
		  state_norm[k] += state_in[j - 1][k] - forward1[k];
		}
		state_norm[k] *= observation[j][k];
		
		norm[j] += state_norm[k];
		break;
	      }
		
		// cas etat markovien
		
	      case MARKOVIAN : {
		if (j == 0) {
		  forward1[k] = initial[k];
		}
		else {
		  forward1[k] = state_in[j - 1][k];
		}
		forward1[k] *= observation[j][k];
		
		norm[j] += forward1[k];
		break;
	      }
	      }
	    }
  

	    if (norm[j] > 0.) {
	      for (k = 0;k < nb_state;k++) {
		switch (state_subtype[k]) {
		case SEMI_MARKOVIAN :
		  state_norm[k] /= norm[j];
		  break;
		case MARKOVIAN :
		  forward1[k] /= norm[j];
		  break;
		}
	      }
	      
	      seq_likelihood += log(norm[j]);
	    }
	    
	    else {
	      seq_likelihood = D_INF;
	      break;
	    }
	    

	    for (k = 0;k < nb_state;k++) {
	      
	      // cas etat semi-markovien
	      
	      if (state_subtype[k] == SEMI_MARKOVIAN) {
		occupancy = nonparametric_process[0]->sojourn_time[k];
		obs_product = 1.;
		forward1[k] = 0.;
		
		if (j < isw_seq.length[i] - 1) {
		  for (m = 1;m <= MIN(j + 1 , occupancy->nb_value - 1);m++) {
		    obs_product *= observation[j - m + 1][k] / norm[j - m + 1];
		    if (obs_product == 0.) {
		      break;
		    }
		    
		    if (m < j + 1) {
		      forward1[k] += obs_product * occupancy->mass[m] * state_in[j - m][k];
		    }
		    
		    else {
		      forward1[k] += obs_product * occupancy->mass[m] * initial[k];
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
		      forward1[k] += obs_product * (1. - occupancy->cumul[m - 1]) * state_in[j - m][k];
		    }
		    
		    else {
                      forward1[k] += obs_product * (1. - occupancy->cumul[m - 1]) * initial[k];
		    }
		  }
		}
	      }

	    }
	    
	    if (j < isw_seq.length[i] - 1) {
	      for (k = 0;k < nb_state;k++) {
		state_in[j][k] = 0.;
 		for (m = 0;m < nb_state;m++) {
		  state_in[j][k] += transition[m][k] * forward1[m];
		}
	      }
	    }

	    for (k = 0;k < isw_seq.nb_variable;k++) {
	      poutput[k]++;
	    }
 	  }

	  if (seq_likelihood != D_INF) {
	    likelihood += seq_likelihood;
	    if (posterior_probability) {
	      posterior_probability[i] = exp(posterior_probability[i] - seq_likelihood);
	    }
	  }
	  
	  else {
	    likelihood = D_INF;
	    break;
	  }
	}
      }

      //       for (i = 0;i < length;i++) {
      // 	delete [] observation[i];
      //       }
      //       delete [] observation;
      
      //       delete [] norm;
      //       delete [] state_norm;
      //       delete [] forward1;
      
      //       for (i = 0;i < length - 1;i++) {
      // 	delete [] state_in[i];
      //       }
      //       delete [] state_in;
      
      //       delete [] poutput;
    }
  } 
 
  return likelihood;
}

/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une chaine de Markov switching
 *  a partir d'un echantillon de sequences par l'algorithme EM.
 *
 *  arguments : reference sur un objet Format_error, stream, chaine de Markov switching initiale,
 *              flags sur le calcul des lois de comptage et sur le calcul des sequences
 *              d'etats optimales, nombre d'iterations.
 *
 *--------------------------------------------------------------*/
Markov_switching* Switching_sequence::markov_switching_estimation(Format_error &error , ostream &os ,
								  const Markov_switching &isw_markov ,
								  int estimator, bool output_R,
								  bool counting_flag , bool state_sequence ,
								  int nb_iter, int mean_computation) const
{
  bool status;
  register int i , j , k , m, n;
  int max_nb_value , iter, nb_likelihood_decrease, offset, nb_value, *occupancy_nb_value,
    *censored_occupancy_nb_value;
  int **poutput;
  double likelihood = D_INF , previous_likelihood , occupancy_likelihood, observation_likelihood , min_likelihood ,
    obs_product, buff, sum, occupancy_mean, **observation, *norm, *state_norm, **forward, **state_in, ***backward, **backward1,
    *auxiliary, *reestim, *ofrequency, *lfrequency, *occupancy_survivor, *censored_occupancy_survivor;
  
  double *temp;

  Chain_reestimation<double> *chain_reestim;
  Parametric *occupancy;
  Reestimation<double> **occupancy_reestim, **length_bias_reestim, **censored_occupancy_reestim;
  Histogram *hoccupancy;
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
    hsw_markov = new Markov_switching(isw_markov, false,I_DEFAULT);// (int)(max_length * SAMPLE_NB_VALUE_COEFF));

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
    for (i = 0;i < hsw_markov->nb_output_process;i++) {
      if ((hsw_markov->sw_process[i + 1]) && (max_nb_value < marginal[i]->nb_value)) {
        max_nb_value = marginal[i]->nb_value;
      }
    }

    poutput = new int*[nb_variable];
    
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
	  poutput[j] = int_sequence[i][j];
	}
  
	// recurrence "forward"
	
	for (j = 0;j < length[i];j++) {
          norm[j] = 0.;

          for (k = 0;k < hsw_markov->nb_state;k++) {

            // calcul des probabilites d'observation

            observation[j][k] = 1.;
            for (m = 0;m < hsw_markov->nb_output_process;m++) {
	      observation[j][k] *= hsw_markov->sw_process[m + 1]->observation[k]->mass_computation(*poutput[m], covar[i][j]);
	    }

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

              //accumulation des quantites de reestimation des probabilites de transition

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

      // estimation des parametres de regression pour chacun des états 
      int s;
      gsl_permutation *p = gsl_permutation_alloc(nb_covariable);
      gsl_matrix *trans = gsl_matrix_alloc(nb_covariable, nb_covariable);
      gsl_matrix *beta = gsl_matrix_alloc(hsw_markov->nb_state, nb_covariable);
      gsl_vector *beta_k = gsl_vector_alloc(nb_covariable);
      gsl_vector *betaprime_k = gsl_vector_alloc(nb_covariable);
      gsl_matrix *inv_temp_add = gsl_matrix_alloc(nb_covariable, nb_covariable);
      gsl_matrix *temp_add = gsl_matrix_calloc(nb_covariable, nb_covariable);
      gsl_matrix *temp_mul2 = gsl_matrix_alloc(nb_covariable, nb_covariable);
      gsl_vector *add_vect = gsl_vector_calloc(nb_covariable); 
      gsl_vector *temp_mul3 = gsl_vector_alloc(nb_covariable);

      for (k = 0; k < hsw_markov->nb_state; k++) {

	for (j = 0; j < nb_covariable; j++) {
	  gsl_vector_set(beta_k, j, hsw_markov->sw_process[1]->observation[k]->regression[j] );
	}
	gsl_vector_set_zero(betaprime_k);
	gsl_matrix_set_zero(temp_add);
	gsl_matrix_set_zero(temp_mul2);
	gsl_vector_set_zero(add_vect);
	gsl_vector_set_zero(temp_mul3);

	// pour l'instant, cet algorithme ne marchera que pour une seule variable j
	for (j = 0 ; j < hsw_markov->nb_output_process; j++) {
	  
	  double sum = 0;
	  double sigma_k = 0;	
	  double *beta_kk;


	  for (i = 0; i < nb_sequence; i++) {
	    
	    gsl_matrix *X_i = gsl_matrix_calloc(length[i], nb_covariable);
	    gsl_matrix *L_ik = gsl_matrix_calloc(length[i], length[i]);
	    gsl_matrix *W_ik = gsl_matrix_calloc(length[i], length[i]);
	    gsl_matrix *temp_mat = gsl_matrix_calloc(length[i], hsw_markov->nb_state);
	    gsl_matrix *temp_mul = gsl_matrix_calloc(nb_covariable, length[i]);
	    gsl_matrix *temp_mul1 = gsl_matrix_calloc(nb_covariable, length[i]);
	    gsl_vector *Y_i = gsl_vector_calloc(length[i]);
	    gsl_vector *mu_k = gsl_vector_calloc(length[i]);
	    gsl_vector *temp_vect = gsl_vector_calloc(length[i]);
	    
	    convert_matrix_double(covar[i], length[i], nb_covariable, X_i);
	    
	    convert_matrix_double(backward[i], length[i], hsw_markov->nb_state, temp_mat);
	    
	    convert_vector_int(int_sequence[i][j], length[i], Y_i);
	    
	    gsl_matrix_get_col(temp_vect, temp_mat, k);
	    convert_matrix_diag_double(temp_vect, length[i], L_ik);

	    for (m = 0; m < length[i]; m++) {
	      gsl_matrix_set(W_ik, m, m, hsw_markov->sw_process[j + 1]->observation[k]->variance_computation(covar[i][m])); 
	    }

	    for (m = 0; m < length[i]; m++) {
	      gsl_vector_set(mu_k, m, hsw_markov->sw_process[j + 1]->observation[k]->mean_computation(covar[i][m])); 
	    }


	    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., X_i, L_ik, 0., temp_mul); // t_XL


	    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., temp_mul, W_ik, 0., temp_mul1); // t_XLW	    
	  


	    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., temp_mul1, X_i, 0., temp_mul2); // t_XLWX
	    gsl_matrix_add(temp_add, temp_mul2);// somme des t_XLWX
	    
	    gsl_vector_sub(Y_i, mu_k); //(Y-mu)

	    // 	    if(i==1){
	    // 	      for(m = 0; m < length[i]; m++){
	    // 		cout<<gsl_vector_get(mu_k,m)<<"  "<< gsl_vector_get(Y_i, m)<< endl;
	    // 	      }
	    // 	    }

	    gsl_blas_dgemv(CblasNoTrans, 1., temp_mul, Y_i, 0., temp_mul3);// t_XL(Y-mu)
	    gsl_vector_add(add_vect, temp_mul3);// somme des t_XL(Y-mu)

	    // 	    for (m = 0; m < nb_covariable; m++){
	    // 	      cout << gsl_vector_get(add_vect, m) <<"   ";
	    // 	    }
	    // 	    cout<<endl<<endl;


	    gsl_vector_free(temp_vect);
	    gsl_vector_free(mu_k);
	    gsl_vector_free(Y_i);
	    gsl_matrix_free(temp_mul1);
	    gsl_matrix_free(temp_mul);
	    gsl_matrix_free(temp_mat);
	    gsl_matrix_free(L_ik);
	    gsl_matrix_free(X_i);
	    
	  }
	  
	  gsl_matrix_memcpy(trans, temp_add);
	  gsl_linalg_LU_decomp(trans, p, &s);
	  gsl_linalg_LU_invert(trans, p, inv_temp_add);

	  gsl_blas_dgemv(CblasNoTrans, 1., inv_temp_add, add_vect, 0., betaprime_k);
	  gsl_vector_add(beta_k, betaprime_k);
	  gsl_matrix_set_row(beta, k, beta_k);

  
	  // mise à jour des lois d'observations
	  
	  beta_kk = new double[nb_covariable];
	  convert_array1D_double (beta_k, nb_covariable, beta_kk);
	  hsw_markov->sw_process[j + 1]->observation[k]->Set_regression(beta_kk);


#   ifdef DEBUG
	  cout<<"beta pour l'état "<<k<<" :  "<<endl;
	  for(i = 0; i<nb_covariable; i++){
	    cout<<beta_kk[i]<<"   ";
	  }
	  cout<<endl;
#endif


	  delete [] beta_kk;

	}

      }

      gsl_vector_free(beta_k);
      gsl_vector_free(betaprime_k);
      gsl_vector_free(add_vect);
      gsl_vector_free(temp_mul3);
      gsl_matrix_free(temp_add);
      gsl_matrix_free(temp_mul2);
      gsl_matrix_free(inv_temp_add);
      gsl_matrix_free(beta);
      gsl_matrix_free(trans);
      gsl_permutation_free(p);


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


//             if(iter==1){
//       	os<< likelihood-previous_likelihood<<", "<<endl;
//             }

      
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
    while ((likelihood != D_INF) && ((likelihood - previous_likelihood) / -likelihood > MARKOV_SWITCHING_LIKELIHOOD_DIFF)
	   && ((nb_iter == I_DEFAULT) && (iter < MARKOV_SWITCHING_NB_ITER) || ((nb_iter != I_DEFAULT) && (iter < nb_iter))));

  
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
	
      for (i = 0;i < hsw_markov->nb_output_process;i++) {
	for (j = 0;j < hsw_markov->nb_state;j++) {
	  hsw_markov->sw_process[i + 1]->observation[j]->parametric_mean_computation();
	  hsw_markov->sw_process[i + 1]->observation[j]->parametric_variance_computation();
	  // hsw_markov->sw_process[i + 1]->observation[j]->min_value_computation();
	  // hsw_markov->sw_process[i + 1]->observation[j]->max_value_computation();
	  hsw_markov->sw_process[i + 1]->observation[j]->nb_value_computation();
	  // hsw_markov->sw_process[i + 1]->observation[j]->mass_computation();
	  hsw_markov->sw_process[i + 1]->observation[j]->cumul_computation();

   
	  // 	    if (observation_likelihood == D_INF) {
	  // 	      min_likelihood = D_INF;
	  // 	    }
	  // 	    else {
	  // 	      hsw_markov->sw_process[i + 1]->observation[j]->computation(OBSERVATION_THRESHOLD);
	  // 	    }
	}
      }
	
    }


#   ifdef DEBUG 
    cout<<"AFFICHAGE DES PARAMETRES DE REGRESSION: "<<endl;
    for(i = 0; i< hsw_markov->nb_output_process; i++){
      for (j = 0; j< hsw_markov->nb_state; j++) {
	for(k=0; k<nb_covariable; k++) {
	  cout << hsw_markov->sw_process[i + 1]->observation[j]->regression[k]<< "   ";
	}
	cout<<endl;
      }
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
	
	for (i = 1;i <= hsw_markov->nb_output_process;i++) {
	  if ((hsw_markov->sw_process[i]) && (hsw_markov_data->characteristics[i])) {
	    delete hsw_markov_data->characteristics[i];
	    hsw_markov_data->characteristics[i] = 0;
	  }
	}

	hsw_markov->create_cumul();	
	
	//hsw_markov->log_computation(); // PROBLEME de memoire

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
							    true, false);// attention a cause de log_computation()
	}
	else{
	  hsw_markov_data->likelihood = hsw_markov->viterbi(*hsw_markov_data, hsw_markov_data->posterior_probability, 
							    false, true);// attention a cause de log_computation()
	}

	hsw_markov->remove_cumul();

	hsw_markov_data->max_value[0] = hsw_markov->nb_state - 1;
	hsw_markov_data->build_marginal_histogram(0);
	hsw_markov_data->build_characteristic(0);

	hsw_markov_data->build_transition_count(*hsw_markov);
	hsw_markov_data->build_observation_histogram(); 


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
	//   << " | " << hsw_markov->Markov_switching::likelihood_computation(*hsw_markov_data) << endl;  // A REVOIR
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
	hsw_markov_data->state_variable_init(INT_VALUE); 

	for (i = 1;i <= hsw_markov->nb_output_process;i++) {
	  if ((hsw_markov->sw_process[i]) && (hsw_markov_data->characteristics[i - 1])) {
	    delete hsw_markov_data->characteristics[i - 1];
	    hsw_markov_data->characteristics[i - 1] = 0;
	  }
	}
      }
    
      // calcul de la vraisemblance et des lois caracteristiques du modele

      hsw_markov_data->hidden_likelihood = hsw_markov->likelihood_computation(*this, hsw_markov_data->posterior_probability,
									      I_DEFAULT, true);
      
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





//--------------------------------------------------------------
// Calcul des sequences d'etats optimale par l'algorithme de Viterbi
//
// attention, pas de log_computation, je fais les logs
// directement dans l'algo
//--------------------------------------------------------------
double Markov_switching::viterbi(const Markov_switching_data &seq , double *posterior_probability,
				 bool fag, bool output_R, int index) const
{
  register int i , j , k , l, m, r;
  int tmp;
  int length , *pstate , **poutput , **input_state, **optimal_state, **optimal_occupancy;
  double likelihood = 0. , obs_product, buff , forward_max , **observation, *forward1 , **state_in;
  Parametric *occupancy;

  // initialisations
  length = (index == I_DEFAULT ? seq.max_length : seq.length[index]);

  observation = new double*[length];
  for (i = 0;i < length;i++) {
    observation[i] = new double[nb_state];
  }

  forward1 = new double[nb_state];

  state_in = new double*[length - 1];
  for (i = 0;i < length - 1;i++) {
    state_in[i] = new double[nb_state];
  }

  input_state = new int*[length - 1];
  for (i = 0;i < length - 1;i++) {
    input_state[i] = new int[nb_state];
  }

  optimal_state = new int*[length];
  for (i = 0;i < length;i++) {
    optimal_state[i] = new int[nb_state];
  }

  optimal_occupancy = new int*[length];
  for (i = 0;i < length;i++) {
    optimal_occupancy[i] = new int[nb_state];
  }

  poutput = new int*[nb_output_process];

  for (i = 0;i < seq.nb_sequence;i++) {
    if ((index == I_DEFAULT) || (index == i)) {
      for (j = 0;j < nb_output_process;j++) {
        poutput[j] = seq.int_sequence[i][j + 1];
      }

      for (j = 0;j < seq.length[i];j++) {
        for (k = 0;k < nb_state;k++) {
	  
          // calcul des probabilites d'observation
	  
          observation[j][k] = 0.;
	  
          for (m = 0; m < nb_output_process;m++) {
	    buff = log(sw_process[m + 1]->observation[k]->mass_computation(*poutput[m], seq.covar[i][j]));
	    
	    if (buff == D_INF) {
              observation[j][k] = D_INF;
              break;
            }
            else {
              observation[j][k] += buff;
            }
          }
	  
	  switch (state_subtype[k]) {
	    
	    // cas etat semi-markovien
	    
          case SEMI_MARKOVIAN : {
            occupancy = nonparametric_process[0]->sojourn_time[k];
            obs_product = 0.;
            forward1[k] = D_INF;
	    
            if (j < seq.length[i] - 1) {
              for (m = 1;m <= MIN(j + 1 , occupancy->nb_value - 1);m++) {
                if (observation[j - m + 1][k] == D_INF) {
                  break;
                }
                else {
                  obs_product += observation[j - m + 1][k];
                }
		
                if (m < j + 1) {
                  buff = obs_product + log(occupancy->mass[m]) + state_in[j - m][k];
                }
		
                else {
		  buff = obs_product + log(occupancy->mass[m]) + log(initial[k]);
		}
		
                if (buff > forward1[k]) {
                  forward1[k] = buff;
                  if (m < j + 1) {
                    optimal_state[j][k] = input_state[j - m][k];
                  }
                  optimal_occupancy[j][k] = m;
                }
              }
            }
	    
            else {
              for (m = 1;m <= MIN(j + 1 , occupancy->nb_value - 1);m++) {
                if (observation[j - m + 1][k] == D_INF) {
                  break;
                }
                else {
                  obs_product += observation[j - m + 1][k];
                }
		
                if (m < j + 1) {
                  buff = obs_product + log(occupancy->cumul[m - 1]) + state_in[j - m][k];
                }
		
                else {
		  buff = obs_product + log(occupancy->cumul[m - 1]) + log(initial[k]);
                }
		
                if (buff > forward1[k]) {
                  forward1[k] = buff;
                  if (m < j + 1) {
                    optimal_state[j][k] = input_state[j - m][k];
                  }
                  optimal_occupancy[j][k] = m;
                }
              }
            }
            break;
          }
	    
	    // cas etat markovien
	    
          case MARKOVIAN : {
            if (j == 0) {
              forward1[k] = log(initial[k]);
            }
            else {
              forward1[k] = state_in[j - 1][k];
              optimal_state[j][k] = input_state[j - 1][k];
            }
            optimal_occupancy[j][k] = 1;
	    
            if (forward1[k] != D_INF) {
              if (observation[j][k] == D_INF) {
                forward1[k] = D_INF;
              }
              else {
                forward1[k] += observation[j][k];
              }
            }
            break;
          }
          }
        }
	
#       ifdef DEBUG
        cout << j << " : ";
        for (k = 0;k < nb_state;k++) {
          cout << forward1[k];
          if (forward1[k] != D_INF) {
            cout << " " << optimal_occupancy[j][k] << " " << optimal_state[j][k];
          }
          cout << " | ";
        }
        cout << endl;
#       endif
	
	
        if (j < seq.length[i] - 1) {
          for (k = 0;k < nb_state;k++) {
            state_in[j][k] = D_INF;
            for (m = 0;m < nb_state;m++) {
              buff = log(transition[m][k]) + forward1[m];
              if (buff > state_in[j][k]) {
                state_in[j][k] = buff;
                input_state[j][k] = m;
              }
            }
          }
        }
	
        for (k = 0;k < nb_output_process;k++) {
          poutput[k]++;
        }
      }

      // extraction de la vraisemblance du chemin optimal
      
      pstate = seq.int_sequence[i][0] + seq.length[i] - 1;
      forward_max = D_INF;
      
      for (j = 0;j < nb_state;j++) {
        if (forward1[j] > forward_max) {
          forward_max = forward1[j];
          *pstate = j;
        }
      }

      if (forward_max != D_INF) {
        likelihood += forward_max;
        if (posterior_probability) {
          posterior_probability[i] = forward_max;
        }
      }

      else {
        likelihood = D_INF;
        if (posterior_probability) {
          posterior_probability[i] = 0.;
        }
        break;
      }

      // restauration

      j = seq.length[i]-1;

      do{
	for (k = 0; k < optimal_occupancy[j][*pstate]-1;k++) {
	  pstate--;
	  *pstate = *(pstate + 1);
	}

	if (j >= optimal_occupancy[j][*pstate]) {
          pstate--;
          *pstate = optimal_state[j][*(pstate + 1)];
          j -= optimal_occupancy[j][*(pstate + 1)];
        }
        else {
          j -= optimal_occupancy[j][*pstate];
        }
      }
      while (j >= 0);
    
      
      //#     ifdef DEBUG
      if(fag){
	//  cout << "\n";
	for ( j = 0; j < seq.length[i]; j++) {
	  cout << seq.int_sequence[i][0][j] << " ";
	}
	cout << endl;
      }
      //#     endif
      
      if(output_R){
	
	tmp = seq.max_length - seq.length[i]; 
	
	for( j = 0; j < tmp ; j++){
	  cout<<"-1 ";
	}
	for (j = tmp; j < seq.max_length; j++){
	  cout << seq.int_sequence[i][0][j-tmp] << " ";
	}
	cout<<endl;
	
      }
      
      //       //#     ifdef DEBUG
      //       cout << "\n";
      //       //    for (j = seq.length[i] - 1;j >= 0;j--) {
      //       for ( j = 0; j < seq.length[i]; j++) {
      //         cout << seq.int_sequence[i][0][j] << " ";
      //       }
      //       cout << endl;
      //       // #     endif

    }
  }

  for (i = 0;i < length;i++) {
    delete [] observation[i];
  }
  delete [] observation;

  delete [] forward1;
  
  for (i = 0;i < length - 1;i++) {
    delete [] state_in[i];
  }
  delete [] state_in;

  for (i = 0;i < length - 1;i++) {
    delete [] input_state[i];
  }
  delete [] input_state;
  
  for (i = 0;i < length;i++) {
    delete [] optimal_state[i];
  }
  delete [] optimal_state;
  
  for (i = 0;i < length;i++) {
    delete [] optimal_occupancy[i];
  }
  delete [] optimal_occupancy;
  
  delete [] poutput;
  
  return likelihood;
}

/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une chaine de Markov cachee a partir
 *  d'un echantillon de sequences.
 *
 *  arguments : reference sur un objet Format_error, stream, type de processus
 *              ('o' : ordinaire), nombre d'etats de la chaine de Markov,
 *              flag sur la nature de la chaine de Markov, flags sur le calcul
 *              des lois de comptage et sur le calcul des sequences d'etats optimales,
 *              probabilite de rester dans un etat initiale, nombre d'iterations.
 *
 *--------------------------------------------------------------*/

Markov_switching* Switching_sequence::markov_switching_estimation(Format_error &error , ostream &os ,
								  int nb_state , bool left_right ,
								  double **regression, double ***regmult,
								  int nb_covariable, int nb_cat, int constant, int ident, int link,
								  bool markov, bool output_R,
								  int order, bool counting_flag, bool state_sequence,
								  double occupancy_mean,
								  double self_transition, int nb_iter,
								  int mean_computation) const
{
  bool status = true;
  register int i,j;
  int nb_value[SEQUENCE_NB_VARIABLE];
  Chain *pchain;
  Discrete_parametric **disc_param;
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

    disc_param = new Discrete_parametric*[nb_state];
    for (i = 0; i < nb_state; i++) {
      disc_param[i] = new Discrete_parametric(NB_VALUE, nb_covariable, nb_cat, 
					      0, regression[i], regmult[i], ident, link, NO_RANDOM); 
    } 
    ihsw_process = new Switching_process(nb_state, disc_param);

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

    ihsw_markov = new Markov_switching(pchain, order, poccupancy, ihsw_process, 20);

    // initialisation des lois d'observation

    for (i = 0;i < ihsw_markov->nb_output_process;i++) {
      ihsw_markov->sw_process[i + 1]->init();
    }


    hsw_markov = new Markov_switching(*ihsw_markov, true);


    if( nb_cat > 0) {
      hsw_markov = markov_switching_estimation_multi(error , os , *ihsw_markov , nb_cat, COMPLETE_LIKELIHOOD, output_R,
						     counting_flag , state_sequence , nb_iter, mean_computation);
    }
    else{
      hsw_markov = markov_switching_estimation(error , os , *ihsw_markov , COMPLETE_LIKELIHOOD, output_R,
					       counting_flag , state_sequence , nb_iter, mean_computation);
    }
    delete ihsw_markov;
  }

  return hsw_markov;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une chaine de Markov cachee
 *  a partir d'un echantillon de sequences par l'algorithme SEM/MCEM.
 *
 *  arguments : reference sur un objet Format_error, stream, chaine de Markov cachee initiale,
 *              parametres pour le nombre de sequences d'etats simulees, flags sur le calcul
 *              des lois de comptage et sur le calcul des sequences d'etats optimales,
 *              nombre d'iterations.
 *
 *--------------------------------------------------------------*/
Markov_switching* Switching_sequence::markov_switching_stochastic_estimation(Format_error &error , std::ostream &os ,
									     const Markov_switching &ihsw_markov ,
									     int min_nb_state_sequence,
									     int max_nb_state_sequence,
									     bool output_R,
									     double parameter,
									     int estimator,
									     bool counting_flag,
									     bool state_sequence,
									     int nb_iter) const
{
  bool status;
  register int i , j , k , l, m , n;
  int max_nb_value , iter , nb_state_sequence , state_occupancy , nb_likelihood_decrease ,
    *occupancy_nb_value , *state_seq , *pstate ;
  int **poutput;
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

  int kk;

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

    hsmarkov = new Markov_switching(ihsw_markov , false, I_DEFAULT);//(int)(max_length*SAMPLE_NB_VALUE_COEFF));

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

    poutput = new int*[nb_variable];

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
          poutput[j] = int_sequence[i][j];
        }


        // recurrence "forward"

	for (j = 0;j < length[i];j++) {
          norm[j] = 0.;

          for (k = 0;k < hsmarkov->nb_state;k++) {

            // calcul des probabilites d'observation

            observation[j][k] = 1.;
            for (m = 0;m < hsmarkov->nb_output_process;m++) {
	      observation[j][k] *= hsmarkov->sw_process[m + 1]->observation[k]->mass_computation(*poutput[m], covar[i][j]);
            }

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
            poutput[m] = int_sequence[i][m] + k;
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

	      kk = k;
	      k -= (state_occupancy - 1);

	      for(m = k; m <= kk; m++){
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
      int s;
      gsl_permutation *p = gsl_permutation_alloc(nb_covariable);
      gsl_matrix *trans = gsl_matrix_alloc(nb_covariable, nb_covariable);
      gsl_matrix *beta = gsl_matrix_alloc(hsmarkov->nb_state, nb_covariable);
      gsl_vector *beta_k = gsl_vector_alloc(nb_covariable);
      gsl_vector *betaprime_k = gsl_vector_alloc(nb_covariable);
      gsl_matrix *inv_temp_add = gsl_matrix_alloc(nb_covariable, nb_covariable);
      gsl_matrix *temp_add = gsl_matrix_calloc(nb_covariable, nb_covariable);
      gsl_matrix *temp_mul2 = gsl_matrix_alloc(nb_covariable, nb_covariable);
      gsl_vector *add_vect = gsl_vector_calloc(nb_covariable); 
      gsl_vector *temp_mul3 = gsl_vector_alloc(nb_covariable);

      for (k = 0; k < hsmarkov->nb_state; k++) {

	for (j = 0; j < nb_covariable; j++) {
	  gsl_vector_set(beta_k, j, hsmarkov->sw_process[1]->observation[k]->regression[j] );
	}
	gsl_vector_set_zero(betaprime_k);
	gsl_matrix_set_zero(temp_add);
	gsl_matrix_set_zero(temp_mul2);
	gsl_vector_set_zero(add_vect);
	gsl_vector_set_zero(temp_mul3);

	// pour l'instant, cet algorithme ne marchera que pour une seule variable j
	for (j = 0 ; j < hsmarkov->nb_output_process; j++) {
	  
	  double sum = 0;
	  double *beta_kk;


	  for (i = 0; i < nb_sequence; i++) {
	    
	    gsl_matrix *X_i = gsl_matrix_calloc(length[i], nb_covariable);
	    gsl_matrix *indic_ik = gsl_matrix_calloc(length[i], length[i]);
	    gsl_matrix *temp_mat = gsl_matrix_calloc(length[i], hsmarkov->nb_state);
	    gsl_matrix *temp_mul = gsl_matrix_calloc(nb_covariable, length[i]);
	    gsl_matrix *temp_mul1 = gsl_matrix_calloc(nb_covariable, length[i]);


	      
	    convert_matrix_double(covar[i], length[i], nb_covariable, X_i);

	    
	    for (m = 0; m < nb_state_sequence; m++){

	      gsl_matrix *W_ik = gsl_matrix_calloc(length[i], length[i]);
	      gsl_vector *temp_vect = gsl_vector_calloc(length[i]);
	      gsl_vector *Y_i = gsl_vector_calloc(length[i]);
	      gsl_vector *mu_k = gsl_vector_calloc(length[i]);

	      convert_vector_int(int_sequence[i][j], length[i], Y_i);

	      for (l = 0; l < length[i]; l++) {
		gsl_matrix_set(W_ik, l, l, hsmarkov->sw_process[j + 1]->observation[k]->variance_computation(covar[i][l])); 
	      }

	      for (l = 0; l < length[i]; l++) {
		gsl_vector_set(mu_k, l, hsmarkov->sw_process[j + 1]->observation[k]->mean_computation(covar[i][l])); 
	      }

	      convert_matrix_double(indic[i][m], length[i], hsmarkov->nb_state, temp_mat);
	      
	      gsl_matrix_get_col(temp_vect, temp_mat, k);
	      convert_matrix_diag_double(temp_vect, length[i], indic_ik);

	      // 	      for (int v=0; v < length[i]; v++){
	      // 		  cout<< gsl_matrix_get(indic_ik,v,v)<<"  ";
	      // 	      }

      
	      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., X_i, indic_ik, 0., temp_mul); // t_XL

	      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., temp_mul, W_ik, 0., temp_mul1); // t_XLW	    
	  

	      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., temp_mul1, X_i, 0., temp_mul2); // t_XLWX
	      gsl_matrix_add(temp_add, temp_mul2);// somme des t_XLWX
	      
	      gsl_vector_sub(Y_i, mu_k); //(Y-mu)
	      gsl_blas_dgemv(CblasNoTrans, 1., temp_mul, Y_i, 0., temp_mul3);// t_XL(Y-mu)
	      gsl_vector_add(add_vect, temp_mul3);// somme des t_XL(Y-mu) 

	      gsl_vector_free(mu_k);
	      gsl_vector_free(Y_i);
	      gsl_vector_free(temp_vect);
	      gsl_matrix_free(W_ik);
	    
	    }


	    gsl_matrix_free(temp_mul1);
	    gsl_matrix_free(temp_mul);
	    gsl_matrix_free(temp_mat);
	    gsl_matrix_free(indic_ik);
	    gsl_matrix_free(X_i);
	    
	  }
	    
	  gsl_matrix_memcpy(trans, temp_add);
	  gsl_linalg_LU_decomp(trans, p, &s);
	  gsl_linalg_LU_invert(trans, p, inv_temp_add);

	  gsl_blas_dgemv(CblasNoTrans, 1., inv_temp_add, add_vect, 0., betaprime_k);
	  gsl_vector_add(beta_k, betaprime_k);
	  gsl_matrix_set_row(beta, k, beta_k);

	  // mise à jour des lois d'observations
	  

	  beta_kk = new double[nb_covariable];
	  convert_array1D_double (beta_k, nb_covariable, beta_kk);
	  hsmarkov->sw_process[j + 1]->observation[k]->Set_regression(beta_kk);


#   ifdef DEBUG
	  cout<<"beta pour l'état "<<k<<" :  "<<endl;
	  for(i = 0; i<nb_covariable; i++){
	    cout<<beta_kk[i]<<"   ";
	  }
	  cout<<endl;
#endif

	  delete [] beta_kk;

	}

      }
      gsl_vector_free(beta_k);
      gsl_vector_free(betaprime_k);
      gsl_vector_free(add_vect);
      gsl_vector_free(temp_mul3);
      gsl_matrix_free(temp_add);
      gsl_matrix_free(temp_mul2);
      gsl_matrix_free(inv_temp_add);
      gsl_matrix_free(beta);
      gsl_matrix_free(trans);
      gsl_permutation_free(p);


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
		complete_run[i]->state_occupancy_estimation(final_run[i] , complete_run[i] ,
							    occupancy_survivor ,
							    censored_occupancy_survivor , false);
	      }

            }

            complete_run[i]->max_computation();
            complete_run[i]->mean_computation();
            complete_run[i]->variance_computation();

            if (iter <= EXPLORATION_NB_ITER) {
              occupancy_likelihood = complete_run[i]->parametric_estimation(occupancy , 1 , true ,
                                                                            OCCUPANCY_THRESHOLD);
            }
            else {
              occupancy_likelihood = complete_run[i]->type_parametric_estimation(occupancy , 1 , true ,
                                                                                 OCCUPANCY_THRESHOLD);
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
    while ((likelihood != D_INF) && ((iter < STOCHASTIC_EXPLORATION_NB_ITER + 2) ||
				     ((nb_iter == I_DEFAULT) && (iter < MARKOV_SWITCHING_NB_ITER) &&
				      (((likelihood - previous_likelihood) / -likelihood > MARKOV_SWITCHING_LIKELIHOOD_DIFF) ||
				       (min_likelihood == D_INF) || (nb_likelihood_decrease == 1))) ||
				     ((nb_iter != I_DEFAULT) && (iter < nb_iter))));


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

      for (i = 0;i < hsmarkov->nb_output_process;i++) {
	for (j = 0;j < hsmarkov->nb_state;j++) {
	
	  hsmarkov->sw_process[i+1]->observation[j]->parametric_variance_computation();
	  // hsmarkov->sw_process[i+1]->observation[j]->min_value_computation();
	  // hsmarkov->sw_process[i+1]->observation[j]->max_value_computation();
	  hsmarkov->sw_process[i+1]->observation[j]->nb_value_computation();
	}
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

        for (i = 1;i <= hsmarkov->nb_output_process;i++) {
          if ((hsmarkov->sw_process[i]) && (seq->characteristics[i])) {
            delete seq->characteristics[i];
            seq->characteristics[i] = 0;
          }
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
	//         hsw_markov_data->min_computation(0);
	//         hsw_markov_data->max_computation(0);

        seq->build_marginal_histogram(0);
        seq->build_characteristic(0);

        seq->build_transition_count(*hsmarkov);
        
#       ifdef MESSAGE
        cout << "\n" << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << ": " << seq->likelihood <<endl;
	//  << " | " << hmarkov->Markov::likelihood_computation(*seq) << endl;
#       endif

	// calcul des lois d'occupation des etats

        for (i = 0;i < hsmarkov->nb_state;i++) {
          if (hsmarkov->state_subtype[i] == SEMI_MARKOVIAN) {
            hsmarkov->nonparametric_process[0]->sojourn_time[i]->computation((seq->characteristics[0] ? seq->characteristics[0]->sojourn_time[i]->nb_value : 1) ,
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
        seq->state_variable_init(INT_VALUE);

        for (i = 1;i <= hsmarkov->nb_output_process;i++) {
          if ((hsmarkov->sw_process[i]) && (seq->characteristics[i - 1])) {
            delete seq->characteristics[i - 1];
            seq->characteristics[i - 1] = 0;
          }
        }
      }

      // calcul de la vraisemblance et des lois caracteristiques du modele

      seq->hidden_likelihood = hsmarkov->likelihood_computation(*this , seq->posterior_probability, I_DEFAULT, true);

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
 *  Estimation des parametres d'une chaine de Markov cachee a partir
 *  d'un echantillon de sequences.
 *
 *  arguments : reference sur un objet Format_error, stream, type de processus
 *              ('o' : ordinaire, 'e' : en equilibre), nombre d'etats de la chaine de Markov,
 *              flag sur la nature de la chaine de Markov, parametres pour nombre de sequences
 *              d'etats simulees, flags sur le calcul des lois de comptage et sur le calcul
 *              des sequences d'etats optimales, probabilite de rester dans un etat initiale,
 *              nombre d'iterations.
 *
 *--------------------------------------------------------------*/
Markov_switching* Switching_sequence::markov_switching_stochastic_estimation(Format_error &error , std::ostream &os ,
									     char type , int nb_state , bool left_right ,
									     double **regression, double ***regmult,
									     int nb_covariable, int nb_cat,
									     int constant, int ident, int link,
									     int order, int min_nb_state_sequence,
									     int max_nb_state_sequence, bool markov, bool output_R,
									     double parameter, bool counting_flag,
									     bool state_sequence, double self_transition,
									     double occupancy_mean, int nb_iter) const
{
  bool status = true;
  register int i,j;
  int nb_value[SEQUENCE_NB_VARIABLE];
  Chain *pchain;
  Discrete_parametric **disc_param;
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

    disc_param = new Discrete_parametric*[nb_state];
    for (i = 0; i < nb_state; i++) {
      disc_param[i] = new Discrete_parametric(NB_VALUE, nb_covariable, nb_cat, 0, regression[i], regmult[i], ident, link, NO_RANDOM);
    }      
    ihsw_process = new Switching_process(nb_state, disc_param);


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

    ihsw_markov = new Markov_switching(pchain, order, poccupancy, ihsw_process, 20);

    // initialisation des lois d'observation

    for (i = 0;i < ihsw_markov->nb_output_process;i++) {
      ihsw_markov->sw_process[i + 1]->init();
    }


    hsw_markov = new Markov_switching(*ihsw_markov, true);

    if (nb_cat >0){
      hsw_markov = markov_switching_stochastic_estimation_multi(error, os, *ihsw_markov, nb_cat, min_nb_state_sequence,
								max_nb_state_sequence, false, parameter,
								COMPLETE_LIKELIHOOD, counting_flag,
								state_sequence, nb_iter);
    }
    else{
      hsw_markov = markov_switching_stochastic_estimation(error , os , *ihsw_markov , min_nb_state_sequence ,
							  max_nb_state_sequence , false, parameter , COMPLETE_LIKELIHOOD,
							  counting_flag , state_sequence , nb_iter);
    }    
    delete ihsw_markov;
  }

  return hsw_markov;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des probabilites a posteriori des etats par l'algorithme forward-backward.
 *
 *  arguments : reference sur un objet Markov_switching_data, indice de la sequence,
 *              stream, type de sortie, format de fichier ('a' : ASCII,
 *              's' : Spreadsheet, 'g' : Gnuplot), references sur l'entropie marginale
 *              maximum et l'entropie (pour la visualisation).
 *
 *--------------------------------------------------------------*/
double Markov_switching::forward_backward(const Markov_switching_data &seq , int index ,
                                          ostream &os , int output , char format,
					  double &max_marginal_entropy , double &entropy1) const

{
  register int i , j , k, m;
  int *pstate , **poutput;
  double seq_likelihood , state_seq_likelihood , obs_product, entropy2 , buff , sum, 
    backward_max , **observation,  *norm , *state_norm, **forward1 , **state_in, 
    **backward ,**backward1, *auxiliary , *occupancy_auxiliary, **backward_output , 
    *transition_predicted , *occupancy_predicted, **state_entropy , **predicted_entropy,
    **transition_entropy , **occupancy_entropy, *partial_entropy , *conditional_entropy , 
    *marginal_entropy;
  Parametric *occupancy;

  // initialisations

  observation = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    observation[i] = new double[nb_state];
  }

  norm = new double[seq.length[index]];
  state_norm = new double[nb_state];

  forward1 = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    forward1[i] = new double[nb_state];
  }

  state_in = new double*[seq.length[index] - 1];
  for (i = 0;i < seq.length[index] - 1;i++) {
    state_in[i] = new double[nb_state];
  }

  backward = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    backward[i] = new double[nb_state];
  }

  backward1 = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    backward1[i] = new double[nb_state];
  }

  auxiliary = new double[nb_state];
  occupancy_auxiliary = new double[seq.length[index] + 1];

  if (output == SSTATE) {
    backward_output = backward;
  }
  else {
    backward_output = new double*[seq.length[index]];
    for (i = 0;i < seq.length[index];i++) {
      backward_output[i] = new double[nb_state];
    }
  }

  transition_predicted = new double[nb_state];
  occupancy_predicted = new double[seq.length[index] + 1];

  state_entropy = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    state_entropy[i] = new double[nb_state];
  }

  predicted_entropy = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    predicted_entropy[i] = new double[nb_state];
  }

  transition_entropy = new double*[nb_state];
  for (i = 0;i < nb_state;i++) {
    transition_entropy[i] = new double[nb_state];
  }

  occupancy_entropy = new double*[nb_state];
  for (i = 0;i < nb_state;i++) {
    switch (state_subtype[i]) {
    case SEMI_MARKOVIAN :
      occupancy = nonparametric_process[0]->sojourn_time[i];
      occupancy_entropy[i] = new double[MIN(seq.length[index] , occupancy->nb_value)];
      break;
    case MARKOVIAN :
      occupancy_entropy[i] = 0;
      break;
    }
  }

  partial_entropy = new double[seq.length[index]];
  conditional_entropy = new double[seq.length[index]];
  marginal_entropy = new double[seq.length[index]];

# ifdef DEBUG
  double *backward0;

  backward0 = new double[nb_state];
# endif

  poutput = new int*[nb_output_process];

  for (i = 0;i < nb_output_process;i++) {
    poutput[i] = seq.int_sequence[index][i + 1];
  }

  // recurrence "forward"


  seq_likelihood = 0.;
  for (i = 0;i < seq.length[index];i++) {
    norm[i] = 0.;

    for (j = 0;j < nb_state;j++) {

      // calcul des probabilites d'observation

      observation[i][j] = 1.;
      for (k = 0;k < nb_output_process;k++) {
	observation[i][j] *= sw_process[k + 1]->observation[j]->mass_computation(*poutput[k], seq.covar[index][i]);
      }

      switch (state_subtype[j]) {
	
	// cas etat semi-markovien 
	
      case SEMI_MARKOVIAN : {
        if (i == 0) {
          state_norm[j] = initial[j];
        }
        else {
          state_norm[j] += state_in[i - 1][j] - forward1[i - 1][j];
        }
        state_norm[j] *= observation[i][j];

        norm[i] += state_norm[j];
        break;
      }

	// cas etat markovien

      case MARKOVIAN : {
        if (i == 0) {
          forward1[i][j] = initial[j];
          state_entropy[i][j] = 0.;
        }
        else {
          forward1[i][j] = state_in[i - 1][j];
          state_entropy[i][j] = predicted_entropy[i - 1][j];
        }
        forward1[i][j] *= observation[i][j];

        norm[i] += forward1[i][j];
        break;
      }
      }
    }

    if (norm[i] > 0.) {
      for (j = 0;j < nb_state;j++) {
        switch (state_subtype[j]) {
        case SEMI_MARKOVIAN :
          state_norm[j] /= norm[i];
          break;
        case MARKOVIAN :
          forward1[i][j] /= norm[i];
          break;
        }
      }

      seq_likelihood += log(norm[i]);
    }

    else {
      seq_likelihood = D_INF;
      break;
    }

    for (j = 0;j < nb_state;j++) {

      // cas etat semi-markovien

      if (state_subtype[j] == SEMI_MARKOVIAN) {
        occupancy = nonparametric_process[0]->sojourn_time[j];
        obs_product = 1.;
        forward1[i][j] = 0.;

        if (i < seq.length[index] - 1) {
          for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
            obs_product *= observation[i - k + 1][j] / norm[i - k + 1];
            if (obs_product == 0.) {
              break;
            }

            if (k < i + 1) {
              occupancy_predicted[k] = obs_product * occupancy->mass[k] * state_in[i - k][j];
	      //              forward1[i][j] += obs_product * occupancy->mass[k] * state_in[i - k][j];
            }

            else {
              switch (type) {
              case 'o' :
                occupancy_predicted[k] = obs_product * occupancy->mass[k] * initial[j];
		//                forward1[i][j] += obs_product * occupancy->mass[k] * initial[j];
                break;
              case 'e' :
                occupancy_predicted[k] = obs_product * forward[j]->mass[k] * initial[j];
		//                forward1[i][j] += obs_product * forward[j]->mass[k] * initial[j];
                break;
              }
            }

            forward1[i][j] += occupancy_predicted[k];
          }
        }

        else {
          for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
            obs_product *= observation[i - k + 1][j] / norm[i - k + 1];
            if (obs_product == 0.) {
              break;
            }

            if (k < i + 1) {
              occupancy_predicted[k] = obs_product * (1. - occupancy->cumul[k - 1]) * state_in[i - k][j];
	      //              forward1[i][j] += obs_product * (1. - occupancy->cumul[k - 1]) * state_in[i - k][j];
            }

            else {
              switch (type) {
              case 'o' :
                occupancy_predicted[k] = obs_product * (1. - occupancy->cumul[k - 1]) * initial[j];
		//                forward1[i][j] += obs_product * (1. - occupancy->cumul[k - 1]) * initial[j];
                break;
              case 'e' :
                occupancy_predicted[k] = obs_product * (1. - forward[j]->cumul[k - 1]) * initial[j];
		//                forward1[i][j] += obs_product * (1. - forward[j]->cumul[k - 1]) * initial[j];
                break;
              }
            }

            forward1[i][j] += occupancy_predicted[k];
          }
        }

        state_entropy[i][j] = 0.;

        if (forward1[i][j] > 0.) {
          for (m = 1;m < k;m++) {
            buff = occupancy_predicted[m] / forward1[i][j];
            if (buff > 0.) {
              if (m < i + 1) {
                state_entropy[i][j] += buff * (predicted_entropy[i - m][j] - log(buff));
              }
              else {
                state_entropy[i][j] -= buff * log(buff);
              }
            }
          }

          if (state_entropy[i][j] < 0.) {
            state_entropy[i][j] = 0.;
          }
        }
      }
    }

    if (i < seq.length[index] - 1) {
      for (j = 0;j < nb_state;j++) {
        state_in[i][j] = 0.;
        for (k = 0;k < nb_state;k++) {
          transition_predicted[k] = transition[k][j] * forward1[i][k];
          state_in[i][j] += transition_predicted[k];

	  //          state_in[i][j] += transition[k][j] * forward1[i][k];
        }

        predicted_entropy[i][j] = 0.;

        if (state_in[i][j] > 0.) {
          for (k = 0;k < nb_state;k++) {
            buff = transition_predicted[k] / state_in[i][j];
            if (buff > 0.) {
              predicted_entropy[i][j] += buff * (state_entropy[i][k] - log(buff));
            }
          }

          if (predicted_entropy[i][j] < 0.) {
            predicted_entropy[i][j] = 0.;
          }
        }
      }
    }

    for (j = 0;j < nb_output_process;j++) {
      poutput[j]++;
    }
  }

  if (seq_likelihood != D_INF) {
    entropy1 = 0.;
    i = seq.length[index] - 1;
    for (j = 0;j < nb_state;j++) {
      if (forward1[i][j] > 0.) {
        entropy1 += forward1[i][j] * (state_entropy[i][j] - log(forward1[i][j]));
      }
    }

    for (i = 0;i < nb_state;i++) {
      if (state_subtype[i] == SEMI_MARKOVIAN) {
        for (j = 0;j < seq.length[index];j++) {
          state_entropy[j][i] = 0.;
        }
      }
    }

    // recurrence "backward"

    for (i = 0;i < nb_output_process;i++) {
      poutput[i]--;
    }

    for (i = 0;i < nb_state;i++) {
      for (j = 0;j < nb_state;j++) {
        transition_entropy[i][j] = 0.;
      }
    }

    for (i = 0;i < nb_state;i++) {
      if (state_subtype[i] == SEMI_MARKOVIAN) {
        occupancy = nonparametric_process[0]->sojourn_time[i];
        for (j = occupancy->offset;j < MIN(seq.length[index] , occupancy->nb_value);j++) {
          occupancy_entropy[i][j] = 0.;
        }
      }
    }

    entropy2 = 0.;

    i = seq.length[index] - 1;
    for (j = 0;j < nb_state;j++) {
      backward[i][j] = forward1[i][j];
      backward1[i][j] = backward[i][j];

      if (output == OUT_STATE) {
        backward_output[i][j] = backward[i][j];
      }

      if (backward[i][j] > 0.) {
        for (k = 0;k < nb_output_process;k++) {
          if (sw_process[k + 1]) {
	    if (sw_process[k + 1]->observation[j]->mass_computation(*poutput[k], seq.covar[index][i]) > 0.) {
	      entropy2 -= backward[i][j] * log(sw_process[k + 1]->observation[j]->mass_computation(*poutput[k], 
												   seq.covar[index][i]));
	    }
	  }
	}
      }
    }

    for (i = seq.length[index] - 2;i >= 0;i--) {
      for (j = 0;j < nb_output_process;j++) {
        poutput[j]--;
      }


      for (j = 0;j < nb_state;j++) {
        auxiliary[j] = 0.;

        switch (state_subtype[j]) {

	  // cas etat semi-markovien

        case SEMI_MARKOVIAN : {
          occupancy = nonparametric_process[0]->sojourn_time[j];
          obs_product = 1.;

          for (k = 1;k < MIN(seq.length[index] - i , occupancy->nb_value);k++) {
            obs_product *= observation[i + k][j] / norm[i + k];
            if (obs_product == 0.) {
              break;
            }

            occupancy_auxiliary[k] = 0.;

            if (backward1[i + k][j] > 0.) {
	      //            if (forward1[i + k][j] > 0.) {
              if (k < seq.length[index] - i - 1) {
                buff = backward1[i + k][j] * obs_product * occupancy->mass[k] /
		  forward1[i + k][j];
                occupancy_auxiliary[k] = buff * state_in[i][j];
                occupancy_entropy[j][k] += occupancy_auxiliary[k];

		/*                if (occupancy->mass[k] > 0.) {
				  entropy2 -= occupancy_auxiliary[k] * log(occupancy->mass[k]);
				  } */
              }

              else {
                buff = obs_product * (1. - occupancy->cumul[k - 1]);
                occupancy_auxiliary[k] = buff * state_in[i][j];
                if (occupancy->cumul[k - 1] < 1.) {
                  entropy2 -= occupancy_auxiliary[k] * log(1. - occupancy->cumul[k - 1]);
                }
              }

              auxiliary[j] += buff;
            }
          }

          sum = 0.;
          for (m = k - 1;m >= 1;m--) {
            sum += occupancy_auxiliary[m];
            if (backward[i + m][j] > 0.) {
              buff = sum / backward[i + m][j];
              if (buff > 0.) {
                state_entropy[i + m][j] += buff * (predicted_entropy[i][j] - log(buff));
              }
            }
          }
          break;
        }

	  // cas etat markovien

        case MARKOVIAN : {
          if (backward1[i + 1][j] > 0.) {
	    //          if (forward1[i + 1][j] > 0.) {
            auxiliary[j] = backward1[i + 1][j] / state_in[i][j];

	    /*            auxiliary[j] = backward1[i + 1][j] * observation[i + 1][j] /
			  (forward1[i + 1][j] * norm[i + 1]); */

            state_entropy[i + 1][j] = predicted_entropy[i][j];
          }
          break;
        }
        }
      }

      for (j = 0;j < nb_state;j++) {
        backward1[i][j] = 0.;

        for (k = 0;k < nb_state;k++) {
          buff = auxiliary[k] * transition[j][k] * forward1[i][j];
          backward1[i][j] += buff;
          transition_entropy[j][k] += buff;

	  /*          if (transition[j][k] > 0.) {
		      entropy2 -= buff * log(transition[j][k]);
		      } */
        }

        switch (state_subtype[j]) {

	  // cas etat semi-markovien

        case SEMI_MARKOVIAN : {
          backward[i][j] = backward[i + 1][j] + backward1[i][j] - auxiliary[j] * state_in[i][j];
          if (backward[i][j] < 0.) {
            backward[i][j] = 0.;
          }
          if (backward[i][j] > 1.) {
            backward[i][j] = 1.;
          }
          break;
        }

	  // cas etat markovien

        case MARKOVIAN : {
          backward[i][j] = backward1[i][j];
          break;
        }
        }

	if (backward[i][j] > 0.) {
	  for (k = 0;k < nb_output_process;k++) {
            if (sw_process[k + 1]) {
	      if (sw_process[k + 1]->observation[j]->mass_computation(*poutput[k], seq.covar[index][i]) > 0.) {
		entropy2 -= backward[i][j] * log(sw_process[k + 1]->observation[j]->mass_computation(*poutput[k], 
												     seq.covar[index][i]));
	      }
	    }
	  }
	}
      }

      switch (output) {

      case IN_STATE : {
        for (j = 0;j < nb_state;j++) {
          switch (state_subtype[j]) {

	    // cas etat semi-markovien

          case SEMI_MARKOVIAN : {
            backward_output[i + 1][j] = auxiliary[j] * state_in[i][j];
            break;
          }

	    // cas etat markovien

          case MARKOVIAN : {
            backward_output[i + 1][j] = 0.;
            for (k = 0;k < nb_state;k++) {
              if (k != j) {
                backward_output[i + 1][j] += transition[k][j] * forward1[i][k];
              }
            }
            backward_output[i + 1][j] *= auxiliary[j];
            break;
          }
          }
        }
        break;
      }

      case OUT_STATE : {
        for (j = 0;j < nb_state;j++) {
          switch (state_subtype[j]) {

	    // cas etat semi-markovien

          case SEMI_MARKOVIAN : {
            backward_output[i][j] = backward1[i][j];
            break;
          }

	    // cas etat markovien

          case MARKOVIAN : {
            backward_output[i][j] = 0.;
            for (k = 0;k < nb_state;k++) {
              if (k != j) {
                backward_output[i][j] += auxiliary[k] * transition[j][k];
              }
            }
            backward_output[i][j] *= forward1[i][j];
            break;
          }
          }
        }
        break;
      }
      }
    }

    if (output == IN_STATE) {
      for (i = 0;i < nb_state;i++) {
        backward_output[0][i] = backward[0][i];
      }
    }

    for (i = 0;i < nb_state;i++) {
      if (initial[i] > 0.) {
        entropy2 -= backward[0][i] * log(initial[i]);
      }
    }

    for (i = 0;i < nb_state;i++) {
      for (j = 0;j < nb_state;j++) {
        if (transition[i][j] > 0.) {
          entropy2 -= transition_entropy[i][j] * log(transition[i][j]);
        }
      }
    }

    for (i = 0;i < nb_state;i++) {
      if (state_subtype[i] == SEMI_MARKOVIAN) {
        occupancy = nonparametric_process[0]->sojourn_time[i];

        if (initial[i] > 0.) {
          obs_product = 1.;

#         ifdef DEBUG
          backward0[i] = 0.;
#         endif

          for (j = 1;j < MIN(seq.length[index] + 1 , occupancy->nb_value);j++) {
            obs_product *= observation[j - 1][i] / norm[j - 1];
            if (obs_product == 0.) {
              break;
            }

            occupancy_auxiliary[j] = 0.;

            if (backward1[j - 1][i] > 0.) {
	      //            if (forward1[j - 1][i] > 0.) {
              if (j < seq.length[index]) {
                switch (type) {

                case 'o' : {
                  occupancy_auxiliary[j] = backward1[j - 1][i] * obs_product * occupancy->mass[j] *
		    initial[i] / forward1[j - 1][i];
                  occupancy_entropy[i][j] += occupancy_auxiliary[j];

		  /*                  if (occupancy->mass[j] > 0.) {
				      entropy2 -= occupancy_auxiliary[j] * log(occupancy->mass[j]);
				      } */
                  break;
                }

                case 'e' : {
                  occupancy_auxiliary[j] = backward1[j - 1][i] * obs_product * forward[i]->mass[j] *
		    initial[i] / forward1[j - 1][i];
                  if (forward[i]->mass[j] > 0.) {
                    entropy2 -= occupancy_auxiliary[j] * log(forward[i]->mass[j]);
                  }
                  break;
                }
                }
              }

              else {
                switch (type) {

                case 'o' : {
                  occupancy_auxiliary[j] = obs_product * (1. - occupancy->cumul[j - 1]) * initial[i];
                  if (occupancy->cumul[j - 1] < 1.) {
                    entropy2 -= occupancy_auxiliary[j] * log(1. - occupancy->cumul[j - 1]);
                  }
                  break;
                }

                case 'e' : {
                  occupancy_auxiliary[j] = obs_product * (1. - forward[i]->cumul[j - 1]) * initial[i];
                  if (forward[i]->cumul[j - 1] < 1.) {
                    entropy2 -= occupancy_auxiliary[j] * log(1. - forward[i]->cumul[j - 1]);
                  }
                  break;
                }
                }
              }

#             ifdef DEBUG
              backward0[i] += occupancy_auxiliary[j];
#             endif

            }
          }

#         ifdef DEBUG
          cout << i << " " << backward[0][i] << " " << backward0[i] << endl;
#         endif

          sum = 0.;
          for (k = j - 1;k >= 1;k--) {
            sum += occupancy_auxiliary[k];
            if (backward[k - 1][i] > 0.) {
              buff = sum / backward[k - 1][i];
              if (buff > 0.) {
                state_entropy[k - 1][i] -= buff * log(buff);
              }
            }
          }
        }

        for (j = occupancy->offset;j < MIN(seq.length[index] , occupancy->nb_value);j++) {
          if (occupancy->mass[j] > 0.) {
            entropy2 -= occupancy_entropy[i][j] * log(occupancy->mass[j]);
          }
        }
      }
    }

    entropy2 += seq_likelihood;

#   ifdef MESSAGE
    if ((entropy2 < entropy1 - DOUBLE_ERROR) || (entropy2 > entropy1 + DOUBLE_ERROR)) {
      cout << "\nERROR: " << entropy1 << " " << entropy2 << endl;
    }
#   endif

    // restauration

    pstate = seq.int_sequence[index][0];

    for (i = 0;i < seq.length[index];i++) {
      backward_max = 0.;
      for (j = 0;j < nb_state;j++) {
        if (backward[i][j] > backward_max) {
          backward_max = backward[i][j];
          *pstate = j;
        }
      }

      pstate++;
    }

    state_seq_likelihood = Markov_switching::likelihood_computation(seq , index);
    for (i = 0;i < seq.length[index];i++) {
      partial_entropy[i] = 0.;
      for (j = 0;j < nb_state;j++) {
        if (state_entropy[i][j] < 0.) {
          state_entropy[i][j] = 0.;
        }
        if (backward[i][j] > 0.) {
          partial_entropy[i] += backward[i][j] * (state_entropy[i][j] - log(backward[i][j]));
        }
      }
      if (partial_entropy[i] < 0.) {
        partial_entropy[i] = 0.;
      }
    }

    conditional_entropy[0] = partial_entropy[0];
    for (i = 1;i < seq.length[index];i++) {
      conditional_entropy[i] = partial_entropy[i] - partial_entropy[i - 1];
    }

    max_marginal_entropy = 0.;
    for (i = 0;i < seq.length[index];i++) {
      marginal_entropy[i] = 0.;
      for (j = 0;j < nb_state;j++) {
        if (backward[i][j] > 0.) {
          marginal_entropy[i] -= backward[i][j] * log(backward[i][j]);
        }
      }
      if (marginal_entropy[i] > max_marginal_entropy) {
        max_marginal_entropy = marginal_entropy[i];
      }
    }

    switch (format) {

    case 'a' : {
      switch (output) {
      case SSTATE :
        os << "\n" << SEQ_label[SEQL_POSTERIOR_STATE_PROBABILITY] << "\n\n";
        break;
      case IN_STATE :
        os << "\n" << SEQ_label[SEQL_POSTERIOR_IN_STATE_PROBABILITY] << "\n\n";
        break;
      case OUT_STATE :
        os << "\n" << SEQ_label[SEQL_POSTERIOR_OUT_STATE_PROBABILITY] << "\n\n";
        break;
      }

      seq.profile_ascii_print(os , index , nb_state , backward_output ,
			      STAT_label[STATL_STATE]);

      //      seq.profile_ascii_print(os , index , nb_state , backward_output , conditional_entropy ,
      //                              marginal_entropy , partial_entropy);

      os << "\n" << STAT_label[STATL_LIKELIHOOD] << ": " << seq_likelihood
         << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_LIKELIHOOD] << ": " << state_seq_likelihood
         << "   (" << exp(state_seq_likelihood - seq_likelihood) << ")" << endl;
      break;
    }

    case 's' : {
      switch (output) {
      case SSTATE :
        os << "\n" << SEQ_label[SEQL_POSTERIOR_STATE_PROBABILITY] << "\n\n";
        break;
      case IN_STATE :
        os << "\n" << SEQ_label[SEQL_POSTERIOR_IN_STATE_PROBABILITY] << "\n\n";
        break;
      case OUT_STATE :
        os << "\n" << SEQ_label[SEQL_POSTERIOR_OUT_STATE_PROBABILITY] << "\n\n";
        break;
      }

      seq.profile_spreadsheet_print(os , index , nb_state , backward_output ,
				    STAT_label[STATL_STATE]);

      //      seq.profile_spreadsheet_print(os , index , nb_state , backward_output , conditional_entropy ,
      //                                    marginal_entropy , partial_entropy) ;

      os << "\n" << STAT_label[STATL_LIKELIHOOD] << "\t" << seq_likelihood
         << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_LIKELIHOOD] << "\t" << state_seq_likelihood
         << "\t" << exp(state_seq_likelihood - seq_likelihood) << endl;
      break;
    }

    case 'g' : {
      //      seq.profile_plot_print(os , index , nb_state , backward_output);
      seq.profile_plot_print(os , index , nb_state , backward_output , conditional_entropy ,
                             marginal_entropy , partial_entropy);
      break;
    }
    }

    if (format != 'g') {
      /*      double gini_index;

      gini_index = 0.;
      for (i = 0;i < seq.length[index];i++) {
      for (j = 0;j < nb_state;j++) {
      gini_index += backward[i][j] * (1. - backward[i][j]);
      }
      } */

      double entropy3 , nb_state_sequence;

      entropy3 = 0.;
      for (i = 0;i < seq.length[index];i++) {
        for (j = 0;j < nb_state;j++) {
          if (backward[i][j] > 0.) {
            entropy3 -= backward[i][j] * log(backward[i][j]);
          }
        }
      }

      // calcul du nombre de sequences d'etats possibles

      for (i = 0;i < nb_output_process;i++) {
        poutput[i] = seq.int_sequence[index][i + 1];
      }

      // recurrence "forward"


      for (i = 0;i < seq.length[index];i++) {
        for (j = 0;j < nb_state;j++) {

          // calcul des probabilites d'observation

          observation[i][j] = 1.;
          for (k = 0;k < nb_output_process;k++) {
	    if (sw_process[k + 1]) {
	      observation[i][j] *= sw_process[k + 1]->observation[j]->mass_computation(*poutput[k], seq.covar[index][i]);
	    }
	  }
	
	  if (observation[i][j] > 0.) {
	    observation[i][j] = 1.;
	  }

	  forward1[i][j] = 0.;

	  switch (state_subtype[j]) {

	    // cas etat semi-markovien

	  case SEMI_MARKOVIAN : {
	    occupancy = nonparametric_process[0]->sojourn_time[j];

	    if (i < seq.length[index] - 1) {
	      for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
		if (observation[i - k + 1][j] == 0.) {
		  break;
		}

		if (k < i + 1) {
		  if (occupancy->mass[k] > 0.) {
		    forward1[i][j] += state_in[i - k][j];
		  }
		}

		else {
		  if (initial[j] > 0.) {
		    switch (type) {

		    case 'o' : {
		      if (occupancy->mass[k] > 0.) {
			forward1[i][j]++;
		      }
		      break;
		    }

		    case 'e' : {
		      if (forward[j]->mass[k] > 0.) {
			forward1[i][j]++;
		      }
		      break;
		    }
		    }
		  }
		}
	      }
	    }

            else {
              for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
                if (observation[i - k + 1][j] == 0.) {
                  break;
                }

                if (k < i + 1) {
                  if (1. - occupancy->cumul[k - 1] > 0.) {
                    forward1[i][j] += state_in[i - k][j];
                  }
                }

                else {
                  if (initial[j] > 0.) {
                    switch (type) {

                    case 'o' : {
                      if (1. - occupancy->cumul[k - 1] > 0.) {
                        forward1[i][j]++;
                      }
                      break;
                    }

                    case 'e' : {
                      if (1. - forward[j]->cumul[k - 1] > 0.) {
                        forward1[i][j]++;
                      }
                      break;
                    }
                    }
                  }
                }
              }
            }
            break;
          }

	    // cas etat markovien

          case MARKOVIAN : {
            if (observation[i][j] == 1.) {
              if (i == 0) {
                if (initial[j] > 0.) {
                  forward1[i][j] = 1.;
                }
              }
              else {
                forward1[i][j] = state_in[i - 1][j];
              }
            }
            break;
          }
          }
        }

        if (i < seq.length[index] - 1) {
          for (j = 0;j < nb_state;j++) {
            state_in[i][j] = 0.;
            for (k = 0;k < nb_state;k++) {
              if (transition[k][j] > 0.) {
                state_in[i][j] += forward1[i][k];
              }
            }
          }
        }

        for (j = 0;j < nb_output_process;j++) {
          poutput[j]++;
        }
      }

      nb_state_sequence = 0.;
      i = seq.length[index] - 1;
      for (j = 0;j < nb_state;j++) {
        nb_state_sequence += forward1[i][j];
      }

      switch (format) {
      case 'a' :
	/*        os << "\n" << SEQ_label[SEQL_GINI_INDEX] << ": " << gini_index << " ("
		  << gini_index / seq.length[index] << ")   " << SEQ_label[SEQL_UPPER_BOUND] << ": "
		  << seq.length[index] * (1. - 1. / nb_state) << " (" << 1. - 1. / nb_state
		  os << ")   " << SEQ_label[SEQL_UPPER_BOUND] << ": "
		  << seq.length[index] * log((double)nb_state) << " (" << log((double)nb_state) */
        os << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_ENTROPY] << ": " << entropy1
           << " (" << entropy1 / seq.length[index] << ")   " << SEQ_label[SEQL_UPPER_BOUND] << ": "
           << log((double)nb_state_sequence) << " ("
           << log((double)nb_state_sequence) / seq.length[index]
           << ")\n" << SEQ_label[SEQL_MARGINAL_ENTROPY_SUM] << ": " << entropy3 << " ("
           << entropy3 / seq.length[index] << ")\n\n"
           << SEQ_label[SEQL_NB_STATE_SEQUENCE] << ": " << nb_state_sequence << endl;
        break;
      case 's' :
	/*        os << "\n" << SEQ_label[SEQL_GINI_INDEX] << "\t" << gini_index << "\t"
		  << gini_index / seq.length[index] << "\t" << SEQ_label[SEQL_UPPER_BOUND] << "\t"
		  << seq.length[index] * (1. - 1. / nb_state) << "\t" << 1. - 1. / nb_state
		  os << "\t" << SEQ_label[SEQL_UPPER_BOUND] << "\t"
		  << seq.length[index] * log((double)nb_state) << "\t" << log((double)nb_state) */
        os << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_ENTROPY] << "\t" << entropy1
           << "\t" << entropy1 / seq.length[index] << "\t" << SEQ_label[SEQL_UPPER_BOUND] << "\t"
           << log((double)nb_state_sequence) << "\t"
           << log((double)nb_state_sequence) / seq.length[index]
           << "\n" << SEQ_label[SEQL_MARGINAL_ENTROPY_SUM] << "\t" << entropy3 << "\t"
           << entropy3 / seq.length[index] << "\n\n"
           << SEQ_label[SEQL_NB_STATE_SEQUENCE] << "\t" << nb_state_sequence << endl;
        break;
      }

#     ifdef DEBUG
      int state;
      double min_nb_state_sequence , smoothed_proba , cumul_smoothed_proba ,
	max_smoothed_proba , **backward2;

      // recurrence "backward"

      min_nb_state_sequence = nb_state_sequence;

      backward2 = new double*[seq.length[index]];
      for (i = 0;i < seq.length[index];i++) {
        backward2[i] = new double[nb_state];
      }

      i = seq.length[index] - 1;
      for (j = 0;j < nb_state;j++) {
        backward2[i][j] = forward1[i][j];
        backward1[i][j] = 1.;
      }

      for (i = seq.length[index] - 2;i >= 0;i--) {
        for (j = 0;j < nb_state;j++) {
          auxiliary[j] = 0.;

          switch (state_subtype[j]) {

	    // cas etat semi-markovien

          case SEMI_MARKOVIAN : {
            occupancy = nonparametric_process[0]->sojourn_time[j];

            for (k = 1;k < MIN(seq.length[index] - i , occupancy->nb_value);k++) {
              if (observation[i + k][j] == 0.) {
                break;
              }

              if (k < seq.length[index] - i - 1) {
                if (occupancy->mass[k] > 0.) {
                  auxiliary[j] += backward1[i + k][j];
                }
              }
              else {
                if (1. - occupancy->cumul[k - 1] > 0.) {
                  auxiliary[j]++;
                }
              }
            }
            break;
          }

	    // cas etat markovien

          case MARKOVIAN : {
            if (observation[i + 1][j] == 1.) {
              auxiliary[j] = backward1[i + 1][j];
            }
            break;
          }
          }
        }

        for (j = 0;j < nb_state;j++) {
          backward1[i][j] = 0.;

          for (k = 0;k < nb_state;k++) {
            if (transition[j][k] > 0.) {
              backward1[i][j] += auxiliary[k];
            }
          }

          switch (state_subtype[j]) {

	    // cas etat semi-markovien

          case SEMI_MARKOVIAN : {

#           ifdef DEBUG
            if ((i == 0) && (initial[j] > 0.)) {
              occupancy = nonparametric_process[0]->sojourn_time[j];
              backward0[j] = 0.;

              for (k = 1;k < MIN(seq.length[index] + 1 , occupancy->nb_value);k++) {
                if (observation[k - 1][j] == 0.) {
                  break;
                }

                if (k < seq.length[index]) {
                  if (occupancy->mass[k] > 0.) {
                    backward0[j] += backward1[k - 1][j];
                  }
                }
                else {
                  if (1. - occupancy->cumul[k - 1] > 0.) {
                    backward0[j]++;
                  }
                }
              }
            }
#           endif

            backward2[i][j] = backward2[i + 1][j] + backward1[i][j] * forward1[i][j] -
	      auxiliary[j] * state_in[i][j];

#           ifdef DEBUG
            if ((i == 0) && (initial[j] > 0.)) {
              cout << j << " " << backward2[i][j] << " " << backward0[j] << endl;
            }
#           endif

            break;
          }

	    // cas etat markovien

          case MARKOVIAN : {
            backward2[i][j] = backward1[i][j] * forward1[i][j];
            break;
          }
          }
        }

        smoothed_proba = 1.1;
        cumul_smoothed_proba = 0.;
        nb_state_sequence = 0;

        for (j = 0;j < nb_state;j++) {
          max_smoothed_proba = 0.;
          for (k = 0;k < nb_state;k++) {
            if ((backward[i][k] > max_smoothed_proba) && (backward[i][k] < smoothed_proba)) {
              max_smoothed_proba = backward[i][k];
              state = k;
            }
          }
          cumul_smoothed_proba += max_smoothed_proba;
          nb_state_sequence += backward2[i][state];

          if (cumul_smoothed_proba < 1. - MIN_SMOOTHED_PROBABILITY) {
            smoothed_proba = max_smoothed_proba;
          }
          else {
            break;
          }
        }

        if (nb_state_sequence < min_nb_state_sequence) {
          min_nb_state_sequence = nb_state_sequence;
        }
      }

      os << SEQ_label[SEQL_NB_STATE_SEQUENCE]
         << " (" << 1. - MIN_SMOOTHED_PROBABILITY << " beam)"
         << ": " << min_nb_state_sequence << endl;

      os << "\n";
      for (i = 0;i < seq.length[index];i++) {
        obs_product = 0.;
        for (j = 0;j < nb_state;j++) {
          os << backward2[i][j] << " (" << backward[i][j] << ")  ";
          obs_product += backward2[i][j];
        }
        os << "| " << obs_product << endl;
      }

      for (i = 0;i < seq.length[index];i++) {
        delete [] backward2[i];
      }
      delete [] backward2;
#     endif

    }
  }

  for (i = 0;i < seq.length[index];i++) {
    delete [] observation[i];
  }
  delete [] observation;

  delete [] norm;
  delete [] state_norm;

  for (i = 0;i < seq.length[index];i++) {
    delete [] forward1[i];
  }
  delete [] forward1;

  for (i = 0;i < seq.length[index] - 1;i++) {
    delete [] state_in[i];
  }
  delete [] state_in;

  for (i = 0;i < seq.length[index];i++) {
    delete [] backward[i];
  }
  delete [] backward;

  for (i = 0;i < seq.length[index];i++) {
    delete [] backward1[i];
  }
  delete [] backward1;

  delete [] auxiliary;
  delete [] occupancy_auxiliary;

  if (output != SSTATE) {
    for (i = 0;i < seq.length[index];i++) {
      delete [] backward_output[i];
    }
    delete [] backward_output;
  }

  delete [] transition_predicted;
  delete [] occupancy_predicted;

  for (i = 0;i < seq.length[index];i++) {
    delete [] state_entropy[i];
  }
  delete [] state_entropy;

  for (i = 0;i < seq.length[index];i++) {
    delete [] predicted_entropy[i];
  }
  delete [] predicted_entropy;

  for (i = 0;i < nb_state;i++) {
    delete [] transition_entropy[i];
  }
  delete [] transition_entropy;

  for (i = 0;i < nb_state;i++) {
    delete [] occupancy_entropy[i];
  }
  delete [] occupancy_entropy;

  delete [] partial_entropy;
  delete [] conditional_entropy;
  delete [] marginal_entropy;

# ifdef DEBUG
  delete [] backward0;
# endif

  delete [] poutput;

  return (seq_likelihood);
}


/*--------------------------------------------------------------*
 *
 *  Simulation de L sequences d'etats correspondant a une sequence observee
 *  par l'algorithme forward-backward.
 *
 *  arguments : reference sur un objet Markov_switching_data, indice de la sequence,
 *              stream, format de fichier ('a' : ASCII, 's' : Spreadsheet),
 *              nombre de sequences d'etats.
 *
 *--------------------------------------------------------------*/
double Markov_switching::forward_backward_sampling(const Markov_switching_data &seq ,
                                                   int index , ostream &os ,
                                                   char format , int nb_state_sequence) const
{
  register int i , j , k;
  int state_occupancy, *pstate, **poutput;
  double seq_likelihood , state_seq_likelihood , obs_product, **observation, 
    * norm , *state_norm, **forward1, **state_in, *backward , *cumul_backward;
  Parametric *occupancy;

# ifdef DEBUG
  register int m;
  double sum;
# endif


  // initialisations

  observation = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    observation[i] = new double[nb_state];
  }

  norm = new double[seq.length[index]];
  state_norm = new double[nb_state];

  forward1 = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    forward1[i] = new double[nb_state];
  }

  state_in = new double*[seq.length[index] - 1];
  for (i = 0;i < seq.length[index] - 1;i++) {
    state_in[i] = new double[nb_state];
  }

  backward = new double[seq.length[index] + 1];
  cumul_backward = new double[seq.length[index] + 1];

  poutput = new int*[nb_output_process];

  for (i = 0;i < nb_output_process;i++) {
    poutput[i] = seq.int_sequence[index][i + 1];
  }

# ifdef DEBUG
  double **state_sequence_probability;


  state_sequence_probability = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    state_sequence_probability[i] = new double[nb_state];
    for (j = 0;j < nb_state;j++) {
      state_sequence_probability[i][j] = 0.;
    }
  }
# endif

  // recurrence "forward"


  seq_likelihood = 0.;
  for (i = 0;i < seq.length[index];i++) {
    norm[i] = 0.;

    for (j = 0;j < nb_state;j++) {

      // calcul des probabilites d'observation

      observation[i][j] = 1.;
      for (k = 0;k < nb_output_process;k++) {
	if (sw_process[k + 1]) {
	  observation[i][j] *= sw_process[k + 1]->observation[j]->mass_computation(*poutput[k], seq.covar[index][i]);
	}
      }


      switch (state_subtype[j]) {

	// cas etat semi-markovien 

      case SEMI_MARKOVIAN : {
        if (i == 0) {
          state_norm[j] = initial[j];
        }
        else {
          state_norm[j] += state_in[i - 1][j] - forward1[i - 1][j];
        }
        state_norm[j] *= observation[i][j];

        norm[i] += state_norm[j];
        break;
      }

	// cas etat markovien

      case MARKOVIAN : {
        if (i == 0) {
          forward1[i][j] = initial[j];
        }
        else {
          forward1[i][j] = state_in[i - 1][j];
        }
        forward1[i][j] *= observation[i][j];

        norm[i] += forward1[i][j];
        break;
      }
      }
    }

    if (norm[i] > 0.) {
      for (j = 0;j < nb_state;j++) {
        switch (state_subtype[j]) {
        case SEMI_MARKOVIAN :
          state_norm[j] /= norm[i];
          break;
        case MARKOVIAN :
          forward1[i][j] /= norm[i];
          break;
        }
      }

      seq_likelihood += log(norm[i]);
    }

    else {
      seq_likelihood = D_INF;
      break;
    }

    for (j = 0;j < nb_state;j++) {

      // cas etat semi-markovien

      if (state_subtype[j] == SEMI_MARKOVIAN) {
        occupancy = nonparametric_process[0]->sojourn_time[j];
        obs_product = 1.;
        forward1[i][j] = 0.;

        if (i < seq.length[index] - 1) {
          for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
            obs_product *= observation[i - k + 1][j] / norm[i - k + 1];
            if (obs_product == 0.) {
              break;
            }

            if (k < i + 1) {
              forward1[i][j] += obs_product * occupancy->mass[k] * state_in[i - k][j];
            }

            else {
              switch (type) {
              case 'o' :
                forward1[i][j] += obs_product * occupancy->mass[k] * initial[j];
                break;
              case 'e' :
                forward1[i][j] += obs_product * forward[j]->mass[k] * initial[j];
                break;
              }
            }
          }
        }

        else {
          for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
            obs_product *= observation[i - k + 1][j] / norm[i - k + 1];
            if (obs_product == 0.) {
              break;
            }

            if (k < i + 1) {
              forward1[i][j] += obs_product * (1. - occupancy->cumul[k - 1]) * state_in[i - k][j];
            }

            else {
              switch (type) {
              case 'o' :
                forward1[i][j] += obs_product * (1. - occupancy->cumul[k - 1]) * initial[j];
                break;
              case 'e' :
                forward1[i][j] += obs_product * (1. - forward[j]->cumul[k - 1]) * initial[j];
                break;
              }
            }
          }
        }
      }
    }

    if (i < seq.length[index] - 1) {
      for (j = 0;j < nb_state;j++) {
        state_in[i][j] = 0.;
        for (k = 0;k < nb_state;k++) {
          state_in[i][j] += transition[k][j] * forward1[i][k];
        }
      }
    }

    for (j = 0;j < nb_output_process;j++) {
      poutput[j]++;
    }
  }

  if (seq_likelihood != D_INF) {

    // passes backward

#   ifdef MESSAGE
    cout << "\n";
#   endif

    for (i = 0;i < nb_state_sequence;i++) {
      j = seq.length[index] - 1;
      pstate = seq.int_sequence[index][0] + j;
      ::cumul_computation(nb_state , forward1[j] , cumul_backward);
      *pstate = cumul_method(nb_state , cumul_backward);

      do {

        // cas etat semi-markovien

        if (state_subtype[*pstate] == SEMI_MARKOVIAN) {
          occupancy = nonparametric_process[0]->sojourn_time[*pstate];
          obs_product = 1.;

          if (j < seq.length[index] - 1) {
            for (k = 1;k <= MIN(j + 1 , occupancy->nb_value - 1);k++) {
              obs_product *= observation[j - k + 1][*pstate] / norm[j - k + 1];
              if (obs_product == 0.) {
                break;
              }

              if (k < j + 1) {
                backward[k] = obs_product * occupancy->mass[k] * state_in[j - k][*pstate] /
		  forward1[j][*pstate];
              }

              else {
                switch (type) {
                case 'o' :
                  backward[k] = obs_product * occupancy->mass[k] * initial[*pstate] /
		    forward1[j][*pstate];
                  break;
                case 'e' :
                  backward[k] = obs_product * forward[*pstate]->mass[k] * initial[*pstate] /
		    forward1[j][*pstate];
                  break;
                }
              }
            }
          }

          else {
            for (k = 1;k <= MIN(j + 1 , occupancy->nb_value - 1);k++) {
              obs_product *= observation[j - k + 1][*pstate] / norm[j - k + 1];
              if (obs_product == 0.) {
                break;
              }

              if (k < j + 1) {
                backward[k] = obs_product * (1. - occupancy->cumul[k - 1]) * state_in[j - k][*pstate] /
		  forward1[j][*pstate];
              }

              else {
                switch (type) {
                case 'o' :
                  backward[k] = obs_product * (1. - occupancy->cumul[k - 1]) * initial[*pstate] /
		    forward1[j][*pstate];
                  break;
                case 'e' :
                  backward[k] = obs_product * (1. - forward[*pstate]->cumul[k - 1]) * initial[*pstate] /
		    forward1[j][*pstate];
                  break;
                }
              }
            }
          }

          ::cumul_computation(k - 1 , backward + 1 , cumul_backward);
          state_occupancy = 1 + cumul_method(k - 1 , cumul_backward);

#         ifdef DEBUG
          sum = 0.;
          for (m = 1;m < k;m++) {
            sum += backward[m];
          }
          if ((sum < 1. - DOUBLE_ERROR) || (sum > 1. + DOUBLE_ERROR)) {
            cout << "\nERROR: " << j << " " << sum << endl;
          }
#         endif

          for (k = 1;k < state_occupancy;k++) {
            pstate--;
            *pstate = *(pstate + 1);
          }
          j -= (state_occupancy - 1);

          if (j == 0) {
            break;
          }
        }

        j--;
        for (k = 0;k < nb_state;k++) {
          backward[k] = transition[k][*pstate] * forward1[j][k] / state_in[j][*pstate];
        }
        ::cumul_computation(nb_state , backward , cumul_backward);
        *--pstate = cumul_method(nb_state , cumul_backward);

#       ifdef DEBUG
        sum = 0.;
        for (k = 0;k < nb_state;k++) {
          sum += backward[k];
        }
        if ((sum < 1. - DOUBLE_ERROR) || (sum > 1. + DOUBLE_ERROR)) {
          cout << "\nERROR: " << j << " " << sum << endl;
        }
#       endif

      }
      while (j > 0);

#     ifdef DEBUG
      pstate = seq.int_sequence[index][0];
      for (j = 0;j < seq.length[index];j++) {
        state_sequence_probability[j][*pstate++]++;
      }
#     endif

#     ifdef MESSAGE
      state_seq_likelihood = Markov_switching::likelihood_computation(seq , index);

      pstate = seq.int_sequence[index][0];

      switch (format) {

      case 'a' : {
        for (j = 0;j < seq.length[index];j++) {
          os << *pstate++ << " ";
        }

        os << "  " << i + 1 << "  " << state_seq_likelihood
           << "   (" << exp(state_seq_likelihood - seq_likelihood) << ")" << endl;
        break;
      }

      case 's' : {
        for (j = 0;j < seq.length[index];j++) {
          os << *pstate++ << "\t";
        }

        os << "\t" << i + 1 << "\t" << state_seq_likelihood
           << "\t" << exp(state_seq_likelihood - seq_likelihood) << endl;
        break;
      }
      }
#     endif

    }

#   ifdef DEBUG
    if (nb_state_sequence >= 1000) {
      for (i = 0;i < seq.length[index];i++) {
        for (j = 0;j < nb_state;j++) {
          state_sequence_probability[i][j] /= nb_state_sequence;
        }
      }

      pstate = seq.int_sequence[index][0];
      for (j = 0;j < seq.length[index];j++) {
        *pstate++ = I_DEFAULT;
      }

      os << "\n" << SEQ_label[SEQL_POSTERIOR_STATE_PROBABILITY] << "\n\n";
      seq.profile_ascii_print(os , index , nb_state , state_sequence_probability ,
                              STAT_label[STATL_STATE]);
    }
#   endif

  }

  for (i = 0;i < seq.length[index];i++) {
    delete [] observation[i];
  }
  delete [] observation;

  delete [] norm;
  delete [] state_norm;

  for (i = 0;i < seq.length[index];i++) {
    delete [] forward1[i];
  }
  delete [] forward1;

  for (i = 0;i < seq.length[index] - 1;i++) {
    delete [] state_in[i];
  }
  delete [] state_in;

  delete [] backward;
  delete [] cumul_backward;

  delete [] poutput;

# ifdef DEBUG
  for (i = 0;i < seq.length[index];i++) {
    delete [] state_sequence_probability[i];
  }
  delete [] state_sequence_probability;
# endif

  return seq_likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des probabilites a posteriori des etats par l'algorithme forward-backward,
 *  des vraisemblances des sequences d'etats optimales par l'algorithme de Viterbi
 *  forward-backward, des L sequences d'etats optimales par l'algorithme de Viterbi
 *  generalise ou l'algorithme forward-backward de simulation et ecriture des resultats.
 *
 *  arguments : reference sur un objet Format_error, stream, sequences,
 *              identificateur de la sequence, type de sortie,
 *              format ('a' : ASCII, 's' : Spreadsheet),
 *              methode de calcul des sequences d'etats (algorithme de Viterbi generalise ou
 *              algorithme forward-backward de simulation), nombre de sequences d'etats.
 *
 *--------------------------------------------------------------*/
bool Markov_switching::state_profile_write(Format_error &error , ostream &os ,
                                           const Markov_switching_data &iseq , int output, 
					   int identifier , char format , int state_sequence ,
                                           int nb_state_sequence) const
{
  bool status = true;
  register int i;
  int index = I_DEFAULT;
  double seq_likelihood , max_marginal_entropy , entropy;
  Markov_switching *hsw_markov1, *hsw_markov2;
  Markov_switching_data *seq;

  error.init();

  if (identifier != I_DEFAULT) {
    for (i = 0;i < iseq.nb_sequence;i++) {
      if (identifier == iseq.identifier[i]) {
        index = i;
        break;
      }
    }

    if (i == iseq.nb_sequence) {
      status = false;
      error.update(SEQ_error[SEQR_SEQUENCE_IDENTIFIER]);
    }
  }

  if (nb_state_sequence < 2) {
    status = false;
    error.update(SEQ_error[SEQR_NB_STATE_SEQUENCE]);
  }

  if (status) {
    if (iseq.type[0] == INT_VALUE) {
      seq = new Markov_switching_data((Switching_sequence&)iseq , 0);
      seq->type[0] = STATE;
    }
    else {
      seq = new Markov_switching_data(iseq , false);
    }

    hsw_markov1 = new Markov_switching(*this , false);

    hsw_markov2 = new Markov_switching(*this, false);
    hsw_markov2->create_cumul();
    hsw_markov2->log_computation();

    for (i = 0;i < seq->nb_sequence;i++) {
      if ((index == I_DEFAULT) || (index == i)) {
        seq_likelihood = hsw_markov1->forward_backward(*seq , i , os , output , format ,
						       max_marginal_entropy , entropy);

	hsw_markov2->viterbi_forward_backward(*seq , i , os , output , format , seq_likelihood);

        switch (state_sequence) {
        case GENERALIZED_VITERBI :
	  hsw_markov2->generalized_viterbi(*seq , i , os , seq_likelihood , format ,
					   nb_state_sequence, false, true);
          break;
        case FORWARD_BACKWARD_SAMPLING :
	  hsw_markov1->forward_backward_sampling(*seq , i , os , format , nb_state_sequence);
          break;
        }
      }
    }

    delete seq;
    delete hsw_markov1;
    delete hsw_markov2;
  }

  return status;
}

/*--------------------------------------------------------------*
 *
 *  Calcul des probabilites a posteriori des etats par l'algorithme forward-backward,
 *  des vraisemblances des sequences d'etats optimales par l'algorithme de Viterbi
 *  forward-backward, des L sequences d'etats optimales par l'algorithme de Viterbi
 *  generalise ou l'algorithme forward-backward de simulation et ecriture des resultats.
 *
 *  arguments : reference sur un objet Format_error, stream,
 *              identificateur de la sequence, type de sortie,
 *              methode de calcul des sequences d'etats (algorithme de Viterbi generalise ou
 *              algorithme forward-backward de simulation), nombre de sequences d'etats.
 *
 *--------------------------------------------------------------*/
bool Markov_switching::state_profile_ascii_write(Format_error &error , ostream &os ,
                                                 int output, int identifier , 
                                                 int state_sequence , int nb_state_sequence) const
{
  bool status;


  error.init();

  if (!sw_markov_data) {
    status = false;
    error.update(STAT_error[STATR_NO_DATA]);
  }
  else {
    status = state_profile_write(error , os , *sw_markov_data , identifier ,
                                 output , 'a' , state_sequence , nb_state_sequence);
  }

  return status;
}

/*--------------------------------------------------------------*
 *
 *  Calcul des probabilites a posteriori des etats par l'algorithme forward-backward,
 *  des vraisemblances des sequences d'etats optimales par l'algorithme de Viterbi
 *  forward-backward, des L sequences d'etats optimales par l'algorithme de Viterbi
 *  generalise ou l'algorithme forward-backward de simulation et ecriture des resultats
 *  dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path, identificateur de la sequence,
 *              type de sortie, format de fichier ('a' : ASCII, 's' : Spreadsheet),
 *              methode de calcul des sequences d'etats (algorithme de Viterbi generalise ou
 *              algorithme forward-backward de simulation), nombre de sequences d'etats.
 *
 *--------------------------------------------------------------*/
bool Markov_switching::state_profile_write(Format_error &error , const char *path ,
                                           int output, int identifier , char format ,
                                           int state_sequence , int nb_state_sequence) const
{
  bool status = true;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }
  if (!sw_markov_data) {
    status = false;
    error.update(STAT_error[STATR_NO_DATA]);
  }

  if (status) {
    status = state_profile_write(error , out_file , *sw_markov_data , identifier ,
                                 output , format , state_sequence , nb_state_sequence);
  }

  return status;
}

/*--------------------------------------------------------------*
 *
 *  Calcul des probabilites a posteriori des etats par l'algorithme forward-backward,
 *  des vraisemblances des sequences d'etats optimales par l'algorithme de Viterbi
 *  forward-backward et affichage des resultats au format Gnuplot.
 *
 *  arguments : reference sur un objet Format_error, prefixe des fichiers, sequences,
 *              identificateur de la sequence, type de sortie, titre des figures.
 *
 *--------------------------------------------------------------*/
bool Markov_switching::state_profile_plot_write(Format_error &error , const char *prefix ,
                                                const Markov_switching_data &iseq , int output ,
                                                int identifier, const char *title) const
{
  bool status = true;
  register int i , j;
  int index;
  double seq_likelihood , max_marginal_entropy , entropy , state_seq_likelihood;
  Markov_switching *hsw_markov;
  Markov_switching_data *seq;
  ostringstream data_file_name[2];
  ofstream *data_out_file;


  error.init();

  for (i = 0;i < iseq.nb_sequence;i++) {
    if (identifier == iseq.identifier[i]) {
      index = i;
      break;
    }
  }

  if (i == iseq.nb_sequence) {
    status = false;
    error.update(SEQ_error[SEQR_SEQUENCE_IDENTIFIER]);
  }

  if (status) {

    // ecriture du fichier de donnees

    data_file_name[0] << prefix << 0 << ".dat";
    data_out_file = new ofstream((data_file_name[0].str()).c_str());

    if (!data_out_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_PREFIX]);
    }

    else {
      if (iseq.type[0] == INT_VALUE) {
        seq = new Markov_switching_data((Switching_sequence&)iseq , 0);
        seq->type[0] = STATE;
      }
      else {
        seq = new Markov_switching_data(iseq , false);
      }

      seq_likelihood = forward_backward(*seq , index , *data_out_file , output , 'g' ,
                                        max_marginal_entropy , entropy);
      data_out_file->close();
      delete data_out_file;

      data_file_name[1] << prefix << 1 << ".dat";
      data_out_file = new ofstream((data_file_name[1].str()).c_str());

      hsw_markov = new Markov_switching(*this , false);

      hsw_markov->create_cumul();
      hsw_markov->log_computation();
      state_seq_likelihood = hsw_markov->viterbi_forward_backward(*seq , index , *data_out_file ,
								  output , 'g' , seq_likelihood);
      data_out_file->close();
      delete data_out_file;

      // ecriture du fichier de commandes et du fichier d'impression

      for (i = 0;i < 2;i++) {
        ostringstream file_name[2];

        switch (i) {
        case 0 :
          file_name[0] << prefix << ".plot";
          break;
        case 1 :
          file_name[0] << prefix << ".print";
          break;
        }

        ofstream out_file((file_name[0].str()).c_str());

        if (i == 1) {
          out_file << "set terminal postscript" << endl;
          file_name[1] << label(prefix) << ".ps";
          out_file << "set output \"" << file_name[1].str() << "\"\n\n";
        }

        out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                 << "set title \"";
        if (title) {
          out_file << title << " - ";
        }
        switch (output) {
        case SSTATE :
          out_file << SEQ_label[SEQL_MAX_POSTERIOR_STATE_PROBABILITY] << "\"\n\n";
          break;
        case IN_STATE :
          out_file << SEQ_label[SEQL_MAX_POSTERIOR_IN_STATE_PROBABILITY] << "\"\n\n";
          break;
        case OUT_STATE :
          out_file << SEQ_label[SEQL_MAX_POSTERIOR_OUT_STATE_PROBABILITY] << "\"\n\n";
          break;
        }

        if (seq->length[index] - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }

        out_file << "plot [0:" << seq->length[index] - 1 << "] [0:1] ";
        for (j = 0;j < nb_state;j++) {
          out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                   << j + 1 << " title \"" << STAT_label[STATL_STATE] << " "
                   << j << "\" with linespoints";
          if (j < nb_state - 1) {
            out_file << ",\\";
          }
          out_file << endl;
        }

        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;

        out_file << "set title \"";
        if (title) {
          out_file << title << " - ";
        }
        switch (output) {
        case SSTATE :
          out_file << SEQ_label[SEQL_MAX_POSTERIOR_STATE_PROBABILITY] << "\"\n\n";
          break;
        case IN_STATE :
          out_file << SEQ_label[SEQL_MAX_POSTERIOR_IN_STATE_PROBABILITY] << "\"\n\n";
          break;
        case OUT_STATE :
          out_file << SEQ_label[SEQL_MAX_POSTERIOR_OUT_STATE_PROBABILITY] << "\"\n\n";
          break;
        }

        out_file << "plot [0:" << seq->length[index] - 1 << "] [0:"
                 << exp(state_seq_likelihood - seq_likelihood) << "] ";
        for (j = 0;j < nb_state;j++) {
          out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using "
                   << j + 1 << " title \"" << STAT_label[STATL_STATE] << " "
                   << j << "\" with linespoints";
          if (j < nb_state - 1) {
            out_file << ",\\";
          }
          out_file << endl;
        }

        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;

        out_file << "set title";
        if (title) {
          out_file << " \"" << title << "\"";
        }
        out_file << "\n\n";

        out_file << "plot [0:" << seq->length[index] - 1 << "] [0:" << max_marginal_entropy << "] "
                 << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                 << nb_state + 1 << " title \"" << "SEQL_CONDITIONAL_ENTROPY:"
                 << "\" with linespoints,\\" << endl;
        out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                 << nb_state + 2 << " title \"" << "SEQL_MARGINAL_ENTROPY:"
                 << "\" with linespoints" << endl;

        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;

        out_file << "set title";
        if (title) {
          out_file << " \"" << title << "\"";
        }
        out_file << "\n\n";

        out_file << "plot [0:" << seq->length[index] - 1 << "] [0:" << entropy << "] "
                 << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                 << nb_state + 3 << " title \"" << "SEQL_PARTIAL_STATE_SEQUENCE_ENTROPY"
                 << "\" with linespoints" << endl;

        if (seq->length[index] - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }

        if (i == 1) {
          out_file << "\nset terminal x11" << endl;
        }

        out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
      }

      delete seq;
      delete hsw_markov;
    }
  }

  return status;
}

/*--------------------------------------------------------------*
 *
 *  Calcul des probabilites a posteriori des etats par l'algorithme forward-backward,
 *  des vraisemblances des sequences d'etats optimales par l'algorithme de Viterbi
 *  forward-backward et affichage des resultats au format Gnuplot.
 *
 *  arguments : reference sur un objet Format_error, prefixe des fichiers,
 *              identificateur de la sequence, type de sortie, titre des figures.
 *
 *--------------------------------------------------------------*/
bool Markov_switching::state_profile_plot_write(Format_error &error , const char *prefix ,
                                                int output, int identifier , const char *title) const
{
  bool status;

  error.init();

  if (!sw_markov_data) {
    status = false;
    error.update(STAT_error[STATR_NO_DATA]);
  }
  else {
    status = state_profile_plot_write(error , prefix , *sw_markov_data , identifier ,
                                      output , title);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des L sequences d'etats optimales par l'algorithme de Viterbi generalise.
 *
 *  arguments : reference sur un objet Markov_switching_data, indice de la sequence,
 *              stream, vraisemblance des donnees, format de fichier
 *              ('a' : ASCII, 's' : Spreadsheet), nombre de sequences d'etats.
 *
 *--------------------------------------------------------------*/
double Markov_switching::generalized_viterbi(const Markov_switching_data &seq , int index ,
                                             ostream &os , double seq_likelihood ,
                                             char format , int inb_state_sequence, 
					     bool output_R, bool posterior_proba) const
{
  bool **active_cell;
  register int i , j , k , m;
  int nb_state_sequence , max_occupancy, brank , previous_rank , nb_cell , *rank , *pstate ,
    ***input_state, ***optimal_state , ***optimal_occupancy, ***input_rank, ***optimal_rank, **poutput;
  double buff , forward_max , state_seq_likelihood , likelihood_cumul , *obs_product, **observation,
    **forward1 , ***state_in;
  Parametric *occupancy;

  //# ifdef MESSAGE
  double entropy = 0.;
  //# endif


  // initialisations

  observation = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    observation[i] = new double[nb_state];
  }

  obs_product = new double[seq.length[index] + 1];

  forward1 = new double*[nb_state];
  for (i = 0;i < nb_state;i++) {
    forward1[i] = new double[inb_state_sequence];
  }

  state_in = new double**[seq.length[index] - 1];
  for (i = 0;i < seq.length[index] - 1;i++) {
    state_in[i] = new double*[nb_state];
    for (j = 0;j < nb_state;j++) {
      state_in[i][j] = new double[inb_state_sequence];
    }
  }

  rank = new int[MAX(seq.length[index] + 1 , nb_state)];

  input_state = new int**[seq.length[index] - 1];
  for (i = 0;i < seq.length[index] - 1;i++) {
    input_state[i] = new int*[nb_state];
    for (j = 0;j < nb_state;j++) {
      input_state[i][j] = new int[inb_state_sequence];
    }
  }

  optimal_state = new int**[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    optimal_state[i] = new int*[nb_state];
    for (j = 0;j < nb_state;j++) {
      optimal_state[i][j] = new int[inb_state_sequence];
    }
  }

  optimal_occupancy = new int**[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    optimal_occupancy[i] = new int*[nb_state];
    for (j = 0;j < nb_state;j++) {
      optimal_occupancy[i][j] = new int[inb_state_sequence];
    }
  }

  input_rank = new int**[seq.length[index] - 1];
  for (i = 0;i < seq.length[index] - 1;i++) {
    input_rank[i] = new int*[nb_state];
    for (j = 0;j < nb_state;j++) {
      input_rank[i][j] = new int[inb_state_sequence];
    }
  }

  optimal_rank = new int**[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    optimal_rank[i] = new int*[nb_state];
    for (j = 0;j < nb_state;j++) {
      optimal_rank[i][j] = new int[inb_state_sequence];
    }
  }

  active_cell = new bool*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    active_cell[i] = new bool[nb_state];
    for (j = 0;j < nb_state;j++) {
      active_cell[i][j] = false;
    }
  }

  poutput = new int*[nb_output_process];

  for (i = 0;i < nb_output_process;i++) {
    poutput[i] = seq.int_sequence[index][i + 1];
  }
  nb_state_sequence = 1;


# ifdef DEBUG
  double  **state_sequence_probability;
  
  state_sequence_probability = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    state_sequence_probability[i] = new double[nb_state];
    for (j = 0;j < nb_state;j++) {
      //      state_sequence_probability[i][j] = 0.;
      state_sequence_probability[i][j] = D_INF;
    }
  }
# endif

  // recurrence "forward"

  for (i = 0;i < seq.length[index];i++) {
    for (j = 0;j < nb_state;j++) {

      // calcul des probabilites d'observation

      observation[i][j] = 0.;
      for (k = 0;k < nb_output_process;k++) {
        if (sw_process[k + 1]) {
	  buff = log(sw_process[k + 1]->observation[j]->mass_computation(*poutput[k], seq.covar[index][i]));
        }

        if (buff == D_INF) {
          observation[i][j] = D_INF;
          break;
        }
        else {
          observation[i][j] += buff;
        }
      }

      switch (state_subtype[j]) {

	// cas etat semi-markovien

      case SEMI_MARKOVIAN : {
        occupancy = nonparametric_process[0]->sojourn_time[j];

        obs_product[0] = 0.;
        for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
          if (observation[i - k + 1][j] == D_INF) {
            break;
          }
          else {
            obs_product[k] = obs_product[k - 1] + observation[i - k + 1][j];
          }
        }
        max_occupancy = k - 1;

        for (k = 1;k <= max_occupancy;k++) {
          rank[k] = 0;
        }

        for (k = 0;k < nb_state_sequence;k++) {
          forward1[j][k] = D_INF;

          if (i < seq.length[index] - 1) {
            for (m = 1;m <= max_occupancy;m++) {
              if (m < i + 1) {
                buff = obs_product[m] + log(occupancy->mass[m]) + state_in[i - m][j][rank[m]];
              }

              else {
                if (rank[i + 1] == 0) {
		  buff = obs_product[m] + log(occupancy->mass[m]) + log(initial[j]);
                }

                else {
                  buff = D_INF;
                }
              }

              if (buff > forward1[j][k]) {
                forward1[j][k] = buff;
                if (m < i + 1) {
                  optimal_state[i][j][k] = input_state[i - m][j][rank[m]];
                  optimal_rank[i][j][k] = input_rank[i - m][j][rank[m]];
                }
                optimal_occupancy[i][j][k] = m;
              }
            }
          }

          else {
            for (m = 1;m <= max_occupancy;m++) {
              if (m < i + 1) {
                buff = obs_product[m] + log(occupancy->cumul[m - 1]) + state_in[i - m][j][rank[m]];
              }

              else {
                if (rank[i + 1] == 0) {
		  buff = obs_product[m] + log(occupancy->cumul[m - 1]) + log(initial[j]);
                }

                else {
                  buff = D_INF;
                }
              }

              if (buff > forward1[j][k]) {
                forward1[j][k] = buff;
                if (m < i + 1) {
                  optimal_state[i][j][k] = input_state[i - m][j][rank[m]];
                  optimal_rank[i][j][k] = input_rank[i - m][j][rank[m]];
                }
                optimal_occupancy[i][j][k] = m;
              }
            }
          }

          if (forward1[j][k] != D_INF) {
            rank[optimal_occupancy[i][j][k]]++;
          }
        }
        break;
      }

	// cas etat markovien

      case MARKOVIAN : {
        for (k = 0;k < nb_state_sequence;k++) {
          if (i == 0) {
            forward1[j][k] = log(initial[j]);
          }
          else {
            forward1[j][k] = state_in[i - 1][j][k];
            optimal_state[i][j][k] = input_state[i - 1][j][k];
            optimal_rank[i][j][k] = input_rank[i - 1][j][k];
          }
          optimal_occupancy[i][j][k] = 1;

          if (forward1[j][k] != D_INF) {
            if (observation[i][j] == D_INF) {
              forward1[j][k] = D_INF;
            }
            else {
              forward1[j][k] += observation[i][j];
            }
          }
        }
        break;
      }
      }

      for (k = nb_state_sequence;k < inb_state_sequence;k++) {
        forward1[j][k] = D_INF;
      }
    }

#   ifdef DEBUG
    cout << i << " : ";
    for (j = 0;j < nb_state;j++) {
      cout << j << " :";
      for (k = 0;k < nb_state_sequence;k++) {
        cout << " " << forward1[j][k];
        if (forward1[j][k] != D_INF) {
          cout << " " << optimal_occupancy[i][j][k];
          if (optimal_occupancy[i][j][k] < i + 1) {
            cout << " " << optimal_state[i][j][k] << " " << optimal_rank[i][j][k];
          }
        }
        cout << " |";
      }
      cout << "| ";
    }
    cout << endl;
#   endif



    if (i < seq.length[index] - 1) {
      if (nb_state_sequence < inb_state_sequence) {
        if (nb_state_sequence * nb_state < inb_state_sequence) {
          nb_state_sequence *= nb_state;
        }
        else {
          nb_state_sequence = inb_state_sequence;
        }
      }

      for (j = 0;j < nb_state;j++) {
        for (k = 0;k < nb_state;k++) {
          rank[k] = 0;
        }

        for (k = 0;k < nb_state_sequence;k++) {
          state_in[i][j][k] = D_INF;
          for (m = 0;m < nb_state;m++) {
            buff = log(transition[m][j]) + forward1[m][rank[m]];
            if (buff > state_in[i][j][k]) {
              state_in[i][j][k] = buff;
              input_state[i][j][k] = m;
              input_rank[i][j][k] = rank[m];
            }
          }

          if (state_in[i][j][k] != D_INF) {
            rank[input_state[i][j][k]]++;
          }
        }

        for (k = nb_state_sequence;k < inb_state_sequence;k++) {
          state_in[i][j][k] = D_INF;
        }
      }

#     ifdef DEBUG
      cout << i << " : ";
      for (j = 0;j < nb_state;j++) {
        cout << j << " :";
        for (k = 0;k < nb_state_sequence;k++) {
          cout << " " << state_in[i][j][k];
          if (state_in[i][j][k] != D_INF) {
            cout << " " << input_state[i][j][k] << " " << input_rank[i][j][k];
          }
          cout << " |";
        }
        cout << "| ";
      }
      cout << endl;
#     endif

    }

    for (j = 0;j < nb_output_process;j++) {
      poutput[j]++;
    }
  }

  // extraction de la vraisemblance du chemin optimal

  for (i = 0;i < nb_state;i++) {
    rank[i] = 0;
  }
  likelihood_cumul = 0.;

  for (i = 0;i < nb_state_sequence;i++) {
    pstate = seq.int_sequence[index][0] + seq.length[index] - 1;
    forward_max = D_INF;

    for (j = 0;j < nb_state;j++) {
      if (forward1[j][rank[j]] > forward_max) {
        forward_max = forward1[j][rank[j]];
        *pstate = j;
      }
    }

    if (i == 0) {
      state_seq_likelihood = forward_max;
    }

    if (forward_max == D_INF) {
      break;
    }

    // restauration

    brank = rank[*pstate];
    rank[*pstate]++;
    j = seq.length[index] - 1;

#   ifdef DEBUG
    cout << "\n" << *pstate << " " << optimal_occupancy[j][*pstate][brank] << " " << brank << " | ";
#   endif

    do {
      for (k = 0;k < optimal_occupancy[j][*pstate][brank];k++) {
        active_cell[j - k][*pstate] = true;
      }

      for (k = 0;k < optimal_occupancy[j][*pstate][brank] - 1;k++) {
        pstate--;
        *pstate = *(pstate + 1);
      }

      if (j >= optimal_occupancy[j][*pstate][brank]) {
        pstate--;
        *pstate = optimal_state[j][*(pstate + 1)][brank];
        previous_rank = optimal_rank[j][*(pstate + 1)][brank];
        j -= optimal_occupancy[j][*(pstate + 1)][brank];
        brank = previous_rank;

#       ifdef DEBUG
        cout << *pstate << " " << optimal_occupancy[j][*pstate][brank] << " " << brank << " | ";
#       endif

      }
      else {
        j -= optimal_occupancy[j][*pstate][brank];
      }
    }
    while (j >= 0);

#   ifdef DEBUG
    cout << endl;
#   endif

    likelihood_cumul += exp(forward_max);

#   ifdef DEBUG
    pstate = seq.int_sequence[index][0];
    for (j = 0;j < seq.length[index];j++) {
      /*      state_sequence_probability[j][*pstate++] += exp(forward_max - seq_likelihood); */

      if (forward_max > state_sequence_probability[j][*pstate]) {
        state_sequence_probability[j][*pstate] = forward_max;
      }
      pstate++;
    }
#   endif

    nb_cell = 0;
    for (j = 0;j < seq.length[index];j++) {
      for (k = 0;k < nb_state;k++) {
        if (active_cell[j][k]) {
          nb_cell++;
        }
      }
    }

    //#   ifdef MESSAGE
    if (i == 0) {
      os << "\n";
    }

    pstate = seq.int_sequence[index][0];



    switch (format) {
      
    case 'a' : {
      if (!output_R){
	for (j = 0;j < seq.length[index];j++) {
	  os << *pstate++ << " ";
	}
      }
      else {
	for (j = 0; j < (seq.max_length - seq.length[index]);j++){
	  os << "-1 " ;
	}
	for (j = 0; j < seq.length[index]; j++){
	  os << *pstate++ << " ";
	}
      }
        
      if(posterior_proba){         
	os << "# Probabilités a posteriori :  " << exp(forward_max-seq_likelihood) ;
      }	
      else {        
	os << 
	  "  " << i + 1 << "  " << forward_max << "   (" << exp(forward_max-seq_likelihood)
	   << "  " << likelihood_cumul / exp(seq_likelihood) << "  " << nb_cell << ")" ;
      }

      if (nb_component == nb_state) {
        os <<" "<<  SEQ_label[SEQL_STATE_BEGIN] << ": ";
	
        pstate = seq.int_sequence[index][0] + 1;
        if (seq.index_parameter) {
          for (j = 1;j < seq.length[index];j++) {
            if (*pstate != *(pstate - 1)) {
              os << seq.index_parameter[index][j] << ", ";
            }
            pstate++;
          }
        }

        else {
          for (j = 1;j < seq.length[index];j++) {
            if (*pstate != *(pstate - 1)) {
              os << j << ", ";
            }
            pstate++;
          }
        }
      }	
      os << endl;
      break;
    }
      
    case 's' : {
      for (j = 0;j < seq.length[index];j++) {
	os << *pstate++ << "\t";
      }
      
      if(posterior_proba) {
	os << "# Probabilités a posteriori :\t" << exp(forward_max - seq_likelihood) ;
      }	
      else {
	os << "\t" << i + 1 << "\t" << forward_max << "\t" << exp(forward_max - seq_likelihood)
	   << "\t" << likelihood_cumul / exp(seq_likelihood) << "\t" << nb_cell ;
      }
      break;
    }
    }
      
    entropy -= exp(forward_max - seq_likelihood) * forward_max;
    //#   endif
  }

# ifdef DEBUG
  os << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_ENTROPY] << ": " << entropy + seq_likelihood << endl;

  if (likelihood_cumul / exp(seq_likelihood) > 0.8) {
    for (i = 0;i < seq.length[index];i++) {
      for (j = 0;j < nb_state;j++) {
        if (state_sequence_probability[i][j] != D_INF) {
          state_sequence_probability[i][j] = exp(state_sequence_probability[i][j] - seq_likelihood);
        }
        else {
          state_sequence_probability[i][j] = 0.;
        }
      }
    }

    pstate = seq.int_sequence[index][0];
    for (j = 0;j < seq.length[index];j++) {
      *pstate++ = I_DEFAULT;
    }

    //    os << "\n" << SEQ_label[SEQL_POSTERIOR_STATE_PROBABILITY] << "\n\n";
    os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_STATE_PROBABILITY] << "\n\n";
    seq.profile_ascii_print(os , index , nb_state , state_sequence_probability ,
                            STAT_label[STATL_STATE]);
  }
# endif

  for (i = 0;i < seq.length[index];i++) {
    delete [] observation[i];
  }
  delete [] observation;

  delete [] obs_product;

  for (i = 0;i < nb_state;i++) {
    delete [] forward1[i];
  }
  delete [] forward1;

  for (i = 0;i < seq.length[index] - 1;i++) {
    for (j = 0;j < nb_state;j++) {
      delete [] state_in[i][j];
    }
    delete [] state_in[i];
  }
  delete [] state_in;

  delete [] rank;

  for (i = 0;i < seq.length[index] - 1;i++) {
    for (j = 0;j < nb_state;j++) {
      delete [] input_state[i][j];
    }
    delete [] input_state[i];
  }
  delete [] input_state;

  for (i = 0;i < seq.length[index];i++) {
    for (j = 0;j < nb_state;j++) {
      delete [] optimal_state[i][j];
    }
    delete [] optimal_state[i];
  }
  delete [] optimal_state;

  for (i = 0;i < seq.length[index];i++) {
    for (j = 0;j < nb_state;j++) {
      delete [] optimal_occupancy[i][j];
    }
    delete [] optimal_occupancy[i];
  }
  delete [] optimal_occupancy;

  for (i = 0;i < seq.length[index] - 1;i++) {
    for (j = 0;j < nb_state;j++) {
      delete [] input_rank[i][j];
    }
    delete [] input_rank[i];
  }
  delete [] input_rank;

  for (i = 0;i < seq.length[index];i++) {
    for (j = 0;j < nb_state;j++) {
      delete [] optimal_rank[i][j];
    }
    delete [] optimal_rank[i];
  }
  delete [] optimal_rank;

  for (i = 0;i < seq.length[index];i++) {
    delete [] active_cell[i];
  }
  delete [] active_cell;

  delete [] poutput;

# ifdef DEBUG
  for (i = 0;i < seq.length[index];i++) {
    delete [] state_sequence_probability[i];
  }
  delete [] state_sequence_probability;
# endif


  return state_seq_likelihood;

}



/*--------------------------------------------------------------*
 *
 *  Calcul des vraisemblances des sequences d'etats optimales
 *  par l'algorithme de Viterbi forward-backward.
 *
 *  arguments : reference sur un objet Markov_switching_data, indice de la sequence,
 *              stream, type de sortie, format de fichier ('a' : ASCII,
 *              's' : Spreadsheet, 'g' : Gnuplot), vraisemblance des donnees.
 *
 *--------------------------------------------------------------*/
double Markov_switching::viterbi_forward_backward(const Markov_switching_data &seq ,
                                                  int index , ostream &os , int output ,
                                                  char format , double seq_likelihood) const
{
  register int i , j , k, m;
  int *pstate, **poutput;
  double obs_product, buff , state_seq_likelihood , backward_max , **observation, 
    **forward1 , **state_in, **backward, **backward1, *auxiliary , *occupancy_auxiliary,
    **backward_output;
  Parametric *occupancy;

  // initialisations
  observation = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    observation[i] = new double[nb_state];
  }

  forward1 = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    forward1[i] = new double[nb_state];
  }

  state_in = new double*[seq.length[index] - 1];
  for (i = 0;i < seq.length[index] - 1;i++) {
    state_in[i] = new double[nb_state];
  }

  backward = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    backward[i] = new double[nb_state];
  }

  backward1 = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    backward1[i] = new double[nb_state];
  }

  auxiliary = new double[nb_state];
  occupancy_auxiliary = new double[seq.length[index] + 1];

  if (output == SSTATE) {
    backward_output = backward;
  }
  else {
    backward_output = new double*[seq.length[index]];
    for (i = 0;i < seq.length[index];i++) {
      backward_output[i] = new double[nb_state];
    }
  }

  poutput = new int*[nb_output_process];

# ifdef MESSAGE
  int *state_sequence , **input_state , **optimal_state , **optimal_forward_occupancy;

  input_state = new int*[seq.length[index] - 1];
  for (i = 0;i < seq.length[index] - 1;i++) {
    input_state[i] = new int[nb_state];
  }

  optimal_state = new int*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    optimal_state[i] = new int[nb_state];
  }

  optimal_forward_occupancy = new int*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    optimal_forward_occupancy[i] = new int[nb_state];
  }

  state_sequence = new int[seq.length[index]];
# endif

  for (i = 0;i < nb_output_process;i++) {
    poutput[i] = seq.int_sequence[index][i + 1];
  }

  // recurrence "forward"


  for (i = 0;i < seq.length[index];i++) {
    for (j = 0;j < nb_state;j++) {

      // calcul des probabilites d'observation

      observation[i][j] = 0.;
      for (k = 0;k < nb_output_process;k++) {
        if (sw_process[k + 1]) {
	  buff = log(sw_process[k + 1]->observation[j]->mass_computation(*poutput[k], seq.covar[index][i]));
	}

        if (buff == D_INF) {
          observation[i][j] = D_INF;
          break;
        }
        else {
          observation[i][j] += buff;
        }
      }

      switch (state_subtype[j]) {

	// cas etat semi-markovien

      case SEMI_MARKOVIAN : {
        occupancy = nonparametric_process[0]->sojourn_time[j];
        obs_product = 0.;
        forward1[i][j] = D_INF;

        if (i < seq.length[index] - 1) {
          for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
            if (observation[i - k + 1][j] == D_INF) {
              break;
            }
            else {
              obs_product += observation[i - k + 1][j];
            }

            if (k < i + 1) {
              buff = obs_product + log(occupancy->mass[k]) + state_in[i - k][j];
            }

            else {
	      buff = obs_product + log(occupancy->mass[k]) + log(initial[j]);
            }

            if (buff > forward1[i][j]) {
              forward1[i][j] = buff;

#             ifdef MESSAGE
              if (k < i + 1) {
                optimal_state[i][j] = input_state[i - k][j];
              }
              optimal_forward_occupancy[i][j] = k;
#             endif

            }
          }
        }

        else {
          for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
            if (observation[i - k + 1][j] == D_INF) {
              break;
            }
            else {
              obs_product += observation[i - k + 1][j];
            }

            if (k < i + 1) {
              buff = obs_product + log(occupancy->cumul[k - 1]) + state_in[i - k][j];
            }

            else {
	      buff = obs_product + log(occupancy->cumul[k - 1]) + log(initial[j]);
            }

            if (buff > forward1[i][j]) {
              forward1[i][j] = buff;

#             ifdef MESSAGE
              if (k < i + 1) {
                optimal_state[i][j] = input_state[i - k][j];
              }
              optimal_forward_occupancy[i][j] = k;
#             endif

            }
          }
        }
        break;
      }

	// cas etat markovien

      case MARKOVIAN : {
        if (i == 0) {
          forward1[i][j] = log(initial[j]);
        }
        else {
          forward1[i][j] = state_in[i - 1][j];

#         ifdef MESSAGE
          optimal_state[i][j] = input_state[i - 1][j];
#         endif

        }

#       ifdef MESSAGE
        optimal_forward_occupancy[i][j] = 1;
#       endif

        if (forward1[i][j] != D_INF) {
          if (observation[i][j] == D_INF) {
            forward1[i][j] = D_INF;
          }
          else {
            forward1[i][j] += observation[i][j];
          }
        }
        break;
      }
      }
    }

    if (i < seq.length[index] - 1) {
      for (j = 0;j < nb_state;j++) {
        state_in[i][j] = D_INF;
        for (k = 0;k < nb_state;k++) {
          buff = log(transition[k][j]) + forward1[i][k];
          if (buff > state_in[i][j]) {
            state_in[i][j] = buff;

#           ifdef MESSAGE
            input_state[i][j] = k;
#           endif

          }
        }
      }
    }

    for (j = 0;j < nb_output_process;j++) {
      poutput[j]++;
    }
  }

  // extraction de la vraisemblance du chemin optimal

# ifdef MESSAGE
  pstate = state_sequence + seq.length[index] - 1;
# endif

  state_seq_likelihood = D_INF;
  i = seq.length[index] - 1;
  for (j = 0;j < nb_state;j++) {
    if (forward1[i][j] > state_seq_likelihood) {
      state_seq_likelihood = forward1[i][j];

#     ifdef MESSAGE
      *pstate = j;
#     endif

    }
  }

  if (state_seq_likelihood != D_INF) {

#   ifdef MESSAGE
    i = seq.length[index] - 1;

    do {
      for (j = 0;j < optimal_forward_occupancy[i][*pstate] - 1;j++) {
        pstate--;
        *pstate = *(pstate + 1);
      }

      if (i >= optimal_forward_occupancy[i][*pstate]) {
        pstate--;
        *pstate = optimal_state[i][*(pstate + 1)];
        i -= optimal_forward_occupancy[i][*(pstate + 1)];
      }
      else {
        i -= optimal_forward_occupancy[i][*pstate];
      }
    }
    while (i >= 0);
#   endif

    // recurrence "backward"

    i = seq.length[index] - 1;
    for (j = 0;j < nb_state;j++) {
      backward1[i][j] = 0.;
      backward[i][j] = forward1[i][j];

      if (output == OUT_STATE) {
        backward_output[i][j] = backward[i][j];
      }
    }

    for (i = seq.length[index] - 2;i >= 0;i--) {
      for (j = 0;j < nb_state;j++) {
        switch (state_subtype[j]) {

	  // cas etat semi-markovien

        case SEMI_MARKOVIAN : {
          occupancy = nonparametric_process[0]->sojourn_time[j];
          obs_product = 0.;

          for (k = 1;k < MIN(seq.length[index] - i , occupancy->nb_value);k++) {
            if (observation[i + k][j] == D_INF) {
              break;
            }
            else {
              obs_product += observation[i + k][j];
            }

            if (k < seq.length[index] - i - 1) {
              occupancy_auxiliary[k] = backward1[i + k][j] + obs_product + occupancy->mass[k];
            }
            else {
              occupancy_auxiliary[k] = obs_product + occupancy->cumul[k - 1];
            }
          }

          auxiliary[j] = D_INF;
          for (m = k - 1;m >= 1;m--) {
            if (occupancy_auxiliary[m] > auxiliary[j]) {
              auxiliary[j] = occupancy_auxiliary[m];
            }

            // transformation des vraisemblances "semi-markoviennes" en vraisemblances "markoviennes"

            if ((auxiliary[j] != D_INF) && (state_in[i][j] != D_INF)) {
              buff = auxiliary[j] + state_in[i][j];
              if (buff > backward[i + m][j]) {
                backward[i + m][j] = buff;
              }
            }
          }
          break;
        }

	  // cas etat markovien

        case MARKOVIAN : {
          if ((backward1[i + 1][j] != D_INF) && (observation[i + 1][j] != D_INF)) {
            auxiliary[j] = backward1[i + 1][j] + observation[i + 1][j];
          }
          else {
            auxiliary[j] = D_INF;
          }
          break;
        }
        }
      }

      for (j = 0;j < nb_state;j++) {
        backward1[i][j] = D_INF;
        for (k = 0;k < nb_state;k++) {
          buff = auxiliary[k] + cumul_transition[j][k];
          if (buff > backward1[i][j]) {
            backward1[i][j] = buff;
          }
        }

        if ((backward1[i][j] != D_INF) && (forward1[i][j] != D_INF)) {
          backward[i][j] = backward1[i][j] + forward1[i][j];
        }
        else {
          backward[i][j] = D_INF;
        }
      }

      switch (output) {

      case IN_STATE : {
        for (j = 0;j < nb_state;j++) {
          switch (state_subtype[j]) {

	    // cas etat semi-markovien

          case SEMI_MARKOVIAN : {
            if ((auxiliary[j] != D_INF) && (state_in[i][j] != D_INF)) {
              backward_output[i + 1][j] = auxiliary[j] + state_in[i][j];
            }
            else {
              backward_output[i + 1][j] = D_INF;
            }
            break;
          }

	    // cas etat markovien

          case MARKOVIAN : {
            backward_output[i + 1][j] = D_INF;

            if (auxiliary[j] != D_INF) {
              for (k = 0;k < nb_state;k++) {
                if (k != j) {
                  buff = cumul_transition[k][j] + forward1[i][k];
                  if (buff > backward_output[i + 1][j]) {
                    backward_output[i + 1][j] = buff;
                  }
                }
              }

              if (backward_output[i + 1][j] != D_INF) {
                backward_output[i + 1][j] += auxiliary[j];
              }
            }
            break;
          }
          }
        }
        break;
      }

      case OUT_STATE : {
        for (j = 0;j < nb_state;j++) {
          switch (state_subtype[j]) {

	    // cas etat semi-markovien

          case SEMI_MARKOVIAN : {
            backward_output[i][j] = backward[i][j];
            break;
          }

	    // cas etat markovien

          case MARKOVIAN : {
            backward_output[i][j] = D_INF;

            if (forward1[i][j] != D_INF) {
              for (k = 0;k < nb_state;k++) {
                if (k != j) {
                  buff = auxiliary[k] + cumul_transition[j][k];
                  if (buff > backward_output[i][j]) {
                    backward_output[i][j] = buff;
                  }
                }
              }

              if (backward_output[i][j] != D_INF) {
                backward_output[i][j] += forward1[i][j];
              }
            }
            break;
          }
          }
        }
        break;
      }
      }
    }

    // cas particulier de rester dans l'etat initial

    for (i = 0;i < nb_state;i++) {
      if ((state_subtype[i] == SEMI_MARKOVIAN) && (cumul_initial[i] != D_INF)) {
        occupancy = nonparametric_process[0]->sojourn_time[i];
        obs_product = 0.;

        for (j = 1;j < MIN(seq.length[index] + 1 , occupancy->nb_value);j++) {
          if (observation[j - 1][i] == D_INF) {
            break;
          }
          else {
            obs_product += observation[j - 1][i];
          }

          if (j < seq.length[index]) {
            switch (type) {
            case 'o' :
              occupancy_auxiliary[j] = backward1[j - 1][i] + obs_product + occupancy->mass[j];
              break;
            case 'e' :
              occupancy_auxiliary[j] = backward1[j - 1][i] + obs_product + forward[i]->mass[j];
              break;
            }
          }

          else {
            switch (type) {
            case 'o' :
              occupancy_auxiliary[j] = obs_product + occupancy->cumul[j - 1];
              break;
            case 'e' :
              occupancy_auxiliary[j] = obs_product + forward[i]->cumul[j - 1];
              break;
            }
          }
        }

        auxiliary[i] = D_INF;
        for (k = j - 1;k >= 1;k--) {
          if (occupancy_auxiliary[k] > auxiliary[i]) {
            auxiliary[i] = occupancy_auxiliary[k];
          }

          // transformation des vraisemblances "semi-markoviennes" en vraisemblances "markoviennes"

          if (auxiliary[i] != D_INF) {
            buff = auxiliary[i] + cumul_initial[i];
            if (buff > backward[k - 1][i]) {
              backward[k - 1][i] = buff;
            }
          }
        }
      }

      if (output == IN_STATE) {
        backward_output[0][i] = backward[0][i];
      }
    }

    // restauration

    pstate = seq.int_sequence[index][0];

    for (i = 0;i < seq.length[index];i++) {
      backward_max = D_INF;
      for (j = 0;j < nb_state;j++) {
        if (backward[i][j] > backward_max) {
          backward_max = backward[i][j];
          *pstate = j;
        }
      }

#     ifdef MESSAGE
      if (*pstate != state_sequence[i]) {
        cout << "\nERROR: " << i << " | " << *pstate << " " << state_sequence[i] << endl;
      }
#     endif

      pstate++;
    }

    //  normalisation

    for (i = 0;i < seq.length[index];i++) {
      for (j = 0;j < nb_state;j++) {
        if (backward_output[i][j] != D_INF) {
          backward_output[i][j] = exp(backward_output[i][j] - seq_likelihood);
	  //          backward_output[i][j] = exp(backward_output[i][j] - state_seq_likelihood);
        }
        else {
          backward_output[i][j] = 0.;
        }
      }
    }

    switch (format) {

    case 'a' : {
      switch (output) {
      case SSTATE :
        os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_STATE_PROBABILITY] << "\n\n";
        break;
      case IN_STATE :
        os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_IN_STATE_PROBABILITY] << "\n\n";
        break;
      case OUT_STATE :
        os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_OUT_STATE_PROBABILITY] << "\n\n";
        break;
      }

      seq.profile_ascii_print(os , index , nb_state , backward_output ,
                              STAT_label[STATL_STATE]);

      os << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_LIKELIHOOD] << ": " << state_seq_likelihood
         << "   (" << exp(state_seq_likelihood - seq_likelihood) << ")" << endl;
      break;
    }

    case 's' : {
      switch (output) {
      case SSTATE :
        os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_STATE_PROBABILITY] << "\n\n";
        break;
      case IN_STATE :
        os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_IN_STATE_PROBABILITY] << "\n\n";
        break;
      case OUT_STATE :
        os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_OUT_STATE_PROBABILITY] << "\n\n";
        break;
      }

      seq.profile_spreadsheet_print(os , index , nb_state , backward_output ,
                                    STAT_label[STATL_STATE]);

      os << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_LIKELIHOOD] << "\t" << state_seq_likelihood
         << "\t" << exp(state_seq_likelihood - seq_likelihood) << endl;
      break;
    }

    case 'g' : {
      seq.profile_plot_print(os , index , nb_state , backward_output);
      break;
    }
    }

#   ifdef DEBUG
    if (format != 'g') {
      double ambiguity = 0.;

      pstate = seq.int_sequence[index][0];
      //      if (output == SSTATE) {
      for (i = 0;i < seq.length[index];i++) {
        for (j = 0;j < nb_state;j++) {
          if (j != *pstate) {
            ambiguity += backward_output[i][j];
          }
        }
        pstate++;
      }
      ambiguity *= exp(seq_likelihood - state_seq_likelihood);
      /*      }

      else {
      for (i = 0;i < seq.length[index];i++) {
      for (j = 0;j < nb_state;j++) {
      if ((backward[i][j] != D_INF) && (j != *pstate)) {
      ambiguity += exp(backward[i][j] - state_seq_likelihood);
      }
      }
      pstate++;
      }
      } */

      switch (format) {
      case 'a' :
        os << "\n" << SEQ_label[SEQL_AMBIGUITY] << ": " << ambiguity
           << " (" << ambiguity / seq.length[index] << ")" << endl;
        break;
      case 's' :
        os << "\n" << SEQ_label[SEQL_AMBIGUITY] << "\t" << ambiguity
           << "\t" << ambiguity / seq.length[index] << "\t" << endl;
        break;
      }
    }
#   endif

  }

  for (i = 0;i < seq.length[index];i++) {
    delete [] observation[i];
  }
  delete [] observation;

  for (i = 0;i < seq.length[index];i++) {
    delete [] forward1[i];
  }
  delete [] forward1;

  for (i = 0;i < seq.length[index] - 1;i++) {
    delete [] state_in[i];
  }
  delete [] state_in;

  for (i = 0;i < seq.length[index];i++) {
    delete [] backward[i];
  }
  delete [] backward;

  for (i = 0;i < seq.length[index];i++) {
    delete [] backward1[i];
  }
  delete [] backward1;

  delete [] auxiliary;
  delete [] occupancy_auxiliary;

  if (output != SSTATE) {
    for (i = 0;i < seq.length[index];i++) {
      delete [] backward_output[i];
    }
    delete [] backward_output;
  }

  delete [] poutput;

# ifdef MESSAGE
  for (i = 0;i < seq.length[index] - 1;i++) {
    delete [] input_state[i];
  }
  delete [] input_state;

  for (i = 0;i < seq.length[index];i++) {
    delete [] optimal_state[i];
  }
  delete [] optimal_state;

  for (i = 0;i < seq.length[index];i++) {
    delete [] optimal_forward_occupancy[i];
  }
  delete [] optimal_forward_occupancy;

  delete [] state_sequence;
# endif

  return state_seq_likelihood;
 
}


//---------------------------------------------------------------
// Calcul des logarithmes des paramètres d'un Markov_switching
//---------------------------------------------------------------
void Markov_switching::log_computation ()
{
  register int i , j;
  double *pcumul;
  Parametric *occupancy;

  Chain::log_computation();

  for (i = 0;i < nb_state;i++) {
    if (state_subtype[i] == SEMI_MARKOVIAN) {
      occupancy = nonparametric_process[0]->sojourn_time[i];

      if (occupancy->mass[occupancy->offset] > 0.) {
        ::log_computation(occupancy->nb_value , occupancy->mass , occupancy->mass);

        pcumul = occupancy->cumul;
        for (j = 0;j < occupancy->nb_value;j++) {
          *pcumul = 1. - *pcumul;
          pcumul++;
        }
        ::log_computation(occupancy->nb_value , occupancy->cumul , occupancy->cumul);

        if (type == 'e') {
          ::log_computation(forward[i]->nb_value , forward[i]->mass , forward[i]->mass);

          pcumul = forward[i]->cumul;
          for (j = 0;j < forward[i]->nb_value;j++) {
            *pcumul = 1. - *pcumul;
            pcumul++;
          }
          ::log_computation(forward[i]->nb_value , forward[i]->cumul , forward[i]->cumul);
        }
      }
    }
  }

  //   for (i = 1;i <= nb_output_process;i++) {
  //     for (j = 0;j < nb_state;j++) {
  //       ::log_computation(sw_process[i]->nb_value , sw_process[i]->observation[j]->density ,
  // 			sw_process[i]->observation[j]->cumul);
  //     }
  //   }

}



/*--------------------------------------------------------------*
 *
 *  Simulation par une chaine de Markov.
 *
 *  arguments : reference sur un objet Format_error,
 *              histogramme des longueurs des sequences,
 *              flag sur le calcul des lois de comptage,
 *              flag calcul d'une divergence de Kullback-Leibler.
 *
 *--------------------------------------------------------------*/
Markov_switching_data* Markov_switching::simulation(gsl_rng *r, Format_error &error , const Histogram &hlength ,
						    int inb_covariable, int iconstant,
						    double **iregression,
						    double ***icovar,
						    bool counting_flag , bool divergence_flag) const
{
  bool status = true;
  register int i , j , k, m;
  int cumul_length , occupancy, *pstate , **poutput;
  double *pcovar, *ccovar;
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

    seq = new Markov_switching_data(nb_output_process + 1 , hlength , inb_covariable, iconstant, false);
    seq->type[0] = STATE;

    seq->sw_markov = new Markov_switching(*this , false);
    sw_markov = seq->sw_markov;

    sw_markov->create_cumul();
    sw_markov->cumul_computation();

    if (sw_markov->nb_output_process > 0) {
      poutput = new int*[sw_markov->nb_output_process];
    }


    for (i = 0; i < sw_markov->nb_state; i++) {
      for( k = 0; k < nb_output_process; k++){
	if (!iconstant){
	  sw_markov->sw_process[k + 1]->observation[i]->nb_covariable = inb_covariable +1;
	  nb_covar = inb_covariable + 1;
	}
	else {
	  sw_markov->sw_process[k + 1]->observation[i]->nb_covariable = inb_covariable;
	  nb_covar = inb_covariable;
	}
	sw_markov->sw_process[k + 1]->observation[i]->Set_regression(iregression[i]);
      }
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
	  //  cout<<seq->covar[i][j][k]<<" ";
	}
      }
    }

    for (i = 0;i < seq->nb_sequence;i++) {
      pstate = seq->int_sequence[i][0];
      *pstate = cumul_method(sw_markov->nb_state , sw_markov->cumul_initial);

      for (j = 0;j < sw_markov->nb_output_process;j++) {
        poutput[j] = seq->int_sequence[i][j + 1];
      }

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
          for (m = 0;m < sw_markov->nb_output_process;m++) {
            if (sw_markov->sw_process[m + 1]) {
              *poutput[m]++ = sw_markov->sw_process[m + 1]->observation[*pstate]->simulation(r,seq->covar[i][j]);
	    }
	  }
	}
      }
      while (j < seq->length[i]);
    }

    sw_markov->remove_cumul();

    if (sw_markov->nb_output_process > 0) {
      delete [] poutput;
    }

    // extraction des caracteristiques des sequences simulees

    for (i = 0;i < seq->nb_variable;i++) {
      seq->max_value_computation(i);
      seq->build_marginal_histogram(i);
    }


    seq->build_transition_count(*sw_markov);
    seq->build_observation_histogram();
    seq->build_characteristic();



    if (!divergence_flag) {
      sw_markov->characteristic_computation(*seq , counting_flag);


      // calcul de la vraisemblance

      //  seq->likelihood = sw_markov->likelihood_computation(*seq); 

#     ifdef MESSAGE
      cout << "\n" << STAT_label[STATL_LIKELIHOOD] << ": " << seq->likelihood
           << " | " << sw_markov->likelihood_computation(*seq, I_DEFAULT, true) << endl;
#     endif

    }

  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par une chaine de Markov .
 *
 *  arguments : reference sur un objet Format_error,
 *              nombre et longueur des sequences,
 *              flag sur le calcul des lois de comptage.
 *
 *--------------------------------------------------------------*/
Markov_switching_data* Markov_switching::simulation(gsl_rng *r, Format_error &error , int nb_sequence ,
						    int length , int inb_covariable, int iconstant, 
						    double **iregression,
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

    seq = simulation(r, error , hlength , inb_covariable, iconstant, 
		     iregression, icovar, counting_flag, false);
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par une chaine de Markov.
 *
 *  arguments : reference sur un objet Format_error, nombre de sequences,
 *              reference sur un objet Markovian_sequences,
 *              flag sur le calcul des lois de comptage.
 *
 *--------------------------------------------------------------*/
Markov_switching_data* Markov_switching::simulation(gsl_rng *r, Format_error &error , int nb_sequence ,
						    const Switching_sequence &iseq ,
						    int inb_covariable, int iconstant, 
						    double **iregression,
						    double ***icovar,
						    bool counting_flag) const
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

    seq = simulation(r, error , *hlength , inb_covariable, iconstant,
		     iregression, icovar, counting_flag, false);
    delete hlength;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par une chaine de Markov cachee.
 *
 *  arguments : reference sur un objet Format_error,
 *              histogramme des longueurs des sequences,
 *              flag sur le calcul des lois de comptage,
 *              flag calcul d'une divergence de Kullback-Leibler.
 *
 *-------------------------------------------------------------*/
Markov_switching_data* Markov_switching::simulation(gsl_rng *r,Format_error &error ,
                                                    const Histogram &hlength ,
						    int inb_covariable, int iconstant, 
						    double **iregression,
						    double ***icovar,
                                                    bool counting_flag ,
                                                    bool divergence_flag, bool hidden) const
{
  Switching_sequence *observ_seq;
  Markov_switching_data *seq;

  seq = Markov_switching::simulation(r, error , hlength , inb_covariable, iconstant,
				     iregression, icovar, counting_flag , divergence_flag);

  if ((seq) && (!divergence_flag)) {
    observ_seq = seq->remove_variable_1();
    seq->hidden_likelihood = likelihood_computation(*observ_seq,0, I_DEFAULT, true);
    delete observ_seq;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par une chaine de Markov cachee.
 *
 *  arguments : reference sur un objet Format_error,
 *              nombre et longueur des sequences.
 *
 *--------------------------------------------------------------*/
Markov_switching_data* Markov_switching::simulation(gsl_rng *r, Format_error &error ,
                                                    int nb_sequence , int length ,
						    int inb_covariable, int iconstant, 
						    double **iregression,
						    double ***icovar,
                                                    bool counting_flag, bool hidden) const
{
  Switching_sequence *observ_seq;
  Markov_switching_data *seq;


  seq = Markov_switching::simulation(r, error , nb_sequence , length , inb_covariable, iconstant,
				     iregression, icovar, counting_flag);

  if (seq) {
    observ_seq = seq->remove_variable_1();
    seq->hidden_likelihood = likelihood_computation(*observ_seq,0, I_DEFAULT, true);
    delete observ_seq;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par une chaine de Markov cachee.
 *
 *  arguments : reference sur un objet Format_error, nombre de sequences,
 *              reference sur un objet Markovian_sequences.
 *
 *--------------------------------------------------------------*/
Markov_switching_data* Markov_switching::simulation(gsl_rng *r, Format_error &error ,
                                                    int nb_sequence ,
                                                    const Switching_sequence &iseq ,
						    int inb_covariable, int iconstant, 
						    double **iregression,
						    double ***icovar,
                                                    bool counting_flag, bool hidden) const

{
  Switching_sequence *observ_seq;
  Markov_switching_data *seq;

  seq = Markov_switching::simulation(r, error , nb_sequence , iseq , inb_covariable, iconstant,
				     iregression, icovar,  counting_flag);

  if (seq) {
    observ_seq = seq->remove_variable_1();
    seq->hidden_likelihood = likelihood_computation(*observ_seq, 0, I_DEFAULT, true);
    delete observ_seq;
  }

  return seq;
}

// //////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////
// ////                      CAS MULTINOMIAL                         ////
// ////////////////////////////////////////////////////////////////////// 
// //////////////////////////////////////////////////////////////////////


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une chaine de Markov switching
 *  a partir d'un echantillon de sequences par l'algorithme EM.
 *
 *  arguments : reference sur un objet Format_error, stream, chaine de Markov switching initiale,
 *              flags sur le calcul des lois de comptage et sur le calcul des sequences
 *              d'etats optimales, nombre d'iterations.
 *
 *--------------------------------------------------------------*/
Markov_switching* Switching_sequence::markov_switching_estimation_multi(Format_error &error , ostream &os ,
									const Markov_switching &isw_markov , int nb_cat,
									int estimator, bool output_R, bool counting_flag, 
									bool state_sequence,
									int nb_iter, int mean_computation) const
{
  bool status;
  register int i , j , k , m, n, cat;
  int max_nb_value , iter, nb_likelihood_decrease, offset, nb_value, *occupancy_nb_value,
    *censored_occupancy_nb_value, **poutput;
  double likelihood = D_INF , previous_likelihood , occupancy_likelihood, observation_likelihood , min_likelihood ,
    obs_product, buff, sum, occupancy_mean, **observation, *norm, *state_norm, **forward, **state_in, ***backward, **backward1,
    *auxiliary, *reestim, *ofrequency, *lfrequency, *occupancy_survivor, *censored_occupancy_survivor;
  
  double *temp;

  Chain_reestimation<double> *chain_reestim;
  Parametric *occupancy;
  Reestimation<double> **occupancy_reestim, **length_bias_reestim, **censored_occupancy_reestim;
  Histogram *hoccupancy;
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

    hsw_markov = new Markov_switching(isw_markov , false, I_DEFAULT);// (int)(max_length * SAMPLE_NB_VALUE_COEFF));

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
    for (i = 0;i < hsw_markov->nb_output_process;i++) {
      if ((hsw_markov->sw_process[i + 1]) && (max_nb_value < marginal[i]->nb_value)) {
        max_nb_value = marginal[i]->nb_value;
      }
    }

    poutput = new int*[nb_variable];
    
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
	  poutput[j] = int_sequence[i][j];
	}

	// recurrence "forward"
	
	for (j = 0;j < length[i];j++) {
          norm[j] = 0.;

          for (k = 0;k < hsw_markov->nb_state;k++) {

            // calcul des probabilites d'observation

            observation[j][k] = 1.;
            for (m = 0;m < hsw_markov->nb_output_process;m++) {
	      observation[j][k] *= hsw_markov->sw_process[m + 1]->observation[k]->mass_computation(*poutput[m], covar[i][j]);
// 	      if((i==40)& (j==0)){
// 		cout<<"etat: "<<k<<", cat: "<<*poutput[m]<<", beta: "<<covar[i][j][0]<<" , "<<covar[i][j][1]<<", masse: "<<
// 		  hsw_markov->sw_process[m + 1]->observation[k]->mass_computation(*poutput[m], covar[i][j])<<endl;
// 	      } 
// 	      if((i==40)& (j==1)){
// 		cout<<"etat: "<<k<<", cat: "<<*poutput[m]<<", beta: "<<covar[i][j][0]<<" , "<<covar[i][j][1]<<", masse: "<<
// 		  hsw_markov->sw_process[m + 1]->observation[k]->mass_computation(*poutput[m], covar[i][j])<<endl;
// 	      } 
// 	      if((i==40)& (j==8)){
// 		cout<<"etat: "<<k<<", cat: "<<*poutput[m]<<", beta: "<<covar[i][j][0]<<" , "<<covar[i][j][1]<<", masse: "<<
// 		  hsw_markov->sw_process[m + 1]->observation[k]->mass_computation(*poutput[m], covar[i][j])<<endl;
// 	      } 
// 	      if((i==40)& (j==6)){
// 		cout<<"etat: "<<k<<", cat: "<<*poutput[m]<<", beta: "<<covar[i][j][0]<<" , "<<covar[i][j][1]<<", masse: "<<
// 		  hsw_markov->sw_process[m + 1]->observation[k]->mass_computation(*poutput[m], covar[i][j])<<endl;
// 	      } 
// 	      // cout<< hsw_markov->sw_process[m + 1]->observation[k]->mass_computation(*poutput[m], covar[i][j])<<endl;
	    }

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







      //      if((iter==0) || ( likelihood - previous_likelihood) / -likelihood > MARKOV_SWITCHING_LIKELIHOOD_DIFF){
	
	// estimation des parametres de regression pour chacun des états 
	
	
	for (k = 0; k < hsw_markov->nb_state; k++){
	  
	  double **beta_kk;
	  gsl_matrix *beta_k = gsl_matrix_calloc(nb_cat, nb_covariable);
	  
	  beta_kk = new double*[nb_cat];
	  for (cat = 0; cat < nb_cat; cat++){
	    beta_kk[cat] = new double[nb_covariable];
	  }
	  
	  	  
	  for(cat = 0;cat< (nb_cat-1); cat++){
	    
	    gsl_vector *add1 = gsl_vector_calloc(nb_covariable);
	    gsl_matrix *add2 = gsl_matrix_calloc(nb_covariable, nb_covariable);
	    
	    int s;
	    gsl_permutation *p = gsl_permutation_alloc(nb_covariable);
	    gsl_matrix *trans = gsl_matrix_alloc(nb_covariable, nb_covariable);
	    gsl_matrix *inv_add2 = gsl_matrix_alloc(nb_covariable, nb_covariable); 
	    
	    gsl_vector *temp6 = gsl_vector_calloc(nb_covariable);
	    
	    gsl_vector *beta_catk = gsl_vector_calloc(nb_covariable);
	    
	    //	  cout<<"paramètres de regression ";
	    for (j = 0; j < nb_covariable; j++){
	      gsl_vector_set(beta_catk, j, hsw_markov->sw_process[1]->observation[k]->regmult[cat][j]);
	      //	    cout<<gsl_vector_get(beta_catk,j)<<" ";
	    }
	    //	  cout<<endl<<endl;
	    
	    
	    
	    for(i = 0; i < nb_sequence; i++){
	      
	      gsl_matrix *V_icatk = gsl_matrix_calloc(length[i],length[i]);
	      gsl_vector *tmp_icatk = gsl_vector_calloc(length[i]);
	      gsl_vector *mu_icatk = gsl_vector_calloc(length[i]);
	      
	      gsl_matrix *X_i = gsl_matrix_calloc(length[i],nb_covariable);
	      
	      gsl_matrix *L_ik = gsl_matrix_calloc(length[i],length[i]);
	      gsl_matrix *temp_mat = gsl_matrix_calloc(length[i], hsw_markov->nb_state);
	      gsl_vector *temp_vect = gsl_vector_calloc(length[i]);
	      
	      gsl_matrix *temp1 = gsl_matrix_calloc(nb_covariable, length[i]);
	      gsl_vector *temp2 = gsl_vector_calloc(nb_covariable);

	      
	      gsl_matrix *temp3 = gsl_matrix_calloc(nb_covariable, length[i]);
	      gsl_matrix *temp4 = gsl_matrix_calloc(nb_covariable, length[i]);
	      gsl_matrix *temp5 = gsl_matrix_calloc(nb_covariable, nb_covariable);
	      

	      
	      for (m=0; m< length[i]; m++){
		gsl_matrix_set(V_icatk, m, m, hsw_markov->sw_process[1]->observation[k]->
			       variance_computation_multi(int_sequence[i][0][m],covar[i][m]));
	      }
	      
	      //cout<<gsl_matrix_get(V_icatk,0,0)<<"  ";
	      
	      for (m=0; m< length[i]; m++){
		gsl_vector_set(tmp_icatk, m, hsw_markov->sw_process[1]->observation[k]->
			       mean_computation_multi(int_sequence[i][0][m],covar[i][m]));
		
		
		if (int_sequence[i][0][m] == cat){
		  gsl_vector_set(mu_icatk, m, 1. - gsl_vector_get(tmp_icatk,m));
		}
		else {
		  gsl_vector_set(mu_icatk, m, 0. - gsl_vector_get(tmp_icatk,m));
		}
		
	      }
	      
	      convert_matrix_double(covar[i], length[i], nb_covariable, X_i);
	      
	      convert_matrix_double(backward[i], length[i], hsw_markov->nb_state, temp_mat);
	      gsl_matrix_get_col(temp_vect, temp_mat, k);
	      convert_matrix_diag_double(temp_vect, length[i], L_ik);
	      
	      //somme des tXLmu
	      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., X_i, L_ik, 0., temp1);
	      gsl_blas_dgemv(CblasNoTrans, 1., temp1, mu_icatk, 0., temp2);
	      gsl_vector_add(add1, temp2); 
	      
	      
	      //somme des tXLVX
	      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., X_i, L_ik, 0., temp3);
	      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., temp3, V_icatk, 0., temp4);
	      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., temp4, X_i, 0., temp5);
	      gsl_matrix_add(add2, temp5); 

	      
	      gsl_matrix_free(temp5);
	      gsl_matrix_free(temp4);
	      gsl_matrix_free(temp3);
	      gsl_vector_free(temp2);
	      gsl_matrix_free(temp1);
	      gsl_vector_free(temp_vect);
	      gsl_matrix_free(temp_mat);
	      gsl_matrix_free(L_ik);
	      gsl_matrix_free(X_i);
	      gsl_vector_free(mu_icatk);
	      gsl_vector_free(tmp_icatk);
	      gsl_matrix_free(V_icatk);
	      
	    }
	    
	    
	    gsl_matrix_memcpy(trans, add2);
	    gsl_linalg_LU_decomp(trans, p, &s);
	    gsl_linalg_LU_invert(trans, p, inv_add2);
	    
	    gsl_blas_dgemv(CblasNoTrans, 1., inv_add2, add1, 0., temp6);
	    
	    
	    //	  cout<<gsl_vector_get(temp6,0)<<"  ";
	    
	    gsl_vector_add(beta_catk, temp6);
	    gsl_matrix_set_row(beta_k, cat, beta_catk);
	    
	    
	    gsl_vector_free(beta_catk);
	    gsl_vector_free(temp6);
	    gsl_matrix_free(inv_add2);
	    gsl_matrix_free(trans);
	    gsl_permutation_free(p);
	    gsl_matrix_free(add2);
	    gsl_vector_free(add1);
	    
	  }
	  
	  
	  
	  for (cat = 0; cat < nb_cat; cat++){
	    for (j = 0; j < nb_covariable; j++){
	      beta_kk[cat][j] = gsl_matrix_get(beta_k, cat, j);
	    }
	  }
	  
	  hsw_markov->sw_process[1]->observation[k]->Set_regmult(beta_kk);
	
	  
#   ifdef DEBUG
	  cout<<"beta pour l'état "<<k<<" :  "<<endl;
	  for (cat = 0; cat < nb_cat; cat++){
	    cout<<"Category "<< cat <<"  ";
	    for(i = 0; i<nb_covariable; i++){
	      cout<<beta_kk[cat][i]<<" |  " << hsw_markov->sw_process[1]->observation[k]->regmult[cat][i]<<"   ";
	    }
	    cout<<endl;
	  }
	  cout<<endl;
#endif
	  
	
	  for (cat = 0; cat < nb_cat; cat++){
	    delete beta_kk[cat];
	  }
	  delete [] beta_kk;
	  
	  gsl_matrix_free(beta_k);
	  
	}
	
	//      }







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


	//cout<<*chain_reestim<<endl;


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


      //if(iter>1){
      //      os<< likelihood-previous_likelihood <<", "<<endl;
      //      }

      
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

//       cout<<previous_likelihood<<endl;
//       cout<<(likelihood-previous_likelihood)/(-likelihood)<<endl;

    } 
    while //(iter<5);
      ((iter==1)|| /*(likelihood != D_INF) &&*/ ((likelihood - previous_likelihood) / -likelihood > MARKOV_SWITCHING_LIKELIHOOD_DIFF)
       && ((nb_iter == I_DEFAULT) && (iter < MARKOV_SWITCHING_NB_ITER) || ((nb_iter != I_DEFAULT) && (iter < nb_iter)))); 
    // revoir conditions d'erret


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
	
      for (i = 0;i < hsw_markov->nb_output_process;i++) {
	for (j = 0;j < hsw_markov->nb_state;j++) {
	  // hsw_markov->sw_process[i + 1]->observation[j]->mean_computation();// A REVOIR
	  hsw_markov->sw_process[i + 1]->observation[j]->parametric_variance_computation();
	  // hsw_markov->sw_process[i + 1]->observation[j]->min_value_computation();
	  //hsw_markov->sw_process[i + 1]->observation[j]->max_value_computation();
	  hsw_markov->sw_process[i + 1]->observation[j]->nb_value_computation();
	  // hsw_markov->sw_process[i + 1]->observation[j]->density_cont_computation();
	  // hsw_markov->sw_process[i + 1]->observation[j]->cumul_cont_computation();
	  
	  
	  // if (observation_likelihood == D_INF) {
	  //   min_likelihood = D_INF;
	  // }
	  // else {
	  //   hsw_markov->sw_process[i + 1]->observation[j]->computation( OBSERVATION_THRESHOLD);
	  // }
	}
      }
       
    }

#   ifdef DEBUG 
    cout<<"AFFICHAGE DES PARAMETRES DE REGRESSION: "<<endl;
    for(i = 0; i< hsw_markov->nb_output_process; i++){
      for (j = 0; j< hsw_markov->nb_state; j++) {
	for(k=0; k<nb_covariable; k++) {
	  cout << hsw_markov->sw_process[i + 1]->observation[j]->regression[k]<< "   ";
	}
	cout<<endl;
      }
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
	
	for (i = 1;i <= hsw_markov->nb_output_process;i++) {
	  if ((hsw_markov->sw_process[i]) && (hsw_markov_data->characteristics[i])) {
	    delete hsw_markov_data->characteristics[i];
	    hsw_markov_data->characteristics[i] = 0;
	  }
	}

	hsw_markov->create_cumul();	
	
	//hsw_markov->log_computation(); // PROBLEME de memoire

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
							    true, false);// attention a cause de log_computation()
	}
	else{
	  hsw_markov_data->likelihood = hsw_markov->viterbi(*hsw_markov_data, hsw_markov_data->posterior_probability, 
							    false, true);// attention a cause de log_computation()
	}

	hsw_markov->remove_cumul();

	hsw_markov_data->max_value[0] = hsw_markov->nb_state - 1;
	hsw_markov_data->build_marginal_histogram(0);
	hsw_markov_data->build_characteristic(0);

	hsw_markov_data->build_transition_count(*hsw_markov);
	//	hsw_markov_data->build_observation_histogram(); // A REVOIR


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
	//   << " | " << hsw_markov->Markov_switching::likelihood_computation(*hsw_markov_data) << endl;  // A REVOIR
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
	hsw_markov_data->state_variable_init(INT_VALUE);
	
	for (i = 1;i <= hsw_markov->nb_output_process;i++) {
	  if ((hsw_markov->sw_process[i]) && (hsw_markov_data->characteristics[i - 1])) {
	    delete hsw_markov_data->characteristics[i - 1];
	    hsw_markov_data->characteristics[i - 1] = 0;
	  }
	}
      }
    
      
      // calcul de la vraisemblance et des lois caracteristiques du modele

      hsw_markov_data->hidden_likelihood = hsw_markov->likelihood_computation(*this, hsw_markov_data->posterior_probability,
									      I_DEFAULT, true);
      
      
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
 *  Estimation des parametres d'une chaine de Markov cachee
 *  a partir d'un echantillon de sequences par l'algorithme SEM/MCEM.
 *
 *  arguments : reference sur un objet Format_error, stream, chaine de Markov cachee initiale,
 *              parametres pour le nombre de sequences d'etats simulees, flags sur le calcul
 *              des lois de comptage et sur le calcul des sequences d'etats optimales,
 *              nombre d'iterations.
 *
 *--------------------------------------------------------------*/
Markov_switching* Switching_sequence::markov_switching_stochastic_estimation_multi(Format_error &error , std::ostream &os ,
										   const Markov_switching &ihsw_markov , int nb_cat,
										   int min_nb_state_sequence,
										   int max_nb_state_sequence,
										   bool output_R, double parameter,
										   int estimator, bool counting_flag,
										   bool state_sequence, int nb_iter) const
{
  bool status;
  register int i , j , k ,l, m , n, cat;
  int max_nb_value , iter , nb_state_sequence , state_occupancy , nb_likelihood_decrease ,
    *occupancy_nb_value , *state_seq , *pstate, **poutput;
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

  int kk;

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


    poutput = new int*[nb_variable];


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
          poutput[j] = int_sequence[i][j];
        }


        // recurrence "forward"

	for (j = 0;j < length[i];j++) {
          norm[j] = 0.;

          for (k = 0;k < hsmarkov->nb_state;k++) {

            // calcul des probabilites d'observation

            observation[j][k] = 1.;
            for (m = 0;m < hsmarkov->nb_output_process;m++) {
	      observation[j][k] *= hsmarkov->sw_process[m + 1]->observation[k]->mass_computation(*poutput[m], covar[i][j]);
            }

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
            poutput[m] = int_sequence[i][m] + k;
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

	      kk = k;
	      k -= (state_occupancy - 1);

	      for(m = k; m <= kk; m++){
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

      if(( likelihood - previous_likelihood) / -likelihood > MARKOV_SWITCHING_LIKELIHOOD_DIFF){

	// estimation des parametres de regression pour chacun des états 
      
	for (k = 0; k < hsmarkov->nb_state; k++) {
	  
	  double **beta_kk;
	  gsl_matrix *beta_k = gsl_matrix_calloc(nb_cat, nb_covariable);
	  
	  beta_kk = new double *[nb_cat];
	  for (cat = 0; cat < nb_cat; cat++){
	    beta_kk[cat] = new double[nb_covariable];
	  }
	  
	  gsl_matrix_set_zero(beta_k);
	  
	  for (cat = 0; cat < (nb_cat-1); cat++){
	    
	    int s;
	    gsl_permutation *p = gsl_permutation_alloc(nb_covariable);
	    gsl_matrix *trans = gsl_matrix_alloc(nb_covariable, nb_covariable);
	    gsl_matrix *inv_temp_add = gsl_matrix_alloc(nb_covariable, nb_covariable);
	    gsl_matrix *temp_add = gsl_matrix_calloc(nb_covariable, nb_covariable);
	    gsl_vector *add_vect = gsl_vector_calloc(nb_covariable); 
	    gsl_vector *beta_temp = gsl_vector_calloc(nb_covariable);
	    gsl_vector *beta_catk = gsl_vector_calloc(nb_covariable);
	    
	    gsl_matrix_set_zero(temp_add);
	    gsl_vector_set_zero(add_vect);
	    
	    for (j = 0; j < nb_covariable; j++) {
	      gsl_vector_set(beta_catk, j, hsmarkov->sw_process[1]->observation[k]->regmult[cat][j]);
	    }
	    
	  
	    // pour l'instant, cet algorithme ne marchera que pour une seule variable j
	    for (j = 0 ; j < hsmarkov->nb_output_process; j++) {
	      
	      for (i = 0; i < nb_sequence; i++) {
	      
		gsl_matrix *X_i = gsl_matrix_calloc(length[i], nb_covariable);
		gsl_matrix *indic_ik = gsl_matrix_calloc(length[i], length[i]);
		gsl_matrix *temp_mat = gsl_matrix_calloc(length[i], hsmarkov->nb_state);
		gsl_matrix *temp_mult2 = gsl_matrix_calloc(nb_covariable, nb_covariable);
		gsl_matrix *temp_mult1 = gsl_matrix_calloc(nb_covariable, length[i]);
		gsl_matrix *temp_mult = gsl_matrix_calloc(nb_covariable, length[i]);
		
		convert_matrix_double(covar[i], length[i], nb_covariable, X_i);
		
		for (m = 0; m < nb_state_sequence; m++){
		
		  gsl_matrix *W_icatk = gsl_matrix_calloc(length[i], length[i]);
		  gsl_vector *mu_icatk = gsl_vector_calloc(length[i]);
		  gsl_vector *delta = gsl_vector_calloc(length[i]);
		  gsl_vector *temp_vect = gsl_vector_calloc(length[i]);
		  gsl_vector *temp_mult3 = gsl_vector_calloc(nb_covariable);
		  
		  
		  for (l = 0; l < length[i]; l++) {
		    gsl_matrix_set(W_icatk, l, l, hsmarkov->sw_process[j + 1]->observation[k]->
				   variance_computation_multi(int_sequence[i][j][l],covar[i][l]));
		    gsl_vector_set(mu_icatk, l, hsmarkov->sw_process[j + 1]->observation[k]->
				   mean_computation_multi(int_sequence[i][j][l], covar[i][l]));
		    
		    if (int_sequence[i][j][l] == cat){
		      gsl_vector_set(delta, l, 1.);
		    }		   
		    else{
		      gsl_vector_set(delta, l, 0.);
		    }
		  }
		  
		  gsl_vector_sub(delta, mu_icatk);
		  
		  convert_matrix_double(indic[i][m], length[i], hsmarkov->nb_state, temp_mat);
	      
		  gsl_matrix_get_col(temp_vect, temp_mat, k);
		  convert_matrix_diag_double(temp_vect, length[i], indic_ik);
		
		  // 	      for (int v=0; v < length[i]; v++){
		  // 		  cout<< gsl_matrix_get(indic_ik,v,v)<<"  ";
		  // 	      }
		  
		  //somme des tXIWX
		  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., X_i, indic_ik, 0., temp_mult); // t_XL
		  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., temp_mult, W_icatk, 0., temp_mult1); // t_XLW	    
		  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., temp_mult1, X_i, 0., temp_mult2); // t_XLWX
		  gsl_matrix_add(temp_add, temp_mult2);// somme des t_XLWX
		  
		  
		  //somme des tXI(delta-mu)
		  gsl_blas_dgemv(CblasNoTrans, 1., temp_mult, delta, 0., temp_mult3);// t_XL(Y-mu)
		  gsl_vector_add(add_vect, temp_mult3);// somme des t_XL(Y-mu) 
		  
		  gsl_vector_free(temp_mult3);
		  gsl_vector_free(temp_vect);
		  gsl_vector_free(delta);
		  gsl_vector_free(mu_icatk);
		  gsl_matrix_free(W_icatk);
		  
		}

		gsl_matrix_free(temp_mult);
		gsl_matrix_free(temp_mult1);
		gsl_matrix_free(temp_mult2);
		gsl_matrix_free(temp_mat);
		gsl_matrix_free(indic_ik);
		gsl_matrix_free(X_i);
	      }
	      
	    }
	    
	    gsl_matrix_memcpy(trans, temp_add);
	    gsl_linalg_LU_decomp(trans, p, &s);
	    gsl_linalg_LU_invert(trans, p, inv_temp_add);
	  
	    gsl_blas_dgemv(CblasNoTrans, 1., inv_temp_add, add_vect, 0., beta_temp);
	    
	    gsl_vector_add(beta_catk, beta_temp);
	    gsl_matrix_set_row(beta_k, cat, beta_catk);

	    gsl_vector_free(beta_catk);
	    gsl_vector_free(beta_temp);
	    gsl_vector_free(add_vect);
	    gsl_matrix_free(temp_add);
	    gsl_matrix_free(inv_temp_add);
	    gsl_matrix_free(trans);
	    gsl_permutation_free(p);
	  }
	  
	  
	  // mise à jour des lois d'observations
	  
	  for (cat = 0; cat < nb_cat; cat++){
	    for(j = 0; j < nb_covariable; j++){
	      beta_kk[cat][j] = gsl_matrix_get(beta_k,cat,j);
	    }
	  }
	  
	  hsmarkov->sw_process[1]->observation[k]->Set_regmult(beta_kk);
	  
	  
#   ifdef DEBUG
	  cout<<"beta pour l'état "<<k<<" :  "<<endl;
	  for (cat = 0; cat < nb_cat; cat++){
	  cout<<"Category "<< cat <<"  ";
	  for(i = 0; i<nb_covariable; i++){
	    cout<<beta_kk[cat][i]<<"   ";
	  }
	  cout<<endl;
	  }
	  cout<<endl;
#endif
	  
	  for (cat = 0; cat < nb_cat; cat++){
	    delete beta_kk[cat];
	  }
	  delete [] beta_kk;

	  gsl_matrix_free(beta_k);
	}
	
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
		complete_run[i]->state_occupancy_estimation(final_run[i] , complete_run[i] ,
							    occupancy_survivor ,
							    censored_occupancy_survivor , false);
	      }

            }

            complete_run[i]->max_computation();
            complete_run[i]->mean_computation();
            complete_run[i]->variance_computation();

            if (iter <= EXPLORATION_NB_ITER) {
              occupancy_likelihood = complete_run[i]->parametric_estimation(occupancy , 1 , true ,
                                                                            OCCUPANCY_THRESHOLD);
            }
            else {
              occupancy_likelihood = complete_run[i]->type_parametric_estimation(occupancy , 1 , true ,
                                                                                 OCCUPANCY_THRESHOLD);
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
    while ((likelihood != D_INF) && ((iter < STOCHASTIC_EXPLORATION_NB_ITER + 2) ||
				     ((nb_iter == I_DEFAULT) && (iter < MARKOV_SWITCHING_NB_ITER) &&
				      (((likelihood - previous_likelihood) / -likelihood > MARKOV_SWITCHING_LIKELIHOOD_DIFF) ||
				       (min_likelihood == D_INF) || (nb_likelihood_decrease == 1))) ||
				     ((nb_iter != I_DEFAULT) && (iter < nb_iter))));


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

      for (i = 0;i < hsmarkov->nb_output_process;i++) {
	for (j = 0;j < hsmarkov->nb_state;j++) {
	
	  hsmarkov->sw_process[i+1]->observation[j]->parametric_variance_computation();
	  // hsmarkov->sw_process[i+1]->observation[j]->min_value_computation();
	  // hsmarkov->sw_process[i+1]->observation[j]->max_value_computation();
	  hsmarkov->sw_process[i+1]->observation[j]->nb_value_computation();
	}
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

        for (i = 1;i <= hsmarkov->nb_output_process;i++) {
          if ((hsmarkov->sw_process[i]) && (seq->characteristics[i])) {
            delete seq->characteristics[i];
            seq->characteristics[i] = 0;
          }
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
//         hsw_markov_data->min_computation(0);
//         hsw_markov_data->max_computation(0);

        seq->build_marginal_histogram(0);
        seq->build_characteristic(0);

        seq->build_transition_count(*hsmarkov);
        
#       ifdef MESSAGE
        cout << "\n" << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << ": " << seq->likelihood <<endl;
	  //  << " | " << hmarkov->Markov::likelihood_computation(*seq) << endl;
#       endif

	// calcul des lois d'occupation des etats

        for (i = 0;i < hsmarkov->nb_state;i++) {
          if (hsmarkov->state_subtype[i] == SEMI_MARKOVIAN) {
            hsmarkov->nonparametric_process[0]->sojourn_time[i]
	      ->computation((seq->characteristics[0] ? seq->characteristics[0]->sojourn_time[i]->nb_value : 1) ,OCCUPANCY_THRESHOLD);
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
        seq->state_variable_init(INT_VALUE);

        for (i = 1;i <= hsmarkov->nb_output_process;i++) {
          if ((hsmarkov->sw_process[i]) && (seq->characteristics[i - 1])) {
            delete seq->characteristics[i - 1];
            seq->characteristics[i - 1] = 0;
          }
        }
      }

      // calcul de la vraisemblance et des lois caracteristiques du modele

      seq->hidden_likelihood = hsmarkov->likelihood_computation(*this , seq->posterior_probability, I_DEFAULT, true);

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

