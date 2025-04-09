/* -*-c++-*-
 *  ----------------------------------------------------------------------------
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
 * Correspond aux fonctions communes à tous les semi-Markov switching
 * (modèle linéaire, modèle linéaire mixte avec effets aléatoires 
 * individuels, modèle linéaire mixte avec effets aléatoires temporels)
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


/*------------------------------------------------------------
 *
 * Tirage aléatoire entre 0 et 1 selon une loi uniforme
 *
 *-----------------------------------------------------------*/

double random_unif() 
{
  return (float)rand() / RAND_MAX;
}


/*--------------------------------------------------------------
 *
 * Calcul de la vraisemblance de sequences pour une semi-chaine 
 * de Markov switching.
 * log f(y|\xi)
 *
 * arguments : référence sur un objet Switching_sequence, 
 *             indice de la séquence
 *
 *--------------------------------------------------------------*/

double Markov_switching::likelihood_computation(const Switching_sequence &isw_seq , int index) const
{
  register int i , j , k, m, l;
  int nb_value , occupancy, *pstate;
  double **poutput;
  double likelihood = 0. , proba;

  // verification de la compatibilite entre le modele et les donnees

  if (isw_seq.nb_variable == 2) {
    for (i = 0;i <= 1;i++) {
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
    poutput = new double*[1];

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

	poutput[0] = isw_seq.real_sequence[i][1];


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
	    switch (sw_process[1]->observation[*pstate]->ident){
	    case NO_RANDOM:
	      proba = sw_process[1]->observation[*pstate]->density_computation(*++poutput[0], isw_seq.covar[i][j]);
	      break;
	    case HETERO_RANDOM:
	      proba = sw_process[1]->observation[*pstate]->density_computation(*++poutput[0], isw_seq.covar[i][j],
									       isw_seq.effect[i][*pstate]);
	      break;
	    case YEAR_RANDOM:
	      proba = sw_process[1]->observation[*pstate]->density_computation(*++poutput[0], isw_seq.covar[i][j], 0,
									       isw_seq.year_effect[isw_seq.index[i][j]-1]);
	      break;
	    default:
	      break;
	    }
	    
	    if (proba > 0.) {
	      likelihood += log(proba);
	    }
	    else {
	      likelihood = D_INF;
	      break;
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


    delete [] poutput;
  }

  poutput = NULL;
  delete poutput;

  pstate = NULL;
  delete pstate;
  
  return likelihood;
}


/*--------------------------------------------------------------
 *
 * Calcul de la vraisemblance de sequences pour une chaine de Markov.
 * log f(y|\xi)
 *
 * arguments : référence sur un objet Markov_switching_data
 *
 * A REVOIR , probleme au niveau du calcul de la vraisemblance 
 * de l'histogramme des observations
 *
 *--------------------------------------------------------------*/

double Markov_switching::likelihood_computation(const Markov_switching_data &isw_markov_data) const
{
  register int i , j;
  int nb_value;
  double buff , likelihood = 0.;
  Histogram **initial_run , **final_run , **single_run;

  // verification de la compatibilite entre le modele et les donnees

  if (isw_markov_data.nb_variable == 2) { 
    for (i = 0;i <= 1;i++) {
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
	
	buff = nonparametric_process[0]->sojourn_time[i]->survivor_likelihood_computation(*(isw_markov_data.characteristics[0]->final_run[i]));
	
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
    }
  }
  
  return likelihood;
}


/*--------------------------------------------------------------
 *
 * Calcul de la vraisemblance de sequences pour une chaine de Markov
 * cachee par l'algorithme forward.
 * log f(y|\xi)
 *
 * arguments : reference sur un objet Switching_sequence,pointeur sur
 *              les probabilites a posteriori des sequences d'etats
 *              les plus probables, indice de la sequence. 
 *
 *--------------------------------------------------------------*/

double Markov_switching::likelihood_computation(const Switching_sequence &isw_seq ,
						double *posterior_probability,
                                                int index) const
{
  register int i , j , k ,l, m, p;
  int nb_value, length;
  double **poutput;
  double likelihood = 0. , seq_likelihood, obs_product, **observation, *norm,
    *state_norm, *forward1 , **state_in;
  Parametric *occupancy;


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
    
    poutput = new double*[isw_seq.nb_variable];
    
    for (i = 0;i < isw_seq.nb_sequence;i++) {
      if ((index == I_DEFAULT) || (index == i)) {
	for (j = 0;j < isw_seq.nb_variable;j++) {
	  poutput[j] = isw_seq.real_sequence[i][j];
	}
	seq_likelihood = 0.;
	
	
	
	for (j = 0; j < isw_seq.length[i];j++){
	  norm[j] = 0.;
	  
	  for (k = 0;k < nb_state;k++) {
	    
	    // calcul des probabilites d'observation
	    observation[j][k] = 1.; 
	    
	    switch (sw_process[1]->observation[k]->ident){
	    case NO_RANDOM:
	      observation[j][k] *=  sw_process[1]->observation[k]->density_computation(*poutput[0], isw_seq.covar[i][j]);
	      break;
	    case HETERO_RANDOM:
	      observation[j][k] *=  sw_process[1]->observation[k]->density_computation(*poutput[0], isw_seq.covar[i][j],
										       isw_seq.effect[i][k]);
	      break;
	    case YEAR_RANDOM:
	      observation[j][k] *=  sw_process[1]->observation[k]->density_computation(*poutput[0], isw_seq.covar[i][j], 0,
										       isw_seq.year_effect[isw_seq.index[i][j]-1]);
	      break;
	    default:
	      break;
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
  }
 
  return likelihood;
  
}


/*--------------------------------------------------------------
 *
 * Calcul des sequences d'etats optimale par l'algorithme de Viterbi
 *
 * arguments : reference sur un objet Markov_switching_data,
 *             pointeur sur les probabilites a posteriori des sequences d'etats
 *             les plus probables,flag d'affichage,
 *              sortie sous format R, indice de la sequence.
 *
 * attention, pas de log_computation, je fais les logs
 * directement dans l'algo
 *
*--------------------------------------------------------------*/

double Markov_switching::viterbi(const Markov_switching_data &seq , double *posterior_probability,
				 bool fag, bool output_R, int index) const
{
  register int i , j , k , l, m, r;
  int length, *pstate, **input_state, **optimal_state, **optimal_occupancy;
  int tmp; 
  double **poutput;
  double likelihood = 0. , obs_product, buff , forward_max , **observation, 
    *forward1 , **state_in;
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

  poutput = new double*[1];


  for (i = 0;i < seq.nb_sequence;i++) {
    if ((index == I_DEFAULT) || (index == i)) {
      poutput[0] = seq.real_sequence[i][1];
      
      for (j = 0;j < seq.length[i];j++) {
        for (k = 0;k < nb_state;k++) {
	  
          // calcul des probabilites d'observation
	  
          observation[j][k] = 0.;

	  switch (sw_process[1]->observation[k]->ident){
	  case NO_RANDOM:{
	    buff = log(sw_process[1]->observation[k]->density_computation(*poutput[0], seq.covar[i][j])) ;
	    
	  }
	    break;
	  case HETERO_RANDOM:
	    buff = log(sw_process[1]->observation[k]->density_computation(*poutput[0], seq.covar[i][j], 
									  seq.effect[i][k])) ;
	    break;
	  case YEAR_RANDOM:
	    buff = log(sw_process[1]->observation[k]->density_computation(*poutput[0], seq.covar[i][j], 0, 
									  seq.year_effect[seq.index[i][j]-1])) ;
	    break;
	  default:
	    break;
	  }
	  
	  if (buff == D_INF) {
	    observation[j][k] = D_INF;
	    break;
	  }
	  else {
	    observation[j][k] += buff;
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

	poutput[0]++;
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
    

    
      if(fag){
	for ( j = 0; j < seq.length[i]; j++) {
	  cout << seq.int_sequence[i][0][j] << " ";
	}
	cout<<endl;
	//cout << "# posterior probability: " << posterior_probability[i]<<endl;
      }
      
      if(output_R){
	tmp = seq.max_length - seq.length[i]; 
	
	for( j = 0; j < tmp ; j++){
	  cout<<"-1 ";
	}
	for (j = tmp; j < seq.max_length; j++){
	  cout << seq.int_sequence[i][0][j-tmp] << " ";
	}
	cout<<endl;
	//cout<<  "# posterior probability: " << posterior_probability[i]<<endl;
      }
      
    }
  }
  
  
  if(seq.year_effect){
    if((fag) || (output_R)){
      cout<<endl;
      
      double *sum1;
      double *sum2;
      double *tmp;
      sum1= new double [nb_state];
      sum2 = new double [nb_state];
      tmp = new double[nb_state];
    
      for( i = 0; i < nb_state; i++){
	sum1[i] = 0.;
	sum2[i] = 0.;
	tmp[i] = 0.;
      }
    
      for (k = 0; k < nb_state; k++){
	for (i = 0; i < seq.nb_sequence; i++){
	  for ( j = 0; j < seq.length[i]; j++){
	    if(seq.int_sequence[i][0][j] == k){
	      sum1[k] = sum1[k] +  seq.year_effect[seq.index[i][j]-1]*seq.year_effect[seq.index[i][j]-1];
	      sum2[k] = sum2[k] +  seq.year_effect[seq.index[i][j]-1];
	      tmp[k] = tmp[k]+1;
	    }
	  }
	}
      
	sum1[k] = sum1[k]/tmp[k];
	sum2[k] = sum2[k]/tmp[k];
      
	cout << "Calcul dans le Viterbi - Empirical year_random_variance, state " << k <<": "
	     << sw_process[1]->observation[k]->random_variance[0] * (sum1[k] - sum2[k]*sum2[k]) <<endl;
      
      }
      cout<<endl;
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
  int *pstate ;
  double **poutput;
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

  poutput = new double*[1];

  poutput[0] = seq.real_sequence[index][1];


  // recurrence "forward"

  seq_likelihood = 0.;
  for (i = 0;i < seq.length[index];i++) {
    norm[i] = 0.;

    for (j = 0;j < nb_state;j++) {

      // calcul des probabilites d'observation

      observation[i][j] = 1.;
      switch(sw_process[1]->observation[j]->ident){
      case NO_RANDOM:
	observation[i][j] *= sw_process[1]->observation[j]->density_computation(*poutput[0], seq.covar[index][i]);
	break;
      case HETERO_RANDOM:
	observation[i][j] *= sw_process[1]->observation[j]->density_computation(*poutput[0], seq.covar[index][i],
										seq.effect[index][j]);
	break;
      case YEAR_RANDOM:
	observation[i][j] *= sw_process[1]->observation[j]->density_computation(*poutput[0], seq.covar[index][i], 0,
										seq.year_effect[seq.index[index][i]-1]);
	break;
      default:
	break;
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
            }

            else {
              switch (type) {
              case 'o' :
                occupancy_predicted[k] = obs_product * occupancy->mass[k] * initial[j];
                break;
              case 'e' :
                occupancy_predicted[k] = obs_product * forward[j]->mass[k] * initial[j];
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
            }

            else {
              switch (type) {
              case 'o' :
                occupancy_predicted[k] = obs_product * (1. - occupancy->cumul[k - 1]) * initial[j];
                break;
              case 'e' :
                occupancy_predicted[k] = obs_product * (1. - forward[j]->cumul[k - 1]) * initial[j];
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

    poutput[0]++;
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

    poutput[0]--;

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
	if (sw_process[1]) {
	  switch (sw_process[1]->observation[j]->ident){
	  case NO_RANDOM:
	    if (sw_process[1]->observation[j]->density_computation(*poutput[0], seq.covar[index][i]) > 0.) {
	      entropy2 -= backward[i][j] * log(sw_process[1]->observation[j]->density_computation(*poutput[0], seq.covar[index][i]));
	    }
	    break;
	  case HETERO_RANDOM:
	    if (sw_process[1]->observation[j]->density_computation(*poutput[0], seq.covar[index][i], seq.effect[index][j]) > 0.) {
	      entropy2 -= backward[i][j] * log(sw_process[1]->observation[j]->density_computation(*poutput[0], seq.covar[index][i],
												  seq.effect[index][j]));
	    }
	    break;
	  case YEAR_RANDOM:
	    if (sw_process[1]->observation[j]->density_computation(*poutput[0], seq.covar[index][i], 
								   0, seq.year_effect[seq.index[index][i]-1]) > 0.) {
	      entropy2 -= backward[i][j] * log(sw_process[1]->observation[j]->density_computation(*poutput[0],seq.covar[index][i],0,
												  seq.year_effect[seq.index[index][i]-1]));
	    }
	    break;
	  default:
	    break;
	  }
	}
      }
    }

    for (i = seq.length[index] - 2;i >= 0;i--) {
      poutput[0]--;
    

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
              if (k < seq.length[index] - i - 1) {
                buff = backward1[i + k][j] * obs_product * occupancy->mass[k] /
		  forward1[i + k][j];
                occupancy_auxiliary[k] = buff * state_in[i][j];
                occupancy_entropy[j][k] += occupancy_auxiliary[k];
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
            auxiliary[j] = backward1[i + 1][j] / state_in[i][j];
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
	  if (sw_process[1]) {
	    switch (sw_process[1]->observation[j]->ident){
	    case NO_RANDOM:
	      if (sw_process[1]->observation[j]->density_computation(*poutput[0], seq.covar[index][i]) > 0.) {
		entropy2 -= backward[i][j] * log(sw_process[1]->observation[j]->density_computation(*poutput[0],seq.covar[index][i]));
	      }
	      break;
	    case HETERO_RANDOM:
	      if (sw_process[1]->observation[j]->density_computation(*poutput[0], seq.covar[index][i], seq.effect[index][j]) > 0.) {
		entropy2 -= backward[i][j] * log(sw_process[1]->observation[j]->density_computation(*poutput[0], seq.covar[index][i],
												    seq.effect[index][j]));
	      }
	      break;
	    case YEAR_RANDOM:
	      if (sw_process[1]->observation[j]->density_computation(*poutput[0], seq.covar[index][i], 
								     0, seq.year_effect[seq.index[index][i]-1]) > 0.) {
		entropy2 -= backward[i][j]*log(sw_process[1]->observation[j]->density_computation(*poutput[0],seq.covar[index][i],0,
												  seq.year_effect[seq.index[index][i]-1]));
	      }
	      break;
	    default:
	      break;
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
              if (j < seq.length[index]) {
                switch (type) {

                case 'o' : {
                  occupancy_auxiliary[j] = backward1[j - 1][i] * obs_product * occupancy->mass[j] *
		    initial[i] / forward1[j - 1][i];
                  occupancy_entropy[i][j] += occupancy_auxiliary[j];
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

      os << "\n" << STAT_label[STATL_LIKELIHOOD] << "\t" << seq_likelihood
         << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_LIKELIHOOD] << "\t" << state_seq_likelihood
         << "\t" << exp(state_seq_likelihood - seq_likelihood) << endl;
      break;
    }

    case 'g' : {
      seq.profile_plot_print(os , index , nb_state , backward_output , conditional_entropy ,
                             marginal_entropy , partial_entropy);
      break;
    }
    }

    if (format != 'g') {

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

      poutput[0] = seq.real_sequence[index][1];

      // recurrence "forward"

      for (i = 0;i < seq.length[index];i++) {
        for (j = 0;j < nb_state;j++) {

          // calcul des probabilites d'observation

          observation[i][j] = 1.;
	  if (sw_process[1]) {
	    switch(sw_process[1]->observation[j]->ident){
	    case NO_RANDOM:
	      observation[i][j] *= sw_process[1]->observation[j]->density_computation(*poutput[0], seq.covar[index][i]);
	      break;
	    case HETERO_RANDOM:
	      observation[i][j] *= sw_process[1]->observation[j]->density_computation(*poutput[0], seq.covar[index][i],
										      seq.effect[index][j]);
	      break;
	    case YEAR_RANDOM:
	      observation[i][j] *= sw_process[1]->observation[j]->density_computation(*poutput[0], seq.covar[index][i], 0,
										      seq.year_effect[seq.index[index][i]-1]);
	      break;
	    default:
	      break;
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

	poutput[0]++;
      }

      nb_state_sequence = 0.;
      i = seq.length[index] - 1;
      for (j = 0;j < nb_state;j++) {
        nb_state_sequence += forward1[i][j];
      }

      switch (format) {
      case 'a' :
        os << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_ENTROPY] << ": " << entropy1
           << " (" << entropy1 / seq.length[index] << ")   " << SEQ_label[SEQL_UPPER_BOUND] << ": "
           << log((double)nb_state_sequence) << " ("
           << log((double)nb_state_sequence) / seq.length[index]
           << ")\n" << SEQ_label[SEQL_MARGINAL_ENTROPY_SUM] << ": " << entropy3 << " ("
           << entropy3 / seq.length[index] << ")\n\n"
           << SEQ_label[SEQL_NB_STATE_SEQUENCE] << ": " << nb_state_sequence << endl;
        break;
      case 's' :
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
  int state_occupancy, *pstate ;
  double **poutput;
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

  poutput = new double*[1];

  poutput[0] = seq.real_sequence[index][1];
  

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
      if (sw_process[1]) {
	switch (sw_process[1]->observation[j]->ident){
	case NO_RANDOM:
	  observation[i][j] *= sw_process[1]->observation[j]->density_computation(*poutput[0], seq.covar[index][i]);
	  break;
	case HETERO_RANDOM:
	  observation[i][j] *= sw_process[1]->observation[j]->density_computation(*poutput[0], seq.covar[index][i],
										  seq.effect[index][j]);
	  break;
	case YEAR_RANDOM:
	  observation[i][j] *= sw_process[1]->observation[j]->density_computation(*poutput[0], seq.covar[index][i], 0,
										  seq.year_effect[seq.index[index][i]-1]);
	  break;
	default:
	  break;
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

    poutput[0]++;
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
    ***input_state, ***optimal_state , ***optimal_occupancy, ***input_rank, ***optimal_rank;
  double **poutput;
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

  poutput = new double*[1];
  poutput[0] = seq.real_sequence[index][1];

  nb_state_sequence = 1;


# ifdef DEBUG
  double  **state_sequence_probability;
  
  state_sequence_probability = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    state_sequence_probability[i] = new double[nb_state];
    for (j = 0;j < nb_state;j++) {
      state_sequence_probability[i][j] = D_INF;
    }
  }
# endif

  // recurrence "forward"

  for (i = 0;i < seq.length[index];i++) {
    for (j = 0;j < nb_state;j++) {

      // calcul des probabilites d'observation

      observation[i][j] = 0.;
      if (sw_process[1]) {
	switch(sw_process[1]->observation[j]->ident){
	case NO_RANDOM:
	  buff = log(sw_process[1]->observation[j]->density_computation(*poutput[0], seq.covar[index][i]));
	  break;
	case HETERO_RANDOM:
	  buff = log(sw_process[1]->observation[j]->density_computation(*poutput[0], seq.covar[index][i],
									seq.effect[index][j]));
	  break;
	case YEAR_RANDOM:
	  buff = log(sw_process[1]->observation[j]->density_computation(*poutput[0], seq.covar[index][i],
									0,seq.year_effect[seq.index[index][i]-1]));
	  break;
	default:
	  break;
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

    poutput[0]++;
    
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
	
        os << endl;
      }
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
  int *pstate ;
  double **poutput;
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

  poutput = new double*[1];

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

  poutput[0] = seq.real_sequence[index][1];


  // recurrence "forward"

  for (i = 0;i < seq.length[index];i++) {
    for (j = 0;j < nb_state;j++) {

      // calcul des probabilites d'observation

      observation[i][j] = 0.;
      if (sw_process[1]) {
	switch(sw_process[1]->observation[j]->ident){
	case NO_RANDOM:
	  buff = log(sw_process[1]->observation[j]->density_computation(*poutput[0], seq.covar[index][i]));
	  break;
	case HETERO_RANDOM:
	  buff = log(sw_process[1]->observation[j]->density_computation(*poutput[0], seq.covar[index][i], seq.effect[index][j]));
	  break;
	case YEAR_RANDOM:
	  buff = log(sw_process[1]->observation[j]->density_computation(*poutput[0], seq.covar[index][i], 0, 
									seq.year_effect[seq.index[index][i]-1]));
	  break;
	default:
	  break;
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

    poutput[0]++;
  
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
      for (i = 0;i < seq.length[index];i++) {
        for (j = 0;j < nb_state;j++) {
          if (j != *pstate) {
            ambiguity += backward_output[i][j];
          }
        }
        pstate++;
      }
      ambiguity *= exp(seq_likelihood - state_seq_likelihood);

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


/*---------------------------------------------------------------
 *
 * Calcul des logarithmes des paramètres d'un Markov_switching
 *
 *---------------------------------------------------------------*/

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

}










