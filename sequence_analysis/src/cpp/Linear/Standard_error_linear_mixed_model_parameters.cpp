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
 *
 *  Estimation des standards errors des parametres de régression 
 *  d'un semi-Markov switching linear model a partir de 
 *  l'algorithme EM.
 *
 *  arguments : reference sur un objet Markov_switching, nombre et
 *              longueur des séquences, observations, nombre de 
 *              covariables, covariables, probabilités de lissage 
 *              calculées dans la passe arrière, état
 *
 *--------------------------------------------------------------*/

double* Switching_sequence::standard_error_regression_parameter(const Markov_switching *hsw_markov, int nb_sequence, 
								int *length, double ***real_sequence, int nb_covariable, 
								double ***covar, double ***backward, int k) const
{
  int i, j, m, s;
  double *sebeta;
  double sigma_k;

  gsl_permutation *p = gsl_permutation_alloc(nb_covariable);
  gsl_matrix *trans = gsl_matrix_alloc(nb_covariable, nb_covariable);
  gsl_matrix *inv_observ = gsl_matrix_alloc(nb_covariable, nb_covariable);

  sebeta = new double[nb_covariable]; 

  gsl_matrix *part1 = gsl_matrix_calloc(nb_covariable, nb_covariable);

  for (i = 0; i < nb_sequence; i++) {
	    
    gsl_matrix *X_i = gsl_matrix_calloc(length[i], nb_covariable);
    gsl_matrix *L_ik = gsl_matrix_calloc(length[i], length[i]);
    gsl_matrix *temp_mat = gsl_matrix_calloc(length[i], hsw_markov->nb_state);
    gsl_matrix *temp_mul = gsl_matrix_calloc(nb_covariable, length[i]);
    gsl_matrix *temp_mul2 = gsl_matrix_calloc(nb_covariable, nb_covariable);
    gsl_vector *temp_vect = gsl_vector_calloc(length[i]);

    convert_matrix_double(covar[i], length[i], nb_covariable, X_i);

    convert_matrix_double(backward[i], length[i], hsw_markov->nb_state, temp_mat);
    gsl_matrix_get_col(temp_vect, temp_mat, k);
    convert_matrix_diag_double(temp_vect, length[i], L_ik);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., X_i, L_ik, 0., temp_mul); // t_XI
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., temp_mul, X_i, 0., temp_mul2); // t_XIX
    gsl_matrix_add(part1, temp_mul2);// somme des t_XIX
	

    gsl_vector_free(temp_vect);
    gsl_matrix_free(temp_mul2);
    gsl_matrix_free(temp_mul);
    gsl_matrix_free(temp_mat);
    gsl_matrix_free(L_ik);
    gsl_matrix_free(X_i);
  }

  sigma_k = 1./ hsw_markov->sw_process[1]->observation[k]->residual_variance;

  gsl_matrix_scale(part1, sigma_k);// (somme des t_XIX)/ sigma_k	  

  gsl_matrix_memcpy(trans, part1);
  gsl_linalg_LU_decomp(trans, p, &s);
  gsl_linalg_LU_invert(trans, p, inv_observ);

  for (i = 0; i < nb_covariable; i++){
    sebeta[i] = sqrt(gsl_matrix_get(inv_observ,i,i));
  }


  gsl_matrix_free(part1);
  gsl_matrix_free(inv_observ);
  gsl_matrix_free(trans);
  gsl_permutation_free(p);

  return sebeta;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des standards errors des variances résiduelles
 *  d'un semi-Markov switching linear model a partir de 
 *  l'algorithme EM.
 *
 *  arguments : reference sur un objet Markov_switching, nombre et
 *              longueur des séquences, observations, nombre de 
 *              covariables, covariables, probabilités de lissage 
 *              calculées dans la passe arrière, état
 *
 *--------------------------------------------------------------*/

double Switching_sequence::standard_error_residual_variance(const Markov_switching *hsw_markov, int nb_sequence, 
									   int *length, double ***real_sequence, int nb_covariable, 
									   double ***covar, double ***backward, int k) const
{
  int i, j, m, s;
  double sesigma;
  double sigma_k;
  double inv_observ;
  double diff;

  double part1 = 0.;
  gsl_vector *beta_k = gsl_vector_alloc(nb_covariable);
	
  for (j = 0; j < nb_covariable; j++){
    gsl_vector_set(beta_k,j, hsw_markov->sw_process[1]->observation[k]->regression[j]); 
  }
	
  sigma_k = hsw_markov->sw_process[1]->observation[k]->residual_variance;

  for (i = 0; i < nb_sequence; i++) {
    gsl_matrix *X_i = gsl_matrix_calloc(length[i], nb_covariable);
    gsl_matrix *L_ik = gsl_matrix_calloc(length[i], length[i]);
    gsl_matrix *temp_mat = gsl_matrix_calloc(length[i], hsw_markov->nb_state);
    gsl_vector *Y_i = gsl_vector_calloc(length[i]);
    double inter = 0.;
    gsl_vector *temp_vect = gsl_vector_calloc(length[i]);
	    
    convert_matrix_double(covar[i], length[i], nb_covariable, X_i);
 
    convert_vector_double(real_sequence[i][0], length[i], Y_i);
    diff = 0.;
    convert_matrix_double(backward[i], length[i], hsw_markov->nb_state, temp_mat);
    gsl_matrix_get_col(temp_vect, temp_mat, k);
    convert_matrix_diag_double(temp_vect, length[i], L_ik);
    gsl_vector_set_zero(temp_vect);
    gsl_blas_dgemv(CblasNoTrans, 1., X_i, beta_k, 0., temp_vect); // X\beta
    gsl_vector_sub(Y_i, temp_vect);// Y-X\beta

    for(j = 0; j < length[i]; j++) {
      inter = (gsl_vector_get(Y_i,j) * gsl_vector_get(Y_i,j)) / (sigma_k*sigma_k*sigma_k) - 1./ (2*sigma_k*sigma_k);
      part1 = part1 + backward[i][j][k]*inter;
    }

    gsl_vector_free(temp_vect);
    gsl_vector_free(Y_i);
    gsl_matrix_free(temp_mat);
    gsl_matrix_free(L_ik);
    gsl_matrix_free(X_i);
  }

  inv_observ = 1./part1;

  sesigma = sqrt(inv_observ);

  gsl_vector_free(beta_k);

  return sesigma;
}

/*--------------------------------------------------------------*
 *
 *  Estimation des standards errors des variances résiduelles
 *  d'un semi-Markov switching linear model a partir de 
 *  l'algorithme EM.
 *
 *  arguments : reference sur un objet Markov_switching, nombre et
 *              longueur des séquences, observations, nombre de 
 *              covariables, covariables, probabilités de lissage 
 *              calculées dans la passe arrière, état
 *
 *--------------------------------------------------------------*/

double Switching_sequence::standard_error_residual_variance_VarCommune(const Markov_switching *hsw_markov, int nb_sequence, 
								       int *length, double ***real_sequence, int nb_covariable, 
								       double ***covar, double ***backward) const
{
  int i, j, k, m, s;
  double sesigma;
  double sigma_k;
  double inv_observ;
  double diff;

  double part1 = 0.;
  gsl_vector *beta_k = gsl_vector_alloc(nb_covariable);

  for (k = 0; k < hsw_markov->nb_state;k++){
    for (j = 0; j < nb_covariable; j++){
      gsl_vector_set(beta_k,j, hsw_markov->sw_process[1]->observation[k]->regression[j]); 
    }
    
    sigma_k = hsw_markov->sw_process[1]->observation[k]->residual_variance;
    
    for (i = 0; i < nb_sequence; i++) {
      gsl_matrix *X_i = gsl_matrix_calloc(length[i], nb_covariable);
      gsl_matrix *L_ik = gsl_matrix_calloc(length[i], length[i]);
      gsl_matrix *temp_mat = gsl_matrix_calloc(length[i], hsw_markov->nb_state);
      gsl_vector *Y_i = gsl_vector_calloc(length[i]);
      double inter = 0.;
      gsl_vector *temp_vect = gsl_vector_calloc(length[i]);
      
      convert_matrix_double(covar[i], length[i], nb_covariable, X_i);
      
      convert_vector_double(real_sequence[i][0], length[i], Y_i);
      diff = 0.;
      convert_matrix_double(backward[i], length[i], hsw_markov->nb_state, temp_mat);
      gsl_matrix_get_col(temp_vect, temp_mat, k);
      convert_matrix_diag_double(temp_vect, length[i], L_ik);
      gsl_vector_set_zero(temp_vect);
      gsl_blas_dgemv(CblasNoTrans, 1., X_i, beta_k, 0., temp_vect); // X\beta
      gsl_vector_sub(Y_i, temp_vect);// Y-X\beta
      
      for(j = 0; j < length[i]; j++) {
	inter = (gsl_vector_get(Y_i,j) * gsl_vector_get(Y_i,j)) / (sigma_k*sigma_k*sigma_k) - 1./ (2*sigma_k*sigma_k);
	part1 = part1 + backward[i][j][k]*inter;
      }

      gsl_vector_free(temp_vect);
      gsl_vector_free(Y_i);
      gsl_matrix_free(temp_mat);
      gsl_matrix_free(L_ik);
      gsl_matrix_free(X_i);
    }
    
  }
  
  inv_observ = 1./part1;

  sesigma = sqrt(inv_observ);

  gsl_vector_free(beta_k);

  return sesigma;
}

/*--------------------------------------------------------------*
 *
 *  Estimation des standards errors des parametres de régression 
 *  d'un semi-Markov switching linear mixed model a partir de 
 *  l'algorithme MCEM avec une étape restauration par simulation-prediction.
 *
 *  arguments : reference sur un objet Markov_switching, nombre et
 *              longueur des séquences, observations, nombre de 
 *              covariables, covariables, nombre d'effets aléatoires, 
 *              effets aléatoires prédits, indicatrice des états 
 *              restaurés, état, nombre de séquences d'état simulées.
 *
 *--------------------------------------------------------------*/

double* Switching_sequence::standard_error_regression_parameter_heterogeneity(const Markov_switching *hsw_markov, int nb_sequence, 
									      int *length, double ***real_sequence, int nb_covariable, 
									      double ***covar, int nb_random, double ***random_predict, 
									      double ****indic, int k,  int nb_state_sequence) const
{
  int i, j, m, s;
  double *sebeta;
  double sigma_k;
  double diff;

  gsl_permutation *p = gsl_permutation_alloc(nb_covariable);
  gsl_matrix *trans = gsl_matrix_alloc(nb_covariable, nb_covariable);
  gsl_matrix *inv_observ = gsl_matrix_alloc(nb_covariable, nb_covariable);

  sebeta = new double[nb_covariable]; 

  // calcul de la partie 1

  gsl_matrix *part1 = gsl_matrix_calloc(nb_covariable, nb_covariable);

  for (i = 0; i < nb_sequence; i++) {
	    
    gsl_matrix *X_i = gsl_matrix_calloc(length[i], nb_covariable);
    gsl_matrix *indic_ik = gsl_matrix_calloc(length[i], length[i]);
    gsl_matrix *temp_mat = gsl_matrix_calloc(length[i], hsw_markov->nb_state);
    gsl_matrix *temp_mul = gsl_matrix_calloc(nb_covariable, length[i]);
    gsl_matrix *temp_mul2 = gsl_matrix_calloc(nb_covariable, nb_covariable);

    convert_matrix_double(covar[i], length[i], nb_covariable, X_i);

    for (m = 0; m < nb_state_sequence; m++){

      gsl_vector *temp_vect = gsl_vector_calloc(length[i]);

      convert_matrix_double(indic[i][m], length[i], hsw_markov->nb_state, temp_mat);
      gsl_matrix_get_col(temp_vect, temp_mat, k);
      convert_matrix_diag_double(temp_vect, length[i], indic_ik);
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., X_i, indic_ik, 0., temp_mul); // t_XI
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., temp_mul, X_i, 0., temp_mul2); // t_XIX
      gsl_matrix_add(part1, temp_mul2);// somme des t_XIX
	
      gsl_vector_free(temp_vect);
    }

    gsl_matrix_free(temp_mul2);
    gsl_matrix_free(temp_mul);
    gsl_matrix_free(temp_mat);
    gsl_matrix_free(indic_ik);
    gsl_matrix_free(X_i);
  }

  sigma_k = 1./ hsw_markov->sw_process[1]->observation[k]->residual_variance;

  gsl_matrix_scale(part1, sigma_k);// (somme des t_XIX)/ sigma_k	  
  gsl_matrix_scale(part1, 1./nb_state_sequence);  // partie 1
	
  // calcul de la partie 2

  gsl_vector *part22 = gsl_vector_calloc(nb_covariable);
  gsl_vector *tau_k = gsl_vector_alloc(nb_random);
  gsl_vector *beta_k = gsl_vector_alloc(nb_covariable);
	
  for (j = 0; j < nb_covariable; j++){
    gsl_vector_set(beta_k,j, hsw_markov->sw_process[1]->observation[k]->regression[j]); 
  }
	
  for (j = 0; j < nb_random; j++){
    gsl_vector_set(tau_k,j, hsw_markov->sw_process[1]->observation[k]->random_variance[j]); 
  }
	    
  for (i = 0; i < nb_sequence; i++) {
	    
    gsl_matrix *X_i = gsl_matrix_calloc(length[i], nb_covariable);
    gsl_matrix *indic_ik = gsl_matrix_calloc(length[i], length[i]);
    gsl_matrix *temp_mat = gsl_matrix_calloc(length[i], hsw_markov->nb_state);
    gsl_vector *temp_mul = gsl_vector_calloc(nb_covariable);
    gsl_vector *Y_i = gsl_vector_calloc(length[i]);
	    
    convert_matrix_double(covar[i], length[i], nb_covariable, X_i);
	 
    for (m = 0; m < nb_state_sequence; m++) {
      gsl_vector *temp_vect = gsl_vector_calloc(length[i]);

      convert_vector_double(real_sequence[i][0], length[i], Y_i);
      diff = 0.;
      convert_matrix_double(indic[i][m], length[i], hsw_markov->nb_state, temp_mat);
      gsl_matrix_get_col(temp_vect, temp_mat, k);
      convert_matrix_diag_double(temp_vect, length[i], indic_ik);
      gsl_vector_set_zero(temp_vect);
      gsl_blas_dgemv(CblasNoTrans, 1., X_i, beta_k, 0., temp_vect); // X\beta
      gsl_vector_sub(Y_i, temp_vect);// Y-X\beta
      diff = sqrt(gsl_vector_get(tau_k,0))*(-random_predict[i][m][k]);
      gsl_vector_add_constant(Y_i,diff); // Y-X\beta-tau\xi
      gsl_vector_set_zero(temp_vect);
      gsl_blas_dgemv(CblasNoTrans, 1., indic_ik, Y_i, 0., temp_vect);// I(Y-X\beta-tau\xi)
      gsl_blas_dgemv(CblasTrans, 1., X_i, temp_vect, 0., temp_mul); // t_XI(Y-X\beta-tau\xi)
      gsl_vector_add(part22, temp_mul);// somme des t_XI(Y-X\beta-tau\xi)
      
      gsl_vector_free(temp_vect);
    }
	  
    gsl_vector_free(Y_i);
    gsl_vector_free(temp_mul);
    gsl_matrix_free(temp_mat);
    gsl_matrix_free(indic_ik);
    gsl_matrix_free(X_i);
  }

  gsl_vector_scale(part22, sigma_k); // ( somme des t_XI(Y-X\beta-tau\xi) ) / sigma_k
  gsl_vector_scale(part22, 1./nb_state_sequence);

  gsl_matrix *part21 = gsl_matrix_calloc(nb_covariable, nb_covariable);
  gsl_matrix *part2 = gsl_matrix_calloc(nb_covariable, nb_covariable);
    
  for (m = 0; m < nb_state_sequence; m++) {
    gsl_vector *part21_inter = gsl_vector_calloc(nb_covariable);

    for (i = 0; i < nb_sequence; i++) {
      gsl_matrix *X_i = gsl_matrix_calloc(length[i], nb_covariable);
      gsl_matrix *indic_ik = gsl_matrix_calloc(length[i], length[i]);
      gsl_matrix *temp_mat = gsl_matrix_calloc(length[i], hsw_markov->nb_state);
      gsl_vector *temp_mul = gsl_vector_calloc(nb_covariable);
      gsl_vector *Y_i = gsl_vector_calloc(length[i]);
      gsl_vector *temp_vect = gsl_vector_calloc(length[i]);	    

      convert_matrix_double(covar[i], length[i], nb_covariable, X_i);
      convert_vector_double(real_sequence[i][0], length[i], Y_i);
      diff = 0.;
      convert_matrix_double(indic[i][m], length[i], hsw_markov->nb_state, temp_mat);
      gsl_matrix_get_col(temp_vect, temp_mat, k);
      convert_matrix_diag_double(temp_vect, length[i], indic_ik);
      gsl_vector_set_zero(temp_vect);
      gsl_blas_dgemv(CblasNoTrans, 1., X_i, beta_k, 0., temp_vect); // X\beta
      gsl_vector_sub(Y_i, temp_vect);// Y-X\beta
      diff = sqrt(gsl_vector_get(tau_k,0))*(-random_predict[i][m][k]);
      gsl_vector_add_constant(Y_i,diff); // Y-X\beta-tau\xi
      gsl_vector_set_zero(temp_vect);
      gsl_blas_dgemv(CblasNoTrans, 1., indic_ik, Y_i, 0., temp_vect);// I(Y-X\beta-tau\xi)
      gsl_blas_dgemv(CblasTrans,  1., X_i, temp_vect, 0., temp_mul); // t_XI(Y-X\beta-tau\xi)
      gsl_vector_add(part21_inter, temp_mul);// somme des t_XI(Y-X\beta-tau\xi)
      
      gsl_vector_free(temp_vect);
      gsl_vector_free(Y_i);
      gsl_vector_free(temp_mul);
      gsl_matrix_free(temp_mat);
      gsl_matrix_free(indic_ik);
      gsl_matrix_free(X_i);
    }
	  
    gsl_vector_scale(part21_inter, sigma_k); // ( somme des t_XI(Y-X\beta-tau\xi) ) /sigma_k
    gsl_vector_sub(part21_inter, part22);

    for (i = 0; i < nb_covariable; i++){
      for (j = 0; j < nb_covariable; j++){
	gsl_matrix_set(part21, i, j, gsl_vector_get(part21_inter, i) * gsl_vector_get(part21_inter,j) );
      }
    }

    gsl_matrix_add(part2, part21);

    gsl_vector_free(part21_inter);
  }

  gsl_matrix_scale(part2, 1./nb_state_sequence);
  gsl_matrix_add(part1,part2);
  gsl_matrix_memcpy(trans, part1);
  gsl_linalg_LU_decomp(trans, p, &s);
  gsl_linalg_LU_invert(trans, p, inv_observ);

  for (i = 0; i < nb_covariable; i++){
    sebeta[i] = sqrt(gsl_matrix_get(inv_observ,i,i));
  }

  gsl_matrix_free(part2);
  gsl_matrix_free(part21);
  gsl_vector_free(beta_k);
  gsl_vector_free(tau_k);
  gsl_vector_free(part22);
  gsl_matrix_free(part1);
  gsl_matrix_free(inv_observ);
  gsl_matrix_free(trans);
  gsl_permutation_free(p);

  return sebeta;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des standards errors des variances résiduelles
 *  d'un semi-Markov switching linear mixed model a partir de 
 *  l'algorithme MCEM avec une étape restauration par simulation-prediction.
 *
 *  arguments : reference sur un objet Markov_switching, nombre et
 *              longueur des séquences, observations, nombre de 
 *              covariables, covariables, nombre d'effets aléatoires, 
 *              effets aléatoires prédits, indicatrice des états 
 *              restaurés, état, nombre de séquences d'état simulées.
 *
 *--------------------------------------------------------------*/

double Switching_sequence::standard_error_residual_variance_heterogeneity(const Markov_switching *hsw_markov, int nb_sequence, 
									   int *length, double ***real_sequence, int nb_covariable, 
									   double ***covar, int nb_random, double ***random_predict, 
									   double ****indic, int k,  int nb_state_sequence) const
{
  int i, j, m, s;
  double sesigma;
  double sigma_k;
  double inv_observ;
  double diff;

  // calcul de la partie 1
	
  double part1 = 0.;
  gsl_vector *tau_k = gsl_vector_alloc(nb_random);
  gsl_vector *beta_k = gsl_vector_alloc(nb_covariable);
	
  for (j = 0; j < nb_covariable; j++){
    gsl_vector_set(beta_k,j, hsw_markov->sw_process[1]->observation[k]->regression[j]); 
  }
	
  for (j = 0; j < nb_random; j++){
    gsl_vector_set(tau_k,j, hsw_markov->sw_process[1]->observation[k]->random_variance[j]); 
  }

  sigma_k = hsw_markov->sw_process[1]->observation[k]->residual_variance;

  for (i = 0; i < nb_sequence; i++) {
    gsl_matrix *X_i = gsl_matrix_calloc(length[i], nb_covariable);
    gsl_matrix *indic_ik = gsl_matrix_calloc(length[i], length[i]);
    gsl_matrix *temp_mat = gsl_matrix_calloc(length[i], hsw_markov->nb_state);
    gsl_vector *Y_i = gsl_vector_calloc(length[i]);
	    
    convert_matrix_double(covar[i], length[i], nb_covariable, X_i);
 
    for (m = 0; m < nb_state_sequence; m++) {
      double inter = 0.;
      gsl_vector *temp_vect = gsl_vector_calloc(length[i]);

      convert_vector_double(real_sequence[i][0], length[i], Y_i);
      diff = 0.;
      convert_matrix_double(indic[i][m], length[i], hsw_markov->nb_state, temp_mat);
      gsl_matrix_get_col(temp_vect, temp_mat, k);
      convert_matrix_diag_double(temp_vect, length[i], indic_ik);
      gsl_vector_set_zero(temp_vect);
      gsl_blas_dgemv(CblasNoTrans, 1., X_i, beta_k, 0., temp_vect); // X\beta
      gsl_vector_sub(Y_i, temp_vect);// Y-X\beta
      diff = sqrt(gsl_vector_get(tau_k,0))*(-random_predict[i][m][k]);
      gsl_vector_add_constant(Y_i,diff); // Y-X\beta-tau\xi

      for(j = 0; j < length[i]; j++) {
	inter = (gsl_vector_get(Y_i,j) * gsl_vector_get(Y_i,j)) / (sigma_k*sigma_k*sigma_k) - 1./ (2*sigma_k*sigma_k);
	part1 = part1 + indic[i][m][j][k]*inter;
      }

      gsl_vector_free(temp_vect);
    }
	  
    gsl_vector_free(Y_i);
    gsl_matrix_free(temp_mat);
    gsl_matrix_free(indic_ik);
    gsl_matrix_free(X_i);
  }

  part1 = part1/nb_state_sequence; 

  // calcul de la partie 2

  double part22 = 0.;
	
  for (i = 0; i < nb_sequence; i++) {
    gsl_matrix *X_i = gsl_matrix_calloc(length[i], nb_covariable);
    gsl_matrix *indic_ik = gsl_matrix_calloc(length[i], length[i]);
    gsl_matrix *temp_mat = gsl_matrix_calloc(length[i], hsw_markov->nb_state);
    gsl_vector *Y_i = gsl_vector_calloc(length[i]);
	    
    convert_matrix_double(covar[i], length[i], nb_covariable, X_i);
	 
    for (m = 0; m < nb_state_sequence; m++) {
      double inter = 0.;
      gsl_vector *temp_vect = gsl_vector_calloc(length[i]);

      convert_vector_double(real_sequence[i][0], length[i], Y_i);
      diff = 0.;
      convert_matrix_double(indic[i][m], length[i], hsw_markov->nb_state, temp_mat);
      gsl_matrix_get_col(temp_vect, temp_mat, k);
      convert_matrix_diag_double(temp_vect, length[i], indic_ik);
      gsl_vector_set_zero(temp_vect);
      gsl_blas_dgemv(CblasNoTrans, 1., X_i, beta_k, 0., temp_vect); // X\beta
      gsl_vector_sub(Y_i, temp_vect);// Y-X\beta
      diff = sqrt(gsl_vector_get(tau_k,0))*(-random_predict[i][m][k]);
      gsl_vector_add_constant(Y_i,diff); // Y-X\beta-tau\xi
      
      for(j = 0; j < length[i]; j++) {
	inter = (gsl_vector_get(Y_i,j) * gsl_vector_get(Y_i,j)) / (2*sigma_k*sigma_k) - 1./ (2*sigma_k);
	part22 = part22 + indic[i][m][j][k]*inter;
      }
    
      gsl_vector_free(temp_vect);
    }
	  
    gsl_vector_free(Y_i);
    gsl_matrix_free(temp_mat);
    gsl_matrix_free(indic_ik);
    gsl_matrix_free(X_i);
  }

  part22 = part22/ nb_state_sequence;

  double part21 = 0.;
  double part2 = 0.;
    
  for (m = 0; m < nb_state_sequence; m++) {
    double part21_inter = 0.;

    for (i = 0; i < nb_sequence; i++) {
      double inter = 0.;
      gsl_matrix *X_i = gsl_matrix_calloc(length[i], nb_covariable);
      gsl_matrix *indic_ik = gsl_matrix_calloc(length[i], length[i]);
      gsl_matrix *temp_mat = gsl_matrix_calloc(length[i], hsw_markov->nb_state);
      gsl_vector *Y_i = gsl_vector_calloc(length[i]);
      gsl_vector *temp_vect = gsl_vector_calloc(length[i]);    
      
      convert_matrix_double(covar[i], length[i], nb_covariable, X_i);
      convert_vector_double(real_sequence[i][0], length[i], Y_i);
      diff = 0.;
      convert_matrix_double(indic[i][m], length[i], hsw_markov->nb_state, temp_mat);
      gsl_matrix_get_col(temp_vect, temp_mat, k);
      convert_matrix_diag_double(temp_vect, length[i], indic_ik);
      gsl_vector_set_zero(temp_vect);
      gsl_blas_dgemv(CblasNoTrans, 1., X_i, beta_k, 0., temp_vect); // X\beta
      gsl_vector_sub(Y_i, temp_vect);// Y-X\beta
      diff = sqrt(gsl_vector_get(tau_k,0))*(-random_predict[i][m][k]);
      gsl_vector_add_constant(Y_i,diff); // Y-X\beta-tau\xi

      for(j = 0; j < length[i]; j++) {
	inter = (gsl_vector_get(Y_i,j) * gsl_vector_get(Y_i,j)) / (2*sigma_k*sigma_k) - 1./ (2*sigma_k);
	part21_inter = part21_inter + indic[i][m][j][k]*inter;
      }
      
      gsl_vector_free(temp_vect);
      gsl_vector_free(Y_i);
      gsl_matrix_free(temp_mat);
      gsl_matrix_free(indic_ik);
      gsl_matrix_free(X_i);

    }
	  
    part21_inter = part21_inter - part22;
    part21 = part21_inter*part21_inter;
    part2 = part2 + part21;
  }

  part2 = part2/ nb_state_sequence;
  part1 = part1 + part2;
  inv_observ = 1./part1;

  sesigma = sqrt(inv_observ);

  gsl_vector_free(beta_k);
  gsl_vector_free(tau_k);

  return sesigma;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des standards errors des écart-types aléatoires
 *  d'un semi-Markov switching linear mixed model a partir de 
 *  l'algorithme MCEM avec une étape restauration par simulation-prediction.
 *
 *  arguments : reference sur un objet Markov_switching, nombre et
 *              longueur des séquences, observations, nombre de 
 *              covariables, covariables, nombre d'effets aléatoires, 
 *              effets aléatoires prédits, indicatrice des états 
 *              restaurés, état, nombre de séquences d'état simulées.
 *
 *--------------------------------------------------------------*/

double* Switching_sequence::standard_error_random_variance_heterogeneity(const Markov_switching *hsw_markov, int nb_sequence, 
									 int *length, double ***real_sequence, int nb_covariable, 
									 double ***covar, int nb_random, double ***random_predict, 
									 double ****indic, int k,  int nb_state_sequence) const
{
  int i, j, m, s;
  double *setau;
  double sigma_k;
  double inv_observ;
  double diff;

  setau = new double[nb_random];

  // calcul de la partie 1
  
  double part1 = 0.;
  
  for (i = 0; i < nb_sequence; i++) {
    for (m = 0; m < nb_state_sequence; m++){
      for(j = 0; j < length[i]; j++) {
	part1 = part1 + indic[i][m][j][k] * random_predict[i][m][k] * random_predict[i][m][k];// somme des I\xi\xi
      }
    }
  }
  sigma_k = hsw_markov->sw_process[1]->observation[k]->residual_variance;
  part1 = part1/(sigma_k*nb_state_sequence); // moyenne de (somme des I\xi\xi)/sigma_k
	
  // calcul de la partie 2

  double part22 = 0.;
  gsl_vector *tau_k = gsl_vector_alloc(nb_random);
  gsl_vector *beta_k = gsl_vector_alloc(nb_covariable);
	
  for (j = 0; j < nb_covariable; j++){
    gsl_vector_set(beta_k,j, hsw_markov->sw_process[1]->observation[k]->regression[j]); 
  }
	
  for (j = 0; j < nb_random; j++){
    gsl_vector_set(tau_k,j, hsw_markov->sw_process[1]->observation[k]->random_variance[j]); 
  }
	    
  for (i = 0; i < nb_sequence; i++) {
	    
    gsl_matrix *X_i = gsl_matrix_calloc(length[i], nb_covariable);
    gsl_matrix *indic_ik = gsl_matrix_calloc(length[i], length[i]);
    gsl_matrix *temp_mat = gsl_matrix_calloc(length[i], hsw_markov->nb_state);
    gsl_vector *Y_i = gsl_vector_calloc(length[i]);
	    
    convert_matrix_double(covar[i], length[i], nb_covariable, X_i);
	 
    for (m = 0; m < nb_state_sequence; m++) {
      gsl_vector *temp_vect = gsl_vector_calloc(length[i]);

      convert_vector_double(real_sequence[i][0], length[i], Y_i);
      diff = 0.;
      convert_matrix_double(indic[i][m], length[i], hsw_markov->nb_state, temp_mat);
      gsl_matrix_get_col(temp_vect, temp_mat, k);
      convert_matrix_diag_double(temp_vect, length[i], indic_ik);
      gsl_vector_set_zero(temp_vect);
      gsl_blas_dgemv(CblasNoTrans, 1., X_i, beta_k, 0., temp_vect); // X\beta
      gsl_vector_sub(Y_i, temp_vect);// Y-X\beta
      diff = sqrt(gsl_vector_get(tau_k,0))*(-random_predict[i][m][k]);
      gsl_vector_add_constant(Y_i,diff); // Y-X\beta-tau\xi

      for(j = 0; j < length[i]; j++) {
	part22 = part22 + indic[i][m][j][k] * random_predict[i][m][k] * gsl_vector_get(Y_i,j);// somme des I\xi(Y-X\beta-tau\xi
      }
    
      gsl_vector_free(temp_vect);
    }
	 
    gsl_vector_free(Y_i);
    gsl_matrix_free(temp_mat);
    gsl_matrix_free(indic_ik);
    gsl_matrix_free(X_i);
  }

  part22 = part22/ (nb_state_sequence*sigma_k); // moyenne de la ( somme des I\xi(Y-X\beta-tau\xi) ) / sigma_k
  
  double part21 = 0.;
  double part2 = 0.;
    
  for (m = 0; m < nb_state_sequence; m++) {
    double part21_inter = 0.;

    for (i = 0; i < nb_sequence; i++) {
      gsl_matrix *X_i = gsl_matrix_calloc(length[i], nb_covariable);
      gsl_matrix *indic_ik = gsl_matrix_calloc(length[i], length[i]);
      gsl_matrix *temp_mat = gsl_matrix_calloc(length[i], hsw_markov->nb_state);
      gsl_vector *Y_i = gsl_vector_calloc(length[i]);
      gsl_vector *temp_vect = gsl_vector_calloc(length[i]);
	           
      convert_matrix_double(covar[i], length[i], nb_covariable, X_i);
      convert_vector_double(real_sequence[i][0], length[i], Y_i);
      diff = 0.;
      convert_matrix_double(indic[i][m], length[i], hsw_markov->nb_state, temp_mat);
      gsl_matrix_get_col(temp_vect, temp_mat, k);
      convert_matrix_diag_double(temp_vect, length[i], indic_ik);
      gsl_vector_set_zero(temp_vect);
      gsl_blas_dgemv(CblasNoTrans, 1., X_i, beta_k, 0., temp_vect); // X\beta
      gsl_vector_sub(Y_i, temp_vect);// Y-X\beta
      diff = sqrt(gsl_vector_get(tau_k,0))*(-random_predict[i][m][k]);
      gsl_vector_add_constant(Y_i,diff); // Y-X\beta-tau\xi
      
      for(j = 0; j < length[i]; j++) {
	part21_inter = part21_inter + indic[i][m][j][k] * random_predict[i][m][k]
	  * gsl_vector_get(Y_i,j);// somme des I\xi(Y-X\beta-tau\xi)
      }
    
      gsl_vector_free(temp_vect);
      gsl_vector_free(Y_i);
      gsl_matrix_free(temp_mat);
      gsl_matrix_free(indic_ik);
      gsl_matrix_free(X_i);
      
    }
    
    part21_inter = part21_inter/ sigma_k; // ( somme des I\xi(Y-X\beta-tau\xi) ) /sigma_k
    part21_inter = part21_inter - part22;
    part21 = part21_inter*part21_inter;
    part2 = part2 + part21;
  }

  part2 = part2/ nb_state_sequence;
  part1=part1+part2;
  inv_observ = 1./part1;

  setau[0] = sqrt(inv_observ); 
  
  gsl_vector_free(beta_k);
  gsl_vector_free(tau_k);

  return setau;
}

