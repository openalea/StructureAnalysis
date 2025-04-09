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
 *  Heterogeneité inter-individuelle 
 *  Un effet aléatoire différent par état
 *  
 *  Espérance conditionnelle des effets aléatoires sachant les 
 *  données observées et les séquences d'états simulées.
 *
 *  arguments : reference sur un objet Markov_switching, longueur 
 *              des séquences, covariables, effets aléatoires, 
 *              observations, indicatrice des états restaurés, 
 *              nombre de séquence d'états simulées, individu.
 *
 *--------------------------------------------------------------*/

double** Switching_sequence::conditional_expectation_heterogeneity_state_wise(const Markov_switching *hsw_markov, int *length, 
									      int nb_covariable, double ***covar, double **effect,
									      double ***real_sequence, double ****indic, 
									      int nb_state_sequence, int indiv) const
{
  int j, k, l;
  double **random_effect;

  gsl_matrix *D= gsl_matrix_calloc(hsw_markov->nb_state, hsw_markov->nb_state);
  gsl_vector *Y_a = gsl_vector_calloc(length[indiv]);
  
  for (k = 0; k < hsw_markov->nb_state; k++) {
    gsl_matrix_set(D, k,k, sqrt(hsw_markov->sw_process[1]->observation[k]->random_variance[0]));
  }
  
  random_effect = new double*[nb_state_sequence];
  
  for (j = 0; j < nb_state_sequence; j++){
    
    int s;
    gsl_permutation *p = gsl_permutation_calloc(length[indiv]);
    gsl_matrix *trans = gsl_matrix_calloc(length[indiv], length[indiv]);
    gsl_vector *beta_p = gsl_vector_calloc(hsw_markov->nb_state);
    gsl_matrix *X_a = gsl_matrix_calloc(length[indiv], length[indiv]);
    gsl_matrix *V_a = gsl_matrix_calloc(length[indiv], length[indiv]);
    gsl_matrix *Sigma_a = gsl_matrix_calloc(length[indiv], length[indiv]);
    gsl_matrix *invSigma_a = gsl_matrix_calloc(length[indiv], length[indiv]);
    gsl_matrix *U_a = gsl_matrix_calloc(length[indiv],hsw_markov->nb_state);
    gsl_matrix *temp_mul = gsl_matrix_calloc(length[indiv], hsw_markov->nb_state);
    gsl_matrix *temp_mul2 = gsl_matrix_calloc(length[indiv], hsw_markov->nb_state);
    gsl_vector *temp_add = gsl_vector_calloc(length[indiv]);
    gsl_vector *temp_add1 = gsl_vector_calloc(length[indiv]);  
    gsl_vector *temp_mul1 = gsl_vector_calloc(hsw_markov->nb_state);
    gsl_vector *effecte = gsl_vector_calloc(hsw_markov->nb_state);	
	
    random_effect[j] = new double[hsw_markov->nb_state];
    	    
    convert_vector_double(real_sequence[indiv][0],length[indiv],Y_a);
	    
    for (k = 0; k < length[indiv]; k++){
	      
      for (l = 0; l < hsw_markov->nb_state; l++) {
	gsl_matrix_set(U_a, k, l, indic[indiv][j][k][l]);
		
	if (indic[indiv][j][k][l] == 1) {
	  gsl_matrix_set(V_a, k, k, hsw_markov->sw_process[1]->observation[l]->residual_variance);
	}
      }
    }
	    
	    
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., U_a, D, 0., temp_mul); 
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., temp_mul, D, 0., temp_mul2); 
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., temp_mul2, U_a, 0., Sigma_a);
    gsl_matrix_add(Sigma_a, V_a);
	    
    //inversion de Sigma_a
    gsl_matrix_memcpy(trans, Sigma_a);
    gsl_linalg_LU_decomp(trans, p, &s);
    gsl_linalg_LU_invert(trans, p, invSigma_a);
	    
	    
    gsl_vector_set_zero(temp_add);
	    
    for (k = 0; k < nb_covariable; k++) {
	      
      gsl_matrix_set_zero(temp_mul);
	      
      for (l = 0; l < length[indiv]; l++) {
	gsl_matrix_set(X_a, l, l, covar[indiv][l][k]);
      }
      for (l = 0; l < hsw_markov->nb_state; l++){
	gsl_vector_set(beta_p,l, hsw_markov->sw_process[1]->observation[l]->regression[k]); 
      }
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., X_a, U_a, 0., temp_mul);
	      
	      
      gsl_blas_dgemv(CblasNoTrans, 1., temp_mul, beta_p, 0., temp_add1);
      gsl_vector_add(temp_add, temp_add1);
    }
	    
    gsl_vector_sub(Y_a, temp_add);
	    
    gsl_vector_set_zero(temp_add1);
    gsl_blas_dgemv(CblasNoTrans, 1., invSigma_a, Y_a, 0., temp_add1); 
    gsl_blas_dgemv(CblasTrans, 1., U_a, temp_add1, 0., temp_mul1);
	    
	    
    gsl_blas_dgemv(CblasTrans, 1., D, temp_mul1, 0., effecte);	  
	  
  
    for (k = 0; k < hsw_markov->nb_state;k++) {
      random_effect[j][k] = gsl_vector_get(effecte,k);
    }
	    
    gsl_vector_free(effecte);
    gsl_vector_free(temp_mul1);
    gsl_vector_free(temp_add1);
    gsl_vector_free(temp_add);
    gsl_matrix_free(temp_mul2);
    gsl_matrix_free(temp_mul);
    gsl_matrix_free(U_a);
    gsl_matrix_free(invSigma_a);
    gsl_matrix_free(Sigma_a);
    gsl_matrix_free(V_a);
    gsl_matrix_free(X_a);
    gsl_vector_free(beta_p);
    gsl_matrix_free(trans);
    gsl_permutation_free(p);

  }

  gsl_vector_free(Y_a);
  gsl_matrix_free(D);

  return random_effect;
}



/*--------------------------------------------------------------*
 *
 *  Heterogeneité inter-individuelle 
 *  Un unique effet aléatoire pour toute la séquence observée
 *  
 *  Espérance conditionnelle des effets aléatoires sachant les 
 *  données observées et les séquences d'états simulées.
 *
 *  arguments : reference sur un objet Markov_switching, longueur 
 *              des séquences, covariables, effets aléatoires, 
 *              observations, indicatrice des états restaurés, 
 *              nombre de séquence d'états simulées, individu.
 *
 *--------------------------------------------------------------*/

double** Switching_sequence::conditional_expectation_heterogeneity_wise(const Markov_switching *hsw_markov, int *length, 
									int nb_covariable, double ***covar, double **effect,
									double ***real_sequence, double ****indic, 
									int nb_state_sequence, int indiv) const
{
  int j, k, l;
  double **random_effect;
	  
  gsl_matrix *D = gsl_matrix_calloc(hsw_markov->nb_state, hsw_markov->nb_state);	
  gsl_vector *tau= gsl_vector_calloc( hsw_markov->nb_state);
  gsl_vector *Y_a = gsl_vector_calloc(length[indiv]);
	  
  for (k = 0; k < hsw_markov->nb_state; k++) {
    gsl_vector_set(tau,k, sqrt(hsw_markov->sw_process[1]->observation[k]->random_variance[0]));
  }
	  
  for (k = 0; k < hsw_markov->nb_state; k++) {
    for (j = 0; j < hsw_markov->nb_state; j++){
      gsl_matrix_set(D,k,j, gsl_vector_get(tau,k)*gsl_vector_get(tau,j));
    }
  }	

  random_effect = new double*[nb_state_sequence];
	  
  for (j = 0; j < nb_state_sequence; j++){
	    
    gsl_permutation *p = gsl_permutation_calloc(length[indiv]);
    gsl_matrix *trans = gsl_matrix_calloc(length[indiv], length[indiv]);
	    
    gsl_vector *beta_p = gsl_vector_calloc(hsw_markov->nb_state);
    gsl_matrix *X_a = gsl_matrix_calloc(length[indiv], length[indiv]);
    gsl_matrix *V_a = gsl_matrix_calloc(length[indiv], length[indiv]);
    gsl_matrix *Sigma_a = gsl_matrix_calloc(length[indiv], length[indiv]);
    gsl_matrix *invSigma_a = gsl_matrix_calloc(length[indiv], length[indiv]);
    gsl_matrix *U_a = gsl_matrix_calloc(length[indiv],hsw_markov->nb_state);
    gsl_matrix *temp_mul = gsl_matrix_calloc(length[indiv], hsw_markov->nb_state);
    gsl_vector *temp_add = gsl_vector_calloc(length[indiv]);
    gsl_vector *temp_add1 = gsl_vector_calloc(length[indiv]);  
    gsl_vector *temp_mul1 = gsl_vector_calloc(hsw_markov->nb_state);
	    
    double effecte = 0;

    random_effect[j] = new double[hsw_markov->nb_state];
    
    convert_vector_double(real_sequence[indiv][0],length[indiv],Y_a);
    
    for (k = 0; k < length[indiv]; k++){
      
      for (l = 0; l < hsw_markov->nb_state; l++) {
	gsl_matrix_set(U_a, k, l, indic[indiv][j][k][l]);
	
	if (indic[indiv][j][k][l] == 1) {
	  gsl_matrix_set(V_a, k, k, hsw_markov->sw_process[1]->observation[l]->residual_variance);
	}
      }
    }
	    
	    
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., U_a, D, 0., temp_mul); 
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., temp_mul, U_a, 0., Sigma_a);
    gsl_matrix_add(Sigma_a, V_a);
	    
    //inversion de Sigma_a
    int s;
    gsl_matrix_memcpy(trans, Sigma_a);
    gsl_linalg_LU_decomp(trans, p, &s);
    gsl_linalg_LU_invert(trans, p, invSigma_a);
	    
    gsl_vector_set_zero(temp_add);
	    
    for (k = 0; k < nb_covariable; k++) {
	      
      gsl_matrix_set_zero(temp_mul);
	      
      for (l = 0; l < length[indiv]; l++) {
	gsl_matrix_set(X_a, l, l, covar[indiv][l][k]);
      }
      for (l = 0; l < hsw_markov->nb_state; l++){
	gsl_vector_set(beta_p,l, hsw_markov->sw_process[1]->observation[l]->regression[k]); 
      }
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., X_a, U_a, 0., temp_mul);
	      
	      
      gsl_blas_dgemv(CblasNoTrans, 1., temp_mul, beta_p, 0., temp_add1);
      gsl_vector_add(temp_add, temp_add1);
    }
	    
    gsl_vector_sub(Y_a, temp_add);
	    
    gsl_vector_set_zero(temp_add1);
    gsl_blas_dgemv(CblasNoTrans, 1., invSigma_a, Y_a, 0., temp_add1); 
    gsl_blas_dgemv(CblasTrans, 1., U_a, temp_add1, 0., temp_mul1);
	  
	    
    for (k = 0; k < hsw_markov->nb_state; k++){
      effecte = effecte + gsl_vector_get(tau,k) * gsl_vector_get(temp_mul1,k);
    }
	    
    for (k = 0; k < hsw_markov->nb_state;k++) {
      random_effect[j][k] = effecte;
    }

    gsl_vector_free(temp_mul1);
    gsl_vector_free(temp_add1);
    gsl_vector_free(temp_add);
    gsl_matrix_free(temp_mul);
    gsl_matrix_free(U_a);
    gsl_matrix_free(invSigma_a);
    gsl_matrix_free(Sigma_a);
    gsl_matrix_free(V_a);
    gsl_matrix_free(X_a);
    gsl_vector_free(beta_p);
    gsl_matrix_free(trans);
    gsl_permutation_free(p);
  }
	  
  
  gsl_vector_free(Y_a);
  gsl_vector_free(tau);
  gsl_matrix_free(D);
	  
  return random_effect;
}

