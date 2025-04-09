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
 *  Estimation des parametres de régression d'un semi-Markov 
 *  switching linear a partir de l'algorithme EM.
 *
 *  arguments : reference sur un objet Markov_switching, nombre et
 *              longueur des séquences, observations, nombre de 
 *              covariables, probabilités de lissage calculées dans la 
 *              passe arrière, état.
 *
 *--------------------------------------------------------------*/

double* Switching_sequence::regression_parameter_no_effect(const Markov_switching *hsw_markov, int nb_sequence, 
							   int *length, double ***real_sequence, int nb_covariable, 
							   double ***covar, double ***backward, int k) const
{
  int i, s;
  double *beta;
  gsl_permutation *p = gsl_permutation_alloc(nb_covariable);
  gsl_matrix *trans = gsl_matrix_alloc(nb_covariable, nb_covariable);
  gsl_vector *beta_k = gsl_vector_calloc(nb_covariable);
  gsl_matrix *inv_temp_add = gsl_matrix_alloc(nb_covariable, nb_covariable);
  gsl_matrix *temp_add = gsl_matrix_calloc(nb_covariable, nb_covariable);
  gsl_matrix *temp_mul2 = gsl_matrix_alloc(nb_covariable, nb_covariable);
  gsl_vector *add_vect = gsl_vector_calloc(nb_covariable); 
  gsl_vector *temp_mul3 = gsl_vector_alloc(nb_covariable);

  beta = new double[nb_covariable];

  for (i = 0; i < nb_sequence; i++) {
	  
    gsl_matrix *X_i = gsl_matrix_calloc(length[i], nb_covariable);
    gsl_matrix *L_ik = gsl_matrix_calloc(length[i], length[i]);
    gsl_matrix *temp_mat = gsl_matrix_calloc(length[i], hsw_markov->nb_state);
    gsl_matrix *temp_mul = gsl_matrix_calloc(nb_covariable, length[i]);
    gsl_vector *Y_i = gsl_vector_calloc(length[i]);
    gsl_vector *temp_vect = gsl_vector_calloc(length[i]);
    
    convert_matrix_double(covar[i], length[i], nb_covariable, X_i);
    convert_matrix_double(backward[i], length[i], hsw_markov->nb_state, temp_mat);
    convert_vector_double(real_sequence[i][0], length[i], Y_i);
    gsl_matrix_get_col(temp_vect, temp_mat, k);
    convert_matrix_diag_double(temp_vect, length[i], L_ik);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., X_i, L_ik, 0., temp_mul); // t_XL
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., temp_mul, X_i, 0., temp_mul2); // t_XLX
    gsl_matrix_add(temp_add, temp_mul2);// somme des t_XL_X
    gsl_blas_dgemv(CblasNoTrans, 1., temp_mul, Y_i, 0., temp_mul3);//t_XLY
    gsl_vector_add(add_vect, temp_mul3);//somme des t_XLY 
    
    gsl_vector_free(temp_vect);
    gsl_vector_free(Y_i);
    gsl_matrix_free(temp_mul);
    gsl_matrix_free(temp_mat);
    gsl_matrix_free(L_ik);
    gsl_matrix_free(X_i);
	    
  }
	  
  gsl_matrix_memcpy(trans, temp_add);
  gsl_linalg_LU_decomp(trans, p, &s);
  gsl_linalg_LU_invert(trans, p, inv_temp_add);
  gsl_blas_dgemv(CblasNoTrans, 1., inv_temp_add, add_vect, 0., beta_k);
  convert_array1D_double (beta_k, nb_covariable, beta);

  gsl_vector_free(temp_mul3);
  gsl_vector_free(add_vect);
  gsl_matrix_free(temp_mul2);
  gsl_matrix_free(temp_add);
  gsl_matrix_free(inv_temp_add);
  gsl_vector_free(beta_k);
  gsl_matrix_free(trans);
  gsl_permutation_free(p);

  return beta;
}


/*--------------------------------------------------------------*
 *
 *  Estimation de la variance résiduelle d'un semi-Markov 
 *  switching linear a partir de l'algorithme EM.
 *
 *  arguments : reference sur un objet Markov_switching, nombre et
 *              longueur des séquences, observations, nombre de 
 *              covariables, probabilités de lissage calculées dans la 
 *              passe arrière, paramètre de régression associé à
 *              l'état k, état.
 *
 *--------------------------------------------------------------*/

double Switching_sequence::residual_variance_no_effect(const Markov_switching *hsw_markov, int nb_sequence, 
						       int *length, double ***real_sequence, int nb_covariable, 
						       double ***covar, double ***backward, double *regression, int k) const
{
  int i, j;
  double sigma_k = 0.;
  double sum = 0.;

  gsl_matrix *temp_add = gsl_matrix_calloc(nb_covariable, nb_covariable);
  gsl_matrix *temp_mul2 = gsl_matrix_calloc(nb_covariable, nb_covariable);
  gsl_vector *add_vect = gsl_vector_calloc(nb_covariable); 
  gsl_vector *temp_mul3 = gsl_vector_alloc(nb_covariable);
  gsl_vector *beta_k = gsl_vector_calloc(nb_covariable);

  for (i = 0; i < nb_covariable; i++){
    gsl_vector_set(beta_k, i, regression[i]);
  }

  for (i = 0; i < nb_sequence; i++) {
    
    gsl_matrix *X_i = gsl_matrix_calloc(length[i], nb_covariable);
    gsl_matrix *L_ik = gsl_matrix_calloc(length[i], length[i]);
    gsl_matrix *temp_mat = gsl_matrix_calloc(length[i], hsw_markov->nb_state);
    gsl_matrix *temp_mul = gsl_matrix_calloc(nb_covariable, length[i]);
    gsl_vector *Y_i = gsl_vector_calloc(length[i]);
    gsl_vector *temp_vect = gsl_vector_calloc(length[i]);
    
    convert_matrix_double(covar[i], length[i], nb_covariable, X_i);
    convert_matrix_double(backward[i], length[i], hsw_markov->nb_state, temp_mat);
    convert_vector_double(real_sequence[i][0], length[i], Y_i);
    gsl_matrix_get_col(temp_vect, temp_mat, k);
    convert_matrix_diag_double(temp_vect, length[i], L_ik);
    gsl_blas_dgemv(CblasNoTrans, 1., X_i, beta_k, 0., temp_vect); //X\beta
    gsl_vector_sub(Y_i, temp_vect);//Y-X\beta
    gsl_blas_dgemv(CblasNoTrans, 1., L_ik, Y_i, 0., temp_vect);//L(Y-X\beta)
	    
    for ( j = 0; j < length[i]; j++) {
      sigma_k = sigma_k + gsl_vector_get(temp_vect, j) * gsl_vector_get(Y_i, j); //(Y-X\beta)L(Y-X\beta)
      sum += backward[i][j][k];//somme des L_ik(s)   
    } 
    
    gsl_vector_free(temp_vect);
    gsl_vector_free(Y_i);
    gsl_matrix_free(temp_mul);
    gsl_matrix_free(temp_mat);
    gsl_matrix_free(L_ik);
    gsl_matrix_free(X_i);
    
  }
  sigma_k = sigma_k / sum;
  
  gsl_vector_free(beta_k);
  gsl_vector_free(temp_mul3);
  gsl_vector_free(add_vect);
  gsl_matrix_free(temp_mul2);
  gsl_matrix_free(temp_add);

  return sigma_k;
}


/*--------------------------------------------------------------*
 *
 *  Estimation de la variance résiduelle (supposée sommune à tous 
 *  les états) d'un semi-Markov  switching linear a partir 
 *  de l'algorithme EM.
 *
 *  arguments : reference sur un objet Markov_switching, nombre et
 *              longueur des séquences, observations, nombre de 
 *              covariables, probabilités de lissage calculées dans la 
 *              passe arrière, paramètres de régression.
 *
 *--------------------------------------------------------------*/

double Switching_sequence::residual_variance_no_effect_VarCommune(const Markov_switching *hsw_markov, int nb_sequence, 
								  int *length, double ***real_sequence, int nb_covariable, 
								  double ***covar, double ***backward, double **regression) const
{
  int i, j, k;
  double sigma = 0.;
  double den = 0.;
  double num = 0.;
  double sigma_k = 0.;
  double sum = 0.;
  gsl_matrix *temp_add = gsl_matrix_calloc(nb_covariable, nb_covariable);
  gsl_matrix *temp_mul2 = gsl_matrix_alloc(nb_covariable, nb_covariable);
  gsl_vector *add_vect = gsl_vector_calloc(nb_covariable); 
  gsl_vector *temp_mul3 = gsl_vector_alloc(nb_covariable);
  
  for (k = 0; k < hsw_markov->nb_state; k++){

    gsl_vector *beta_k = gsl_vector_calloc(nb_covariable);

    for (i = 0; i < nb_covariable; i++){
      gsl_vector_set(beta_k, i, regression[k][i]);
    }

    for (i = 0; i < nb_sequence; i++) {
      gsl_matrix *X_i = gsl_matrix_calloc(length[i], nb_covariable);
      gsl_matrix *L_ik = gsl_matrix_calloc(length[i], length[i]);
      gsl_matrix *temp_mat = gsl_matrix_calloc(length[i], hsw_markov->nb_state);
      gsl_matrix *temp_mul = gsl_matrix_calloc(nb_covariable, length[i]);
      gsl_vector *Y_i = gsl_vector_calloc(length[i]);
      gsl_vector *temp_vect = gsl_vector_calloc(length[i]);
      
      convert_matrix_double(covar[i], length[i], nb_covariable, X_i);
      convert_matrix_double(backward[i], length[i], hsw_markov->nb_state, temp_mat);
      convert_vector_double(real_sequence[i][0], length[i], Y_i);
      gsl_matrix_get_col(temp_vect, temp_mat, k);
      convert_matrix_diag_double(temp_vect, length[i], L_ik);
      gsl_blas_dgemv(CblasNoTrans, 1., X_i, beta_k, 0., temp_vect); //X\beta
      gsl_vector_sub(Y_i, temp_vect);//Y-X\beta
      gsl_blas_dgemv(CblasNoTrans, 1., L_ik, Y_i, 0., temp_vect);//L(Y-X\beta)
      
      for ( j = 0; j < length[i]; j++) {
	sigma_k = sigma_k + gsl_vector_get(temp_vect, j) * gsl_vector_get(Y_i, j); //(Y-X\beta)L(Y-X\beta)
	sum += backward[i][j][k];//somme des L_ik(s)   
      }
	   
      gsl_vector_free(temp_vect);
      gsl_vector_free(Y_i);
      gsl_matrix_free(temp_mul);
      gsl_matrix_free(temp_mat);
      gsl_matrix_free(L_ik);
      gsl_matrix_free(X_i);
      
    }

    den = den + sigma_k;
    num = num + sum;

  }
    
  sigma = den/num;
       
  gsl_vector_free(temp_mul3);
  gsl_vector_free(add_vect);
  gsl_matrix_free(temp_mul2);
  gsl_matrix_free(temp_add);

  return sigma;

}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres de régression d'un semi-Markov 
 *  switching linear a partir de l'algorithme MCEM.
 *
 *  arguments : reference sur un objet Markov_switching, nombre et
 *              longueur des séquences, observations, nombre de 
 *              covariables, covariables,  indicatrice des états restaurés, 
 *              état, nombre de séquences d'état simulées.
 *
 *--------------------------------------------------------------*/

double* Switching_sequence::regression_parameter_no_effect_stochastic(const Markov_switching *hsmarkov, int nb_sequence, 
								      int *length, double ***real_sequence, int nb_covariable, 
								      double ***covar, double ****indic, int k,
								      int nb_state_sequence) const
{
  int i, m, s;
  double *beta;

  gsl_permutation *p = gsl_permutation_alloc(nb_covariable);
  gsl_matrix *trans = gsl_matrix_alloc(nb_covariable, nb_covariable);
  gsl_vector *beta_k = gsl_vector_calloc(nb_covariable);
  gsl_matrix *inv_temp_add = gsl_matrix_alloc(nb_covariable, nb_covariable);
  gsl_matrix *temp_add = gsl_matrix_calloc(nb_covariable, nb_covariable);
  gsl_matrix *temp_mul2 = gsl_matrix_calloc(nb_covariable, nb_covariable);
  gsl_vector *add_vect = gsl_vector_calloc(nb_covariable); 
  gsl_vector *temp_mul3 = gsl_vector_calloc(nb_covariable);

  beta = new double[nb_covariable];
	  
  for (i = 0; i < nb_sequence; i++) {
    
    gsl_matrix *X_i = gsl_matrix_calloc(length[i], nb_covariable);
    gsl_matrix *indic_ik = gsl_matrix_calloc(length[i], length[i]);
    gsl_matrix *temp_mat = gsl_matrix_calloc(length[i], hsmarkov->nb_state);
    gsl_matrix *temp_mul = gsl_matrix_calloc(nb_covariable, length[i]);
    gsl_vector *Y_i = gsl_vector_calloc(length[i]);
    
    convert_matrix_double(covar[i], length[i], nb_covariable, X_i);
    convert_vector_double(real_sequence[i][0], length[i], Y_i);
	    
    for (m = 0; m < nb_state_sequence; m++){
      
      gsl_vector *temp_vect = gsl_vector_calloc(length[i]);

      convert_matrix_double(indic[i][m], length[i], hsmarkov->nb_state, temp_mat);
      gsl_matrix_get_col(temp_vect, temp_mat, k);
      convert_matrix_diag_double(temp_vect, length[i], indic_ik);
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., X_i, indic_ik, 0., temp_mul); // t_XI
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., temp_mul, X_i, 0., temp_mul2); // t_XIX
      gsl_matrix_add(temp_add, temp_mul2);// somme des t_XI_X
      gsl_blas_dgemv(CblasNoTrans, 1., temp_mul, Y_i, 0., temp_mul3);// t_XIY
      gsl_vector_add(add_vect, temp_mul3);// somme des t_XIY 
      
      gsl_vector_free(temp_vect);
    }

    gsl_vector_free(Y_i);
    gsl_matrix_free(temp_mul);
    gsl_matrix_free(temp_mat);
    gsl_matrix_free(indic_ik);
    gsl_matrix_free(X_i);
	    
  }
	    
  gsl_matrix_memcpy(trans, temp_add);
  gsl_linalg_LU_decomp(trans, p, &s);
  gsl_linalg_LU_invert(trans, p, inv_temp_add);
  gsl_blas_dgemv(CblasNoTrans, 1., inv_temp_add, add_vect, 0., beta_k);
  convert_array1D_double (beta_k, nb_covariable, beta);

  gsl_vector_free(temp_mul3);
  gsl_vector_free(add_vect);
  gsl_matrix_free(temp_mul2);
  gsl_matrix_free(temp_add);
  gsl_matrix_free(inv_temp_add);
  gsl_vector_free(beta_k);
  gsl_matrix_free(trans);
  gsl_permutation_free(p);
  
  return beta;
}


/*--------------------------------------------------------------*
 *
 *  Estimation de la variance résiduelle d'un semi-Markov 
 *  switching linear a partir de l'algorithme MCEM.
 *
 *  arguments : reference sur un objet Markov_switching, nombre et
 *              longueur des séquences, observations, nombre de 
 *              covariables, covariables, indicatrices des états 
 *              restaurés, paramètre de régression associé à
 *              l'état k, état, nombre de séquences d'états simulées.
 *
 *--------------------------------------------------------------*/

double Switching_sequence::residual_variance_no_effect_stochastic(const Markov_switching *hsmarkov, int nb_sequence, 
								  int *length, double ***real_sequence, int nb_covariable, 
								  double ***covar, double ****indic, double *regression, 
								  int k, int nb_state_sequence) const
{
  int i, j, m;
  double sigma_k = 0.;
  double sum = 0.;

  gsl_matrix *temp_add = gsl_matrix_calloc(nb_covariable, nb_covariable);
  gsl_matrix *temp_mul2 = gsl_matrix_calloc(nb_covariable, nb_covariable);
  gsl_vector *add_vect = gsl_vector_calloc(nb_covariable); 
  gsl_vector *temp_mul3 = gsl_vector_calloc(nb_covariable);
  gsl_vector *beta_k = gsl_vector_calloc(nb_covariable);

  for (i = 0; i < nb_covariable; i++){
    gsl_vector_set(beta_k, i, regression[i]);
  }

  for (i = 0; i < nb_sequence; i++) {
	    
    gsl_matrix *X_i = gsl_matrix_calloc(length[i], nb_covariable);
    gsl_matrix *indic_ik = gsl_matrix_calloc(length[i], length[i]);
    gsl_matrix *temp_mat = gsl_matrix_calloc(length[i], hsmarkov->nb_state);
    gsl_matrix *temp_mul = gsl_matrix_calloc(nb_covariable, length[i]);
    
    convert_matrix_double(covar[i], length[i], nb_covariable, X_i);
    
    for (m = 0; m < nb_state_sequence; m++) {
      gsl_vector *Y_i = gsl_vector_calloc(length[i]);
      gsl_vector *temp_vect = gsl_vector_calloc(length[i]);
      
      convert_vector_double(real_sequence[i][0], length[i], Y_i);
      convert_matrix_double(indic[i][m], length[i], hsmarkov->nb_state, temp_mat);
      gsl_matrix_get_col(temp_vect, temp_mat, k);
      convert_matrix_diag_double(temp_vect, length[i], indic_ik);
      gsl_blas_dgemv(CblasNoTrans, 1., X_i, beta_k, 0., temp_vect); // X\beta
      gsl_vector_sub(Y_i, temp_vect);// Y-X\beta
      gsl_blas_dgemv(CblasNoTrans, 1., indic_ik, Y_i, 0., temp_vect);// I(Y-X\beta)
	      
      for ( j = 0; j < length[i]; j++) {
	sigma_k = sigma_k + gsl_vector_get(temp_vect, j) * gsl_vector_get(Y_i, j); // (Y-X\beta)I(Y-X\beta)
	sum += indic[i][m][j][k];	// somme des I_ik(s)   
      } 

      gsl_vector_free(temp_vect);
      gsl_vector_free(Y_i);
    }

    gsl_matrix_free(temp_mul);
    gsl_matrix_free(temp_mat);
    gsl_matrix_free(indic_ik);
    gsl_matrix_free(X_i);
    
  }

  sigma_k = sigma_k / sum;


  gsl_vector_free(beta_k);
  gsl_vector_free(temp_mul3);
  gsl_vector_free(add_vect);
  gsl_matrix_free(temp_mul2);
  gsl_matrix_free(temp_add);
  
  return sigma_k;
}


/*--------------------------------------------------------------*
 *
 *  Estimation de la variance résiduelle (supposée commune à tous 
 *  les états) d'un semi-Markov switching linear a partir 
 *  de l'algorithme MCEM.
 *
 *  arguments : reference sur un objet Markov_switching, nombre et
 *              longueur des séquences, observations, nombre de 
 *              covariables, covariables, indicatrices des états 
 *              restaurés, paramètres de régression, nombre de 
 *              sequence d'états simulées.
 *
 *--------------------------------------------------------------*/

double Switching_sequence::residual_variance_no_effect_VarCommune_stochastic(const Markov_switching *hsmarkov, int nb_sequence, 
									     int *length, double ***real_sequence, int nb_covariable, 
									     double ***covar, double ****indic, double **regression,
									     int nb_state_sequence) const
{
  int i, j, k,m;
  double sigma = 0.;
  double sum = 0.;

  for (k = 0; k < hsmarkov->nb_state; k++){

    gsl_matrix *temp_add = gsl_matrix_calloc(nb_covariable, nb_covariable);
    gsl_matrix *temp_mul2 = gsl_matrix_calloc(nb_covariable, nb_covariable);
    gsl_vector *add_vect = gsl_vector_calloc(nb_covariable); 
    gsl_vector *temp_mul3 = gsl_vector_calloc(nb_covariable);
    gsl_vector *beta_k = gsl_vector_calloc(nb_covariable);
    
    for (i = 0; i < nb_covariable; i++){
      gsl_vector_set(beta_k, i, regression[k][i]);
    }
    
    for (i = 0; i < nb_sequence; i++) {
      
      gsl_matrix *X_i = gsl_matrix_calloc(length[i], nb_covariable);
      gsl_matrix *indic_ik = gsl_matrix_calloc(length[i], length[i]);
      gsl_matrix *temp_mat = gsl_matrix_calloc(length[i], hsmarkov->nb_state);
      gsl_matrix *temp_mul = gsl_matrix_calloc(nb_covariable, length[i]);
      
      convert_matrix_double(covar[i], length[i], nb_covariable, X_i);
      
      for (m = 0; m < nb_state_sequence; m++) {
	gsl_vector *Y_i = gsl_vector_calloc(length[i]);
	gsl_vector *temp_vect = gsl_vector_calloc(length[i]);
	
	convert_vector_double(real_sequence[i][0], length[i], Y_i);
	convert_matrix_double(indic[i][m], length[i], hsmarkov->nb_state, temp_mat);
	gsl_matrix_get_col(temp_vect, temp_mat, k);
	convert_matrix_diag_double(temp_vect, length[i], indic_ik);
	gsl_blas_dgemv(CblasNoTrans, 1., X_i, beta_k, 0., temp_vect); // X\beta
	gsl_vector_sub(Y_i, temp_vect);// Y-X\beta
	gsl_blas_dgemv(CblasNoTrans, 1., indic_ik, Y_i, 0., temp_vect);// I(Y-X\beta)
	
	for ( j = 0; j < length[i]; j++) {
	  sigma = sigma + gsl_vector_get(temp_vect, j) * gsl_vector_get(Y_i, j); // (Y-X\beta)I(Y-X\beta)
	  sum += indic[i][m][j][k];	// somme des I_ik(s)   
	} 
	
	gsl_vector_free(temp_vect);
	gsl_vector_free(Y_i);
      }
      
      gsl_matrix_free(temp_mul);
      gsl_matrix_free(temp_mat);
      gsl_matrix_free(indic_ik);
      gsl_matrix_free(X_i);
      
    }

    gsl_vector_free(beta_k);
    gsl_vector_free(temp_mul3);
    gsl_vector_free(add_vect);
    gsl_matrix_free(temp_mul2);
    gsl_matrix_free(temp_add);
    
  }

  sigma = sigma/sum;

  return sigma;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres de régression d'un semi-Markov 
 *  switching linear mixed model a partir de l'algorithme MCEM 
 *  avec une étape restauration par simulation-prediction.
 *
 *  arguments : reference sur un objet Markov_switching, nombre et
 *              longueur des séquences, observations, nombre de 
 *              covariables, covariables, effets aléatoires prédits, 
 *              indicatrice des états restaurés, état, nombre de 
 *              séquences d'état simulées.
 *
 *--------------------------------------------------------------*/

double* Switching_sequence::regression_parameter_heterogeneity_stochastic(const Markov_switching *hsw_markov, int nb_sequence, 
									  int *length, double ***real_sequence, int nb_covariable, 
									  double ***covar, double ***random_predict, double ****indic,  
									  int k,  int nb_state_sequence) const
{
  int i, m, s;
  double *beta;
  double diff;

  gsl_permutation *p = gsl_permutation_alloc(nb_covariable);
  gsl_matrix *trans = gsl_matrix_alloc(nb_covariable, nb_covariable);
  gsl_matrix *temp_add = gsl_matrix_calloc(nb_covariable, nb_covariable);
  gsl_vector *add_vect = gsl_vector_calloc(nb_covariable);
  gsl_vector *beta_k = gsl_vector_calloc(nb_covariable);
  gsl_matrix *inv_temp_add = gsl_matrix_alloc(nb_covariable, nb_covariable);

  beta = new double[nb_covariable];

  for (i = 0; i < nb_sequence; i++) {
	    
    gsl_matrix *X_i = gsl_matrix_calloc(length[i], nb_covariable);
    gsl_matrix *indic_ik = gsl_matrix_calloc(length[i], length[i]);
    gsl_matrix *temp_mat = gsl_matrix_calloc(length[i], hsw_markov->nb_state);
    gsl_matrix *temp_mul = gsl_matrix_calloc(nb_covariable, length[i]);
    gsl_matrix *temp_mul2 = gsl_matrix_calloc(nb_covariable, nb_covariable);
    gsl_vector *temp_mul3 = gsl_vector_calloc(nb_covariable);
    gsl_vector *Y_i = gsl_vector_calloc(length[i]);

    convert_matrix_double(covar[i], length[i], nb_covariable, X_i);
	    
    for (m = 0; m < nb_state_sequence; m++){

      gsl_vector *temp_vect = gsl_vector_calloc(length[i]);

      convert_vector_double(real_sequence[i][0], length[i], Y_i);
      diff = 0.;
      convert_matrix_double(indic[i][m], length[i], hsw_markov->nb_state, temp_mat);
      gsl_matrix_get_col(temp_vect, temp_mat, k);
      convert_matrix_diag_double(temp_vect, length[i], indic_ik);
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., X_i, indic_ik, 0., temp_mul); // t_XIL
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., temp_mul, X_i, 0., temp_mul2); // t_XIX
      gsl_matrix_add(temp_add, temp_mul2);// somme des t_XIX
      diff = sqrt( hsw_markov->sw_process[1]->observation[k]->random_variance[0])*(-random_predict[i][m][k]);
      gsl_vector_add_constant(Y_i,diff); // Y-\tau\xi
      gsl_blas_dgemv(CblasNoTrans, 1., temp_mul, Y_i, 0., temp_mul3);// t_XI(Y-\tau\xi)
      gsl_vector_add(add_vect, temp_mul3);// somme des t_XI(Y-\tau\xi) 
	   
      gsl_vector_free(temp_vect);
    }
	    
    gsl_vector_free(Y_i);
    gsl_vector_free(temp_mul3);
    gsl_matrix_free(temp_mul2);
    gsl_matrix_free(temp_mul);
    gsl_matrix_free(temp_mat);
    gsl_matrix_free(indic_ik);
    gsl_matrix_free(X_i);
	    
  }
	  
  gsl_matrix_memcpy(trans, temp_add);
  gsl_linalg_LU_decomp(trans, p, &s);
  gsl_linalg_LU_invert(trans, p, inv_temp_add);
  gsl_blas_dgemv(CblasNoTrans, 1., inv_temp_add, add_vect, 0., beta_k);
  convert_array1D_double (beta_k, nb_covariable, beta);

  gsl_matrix_free(inv_temp_add);
  gsl_vector_free(beta_k);
  gsl_vector_free(add_vect);
  gsl_matrix_free(temp_add);
  gsl_matrix_free(trans);
  gsl_permutation_free(p);
   
  return beta;
}


/*--------------------------------------------------------------*
 *
 *  Estimation de la variance résiduelle d'un semi-Markov 
 *  switching linear mixed model a partir de l'algorithme MCEM 
 *  avec une étape restauration par simulation-prediction.
 *
 *  arguments : reference sur un objet Markov_switching, nombre et
 *              longueur des séquences, observations, nombre de 
 *              covariables, covariables, effets aléatoires prédits,
 *              indicatrice des états restaurés, paramètre de 
 *              regression associé à l'état k, état, nombre de 
 *              séquences d'état simulées.
 *
 *--------------------------------------------------------------*/

double Switching_sequence::residual_variance_heterogeneity_stochastic(const Markov_switching *hsw_markov, int nb_sequence, 
								      int *length, double ***real_sequence, int nb_covariable, 
								      double ***covar, double ***random_predict, double ****indic, 
								      double *regression, int k, int nb_state_sequence) const
{
  int i, j, m;
  double sigma_k = 0.;
  double sum = 0.;
  double diff;
  double tau_k;

  gsl_vector *beta_k = gsl_vector_calloc(nb_covariable);

  for (i = 0; i < nb_covariable; i++){
    gsl_vector_set(beta_k, i, regression[i]);
  }
  
  tau_k =  hsw_markov->sw_process[1]->observation[k]->random_variance[0];  
  
  for (i = 0; i < nb_sequence; i++) {
	    
    gsl_matrix *X_i = gsl_matrix_calloc(length[i], nb_covariable);
    gsl_matrix *indic_ik = gsl_matrix_calloc(length[i], length[i]);
    gsl_matrix *temp_mat = gsl_matrix_calloc(length[i], hsw_markov->nb_state);
    gsl_matrix *temp_mul = gsl_matrix_calloc(nb_covariable, length[i]);
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
      diff = sqrt(tau_k)*(-random_predict[i][m][k]);
      gsl_vector_add_constant(Y_i,diff);
      gsl_vector_set_zero(temp_vect);
      gsl_blas_dgemv(CblasNoTrans, 1., indic_ik, Y_i, 0., temp_vect);// I(Y-X\beta-tau\xi)
	      
      for (j = 0; j < length[i]; j++) {
	sigma_k = sigma_k + gsl_vector_get(temp_vect, j) * gsl_vector_get(Y_i, j); // (Y-X\beta-tau\xi)I(Y-X\beta-tau\xi)
	sum += indic[i][m][j][k];	// somme des I_ik(s)   
      } 
	   
      gsl_vector_free(temp_vect);
	      
    }
	    
    gsl_vector_free(Y_i);
    gsl_matrix_free(temp_mul);
    gsl_matrix_free(temp_mat);
    gsl_matrix_free(indic_ik);
    gsl_matrix_free(X_i);
  }

  gsl_vector_free(beta_k);

  sigma_k = sigma_k / sum;	  

  return sigma_k;
}


/*--------------------------------------------------------------*
 *
 *  Estimation de la variance aléatoire (tau*tau) associé aux effets 
 *  aléatoires hétérogénéité d'un semi-Markov switching linear mixed
 *  model a partir de l'algorithme MCEM avec une étape restauration 
 *  par simulation-prediction.
 *
 *  arguments : reference sur un objet Markov_switching, nombre et
 *              longueur des séquences, observations, nombre de 
 *              covariables, covariables, effets aléatoires prédits, 
 *              indicatrice des états restaurés, état, nombre de 
 *              séquences d'état simulées.
 *
 *--------------------------------------------------------------*/

double* Switching_sequence::random_variance_heterogeneity_stochastic(const Markov_switching *hsw_markov, int nb_sequence, 
								     int *length, double ***real_sequence, int nb_covariable, 
								     double ***covar, double ***random_predict, double ****indic,
								     double *regression, int k, int nb_state_sequence) const
{
  int i, m, j, s;
  double *tau;
  double num = 0.;
  double den = 0.;
  gsl_vector *beta_k = gsl_vector_calloc(nb_covariable);
  gsl_vector *tau_k = gsl_vector_calloc(1); // nb_random = 1 (un seul effet aléatoire)

  for (i = 0; i < nb_covariable; i++){
    gsl_vector_set(beta_k, i, regression[i]);
  }

  tau = new double[nb_random];
	  
  for (i = 0; i < nb_sequence; i++) {
	    
    gsl_matrix *X_i = gsl_matrix_calloc(length[i], nb_covariable);
    gsl_matrix *indic_ik = gsl_matrix_calloc(length[i], length[i]);
    gsl_matrix *temp_mat = gsl_matrix_calloc(length[i], hsw_markov->nb_state);
    gsl_matrix *temp_mul = gsl_matrix_calloc(nb_covariable, length[i]);
    gsl_vector *Y_i = gsl_vector_calloc(length[i]);
	    	  
    convert_matrix_double(covar[i], length[i], nb_covariable, X_i);
	    	    
    for (m = 0; m < nb_state_sequence; m++){
      gsl_vector *temp_vect = gsl_vector_calloc(length[i]);
     
      convert_vector_double(real_sequence[i][0], length[i], Y_i);
      convert_matrix_double(indic[i][m], length[i], hsw_markov->nb_state, temp_mat);
      gsl_matrix_get_col(temp_vect, temp_mat, k);
      convert_matrix_diag_double(temp_vect, length[i], indic_ik);
      gsl_vector_set_zero(temp_vect);
      gsl_blas_dgemv(CblasNoTrans, 1., X_i, beta_k, 0., temp_vect); // X\beta
      gsl_vector_sub(Y_i, temp_vect);// Y - X\beta 
      gsl_vector_set_zero(temp_vect);
      gsl_blas_dgemv(CblasNoTrans, 1., indic_ik, Y_i, 0., temp_vect); // I(Y- X\beta)	   
	      
      for(j = 0; j < length[i]; j++) {
	den = den + gsl_vector_get(temp_vect, j) * random_predict[i][m][k];
	num = num + indic[i][m][j][k] * random_predict[i][m][k] * random_predict[i][m][k];
      }

      gsl_vector_free(temp_vect);
    }
	    	    
    gsl_vector_free(Y_i);
    gsl_matrix_free(temp_mul);
    gsl_matrix_free(temp_mat);
    gsl_matrix_free(indic_ik);
    gsl_matrix_free(X_i);
	    
  }

  gsl_vector_set(tau_k,0,den/num);
  gsl_vector_set(tau_k,0,gsl_vector_get(tau_k,0)*gsl_vector_get(tau_k,0));
  convert_array1D_double(tau_k, nb_random, tau);

  gsl_vector_free(tau_k);
  gsl_vector_free(beta_k);
   
  return tau;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres de régression d'un semi-Markov 
 *  switching linear mixed model a partir de l'algorithme MCEM 
 *  avec une étape restauration probabiliste-simulation.
 *
 *  arguments : reference sur un objet Markov_switching, nombre et
 *              longueur des séquences, observations, nombre de 
 *              covariables, covariables, effets aléatoires temporels, 
 *              probabilités de lissage calculées dans la passe 
 *              arrière, paramètres d'index, nombre d'effets 
 *              aléatoires, état
 *
 *--------------------------------------------------------------*/

double* Switching_sequence::regression_parameter_temporal_stochastic(const Markov_switching *hsw_markov, int nb_sequence, 
								     int *length, double ***real_sequence, int nb_covariable, 
								     double ***covar, double *year_effect, double ***backward,
								     int **index, int nb_random, int k) const
{
  int i, m, s;
  double *beta;
  double rand_effect;

  gsl_permutation *p = gsl_permutation_alloc(nb_covariable);
  gsl_matrix *trans = gsl_matrix_alloc(nb_covariable, nb_covariable);
  gsl_vector *beta_k = gsl_vector_alloc(nb_covariable);
  gsl_matrix *inv_temp_add = gsl_matrix_alloc(nb_covariable, nb_covariable);
  gsl_matrix *temp_add = gsl_matrix_calloc(nb_covariable, nb_covariable);
  gsl_vector *add_vect = gsl_vector_calloc(nb_covariable); 
  
  beta = new double[nb_covariable];

  rand_effect = sqrt(hsw_markov->sw_process[1]->observation[k]->random_variance[nb_random-1]);// écart type de l'effet année

  for (i = 0; i < nb_sequence; i++) {
    
    gsl_matrix *X_i = gsl_matrix_calloc(length[i], nb_covariable);
    gsl_matrix *L_ik = gsl_matrix_calloc(length[i], length[i]);
    gsl_matrix *temp_mat = gsl_matrix_calloc(length[i], hsw_markov->nb_state);
    gsl_matrix *temp_mul = gsl_matrix_calloc(nb_covariable, length[i]);
    gsl_matrix *temp_mul2 = gsl_matrix_alloc(nb_covariable, nb_covariable);
    gsl_vector *temp_mul3 = gsl_vector_alloc(nb_covariable);
    gsl_vector *Y_i = gsl_vector_calloc(length[i]);
    gsl_vector *temp_vect = gsl_vector_calloc(length[i]);
    gsl_vector *year_i = gsl_vector_calloc(length[i]);
    
    for (m = 0; m < length[i]; m++){
      gsl_vector_set(year_i, m, year_effect[index[i][m]-1]);
    }
	      
    gsl_vector_scale(year_i, rand_effect);//\tau\xi
    convert_matrix_double(covar[i], length[i], nb_covariable, X_i);
    convert_matrix_double(backward[i], length[i], hsw_markov->nb_state, temp_mat);
    convert_vector_double(real_sequence[i][0], length[i], Y_i);
    
    gsl_matrix_get_col(temp_vect, temp_mat, k);
    convert_matrix_diag_double(temp_vect, length[i], L_ik);	     
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., X_i, L_ik, 0., temp_mul); // t_XL
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., temp_mul, X_i, 0., temp_mul2); // t_XLX
    gsl_matrix_add(temp_add, temp_mul2);// somme des t_XL_X
    gsl_vector_sub(Y_i, year_i);//Y - \tau\xi
    gsl_blas_dgemv(CblasNoTrans, 1., temp_mul, Y_i, 0., temp_mul3);// t_XL(Y - \tau\xi)
    gsl_vector_add(add_vect, temp_mul3);// somme des t_XL(Y - \tau\xi) 
	
    gsl_vector_free(year_i);
    gsl_vector_free(temp_vect);
    gsl_vector_free(Y_i);
    gsl_vector_free(temp_mul3);
    gsl_matrix_free(temp_mul2);
    gsl_matrix_free(temp_mul);
    gsl_matrix_free(temp_mat);
    gsl_matrix_free(L_ik);
    gsl_matrix_free(X_i);
    
  }
	  
  gsl_matrix_memcpy(trans, temp_add);
  gsl_linalg_LU_decomp(trans, p, &s);
  gsl_linalg_LU_invert(trans, p, inv_temp_add);
  gsl_blas_dgemv(CblasNoTrans, 1., inv_temp_add, add_vect, 0., beta_k);
  convert_array1D_double (beta_k, nb_covariable, beta);

  gsl_vector_free(add_vect);
  gsl_matrix_free(temp_add);
  gsl_matrix_free(inv_temp_add);
  gsl_vector_free(beta_k);
  gsl_matrix_free(trans);
  gsl_permutation_free(p);
	
   
  return beta;
}


/*--------------------------------------------------------------*
 *
 *  Estimation de la variance résiduelle d'un semi-Markov switching
 *  linear mixed model a  partir de l'algorithme MCEM.
 *
 *  arguments : reference sur un objet Markov_switching, nombre et
 *              longueur des séquences, observations, nombre de 
 *              covariables, covariables, effets aléatoires temporels, 
 *              probabilités de lissage calculées dans la passe 
 *              arrière, paramètres d'index, paramètres de régression, 
 *              nombre d'effets aléatoires, état
 *
 *--------------------------------------------------------------*/

double Switching_sequence::residual_variance_temporal_stochastic(const Markov_switching *hsw_markov, int nb_sequence, 
								 int *length, double ***real_sequence, int nb_covariable, 
								 double ***covar, double *year_effect, double ***backward,
								 int **index, double *regression, int nb_random, int k) const
{
  int i, j, m;
  double sigma = 0.;
  double sum = 0.;
  double rand_effect;
  gsl_vector *beta_k = gsl_vector_alloc(nb_covariable);
	
  for(j = 0; j < nb_covariable; j++){
    gsl_vector_set(beta_k, j, regression[j]);
  }
	
  rand_effect = sqrt(hsw_markov->sw_process[1]->observation[k]->random_variance[nb_random-1]);

  for (i = 0; i < nb_sequence; i++) {
    
    gsl_matrix *X_i = gsl_matrix_calloc(length[i], nb_covariable);
    gsl_matrix *L_ik = gsl_matrix_calloc(length[i], length[i]);
    gsl_matrix *temp_mat = gsl_matrix_calloc(length[i], hsw_markov->nb_state);
    gsl_matrix *temp_mul = gsl_matrix_calloc(nb_covariable, length[i]);
    gsl_vector *Y_i = gsl_vector_calloc(length[i]);
    gsl_vector *temp_vect = gsl_vector_calloc(length[i]);
    gsl_vector *year_i = gsl_vector_calloc(length[i]);

    for (m = 0; m < length[i]; m++){
      gsl_vector_set(year_i, m, year_effect[index[i][m]-1]);
    }
    
    gsl_vector_scale(year_i, rand_effect);//\tau\xi
    convert_matrix_double(covar[i], length[i], nb_covariable, X_i);
    convert_matrix_double(backward[i], length[i], hsw_markov->nb_state, temp_mat);
    convert_vector_double(real_sequence[i][0], length[i], Y_i);
    gsl_matrix_get_col(temp_vect, temp_mat, k);
    convert_matrix_diag_double(temp_vect, length[i], L_ik);
    gsl_blas_dgemv(CblasNoTrans, 1., X_i, beta_k, 0., temp_vect); // X\beta
    gsl_vector_sub(Y_i, temp_vect);// Y-X\beta
    gsl_vector_sub(Y_i, year_i);//Y- X\beta - \tau\xi
    gsl_vector_set_zero(temp_vect);
    gsl_blas_dgemv(CblasNoTrans, 1., L_ik, Y_i, 0., temp_vect);// L(Y-X\beta - \tau\xi)
	    
    for ( j = 0; j < length[i]; j++) {
      sigma = sigma + gsl_vector_get(temp_vect, j) * gsl_vector_get(Y_i, j); // (Y-X\beta - \tau\xi)L(Y-X\beta-\tau\xi)
      sum += backward[i][j][k];	// somme des L_ik(s)   
    } 
	   
    gsl_vector_free(year_i);
    gsl_vector_free(temp_vect);
    gsl_vector_free(Y_i);
    gsl_matrix_free(temp_mul);
    gsl_matrix_free(temp_mat);
    gsl_matrix_free(L_ik);
    gsl_matrix_free(X_i);
 
  }
	  
  sigma = sigma / sum;

  gsl_vector_free(beta_k);
      
  return sigma;
}


/*--------------------------------------------------------------*
 *
 *  Estimation de la variance résiduelle (supposée commune à tous 
 *  les états) d'un semi-Markov switching linear mixed model a 
 *  partir de l'algorithme MCEM.
 *
 *  arguments : reference sur un objet Markov_switching, nombre et
 *              longueur des séquences, observations, nombre de 
 *              covariables, covariables, effets aléatoires temporels, 
 *              probabilités de lissage calculées dans la passe 
 *              arrière, paramètres d'index, paramètres de régression, 
 *              nombre d'effets aléatoires
 *
 *--------------------------------------------------------------*/

double Switching_sequence::residual_variance_temporal_VarCommune_stochastic(const Markov_switching *hsw_markov, int nb_sequence, 
									    int *length, double ***real_sequence, int nb_covariable, 
									    double ***covar, double *year_effect, double ***backward,
									    int **index, double **regression, int nb_random) const
{
  int i, j, k,m;
  double sigma = 0.;
  double den = 0.;
  double num = 0.;
  double sum, sigma_k, rand_effect;

  for (k = 0; k < hsw_markov->nb_state; k++){
	
    gsl_vector *beta_k = gsl_vector_alloc(nb_covariable);
	
    for(j = 0; j < nb_covariable; j++){
      gsl_vector_set(beta_k, j, regression[k][j]);
    }
	
    rand_effect = sqrt(hsw_markov->sw_process[1]->observation[k]->random_variance[nb_random-1]);
    sum = 0.;
    sigma_k = 0.;

    for (i = 0; i < nb_sequence; i++) {
	    
      gsl_matrix *X_i = gsl_matrix_calloc(length[i], nb_covariable);
      gsl_matrix *L_ik = gsl_matrix_calloc(length[i], length[i]);
      gsl_matrix *temp_mat = gsl_matrix_calloc(length[i], hsw_markov->nb_state);
      gsl_matrix *temp_mul = gsl_matrix_calloc(nb_covariable, length[i]);
      gsl_vector *Y_i = gsl_vector_calloc(length[i]);
      gsl_vector *temp_vect = gsl_vector_calloc(length[i]);
      gsl_vector *year_i = gsl_vector_calloc(length[i]);

      for (m = 0; m < length[i]; m++){
	gsl_vector_set(year_i, m, year_effect[index[i][m]-1]);
      }

      gsl_vector_scale(year_i, rand_effect);//\tau\xi
      convert_matrix_double(covar[i], length[i], nb_covariable, X_i);
      convert_matrix_double(backward[i], length[i], hsw_markov->nb_state, temp_mat);
      convert_vector_double(real_sequence[i][0], length[i], Y_i);
      gsl_matrix_get_col(temp_vect, temp_mat, k);
      convert_matrix_diag_double(temp_vect, length[i], L_ik);
      gsl_blas_dgemv(CblasNoTrans, 1., X_i, beta_k, 0., temp_vect); // X\beta
      gsl_vector_sub(Y_i, temp_vect);// Y-X\beta
      gsl_vector_sub(Y_i, year_i);//Y- X\beta - \tau\xi
      gsl_vector_set_zero(temp_vect);
      gsl_blas_dgemv(CblasNoTrans, 1., L_ik, Y_i, 0., temp_vect);// L(Y-X\beta - \tau\xi)
	    
      for ( j = 0; j < length[i]; j++) {
	sigma_k = sigma_k + gsl_vector_get(temp_vect, j) * gsl_vector_get(Y_i, j); // (Y-X\beta - \tau\xi)L(Y-X\beta-\tau\xi)
	sum += backward[i][j][k];	// somme des L_ik(s)   
      } 	   

      gsl_vector_free(year_i);
      gsl_vector_free(temp_vect);
      gsl_vector_free(Y_i);
      gsl_matrix_free(temp_mul);
      gsl_matrix_free(temp_mat);
      gsl_matrix_free(L_ik);
      gsl_matrix_free(X_i);
 
    }
	  
    den = den + sigma_k;
    num = num + sum;

    gsl_vector_free(beta_k);

  }

  sigma = den / num;
      
  return sigma;
}


/*--------------------------------------------------------------*
 *
 *  Estimation de la variance aléatoire d'un semi-Markov switching 
 *  linear mixed model a partir de l'algorithme MCEM.
 *
 *  arguments : reference sur un objet Markov_switching, nombre et
 *              longueur des séquences, observations, nombre de 
 *              covariables, covariables, effets aléatoires temporels, 
 *              probabilités de lissage calculées dans la passe 
 *              arrière, paramètres d'index, paramètres de régression, 
 *              nombre d'effets aléatoires, état
 *
 *--------------------------------------------------------------*/

double* Switching_sequence::random_variance_temporal_stochastic(const Markov_switching *hsw_markov, int nb_sequence, 
							       int *length, double ***real_sequence, int nb_covariable, 
							       double ***covar, double *year_effect, double ***backward,
							       int **index, double *regression, int nb_random, int k) const
{
  int i, j, m;
  double *tau;
  double den = 0.;
  double num = 0.;
  gsl_vector *tau_k = gsl_vector_alloc(nb_random);
  gsl_vector *beta_k = gsl_vector_alloc(nb_covariable);
  
  for (j = 0; j < nb_covariable; j++){
    gsl_vector_set(beta_k,j, regression[j]); 
  }
	
  tau = new double[nb_random];
  
  for (i = 0; i < nb_sequence; i++) {
	      
    gsl_vector *temp_vect = gsl_vector_calloc(length[i]);
    gsl_vector *temp_mul2 = gsl_vector_calloc(length[i]);
    gsl_matrix *X_i = gsl_matrix_calloc(length[i], nb_covariable);
    gsl_matrix *L_ik = gsl_matrix_calloc(length[i], length[i]);
    gsl_matrix *temp_mat = gsl_matrix_calloc(length[i], hsw_markov->nb_state);
    gsl_matrix *temp_mul = gsl_matrix_calloc(nb_covariable, length[i]);
    gsl_vector *Y_i = gsl_vector_calloc(length[i]);
    gsl_vector *year_i = gsl_vector_calloc(length[i]);
	      
    for (m = 0; m < length[i]; m++){
      gsl_vector_set(year_i, m, year_effect[index[i][m]-1]);
    }
	      
    convert_matrix_double(covar[i], length[i], nb_covariable, X_i);
    convert_vector_double(real_sequence[i][0], length[i], Y_i);
    convert_matrix_double(backward[i], length[i], hsw_markov->nb_state, temp_mat);
    gsl_matrix_get_col(temp_vect, temp_mat, k);
    convert_matrix_diag_double(temp_vect, length[i], L_ik);
    gsl_vector_set_zero(temp_vect);
    gsl_blas_dgemv(CblasNoTrans, 1., X_i, beta_k, 0., temp_vect); // X\beta
    gsl_vector_sub(Y_i, temp_vect);// Y - X\beta 
    gsl_vector_set_zero(temp_vect);
    gsl_blas_dgemv(CblasNoTrans, 1., L_ik, Y_i, 0., temp_vect); // I(Y- X\beta)	
    gsl_blas_dgemv(CblasNoTrans, 1., L_ik, year_i, 0., temp_mul2);//I\xi   
		
    for ( j = 0; j < length[i]; j++) {
      num = num + gsl_vector_get(year_i,j) * gsl_vector_get(temp_mul2,j); // sommme des \xiX\xi
      den = den + gsl_vector_get(year_i,j) * gsl_vector_get(temp_vect,j); //somme des \xiI(Y-X\beta)    
    } 
		
    gsl_vector_free(year_i);	    
    gsl_vector_free(Y_i);
    gsl_matrix_free(temp_mul);
    gsl_matrix_free(temp_mat);
    gsl_matrix_free(L_ik);
    gsl_matrix_free(X_i);
    gsl_vector_free(temp_mul2);
    gsl_vector_free(temp_vect);
  }

  gsl_vector_set(tau_k,0, den/num);
  gsl_vector_set(tau_k,0,gsl_vector_get(tau_k,0)*gsl_vector_get(tau_k,0));
  convert_array1D_double(tau_k, nb_random, tau);

  gsl_vector_free(beta_k);
  gsl_vector_free(tau_k);
   
  return tau;
}


/*--------------------------------------------------------------*
 *
 *  Estimation de la variance aléatoire (variance commune à tous les 
 *  états) d'un semi-Markov switching linear mixed model a partir 
 *  de l'algorithme MCEM.
 *
 *  arguments : reference sur un objet Markov_switching, nombre et
 *              longueur des séquences, observations, nombre de 
 *              covariables, covariables, effets aléatoires temporels, 
 *              probabilités de lissage calculées dans la passe 
 *              arrière, paramètres d'index, paramètres de régression, 
 *              nombre d'effets aléatoires
 *
 *--------------------------------------------------------------*/

double* Switching_sequence::random_variance_temporal_VarCommune_stochastic(const Markov_switching *hsw_markov, int nb_sequence, 
									   int *length, double ***real_sequence, int nb_covariable, 
									   double ***covar, double *year_effect, double ***backward,
									   int **index, double **regression, int nb_random) const
{
  int i, j, k, m;
  double *tau;
  double den = 0.;
  double num = 0.;
  double denk, numk;
  double tau0;

  gsl_vector *tau_1 = gsl_vector_alloc(nb_random);

  tau = new double[nb_random];
  
  for (k = 0; k < hsw_markov->nb_state; k++) {
    
    gsl_vector *tau_k = gsl_vector_alloc(nb_random);
    gsl_vector *beta_k = gsl_vector_alloc(nb_covariable);

    for (j = 0; j < nb_covariable; j++){
      gsl_vector_set(beta_k,j, regression[k][j]); 
    }
	
    denk = 0.;
    numk = 0.;
	  
    for (i = 0; i < nb_sequence; i++) {
	    
      gsl_vector *temp_vect = gsl_vector_calloc(length[i]);
      gsl_vector *temp_mul2 = gsl_vector_calloc(length[i]);
      gsl_matrix *X_i = gsl_matrix_calloc(length[i], nb_covariable);
      gsl_matrix *L_ik = gsl_matrix_calloc(length[i], length[i]);
      gsl_matrix *temp_mat = gsl_matrix_calloc(length[i], hsw_markov->nb_state);
      gsl_matrix *temp_mul = gsl_matrix_calloc(nb_covariable, length[i]);
      gsl_vector *Y_i = gsl_vector_calloc(length[i]);
      gsl_vector *year_i = gsl_vector_calloc(length[i]);
	    
      for (m = 0; m < length[i]; m++){
	gsl_vector_set(year_i, m, year_effect[index[i][m]-1]);
      }
	    
      convert_matrix_double(covar[i], length[i], nb_covariable, X_i);
      convert_vector_double(real_sequence[i][0], length[i], Y_i);
      convert_matrix_double(backward[i], length[i], hsw_markov->nb_state, temp_mat);
      gsl_matrix_get_col(temp_vect, temp_mat, k);
      convert_matrix_diag_double(temp_vect, length[i], L_ik);
      gsl_vector_set_zero(temp_vect);
      gsl_blas_dgemv(CblasNoTrans, 1., X_i, beta_k, 0., temp_vect); // X\beta
      gsl_vector_sub(Y_i, temp_vect);// Y - X\beta 
      gsl_vector_set_zero(temp_vect);
      gsl_blas_dgemv(CblasNoTrans, 1., L_ik, Y_i, 0., temp_vect); // I(Y- X\beta)	
      gsl_blas_dgemv(CblasNoTrans, 1., L_ik, year_i, 0., temp_mul2);//I\xi   
	    
      for (j = 0; j < length[i]; j++) {
	numk = numk + gsl_vector_get(year_i,j) * gsl_vector_get(temp_mul2,j); // sommme des \xiX\xi
	denk = denk + gsl_vector_get(year_i,j) * gsl_vector_get(temp_vect,j); //somme des \xiI(Y-X\beta)    
      } 
			    
      gsl_vector_free(year_i);	    
      gsl_vector_free(Y_i);
      gsl_matrix_free(temp_mul);
      gsl_matrix_free(temp_mat);
      gsl_matrix_free(L_ik);
      gsl_matrix_free(X_i);
      gsl_vector_free(temp_mul2);
      gsl_vector_free(temp_vect);
    }
	 
    den = den + denk;
    num = num + numk; 
	
    gsl_vector_free(beta_k);
    gsl_vector_free(tau_k);
  }
	
  tau0 = den / num;

  gsl_vector_set(tau_1,0,tau0*tau0);
  convert_array1D_double(tau_1, nb_random, tau);

  gsl_vector_free(tau_1);
     
  return tau;
}
