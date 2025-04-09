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

#include <iostream>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

using namespace std;


/*---------------------------------------------------------------
 *
 * Transformation d'une matrice d'entiers en objet de type
 * gsl_matrix
 *
 * arguments : matrice d'entiers, nombre de lignes, nombre de colonnes,
 *             objet de type gsl_matrix en sortie 
 *
 *---------------------------------------------------------------*/

void convert_matrix_int (int **matrix, int nb_row, int nb_col, gsl_matrix* convert_matrix)
{
  register int i, j;

  for (i = 0; i < nb_row; i++) {
    for (j = 0; j < nb_col; j++) {
      gsl_matrix_set(convert_matrix , i, j , matrix[i][j]);  
    }
  }

}


/*---------------------------------------------------------------
 *
 * Transformation d'une matrice diagonales d'entiers en objet de type
 * gsl_matrix
 *
 * arguments : matrice diagonale d'entiers, nombre de lignes,
 *             objet de type gsl_matrix en sortie 
 *
 *---------------------------------------------------------------*/

void convert_matrix_diag_int(int *matrix, int nb_row, gsl_matrix* convert_matrix)
{
  register int i;
  for (i = 0; i < nb_row; i++) {
    gsl_matrix_set(convert_matrix, i, i, matrix[i]);
  }
}


/*---------------------------------------------------------------
 *
 * Transformation d'un vecteur d'entiers en objet de type
 * gsl_vector
 *
 * arguments : vecteur d'entiers, nombre d'éléments du vecteur,
 *             objet de type gsl_vector en sortie 
 *
 *---------------------------------------------------------------*/

void convert_vector_int (int *vect, int length, gsl_vector* convert_vect)
{
  register int i;

  for (i = 0; i < length; i++) {
    gsl_vector_set(convert_vect, i, vect[i]);
  }
} 


/*---------------------------------------------------------------
 *
 * Transformation d'un objet de type gsl_matrix en matrice d'entiers
 *
 * arguments : objet de type gsl_matrix, nombre de lignes, 
 *             nombre de colonnes, matrice d'entiers en sortie 
 *
 *---------------------------------------------------------------*/

void convert_array2D_int (gsl_matrix *tab, int nb_row, int nb_col, int **convert_tab)
{
  register int i, j;

  for (i = 0; i < nb_row ; i++) {
    for (j=0; j < nb_col; j++) {
      convert_tab[i][j] = gsl_matrix_get(tab, i, j);
    }
  }
}


/*---------------------------------------------------------------
 *
 * Transformation d'un objet de type gsl_matrix diagonale 
 * en vecteur d'entiers
 *
 * arguments : objet de type gsl_matrix diagonale, nombre de lignes, 
 *             vecteur d'entiers en sortie 
 *
 *---------------------------------------------------------------*/

void convert_array1D_diag_int (gsl_matrix *tab, int length, int *convert_tab)
{
  register int i;
  
  for (i = 0; i < length; i++) {
    convert_tab[i] = gsl_matrix_get( tab, i, i);
  }
}


/*---------------------------------------------------------------
 *
 * Transformation d'un objet de type gsl_vector en vecteur d'entiers
 *
 * arguments : objet de type gsl_vector, nombre d'éléments du vecteur,
 *             vecteur d'entiers en sortie 
 *
 *---------------------------------------------------------------*/

void convert_array1D_int (gsl_vector *tab, int length, int *convert_tab)
{
  register int i;

  for (i = 0; i < length; i++) {
    convert_tab[i] = gsl_vector_get(tab, i);
  }
} 


/*---------------------------------------------------------------
 *
 * Transformation d'une matrice en objet de type gsl_matrix
 *
 * arguments : matrice, nombre de lignes, nombre de colonnes,
 *             objet de type gsl_matrix en sortie 
 *
 *---------------------------------------------------------------*/

void convert_matrix_double (double **matrix, int nb_row, int nb_col, gsl_matrix* convert_matrix)
{
  register int i, j;

  for (i = 0; i < nb_row; i++) {
    for (j = 0; j < nb_col; j++) {
      gsl_matrix_set(convert_matrix , i, j , matrix[i][j]);  
    }
  }
}


/*---------------------------------------------------------------
 *
 * Transformation d'une matrice diagonale en objet de type gsl_matrix
 *
 * arguments : matrice diagonale, nombre de lignes, nombre de colonnes,
 *             objet de type gsl_matrix en sortie 
 *
 *---------------------------------------------------------------*/

void convert_matrix_diag_double(gsl_vector *matrix, int nb_row, gsl_matrix* convert_matrix)
{
  register int i;
  for (i = 0; i < nb_row; i++) {
    gsl_matrix_set(convert_matrix, i, i, gsl_vector_get(matrix,i));
  }
}


/*---------------------------------------------------------------
 *
 * Transformation d'un vecteur en objet de type gsl_vector
 *
 * arguments : vecteur, nombre de lignes, nombre de colonnes,
 *             objet de type gsl_vector en sortie 
 *
 *---------------------------------------------------------------*/

void convert_vector_double (double *vect, int length, gsl_vector* convert_vect)
{
  register int i;

  for (i = 0; i < length; i++) {
    gsl_vector_set(convert_vect, i, vect[i]);
  }
} 


/*---------------------------------------------------------------
 *
 * Transformation d'un objet de type gsl_matrix en matrice
 *
 * arguments : objet de type gsl_matrix, nombre de lignes, 
 *             nombre de colonnes, matrice en sortie 
 *
 *---------------------------------------------------------------*/

void convert_array2D_double (gsl_matrix *tab, int nb_row, int nb_col, double **convert_tab)
{
  register int i, j;

  for (i = 0; i < nb_row ; i++) {
    for (j=0; j <nb_col; j++) {
      convert_tab[i][j] = gsl_matrix_get(tab, i, j);
    }
  }
}


/*---------------------------------------------------------------
 *
 * Transformation d'un objet de type gsl_matrix diagonale en vecteur
 *
 * arguments : objet de type gsl_matrix diagonale, nombre de lignes, 
 *             vecteur en sortie 
 *
 *---------------------------------------------------------------*/

void convert_array1D_diag_double (gsl_matrix *tab, int length, double *convert_tab)
{
  register int i;
  
  for (i = 0; i < length; i++) {
    convert_tab[i] = gsl_matrix_get( tab, i, i);
  }
}


/*---------------------------------------------------------------
 *
 * Transformation d'un objet de type gsl_vector en vecteur
 *
 * arguments : objet de type gsl_vector, nombre d'éléments du vecteur,
 *             vecteurs en sortie 
 *
 *---------------------------------------------------------------*/

void convert_array1D_double (gsl_vector *tab, int length, double *convert_tab)
{
  register int i;

  for (i = 0; i < length; i++) {
    convert_tab[i] = gsl_vector_get(tab, i);
  }
} 



