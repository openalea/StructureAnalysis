#include <iostream>
#include <math.h>
//#include <lapackpp/blas1.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

using namespace std;

void convert_matrix_int (int **essai, int nb_row, int nb_col, gsl_matrix* convert_essai)
{
  register int i, j;

  for (i = 0; i < nb_row; i++) {
    for (j = 0; j < nb_col; j++) {
      gsl_matrix_set(convert_essai , i, j , essai[i][j]);  
      //    cout << gsl_matrix_get(convert_essai, i, j ) << "  ";
    }
    //   cout<<endl;
  }

}

void convert_matrix_diag_int(int *essai, int nb_row, gsl_matrix* convert_essai)
{
  register int i;
  for (i = 0; i < nb_row; i++) {
    gsl_matrix_set(convert_essai, i, i, essai[i]);
    //   cout << gsl_matrix_get(convert_essai, i, i) << "  ";
  }
  //  cout << endl;
}

void convert_vector_int (int *essai, int length, gsl_vector* essai_convert)
{
  register int i;

  for (i = 0; i < length; i++) {
    gsl_vector_set(essai_convert, i, essai[i]);
    //  cout << gsl_vector_get(essai_convert, i) << "  ";
  }
  // cout << endl ;
} 

void convert_array2D_int (gsl_matrix *essai, int nb_row, int nb_col, int **convert_essai)
{
  register int i, j;

  for (i = 0; i < nb_row ; i++) {
    for (j=0; j <nb_col; j++) {
      convert_essai[i][j] = gsl_matrix_get(essai, i, j);
      //    cout << convert_essai[i][j] << "  ";
    }
    //  cout<<endl;
  }
}

void convert_array1D_diag_int (gsl_matrix *essai, int length, int *essai_convert)
{
  register int i;
  
  for (i = 0; i < length; i++) {
    essai_convert[i] = gsl_matrix_get( essai, i, i);
    //    cout << essai_convert[i] << "  ";
  }
  //  cout << endl;
}

void convert_array1D_int (gsl_vector *essai, int length, int *essai_convert)
{
  register int i;

  for (i = 0; i < length; i++) {
    essai_convert[i] = gsl_vector_get(essai, i);
    //   cout << essai_convert[i] << "  ";
  }
  // cout << endl ;
} 

void convert_matrix_double (double **essai, int nb_row, int nb_col, gsl_matrix* convert_essai)
{
  register int i, j;

  for (i = 0; i < nb_row; i++) {
    for (j = 0; j < nb_col; j++) {
      gsl_matrix_set(convert_essai , i, j , essai[i][j]);  
      //    cout << gsl_matrix_get(convert_essai, i, j ) << "  ";
    }
    // cout<<endl;
  }
}

void convert_matrix_diag_double(gsl_vector *essai, int nb_row, gsl_matrix* convert_essai)
{
  register int i;
  for (i = 0; i < nb_row; i++) {
    gsl_matrix_set(convert_essai, i, i, gsl_vector_get(essai,i));
    // cout << gsl_matrix_get(convert_essai, i, i) << "  ";
  }
  // cout << endl;
}

void convert_vector_double (double *essai, int length, gsl_vector* essai_convert)
{
  register int i;

  for (i = 0; i < length; i++) {
    gsl_vector_set(essai_convert, i, essai[i]);
    //  cout << gsl_vector_get(essai_convert, i) << "  ";
  }
  //cout << endl ;
} 

void convert_array2D_double (gsl_matrix *essai, int nb_row, int nb_col, double **convert_essai)
{
  register int i, j;

  for (i = 0; i < nb_row ; i++) {
    for (j=0; j <nb_col; j++) {
      convert_essai[i][j] = gsl_matrix_get(essai, i, j);
      //   cout << convert_essai[i][j] << "  ";
    }
    // cout<<endl;
  }
}

void convert_array1D_diag_double (gsl_matrix *essai, int length, double *essai_convert)
{
  register int i;
  
  for (i = 0; i < length; i++) {
    essai_convert[i] = gsl_matrix_get( essai, i, i);
    //  cout << essai_convert[i] << "  ";
  }
  // cout << endl;
}

void convert_array1D_double (gsl_vector *essai, int length, double *essai_convert)
{
  register int i;

  for (i = 0; i < length; i++) {
    essai_convert[i] = gsl_vector_get(essai, i);
    // cout << essai_convert[i] << "  ";
  }
  //cout << endl ;
} 



