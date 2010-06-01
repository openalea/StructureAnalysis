/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2010 CIRAD/INRIA Virtual Plants
 *
 *       File author(s): Y. Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id$
 *
 *       Forum for V-Plants developers:
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



#include "tool/util_math.h"
#include "stat_tools.h"
#include "stat_label.h"

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/students_t.hpp>

using namespace std;
using namespace boost::math;


static const double  CRITICAL_PROBABILITY_FACTOR = 1.2;

const static double correct_term[15] = {-0.0118 , -0.0067 , -0.0033 , -0.0010 ,
                                        0.0001 , 0.0006 , 0.0006 , 0.0002 ,
                                        -0.0003 , -0.0006 , -0.0005 , 0.0002 ,
                                        0.0017 , 0.0043 , 0.0082};

const static double a_term[5] = {0.3183 , 0.4991 , 1.1094 , 3.0941 , 9.948};

const static double b_term[5] = {0.0 , 0.0518 , -0.0460 , -2.756 , -14.05};



/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Test.
 *
 *  arguments : identificateur, unilateral/bilateral.
 *
 *--------------------------------------------------------------*/

Test::Test(int iident , bool ione_side)

{
  ident = iident;
  one_side = ione_side;
  df1 = I_DEFAULT;
  df2 = I_DEFAULT;
  value = -D_INF;
  critical_probability = D_DEFAULT;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Test.
 *
 *  arguments : identificateur, unilateral/bilateral, nombres de degres de liberte ,
 *              valeur.
 *
 *--------------------------------------------------------------*/

Test::Test(int iident , bool ione_side , int idf1 , int idf2 , double ivalue)

{
  ident = iident;
  one_side = ione_side;
  df1 = idf1;
  df2 = idf2;
  value = ivalue;
  critical_probability = D_DEFAULT;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Test.
 *
 *  arguments : identificateur, unilateral/bilateral, nombres de degres de liberte,
 *              valeur, probabilite critique.
 *
 *--------------------------------------------------------------*/

Test::Test(int iident , bool ione_side , int idf1 , int idf2 , double ivalue ,
           double icritical_probability)

{
/*
  ident = iident;
  one_side = ione_side;
  df1 = idf1;
  df2 = idf2;

  value = ivalue;

  switch (ident) {
  case STANDARD_NORMAL :
    standard_normal_critical_probability_computation();
    break;
  case CHI2 :
    chi2_critical_probability_computation();
    break;
  case FISHER :
    F_critical_probability_computation();
    break;
  case STUDENT :
    t_critical_probability_computation();
    break;
  }

# ifdef DEBUG
  ascii_print(cout , false , false);
# endif

  critical_probability = icritical_probability;

  switch (ident) {
  case STANDARD_NORMAL :
    standard_normal_value_computation();
    break;
  case CHI2 :
    chi2_value_computation();
    break;
  case FISHER :
    F_value_computation();
    break;
  case STUDENT :
    t_value_computation();
    break;
  }

# ifdef DEBUG
  cout << STAT_label[STATL_REFERENCE] << " ";

  switch (ident) {
  case STANDARD_NORMAL :
    cout << STAT_label[STATL_STANDARD_NORMAL_VALUE];
    break;
  case CHI2 :
    cout << STAT_label[STATL_CHI2_VALUE];
    break;
  case FISHER :
    cout << STAT_label[STATL_F_VALUE];
    break;
  case STUDENT :
    cout << STAT_label[STATL_T_VALUE];
    break;
  }

  cout << ": " << value << "   " << STAT_label[STATL_REFERENCE] << " "
       << STAT_label[STATL_CRITICAL_PROBABILITY] << ": " << critical_probability << endl;
# endif
*/
}


/*--------------------------------------------------------------*
 *
 *  Construction par copie d'un objet Test.
 *
 *  argument : reference sur un objet Test.
 *
 *--------------------------------------------------------------*/

void Test::copy(const Test &test)

{
  ident = test.ident;
  one_side = test.one_side;
  df1 = test.df1;
  df2 = test.df2;
  value = test.value;
  critical_probability = test.critical_probability;
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Test.
 *
 *  argument : reference sur un objet Test.
 *
 *--------------------------------------------------------------*/

Test& Test::operator=(const Test &test)

{
  if (&test != this) {
    copy(test);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des resultats d'un test parametrique.
 *
 *  arguments : stream, flag commentaire, flag resultat de reference.
 *
 *--------------------------------------------------------------*/

ostream& Test::ascii_print(ostream &os , bool comment_flag , bool reference_flag) const

{
  if (value != 0.) {
    if (comment_flag) {
      os << "# ";
    }

    if (ident == STUDENT) {
      switch (one_side) {
      case true :
        os << STAT_label[STATL_ONE_SIDED] << " ";
        break;
      case false :
        os << STAT_label[STATL_TWO_SIDED] << " ";
        break;
      }
    }

    switch (ident) {
    case CHI2 :
      os << STAT_label[STATL_CHI2_TEST];
      break;
    case FISHER :
      os << STAT_label[STATL_F_TEST];
      break;
    case STUDENT :
      os << STAT_label[STATL_T_TEST];
      break;
    }

    if (df1 > 0) {
      os << " (" << df1 << " " << STAT_label[df1 == 1 ? STATL_FREEDOM_DEGREE : STATL_FREEDOM_DEGREES];
      if (df2 > 0) {
        os << ", " << df2 << " " << STAT_label[df2 == 1 ? STATL_FREEDOM_DEGREE : STATL_FREEDOM_DEGREES];
      }
      os << ")";
    }
    os << endl;

    if (comment_flag) {
      os << "# ";
    }
    switch (ident) {
    case STANDARD_NORMAL :
      os << STAT_label[STATL_STANDARD_NORMAL_VALUE];
      break;
    case CHI2 :
      os << STAT_label[STATL_CHI2_VALUE];
      break;
    case FISHER :
      os << STAT_label[STATL_F_VALUE];
      break;
    case STUDENT :
      os << STAT_label[STATL_T_VALUE];
      break;
    }
    os << ": " << value << "   "
       << STAT_label[STATL_CRITICAL_PROBABILITY] << ": " << critical_probability << endl;

    if (reference_flag) {
      register int i;
      Test *test;


      test = new Test(*this);

      for (i = 0;i < NB_CRITICAL_PROBABILITY;i++) {
        test->critical_probability = ref_critical_probability[i];

/*        switch (test->one_side) {
        case true :
          test->critical_probability = ref_critical_probability[i];
          break;
        case false :
          test->critical_probability = 2 * ref_critical_probability[i];
          break;
        } */

        switch (test->ident) {
        case STANDARD_NORMAL :
          test->standard_normal_value_computation();
          break;
        case CHI2 :
          test->chi2_value_computation();
          break;
        case FISHER :
          test->F_value_computation();
          break;
        case STUDENT :
          test->t_value_computation();
          break;
        }

        if (comment_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_REFERENCE] << " ";
        switch (test->ident) {
        case STANDARD_NORMAL :
          os << STAT_label[STATL_STANDARD_NORMAL_VALUE];
          break;
        case CHI2 :
          os << STAT_label[STATL_CHI2_VALUE];
          break;
        case FISHER :
          os << STAT_label[STATL_F_VALUE];
          break;
        case STUDENT :
          os << STAT_label[STATL_T_VALUE];
          break;
        }
        os << ": " << test->value << "   " << STAT_label[STATL_REFERENCE] << " "
           << STAT_label[STATL_CRITICAL_PROBABILITY] << ": " << test->critical_probability << endl;

        if (critical_probability > test->critical_probability * CRITICAL_PROBABILITY_FACTOR) {
          break;
        }
      }

      delete test;
    }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des resultats d'un test parametrique au format tableur.
 *
 *  arguments : stream, flag resultat de reference.
 *
 *--------------------------------------------------------------*/

ostream& Test::spreadsheet_print(ostream &os , bool reference_flag) const

{
  if (value != 0.) {
    if (ident == STUDENT) {
      switch (one_side) {
      case true :
        os << STAT_label[STATL_ONE_SIDED] << " ";
        break;
      case false :
        os << STAT_label[STATL_TWO_SIDED] << " ";
        break;
      }
    }

    switch (ident) {
    case CHI2 :
      os << STAT_label[STATL_CHI2_TEST];
      break;
    case FISHER :
      os << STAT_label[STATL_F_TEST];
      break;
    case STUDENT :
      os << STAT_label[STATL_T_TEST];
      break;
    }

    if (df1 > 0) {
      os << "\t" << df1 << "\t" << STAT_label[df1 == 1 ? STATL_FREEDOM_DEGREE : STATL_FREEDOM_DEGREES];
      if (df2 > 0) {
        os << "\t\t" << df2 << "\t" << STAT_label[df2 == 1 ? STATL_FREEDOM_DEGREE : STATL_FREEDOM_DEGREES];
      }
    }
    os << endl;

    switch (ident) {
    case STANDARD_NORMAL :
      os << STAT_label[STATL_STANDARD_NORMAL_VALUE];
      break;
    case CHI2 :
      os << STAT_label[STATL_CHI2_VALUE];
      break;
    case FISHER :
      os << STAT_label[STATL_F_VALUE];
      break;
    case STUDENT :
      os << STAT_label[STATL_T_VALUE];
      break;
    }
    os << "\t" << value << "\t\t"
       << STAT_label[STATL_CRITICAL_PROBABILITY] << "\t" << critical_probability << endl;

    if (reference_flag) {
      register int i;
      Test *test;


      test = new Test(*this);

      for (i = 0;i < NB_CRITICAL_PROBABILITY;i++) {
        test->critical_probability = ref_critical_probability[i];

/*        switch (test->one_side) {
        case true :
          test->critical_probability = ref_critical_probability[i];
          break;
        case false :
          test->critical_probability = 2 * ref_critical_probability[i];
          break;
        } */

        switch (test->ident) {
        case STANDARD_NORMAL :
          test->standard_normal_value_computation();
          break;
        case CHI2 :
          test->chi2_value_computation();
          break;
        case FISHER :
          test->F_value_computation();
          break;
        case STUDENT :
          test->t_value_computation();
          break;
        }

        os << STAT_label[STATL_REFERENCE] << " ";
        switch (test->ident) {
          case STANDARD_NORMAL :
          os << STAT_label[STATL_STANDARD_NORMAL_VALUE];
          break;
        case CHI2 :
          os << STAT_label[STATL_CHI2_VALUE];
          break;
        case FISHER :
          os << STAT_label[STATL_F_VALUE];
          break;
        case STUDENT :
          os << STAT_label[STATL_T_VALUE];
          break;
        }
        os << "\t" << test->value << "\t\t" << STAT_label[STATL_REFERENCE] << " "
           << STAT_label[STATL_CRITICAL_PROBABILITY] << "\t" << test->critical_probability << endl;

        if (critical_probability > test->critical_probability * CRITICAL_PROBABILITY_FACTOR) {
          break;
        }
      }

      delete test;
    }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la probabilite critique a partir de la valeur prise
 *  par une variable normale centree reduite (formule 26.2.19,
 *  Handbook of Mathemetical Functions (M. Abramowitz & I.A. Stegum)).
 *
 *  argument : valeur prise par une variable normale centree reduite.
 *
 *--------------------------------------------------------------*/

double standard_normal_critical_probability_computation(double value)

{
  register int i;
  double critical_probability , term , var[7];


  var[1] = fabs(value);
  for (i = 2;i <= 6;i++) {
    var[i] = var[i - 1] * var[1];
  }

  term = 1. + 0.0498673470 * var[1] + 0.0211410061 * var[2] +
         0.0032776263 * var[3] + 0.0000380036 * var[4] +
         0.0000488906 * var[5] + 0.0000053830 * var[6];

  var[0] = 1.;
  for (i = 0;i < 16;i++) {
    var[0] *= term;
  }
  critical_probability = 0.5 / var[0];

  if (value < 0.) {
    critical_probability = 1. - critical_probability;
  }

  return critical_probability;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la probabilite critique a partir de la valeur
 *  prise par une variable normale centree reduite.
 *
 *--------------------------------------------------------------*/

void Test::standard_normal_critical_probability_computation()

{
  critical_probability = ::standard_normal_critical_probability_computation(value);

# ifdef DEBUG
  normal dist;

  cout << "\nTEST Gaussian distribution: " << cdf(complement(dist , value))
       << " | " <<  critical_probability << endl;
# endif

  if (!one_side) {
    critical_probability *= 2.;
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la valeur prise par une variable normale centree reduite
 *  a partir de la probabilite critique (formule 26.2.23,
 *  Handbook of Mathemetical Functions (M. Abramowitz & I.A. Stegum)).
 *
 *  argument : "queue" de la fonction de repartition.
 *
 *--------------------------------------------------------------*/

double standard_normal_value_computation(double critical_probability)

{
  register int i;
  double value , var[4];


  if (critical_probability <= 0.5) {
    var[1] = sqrt(log(1. / (critical_probability * critical_probability)));
  }
  else {
    var[1] = sqrt(log(1. / ((1. - critical_probability) * (1. - critical_probability))));
  }
  for (i = 2;i <= 3;i++) {
    var[i] = var[i - 1] * var[1];
  }

  value = var[1] - (2.515517 + 0.802853 * var[1] + 0.010328 * var[2]) /
                   (1. + 1.432788 * var[1] + 0.189269 * var[2] + 0.001308 * var[3]);

  if (critical_probability > 0.5) {
    value = -value;
  }

  return value;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la valeur prise par une variable normale centree reduite
 *  a partir de la probabilite critique.
 *
 *--------------------------------------------------------------*/

void Test::standard_normal_value_computation()

{
  value = ::standard_normal_value_computation(one_side ? critical_probability : critical_probability / 2.);

# ifdef DEBUG
  normal dist;

  cout << "\nTEST Gaussian distribution: " << quantile(complement(dist , (one_side ? critical_probability : critical_probability / 2.)))
       << " | " << value << endl;
# endif

}


/*--------------------------------------------------------------*
 *
 *  Calcul de la probabilite critique a partir de la valeur prise
 *  par une variable du chi2 (d'apres Handbook of Mathemetical Functions
 *  (M. Abramowitz & I.A. Stegum)).
 *
 *--------------------------------------------------------------*/

void Test::chi2_critical_probability_computation()

{
  register int i;
  double normal_var , correct;


  if (df1 > 0) {

    // calcul de la variable normale centree reduite correspondant
    // a la variable du chi2 (formule 26.4.14)

#   ifdef _WIN32
    normal_var = (pow((value / df1) , 1. / 3.) - (1. - 2. / (9. * df1))) / sqrt(2. / (9. * df1));
#   else
    normal_var = (cbrt(value / df1) - (1. - 2. / (9. * df1))) / sqrt(2. / (9. * df1));
#   endif

    // calcul du terme de correction par interpolation lineaire (table 26.4.15)

    i = (int)(normal_var * 2 + 7);
    if (i < 0) {
      correct = correct_term[0];
    }
    else if (i >= 14) {
      correct = correct_term[14];
    }
    else {
      correct = correct_term[i] * (normal_var * 2 + 7 - i) +
                correct_term[i + 1] * (i + 1 - (normal_var * 2 + 7));
    }
    correct *= 60. / (double)df1;
 
    critical_probability = ::standard_normal_critical_probability_computation(normal_var + correct);

#   ifdef DEBUG
    chi_squared dist(df1);

    cout << "\nTEST Chi2 distribution: " << cdf(complement(dist , value))
         << " | " <<  critical_probability << endl;
#   endif

  }

  else {
    critical_probability = D_DEFAULT;
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la valeur prise par une variable du chi2 a partir
 *  de la probabilite critique (d'apres Handbook of Mathemetical Functions
 *  (M. Abramowitz & I.A. Stegum)).
 *
 *--------------------------------------------------------------*/

void Test::chi2_value_computation()

{
  register int i;
  double normal_var , correct , var;


  if ((df1 > 0) && (critical_probability > 0.)) {
    normal_var = ::standard_normal_value_computation(critical_probability);

    // calcul du terme de correction par interpolation lineaire (table 26.4.15)

    i = (int)(normal_var * 2 + 7);
    if (i < 0) {
      correct = correct_term[0];
    }
    else if (i >= 14) {
      correct = correct_term[14];
    }
    else {
      correct = correct_term[i] * (normal_var * 2 + 7 - i) +
                correct_term[i + 1] * (i + 1 - (normal_var * 2 + 7));
    }
    correct *= 60. / (double)df1;

    // approximation cubique finale (formule 26.4.18)

    var = 1. - 2. / (9. * df1) +
          (normal_var - correct) * sqrt(2. / (9. * df1));
    value = df1 * var * var * var;

#   ifdef DEBUG
    chi_squared dist(df1);

    cout << "\nTEST Chi2 distribution: " << quantile(complement(dist , critical_probability))
         << " | " << value << endl;
#   endif

  }

  else {
    value = -D_INF;
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la probabilite critique a partir de la valeur prise
 *  par une variable F (d'apres Handbook of Mathemetical Functions
 *  (M. Abramowitz & I.A. Stegum)).
 *
 *--------------------------------------------------------------*/

void Test::F_critical_probability_computation()

{
  double normal_var , cbrt_value;


  if ((df1 > 0) && (df2 > 0)) {

    // calcul de la variable normale centree reduite correspondant
    // a la variable F (formule 26.6.15)

#   ifdef _WIN32
    cbrt_value = pow(value , 1. / 3.);
#   else
    cbrt_value = cbrt(value);
#   endif

    normal_var = (cbrt_value * (1. - 2. / (9. * df2)) - (1. - 2. / (9. * df1))) /
                 sqrt(2. / (9. * df1) + cbrt_value * cbrt_value * 2. / (9. * df2));

    critical_probability = ::standard_normal_critical_probability_computation(normal_var);

#   ifdef DEBUG
    fisher_f dist(df1 , df2);

    cout << "\nTEST F-distribution: " << cdf(complement(dist , value))
         << " | " <<  critical_probability << endl;
#   endif

  }

  else {
    critical_probability = D_DEFAULT;
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la valeur prise par une variable F a partir de
 *  la probabilite critique (d'apres Handbook of Mathemetical Functions
 *  (M. Abramowitz & I.A. Stegum)).
 *
 *--------------------------------------------------------------*/

void Test::F_value_computation()

{
  double normal_var , h , lambda , w;


  if ((df1 > 1) && (df2 > 1) && (critical_probability > 0.)) {
    normal_var = ::standard_normal_value_computation(critical_probability);

    // approximation par la fonction Beta incomplete (formule 26.5.22)

    h = 2. / (1. / (double)(df2 - 1) + 1. / (double)(df1 - 1));
    lambda = (normal_var * normal_var - 3.) / 6.;
    w = (normal_var * sqrt(h + lambda)) / h -
        (1. / (double)(df1 - 1) - 1. / (double)(df2 - 1)) *
        (lambda + 5. / 6. - 2. / (3. * h));

    // approximation finale (formule 26.6.16)

    value = exp(2. * w);

#   ifdef DEBUG
    fisher_f dist(df1 , df2);

    cout << "\nTEST F-distribution: " << quantile(complement(dist , critical_probability))
         << " | " << value << endl;
#   endif

  }

  else {
    value = -D_INF;
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la probabilite critique a partir de la valeur prise
 *  par une variable t (d'apres Handbook of Mathemetical Functions
 *  (M. Abramowitz & I.A. Stegum)).
 *
 *--------------------------------------------------------------*/

void Test::t_critical_probability_computation()

{
  register int i;
  double abs_value , buff , normal_var;


  if (df1 > 0) {
    abs_value = fabs(value);

    if ((df1 <= 5) && (abs_value > MIN_T_VALUE)) {

      // calcul directe de la queue de la fonction de repartition (formule 26.7.7)

      buff = 1.;
      for (i = 0;i < df1;i++) {
        buff *= abs_value;
      }

      critical_probability = a_term[df1 - 1] / buff + b_term[df1 - 1] / (buff * abs_value);
    }

    else {

      // calcul de la variable normale centree reduite correspondant
      // a la variable t (formule 26.7.8)

      normal_var = (abs_value * (1. - 1. / (4. * df1))) /
                   sqrt(1. + (abs_value * abs_value) / (2. * df1));

      critical_probability = ::standard_normal_critical_probability_computation(normal_var);
    }

#   ifdef DEBUG
    students_t dist(df1);

    cout << "\nTEST t-distribution: " << cdf(complement(dist , value))
         << " | " <<  critical_probability << endl;
#   endif

    switch (one_side) {

    case true : {
      if (value < 0.) {
        critical_probability = 1. - critical_probability;
      }
      break;
    }

    case false : {
      critical_probability *= 2.;
      break;
    }
    }
  }

  else {
    critical_probability = D_DEFAULT;
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la valeur prise par une variable t a partir de
 *  la probabilite critique (d'apres Handbook of Mathemetical Functions
 *  (M. Abramowitz & I.A. Stegum)).
 *
 *  arguments : unilateral/bilateral, nombres de degres de liberte, probabilite critique.
 *
 *--------------------------------------------------------------*/

double t_value_computation(bool one_side , int df , double critical_probability)

{
  register int i;
  double normal_var , value , var[10];


  if ((df > 0) && (critical_probability > 0.)) {
    normal_var = ::standard_normal_value_computation(one_side ? critical_probability : critical_probability / 2.);

    // (formule 26.7.5)

    var[1] = normal_var;
    for (i = 2;i <= 9;i++) {
      var[i] = var[i - 1] * var[1];
    }

    value = var[1] + (var[3] + var[1]) / (4 * (double)df) +
            (5 * var[5] + 16 * var[3] + 3 * var[1]) / (96 * (double)df * (double)df) +
            (3 * var[7] + 19 * var[5] + 17 * var[3] - 15 * var[1]) /
            (384 * (double)df * (double)df * (double)df) +
            (79 * var[9] + 776 * var[7] + 1482 * var[5] - 1920 * var[3] - 945 * var[1]) /
            (92160 * (double)df * (double)df * (double)df * (double)df);
  }

  else {
    value = -D_INF;
  }

  return value;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la valeur prise par une variable t a partir de
 *  la probabilite critique.
 *
 *--------------------------------------------------------------*/

void Test::t_value_computation()

{
  value = ::t_value_computation(one_side , df1 , critical_probability);

# ifdef DEBUG
  students_t dist(df1);

  cout << "\nTEST t-distribution: " << quantile(complement(dist , critical_probability))
       << " | " << value << endl;
# endif

}
