/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2015 CIRAD/INRA/Inria Virtual Plants
 *
 *       File author(s): Yann Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id: test.cpp 18021 2015-04-23 07:07:14Z guedon $
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



#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/students_t.hpp>

#include "stat_tools.h"
#include "stat_label.h"

using namespace std;
using namespace boost::math;


namespace stat_tool {


static const double  CRITICAL_PROBABILITY_FACTOR = 1.2;



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
 *  Calcul de la probabilite critique a partir de la valeur
 *  prise par une variable normale centree reduite.
 *
 *--------------------------------------------------------------*/

void Test::standard_normal_critical_probability_computation()

{
  normal dist;


  critical_probability = cdf(complement(dist , value));

  if (!one_side) {
    critical_probability *= 2.;
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la valeur prise par une variable normale centree reduite
 *  a partir de la probabilite critique.
 *
 *--------------------------------------------------------------*/

void Test::standard_normal_value_computation()

{
  normal dist;


  value = quantile(complement(dist , (one_side ? critical_probability : critical_probability / 2.)));
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la probabilite critique a partir de la valeur prise
 *  par une variable du chi2.
 *
 *--------------------------------------------------------------*/

void Test::chi2_critical_probability_computation()

{
  if (df1 > 0) {
    chi_squared dist(df1);


    critical_probability = cdf(complement(dist , value));
  }

  else {
    critical_probability = D_DEFAULT;
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la valeur prise par une variable du chi2 a partir
 *  de la probabilite critique.
 *
 *--------------------------------------------------------------*/

void Test::chi2_value_computation()

{
  if ((df1 > 0) && (critical_probability > 0.)) {
    chi_squared dist(df1);


    value = quantile(complement(dist , critical_probability));
  }

  else {
    value = -D_INF;
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la probabilite critique a partir de la valeur prise
 *  par une variable F.
 *
 *--------------------------------------------------------------*/

void Test::F_critical_probability_computation()

{
  if ((df1 > 0) && (df2 > 0)) {
    fisher_f dist(df1 , df2);


    critical_probability = cdf(complement(dist , value));
  }

  else {
    critical_probability = D_DEFAULT;
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la valeur prise par une variable F a partir de
 *  la probabilite critique.
 *
 *--------------------------------------------------------------*/

void Test::F_value_computation()

{
  if ((df1 > 1) && (df2 > 1) && (critical_probability > 0.)) {
    fisher_f dist(df1 , df2);


    value = quantile(complement(dist , critical_probability));
  }

  else {
    value = -D_INF;
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la probabilite critique a partir de la valeur prise
 *  par une variable t.
 *
 *--------------------------------------------------------------*/

void Test::t_critical_probability_computation()

{
  if (df1 > 0) {
    students_t dist(df1);


    critical_probability = cdf(complement(dist , value));

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
 *  la probabilite critique.
 *
 *--------------------------------------------------------------*/

void Test::t_value_computation()

{
  if (df1 > 0) {
    students_t dist(df1);


    value = quantile(complement(dist , critical_probability));
  }

  else {
    value = -D_INF;
  }
}


};  // namespace stat_tool
