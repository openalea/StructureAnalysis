/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2017 CIRAD/INRA/Inria Virtual Plants
 *
 *       File author(s): Yann Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id: test.cpp 18469 2015-07-29 11:27:56Z guedon $
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



/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the Test class.
 *
 *  \param[in] iident    identifier,
 *  \param[in] ione_side flag one-sided/two-sided.
 */
/*--------------------------------------------------------------*/

Test::Test(test_distribution iident , bool ione_side)

{
  ident = iident;
  one_side = ione_side;
  df1 = I_DEFAULT;
  df2 = I_DEFAULT;
  value = -D_INF;
  critical_probability = D_DEFAULT;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the Test class.
 *
 *  \param[in] iident    identifier,
 *  \param[in] ione_side flag one-sided/two-sided,
 *  \param[in] idf1      degrees of freedom,
 *  \param[in] idf2      degrees of freedom,
 *  \param[in] ivalue    value.
 */
/*--------------------------------------------------------------*/

Test::Test(test_distribution iident , bool ione_side , int idf1 , int idf2 , double ivalue)

{
  ident = iident;
  one_side = ione_side;
  df1 = idf1;
  df2 = idf2;
  value = ivalue;
  critical_probability = D_DEFAULT;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the Test class.
 *
 *  \param[in] iident                identifier,
 *  \param[in] ione_side             flag one-sided/two-sided,
 *  \param[in] idf1                  degrees of freedom,
 *  \param[in] idf2                  degrees of freedom,
 *  \param[in] ivalue                value,
 *  \param[in] icritical_probability critical probability.
 */
/*--------------------------------------------------------------*/

Test::Test(test_distribution iident , bool ione_side , int idf1 , int idf2 , double ivalue ,
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


/*--------------------------------------------------------------*/
/**
 *  \brief Construction by copy of a Test object.
 *
 *  \param[in] test reference on a Test object.
 */
/*--------------------------------------------------------------*/

void Test::copy(const Test &test)

{
  ident = test.ident;
  one_side = test.one_side;
  df1 = test.df1;
  df2 = test.df2;
  value = test.value;
  critical_probability = test.critical_probability;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the Test class.
 *
 *  \param[in] test reference on a Test object.
 *
 *  \return         Test object.
 */
/*--------------------------------------------------------------*/

Test& Test::operator=(const Test &test)

{
  if (&test != this) {
    copy(test);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the results of a test.
 *
 *  \param[in,out] os             stream,
 *  \param[in]     comment_flag   flag comment,
 *  \param[in]     reference_flag flag reference result.
 */
/*--------------------------------------------------------------*/

ostream& Test::ascii_print(ostream &os , bool comment_flag , bool reference_flag) const

{
  if (value != 0.) {
    if (comment_flag) {
      os << "# ";
    }

    if (ident == STUDENT) {
      os << (one_side ? STAT_label[STATL_ONE_SIDED] : STAT_label[STATL_TWO_SIDED]) << " ";
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
      int i;
      Test *test;


      test = new Test(*this);

      for (i = 0;i < NB_CRITICAL_PROBABILITY;i++) {
        test->critical_probability = ref_critical_probability[i];

/*        if (test->one_side) {
          test->critical_probability = ref_critical_probability[i];
        }
        else {
          test->critical_probability = 2 * ref_critical_probability[i];
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


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the results of a test at the spreadsheet format.
 *
 *  \param[in,out] os             stream,
 *  \param[in]     reference_flag flag reference result.
 */
/*--------------------------------------------------------------*/

ostream& Test::spreadsheet_print(ostream &os , bool reference_flag) const

{
  if (value != 0.) {
    if (ident == STUDENT) {
      os << (one_side ? STAT_label[STATL_ONE_SIDED] : STAT_label[STATL_TWO_SIDED]) << " ";
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
      int i;
      Test *test;


      test = new Test(*this);

      for (i = 0;i < NB_CRITICAL_PROBABILITY;i++) {
        test->critical_probability = ref_critical_probability[i];

/*        if (test->one_side) {
          test->critical_probability = ref_critical_probability[i];
        }
        else {
          test->critical_probability = 2 * ref_critical_probability[i];
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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the critical probability from the value
 *         taken by a standard Gaussian random variable.
 */
/*--------------------------------------------------------------*/

void Test::standard_normal_critical_probability_computation()

{
  normal dist;


  if (one_side) {
    critical_probability = cdf(complement(dist , value));
  }
  else {
    critical_probability = 2 * cdf(complement(dist , fabs(value)));
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the value taken by a standard Gaussian random variable
 *         from the critical probability.
 */
/*--------------------------------------------------------------*/

void Test::standard_normal_value_computation()

{
  normal dist;


  value = quantile(complement(dist , (one_side ? critical_probability : critical_probability / 2)));
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of la critical probability from the value taken by
 *         a Chi2 random variable.
 */
/*--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the value taken by a Chi2 random variable from
 *         the critical probability.
 */
/*--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the critical probability from the value taken by
 *         a F random variable.
 */
/*--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the value taken by a F random variable from
 *         the critical probability.
 */
/*--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the critical probability from the value taken by
 *         a Student's t-random variable.
 */
/*--------------------------------------------------------------*/

void Test::t_critical_probability_computation()

{
  if (df1 > 0) {
    students_t dist(df1);


    if (one_side) {
      critical_probability = cdf(complement(dist , value));
    }
    else {
      critical_probability = 2 * cdf(complement(dist , fabs(value)));
    }
  }

  else {
    critical_probability = D_DEFAULT;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the value taken by a Student's t-random variable from
 *         the critical probability.
 */
/*--------------------------------------------------------------*/

void Test::t_value_computation()

{
  if (df1 > 0) {
    students_t dist(df1);


    value = quantile(complement(dist , (one_side ? critical_probability : critical_probability / 2)));
  }

  else {
    value = -D_INF;
  }
}


};  // namespace stat_tool
