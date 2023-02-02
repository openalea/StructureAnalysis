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



#include <limits.h>
#include <math.h>

#include <string>
#include <vector>
#include <sstream>
#include <iomanip>

#include <boost/tokenizer.hpp>

#include "quantile_computation.hpp"

#include "vectors.h"
#include "stat_label.h"

// #include "quantile_computation.h"   problem compiler C++ Windows

using namespace std;
using namespace boost;


namespace stat_tool {



/*--------------------------------------------------------------*/
/**
 *  \brief Default constructor of the Vectors class.
 */
/*--------------------------------------------------------------*/

Vectors::Vectors()

{
  nb_vector = 0;
  identifier = NULL;

  nb_variable = 0;

  type = NULL;
  min_value = NULL;
  max_value = NULL;
  min_interval = NULL;
  marginal_distribution = NULL;
  marginal_histogram = NULL;

  mean = NULL;
  covariance = NULL;

  int_vector = NULL;
  real_vector = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Initialization of a Vectors object.
 *
 *  \param[in] inb_vector   number of individuals,
 *  \param[in] iidentifier  individual identifiers,
 *  \param[in] inb_variable number of variables,
 *  \param[in] itype        variable types,
 *  \param[in] init_flag    flag initialization.
 */
/*--------------------------------------------------------------*/

void Vectors::init(int inb_vector , int *iidentifier , int inb_variable ,
                   variable_nature *itype , bool init_flag)

{
  int i , j;


  nb_vector = inb_vector;

  identifier = new int[nb_vector];

  if (iidentifier) {
    for (i = 0;i < nb_vector;i++) {
      identifier[i] = iidentifier[i];
    }
  }

  else {
    for (i = 0;i < nb_vector;i++) {
      identifier[i] = i + 1;
    }
  }

  nb_variable = inb_variable;

  type = new variable_nature[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  min_interval = new double[nb_variable];
  marginal_distribution = new FrequencyDistribution*[nb_variable];
  marginal_histogram = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    type[i] = itype[i];
    min_value[i] = 0.;
    max_value[i] = 0.;
    min_interval[i] = 0.;
    marginal_distribution[i] = NULL;
    marginal_histogram[i] = NULL;
  }

  mean = new double[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    mean[i] = D_INF;
  }

  covariance = new double*[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    covariance[i] = new double[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      covariance[i][j] = D_DEFAULT;
    }
  }

  int_vector = new int*[nb_vector];
  real_vector = new double*[nb_vector];
  for (i = 0;i < nb_vector;i++) {
    int_vector[i] = new int[nb_variable];
    real_vector[i] = new double[nb_variable];

    if (init_flag) {
      for (j = 0;j < nb_variable;j++) {
        int_vector[i][j] = 0;
        real_vector[i][j] = 0.;
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the Vectors class.
 *
 *  \param[in] inb_vector   number of individuals,
 *  \param[in] iidentifier  individual identifiers,
 *  \param[in] inb_variable number of variables,
 *  \param[in] iint_vector  integer-valued vectors.
 */
/*--------------------------------------------------------------*/

Vectors::Vectors(int inb_vector , int *iidentifier , int inb_variable ,
                 int **iint_vector)

{
  int i , j;
  variable_nature *itype;


  itype = new variable_nature[inb_variable];
  for (i = 0;i < inb_variable;i++) {
    itype[i] = INT_VALUE;
  }

  init(inb_vector , iidentifier , inb_variable , itype , false);
  delete [] itype;

  for (i = 0;i < nb_vector;i++) {
    for (j = 0;j < nb_variable;j++) {
      int_vector[i][j] = iint_vector[i][j];
    }
  }

  for (i = 0;i < nb_variable;i++) {
    min_value_computation(i);
    max_value_computation(i);

    build_marginal_frequency_distribution(i);
    min_interval_computation(i);
  }

  covariance_computation();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the Vectors class.
 *
 *  \param[in] inb_vector   number of individuals,
 *  \param[in] iidentifier  individual identifiers,
 *  \param[in] inb_variable number of variables,
 *  \param[in] ireal_vector real-valued vectors.
 */
/*--------------------------------------------------------------*/

Vectors::Vectors(int inb_vector , int *iidentifier , int inb_variable ,
                 double **ireal_vector)

{
  int i , j;
  variable_nature *itype;


  itype = new variable_nature[inb_variable];
  for (i = 0;i < inb_variable;i++) {
    itype[i] = REAL_VALUE;
  }

  init(inb_vector , iidentifier , inb_variable , itype , false);
  delete [] itype;

  for (i = 0;i < nb_vector;i++) {
    for (j = 0;j < nb_variable;j++) {
      real_vector[i][j] = ireal_vector[i][j];
    }
  }

  for (i = 0;i < nb_variable;i++) {
    min_value_computation(i);
    max_value_computation(i);

    build_marginal_histogram(i);
    min_interval_computation(i);

    mean_computation(i);
    variance_computation(i);
  }

  covariance_computation();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the Vectors class.
 *
 *  \param[in] inb_vector     number of individuals,
 *  \param[in] iidentifier    individual identifiers,
 *  \param[in] inb_variable   number of variables,
 *  \param[in] itype          variable types,
 *  \param[in] iint_vector    integer variables,
 *  \param[in] ireal_vector   real variables,
 *  \param[in] variable_index variable indexing.
 */
/*--------------------------------------------------------------*/

Vectors::Vectors(int inb_vector , int *iidentifier , int inb_variable , variable_nature *itype ,
                 int **iint_vector , double **ireal_vector , bool variable_index)

{
  int i , j , k , m;


  init(inb_vector , iidentifier , inb_variable , itype , false);

  if (variable_index) {
    for (i = 0;i < nb_variable;i++) {
      if (type[i] != REAL_VALUE) {
        for (j = 0;j < nb_vector;j++) {
          int_vector[j][i] = iint_vector[j][i];
        }
      }

      else {
        for (j = 0;j < nb_vector;j++) {
          real_vector[j][i] = ireal_vector[j][i];
        }
      }
    }
  }

  else {
    i = 0;
    j = 0;
    for (k = 0;k < nb_variable;k++) {
      if (type[k] != REAL_VALUE) {
        for (m = 0;m < nb_vector;m++) {
          int_vector[m][k] = iint_vector[m][i];
        }
        i++;
      }

      else {
        for (m = 0;m < nb_vector;m++) {
          real_vector[m][k] = ireal_vector[m][j];
        }
        j++;
      }
    }
  }

  for (i = 0;i < nb_variable;i++) {
    min_value_computation(i);
    max_value_computation(i);

    if (type[i] != REAL_VALUE) {
      build_marginal_frequency_distribution(i);
    }

    else {
      build_marginal_histogram(i);

      mean_computation(i);
      variance_computation(i);
    }

    min_interval_computation(i);
  }

  covariance_computation();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the Vectors class.
 *
 *  \param[in] vec      reference on a Vectors object,
 *  \param[in] variable variable index,
 *  \param[in] itype    variable type.
 */
/*--------------------------------------------------------------*/

Vectors::Vectors(const Vectors &vec , int variable , variable_nature itype)

{
  int i , j;


  nb_vector = vec.nb_vector;

  identifier = new int[nb_vector];
  for (i = 0;i < nb_vector;i++) {
    identifier[i] = vec.identifier[i];
  }

  nb_variable = vec.nb_variable;

  type = new variable_nature[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  min_interval = new double[nb_variable];
  marginal_distribution = new FrequencyDistribution*[nb_variable];
  marginal_histogram = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    if (i != variable) {
      type[i] = vec.type[i];
      min_value[i] = vec.min_value[i];
      max_value[i] = vec.max_value[i];
      min_interval[i] = vec.min_interval[i];

      if (vec.marginal_distribution[i]) {
        marginal_distribution[i] = new FrequencyDistribution(*(vec.marginal_distribution[i]));
      }
      else {
        marginal_distribution[i] = NULL;
      }

      if (vec.marginal_histogram[i]) {
        marginal_histogram[i] = new Histogram(*(vec.marginal_histogram[i]));
      }
      else {
        marginal_histogram[i] = NULL;
      }
    }

    else {
      type[i] = itype;
      min_value[i] = 0.;
      max_value[i] = 0.;
      min_interval[i] = 0.;
      marginal_distribution[i] = NULL;
      marginal_histogram[i] = NULL;
    }
  }

  mean = new double[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    if (i != variable) {
      mean[i] = vec.mean[i];
    }
    else {
      mean[i] = D_INF;
    }
  }

  covariance = new double*[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    covariance[i] = new double[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      if ((i != variable) && (j != variable)) {
        covariance[i][j] = vec.covariance[i][j];
      }
      else {
        covariance[i][j] = D_DEFAULT;
      }
    }
  }

  int_vector = new int*[nb_vector];
  real_vector = new double*[nb_vector];
  for (i = 0;i < nb_vector;i++) {
    int_vector[i] = new int[nb_variable];
    real_vector[i] = new double[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      if (j != variable) {
        int_vector[i][j] = vec.int_vector[i][j];
        real_vector[i][j] = vec.real_vector[i][j];
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the Vectors class.
 *
 *  \param[in] vec        reference on a Vectors object,
 *  \param[in] inb_vector number of individuals,
 *  \param[in] index      indices of the selected individuals.
 */
/*--------------------------------------------------------------*/

Vectors::Vectors(const Vectors &vec , int inb_vector , int *index)

{
  int i , j;
  int *pivector , *civector;


  nb_vector = inb_vector;

  identifier = new int[nb_vector];
  for (i = 0;i < nb_vector;i++) {
    identifier[i] = vec.identifier[index[i]];
  }

  nb_variable = vec.nb_variable;

  type = new variable_nature[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  min_interval = new double[nb_variable];
  marginal_distribution = new FrequencyDistribution*[nb_variable];
  marginal_histogram = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    type[i] = vec.type[i];
    marginal_distribution[i] = NULL;
    marginal_histogram[i] = NULL;
  }

  mean = new double[nb_variable];

  covariance = new double*[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    covariance[i] = new double[nb_variable];
  }

  int_vector = new int*[nb_vector];
  real_vector = new double*[nb_vector];
  for (i = 0;i < nb_vector;i++) {
    int_vector[i] = new int[nb_variable];
    real_vector[i] = new double[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      int_vector[i][j] = vec.int_vector[index[i]][j];
      real_vector[i][j] = vec.real_vector[index[i]][j];
    }
  }

  for (i = 0;i < nb_variable;i++) {
    min_value_computation(i);
    max_value_computation(i);

    if (type[i] != REAL_VALUE) {
      build_marginal_frequency_distribution(i);
    }

    else  {
      build_marginal_histogram(i , vec.marginal_histogram[i]->bin_width);

      mean_computation(i);
      variance_computation(i);
    }

    min_interval_computation(i);
  }

  covariance_computation();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a Vectors object.
 *
 *  \param[in] vec reference on a Vectors object.
 */
/*--------------------------------------------------------------*/

void Vectors::copy(const Vectors &vec)

{
  int i , j;


  nb_vector = vec.nb_vector;

  identifier = new int[nb_vector];
  for (i = 0;i < nb_vector;i++) {
    identifier[i] = vec.identifier[i];
  }

  nb_variable = vec.nb_variable;

  type = new variable_nature[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  min_interval = new double[nb_variable];
  marginal_distribution = new FrequencyDistribution*[nb_variable];
  marginal_histogram = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    type[i] = vec.type[i];
    min_value[i] = vec.min_value[i];
    max_value[i] = vec.max_value[i];
    min_interval[i] = vec.min_interval[i];

    if (vec.marginal_distribution[i]) {
      marginal_distribution[i] = new FrequencyDistribution(*(vec.marginal_distribution[i]));
    }
    else {
      marginal_distribution[i] = NULL;
    }

    if (vec.marginal_histogram[i]) {
      marginal_histogram[i] = new Histogram(*(vec.marginal_histogram[i]));
    }
    else {
      marginal_histogram[i] = NULL;
    }
  }

  mean = new double[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    mean[i] = vec.mean[i];
  }

  covariance = new double*[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    covariance[i] = new double[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      covariance[i][j] = vec.covariance[i][j];
    }
  }

  int_vector = new int*[nb_vector];
  real_vector = new double*[nb_vector];
  for (i = 0;i < nb_vector;i++) {
    int_vector[i] = new int[nb_variable];
    real_vector[i] = new double[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      int_vector[i][j] = vec.int_vector[i][j];
      real_vector[i][j] = vec.real_vector[i][j];
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a Vectors object adding a state variable.
 *
 *  \param[in] vec reference on a Vectors object.
 */
/*--------------------------------------------------------------*/

void Vectors::add_state_variable(const Vectors &vec)

{
  int i , j;


  nb_vector = vec.nb_vector;

  identifier = new int[nb_vector];
  for (i = 0;i < nb_vector;i++) {
    identifier[i] = vec.identifier[i];
  }

  nb_variable = vec.nb_variable + 1;

  type = new variable_nature[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  min_interval = new double[nb_variable];
  marginal_distribution = new FrequencyDistribution*[nb_variable];
  marginal_histogram = new Histogram*[nb_variable];

  type[0] = STATE;
  min_value[0] = 0.;
  max_value[0] = 0.;
  min_interval[0] = 0.; 
  marginal_distribution[0] = NULL;
  marginal_histogram[0] = NULL;

  for (i = 0;i < vec.nb_variable;i++) {
    type[i + 1] = vec.type[i];
    min_value[i + 1] = vec.min_value[i];
    max_value[i + 1] = vec.max_value[i];
    min_interval[i + 1] = vec.min_interval[i];

    if (vec.marginal_distribution[i]) {
      marginal_distribution[i + 1] = new FrequencyDistribution(*(vec.marginal_distribution[i]));
    }
    else {
      marginal_distribution[i + 1] = NULL;
    }

    if (vec.marginal_histogram[i]) {
      marginal_histogram[i + 1] = new Histogram(*(vec.marginal_histogram[i]));
    }
    else {
      marginal_histogram[i + 1] = NULL;
    }
  }

  mean = new double[nb_variable];
  for (i = 0;i < vec.nb_variable;i++) {
    mean[i + 1] = vec.mean[i];
  }

  covariance = new double*[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    covariance[i] = new double[nb_variable];
  }

  for (i = 0;i < vec.nb_variable;i++) {
    for (j = 0;j < vec.nb_variable;j++) {
      covariance[i + 1][j + 1] = vec.covariance[i][j];
    }
  }

  int_vector = new int*[nb_vector];
  real_vector = new double*[nb_vector];
  for (i = 0;i < nb_vector;i++) {
    int_vector[i] = new int[nb_variable];
    real_vector[i] = new double[nb_variable];
    for (j = 0;j < vec.nb_variable;j++) {
      int_vector[i][j + 1] = vec.int_vector[i][j];
      real_vector[i][j + 1] = vec.real_vector[i][j];
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor by copy of the Vectors class.
 *
 *  \param[in] vec       reference on a Vectors object,
 *  \param[in] transform type of transform.
 */
/*--------------------------------------------------------------*/

Vectors::Vectors(const Vectors &vec , vector_transformation transform)

{
  switch (transform) {
  case ADD_COMPONENT_VARIABLE :
    Vectors::add_state_variable(vec);
    break;
  default :
    Vectors::copy(vec);
    break;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destruction of the data members of a Vectors object.
 */
/*--------------------------------------------------------------*/

void Vectors::remove()

{
  int i;


  delete [] identifier;

  delete [] type;

  delete [] min_value;
  delete [] max_value;
  delete [] min_interval;

  if (marginal_distribution) {
    for (i = 0;i < nb_variable;i++) {
      delete marginal_distribution[i];
    }
    delete [] marginal_distribution;
  }

  if (marginal_histogram) {
    for (i = 0;i < nb_variable;i++) {
      delete marginal_histogram[i];
    }
    delete [] marginal_histogram;
  }

  delete [] mean;

  if (covariance) {
    for (i = 0;i < nb_variable;i++) {
      delete [] covariance[i];
    }
    delete [] covariance;
  }

  if (int_vector) {
    for (i = 0;i < nb_vector;i++) {
      delete [] int_vector[i];
    }
    delete [] int_vector;
  }

  if (real_vector) {
    for (i = 0;i < nb_vector;i++) {
      delete [] real_vector[i];
    }
    delete [] real_vector;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the Vectors class.
 */
/*--------------------------------------------------------------*/

Vectors::~Vectors()

{
  remove();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the Vectors class.
 *
 *  \param[in] vec reference on a Vectors object.
 *
 *  \return        Vectors object.
 */
/*--------------------------------------------------------------*/

Vectors& Vectors::operator=(const Vectors &vec)

{
  if (&vec != this) {
    remove();
    copy(vec);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of an integer-valued variable as a real-valued variable.
 *
 *  \param[in] variable variable index.
 */
/*--------------------------------------------------------------*/

void Vectors::build_real_vector(int variable)

{
  int i , j;


  for (i = 0;i < nb_variable;i++) {
    if (((variable == I_DEFAULT) || (variable == i)) && (type[i] == INT_VALUE)) {
      for (j = 0;j < nb_vector;j++) {
        real_vector[j][i] = int_vector[j][i];
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of the marginal frequency distribution for a positive integer-valued variable.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] variable variable index.
 *
 *  \return             DiscreteDistributionData object.
 */
/*--------------------------------------------------------------*/

DiscreteDistributionData* Vectors::extract(StatError &error , int variable) const

{
  bool status = true;
  DiscreteDistributionData *histo;


  histo = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if ((type[variable] != INT_VALUE) && (type[variable] != STATE)) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[STATE];
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }

    else if (!marginal_distribution[variable]) {
      status = false;
      error.update(STAT_error[STATR_MARGINAL_FREQUENCY_DISTRIBUTION]);
    }
  }

  if (status) {
    histo = new DiscreteDistributionData(*marginal_distribution[variable]);
  }

  return histo;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Checking of a list of individual identifiers (unique identifiers).
 *
 *  \param[in] error         reference on a StatError object,
 *  \param[in] nb_individual number of individuals,
 *  \param[in] identifier    individual identifiers.
 *
 *  \return                  error status.
 */
/*--------------------------------------------------------------*/

bool identifier_checking(StatError &error , int nb_individual , int *identifier)

{
  bool status = true , *selected_identifier;
  int i;
  int max_identifier;


  max_identifier = identifier[0];
  for (i = 1;i < nb_individual;i++) {
    if (identifier[i] > max_identifier) {
      max_identifier = identifier[i];
    }
  }

  selected_identifier = new bool[max_identifier + 1];
  for (i = 0;i <= max_identifier;i++) {
    selected_identifier[i] = false;
  }

  for (i = 0;i < nb_individual;i++) {
    if (selected_identifier[identifier[i]]) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_IDENTIFIER] << " " << identifier[i] << " "
                    << STAT_error[STATR_ALREADY_USED];
      error.update((error_message.str()).c_str());
    }
    else {
      selected_identifier[identifier[i]] = true;
    }
  }

  delete [] selected_identifier;

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Checking of the identifiers of a Vectors object.
 *
 *  \param[in] error reference on a StatError object.
 *
 *  \return          error status.
 */
/*--------------------------------------------------------------*/

bool Vectors::check(StatError &error)

{
  bool status = true , lstatus;


  error.init();

  if (nb_variable > VECTOR_NB_VARIABLE) {
    status = false;
    error.update(STAT_error[STATR_NB_VARIABLE]);
  }

  lstatus = identifier_checking(error , nb_vector , identifier);
  if (!lstatus) {
    status = false;
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Merging of Vectors objects.
 *
 *  \param[in] error     reference on a StatError object,
 *  \param[in] nb_sample number of Vectors objects,
 *  \param[in] ivec      pointer on the Vectors objects.
 *
 *  \return              Vectors object.
 */
/*--------------------------------------------------------------*/

Vectors* Vectors::merge(StatError &error , int nb_sample , const Vectors **ivec) const

{
  bool status = true;
  int i , j , k , m;
  int inb_vector , cumul_nb_vector , *iidentifier;
  const FrequencyDistribution **phisto;
  Vectors *vec;
  const Vectors **pvec;


  vec = NULL;
  error.init();

  for (i = 0;i < nb_sample;i++) {
    if (ivec[i]->nb_variable != nb_variable) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                    << STAT_error[STATR_NB_VARIABLE];
      error.correction_update((error_message.str()).c_str() , nb_variable);
    }

    else {
      for (j = 0;j < nb_variable;j++) {
        if (ivec[i]->type[j] != type[j]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                        << STAT_label[STATL_VARIABLE] << " " << j + 1 << ": "
                        << STAT_error[STATR_VARIABLE_TYPE];
          error.correction_update((error_message.str()).c_str() , STAT_variable_word[type[j]]);
        }
      }
    }
  }

  if (status) {
    nb_sample++;
    pvec = new const Vectors*[nb_sample];

    pvec[0] = this;
    for (i = 1;i < nb_sample;i++) {
      pvec[i] = ivec[i - 1];
    }

    // computation of the number of individuals

    inb_vector = 0;
    for (i = 0;i < nb_sample;i++) {
      inb_vector += pvec[i]->nb_vector;
    }

    // comparaison of individual identifiers

    iidentifier = new int[inb_vector];

    cumul_nb_vector = 0;
    i = 0;
    for (j = 0;j < nb_sample;j++) {
      for (k = 0;k < pvec[j]->nb_vector;k++) {
        iidentifier[i] = pvec[j]->identifier[k];

        for (m = 0;m < cumul_nb_vector;m++) {
          if (iidentifier[i] == iidentifier[m]) {
            delete [] iidentifier;
            iidentifier = NULL;
            break;
          }
        }

        if (!iidentifier) {
          break;
        }
        i++;
      }

      if (!iidentifier) {
        break;
      }
      cumul_nb_vector += pvec[j]->nb_vector;
    }

    vec = new Vectors(inb_vector , iidentifier , nb_variable , type);
    delete [] iidentifier;

    // copy of vectors

    i = 0;
    for (j = 0;j < nb_sample;j++) {
      for (k = 0;k < pvec[j]->nb_vector;k++) {
        for (m = 0;m < pvec[j]->nb_variable;m++) {
          vec->int_vector[i][m] = pvec[j]->int_vector[k][m];
          vec->real_vector[i][m] = pvec[j]->real_vector[k][m];
        }
        i++;
      }
    }

    phisto = new const FrequencyDistribution*[nb_sample];

    for (i = 0;i < vec->nb_variable;i++) {
      vec->min_value[i] = pvec[0]->min_value[i];
      vec->max_value[i] = pvec[0]->max_value[i];
      for (j = 1;j < nb_sample;j++) {
        if (pvec[j]->min_value[i] < vec->min_value[i]) {
          vec->min_value[i] = pvec[j]->min_value[i];
        }
        if (pvec[j]->max_value[i] > vec->max_value[i]) {
          vec->max_value[i] = pvec[j]->max_value[i];
        }
      }

      vec->min_interval[i] = pvec[0]->min_interval[i];
      for (j = 1;j < nb_sample;j++) {
        if (pvec[j]->min_interval[i] < vec->min_interval[i]) {
          vec->min_interval[i] = pvec[j]->min_interval[i];
        }
      }

      for (j = 0;j < nb_sample;j++) {
        if (pvec[j]->marginal_distribution[i]) {
          phisto[j] = pvec[j]->marginal_distribution[i];
        }
        else {
          break;
        }
      }

      if (j == nb_sample) {
        vec->marginal_distribution[i] = new FrequencyDistribution(nb_sample , phisto);

        vec->mean[i] = vec->marginal_distribution[i]->mean;
        vec->covariance[i][i] = vec->marginal_distribution[i]->variance;
      }

      else {
        vec->build_marginal_histogram(i);

        vec->mean_computation(i);
        vec->variance_computation(i);
      }
    }

    vec->covariance_computation();

    delete [] phisto;
    delete [] pvec;
  }

  return vec;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Merging of Vectors objects.
 *
 *  \param[in] error     reference on a StatError object,
 *  \param[in] nb_sample number of Vectors objects,
 *  \param[in] ivec      pointer on the Vectors objects.
 *
 *  \return              Vectors object.
 */
/*--------------------------------------------------------------*/

Vectors* Vectors::merge(StatError &error , int nb_sample , const vector<Vectors> ivec) const

{
  int i;
  Vectors *vec;
  const Vectors **pvec;


  pvec = new const Vectors*[nb_sample];
  for (i = 0;i < nb_sample;i++) {
    pvec[i] = new Vectors(ivec[i]);
  }

  vec = merge(error , nb_sample , pvec);

  for (i = 0;i < nb_sample;i++) {
    delete pvec[i];
  }
  delete [] pvec;

  return vec;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Shifting of values of a variable.
 *
 *  \param[in] error       reference on a StatError object,
 *  \param[in] variable    variable index,
 *  \param[in] shift_param integer shifting parameter.
 *
 *  \return                Vectors object.
 */
/*--------------------------------------------------------------*/

Vectors* Vectors::shift(StatError &error , int variable , int shift_param) const

{
  bool status = true;
  int i;
  Vectors *vec;


  vec = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if (type[variable] == INT_VALUE) {
      if (shift_param + min_value[variable] < INT_MIN) {
        status = false;
        ostringstream correction_message;
        correction_message << STAT_error[STATR_GREATER_THAN] << " "
                           << INT_MIN - min_value[variable];
        error.correction_update(STAT_error[STATR_SHIFT_VALUE] , (correction_message.str()).c_str());
      }

      if (shift_param + max_value[variable] > INT_MAX) {
        status = false;
        ostringstream correction_message;
        correction_message << STAT_error[STATR_SMALLER_THAN] << " "
                           << INT_MAX - max_value[variable];
        error.correction_update(STAT_error[STATR_SHIFT_VALUE] , (correction_message.str()).c_str());
      }
    }

    else if (type[variable] != REAL_VALUE) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[REAL_VALUE];
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }
  }

  if (status) {
    vec = new Vectors(*this , variable , type[variable]);

    switch (vec->type[variable]) {

    // shifting of integer values

    case INT_VALUE : {
      for (i = 0;i < vec->nb_vector;i++) {
        vec->int_vector[i][variable] = int_vector[i][variable] + shift_param;
      }
      break;
    }

    // shifting of real values

    case REAL_VALUE : {
      for (i = 0;i < vec->nb_vector;i++) {
        vec->real_vector[i][variable] = real_vector[i][variable] + shift_param;
      }
      break;
    }
    }

    vec->min_value[variable] = min_value[variable] + shift_param;
    vec->max_value[variable] = max_value[variable] + shift_param;
    vec->min_interval[variable] = min_interval[variable];

    if ((vec->type[variable] == INT_VALUE) && (vec->min_value[variable] >= 0) &&
        (vec->max_value[variable] <= MARGINAL_DISTRIBUTION_MAX_VALUE)) {
      if (marginal_distribution[variable]) {
        vec->marginal_distribution[variable] = new FrequencyDistribution(*marginal_distribution[variable] ,
                                                                         SHIFT , shift_param);

        vec->mean[variable] = vec->marginal_distribution[variable]->mean;
      }

      else {
        vec->build_marginal_frequency_distribution(variable);
      }
    }

    else {
      vec->build_marginal_histogram(variable);

      vec->mean[variable] = mean[variable] + shift_param;
    }

    for (i = 0;i < vec->nb_variable;i++) {
      vec->covariance[i][variable] = covariance[i][variable];
      vec->covariance[variable][i] = covariance[variable][i];
    }
  }

  return vec;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Shifting of values of a real-valued variable.
 *
 *  \param[in] error       reference on a StatError object,
 *  \param[in] variable    variable index,
 *  \param[in] shift_param real shifting parameter.
 *
 *  \return                Vectors object.
 */
/*--------------------------------------------------------------*/

Vectors* Vectors::shift(StatError &error , int variable , double shift_param) const

{
  bool status = true;
  int i;
  Vectors *vec;


  vec = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if (type[variable] != REAL_VALUE) {
      status = false;
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , STAT_variable_word[REAL_VALUE]);
    }
  }

  if (status) {
    vec = new Vectors(*this , variable , type[variable]);

    // shifting of real values

    for (i = 0;i < vec->nb_vector;i++) {
      vec->real_vector[i][variable] = real_vector[i][variable] + shift_param;
    }

    vec->min_value[variable] = min_value[variable] + shift_param;
    vec->max_value[variable] = max_value[variable] + shift_param;
    vec->min_interval[variable] = min_interval[variable];

    vec->build_marginal_histogram(variable);

    vec->mean[variable] = mean[variable] + shift_param;

    for (i = 0;i < vec->nb_variable;i++) {
      vec->covariance[i][variable] = covariance[i][variable];
      vec->covariance[variable][i] = covariance[variable][i];
    }
  }

  return vec;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Thresholding of values of a variable.
 *
 *  \param[in] error     reference on a StatError object,
 *  \param[in] variable  variable index,
 *  \param[in] threshold integer threshold,
 *  \param[in] mode      mode (ABOVE/BELOW).
 *
 *  \return              Vectors object.
 */
/*--------------------------------------------------------------*/

Vectors* Vectors::thresholding(StatError &error , int variable , int threshold ,
                               threshold_direction mode) const

{
  bool status = true;
  int i;
  Vectors *vec;


  vec = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if ((type[variable] != INT_VALUE) && (type[variable] != REAL_VALUE)) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[REAL_VALUE];
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }

    if (threshold <= min_value[variable]) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_error[STATR_GREATER_THAN] << " " << min_value[variable];
      error.correction_update(STAT_error[STATR_THRESHOLD_VALUE] , (correction_message.str()).c_str());
    }

    if (threshold >= max_value[variable]) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_error[STATR_SMALLER_THAN] << " " << max_value[variable];
      error.correction_update(STAT_error[STATR_THRESHOLD_VALUE] , (correction_message.str()).c_str());
    }
  }

  if (status) {
    vec = new Vectors(*this , variable , type[variable]);

    switch (vec->type[variable]) {

    // thresholding of integer values

    case INT_VALUE : {
      switch (mode) {

      case ABOVE : {
        for (i = 0;i < vec->nb_vector;i++) {
          if (int_vector[i][variable] > threshold) {
            vec->int_vector[i][variable] = threshold;
          }
          else {
            vec->int_vector[i][variable] = int_vector[i][variable];
          }
        }
        break;
      }

      case BELOW : {
        for (i = 0;i < vec->nb_vector;i++) {
          if (int_vector[i][variable] < threshold) {
            vec->int_vector[i][variable] = threshold;
          }
          else {
            vec->int_vector[i][variable] = int_vector[i][variable];
          }
        }
        break;
      }
      }

      break;
    }

    // thresholding of real values

    case REAL_VALUE : {
      switch (mode) {

      case ABOVE : {
        for (i = 0;i < vec->nb_vector;i++) {
          if (real_vector[i][variable] > threshold) {
            vec->real_vector[i][variable] = threshold;
          }
          else {
            vec->real_vector[i][variable] = real_vector[i][variable];
          }
        }
        break;
      }

      case BELOW : {
        for (i = 0;i < vec->nb_vector;i++) {
          if (real_vector[i][variable] < threshold) {
            vec->real_vector[i][variable] = threshold;
          }
          else {
            vec->real_vector[i][variable] = real_vector[i][variable];
          }
        }
        break;
      }
      }

      break;
    }
    }

    switch (mode) {
    case ABOVE :
      vec->min_value[variable] = min_value[variable];
      vec->max_value[variable] = threshold;
      break;
    case BELOW :
      vec->min_value[variable] = threshold;
      vec->max_value[variable] = max_value[variable];
      break;
    }

    if ((vec->type[variable] == INT_VALUE) && (vec->min_value[variable] >= 0) &&
        (vec->max_value[variable] <= MARGINAL_DISTRIBUTION_MAX_VALUE)) {
      vec->build_marginal_frequency_distribution(variable);
    }
    else {
      vec->build_marginal_histogram(variable);
      vec->mean_computation(variable);
    }

    vec->min_interval_computation(variable);

    vec->covariance_computation(variable);
  }

  return vec;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Thresholding of values of a real-valued variable.
 *
 *  \param[in] error     reference on a StatError object,
 *  \param[in] variable  variable index,
 *  \param[in] threshold real threshold,
 *  \param[in] mode      mode (ABOVE/BELOW).
 *
 *  \return              Vectors object.
 */
/*--------------------------------------------------------------*/

Vectors* Vectors::thresholding(StatError &error , int variable , double threshold ,
                               threshold_direction mode) const

{
  bool status = true;
  int i;
  Vectors *vec;


  vec = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if (type[variable] != REAL_VALUE) {
      status = false;
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , STAT_variable_word[REAL_VALUE]);
    }

    if (threshold <= min_value[variable]) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_error[STATR_GREATER_THAN] << " " << min_value[variable];
      error.correction_update(STAT_error[STATR_THRESHOLD_VALUE] , (correction_message.str()).c_str());
    }

    if (threshold >= max_value[variable]) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_error[STATR_SMALLER_THAN] << " " << max_value[variable];
      error.correction_update(STAT_error[STATR_THRESHOLD_VALUE] , (correction_message.str()).c_str());
    }
  }

  if (status) {
    vec = new Vectors(*this , variable , type[variable]);

    // thresholding of real values

    switch (mode) {

    case ABOVE : {
      for (i = 0;i < vec->nb_vector;i++) {
        if (real_vector[i][variable] > threshold) {
          vec->real_vector[i][variable] = threshold;
        }
        else {
          vec->real_vector[i][variable] = real_vector[i][variable];
        }
      }

      vec->min_value[variable] = min_value[variable];
      vec->max_value[variable] = threshold;
      break;
    }

    case BELOW : {
      for (i = 0;i < vec->nb_vector;i++) {
        if (real_vector[i][variable] < threshold) {
          vec->real_vector[i][variable] = threshold;
        }
        else {
          vec->real_vector[i][variable] = real_vector[i][variable];
        }
      }

      vec->min_value[variable] = threshold;
      vec->max_value[variable] = max_value[variable];
      break;
    }
    }

    vec->build_marginal_histogram(variable);

    vec->min_interval_computation(variable);

    vec->mean_computation(variable);
    vec->covariance_computation(variable);
  }

  return vec;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Clustering of values of a variable.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] variable variable index,
 *  \param[in] step     clustering step,
 *  \param[in] mode     mode (FLOOR/ROUND/CEIL).
 *
 *  \return             Vectors object.
 */
/*--------------------------------------------------------------*/

Vectors* Vectors::cluster(StatError &error , int variable , int step , rounding mode) const

{
  bool status = true;
  int i;
  Vectors *vec;


  vec = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if ((type[variable] != INT_VALUE) && (type[variable] != REAL_VALUE)) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[REAL_VALUE];
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }
  }

  if (step < 1) {
    status = false;
    error.update(STAT_error[STATR_CLUSTERING_STEP]);
  }

  if (status) {
    vec = new Vectors(*this , variable , type[variable]);

    switch (vec->type[variable]) {

    // clustering of integer values

    case INT_VALUE : {
      switch (mode) {

      case FLOOR : {
        for (i = 0;i < vec->nb_vector;i++) {
          vec->int_vector[i][variable] = int_vector[i][variable] / step;
        }

        vec->min_value[variable] = (int)min_value[variable] / step;
        vec->max_value[variable] = (int)max_value[variable] / step;
//        vec->min_value[variable] = floor(min_value[variable] / step);
//        vec->max_value[variable] = floor(max_value[variable] / step);
        break;
      }

      case ROUND : {
        for (i = 0;i < vec->nb_vector;i++) {
          vec->int_vector[i][variable] = (int_vector[i][variable] + step / 2) / step;
//          vec->int_vector[i][variable] = (int)::round((double)int_vector[i][variable] / (double)step);
        }

        vec->min_value[variable] = ((int)min_value[variable] + step / 2) / step;
        vec->max_value[variable] = ((int)max_value[variable] + step / 2) / step;
//        vec->min_value[variable] = ::round(min_value[variable] / step);
//        vec->max_value[variable] = ::round(max_value[variable] / step);
        break;
      }

      case CEIL : {
        for (i = 0;i < vec->nb_vector;i++) {
          vec->int_vector[i][variable] = (int_vector[i][variable] + step - 1) / step;
//          vec->int_vector[i][variable] = (int)ceil((double)int_vector[i][variable] / (double)step);
        }

        vec->min_value[variable] = ((int)min_value[variable] + step - 1) / step;
        vec->max_value[variable] = ((int)max_value[variable] + step - 1) / step;
//        vec->min_value[variable] = ceil(min_value[variable] / step);
//        vec->max_value[variable] = ceil(max_value[variable] / step);
        break;
      }
      }

      if (marginal_distribution[variable]) {
        vec->marginal_distribution[variable] = new FrequencyDistribution(*marginal_distribution[variable] ,
                                                                         CLUSTER , step , mode);

        vec->mean[variable] = vec->marginal_distribution[variable]->mean;
        vec->covariance[variable][variable] = vec->marginal_distribution[variable]->variance;
      }

      else {
        vec->build_marginal_frequency_distribution(variable);
      }

      vec->min_interval_computation(variable);
      break;
    }

    // clustering of real values

    case REAL_VALUE : {
      for (i = 0;i < vec->nb_vector;i++) {
        vec->real_vector[i][variable] = real_vector[i][variable] / step;
      }

      vec->min_value[variable] = min_value[variable] / step;
      vec->max_value[variable] = max_value[variable] / step;

      vec->build_marginal_histogram(variable , marginal_histogram[variable]->bin_width / step);

      vec->min_interval_computation(variable);

      vec->mean[variable] = mean[variable] / step;
      vec->covariance[variable][variable] = covariance[variable][variable] / (step * step);
      break;
    }
    }

    vec->covariance_computation(variable);
  }

  return vec;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Transcoding of categories of an integer-valued variable.
 *
 *  \param[in] vec          reference on a Vectors object,
 *  \param[in] variable     variable index,
 *  \param[in] min_category lowest category,
 *  \param[in] max_category highest category,
 *  \param[in] category     transcoding table.
 */
/*--------------------------------------------------------------*/

void Vectors::transcode(const Vectors &vec , int variable , int min_category ,
                        int max_category , int *category)

{
  int i;


  for (i = 0;i < nb_vector;i++) {
    int_vector[i][variable] = category[vec.int_vector[i][variable] -
                                       (int)vec.min_value[variable]] + min_category;
  }

  min_value[variable] = min_category;
  max_value[variable] = max_category;

  build_marginal_frequency_distribution(variable);
  min_interval_computation(variable);

  covariance_computation(variable);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Transcoding of categories of an integer-valued variable.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] variable variable index,
 *  \param[in] category transcoding table.
 *
 *  \return             Vectors object.
 */
/*--------------------------------------------------------------*/

Vectors* Vectors::transcode(StatError &error , int variable , int *category) const

{
  bool status = true , *presence;
  int i;
  int min_category , max_category;
  Vectors *vec;


  vec = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if ((type[variable] != INT_VALUE) && (type[variable] != STATE)) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[STATE];
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }

    else {
      min_category = category[0];
      max_category = category[0];

      for (i = 1;i <= (int)(max_value[variable] - min_value[variable]);i++) {
        if (category[i] < min_category) {
          min_category = category[i];
        }
        if (category[i] > max_category) {
          max_category = category[i];
        }
      }

      if (max_category - min_category == 0) {
        status = false;
        error.update(STAT_error[STATR_NB_CATEGORY]);
      }

      if (max_category - min_category > (int)(max_value[variable] - min_value[variable])) {
        status = false;
        error.update(STAT_error[STATR_NON_CONSECUTIVE_CATEGORIES]);
      }
    }

    if (status) {
      presence = new bool[max_category - min_category + 1];
      for (i = 0;i <= max_category - min_category;i++) {
        presence[i] = false;
      }

      for (i = 0;i <= (int)(max_value[variable] - min_value[variable]);i++) {
        presence[category[i] - min_category] = true;
      }

      for (i = 0;i <= max_category - min_category;i++) {
        if (!presence[i]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_error[STATR_MISSING_CATEGORY] << " " << i + min_category;
          error.update((error_message.str()).c_str());
        }
      }

      delete [] presence;
    }

    if (status) {
      for (i = 0;i <= (int)(max_value[variable] - min_value[variable]);i++) {
        category[i] -= min_category;
      }

      vec = new Vectors(*this , variable , type[variable]);
      vec->transcode(*this , variable , min_category , max_category , category);
    }
  }

  return vec;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Transcoding of categories of an integer-valued variable.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] variable variable index,
 *  \param[in] category transcoding table.
 *
 *  \return             Vectors object.
 */
/*--------------------------------------------------------------*/

Vectors* Vectors::transcode(StatError &error , int variable , vector<int> category) const

{
  return transcode(error , variable , category.data());
}


/*--------------------------------------------------------------*/
/**
 *  \brief Partitioning of values of a variable.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] variable variable index,
 *  \param[in] nb_class number of classes,
 *  \param[in] ilimit   integer limits between classes (beginning of classes).
 *
 *  \return             Vectors object.
 */
/*--------------------------------------------------------------*/

Vectors* Vectors::cluster(StatError &error , int variable ,
                          int nb_class , int *ilimit) const

{
  bool status = true;
  int i , j , k;
  int *int_limit , *category;
  double *real_limit;
  Vectors *vec;


  vec = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if ((nb_class < 2) || (nb_class >= (int)(max_value[variable] - min_value[variable]) + 1)) {
      status = false;
      error.update(STAT_error[STATR_NB_CLASS]);
    }
  }

  if (status) {
    if (type[variable] != REAL_VALUE) {
      int_limit = new int[nb_class + 1];
      int_limit[0] = (int)min_value[variable];
      for (i = 1;i < nb_class;i++) {
        int_limit[i] = ilimit[i - 1];
      }
      int_limit[nb_class] = (int)max_value[variable] + 1;

      for (i = 0;i < nb_class;i++) {
        if (int_limit[i] >= int_limit[i + 1]) {
          status = false;
          error.update(STAT_error[STATR_CLUSTER_LIMIT]);
        }
      }

      if (status) {
        category = new int[(int)(max_value[variable] - min_value[variable]) + 1];

        i = 0;
        for (j = 0;j < nb_class;j++) {
          for (k = int_limit[j];k < int_limit[j + 1];k++) {
            category[i++] = j;
          }
        }

        vec = new Vectors(*this , variable , type[variable]);
        vec->transcode(*this , variable , 0 , nb_class - 1 , category);

        delete [] category;
      }

      delete [] int_limit;
    }

    else {
      real_limit = new double[nb_class + 1];
      real_limit[0] = min_value[variable];
      for (i = 1;i < nb_class;i++) {
        real_limit[i] = ilimit[i - 1];
      }
      real_limit[nb_class] = max_value[variable] + 1;

      for (i = 0;i < nb_class;i++) {
        if (real_limit[i] >= real_limit[i + 1]) {
          status = false;
          error.update(STAT_error[STATR_CLUSTER_LIMIT]);
        }
      }

      if (status) {
        vec = new Vectors(*this , variable , INT_VALUE);
        vec->cluster(*this , variable , nb_class , real_limit);
      }

      delete [] real_limit;
    }
  }

  return vec;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Partitioning of values of a variable.
 *
 *  \param[in] error     reference on a StatError object,
 *  \param[in] variable  variable index,
 *  \param[in] nb_class  number of classes,
 *  \param[in] ilimit    integer limits between classes (beginning of classes).
 *
 *  \return              Vectors object.
 */
/*--------------------------------------------------------------*/

Vectors* Vectors::cluster(StatError &error , int variable ,
                          int nb_class , vector<int> ilimit) const

{
  return cluster(error , variable , nb_class , ilimit.data());
}


/*--------------------------------------------------------------*/
/**
 *  \brief Partitioning of values of a real-valued variable.
 *
 *  \param[in] vec      reference on a Vectors object,
 *  \param[in] variable variable index,
 *  \param[in] nb_class number of classes,
 *  \param[in] limit    real limits between classes (beginning of classes).
 */
/*--------------------------------------------------------------*/

void Vectors::cluster(const Vectors &vec , int variable , int nb_class , double *limit)

{
  int i , j;


  for (i = 0;i < nb_vector;i++) {
    for (j = 0;j < nb_class;j++) {
      if (vec.real_vector[i][variable] < limit[j + 1]) {
        int_vector[i][variable] = j;
        break;
      }
    }
  }

  min_value_computation(variable);
  max_value_computation(variable);

  build_marginal_frequency_distribution(variable);
  min_interval_computation(variable);

  covariance_computation(variable);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Partitioning of values of a real-valued variable.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] variable variable index,
 *  \param[in] nb_class number of classes,
 *  \param[in] ilimit   real limits between classes (beginning of classes).
 *
 *  \return             Vectors object.
 */
/*--------------------------------------------------------------*/

Vectors* Vectors::cluster(StatError &error , int variable ,
                          int nb_class , double *ilimit) const

{
  bool status = true;
  int i;
  double *limit;
  Vectors *vec;


  vec = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if (type[variable] != REAL_VALUE) {
      status = false;
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , STAT_variable_word[REAL_VALUE]);
    }
  }

  if (nb_class < 2) {
    status = false;
    error.update(STAT_error[STATR_NB_CLASS]);
  }

  if (status) {
    limit = new double[nb_class + 1];
    limit[0] = min_value[variable];
    for (i = 1;i < nb_class;i++) {
      limit[i] = ilimit[i - 1];
    }
    limit[nb_class] = max_value[variable] + DOUBLE_ERROR;

    for (i = 0;i < nb_class;i++) {
      if (limit[i] >= limit[i + 1]) {
        status = false;
        error.update(STAT_error[STATR_CLUSTER_LIMIT]);
      }
    }

    if (status) {
      vec = new Vectors(*this , variable , INT_VALUE);
      vec->cluster(*this , variable , nb_class , limit);
    }

    delete [] limit;
  }

  return vec;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Partitioning of values of a real-valued variable.
 *
 *  \param[in] error     reference on a StatError object,
 *  \param[in] variable  variable index,
 *  \param[in] nb_class  number of classes,
 *  \param[in] ilimit    real limits between classes (beginning of classes).
 *
 *  \return              Vectors object.
 */
/*--------------------------------------------------------------*/

Vectors* Vectors::cluster(StatError &error , int variable ,
                          int nb_class , vector<double> ilimit) const

{
  return cluster(error , variable , nb_class , ilimit.data());
}


/*--------------------------------------------------------------*/
/**
 *  \brief Scaling of a variable.
 *
 *  \param[in] error         reference on a StatError object,
 *  \param[in] variable      variable index,
 *  \param[in] scaling_coeff integer scaling factor.
 *
 *  \return                  Vectors object.
 */
/*--------------------------------------------------------------*/

Vectors* Vectors::scaling(StatError &error , int variable , int scaling_coeff) const

{
  bool status = true;
  int i;
  Vectors *vec;


  vec = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if (type[variable] == INT_VALUE) {
      if ((min_value[variable] * scaling_coeff < INT_MIN) ||
          (max_value[variable] * scaling_coeff > INT_MAX)) {
        status = false;
        error.update(STAT_error[STATR_SCALING_COEFF]);
      }
    }

    else if (type[variable] != REAL_VALUE) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[REAL_VALUE];
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }
  }

  if (scaling_coeff <= 1) {
    status = false;
    error.update(STAT_error[STATR_SCALING_COEFF]);
  }

  if (status) {
    vec = new Vectors(*this , variable , type[variable]);

    switch (vec->type[variable]) {

    // scaling of integer values

    case INT_VALUE : {
      for (i = 0;i < vec->nb_vector;i++) {
        vec->int_vector[i][variable] = int_vector[i][variable] * scaling_coeff;
      }
      break;
    }

    // scaling of real values

    case REAL_VALUE : {
      for (i = 0;i < vec->nb_vector;i++) {
        vec->real_vector[i][variable] = real_vector[i][variable] * scaling_coeff;
      }
      break;
    }
    }

    vec->min_value[variable] = min_value[variable] * scaling_coeff;
    vec->max_value[variable] = max_value[variable] * scaling_coeff;
    vec->min_interval[variable] = min_interval[variable] * scaling_coeff;

    vec->build_marginal_frequency_distribution(variable);

    vec->covariance_computation(variable);
  }

  return vec;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Scaling of a variable.
 *
 *  \param[in] error         reference on a StatError object,
 *  \param[in] variable      variable index,
 *  \param[in] scaling_coeff real scaling factor.
 *
 *  \return                  Vectors object.
 */
/*--------------------------------------------------------------*/

Vectors* Vectors::scaling(StatError &error , int variable , double scaling_coeff) const

{
  bool status = true;
  int i;
  Vectors *vec;


  vec = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if ((type[variable] != INT_VALUE) && (type[variable] != REAL_VALUE)) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[REAL_VALUE];
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }
  }

  if (scaling_coeff <= 0) {
    status = false;
    error.update(STAT_error[STATR_SCALING_COEFF]);
  }

  if (status) {
    vec = new Vectors(*this , variable , REAL_VALUE);

    switch (type[variable]) {

    // scaling of integer values

    case INT_VALUE : {
      for (i = 0;i < vec->nb_vector;i++) {
        vec->real_vector[i][variable] = int_vector[i][variable] * scaling_coeff;
      }
      break;
    }

    // scaling of real values

    case REAL_VALUE : {
      for (i = 0;i < vec->nb_vector;i++) {
        vec->real_vector[i][variable] = real_vector[i][variable] * scaling_coeff;
      }
      break;
    }
    }

    vec->min_value[variable] = min_value[variable] * scaling_coeff;
    vec->max_value[variable] = max_value[variable] * scaling_coeff;
    vec->min_interval[variable] = min_interval[variable] * scaling_coeff;

    vec->build_marginal_histogram(variable);

    vec->mean[variable] = mean[variable] * scaling_coeff;
    vec->covariance_computation(variable);
  }

  return vec;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Rounding of values of a real-valued variable.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] variable variable index,
 *  \param[in] mode     mode (FLOOR/ROUND/CEIL).
 *
 *  \return             Vectors object.
 */
/*--------------------------------------------------------------*/

Vectors* Vectors::round(StatError &error , int variable , rounding mode) const

{
  bool status = true;
  int i , j;
  variable_nature *itype;
  Vectors *vec;


  vec = NULL;
  error.init();

  if (variable != I_DEFAULT) {
    if ((variable < 1) || (variable > nb_variable)) {
      status = false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
    }

    else {
      variable--;

      if (type[variable] != REAL_VALUE) {
        status = false;
        error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , STAT_variable_word[REAL_VALUE]);
      }
    }
  }

  if (status) {
    for (i = 0;i < nb_variable;i++) {
      if (((variable == I_DEFAULT) && (type[i] == REAL_VALUE)) || (variable == i)) {
        if (min_value[i] < INT_MIN) {
          status = false;
          ostringstream error_message , correction_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                        << STAT_error[STATR_ROUNDED_VALUE];
          correction_message << STAT_error[STATR_GREATER_THAN] << " " << INT_MIN;
          error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
        }

        if (max_value[i] > INT_MAX) {
          status = false;
          ostringstream error_message , correction_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                        << STAT_error[STATR_ROUNDED_VALUE];
          correction_message << STAT_error[STATR_SMALLER_THAN] << " " << INT_MAX;
          error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
        }
      }
    }
  }

  if (status) {
    itype = new variable_nature[nb_variable];

    if (variable == I_DEFAULT) {
      for (i = 0;i < nb_variable;i++) {
        itype[i] = INT_VALUE;
      }
    }

    else {
      for (i = 0;i < nb_variable;i++) {
        itype[i] = type[i];
      }
      itype[variable] = INT_VALUE;
    }

    vec = new Vectors(nb_vector , identifier , nb_variable , itype);
    delete [] itype;

    for (i = 0;i < vec->nb_vector;i++) {
      for (j = 0;j < vec->nb_variable;j++) {

        // copy of integer values

        if (type[j] != REAL_VALUE) {
          vec->int_vector[i][j] = int_vector[i][j];
        }

        else {

          // rounding of real values

          if ((variable == I_DEFAULT) || (variable == j)) {
            switch (mode) {
            case FLOOR :
              vec->int_vector[i][j] = (int)floor(real_vector[i][j]);
              break;
            case ROUND :
              vec->int_vector[i][j] = (int)::round(real_vector[i][j]);
              break;
            case CEIL :
              vec->int_vector[i][j] = (int)ceil(real_vector[i][j]);
              break;
            }
          }

          // copy of real values

          else {
            vec->real_vector[i][j] = real_vector[i][j];
          }
        }
      }
    }

    for (i = 0;i < vec->nb_variable;i++) {
      if ((variable == I_DEFAULT) || (variable == i)) {
        switch (mode) {
        case FLOOR :
          vec->min_value[i] = floor(min_value[i]);
          vec->max_value[i] = floor(max_value[i]);
          break;
        case ROUND :
          vec->min_value[i] = ::round(min_value[i]);
          vec->max_value[i] = ::round(max_value[i]);
          break;
        case CEIL :
          vec->min_value[i] = ceil(min_value[i]);
          vec->max_value[i] = ceil(max_value[i]);
          break;
        }

        vec->build_marginal_frequency_distribution(i);
        vec->min_interval_computation(i);
      }

      else {
        vec->min_value[i] = min_value[i];
        vec->max_value[i] = max_value[i];
        vec->min_interval[i] = min_interval[i];

        if (marginal_distribution[i]) {
          vec->marginal_distribution[i] = new FrequencyDistribution(*marginal_distribution[i]);
        }
        if (marginal_histogram[i]) {
          vec->marginal_histogram[i] = new Histogram(*marginal_histogram[i]);
        }

        vec->mean[i] = mean[i];
        for (j = 0;j < vec->nb_variable;j++) {
          vec->covariance[i][j] = covariance[i][j];
        }
      }
    }

    vec->covariance_computation(variable);
  }

  return vec;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Log-transform of values.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] variable variable index,
 *  \param[in] base     base of the logarithm.
 *
 *  \return             Vectors object.
 */
/*--------------------------------------------------------------*/

Vectors* Vectors::log_transform(StatError &error , int variable , log_base base) const

{
  bool status = true;
  int i , j , k;
  int offset , inb_variable;
  variable_nature *itype;
  Vectors *vec;


  vec = NULL;
  error.init();

  if (type[0] == STATE) {
    offset = 1;
    if (nb_variable == 1) {
      status = false;
      error.update(STAT_error[STATR_NB_VARIABLE]);
    }
  }
  else {
    offset = 0;
  }

  if (variable != I_DEFAULT) {
    if ((variable < offset + 1) || (variable > nb_variable)) {
      status = false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
    }
    else {
      variable--;
    }
  }

  for (i = offset;i < nb_variable;i++) {
    if ((variable == I_DEFAULT) || (variable == i)) {
      if ((type[i] != INT_VALUE) && (type[i] != REAL_VALUE)) {
        status = false;
        ostringstream error_message , correction_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                      << STAT_error[STATR_VARIABLE_TYPE];
        correction_message << STAT_variable_word[INT_VALUE] << " or "
                           << STAT_variable_word[REAL_VALUE];
        error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
      }

      else if (min_value[i] <= 0.) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                      << STAT_error[STATR_POSITIVE_MIN_VALUE];
        error.update((error_message.str()).c_str());
      }
    }
  }

  if (status) {
    if (variable == I_DEFAULT) {
      inb_variable = nb_variable;
    }
    else {
      inb_variable = offset + 1;
    }

    itype = new variable_nature[inb_variable];
    if (type[0] == STATE) {
      itype[0] = type[0];
    }
    for (i = offset;i < inb_variable;i++) {
      itype[i] = REAL_VALUE;
    }

    vec = new Vectors(nb_vector , identifier , inb_variable , itype);
    delete [] itype;

    // copy of the state variable

    if (type[0] == STATE) {
      for (i = 0;i < vec->nb_vector;i++) {
        vec->int_vector[i][0] = int_vector[i][0];
      }

      vec->min_value[0] = min_value[0];
      vec->max_value[0] = max_value[0];

      vec->marginal_distribution[0] = new FrequencyDistribution(*marginal_distribution[0]);
    }

    // log-transform of values

    switch (base) {

    case NATURAL : {
      for (i = 0;i < nb_vector;i++) {
        j = offset;
        for (k = offset;k < nb_variable;k++) {
          if ((variable == I_DEFAULT) || (variable == k)) {
            if (type[k] != REAL_VALUE) {
              vec->real_vector[i][j] = log(int_vector[i][k]);
            }

            else {
              vec->real_vector[i][j] = log(real_vector[i][k]);
            }

            j++;
          }
        }
      }
      break;
    }

    case TWO : {
      for (i = 0;i < nb_vector;i++) {
        j = offset;
        for (k = offset;k < nb_variable;k++) {
          if ((variable == I_DEFAULT) || (variable == k)) {
            if (type[k] != REAL_VALUE) {
              vec->real_vector[i][j] = log2(int_vector[i][k]);
            }

            else {
              vec->real_vector[i][j] = log2(real_vector[i][k]);
            }

            j++;
          }
        }
      }
      break;
    }

    case TEN : {
       for (i = 0;i < nb_vector;i++) {
        j = offset;
        for (k = offset;k < nb_variable;k++) {
          if ((variable == I_DEFAULT) || (variable == k)) {
            if (type[k] != REAL_VALUE) {
              vec->real_vector[i][j] = log10(int_vector[i][k]);
            }

            else {
              vec->real_vector[i][j] = log10(real_vector[i][k]);
            }

            j++;
          }
        }
      }
     break;
    }
    }

    for (i = offset;i < vec->nb_variable;i++) {
      vec->min_value_computation(i);
      vec->max_value_computation(i);

      vec->build_marginal_histogram(i);
      vec->min_interval_computation(i);

      vec->mean_computation(i);
      vec->variance_computation(i);
    }

    vec->covariance_computation();
 }

  return vec;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Selection of individuals taking values in a given range for a variable.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] display    flag for displaying the selected individuals,
 *  \param[in] variable   variable index,
 *  \param[in] imin_value lowest integer value,
 *  \param[in] imax_value highest integer value,
 *  \param[in] keep       flag for keeping or rejecting the selected individuals.
 *
 *  \return               Vectors object.
 */
/*--------------------------------------------------------------*/

Vectors* Vectors::value_select(StatError &error , bool display , int variable ,
                               int imin_value , int imax_value , bool keep) const

{
  bool status = true;
  int i;
  int inb_vector , *index , *iidentifier;
  Vectors *vec;


  vec = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if ((imin_value > max_value[variable]) || (imin_value > imax_value)) {
      status = false;
      error.update(STAT_error[STATR_MIN_VALUE]);
    }
    if ((imax_value < min_value[variable]) || (imax_value < imin_value)) {
      status = false;
      error.update(STAT_error[STATR_MAX_VALUE]);
    }
  }

  if (status) {

    // selection of individuals

    iidentifier = new int[nb_vector];
    index = new int[nb_vector];
    inb_vector = 0;

    if (type[variable] != REAL_VALUE) {
      for (i = 0;i < nb_vector;i++) {
        if ((int_vector[i][variable] >= imin_value) && (int_vector[i][variable] <= imax_value)) {
          if (keep) {
            iidentifier[inb_vector] = identifier[i];
            index[inb_vector++] = i;
          }
        }

        else if (!keep) {
          iidentifier[inb_vector] = identifier[i];
          index[inb_vector++] = i;
        }
      }
    }

    else {
      for (i = 0;i < nb_vector;i++) {
        if ((real_vector[i][variable] >= imin_value) && (real_vector[i][variable] <= imax_value)) {
          if (keep) {
            iidentifier[inb_vector] = identifier[i];
            index[inb_vector++] = i;
          }
        }

        else if (!keep) {
          iidentifier[inb_vector] = identifier[i];
          index[inb_vector++] = i;
        }
      }
    }

    if (inb_vector == 0) {
      status = false;
      error.update(STAT_error[STATR_EMPTY_SAMPLE]);
    }

    // copy of vectors

    if (status) {
      if ((display) && (inb_vector <= DISPLAY_NB_INDIVIDUAL)) {
        cout << "\n" << STAT_label[inb_vector == 1 ? STATL_VECTOR : STATL_VECTORS] << ": ";
        for (i = 0;i < inb_vector;i++) {
          cout << iidentifier[i] << ", ";
        }
        cout << endl;
      }

      vec = new Vectors(*this , inb_vector , index);
    }

    delete [] iidentifier;
    delete [] index;
  }

  return vec;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Selection of individuals taking values in a given range for a real-valued variable.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] display    flag for displaying the selected individuals,
 *  \param[in] variable   variable index,
 *  \param[in] imin_value lowest real value,
 *  \param[in] imax_value highest real value,
 *  \param[in] keep       flag for keeping or rejecting the selected individuals.
 *
 *  \return               Vectors object.
 */
/*--------------------------------------------------------------*/

Vectors* Vectors::value_select(StatError &error , bool display , int variable ,
                               double imin_value , double imax_value , bool keep) const

{
  bool status = true;
  int i;
  int inb_vector , *index , *iidentifier;
  Vectors *vec;


  vec = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if (type[variable] != REAL_VALUE) {
      status = false;
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , STAT_variable_word[REAL_VALUE]);
    }
    
    if ((imin_value > max_value[variable]) || (imin_value > imax_value)) {
      status = false;
      error.update(STAT_error[STATR_MIN_VALUE]);
    }
    if ((imax_value < min_value[variable]) || (imax_value < imin_value)) {
      status = false;
      error.update(STAT_error[STATR_MAX_VALUE]);
    }
  }

  if (status) {

    // selection of individuals

    iidentifier = new int[nb_vector];
    index = new int[nb_vector];
    inb_vector = 0;

    for (i = 0;i < nb_vector;i++) {
      if ((real_vector[i][variable] >= imin_value) && (real_vector[i][variable] <= imax_value)) {
        if (keep) {
          iidentifier[inb_vector] = identifier[i];
          index[inb_vector++] = i;
        }
      }

      else if (!keep) {
        iidentifier[inb_vector] = identifier[i];
        index[inb_vector++] = i;
      }
    }

    if (inb_vector == 0) {
      status = false;
      error.update(STAT_error[STATR_EMPTY_SAMPLE]);
    }

    // copy of vectors

    if (status) {
      if ((display) && (inb_vector <= DISPLAY_NB_INDIVIDUAL)) {
        cout << "\n" << STAT_label[inb_vector == 1 ? STATL_VECTOR : STATL_VECTORS] << ": ";
        for (i = 0;i < inb_vector;i++) {
          cout << iidentifier[i] << ", ";
        }
        cout << endl;
      }

      vec = new Vectors(*this , inb_vector , index);
    }

    delete [] iidentifier;
    delete [] index;
  }

  return vec;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Checking of an identifier list.
 *
 *  \param[in] error                  reference on a StatError object,
 *  \param[in] nb_individual          number of individuals,
 *  \param[in] identifier             individual identifiers,
 *  \param[in] nb_selected_individual number of selected individuals,
 *  \param[in] selected_identifier    selected individual identifiers,
 *  \param[in] data_label             data type label.
 *
 *  \return                           error status.
 */
/*--------------------------------------------------------------*/

bool selected_identifier_checking(StatError &error , int nb_individual , int *identifier ,
                                  int nb_selected_individual , int *selected_identifier ,
                                  const char *data_label)

{
  bool status = true , *selected_individual;
  int i , j;
  int max_identifier;


  max_identifier = selected_identifier[0];
  for (i = 1;i < nb_selected_individual;i++) {
    if (selected_identifier[i] > max_identifier) {
      max_identifier = selected_identifier[i];
    }
  }

  selected_individual = new bool[max_identifier + 1];
  for (i = 0;i <= max_identifier;i++) {
    selected_individual[i] = false;
  }

  for (i = 0;i < nb_selected_individual;i++) {
    for (j = 0;j < nb_individual;j++) {
      if (selected_identifier[i] == identifier[j]) {
        break;
      }
    }

    if (j == nb_individual) {
      status = false;
      ostringstream error_message;
      error_message << selected_identifier[i] << ": " << STAT_error[STATR_BAD] << " "
                    << data_label << " " << STAT_label[STATL_IDENTIFIER];
      error.update((error_message.str()).c_str());
    }

    else if (selected_individual[selected_identifier[i]]) {
      status = false;
      ostringstream error_message;
      error_message << data_label << " " << selected_identifier[i] << " "
                    << STAT_error[STATR_ALREADY_SELECTED];
      error.update((error_message.str()).c_str());
    }
    else {
      selected_individual[selected_identifier[i]] = true;
    }
  }

  delete [] selected_individual;

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Selection of individuals by their identifiers.
 *
 *  \param[in] nb_individual          number of individuals,
 *  \param[in] identifier             individual identifiers,
 *  \param[in] nb_selected_individual number of selected individuals,
 *  \param[in] selected_identifier    selected individual identifiers,
 *  \param[in] keep                   flag for keeping or rejecting the selected individuals.
 *
 *  \return                           selected individual indices.
 */
/*--------------------------------------------------------------*/

int* identifier_select(int nb_individual , int *identifier , int nb_selected_individual ,
                       int *selected_identifier , bool keep)

{
  int i , j , k;
  int *index;


  if (keep) {
    index = new int[nb_selected_individual];

    for (i = 0;i < nb_selected_individual;i++) {
      for (j = 0;j < nb_individual;j++) {
        if (selected_identifier[i] == identifier[j]) {
          index[i] = j;
          break;
        }
      }
    }
  }

  else {
    index = new int[nb_individual - nb_selected_individual];

    i = 0;
    for (j = 0;j < nb_individual;j++) {
      for (k = 0;k < nb_selected_individual;k++) {
        if (selected_identifier[k] == identifier[j]) {
          break;
        }
      }

      if (k == nb_selected_individual) {
        index[i++] = j;
      }
    }
  }

  return index;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Selection of individuals by their identifiers.
 *
 *  \param[in] error       reference on a StatError object,
 *  \param[in] inb_vector  number of individuals,
 *  \param[in] iidentifier identifiers,
 *  \param[in] keep        flag for keeping or rejecting the selected individuals.
 *
 *  \return                Vectors object.
 */
/*--------------------------------------------------------------*/

Vectors* Vectors::select_individual(StatError &error , int inb_vector ,
                                    int *iidentifier , bool keep) const

{
  bool status = true;
  int *index;
  Vectors *vec;


  vec = NULL;
  error.init();

  if ((inb_vector < 1) || (inb_vector > (keep ? nb_vector : nb_vector - 1))) {
    status = false;
    error.update(STAT_error[STATR_NB_VECTOR]);
  }

  else {
    status = selected_identifier_checking(error , nb_vector , identifier , inb_vector ,
                                          iidentifier , STAT_label[STATL_VECTOR]);
  }

  if (status) {
    index = identifier_select(nb_vector , identifier , inb_vector , iidentifier , keep);

    vec = new Vectors(*this , (keep ? inb_vector : nb_vector - inb_vector) , index);

    delete [] index;
  }

  return vec;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Selection of individuals by their identifiers.
 *
 *  \param[in] error          reference on a StatError object,
 *  \param[in] inb_vector     number of individuals,
 *  \param[in] vec_identifier identifiers,
 *  \param[in] keep           flag for keeping or rejecting the selected individuals.
 *
 *  \return                   Vectors object.
 */
/*--------------------------------------------------------------*/

Vectors* Vectors::select_individual(StatError &error , int inb_vector ,
                                    vector<int> iidentifier , bool keep) const

{
  return select_individual(error , inb_vector , iidentifier.data() , keep);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Selection of variables.
 *
 *  \param[in] vec      reference on a Vectors object,
 *  \param[in] variable selected variable indices.
 */
/*--------------------------------------------------------------*/

void Vectors::select_variable(const Vectors &vec , int *variable)

{
  int i , j;


  for (i = 0;i < nb_vector;i++) {
    for (j = 0;j < nb_variable;j++) {
      int_vector[i][j] = vec.int_vector[i][variable[j]];
      real_vector[i][j] = vec.real_vector[i][variable[j]];
    }
  }

  for (i = 0;i < nb_variable;i++) {
    min_value[i] = vec.min_value[variable[i]];
    max_value[i] = vec.max_value[variable[i]];
    min_interval[i] = vec.min_interval[variable[i]];

    if (vec.marginal_distribution[variable[i]]) {
      marginal_distribution[i] = new FrequencyDistribution(*(vec.marginal_distribution[variable[i]]));
    }
    if (vec.marginal_histogram[variable[i]]) {
      marginal_histogram[i] = new Histogram(*(vec.marginal_histogram[variable[i]]));
    }

    mean[i] = vec.mean[variable[i]];
    for (j = 0;j < nb_variable;j++) {
      covariance[i][j] = vec.covariance[variable[i]][variable[j]];
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Selection of 2 variables supposed to be an explanatory variable and a response variable.
 *
 *  \param[in] explanatory_variable explanatory variable,
 *  \param[in] response_variable    response variable.
 *
 *  \return                         Vectors object.
 */
/*--------------------------------------------------------------*/

Vectors* Vectors::select_variable(int explanatory_variable , int response_variable) const

{
  int variable[2];
  variable_nature itype[2];
  Vectors *vec;


  variable[0] = explanatory_variable;
  variable[1] = response_variable;

  itype[0] = type[explanatory_variable];
  itype[1] = type[response_variable];

  vec = new Vectors(nb_vector , identifier , 2 , itype);

  vec->select_variable(*this , variable);

  return vec;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Selection of variables.
 *
 *  \param[in] nb_variable          number of variables,
 *  \param[in] nb_selected_variable number of selected variables,
 *  \param[in] selected_variable    selected variable indices,
 *  \param[in] keep                 flag for keeping or rejecting the selected variables.
 *
 *  \return    variable             variable indices.
 */
/*--------------------------------------------------------------*/

int* select_variable(int nb_variable , int nb_selected_variable ,
                     int *selected_variable , bool keep)

{
  int i , j , k;
  int *variable;


  for (i = 0;i < nb_selected_variable;i++) {
    selected_variable[i]--;
  }

  if (keep) {
    variable = new int[nb_selected_variable];

    for (i = 0;i < nb_selected_variable;i++) {
      variable[i] = selected_variable[i];
    }
  }

  else {
    variable = new int[nb_variable - nb_selected_variable];

    i = 0;
    for (j = 0;j < nb_variable;j++) {
      for (k = 0;k < nb_selected_variable;k++) {
        if (selected_variable[k] == j) {
          break;
        }
      }

      if (k == nb_selected_variable) {
        variable[i++] = j;
      }
    }
  }

  return variable;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Selection of variables.
 *
 *  \param[in] error        reference on a StatError object,
 *  \param[in] inb_variable number of variables,
 *  \param[in] ivariable    variable indices,
 *  \param[in] keep         flag for keeping or rejecting the selected variables.
 *
 *  \return                 Vectors object.
 */
/*--------------------------------------------------------------*/

Vectors* Vectors::select_variable(StatError &error , int inb_variable ,
                                  int *ivariable , bool keep) const

{
  bool status = true , *selected_variable;
  int i;
  int bnb_variable , *variable;
  variable_nature *itype;
  Vectors *vec;


  vec = NULL;
  error.init();

  if ((inb_variable < 1) || (inb_variable > (keep ? nb_variable : nb_variable - 1))) {
    status = false;
    error.update(STAT_error[STATR_NB_SELECTED_VARIABLE]);
  }

  else {
    selected_variable = new bool[nb_variable + 1];
    for (i = 1;i <= nb_variable;i++) {
      selected_variable[i] = false;
    }

    for (i = 0;i < inb_variable;i++) {
      if ((ivariable[i] < 1) || (ivariable[i] > nb_variable)) {
        status = false;
        ostringstream error_message;
        error_message << ivariable[i] << ": " << STAT_error[STATR_VARIABLE_INDEX];
        error.update((error_message.str()).c_str());
      }

      else if (selected_variable[ivariable[i]]) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << ivariable[i] << " "
                      << STAT_error[STATR_ALREADY_SELECTED];
        error.update((error_message.str()).c_str());
      }
      else {
        selected_variable[ivariable[i]] = true;
      }
    }

    delete [] selected_variable;
  }

  if (status) {
    variable = stat_tool::select_variable(nb_variable , inb_variable , ivariable , keep);

    bnb_variable = (keep ? inb_variable : nb_variable - inb_variable);

    itype = new variable_nature[bnb_variable];
    for (i = 0;i < bnb_variable;i++) {
      itype[i] = type[variable[i]];
    }

    vec = new Vectors(nb_vector , identifier , bnb_variable , itype);
    vec->select_variable(*this , variable);

    delete [] variable;
    delete [] itype;
  }

  return vec;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Selection of variables.
 *
 *  \param[in] error        reference on a StatError object,
 *  \param[in] inb_variable number of variables,
 *  \param[in] ivariable    variable indices,
 *  \param[in] keep         flag for keeping or rejecting the selected variables.
 *
 *  \return                 Vectors object.
 */
/*--------------------------------------------------------------*/

Vectors* Vectors::select_variable(StatError &error , int inb_variable ,
                                  vector<int>ivariable , bool keep) const

{
  return select_variable(error , inb_variable , ivariable.data() , keep);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Removing of the 1st variable (supposed to be a state variable).
 *
 *  \return Vectors object.
 */
/*--------------------------------------------------------------*/

Vectors* Vectors::remove_variable_1() const

{
  int i;
  int *variable;
  variable_nature *itype;
  Vectors *vec;


  variable = new int[nb_variable - 1];
  itype = new variable_nature[nb_variable - 1];
  for (i = 0;i < nb_variable - 1;i++) {
    variable[i] = i + 1;
    itype[i] = type[i + 1];
  }

  vec = new Vectors(nb_vector , identifier , nb_variable - 1 , itype);
  vec->select_variable(*this , variable);

  delete [] variable;
  delete [] itype;

  return vec;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Summation of variables.
 *
 *  \param[in] error              reference on a StatError object,
 *  \param[in] nb_summed_variable number of variables to be summed,
 *  \param[in] variable           variable indices.
 *
 *  \return                       Vectors object.
 */
/*--------------------------------------------------------------*/

Vectors* Vectors::sum_variable(StatError &error , int nb_summed_variable , int *ivariable) const

{
  bool status = true , *selected_variable;
  int i , j , k , m;
  int inb_variable , *copied_variable , *summed_variable;
  variable_nature *itype;
  Vectors *vec;


  vec = NULL;
  error.init();

  if (nb_variable == 1) {
    status = false;
    error.update(STAT_error[STATR_NB_VARIABLE]);
  }

  if ((nb_summed_variable < 2) || (nb_summed_variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_NB_SUMMED_VARIABLE]);
  }

  else {
    selected_variable = new bool[nb_variable + 1];
    for (i = 1;i <= nb_variable;i++) {
      selected_variable[i] = false;
    }

    for (i = 0;i < nb_summed_variable;i++) {
      if ((ivariable[i] < 1) || (ivariable[i] > nb_variable)) {
        status = false;
        ostringstream error_message;
        error_message << ivariable[i] << ": " << STAT_error[STATR_VARIABLE_INDEX];
        error.update((error_message.str()).c_str());
      }

      else {
        if (selected_variable[ivariable[i]]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << ivariable[i] << " "
                        << STAT_error[STATR_ALREADY_SELECTED];
          error.update((error_message.str()).c_str());
        }

        else {
          selected_variable[ivariable[i]] = true;

          if ((type[ivariable[i] - 1] != INT_VALUE) && (type[ivariable[i] - 1] != REAL_VALUE)) {
            status = false;
            ostringstream error_message , correction_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << ivariable[i] << ": "
                          << STAT_error[STATR_VARIABLE_TYPE];
            correction_message << STAT_variable_word[INT_VALUE] << " or "
                               << STAT_variable_word[REAL_VALUE];
            error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
          }
        }
      }
    }

    delete [] selected_variable;
  }

  if (status) {
    inb_variable = nb_variable - nb_summed_variable + 1;
    for (i = 0;i < nb_summed_variable;i++) {
      ivariable[i]--;
    }

    summed_variable = new int[nb_summed_variable];
    copied_variable = new int[nb_variable - nb_summed_variable];
    i = 0;
    j = 0;
    for (k = 0;k < nb_variable;k++) {
      for (m = 0;m < nb_summed_variable;m++) {
        if (ivariable[m] == k) {
          summed_variable[i++] = k;
          break;
        }
      }
      if (m == nb_summed_variable) {
        copied_variable[j++] = k;
      }
    }

    itype = new variable_nature[inb_variable];
    i = 0;
    for (j = 0;j < inb_variable;j++) {
      if (j == summed_variable[0]) {
        itype[j] = type[j];
        for (k = 1;k < nb_summed_variable;k++) {
          if (type[summed_variable[k]] == REAL_VALUE) {
           itype[j] = REAL_VALUE;
          }
        }
      }
      else {
        itype[j] = type[copied_variable[i++]];
      }
    }

    vec = new Vectors(nb_vector , identifier , inb_variable , itype);
    delete [] itype;

    for (i = 0;i < nb_vector;i++) {
      j = 0;
      for (k = 0;k < inb_variable;k++) {

        // summation of values

        if (k == summed_variable[0]) {
          switch (vec->type[k]) {

          case INT_VALUE : {
            vec->int_vector[i][k] = 0.;

            for (m = 0;m < nb_summed_variable;m++) {
              switch (type[summed_variable[m]]) {
              case INT_VALUE :
                vec->int_vector[i][k] += int_vector[i][summed_variable[m]];
                break;
              case REAL_VALUE :
                vec->int_vector[i][k] += real_vector[i][summed_variable[m]];
                break;
              }
            }
            break;
          }

          case REAL_VALUE : {
            vec->real_vector[i][k] = 0.;

            for (m = 0;m < nb_summed_variable;m++) {
              switch (type[summed_variable[m]]) {
              case INT_VALUE :
                vec->real_vector[i][k] += int_vector[i][summed_variable[m]];
                break;
              case REAL_VALUE :
                vec->real_vector[i][k] += real_vector[i][summed_variable[m]];
                break;
              }
            }
            break;
          }
          }
        }

        // copy of values

        else {
          switch (vec->type[k]) {
          case INT_VALUE :
            vec->int_vector[i][k] = int_vector[i][copied_variable[j++]];
            break;
          case REAL_VALUE :
            vec->real_vector[i][k] = real_vector[i][copied_variable[j++]];
            break;
          }
        }
      }
    }

    i = 0;
    for (j = 0;j < inb_variable;j++) {
      if (j == summed_variable[0]) {
        vec->min_value_computation(j);
        vec->max_value_computation(j);

        vec->build_marginal_frequency_distribution(j);
      }

      else {
        vec->min_value[j] = min_value[copied_variable[i]];
        vec->max_value[j] = max_value[copied_variable[i]];
        vec->min_interval[j] = min_interval[copied_variable[i]];

        if (marginal_distribution[copied_variable[i]]) {
          vec->marginal_distribution[j] = new FrequencyDistribution(*(marginal_distribution[copied_variable[i]]));
        }
        if (marginal_histogram[copied_variable[i]]) {
          vec->marginal_histogram[j] = new Histogram(*(marginal_histogram[copied_variable[i]]));
        }

        vec->mean[j] = mean[copied_variable[i]];

        k = 0;
        for (m = 0;m < inb_variable;m++) {
          if (m != summed_variable[0]) {
            vec->covariance[j][m] = covariance[copied_variable[i]][copied_variable[k++]];
            if (m != j) {
              vec->covariance[m][j] = vec->covariance[j][m];
            }
          }
        }

        i++;
      }

      vec->covariance_computation(summed_variable[0]);
    }

    delete [] summed_variable;
    delete [] copied_variable;
  }

  return vec;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Summation of variables.
 *
 *  \param[in] error              reference on a StatError object,
 *  \param[in] nb_summed_variable number of variables to be summed,
 *  \param[in] variable           variable indices.
 *
 *  \return                       Vectors object.
 */
/*--------------------------------------------------------------*/

Vectors* Vectors::sum_variable(StatError &error , int nb_summed_variable , vector<int>ivariable) const

{
  return sum_variable(error , nb_summed_variable , ivariable.data());
}


/*--------------------------------------------------------------*/
/**
 *  \brief Merging of variables of Vectors objects.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] nb_sample  number of Vectors objects,
 *  \param[in] ivec       pointer on the Vectors objects,
 *  \param[in] ref_sample reference Vectors object for the identifiers.
 *
 *  \return               Vectors object.
 */
/*--------------------------------------------------------------*/

Vectors* Vectors::merge_variable(StatError &error , int nb_sample ,
                                 const Vectors **ivec , int ref_sample) const

{
  bool status = true;
  int i , j , k;
  int inb_variable , *iidentifier;
  variable_nature *itype;
  Vectors *vec;
  const Vectors **pvec;


  vec = NULL;
  error.init();

  for (i = 0;i < nb_sample;i++) {
    if (ivec[i]->nb_vector != nb_vector) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                    << STAT_error[STATR_NB_VECTOR];
      error.update((error_message.str()).c_str());
    }
  }

  if ((ref_sample != I_DEFAULT) && ((ref_sample < 1) || (ref_sample > nb_sample + 1))) {
    status = false;
    error.update(STAT_error[STATR_SAMPLE_INDEX]);
  }

  if (status) {
    nb_sample++;
    pvec = new const Vectors*[nb_sample];

    pvec[0] = this;
    inb_variable = nb_variable;
    for (i = 1;i < nb_sample;i++) {
      pvec[i] = ivec[i - 1];
      inb_variable += ivec[i - 1]->nb_variable;
    }

    // comparaison of individual identifiers

    if (ref_sample == I_DEFAULT) {
      for (i = 0;i < nb_vector;i++) {
        for (j = 1;j < nb_sample;j++) {
          if (pvec[j]->identifier[i] != pvec[0]->identifier[i]) {
            break;
          }
        }
        if (j < nb_sample) {
          break;
        }
      }

      if (i < nb_vector) {
        iidentifier = NULL;
      }
      else {
        iidentifier = pvec[0]->identifier;
      }
    }

    else {
      ref_sample--;
      iidentifier = pvec[ref_sample]->identifier;
    }

    itype = new variable_nature[inb_variable];
    inb_variable = 0;
    for (i = 0;i < nb_sample;i++) {
      for (j = 0;j < pvec[i]->nb_variable;j++) {
        itype[inb_variable] = pvec[i]->type[j];
        if ((inb_variable > 0) && (itype[inb_variable] == STATE)) {
          itype[inb_variable] = INT_VALUE;
        }
        inb_variable++;
      }
    }

    vec = new Vectors(nb_vector , iidentifier , inb_variable , itype);
    delete [] itype;

    // copy of vectors

    for (i = 0;i < nb_vector;i++) {
      inb_variable = 0;
      for (j = 0;j < nb_sample;j++) {
        for (k = 0;k < pvec[j]->nb_variable;k++) {
          vec->int_vector[i][inb_variable] = pvec[j]->int_vector[i][k];
          vec->real_vector[i][inb_variable++] = pvec[j]->real_vector[i][k];
        }
      }
    }

    inb_variable = 0;
    for (i = 0;i < nb_sample;i++) {
      for (j = 0;j < pvec[i]->nb_variable;j++) {
        vec->min_value[inb_variable] = pvec[i]->min_value[j];
        vec->max_value[inb_variable] = pvec[i]->max_value[j];
        vec->min_interval[inb_variable] = pvec[i]->min_interval[j];

        if (pvec[i]->marginal_distribution[j]) {
          vec->marginal_distribution[inb_variable] = new FrequencyDistribution(*(pvec[i]->marginal_distribution[j]));
        }
        if (pvec[i]->marginal_histogram[j]) {
          vec->marginal_histogram[inb_variable] = new Histogram(*(pvec[i]->marginal_histogram[j]));
        }

        vec->mean[inb_variable] = pvec[i]->mean[j];
        vec->covariance[inb_variable][inb_variable] = pvec[i]->covariance[j][j];

        inb_variable++;
      }
    }

    vec->covariance_computation();

    delete [] pvec;
  }

  return vec;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Merging of variables of Vectors objects.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] nb_sample  number of Vectors objects,
 *  \param[in] ivec       pointer on the Vectors objects,
 *  \param[in] ref_sample reference Vectors object for the identifiers.
 *
 *  \return               Vectors object.
 */
/*--------------------------------------------------------------*/

Vectors* Vectors::merge_variable(StatError &error , int nb_sample ,
                                 const vector<Vectors> ivec , int ref_sample) const

{
  int i;
  Vectors *vec;
  const Vectors **pvec;


  pvec = new const Vectors*[nb_sample];
  for (i = 0;i < nb_sample;i++) {
    pvec[i] = new Vectors(ivec[i]);
  }

  vec = merge_variable(error , nb_sample , pvec , ref_sample);

  for (i = 0;i < nb_sample;i++) {
    delete pvec[i];
  }
  delete [] pvec;

  return vec;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a Vectors object from a file.
 *         Format: rows: individuals | columns: variables.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 *
 *  \return          Vectors object.
 */
/*--------------------------------------------------------------*/

Vectors* Vectors::ascii_read(StatError &error , const string path)

{
  string buffer;
  size_t position;
  typedef tokenizer<char_separator<char>> tokenizer;
  char_separator<char> separator(" \t");
  bool status , lstatus;
  int i , j;
  int line , read_line , initial_nb_line , nb_variable = 0 , nb_vector , index , int_value;
  variable_nature *type;
  double real_value;
  Vectors *vec;
  ifstream in_file(path.c_str());


  vec = NULL;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {

    // 1st pass: analysis of the line defining the number of variables

    status = true;
    line = 0;
    read_line = 0;

    while (getline(in_file , buffer)) {
      line++;

      position = buffer.find('#');
      if (position != string::npos) {
        buffer.erase(position);
      }
      i = 0;

      tokenizer tok_buffer(buffer , separator);

      for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
        switch (i) {

        // test number of variables

        case 0 : {
          lstatus = true;

/*          try {
            int_value = stoi(*token);   in C++ 11
          }
          catch(invalid_argument &arg) {
            lstatus = false;
          } */
          int_value = atoi(token->c_str());

          if (lstatus) {
            if ((int_value < 1) || (int_value > VECTOR_NB_VARIABLE)) {
              lstatus = false;
            }
            else {
              nb_variable = int_value;
            }
          }

          if (!lstatus) {
            status = false;
            error.update(STAT_parsing[STATP_NB_VARIABLE] , line , i + 1);
          }
          break;
        }

        // test VARIABLE(S) keyword

        case 1 : {
          if (*token != STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES]) {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEYWORD] ,
                                    STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] , line , i + 1);
          }
          break;
        }
        }

        i++;
      }

      if (i > 0) {
        if (i != 2) {
          status = false;
          error.update(STAT_parsing[STATP_FORMAT] , line);
        }

        read_line++;
        break;
      }
    }

    if (read_line < 1) {
      status = false;
      error.update(STAT_parsing[STATP_FORMAT]);
    }

    // analysis of the lines defining the variable types

    if (status) {
      type = new variable_nature[nb_variable];
      for (i = 0;i < nb_variable;i++) {
        type[i] = AUXILIARY;
      }

      read_line = 0;

      while (getline(in_file , buffer)) {
        line++;

        position = buffer.find('#');
        if (position != string::npos) {
          buffer.erase(position);
        }
        i = 0;

        tokenizer tok_buffer(buffer , separator);

        for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
          switch (i) {

          // test VARIABLE keyword

          case 0 : {
            if (*token != STAT_word[STATW_VARIABLE]) {
              status = false;
              error.correction_update(STAT_parsing[STATP_KEYWORD] , STAT_word[STATW_VARIABLE] , line , i + 1);
            }
            break;
          }

          // test variable index

          case 1 : {
            lstatus = true;

/*            try {
              int_value = stoi(*token);   in C++ 11
            }
            catch(invalid_argument &arg) {
              lstatus = false;
            } */
            int_value = atoi(token->c_str());

            if ((lstatus) && (int_value != read_line + 1)) {
              lstatus = false;
            }

            if (!lstatus) {
              status = false;
              error.correction_update(STAT_parsing[STATP_VARIABLE_INDEX] , read_line + 1 , line , i + 1);
            }
            break;
          }

          // test separator

          case 2 : {
            if (*token != ":") {
              status = false;
              error.update(STAT_parsing[STATP_SEPARATOR] , line , i + 1);
            }
            break;
          }

          // test keyword defining the variable type

          case 3 : {
            for (j = INT_VALUE;j <= STATE;j++) {
              if (*token == STAT_variable_word[j]) {
                type[read_line] = (variable_nature)j;
                break;
              }
            }

            if (j == REAL_VALUE + 1) {
              status = false;
              error.update(STAT_parsing[STATP_KEYWORD] , line , i + 1);
            }
            break;
          }
          }

          i++;
        }

        if (i > 0) {
          if (i != 4) {
            status = false;
            error.update(STAT_parsing[STATP_FORMAT] , line);
          }

          read_line++;
          if (read_line == nb_variable) {
            break;
          }
        }
      }

      if (read_line < nb_variable) {
        status = false;
        error.update(STAT_parsing[STATP_FORMAT]);
      }

      initial_nb_line = line;
    }

    if (status) {

      // analysis of the lines and determination of the number of individuals

      line = 0;
      nb_vector = 0;

      while (getline(in_file , buffer)) {
        line++;

#       ifdef DEBUG
        cout << line << "  " << buffer << endl;
#       endif

        position = buffer.find('#');
        if (position != string::npos) {
          buffer.erase(position);
        }
        i = 0;

        tokenizer tok_buffer(buffer , separator);

        for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
          if (i <= nb_variable) {
            lstatus = true;

            if (type[i] != REAL_VALUE) {
/*              try {
                int_value = stoi(*token);   in C++ 11
              }
              catch(invalid_argument &arg) {
                lstatus = false;
              } */
              int_value = atoi(token->c_str());
            }

            else {
/*              try {
                real_value = stod(*token);   in C++ 11
              }
              catch(invalid_argument &arg) {
                lstatus = false;
              } */
              real_value = atof(token->c_str());
            }

            if (!lstatus) {
              status = false;
              error.update(STAT_parsing[STATP_DATA_TYPE] , line , i + 1);
            }
          }

          i++;
        }

        // test constant number of items per line

        if (i > 0) {
          if (i != nb_variable) {
            status = false;
            error.correction_update(STAT_parsing[STATP_NB_TOKEN] , nb_variable , line);
          }

          nb_vector++;
        }
      }

      if (nb_vector == 0) {
        status = false;
        error.update(STAT_parsing[STATP_EMPTY_SAMPLE]);
      }
    }

    // 2nd pass: copy of vectors

    if (status) {
//      in_file.close();
//      in_file.open(path.c_str() , ios::in);

      in_file.clear();
      in_file.seekg(0 , ios::beg);

      line = 0;

      do {
        getline(in_file , buffer);
        line++;

#       ifdef DEBUG
        cout << line << "  " << buffer << endl;
#       endif

      }
      while (line < initial_nb_line);

      vec = new Vectors(nb_vector , NULL , nb_variable , type);

      index = 0;

      while (getline(in_file , buffer)) {
        position = buffer.find('#');
        if (position != string::npos) {
          buffer.erase(position);
        }
        i = 0;

        tokenizer tok_buffer(buffer , separator);

        for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
          if (type[i] != REAL_VALUE) {
//            vec->int_vector[index][i++] = stoi(*token);   in C++ 11
            vec->int_vector[index][i++] = atoi(token->c_str());
          }

          else {
//            vec->real_vector[index][i++] = stod(*token);   in C++ 11
            vec->real_vector[index][i++] = atof(token->c_str());
          }
        }

        if (i > 0) {
          index++;
        }
      }

      for (i = 0;i < vec->nb_variable;i++) {
        vec->min_value_computation(i);
        vec->max_value_computation(i);

        vec->build_marginal_frequency_distribution(i);
        vec->min_interval_computation(i);
      }

      vec->covariance_computation();
    }
  }

  return vec;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing on a single line of a Vectors object.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& Vectors::line_write(ostream &os) const

{
  os << nb_vector << " " << STAT_label[nb_vector == 1 ? STATL_VECTOR : STATL_VECTORS] << "   "
     << nb_variable << " " << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES];

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a Vectors object.
 *
 *  \param[in,out] os           stream,
 *  \param[in]     exhaustive   flag detail level,
 *  \param[in]     comment_flag flag comments.
 */
/*--------------------------------------------------------------*/

ostream& Vectors::ascii_write(ostream &os , bool exhaustive , bool comment_flag) const

{
  int i , j;
  int buff , width[2] , *int_value;
  double median , lower_quartile , upper_quartile , *real_value , **correlation;
  Test *test;
  ios_base::fmtflags format_flags;


  format_flags = os.setf(ios::right , ios::adjustfield);

  if (comment_flag) {
    os << "# ";
  }
  os << nb_vector << " " << STAT_label[nb_vector == 1 ? STATL_VECTOR : STATL_VECTORS] << endl;

  os << "\n" << nb_variable << " " << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] << endl;

  for (i = 0;i < nb_variable;i++) {
    os << "\n" << STAT_word[STATW_VARIABLE] << " " << i + 1 << " : "
       << STAT_variable_word[type[i]];

    os << "   ";
    if (comment_flag) {
      os << "# ";
    }
    os << "("  << STAT_label[STATL_MIN_VALUE] << ": " << min_value[i] << ", "
       << STAT_label[STATL_MAX_VALUE] << ": " << max_value[i] << ")" << endl;

    if (marginal_distribution[i]) {
      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";

      marginal_distribution[i]->ascii_characteristic_print(os , exhaustive , comment_flag);

      if ((marginal_distribution[i]->nb_value <= ASCII_NB_VALUE) || (exhaustive)) {
        os << "\n";
        if (comment_flag) {
          os << "# ";
        }
        os << "   | " << STAT_label[STATL_FREQUENCY] << endl;
        marginal_distribution[i]->ascii_print(os , comment_flag);
      }
    }

    else {
      if (covariance[i][i] > 0.) {
        switch (type[i]) {

        case INT_VALUE : {
          int_value = new int[nb_vector];
          for (j = 0;j < nb_vector;j++) {
            int_value[j] = int_vector[j][i];
          }

          lower_quartile = quantile_computation(nb_vector , int_value , 0.25);
          median = quantile_computation(nb_vector , int_value , 0.5);
          upper_quartile = quantile_computation(nb_vector , int_value , 0.75);

          delete [] int_value;
          break;
        }

        case REAL_VALUE : {
          real_value = new double[nb_vector];
          for (j = 0;j < nb_vector;j++) {
            real_value[j] = real_vector[j][i];
          }

          lower_quartile = quantile_computation(nb_vector , real_value , 0.25);
          median = quantile_computation(nb_vector , real_value , 0.5);
          upper_quartile = quantile_computation(nb_vector , real_value , 0.75);

          delete [] real_value;
          break;
        }
        }
      }

      else {
        median = mean[i];
      }

      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_SAMPLE_SIZE] << ": " << nb_vector << endl;

      if (comment_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_MEAN] << ": " << mean[i] << "   "
         << STAT_label[STATL_MEDIAN] << ": " << median << endl;

      if (comment_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_VARIANCE] << ": " << covariance[i][i] << "   "
         << STAT_label[STATL_STANDARD_DEVIATION] << ": " << sqrt(covariance[i][i]);
      if (covariance[i][i] > 0.) {
        os << "   " << STAT_label[STATL_LOWER_QUARTILE] << ": " << lower_quartile
           << "   " << STAT_label[STATL_UPPER_QUARTILE] << ": " << upper_quartile;
      }
      os << endl;

      if ((exhaustive) && (covariance[i][i] > 0.)) {
        if (comment_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_SKEWNESS_COEFF] << ": " << skewness_computation(i) << "   "
           << STAT_label[STATL_KURTOSIS_COEFF] << ": " << kurtosis_computation(i) << endl;
      }

      if (exhaustive) {
        os << "\n";
        if (comment_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_HISTOGRAM] << endl;

        os << "\n";
        if (comment_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_VALUE] << " | " << STAT_label[STATL_FREQUENCY] << endl;
        marginal_histogram[i]->ascii_print(os , comment_flag);
      }
    }
  }

  width[0] = column_width(nb_variable);

  if (exhaustive) {

    // computation of the column width

    width[1] = 0;
    for (i = 0;i < nb_variable;i++) {
      buff = column_width(nb_variable , covariance[i]);
      if (buff > width[1]) {
        width[1] = buff;
      }
    }
    width[1] += ASCII_SPACE;

    // writing of the variance-covariance matrix

    os << "\n";
    if (comment_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_VARIANCE_COVARIANCE_MATRIX] << endl;

    os << "\n";
    if (comment_flag) {
      os << "# ";
    }
    os << setw(width[0] + width[1]) << 1;
    for (i = 1;i < nb_variable;i++) {
      os << setw(width[1]) << i + 1;
    }
    for (i = 0;i < nb_variable;i++) {
      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << setw(width[0]) << i + 1;
      for (j = 0;j < nb_variable;j++) {
        os << setw(width[1]) << covariance[i][j];
      }
    }
    os << endl;
  }

  correlation = correlation_computation();

  // computation of the column width

  width[1] = 0;
  for (i = 0;i < nb_variable;i++) {
    buff = column_width(nb_variable , correlation[i]);
    if (buff > width[1]) {
      width[1] = buff;
    }
  }
  width[1] += ASCII_SPACE;

  // writing of the correlation matrix

  os << "\n";
  if (comment_flag) {
    os << "# ";
  }
  os << STAT_label[STATL_CORRELATION_MATRIX] << endl;

  os << "\n";
  if (comment_flag) {
    os << "# ";
  }
  os << setw(width[0] + width[1]) << 1;
  for (i = 1;i < nb_variable;i++) {
    os << setw(width[1]) << i + 1;
  }
  for (i = 0;i < nb_variable;i++) {
    os << "\n";
    if (comment_flag) {
      os << "# ";
    }
    os << setw(width[0]) << i + 1;
    for (j = 0;j < nb_variable;j++) {
      os << setw(width[1]) << correlation[i][j];
    }
  }
  os << endl;

  // test significant character of the correlation coefficients

  if (nb_vector > 2) {
#   ifdef DEBUG
    double value = fabs(correlation[0][1]) * sqrt(nb_vector - 2) / sqrt(1 - correlation[0][1] * correlation[0][1]);

    test = new Test(STUDENT , false , nb_vector - 2 , I_DEFAULT , value);
    test->t_critical_probability_computation();

    os << "\n" << test->value << "   " << STAT_label[STATL_CRITICAL_PROBABILITY]
       << ": " << test->critical_probability << endl;

    delete test;
#   endif

    test = new Test(STUDENT , false , nb_vector - 2 , I_DEFAULT , D_DEFAULT);

    for (i = 0;i < NB_CRITICAL_PROBABILITY;i++) {
      test->critical_probability = ref_critical_probability[i];
      test->t_value_computation();

      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_REFERENCE] << " " << STAT_label[STATL_T_VALUE] << ": "
         << test->value << "   " << STAT_label[STATL_REFERENCE] << " "
         << STAT_label[STATL_CRITICAL_PROBABILITY] << ": " << test->critical_probability << endl;

      if (comment_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_LIMIT_CORRELATION_COEFF] << ": "
         << test->value / sqrt(test->value * test->value + nb_vector - 2) << endl;
    }

    delete test;
  }

  for (i = 0;i < nb_variable;i++) {
    delete [] correlation[i];
  }
  delete [] correlation;

  os.setf(format_flags , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a Vectors object.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     exhaustive flag detail level.
 */
/*--------------------------------------------------------------*/

ostream& Vectors::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os, exhaustive, false);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a Vectors object in a file.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] path       file path,
 *  \param[in] exhaustive flag detail level.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool Vectors::ascii_write(StatError &error , const string path , bool exhaustive) const

{
  bool status;
  ofstream out_file(path.c_str());


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    ascii_write(out_file , exhaustive , false);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of vectors.
 *
 *  \param[in,out] os                    stream,
 *  \param[in]     comment_flag          flag comments,
 *  \param[in]     posterior_probability posterior probabilities of the optimal assignments.
 *  \param[in]     entropy               assignment entropies.
 */
/*--------------------------------------------------------------*/

ostream& Vectors::ascii_print(ostream &os , bool comment_flag ,
                              double *posterior_probability , double *entropy) const

{
  int i , j;
  int buff , width[2];
  ios_base::fmtflags format_flags;


  format_flags = os.setf(ios::right , ios::adjustfield);

  width[0] = 0;
  for (i = 0;i < nb_variable;i++) {
    if (type[i] != REAL_VALUE) {
      buff = column_width((int)min_value[i] , (int)max_value[i]);
      if (buff > width[0]) {
        width[0] = buff;
      }
    }

    else {
/*      buff = column_width(min_value[i] , max_value[i]);
      if (buff > width[0]) {
        width[0] = buff;
      } */

      for (j = 0;j < nb_vector;j++) {
        buff = column_width(1 , real_vector[j] + i);
        if (buff > width[0]) {
          width[0] = buff;
        }
      }
    }
  }
  width[1] = width[0] + ASCII_SPACE;

  os << "\n";
  for (i = 0;i < nb_vector;i++) {
    if (type[0] != REAL_VALUE) {
      os << setw(width[0]) << int_vector[i][0];
    }
    else {
      os << setw(width[0]) << real_vector[i][0];
    }

    for (j = 1;j < nb_variable;j++) {
      if (type[j] != REAL_VALUE) {
        os << setw(width[1]) << int_vector[i][j];
      }
      else {
        os << setw(width[1]) << real_vector[i][j];
      }
    }

    os << "   ";
    if (comment_flag) {
      os << "# ";
    }
    os << "(" << identifier[i] << ")";

    if ((posterior_probability) && (entropy)) {
      os << "   " << STAT_label[STATL_POSTERIOR_ASSIGNMENT_PROBABILITY]
         << ": " << posterior_probability[i] << "   "
         << STAT_label[STATL_ASSIGNMENT_ENTROPY] << ": " << entropy[i];
    }
    os << endl;
  }

  os.setf(format_flags , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a Vectors object.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     exhaustive flag detail level.
 */
/*--------------------------------------------------------------*/

ostream& Vectors::ascii_data_write(ostream &os , bool exhaustive) const

{
  ascii_write(os , exhaustive , false);
  ascii_print(os , false);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a Vectors object in a file.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] path       file path,
 *  \param[in] exhaustive flag detail level.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool Vectors::ascii_data_write(StatError &error , const string path , bool exhaustive) const

{
  bool status = false;
  ofstream out_file(path.c_str());


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    ascii_write(out_file , exhaustive , true);
    ascii_print(out_file , true);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a Vectors object in a file at the spreadsheet format.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 *
 *  \return          error status.
 */
/*--------------------------------------------------------------*/

bool Vectors::spreadsheet_write(StatError &error , const string path) const

{
  bool status;
  int i , j;
  int *int_value;
  double median , lower_quartile , upper_quartile , *real_value , **correlation;
  Test *test;
  ofstream out_file(path.c_str());


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    out_file << nb_vector << "\t" << STAT_label[nb_vector == 1 ? STATL_VECTOR : STATL_VECTORS] << endl;

    out_file << "\n" << nb_variable << "\t" << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] << endl;

    for (i = 0;i < nb_variable;i++) {
      out_file << "\n" << STAT_word[STATW_VARIABLE] << "\t" << i + 1<< "\t"
               << STAT_variable_word[type[i]];

      out_file << "\t\t" << STAT_label[STATL_MIN_VALUE] << "\t" << min_value[i]
               << "\t\t" << STAT_label[STATL_MAX_VALUE] << "\t" << max_value[i] << endl;

      if (marginal_distribution[i]) {
        out_file << "\n" << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
        marginal_distribution[i]->spreadsheet_characteristic_print(out_file);

        out_file << "\n\t" << STAT_label[STATL_FREQUENCY] << endl;
        marginal_distribution[i]->spreadsheet_print(out_file);
      }

      else {
        if (covariance[i][i] > 0.) {
          switch (type[i]) {

          case INT_VALUE : {
            int_value = new int[nb_vector];
            for (j = 0;j < nb_vector;j++) {
              int_value[j] = int_vector[j][i];
            }

            lower_quartile = quantile_computation(nb_vector , int_value , 0.25);
            median = quantile_computation(nb_vector , int_value , 0.5);
            upper_quartile = quantile_computation(nb_vector , int_value , 0.75);

            delete [] int_value;
            break;
          }

          case REAL_VALUE : {
            real_value = new double[nb_vector];
            for (j = 0;j < nb_vector;j++) {
              real_value[j] = real_vector[j][i];
            }

            lower_quartile = quantile_computation(nb_vector , real_value , 0.25);
            median = quantile_computation(nb_vector , real_value , 0.5);
            upper_quartile = quantile_computation(nb_vector , real_value , 0.75);

            delete [] real_value;
            break;
          }
          }
        }

        else {
          median = mean[i];
        }

        out_file << "\n" << STAT_label[STATL_SAMPLE_SIZE] << "\t" << nb_vector << endl;

        out_file << STAT_label[STATL_MEAN] << "\t" << mean[i] << "\t\t"
                 << STAT_label[STATL_MEDIAN] << "\t" << median << endl;

        out_file << STAT_label[STATL_VARIANCE] << "\t" << covariance[i][i] << "\t\t"
                 << STAT_label[STATL_STANDARD_DEVIATION] << "\t" << sqrt(covariance[i][i]);
        if (covariance[i][i] > 0.) {
          out_file << "\t\t" << STAT_label[STATL_LOWER_QUARTILE] << "\t" << lower_quartile
                   << "\t\t" << STAT_label[STATL_UPPER_QUARTILE] << "\t" << upper_quartile;
        }
        out_file << endl;

        if (covariance[i][i] > 0.) {
          out_file << STAT_label[STATL_SKEWNESS_COEFF] << "\t" << skewness_computation(i) << "\t\t"
                   << STAT_label[STATL_KURTOSIS_COEFF] << "\t" << kurtosis_computation(i) << endl;
        }

        out_file << "\n" << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_HISTOGRAM] << endl;
        out_file << "\n" << STAT_label[STATL_VALUE] << "\t" << STAT_label[STATL_FREQUENCY] << endl;
        marginal_histogram[i]->spreadsheet_print(out_file);
      }
    }

    // writing of the variance-covariance matrix

    out_file << "\n" << STAT_label[STATL_VARIANCE_COVARIANCE_MATRIX] << endl;

    out_file << "\n";
    for (i = 0;i < nb_variable;i++) {
      out_file << "\t" << i + 1;
    }
    for (i = 0;i < nb_variable;i++) {
      out_file << "\n" << i + 1;
      for (j = 0;j < nb_variable;j++) {
        out_file << "\t" << covariance[i][j];
      }
    }
    out_file << endl;

    // writing of the correlation matrix

    correlation = correlation_computation();

    out_file << "\n" << STAT_label[STATL_CORRELATION_MATRIX] << endl;

    out_file << "\n";
    for (i = 0;i < nb_variable;i++) {
      out_file << "\t" << i + 1;
    }
    for (i = 0;i < nb_variable;i++) {
      out_file << "\n" << i + 1;
      for (j = 0;j < nb_variable;j++) {
        out_file << "\t" << correlation[i][j];
      }
    }
    out_file << endl;

    // test significant character of the correlation coefficients

    if (nb_vector > 2) {
      test = new Test(STUDENT , false , nb_vector - 2 , I_DEFAULT , D_DEFAULT);

      for (i = 0;i < NB_CRITICAL_PROBABILITY;i++) {
        test->critical_probability = ref_critical_probability[i];
        test->t_value_computation();

        out_file << "\n" << STAT_label[STATL_REFERENCE] << " " << STAT_label[STATL_T_VALUE] << "\t"
                 << test->value << "\t" << STAT_label[STATL_REFERENCE] << " "
                 << STAT_label[STATL_CRITICAL_PROBABILITY] << "\t" << test->critical_probability << endl;

        out_file << STAT_label[STATL_LIMIT_CORRELATION_COEFF] << "\t"
                 << test->value / sqrt(test->value * test->value + nb_vector - 2) << endl;
      }

      delete test;
    }

    for (i = 0;i < nb_variable;i++) {
      delete [] correlation[i];
    }
    delete [] correlation;

    // writing of vectors

    for (i = 0;i < nb_vector;i++) {
      out_file << "\n";
      for (j = 0;j < nb_variable;j++) {
        if (type[j] != REAL_VALUE) {
          out_file << int_vector[i][j] << "\t";
        }
        else {
          out_file << real_vector[i][j] << "\t";
        }
      }
      out_file << "\t" << identifier[i];
    }
    out_file << endl;
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of vectors at the Gnuplot format.
 *
 *  \param[in] path              file path,
 *  \param[in] standard_residual standardized residuals.
 *
 *  \return                      error status.
 */
/*--------------------------------------------------------------*/

bool Vectors::plot_print(const char *path , double *standard_residual) const

{
  bool status = false;
  int i , j;
  ofstream out_file(path);


  if (out_file) {
    status = true;

    for (i = 0;i < nb_vector;i++) {
      for (j = 0;j < nb_variable;j++) {
        if (type[j] != REAL_VALUE) {
          out_file << int_vector[i][j] << " ";
        }
        else {
          out_file << real_vector[i][j] << " ";
        }
      }

      if (standard_residual) {
        out_file << standard_residual[i];
      }
      out_file << endl;
    }
  }

  return status;
}



/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a of a Vectors object using Gnuplot.
 *
 *  \param[in] error  reference on a StatError object,
 *  \param[in] prefix file prefix,
 *  \param[in] title  figure title.
 *
 *  \return           error status.
 */
/*--------------------------------------------------------------*/

bool Vectors::plot_write(StatError &error , const char *prefix ,
                         const char *title) const

{
  bool status;
  int i , j , k , m , n;
  int **frequency;
  ostringstream *data_file_name;


  error.init();

  // writing of the data files

  data_file_name = new ostringstream[nb_variable + 1];

  data_file_name[0] << prefix << 0 << ".dat";
  status = plot_print((data_file_name[0].str()).c_str());

  if (!status) {
    error.update(STAT_error[STATR_FILE_PREFIX]);
  }

  else {
    for (i = 0;i < nb_variable;i++) {
      data_file_name[i + 1] << prefix << i + 1 << ".dat";

      if (marginal_distribution[i]) {
        marginal_distribution[i]->plot_print((data_file_name[i + 1].str()).c_str());
      }
      else {
        marginal_histogram[i]->plot_print((data_file_name[i + 1].str()).c_str());
      }
    }

    // writing of the script files

    for (i = 0;i < 2;i++) {
      for (j = 0;j < nb_variable;j++) {
        ostringstream file_name[2];

        switch (i) {

        case 0 : {
          if (nb_variable == 1) {
            file_name[0] << prefix << ".plot";
          }
          else {
            file_name[0] << prefix << j + 1 << ".plot";
          }
          break;
        }

        case 1 : {
          if (nb_variable == 1) {
            file_name[0] << prefix << ".print";
          }
          else {
            file_name[0] << prefix << j + 1 << ".print";
          }
          break;
        }
        }

        ofstream out_file((file_name[0].str()).c_str());

        if (i == 1) {
          out_file << "set terminal postscript" << endl;

          if (nb_variable == 1) {
            file_name[1] << label(prefix) << ".ps";
          }
          else {
            file_name[1] << label(prefix) << j + 1 << ".ps";
          }
          out_file << "set output \"" << file_name[1].str() << "\"\n\n";
        }

        out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                 << "set title";
        if (title) {
          out_file << " \"" << title << "\"";
        }
        out_file << "\n\n";

        if (marginal_distribution[j]) {
          if (marginal_distribution[j]->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }
          if ((int)(marginal_distribution[j]->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics 0,1" << endl;
          }

          out_file << "plot [0:" << MAX(marginal_distribution[j]->nb_value - 1 , 1) << "] [0:"
                   << (int)(marginal_distribution[j]->max * YSCALE) + 1 << "] \""
                   << label((data_file_name[j + 1].str()).c_str()) << "\" using 1 title \""
                   << STAT_label[STATL_VARIABLE] << " " << j + 1 << " "
                   << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                   << "\" with impulses" << endl;

          if (marginal_distribution[j]->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }
          if ((int)(marginal_distribution[j]->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics autofreq" << endl;
          }
        }

        else {
          if ((int)(marginal_histogram[j]->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics 0,1" << endl;
          }

          out_file << "plot [" << marginal_histogram[j]->min_value << ":"
                   << marginal_histogram[j]->max_value << "] [0:"
                   << (int)(marginal_histogram[j]->max * YSCALE) + 1 << "] \""
                   << label((data_file_name[j + 1].str()).c_str()) << "\" using 1:2 title \""
                   << STAT_label[STATL_VARIABLE] << " " << j + 1 << " "
                   << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_HISTOGRAM]
                   << "\" with histeps" << endl;

          if ((int)(marginal_histogram[j]->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics autofreq" << endl;
          }
        }

        if ((i == 0) && (nb_variable > 1)) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;

        for (k = 0;k < nb_variable;k++) {
          if (k != j) {
            out_file << "set xlabel \"" << STAT_label[STATL_VARIABLE] << " " << j + 1 << "\"" << endl;
            out_file << "set ylabel \"" << STAT_label[STATL_VARIABLE] << " " << k + 1 << "\"" << endl;

            if (max_value[j] - min_value[j] < TIC_THRESHOLD) {
              out_file << "set xtics " << MIN(min_value[j] , 0) << ",1" << endl;
            }
            if (max_value[k] - min_value[k] < TIC_THRESHOLD) {
              out_file << "set ytics " << MIN(min_value[k] , 0) << ",1" << endl;
            }

            if (((marginal_distribution[j]) &&
                 (marginal_distribution[j]->nb_value <= PLOT_NB_VALUE)) &&
                ((marginal_distribution[k]) &&
                 (marginal_distribution[k]->nb_value <= PLOT_NB_VALUE))) {
              frequency = joint_frequency_computation(j , k);

              for (m = marginal_distribution[j]->offset;m < marginal_distribution[j]->nb_value;m++) {
                for (n = marginal_distribution[k]->offset;n < marginal_distribution[k]->nb_value;n++) {
                  if (frequency[m][n] > 0) {
                    out_file << "set label \"" << frequency[m][n] << "\" at " << m << ","
                             << n << endl;
                  }
                }
              }

              for (m = 0;m < marginal_distribution[j]->nb_value;m++) {
                delete [] frequency[m];
              }
              delete [] frequency;
            }

            if ((min_value[j] >= 0.) && (max_value[j] - min_value[j] > min_value[j] * PLOT_RANGE_RATIO)) {
              out_file << "plot [" << 0;
            }
            else {
              out_file << "plot [" << min_value[j];
            }
            out_file << ":" << max_value[j] * YSCALE << "] [";
            if ((min_value[k] >= 0.) && (max_value[k] - min_value[k] > min_value[k] * PLOT_RANGE_RATIO)) {
              out_file << 0;
            }
            else {
              out_file << min_value[k];
            }
            out_file << ":" << max_value[k] * YSCALE << "] \""
                     << label((data_file_name[0].str()).c_str()) << "\" using "
                     << j + 1 << ":" << k + 1 << " notitle with points" << endl;

            out_file << "unset label" << endl;

            out_file << "set xlabel" << endl;
            out_file << "set ylabel" << endl;

            if (max_value[j] - min_value[j] < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }
            if (max_value[k] - min_value[k] < TIC_THRESHOLD) {
              out_file << "set ytics autofreq" << endl;
            }

            if ((i == 0) && (((j < nb_variable - 1) && (k < nb_variable - 1)) ||
                 ((j == nb_variable - 1) && (k < nb_variable - 2)))) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;
          }
        }

        if (i == 1) {
          out_file << "\nset terminal x11" << endl;
        }

        out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
      }
    }
  }

  delete [] data_file_name;

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of 2-dimensional vectors for plotting.
 *
 *  \param[in] plot      reference on a SinglePlot object,
 *  \param[in] variable1 variable 1 index.
 *  \param[in] variable2 variable 2 index.
 */
/*--------------------------------------------------------------*/

void Vectors::plotable_write(SinglePlot &plot , int variable1 , int variable2) const

{
  int i;


  if (type[variable1] != REAL_VALUE) {
    if (type[variable2] != REAL_VALUE) {
      for (i = 0;i < nb_vector;i++) {
        plot.add_point(int_vector[i][variable1] , int_vector[i][variable2]);
      }
    }
    else {
      for (i = 0;i < nb_vector;i++) {
        plot.add_point(int_vector[i][variable1] , real_vector[i][variable2]);
      }
    }
  }

  else {
    if (type[variable2] != REAL_VALUE) {
      for (i = 0;i < nb_vector;i++) {
        plot.add_point(real_vector[i][variable1] , int_vector[i][variable2]);
      }
    }
    else {
      for (i = 0;i < nb_vector;i++) {
        plot.add_point(real_vector[i][variable1] , real_vector[i][variable2]);
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of frequencies associated with 2-dimensional vectors.
 *
 *  \param[in] plot      reference on a SinglePlot object,
 *  \param[in] variable1 variable 1 index.
 *  \param[in] variable2 variable 2 index.
 */
/*--------------------------------------------------------------*/

void Vectors::plotable_frequency_write(SinglePlot &plot , int variable1 , int variable2) const

{
  int i , j;
  int **frequency;
  ostringstream label;


  frequency = joint_frequency_computation(variable1 , variable2);

  for (i = marginal_distribution[variable1]->offset;i < marginal_distribution[variable1]->nb_value;i++) {
    for (j = marginal_distribution[variable2]->offset;j < marginal_distribution[variable2]->nb_value;j++) {
      if (frequency[i][j] > 0) {
        label.str("");
        label << frequency[i][j];

        plot.add_text(i , j , label.str());
      }
    }
  }

  for (i = 0;i < marginal_distribution[variable1]->nb_value;i++) {
    delete [] frequency[i];
  }
  delete [] frequency;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a Vectors object.
 *
 *  \return MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* Vectors::get_plotable() const

{
  int i , j , k;
  int nb_plot_set;
  double xmin , ymin;
  ostringstream legend , label;
  MultiPlotSet *plot_set;


  nb_plot_set = nb_variable * nb_variable;

  plot_set = new MultiPlotSet(nb_plot_set , nb_variable);
  MultiPlotSet &plot = *plot_set;

  plot.border = "15 lw 0";

  i = 0;
  for (j = 0;j < nb_variable;j++) {
    plot.variable_nb_viewpoint[j] = 1;
    plot.variable[i] = j;

    if (marginal_distribution[j]) {

      // marginal frequency distribution

      plot[i].xrange = Range(0 , MAX(marginal_distribution[j]->nb_value - 1 , 1));
      plot[i].yrange = Range(0 , ceil(marginal_distribution[j]->max * YSCALE));

      if (marginal_distribution[j]->nb_value - 1 < TIC_THRESHOLD) {
        plot[i].xtics = 1;
      }
      if (ceil(marginal_distribution[j]->max * YSCALE) < TIC_THRESHOLD) {
        plot[i].ytics = 1;
      }

      plot[i].resize(1);

      legend.str("");
      legend << STAT_label[STATL_VARIABLE] << " " << j + 1 << " "
             << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      plot[i][0].legend = legend.str();

      plot[i][0].style = "impulses";

      marginal_distribution[j]->plotable_frequency_write(plot[i][0]);
    }

    else {

      // marginal histogram

      plot[i].xrange = Range(marginal_histogram[j]->min_value , marginal_histogram[j]->max_value);
      plot[i].yrange = Range(0 , ceil(marginal_histogram[j]->max * YSCALE));

      if (ceil(marginal_histogram[j]->max * YSCALE) < TIC_THRESHOLD) {
        plot[i].ytics = 1;
      }

      plot[i].resize(1);

      legend.str("");
      legend << STAT_label[STATL_VARIABLE] << " " << j + 1 << " "
             << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_HISTOGRAM];
      plot[i][0].legend = legend.str();

      plot[i][0].style = "histeps";

      marginal_histogram[j]->plotable_write(plot[i][0]);
    }
    i++;

    for (k = 0;k < nb_variable;k++) {
      if (k != j) {
        plot.variable[i] = j;

         // joint frequency distribution of 2 variables

        if ((min_value[j] >= 0.) && (max_value[j] - min_value[j] > min_value[j] * PLOT_RANGE_RATIO)) {
          xmin = 0.;
        }
        else {
          xmin = min_value[j];
        }
        plot[i].xrange = Range(xmin , max_value[j] * YSCALE);

        if ((min_value[k] >= 0.) && (max_value[k] - min_value[k] > min_value[k] * PLOT_RANGE_RATIO)) {
          ymin = 0.;
        }
        else {
          ymin = min_value[k];
        }
        plot[i].yrange = Range(ymin , max_value[k] * YSCALE);

        if (max_value[j] - min_value[j] < TIC_THRESHOLD) {
          plot[i].xtics = 1;
        }
        if (max_value[k] - min_value[k] < TIC_THRESHOLD) {
          plot[i].ytics = 1;
        }

        label.str("");
        label << STAT_label[STATL_VARIABLE] << " " << j + 1;
        plot[i].xlabel = label.str();

        label.str("");
        label << STAT_label[STATL_VARIABLE] << " " << k + 1;
        plot[i].ylabel = label.str();

        if (((marginal_distribution[j]) &&
             (marginal_distribution[j]->nb_value <= PLOT_NB_VALUE)) &&
            ((marginal_distribution[k]) &&
             (marginal_distribution[k]->nb_value <= PLOT_NB_VALUE))) {
          plot[i].resize(2);
        }
        else {
          plot[i].resize(1);
        }

        plot[i][0].style = "points";

        plotable_write(plot[i][0] , j , k);

        if (((marginal_distribution[j]) &&
             (marginal_distribution[j]->nb_value <= PLOT_NB_VALUE)) &&
            ((marginal_distribution[k]) &&
             (marginal_distribution[k]->nb_value <= PLOT_NB_VALUE))) {
          plot[i][1].label = "true";

          plotable_frequency_write(plot[i][1] , j , k);
        }

        i++;
      }
    }
  }

  return plot_set;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the minimum value taken by a variable.
 *
 *  \param[in] variable variable index.
 */
/*--------------------------------------------------------------*/

void Vectors::min_value_computation(int variable)

{
  int i;


  if (type[variable] != REAL_VALUE) {
    min_value[variable] = int_vector[0][variable];

    for (i = 1;i < nb_vector;i++) {
      if (int_vector[i][variable] < min_value[variable]) {
        min_value[variable] = int_vector[i][variable];
      }
    }
  }

  else {
    min_value[variable] = real_vector[0][variable];

    for (i = 1;i < nb_vector;i++) {
      if (real_vector[i][variable] < min_value[variable]) {
        min_value[variable] = real_vector[i][variable];
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the maximum taken by a variable.
 *
 *  \param[in] variable variable index.
 */
/*--------------------------------------------------------------*/

void Vectors::max_value_computation(int variable)

{
  int i;


  if (type[variable] != REAL_VALUE) {
    max_value[variable] = int_vector[0][variable];

    for (i = 1;i < nb_vector;i++) {
      if (int_vector[i][variable] > max_value[variable]) {
        max_value[variable] = int_vector[i][variable];
      }
    }
  }

  else {
    max_value[variable] = real_vector[0][variable];

    for (i = 1;i < nb_vector;i++) {
      if (real_vector[i][variable] > max_value[variable]) {
        max_value[variable] = real_vector[i][variable];
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of the marginal frequency distribution for
 *         a positive integer-valued variable.
 *
 *  \param[in] variable variable index.
 */
/*--------------------------------------------------------------*/

void Vectors::build_marginal_frequency_distribution(int variable)

{
  if ((type[variable] != REAL_VALUE) && (min_value[variable] >= 0) &&
      (max_value[variable] <= MARGINAL_DISTRIBUTION_MAX_VALUE)) {
    int i;


    marginal_distribution[variable] = new FrequencyDistribution((int)max_value[variable] + 1);

    for (i = 0;i < nb_vector;i++) {
      (marginal_distribution[variable]->frequency[int_vector[i][variable]])++;
    }

    marginal_distribution[variable]->offset = (int)min_value[variable];
    marginal_distribution[variable]->nb_element = nb_vector;
    marginal_distribution[variable]->max_computation();
    marginal_distribution[variable]->mean_computation();
    marginal_distribution[variable]->variance_computation();

    mean[variable] = marginal_distribution[variable]->mean;
    covariance[variable][variable] = marginal_distribution[variable]->variance;
  }

  else {
    build_marginal_histogram(variable);

    mean_computation(variable);
    variance_computation(variable);
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of the marginal histogram for a variable.
 *
 *  \param[in] variable   variable index,
 *  \param[in] bin_width  bin width,
 *  \param[in] imin_value minimum value.
 */
/*--------------------------------------------------------------*/

void Vectors::build_marginal_histogram(int variable , double bin_width , double imin_value)

{
  if ((!marginal_histogram[variable]) || (bin_width != marginal_histogram[variable]->bin_width) ||
      (imin_value != D_INF)) {
    int i;


    // construction of the histogram

    if (bin_width == D_DEFAULT) {
      bin_width = MAX(::round((max_value[variable] - min_value[variable]) * HISTOGRAM_FREQUENCY / nb_vector) , 1);

#     ifdef MESSAGE
      cout << "\n" << STAT_label[STATL_VARIABLE] << " " << variable + 1 << " - "
           << STAT_label[STATL_BIN_WIDTH] << ": " << bin_width << endl;
//           << " (" << min_value[variable] << ", " << max_value[variable] << ")"
#     endif

    }

    if (imin_value == D_INF) {
      imin_value = floor(min_value[variable] / bin_width) * bin_width;
    }

    if (marginal_histogram[variable]) {
      marginal_histogram[variable]->nb_bin = (int)floor((max_value[variable] - imin_value) / bin_width) + 1;

      delete [] marginal_histogram[variable]->frequency;
      marginal_histogram[variable]->frequency = new int[marginal_histogram[variable]->nb_bin];
    }

    else {
      marginal_histogram[variable] = new Histogram((int)floor((max_value[variable] - imin_value) / bin_width) + 1 , false);

      marginal_histogram[variable]->nb_element = nb_vector;
      marginal_histogram[variable]->type = type[variable];
    }

    marginal_histogram[variable]->bin_width = bin_width;
    marginal_histogram[variable]->min_value = imin_value;
    marginal_histogram[variable]->max_value = ceil(max_value[variable] / bin_width) * bin_width;

#   ifdef DEBUG
    cout << "\nTEST: " << marginal_histogram[variable]->min_value << " " << marginal_histogram[variable]->max_value
         << " | " << marginal_histogram[variable]->nb_bin 
        << " " << (marginal_histogram[variable]->max_value - marginal_histogram[variable]->min_value) / marginal_histogram[variable]->bin_width << endl;
#    endif

    // computation of bin frequencies

    for (i = 0;i < marginal_histogram[variable]->nb_bin;i++) {
      marginal_histogram[variable]->frequency[i] = 0;
    }

    if (type[variable] != REAL_VALUE) {
      for (i = 0;i < nb_vector;i++) {
//        (marginal_histogram[variable]->frequency[(int)((int_vector[i][variable] - imin_value) / bin_width)])++;
        (marginal_histogram[variable]->frequency[(int)floor((int_vector[i][variable] - imin_value) / bin_width)])++;
      }
    }

    else {
      for (i = 0;i < nb_vector;i++) {
//        (marginal_histogram[variable]->frequency[(int)((real_vector[i][variable] - imin_value) / bin_width)])++;
        (marginal_histogram[variable]->frequency[(int)floor((real_vector[i][variable] - imin_value) / bin_width)])++;
      }
    }

    marginal_histogram[variable]->max_computation();
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Change of the bin width of the marginal histogram for a variable.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] variable   variable index,
 *  \param[in] bin_width  bin width,
 *  \param[in] imin_value minimum value,
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool Vectors::select_bin_width(StatError &error , int variable , double bin_width ,
                               double imin_value)

{
  bool status = true;


  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if (!marginal_histogram[variable]) {
      status = false;
      error.update(STAT_error[STATR_MARGINAL_HISTOGRAM]);
    }
    if ((bin_width <= 0.) || ((type[variable] != REAL_VALUE) && ((int)bin_width != bin_width))) {
      status = false;
      error.update(STAT_error[STATR_HISTOGRAM_BIN_WIDTH]);
    }
    if ((imin_value != D_INF) && ((imin_value <= min_value[variable] - bin_width) ||
         (imin_value > min_value[variable]) || ((type[variable] != REAL_VALUE) &&
          ((int)imin_value != imin_value)))) {
      status = false;
      error.update(STAT_error[STATR_HISTOGRAM_MIN_VALUE]);
    }
  }

  if (status) {
    build_marginal_histogram(variable , bin_width , imin_value);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of an order on the individuals on the basis of
 *         the values taken by a variable.
 *
 *  \param[in] variable variable index.
 *
 *  \return             ordered individual indices.
 */
/*--------------------------------------------------------------*/

int* Vectors::order_computation(int variable) const

{
  int i , j;
  int int_min , int_value , *index;
  double real_min , real_value;


  index = new int[nb_vector];

  i = 0;

  if (type[variable] != REAL_VALUE) {
    do {

      // determination of the current minimum value

      if (i == 0) {
        int_value = (int)min_value[variable];
      }

      else {
        int_min = (int)max_value[variable] + 1;
        for (j = 0;j < nb_vector;j++) {
          if ((int_vector[j][variable] > int_value) && (int_vector[j][variable] < int_min)) {
            int_min = int_vector[j][variable];
          }
        }
        int_value = int_min;
      }

      // selection of the individuals taken the current minimum value for the selected variable

      for (j = 0;j < nb_vector;j++) {
        if (int_vector[j][variable] == int_value) {
          index[i++] = j;
        }
      }
    }
    while (i < nb_vector);
  }

  else {
    do {

      // determination of the current minimum value

      if (i == 0) {
        real_value = min_value[variable];
      }

      else {
        real_min = max_value[variable] + 1;
        for (j = 0;j < nb_vector;j++) {
          if ((real_vector[j][variable] > real_value) && (real_vector[j][variable] < real_min)) {
            real_min = real_vector[j][variable];
          }
        }
        real_value = real_min;
      }

      // selection of the individuals taken the current minimum value for the selected variable

      for (j = 0;j < nb_vector;j++) {
        if (real_vector[j][variable] == real_value) {
          index[i++] = j;
        }
      }
    }
    while (i < nb_vector);
  }

  return index;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the cumulative frequency distribution function for a variable.
 *
 *  \param[in] variable variable index,
 *  \param[in] cdf      (values, cumulative distribution function).
 *
 *  \return             cumulative frequency distribution function.
 */
/*--------------------------------------------------------------*/

int Vectors::cumulative_distribution_function_computation(int variable , double **cdf) const

{
  int i , j;
  int cumul , int_min , int_value , frequency;
  double real_min , real_value;


  if (marginal_distribution[variable]) {
    i = marginal_distribution[variable]->cumulative_distribution_function_computation(cdf);
  }

  else {
    cdf[0] = new double[nb_vector];
    cdf[1] = new double[nb_vector];

    cumul = 0;
    i = 0;

    if (type[variable] != REAL_VALUE) {

      do {

        // determination of the current minimum value

        if (cumul == 0) {
          int_value = (int)min_value[variable];
        }

        else {
          int_min = (int)max_value[variable] + 1;
          for (j = 0;j < nb_vector;j++) {
            if ((int_vector[j][variable] > int_value) &&
                (int_vector[j][variable] < int_min)) {
              int_min = int_vector[j][variable];
            }
          }
          int_value = int_min;
        }

        // computation of the number of individuals taken the current minimum value for the selected variable

        frequency = 0;
        for (j = 0;j < nb_vector;j++) {
          if (int_vector[j][variable] == int_value) {
            frequency++;
          }
        }

        cdf[0][i] = int_value;
        cdf[1][i] = (cumul + (double)(frequency + 1) / 2.) / (double)nb_vector;
        cumul += frequency;
        i++;
      }
      while (cumul < nb_vector);
    }

    else {
      do {

        // determination of the current minimum value

        if (cumul == 0) {
          real_value = min_value[variable];
        }

        else {
          real_min = max_value[variable] + 1;
          for (j = 0;j < nb_vector;j++) {
            if ((real_vector[j][variable] > real_value) &&
                (real_vector[j][variable] < real_min)) {
              real_min = real_vector[j][variable];
            }
          }
          real_value = real_min;
        }

        // computation of the number of individuals taken the current minimum value for the selected variable

        frequency = 0;
        for (j = 0;j < nb_vector;j++) {
          if (real_vector[j][variable] == real_value) {
            frequency++;
          }
        }

        cdf[0][i] = real_value;
        cdf[1][i] = (cumul + (double)(frequency + 1) / 2.) / (double)nb_vector;
        cumul += frequency;
        i++;
      }
      while (cumul < nb_vector);
    }
  }

  return i;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the shortest interval between 2 successive values for a variable.
 *
 *  \param[in] variable variable index.
 */
/*--------------------------------------------------------------*/

void Vectors::min_interval_computation(int variable)

{
  if (marginal_distribution[variable]) {
    min_interval[variable] = marginal_distribution[variable]->min_interval_computation();
  }

  else {
    int i , j;
    int int_min , int_value;
    double real_min , real_value;


    min_interval[variable] = max_value[variable] - min_value[variable];
    i = 0;

    if (type[variable] != REAL_VALUE) {
      do {

        // determination of the current minimum value

        if (i == 0) {
          int_value = (int)min_value[variable];
        }

        else {
          int_min = (int)max_value[variable] + 1;
          for (j = 0;j < nb_vector;j++) {
            if ((int_vector[j][variable] > int_value) &&
                (int_vector[j][variable] < int_min)) {
              int_min = int_vector[j][variable];
            }
          }

          if (int_min - int_value < min_interval[variable]) {
            min_interval[variable] = int_min - int_value;
          }
          int_value = int_min;
        }

        // computation of the number of individuals taken the current minimum value for the selected variable

        for (j = 0;j < nb_vector;j++) {
          if (int_vector[j][variable] == int_value) {
            i++;
          }
        }
      }
      while (i < nb_vector);
    }

    else {
      do {

        // determination of the current minimum value

        if (i == 0) {
          real_value = min_value[variable];
        }

        else {
          real_min = max_value[variable] + 1;
          for (j = 0;j < nb_vector;j++) {
            if ((real_vector[j][variable] > real_value) &&
                (real_vector[j][variable] < real_min)) {
              real_min = real_vector[j][variable];
            }
          }

          if (real_min - real_value < min_interval[variable]) {
            min_interval[variable] = real_min - real_value;
          }
          real_value = real_min;
        }

        // computation of the number of individuals taken the current minimum value for the selected variable

        for (j = 0;j < nb_vector;j++) {
          if (real_vector[j][variable] == real_value) {
            i++;
          }
        }
      }
      while (i < nb_vector);
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Mean computation for a variable.
 *
 *  \param[in] variable variable index.
 */
/*--------------------------------------------------------------*/

void Vectors::mean_computation(int variable)

{
  int i;


  mean[variable] = 0.;

  if (type[variable] != REAL_VALUE) {
    for (i = 0;i < nb_vector;i++) {
      mean[variable] += int_vector[i][variable];
    }
  }

  else {
    for (i = 0;i < nb_vector;i++) {
      mean[variable] += real_vector[i][variable];
    }
  }

  mean[variable] /= nb_vector;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Variance computation for a variable.
 *
 *  \param[in] variable variable index.
 */
/*--------------------------------------------------------------*/

void Vectors::variance_computation(int variable)

{
  if (mean[variable] != D_INF) {
    if (nb_vector > 1) {
      int i;
      double diff;
      long double square_sum;


      square_sum = 0.;

      if (type[variable] != REAL_VALUE) {
        for (i = 0;i < nb_vector;i++) {
          diff = int_vector[i][variable] - mean[variable];
          square_sum += diff * diff;
        }
      }

      else{
        for (i = 0;i < nb_vector;i++) {
          diff = real_vector[i][variable] - mean[variable];
          square_sum += diff * diff;
        }
      }

      covariance[variable][variable] = square_sum / (nb_vector - 1);
    }

    else {
      covariance[variable][variable] = 0.;
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the variance-covariance matrix.
 *
 *  \param[in] variable variable index.
 */
/*--------------------------------------------------------------*/

void Vectors::covariance_computation(int variable)

{
  if (mean[0] != D_INF) {
    int i , j , k;
    long double square_sum;


    for (i = 0;i < nb_variable;i++) {
      if ((variable == I_DEFAULT) || (i == variable)) {
        for (j = (variable == I_DEFAULT ? i + 1 : 0);j < nb_variable;j++) {
          if (nb_vector > 1) {
            square_sum = 0.;

            if ((type[i] != REAL_VALUE) && (type[j] != REAL_VALUE)) {
              for (k = 0;k < nb_vector;k++) {
                square_sum += (int_vector[k][i] - mean[i]) * (int_vector[k][j] - mean[j]);
              }
            }
            else if ((type[i] != REAL_VALUE) && (type[j] == REAL_VALUE)) {
              for (k = 0;k < nb_vector;k++) {
                square_sum += (int_vector[k][i] - mean[i]) * (real_vector[k][j] - mean[j]);
              }
            }
            else if ((type[i] == REAL_VALUE) && (type[j] != REAL_VALUE)) {
              for (k = 0;k < nb_vector;k++) {
                square_sum += (real_vector[k][i] - mean[i]) * (int_vector[k][j] - mean[j]);
              }
            }
//            else if ((type[i] == REAL_VALUE) && (type[j] == REAL_VALUE)) {
            else {
              for (k = 0;k < nb_vector;k++) {
                square_sum += (real_vector[k][i] - mean[i]) * (real_vector[k][j] - mean[j]);
              }
            }

            covariance[i][j] = square_sum / (nb_vector - 1);
          }

          else {
            covariance[i][j] = 0.;
          }

          if (j != i) {
            covariance[j][i] = covariance[i][j];
          }
        }
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the mean absolute deviation for a variable.
 *
 *  \param[in] variable variable index,
 *  \param[in] location location measure (e.g. mean or median).
 *
 *  \return             mean absolute deviation.
 */
/*--------------------------------------------------------------*/

double Vectors::mean_absolute_deviation_computation(int variable , double location) const

{
  int i;
  double mean_absolute_deviation;


  if (marginal_distribution[variable]) {
    mean_absolute_deviation = marginal_distribution[variable]->mean_absolute_deviation_computation(location);
  }

  else {
    mean_absolute_deviation = 0.;

    if (type[variable] != REAL_VALUE) {
      for (i = 0;i < nb_vector;i++) {
        mean_absolute_deviation += fabs(int_vector[i][variable] - location);
      }
    }

    else {
      for (i = 0;i < nb_vector;i++) {
        mean_absolute_deviation += fabs(real_vector[i][variable] - location);
      }
    }

    mean_absolute_deviation /= nb_vector;
  }

  return mean_absolute_deviation;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the mean absolute difference for a variable.
 *
 *  \param[in] variable variable index.
 *
 *  \return             mean absolute difference.
 */
/*--------------------------------------------------------------*/

double Vectors::mean_absolute_difference_computation(int variable) const

{
  int i , j;
  double mean_absolute_difference;


  mean_absolute_difference = 0.;

  if (nb_vector > 1) {
    if (type[variable] != REAL_VALUE) {
      for (i = 0;i < nb_vector;i++) {
        for (j = i + 1;j < nb_vector;j++) {
          mean_absolute_difference += abs(int_vector[i][variable] - int_vector[j][variable]);
        }
      }
    }

    else {
      for (i = 0;i < nb_vector;i++) {
        for (j = i + 1;j < nb_vector;j++) {
          mean_absolute_difference += fabs(real_vector[i][variable] - real_vector[j][variable]);
        }
      }
    }

    mean_absolute_difference = 2 * mean_absolute_difference / (nb_vector * (nb_vector - 1));
  }

  return mean_absolute_difference;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the coefficient of skewness for a variable.
 *
 *  \param[in] variable variable index.
 *
 *  \return             coefficient of skewness.
 */
/*--------------------------------------------------------------*/

double Vectors::skewness_computation(int variable) const

{
  int i;
  double skewness = D_INF , diff;
  long double cube_sum;


  if ((mean[variable] != D_INF) && (covariance[variable][variable] != D_DEFAULT)) {
    if ((nb_vector > 2) && (covariance[variable][variable] > 0.)) {
      cube_sum = 0.;

      if (type[variable] != REAL_VALUE) {
        for (i = 0;i < nb_vector;i++) {
          diff = int_vector[i][variable] - mean[variable];
          cube_sum += diff * diff * diff;
        }
      }

      else {
        for (i = 0;i < nb_vector;i++) {
          diff = real_vector[i][variable] - mean[variable];
          cube_sum += diff * diff * diff;
        }
      }

      skewness = cube_sum * nb_vector / ((nb_vector - 1) * (nb_vector - 2) *
                  pow(covariance[variable][variable] , 1.5));
    }

    else {
      skewness = 0.;
    }
  }

  return skewness;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the excess kurtosis for a variable:
 *         excess kurtosis = coefficient of kurtosis - 3.
 *
 *  \param[in] variable variable index.
 *
 *  \return             excess kurtosis.
 */
/*--------------------------------------------------------------*/

double Vectors::kurtosis_computation(int variable) const

{
  int i;
  double kurtosis = D_INF , diff;
  long double power_sum;


  if ((mean[variable] != D_INF) && (covariance[variable][variable] != D_DEFAULT)) {
    if (covariance[variable][variable] > 0.) {
      power_sum = 0.;

      if (type[variable] != REAL_VALUE) {
        for (i = 0;i < nb_vector;i++) {
          diff = int_vector[i][variable] - mean[variable];
          power_sum += diff * diff * diff * diff;
        }
      }

      else {
        for (i = 0;i < nb_vector;i++) {
          diff = real_vector[i][variable] - mean[variable];
          power_sum += diff * diff * diff * diff;
        }
      }

      kurtosis = power_sum / ((nb_vector - 1) * covariance[variable][variable] *
                  covariance[variable][variable]) - 3.;
    }

    else {
      kurtosis = -2.;
    }
  }

  return kurtosis;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the correlation matrix from
 *         the variance-covariance matrix.
 *
 *  \return correlation matrix.
 */
/*--------------------------------------------------------------*/

double** Vectors::correlation_computation() const

{
  int i , j;
  double norm , **correlation = NULL;


  if (covariance[0][0] != D_DEFAULT) {
    correlation = new double*[nb_variable];
    for (i = 0;i < nb_variable;i++) {
      correlation[i] = new double[nb_variable];
    }

    for (i = 0;i < nb_variable;i++) {
      correlation[i][i] = 1.;
      for (j = i + 1;j < nb_variable;j++) {
        correlation[i][j] = covariance[i][j];
        norm = covariance[i][i] * covariance[j][j];
        if (norm > 0.) {
          correlation[i][j] /= sqrt(norm);
        }
        correlation[j][i] = correlation[i][j];
      }
    }
  }

  return correlation;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the mean direction for a circular variable.
 *
 *  \param[in] variable variable index,
 *  \param[in] unit     unit (DEGREE/RADIAN).
 *
 *  \return             mean direction.
 */
/*--------------------------------------------------------------*/

double* Vectors::mean_direction_computation(int variable , angle_unit unit) const

{
  int i;
  double *mean_direction;


  mean_direction = new double[4];
//  mean_direction = new double[2];

  mean_direction[0] = 0.;
  mean_direction[1] = 0.;

  switch (type[variable]) {

  case INT_VALUE : {
    switch (unit) {

    case DEGREE : {
      for (i = 0;i < nb_vector;i++) {
        mean_direction[0] += cos(int_vector[i][variable] * M_PI / 180);
        mean_direction[1] += sin(int_vector[i][variable] * M_PI / 180);
      }
      break;
    }

    case RADIAN : {
      for (i = 0;i < nb_vector;i++) {
        mean_direction[0] += cos(int_vector[i][variable]);
        mean_direction[1] += sin(int_vector[i][variable]);
      }
      break;
    }
    }
    break;
  }

  case REAL_VALUE : {
    switch (unit) {

    case DEGREE : {
      for (i = 0;i < nb_vector;i++) {
        mean_direction[0] += cos(real_vector[i][variable] * M_PI / 180);
        mean_direction[1] += sin(real_vector[i][variable] * M_PI / 180);
      }
      break;
    }

    case RADIAN : {
      for (i = 0;i < nb_vector;i++) {
        mean_direction[0] += cos(real_vector[i][variable]);
        mean_direction[1] += sin(real_vector[i][variable]);
      }
      break;
    }
    }
    break;
  }
  }

  mean_direction[0] /= nb_vector;
  mean_direction[1] /= nb_vector;

  mean_direction[2] = sqrt(mean_direction[0] * mean_direction[0] +
                           mean_direction[1] * mean_direction[1]);

  if (mean_direction[2] > 0.) {
    mean_direction[3] = atan(mean_direction[1] / mean_direction[0]);

    if (mean_direction[0] < 0.) {
      mean_direction[3] += M_PI;
    }
    if (unit == DEGREE) {
      mean_direction[3] *= (180 / M_PI);
    }
  }

  else {
    mean_direction[3] = D_DEFAULT;
  }

  return mean_direction;
}


};  // namespace stat_tool
