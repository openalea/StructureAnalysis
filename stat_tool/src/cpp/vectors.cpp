/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2002 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): Y. Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id$
 *
 *       Forum for AMAPmod developers: amldevlp@cirad.fr
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
#include <sstream>
#include <iomanip>

#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"
#include "tool/config.h"

// #include <rw/vstream.h>
// #include <rw/rwfile.h>
#include "stat_tools.h"
#include "distribution.h"
#include "vectors.h"
#include "stat_label.h"

using namespace std;

extern int column_width(int value);
extern int column_width(int min_value , int max_value);
extern int column_width(int nb_value , const double *value , double scale = 1.);
extern char* label(const char *file_name);



/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe Vectors.
 *
 *--------------------------------------------------------------*/

Vectors::Vectors()

{
  nb_variable = 0;

  min_value = 0;
  max_value = 0;
  marginal = 0;

  mean = 0;
  covariance = 0;

  nb_vector = 0;
  identifier = 0;
  vector = 0;
}


/*--------------------------------------------------------------*
 *
 *  Initialisation d'un objet Vectors.
 *
 *  arguments : nombre de variables, nombre de vecteurs, identificateurs des vecteurs,
 *              flag initialisation.
 *
 *--------------------------------------------------------------*/

void Vectors::init(int inb_variable , int inb_vector , int *iidentifier , bool init_flag)

{
  register int i , j;
  int *pvector;


  nb_variable = inb_variable;

  min_value = new int[nb_variable];
  max_value = new int[nb_variable];
  marginal = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    min_value[i] = 0;
    max_value[i] = 0;
    marginal[i] = 0;
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

  vector = new int*[nb_vector];
  for (i = 0;i < nb_vector;i++) {
    vector[i] = new int[nb_variable];

    if (init_flag) {
      pvector = vector[i];
      for (j = 0;j < nb_variable;j++) {
        *pvector++ = 0;
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Vectors.
 *
 *  arguments : nombre de variables, nombre de vecteurs, vecteurs,
 *              identificateurs des vecteurs.
 *
 *--------------------------------------------------------------*/

Vectors::Vectors(int inb_variable , int inb_vector , int **ivector , int *iidentifier)

{
  register int i , j;
  int *pvector , *cvector;


  init(inb_variable , inb_vector , iidentifier , false);

  for (i = 0;i < nb_vector;i++) {
    pvector = vector[i];
    cvector = ivector[i];
    for (j = 0;j < nb_variable;j++) {
      *pvector++ = *cvector++;
    }
  }

  for (i = 0;i < nb_variable;i++) {
    min_value_computation(i);
    max_value_computation(i);
    build_marginal_histogram(i);
  }

  covariance_computation();
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Vectors.
 *
 *  arguments : reference sur un objet Vectors, nombre de vecteurs,
 *              index des vecteurs selectionnes.
 *
 *--------------------------------------------------------------*/

Vectors::Vectors(const Vectors &vec , int inb_vector , int *index)

{
  register int i , j;
  int *pvector , *cvector;


  nb_variable = vec.nb_variable;

  min_value = new int[nb_variable];
  max_value = new int[nb_variable];
  marginal = new Histogram*[nb_variable];

  mean = new double[nb_variable];

  covariance = new double*[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    covariance[i] = new double[nb_variable];
  }

  nb_vector = inb_vector;

  identifier = new int[nb_vector];
  for (i = 0;i < nb_vector;i++) {
    identifier[i] = vec.identifier[index[i]];
  }

  vector = new int*[nb_vector];
  for (i = 0;i < nb_vector;i++) {
    vector[i] = new int[nb_variable];
    pvector = vector[i];
    cvector = vec.vector[index[i]];
    for (j = 0;j < nb_variable;j++) {
      *pvector++ = *cvector++;
    }
  }

  for (i = 0;i < nb_variable;i++) {
    min_value_computation(i);
    max_value_computation(i);
    marginal[i] = 0;
    build_marginal_histogram(i);
  }

  covariance_computation();
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Vectors.
 *
 *  argument : reference sur un objet Vectors.
 *
 *--------------------------------------------------------------*/

void Vectors::copy(const Vectors &vec)

{
  register int i , j;
  int *pvector , *cvector;


  nb_variable = vec.nb_variable;

  min_value = new int[nb_variable];
  max_value = new int[nb_variable];
  marginal = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    min_value[i] = vec.min_value[i];
    max_value[i] = vec.max_value[i];

    if (vec.marginal[i]) {
      marginal[i] = new Histogram(*(vec.marginal[i]));
    }
    else {
      marginal[i] = 0;
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

  nb_vector = vec.nb_vector;

  identifier = new int[nb_vector];
  for (i = 0;i < nb_vector;i++) {
    identifier[i] = vec.identifier[i];
  }

  vector = new int*[nb_vector];
  for (i = 0;i < nb_vector;i++) {
    vector[i] = new int[nb_variable];
    pvector = vector[i];
    cvector = vec.vector[i];
    for (j = 0;j < nb_variable;j++) {
      *pvector++ = *cvector++;
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet Vectors.
 *
 *--------------------------------------------------------------*/

void Vectors::remove()

{
  register int i;


  delete [] min_value;
  delete [] max_value;

  if (marginal) {
    for (i = 0;i < nb_variable;i++) {
      delete marginal[i];
    }
    delete [] marginal;
  }

  delete [] mean;

  if (covariance) {
    for (i = 0;i < nb_variable;i++) {
      delete [] covariance[i];
    }
    delete [] covariance;
  }

  delete [] identifier;

  if (vector) {
    for (i = 0;i < nb_vector;i++) {
      delete [] vector[i];
    }
    delete [] vector;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Vectors.
 *
 *--------------------------------------------------------------*/

Vectors::~Vectors()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Vectors.
 *
 *  argument : reference sur un objet Vectors.
 *
 *--------------------------------------------------------------*/

Vectors& Vectors::operator=(const Vectors &vec)

{
  if (&vec != this) {
    remove();
    copy(vec);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Extraction des histogrammes correspondant aux lois marginales.
 *
 *  arguments : reference sur un objet Format_error, variable.
 *
 *--------------------------------------------------------------*/

Distribution_data* Vectors::extract(Format_error &error , int variable) const

{
  Distribution_data *histo;


  histo = 0;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if (marginal[variable]) {
      histo = new Distribution_data(*(marginal[variable]));
    }
    else {
      error.update(STAT_error[STATR_MARGINAL_HISTOGRAM]);
    }
  }

  return histo;
}


/*--------------------------------------------------------------*
 *
 *  Verification d'une liste d'identificateurs.
 *
 *  arguments : reference sur un objet Format_error, nombre de formes,
 *              identificateurs des formes.
 *
 *--------------------------------------------------------------*/

bool identifier_checking(Format_error &error , int nb_pattern , int *identifier)

{
  bool status = true , *selected_identifier;
  register int i;
  int max_identifier;


  max_identifier = 1;
  for (i = 0;i < nb_pattern;i++) {
    if (identifier[i] > max_identifier) {
      max_identifier = identifier[i];
    }
  }

  selected_identifier = new bool[max_identifier + 1];
  for (i = 0;i <= max_identifier;i++) {
    selected_identifier[i] = false;
  }

  for (i = 0;i < nb_pattern;i++) {
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


/*--------------------------------------------------------------*
 *
 *  Verification d'un objet Vectors.
 *
 *  argument : reference sur un objet Format_error.
 *
 *--------------------------------------------------------------*/

bool Vectors::check(Format_error &error)

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


/*--------------------------------------------------------------*
 *
 *  Fusion d'objets Vectors.
 *
 *  argument : reference sur un objet Format_error, nombre d'objets Vectors,
 *             pointeurs sur les objets Vectors.
 *
 *--------------------------------------------------------------*/

Vectors* Vectors::merge(Format_error &error , int nb_sample , const Vectors **ivec) const

{
  bool status = true;
  register int i , j , k , m;
  int inb_vector , *pvector , *cvector;
  const Histogram **phisto;
  Vectors *vec;
  const Vectors **pvec;


  vec = 0;
  error.init();

  for (i = 0;i < nb_sample;i++) {
    if (ivec[i]->nb_variable != nb_variable) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                    << STAT_error[STATR_NB_VARIABLE];
      error.correction_update((error_message.str()).c_str() , nb_variable);
    }
  }

  if (status) {
    nb_sample++;
    pvec = new const Vectors*[nb_sample];

    pvec[0] = this;
    for (i = 1;i < nb_sample;i++) {
      pvec[i] = ivec[i - 1];
    }

    inb_vector = 0;
    for (i = 0;i < nb_sample;i++) {
      inb_vector += pvec[i]->nb_vector;
    }

    vec = new Vectors(nb_variable , inb_vector , 0 , false);

    // copie des vecteurs

    i = 0;
    for (j = 0;j < nb_sample;j++) {
      for (k = 0;k < pvec[j]->nb_vector;k++) {
        pvector = vec->vector[i++];
        cvector = pvec[j]->vector[k];
        for (m = 0;m < vec->nb_variable;m++) {
          *pvector++ = *cvector++;
        }
      }
    }

    phisto = new const Histogram*[nb_sample];

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

      for (j = 0;j < nb_sample;j++) {
        if (pvec[j]->marginal[i]) {
          phisto[j] = pvec[j]->marginal[i];
        }
        else {
          break;
        }
      }

      if (j == nb_sample) {
        vec->marginal[i] = new Histogram(nb_sample , phisto);
        vec->mean[i] = vec->marginal[i]->mean;
        vec->covariance[i][i] = vec->marginal[i]->variance;
      }

      else {
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


/*--------------------------------------------------------------*
 *
 *  Translation des valeurs d'une variable donnee.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable,
 *              parametre de translation.
 *
 *--------------------------------------------------------------*/

Vectors* Vectors::shift(Format_error &error , int variable , int shift_param) const

{
  bool status = true;
  register int i , j;
  int *pvector , *cvector;
  Vectors *vec;


  vec = 0;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

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

  if (status) {
    vec = new Vectors(nb_variable , nb_vector , identifier , false);

    for (i = 0;i < vec->nb_vector;i++) {
      pvector = vec->vector[i];
      cvector = vector[i];
      for (j = 0;j < vec->nb_variable;j++) {
        if (j == variable) {
          *pvector++ = *cvector++ + shift_param;
        }
        else {
          *pvector++ = *cvector++;
        }
      }
    }

    for (i = 0;i < vec->nb_variable;i++) {
      if (i == variable) {
        vec->min_value[i] = min_value[i] + shift_param;
        vec->max_value[i] = max_value[i] + shift_param;

        if ((marginal[i]) && (vec->min_value[i] >= 0) && (vec->max_value[i] <= MARGINAL_MAX_VALUE)) {
          vec->marginal[i] = new Histogram(*(marginal[i]) , 's' , shift_param);
        }
        else {
          vec->build_marginal_histogram(i);
        }

        vec->mean[i] = mean[i] + shift_param;
      }

      else {
        vec->min_value[i] = min_value[i];
        vec->max_value[i] = max_value[i];

        if (marginal[i]) {
          vec->marginal[i] = new Histogram(*(marginal[i]));
        }
        vec->mean[i] = mean[i];
      }

      for (j = 0;j < vec->nb_variable;j++) {
        vec->covariance[i][j] = covariance[i][j];
      }
    }
  }

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Regroupement des valeurs d'une variable donnee.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable,
 *              pas de regroupement.
 *
 *--------------------------------------------------------------*/

Vectors* Vectors::cluster(Format_error &error , int variable , int step) const

{
  bool status = true;
  register int i , j;
  int *pvector , *cvector;
  Vectors *vec;


  vec = 0;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }
  if (step < 1) {
    status = false;
    error.update(STAT_error[STATR_CLUSTERING_STEP]);
  }

  if (status) {
    variable--;

    vec = new Vectors(nb_variable , nb_vector , identifier , false);

    for (i = 0;i < vec->nb_vector;i++) {
      pvector = vec->vector[i];
      cvector = vector[i];
      for (j = 0;j < vec->nb_variable;j++) {
        if (j == variable) {
          *pvector++ = *cvector++ / step;
        }
        else {
          *pvector++ = *cvector++;
        }
      }
    }

    for (i = 0;i < vec->nb_variable;i++) {
      if (i == variable) {
        vec->min_value[i] = min_value[i] / step;
        vec->max_value[i] = max_value[i] / step;

        if (marginal[i]) {
          vec->marginal[i] = new Histogram(*(marginal[i]) , 'c' , step);
          vec->mean[i] = vec->marginal[i]->mean;
          vec->covariance[i][i] = vec->marginal[i]->variance;
        }
        else {
          vec->build_marginal_histogram(i);
        }
      }

      else {
        vec->min_value[i] = min_value[i];
        vec->max_value[i] = max_value[i];

        if (marginal[i]) {
          vec->marginal[i] = new Histogram(*(marginal[i]));
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


/*--------------------------------------------------------------*
 *
 *  Transcodage des symboles d'une variable donnee.
 *
 *  arguments : reference sur un objet Vectors, indice de la variable,
 *              plus petit et plus grand symboles, table de transcodage des symboles.
 *
 *--------------------------------------------------------------*/

void Vectors::transcode(const Vectors &vec , int variable , int min_symbol ,
                        int max_symbol , int *symbol)

{
  register int i , j;
  int *pvector , *cvector;


  for (i = 0;i < nb_vector;i++) {
    pvector = vector[i];
    cvector = vec.vector[i];
    for (j = 0;j < nb_variable;j++) {
      if (j == variable) {
        *pvector++ = symbol[*cvector++ - vec.min_value[variable]] + min_symbol;
      }
      else {
        *pvector++ = *cvector++;
      }
    }
  }

  for (i = 0;i < nb_variable;i++) {
    if (i == variable) {
      min_value[i] = min_symbol;
      max_value[i] = max_symbol;
      build_marginal_histogram(i);
    }

    else {
      min_value[i] = vec.min_value[i];
      max_value[i] = vec.max_value[i];

      if (vec.marginal[i]) {
        marginal[i] = new Histogram(*(vec.marginal[i]));
      }
      mean[i] = vec.mean[i];
      for (j = 0;j < nb_variable;j++) {
        covariance[i][j] = vec.covariance[i][j];
      }
    }
  }

  covariance_computation(variable);
}


/*--------------------------------------------------------------*
 *
 *  Transcodage des symboles d'une variable donnee.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable,
 *              table de transcodage des symboles.
 *
 *--------------------------------------------------------------*/

Vectors* Vectors::transcode(Format_error &error , int variable , int *symbol) const

{
  bool status = true , *presence;
  register int i;
  int min_symbol , max_symbol;
  Vectors *vec;


  vec = 0;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    min_symbol = INT_MAX;
    max_symbol = INT_MIN;

    for (i = 0;i <= max_value[variable] - min_value[variable];i++) {
      if (symbol[i] < min_symbol) {
        min_symbol = symbol[i];
      }
      if (symbol[i] > max_symbol) {
        max_symbol = symbol[i];
      }
    }

    if (max_symbol - min_symbol == 0) {
      status = false;
      error.update(STAT_error[STATR_NB_SYMBOL]);
    }

    if (max_symbol - min_symbol > max_value[variable] - min_value[variable]) {
      status = false;
      error.update(STAT_error[STATR_NON_CONSECUTIVE_SYMBOLS]);
    }

    if (status) {
      presence = new bool[max_symbol - min_symbol + 1];
      for (i = 0;i <= max_symbol - min_symbol;i++) {
        presence[i] = false;
      }

      for (i = 0;i <= max_value[variable] - min_value[variable];i++) {
        presence[symbol[i] - min_symbol] = true;
      }

      for (i = 0;i <= max_symbol - min_symbol;i++) {
        if (!presence[i]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_error[STATR_MISSING_SYMBOL] << " " << i + min_symbol;
          error.update((error_message.str()).c_str());
        }
      }

      delete [] presence;
    }

    if (status) {
      for (i = 0;i <= max_value[variable] - min_value[variable];i++) {
        symbol[i] -= min_symbol;
      }

      vec = new Vectors(nb_variable , nb_vector , identifier , false);
      vec->transcode(*this , variable , min_symbol , max_symbol , symbol);
    }
  }

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Regroupement des symboles d'une variable donnee.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable,
 *              nombre de classes, bornes pour regrouper les symboles.
 *
 *--------------------------------------------------------------*/

Vectors* Vectors::cluster(Format_error &error , int variable , int nb_class ,
                          int *ilimit) const

{
  bool status = true;
  register int i , j , k;
  int *limit , *symbol;
  Vectors *vec;


  vec = 0;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if ((nb_class < 2) || (nb_class >= max_value[variable] - min_value[variable] + 1)) {
      status = false;
      error.update(STAT_error[STATR_NB_CLASS]);
    }
  }

  if (status) {
    limit = new int[nb_class + 1];
    limit[0] = min_value[variable];
    for (i = 1;i < nb_class;i++) {
      limit[i] = ilimit[i - 1];
    }
    limit[nb_class] = max_value[variable] + 1;

    for (i = 0;i < nb_class;i++) {
      if (limit[i] >= limit[i + 1]) {
        status = false;
        error.update(STAT_error[STATR_CLUSTER_LIMIT]);
      }
    }

    if (status) {
      symbol = new int[max_value[variable] - min_value[variable] + 1];

      i = 0;
      for (j = 0;j < nb_class;j++) {
        for (k = limit[j];k < limit[j + 1];k++) {
          symbol[i++] = j;
        }
      }

      vec = new Vectors(nb_variable , nb_vector , identifier , false);
      vec->transcode(*this , variable , 0 , nb_class - 1 , symbol);

      delete [] symbol;
    }

    delete [] limit;
  }

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Selection de vecteurs sur la valeur prise par une variable donnee.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable,
 *              bornes sur les valeurs, flag pour conserver ou rejeter
 *              les vecteurs selectionnes.
 *
 *--------------------------------------------------------------*/

Vectors* Vectors::value_select(Format_error &error , int variable , int imin_value ,
                               int imax_value , bool keep) const

{
  bool status = true;
  register int i;
  int inb_vector , *index;
  Vectors *vec;


  vec = 0;
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

    // selection des vecteurs

    index = new int[nb_vector];
    inb_vector = 0;

    for (i = 0;i < nb_vector;i++) {
      if ((vector[i][variable] >= imin_value) && (vector[i][variable] <= imax_value)) {
        if (keep) {
          index[inb_vector++] = i;
        }
      }

      else if (!keep) {
        index[inb_vector++] = i;
      }
    }

    if (inb_vector == 0) {
      status = false;
      error.update(STAT_error[STATR_EMPTY_SAMPLE]);
    }

    // copie des vecteurs

    if (status) {
      vec = new Vectors(*this , inb_vector , index);
    }

    delete [] index;
  }

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Selection de formes par l'identificateur.
 *
 *  arguments : nombre de formes, identificateurs des formes,
 *              nombre de formes selectionnees, identificateurs des formes selectionnees,
 *              flag pour conserver ou rejeter les formes selectionnees.
 *
 *--------------------------------------------------------------*/

int* identifier_select(int nb_pattern , int *pattern_identifier , int selected_nb_pattern ,
                       int *selected_identifier , bool keep)

{
  register int i , j , k;
  int *index;


  switch (keep) {

  case false : {
    index = new int[nb_pattern - selected_nb_pattern];

    i = 0;
    for (j = 0;j < nb_pattern;j++) {
      for (k = 0;k < selected_nb_pattern;k++) {
        if (selected_identifier[k] == pattern_identifier[j]) {
          break;
        }
      }

      if (k == selected_nb_pattern) {
        index[i++] = j;
      }
    }
    break;
  }

  case true : {
    index = new int[selected_nb_pattern];

    for (i = 0;i < selected_nb_pattern;i++) {
      for (j = 0;j < nb_pattern;j++) {
        if (selected_identifier[i] == pattern_identifier[j]) {
          index[i] = j;
          break;
        }
      }
    }
    break;
  }
  }

  return index;
}


/*--------------------------------------------------------------*
 *
 *  Selection de vecteurs par l'identificateur.
 *
 *  arguments : reference sur un objet Format_error, nombre de vecteurs,
 *              identificateurs des vecteurs, flag pour conserver ou rejeter
 *              les vecteurs selectionnees.
 *
 *--------------------------------------------------------------*/

Vectors* Vectors::select_individual(Format_error &error , int inb_vector ,
                                    int *iidentifier , bool keep) const

{
  bool status = true , *selected_vector;
  register int i , j;
  int max_identifier , *index;
  Vectors *vec;


  vec = 0;
  error.init();

  if ((inb_vector < 1) || (inb_vector > (keep ? nb_vector : nb_vector - 1))) {
    status = false;
    error.update(STAT_error[STATR_NB_VECTOR]);
  }

  else {
    max_identifier = 1;
    for (i = 0;i < inb_vector;i++) {
      if (iidentifier[i] > max_identifier) {
        max_identifier = iidentifier[i];
      }
    }

    selected_vector = new bool[max_identifier + 1];
    for (i = 0;i <= max_identifier;i++) {
      selected_vector[i] = false;
    }

    for (i = 0;i < inb_vector;i++) {
      for (j = 0;j < nb_vector;j++) {
        if (iidentifier[i] == identifier[j]) {
          break;
        }
      }

      if (j == nb_vector) {
        status = false;
        ostringstream error_message;
        error_message << iidentifier[i] << ": " << STAT_error[STATR_VECTOR_IDENTIFIER];
        error.update((error_message.str()).c_str());
      }

      else if (selected_vector[iidentifier[i]]) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VECTOR] << " " << iidentifier[i] << " "
                      << STAT_error[STATR_ALREADY_SELECTED];
        error.update((error_message.str()).c_str());
      }
      else {
        selected_vector[iidentifier[i]] = true;
      }
    }

    delete [] selected_vector;
  }

  if (status) {
    index = identifier_select(nb_vector , identifier , inb_vector , iidentifier , keep);

    vec = new Vectors(*this , (keep ? inb_vector : nb_vector - inb_vector) , index);

    delete [] index;
  }

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Selection de variables.
 *
 *  arguments : reference sur un objet Vectors, indices des variables selectionnees.
 *
 *--------------------------------------------------------------*/

void Vectors::select_variable(const Vectors &vec , int *variable)

{
  register int i , j;
  int *pvector;


  for (i = 0;i < nb_vector;i++) {
    pvector = vector[i];
    for (j = 0;j < nb_variable;j++) {
      *pvector++ = vec.vector[i][variable[j]];
    }
  }

  for (i = 0;i < nb_variable;i++) {
    min_value[i] = vec.min_value[variable[i]];
    max_value[i] = vec.max_value[variable[i]];

    if (vec.marginal[variable[i]]) {
      marginal[i] = new Histogram(*(vec.marginal[variable[i]]));
    }
    mean[i] = vec.mean[variable[i]];
    for (j = 0;j < nb_variable;j++) {
      covariance[i][j] = vec.covariance[variable[i]][variable[j]];
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Selection de 2 variables.
 *
 *  arguments : variable explicative, variable expliquee.
 *
 *--------------------------------------------------------------*/

Vectors* Vectors::select_variable(int explanatory_variable , int response_variable) const

{
  int variable[2];
  Vectors *vec;


  variable[0] = explanatory_variable;
  variable[1] = response_variable;

  vec = new Vectors(2 , nb_vector , identifier , false);

  vec->select_variable(*this , variable);

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Selection de variables.
 *
 *  arguments : nombre de variables, nombre de variables selectionnees,
 *              indices des variables selectionnees, flag pour conserver ou
 *              rejeter les variables selectionnees.
 *
 *--------------------------------------------------------------*/

int* select_variable(int nb_variable , int selected_nb_variable ,
                     int *selected_variable , bool keep)

{
  register int i , j , k;
  int *variable;


  for (i = 0;i < selected_nb_variable;i++) {
    selected_variable[i]--;
  }

  switch (keep) {

  case false : {
    variable = new int[nb_variable - selected_nb_variable];

    i = 0;
    for (j = 0;j < nb_variable;j++) {
      for (k = 0;k < selected_nb_variable;k++) {
        if (selected_variable[k] == j) {
          break;
        }
      }

      if (k == selected_nb_variable) {
        variable[i++] = j;
      }
    }
    break;
  }

  case true : {
    variable = new int[selected_nb_variable];

    for (i = 0;i < selected_nb_variable;i++) {
      variable[i] = selected_variable[i];
    }
    break;
  }
  }

  return variable;
}


/*--------------------------------------------------------------*
 *
 *  Selection de variables.
 *
 *  arguments : reference sur un objet Format_error, nombre de variables,
 *              indices des variables, flag pour conserver ou rejeter
 *              les variables selectionnees.
 *
 *--------------------------------------------------------------*/

Vectors* Vectors::select_variable(Format_error &error , int inb_variable ,
                                  int *ivariable , bool keep) const

{
  bool status = true , *selected_variable;
  register int i;
  int *variable;
  Vectors *vec;


  vec = 0;
  error.init();

  if ((inb_variable < 1) || (inb_variable > (keep ? nb_variable : nb_variable - 1))) {
    status = false;
    error.update(STAT_error[STATR_NB_SELECTED_VARIABLE]);
  }

  if (status) {
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
    variable = ::select_variable(nb_variable , inb_variable , ivariable , keep);

    vec = new Vectors((keep ? inb_variable : nb_variable - inb_variable) ,
                      nb_vector , identifier , false);

    vec->select_variable(*this , variable);
    delete [] variable;
  }

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Concatenation des variables d'objets Vectors.
 *
 *  arguments : reference sur un objet Format_error, nombre d'objets Vectors,
 *              pointeurs sur les objets Vectors, echantillon de reference pour les identificateurs.
 *
 *--------------------------------------------------------------*/

Vectors* Vectors::merge_variable(Format_error &error , int nb_sample ,
                                 const Vectors **ivec , int ref_sample) const

{
  bool status = true;
  register int i , j , k;
  int inb_variable , *iidentifier , *pvector , *cvector;
  Vectors *vec;
  const Vectors **pvec;


  vec = 0;
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

    // calcul des identificateurs des vecteurs

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
        iidentifier = 0;
      }
      else {
        iidentifier = pvec[0]->identifier;
      }
    }

    else {
      iidentifier = pvec[--ref_sample]->identifier;
    }

    vec = new Vectors(inb_variable , nb_vector , iidentifier , false);

    // copie des vecteurs

    for (i = 0;i < nb_vector;i++) {
      pvector = vec->vector[i];
      for (j = 0;j < nb_sample;j++) {
        cvector = pvec[j]->vector[i];
        for (k = 0;k < pvec[j]->nb_variable;k++) {
          *pvector++ = *cvector++;
        }
      }
    }

    inb_variable = 0;
    for (i = 0;i < nb_sample;i++) {
      for (j = 0;j < pvec[i]->nb_variable;j++) {
        vec->min_value[inb_variable] = pvec[i]->min_value[j];
        vec->max_value[inb_variable] = pvec[i]->max_value[j];

        if (pvec[i]->marginal[j]) {
          vec->marginal[inb_variable] = new Histogram(*(pvec[i]->marginal[j]));
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


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Vectors a partir d'un fichier.
 *  Format : chaque ligne represente un individu.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

Vectors* vectors_ascii_read(Format_error &error , const char *path)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  bool status , lstatus;
  register int i;
  int line , nb_variable , nb_vector , index;
  long value;
  Vectors *vec;
  ifstream in_file(path);


  vec = 0;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {

    // 1ere passe : analyse des lignes et recherche du nombre de variable et
    // du nomber de vecteurs

    status = true;
    line = 0;
    nb_vector = 0;

    while (buffer.readLine(in_file , false)) {
      line++;

#     ifdef DEBUG
      cout << line << "  " << buffer << endl;
#     endif

      position = buffer.first('#');
      if (position != RW_NPOS) {
        buffer.remove(position);
      }
      i = 0;

      RWCTokenizer next(buffer);

      while (!((token = next()).isNull())) {
        if ((nb_vector == 0) || (i <= nb_variable)) {
          lstatus = locale.stringToNum(token , &value);
          if (!lstatus) {
            status = false;
            error.update(STAT_parsing[STATP_DATA_TYPE] , line , i + 1);
          }
        }

        i++;
      }

      // test nombre de valeurs par ligne constant

      if (i > 0) {
        if (nb_vector == 0) {
          nb_variable = i;
        }
        else if (i != nb_variable) {
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

    // 2eme passe : copie des vecteurs

    if (status) {
//      in_file.close();
//      in_file.open(path , ios::in);

      in_file.clear();
      in_file.seekg(0 , ios::beg);

      vec = new Vectors(nb_variable , nb_vector , 0 , false);

      index = 0;

      while (buffer.readLine(in_file , false)) {
        position = buffer.first('#');
        if (position != RW_NPOS) {
          buffer.remove(position);
        }
        i = 0;

        RWCTokenizer next(buffer);

        while (!((token = next()).isNull())) {
          locale.stringToNum(token , &value);
          vec->vector[index][i++] = value;
        }

        if (i > 0) {
          index++;
        }
      }

      for (i = 0;i < vec->nb_variable;i++) {
        vec->min_value_computation(i);
        vec->max_value_computation(i);
        vec->build_marginal_histogram(i);
      }

      vec->covariance_computation();
    }
  }

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet Vectors.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& Vectors::line_write(ostream &os) const

{
  os << nb_variable << " " << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] << "   "
     << nb_vector << " " << STAT_label[nb_vector == 1 ? STATL_VECTOR : STATL_VECTORS];

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Vectors.
 *
 *  arguments : stream, flag niveau de detail, flag commentaire.
 *
 *--------------------------------------------------------------*/

ostream& Vectors::ascii_write(ostream &os , bool exhaustive , bool comment_flag) const

{
  register int i , j;
  int buff , width[2];
  long old_adjust;
  double **correlation;
  Test *test;


  old_adjust = os.setf(ios::right , ios::adjustfield);

  if (comment_flag) {
    os << "# ";
  }
  os << nb_variable << " " << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] << endl;

  if (comment_flag) {
    os << "# ";
  }
  os << nb_vector << " " << STAT_label[nb_vector == 1 ? STATL_VECTOR : STATL_VECTORS] << endl;

  for (i = 0;i < nb_variable;i++) {
    os << "\n";
    if (comment_flag) {
      os << "# ";
    }
    os << STAT_word[STATW_VARIABLE] << " " << i + 1 << "   (" << STAT_label[STATL_MIN_VALUE] << ": "
       << min_value[i] << ", " << STAT_label[STATL_MAX_VALUE] << ": " << max_value[i] << ")" << endl;

    os << "\n";
    if (comment_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_HISTOGRAM] << " - ";

    if (marginal[i]) {
      marginal[i]->ascii_characteristic_print(os , exhaustive , comment_flag);

      if ((marginal[i]->nb_value <= ASCII_NB_VALUE) || (exhaustive)) {
        os << "\n";
        if (comment_flag) {
          os << "# ";
        }
        os << "   | " << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_HISTOGRAM] << endl;
        marginal[i]->ascii_print(os , comment_flag);
      }
    }

    else {
      os << STAT_label[STATL_SAMPLE_SIZE] << ": " << nb_vector << endl;

      if (comment_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_MEAN] << ": " << mean[i] << "   "
         << STAT_label[STATL_VARIANCE] << ": " << covariance[i][i] << "   "
         << STAT_label[STATL_STANDARD_DEVIATION] << ": " << sqrt(covariance[i][i]) << endl;

      if ((exhaustive) && (covariance[i][i] > 0.)) {
        if (comment_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_SKEWNESS_COEFF] << ": " << skewness_computation(i) << "   "
           << STAT_label[STATL_KURTOSIS_COEFF] << ": " << kurtosis_computation(i) << endl;
      }
    }
  }

  width[0] = column_width(nb_variable);

  if (exhaustive) {

    // calcul des largeurs des colonnes

    width[1] = 0;
    for (i = 0;i < nb_variable;i++) {
      buff = column_width(nb_variable , covariance[i]);
      if (buff > width[1]) {
        width[1] = buff;
      }
    }
    width[1] += ASCII_SPACE;

    // ecriture de la matrice de variance-covariance

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

  // calcul des largeurs des colonnes

  width[1] = 0;
  for (i = 0;i < nb_variable;i++) {
    buff = column_width(nb_variable , correlation[i]);
    if (buff > width[1]) {
      width[1] = buff;
    }
  }
  width[1] += ASCII_SPACE;

  // ecriture de la matrice des coefficients de correlation

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

  // test du caractere significatif des coefficients de correlation

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

  for (i = 0;i < nb_variable;i++) {
    delete [] correlation[i];
  }
  delete [] correlation;

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Vectors.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Vectors::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , exhaustive , false);
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Vectors dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Vectors::ascii_write(Format_error &error , const char *path ,
                          bool exhaustive) const

{
  bool status;
  ofstream out_file(path);


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


/*--------------------------------------------------------------*
 *
 *  Ecriture des vecteurs discrets.
 *
 *  arguments : stream, flag commentaire.
 *
 *--------------------------------------------------------------*/

ostream& Vectors::ascii_print(ostream &os , bool comment_flag) const

{
  register int i , j;
  int buff , *pvector , width[2];
  long old_adjust;


  old_adjust = os.setf(ios::right , ios::adjustfield);

  width[0] = 0;
  for (i = 0;i < nb_variable;i++) {
    buff = column_width(min_value[i] , max_value[i]);
    if (buff > width[0]) {
      width[0] = buff;
    }
  }
  width[1] = width[0] + ASCII_SPACE;

  os << "\n";
  for (i = 0;i < nb_vector;i++) {
    pvector = vector[i];
    os << setw(width[0]) << *pvector++;
    for (j = 1;j < nb_variable;j++) {
      os << setw(width[1]) << *pvector++;
    }

    os << "   ";
    if (comment_flag) {
      os << "# ";
    }
    os << "(" << identifier[i] << ")" << endl;
  }

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Vectors.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Vectors::ascii_data_write(ostream &os , bool exhaustive) const

{
  ascii_write(os , exhaustive , false);
  ascii_print(os , false);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Vectors dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Vectors::ascii_data_write(Format_error &error , const char *path , bool exhaustive) const

{
  bool status = false;
  ofstream out_file(path);


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


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Vectors dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

bool Vectors::spreadsheet_write(Format_error &error , const char *path) const

{
  bool status;
  register int i , j;
  int *pvector;
  double **correlation;
  Test *test;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    out_file << nb_variable << "\t" << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] << endl;
    out_file << nb_vector << "\t" << STAT_label[nb_vector == 1 ? STATL_VECTOR : STATL_VECTORS] << endl;

    for (i = 0;i < nb_variable;i++) {
      out_file << "\n" << STAT_word[STATW_VARIABLE] << "\t" << i + 1 << "\t\t" << STAT_label[STATL_MIN_VALUE]
               << "\t" << min_value[i] << "\t\t" << STAT_label[STATL_MAX_VALUE] << "\t" << max_value[i] << endl;

      out_file << "\n" << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_HISTOGRAM] << "\t";

      if (marginal[i]) {
        marginal[i]->spreadsheet_characteristic_print(out_file);

        out_file << "\n\t" << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_HISTOGRAM] << endl;
        marginal[i]->spreadsheet_print(out_file);
      }

      else {
        out_file << STAT_label[STATL_SAMPLE_SIZE] << "\t" << nb_vector << endl;

        out_file << STAT_label[STATL_MEAN] << "\t" << mean[i] << "\t\t"
                 << STAT_label[STATL_VARIANCE] << "\t" << covariance[i][i] << "\t\t"
                 << STAT_label[STATL_STANDARD_DEVIATION] << "\t" << sqrt(covariance[i][i]) << endl;

        if (covariance[i][i] > 0.) {
          out_file << STAT_label[STATL_SKEWNESS_COEFF] << "\t" << skewness_computation(i) << "\t\t"
                   << STAT_label[STATL_KURTOSIS_COEFF] << "\t" << kurtosis_computation(i) << endl;
        }
      }
    }

    // ecriture de la matrice de variance-covariance

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

    // ecriture de la matrice des coefficients de correlation

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

    // test du caractere significatif des coefficients de correlation

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

    for (i = 0;i < nb_variable;i++) {
      delete [] correlation[i];
    }
    delete [] correlation;

    // ecriture des vecteurs

    for (i = 0;i < nb_vector;i++) {
      pvector = vector[i];
      out_file << "\n";
      for (j = 0;j < nb_variable;j++) {
        out_file << *pvector++ << "\t";
      }
      out_file << "\t" << identifier[i];
    }
    out_file << endl;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des vecteurs discrets au format Gnuplot.
 *
 *  arguments : path, residus reduits.
 *
 *--------------------------------------------------------------*/

bool Vectors::plot_print(const char *path , double *standard_residual) const

{
  bool status = false;
  register int i , j;
  int *pvector;
  ofstream out_file(path);


  if (out_file) {
    status = true;

    for (i = 0;i < nb_vector;i++) {
      pvector = vector[i];
      for (j = 0;j < nb_variable;j++) {
        out_file << *pvector++ << " ";
      }
      if (standard_residual) {
        out_file << standard_residual[i];
      }
      out_file << endl;
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet Vectors.
 *
 *  arguments : reference sur un objet Format_error, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Vectors::plot_write(Format_error &error , const char *prefix ,
                         const char *title) const

{
  bool status;
  register int i , j , k , m , n;
  int nb_histo , histo_index , **frequency;
  const Histogram **phisto;
  ostringstream data_file_name[2];


  error.init();

  // ecriture des fichiers de donnees

  data_file_name[0] << prefix << 0 << ".dat";
  status = plot_print((data_file_name[0].str()).c_str());

  if (!status) {
    error.update(STAT_error[STATR_FILE_PREFIX]);
  }

  else {
    phisto = new const Histogram*[nb_variable];

    nb_histo = 0;
    for (i = 0;i < nb_variable;i++) {
      if (marginal[i]) {
        phisto[nb_histo++] = marginal[i];
      }
    }

    if (nb_histo > 0) {
      data_file_name[1] << prefix << 1 << ".dat";
      phisto[0]->plot_print((data_file_name[1].str()).c_str() , nb_histo - 1 , phisto + 1);
    }

    delete [] phisto;

    // ecriture des fichiers de commandes et des fichier d'impression

    histo_index = 1;

    for (i = 0;i < nb_variable;i++) {
      for (j = 0;j < 2;j++) {
        ostringstream file_name[2];

        switch (j) {

        case 0 : {
          if (nb_variable == 1) {
            file_name[0] << prefix << ".plot";
          }
          else {
            file_name[0] << prefix << i + 1 << ".plot";
          }
          break;
        }

        case 1 : {
          if (nb_variable == 1) {
            file_name[0] << prefix << ".print";
          }
          else {
            file_name[0] << prefix << i + 1 << ".print";
          }
          break;
        }
        }

        ofstream out_file((file_name[0].str()).c_str());

        if (j == 1) {
          out_file << "set terminal postscript" << endl;

          if (nb_variable == 1) {
            file_name[1] << label(prefix) << ".ps";
          }
          else {
            file_name[1] << label(prefix) << i + 1 << ".ps";
          }
          out_file << "set output \"" << file_name[1].str() << "\"\n\n";
        }

        out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                 << "set title";
        if (title) {
          out_file << " \"" << title << "\"";
        }
        out_file << "\n\n";

        if (marginal[i]) {
          if (marginal[i]->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }
          if ((int)(marginal[i]->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics 0,1" << endl;
          }

          out_file << "plot [0:" << MAX(marginal[i]->nb_value - 1 , 1) << "] [0:"
                   << (int)(marginal[i]->max * YSCALE) + 1 << "] \""
                   << label((data_file_name[1].str()).c_str()) << "\" using " << histo_index
                   << " title \"" << STAT_label[STATL_VARIABLE] << " " << i + 1 << " "
                   << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_HISTOGRAM]
                   << "\" with impulses" << endl;

          if (marginal[i]->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }
          if ((int)(marginal[i]->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics autofreq" << endl;
          }

          if ((j == 0) && (nb_variable > 1)) {
            out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
          out_file << endl;
        }

        for (k = 0;k < nb_variable;k++) {
          if (k != i) {
            out_file << "set xlabel \"" << STAT_label[STATL_VARIABLE] << " " << i + 1 << "\"" << endl;
            out_file << "set ylabel \"" << STAT_label[STATL_VARIABLE] << " " << k + 1 << "\"" << endl;

            if (max_value[i] - min_value[i] < TIC_THRESHOLD) {
              out_file << "set xtics " << MIN(min_value[i] , 0) << ",1" << endl;
            }
            if (max_value[k] - min_value[k] < TIC_THRESHOLD) {
              out_file << "set ytics " << MIN(min_value[k] , 0) << ",1" << endl;
            }

            if (((marginal[i]) && (marginal[i]->nb_value <= PLOT_NB_VALUE)) &&
                ((marginal[k]) && (marginal[k]->nb_value <= PLOT_NB_VALUE))) {
              frequency = joint_frequency_computation(i , k);

              for (m = marginal[i]->offset;m < marginal[i]->nb_value;m++) {
                for (n = marginal[k]->offset;n < marginal[k]->nb_value;n++) {
                  if (frequency[m][n] > 0) {
                    out_file << "set label \"" << frequency[m][n] << "\" at " << m << ","
                             << n << endl;
                  }
                }
              }

              for (m = 0;m < marginal[i]->nb_value;m++) {
                delete [] frequency[m];
              }
              delete [] frequency;
            }

            if ((min_value[i] >= 0) && (max_value[i] - min_value[i] > min_value[i] * PLOT_RANGE_RATIO)) {
              out_file << "plot [" << 0;
            }
            else {
              out_file << "plot [" << min_value[i];
            }
            out_file << ":" << MAX(max_value[i] , min_value[i] + 1) << "] [";
            if ((min_value[k] >= 0) && (max_value[k] - min_value[k] > min_value[k] * PLOT_RANGE_RATIO)) {
              out_file << 0;
            }
            else {
              out_file << min_value[k];
            }
            out_file << ":" << MAX(max_value[k] , min_value[k] + 1) << "] \""
                     << label((data_file_name[0].str()).c_str()) << "\" using "
                     << i + 1 << ":" << k + 1 << " notitle with points" << endl;

            out_file << "unset label" << endl;

            out_file << "set xlabel" << endl;
            out_file << "set ylabel" << endl;

            if (max_value[i] - min_value[i] < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }
            if (max_value[k] - min_value[k] < TIC_THRESHOLD) {
              out_file << "set ytics autofreq" << endl;
            }

            if ((j == 0) && (((i < nb_variable - 1) && (k < nb_variable - 1)) ||
                 ((i == nb_variable - 1) && (k < nb_variable - 2)))) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;
          }
        }

        if (j == 1) {
          out_file << "\nset terminal x11" << endl;
        }

        out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
      }

      if (marginal[i]) {
        histo_index++;
      }
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Fonctions pour la persistance.
 *
 *--------------------------------------------------------------*/

/* RWDEFINE_COLLECTABLE(Vectors , STATI_VECTORS);


RWspace Vectors::binaryStoreSize() const

{
  register int i;
  RWspace size;


  size = sizeof(nb_variable) + sizeof(min_value) * nb_variable + sizeof(max_value) * nb_variable;
  for (i = 0;i < nb_variable;i++) {
    size += sizeof(true);
    if (marginal[i]) {
      size += marginal[i]->binaryStoreSize();
    }
  }

  size += sizeof(mean) * nb_variable + sizeof(covariance) * nb_variable * nb_variable +
          sizeof(nb_vector) + sizeof(identifier) * nb_vector +
          sizeof(vector) * nb_vector * nb_variable;

  return size;
}


void Vectors::restoreGuts(RWvistream &is)

{
  bool status;
  register int i , j;
  int *pvector;


  remove();

  is >> nb_variable;

  min_value = new int[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    is >> min_value[i];
  }
  max_value = new int[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    is >> max_value[i];
  }

  marginal = new Histogram*[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    is >> status;
    if (status) {
      marginal[i] = new Histogram();
      marginal[i]->restoreGuts(is);
    }
    else {
      marginal[i] = 0;
    }
  }

  mean = new double[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    is >> mean[i];
  }

  covariance = new double*[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    covariance[i] = new double[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      is >> covariance[i][j];
    }
  }

  is >> nb_vector;

  identifier = new int[nb_vector];
  for (i = 0;i < nb_vector;i++) {
    is >> identifier[i];
  }

  vector = new int*[nb_vector];
  for (i = 0;i < nb_vector;i++) {
    vector[i] = new int[nb_variable];
    pvector = vector[i];
    for (j = 0;j < nb_variable;j++) {
      is >> *pvector++;
    }
  }
}


void Vectors::restoreGuts(RWFile &file)

{
  bool status;
  register int i;


  remove();

  file.Read(nb_variable);

  min_value = new int[nb_variable];
  file.Read(min_value , nb_variable);
  max_value = new int[nb_variable];
  file.Read(max_value , nb_variable);

  marginal = new Histogram*[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    file.Read(status);
    if (status) {
      marginal[i] = new Histogram();
      marginal[i]->restoreGuts(file);
    }
    else {
      marginal[i] = 0;
    }
  }

  mean = new double[nb_variable];
  file.Read(mean , nb_variable);

  covariance = new double*[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    covariance[i] = new double[nb_variable];
    file.Read(covariance[i] , nb_variable);
  }

  file.Read(nb_vector);

  identifier = new int[nb_vector];
  file.Read(identifier , nb_vector);

  vector = new int*[nb_vector];
  for (i = 0;i < nb_vector;i++) {
    vector[i] = new int[nb_variable];
    file.Read(vector[i] , nb_variable);
  }
}


void Vectors::saveGuts(RWvostream &os) const

{
  register int i , j;
  int *pvector;


  os << nb_variable;

  for (i = 0;i < nb_variable;i++) {
    os << min_value[i];
  }
  for (i = 0;i < nb_variable;i++) {
    os << max_value[i];
  }

  for (i = 0;i < nb_variable;i++) {
    if (marginal[i]) {
      os << true;
      marginal[i]->saveGuts(os);
    }
    else {
      os << false;
    }
  }

  for (i = 0;i < nb_variable;i++) {
    os << mean[i];
  }
  for (i = 0;i < nb_variable;i++) {
    for (j = 0;j < nb_variable;j++) {
      os << covariance[i][j];
    }
  }

  os << nb_vector;

  for (i = 0;i < nb_vector;i++) {
    os << identifier[i];
  }

  for (i = 0;i < nb_vector;i++) {
    pvector = vector[i];
    for (j = 0;j < nb_variable;j++) {
      os << *pvector++;
    }
  }
}


void Vectors::saveGuts(RWFile &file) const

{
  register int i;


  file.Write(nb_variable);

  file.Write(min_value , nb_variable);
  file.Write(max_value , nb_variable);

  for (i = 0;i < nb_variable;i++) {
    if (marginal[i]) {
      file.Write(true);
      marginal[i]->saveGuts(file);
    }
    else {
      file.Write(false);
    }
  }

  file.Write(mean , nb_variable);
  for (i = 0;i < nb_variable;i++) {
    file.Write(covariance[i] , nb_variable);
  }

  file.Write(nb_vector);

  file.Write(identifier , nb_vector);

  for (i = 0;i < nb_vector;i++) {
    file.Write(vector[i] , nb_variable);
  }
} */


/*--------------------------------------------------------------*
 *
 *  Calcul de la valeur minimum prise par une variable.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void Vectors::min_value_computation(int variable)

{
  register int i;


  min_value[variable] = vector[0][variable];

  for (i = 1;i < nb_vector;i++) {
    if (vector[i][variable] < min_value[variable]) {
      min_value[variable] = vector[i][variable];
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la valeur maximum prise par une variable.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void Vectors::max_value_computation(int variable)

{
  register int i;


  max_value[variable] = vector[0][variable];

  for (i = 1;i < nb_vector;i++) {
    if (vector[i][variable] > max_value[variable]) {
      max_value[variable] = vector[i][variable];
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Construction des echantillons correspondants au lois marginales.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void Vectors::build_marginal_histogram(int variable)

{
  if ((min_value[variable] >= 0) && (max_value[variable] <= MARGINAL_MAX_VALUE)) {
    register int i;


    marginal[variable] = new Histogram(max_value[variable] + 1);

    for (i = 0;i < nb_vector;i++) {
      (marginal[variable]->frequency[vector[i][variable]])++;
    }

    marginal[variable]->offset = min_value[variable];
    marginal[variable]->nb_element = nb_vector;
    marginal[variable]->max_computation();
    marginal[variable]->mean_computation();
    marginal[variable]->variance_computation();

    mean[variable] = marginal[variable]->mean;
    covariance[variable][variable] = marginal[variable]->variance;
  }

  else {
    mean_computation(variable);
    variance_computation(variable);
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul d'un ordre sur les vecteurs a partir des valeurs prises
 *  par une variable donnee.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

int* Vectors::order_computation(int variable) const

{
  register int i , j;
  int min , value , *index;


  index = new int[nb_vector];

  i = 0;
  do {

    // recherche de la valeur minimum courante

    if (i == 0) {
      value = min_value[variable];
    }

    else {
      min = max_value[variable] + 1;
      for (j = 0;j < nb_vector;j++) {
        if ((vector[j][variable] > value) && (vector[j][variable] < min)) {
          min = vector[j][variable];
        }
      }
      value = min;
    }

    // recherche des vecteurs prenant pour la variable selectionnee
    // la valeur minimum courante

    for (j = 0;j < nb_vector;j++) {
      if (vector[j][variable] == value) {
        index[i++] = j;
      }
    }
  }
  while (i < nb_vector);

  return index;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la moyenne d'une variable.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void Vectors::mean_computation(int variable)

{
  register int i;


  mean[variable] = 0.;

  for (i = 0;i < nb_vector;i++) {
    mean[variable] += vector[i][variable];
  }

  mean[variable] /= nb_vector;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la variance d'une variable.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void Vectors::variance_computation(int variable)

{
  if (mean[variable] != D_INF) {
    covariance[variable][variable] = 0.;

    if (nb_vector > 1) {
      register int i;
      double diff;


      for (i = 0;i < nb_vector;i++) {
        diff = vector[i][variable] - mean[variable];
        covariance[variable][variable] += diff * diff;
      }

      covariance[variable][variable] /= (nb_vector - 1);
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la matrice de variance-covariance.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void Vectors::covariance_computation(int variable)

{
  if (mean[0] != D_INF) {
    register int i , j , k;


    for (i = 0;i < nb_variable;i++) {
      if ((variable == I_DEFAULT) || (i == variable)) {
        for (j = (variable == I_DEFAULT ? i + 1 : 0);j < nb_variable;j++) {
          covariance[i][j] = 0.;

          if (nb_vector > 1) {
            for (k = 0;k < nb_vector;k++) {
              covariance[i][j] += (vector[k][i] - mean[i]) * (vector[k][j] - mean[j]);
            }
            covariance[i][j] /= (nb_vector - 1);
          }

          covariance[j][i] = covariance[i][j];
        }
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de l'ecart absolu moyen d'une variable.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

double Vectors::mean_absolute_deviation_computation(int variable) const

{
  register int i;
  double mean_absolute_deviation = D_DEFAULT;


  if (mean[variable] != D_INF) {
    mean_absolute_deviation = 0.;

    for (i = 0;i < nb_vector;i++) {
      mean_absolute_deviation += fabs(vector[i][variable] - mean[variable]);
    }

    mean_absolute_deviation /= nb_vector;
  }

  return mean_absolute_deviation;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la difference absolue moyenne pour une variable.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

double Vectors::mean_absolute_difference_computation(int variable) const

{
  register int i , j;
  double mean_absolute_difference;


  mean_absolute_difference = 0.;

  if (nb_vector > 1) {
    for (i = 0;i < nb_vector;i++) {
      for (j = i + 1;j < nb_vector;j++) {
        mean_absolute_difference += abs(vector[i][variable] - vector[j][variable]);
      }
    }

    mean_absolute_difference = 2 * mean_absolute_difference / (nb_vector * (nb_vector - 1));
  }

  return mean_absolute_difference;
}


/*--------------------------------------------------------------*
 *
 *  Calcul du coefficient d'asymetrie d'une variable.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

double Vectors::skewness_computation(int variable) const

{
  register int i;
  double skewness = D_INF , diff;


  if ((mean[variable] != D_INF) && (covariance[variable][variable] != D_DEFAULT)) {
    skewness = 0.;

    if ((nb_vector > 2) && (covariance[variable][variable] > 0.)) {
      for (i = 0;i < nb_vector;i++) {
        diff = vector[i][variable] - mean[variable];
        skewness += diff * diff * diff;
      }

      skewness = skewness * nb_vector / ((nb_vector - 1) * (nb_vector - 2) *
                  pow(covariance[variable][variable] , 1.5));
    }
  }

  return skewness;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de l'exces d'applatissement d'une variable :
 *  exces d'applatissement = coefficient d'applatissement - 3.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

double Vectors::kurtosis_computation(int variable) const

{
  register int i;
  double kurtosis = D_INF , diff;


  if ((mean[variable] != D_INF) && (covariance[variable][variable] != D_DEFAULT)) {
    if (covariance[variable][variable] == 0.) {
      kurtosis = -2.;
    }

    else {
      kurtosis = 0.;

      for (i = 0;i < nb_vector;i++) {
        diff = vector[i][variable] - mean[variable];
        kurtosis += diff * diff * diff * diff;
      }

      kurtosis = kurtosis / ((nb_vector - 1) * covariance[variable][variable] *
                  covariance[variable][variable]) - 3.;
    }
  }

  return kurtosis;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la matrice des coefficients de correlation a partir
 *  de la matrice de variance-covariance.
 *
 *--------------------------------------------------------------*/

double** Vectors::correlation_computation() const

{
  register int i , j;
  double norm , **correlation = 0;


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
