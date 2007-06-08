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
#include "stat_tool/stat_tools.h"
#include "stat_tool/distribution.h"
#include "stat_tool/vectors.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/stat_label.h"
#include "renewal.h"
#include "sequences.h"
#include "tops.h"
#include "sequence_label.h"

using namespace std;


extern bool identifier_checking(Format_error &error , int nb_pattern , int *pattern_identifier);
extern int* identifier_select(int nb_pattern , int *pattern_identifier , int selected_nb_pattern ,
                              int *selected_identifier , bool keep);
extern int* select_variable(int nb_variable , int selected_nb_variable ,
                            int *selected_variable , bool keep);



/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe Sequences.
 *
 *--------------------------------------------------------------*/

Sequences::Sequences()

{
  nb_variable = 0;

  type = 0;
  min_value = 0;
  max_value = 0;
  marginal = 0;

  nb_sequence = 0;
  identifier = 0;

  max_length = 0;
  cumul_length = 0;
  length = 0;
  hlength = 0;

  index_interval = 0;

  sequence = 0;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Sequences.
 *
 *  arguments : nombre de variables, nombre de sequences.
 *
 *--------------------------------------------------------------*/

Sequences::Sequences(int inb_variable , int inb_sequence)

{
  register int i , j;


  nb_variable = inb_variable;

  type = new int[nb_variable];
  min_value = new int[nb_variable];
  max_value = new int[nb_variable];
  marginal = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    type[i] = INT_VALUE;
    min_value[i] = 0;
    max_value[i] = 0;
    marginal[i] = 0;
  }

  nb_sequence = inb_sequence;

  identifier = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    identifier[i] = i + 1;
  }

  length = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    length[i] = 0;
  }

  max_length = 0;
  cumul_length = 0;
  hlength = 0;

  index_interval = 0;

  sequence = new int**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    sequence[i] = new int*[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      sequence[i][j] = 0;
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Initialisation d'un objet Sequences.
 *
 *  arguments : nombre de variables, types de chaque variable,
 *              nombre de sequences, identificateurs des sequences,
 *              longueurs des sequences, flag initialisation.
 *
 *--------------------------------------------------------------*/

void Sequences::init(int inb_variable , int *itype , int inb_sequence ,
                     int *iidentifier , int *ilength , bool init_flag)

{
  register int i , j , k;
  int blength , *psequence;


  nb_variable = inb_variable;

  type = new int[nb_variable];
  min_value = new int[nb_variable];
  max_value = new int[nb_variable];
  marginal = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    type[i] = itype[i];
    min_value[i] = 0;
    max_value[i] = 0;
    marginal[i] = 0;
  }

  nb_sequence = inb_sequence;

  identifier = new int[nb_sequence];
  if (iidentifier) {
    for (i = 0;i < nb_sequence;i++) {
      identifier[i] = iidentifier[i];
    }
  }
  else {
    for (i = 0;i < nb_sequence;i++) {
      identifier[i] = i + 1;
    }
  }

  length = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    length[i] = ilength[i];
  }

  max_length_computation();
  cumul_length_computation();
  build_length_histogram();

  index_interval = 0;

  sequence = new int**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    sequence[i] = new int*[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      blength = ((type[j] == POSITION) || (type[j] == POSITION_INTERVAL) ? length[i] + 1 : length[i]);
      sequence[i][j] = new int[blength];

      if (init_flag) {
        psequence = sequence[i][j];
        for (k = 0;k < blength;k++) {
          *psequence++ = 0;
        }
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Sequences.
 *
 *  arguments : nombre de variables, types de chaque variable,
 *              nombre de sequences, longueurs des sequences,
 *              sequences, identificateurs des sequences.
 *
 *--------------------------------------------------------------*/

Sequences::Sequences(int inb_variable , int *itype , int inb_sequence ,
                     int *ilength , int ***isequence , int *iidentifier)

{
  register int i , j , k;
  int *psequence , *csequence;


  init(inb_variable , itype , inb_sequence , iidentifier , ilength , false);

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < nb_variable;j++) {
      psequence = sequence[i][j];
      csequence = isequence[i][j];
      for (k = 0;k < (type[j] == POSITION ? length[i] + 1 : length[i]);k++) {
        *psequence++ = *csequence++;
      }
    }
  }

  for (i = 0;i < nb_variable;i++) {
    if (type[i] != POSITION) {
      min_value_computation(i);
      max_value_computation(i);
      build_marginal_histogram(i);
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Initialisation d'un objet Sequences.
 *
 *  arguments : nombre de variables, nombre de sequences, identificateurs des sequences,
 *              longueurs des sequences, flag initialisation.
 *
 *--------------------------------------------------------------*/

void Sequences::init(int inb_variable , int inb_sequence , int *iidentifier ,
                     int *ilength , bool init_flag)

{
  register int i , j , k;
  int *psequence;


  nb_variable = inb_variable;

  type = new int[nb_variable];
  min_value = new int[nb_variable];
  max_value = new int[nb_variable];
  marginal = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    type[i] = INT_VALUE;
    min_value[i] = 0;
    max_value[i] = 0;
    marginal[i] = 0;
  }

  nb_sequence = inb_sequence;

  identifier = new int[nb_sequence];
  if (iidentifier) {
    for (i = 0;i < nb_sequence;i++) {
      identifier[i] = iidentifier[i];
    }
  }
  else {
    for (i = 0;i < nb_sequence;i++) {
      identifier[i] = i + 1;
    }
  }

  length = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    length[i] = ilength[i];
  }

  max_length_computation();
  cumul_length_computation();
  build_length_histogram();

  index_interval = 0;

  sequence = new int**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    sequence[i] = new int*[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      sequence[i][j] = new int[length[i]];

      if (init_flag) {
        psequence = sequence[i][j];
        for (k = 0;k < length[i];k++) {
          *psequence++ = 0;
        }
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Sequences.
 *
 *  arguments : nombre de variables, nombre de sequences, longueurs des sequences,
 *              sequences, identificateurs des sequences.
 *
 *--------------------------------------------------------------*/

Sequences::Sequences(int inb_variable , int inb_sequence , int *ilength ,
                     int ***isequence , int *iidentifier)

{
  register int i , j , k;
  int *psequence , *csequence;


  init(inb_variable , inb_sequence , iidentifier , ilength , false);

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < nb_variable;j++) {
      psequence = sequence[i][j];
      csequence = isequence[i][j];
      for (k = 0;k < length[i];k++) {
        *psequence++ = *csequence++;
      }
    }
  }

  for (i = 0;i < nb_variable;i++) {
    min_value_computation(i);
    max_value_computation(i);
    build_marginal_histogram(i);
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Sequences.
 *
 *  arguments : nombre de variables, histogramme des longueurs des sequences,
 *              flag initialisation.
 *
 *--------------------------------------------------------------*/

Sequences::Sequences(int inb_variable , const Histogram &ihlength , bool init_flag)

{
  register int i , j , k;
  int *plength , *psequence;


  nb_variable = inb_variable;

  type = new int[nb_variable];
  min_value = new int[nb_variable];
  max_value = new int[nb_variable];
  marginal = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    type[i] = INT_VALUE;
    min_value[i] = 0;
    max_value[i] = 0;
    marginal[i] = 0;
  }

  nb_sequence = ihlength.nb_element;

  identifier = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    identifier[i] = i + 1;
  }

  length = new int[nb_sequence];
  plength = length;
  for (i = ihlength.offset;i < ihlength.nb_value;i++) {
    for (j = 0;j < ihlength.frequency[i];j++) {
      *plength++ = i;
    }
  }

  max_length = ihlength.nb_value - 1;
  cumul_length_computation();
  hlength = new Histogram(ihlength);

  index_interval = 0;

  sequence = new int**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    sequence[i] = new int*[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      sequence[i][j] = new int[length[i]];

      if (init_flag) {
        psequence = sequence[i][j];
        for (k = 0;k < length[i];k++) {
          *psequence++ = 0;
        }
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Sequences a partir d'un objet Renewal_data.
 *
 *  argument : reference sur un objet Renewal_data.
 *
 *--------------------------------------------------------------*/

Sequences::Sequences(const Renewal_data &timev)

{
  register int i , j;
  int *psequence , *csequence;


  nb_variable = 1;

  type = new int[nb_variable];
  type[0] = INT_VALUE;

  min_value = new int[nb_variable];
  max_value = new int[nb_variable];
  marginal = new Histogram*[nb_variable];

  nb_sequence = timev.nb_element;

  identifier = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    identifier[i] = i + 1;
  }

  length = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    length[i] = timev.length[i];
  }

  max_length_computation();
  cumul_length_computation();
  build_length_histogram();

  index_interval = 0;

  sequence = new int**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    sequence[i] = new int*[nb_variable];
    sequence[i][0] = new int[length[i]];

    psequence = sequence[i][0];
    csequence = timev.sequence[i];
    for (j = 0;j < length[i];j++) {
      *psequence++ = *csequence++;
    }
  }

  min_value_computation(0);
  max_value_computation(0);
  build_marginal_histogram(0);
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Sequences.
 *
 *  arguments : reference sur un objet Sequences, nombre de sequences,
 *              indices des sequences selectionnees.
 *
 *--------------------------------------------------------------*/

Sequences::Sequences(const Sequences &seq , int inb_sequence , int *index)

{
  register int i , j , k;
  int blength , *psequence , *csequence;


  nb_variable = seq.nb_variable;

  type = new int[nb_variable];
  min_value = new int[nb_variable];
  max_value = new int[nb_variable];
  marginal = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    type[i] = seq.type[i];
    min_value[i] = 0;
    max_value[i] = 0;
    marginal[i] = 0;
  }

  nb_sequence = inb_sequence;

  identifier = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    identifier[i] = seq.identifier[index[i]];
  }

  length = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    length[i] = seq.length[index[i]];
  }

  max_length_computation();
  cumul_length_computation();
  build_length_histogram();

  index_interval = 0;

  sequence = new int**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    sequence[i] = new int*[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      blength = (type[j] == POSITION ? length[i] + 1 : length[i]);
      sequence[i][j] = new int[blength];

      psequence = sequence[i][j];
      csequence = seq.sequence[index[i]][j];
      for (k = 0;k < blength;k++) {
        *psequence++ = *csequence++;
      }
    }
  }

  for (i = 0;i < nb_variable;i++) {
    if (type[i] != POSITION) {
      min_value_computation(i);
      max_value_computation(i);
      build_marginal_histogram(i);
    }
  }

  if ((type[0] == TIME) || ((type[0] == POSITION) && (type[1] == INT_VALUE))) {
    index_interval_computation();
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Sequences.
 *
 *  arguments : reference sur un objet Sequences, flag inversion.
 *
 *--------------------------------------------------------------*/

void Sequences::copy(const Sequences &seq , bool reverse_flag)

{
  register int i , j , k;
  int blength , end_position , *psequence , *csequence;


  nb_variable = seq.nb_variable;

  type = new int[nb_variable];
  min_value = new int[nb_variable];
  max_value = new int[nb_variable];
  marginal = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    type[i] = seq.type[i];
    min_value[i] = seq.min_value[i];
    max_value[i] = seq.max_value[i];

    if (seq.marginal[i]) {
      marginal[i] = new Histogram(*(seq.marginal[i]));
    }
    else {
      marginal[i] = 0;
    }
  }

  nb_sequence = seq.nb_sequence;

  identifier = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    identifier[i] = seq.identifier[i];
  }

  max_length = seq.max_length;
  cumul_length = seq.cumul_length;

  length = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    length[i] = seq.length[i];
  }

  hlength = new Histogram(*(seq.hlength));

  if (seq.index_interval) {
    index_interval = new Histogram(*(seq.index_interval));
  }
  else {
    index_interval = 0;
  }

  sequence = new int**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    sequence[i] = new int*[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      blength = (type[j] == POSITION ? length[i] + 1 : length[i]);
      sequence[i][j] = new int[blength];
      psequence = sequence[i][j];

      switch (reverse_flag) {

      case false : {
        csequence = seq.sequence[i][j];
        for (k = 0;k < blength;k++) {
          *psequence++ = *csequence++;
        }
        break;
      }

      case true : {
        csequence = seq.sequence[i][j] + length[i];

        if (type[j] == POSITION) {
          end_position = *csequence;
          for (k = 0;k < length[i];k++) {
            *psequence++ = end_position - *--csequence;
          }
          *psequence = end_position;
        }

        else {
          for (k = 0;k < length[i];k++) {
            *psequence++ = *--csequence;
          }
        }
        break;
      }
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Sequences avec ajout d'une variable.
 *
 *  arguments : reference sur un objet Sequences, indice de la variable.
 *
 *--------------------------------------------------------------*/

void Sequences::add_variable(const Sequences &seq , int variable)

{
  register int i , j , k , m;
  int blength , *psequence , *csequence;


  nb_variable = seq.nb_variable + 1;

  type = new int[nb_variable];
  min_value = new int[nb_variable];
  max_value = new int[nb_variable];
  marginal = new Histogram*[nb_variable];

  i = 0;
  for (j = 0;j < nb_variable;j++) {
    if (j != variable) {
      type[j] = seq.type[i];
      min_value[j] = seq.min_value[i];
      max_value[j] = seq.max_value[i];

      if (seq.marginal[i]) {
        marginal[j] = new Histogram(*(seq.marginal[i]));
      }
      else {
        marginal[j] = 0;
      }

      i++;
    }

    else {
      type[j] = INT_VALUE;
      min_value[j] = 0;
      max_value[j] = 0;
      marginal[j] = 0;
    }
  }

  nb_sequence = seq.nb_sequence;

  identifier = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    identifier[i] = seq.identifier[i];
  }

  max_length = seq.max_length;
  cumul_length = seq.cumul_length;

  length = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    length[i] = seq.length[i];
  }

  hlength = new Histogram(*(seq.hlength));

  if (seq.index_interval) {
    index_interval = new Histogram(*(seq.index_interval));
  }
  else {
    index_interval = 0;
  }

  sequence = new int**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    sequence[i] = new int*[nb_variable];
    j = 0;
    for (k = 0;k < nb_variable;k++) {
      blength = (type[k] == POSITION ? length[i] + 1 : length[i]);
      sequence[i][k] = new int[blength];

      if (k != variable) {
        psequence = sequence[i][k];
        csequence = seq.sequence[i][j++];
        for (m = 0;m < blength;m++) {
          *psequence++ = *csequence++;
        }
      }

      else {
        psequence = sequence[i][k];
        for (m = 0;m < blength;m++) {
          *psequence++ = 0;
        }
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par copie de la classe Sequences.
 *
 *  arguments : reference sur un objet Sequences, type de transformation
 *              ('c' : copie, 'a' : addition d'une variable),
 *              inversion / indice de la variable ('a').
 *
 *--------------------------------------------------------------*/

Sequences::Sequences(const Sequences &seq , char transform , int param)

{
  switch (transform) {
  case 'c' :
    Sequences::copy(seq , (param == REVERSE ? true : false));
    break;
  case 'a' :
    Sequences::add_variable(seq , param);
    break;
  default :
    Sequences::copy(seq);
    break;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet Sequences.
 *
 *--------------------------------------------------------------*/

void Sequences::remove()

{
  register int i , j;


  delete [] type;
  delete [] min_value;
  delete [] max_value;

  if (marginal) {
    for (i = 0;i < nb_variable;i++) {
      delete marginal[i];
    }
    delete [] marginal;
  }

  delete [] identifier;

  delete [] length;
  delete hlength;

  delete index_interval;

  if (sequence) {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < nb_variable;j++) {
        delete [] sequence[i][j];
      }
      delete [] sequence[i];
    }
    delete [] sequence;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Sequences.
 *
 *--------------------------------------------------------------*/

Sequences::~Sequences()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Sequences.
 *
 *  argument : reference sur un objet Sequences.
 *
 *--------------------------------------------------------------*/

Sequences& Sequences::operator=(const Sequences &seq)

{
  if (&seq != this) {
    remove();
    copy(seq);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Extraction d'un histogramme du nombre d'occurrences des valeurs.
 *
 *  arguments : reference sur un objet Format_error, variable.
 *
 *--------------------------------------------------------------*/

Distribution_data* Sequences::extract(Format_error &error , int variable) const

{
  bool status = true;
  Distribution_data *histo;


  histo = 0;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if ((type[variable] != INT_VALUE) && (type[variable] != STATE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                    << SEQ_error[SEQR_VARIABLE_TYPE];
      correction_message << STAT_sequence_word[INT_VALUE] << " or " << STAT_sequence_word[STATE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }

    else if (!marginal[variable]) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                    << STAT_error[STATR_MARGINAL_HISTOGRAM];
      error.update((error_message.str()).c_str());
    }
  }

  if (status) {
    histo = new Distribution_data(*(marginal[variable]));
  }

  return histo;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un object Vectors a partir d'un objet Sequences.
 *
 *  argument : flag variable d'index.
 *
 *--------------------------------------------------------------*/

Vectors* Sequences::build_vectors(bool index_variable) const

{
  register int i , j , k , m;
  int offset;
  Vectors *vec;


  if ((type[0] == TIME) || (type[0] == POSITION)) {
    index_variable = false;
  }

  switch (index_variable) {
  case false :
    offset = 0;
    break;
  case true :
    offset = 1;
    break;
  }

  vec = new Vectors(nb_variable + offset , cumul_length , 0 , false);

  i = 0;
  for (j = 0;j < nb_sequence;j++) {
    for (k = 0;k < length[j];k++) {
      if (index_variable) {
        vec->vector[i][0] = k;
      }
      for (m = 0;m < nb_variable;m++) {
        vec->vector[i][m + offset] = sequence[j][m][k];
      }
      i++;
    }
  }

  if (index_variable) {
    vec->min_value_computation(0);
    vec->max_value_computation(0);
    vec->build_marginal_histogram(0);
  }

  for (i = 0;i < nb_variable;i++) {
    if (type[i] != POSITION) {
      vec->min_value[i + offset] = min_value[i];
      vec->max_value[i + offset] = max_value[i];
    }
    else {
      vec->min_value_computation(i + offset);
      vec->max_value_computation(i + offset);
    }

    if (marginal[i]) {
      vec->marginal[i + offset] = new Histogram(*marginal[i]);
      vec->mean[i + offset] = vec->marginal[i + offset]->mean;
      vec->covariance[i + offset][i + offset] = vec->marginal[i + offset]->variance;
    }
    else {
      vec->build_marginal_histogram(i + offset);
    }
  }

  vec->covariance_computation();

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Extraction de mesures globales (longueur, nombre de series/d'occurrences d'une valeur)
 *  par sequence.
 *
 *  arguments : reference sur un objet Format_error, type, variable, valeur.
 *
 *--------------------------------------------------------------*/

Vectors* Sequences::extract_vectors(Format_error &error , int type , int variable , int value) const

{
  bool status = true;
  register int i , j;
  int count , *psequence;
  Vectors *vec;


  vec = 0;
  error.init();

  if (variable != I_DEFAULT) {
    if ((variable < 1) || (variable > nb_variable)) {
      status = false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
    }

    else {
      variable--;

      if ((value < min_value[variable]) || (value > max_value[variable]) ||
          ((marginal[variable]) && (marginal[variable]->frequency[value] == 0))) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VALUE] << " " << value << " "
                      << SEQ_error[SEQR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }
    }
  }

  if (status) {
    vec = new Vectors(1 , nb_sequence , identifier , false);

    switch (type) {

    case LENGTH : {
      for (i = 0;i < nb_sequence;i++) {
        vec->vector[i][0] = length[i];
      }
      break;
    }

    case NB_RUN : {
      for (i = 0;i < nb_sequence;i++) {
        psequence = sequence[i][variable];
        count = 0;
        if (*psequence++ == value) {
          count++;
        }
        for (j = 1;j < length[i];j++) {
          if ((*psequence != *(psequence - 1)) && (*psequence == value)) {
            count++;
          }
          psequence++;
        }

        vec->vector[i][0] = count;
      }
      break;
    }

    case NB_OCCURRENCE : {
      for (i = 0;i < nb_sequence;i++) {
        psequence = sequence[i][variable];
        count = 0;
        for (j = 0;j < length[i];j++) {
          if (*psequence++ == value) {
            count++;
          }
        }

        vec->vector[i][0] = count;
      }
      break;
    }
    }

    vec->min_value_computation(0);
    vec->max_value_computation(0);
    vec->build_marginal_histogram(0);
  }

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Markovian_sequences a partir d'un objet Sequences.
 *
 *  argument : reference sur un objet Format_error.
 *
 *--------------------------------------------------------------*/

Markovian_sequences* Sequences::markovian_sequences(Format_error &error) const

{
  bool status = true;
  register int i;
  Markovian_sequences *seq;


  seq = 0;
  error.init();

  for (i = 0;i < nb_variable;i++) {
    if ((type[i] != INT_VALUE) && (type[i] != STATE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << SEQ_error[SEQR_VARIABLE_TYPE];
      correction_message << STAT_sequence_word[INT_VALUE] << " or " << STAT_sequence_word[STATE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }

    else {
      if (min_value[i] < 0) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                      << SEQ_error[SEQR_POSITIVE_MIN_VALUE];
        error.update((error_message.str()).c_str());
      }

      if (max_value[i] == min_value[i]) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                      << STAT_error[STATR_NB_VALUE];
        error.update((error_message.str()).c_str());
      }

      if (!marginal[i]) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                      << STAT_error[STATR_MARGINAL_HISTOGRAM];
        error.update((error_message.str()).c_str());
      }
    }
  }

  if (status) {
    seq = new Markovian_sequences(*this);
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Tops a partir d'un objet Sequences.
 *
 *  argument : reference sur un objet Format_error.
 *
 *--------------------------------------------------------------*/

Tops* Sequences::tops(Format_error &error) const

{
  bool status = true;
  register int i;
  Tops *tops;


  tops = 0;
  error.init();

  if (nb_variable != 2) {
    status = false;
    error.correction_update(STAT_error[STATR_NB_VARIABLE] , 2);
  }

  if (type[0] != POSITION) {
    status = false;
    error.correction_update(SEQ_error[SEQR_VARIABLE_1_TYPE] , STAT_sequence_word[POSITION]);
  }

  if (nb_variable >= 2) {
    if (type[1] != NB_INTERNODE) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << 2 << ": "
                    << SEQ_error[SEQR_VARIABLE_TYPE];
      error.correction_update((error_message.str()).c_str() , STAT_sequence_word[NB_INTERNODE]);
    }

    else if (!marginal[1]) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << 2 << ": "
                    << STAT_error[STATR_MARGINAL_HISTOGRAM];
      error.update((error_message.str()).c_str());
    }
  }

  if (status) {
    for (i = 0;i < nb_sequence;i++) {
      if (sequence[i][0][0] == 0) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_SEQUENCE] << " " << i + 1 << ": "
                      << SEQ_parsing[SEQP_POSITION];
        error.update((error_message.str()).c_str());
      }
    }
  }

  if (status) {
    tops = new Tops(*this);
  }

  return tops;
}


/*--------------------------------------------------------------*
 *
 *  Verification du caractere (strictement) croissant des sequences
 *  pour une variable donnee.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable,
 *              flag croissance strict ou non, labels de l'objet et de la variable.
 *
 *--------------------------------------------------------------*/

bool Sequences::increasing_sequence_checking(Format_error &error , int variable , int strict ,
                                             const char *pattern_label , const char *variable_label) const

{
  bool status = true;
  register int i , j;
  int *psequence;


  for (i = 0;i < nb_sequence;i++) {
    psequence = sequence[i][variable] + 1;

    for (j = 1;j < (type[variable] == POSITION ? length[i] + 1 : length[i]);j++) {
      if ((((!strict) || (j == length[i])) && (*psequence < *(psequence - 1))) ||
          ((strict) && (j < length[i]) && (*psequence <= *(psequence - 1)))) {
        status = false;
        ostringstream error_message;

        if (type[variable] == INT_VALUE) {
          error_message << pattern_label << " " << i + 1 << ": " << STAT_label[STATL_VARIABLE] << " "
                        << variable + 1 << ": " << variable_label << " " << *psequence << " "
                        << STAT_error[STATR_NOT_ALLOWED];
        }

        else {
          error_message << pattern_label << " " << i + 1 << ": " << variable_label << " "
                        << *psequence << " " << STAT_error[STATR_NOT_ALLOWED];
        }

        error.update((error_message.str()).c_str());
      }

      psequence++;
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Verification d'un objet Sequences.
 *
 *  arguments : reference sur un objet Format_error, label.
 *
 *--------------------------------------------------------------*/

bool Sequences::check(Format_error &error , const char *pattern_label)

{
  bool status = true , lstatus;


  error.init();

  if (nb_variable > SEQUENCE_NB_VARIABLE) {
    status = false;
    error.update(STAT_error[STATR_NB_VARIABLE]);
  }

  if (max_length <= 1) {
    status = false;
    error.update(SEQ_parsing[SEQP_MAX_SEQUENCE_LENGTH]);
  }

  lstatus = identifier_checking(error , nb_sequence , identifier);
  if (!lstatus) {
    status = false;
  }

  if ((type[0] == TIME) || (type[0] == POSITION)) {
    lstatus = increasing_sequence_checking(error , 0 , (type[0] == POSITION ? false : true) ,
                                           pattern_label , STAT_sequence_word[type[0]]);

    if (!lstatus) {
      status = false;
    }
  }

  if ((status) && ((type[0] == TIME) || ((type[0] == POSITION) && (type[1] == INT_VALUE)))) {
    index_interval_computation();
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Extraction d'un objet Time_events a partir d'un objet Sequences.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable,
 *              nombre d'evenements, dates de debut et de fin,
 *              dates precedente et suivante.
 *
 *--------------------------------------------------------------*/

Time_events* Sequences::extract_time_events(Format_error &error , int variable ,
                                            int begin_date , int end_date ,
                                            int previous_date , int next_date) const

{
  bool status = true , lstatus;
  register int i , j;
  int nb_element , previous , begin , end , next , *time , *nb_event , *pdate;
  Time_events *timev;


  timev = 0;
  error.init();

  if (type[0] != TIME) {
    status = false;
    error.correction_update(SEQ_error[SEQR_VARIABLE_1_TYPE] , STAT_sequence_word[TIME]);
  }

  if ((variable < 2) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    lstatus = increasing_sequence_checking(error , variable , false , SEQ_label[SEQL_SEQUENCE] ,
                                           STAT_label[STATL_VALUE]);
    if (!lstatus) {
      status = false;
    }
  }

  if (begin_date >= end_date) {
    status = false;
    error.update(SEQ_error[SEQR_DATE_ORDER]);
  }

  if (previous_date != I_DEFAULT) {
    if (previous_date > begin_date) {
      status = false;
      error.update(SEQ_error[SEQR_DATE_ORDER]);
    }
  }
  else {
    previous_date = begin_date;
  }

  if (next_date != I_DEFAULT) {
    if (next_date < end_date) {
      status = false;
      error.update(SEQ_error[SEQR_DATE_ORDER]);
    }
  }
  else {
    next_date = end_date;
  }

  if (status) {
    time = new int[nb_sequence];
    nb_event = new int[nb_sequence];
    nb_element = 0;

    for (i = 0;i < nb_sequence;i++) {
      pdate = sequence[i][0];
      previous = I_DEFAULT;
      begin = I_DEFAULT;
      end = I_DEFAULT;
      next = I_DEFAULT;

      for (j = 0;j < length[i];j++) {
        if (*pdate == previous_date) {
          previous = j;
        }
        if (*pdate == begin_date) {
          begin = j;
        }
        if (*pdate == end_date) {
          end = j;
        }
        if (*pdate == next_date) {
          next = j;
          break;
        }
        pdate++;
      }

      if ((previous != I_DEFAULT) && (begin != I_DEFAULT) && (end != I_DEFAULT) &&
          (next != I_DEFAULT) && ((previous == begin) || ((previous <  begin) &&
            (sequence[i][variable][previous] < sequence[i][variable][begin]))) && ((end == next) ||
           ((end < next) && (sequence[i][variable][end] < sequence[i][variable][next])))) {
        time[nb_element] = end_date - begin_date;
        nb_event[nb_element++] = sequence[i][variable][end] - sequence[i][variable][begin];
      }
    }

    if (nb_element == 0) {
      status = false;
      error.update(STAT_error[STATR_EMPTY_SAMPLE]);
    }

    else {
      timev = new Time_events(nb_element , time , nb_event);
    }

    delete [] time;
    delete [] nb_event;
  }

  return timev;
}


/*--------------------------------------------------------------*
 *
 *  Extraction d'un objet Renewal_data a partir d'un objet Sequences.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable,
 *              nombre d'evenements, index de debut et de fin.
 *
 *--------------------------------------------------------------*/

Renewal_data* Sequences::extract_renewal_data(Format_error &error , int variable ,
                                              int begin_index , int end_index) const

{
  bool status = true , lstatus;
  register int i , j;
  int nb_element , index , *ptime , *pnb_event , *psequence , *csequence;
  Renewal_data *timev;


  timev = 0;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if (type[variable] != INT_VALUE) {
      status = false;
      error.correction_update(SEQ_error[SEQR_VARIABLE_TYPE] , STAT_sequence_word[INT_VALUE]);
    }

    else {
      lstatus = increasing_sequence_checking(error , variable , false , SEQ_label[SEQL_SEQUENCE] ,
                                             STAT_label[STATL_VALUE]);
      if (!lstatus) {
        status = false;
      }
    }
  }

  if ((begin_index < 0) || (begin_index + 1 >= max_length) || (begin_index > end_index)) {
    status = false;
    error.update(SEQ_error[SEQR_BEGIN_INDEX]);
  }
  if ((end_index < 0) || (end_index + 1 >= max_length) || (end_index < begin_index)) {
    status = false;
    error.update(SEQ_error[SEQR_END_INDEX]);
  }

  if (status) {
    timev = new Renewal_data(nb_sequence , end_index + 1 - begin_index);

    ptime = new int[nb_sequence];
    pnb_event = new int[nb_sequence];

    nb_element = 0;
    for (i = 0;i < nb_sequence;i++) {
      if (end_index + 1 < length[i]) {
        *ptime++ = end_index + 1 - begin_index;
        *pnb_event++ = sequence[i][variable][end_index + 1] -
                       sequence[i][variable][begin_index];

        timev->length[nb_element] = end_index + 1 - begin_index;
        timev->sequence[nb_element] = new int[timev->length[nb_element]];

        psequence = timev->sequence[nb_element++];
        csequence = sequence[i][variable] + begin_index;
        index = begin_index;
        for (j = begin_index + 1;j <= end_index + 1;j++) {
          *psequence = *(csequence + 1) - *csequence;
          csequence++;

          if (*psequence > 0) {
            if (index == begin_index) {
              (timev->forward->frequency[j - index])++;
             }
             else {
               (timev->within->frequency[j - index])++;
             }
             index = j;
          }
          psequence++;
        }

        if (index > begin_index) {
          (timev->backward->frequency[end_index + 1 - index])++;
        }
      }
    }

    // construction des echantillons {temps, nombre d'evenements, frequence} et
    // des histogrammes du temps d'observation et du nombre d'evenements

    ptime -= nb_element;
    pnb_event -= nb_element;

    timev->build(nb_element , ptime , pnb_event);
    delete [] ptime;
    delete [] pnb_event;

    // extraction des caracteristiques des histogrammes des intervalles de temps entre 2 evenements,
    // des intervalles de temps entre 2 evenements compris entre les 2 dates d'observation,
    // des intervalles de temps apres le dernier evenement, des intervalles de temps residuel

    timev->within->nb_value_computation();
    timev->within->offset_computation();
    timev->within->nb_element_computation();
    timev->within->max_computation();
    timev->within->mean_computation();
    timev->within->variance_computation();

    timev->backward->nb_value_computation();
    timev->backward->offset_computation();
    timev->backward->nb_element_computation();
    timev->backward->max_computation();
    timev->backward->mean_computation();
    timev->backward->variance_computation();

    timev->forward->nb_value_computation();
    timev->forward->offset_computation();
    timev->forward->nb_element_computation();
    timev->forward->max_computation();
    timev->forward->mean_computation();
    timev->forward->variance_computation();

    timev->build_index_event(1);

    if ((timev->backward->nb_element == 0) && (timev->forward->nb_element == 0)) {
      delete timev;
      timev = 0;
      error.update(SEQ_error[SEQR_BOTH_END_CENSORED_INTERVAL]);
    }
  }

  return timev;
}


/*--------------------------------------------------------------*
 *
 *  Fusion d'objets Sequences.
 *
 *  arguments : reference sur un objet Format_error, nombre d'objets Sequences,
 *              pointeurs sur les objets Sequences.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::merge(Format_error &error , int nb_sample , const Sequences **iseq) const

{
  bool status = true;
  register int i , j , k , m , n;
  int inb_sequence , *ilength , *psequence , *csequence;
  const Histogram **phisto;
  Sequences *seq;
  const Sequences **pseq;


  seq = 0;
  error.init();

  for (i = 0;i < nb_sample;i++) {
    if (iseq[i]->nb_variable != nb_variable) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                    << STAT_error[STATR_NB_VARIABLE];
      error.correction_update((error_message.str()).c_str() , nb_variable);
    }

    else {
      for (j = 0;j < nb_variable;j++) {
        if (iseq[i]->type[j] != type[j]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                        << STAT_label[STATL_VARIABLE] << " " << j + 1 << ": "
                        << SEQ_error[SEQR_VARIABLE_TYPE];
          error.correction_update((error_message.str()).c_str() , STAT_sequence_word[type[j]]);
        }
      }
    }
  }

  if (status) {
    nb_sample++;
    pseq = new const Sequences*[nb_sample];

    pseq[0] = this;
    for (i = 1;i < nb_sample;i++) {
      pseq[i] = iseq[i - 1];
    }

    // calcul du nombre et de la longueur des sequences

    inb_sequence = 0;
    for (i = 0;i < nb_sample;i++) {
      inb_sequence += pseq[i]->nb_sequence;
    }

    ilength = new int[inb_sequence];

    i = 0;
    for (j = 0;j < nb_sample;j++) {
      for (k = 0;k < pseq[j]->nb_sequence;k++) {
        ilength[i++] = pseq[j]->length[k];
      }
    }

    seq = new Sequences(nb_variable , type , inb_sequence , 0 , ilength , false);
    delete [] ilength;

    // copie des sequences

    i = 0;
    for (j = 0;j < nb_sample;j++) {
      for (k = 0;k < pseq[j]->nb_sequence;k++) {
        for (m = 0;m < pseq[j]->nb_variable;m++) {
          psequence = seq->sequence[i][m];
          csequence = pseq[j]->sequence[k][m];
          for (n = 0;n < (pseq[j]->type[m] == POSITION ? pseq[j]->length[k] + 1 : pseq[j]->length[k]);n++) {
            *psequence++ = *csequence++;
          }
        }
        i++;
      }
    }

    phisto = new const Histogram*[nb_sample];

    for (i = 0;i < seq->nb_variable;i++) {
      if (seq->type[i] != POSITION) {
        seq->min_value[i] = pseq[0]->min_value[i];
        seq->max_value[i] = pseq[0]->max_value[i];
        for (j = 1;j < nb_sample;j++) {
          if (pseq[j]->min_value[i] < seq->min_value[i]) {
            seq->min_value[i] = pseq[j]->min_value[i];
          }
          if (pseq[j]->max_value[i] > seq->max_value[i]) {
            seq->max_value[i] = pseq[j]->max_value[i];
          }
        }

        for (j = 0;j < nb_sample;j++) {
          if (pseq[j]->marginal[i]) {
            phisto[j] = pseq[j]->marginal[i];
          }
          else {
            break;
          }
        }

        if (j == nb_sample) {
          seq->marginal[i] = new Histogram(nb_sample , phisto);
        }
      }
    }

    if ((seq->type[0] == TIME) || ((seq->type[0] == POSITION) && (seq->type[1] == INT_VALUE))) {
      for (i = 0;i < nb_sample;i++) {
        phisto[i] = pseq[i]->index_interval;
      }
      seq->index_interval = new Histogram(nb_sample , phisto);
    }

    delete [] phisto;
    delete [] pseq;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Translation des valeurs d'une variable donnee.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable,
 *              parametre de translation.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::shift(Format_error &error , int variable , int shift_param) const

{
  bool status = true;
  register int i , j , k;
  int *psequence , *csequence;
  Sequences *seq;


  seq = 0;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if ((type[variable] != INT_VALUE) && (type[variable] != TIME)) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_sequence_word[INT_VALUE] << " or " << STAT_sequence_word[TIME];
      error.correction_update(SEQ_error[SEQR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }

    else {
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
  }

  if (status) {
    seq = new Sequences(nb_variable , type , nb_sequence , identifier , length , false);

    for (i = 0;i < seq->nb_sequence;i++) {
      for (j = 0;j < seq->nb_variable;j++) {
        psequence = seq->sequence[i][j];
        csequence = sequence[i][j];

        // translation des valeurs

        if (j == variable) {
          for (k = 0;k < seq->length[i];k++) {
            *psequence++ = *csequence++ + shift_param;
          }
        }

        // copie des valeurs

        else {
          for (k = 0;k < (seq->type[j] == POSITION ? seq->length[i] + 1 : seq->length[i]);k++) {
            *psequence++ = *csequence++;
          }
        }
      }
    }

    for (i = 0;i < seq->nb_variable;i++) {
      if (i == variable) {
        seq->min_value[i] = min_value[i] + shift_param;
        seq->max_value[i] = max_value[i] + shift_param;

        if ((seq->min_value[i] >= 0) && (seq->max_value[i] <= MARGINAL_MAX_VALUE)) {
          if (marginal[i])  {
            seq->marginal[i] = new Histogram(*(marginal[i]) , 's' , shift_param);
          }
          else {
            seq->build_marginal_histogram(i);
          }
        }
      }

      else if (seq->type[i] != POSITION) {
        seq->min_value[i] = min_value[i];
        seq->max_value[i] = max_value[i];
        if (marginal[i]) {
          seq->marginal[i] = new Histogram(*(marginal[i]));
        }
      }
    }

    if ((seq->type[0] == TIME) || ((seq->type[0] == POSITION) && (seq->type[1] == INT_VALUE))) {
      seq->index_interval = new Histogram(*(index_interval));
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Regroupement des valeurs d'une variable donnee.
 *
 *  arguments : reference sur un objet Sequences, indice de la variable,
 *              pas de regroupement, flag pour ajouter une variable.
 *
 *--------------------------------------------------------------*/

void Sequences::cluster(const Sequences &seq , int ivariable , int step , bool add_flag)

{
  register int i , j , k;
  int variable , offset , *psequence , *csequence;


  switch (add_flag) {
  case false :
    variable = ivariable;
    offset = 0;
    break;
  case true :
    variable = 0;
    offset = 1;
    break;
  }

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < nb_variable;j++) {
      psequence = sequence[i][j];

      // regroupement des valeurs

      if (j == variable) {
        csequence = seq.sequence[i][ivariable];
        for (k = 0;k < length[i];k++) {
          *psequence++ = *csequence++ / step;
        }
      }

      // copie des valeurs

      else {
        csequence = seq.sequence[i][j - offset];
        for (k = 0;k < (type[j] == POSITION ? length[i] + 1 : length[i]);k++) {
          *psequence++ = *csequence++;
        }
      }
    }
  }

  for (i = 0;i < nb_variable;i++) {
    if (i == variable) {
      min_value[i] = seq.min_value[i - offset] / step;
      max_value[i] = seq.max_value[i - offset] / step;
      if (seq.marginal[i - offset]) {
        marginal[i] = new Histogram(*(seq.marginal[i - offset]) , 'c' , step);
      }
      else {
        build_marginal_histogram(i);
      }
    }

    else if (type[i] != POSITION) {
      min_value[i] = seq.min_value[i - offset];
      max_value[i] = seq.max_value[i - offset];
      if (seq.marginal[i - offset]) {
        marginal[i] = new Histogram(*(seq.marginal[i - offset]));
      }
    }
  }

  if ((type[0] == TIME) || ((type[0] == POSITION) && (type[1] == INT_VALUE))) {
    index_interval = new Histogram(*(seq.index_interval));
  }
}


/*--------------------------------------------------------------*
 *
 *  Regroupement des valeurs d'une variable donnee.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable,
 *              pas de regroupement.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::cluster(Format_error &error , int variable , int step) const

{
  bool status = true;
  Sequences *seq;


  seq = 0;
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
      correction_message << STAT_sequence_word[INT_VALUE] << " or " << STAT_sequence_word[STATE];
      error.correction_update(SEQ_error[SEQR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }
  }

  if (step < 1) {
    status = false;
    error.update(STAT_error[STATR_CLUSTERING_STEP]);
  }

  if (status) {
    seq = new Sequences(nb_variable , type , nb_sequence , identifier , length , false);
    seq->cluster(*this , variable , step);
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Transcodage des symboles d'une variable donnee.
 *
 *  arguments : reference sur un objet Sequences, indice de la variable,
 *              plus petit et plus grand symboles, table de transcodage des symboles,
 *              flag pour ajouter une variable.
 *
 *--------------------------------------------------------------*/

void Sequences::transcode(const Sequences &seq , int ivariable , int min_symbol ,
                          int max_symbol , int *symbol , bool add_flag)

{
  register int i , j , k;
  int variable , offset , *psequence , *csequence;


  switch (add_flag) {
  case false :
    variable = ivariable;
    offset = 0;
    break;
  case true :
    variable = 0;
    offset = 1;
    break;
  }

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < nb_variable;j++) {
      psequence = sequence[i][j];

      // transcodage des symboles

      if (j == variable) {
        csequence = seq.sequence[i][ivariable];
        for (k = 0;k < length[i];k++) {
          *psequence++ = symbol[*csequence++ - seq.min_value[variable]] + min_symbol;
        }
      }

      // copie des valeurs

      else {
        csequence = seq.sequence[i][j - offset];
        for (k = 0;k < (type[j] == POSITION ? length[i] + 1 : length[i]);k++) {
          *psequence++ = *csequence++;
        }
      }
    }
  }

  for (i = 0;i < nb_variable;i++) {
    if (i == variable) {
      min_value[i] = min_symbol;
      max_value[i] = max_symbol;
      build_marginal_histogram(i);
    }

    else if (type[i] != POSITION) {
      min_value[i] = seq.min_value[i - offset];
      max_value[i] = seq.max_value[i - offset];
      if (seq.marginal[i - offset]) {
        marginal[i] = new Histogram(*(seq.marginal[i - offset]));
      }
    }
  }

  if ((type[0] == TIME) || ((type[0] == POSITION) && (type[1] == INT_VALUE))) {
    index_interval = new Histogram(*(seq.index_interval));
  }
}


/*--------------------------------------------------------------*
 *
 *  Transcodage des symboles d'une variable donnee.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable,
 *              table de transcodage des symboles.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::transcode(Format_error &error , int variable , int *symbol) const

{
  bool status = true , *presence;
  register int i;
  int min_symbol , max_symbol;
  Sequences *seq;


  seq = 0;
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
      correction_message << STAT_sequence_word[INT_VALUE] << " or " << STAT_sequence_word[STATE];
      error.correction_update(SEQ_error[SEQR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }

    else {
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

      seq = new Sequences(nb_variable , type , nb_sequence , identifier , length , false);
      seq->transcode(*this , variable , min_symbol , max_symbol , symbol);
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Regroupement des symboles d'une variable donnee.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable,
 *              nombre de classes, bornes pour regrouper les symboles.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::cluster(Format_error &error , int variable ,
                              int nb_class , int *ilimit) const

{
  bool status = true;
  register int i , j , k;
  int *limit , *symbol;
  Sequences *seq;


  seq = 0;
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
      correction_message << STAT_sequence_word[INT_VALUE] << " or " << STAT_sequence_word[STATE];
      error.correction_update(SEQ_error[SEQR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }
    else if ((nb_class < 2) || (nb_class >= max_value[variable] - min_value[variable] + 1)) {
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

      seq = new Sequences(nb_variable , type , nb_sequence , identifier , length , false);
      seq->transcode(*this , variable , 0 , nb_class - 1 , symbol);

      delete [] symbol;
    }

    delete [] limit;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Selection de sequences sur les valeurs prises pour une variable donnee.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable,
 *              bornes sur les valeurs, flag pour conserver ou rejeter
 *              les sequences selectionnees.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::value_select(Format_error &error , int variable , int imin_value ,
                                   int imax_value , bool keep) const

{
  bool status = true;
  register int i , j;
  int inb_sequence , *index , *psequence;
  Sequences *seq;


  seq = 0;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if (type[variable] != POSITION) {
      if ((imin_value > max_value[variable]) || (imin_value > imax_value)) {
        status = false;
        error.update(STAT_error[STATR_MIN_VALUE]);
      }
      if ((imax_value < min_value[variable]) || (imax_value < imin_value)) {
        status = false;
        error.update(STAT_error[STATR_MAX_VALUE]);
      }
    }
  }

  if (status) {

    // selection des sequences

    index = new int[nb_sequence];
    inb_sequence = 0;

    for (i = 0;i < nb_sequence;i++) {
      psequence = sequence[i][variable];
      for (j = 0;j < (type[variable] == POSITION ? length[i] + 1 : length[i]);j++) {
        if ((*psequence >= imin_value) && (*psequence <= imax_value)) {
          if (keep) {
            index[inb_sequence++] = i;
          }
          break;
        }

        psequence++;
      }

      if ((!keep) && (j == length[i])) {
        index[inb_sequence++] = i;
      }
    }

    if (inb_sequence == 0) {
      status = false;
      error.update(STAT_error[STATR_EMPTY_SAMPLE]);
    }

    // copie des sequences

    if (status) {
      seq = new Sequences(*this , inb_sequence , index);
    }

    delete [] index;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Selection de sequences par l'identificateur.
 *
 *  arguments : reference sur un objet Format_error, nombre de sequences,
 *              identificateur des sequences, flag pour conserver ou rejeter
 *              les sequences selectionnees.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::select_individual(Format_error &error , int inb_sequence ,
                                        int *iidentifier , bool keep) const

{
  bool status = true , *selected_sequence;
  register int i , j;
  int max_identifier , *index;
  Sequences *seq;


  seq = 0;
  error.init();

  if ((inb_sequence < 1) || (inb_sequence > (keep ? nb_sequence : nb_sequence - 1))) {
    status = false;
    error.update(SEQ_error[SEQR_NB_SEQUENCE]);
  }

  else {
    max_identifier = 1;
    for (i = 0;i < inb_sequence;i++) {
      if (iidentifier[i] > max_identifier) {
        max_identifier = iidentifier[i];
      }
    }

    selected_sequence = new bool[max_identifier + 1];
    for (i = 0;i <= max_identifier;i++) {
      selected_sequence[i] = false;
    }

    for (i = 0;i < inb_sequence;i++) {
      for (j = 0;j < nb_sequence;j++) {
        if (iidentifier[i] == identifier[j]) {
          break;
        }
      }

      if (j == nb_sequence) {
        status = false;
        ostringstream error_message;
        error_message << iidentifier[i] << ": " << SEQ_error[SEQR_SEQUENCE_IDENTIFIER];
        error.update((error_message.str()).c_str());
      }

      else if (selected_sequence[iidentifier[i]]) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_SEQUENCE] << " " << iidentifier[i] << " "
                      << STAT_error[STATR_ALREADY_SELECTED];
        error.update((error_message.str()).c_str());
      }
      else {
        selected_sequence[iidentifier[i]] = true;
      }
    }

    delete [] selected_sequence;
  }

  if (status) {
    index = identifier_select(nb_sequence , identifier , inb_sequence , iidentifier , keep);

    seq = new Sequences(*this , (keep ? inb_sequence : nb_sequence - inb_sequence) , index);

    delete [] index;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Selection de variables.
 *
 *  arguments : reference sur un objet Sequences, indices des variables.
 *
 *--------------------------------------------------------------*/

void Sequences::select_variable(const Sequences &seq , int *variable)

{
  register int i , j , k;
  int *psequence , *csequence;


  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < nb_variable;j++) {
      psequence = sequence[i][j];
      csequence = seq.sequence[i][variable[j]];
      for (k = 0;k < (type[j] == POSITION ? length[i] + 1 : length[i]);k++) {
        *psequence++ = *csequence++;
      }
    }
  }

  for (i = 0;i < nb_variable;i++) {
    if (type[i] != POSITION) {
      min_value[i] = seq.min_value[variable[i]];
      max_value[i] = seq.max_value[variable[i]];
      if (seq.marginal[variable[i]]) {
        marginal[i] = new Histogram(*(seq.marginal[variable[i]]));
      }
    }
  }

  if ((type[0] == TIME) || ((type[0] == POSITION) && (type[1] == INT_VALUE))) {
    index_interval = new Histogram(*(seq.index_interval));
  }
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

Sequences* Sequences::select_variable(Format_error &error , int inb_variable ,
                                      int *ivariable , bool keep) const

{
  bool status = true , *selected_variable;
  register int i;
  int bnb_variable , *variable , *itype;
  Sequences *seq;


  seq = 0;
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

    bnb_variable = (keep ? inb_variable : nb_variable - inb_variable);

    itype = new int[bnb_variable];
    for (i = 0;i < bnb_variable;i++) {
      itype[i] = type[variable[i]];
    }

    if ((bnb_variable == 1) && ((itype[0] == TIME) || (itype[0] == POSITION))) {
      status = false;
      error.update(STAT_error[STATR_NB_SELECTED_VARIABLE]);
    }

    if (status) {
      seq = new Sequences(bnb_variable , itype , nb_sequence , identifier , length , false);
      seq->select_variable(*this , variable);
    }

    delete [] variable;
    delete [] itype;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Concatenation des variables d'objets Sequences.
 *
 *  arguments : reference sur un objet Format_error, nombre d'objets Sequences,
 *              pointeurs sur les objets Sequences, echantillon de reference pour les identificateurs.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::merge_variable(Format_error &error , int nb_sample ,
                                     const Sequences **iseq , int ref_sample) const

{
  bool status = true;
  register int i , j , k , m;
  int inb_variable , *itype , *iidentifier , *psequence , *csequence;
  Sequences *seq;
  const Sequences **pseq;


  seq = 0;
  error.init();

  if ((type[0] != INT_VALUE) && (type[0] != STATE) && (type[0] != TIME)) {
    status = false;
    ostringstream error_message , correction_message;
    error_message << STAT_label[STATL_SAMPLE] << " " << 1 << ": "
                  << SEQ_error[SEQR_VARIABLE_1_TYPE];
    correction_message << STAT_sequence_word[INT_VALUE] << " or " << STAT_sequence_word[STATE]
                       << " or " << STAT_sequence_word[TIME];
    error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
  }

  for (i = 0;i < nb_sample;i++) {
    if ((iseq[i]->type[0] != INT_VALUE) && (iseq[i]->type[0] != STATE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                    << SEQ_error[SEQR_VARIABLE_1_TYPE];
      correction_message << STAT_sequence_word[INT_VALUE] << " or " << STAT_sequence_word[STATE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }
  }

  for (i = 0;i < nb_sample;i++) {
    if (iseq[i]->nb_sequence != nb_sequence) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                    << SEQ_error[SEQR_NB_SEQUENCE];
      error.update((error_message.str()).c_str());
    }

    else {
      for (j = 0;j < nb_sequence;j++) {
        if (iseq[i]->length[j] != length[j]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                        << SEQ_label[SEQL_SEQUENCE] << " " << j + 1 << ": "
                        << SEQ_error[SEQR_SEQUENCE_LENGTH];
          error.update((error_message.str()).c_str());
        }
      }
    }
  }

  if ((ref_sample != I_DEFAULT) && ((ref_sample < 1) || (ref_sample > nb_sample + 1))) {
    status = false;
    error.update(STAT_error[STATR_SAMPLE_INDEX]);
  }

  if (status) {
    nb_sample++;
    pseq = new const Sequences*[nb_sample];

    pseq[0] = this;
    inb_variable = nb_variable;
    for (i = 1;i < nb_sample;i++) {
      pseq[i] = iseq[i - 1];
      inb_variable += iseq[i - 1]->nb_variable;
    }

    itype = new int[inb_variable];
    itype[0] = type[0];
    for (i = 1;i < inb_variable;i++) {
      itype[i] = INT_VALUE;
    }

    // calcul des identificateurs des sequences

    if (ref_sample == I_DEFAULT) {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 1;j < nb_sample;j++) {
          if (pseq[j]->identifier[i] != pseq[0]->identifier[i]) {
            break;
          }
        }
        if (j < nb_sample) {
          break;
        }
      }

      if (i < nb_sequence) {
        iidentifier = 0;
      }
      else {
        iidentifier = pseq[0]->identifier;
      }
    }

    else {
      iidentifier = pseq[--ref_sample]->identifier;
    }

    seq = new Sequences(inb_variable , itype , nb_sequence , iidentifier , length , false);

    // copie des sequences

    for (i = 0;i < nb_sequence;i++) {
      inb_variable = 0;
      for (j = 0;j < nb_sample;j++) {
        for (k = 0;k < pseq[j]->nb_variable;k++) {
          psequence = seq->sequence[i][inb_variable++];
          csequence = pseq[j]->sequence[i][k];
          for (m = 0;m < length[i];m++) {
            *psequence++ = *csequence++;
          }
        }
      }
    }

    inb_variable = 0;
    for (i = 0;i < nb_sample;i++) {
      for (j = 0;j < pseq[i]->nb_variable;j++) {
        seq->min_value[inb_variable] = pseq[i]->min_value[j];
        seq->max_value[inb_variable] = pseq[i]->max_value[j];
        if (pseq[i]->marginal[j]) {
          seq->marginal[inb_variable] = new Histogram(*(pseq[i]->marginal[j]));
        }
        inb_variable++;
      }
    }

    if (index_interval) {
      seq->index_interval = new Histogram(*index_interval);
    }

    delete [] pseq;
    delete [] itype;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Inversion du sens de parcours des sequences.
 *
 *  argument : reference sur un objet Format_error.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::reverse(Format_error &error) const

{
  Sequences *seq;


  error.init();

  if ((type[0] != INT_VALUE) && (type[0] != STATE) && (type[0] != POSITION)) {
    seq = 0;
    ostringstream correction_message;
    correction_message << STAT_sequence_word[INT_VALUE] << " or " << STAT_sequence_word[STATE]
                       << " or " << STAT_sequence_word[POSITION];
    error.correction_update(SEQ_error[SEQR_VARIABLE_1_TYPE] , (correction_message.str()).c_str());
  }

  else {
    seq = new Sequences(*this , 'c' , REVERSE);
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Selection des sequences sur un critere de longueur.
 *
 *  arguments : reference sur un objet Format_error, bornes sur la longueur,
 *              flag pour conserver ou rejeter les sequences selectionnees.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::length_select(Format_error &error , int min_length ,
                                    int imax_length , bool keep) const

{
  bool status = true;
  register int i;
  int inb_sequence , *index;
  Sequences *seq;


  seq = 0;
  error.init();

  if ((min_length < 1) || (min_length >= hlength->nb_value) || (min_length > imax_length)) {
    status = false;
    error.update(SEQ_error[SEQR_MIN_SEQUENCE_LENGTH]);
  }
  if ((imax_length < hlength->offset) || (imax_length < min_length)) {
    status = false;
    error.update(SEQ_error[SEQR_MAX_SEQUENCE_LENGTH]);
  }

  if (status) {

    // selection des sequences

    index = new int[nb_sequence];
    inb_sequence = 0;

    for (i = 0;i < nb_sequence;i++) {
      if ((length[i] >= min_length) && (length[i] <= imax_length)) {
        if (keep) {
          index[inb_sequence++] = i;
        }
      }

      else if (!keep) {
        index[inb_sequence++] = i;
      }
    }

    if (inb_sequence == 0) {
      status = false;
      error.update(STAT_error[STATR_EMPTY_SAMPLE]);
    }

    // copie des sequences selectionnees

    if (status) {
      seq = new Sequences(*this , inb_sequence , index);
    }

    delete [] index;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Suppression des premieres/dernieres series d'une valeur donne.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable,
 *              valeur, position ('b' : begin, 'e' : end),
 *              longueur maximum de la serie supprimee.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::remove_run(Format_error &error , int variable , int ivalue ,
                                 char position , int max_run_length) const

{
  bool status = true;
  register int i , j , k;
  int smax_length , inb_sequence , *iidentifier , *ilength , *index , *psequence , *csequence;
  Sequences *seq;


  seq = 0;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }
  if ((type[0] != INT_VALUE) && (type[0] != STATE) && (type[0] != TIME)) {
    status = false;
    ostringstream correction_message;
    correction_message << STAT_sequence_word[INT_VALUE] << " or " << STAT_sequence_word[STATE]
                       << " or " << STAT_sequence_word[TIME];
    error.correction_update(SEQ_error[SEQR_VARIABLE_1_TYPE] , (correction_message.str()).c_str());
  }

  if (status) {
    variable--;

    if ((type[variable] != INT_VALUE) && (type[variable] != STATE)) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_sequence_word[INT_VALUE] << " or " << STAT_sequence_word[STATE];
      error.correction_update(SEQ_error[SEQR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }

    else {
      if (!marginal[variable]) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                      << STAT_error[STATR_MARGINAL_HISTOGRAM];
        error.update((error_message.str()).c_str());
      }

      else if ((ivalue < marginal[variable]->offset) || (ivalue >= marginal[variable]->nb_value) ||
               (marginal[variable]->frequency[ivalue] == 0)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                      << STAT_label[STATL_VALUE] << " " << ivalue << " "
                      << SEQ_error[SEQR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }
    }
  }

  if ((max_run_length < I_DEFAULT) || (max_run_length == 0)) {
    status = false;
    error.update(SEQ_error[SEQR_MAX_RUN_LENGTH]);
  }

  if (status) {

    // calcul de la longueur des sequences

    iidentifier = new int[nb_sequence];
    ilength = new int[nb_sequence];
    index = new int[nb_sequence];
    inb_sequence = 0;

    for (i = 0;i < nb_sequence;i++) {
      if (max_run_length == I_DEFAULT) {
        smax_length = length[i];
      }
      else {
        smax_length = MIN(max_run_length , length[i]);
      }

      switch (position) {

      case 'b' : {
        csequence = sequence[i][variable];
        for (j = 0;j < smax_length;j++) {
          if (*csequence++ != ivalue) {
            break;
          }
        }
        break;
      }

      case 'e' : {
        csequence = sequence[i][variable] + length[i];
        for (j = 0;j < smax_length;j++) {
          if (*--csequence != ivalue) {
            break;
          }
        }
        break;
      }
      }

      if (j < length[i]) {
        iidentifier[inb_sequence] = identifier[i];
        ilength[inb_sequence] = length[i] - j;
        index[inb_sequence++] = i;
      }
    }

    // copie des sequences

    seq = new Sequences(nb_variable , type , inb_sequence , iidentifier , ilength , false);

    for (i = 0;i < seq->nb_sequence;i++) {
      for (j = 0;j < seq->nb_variable;j++) {
        psequence = seq->sequence[i][j];

        switch (position) {
        case 'b' :
          csequence = sequence[index[i]][j] + length[index[i]] - seq->length[i];
          break;
        case 'e' :
          csequence = sequence[index[i]][j];
          break;
        }

        for (k = 0;k < seq->length[i];k++) {
          *psequence++ = *csequence++;
        }
      }
    }

    for (i = 0;i < seq->nb_variable;i++) {
      seq->min_value_computation(i);
      seq->max_value_computation(i);
      seq->build_marginal_histogram(i);
    }

    if (seq->type[0] == TIME) {
      seq->index_interval_computation();
    }

    delete [] iidentifier;
    delete [] ilength;
    delete [] index;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Extraction de sous-sequences.
 *
 *  arguments : reference sur un objet Format_error, index minimum et
 *              maximum dans la sequence.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::index_extract(Format_error &error , int min_index , int max_index) const

{
  bool status = true;
  register int i , j , k;
  int inb_sequence , *iidentifier , *ilength , *index , *first_index , *pindex ,
      *psequence , *csequence;
  Sequences *seq;


  seq = 0;
  error.init();

  if ((min_index < 0) || (((type[0] == INT_VALUE) || (type[0] == STATE)) && (min_index >= max_length)) ||
      ((max_index != I_DEFAULT) && (min_index > max_index))) {
    status = false;
    error.update(SEQ_error[SEQR_MIN_INDEX]);
  }
  if ((max_index != I_DEFAULT) && ((max_index < 0) || (((type[0] == INT_VALUE) || (type[0] == STATE)) &&
        (max_index >= max_length)) || (max_index < min_index))) {
    status = false;
    error.update(SEQ_error[SEQR_MAX_INDEX]);
  }

  if (status) {

    // selection des sequences

    iidentifier = new int[nb_sequence];
    ilength = new int[nb_sequence];
    index = new int[nb_sequence];
    inb_sequence = 0;

    // index implicite

    if ((type[0] == INT_VALUE) || (type[0] == STATE)) {
      if (max_index == I_DEFAULT) {
        for (i = 0;i < nb_sequence;i++) {
          if (length[i] > min_index) {
            iidentifier[inb_sequence] = identifier[i];
            ilength[inb_sequence] = length[i] - min_index;
            index[inb_sequence++] = i;
          }
        }
      }

      else {
        for (i = 0;i < nb_sequence;i++) {
          if (length[i] > max_index) {
            iidentifier[inb_sequence] = identifier[i];
            ilength[inb_sequence] = max_index - min_index + 1;
            index[inb_sequence++] = i;
          }
        }
      }
    }

    // index explicite

    else {
      first_index = new int[nb_sequence];

      if (max_index == I_DEFAULT) {
        for (i = 0;i < nb_sequence;i++) {
          if (sequence[i][0][length[i] - 1] >= min_index) {
            iidentifier[inb_sequence] = identifier[i];

            pindex = sequence[i][0];
            for (j = 0;j < length[i];j++) {
              if (*pindex++ >= min_index) {
                break;
              }
            }
            first_index[inb_sequence] = j;
            ilength[inb_sequence] = length[i] - j;

            index[inb_sequence++] = i;
          }
        }
      }

      else {
        for (i = 0;i < nb_sequence;i++) {
          if (sequence[i][0][length[i] - 1] >= max_index) {
            iidentifier[inb_sequence] = identifier[i];

            pindex = sequence[i][0];
            for (j = 0;j < length[i];j++) {
              if (*pindex++ >= min_index) {
                break;
              }
            }
            first_index[inb_sequence] = j;
            ilength[inb_sequence] = length[i] - j;

            pindex = sequence[i][0] + length[i];
            for (j = 0;j < length[i];j++) {
              if (*--pindex <= max_index) {
                break;
              }
            }
            ilength[inb_sequence] -= j;

            index[inb_sequence++] = i;
          }
        }
      }
    }

    if (inb_sequence == 0) {
      status = false;
      error.update(STAT_error[STATR_EMPTY_SAMPLE]);
    }

    // extraction des sous-sequences

    if (status) {
      seq = new Sequences(nb_variable , type , inb_sequence , iidentifier , ilength , false);

      // index implicite

      if ((type[0] == INT_VALUE) || (type[0] == STATE)) {
        for (i = 0;i < seq->nb_sequence;i++) {
          for (j = 0;j < seq->nb_variable;j++) {
            psequence = seq->sequence[i][j];
            csequence = sequence[index[i]][j] + min_index;
            for (k = 0;k < seq->length[i];k++) {
              *psequence++ = *csequence++;
            }
          }
        }
      }

      // index explicite

      else {
        for (i = 0;i < seq->nb_sequence;i++) {
          for (j = 0;j < seq->nb_variable;j++) {
            psequence = seq->sequence[i][j];
            csequence = sequence[index[i]][j] + first_index[i];
            for (k = 0;k < seq->length[i];k++) {
              *psequence++ = *csequence++;
            }

            if (seq->type[j] == POSITION) {
              if (max_index == I_DEFAULT) {
                *psequence = *csequence;
              }
              else {
                *psequence = max_index;
              }
            }
          }
        }
      }

      for (i = 0;i < seq->nb_variable;i++) {
        if (seq->type[i] != POSITION) {
          seq->min_value_computation(i);
          seq->max_value_computation(i);
          seq->build_marginal_histogram(i);
        }
      }

      if ((seq->type[0] == TIME) || ((seq->type[0] == POSITION) && (seq->type[1] == INT_VALUE))) {
        seq->index_interval_computation();
      }

      delete [] iidentifier;
      delete [] ilength;
      delete [] index;
      if ((type[0] != INT_VALUE) && (type[0] != STATE)) {
        delete [] first_index;
      }
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Extraction par segmentation d'un objet Sequences.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable,
 *              nombre de valeurs, valeurs, flag zones correspondant aux valeurs
 *              extraites/pas extraites.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::segmentation_extract(Format_error &error , int variable ,
                                           int nb_value , int *ivalue , bool keep) const

{
  bool status = true;
  register int i , j , k , m;
  int nb_present_value , nb_selected_value , *pfrequency , *selected_value , *itype ,
      *psequence , *csequence , *zone_length , **zone_start;
  Sequences *seq;


  seq = 0;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }
  if ((type[0] != INT_VALUE) && (type[0] != STATE) && (type[0] != TIME)) {
    status = false;
    ostringstream correction_message;
    correction_message << STAT_sequence_word[INT_VALUE] << " or " << STAT_sequence_word[STATE]
                       << " or " << STAT_sequence_word[TIME];
    error.correction_update(SEQ_error[SEQR_VARIABLE_1_TYPE] , (correction_message.str()).c_str());
  }
  if (((type[0] == TIME) && (nb_variable == 2)) || (nb_variable == 1)) {
    status = false;
    error.update(STAT_error[STATR_NB_VARIABLE]);
  }

  if (status) {
    variable--;

    if ((type[variable] != INT_VALUE) && (type[variable] != STATE)) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_sequence_word[INT_VALUE] << " or " << STAT_sequence_word[STATE];
      error.correction_update(SEQ_error[SEQR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }

    else {
      if (!marginal[variable]) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                      << STAT_error[STATR_MARGINAL_HISTOGRAM];
        error.update((error_message.str()).c_str());
      }

      else {
        pfrequency = marginal[variable]->frequency + marginal[variable]->offset;
        nb_present_value = 0;
        for (i = marginal[variable]->offset;i < marginal[variable]->nb_value;i++) {
          if (*pfrequency > 0) {
            nb_present_value++;
          }
          pfrequency++;
        }

        if ((nb_value < 1) || (nb_value > (keep ? nb_present_value : nb_present_value - 1))) {
          status = false;
          error.update(SEQ_error[SEQR_NB_SELECTED_VALUE]);
        }

        else {
          selected_value = new int[marginal[variable]->nb_value];
          for (i = marginal[variable]->offset;i < marginal[variable]->nb_value;i++) {
            selected_value[i] = false;
          }

          for (i = 0;i < nb_value;i++) {
            if ((ivalue[i] < marginal[variable]->offset) || (ivalue[i] >= marginal[variable]->nb_value) ||
                (marginal[variable]->frequency[ivalue[i]] == 0)) {
              status = false;
              ostringstream error_message;
              error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                            << STAT_label[STATL_VALUE] << " " << ivalue[i] << " "
                            << SEQ_error[SEQR_NOT_PRESENT];
              error.update((error_message.str()).c_str());
            }

            else if (selected_value[ivalue[i]]) {
              status = false;
              ostringstream error_message;
              error_message << STAT_label[STATL_VALUE] << " " << ivalue[i] << " "
                            << STAT_error[STATR_ALREADY_SELECTED];
              error.update((error_message.str()).c_str());
            }
            else {
              selected_value[ivalue[i]] = true;
            }
          }

          delete [] selected_value;
        }
      }
    }
  }

  if (status) {
    itype = new int[nb_variable - 1];
    if (type[0] == TIME) {
      itype[0] = type[0];
    }
    for (i = (type[0] == TIME ? 1 : 0);i < nb_variable - 1;i++) {
      itype[i] = INT_VALUE;
    }

    switch (keep) {

    case false : {
      nb_selected_value = nb_present_value - nb_value;
      selected_value = new int[nb_selected_value];

      pfrequency = marginal[variable]->frequency + marginal[variable]->offset;
      i = 0;

      for (j = marginal[variable]->offset;j < marginal[variable]->nb_value;j++) {
        if (*pfrequency > 0) {
          for (k = 0;k < nb_value;k++) {
            if (ivalue[k] == j) {
              break;
            }
          }

          if (k == nb_value) {
            selected_value[i++] = j;
          }
        }

        pfrequency++;
      }
      break;
    }

    case true : {
      nb_selected_value = nb_value;
      selected_value = ivalue;
      break;
    }
    }

    // initialisations

    zone_length = new int[cumul_length];

    zone_start = new int*[cumul_length];
    for (i = 0;i < cumul_length;i++) {
      zone_start[i] = 0;
    }

    // recherche des sequences a extraire

    i = -1;

    for (j = 0;j < nb_sequence;j++) {
      psequence = sequence[j][variable];

      for (k = 0;k < nb_selected_value;k++) {
        if (*psequence == selected_value[k]) {
          break;
        }
      }

      if (k < nb_selected_value) {
        i++;
        zone_start[i] = new int[2];

        zone_start[i][0] = j;
        zone_start[i][1] = 0;
        zone_length[i] = 1;
      }

      for (k = 1;k < length[j];k++) {
        psequence++;

        for (m = 0;m < nb_selected_value;m++) {
          if (*psequence == selected_value[m]) {
            break;
          }
        }

        if (m  < nb_selected_value) {
          if (*psequence != *(psequence - 1)) {
            i++;
            zone_start[i] = new int[2];

            zone_start[i][0] = j;
            zone_start[i][1] = k;
            zone_length[i] = 1;
          }
          else {
            zone_length[i]++;
          }
        }
      }
    }

    // creation de l'objet Sequences

    seq = new Sequences(nb_variable - 1 , itype , i + 1 , 0 , zone_length , false);

    // copie des sequences

    for (i = 0;i < seq->nb_sequence;i++) {
      j = 0;
      for (k = 0;k < nb_variable;k++) {
        if (k != variable) {
          psequence = seq->sequence[i][j];
          csequence = sequence[zone_start[i][0]][k] + zone_start[i][1];
          for (m = 0;m < seq->length[i];m++) {
            *psequence++ = *csequence++;
          }
          j++;
        }
      }
    }

    for (i = 0;i < seq->nb_variable;i++) {
      seq->min_value_computation(i);
      seq->max_value_computation(i);
      seq->build_marginal_histogram(i);
    }

    if (seq->type[0] == TIME) {
      seq->index_interval_computation();
    }

    delete [] itype;

    if (!keep) {
      delete [] selected_value;
    }

    delete [] zone_length;
    for (i = 0;i < seq->nb_sequence;i++) {
      delete [] zone_start[i];
    }
    delete [] zone_start;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Cumul des valeurs successives des sequences.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::cumulate(Format_error &error , int variable) const

{
  bool status = true;
  register int i , j , k , m;
  int inb_variable , *itype , *psequence , *csequence;
  Sequences *seq;


  seq = 0;
  error.init();

  if (variable != I_DEFAULT) {
    if ((variable < 1) || (variable > nb_variable)) {
      status = false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
    }

    else {
      variable--;

      if ((type[variable] != INT_VALUE) && (type[variable] != STATE)) {
        status = false;
        ostringstream correction_message;
        correction_message << STAT_sequence_word[INT_VALUE] << " or " << STAT_sequence_word[STATE];
        error.correction_update(SEQ_error[SEQR_VARIABLE_TYPE] , (correction_message.str()).c_str());
      }
    }
  }

  if (status) {
    if (variable == I_DEFAULT) {
      inb_variable = nb_variable;
    }
    else {
      inb_variable = ((type[0] == TIME) || (type[0] == POSITION) ? 2 : 1);
    }

    itype = new int[inb_variable];

    i = 0;
    for (j = 0;j < nb_variable;j++) {
      if ((variable == I_DEFAULT) || (variable == j) || (type[j] == TIME) || (type[j] == POSITION)) {
        itype[i++] = type[j];
      }
    }

    seq = new Sequences(inb_variable , itype , nb_sequence , identifier , length , false);
    delete [] itype;

    for (i = 0;i < nb_sequence;i++) {
      j = 0;
      for (k = 0;k < nb_variable;k++) {
        if ((variable == I_DEFAULT) || (variable == k) || (type[k] == TIME) || (type[k] == POSITION)) {
          psequence = seq->sequence[i][j++];
          csequence = sequence[i][k];

          // copie des index

          if ((type[k] == TIME) || (type[k] == POSITION)) {
            for (m = 0;m < (type[k] == POSITION ? length[i] + 1 : length[i]);m++) {
              *psequence++ = *csequence++;
            }
          }

          // cumul des valeurs

          else {
            *psequence++ = *csequence++;
            for (m = 1;m < length[i];m++) {
              *psequence = *(psequence - 1) + *csequence++;
              psequence++;
            }
          }
        }
      }
    }

    for (i = 0;i < seq->nb_variable;i++) {
      if (seq->type[i] != POSITION) {
        seq->min_value_computation(i);
        seq->max_value_computation(i);
        seq->build_marginal_histogram(i);
      }
    }

    if ((seq->type[0] == TIME) || ((seq->type[0] == POSITION) && (seq->type[1] == INT_VALUE))) {
      seq->index_interval = new Histogram(*index_interval);
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture de sequences reelles.
 *
 *  arguments : stream, sequences reelles, format ('c' : colonne /
 *              'l' : ligne / 'v' : vecteur), indice de la sequence.
 *
 *--------------------------------------------------------------*/

ostream& Sequences::ascii_print(ostream &os , double ***real_sequence ,
                                char format , int index) const

{
  register int i , j , k;
  int *psequence;
  double *rsequence;


  os << nb_variable << " " << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] << endl;

  for (i = 0;i < nb_variable;i++) {
    os << "\n" << STAT_word[STATW_VARIABLE] << " " << i + 1 << " : " << STAT_sequence_word[type[i]];
  }
  os << endl;

  switch (format) {

  case 'c' : {
    for (i = 0;i < nb_sequence;i++) {
      if ((index == I_DEFAULT) || (index == i)) {
        os << "\n";
        for (j = 0;j < length[i];j++) {
          if ((type[0] == TIME) || (type[0] == POSITION)) {
            os << sequence[i][0][j] << " ";
          }
          for (k = 0;k < ((type[0] == TIME) || (type[0] == POSITION) ? nb_variable - 1 : nb_variable);k++) {
            os << real_sequence[i][k][j] << " ";
          }

          if (j < length[i] - 1) {
            if (nb_variable > 1) {
              os << "| ";
            }
          }
        }

        if (type[0] == POSITION) {
          os << "| " << sequence[i][0][length[i]];
        }

        os << "   # (" << identifier[i] << ")" << endl;
      }
    }
    break;
  }

  case 'l' : {
    for (i = 0;i < nb_sequence;i++) {
      if ((index == I_DEFAULT) || (index == i)) {
        os << "\n";
        if ((type[0] == TIME) || (type[0] == POSITION)) {
          psequence = sequence[i][0];
          for (j = 0;j < (type[0] == POSITION ? length[i] + 1 : length[i]);j++) {
            os << *psequence++ << " ";
          }
          os << endl;
        }

        for (j = 0;j < ((type[0] == TIME) || (type[0] == POSITION) ? nb_variable - 1 : nb_variable);j++) {
          rsequence = real_sequence[i][j];
          for (k = 0;k < length[i];k++) {
            os << *rsequence++ << " ";
          }
          if (j < ((type[0] == TIME) || (type[0] == POSITION) ? nb_variable - 2 : nb_variable - 1)) {
            os << endl;
          }
        }

        os << "   # (" << identifier[i] << ")" << endl;
      }
    }
    break;
  }

  case 'v' : {
    for (i = 0;i < nb_sequence;i++) {
      if ((index == I_DEFAULT) || (index == i)) {
        os << "\n";
        for (j = 0;j < length[i];j++) {
          if ((type[0] == TIME) || (type[0] == POSITION)) {
            os << sequence[i][0][j] << " ";
          }
          for (k = 0;k < ((type[0] == TIME) || (type[0] == POSITION) ? nb_variable - 1 : nb_variable);k++) {
            os << real_sequence[i][k][j] << " ";
          }
          os << identifier[i] << endl;
        }
      }
    }
    break;
  }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture de sequences reelles dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path, sequences reelles,
 *              format ('c' : colonne / 'l' : ligne / 'v' : vecteur), indice de la sequence.
 *
 *--------------------------------------------------------------*/

bool Sequences::ascii_print(Format_error &error , const char *path ,
                            double ***real_sequence , char format , int index) const

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
    ascii_print(out_file , real_sequence , format , index);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Mise a l'echelle de sequences reelles.
 *
 *  arguments : indice de la variable, sequences reelles.
 *
 *--------------------------------------------------------------*/

int Sequences::scaling_coefficient(int variable , double ***real_sequence)

{
  register int i , j;
  int scaling_coeff;
  double absolute_value_mean , *psequence;


  absolute_value_mean = 0.;
  for (i = 0;i < nb_sequence;i++) {
    psequence = real_sequence[i][variable];
    for (j = 0;j < length[i];j++) {
      absolute_value_mean += fabs(*psequence);
      psequence++;
    }
  }
  absolute_value_mean /= cumul_length;

  scaling_coeff = 1;
  while (absolute_value_mean < ABSOLUTE_VALUE_MEAN) {
    absolute_value_mean *= 10;
    scaling_coeff *= 10;
  }

  return scaling_coeff;
}


/*--------------------------------------------------------------*
 *
 *  Differenciation au 1er ordre des sequences.
 *
 *  arguments : reference sur un objet Format_error, stream, indice de la variable,
 *              premier element de la sequence garde ou pas, path.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::difference(Format_error &error , ostream &os , int variable ,
                                 bool first_element , const char *path) const

{
  bool status = true;
  register int i , j , k , m;
  int inb_variable , *itype , *ilength , *psequence , *csequence;
  Sequences *seq;


  seq = 0;
  error.init();

  if (type[0] == POSITION) {
    status = increasing_sequence_checking(error , 0 , true , SEQ_label[SEQL_SEQUENCE] ,
                                          STAT_sequence_word[POSITION]);
  }

  if (variable != I_DEFAULT) {
    if ((variable < 1) || (variable > nb_variable)) {
      status = false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
    }

    else {
      variable--;

      if ((type[variable] != INT_VALUE) && (type[variable] != STATE)) {
        status = false;
        ostringstream correction_message;
        correction_message << STAT_sequence_word[INT_VALUE] << " or " << STAT_sequence_word[STATE];
        error.correction_update(SEQ_error[SEQR_VARIABLE_TYPE] , (correction_message.str()).c_str());
      }
    }
  }

  if ((!first_element) && (hlength->offset < 2)) {
    status = false;
    ostringstream correction_message;
    correction_message << STAT_error[STATR_GREATER_THAN] << " " << 1;
    error.correction_update(SEQ_error[SEQR_MIN_SEQUENCE_LENGTH] , (correction_message.str()).c_str());
  }

  if (status) {
    if (variable == I_DEFAULT) {
      inb_variable = nb_variable;
    }
    else {
      inb_variable = ((type[0] == TIME) || (type[0] == POSITION) ? 2 : 1);
    }

    itype = new int[inb_variable];

    i = 0;
    for (j = 0;j < nb_variable;j++) {
      if ((variable == I_DEFAULT) || (variable == j) || (type[j] == TIME) || (type[j] == POSITION)) {
        itype[i++] = type[j];
      }
    }

    if (first_element) {
      ilength = length;
    }
    else {
      ilength = new int[nb_sequence];
      for (i = 0;i < nb_sequence;i++) {
        ilength[i] = length[i] - 1;
      }
    }

    seq = new Sequences(inb_variable , itype , nb_sequence , identifier , ilength , false);

    delete [] itype;
    if (!first_element) {
      delete [] ilength;
    }

    // copie des index

    if ((seq->type[0] == TIME) || (seq->type[0] == POSITION)) {
      for (i = 0;i < seq->nb_sequence;i++) {
        psequence = seq->sequence[i][0];
        if (first_element) {
          csequence = sequence[i][0];
        }
        else {
          csequence = sequence[i][0] + 1;
        }

        for (j = 0;j < (seq->type[0] == POSITION ? seq->length[i] + 1 : seq->length[i]);j++) {
          *psequence++ = *csequence++;
        }
      }
    }

    // differenciation des sequences

    if (((type[0] == TIME) || (type[0] == POSITION)) && ((index_interval->mean != 1.) ||
         (index_interval->variance > 0.))) {
      int scaling_coeff , *pindex;
      double *rsequence , ***real_sequence;


      real_sequence = new double**[seq->nb_sequence];
      for (i = 0;i < seq->nb_sequence;i++) {
        real_sequence[i] = new double*[seq->nb_variable];
        for (j = 0;j < seq->nb_variable;j++) {
          real_sequence[i][j] = new double[seq->length[i]];
        }
      }

      for (i = 0;i < nb_sequence;i++) {
        j = 0;
        for (k = 1;k < nb_variable;k++) {
          if ((type[k] == INT_VALUE) && ((variable == I_DEFAULT) || (variable == k))) {
            rsequence = real_sequence[i][j];
            pindex = sequence[i][0];
            csequence = sequence[i][k];

            if (first_element) {
              *rsequence++ = (double)*csequence / (double)*pindex;
            }
            for (m = 0;m < length[i] - 1;m++) {
              *rsequence++ = (double)(*(csequence + 1) - *csequence) / (double)(*(pindex + 1) - *pindex);
              pindex++;
              csequence++;
            }
            j++;
          }
        }
      }

      // ecriture des sequences reelles

      if (path) {
        status = seq->ascii_print(error , path , real_sequence , 'l');

        if (!status) {

#         ifdef MESSAGE
          os << error;
#         endif

        }
      }

      // mise a l'echelle des valeurs

      for (i = 1;i < seq->nb_variable;i++) {
        scaling_coeff = seq->scaling_coefficient(i , real_sequence);

        for (j = 0;j < seq->nb_sequence;j++) {
          psequence = seq->sequence[j][i];
          rsequence = real_sequence[j][i - 1];
          for (k = 0;k < seq->length[j];k++) {
            *psequence++ = (int)round(*rsequence++ * scaling_coeff);
          }
        }
      }

      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < seq->nb_variable;j++) {
          delete [] real_sequence[i][j];
        }
        delete [] real_sequence[i];
      }
      delete [] real_sequence;
    }

    else {
      for (i = 0;i < nb_sequence;i++) {
        j = 0;
        for (k = 0;k < nb_variable;k++) {
          if ((variable == I_DEFAULT) || (variable == k) || (type[k] == TIME) || (type[k] == POSITION)) {
            if ((type[k] == INT_VALUE) || (type[k] == STATE)) {
              psequence = seq->sequence[i][j];
              csequence = sequence[i][k];

              if (first_element) {
                *psequence++ = *csequence;
              }
              for (m = 0;m < length[i] - 1;m++) {
                *psequence++ = *(csequence + 1) - *csequence;
                csequence++;
              }
            }

            j++;
          }
        }
      }
    }

    for (i = 0;i < seq->nb_variable;i++) {
      if (seq->type[i] != POSITION) {
        seq->min_value_computation(i);
        seq->max_value_computation(i);
        seq->build_marginal_histogram(i);
      }
    }

    if ((seq->type[0] == TIME) || ((seq->type[0] == POSITION) && (seq->type[1] == INT_VALUE))) {
      seq->index_interval_computation();
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Filtrage de type moyenne mobile des sequences.
 *
 *  arguments : reference sur un objet Format_error, stream, demi-largeur du filtre,
 *              filtre, indice de la variable, debut/fin garde ou pas,
 *              tendance ou residus (par soustraction ou par division), path,
 *              format ('s' : sequence / 'v' : vecteur).
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::moving_average(Format_error &error , ostream &os , int nb_point ,
                                     double *filter , int variable , bool begin_end ,
                                     int output , const char *path , char format) const

{
  bool status = true;
  register int i , j , k , m , n;
  int inb_variable , *itype , *ilength , *psequence , *csequence;
  double *rsequence , *ppoint , ***real_sequence;
  Sequences *seq;


  seq = 0;
  error.init();

  if (variable == I_DEFAULT) {
    if ((index_interval) && (index_interval->variance > 0.)) {
      status = false;
      error.update(SEQ_error[SEQR_UNEQUAL_INDEX_INTERVALS]);
    }
  }

  else {
    if ((variable < 1) || (variable > nb_variable)) {
      status = false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
    }

    else {
      variable--;

      if ((type[variable] != INT_VALUE) && (type[variable] != STATE)) {
        status = false;
        ostringstream correction_message;
        correction_message << STAT_sequence_word[INT_VALUE] << " or " << STAT_sequence_word[STATE];
        error.correction_update(SEQ_error[SEQR_VARIABLE_TYPE] , (correction_message.str()).c_str());
      }
    }
  }

  if ((!begin_end) && (hlength->offset < 2 * nb_point + 1)) {
    status = false;
    ostringstream correction_message;
    correction_message << STAT_error[STATR_GREATER_THAN] << " " << 2 * nb_point;
    error.correction_update(SEQ_error[SEQR_MIN_SEQUENCE_LENGTH] , (correction_message.str()).c_str());
  }

  if (status) {
    if (variable == I_DEFAULT) {
      inb_variable = nb_variable;
    }
    else {
      inb_variable = ((type[0] == TIME) || (type[0] == POSITION) ? 2 : 1);
    }

    itype = new int[inb_variable];

    i = 0;
    for (j = 0;j < nb_variable;j++) {
      if ((variable == I_DEFAULT) || (variable == j) || (type[j] == TIME) || (type[j] == POSITION)) {
        itype[i++] = type[j];
      }
    }

    if (begin_end) {
      ilength = length;
    }
    else {
      ilength = new int[nb_sequence];
      for (i = 0;i < nb_sequence;i++) {
        ilength[i] = length[i] - 2 * nb_point;
      }
    }

    seq = new Sequences(inb_variable , itype , nb_sequence , identifier , ilength , false);

    delete [] itype;
    if (!begin_end) {
      delete [] ilength;
    }

    // copie des index

    if ((seq->type[0] == TIME) || (seq->type[0] == POSITION)) {
      for (i = 0;i < seq->nb_sequence;i++) {
        psequence = seq->sequence[i][0];
        if (begin_end) {
          csequence = sequence[i][0];
        }
        else {
          csequence = sequence[i][0] + nb_point;
        }

        for (j = 0;j < (seq->type[0] == POSITION ? seq->length[i] + 1 : seq->length[i]);j++) {
          *psequence++ = *csequence++;
        }
      }
    }

    // filtrage de type moyenne mobile

    real_sequence = new double**[seq->nb_sequence];
    for (i = 0;i < seq->nb_sequence;i++) {
      real_sequence[i] = new double*[seq->nb_variable];
      for (j = 0;j < seq->nb_variable;j++) {
        real_sequence[i][j] = new double[seq->length[i]];
      }
    }

    for (i = 0;i < nb_sequence;i++) {
      j = 0;
      for (k = 0;k < nb_variable;k++) {
        if (((type[k] == INT_VALUE) || (type[k] == STATE)) &&
            ((variable == I_DEFAULT) || (variable == k))) {
          rsequence = real_sequence[i][j];

          if (begin_end) {
            for (m = 0;m < MIN(nb_point , length[i]);m++) {
              csequence = sequence[i][k];
              ppoint = filter;
              *rsequence = 0.;
              for (n = 0;n < 2 * nb_point + 1;n++) {
                *rsequence += *csequence * *ppoint++;
                if ((m - nb_point + n >= 0) && (m - nb_point + n < length[i] - 1)) {
                  csequence++;
                }
              }
              rsequence++;
            }
          }

          for (m = nb_point;m < length[i] - nb_point;m++) {
            csequence = sequence[i][k] + m - nb_point;
            ppoint = filter;
            *rsequence = 0.;
            for (n = 0;n < 2 * nb_point + 1;n++) {
              *rsequence += *csequence++ * *ppoint++;
            }
            rsequence++;
          }

          if (begin_end) {
            for (m = MAX(length[i] - nb_point , nb_point);m < length[i];m++) {
              csequence = sequence[i][k] + m - nb_point;
              ppoint = filter;
              *rsequence = 0.;
              for (n = 0;n < 2 * nb_point + 1;n++) {
                *rsequence += *csequence * *ppoint++;
                if ((m - nb_point + n >= 0) && (m - nb_point + n < length[i] - 1)) {
                  csequence++;
                }
              }
              rsequence++;
            }
          }

          if (output != TREND) {
            rsequence = real_sequence[i][j];
            if (begin_end) {
              csequence = sequence[i][k];
            }
            else {
              csequence = sequence[i][k] + nb_point;
            }

            switch (output) {

            case SUBTRACTION_RESIDUAL : {
              for (m = 0;m < seq->length[i];m++) {
                *rsequence = *csequence++ - *rsequence;
                rsequence++;
              }
              break;
            }

            case DIVISION_RESIDUAL : {
              for (m = 0;m < seq->length[i];m++) {
                *rsequence = *csequence++ / *rsequence;
                rsequence++;
              }
              break;
            }
            }
          }

          j++;
        }
      }
    }

    // ecriture des sequences reelles

    if (path) {
      status = seq->ascii_print(error , path , real_sequence , (format == 's' ? 'l' : format));

      if (!status) {

#       ifdef MESSAGE
        os << error;
#       endif

      }
    }

    // arrondi des valeurs

    i = 0;
    for (j = 0;j < seq->nb_variable;j++) {
      if ((seq->type[j] == INT_VALUE) || (seq->type[j] == STATE)) {
        for (k = 0;k < seq->nb_sequence;k++) {
          psequence = seq->sequence[k][j];
          rsequence = real_sequence[k][i];

          if (output != DIVISION_RESIDUAL) {
            for (m = 0;m < seq->length[k];m++) {
              *psequence++ = (int)round(*rsequence++);
            }
          }

          else {
            for (m = 0;m < seq->length[k];m++) {
              *psequence++ = (int)round(RESIDUAL_SCALE * *rsequence++);
            }
          }
        }

        i++;
      }
    }

    for (i = 0;i < seq->nb_sequence;i++) {
      for (j = 0;j < seq->nb_variable;j++) {
        delete [] real_sequence[i][j];
      }
      delete [] real_sequence[i];
    }
    delete [] real_sequence;

    for (i = 0;i < seq->nb_variable;i++) {
      if (seq->type[i] != POSITION) {
        seq->min_value_computation(i);
        seq->max_value_computation(i);
        seq->build_marginal_histogram(i);
      }
    }

    if ((seq->type[0] == TIME) || ((seq->type[0] == POSITION) && (seq->type[1] == INT_VALUE))) {
      seq->index_interval_computation();
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Filtrage de type moyenne mobile des sequences.
 *
 *  arguments : reference sur un objet Format_error, stream, loi symmetrique,
 *              indice de la variable, debut/fin supprime ou pas,
 *              tendance ou residus (par soustraction ou par division), path,
 *              format ('s' : sequence / 'v' : vecteur).
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::moving_average(Format_error &error , ostream &os , const Distribution &dist ,
                                     int variable , bool begin_end , int output , const char *path,
                                     char format) const

{
  bool status = true;
  Sequences *seq;


  seq = 0;
  error.init();

  if ((dist.offset != 0) || ((dist.nb_value - dist.offset) % 2 == 0)) {
    status = false;
    error.correction_update(STAT_error[STATR_NB_VALUE] , STAT_error[STATR_ODD]);
  }
  if (fabs(dist.skewness_computation()) > SKEWNESS_ROUNDNESS) {
    status = false;
    error.update(STAT_error[STATR_NON_SYMMETRICAL_DISTRIBUTION]);
  }
  if (dist.complement > 0.) {
    status = false;
    error.update(STAT_error[STATR_UNPROPER_DISTRIBUTION]);
  }

  if (status) {
    seq = moving_average(error , os , dist.nb_value / 2 , dist.mass ,
                         variable , begin_end , output , path , format);
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des sequences des temps de retour pour une valeur prise
 *  par une variable donnee.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable, valeur.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::recurrence_time_sequences(Format_error &error , int variable , int value) const

{
  bool status = true;
  register int i , j;
  int inb_sequence , ilength , previous_index , *psequence , *rsequence;
  Sequences *seq;


  seq = 0;
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
      correction_message << STAT_sequence_word[INT_VALUE] << " or " << STAT_sequence_word[STATE];
      error.correction_update(SEQ_error[SEQR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }

    else {
      if (!marginal[variable]) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                      << STAT_error[STATR_MARGINAL_HISTOGRAM];
        error.update((error_message.str()).c_str());
      }

      else if ((value < marginal[variable]->offset) || (value >= marginal[variable]->nb_value) ||
               (marginal[variable]->frequency[value] == 0)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                      << STAT_label[STATL_VALUE] << " " << value << " "
                      << SEQ_error[SEQR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }
    }
  }

  if (status) {
    seq = new Sequences(1 , nb_sequence);

    // calcul des sequences des temps de retour

    inb_sequence = 0;

    for (i = 0;i < nb_sequence;i++) {
      psequence = sequence[i][variable];
      previous_index = 0;
      ilength = 0;
      for (j = 0;j < length[i];j++) {
        if (*psequence == value) {
          if (ilength == 0) {
            seq->sequence[inb_sequence][0] = new int[length[i]];
            rsequence = seq->sequence[inb_sequence][0];
          }

          *rsequence++ = j - previous_index;
          previous_index = j;
          ilength++;
        }

        psequence++;
      }

      if (ilength > 0) {
        seq->length[inb_sequence] = ilength;
        seq->identifier[inb_sequence++] = identifier[i];
      }
    }

    seq->nb_sequence = inb_sequence;

    seq->max_length_computation();
    seq->cumul_length_computation();
    seq->build_length_histogram();

    seq->min_value_computation(0);
    seq->max_value_computation(0);
    seq->build_marginal_histogram(0);
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Discretisation des positions.
 *
 *  arguments : reference sur un objet Format_error,
 *              pas de discretisation.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::transform_position(Format_error &error , int step) const

{
  bool status = true;
  register int i , j , k , m;
  int inter_pos , nb_unit , *ilength , *position , **psequence;
  Sequences *seq;


  seq = 0;
  error.init();

  if (type[0] != POSITION) {
    status = false;
    error.correction_update(SEQ_error[SEQR_VARIABLE_1_TYPE] , STAT_sequence_word[POSITION]);
  }

  if ((step < 1) || (step > index_interval->mean))  {
    status = false;
    error.update(SEQ_error[SEQR_POSITION_STEP]);
  }

  if (status) {
    ilength = new int[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      ilength[i] = sequence[i][0][length[i]] / step + 1 + length[i];
    }

    seq = new Sequences(nb_variable - 1 , nb_sequence , identifier , ilength , false);
    delete [] ilength;

    // extraction des sequences

    psequence = new int*[nb_variable];

    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < nb_variable - 1;j++) {
        psequence[j] = seq->sequence[i][j];
      }
      position = sequence[i][0];
      seq->length[i] = 0;

      for (j = 0;j <= length[i];j++) {
        if (j == 0) {
          inter_pos = *position;
        }
        else {
          inter_pos = *position - *(position - 1);
        }
        position++;

        nb_unit = (inter_pos % step == 0 ? inter_pos / step : inter_pos / step + 1);
        if (j < length[i]) {
          nb_unit -= MIN(nb_unit , 1);
        }

        for (k = 0;k < nb_variable - 1;k++) {
          for (m = 0;m < nb_unit;m++) {
            *psequence[k]++ = min_value[k + 1] - 1;
          }
          if (j < length[i]) {
            *psequence[k]++ = sequence[i][k + 1][j];
          }
        }

        if (j < length[i]) {
          seq->length[i] += nb_unit + 1;
        }
        else {
          seq->length[i] += nb_unit;
        }
      }
    }

    seq->max_length_computation();
    seq->cumul_length_computation();
    delete seq->hlength;
    seq->build_length_histogram();

    for (i = 0;i < seq->nb_variable;i++) {
      seq->min_value[i] = min_value[i + 1] - 1;
      seq->max_value[i] = max_value[i + 1];
      seq->build_marginal_histogram(i);
    }

    delete [] psequence;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Changement d'unite d'une variable.
 *
 *  arguments : reference sur un objet Format_error, variable, facteur d'echelle.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::scaling(Format_error &error , int variable , int scaling_coeff) const

{
  bool status = true;
  register int i , j , k;
  int *psequence , *csequence;
  Sequences *seq;


  seq = 0;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if (type[variable] != INT_VALUE) {
      status = false;
      error.correction_update(SEQ_error[SEQR_VARIABLE_TYPE] , STAT_sequence_word[INT_VALUE]);
    }
  }

  if (scaling_coeff < 1) {
    status = false;
    error.update(SEQ_error[SEQR_SCALING_COEFF]);
  }

  if (status) {
    seq = new Sequences(nb_variable , type , nb_sequence , identifier , length , false);

    for (i = 0;i < seq->nb_sequence;i++) {
      for (j = 0;j < seq->nb_variable;j++) {
        psequence = seq->sequence[i][j];
        csequence = sequence[i][j];

        // mise a l'echelle des valeurs

        if (j == variable) {
          for (k = 0;k < (seq->type[j] == POSITION ? seq->length[i] + 1 : seq->length[i]);k++) {
            *psequence++ = *csequence++ * scaling_coeff;
          }
        }

        // copie des valeurs

        else {
          for (k = 0;k < (seq->type[j] == POSITION ? seq->length[i] + 1 : seq->length[i]);k++) {
            *psequence++ = *csequence++;
          }
        }
      }
    }

    for (i = 0;i < seq->nb_variable;i++) {
      if (i == variable) {
        seq->min_value[i] = min_value[i] * scaling_coeff;
        seq->max_value[i] = max_value[i] * scaling_coeff;
        seq->build_marginal_histogram(i);
      }

      else if (seq->type[i] != POSITION) {
        seq->min_value[i] = min_value[i];
        seq->max_value[i] = max_value[i];
        if (marginal[i]) {
          seq->marginal[i] = new Histogram(*(marginal[i]));
        }
      }
    }

    if ((seq->type[0] == TIME) || ((seq->type[0] == POSITION) && (seq->type[1] == INT_VALUE))) {
      if (variable == 0) {
        seq->index_interval_computation();
      }
      else {
        seq->index_interval = new Histogram(*(index_interval));
      }
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Croisement des sequences.
 *
 *  argument : reference sur un objet Format_error.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::cross(Format_error &error) const

{
  bool status = true;
  register int i , j , k , m;
  int sense = 0 , *ilength;
  Sequences *seq;


  seq = 0;
  error.init();

  if (type[0] != INT_VALUE) {
    status = false;
    error.correction_update(SEQ_error[SEQR_VARIABLE_1_TYPE] , STAT_sequence_word[INT_VALUE]);
  }

  for (i = 1;i < nb_sequence;i++) {
    if (length[i] > length[i - 1]) {
      if (sense == 0) {
        sense++;
      }
      else if (sense == -1) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_SEQUENCE] << " " << i + 1 << ": "
                      << SEQ_error[SEQR_LENGTH];
        error.update((error_message.str()).c_str());
      }
    }

    else if (length[i] < length[i - 1]) {
      if (sense == 0) {
        sense--;
      }
      else if (sense == 1) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_SEQUENCE] << " " << i + 1 << ": "
                      << SEQ_error[SEQR_LENGTH];
        error.update((error_message.str()).c_str());
      }
    }
  }

  if (status) {
    ilength = new int[max_length];
    for (i = 0;i < max_length;i++) {
      ilength[i] = nb_sequence;
    }

    seq = new Sequences(nb_variable , max_length , 0 , ilength , false);
    delete [] ilength;

    // constitution des sequences croisees

    for (i = 0;i < seq->nb_sequence;i++) {
      j = 0;
      while (length[j] <= i) {
        j++;
      }

      k = 0;
      do {
        for (m = 0;m < seq->nb_variable;m++) {
          seq->sequence[i][m][k] = sequence[j][m][i];
        }
        j++;
        k++;
      }
      while ((j < nb_sequence) && (length[j] > i));

      seq->length[i] = k;
    }

    seq->max_length = nb_sequence;
    seq->cumul_length = cumul_length;
    delete seq->hlength;
    seq->build_length_histogram();

    for (i = 0;i < seq->nb_variable;i++) {
      seq->min_value[i] = min_value[i];
      seq->max_value[i] = max_value[i];
      if (marginal[i]) {
        seq->marginal[i] = new Histogram(*(marginal[i]));
      }
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la valeur minimum prise par une variable.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void Sequences::min_value_computation(int variable)

{
  if (type[variable] != POSITION) {
    register int i , j;
    int *psequence;


    min_value[variable] = INT_MAX;

    for (i = 0;i < nb_sequence;i++) {
      psequence = sequence[i][variable];
      for (j = 0;j < length[i];j++) {
        if (*psequence < min_value[variable]) {
          min_value[variable] = *psequence;
        }
        psequence++;
      }
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

void Sequences::max_value_computation(int variable)

{
  if (type[variable] != POSITION) {
    register int i , j;
    int *psequence;


    max_value[variable] = INT_MIN;

    for (i = 0;i < nb_sequence;i++) {
      psequence = sequence[i][variable];
      for (j = 0;j < length[i];j++) {
        if (*psequence > max_value[variable]) {
          max_value[variable] = *psequence;
        }
        psequence++;
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de l'histogramme du nombre d'occurrences des valeurs
 *  (pour une variable donnee).
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void Sequences::marginal_histogram_computation(int variable)

{
  register int i , j;
  int *pfrequency , *psequence;


  pfrequency = marginal[variable]->frequency;
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    *pfrequency++ = 0;
  }

  for (i = 0;i < nb_sequence;i++) {
    psequence = sequence[i][variable];
    for (j = 0;j < length[i];j++) {
      (marginal[variable]->frequency[*psequence++])++;
    }
  }

  marginal[variable]->offset = min_value[variable];
  marginal[variable]->nb_element_computation();
  marginal[variable]->max_computation();
  marginal[variable]->mean_computation();
  marginal[variable]->variance_computation();
}


/*--------------------------------------------------------------*
 *
 *  Construction de l'histogramme du nombre d'occurrences des valeurs
 *  (pour une variable donnee).
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void Sequences::build_marginal_histogram(int variable)

{
  if ((type[variable] != POSITION) &&
      (min_value[variable] >= 0) && (max_value[variable] <= MARGINAL_MAX_VALUE)) {
    marginal[variable] = new Histogram(max_value[variable] + 1);
    marginal_histogram_computation(variable);
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la moyenne d'une variable.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

double Sequences::mean_computation(int variable) const

{
  register int i , j;
  int *psequence;
  double mean = D_INF;


  if (marginal[variable]) {
    mean = marginal[variable]->mean;
  }

  else if ((type[variable] == INT_VALUE) || (type[variable] == NB_INTERNODE)) {
    mean = 0.;

    for (i = 0;i < nb_sequence;i++) {
      psequence = sequence[i][variable];
      for (j = 0;j < length[i];j++) {
        mean += *psequence++;
      }
    }

    mean /= cumul_length;
  }

  return mean;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la variance d'une variable.
 *
 *  arguments : indice de la variable, moyenne.
 *
 *--------------------------------------------------------------*/

double Sequences::variance_computation(int variable , double mean) const

{
  register int i , j;
  int *psequence;
  double variance = D_DEFAULT , diff;


  if (marginal[variable]) {
    variance = marginal[variable]->variance;
  }

  else if ((type[variable] == INT_VALUE) || (type[variable] == NB_INTERNODE)) {
    variance = 0.;

    if (cumul_length > 1) {
      for (i = 0;i < nb_sequence;i++) {
        psequence = sequence[i][variable];
        for (j = 0;j < length[i];j++) {
          diff = *psequence++ - mean;
          variance += diff * diff;
        }
      }

      variance /= (cumul_length - 1);
    }
  }

  return variance;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de l'ecart absolu moyen d'une variable.
 *
 *  arguments : indice de la variable, moyenne.
 *
 *--------------------------------------------------------------*/

double Sequences::mean_absolute_deviation_computation(int variable , double mean) const

{
  register int i , j;
  int *psequence;
  double mean_absolute_deviation = D_DEFAULT;


  if (marginal[variable]) {
    mean_absolute_deviation = marginal[variable]->mean_absolute_deviation_computation();
  }

  else if ((type[variable] == INT_VALUE) || (type[variable] == NB_INTERNODE)) {
    mean_absolute_deviation = 0.;

    for (i = 0;i < nb_sequence;i++) {
      psequence = sequence[i][variable];
      for (j = 0;j < length[i];j++) {
        mean_absolute_deviation += fabs(*psequence - mean);
        psequence++;
      }
    }

    mean_absolute_deviation /= cumul_length;
  }

  return mean_absolute_deviation;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la difference absolue moyenne d'une variable.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

double Sequences::mean_absolute_difference_computation(int variable) const

{
  register int i , j , k , m;
  int *psequence1 , *psequence2;
  double mean_absolute_difference = D_DEFAULT;


  if ((type[variable] == INT_VALUE) || (type[variable] == NB_INTERNODE)) {
    mean_absolute_difference = 0.;

    for (i = 0;i < nb_sequence;i++) {
      psequence1 = sequence[i][variable];
      for (j = 0;j < length[i];j++) {
        psequence2 = sequence[i][variable] + j + 1;
        for (k = j + 1;k < length[i];k++) {
          mean_absolute_difference += abs(*psequence1 - *psequence2);
          psequence2++;
        }
        psequence1++;
      }

      for (j = i + 1;j < nb_sequence;j++) {
        psequence1 = sequence[i][variable];
        for (k = 0;k < length[i];k++) {
          psequence2 = sequence[j][variable];
          for (m = 0;m < length[j];m++) {
            mean_absolute_difference += abs(*psequence1 - *psequence2);
            psequence2++;
          }
          psequence1++;
        }
      }
    }

    mean_absolute_difference = 2 * mean_absolute_difference / (cumul_length * (double)(cumul_length - 1));
  }

  return mean_absolute_difference;
}


/*--------------------------------------------------------------*
 *
 *  Calcul du coefficient d'asymetrie d'une variable.
 *
 *  arguments : indice de la variable, moyenne, variance.
 *
 *--------------------------------------------------------------*/

double Sequences::skewness_computation(int variable , double mean , double variance) const

{
  register int i , j;
  int *psequence;
  double skewness = D_INF , diff;


  if (marginal[variable]) {
    skewness = marginal[variable]->skewness_computation();
  }

  else if ((type[variable] == INT_VALUE) || (type[variable] == NB_INTERNODE)) {
    skewness = 0.;

    if ((cumul_length > 2) && (variance > 0.)) {
      for (i = 0;i < nb_sequence;i++) {
        psequence = sequence[i][variable];
        for (j = 0;j < length[i];j++) {
          diff = *psequence++ - mean;
          skewness += diff * diff * diff;
        }
      }

      skewness = skewness * cumul_length / ((cumul_length - 1) * (double)(cumul_length - 2) *
                  pow(variance , 1.5));
    }
  }

  return skewness;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de l'exces d'applatissement d'une variable :
 *  exces d'applatissement = coefficient d'applatissement - 3.
 *
 *  arguments : indice de la variable, moyenne, variance.
 *
 *--------------------------------------------------------------*/

double Sequences::kurtosis_computation(int variable , double mean , double variance) const

{
  register int i , j;
  int *psequence;
  double kurtosis = D_INF , diff;


  if (marginal[variable]) {
    kurtosis = marginal[variable]->kurtosis_computation();
  }

  else if ((type[variable] == INT_VALUE) || (type[variable] == NB_INTERNODE)) {
    if (variance == 0.) {
      kurtosis = -2.;
    }

    else {
      kurtosis = 0.;

      for (i = 0;i < nb_sequence;i++) {
        psequence = sequence[i][variable];
        for (j = 0;j < length[i];j++) {
          diff = *psequence++ - mean;
          kurtosis += diff * diff * diff * diff;
        }
      }

      kurtosis = kurtosis / ((cumul_length - 1) * variance * variance) - 3.;
    }
  }

  return kurtosis;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la longueur maximum des sequences.
 *
 *--------------------------------------------------------------*/

void Sequences::max_length_computation()

{
  register int i;
  int *plength;


  max_length = 0;
  plength = length;
  for (i = 0;i < nb_sequence;i++) {
    if (*plength > max_length) {
      max_length = *plength;
    }
    plength++;
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la longueur cumulee des sequences.
 *
 *--------------------------------------------------------------*/

void Sequences::cumul_length_computation()

{
  register int i;
  int *plength;


  cumul_length = 0;
  plength = length;
  for (i = 0;i < nb_sequence;i++) {
    cumul_length += *plength++;
  }
}


/*--------------------------------------------------------------*
 *
 *  Construction de l'histogramme des longueurs des sequences.
 *
 *--------------------------------------------------------------*/

void Sequences::build_length_histogram()

{
  register int i;
  int *plength;


  hlength = new Histogram(max_length + 1);

  plength = length;
  hlength->nb_element = nb_sequence;
  for (i = 0;i < nb_sequence;i++) {
    (hlength->frequency[*plength++])++;
  }

  hlength->nb_value_computation();
  hlength->offset_computation();
  hlength->max_computation();
  hlength->mean_computation();
  hlength->variance_computation();
}


/*--------------------------------------------------------------*
 *
 *  Calcul des dates a partir des intervalles de temps /
 *  calcul des positions a partir des intervalles inter-positions.
 *
 *--------------------------------------------------------------*/

void Sequences::index_computation()

{
  if ((type[0] == TIME_INTERVAL) || (type[0] == POSITION_INTERVAL)) {
    register int i , j;
    int *psequence;


    switch (type[0]) {
    case TIME_INTERVAL :
      type[0] = TIME;
      break;
    case POSITION_INTERVAL :
      type[0] = POSITION;
      break;
    }

    for (i = 0;i < nb_sequence;i++) {
      psequence = sequence[i][0] + 1;
      for (j = 1;j < (type[0] == POSITION ? length[i] + 1 : length[i]);j++) {
        *psequence += *(psequence - 1);
        psequence++;
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Extraction de l'histogramme des intervalles entre index successifs.
 *
 *--------------------------------------------------------------*/

void Sequences::index_interval_computation()

{
  if ((type[0] == TIME) || ((type[0] == POSITION) && (type[1] == INT_VALUE))) {
    register int i , j;
    int max_index , *psequence;


    max_index = 0;
    for (i = 0;i < nb_sequence;i++) {
      if (sequence[i][0][length[i] - 1] > max_index) {
        max_index = sequence[i][0][length[i] - 1];
      }
    }

    index_interval = new Histogram(max_index + 1);

    // constitution de l'histogramme des intervalles entre index successifs

    for (i = 0;i < nb_sequence;i++) {
      psequence = sequence[i][0] + 1;
      for (j = 1;j < length[i];j++) {
        (index_interval->frequency[*psequence - *(psequence - 1)])++;
        psequence++;
      }
    }

    index_interval->nb_value_computation();
    index_interval->offset_computation();
    index_interval->nb_element_computation();
    index_interval->max_computation();
    index_interval->mean_computation();
    index_interval->variance_computation();
  }
}


/*--------------------------------------------------------------*
 *
 *  Extraction de l'histogramme des intervalles entre index successifs
 *  pour une valeur donnee d'une variable donnee.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable, valeur.
 *
 *--------------------------------------------------------------*/

Histogram* Sequences::value_index_interval_computation(Format_error &error , int variable ,
                                                       int value) const

{
  bool status = true;
  register int i , j;
  int max_index , previous_index , *pindex , *psequence;
  Histogram *value_index_interval;


  value_index_interval = 0;
  error.init();

  if ((type[0] != TIME) && (type[0] != POSITION)) {
    status = false;
    ostringstream correction_message;
    correction_message << STAT_sequence_word[TIME] << " or "
                       << STAT_sequence_word[POSITION];
    error.correction_update(SEQ_error[SEQR_VARIABLE_1_TYPE] , (correction_message.str()).c_str());
  }

  if ((variable < 2) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if (!marginal[variable]) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                    << STAT_error[STATR_MARGINAL_HISTOGRAM];
      error.update((error_message.str()).c_str());
    }

    else if ((value < marginal[variable]->offset) || (value >= marginal[variable]->nb_value) ||
             (marginal[variable]->frequency[value] <= 1)) {
      status = false;
      error.update(SEQ_error[SEQR_VALUE]);
    }
  }

  if (status) {
    max_index = 0;
    for (i = 0;i < nb_sequence;i++) {
     if (sequence[i][0][length[i] - 1] > max_index) {
        max_index = sequence[i][0][length[i] - 1];
      }
    }

    value_index_interval = new Histogram(max_index + 1);

    for (i = 0;i < nb_sequence;i++) {
      pindex = sequence[i][0];
      psequence = sequence[i][variable];
      previous_index = I_DEFAULT;

      for (j = 0;j < length[i];j++) {
        if (*psequence == value) {
          if (previous_index != I_DEFAULT) {
            (value_index_interval->frequency[*pindex - previous_index])++;
          }
          previous_index = *pindex;
        }

        pindex++;
        psequence++;
      }
    }

    // extraction des caracteristiques de l'histogramme

    value_index_interval->nb_value_computation();
    value_index_interval->offset_computation();
    value_index_interval->nb_element_computation();

    if (value_index_interval->nb_element == 0) {
      delete value_index_interval;
      value_index_interval = 0;
      error.update(STAT_error[STATR_EMPTY_HISTOGRAM]);
    }

    else {
      value_index_interval->max_computation();
      value_index_interval->mean_computation();
      value_index_interval->variance_computation();
    }
  }

  return value_index_interval;
}
