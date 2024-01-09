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
 *       $Id: sequence_characteristics.cpp 18069 2015-04-23 10:52:12Z guedon $
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



#include <sstream>

#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/distribution.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/stat_label.h"

#include "sequences.h"
#include "sequence_label.h"

using namespace std;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe SequenceCharacteristics.
 *
 *  argument : nombre de valeurs.
 *
 *--------------------------------------------------------------*/

SequenceCharacteristics::SequenceCharacteristics(int inb_value)

{
  nb_value = inb_value;

  index_value = NULL;
  explicit_index_value = NULL;

  first_occurrence = NULL;
  recurrence_time = NULL;
  sojourn_time = NULL;
  initial_run = NULL;
  final_run = NULL;

  nb_run = NULL;
  nb_occurrence = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe SequenceCharacteristics avec ajout/suppression
 *  des lois empiriques de temps de sejour initial.
 *
 *  arguments : reference sur un objet SequenceCharacteristics,
 *              flag construction des lois empiriques de temps de sejour initial.
 *
 *--------------------------------------------------------------*/

SequenceCharacteristics::SequenceCharacteristics(const SequenceCharacteristics &characteristics ,
                                                 bool initial_run_flag)

{
  register int i;


  nb_value = characteristics.nb_value;

  index_value = new Curves(*(characteristics.index_value));

  if (characteristics.explicit_index_value) {
    explicit_index_value = new Curves(*(characteristics.explicit_index_value));
  }
  else {
    explicit_index_value = NULL;
  }

  first_occurrence = new FrequencyDistribution*[nb_value];
  for (i = 0;i < nb_value;i++) {
    first_occurrence[i] = new FrequencyDistribution(*(characteristics.first_occurrence[i]));
  }

  recurrence_time = new FrequencyDistribution*[nb_value];
  for (i = 0;i < nb_value;i++) {
    recurrence_time[i] = new FrequencyDistribution(*(characteristics.recurrence_time[i]));
  }

  sojourn_time = new FrequencyDistribution*[nb_value];

  if (initial_run_flag) {
    initial_run = new FrequencyDistribution*[nb_value];
  }
  else {
    initial_run = NULL;
  }

  final_run = new FrequencyDistribution*[nb_value];

  for (i = 0;i < nb_value;i++) {
    if (((characteristics.initial_run) && (initial_run_flag)) ||
        ((!(characteristics.initial_run)) && (!initial_run_flag))) {
      sojourn_time[i] = new FrequencyDistribution(*(characteristics.sojourn_time[i]));
      final_run[i] = new FrequencyDistribution(*(characteristics.final_run[i]));
    }
    else {
      sojourn_time[i] = NULL;
      final_run[i] = NULL;
    }

    if ((characteristics.initial_run) && (initial_run_flag)) {
      initial_run[i] = new FrequencyDistribution(*(characteristics.initial_run[i]));
    }
  }

  if (characteristics.nb_run) {
    nb_run = new FrequencyDistribution*[nb_value];
    for (i = 0;i < nb_value;i++) {
      nb_run[i] = new FrequencyDistribution(*(characteristics.nb_run[i]));
    }
  }
  else {
    nb_run = NULL;
  }

  if (characteristics.nb_occurrence) {
    nb_occurrence = new FrequencyDistribution*[nb_value];
    for (i = 0;i < nb_value;i++) {
      nb_occurrence[i] = new FrequencyDistribution(*(characteristics.nb_occurrence[i]));
    }
  }
  else {
    nb_occurrence = NULL;
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie des caracteristiques d'un echantillon de sequences
 *  pour une variable donnee.
 *
 *  argument : reference sur un objet SequenceCharacteristics.
 *
 *--------------------------------------------------------------*/

void SequenceCharacteristics::copy(const SequenceCharacteristics &characteristics)

{
  register int i;


  nb_value = characteristics.nb_value;

  index_value = new Curves(*(characteristics.index_value));

  if (characteristics.explicit_index_value) {
    explicit_index_value = new Curves(*(characteristics.explicit_index_value));
  }
  else {
    explicit_index_value = NULL;
  }

  first_occurrence = new FrequencyDistribution*[nb_value];
  for (i = 0;i < nb_value;i++) {
    first_occurrence[i] = new FrequencyDistribution(*(characteristics.first_occurrence[i]));
  }

  recurrence_time = new FrequencyDistribution*[nb_value];
  for (i = 0;i < nb_value;i++) {
    recurrence_time[i] = new FrequencyDistribution(*(characteristics.recurrence_time[i]));
  }

  sojourn_time = new FrequencyDistribution*[nb_value];

  if (characteristics.initial_run) {
    initial_run = new FrequencyDistribution*[nb_value];
  }
  else {
    initial_run = NULL;
  }

  final_run = new FrequencyDistribution*[nb_value];

  for (i = 0;i < nb_value;i++) {
    sojourn_time[i] = new FrequencyDistribution(*(characteristics.sojourn_time[i]));
    if (characteristics.initial_run) {
      initial_run[i] = new FrequencyDistribution(*(characteristics.initial_run[i]));
    }
    final_run[i] = new FrequencyDistribution(*(characteristics.final_run[i]));
  }

  if (characteristics.nb_run) {
    nb_run = new FrequencyDistribution*[nb_value];
    for (i = 0;i < nb_value;i++) {
      nb_run[i] = new FrequencyDistribution(*(characteristics.nb_run[i]));
    }
  }
  else {
    nb_run = NULL;
  }

  if (characteristics.nb_occurrence) {
    nb_occurrence = new FrequencyDistribution*[nb_value];
    for (i = 0;i < nb_value;i++) {
      nb_occurrence[i] = new FrequencyDistribution(*(characteristics.nb_occurrence[i]));
    }
  }
  else {
    nb_occurrence = NULL;
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie des caracteristiques inchangees d'un echantillon de sequences
 *  pour une variable donnee dans le cas d'une inversion des sequences.
 *
 *  argument : reference sur un objet SequenceCharacteristics.
 *
 *--------------------------------------------------------------*/

void SequenceCharacteristics::reverse(const SequenceCharacteristics &characteristics)

{
  register int i;


  nb_value = characteristics.nb_value;

  index_value = NULL;
  explicit_index_value = NULL;
  first_occurrence = NULL;

  recurrence_time = new FrequencyDistribution*[nb_value];
  for (i = 0;i < nb_value;i++) {
    recurrence_time[i] = new FrequencyDistribution(*(characteristics.recurrence_time[i]));
  }

  if (characteristics.initial_run) {
    sojourn_time = new FrequencyDistribution*[nb_value];
    initial_run = new FrequencyDistribution*[nb_value];
    final_run = new FrequencyDistribution*[nb_value];

    for (i = 0;i < nb_value;i++) {
      sojourn_time[i] = new FrequencyDistribution(*(characteristics.sojourn_time[i]));
      initial_run[i] = new FrequencyDistribution(*(characteristics.final_run[i]));
      final_run[i] = new FrequencyDistribution(*(characteristics.initial_run[i]));
    }
  }

  else {
    sojourn_time = NULL;
    initial_run = NULL;
    final_run = NULL;
  }

  if (characteristics.nb_run) {
    nb_run = new FrequencyDistribution*[nb_value];
    for (i = 0;i < nb_value;i++) {
      nb_run[i] = new FrequencyDistribution(*(characteristics.nb_run[i]));
    }
  }
  else {
    nb_run = NULL;
  }

  if (characteristics.nb_occurrence) {
    nb_occurrence = new FrequencyDistribution*[nb_value];
    for (i = 0;i < nb_value;i++) {
      nb_occurrence[i] = new FrequencyDistribution(*(characteristics.nb_occurrence[i]));
    }
  }
  else {
    nb_occurrence = NULL;
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par copie de la classe SequenceCharacteristics.
 *
 *  arguments : reference sur un objet SequenceCharacteristics,
 *              type de transformation ('r' : inversion).
 *
 *--------------------------------------------------------------*/

SequenceCharacteristics::SequenceCharacteristics(const SequenceCharacteristics &characteristics ,
                                                 char transform)

{
  switch (transform) {
  case 'r' :
    reverse(characteristics);
    break;
  default :
    copy(characteristics);
    break;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur des champs de la classe SequenceCharacteristics.
 *
 *--------------------------------------------------------------*/

void SequenceCharacteristics::remove()

{
  register int i;


  delete index_value;
  delete explicit_index_value;

  if (first_occurrence) {
    for (i = 0;i < nb_value;i++) {
      delete first_occurrence[i];
    }
    delete [] first_occurrence;
  }

  if (recurrence_time) {
    for (i = 0;i < nb_value;i++) {
      delete recurrence_time[i];
    }
    delete [] recurrence_time;
  }

  if (sojourn_time) {
    for (i = 0;i < nb_value;i++) {
      delete sojourn_time[i];
    }
    delete [] sojourn_time;
  }

  if (initial_run) {
    for (i = 0;i < nb_value;i++) {
      delete initial_run[i];
    }
    delete [] initial_run;
  }

  if (final_run) {
    for (i = 0;i < nb_value;i++) {
      delete final_run[i];
    }
    delete [] final_run;
  }

  if (nb_run) {
    for (i = 0;i < nb_value;i++) {
      delete nb_run[i];
    }
    delete [] nb_run;
  }

  if (nb_occurrence) {
    for (i = 0;i < nb_value;i++) {
      delete nb_occurrence[i];
    }
    delete [] nb_occurrence;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe SequenceCharacteristics.
 *
 *--------------------------------------------------------------*/

SequenceCharacteristics::~SequenceCharacteristics()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe SequenceCharacteristics.
 *
 *  argument : reference sur un objet SequenceCharacteristics.
 *
 *--------------------------------------------------------------*/

SequenceCharacteristics& SequenceCharacteristics::operator=(const SequenceCharacteristics &characteristics)

{
  if (&characteristics != this) {
    remove();
    copy(characteristics);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Creation des lois empiriques de temps de sejour dans une valeur
 *  pour une variable donnee.
 *
 *  arguments : longueur maximum des sequences, flag sur la creation
 *              des lois empiriques de temps de sejour initial.
 *
 *--------------------------------------------------------------*/

void SequenceCharacteristics::create_sojourn_time_frequency_distribution(int max_length , int initial_run_flag)

{
  register int i;


  sojourn_time = new FrequencyDistribution*[nb_value];
  for (i = 0;i < nb_value;i++) {
    sojourn_time[i] = new FrequencyDistribution(max_length + 1);
  }

  if (initial_run_flag) {
    initial_run = new FrequencyDistribution*[nb_value];
    for (i = 0;i < nb_value;i++) {
      initial_run[i] = new FrequencyDistribution(max_length + 1);
    }
  }

  final_run = new FrequencyDistribution*[nb_value];
  for (i = 0;i < nb_value;i++) {
    final_run[i] = new FrequencyDistribution(max_length + 1);
  }
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet SequenceCharacteristics.
 *
 *  arguments : stream, type de la variable, loi empirique
 *              des longueurs des sequences, flag niveau de detail, flag fichier.
 *
 *--------------------------------------------------------------*/

ostream& SequenceCharacteristics::ascii_print(ostream &os , int type ,
                                              const FrequencyDistribution &length_distribution ,
                                              bool exhaustive , bool comment_flag) const

{
  register int i;


  if (exhaustive) {
    os << "\n";
    if (comment_flag) {
      os << "# ";
    }
    os << "  ";
    for (i = 0;i < nb_value;i++) {
      os << " | " << SEQ_label[SEQL_OBSERVED] << " "
         << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << i;
    }
    os << " | " << STAT_label[STATL_FREQUENCY] << endl;

    index_value->ascii_print(os , comment_flag);

    if (explicit_index_value) {
      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_INDEX_PARAMETER];
      for (i = 0;i < nb_value;i++) {
        os << " | " << SEQ_label[SEQL_OBSERVED] << " "
           << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << i;
      }
      os << " | " << STAT_label[STATL_FREQUENCY] << endl;

      explicit_index_value->ascii_print(os , comment_flag);
    }
  }

  for (i = 0;i < nb_value;i++) {
    os << "\n";
    if (comment_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_FIRST_OCCURRENCE_OF] << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE]
       << " " << i << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    first_occurrence[i]->ascii_characteristic_print(os , false , comment_flag);

    if ((first_occurrence[i]->nb_element > 0) && (exhaustive)) {
      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << "   | " << SEQ_label[SEQL_FIRST_OCCURRENCE_OF] << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE]
         << " " << i << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
      first_occurrence[i]->ascii_print(os , comment_flag);
    }
  }

  for (i = 0;i < nb_value;i++) {
    os << "\n";
    if (comment_flag) {
      os << "# ";
    }
    os << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << i << " "
       << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    recurrence_time[i]->ascii_characteristic_print(os , false , comment_flag);

    if ((recurrence_time[i]->nb_element > 0) && (exhaustive)) {
      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << "   | " << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << i << " "
         << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
      recurrence_time[i]->ascii_print(os , comment_flag);
    }
  }

  for (i = 0;i < nb_value;i++) {
    os << "\n";
    if (comment_flag) {
      os << "# ";
    }
    os << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << i << " "
       << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    sojourn_time[i]->ascii_characteristic_print(os , false , comment_flag);

    if ((sojourn_time[i]->nb_element > 0) && (exhaustive)) {
      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << "   | " << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << i << " "
         << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
      sojourn_time[i]->ascii_print(os , comment_flag);
    }

    if (initial_run) {
      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_INITIAL_RUN] << " - "
         << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << i << " "
         << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
      initial_run[i]->ascii_characteristic_print(os , false , comment_flag);

      if ((initial_run[i]->nb_element > 0) && (exhaustive)) {
        os << "\n";
        if (comment_flag) {
          os << "# ";
        }
        os << "   | " << SEQ_label[SEQL_INITIAL_RUN] << " - "
           << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << i << " "
           << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
        initial_run[i]->ascii_print(os , comment_flag);
      }
    }

    os << "\n";
    if (comment_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_FINAL_RUN] << " - "
       << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << i << " "
       << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    final_run[i]->ascii_characteristic_print(os , false , comment_flag);

    if ((final_run[i]->nb_element > 0) && (exhaustive)) {
      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << "   | " << SEQ_label[SEQL_FINAL_RUN] << " - "
         << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << i << " "
         << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
      final_run[i]->ascii_print(os , comment_flag);
    }
  }

  if (nb_run) {
    for (i = 0;i < nb_value;i++) {
      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_NB_RUN_OF] << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE]
         << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
      nb_run[i]->ascii_characteristic_print(os , (length_distribution.variance > 0. ? false : true) , comment_flag);

      if ((nb_run[i]->nb_element > 0) && (exhaustive)) {
        os << "\n";
        if (comment_flag) {
          os << "# ";
        }
        os << "   | " << SEQ_label[SEQL_NB_RUN_OF] << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE]
           << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
        nb_run[i]->ascii_print(os , comment_flag);
      }
    }
  }

  if (nb_occurrence) {
    for (i = 0;i < nb_value;i++) {
      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_NB_OCCURRENCE_OF] << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE]
         << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
      nb_occurrence[i]->ascii_characteristic_print(os , (length_distribution.variance > 0. ? false : true) , comment_flag);

      if ((nb_occurrence[i]->nb_element > 0) && (exhaustive)) {
        os << "\n";
        if (comment_flag) {
          os << "# ";
        }
        os << "   | " << SEQ_label[SEQL_NB_OCCURRENCE_OF] << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE]
           << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
        nb_occurrence[i]->ascii_print(os , comment_flag);
      }
    }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet SequenceCharacteristics au format tableur.
 *
 *  arguments : stream, type de la variable, loi empirique
 *              des longueurs des sequences.
 *
 *--------------------------------------------------------------*/

ostream& SequenceCharacteristics::spreadsheet_print(ostream &os , int type ,
                                                    const FrequencyDistribution &length_distribution) const

{
  register int i;
  Curves *smoothed_curves;


  os << "\n";
  for (i = 0;i < nb_value;i++) {
    os << "\t" << SEQ_label[SEQL_OBSERVED] << " "
       << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << i;
  }
  os << "\t" << STAT_label[STATL_FREQUENCY] << endl;
  index_value->spreadsheet_print(os);

  smoothed_curves = new Curves(*index_value , 's');

  os << "\n" << SEQ_label[SEQL_SMOOTHED_OBSERVED_PROBABILITIES] << endl;
  for (i = 0;i < nb_value;i++) {
    os << "\t" << SEQ_label[SEQL_OBSERVED] << " "
       << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << i;
  }
  os << "\t" << STAT_label[STATL_FREQUENCY] << endl;
  smoothed_curves->spreadsheet_print(os);

  delete smoothed_curves;

  if (explicit_index_value) {
    os << "\n" << SEQ_label[SEQL_INDEX_PARAMETER];
    for (i = 0;i < nb_value;i++) {
      os << "\t" << SEQ_label[SEQL_OBSERVED] << " "
         << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << i;
    }
    os << "\t" << STAT_label[STATL_FREQUENCY] << endl;
    explicit_index_value->spreadsheet_print(os);

    smoothed_curves = new Curves(*explicit_index_value , 's');

    os << "\n" << SEQ_label[SEQL_SMOOTHED_OBSERVED_PROBABILITIES] << endl;
    os << SEQ_label[SEQL_INDEX_PARAMETER];
    for (i = 0;i < nb_value;i++) {
      os << "\t" << SEQ_label[SEQL_OBSERVED] << " "
         << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << i;
    }
    os << "\t" << STAT_label[STATL_FREQUENCY] << endl;
    smoothed_curves->spreadsheet_print(os);

    delete smoothed_curves;
  }

  for (i = 0;i < nb_value;i++) {
    os << "\n" << SEQ_label[SEQL_FIRST_OCCURRENCE_OF] << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE]
       << " " << i << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
    first_occurrence[i]->spreadsheet_characteristic_print(os);

    if (first_occurrence[i]->nb_element > 0) {
      os << "\n\t" << SEQ_label[SEQL_FIRST_OCCURRENCE_OF] << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE]
         << " " << i << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
      first_occurrence[i]->spreadsheet_print(os);
    }
  }

  for (i = 0;i < nb_value;i++) {
    os << "\n" << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << i << " "
       << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
    recurrence_time[i]->spreadsheet_characteristic_print(os);

    if (recurrence_time[i]->nb_element > 0) {
      os << "\n\t" << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << i << " "
         << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
      recurrence_time[i]->spreadsheet_print(os);
    }
  }

  for (i = 0;i < nb_value;i++) {
    os << "\n" << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << i << " "
       << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
    sojourn_time[i]->spreadsheet_characteristic_print(os);

    if (sojourn_time[i]->nb_element > 0) {
      os << "\n\t" << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << i << " "
         << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
      sojourn_time[i]->spreadsheet_print(os);
    }

    if (initial_run) {
      os << "\n" << SEQ_label[SEQL_INITIAL_RUN] << " - "
         << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << i << " "
         << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
      initial_run[i]->spreadsheet_characteristic_print(os);

      if (initial_run[i]->nb_element > 0) {
        os << "\n\t" << SEQ_label[SEQL_INITIAL_RUN] << " - "
           << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << i << " "
           << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
        initial_run[i]->spreadsheet_print(os);
      }
    }

    os << "\n" << SEQ_label[SEQL_FINAL_RUN] << " - "
       << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << i << " "
       << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
    final_run[i]->spreadsheet_characteristic_print(os);

    if (final_run[i]->nb_element > 0) {
      os << "\n\t" << SEQ_label[SEQL_FINAL_RUN] << " - "
         << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << i << " "
         << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
      final_run[i]->spreadsheet_print(os);
    }
  }

  if (nb_run) {
    for (i = 0;i < nb_value;i++) {
      os << "\n" << SEQ_label[SEQL_NB_RUN_OF] << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE]
         << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
      nb_run[i]->spreadsheet_characteristic_print(os , (length_distribution.variance > 0. ? false : true));

      if (nb_run[i]->nb_element > 0) {
        os << "\n\t" << SEQ_label[SEQL_NB_RUN_OF] << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE]
           << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
        nb_run[i]->spreadsheet_print(os);
      }
    }
  }

  if (nb_occurrence) {
    for (i = 0;i < nb_value;i++) {
      os << "\n" << SEQ_label[SEQL_NB_OCCURRENCE_OF] << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE]
         << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
      nb_occurrence[i]->spreadsheet_characteristic_print(os , (length_distribution.variance > 0. ? false : true));

      if (nb_occurrence[i]->nb_element > 0) {
        os << "\n\t" << SEQ_label[SEQL_NB_OCCURRENCE_OF] << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE]
           << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
        nb_occurrence[i]->spreadsheet_print(os);
      }
    }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet Sequences_characteristics.
 *
 *  arguments : prefixe des fichiers, titre des figures, indice de la variable,
 *              nombre de variables, type de la variable,
 *              loi empirique des longueurs des sequences.
 *
 *--------------------------------------------------------------*/

bool SequenceCharacteristics::plot_print(const char *prefix , const char *title ,
                                         int variable , int nb_variable , int type ,
                                         const FrequencyDistribution &length_distribution) const

{
  bool status , start;
  register int i , j , k;
  int index_length , nb_histo , histo_index;
  Curves *smoothed_curves;
  const FrequencyDistribution **phisto;
  ostringstream data_file_name[3];


  // ecriture des fichiers de donnees

  data_file_name[0] << prefix << variable + 1 << 0 << ".dat";

  index_length = index_value->plot_length_computation();

  if (index_value->frequency[index_length - 1] < MAX_FREQUENCY) {
    smoothed_curves = new Curves(*index_value , 's');
  }
  else {
    smoothed_curves = NULL;
  }

  status = index_value->plot_print((data_file_name[0].str()).c_str() ,
                                   index_length , smoothed_curves);
  delete smoothed_curves;

  if (explicit_index_value) {
    data_file_name[2] << prefix << variable + 1 << 2 << ".dat";
    status = explicit_index_value->plot_print((data_file_name[2].str()).c_str());
  }

  if (status) {
    phisto = new const FrequencyDistribution*[1 + NB_OUTPUT * 6];

    data_file_name[1] << prefix << variable + 1 << 1 << ".dat";

    nb_histo = 0;
    for (i = 0;i < nb_value;i++) {
      if (first_occurrence[i]->nb_element > 0) {
        phisto[nb_histo++] = first_occurrence[i];
      }
    }

    for (i = 0;i < nb_value;i++) {
      if (recurrence_time[i]->nb_element > 0) {
        phisto[nb_histo++] = recurrence_time[i];
      }
    }

    for (i = 0;i < nb_value;i++) {
      if (sojourn_time[i]->nb_element > 0) {
        phisto[nb_histo++] = sojourn_time[i];
      }
      if ((initial_run) && (initial_run[i]->nb_element > 0)) {
        phisto[nb_histo++] = initial_run[i];
      }
      if (final_run[i]->nb_element > 0) {
        phisto[nb_histo++] = final_run[i];
      }
    }

    if ((nb_run) && (nb_occurrence)) {
      for (i = 0;i < nb_value;i++) {
        if ((nb_run[i]->nb_element > 0) && (nb_occurrence[i]->nb_element > 0)) {
          phisto[nb_histo++] = nb_run[i];
          phisto[nb_histo++] = nb_occurrence[i];
        }
      }
    }

    length_distribution.plot_print((data_file_name[1].str()).c_str() , nb_histo , phisto);

    // ecriture des fichiers de commandes et des fichiers d'impression

    for (i = 0;i < 2;i++) {
      ostringstream file_name[2];

      switch (i) {
      case 0 :
        file_name[0] << prefix << variable + 1 << 1 << ".plot";
        break;
      case 1 :
        file_name[0] << prefix << variable + 1 << 1 << ".print";
        break;
      }

      ofstream out_file((file_name[0].str()).c_str());

      if (i == 1) {
        out_file << "set terminal postscript" << endl;
        file_name[1] << label(prefix) << variable + 1 << 1 << ".ps";
        out_file << "set output \"" << file_name[1].str() << "\"\n\n";
      }

      out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n";

      if (index_value->frequency[index_length - 1] < MAX_FREQUENCY) {
        out_file << "set title" << " \"";
        if (title) {
          out_file << title << " - ";
        }
        if (nb_variable > 1) {
          out_file << STAT_label[STATL_VARIABLE] << " " << variable + 1 << " - ";
        }
        out_file << SEQ_label[SEQL_SMOOTHED_OBSERVED_PROBABILITIES] << "\"\n\n";

        if (index_length - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }

        j = nb_value + 1;

        out_file << "plot [0:" << index_length - 1 << "] [0:1] ";
        for (k = 0;k < nb_value;k++) {
          out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                   << j++ << " title \"" << SEQ_label[SEQL_OBSERVED] << " "
                   << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " "
                   << k << "\" with linespoints";
          if (k < nb_value - 1) {
            out_file << ",\\";
          }
          out_file << endl;
        }

        if (index_length - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }

        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;
      }

      out_file << "set title \"";
      if (title) {
        out_file << title;
        if (nb_variable > 1) {
          out_file << " - ";
        }
      }
      if (nb_variable > 1) {
        out_file << STAT_label[STATL_VARIABLE] << " " << variable + 1;
      }
      out_file << "\"\n\n";

      if (index_length - 1 < TIC_THRESHOLD) {
        out_file << "set xtics 0,1" << endl;
      }

      out_file << "plot [0:" << index_length - 1 << "] [0:1] ";
      for (j = 0;j < nb_value;j++) {
        out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                 << j + 1 << " title \"" << SEQ_label[SEQL_OBSERVED] << " "
                 << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " "
                 << j << "\" with linespoints";
        if (j < nb_value - 1) {
          out_file << ",\\";
        }
        out_file << endl;
      }

      if (index_length - 1 < TIC_THRESHOLD) {
        out_file << "set xtics autofreq" << endl;
      }

      if (i == 0) {
        out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
      }
      out_file << endl;

      if (explicit_index_value) {
        out_file << "set title \"";
        if (title) {
          out_file << title;
          if (nb_variable > 1) {
            out_file << " - ";
          }
        }
        if (nb_variable > 1) {
          out_file << STAT_label[STATL_VARIABLE] << " " << variable + 1;
        }
        out_file << "\"\n\n";

        out_file << "set xlabel \"" << SEQ_label[SEQL_INDEX] << "\"" << endl;
        if (explicit_index_value->index_parameter[explicit_index_value->length - 1] - explicit_index_value->index_parameter[0] < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }

        out_file << "plot [" << explicit_index_value->index_parameter[0] << ":"
                 << explicit_index_value->index_parameter[explicit_index_value->length - 1] << "] [0:1] ";
        for (j = 0;j < nb_value;j++) {
          out_file << "\"" << label((data_file_name[2].str()).c_str()) << "\" using 1:"
                   << j + 2 << " title \"" << SEQ_label[SEQL_OBSERVED] << " "
                   << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " "
                   << j << "\" with linespoints";
          if (j < nb_value - 1) {
            out_file << ",\\";
          }
          out_file << endl;
        }

        if (explicit_index_value->index_parameter[explicit_index_value->length - 1] - explicit_index_value->index_parameter[0] < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }
        out_file << "set xlabel" << endl;

        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;
      }

      if (length_distribution.nb_value - 1 < TIC_THRESHOLD) {
        out_file << "set xtics 0,1" << endl;
      }
      if ((int)(length_distribution.max * YSCALE) + 1 < TIC_THRESHOLD) {
        out_file << "set ytics 0,1" << endl;
      }

      out_file << "plot [0:" << length_distribution.nb_value - 1 << "] [0:"
               << (int)(length_distribution.max * YSCALE) + 1 << "] \"" 
               << label((data_file_name[1].str()).c_str()) << "\" using 1 title \""
               << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
               << "\" with impulses" << endl;

      if (length_distribution.nb_value - 1 < TIC_THRESHOLD) {
        out_file << "set xtics autofreq" << endl;
      }
      if ((int)(length_distribution.max * YSCALE) + 1 < TIC_THRESHOLD) {
        out_file << "set ytics autofreq" << endl;
      }

      if (i == 1) {
        out_file << "\nset terminal x11" << endl;
      }

      out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
    }

    histo_index = 2;

    for (i = 0;i < 2;i++) {
      ostringstream file_name[2];

      switch (i) {
      case 0 :
        file_name[0] << prefix << variable + 1 << 2 << ".plot";
        break;
      case 1 :
        file_name[0] << prefix << variable + 1 << 2 << ".print";
        break;
      }

      ofstream out_file((file_name[0].str()).c_str());

      if (i == 1) {
        out_file << "set terminal postscript" << endl;
        file_name[1] << label(prefix) << variable + 1 << 2 << ".ps";
        out_file << "set output \"" << file_name[1].str() << "\"\n\n";
      }

      out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
               << "set title";
      if ((title) || (nb_variable > 1)) {
        out_file << " \"";
        if (title) {
          out_file << title;
          if (nb_variable > 1) {
            out_file << " - ";
          }
        }
        if (nb_variable > 1) {
          out_file << STAT_label[STATL_VARIABLE] << " " << variable + 1;
        }
        out_file << "\"";
      }
      out_file << "\n\n";

      j = histo_index;

      start = true;
      for (k = 0;k < nb_value;k++) {
        if (first_occurrence[k]->nb_element > 0) {
          if (!start) {
            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;
          }
          else {
            start = false;
          }

          if (MAX(1 , first_occurrence[k]->nb_value - 1) < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }
          if ((int)(first_occurrence[k]->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics 0,1" << endl;
          }

          out_file << "plot [0:" << MAX(first_occurrence[k]->nb_value - 1 , 1) << "] [0:"
                   << (int)(first_occurrence[k]->max * YSCALE) + 1 << "] \""
                   << label((data_file_name[1].str()).c_str()) << "\" using " << j++
                   << " title \"" << SEQ_label[SEQL_FIRST_OCCURRENCE_OF]
                   << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << k
                   << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses" << endl;

          if (MAX(1 , first_occurrence[k]->nb_value - 1) < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }
          if ((int)(first_occurrence[k]->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics autofreq" << endl;
          }
        }
      }

      if (i == 1) {
        out_file << "\nset terminal x11" << endl;
      }

      out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
    }

    histo_index = j;

    for (i = 0;i < 2;i++) {
      ostringstream file_name[2];

      switch (i) {
      case 0 :
        file_name[0] << prefix << variable + 1 << 3 << ".plot";
        break;
      case 1 :
        file_name[0] << prefix << variable + 1 << 3 << ".print";
        break;
      }

      ofstream out_file((file_name[0].str()).c_str());

      if (i == 1) {
        out_file << "set terminal postscript" << endl;
        file_name[1] << label(prefix) << variable + 1 << 3 << ".ps";
        out_file << "set output \"" << file_name[1].str() << "\"\n\n";
      }

      out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
               << "set title";
      if ((title) || (nb_variable > 1)) {
        out_file << " \"";
        if (title) {
          out_file << title;
          if (nb_variable > 1) {
            out_file << " - ";
          }
        }
        if (nb_variable > 1) {
          out_file << STAT_label[STATL_VARIABLE] << " " << variable + 1;
        }
        out_file << "\"";
      }
      out_file << "\n\n";

      j = histo_index;

      start = true;
      for (k = 0;k < nb_value;k++) {
        if (recurrence_time[k]->nb_element > 0) {
          if (!start) {
            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;
          }
          else {
            start = false;
          }

          if (recurrence_time[k]->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }
          if ((int)(recurrence_time[k]->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics 0,1" << endl;
          }

          out_file << "plot [0:" << recurrence_time[k]->nb_value - 1 << "] [0:"
                   << (int)(recurrence_time[k]->max * YSCALE) + 1 << "] \""
                   << label((data_file_name[1].str()).c_str()) << "\" using " << j++
                   << " title \"" << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE]
                   << " " << k << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                   << "\" with impulses" << endl;

          if (recurrence_time[k]->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }
          if ((int)(recurrence_time[k]->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics autofreq" << endl;
          }
        }
      }

      if (i == 1) {
        out_file << "\nset terminal x11" << endl;
      }

      out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
    }

    histo_index = j;

    for (i = 0;i < 2;i++) {
      ostringstream file_name[2];

      switch (i) {
      case 0 :
        file_name[0] << prefix << variable + 1 << 4 << ".plot";
        break;
      case 1 :
        file_name[0] << prefix << variable + 1 << 4 << ".print";
        break;
      }

      ofstream out_file((file_name[0].str()).c_str());

      if (i == 1) {
        out_file << "set terminal postscript" << endl;
        file_name[1] << label(prefix) << variable + 1 << 4 << ".ps";
        out_file << "set output \"" << file_name[1].str() << "\"\n\n";
      }

      out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
               << "set title";
      if ((title) || (nb_variable > 1)) {
        out_file << " \"";
        if (title) {
          out_file << title;
          if (nb_variable > 1) {
            out_file << " - ";
          }
        }
        if (nb_variable > 1) {
          out_file << STAT_label[STATL_VARIABLE] << " " << variable + 1;
        }
        out_file << "\"";
      }
      out_file << "\n\n";

      j = histo_index;

      start = true;
      for (k = 0;k < nb_value;k++) {
        if (sojourn_time[k]->nb_element > 0) {
          if (!start) {
            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;
          }
          else {
            start = false;
          }

          if (sojourn_time[k]->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }
          if ((int)(sojourn_time[k]->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics 0,1" << endl;
          }

          out_file << "plot [0:" << sojourn_time[k]->nb_value - 1 << "] [0:"
                   << (int)(sojourn_time[k]->max * YSCALE) + 1 << "] \""
                   << label((data_file_name[1].str()).c_str()) << "\" using " << j++
                   << " title \"" << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << k
                   << " " << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                   << "\" with impulses" << endl;

          if (sojourn_time[k]->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }
          if ((int)(sojourn_time[k]->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics autofreq" << endl;
          }
        }

        if ((initial_run) && (initial_run[k]->nb_element > 0)) {
          if (!start) {
            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;
          }
          else {
            start = false;
          }

          if (initial_run[k]->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }
          if ((int)(initial_run[k]->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics 0,1" << endl;
          }

          out_file << "plot [0:" << initial_run[k]->nb_value - 1 << "] [0:"
                   << (int)(initial_run[k]->max * YSCALE) + 1 << "] \""
                   << label((data_file_name[1].str()).c_str()) << "\" using " << j++
                   << " title \"" << SEQ_label[SEQL_INITIAL_RUN] << " - "
                   << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << k
                   << " " << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                   << "\" with impulses" << endl;

          if (initial_run[k]->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }
          if ((int)(initial_run[k]->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics autofreq" << endl;
          }
        }

        if (final_run[k]->nb_element > 0) {
          if (!start) {
            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;
          }
          else {
            start = false;
          }

          if (final_run[k]->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }
          if ((int)(final_run[k]->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics 0,1" << endl;
          }

          out_file << "plot [0:" << final_run[k]->nb_value - 1 << "] [0:"
                   << (int)(final_run[k]->max * YSCALE) + 1 << "] \""
                   << label((data_file_name[1].str()).c_str()) << "\" using " << j++
                   << " title \"" << SEQ_label[SEQL_FINAL_RUN] << " - "
                   << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << k
                   << " " << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                   << "\" with impulses" << endl;

          if (final_run[k]->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }
          if ((int)(final_run[k]->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics autofreq" << endl;
          }
        }
      }

      if (i == 1) {
        out_file << "\nset terminal x11" << endl;
      }

      out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
    }

    histo_index = j;

    if ((nb_run) && (nb_occurrence)) {
      for (i = 0;i < 2;i++) {
        ostringstream file_name[2];

        switch (i) {
        case 0 :
          file_name[0] << prefix << variable + 1 << 5 << ".plot";
          break;
        case 1 :
          file_name[0] << prefix << variable + 1 << 5 << ".print";
          break;
        }

        ofstream out_file((file_name[0].str()).c_str());

        if (i == 1) {
          out_file << "set terminal postscript" << endl;
          file_name[1] << label(prefix) << variable + 1 << 5 << ".ps";
          out_file << "set output \"" << file_name[1].str() << "\"\n\n";
        }

        out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                 << "set title";
        if ((title) || (nb_variable > 1)) {
          out_file << " \"";
          if (title) {
            out_file << title;
            if (nb_variable > 1) {
              out_file << " - ";
            }
          }
          if (nb_variable > 1) {
            out_file << STAT_label[STATL_VARIABLE] << " " << variable + 1;
          }
          out_file << "\"";
        }
        out_file << "\n\n";

        j = histo_index;

        start = true;
        for (k = 0;k < nb_value;k++) {
          if ((nb_run[k]->nb_element > 0) && (nb_occurrence[k]->nb_element > 0)) {
            if (!start) {
              if (i == 0) {
                out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
              }
              out_file << endl;
            }
            else {
              start = false;
            }

            if (nb_run[k]->nb_value - 1 < TIC_THRESHOLD) {
              out_file << "set xtics 0,1" << endl;
            }
            if ((int)(nb_run[k]->max * YSCALE) + 1 < TIC_THRESHOLD) {
              out_file << "set ytics 0,1" << endl;
            }

            out_file << "plot [0:" << nb_run[k]->nb_value - 1 << "] [0:"
                     << (int)(nb_run[k]->max * YSCALE) + 1 << "] \""
                     << label((data_file_name[1].str()).c_str()) << "\" using " << j++
                     << " title \"" << SEQ_label[SEQL_NB_RUN_OF]
                     << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << k
                     << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                     << "\" with impulses" << endl;

            if (nb_run[k]->nb_value - 1 < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }
            if ((int)(nb_run[k]->max * YSCALE) + 1 < TIC_THRESHOLD) {
              out_file << "set ytics autofreq" << endl;
            }

            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;

            if (nb_occurrence[k]->nb_value - 1 < TIC_THRESHOLD) {
              out_file << "set xtics 0,1" << endl;
            }
            if ((int)(nb_occurrence[k]->max * YSCALE) + 1 < TIC_THRESHOLD) {
              out_file << "set ytics 0,1" << endl;
            }

            out_file << "plot [0:" << nb_occurrence[k]->nb_value - 1 << "] [0:"
                     << (int)(nb_occurrence[k]->max * YSCALE) + 1 << "] \""
                     << label((data_file_name[1].str()).c_str()) << "\" using " << j++
                     << " title \"" << SEQ_label[SEQL_NB_OCCURRENCE_OF]
                     << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << k
                     << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                     << "\" with impulses" << endl;

            if (nb_occurrence[k]->nb_value - 1 < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }
            if ((int)(nb_occurrence[k]->max * YSCALE) + 1 < TIC_THRESHOLD) {
              out_file << "set ytics autofreq" << endl;
            }
          }
        }

        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;

        if (length_distribution.nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }
        if ((int)(length_distribution.max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics 0,1" << endl;
        }

        out_file << "plot [0:" << length_distribution.nb_value - 1 << "] [0:"
                 << (int)(length_distribution.max * YSCALE) + 1 << "] \"" 
                 << label((data_file_name[1].str()).c_str()) << "\" using 1 title \""
                 << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                 << "\" with impulses" << endl;

        if (length_distribution.nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }
        if ((int)(length_distribution.max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics autofreq" << endl;
        }

        if (i == 1) {
          out_file << "\nset terminal x11" << endl;
        }

        out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
      }
    }

    delete [] phisto;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet Sequences_characteristics.
 *
 *  arguments : reference sur un objet MultiPlotSet, indice du MultiPlot,
 *              indice et type de la variable,
 *              loi empirique des longueurs des sequences.
 *
 *--------------------------------------------------------------*/

void SequenceCharacteristics::plotable_write(MultiPlotSet &plot , int &index ,
                                             int variable , int type ,
                                             const FrequencyDistribution &length_distribution) const

{
  register int i , j , k;
  int index_length , nb_histo , max_nb_value , max_frequency;
  double shift;
  Curves *smoothed_curves;
  ostringstream title , legend;


  index_length = index_value->plot_length_computation();

  // calcul du nombre de vues

  /*  nb_plot_set = 2;
  if (index_value->frequency[index_length - 1] < MAX_FREQUENCY) {
    nb_plot_set++;
  }

  nb_plot_set++;
  for (i = 0;i < nb_value;i++) {
    if (first_occurrence[i]->nb_element > 0) {
      nb_plot_set++;
    }
  }

  nb_plot_set++;
  for (i = 0;i < nb_value;i++) {
    if (recurrence_time[i]->nb_element > 0) {
      nb_plot_set++;
    }
  }

  nb_plot_set++;
  for (i = 0;i < nb_value;i++) {
    if (sojourn_time[i]->nb_element > 0) {
      nb_plot_set++;
    }
    if ((initial_run) && (initial_run[i]->nb_element > 0)) {
      nb_plot_set++;
    }
    if (final_run[i]->nb_element > 0) {
      nb_plot_set++;
    }
  }

  if ((nb_run) && (nb_occurrence)) {
    nb_plot_set += 3;
    for (i = 0;i < nb_value;i++) {
      if ((nb_run[i]->nb_element > 0) && (nb_occurrence[i]->nb_element > 0)) {
        nb_plot_set += 2;
      }
    }
  } */

  plot.variable_nb_viewpoint[variable] += 4;
  if ((nb_run) && (nb_occurrence)) {
    plot.variable_nb_viewpoint[variable]++;
  }

  if (index_value->frequency[index_length - 1] < MAX_FREQUENCY) {

    // vue : intensite lissee

    plot.variable[index] = variable;
    plot.viewpoint[index] = INTENSITY;

    smoothed_curves = new Curves(*index_value , 's');

    title.str("");
    if (plot.nb_variable > 1) {
      title << STAT_label[STATL_VARIABLE] << " " << variable + 1 << " - ";
    }
    title << SEQ_label[SEQL_SMOOTHED_OBSERVED_PROBABILITIES];
    plot[index].title = title.str();

    plot[index].xrange = Range(0 , index_length - 1);
    plot[index].yrange = Range(0. , 1.);

    if (index_length - 1 < TIC_THRESHOLD) {
      plot[index].xtics = 1;
    }

    plot[index].resize(nb_value);

    for (i = 0;i < nb_value;i++) {
      legend.str("");
      legend << SEQ_label[SEQL_OBSERVED] << " "
             << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << i;
      plot[index][i].legend = legend.str();

      plot[index][i].style = "linespoints";
    }

    smoothed_curves->plotable_write(plot[index]);

    delete smoothed_curves;
    index++;
  }

  // vue : intensite

  plot.variable[index] = variable;
  plot.viewpoint[index] = INTENSITY;

  if (plot.nb_variable > 1) {
    title.str("");
    title << STAT_label[STATL_VARIABLE] << " " << variable + 1;
    plot[index].title = title.str();
  }

  plot[index].xrange = Range(0 , index_length - 1);
  plot[index].yrange = Range(0. , 1.);

  if (index_length - 1 < TIC_THRESHOLD) {
    plot[index].xtics = 1;
  }

  plot[index].resize(nb_value);

  for (i = 0;i < nb_value;i++) {
    legend.str("");
    legend << SEQ_label[SEQL_OBSERVED] << " "
           << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << i;
    plot[index][i].legend = legend.str();

    plot[index][i].style = "linespoints";
  }

  index_value->plotable_write(plot[index]);
  index++;

  if (explicit_index_value) {

    // vue : intensite fonction d'un parametre d'index explicite

    plot.variable[index] = variable;
    plot.viewpoint[index] = INTENSITY;

    if (plot.nb_variable > 1) {
      title.str("");
      title << STAT_label[STATL_VARIABLE] << " " << variable + 1;
      plot[index].title = title.str();
    }

    plot[index].xrange = Range(explicit_index_value->index_parameter[0] ,
                               explicit_index_value->index_parameter[explicit_index_value->length - 1]);
    plot[index].yrange = Range(0. , 1.);

    if (explicit_index_value->index_parameter[explicit_index_value->length - 1] - explicit_index_value->index_parameter[0] < TIC_THRESHOLD) {
      plot[index].xtics = 1;
    }
    plot[index].xlabel = SEQ_label[SEQL_INDEX];

    plot[index].resize(nb_value);

    for (i = 0;i < nb_value;i++) {
      legend.str("");
      legend << SEQ_label[SEQL_OBSERVED] << " "
             << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << i;
      plot[index][i].legend = legend.str();

      plot[index][i].style = "linespoints";
    }

    explicit_index_value->plotable_write(plot[index]);
    index++;
  }

  // vue : loi empirique des longueurs des sequences

  plot.variable[index] = variable;
  plot.viewpoint[index] = INTENSITY;

  plot[index].xrange = Range(0 , length_distribution.nb_value - 1);
  plot[index].yrange = Range(0 , ceil(length_distribution.max * YSCALE));

  if (length_distribution.nb_value - 1 < TIC_THRESHOLD) {
    plot[index].xtics = 1;
  }
  if (ceil(length_distribution.max * YSCALE) < TIC_THRESHOLD) {
    plot[index].ytics = 1;
  }

  plot[index].resize(1);

  legend.str("");
  legend << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
  plot[index][0].legend = legend.str();

  plot[index][0].style = "impulses";

  length_distribution.plotable_frequency_write(plot[index][0]);
  index++;

  // vue : lois empiriques du temps avant 1ere occurrence d'une observation

  plot.variable[index] = variable;
  plot.viewpoint[index] = FIRST_OCCURRENCE;

  title.str("");
  if (plot.nb_variable > 1) {
    title << STAT_label[STATL_VARIABLE] << " " << variable + 1 << " - ";
  }
  title << SEQ_label[SEQL_FIRST_OCCURRENCE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTIONS];
  plot[index].title = title.str();

  // calcul du temps maximum et de la frequence maximum

  nb_histo = 0;
  max_nb_value = 0;
  max_frequency = 0;

  for (i = 0;i < nb_value;i++) {
    if (first_occurrence[i]->nb_element > 0) {
      nb_histo++;

      if (first_occurrence[i]->nb_value > max_nb_value) {
        max_nb_value = first_occurrence[i]->nb_value;
      }
      if (first_occurrence[i]->max > max_frequency) {
        max_frequency = first_occurrence[i]->max;
      }
    }
  }

  plot[index].xrange = Range(0 , max_nb_value);
  plot[index].yrange = Range(0 , ceil(max_frequency * YSCALE));

  if (max_nb_value < TIC_THRESHOLD) {
    plot[index].xtics = 1;
  }
  if (ceil(max_frequency * YSCALE) < TIC_THRESHOLD) {
    plot[index].ytics = 1;
  }

  plot[index].resize(nb_histo);

  i = 0;
  shift = 0.;

  for (j = 0;j < nb_value;j++) {
    if (first_occurrence[j]->nb_element > 0) {
      legend.str("");
      legend << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << j;
      plot[index][i].legend = legend.str();

      plot[index][i].style = "impulses";

      for (k = first_occurrence[j]->offset;k < first_occurrence[j]->nb_value;k++) {
        if (first_occurrence[j]->frequency[k] > 0) {
          plot[index][i].add_point(k + shift , first_occurrence[j]->frequency[k]);
        }
      }

      if (PLOT_SHIFT * (nb_histo - 1) < PLOT_MAX_SHIFT) {
        shift += PLOT_SHIFT;
      }
      else {
        shift += PLOT_MAX_SHIFT / (nb_histo - 1);
      }

      i++;
    }
  }
  index++;

  for (i = 0;i < nb_value;i++) {
    if (first_occurrence[i]->nb_element > 0) {

      // vue : loi empirique du temps avant la 1ere occurrence d'une observation

      plot.variable[index] = variable;
      plot.viewpoint[index] = FIRST_OCCURRENCE;

      if (plot.nb_variable > 1) {
        title.str("");
        title << STAT_label[STATL_VARIABLE] << " " << variable + 1;
        plot[index].title = title.str();
      }

      plot[index].xrange = Range(0 , MAX(first_occurrence[i]->nb_value - 1 , 1));
      plot[index].yrange = Range(0 , ceil(first_occurrence[i]->max * YSCALE));

      if (MAX(first_occurrence[i]->nb_value - 1 , 1) < TIC_THRESHOLD) {
        plot[index].xtics = 1;
      }
      if (ceil(first_occurrence[i]->max * YSCALE) < TIC_THRESHOLD) {
        plot[index].ytics = 1;
      }

      plot[index].resize(1);

      legend.str("");
      legend << SEQ_label[SEQL_FIRST_OCCURRENCE_OF] << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE]
             << " " << i << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      plot[index][0].legend = legend.str();

      plot[index][0].style = "impulses";

      first_occurrence[i]->plotable_frequency_write(plot[index][0]);
      index++;
    }
  }

  // vue : lois empiriques du temps de retour dans une observation

  plot.variable[index] = variable;
  plot.viewpoint[index] = RECURRENCE_TIME;

  title.str("");
  if (plot.nb_variable > 1) {
    title << STAT_label[STATL_VARIABLE] << " " << variable + 1 << " - ";
  }
  title << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTIONS];
  plot[index].title = title.str();

  // calcul du temps maximum et de la frequence maximum

  nb_histo = 0;
  max_nb_value = 0;
  max_frequency = 0;

  for (i = 0;i < nb_value;i++) {
    if (recurrence_time[i]->nb_element > 0) {
      nb_histo++;

      if (recurrence_time[i]->nb_value > max_nb_value) {
        max_nb_value = recurrence_time[i]->nb_value;
      }
      if (recurrence_time[i]->max > max_frequency) {
        max_frequency = recurrence_time[i]->max;
      }
    }
  }

  plot[index].xrange = Range(0 , max_nb_value);
  plot[index].yrange = Range(0 , ceil(max_frequency * YSCALE));

  if (max_nb_value < TIC_THRESHOLD) {
    plot[index].xtics = 1;
  }
  if (ceil(max_frequency * YSCALE) < TIC_THRESHOLD) {
    plot[index].ytics = 1;
  }

  plot[index].resize(nb_histo);

  i = 0;
  shift = 0.;

  for (j = 0;j < nb_value;j++) {
    if (recurrence_time[j]->nb_element > 0) {
      legend.str("");
      legend << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << j;
      plot[index][i].legend = legend.str();

      plot[index][i].style = "impulses";

      for (k = recurrence_time[j]->offset;k < recurrence_time[j]->nb_value;k++) {
        if (recurrence_time[j]->frequency[k] > 0) {
          plot[index][i].add_point(k + shift , recurrence_time[j]->frequency[k]);
        }
      }

      if (PLOT_SHIFT * (nb_histo - 1) < PLOT_MAX_SHIFT) {
        shift += PLOT_SHIFT;
      }
      else {
        shift += PLOT_MAX_SHIFT / (nb_histo - 1);
      }

      i++;
    }
  }
  index++;

  for (i = 0;i < nb_value;i++) {
    if (recurrence_time[i]->nb_element > 0) {

      // vue : loi empirique du temps de retour dans une observation

      plot.variable[index] = variable;
      plot.viewpoint[index] = RECURRENCE_TIME;

      if (plot.nb_variable > 1) {
        title.str("");
        title << STAT_label[STATL_VARIABLE] << " " << variable + 1;
        plot[index].title = title.str();
      }

      plot[index].xrange = Range(0 , recurrence_time[i]->nb_value - 1);
      plot[index].yrange = Range(0 , ceil(recurrence_time[i]->max * YSCALE));

      if (recurrence_time[i]->nb_value - 1 < TIC_THRESHOLD) {
        plot[index].xtics = 1;
      }
      if (ceil(recurrence_time[i]->max * YSCALE) < TIC_THRESHOLD) {
        plot[index].ytics = 1;
      }

      plot[index].resize(1);

      legend.str("");
      legend << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " "
             << i << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
             << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      plot[index][0].legend = legend.str();

      plot[index][0].style = "impulses";

      recurrence_time[i]->plotable_frequency_write(plot[index][0]);
      index++;
    }
  }

  // vue : lois empiriques du temps de sejour dans une observation

  plot.variable[index] = variable;
  plot.viewpoint[index] = SOJOURN_TIME;

  title.str("");
  if (plot.nb_variable > 1) {
    title << STAT_label[STATL_VARIABLE] << " " << variable + 1 << " - ";
  }
  title << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTIONS];
  plot[index].title = title.str();

  // calcul du temps maximum et de la frequence maximum

  nb_histo = 0;
  max_nb_value = 0;
  max_frequency = 0;

  for (i = 0;i < nb_value;i++) {
    if (sojourn_time[i]->nb_element > 0) {
      nb_histo++;

      if (sojourn_time[i]->nb_value > max_nb_value) {
        max_nb_value = sojourn_time[i]->nb_value;
      }
      if (sojourn_time[i]->max > max_frequency) {
        max_frequency = sojourn_time[i]->max;
      }
    }
  }

  plot[index].xrange = Range(0 , max_nb_value);
  plot[index].yrange = Range(0 , ceil(max_frequency * YSCALE));

  if (max_nb_value < TIC_THRESHOLD) {
    plot[index].xtics = 1;
  }
  if (ceil(max_frequency * YSCALE) < TIC_THRESHOLD) {
    plot[index].ytics = 1;
  }

  plot[index].resize(nb_histo);

  i = 0;
  shift = 0.;

  for (j = 0;j < nb_value;j++) {
    if (sojourn_time[j]->nb_element > 0) {
      legend.str("");
      legend << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << j;
      plot[index][i].legend = legend.str();

      plot[index][i].style = "impulses";

      for (k = sojourn_time[j]->offset;k < sojourn_time[j]->nb_value;k++) {
        if (sojourn_time[j]->frequency[k] > 0) {
          plot[index][i].add_point(k + shift , sojourn_time[j]->frequency[k]);
        }
      }

      if (PLOT_SHIFT * (nb_histo - 1) < PLOT_MAX_SHIFT) {
        shift += PLOT_SHIFT;
      }
      else {
        shift += PLOT_MAX_SHIFT / (nb_histo - 1);
      }

      i++;
    }
  }
  index++;

  for (i = 0;i < nb_value;i++) {
    if (sojourn_time[i]->nb_element > 0) {

      // vue : loi empirique du temps de sejour dans une observation

      plot.variable[index] = variable;
      plot.viewpoint[index] = SOJOURN_TIME;

      if (plot.nb_variable > 1) {
        title.str("");
        title << STAT_label[STATL_VARIABLE] << " " << variable + 1;
        plot[index].title = title.str();
      }

      plot[index].xrange = Range(0 , sojourn_time[i]->nb_value - 1);
      plot[index].yrange = Range(0 , ceil(sojourn_time[i]->max * YSCALE));

      if (sojourn_time[i]->nb_value - 1 < TIC_THRESHOLD) {
        plot[index].xtics = 1;
      }
      if (ceil(sojourn_time[i]->max * YSCALE) < TIC_THRESHOLD) {
        plot[index].ytics = 1;
      }

      plot[index].resize(1);

      legend.str("");
      legend << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " "
             << i << " " << SEQ_label[SEQL_SOJOURN_TIME] << " "
             << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      plot[index][0].legend = legend.str();

      plot[index][0].style = "impulses";

      sojourn_time[i]->plotable_frequency_write(plot[index][0]);
      index++;
    }

    if ((initial_run) && (initial_run[i]->nb_element > 0)) {

      // vue : loi empirique du temps de sejour dans la 1ere observation rencontree

      plot.variable[index] = variable;
      plot.viewpoint[index] = SOJOURN_TIME;

      if (plot.nb_variable > 1) {
        title.str("");
        title << STAT_label[STATL_VARIABLE] << " " << variable + 1;
        plot[index].title = title.str();
      }

      plot[index].xrange = Range(0 , initial_run[i]->nb_value - 1);
      plot[index].yrange = Range(0 , ceil(initial_run[i]->max * YSCALE));

      if (initial_run[i]->nb_value - 1 < TIC_THRESHOLD) {
        plot[index].xtics = 1;
      }
      if (ceil(initial_run[i]->max * YSCALE) < TIC_THRESHOLD) {
        plot[index].ytics = 1;
      }

      plot[index].resize(1);

      legend.str("");
      legend << SEQ_label[SEQL_INITIAL_RUN] << " - "
             << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " "
             << i << " " << SEQ_label[SEQL_SOJOURN_TIME] << " "
             << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      plot[index][0].legend = legend.str();

      plot[index][0].style = "impulses";

      initial_run[i]->plotable_frequency_write(plot[index][0]);
      index++;
    }

    if (final_run[i]->nb_element > 0) {

      // vue : loi empirique du temps de sejour dans la derniere observation rencontree

      plot.variable[index] = variable;
      plot.viewpoint[index] = SOJOURN_TIME;

      if (plot.nb_variable > 1) {
        title.str("");
        title << STAT_label[STATL_VARIABLE] << " " << variable + 1;
        plot[index].title = title.str();
      }

      plot[index].xrange = Range(0 , final_run[i]->nb_value - 1);
      plot[index].yrange = Range(0 , ceil(final_run[i]->max * YSCALE));

      if (final_run[i]->nb_value - 1 < TIC_THRESHOLD) {
        plot[index].xtics = 1;
      }
      if (ceil(final_run[i]->max * YSCALE) < TIC_THRESHOLD) {
        plot[index].ytics = 1;
      }

      legend.str("");
      legend << SEQ_label[SEQL_FINAL_RUN] << " - "
             << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " "
             << i << " " << SEQ_label[SEQL_SOJOURN_TIME] << " "
             << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      plot[index][0].legend = legend.str();

      plot[index][0].style = "impulses";

      final_run[i]->plotable_frequency_write(plot[index][0]);
      index++;
    }
  }

  if ((nb_run) && (nb_occurrence)) {

    // vue : lois empiriques du nombre de series d'une observation

    plot.variable[index] = variable;
    plot.viewpoint[index] = COUNTING;

    title.str("");
    if (plot.nb_variable > 1) {
      title << STAT_label[STATL_VARIABLE] << " " << variable + 1 << " - ";
    }
    title << SEQ_label[SEQL_NB_RUN] << " " << SEQ_label[SEQL_PER_SEQUENCE] << " "
          << STAT_label[STATL_FREQUENCY_DISTRIBUTIONS];
    plot[index].title = title.str();

    // calcul du nombre de series maximum et de la frequence maximum

    nb_histo = 0;
    max_nb_value = 0;
    max_frequency = 0;

    for (i = 0;i < nb_value;i++) {
      if (nb_run[i]->nb_element > 0) {
        nb_histo++;

        if (nb_run[i]->nb_value > max_nb_value) {
          max_nb_value = nb_run[i]->nb_value;
        }
        if (nb_run[i]->max > max_frequency) {
          max_frequency = nb_run[i]->max;
        }
      }
    }

    plot[index].xrange = Range(0 , max_nb_value);
    plot[index].yrange = Range(0 , ceil(max_frequency * YSCALE));

    if (max_nb_value < TIC_THRESHOLD) {
      plot[index].xtics = 1;
    }
    if (ceil(max_frequency * YSCALE) < TIC_THRESHOLD) {
      plot[index].ytics = 1;
    }

    plot[index].resize(nb_histo);

    i = 0;
    shift = 0.;

    for (j = 0;j < nb_value;j++) {
      if (nb_run[j]->nb_element > 0) {
        legend.str("");
        legend << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << j;
        plot[index][i].legend = legend.str();

        plot[index][i].style = "impulses";

        for (k = nb_run[j]->offset;k < nb_run[j]->nb_value;k++) {
          if (nb_run[j]->frequency[k] > 0) {
            plot[index][i].add_point(k + shift , nb_run[j]->frequency[k]);
          }
        }

        if (PLOT_SHIFT * (nb_histo - 1) < PLOT_MAX_SHIFT) {
          shift += PLOT_SHIFT;
        }
        else {
          shift += PLOT_MAX_SHIFT / (nb_histo - 1);
        }

        i++;
      }
    }
    index++;

    // vue : lois empiriques du nombre d'occurrences d'une observation

    plot.variable[index] = variable;
    plot.viewpoint[index] = COUNTING;

    title.str("");
    if (plot.nb_variable > 1) {
      title << STAT_label[STATL_VARIABLE] << " " << variable + 1 << " - ";
    }
    title << SEQ_label[SEQL_NB_OCCURRENCE] << " " << SEQ_label[SEQL_PER_SEQUENCE] << " "
          << STAT_label[STATL_FREQUENCY_DISTRIBUTIONS];
    plot[index].title = title.str();

    // calcul du nombre d'occurrences maximum et de la frequence maximum

    nb_histo = 0;
    max_nb_value = 0;
    max_frequency = 0;

    for (i = 0;i < nb_value;i++) {
      if (nb_occurrence[i]->nb_element > 0) {
        nb_histo++;

        if (nb_occurrence[i]->nb_value > max_nb_value) {
          max_nb_value = nb_occurrence[i]->nb_value;
        }
        if (nb_occurrence[i]->max > max_frequency) {
          max_frequency = nb_occurrence[i]->max;
        }
      }
    }

    plot[index].xrange = Range(0 , max_nb_value);
    plot[index].yrange = Range(0 , ceil(max_frequency * YSCALE));

    if (max_nb_value < TIC_THRESHOLD) {
      plot[index].xtics = 1;
    }
    if (ceil(max_frequency * YSCALE) < TIC_THRESHOLD) {
      plot[index].ytics = 1;
    }

    plot[index].resize(nb_histo);

    i = 0;
    shift = 0.;

    for (j = 0;j < nb_value;j++) {
      if (nb_occurrence[j]->nb_element > 0) {
        legend.str("");
        legend << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << j;
        plot[index][i].legend = legend.str();

        plot[index][i].style = "impulses";

        for (k = nb_occurrence[j]->offset;k < nb_occurrence[j]->nb_value;k++) {
          if (nb_occurrence[j]->frequency[k] > 0) {
            plot[index][i].add_point(k + shift , nb_occurrence[j]->frequency[k]);
          }
        }

        if (PLOT_SHIFT * (nb_histo - 1) < PLOT_MAX_SHIFT) {
          shift += PLOT_SHIFT;
        }
        else {
          shift += PLOT_MAX_SHIFT / (nb_histo - 1);
        }

        i++;
      }
    }
    index++;

    for (i = 0;i < nb_value;i++) {
      if ((nb_run[i]->nb_element > 0) && (nb_occurrence[i]->nb_element > 0)) {

        // vue : loi empirique du nombre de series d'une observation par sequence

        plot.variable[index] = variable;
        plot.viewpoint[index] = COUNTING;

        if (plot.nb_variable > 1) {
          title.str("");
          title << STAT_label[STATL_VARIABLE] << " " << variable + 1;
          plot[index].title = title.str();
        }

        plot[index].xrange = Range(0 , nb_run[i]->nb_value - 1);
        plot[index].yrange = Range(0 , ceil(nb_run[i]->max * YSCALE));

        if (nb_run[i]->nb_value - 1 < TIC_THRESHOLD) {
          plot[index].xtics = 1;
        }
        if (ceil(nb_run[i]->max * YSCALE) + 1 < TIC_THRESHOLD) {
          plot[index].ytics = 1;
        }

        plot[index].resize(1);

        legend.str("");
        legend << SEQ_label[SEQL_NB_RUN_OF]
               << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << i << " "
               << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        plot[index][0].legend = legend.str();

        plot[index][0].style = "impulses";

        nb_run[i]->plotable_frequency_write(plot[index][0]);
        index++;

        // vue : loi empirique du nombre d'occurrences d'une observation par sequence

        plot.variable[index] = variable;
        plot.viewpoint[index] = COUNTING;

        if (plot.nb_variable > 1) {
          title.str("");
          title << STAT_label[STATL_VARIABLE] << " " << variable + 1;
          plot[index].title = title.str();
        }

        plot[index].xrange = Range(0 , nb_occurrence[i]->nb_value - 1);
        plot[index].yrange = Range(0 , ceil(nb_occurrence[i]->max * YSCALE));

        if (nb_occurrence[i]->nb_value - 1 < TIC_THRESHOLD) {
          plot[index].xtics = 1;
        }
        if (ceil(nb_occurrence[i]->max * YSCALE) < TIC_THRESHOLD) {
          plot[index].ytics = 1;
        }

        plot[index].resize(1);

        legend.str("");
        legend << SEQ_label[SEQL_NB_OCCURRENCE_OF]
               << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << i << " "
               << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        plot[index][0].legend = legend.str();

        plot[index][0].style = "impulses";

        nb_occurrence[i]->plotable_frequency_write(plot[index][0]);
        index++;
      }
    }

    // vue : loi empirique des longueurs des sequences

    plot.variable[index] = variable;
    plot.viewpoint[index] = COUNTING;

    plot[index].xrange = Range(0 , length_distribution.nb_value - 1);
    plot[index].yrange = Range(0 , ceil(length_distribution.max * YSCALE));

    if (length_distribution.nb_value - 1 < TIC_THRESHOLD) {
      plot[index].xtics = 1;
    }
    if (ceil(length_distribution.max * YSCALE) < TIC_THRESHOLD) {
      plot[index].ytics = 1;
    }

    plot[index].resize(1);

    legend.str("");
    legend << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    plot[index][0].legend = legend.str();

    plot[index][0].style = "impulses";

    length_distribution.plotable_frequency_write(plot[index][0]);
    index++;
  }
}


};  // namespace sequence_analysis
