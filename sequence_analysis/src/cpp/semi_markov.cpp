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
 *       $Id: semi_markov.cpp 18067 2015-04-23 10:51:22Z guedon $
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
#include <iomanip>

#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"
#include "tool/config.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/distribution.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/stat_label.h"

#include "sequences.h"
#include "semi_markov.h"
#include "sequence_label.h"

using namespace std;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe SemiMarkovChain.
 *
 *--------------------------------------------------------------*/

SemiMarkovChain::SemiMarkovChain()

{
  state_subtype = NULL;
  state_process = NULL;
  forward = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe SemiMarkovChain.
 *
 *  arguments : type, nombre d'etats.
 *
 *--------------------------------------------------------------*/

SemiMarkovChain::SemiMarkovChain(char itype , int inb_state)
:Chain(itype , inb_state)

{
  state_subtype = NULL;
  state_process = new CategoricalSequenceProcess(nb_state , nb_state , false);
  forward = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe SemiMarkovChain.
 *
 *  arguments : pointeur sur un objet Chain et sur un objet CategoricalSequenceProcess.
 *
 *--------------------------------------------------------------*/

SemiMarkovChain::SemiMarkovChain(const Chain *pchain , const CategoricalSequenceProcess *poccupancy)
:Chain(*pchain)

{
  register int i;


  state_subtype = new int[nb_state];

  state_process = new CategoricalSequenceProcess(*poccupancy);
  for (i = 0;i < nb_state;i++) {
    if (transition[i][i] < 1.) {
      state_process->absorption[i] = 0.;
    }
    else {
      state_process->absorption[i] = 1.;
    }
  }

  forward = new Forward*[nb_state];

  for (i = 0;i < nb_state;i++) {
    state_subtype[i] = (state_process->sojourn_time[i] ? SEMI_MARKOVIAN : MARKOVIAN);

    if ((state_subtype[i] == SEMI_MARKOVIAN) && (state_type[i] == 'r')) {
      forward[i] = new Forward(*(state_process->sojourn_time[i]));
    }
    else {
      forward[i] = NULL;
    }
  }

  if (type == 'e') {
    for (i = 0;i < nb_state;i++) {
      initial[i] = 1. / (double)nb_state;
    }
    initial_probability_computation();
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet SemiMarkovChain.
 *
 *  arguments : reference sur un objet SemiMarkovChain,
 *              parametre (si strictement positif :
 *              nombre de valeurs allouees pour les lois d'occupation des etats).
 *
 *--------------------------------------------------------------*/

void SemiMarkovChain::copy(const SemiMarkovChain &smarkov , int param)

{
  register int i;


  state_subtype = new int[nb_state];
  for (i = 0;i < nb_state;i++) {
    state_subtype[i] = smarkov.state_subtype[i];
  }

  forward = new Forward*[nb_state];

  for (i = 0;i < nb_state;i++) {
    if (smarkov.forward[i]) {
      forward[i] = new Forward(*(smarkov.forward[i]) , param);
    }
    else {
      forward[i] = NULL;
    }
  }

  switch (param) {
  case I_DEFAULT :
    state_process = new CategoricalSequenceProcess(*(smarkov.state_process));
    break;
  case 0 :
    state_process = new CategoricalSequenceProcess(*(smarkov.state_process) , 'o' , I_DEFAULT);
    break;
  default :
    state_process = new CategoricalSequenceProcess(*(smarkov.state_process) , 'o' , param);
    break;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet SemiMarkovChain.
 *
 *--------------------------------------------------------------*/

void SemiMarkovChain::remove()

{
  register int i;


  delete [] state_subtype;

  delete state_process;

  if (forward) {
    for (i = 0;i < nb_state;i++) {
      delete forward[i];
    }
    delete [] forward;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe SemiMarkovChain.
 *
 *--------------------------------------------------------------*/

SemiMarkovChain::~SemiMarkovChain()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe SemiMarkovChain.
 *
 *  argument : reference sur un objet SemiMarkovChain.
 *
 *--------------------------------------------------------------*/

SemiMarkovChain& SemiMarkovChain::operator=(const SemiMarkovChain &smarkov)

{
  if (&smarkov != this) {
    remove();
    Chain::remove();

    Chain::copy(smarkov);
    copy(smarkov);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Calcul du nombre de parametres independants d'un objet SemiMarkovChain.
 *
 *  argument : probabilite minimum.
 *
 *--------------------------------------------------------------*/

int SemiMarkovChain::nb_parameter_computation(double min_probability) const

{
  register int i;
  int nb_parameter = Chain::nb_parameter_computation(min_probability);


  for (i = 0;i < nb_state;i++) {
    if (state_subtype[i] == SEMI_MARKOVIAN) {
      nb_parameter += state_process->sojourn_time[i]->nb_parameter_computation();
      if (state_process->sojourn_time[i]->inf_bound == 1) {
        nb_parameter--;
      }
    }
  }

  return nb_parameter;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe SemiMarkov.
 *
 *--------------------------------------------------------------*/

SemiMarkov::SemiMarkov()

{
  nb_iterator = 0;
  semi_markov_data = NULL;

  nb_output_process = 0;
  categorical_process = NULL;
  discrete_parametric_process = NULL;
  continuous_parametric_process = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe SemiMarkov.
 *
 *  arguments : type, nombre d'etats, nombre de processus d'observation,
 *              nombre de valeurs observees par processus.
 *
 *--------------------------------------------------------------*/

SemiMarkov::SemiMarkov(char itype , int inb_state , int inb_output_process , int *nb_value)
:SemiMarkovChain(itype , inb_state)

{
  register int i;


  nb_iterator = 0;
  semi_markov_data = NULL;

  nb_output_process = inb_output_process;

  categorical_process = new CategoricalSequenceProcess*[nb_output_process];
  discrete_parametric_process = new DiscreteParametricProcess*[nb_output_process];
  continuous_parametric_process = new ContinuousParametricProcess*[nb_output_process];

  for (i = 0;i < nb_output_process;i++) {
    if (nb_value[i] == I_DEFAULT) {
      categorical_process[i] = NULL;
      discrete_parametric_process[i] = NULL;
      continuous_parametric_process[i] = new ContinuousParametricProcess(nb_state);
    }

    else if (nb_value[i] <= NB_OUTPUT) {
      categorical_process[i] = new CategoricalSequenceProcess(nb_state , nb_value[i] , true);
      discrete_parametric_process[i] = NULL;
      continuous_parametric_process[i] = NULL;
    }

    else {
      categorical_process[i] = NULL;
      discrete_parametric_process[i] = new DiscreteParametricProcess(nb_state , (int)(nb_value[i] * SAMPLE_NB_VALUE_COEFF));
      continuous_parametric_process[i] = NULL;
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe SemiMarkov.
 *
 *  arguments : pointeur sur un objet Chain, sur un objet CategoricalSequenceProcess et
 *              sur un objet CategoricalProcess, longueur des sequences,
 *              flag sur le calcul des lois de comptage.
 *
 *--------------------------------------------------------------*/

SemiMarkov::SemiMarkov(const Chain *pchain , const CategoricalSequenceProcess *poccupancy ,
                       const CategoricalProcess *pobservation ,
                       int length , bool counting_flag)
:SemiMarkovChain(pchain , poccupancy)

{
  register int i;


  nb_iterator = 0;
  semi_markov_data = NULL;

  nb_output_process = (pobservation ? 1 : 0);

  if (nb_output_process == 1) {
    categorical_process = new CategoricalSequenceProcess*[nb_output_process];
    categorical_process[0] = new CategoricalSequenceProcess(*pobservation);
  }
  else {
    categorical_process = NULL;
  }

  discrete_parametric_process = NULL;
  continuous_parametric_process = NULL;

  if (length > COUNTING_MAX_LENGTH) {
    counting_flag = false;
  }
  characteristic_computation(length , counting_flag);
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet SemiMarkov.
 *
 *  arguments : reference sur un objet SemiMarkov,
 *              flag copie de l'objet SemiMarkovData, parametre.
 *
 *--------------------------------------------------------------*/

void SemiMarkov::copy(const SemiMarkov &smarkov , bool data_flag , int param)

{
  register int i;


  nb_iterator = 0;

  if ((data_flag) && (smarkov.semi_markov_data)) {
    semi_markov_data = new SemiMarkovData(*(smarkov.semi_markov_data) , false);
  }
  else {
    semi_markov_data = NULL;
  }

  nb_output_process = smarkov.nb_output_process;

  if (smarkov.categorical_process) {
    categorical_process = new CategoricalSequenceProcess*[nb_output_process];

    switch (param) {

    case I_DEFAULT : {
      for (i = 0;i < nb_output_process;i++) {
        if (smarkov.categorical_process[i]) {
          categorical_process[i] = new CategoricalSequenceProcess(*(smarkov.categorical_process[i]));
        }
        else {
          categorical_process[i] = NULL;
        }
      }
      break;
    }

    default : {
      for (i = 0;i < nb_output_process;i++) {
        if (smarkov.categorical_process[i]) {
          categorical_process[i] = new CategoricalSequenceProcess(*(smarkov.categorical_process[i]) ,
                                                                  'c' , false);
        }
        else {
          categorical_process[i] = NULL;
        }
      }
      break;
    }
    }
  }

  else {
    categorical_process = NULL;
  }

  if (smarkov.discrete_parametric_process) {
    discrete_parametric_process = new DiscreteParametricProcess*[nb_output_process];

    for (i = 0;i < nb_output_process;i++) {
      if (smarkov.discrete_parametric_process[i]) {
        discrete_parametric_process[i] = new DiscreteParametricProcess(*(smarkov.discrete_parametric_process[i]));
      }
      else {
        discrete_parametric_process[i] = NULL;
      }
    }
  }

  else {
    discrete_parametric_process = NULL;
  }

  if (smarkov.continuous_parametric_process) {
    continuous_parametric_process = new ContinuousParametricProcess*[nb_output_process];

    for (i = 0;i < nb_output_process;i++) {
      if (smarkov.continuous_parametric_process[i]) {
        continuous_parametric_process[i] = new ContinuousParametricProcess(*(smarkov.continuous_parametric_process[i]));
      }
      else {
        continuous_parametric_process[i] = NULL;
      }
    }
  }

  else {
    continuous_parametric_process = NULL;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet SemiMarkov.
 *
 *--------------------------------------------------------------*/

void SemiMarkov::remove()

{
  register int i;


  delete semi_markov_data;

  if (categorical_process) {
    for (i = 0;i < nb_output_process;i++) {
      delete categorical_process[i];
    }
    delete [] categorical_process;
  }

  if (discrete_parametric_process) {
    for (i = 0;i < nb_output_process;i++) {
      delete discrete_parametric_process[i];
    }
    delete [] discrete_parametric_process;
  }

  if (continuous_parametric_process) {
    for (i = 0;i < nb_output_process;i++) {
      delete continuous_parametric_process[i];
    }
    delete [] continuous_parametric_process;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe SemiMarkov.
 *
 *--------------------------------------------------------------*/

SemiMarkov::~SemiMarkov()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Destruction d'un objet SemiMarkov en tenant compte du nombre
 *  d'iterateurs pointant dessus.
 *
 *--------------------------------------------------------------*/

void SemiMarkov::conditional_delete()

{
  if (nb_iterator == 0) {
    delete this;
  }
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe SemiMarkov.
 *
 *  argument : reference sur un objet SemiMarkov.
 *
 *--------------------------------------------------------------*/

SemiMarkov& SemiMarkov::operator=(const SemiMarkov &smarkov)

{
  if ((&smarkov != this) && (nb_iterator == 0)) {
    remove();
    SemiMarkovChain::remove();
    Chain::remove();

    Chain::copy(smarkov);
    SemiMarkovChain::copy(smarkov);
    copy(smarkov);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Extraction d'une loi.
 *
 *  arguments : reference sur un objet StatError, type de loi,
 *              variable, etat ou observation.
 *
 *--------------------------------------------------------------*/

DiscreteParametricModel* SemiMarkov::extract(StatError &error , int type ,
                                             int variable , int value) const

{
  bool status = true;
  int hvariable;
  Distribution *pdist;
  DiscreteParametric *pparam;
  DiscreteParametricModel *dist;
  FrequencyDistribution *phisto;
  CategoricalSequenceProcess *process;


  dist = NULL;
  error.init();

  pdist = NULL;
  pparam = NULL;

  if (type == OBSERVATION) {
    if ((variable < 1) || (variable > nb_output_process)) {
      status = false;
      error.update(STAT_error[STATR_OUTPUT_PROCESS_INDEX]);
    }

    else {
      if ((value < 0) || (value >= nb_state)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_STATE] << " " << value << " "
                      << STAT_error[STATR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }

      else {
        if (categorical_process[variable - 1]) {
          pdist = categorical_process[variable - 1]->observation[value];
        }
        else if (discrete_parametric_process[variable - 1]) {
          pparam = discrete_parametric_process[variable - 1]->observation[value];
        }
        else {
          status = false;
          ostringstream correction_message;
          correction_message << STAT_label[STATL_CATEGORICAL] << " or "
                             << STAT_label[STATL_DISCRETE_PARAMETRIC];
          error.correction_update(STAT_error[STATR_OUTPUT_PROCESS_TYPE] , (correction_message.str()).c_str());
        }
      }
    }
  }

  else {
    if ((variable < 0) || (variable > nb_output_process)) {
      status = false;
      error.update(STAT_error[STATR_OUTPUT_PROCESS_INDEX]);
    }

    else {
      if (variable == 0) {
        process = state_process;
      }

      else {
        process = categorical_process[variable - 1];

        if (!process) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << variable << ": " 
                        << SEQ_error[SEQR_CHARACTERISTICS_NOT_COMPUTED];
          error.update((error_message.str()).c_str());
        }
      }

      if ((process) && ((value < 0) || (value >= process->nb_value))) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[variable == 0 ? STATL_STATE : STATL_OUTPUT] << " "
                      << value << " " << STAT_error[STATR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }

      if (status) {
        switch (type) {
        case FIRST_OCCURRENCE :
          pdist = process->first_occurrence[value];
          break;
        case RECURRENCE_TIME :
          pdist = process->recurrence_time[value];
          break;
        case SOJOURN_TIME :
          pparam = process->sojourn_time[value];
          break;
        case NB_RUN :
          pdist = process->nb_run[value];
          break;
        case NB_OCCURRENCE :
          pdist = process->nb_occurrence[value];
          break;
        }

        if ((!pdist) && (!pparam)) {
          status = false;
          error.update(SEQ_error[SEQR_NON_EXISTING_CHARACTERISTIC_DISTRIBUTION]);
        }
      }
    }
  }

  if (status) {
    phisto = NULL;

    if (semi_markov_data) {
      switch (semi_markov_data->type[0]) {
      case STATE :
        hvariable = variable;
        break;
      case INT_VALUE :
        hvariable = variable - 1;
        break;
      }

      if (hvariable >= 0) {
        switch (type) {

        case OBSERVATION : {
          if ((semi_markov_data->observation_distribution) &&
              (semi_markov_data->observation_distribution[hvariable])) {
            phisto = semi_markov_data->observation_distribution[hvariable][value];
          }
          break;
        }

        case FIRST_OCCURRENCE : {
          phisto = semi_markov_data->characteristics[hvariable]->first_occurrence[value];
          break;
        }

        case RECURRENCE_TIME : {
          if (semi_markov_data->characteristics[hvariable]->recurrence_time[value]->nb_element > 0) {
            phisto = semi_markov_data->characteristics[hvariable]->recurrence_time[value];
          }
          break;
        }

        case SOJOURN_TIME : {
          if (semi_markov_data->characteristics[hvariable]->sojourn_time[value]->nb_element > 0) {
            phisto = semi_markov_data->characteristics[hvariable]->sojourn_time[value];
          }
          break;
        }

        case NB_RUN : {
          phisto = semi_markov_data->characteristics[hvariable]->nb_run[value];
          break;
        }

        case NB_OCCURRENCE : {
          phisto = semi_markov_data->characteristics[hvariable]->nb_occurrence[value];
          break;
        }
        }
      }
    }

    if (pdist) {
      dist = new DiscreteParametricModel(*pdist , phisto);
    }
    else if (pparam) {
      dist = new DiscreteParametricModel(*pparam , phisto);
    }
  }

  return dist;
}


/*--------------------------------------------------------------*
 *
 *  Extraction d'une loi de l'intervalle de temps residuel.
 *
 *  arguments : reference sur un objet StatError, etat,
 *              type de loi empirique associe.
 *
 *--------------------------------------------------------------*/

DiscreteParametricModel* SemiMarkov::extract(StatError &error , int state ,
                                             int frequency_distribution_type) const

{
  bool status = true;
  Distribution *pdist;
  DiscreteParametricModel *dist;
  FrequencyDistribution *phisto;


  dist = NULL;
  error.init();

  pdist = NULL;

  if ((state < 0) || (state >= nb_state)) {
    status = false;
    ostringstream error_message;
    error_message << STAT_label[STATL_STATE] << " " << state << " "
                  << STAT_error[STATR_NOT_PRESENT];
    error.update((error_message.str()).c_str());
  }

  else {
    pdist = forward[state];

    if (!pdist) {
      status = false;
      error.update(SEQ_error[SEQR_NON_EXISTING_FORWARD_DISTRIBUTION]);
    }

    else {
      phisto = NULL;

      if ((semi_markov_data) && (semi_markov_data->type[0] == STATE)) {
        switch (frequency_distribution_type) {

        case INITIAL_RUN : {
          if ((semi_markov_data->characteristics[0]->initial_run) &&
              (semi_markov_data->characteristics[0]->initial_run[state]->nb_element > 0)) {
            phisto = semi_markov_data->characteristics[0]->initial_run[state];
          }
          break;
        }

        case FINAL_RUN : {
          if (semi_markov_data->characteristics[0]->final_run[state]->nb_element > 0) {
            phisto = semi_markov_data->characteristics[0]->final_run[state];
          }
          break;
        }
        }
      }

      dist = new DiscreteParametricModel(*pdist , phisto);
    }
  }

  return dist;
}


/*--------------------------------------------------------------*
 *
 *  Extraction de la partie "donnees" d'un objet SemiMarkov.
 *
 *  argument : reference sur un objet StatError.
 *
 *--------------------------------------------------------------*/

SemiMarkovData* SemiMarkov::extract_data(StatError &error) const

{
  bool status = true;
  SemiMarkovData *seq;


  seq = NULL;
  error.init();

  if (!semi_markov_data) {
    status = false;
    error.update(STAT_error[STATR_NO_DATA]);
  }
  else if (nb_output_process + 1 != semi_markov_data->nb_variable) {
    status = false;
    error.update(SEQ_error[SEQR_STATE_SEQUENCES]);
  }

  if (status) {
    seq = new SemiMarkovData(*semi_markov_data);
    seq->semi_markov = new SemiMarkov(*this , false);
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Application d'un seuil sur les parametres d'une semi-chaine de Markov.
 *
 *  argument : probabilite minimum.
 *
 *--------------------------------------------------------------*/

SemiMarkov* SemiMarkov::thresholding(double min_probability) const

{
  register int i;
  SemiMarkov *smarkov;


  smarkov = new SemiMarkov(*this , false , 0);
  smarkov->Chain::thresholding(min_probability , true);

  for (i = 0;i < smarkov->nb_output_process;i++) {
    if (smarkov->categorical_process[i]) {
      smarkov->categorical_process[i]->thresholding(min_probability);
    }
  }

  return smarkov;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet SemiMarkov a partir d'un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              longueur des sequences, flag sur le calcul des lois de comptage,
 *              seuil sur les fonctions de repartition des lois d'occupation des etats.
 *
 *--------------------------------------------------------------*/

SemiMarkov* semi_markov_ascii_read(StatError &error , const char *path , int length ,
                                   bool counting_flag , double cumul_threshold)

{
  RWCString buffer , token;
  size_t position;
  char type = 'v';
  bool status;
  register int i;
  int line;
  const Chain *chain;
  const CategoricalSequenceProcess *occupancy;
  const CategoricalProcess *observation;
  SemiMarkov *smarkov;
  ifstream in_file(path);


  smarkov = NULL;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    line = 0;

    if (length < 2) {
      status = false;
      error.update(SEQ_error[SEQR_SHORT_SEQUENCE_LENGTH]);
    }
    if (length > MAX_LENGTH) {
      status = false;
      error.update(SEQ_error[SEQR_LONG_SEQUENCE_LENGTH]);
    }

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

        // test mot cle (EQUILIBRIUM) SEMI-MARKOV_CHAIN

        if (i == 0) {
          if (token == SEQ_word[SEQW_SEMI_MARKOV_CHAIN]) {
            type = 'o';
          }
          else if (token == SEQ_word[SEQW_EQUILIBRIUM_SEMI_MARKOV_CHAIN]) {
            type = 'e';
          }
          else {
            status = false;
            ostringstream correction_message;
            correction_message << SEQ_word[SEQW_SEMI_MARKOV_CHAIN] << " or "
                               << SEQ_word[SEQW_EQUILIBRIUM_SEMI_MARKOV_CHAIN];
            error.correction_update(STAT_parsing[STATP_KEY_WORD] ,
                                    (correction_message.str()).c_str() , line);
          }
        }

        i++;
      }

      if (i > 0) {
        if (i != 1) {
          status = false;
          error.update(STAT_parsing[STATP_FORMAT] , line);
        }
        break;
      }
    }

    if (type != 'v') {
      // analyse du format et lecture de la chaine de Markov

      chain = chain_parsing(error , in_file , line , type);

      if (chain) {

        // analyse du format et lecture des lois d'occupation de etats

        occupancy = occupancy_parsing(error , in_file , line , *chain , cumul_threshold);
        if (!occupancy) {
          status = false;
        }
        // analyse du format et lecture des lois d'observation

        observation = NULL;

        while (buffer.readLine(in_file , false)) {
          line++;

#         ifdef DEBUG
          cout << line << "  " << buffer << endl;
#         endif

          position = buffer.first('#');
          if (position != RW_NPOS) {
            buffer.remove(position);
          }
          i = 0;

          RWCTokenizer next(buffer);

          while (!((token = next()).isNull())) {

            // test mot cle OUTPUT_PROCESS

            if (i == 0) {
              if (token != STAT_word[STATW_OUTPUT_PROCESS]) {
                status = false;
                error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_OUTPUT_PROCESS] , line);
              }
            }
            i++;

          }

          if (i > 0) {
            if (i != 1) {
              status = false;
              error.update(STAT_parsing[STATP_FORMAT] , line);
            }
            // analyse du format et lecture des lois d'observation

            observation = categorical_observation_parsing(error , in_file , line , chain->nb_state ,
                                                          HIDDEN_MARKOV , false);
            if (!observation) {
              status = false;
            }
            break;
          }
        }

        while (buffer.readLine(in_file , false)) {
          line++;
#         ifdef DEBUG
          cout << line << "  " << buffer << endl;
#         endif

          position = buffer.first('#');
          if (position != RW_NPOS) {
            buffer.remove(position);
          }
          if (!(buffer.isNull())) {
            status = false;
            error.update(STAT_parsing[STATP_FORMAT] , line);
          }
        }

        if (status) {
          smarkov = new SemiMarkov(chain , occupancy , observation , length , counting_flag);
        }

        delete chain;
        delete occupancy;
        delete observation;
      }
    }
  }

  return smarkov;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet SemiMarkov.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& SemiMarkov::line_write(ostream &os) const

{
  os << nb_state << " " << STAT_word[STATW_STATES];

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet SemiMarkov et de la structure de donnees associee.
 *
 *  arguments : stream, pointeur sur un objet SemiMarkovData,
 *              flag niveau de detail, flag fichier, flag semi-Markov cache.
 *
 *--------------------------------------------------------------*/

ostream& SemiMarkov::ascii_write(ostream &os , const SemiMarkovData *seq ,
                                 bool exhaustive , bool file_flag , bool hidden) const

{
  register int i , j , k;
  int buff , width , variable;
  double **distance;
  FrequencyDistribution *marginal_dist = NULL , **observation_dist = NULL;
  Histogram *marginal_histo = NULL , **observation_histo = NULL;
  SequenceCharacteristics *characteristics = NULL;
  long old_adjust;


  old_adjust = os.setf(ios::left , ios::adjustfield);

  switch (hidden) {

  case false : {
    switch (type) {
    case 'o' :
      os << SEQ_word[SEQW_SEMI_MARKOV_CHAIN] << endl;
      break;
    case 'e' :
      os << SEQ_word[SEQW_EQUILIBRIUM_SEMI_MARKOV_CHAIN] << endl;
      break;
    }
    break;
  }

  case true : {
    switch (type) {
    case 'o' :
      os << SEQ_word[SEQW_HIDDEN_SEMI_MARKOV_CHAIN] << endl;
      break;
    case 'e' :
      os << SEQ_word[SEQW_EQUILIBRIUM_HIDDEN_SEMI_MARKOV_CHAIN] << endl;
      break;
    }
    break;
  }
  }

  // ecriture des parametres de la chaine de Markov

  ascii_print(os , file_flag);

  // ecriture des lois d'occupation des etats

  if ((seq) && (seq->type[0] == STATE)) {
    characteristics = seq->characteristics[0];
  }
  else {
    characteristics = NULL;
  }

  state_process->ascii_print(os , 0 , NULL , NULL , characteristics ,
                             exhaustive , file_flag , forward);

  if (hidden) {
    for (i = 0;i < nb_output_process;i++) {
      if (discrete_parametric_process[i]) {
        if (discrete_parametric_process[i]->weight) {
          width = column_width(nb_state , discrete_parametric_process[i]->weight->mass);
        }
        else {
          width = 0;
        }
        if (discrete_parametric_process[i]->restoration_weight) {
          buff = column_width(nb_state , discrete_parametric_process[i]->restoration_weight->mass);
          if (buff > width) {
            width = buff;
          }
        }
        width++;

        if (discrete_parametric_process[i]->weight) {
          os << "\n";
          if (file_flag) {
            os << "# ";
          }
          os << STAT_label[STATL_THEORETICAL] << " " << SEQ_label[SEQL_STATE_PROBABILITY] << ": ";

          for (j = 0;j < nb_state;j++) {
            os << setw(width) << discrete_parametric_process[i]->weight->mass[j];
          }
          os << endl;
        }

        if (discrete_parametric_process[i]->restoration_weight) {
          if (file_flag) {
            os << "# ";
          }
          os << STAT_label[STATL_RESTORATION] << " " << SEQ_label[SEQL_STATE_PROBABILITY] << ": ";

          for (j = 0;j < nb_state;j++) {
            os << setw(width) << discrete_parametric_process[i]->restoration_weight->mass[j];
          }
          os << endl;
        }

        break;
      }

      else if (continuous_parametric_process[i]) {
        if (continuous_parametric_process[i]->weight) {
          width = column_width(nb_state , continuous_parametric_process[i]->weight->mass);
        }
        else {
          width = 0;
        }
        if (continuous_parametric_process[i]->restoration_weight) {
          buff = column_width(nb_state , continuous_parametric_process[i]->restoration_weight->mass);
          if (buff > width) {
            width = buff;
          }
        }
        width++;

        if (continuous_parametric_process[i]->weight) {
          os << "\n";
          if (file_flag) {
            os << "# ";
          }
          os << STAT_label[STATL_THEORETICAL] << " " << SEQ_label[SEQL_STATE_PROBABILITY] << ": ";

          for (j = 0;j < nb_state;j++) {
            os << setw(width) << continuous_parametric_process[i]->weight->mass[j];
          }
          os << endl;
        }

        if (continuous_parametric_process[i]->restoration_weight) {
          if (file_flag) {
            os << "# ";
          }
          os << STAT_label[STATL_RESTORATION] << " " << SEQ_label[SEQL_STATE_PROBABILITY] << ": ";

          for (j = 0;j < nb_state;j++) {
            os << setw(width) << continuous_parametric_process[i]->restoration_weight->mass[j];
          }
          os << endl;
        }

        break;
      }
    }

    os << "\n" << nb_output_process << " "
       << STAT_word[nb_output_process == 1 ? STATW_OUTPUT_PROCESS : STATW_OUTPUT_PROCESSES] << endl;
  }

  // ecriture des lois associees a chaque processus d'observation

  if (hidden) {
    distance = new double*[nb_state];
    for (i = 0;i < nb_state;i++) {
      distance[i] = new double[nb_state];
    }
  }

  for (i = 0;i < nb_output_process;i++) {
    os << "\n" << STAT_word[STATW_OUTPUT_PROCESS];

    if (hidden) {
      os << " " << i + 1;

      if (categorical_process[i]) {
        os << " : " << STAT_word[STATW_CATEGORICAL];
      }
      else if (discrete_parametric_process[i]) {
        os << " : " << STAT_word[STATW_DISCRETE_PARAMETRIC];
      }
      else {
        os << " : " << STAT_word[STATW_CONTINUOUS_PARAMETRIC];
      }
    }
    os << endl;
    if ((continuous_parametric_process[i]) && (continuous_parametric_process[i]->ident == LINEAR_MODEL)) {
      for (j = 0;j < nb_state;j++) {
        os << "\n" << STAT_word[STATW_STATE] << " " << j << " "
           << STAT_word[STATW_OBSERVATION_MODEL] << endl;
        continuous_parametric_process[i]->observation[j]->ascii_parameter_print(os , file_flag);
      }
    }

    else {
      if (seq) {
        switch (seq->type[0]) {
        case STATE :
          variable = i + 1;
          break;
        default :
          variable = i;
          break;
        }

        if (seq->observation_distribution) {
          observation_dist = seq->observation_distribution[variable];
        }
        marginal_dist = seq->marginal_distribution[variable];

        if (seq->observation_histogram) {
          observation_histo = seq->observation_histogram[variable];
        }
        marginal_histo = seq->marginal_histogram[variable];

        characteristics = seq->characteristics[variable];
      }

      if (categorical_process[i]) {
        categorical_process[i]->ascii_print(os , i + 1 , observation_dist , marginal_dist ,
                                            characteristics , exhaustive , file_flag);

        if (hidden) {
          for (j = 0;j < nb_state;j++) {
            distance[j][j] = 0.;

            for (k = j + 1;k < nb_state;k++) {
              if ((transition[j][k] > MIN_PROBABILITY) || (transition[k][j] > MIN_PROBABILITY)) {
                distance[j][k] = categorical_process[i]->observation[j]->overlap_distance_computation(*(categorical_process[i]->observation[k]));
              }
              else {
                distance[j][k] = 1.;
              }

              distance[k][j] = distance[j][k];
            }
          }
        }
      }

      else if (discrete_parametric_process[i]) {
        discrete_parametric_process[i]->ascii_print(os , observation_dist , marginal_dist ,
                                                    exhaustive , file_flag);

        if (hidden) {
          for (j = 0;j < nb_state;j++) {
            distance[j][j] = 0.;

            for (k = j + 1;k < nb_state;k++) {
              if ((transition[j][k] > MIN_PROBABILITY) || (transition[k][j] > MIN_PROBABILITY)) {
                distance[j][k] = discrete_parametric_process[i]->observation[j]->sup_norm_distance_computation(*(discrete_parametric_process[i]->observation[k]));
              }
              else {
                distance[j][k] = 1.;
              }

              distance[k][j] = distance[j][k];
            }
          }
        }
      }

      else {
        continuous_parametric_process[i]->ascii_print(os , observation_histo , observation_dist ,
                                                      marginal_histo , marginal_dist ,
                                                      exhaustive , file_flag);

        if (hidden) {
          for (j = 0;j < nb_state;j++) {
            distance[j][j] = 0.;

            for (k = j + 1;k < nb_state;k++) {
              if ((transition[j][k] > MIN_PROBABILITY) || (transition[k][j] > MIN_PROBABILITY)) {
                distance[j][k] = continuous_parametric_process[i]->observation[j]->sup_norm_distance_computation(*(continuous_parametric_process[i]->observation[k]));
              }
              else {
                distance[j][k] = 1.;
              }

              distance[k][j] = distance[j][k];
            }
          }
        }
      }

      if (hidden) {
        width = column_width(nb_state , distance[0]);
        for (j = 1;j < nb_state;j++) {
          buff = column_width(nb_state , distance[j]);
          if (buff > width) {
            width = buff;
          }
        }
        width += ASCII_SPACE;

        os.setf(ios::left , ios::adjustfield);

        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_CONSECUTIVE_STATE_OBSERVATION_DISTRIBUTION_DISTANCE] << endl;

        for (j = 0;j < nb_state;j++) {
          if (file_flag) {
            os << "# ";
          }
          for (k = 0;k < nb_state;k++) {
            if ((k != j) && (transition[j][k] > MIN_PROBABILITY)) {
              os << setw(width) << distance[j][k];
            }
            else {
              os << setw(width) << "_";
            }
          }
          os << endl;
        }
      }
    }
  }

  if (hidden) {
    for (i = 0;i < nb_state;i++) {
      delete [] distance[i];
    }
    delete [] distance;
  }

  if (seq) {
    int nb_parameter = nb_parameter_computation(hidden ? MIN_PROBABILITY : 0.);
    double information;


    // ecriture de la loi empirique des longueurs des sequences

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    seq->length_distribution->ascii_characteristic_print(os , false , file_flag);

    if (exhaustive) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << "   | " << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
      seq->length_distribution->ascii_print(os , file_flag);
    }

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_CUMUL_LENGTH] << ": " << seq->cumul_length << endl;

    // ecriture de la quantite d'information des sequences dans le cas iid

    for (i = 0;i < seq->nb_variable;i++) {
      if (seq->type[i] == REAL_VALUE) {
        break;
      }
    }

    if (i == seq->nb_variable) {
      information = seq->iid_information_computation();

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_IID_INFORMATION] << ": " << information << " ("
         << information / seq->cumul_length << ")" << endl;
    }

    // ecriture des vraisemblances des sequences

    switch (hidden) {

    case false : {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_LIKELIHOOD] << ": " << seq->likelihood << "   ("
         << STAT_label[STATL_NORMALIZED] << ": " << seq->likelihood / seq->cumul_length << ")" << endl;
      break;
    }

    case true : {
      if (seq->restoration_likelihood != D_INF) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << ": " << seq->restoration_likelihood << "   ("
           << STAT_label[STATL_NORMALIZED] << ": " << seq->restoration_likelihood / seq->cumul_length << ")" << endl;
      }

      if (seq->sample_entropy != D_DEFAULT) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_STATE_SEQUENCE_ENTROPY] << ": " << seq->sample_entropy << "   ("
           << STAT_label[STATL_NORMALIZED] << ": " << seq->sample_entropy / seq->cumul_length << ")" << endl;
      }

      if (seq->likelihood != D_INF) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_OBSERVED_SEQUENCES_LIKELIHOOD] << ": " << seq->likelihood << "   ("
           << STAT_label[STATL_NORMALIZED] << ": " << seq->likelihood / seq->cumul_length << ")" << endl;
      }
      break;
    }
    }

    if (seq->likelihood != D_INF) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AIC] << "): "
         << 2 * (seq->likelihood - nb_parameter) << endl;

      if (nb_parameter < seq->cumul_length - 1) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
           << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AICc] << "): "
           << 2 * (seq->likelihood - (double)(nb_parameter * seq->cumul_length) /
             (double)(seq->cumul_length - nb_parameter - 1)) << endl;
      }

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BIC] << "): "
         << 2 * seq->likelihood - nb_parameter * log((double)seq->cumul_length) << endl;

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BICc] << "): "
         << 2 * seq->likelihood - penalty_computation(hidden , (hidden ? MIN_PROBABILITY : 0.)) << endl;
    }

//    if ((hidden) && (seq->restoration_likelihood != D_INF)) {
    if ((hidden) && (seq->likelihood != D_INF)) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[ICL] << "): "
//         << 2 * seq->restoration_likelihood - nb_parameter * log((double)seq->cumul_length) << endl;
         << 2 * (seq->likelihood - seq->sample_entropy) - nb_parameter * log((double)seq->cumul_length) << endl;

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[ICLc] << "): "
//         << 2 * seq->restoration_likelihood - penalty_computation(hidden , MIN_PROBABILITY) << endl;
         << 2 * (seq->likelihood - seq->sample_entropy) - penalty_computation(hidden , MIN_PROBABILITY) << endl;
    }
  }

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet SemiMarkov.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& SemiMarkov::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , semi_markov_data , exhaustive , false);
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet SemiMarkov dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool SemiMarkov::ascii_write(StatError &error , const char *path ,
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
    ascii_write(out_file , semi_markov_data , exhaustive , true);
  }

  return status;
}

/*--------------------------------------------------------------*
 *  B.O.
 *  Write an initialization file for a Hidden Semi Markov Model
 *
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool SemiMarkov::write_hidden_semi_markov_init_file(StatError &error , const char *path) const {
    bool status;
    register int i , j , k;
    int buff , width , variable;
    double sum_probabilities, current_probability;
    FrequencyDistribution *marginal_dist = NULL , **observation_dist = NULL;
    SequenceCharacteristics *characteristics = NULL;
    //SemiMarkovData* seq = semi_markov_data;

    ofstream out_file(path);
    error.init();
    if (!out_file) {
        status = false;
        error.update(STAT_error[STATR_FILE_NAME]);
    }
    else {
        status = true;
        // FILE HEADER
        out_file << SEQ_word[SEQW_HIDDEN_SEMI_MARKOV_CHAIN] << endl;
        out_file << "\n" << nb_state << " " << STAT_word[STATW_STATES] << endl;
        out_file << "\n" << STAT_word[STATW_INITIAL_PROBABILITIES] << endl;
        width = column_width(nb_state , initial);
        for (i = 0;i < nb_row;i++) {
            buff = column_width(nb_state , transition[i]);
            if (buff > width) {
                width = buff;
            }
        }
        width += ASCII_SPACE;

        // INITIAL_PROBABILITIES
        sum_probabilities = 0.0;
        current_probability = 0.0;
        for (i = 0;i < nb_state - 1;i++) {
            current_probability = (int)initial[i] + floor((initial[i] - (int)initial[i])*10000)/10000;
            sum_probabilities += current_probability;
            out_file << current_probability << setw(width);
        }
        out_file << setw(width) << 1 - sum_probabilities;
        out_file << endl;

        // TRANSITION_PROBABILITIES
        out_file << "\n" << STAT_word[STATW_TRANSITION_PROBABILITIES] << endl;

        for (i = 0;i < nb_row-1;i++) {
            sum_probabilities = 0.0;
            current_probability = 0.0;
            for (j = 0;j < nb_state - 1;j++) {
                current_probability = (int)transition[i][j] + floor((transition[i][j] - (int)transition[i][j])*10000)/10000;
                sum_probabilities += current_probability;
                out_file << current_probability << setw(width);
            }
            out_file << setw(width) << 1 - sum_probabilities;
            out_file << endl;
        }
        // sum probs last line debug
        sum_probabilities = 0.0;
        current_probability = 0.0;
        for (j = 0;j < nb_state - 1;j++) {
            current_probability = (int)transition[i][j] + floor((transition[i][j] - (int)transition[i][j])*10000)/10000;
            sum_probabilities += current_probability;
        }
        current_probability = (int)transition[i][0] + floor((transition[i][0] - (int)transition[i][0])*10000)/10000;
        current_probability += 1 - sum_probabilities;
        out_file << current_probability << setw(width);
        for (j = 1;j < nb_state;j++) {
            current_probability = (int)transition[i][j] + floor((transition[i][j] - (int)transition[i][j])*10000)/10000;
            out_file << current_probability << setw(width);
        }
        out_file << endl;


        // OCCUPANCY_DISTRIBUTION
        if ((semi_markov_data) && (semi_markov_data->type[0] == STATE)) {
            characteristics = semi_markov_data->characteristics[0];
        }
        else {
            characteristics = NULL;
        }
        state_process->ascii_print(out_file, 0, NULL, NULL, characteristics, false, true, forward);

        // OUTPUT_PROCESS
        out_file << "\n" << nb_output_process << " " << STAT_word[nb_output_process == 1 ? STATW_OUTPUT_PROCESS : STATW_OUTPUT_PROCESSES] << endl;
        for (i = 0;i < nb_output_process;i++) {
            out_file << "\n" << STAT_word[STATW_OUTPUT_PROCESS] << " " << i+1;
            out_file << " : " << STAT_word[STATW_CATEGORICAL] << endl;

            if (semi_markov_data->observation_distribution) {
              observation_dist = semi_markov_data->observation_distribution[variable];
            }
            marginal_dist = semi_markov_data->marginal_distribution[variable];
            characteristics = semi_markov_data->characteristics[variable];
            categorical_process[i]->ascii_print(out_file, i+1, observation_dist, marginal_dist, characteristics, false, true);
        }
    }
    return status;
}

/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet SemiMarkov et de la structure de donnees associee
 *  dans un fichier au format tableur.
 *
 *  arguments : stream, pointeur sur les sequences observees,
 *              flag semi-Markov cache.
 *
 *--------------------------------------------------------------*/

ostream& SemiMarkov::spreadsheet_write(ostream &os , const SemiMarkovData *seq ,
                                       bool hidden) const

{
  register int i , j , k;
  int variable;
  double **distance;
  FrequencyDistribution *marginal_dist = NULL , **observation_dist = NULL;
  Histogram *marginal_histo = NULL , **observation_histo = NULL;
  SequenceCharacteristics *characteristics = NULL;


  switch (hidden) {

  case false : {
    switch (type) {
    case 'o' :
      os << SEQ_word[SEQW_SEMI_MARKOV_CHAIN] << endl;
      break;
    case 'e' :
      os << SEQ_word[SEQW_EQUILIBRIUM_SEMI_MARKOV_CHAIN] << endl;
      break;
    }
    break;
  }

  case true : {
    switch (type) {
    case 'o' :
      os << SEQ_word[SEQW_HIDDEN_SEMI_MARKOV_CHAIN] << endl;
      break;
    case 'e' :
      os << SEQ_word[SEQW_EQUILIBRIUM_HIDDEN_SEMI_MARKOV_CHAIN] << endl;
      break;
    }
    break;
  }
  }

  // ecriture des parametres de la chaine de Markov

  spreadsheet_print(os);

  // ecriture des lois d'occupation des etats

  if ((seq) && (seq->type[0] == STATE)) {
    characteristics = seq->characteristics[0];
  }
  else {
    characteristics = NULL;
  }

  state_process->spreadsheet_print(os , 0 , NULL , NULL , characteristics , forward);

  // ecriture des lois associees a chaque processus d'observation

  if (hidden) {
    os << "\n" << nb_output_process << "\t"
       << STAT_word[nb_output_process == 1 ? STATW_OUTPUT_PROCESS : STATW_OUTPUT_PROCESSES] << endl;
  }

  if (hidden) {
    distance = new double*[nb_state];
    for (i = 0;i < nb_state;i++) {
      distance[i] = new double[nb_state];
    }
  }

  for (i = 0;i < nb_output_process;i++) {
    os << "\n" << STAT_word[STATW_OUTPUT_PROCESS];

    if (hidden) {
      os << "\t" << i + 1;

      if (categorical_process[i]) {
        os << "\t" << STAT_word[STATW_CATEGORICAL];
      }
      else if (discrete_parametric_process[i]) {
        os << "\t" << STAT_word[STATW_DISCRETE_PARAMETRIC];
      }
      else {
        os << "\t" << STAT_word[STATW_CONTINUOUS_PARAMETRIC];
      }
    }
    os << endl;

    if ((continuous_parametric_process[i]) && (continuous_parametric_process[i]->ident == LINEAR_MODEL)) {
      for (j = 0;j < nb_state;j++) {
        os << "\n" << STAT_word[STATW_STATE] << " " << j << "\t"
           << STAT_word[STATW_OBSERVATION_MODEL] << endl;
        continuous_parametric_process[i]->observation[j]->spreadsheet_parameter_print(os);
      }
    }

    if (seq) {
      switch (seq->type[0]) {
      case STATE :
        variable = i + 1;
        break;
      default :
        variable = i;
        break;
      }

      if ((continuous_parametric_process[i]) && (continuous_parametric_process[i]->ident == LINEAR_MODEL)) {
        seq->linear_model_spreadsheet_print(os , variable , continuous_parametric_process[i]);
      }

      else {
        if (seq->observation_distribution) {
          observation_dist = seq->observation_distribution[variable];
        }
        marginal_dist = seq->marginal_distribution[variable];

        if (seq->observation_histogram) {
          observation_histo = seq->observation_histogram[variable];
        }
        marginal_histo = seq->marginal_histogram[variable];

        characteristics = seq->characteristics[variable];

        if (categorical_process[i]) {
          categorical_process[i]->spreadsheet_print(os , i + 1 , observation_dist , marginal_dist ,
                                                    characteristics);

          if (hidden) {
            for (j = 0;j < nb_state;j++) {
              for (k = j + 1;k < nb_state;k++) {
                if ((transition[j][k] > MIN_PROBABILITY) || (transition[k][j] > MIN_PROBABILITY)) {
                  distance[j][k] = categorical_process[i]->observation[j]->overlap_distance_computation(*(categorical_process[i]->observation[k]));
                  distance[k][j] = distance[j][k];
                }
              }
            }
          }
        }

        else if (discrete_parametric_process[i]) {
          discrete_parametric_process[i]->spreadsheet_print(os , observation_dist , marginal_dist);

          if (hidden) {
            for (j = 0;j < nb_state;j++) {
              for (k = j + 1;k < nb_state;k++) {
                if ((transition[j][k] > MIN_PROBABILITY) || (transition[k][j] > MIN_PROBABILITY)) {
                  distance[j][k] = discrete_parametric_process[i]->observation[j]->sup_norm_distance_computation(*(discrete_parametric_process[i]->observation[k]));
                  distance[k][j] = distance[j][k];
                }
              }
            }
          }
        }

        else {
          continuous_parametric_process[i]->spreadsheet_print(os , observation_histo , observation_dist ,
                                                              marginal_histo , marginal_dist);

          if (hidden) {
            for (j = 0;j < nb_state;j++) {
              for (k = j + 1;k < nb_state;k++) {
                if ((transition[j][k] > MIN_PROBABILITY) || (transition[k][j] > MIN_PROBABILITY)) {
                  distance[j][k] = continuous_parametric_process[i]->observation[j]->sup_norm_distance_computation(*(continuous_parametric_process[i]->observation[k]));
                  distance[k][j] = distance[j][k];
                }
              }
            }
          }
        }

        if (hidden) {
          os << "\n" << STAT_label[STATL_CONSECUTIVE_STATE_OBSERVATION_DISTRIBUTION_DISTANCE] << endl;

          for (j = 0;j < nb_state;j++) {
            for (k = 0;k < nb_state;k++) {
              if ((k != j) && (transition[j][k] > MIN_PROBABILITY)) {
                os << distance[j][k];
              }
              os << "\t";
            }
            os << endl;
          }
        }
      }
    }
  }

  if (hidden) {
    for (i = 0;i < nb_state;i++) {
      delete [] distance[i];
    }
    delete [] distance;
  }

  if (seq) {
    int nb_parameter = nb_parameter_computation(hidden ? MIN_PROBABILITY : 0.);
    double information;


    // ecriture de la loi empirique des longueurs des sequences

    os << "\n" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
    seq->length_distribution->spreadsheet_characteristic_print(os);

    os << "\n\t" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
    seq->length_distribution->spreadsheet_print(os);

    os << "\n" << SEQ_label[SEQL_CUMUL_LENGTH] << "\t" << seq->cumul_length << endl;

    // ecriture de la quantite d'information des sequences dans le cas iid

    for (i = 0;i < seq->nb_variable;i++) {
      if (seq->type[i] == REAL_VALUE) {
        break;
      }
    }

    if (i == seq->nb_variable) {
      information = seq->iid_information_computation();

      os << "\n" << SEQ_label[SEQL_IID_INFORMATION] << "\t" << information << "\t"
         << information / seq->cumul_length << endl;
    }

    // ecriture des vraisemblances des sequences

    switch (hidden) {

    case false : {
      os << "\n" << STAT_label[STATL_LIKELIHOOD] << "\t" << seq->likelihood << "\t"
         << STAT_label[STATL_NORMALIZED] << "\t" << seq->likelihood / seq->cumul_length << endl;
      break;
    }

    case true : {
      if (seq->restoration_likelihood != D_INF) {
        os << "\n" << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << "\t" << seq->restoration_likelihood << "\t"
           << STAT_label[STATL_NORMALIZED] << "\t" << seq->restoration_likelihood / seq->cumul_length << endl;
      }

      if (seq->sample_entropy != D_DEFAULT) {
        os << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_ENTROPY] << "\t" << seq->sample_entropy << "\t"
           << STAT_label[STATL_NORMALIZED] << "\t" << seq->sample_entropy / seq->cumul_length << endl;
      }

      if (seq->likelihood != D_INF) {
        os << "\n" << SEQ_label[SEQL_OBSERVED_SEQUENCES_LIKELIHOOD] << "\t" << seq->likelihood << "\t"
           << STAT_label[STATL_NORMALIZED] << "\t" << seq->likelihood / seq->cumul_length << endl;
      }
      break;
    }
    }

    if (seq->likelihood != D_INF) {
      os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AIC] << ")\t"
         << 2 * (seq->likelihood - nb_parameter) << endl;

      if (nb_parameter < seq->cumul_length - 1) {
        os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
           << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AICc] << ")\t"
           << 2 * (seq->likelihood - (double)(nb_parameter * seq->cumul_length) /
              (double)(seq->cumul_length - nb_parameter - 1)) << endl;
      }

      os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BIC] << ")\t"
         << 2 * seq->likelihood - nb_parameter * log((double)seq->cumul_length) << endl;

      os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BICc] << ")\t"
         << 2 * seq->likelihood - penalty_computation(hidden , (hidden ? MIN_PROBABILITY : 0.)) << endl;
    }

//    if ((hidden) && (seq->restoration_likelihood != D_INF)) {
    if ((hidden) && (seq->likelihood != D_INF)) {
      os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[ICL] << ")\t"
//         << 2 * seq->restoration_likelihood - nb_parameter * log((double)seq->cumul_length) << endl;
         << 2 * (seq->likelihood - seq->sample_entropy) - nb_parameter * log((double)seq->cumul_length) << endl;

      os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[ICLc] << ")\t"
//         << 2 * seq->restoration_likelihood - penalty_computation(hidden , MIN_PROBABILITY) << endl;
         << 2 * (seq->likelihood - seq->sample_entropy) - penalty_computation(hidden , MIN_PROBABILITY) << endl;
    }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet SemiMarkov dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool SemiMarkov::spreadsheet_write(StatError &error , const char *path) const

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
    spreadsheet_write(out_file , semi_markov_data);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet SemiMarkov et de la structure
 *  de donnees associee.
 *
 *  arguments : prefixe des fichiers, titre des figures,
 *              pointeur sur les sequences observees.
 *
 *--------------------------------------------------------------*/

bool SemiMarkov::plot_write(const char *prefix , const char *title ,
                            const SemiMarkovData *seq) const

{
  bool status;
  register int i;
  int variable , nb_value = I_DEFAULT;
  double *empirical_cdf[2];
  FrequencyDistribution *length_distribution = NULL , *marginal_dist = NULL , **observation_dist = NULL;
  Histogram *marginal_histo = NULL , **observation_histo = NULL;
  SequenceCharacteristics *characteristics = NULL;


  if ((seq) && (seq->type[0] == STATE)) {
    characteristics = seq->characteristics[0];
    length_distribution = seq->length_distribution;
  }
  else {
    characteristics = NULL;
  }

  status = state_process->plot_print(prefix , title , 0 , NULL , NULL ,
                                     characteristics , length_distribution , forward);

  if (status) {
    if (seq) {
      length_distribution = seq->length_distribution;
    }

    for (i = 0;i < nb_output_process;i++) {
      if (seq) {
        switch (seq->type[0]) {
        case STATE :
          variable = i + 1;
          break;
        default :
          variable = i;
          break;
        }

        if ((continuous_parametric_process[i]) && (continuous_parametric_process[i]->ident == LINEAR_MODEL)) {
          seq->linear_model_plot_print(prefix , title , variable , continuous_parametric_process[i]);
        }

        else {
          if (seq->observation_distribution) {
            observation_dist = seq->observation_distribution[variable];
          }
          marginal_dist = seq->marginal_distribution[variable];

          if (seq->observation_histogram) {
            observation_histo = seq->observation_histogram[variable];
          }
          marginal_histo = seq->marginal_histogram[variable];

          characteristics = seq->characteristics[variable];

          if (continuous_parametric_process[i]) {
            nb_value = seq->cumulative_distribution_function_computation(variable , empirical_cdf);
          }
        }
      }

      if (categorical_process[i]) {
        categorical_process[i]->plot_print(prefix , title , i + 1 , observation_dist ,
                                           marginal_dist , characteristics ,
                                           length_distribution);
      }
      else if (discrete_parametric_process[i]) {
        discrete_parametric_process[i]->plot_print(prefix , title , i + 1 , observation_dist ,
                                                   marginal_dist);
      }
      else if (continuous_parametric_process[i]->ident != LINEAR_MODEL) {
        continuous_parametric_process[i]->plot_print(prefix , title , i + 1 ,
                                                     observation_histo , observation_dist ,
                                                     marginal_histo , marginal_dist ,
                                                     nb_value , (seq ? empirical_cdf : NULL));
        if (seq) {
          delete [] empirical_cdf[0];
          delete [] empirical_cdf[1];
        }
      }
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet SemiMarkov.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool SemiMarkov::plot_write(StatError &error , const char *prefix ,
                            const char *title) const

{
  bool status = plot_write(prefix , title , semi_markov_data);

  error.init();

  if (!status) {
    error.update(STAT_error[STATR_FILE_PREFIX]);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet SemiMarkov et de la structure
 *  de donnees associee.
 *
 *  argument : pointeur sur les sequences observees.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* SemiMarkov::get_plotable(const SemiMarkovData *seq) const

{
  register int i , j;
  int nb_plot_set , index_length , index , variable;
  FrequencyDistribution *length_distribution = NULL , *marginal_dist = NULL , **observation_dist = NULL;
  Histogram *marginal_histo = NULL , **observation_histo = NULL;
  SequenceCharacteristics *characteristics = NULL;
  MultiPlotSet *plot_set;


  if ((seq) && (seq->type[0] == STATE)) {
    characteristics = seq->characteristics[0];
  }
  else {
    characteristics = NULL;
  }

  // calcul du nombre de vues

  nb_plot_set = 0;

  if ((state_process->index_value) || (characteristics)) {
    nb_plot_set++;

    if (characteristics) {
      index_length = characteristics->index_value->plot_length_computation();

      if (characteristics->index_value->frequency[index_length - 1] < MAX_FREQUENCY) {
        nb_plot_set++;
      }
      nb_plot_set++;
    }
  }

  if ((state_process->first_occurrence) || (characteristics)) {
    for (i = 0;i < nb_state;i++) {
      if ((state_process->first_occurrence) &&
          (state_process->first_occurrence[i])) {
        nb_plot_set++;
      }
      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->first_occurrence[i]->nb_element > 0)) {
        nb_plot_set++;
      }
    }
  }

  if ((state_process->recurrence_time) || (characteristics)) {
    for (i = 0;i < nb_state;i++) {
      if ((state_process->recurrence_time) &&
          (state_process->recurrence_time[i])) {
        nb_plot_set++;
      }
      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->recurrence_time[i]->nb_element > 0)) {
        nb_plot_set++;
      }
    }
  }

  if ((state_process->sojourn_time) || (characteristics)) {
    for (i = 0;i < nb_state;i++) {
      if ((state_process->sojourn_time) &&
          (state_process->sojourn_time[i])) {
        nb_plot_set++;
      }
      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->sojourn_time[i]->nb_element > 0)) {
        nb_plot_set++;
      }

      if ((characteristics) && (i < characteristics->nb_value) &&
          (characteristics->initial_run) &&
          (characteristics->initial_run[i]->nb_element > 0)) {
        nb_plot_set++;
      }

      if ((forward) && (forward[i])) {
        nb_plot_set++;
      }
      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->final_run[i]->nb_element > 0)) {
        nb_plot_set++;
      }
    }
  }

  if ((state_process->nb_run) || (state_process->nb_occurrence) ||
      ((characteristics) && (characteristics->nb_run) && (characteristics->nb_occurrence))) {
    for (i = 0;i < nb_state;i++) {
      if (state_process->nb_run) {
        nb_plot_set++;
      }
      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->nb_run) && (characteristics->nb_run[i]->nb_element > 0)) {
        nb_plot_set++;
      }

      if (state_process->nb_occurrence) {
        nb_plot_set++;
      }
      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->nb_occurrence) &&
               (characteristics->nb_occurrence[i]->nb_element > 0)) {
        nb_plot_set++;
      }
    }

    if ((characteristics) && (characteristics->nb_run) && (characteristics->nb_occurrence)) {
      nb_plot_set++;
    }
  }

  for (i = 0;i < nb_output_process;i++) {
    if (seq) {
      switch (seq->type[0]) {
      case STATE :
        variable = i + 1;
        break;
      default :
        variable = i;
        break;
      }

      characteristics = seq->characteristics[variable];
    }

    if (categorical_process[i]) {
      if ((categorical_process[i]->index_value) || (characteristics)) {
        nb_plot_set++;

        if (characteristics) {
          index_length = characteristics->index_value->plot_length_computation();

          if (characteristics->index_value->frequency[index_length - 1] < MAX_FREQUENCY) {
            nb_plot_set++;
          }
          nb_plot_set++;
        }
      }

      if ((categorical_process[i]->first_occurrence) || (characteristics)) {
        for (j = 0;j < categorical_process[i]->nb_value;j++) {
          if ((categorical_process[i]->first_occurrence) &&
              (categorical_process[i]->first_occurrence[j])) {
            nb_plot_set++;
          }
          else if ((characteristics) && (j < characteristics->nb_value) &&
                   (characteristics->first_occurrence[j]->nb_element > 0)) {
            nb_plot_set++;
          }
        }
      }

      if ((categorical_process[i]->recurrence_time) || (characteristics)) {
        for (j = 0;j < categorical_process[i]->nb_value;j++) {
          if ((categorical_process[i]->recurrence_time) &&
              (categorical_process[i]->recurrence_time[j])) {
            nb_plot_set++;
          }
          else if ((characteristics) && (i < characteristics->nb_value) &&
                   (characteristics->recurrence_time[j]->nb_element > 0)) {
            nb_plot_set++;
          }
        }
      }

      if ((categorical_process[i]->sojourn_time) || (characteristics)) {
        for (j = 0;j < categorical_process[i]->nb_value;j++) {
          if ((categorical_process[i]->sojourn_time) &&
              (categorical_process[i]->sojourn_time[j])) {
            nb_plot_set++;
          }
          else if ((characteristics) && (i < characteristics->nb_value) &&
                   (characteristics->sojourn_time[j]->nb_element > 0)) {
            nb_plot_set++;
          }

/*          if ((characteristics) && (j < characteristics->nb_value) &&
              (characteristics->initial_run) &&
              (characteristics->initial_run[j]->nb_element > 0)) {
            nb_plot_set++;
          } */

          if ((characteristics) && (j < characteristics->nb_value) &&
              (characteristics->final_run[j]->nb_element > 0)) {
            nb_plot_set++;
          }
        }
      }

      if ((categorical_process[i]->nb_run) || (categorical_process[i]->nb_occurrence) ||
          ((characteristics) && (characteristics->nb_run) && (characteristics->nb_occurrence))) {
        for (j = 0;j < categorical_process[i]->nb_value;j++) {
          if (categorical_process[i]->nb_run) {
            nb_plot_set++;
          }
          else if ((characteristics) && (j < characteristics->nb_value) &&
                   (characteristics->nb_run) && (characteristics->nb_run[j]->nb_element > 0)) {
            nb_plot_set++;
          }

          if (categorical_process[i]->nb_occurrence) {
            nb_plot_set++;
          }
          else if ((characteristics) && (j < characteristics->nb_value) &&
                   (characteristics->nb_occurrence) &&
                   (characteristics->nb_occurrence[j]->nb_element > 0)) {
            nb_plot_set++;
          }
        }

        if ((characteristics) && (characteristics->nb_run) && (characteristics->nb_occurrence)) {
          nb_plot_set++;
        }
      }
    }

    if ((seq->observation_distribution) || (seq->observation_histogram)) {
      nb_plot_set += nb_state;
    }
    else {
      nb_plot_set++;
    }

    if ((categorical_process[i]) && (seq->marginal_distribution[variable])) {
      if ((categorical_process[i]->weight) &&
          (categorical_process[i]->mixture)) {
        nb_plot_set++;
      }
      if ((categorical_process[i]->restoration_weight) &&
          (categorical_process[i]->restoration_mixture)) {
        nb_plot_set++;
      }
    }

    if ((discrete_parametric_process[i]) && (seq->marginal_distribution[variable])) {
      if ((discrete_parametric_process[i]->weight) &&
          (discrete_parametric_process[i]->mixture)) {
        nb_plot_set += 2;
      }
      if ((discrete_parametric_process[i]->restoration_weight) &&
          (discrete_parametric_process[i]->restoration_mixture)) {
        nb_plot_set += 2;
      }
    }

    if ((continuous_parametric_process[i]) && ((seq->marginal_histogram[variable]) ||
         (seq->marginal_distribution[variable]))) {
      if (continuous_parametric_process[i]->weight) {
        nb_plot_set += 2;
      }
      if (continuous_parametric_process[i]->restoration_weight) {
        nb_plot_set += 2;
      }
    }
  }

  plot_set = new MultiPlotSet(nb_plot_set , nb_output_process + 1);
  plot_set->border = "15 lw 0";

  if ((seq) && (seq->type[0] == STATE)) {
    characteristics = seq->characteristics[0];
    length_distribution = seq->length_distribution;
  }
  else {
    characteristics = NULL;
  }

  index = 0;
  plot_set->variable_nb_viewpoint[0] = 0;
  state_process->plotable_write(*plot_set , index , 0 , NULL , NULL , characteristics ,
                                length_distribution , forward);

  if (seq) {
    length_distribution = seq->length_distribution;
  }

  for (i = 0;i < nb_output_process;i++) {
    if (seq) {
      switch (seq->type[0]) {
      case STATE :
        variable = i + 1;
        break;
      default :
        variable = i;
        break;
      }

      if (seq->observation_distribution) {
        observation_dist = seq->observation_distribution[variable];
      }
      marginal_dist = seq->marginal_distribution[variable];

      if (seq->observation_histogram) {
        observation_histo = seq->observation_histogram[variable];
      }
      marginal_histo = seq->marginal_histogram[variable];

      characteristics = seq->characteristics[variable];
    }

    if (categorical_process[i]) {
      plot_set->variable_nb_viewpoint[i] = 0;
      categorical_process[i]->plotable_write(*plot_set , index , i + 1 , observation_dist ,
                                             marginal_dist , characteristics ,
                                             length_distribution);
    }
    else if (discrete_parametric_process[i]) {
      discrete_parametric_process[i]->plotable_write(*plot_set , index , i + 1 , observation_dist ,
                                                     marginal_dist);
    }
    else {
      continuous_parametric_process[i]->plotable_write(*plot_set , index , i + 1 ,
                                                       observation_histo , observation_dist ,
                                                       marginal_histo , marginal_dist);
    }
  }

  return plot_set;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet SemiMarkov.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* SemiMarkov::get_plotable() const

{
  return get_plotable(semi_markov_data);
}


/*--------------------------------------------------------------*
 *
 *  Calcul du nombre de parametres independants d'un objet SemiMarkov.
 *
 *  argument : probabilite minimum.
 *
 *--------------------------------------------------------------*/

int SemiMarkov::nb_parameter_computation(double min_probability) const

{
  register int i;
  int nb_parameter = SemiMarkovChain::nb_parameter_computation(min_probability);


  for (i = 0;i < nb_output_process;i++) {
    if (categorical_process[i]) {
      nb_parameter += categorical_process[i]->nb_parameter_computation(min_probability);
    }
    else if (discrete_parametric_process[i]) {
      nb_parameter += discrete_parametric_process[i]->nb_parameter_computation();
    }
    else {
      nb_parameter += continuous_parametric_process[i]->nb_parameter_computation();
    }
  }

  return nb_parameter;
}


/*--------------------------------------------------------------*
 *
 *  Calcul d'une penalite adaptative.
 *
 *  arguments : flag semi-Markov cache, probabilite minimum.
 *  la proba minimum correspond au seuil a partir duquel on considere que le param vaut 0. e.g. 1e-5
 *
 *--------------------------------------------------------------*/

double SemiMarkov::penalty_computation(bool hidden , double min_probability) const

{
  register int i , j , k;
  int nb_parameter , sample_size;
  double sum , *memory , *state_marginal;
  double penalty = 0.;


  if (semi_markov_data) {
    if (hidden) {
      memory = memory_computation();

      state_marginal = new double[nb_state];

      switch (type) {

      case 'o' : {
        sum = 0.;
        for (i = 0;i < state_process->length->nb_value - 2;i++) {
          sum += (1. - state_process->length->cumul[i + 1]);
        }
        for (i = 0;i < nb_state;i++) {
          memory[i] /= sum;
        }

        for (i = 0;i < nb_state;i++) {
          state_marginal[i] = 0.;
        }
        for (i = 0;i < state_process->length->nb_value - 1;i++) {
          for (j = 0;j < nb_state;j++) {
            state_marginal[j] += state_process->index_value->point[j][i] *
                                 (1. - state_process->length->cumul[i]);
          }
        }

        sum = 0.;
        for (i = 0;i < nb_state;i++) {
          sum += state_marginal[i];
        }
        for (i = 0;i < nb_state;i++) {
          state_marginal[i] /= sum;
        }
        break;
      }

      case 'e' : {
        for (i = 0;i < nb_state;i++) {
          state_marginal[i] = initial[i];
        }
        break;
      }
      }
    }

    for (i = 0;i < nb_state;i++) {
      nb_parameter = 0;
      if (!hidden) {
        sample_size = 0;
      }
      for (j = 0;j < nb_state;j++) {
        if (transition[i][j] > min_probability) {
          nb_parameter++;
          if (!hidden) {
            sample_size += semi_markov_data->chain_data->transition[i][j];
          }
        }
      }

      nb_parameter--;

      if (nb_parameter > 0) {
        switch (hidden) {

        case false : {
          if (sample_size > 0) {
            penalty += nb_parameter * log((double)sample_size);
          }
          break;
        }

        case true : {
          if (memory[i] > 0.) {
            penalty += nb_parameter * log(memory[i] * semi_markov_data->cumul_length);
          }
          break;
        }
        }
      }
    }

    for (i = 0;i < nb_state;i++) {
      if (state_subtype[i] == SEMI_MARKOVIAN) {
        nb_parameter = state_process->sojourn_time[i]->nb_parameter_computation();
        if (state_process->sojourn_time[i]->inf_bound == 1) {
          nb_parameter--;
        }

        switch (hidden) {
        case false :
          penalty += nb_parameter *
                     log((double)semi_markov_data->marginal_distribution[0]->frequency[i]);
          break;
        case true :
          penalty += nb_parameter * log(state_marginal[i] * semi_markov_data->cumul_length);
          break;
        }
      }
    }

    for (i = 0;i < nb_output_process;i++) {
      if (categorical_process[i]) {
        for (j = 0;j < nb_state;j++) {
          nb_parameter = 0;
          for (k = 0;k < categorical_process[i]->nb_value;k++) {
            if (categorical_process[i]->observation[j]->mass[k] > min_probability) {
              nb_parameter++;
            }
          }

          nb_parameter--;

          if (nb_parameter > 0) {
            switch (hidden) {
              case false :
              penalty += nb_parameter *
                         log((double)semi_markov_data->marginal_distribution[0]->frequency[j]);
              break;
            case true :
              penalty += nb_parameter * log(state_marginal[j] * semi_markov_data->cumul_length);
              break;
            }
          }
        }
      }

      else if (discrete_parametric_process[i]) {
        for (j = 0;j < nb_state;j++) {
          nb_parameter = discrete_parametric_process[i]->observation[j]->nb_parameter_computation();

          switch (hidden) {
            case false :
            penalty += nb_parameter *
                       log((double)semi_markov_data->marginal_distribution[0]->frequency[j]);
            break;
          case true :
            penalty += nb_parameter * log(state_marginal[j] * semi_markov_data->cumul_length);
            break;
          }
        }
      }

      else {
        for (j = 0;j < nb_state;j++) {
          nb_parameter = continuous_parametric_process[i]->observation[j]->nb_parameter_computation();

          switch (hidden) {
            case false :
            penalty += nb_parameter *
                       log((double)semi_markov_data->marginal_distribution[0]->frequency[j]);
            break;
          case true :
            penalty += nb_parameter * log(state_marginal[j] * semi_markov_data->cumul_length);
            break;
          }
        }
      }
    }

    if (hidden) {
      delete [] memory;
      delete [] state_marginal;
    }
  }

  return penalty;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe SemiMarkovData.
 *
 *--------------------------------------------------------------*/

SemiMarkovData::SemiMarkovData()

{
  semi_markov = NULL;
  chain_data = NULL;

  likelihood = D_INF;
  restoration_likelihood = D_INF;
  sample_entropy = D_DEFAULT;

  posterior_probability = NULL;
  entropy = NULL;
  nb_state_sequence = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe SemiMarkovData.
 *
 *  arguments : loi empirique des longueurs des sequences, nombre de variables,
 *              type de chaque variable, flag initialisation.
 *
 *--------------------------------------------------------------*/

SemiMarkovData::SemiMarkovData(const FrequencyDistribution &ilength_distribution , int inb_variable ,
                               int *itype , bool init_flag)
:MarkovianSequences(ilength_distribution , inb_variable , itype , init_flag)

{
  semi_markov = NULL;
  chain_data = NULL;

  likelihood = D_INF;
  restoration_likelihood = D_INF;
  sample_entropy = D_DEFAULT;

  posterior_probability = NULL;
  entropy = NULL;
  nb_state_sequence = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet SemiMarkovData a partir d'un objet MarkovianSequences
 *  avec ajout d'une variable d'etat.
 *
 *  argument : reference sur un objet MarkovianSequences.
 *
 *--------------------------------------------------------------*/

SemiMarkovData::SemiMarkovData(const MarkovianSequences &seq)
:MarkovianSequences(seq , 'a' , UNCHANGED)

{
  semi_markov = NULL;
  chain_data = NULL;

  likelihood = D_INF;
  restoration_likelihood = D_INF;
  sample_entropy = D_DEFAULT;

  posterior_probability = NULL;
  entropy = NULL;
  nb_state_sequence = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet SemiMarkovData a partir d'un objet MarkovianSequences.
 *
 *  arguments : reference sur un objet MarkovianSequences, type de transformation
 *              ('c' : copie, 'a' : ajout d'une variable d'etat),
 *              ajout/suppression des lois empiriques de temps de sejour initial.
 *
 *--------------------------------------------------------------*/

SemiMarkovData::SemiMarkovData(const MarkovianSequences &seq , char transform ,
                               bool initial_run_flag)
:MarkovianSequences(seq , transform , (initial_run_flag ? ADD_INITIAL_RUN : REMOVE_INITIAL_RUN))

{
  semi_markov = NULL;
  chain_data = NULL;

  likelihood = D_INF;
  restoration_likelihood = D_INF;
  sample_entropy = D_DEFAULT;

  posterior_probability = NULL;
  entropy = NULL;
  nb_state_sequence = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet SemiMarkovData.
 *
 *  arguments : reference sur un objet SemiMarkovData,
 *              flag copie de l'objet SemiMarkov.
 *
 *--------------------------------------------------------------*/

void SemiMarkovData::copy(const SemiMarkovData &seq , bool model_flag)

{
  register int i;


  if ((model_flag) && (seq.semi_markov)) {
    semi_markov = new SemiMarkov(*(seq.semi_markov) , false);
  }
  else {
    semi_markov = NULL;
  }

  if (seq.chain_data) {
    chain_data = new ChainData(*(seq.chain_data));
  }
  else {
    chain_data = NULL;
  }

  likelihood = seq.likelihood;
  restoration_likelihood = seq.restoration_likelihood;
  sample_entropy = seq.sample_entropy;

  if (seq.posterior_probability) {
    posterior_probability = new double[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      posterior_probability[i] = seq.posterior_probability[i];
    }
  }
  else {
    posterior_probability = NULL;
  }

  if (seq.entropy) {
    entropy = new double[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      entropy[i] = seq.entropy[i];
    }
  }
  else {
    entropy = NULL;
  }

  if (seq.nb_state_sequence) {
    nb_state_sequence = new double[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      nb_state_sequence[i] = seq.nb_state_sequence[i];
    }
  }
  else {
    nb_state_sequence = NULL;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe SemiMarkovData.
 *
 *--------------------------------------------------------------*/

SemiMarkovData::~SemiMarkovData()

{
  delete semi_markov;
  delete chain_data;

  delete [] posterior_probability;
  delete [] entropy;
  delete [] nb_state_sequence;
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe SemiMarkovData.
 *
 *  argument : reference sur un objet SemiMarkovData.
 *
 *--------------------------------------------------------------*/

SemiMarkovData& SemiMarkovData::operator=(const SemiMarkovData &seq)

{
  if (&seq != this) {
    delete semi_markov;
    delete chain_data;

    delete [] posterior_probability;
    delete [] entropy;
    delete [] nb_state_sequence;

    remove();
    Sequences::remove();

    Sequences::copy(seq);
    MarkovianSequences::copy(seq);
    copy(seq);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Extraction d'une loi empirique.
 *
 *  arguments : reference sur un objet StatError, type de loi empirique,
 *              variable, etat ou observation.
 *
 *--------------------------------------------------------------*/

DiscreteDistributionData* SemiMarkovData::extract(StatError &error , int type ,
                                                  int variable , int value) const

{
  bool status = true;
  Distribution *pdist;
  DiscreteParametric *pparam;
  FrequencyDistribution *phisto;
  DiscreteDistributionData *histo;
  CategoricalSequenceProcess *process;


  histo = NULL;
  error.init();

  phisto = NULL;

  if (type == OBSERVATION) {
    if ((variable < 2) || (variable > nb_variable)) {
      status = false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
    }

    else {
      variable--;

      if ((value < 0) || (value >= marginal_distribution[0]->nb_value)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_STATE] << " " << value << " "
                      << STAT_error[STATR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }

      else {
        if (!observation_distribution[variable]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                        << STAT_error[STATR_VARIABLE_TYPE];
          error.correction_update((error_message.str()).c_str() , STAT_variable_word[INT_VALUE]);
        }

        else {
          phisto = observation_distribution[variable][value];

          if (phisto->nb_element == 0) {
            status = false;
            error.update(STAT_error[STATR_EMPTY_SAMPLE]);
          }
        }
      }
    }
  }

  else {
    if ((variable < 1) || (variable > nb_variable)) {
      status = false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
    }

    else {
      variable--;

      if (!characteristics[variable]) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                      << SEQ_error[SEQR_CHARACTERISTICS_NOT_COMPUTED];
        error.update((error_message.str()).c_str());
      }

      else if ((value < 0) || (value >= marginal_distribution[variable]->nb_value)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[variable == 0 ? STATL_STATE : STATL_OUTPUT] << " "
                      << value << " " << STAT_error[STATR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }

      if (status) {
        switch (type) {

        case FIRST_OCCURRENCE : {
          phisto = characteristics[variable]->first_occurrence[value];
          break;
        }

        case RECURRENCE_TIME : {
          phisto = characteristics[variable]->recurrence_time[value];
          break;
        }

        case SOJOURN_TIME : {
          phisto = characteristics[variable]->sojourn_time[value];
          break;
        }

        case INITIAL_RUN : {
          if (characteristics[variable]->initial_run) {
            phisto = characteristics[variable]->initial_run[value];
          }
          else {
            status = false;
            error.update(STAT_error[STATR_NON_EXISTING_FREQUENCY_DISTRIBUTION]);
          }
          break;
        }

        case FINAL_RUN : {
          phisto = characteristics[variable]->final_run[value];
          break;
        }

        case NB_RUN : {
          phisto = characteristics[variable]->nb_run[value];
          break;
        }

        case NB_OCCURRENCE : {
          phisto = characteristics[variable]->nb_occurrence[value];
          break;
        }
        }

        if ((phisto) && (phisto->nb_element == 0)) {
          status = false;
          error.update(STAT_error[STATR_EMPTY_SAMPLE]);
        }
      }
    }
  }

  if (status) {
    if (variable == 0) {
      process = semi_markov->state_process;
    }
    else {
      process = semi_markov->categorical_process[variable - 1];
    }

    pdist = NULL;
    pparam = NULL;

    switch (type) {

    case OBSERVATION : {
      if (semi_markov->categorical_process[variable - 1]) {
        pdist = semi_markov->categorical_process[variable - 1]->observation[value];
      }
      else if (semi_markov->discrete_parametric_process[variable - 1]) {
        pparam = semi_markov->discrete_parametric_process[variable - 1]->observation[value];
      }
      break;
    }

    case FIRST_OCCURRENCE : {
      pdist = process->first_occurrence[value];
      break;
    }

    case RECURRENCE_TIME : {
      pdist = process->recurrence_time[value];
      break;
    }

    case SOJOURN_TIME : {
      pparam = process->sojourn_time[value];
      break;
    }

    case INITIAL_RUN : {
      if ((variable == 0) && (semi_markov->forward)) {
        pdist = semi_markov->forward[value];
      }
      break;
    }

    case FINAL_RUN : {
      if ((variable == 0) && (semi_markov->forward)) {
        pdist = semi_markov->forward[value];
      }
      break;
    }

    case NB_RUN : {
      pdist = process->nb_run[value];
      break;
    }

    case NB_OCCURRENCE : {
      pdist = process->nb_occurrence[value];
      break;
    }
    }

    if (pdist) {
      histo = new DiscreteDistributionData(*phisto , pdist);
    }
    else {
      histo = new DiscreteDistributionData(*phisto , pparam);
    }
  }

  return histo;
}


/*--------------------------------------------------------------*
 *
 *  Suppression du parametre d'index.
 *
 *  argument : reference sur un objet StatError.
 *
 *--------------------------------------------------------------*/

SemiMarkovData* SemiMarkovData::remove_index_parameter(StatError &error) const

{
  SemiMarkovData *seq;


  error.init();

  if (!index_parameter) {
    seq = NULL;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }
  else {
    seq = new SemiMarkovData(*this , true , 'm');
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet SemiMarkovData avec transformation du parametre d'index implicite
 *  en parametre d'index explicite.
 *
 *  argument : reference sur un objet StatError.
 *
 *--------------------------------------------------------------*/

SemiMarkovData* SemiMarkovData::explicit_index_parameter(StatError &error) const

{
  SemiMarkovData *seq;


  error.init();

  if (index_parameter) {
    seq = NULL;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }
  else {
    seq = new SemiMarkovData(*this , true , 'e');
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Construction des variables auxilliaires correspondant a
 *  la restauration des sequences d'etats optimales.
 *
 *  argument : reference sur un objet StatError.
 *
 *--------------------------------------------------------------*/

MarkovianSequences* SemiMarkovData::build_auxiliary_variable(StatError &error) const

{
  bool status = true;
  register int i;
  MarkovianSequences *seq;


  seq = NULL;
  error.init();

  if (type[0] != STATE) {
    status = false;
    ostringstream error_message;
    error_message << STAT_label[STATL_VARIABLE] << " 1: "
                  << STAT_error[STATR_VARIABLE_TYPE];
    error.correction_update((error_message.str()).c_str() , STAT_variable_word[STATE]);
  }

  for (i = 0;i < semi_markov->nb_output_process;i++) {
    if (((semi_markov->discrete_parametric_process) && (semi_markov->discrete_parametric_process[i])) ||
        ((semi_markov->continuous_parametric_process) && (semi_markov->continuous_parametric_process[i]))) {
      break;
    }
  }

  if (i == semi_markov->nb_output_process) {
    status = false;
    error.update(SEQ_error[SEQR_PARAMETRIC_PROCESS]);
  }

  if (status) {
    seq = MarkovianSequences::build_auxiliary_variable(semi_markov->discrete_parametric_process ,
                                                       semi_markov->continuous_parametric_process);
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet SemiMarkovData.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& SemiMarkovData::ascii_write(ostream &os , bool exhaustive) const

{
  if (semi_markov) {
    semi_markov->ascii_write(os , this , exhaustive , false ,
                             sequence_analysis::test_hidden(semi_markov->nb_output_process , semi_markov->categorical_process));
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet SemiMarkovData dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool SemiMarkovData::ascii_write(StatError &error , const char *path ,
                                 bool exhaustive) const

{
  bool status = false;


  if (semi_markov) {
    ofstream out_file(path);

    error.init();

    if (!out_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_NAME]);
    }

    else {
      status = true;
      semi_markov->ascii_write(out_file , this , exhaustive , true ,
                               sequence_analysis::test_hidden(semi_markov->nb_output_process , semi_markov->categorical_process));
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet SemiMarkovData.
 *
 *  arguments : stream, format (ligne/colonne), flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& SemiMarkovData::ascii_data_write(ostream &os , char format , bool exhaustive) const

{
  MarkovianSequences::ascii_write(os , exhaustive , false);
  ascii_print(os , format , false , posterior_probability , entropy , nb_state_sequence);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet SemiMarkovData dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              format (ligne/colonne), flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool SemiMarkovData::ascii_data_write(StatError &error , const char *path ,
                                      char format , bool exhaustive) const

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
    if (format != 'a') {
      MarkovianSequences::ascii_write(out_file , exhaustive , true);
    }
    ascii_print(out_file , format , true , posterior_probability , entropy , nb_state_sequence);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet SemiMarkovData dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool SemiMarkovData::spreadsheet_write(StatError &error , const char *path) const

{
  bool status = false;


  if (semi_markov) {
    ofstream out_file(path);

    error.init();

    if (!out_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_NAME]);
    }

    else {
      status = true;
      semi_markov->spreadsheet_write(out_file , this ,
                                     sequence_analysis::test_hidden(semi_markov->nb_output_process , semi_markov->categorical_process));
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet SemiMarkovData.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool SemiMarkovData::plot_write(StatError &error , const char *prefix ,
                                const char *title) const

{
  bool status = false;


  if (semi_markov) {
    status = semi_markov->plot_write(prefix , title , this);

    error.init();

    if (!status) {
      error.update(STAT_error[STATR_FILE_PREFIX]);
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet SemiMarkovData.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* SemiMarkovData::get_plotable() const

{
  MultiPlotSet *plot_set;


  if (semi_markov) {
    plot_set = semi_markov->get_plotable(this);
  }
  else {
    plot_set = NULL;
  }

  return plot_set;
}


};  // namespace sequence_analysis
