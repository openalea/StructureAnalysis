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



#include <math.h>
#include <sstream>
#include <iomanip>
#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"
#include "stat_tool/stat_tools.h"
#include "stat_tool/distribution.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/stat_label.h"
#include "sequences.h"
#include "variable_order_markov.h"
#include "sequence_label.h"
// #include "tool/config.h"

using namespace std;


extern int column_width(int value);
extern int column_width(int nb_value , const double *value , double scale = 1.);
extern char* label(const char *file_name);

extern double standard_normal_value_computation(double critical_probability);



/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe VariableOrderMarkov.
 *
 *--------------------------------------------------------------*/

VariableOrderMarkov::VariableOrderMarkov()

{
  nb_iterator = 0;
  markov_data = NULL;

  max_order = 0;

  memory_type = NULL;
  order = NULL;
  state = NULL;

  parent = NULL;
  child = NULL;

  next = NULL;
  nb_memory = NULL;
  previous = NULL;

  nb_output_process = 0;
  nonparametric_process = NULL;
  parametric_process = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe VariableOrderMarkov.
 *
 *  arguments : type, nombre d'etats, nombres de memoires.
 *
 *--------------------------------------------------------------*/

VariableOrderMarkov::VariableOrderMarkov(char itype , int inb_state , int inb_row)
:Chain(itype , inb_state , inb_row , true)

{
  register int i;


  nb_iterator = 0;
  markov_data = NULL;

  max_order = 0;

  memory_type = new int[nb_row];
  order = new int[nb_row];
  state = new int*[nb_row];
  parent = new int[nb_row];
  child = new int*[nb_row];

  for (i = 0;i < nb_row;i++) {
    state[i] = NULL;
    child[i] = NULL;
  }

  next = NULL;
  nb_memory = NULL;
  previous = NULL;

  nb_output_process = 0;
  nonparametric_process = new NonparametricSequenceProcess*[1];
  nonparametric_process[0] = new NonparametricSequenceProcess(nb_state , nb_state);

  parametric_process = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe VariableOrderMarkov.
 *
 *  arguments : type, nombre d'etats, nombres de memoires, ordre maximum.
 *
 *--------------------------------------------------------------*/

VariableOrderMarkov::VariableOrderMarkov(char itype , int inb_state ,
                                         int inb_row , int imax_order)
:Chain(itype , inb_state , inb_row , true)

{
  register int i;


  nb_iterator = 0;
  markov_data = NULL;

  max_order = imax_order;

  memory_type = new int[nb_row];
  order = new int[nb_row];
  state = new int*[nb_row];
  parent = new int[nb_row];
  child = new int*[nb_row];

  for (i = 0;i < nb_row;i++) {
    state[i] = new int[max_order];
    child[i] = NULL;
  }

  next = NULL;
  nb_memory = NULL;
  previous = NULL;

  nb_output_process = 0;
  nonparametric_process = new NonparametricSequenceProcess*[1];
  nonparametric_process[0] = new NonparametricSequenceProcess(nb_state , nb_state);

  parametric_process = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet VariableOrderMarkov d'ordre fixe.
 *
 *  arguments : type, nombre d'etats, ordre, flag initialisation,
 *              nombre de processus d'observation,
 *              nombre de valeurs observees par processus.
 *
 *--------------------------------------------------------------*/

VariableOrderMarkov::VariableOrderMarkov(char itype , int inb_state ,
                                         int iorder , bool init_flag ,
                                         int inb_output_process , int nb_value)
:Chain(itype , inb_state , (int)(pow((double)inb_state , iorder + 1) - 1) / (inb_state - 1) , init_flag)

{
  register int i , j;


  nb_iterator = 0;
  markov_data = NULL;

  max_order = iorder;

  memory_type = new int[nb_row];
  order = new int[nb_row];
  state = new int*[nb_row];
  parent = new int[nb_row];
  child = new int*[nb_row];

  // racine (ordre 0)

  memory_type[0] = NON_TERMINAL;
  order[0] = 0;
  state[0] = NULL;
  parent[0] = -1;
  child[0] = new int[nb_state];

  for (i = 1;i < nb_row;i++) {

    // cas augmentation de l'ordre

    if (order[i - 1] < max_order) {
      order[i] = order[i - 1] + 1;
      if (order[i] < max_order) {
        memory_type[i] = NON_TERMINAL;
      }
      else {
        memory_type[i] = TERMINAL;
      }

      state[i] = new int[order[i]];
      for (j = 0;j < order[i - 1];j++) {
        state[i][j] = state[i - 1][j];
      }
      state[i][order[i] - 1] = 0;

      parent[i] = i - 1;
      child[i - 1][0] = i;
    }

    else {

      // cas ordre (maximum) stable

      if (state[i - 1][order[i - 1] - 1] < nb_state - 1) {
        memory_type[i] = TERMINAL;
        order[i] = max_order;
        state[i] = new int[order[i]];
        for (j = 0;j < order[i] - 1;j++) {
          state[i][j] = state[i - 1][j];
        }
        state[i][order[i] - 1] = state[i - 1][order[i] - 1] + 1;

        parent[i] = i - state[i][order[i] - 1] - 1;
        child[parent[i]][state[i][order[i] - 1]] = i;
      }

      // cas diminution de l'ordre

      else {
        memory_type[i] = NON_TERMINAL;

        for (j = order[i - 1] - 2;j >= 0;j--) {
          if (state[i - 1][j] != nb_state - 1) {
            break;
          }
        }
        order[i] = j + 1;

        state[i] = new int[order[i]];
        for (j = 0;j < order[i] - 1;j++) {
          state[i][j] = state[i - 1][j];
        }
        state[i][order[i] - 1] = state[i - 1][order[i] - 1] + 1;

        // recherche du noeud pere

        find_parent_memory(i);
      }
    }

    if (memory_type[i] == NON_TERMINAL) {
      child[i] = new int[nb_state];
    }
    else {
      child[i] = NULL;
    }
  }

  // calcul des transitions entre memoires terminales

  next = NULL;
  nb_memory = NULL;
  previous = NULL;

  build_memory_transition();

  nb_output_process = inb_output_process;

  nonparametric_process = new NonparametricSequenceProcess*[nb_output_process + 1];
  nonparametric_process[0] = new NonparametricSequenceProcess(nb_state , nb_state , false);
  if (nb_output_process == 1) {
    nonparametric_process[1] = new NonparametricSequenceProcess(nb_state , nb_value , true);
  }

  parametric_process = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Completion de l'arborescence des memoires.
 *
 *  argument : reference sur un objet VariableOrderMarkov.
 *
 *--------------------------------------------------------------*/

void VariableOrderMarkov::memory_tree_completion(const VariableOrderMarkov &markov)

{
  bool prefix;
  register int i , j , k , m;
  int bnb_memory , border , *markov_next , *completion_next;
  VariableOrderMarkov *completion;


  bnb_memory = markov.nb_row;
  for (i = 1;i < markov.nb_row;i++) {
    if (markov.order[i] == markov.max_order) {
      bnb_memory--;
    }
  }

  completion = new VariableOrderMarkov(markov.type , markov.nb_state ,
                                       (int)(pow((double)markov.nb_state , markov.max_order) - 1) / (markov.nb_state - 1) - bnb_memory);
  completion->nb_row = 0;

  for (i = 1;i < markov.nb_row;i++) {
    if (markov.order[i] > 1) {

      // recherche du noeud correspondant au plus grand prefixe propre

      prefix = false;

      for (j = 0;j < markov.nb_row;j++) {
        if (markov.order[j] == markov.order[i] - 1) {
          for (k = 0;k < markov.order[j];k++) {
            if (markov.state[j][k] != markov.state[i][k + 1]) {
              break;
            }
          }

          if (k == markov.order[j]) {
            prefix = true;
            break;
          }
        }
      }

      if (!prefix) {
        for (j = 0;j < completion->nb_row;j++) {
          if (completion->order[j] == markov.order[i] - 1) {
            for (k = 0;k < completion->order[j];k++) {
              if (completion->state[j][k] != markov.state[i][k + 1]) {
                break;
              }
            }

            if (k == completion->order[j]) {
              prefix = true;
              break;
            }
          }
        }

        if (!prefix) {

          // construction du noeud correspondant au plus grand prefixe propre

          completion->order[completion->nb_row] = markov.order[i] - 1;
          completion->state[completion->nb_row] = new int[completion->order[completion->nb_row]];
          for (j = 0;j < completion->order[completion->nb_row];j++) {
            completion->state[completion->nb_row][j] = markov.state[i][j + 1];
          }
          for (j = 0;j < markov.nb_state;j++) {
            completion->transition[completion->nb_row][j] = markov.transition[i][j];
          }
          (completion->nb_row)++;

          for (j = completion->order[completion->nb_row - 1] - 1;j >= 2;j--) {
            prefix = false;

            for (k = 0;k < markov.nb_row;k++) {
              if (markov.order[k] == j) {
                for (m = 0;m < markov.order[k];m++) {
                  if (markov.state[k][m] != completion->state[completion->nb_row - 1][m + 1]) {
                    break;
                  }
                }

                if (m == markov.order[k]) {
                  prefix = true;
                  break;
                }
              }
            }

            if (!prefix) {
              for (k = 0;k < completion->nb_row - 1;k++) {
                if (completion->order[k] == j) {
                  for (m = 0;m < completion->order[k];m++) {
                    if (completion->state[k][m] != completion->state[completion->nb_row - 1][m + 1]) {
                      break;
                    }
                  }

                  if (m == completion->order[k]) {
                    prefix = true;
                    break;
                  }
                }
              }

              if (!prefix) {

                // construction du noeud correspondant au plus grand prefixe propre

                completion->order[completion->nb_row] = j;
                completion->state[completion->nb_row] = new int[completion->order[completion->nb_row]];
                for (k = 0;k < completion->order[completion->nb_row];k++) {
                  completion->state[completion->nb_row][k] = completion->state[completion->nb_row - 1][k + 1];
                }
                for (k = 0;k < markov.nb_state;k++) {
                  completion->transition[completion->nb_row][k] = completion->transition[completion->nb_row - 1][k];
                }
                (completion->nb_row)++;
              }

              else {
                break;
              }
            }
          }
        }
      }
    }
  }

# ifdef DEBUG
  {
    cout << "\n";
    for (i = 0;i < completion->nb_row;i++) {
      cout << completion->order[i] << " | ";
      for (j = completion->order[i] - 1;j >= 0;j--) {
        cout << completion->state[i][j] << " ";
      }
      cout << endl;
    }
  }
# endif

  // recherche du noeud suivant dans l'ordonancement de l'arborescence completee

  markov_next = new int[markov.nb_row];
  completion_next = new int[completion->nb_row];

  for (i = 0;i < markov.nb_row;i++) {
    markov_next[i] = I_DEFAULT;

    if ((markov.memory_type[i] == TERMINAL) && (markov.order[i] < markov.max_order - 1)) {

      // recherche du 1er fils (etat 0) du noeud terminal

      for (j = 0;j < completion->nb_row;j++) {
        if (completion->order[j] == markov.order[i] + 1) {
          for (k = 0;k < markov.order[i];k++) {
            if (completion->state[j][k] != markov.state[i][k]) {
              break;
            }
          }

          if ((k == markov.order[i]) && (completion->state[j][completion->order[j] - 1] == 0)) {
            markov_next[i] = j;
            break;
          }
        }
      }
    }
  }

  for (i = 0;i < completion->nb_row;i++) {
    completion_next[i] = I_DEFAULT;

    // recherche du 1er fils (etat 0) du noeud cree

    if (completion->order[i] < markov.max_order - 1) {
      for (j = 0;j < completion->nb_row;j++) {
        if (completion->order[j] == completion->order[i] + 1) {
          for (k = 0;k < completion->order[i];k++) {
            if (completion->state[j][k] != completion->state[i][k]) {
              break;
            }
          }

          if ((k == completion->order[i]) && (completion->state[j][completion->order[j] - 1] == 0)) {
            completion_next[i] = j;
            break;
          }
        }
      }
    }

    if ((completion->order[i] == markov.max_order - 1) || (j == completion->nb_row)) {
      if (completion->state[i][completion->order[i] - 1] < markov.nb_state - 1) {

        // recherche du frere (etat suivant) du noeud cree

        for (j = 0;j < completion->nb_row;j++) {
          if (completion->order[j] == completion->order[i]) {
            for (k = 0;k < completion->order[i] - 1;k++) {
              if (completion->state[j][k] != completion->state[i][k]) {
                break;
              }
            }

            if ((k == completion->order[i] - 1) &&
                (completion->state[j][completion->order[j] - 1] == completion->state[i][completion->order[i] - 1] + 1)) {
              completion_next[i] = j;
              break;
            }
          }
        }
      }

      else {
        for (j = completion->order[i] - 2;j >= 0;j--) {
          if (completion->state[i][j] != markov.nb_state - 1) {
            break;
          }
        }

        border = j + 1;
        for (j = 0;j < completion->nb_row;j++) {
          if (completion->order[j] == border) {
            for (k = 0;k < border - 1;k++) {
              if (completion->state[j][k] != completion->state[i][k]) {
                break;
              }
            }

            if ((k == border - 1) &&
                (completion->state[j][border - 1] == completion->state[i][border - 1] + 1)) {
              completion_next[i] = j;
              break;
            }
          }
        }
      }
    }
  }

# ifdef DEBUG
  {
    cout << "\n";
    for (i = 0;i < markov.nb_row;i++) {
      for (j = markov.order[i] - 1;j >= 0;j--) {
        cout << markov.state[i][j] << " ";
      }
      if (markov_next[i] != I_DEFAULT) {
        cout << " -> ";
        for (j = completion->order[markov_next[i]] - 1;j >= 0;j--) {
          cout << completion->state[markov_next[i]][j] << " ";
        }
      }
      cout << endl;
    }

    cout << "\n";
    for (i = 0;i < completion->nb_row;i++) {
      for (j = completion->order[i] - 1;j >= 0;j--) {
        cout << completion->state[i][j] << " ";
      }
      if (completion_next[i] != I_DEFAULT) {
        cout << " -> ";
        for (j = completion->order[completion_next[i]] - 1;j >= 0;j--) {
          cout << completion->state[completion_next[i]][j] << " ";
        }
      }
      cout << endl;
    }
  }
# endif

  // copie des parametres

  type = markov.type;
  nb_state = markov.nb_state;
  nb_row = markov.nb_row + completion->nb_row;

  if (markov.nb_component > 0) {
    accessibility = new bool*[nb_state];
    for (i = 0;i < nb_state;i++) {
      accessibility[i] = new bool[nb_state];
      for (j = 0;j < nb_state;j++) {
        accessibility[i][j] = markov.accessibility[i][j];
      }
    }

    nb_component = markov.nb_component;
    component_nb_state = new int[nb_component];
    component = new int*[nb_component];

    for (i = 0;i < nb_component;i++) {
      component_nb_state[i] = markov.component_nb_state[i];
      component[i] = new int[component_nb_state[i]];
      for (j = 0;j < component_nb_state[i];j++) {
        component[i][j] = markov.component[i][j];
      }
    }

    state_type = new char[nb_state];
    for (i = 0;i < nb_state;i++) {
      state_type[i] = markov.state_type[i];
    }
  }

  else {
    accessibility = NULL;
    nb_component = 0;
    component_nb_state = NULL;
    component = NULL;
    state_type = NULL;
  }

  initial = new double[type == 'o' ? nb_state : nb_row];

  if (type == 'o') {
    for (i = 0;i < nb_state;i++) {
      initial[i] = markov.initial[i];
    }
  }

  transition = new double*[nb_row];
  for (i = 0;i < nb_row;i++) {
    transition[i] = new double[nb_state];
  }

  cumul_initial = NULL;
  cumul_transition = NULL;

  nb_iterator = 0;

  if (markov.markov_data) {
    markov_data = new VariableOrderMarkovData(*(markov.markov_data) , false);
  }
  else {
    markov_data = NULL;
  }

  max_order = markov.max_order;

  memory_type = new int[nb_row];
  order = new int[nb_row];
  state = new int*[nb_row];
  parent = new int[nb_row];
  child = new int*[nb_row];

  next = NULL;
  nb_memory = NULL;
  previous = NULL;

  // insertion des memoires dans l'arborescence

  i = 0;
  for (j = 0;j < markov.nb_row;j++) {
    for (k = 0;k < nb_state;k++) {
      transition[i][k] = markov.transition[j][k];
    }

    memory_type[i] = markov.memory_type[j];

    order[i] = markov.order[j];
    state[i] = new int[order[i]];
    for (k = 0;k < order[i];k++) {
      state[i][k] = markov.state[j][k];
    }

    parent[i] = markov.parent[j];

    if ((markov.memory_type[j] == NON_TERMINAL) || (markov_next[j] != I_DEFAULT)) {
      child[i] = new int[nb_state];

      if (markov.memory_type[j] == NON_TERMINAL) {
        for (k = 0;k < nb_state;k++) {
          child[i][k] = markov.child[j][k];
        }
      }
    }

    else {
      child[i] = NULL;
    }

    i++;

    if (markov_next[j] != I_DEFAULT) {
      k = markov_next[j];
      bnb_memory = i;

      do {
        memory_type[i] = COMPLETION;

        order[i] = completion->order[k];
        state[i] = new int[order[i]];
        for (m = 0;m < order[i];m++) {
          state[i][m] = completion->state[k][m];
        }

        find_parent_memory(i);

        for (m = 0;m < nb_state;m++) {
          transition[i][m] = transition[parent[i]][m];
        }

        if ((completion_next[k] != I_DEFAULT) &&
            (completion->order[completion_next[k]] == completion->order[k] + 1)) {
          child[i] = new int[nb_state];
        }
        else {
          child[i] = NULL;
        }

        i++;
        k = completion_next[k];
      }
      while (k != I_DEFAULT);

      // mise a jour des relations parent/enfants

      for (k = 0;k < bnb_memory;k++) {
        if (memory_type[k] == NON_TERMINAL) {
          for (m = 0;m < nb_state;m++) {
            if (child[k][m] >= bnb_memory) {
              child[k][m] += i - bnb_memory;
            }
          }
        }
      }

      for (k = j + 1;k < markov.nb_row;k++) {
        if (markov.parent[k] >= bnb_memory) {
          markov.parent[k] += i - bnb_memory;
        }
        if (markov.memory_type[k] == NON_TERMINAL) {
          for (m = 0;m < nb_state;m++) {
            if (markov.child[k][m] >= bnb_memory) {
              markov.child[k][m] += i - bnb_memory;
            }
          }
        }
      }
    }
  }

# ifdef DEBUG
  {
    cout << "\n";
    cout << "Suffix free? " << (check_free_suffix() ? "True" : "False") << endl;

    for (i = 0;i < nb_row;i++) {
      cout << i << "  ";
      for (j = max_order - 1;j >= order[i];j--) {
        cout << "  ";
      }
      for (j = order[i] - 1;j >= 0;j--) {
        cout << state[i][j] << " ";
      }

      switch (memory_type[i]) {
      case NON_TERMINAL :
        cout << "  " << SEQ_label[SEQL_NON_TERMINAL];
        break;
      case TERMINAL :
        cout << "      " << SEQ_label[SEQL_TERMINAL];
        break;
      case COMPLETION :
        cout << "    " << SEQ_label[SEQL_COMPLETION];
        break;
      }

      cout << " | " << parent[i];

      if (child[i]) {
        cout << " |";
        for (j = 0;j < nb_state;j++) {
          cout << " " << child[i][j];
        }
      }

      else {
        cout << "  ";
        for (j = 0;j < nb_state;j++) {
          cout << "  ";
        }
      }

      cout << endl;
    }
  }
# endif

  delete completion;

  delete [] markov_next;
  delete [] completion_next;

  build_memory_transition();
  build_previous_memory();
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe VariableOrderMarkov avec completion
 *  de l'arborescence des memoires.
 *
 *  arguments : reference sur un objet VariableOrderMarkov,
 *              nombre de processus d'observation,
 *              nombre de valeurs observees par processus.
 *
 *--------------------------------------------------------------*/

VariableOrderMarkov::VariableOrderMarkov(const VariableOrderMarkov &markov ,
                                         int inb_output_process , int nb_value)

{
  memory_tree_completion(markov);

  nb_output_process = inb_output_process;

  nonparametric_process = new NonparametricSequenceProcess*[nb_output_process + 1];
  nonparametric_process[0] = new NonparametricSequenceProcess(nb_state , nb_state , false);
  if (nb_output_process == 1) {
    nonparametric_process[1] = new NonparametricSequenceProcess(nb_state , nb_value , true);
  }

  parametric_process = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe VariableOrderMarkov avec completion
 *  de l'arborescence des memoires.
 *
 *  arguments : reference sur un objet VariableOrderMarkov,
 *              nombre de processus d'observation,
 *              nombre de valeurs observees par processus.
 *
 *--------------------------------------------------------------*/

/* VariableOrderMarkov::VariableOrderMarkov(const VariableOrderMarkov &markov ,
                                         int inb_output_process , int *nb_value)

{
  register int i;


  memory_tree_completion(markov);

  nb_output_process = inb_output_process;

  nonparametric_process = new NonparametricSequenceProcess*[nb_output_process + 1];
  nonparametric_process[0] = new NonparametricSequenceProcess(nb_state , nb_state , false);
  parametric_process = new DiscreteParametricProcess*[nb_output_process + 1];
  parametric_process[0] = NULL;

  for (i = 1;i <= nb_output_process;i++) {
    if (*nb_value <= NB_OUTPUT) {
      nonparametric_process[i] = new NonparametricSequenceProcess(nb_state , *nb_value++ , true);
      parametric_process[i] = NULL;
    }
    else {
      nonparametric_process[i] = NULL;
      parametric_process[i] = new DiscreteParametricProcess(nb_state , (int)(*nb_value++ * SAMPLE_NB_VALUE_COEFF));
    }
  }
} */


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe VariableOrderMarkov avec completion
 *  de l'arborescence des memoires.
 *
 *  arguments : pointeur sur un objet VariableOrderMarkov,
 *              pointeur sur un objet NonparametricSequenceProcess,
 *              longueur des sequences.
 *
 *--------------------------------------------------------------*/

VariableOrderMarkov::VariableOrderMarkov(const VariableOrderMarkov *pmarkov ,
                                         const NonparametricProcess *pobservation , int length)

{
  register int i;
  int nb_terminal;


  memory_tree_completion(*pmarkov);

  switch (type) {

  case 'o' : {
    non_terminal_transition_probability_computation();
    break;
  }

  case 'e' : {
    nb_terminal = (nb_row - 1) * (nb_state - 1) / nb_state + 1;

    for (i = 1;i < nb_row;i++) {
      if (!child[i]) {
        initial[i] = 1. / (double)nb_terminal;
      }
      else {
        initial[i] = 0.;
      }
    }

    initial_probability_computation();
    break;
  }
  }

  nb_output_process = (pobservation ? 1 : 0);
  nonparametric_process = new NonparametricSequenceProcess*[nb_output_process + 1];
  nonparametric_process[0] = new NonparametricSequenceProcess(nb_state , nb_state);
  if (pobservation) {
    nonparametric_process[1] = new NonparametricSequenceProcess(*pobservation);
  }

  parametric_process = NULL;

  characteristic_computation(length , true);
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet VariableOrderMarkov.
 *
 *  arguments : reference sur un objet VariableOrderMarkov.
 *              flag copie de l'objet VariableOrderMarkovData.
 *
 *--------------------------------------------------------------*/

void VariableOrderMarkov::copy(const VariableOrderMarkov &markov , bool data_flag)

{
  register int i , j;


  nb_iterator = 0;

  if ((data_flag) && (markov.markov_data)) {
    markov_data = new VariableOrderMarkovData(*(markov.markov_data) , false);
  }
  else {
    markov_data = NULL;
  }

  memory_type = new int[nb_row];
  for (i = 0;i < nb_row;i++) {
    memory_type[i] = markov.memory_type[i];
  }

  order = new int[nb_row];
  for (i = 0;i < nb_row;i++) {
    order[i] = markov.order[i];
  }

  max_order = markov.max_order;

  state = new int*[nb_row];
  for (i = 0;i < nb_row;i++) {
    state[i] = new int[order[i]];
    for (j = 0;j < order[i];j++) {
      state[i][j] = markov.state[i][j];
    }
  }

  parent = new int[nb_row];
  for (i = 0;i < nb_row;i++) {
    parent[i] = markov.parent[i];
  }

  child = new int*[nb_row];
  for (i = 0;i < nb_row;i++) {
    if (markov.child[i]) {
      child[i] = new int[nb_state];
      for (j = 0;j < nb_state;j++) {
        child[i][j] = markov.child[i][j];
      }
    }
    else {
      child[i] = NULL;
    }
  }

  if (markov.next) {
    next = new int*[nb_row];
    for (i = 0;i < nb_row;i++) {
      if (markov.next[i]) {
        next[i] = new int[nb_state];
        for (j = 0;j < nb_state;j++) {
          next[i][j] = markov.next[i][j];
        }
      }
      else {
        next[i] = NULL;
      }
    }
  }
  else {
    next = NULL;
  }

  if (markov.nb_memory) {
    nb_memory = new int[nb_row];
    for (i = 0;i < nb_row;i++) {
      nb_memory[i] = markov.nb_memory[i];
    }
  }
  else {
    nb_memory = NULL;
  }

  if (markov.previous) {
    previous = new int*[nb_row];
    for (i = 0;i < nb_row;i++) {
      if (markov.previous[i]) {
        previous[i] = new int[nb_memory[i]];
        for (j = 0;j < nb_memory[i];j++) {
          previous[i][j] = markov.previous[i][j];
        }
      }
      else {
        previous[i] = NULL;
      }
    }
  }
  else {
    previous = NULL;
  }

  nb_output_process = markov.nb_output_process;

  nonparametric_process = new NonparametricSequenceProcess*[nb_output_process + 1];

  if (markov.parametric_process) {
    parametric_process = new DiscreteParametricProcess*[nb_output_process + 1];
    parametric_process[0] = NULL;
  }
  else {
    parametric_process = NULL;
  }

  nonparametric_process[0] = new NonparametricSequenceProcess(*(markov.nonparametric_process[0]));

  for (i = 1;i <= nb_output_process;i++) {
    if (markov.nonparametric_process[i]) {
      nonparametric_process[i] = new NonparametricSequenceProcess(*(markov.nonparametric_process[i]));
      if (markov.parametric_process) {
        parametric_process[i] = NULL;
      }
    }
    else {
      nonparametric_process[i] = NULL;
      parametric_process[i] = new DiscreteParametricProcess(*(markov.parametric_process[i]));
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet VariableOrderMarkov.
 *
 *--------------------------------------------------------------*/

void VariableOrderMarkov::remove()

{
  register int i;


  delete markov_data;

  delete [] memory_type;
  delete [] order;

  if (state) {
    for (i = 0;i < nb_row;i++) {
      delete [] state[i];
    }
    delete [] state;
  }

  delete [] parent;

  if (child) {
    for (i = 0;i < nb_row;i++) {
      delete [] child[i];
    }
    delete [] child;
  }

  if (next) {
    for (i = 1;i < nb_row;i++) {
      delete [] next[i];
    }
    delete [] next;
  }

  delete [] nb_memory;

  if (previous) {
    for (i = 1;i < nb_row;i++) {
      delete [] previous[i];
    }
    delete [] previous;
  }

  if (nonparametric_process) {
    for (i = 0;i <= nb_output_process;i++) {
      delete nonparametric_process[i];
    }
    delete [] nonparametric_process;
  }

  if (parametric_process) {
    for (i = 1;i <= nb_output_process;i++) {
      delete parametric_process[i];
    }
    delete [] parametric_process;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe VariableOrderMarkov.
 *
 *--------------------------------------------------------------*/

VariableOrderMarkov::~VariableOrderMarkov()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe VariableOrderMarkov.
 *
 *  argument : reference sur un objet VariableOrderMarkov.
 *
 *--------------------------------------------------------------*/

VariableOrderMarkov& VariableOrderMarkov::operator=(const VariableOrderMarkov &markov)

{
  if (&markov != this) {
    remove();
    Chain::remove();

    Chain::copy(markov);
    copy(markov);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  recherche de la memoire pere.
 *
 *  argument : indice de la memoire.
 *
 *--------------------------------------------------------------*/

void VariableOrderMarkov::find_parent_memory(int index)

{
  register int i;


  for (i = index - 1;i >= 0;i--) {
    if (order[i] == order[index] - 1) {
      parent[index] = i;
      child[i][state[index][order[index] - 1]] = index;
      break;
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Construction des transitions entre memoires.
 *
 *--------------------------------------------------------------*/

void VariableOrderMarkov::build_memory_transition()

{
  if (!next) {
    register int i , j , k;
    int bnb_memory;


#   ifdef DEBUG
    cout << "\n";
#   endif

    next = new int*[nb_row];
    next[0] = NULL;

    for (i = 1;i < nb_row;i++) {
      if ((type == 'o') || (!child[i])) {
        next[i] = new int[nb_state];

#       ifdef DEBUG
        for (j = order[i] - 1;j >= 0;j--) {
          cout << state[i][j] << " ";
        }
#       endif

        bnb_memory = 0;
        for (j = 1;j < nb_row;j++) {
          if (((child[i]) && (child[j]) && (order[j] == order[i] + 1)) ||
              ((!child[j]) && (order[j] <= order[i] + 1))) {
            for (k = 0;k < order[j] - 1;k++) {
              if (state[j][k + 1] != state[i][k]) {
                break;
              }
            }

            if ((order[j] == 1) || (k == order[j] - 1)) {
//              if ((memory_type[i] == NON_TERMINAL) || (transition[i][state[j][0]] > 0.)) {
                next[i][state[j][0]] = j;
/*              }
              else {
                next[i][state[j][0]] = I_DEFAULT;
              } */

#             ifdef DEBUG
              cout << "| ";
              for (k = order[j] - 1;k >= 0;k--) {
                cout << state[j][k] << " ";
              }
#             endif

              bnb_memory++;
              if (bnb_memory == nb_state) {
                break;
              }
            }
          }
        }

#       ifdef DEBUG
        cout << endl;
#       endif

      }

      else {
        next[i] = NULL;
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Extraction des memoires precedentes.
 *
 *--------------------------------------------------------------*/

void VariableOrderMarkov::build_previous_memory()

{
  if ((next) && (!nb_memory) && (!previous)) {
    register int i , j;
    int *buffer;


    nb_memory = new int[nb_row];
    previous = new int*[nb_row];
    nb_memory[0] = 0;
    previous[0] = NULL;

    buffer = new int[nb_row - 1];
    for (i = 1;i < nb_row;i++) {
      nb_memory[i] = 0;

//      if (next[i]) {
      if ((type == 'o') || (!child[i])) {
        for (j = 1;j < nb_row;j++) {
          if ((next[j]) && (next[j][state[i][0]] == i)) {
            buffer[nb_memory[i]] = j;
            nb_memory[i]++;
          }
        }
      }

      if (nb_memory[i] > 0) {
        previous[i] = new int[nb_memory[i]];
        for (j = 0;j < nb_memory[i];j++) {
          previous[i][j] = buffer[j];
        }
      }
      else {
        previous[i] = NULL;
      }
    }

    delete [] buffer;
  }
}


/*--------------------------------------------------------------*
 *
 *  Verification que les memoires terminales est un ensemble libre de suffixe.
 *
 *--------------------------------------------------------------*/

bool VariableOrderMarkov::check_free_suffix() const

{
  bool free_suffix = true;
  register int i , j , k;


  for (i = 1;i < nb_row;i++) {
    if ((!child[i]) && (order[i] >= 2)) {
      for (j = 1;j < nb_row;j++) {
        if ((!child[j]) && (order[j] < order[i])) {
          for (k = 0;k < order[j];k++) {
            if (state[i][k] != state[j][k]) {
              break;
            }
          }

          if (k == order[j]) {
            free_suffix = false;
            break;
          }
        }
      }

      if (j < nb_row) {
        break;
      }
    }
  }

  return free_suffix;
}


/*--------------------------------------------------------------*
 *
 *  Construction de la matrice des transitions possibles entre etats.
 *
 *--------------------------------------------------------------*/

bool** VariableOrderMarkov::logic_transition_computation() const

{
  bool **logic_transition;
  register int i , j;
  double **order1_transition;


  logic_transition = new bool*[nb_state];
  for (i = 0;i < nb_state;i++) {
    logic_transition[i] = new bool[nb_state];
    logic_transition[i][i] = false;
  }

  order1_transition = new double*[nb_state];
  for (i = 0;i < nb_state;i++) {
    order1_transition[i] = new double[nb_state];
    for (j = 0;j < nb_state;j++) {
      order1_transition[i][j] = 0.;
    }
  }

  for (i = 1;i < nb_row;i++) {
    if (memory_type[i] == TERMINAL) {
      for (j = 0;j < nb_state;j++) {
        order1_transition[state[i][0]][j] += transition[i][j];
      }
    }
  }

  for (i = 0;i < nb_state;i++) {
    for (j = 0;j < nb_state;j++) {
      if (j != i) {
        logic_transition[i][j] = (order1_transition[i][j] == 0. ? false : true);
      }
    }
  }

  for (i = 0;i < nb_state;i++) {
    delete [] order1_transition[i];
  }
  delete [] order1_transition;

  return logic_transition;
}


/*--------------------------------------------------------------*
 *
 *  Extraction des classes d'une chaine de Markov d'ordre variable
 *  a partir de l'accessibilite des etats.
 *
 *--------------------------------------------------------------*/

void VariableOrderMarkov::component_computation()

{
  bool **logic_transition;
  register int i;


  logic_transition = logic_transition_computation();

  Chain::component_computation(logic_transition);

  for (i = 0;i < nb_state;i++) {
    delete [] logic_transition[i];
  }
  delete [] logic_transition;
}


/*--------------------------------------------------------------*
 *
 *  Construction des memoires non-terminales.
 *
 *--------------------------------------------------------------*/

void VariableOrderMarkov::build_non_terminal()

{
  register int i , j , k , m;
  int nb_non_terminal , nb_terminal , bnb_memory;


  nb_non_terminal = (nb_row - 1) / nb_state;
  nb_terminal = (nb_row - 1) * (nb_state - 1) / nb_state + 1;
//  nb_terminal = nb_non_terminal * (nb_state - 1) + 1;

  parent[0] = -1;

  i = 0;
  for (j = nb_non_terminal;j < nb_non_terminal + nb_terminal;j++) {
    for (k = order[j] - 1;k >= 0;k--) {
      if (state[j][k] != 0) {
        break;
      }
    }
    bnb_memory = order[j] - 1 - k;

    // insertion des memoires non-terminales

    for (k = order[j] - bnb_memory;k < order[j];k++) {
      memory_type[i] = NON_TERMINAL;

      for (m = 0;m < nb_state;m++) {
        transition[i][m] = 0.;
      }

      order[i] = k;
      for (m = 0;m < order[i];m++) {
        state[i][m] = state[j][m];
      }

      child[i] = new int[nb_state];
      i++;
    }

    // copie de la memoire terminale

    memory_type[i] = TERMINAL;
    if (i < j) {
      for (k = 0;k < nb_state;k++) {
        transition[i][k] = transition[j][k];
      }
      order[i] = order[j];
      for (k = 0;k < order[i];k++) {
        state[i][k] = state[j][k];
      }
    }

    // relations parent-enfant

    for (k = 0;k < bnb_memory;k++) {
      parent[i - k] = i - k - 1;
      child[i - k - 1][state[i - k][order[i - k] - 1]] = i - k;
    }

    find_parent_memory(i - bnb_memory);
    i++;
  }

# ifdef DEBUG
  {
    cout << "\n";
    cout << "Suffix free? " << (check_free_suffix() ? "True" : "False") << endl;

    for (i = 0;i < nb_row;i++) {
      for (j = max_order - 1;j >= order[i];j--) {
        cout << "  ";
      }
      for (j = order[i] - 1;j >= 0;j--) {
        cout << state[i][j] << " ";
      }

      switch (memory_type[i]) {
      case NON_TERMINAL :
        cout << "  " << SEQ_label[SEQL_NON_TERMINAL];
        break;
      case TERMINAL :
        cout << "      " << SEQ_label[SEQL_TERMINAL];
        break;
      }

      cout << " | " << parent[i];

      if (child[i]) {
        cout << " |";
        for (j = 0;j < nb_state;j++) {
          cout << " " << child[i][j];
        }
      }

      else {
        cout << "  ";
        for (j = 0;j < nb_state;j++) {
          cout << "  ";
        }
      }

      cout << endl;
    }
  }
# endif

}


/*--------------------------------------------------------------*
 *
 *  Extraction d'une loi.
 *
 *  arguments : reference sur un objet StatError, type de loi,
 *              variable, etat ou observation.
 *
 *--------------------------------------------------------------*/

DiscreteParametricModel* VariableOrderMarkov::extract(StatError &error , int type ,
                                                      int variable , int value) const

{
  bool status = true;
  int hvariable;
  Distribution *pdist;
  DiscreteParametric *pparam;
  DiscreteParametricModel *dist;
  FrequencyDistribution *phisto;


  dist = NULL;
  error.init();

  pdist = NULL;
  pparam = NULL;

  if (type == OBSERVATION) {
    if ((variable < 1) || (variable > nb_output_process)) {
      status = false;
      error.update(SEQ_error[SEQR_OUTPUT_PROCESS_INDEX]);
    }

    else {
      if ((value < 0) || (value >= nb_state)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_STATE] << " " << value << " "
                      << SEQ_error[SEQR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }

      else {
        if (nonparametric_process[variable]) {
          pdist = nonparametric_process[variable]->observation[value];
        }
        else {
          pparam = parametric_process[variable]->observation[value];
        }
      }
    }
  }

  else {
    if ((variable < 0) || (variable > nb_output_process)) {
      status = false;
      error.update(SEQ_error[SEQR_OUTPUT_PROCESS_INDEX]);
    }

    else {
      if (!nonparametric_process[variable]) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << variable << ": " 
                      << SEQ_error[SEQR_CHARACTERISTICS_NOT_COMPUTED];
        error.update((error_message.str()).c_str());
      }

      else if ((value < 0) || (value >= nonparametric_process[variable]->nb_value)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[variable == 0 ? STATL_STATE : STATL_OUTPUT] << " "
                      << value << " " << SEQ_error[SEQR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }

      if (status) {
        switch (type) {
        case FIRST_OCCURRENCE :
          pdist = nonparametric_process[variable]->first_occurrence[value];
          break;
        case RECURRENCE_TIME :
          pdist = nonparametric_process[variable]->recurrence_time[value];
          break;
        case SOJOURN_TIME :
          pparam = nonparametric_process[variable]->sojourn_time[value];
          break;
        case NB_RUN :
          pdist = nonparametric_process[variable]->nb_run[value];
          break;
        case NB_OCCURRENCE :
          pdist = nonparametric_process[variable]->nb_occurrence[value];
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

    if (markov_data) {
      switch (markov_data->type[0]) {
      case INT_VALUE :
        hvariable = variable - 1;
        break;
      case STATE :
        hvariable = variable;
        break;
      }

      if (hvariable >= 0) {
        switch (type) {

        case OBSERVATION : {
          if ((markov_data->observation) && (markov_data->observation[hvariable])) {
            phisto = markov_data->observation[hvariable][value];
          }
          break;
        }

        case FIRST_OCCURRENCE : {
          phisto = markov_data->characteristics[hvariable]->first_occurrence[value];
          break;
        }

        case RECURRENCE_TIME : {
          if (markov_data->characteristics[hvariable]->recurrence_time[value]->nb_element > 0) {
            phisto = markov_data->characteristics[hvariable]->recurrence_time[value];
          }
          break;
        }

        case SOJOURN_TIME : {
          if (markov_data->characteristics[hvariable]->sojourn_time[value]->nb_element > 0) {
            phisto = markov_data->characteristics[hvariable]->sojourn_time[value];
          }
          break;
        }

        case NB_RUN : {
          phisto = markov_data->characteristics[hvariable]->nb_run[value];
          break;
        }

        case NB_OCCURRENCE : {
          phisto = markov_data->characteristics[hvariable]->nb_occurrence[value];
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
 *  Extraction de la partie "donnees" d'un objet VariableOrderMarkov.
 *
 *  argument : reference sur un objet StatError.
 *
 *--------------------------------------------------------------*/

VariableOrderMarkovData* VariableOrderMarkov::extract_data(StatError &error) const

{
  bool status = true;
  VariableOrderMarkovData *seq;


  seq = NULL;
  error.init();

  if (!markov_data) {
    status = false;
    error.update(STAT_error[STATR_NO_DATA]);
  }
  else if (nb_output_process + 1 != markov_data->nb_variable) {
    status = false;
    error.update(SEQ_error[SEQR_STATE_SEQUENCES]);
  }

  if (status) {
    seq = new VariableOrderMarkovData(*markov_data);
    seq->markov = new VariableOrderMarkov(*this , false);
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Application d'un seuil sur les parametres d'une chaine de Markov d'ordre variable.
 *
 *  argument : probabilite minimum.
 *
 *--------------------------------------------------------------*/

void VariableOrderMarkov::threshold_application(double min_probability)

{
  bool stop;
  register int i , j;
  int nb_correction;
  double norm;


  if (min_probability > THRESHOLDING_FACTOR / (double)nb_state) {
    min_probability = THRESHOLDING_FACTOR / (double)nb_state;
  }

  if (type == 'o') {
    do {
      stop = true;
      nb_correction = 0;
      norm = 0.;

      for (i = 0;i < nb_state;i++) {
        if (initial[i] <= min_probability) {
          nb_correction++;
          initial[i] = min_probability;
        }
        else {
          norm += initial[i];
        }
      }

      if (nb_correction > 0) {
        for (i = 0;i < nb_state;i++) {
          if (initial[i] > min_probability) {
            initial[i] *= (1. - nb_correction * min_probability) / norm;
            if (initial[i] < min_probability) {
              stop = false;
            }
          }
        }
      }
    }
    while (!stop);
  }

  for (i = 1;i < nb_row;i++) {
    if ((memory_type[i] == TERMINAL) || ((type == 'o') &&
         (memory_type[i] == NON_TERMINAL))) {
      do {
        stop = true;
        nb_correction = 0;
        norm = 0.;

        for (j = 0;j < nb_state;j++) {
          if (transition[i][j] <= min_probability) {
            nb_correction++;
            transition[i][j] = min_probability;
          }
          else {
            norm += transition[i][j];
          }
        }

        if (nb_correction > 0) {
          for (j = 0;j < nb_state;j++) {
            if (transition[i][j] > min_probability) {
              transition[i][j] *= (1. - nb_correction * min_probability) / norm;
              if (transition[i][j] < min_probability) {
                stop = false;
              }
            }
          }
        }
      }
      while (!stop);
    }

    else if (memory_type[i] == COMPLETION) {
      for (j = 0;j < nb_state;j++) {
        transition[i][j] = transition[parent[i]][j];
      }
    }
  }

  if (accessibility) {
    for (i = 0;i < nb_state;i++) {
      delete [] accessibility[i];
    }
    delete [] accessibility;
  }
  accessibility = NULL;

  delete [] component_nb_state;

  if (component) {
    for (i = 0;i < nb_component;i++) {
      delete [] component[i];
    }
    delete [] component;
  }
  nb_component = 0;

  component_computation();

  if (next) {
    for (i = 1;i < nb_row;i++) {
      delete [] next[i];
    }
    delete [] next;
  }
  next = NULL;

  delete [] nb_memory;
  nb_memory = NULL;

  if (previous) {
    for (i = 1;i < nb_row;i++) {
      delete [] previous[i];
    }
    delete [] previous;
  }
  previous = NULL;

  build_memory_transition();
  build_previous_memory();
}


/*--------------------------------------------------------------*
 *
 *  Application d'un seuil sur les parametres d'une chaine de Markov d'ordre variable.
 *
 *  argument : probabilite minimum.
 *
 *--------------------------------------------------------------*/

VariableOrderMarkov* VariableOrderMarkov::thresholding(double min_probability) const

{
  register int i;
  VariableOrderMarkov *markov;


  markov = new VariableOrderMarkov(*this , false);
  markov->threshold_application(min_probability);

  for (i = 1;i <= markov->nb_output_process;i++) {
    if (markov->nonparametric_process[i]) {
      markov->nonparametric_process[i]->thresholding(min_probability);
    }
  }

  return markov;
}


/*--------------------------------------------------------------*
 *
 *  Analyse du format d'un objet VariableOrderMarkov.
 *
 *  arguments : reference sur un objet StatError, stream,
 *              reference sur l'indice de la ligne lue, type du processus
 *              ('o' : ordinaire, 'e' : en equilibre).
 *
 *--------------------------------------------------------------*/

VariableOrderMarkov* variable_order_markov_parsing(StatError &error , ifstream &in_file ,
                                                   int &line , char type)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  streampos transition_line;
  bool status = true , lstatus , increase , **logic_transition;
  register int i , j;
  int read_line , tline , nb_state = 0 , order , previous_order , max_order = 0 , buff ,
      nb_terminal , nb_non_terminal , memory , state[ORDER] , previous_state[ORDER];
  long value;
  double proba , cumul , *initial;
  VariableOrderMarkov *markov;


  markov = NULL;

  // analyse ligne definissant le nombre d'etats

  while (buffer.readLine(in_file , false)) {
    line++;

#   ifdef DEBUG
    cout << line << "  " << buffer << endl;
#   endif

    position = buffer.first('#');
    if (position != RW_NPOS) {
      buffer.remove(position);
    }
    i = 0;

    RWCTokenizer next(buffer);

    while (!((token = next()).isNull())) {
      switch (i) {

      // test valeur nombre d'etats

      case 0 : {
        lstatus = locale.stringToNum(token , &value);
        if (lstatus) {
          if ((value < 2) || (value > NB_STATE)) {
            lstatus = false;
          }
          else {
            nb_state = value;
          }
        }

        if (!lstatus) {
          status = false;
          error.update(STAT_parsing[STATP_NB_STATE] , line , i + 1);
        }
        break;
      }

      // test mot cle STATES

      case 1 : {
        if (token != STAT_word[STATW_STATES]) {
          status = false;
          error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_STATES] , line , i + 1);
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
      break;
    }
  }

  if (nb_state == 0) {
    status = false;
    error.update(STAT_parsing[STATP_FORMAT] , line);
  }

  if (status) {
    initial = new double[nb_state];

    // 1ere passe : recherche du nombre de memoires terminales et de l'ordre maximum,
    // analyse probabilites initiales, probabilites de transition et memoires

    read_line = 0;
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

      if ((read_line == 0) || ((type == 'o') && (read_line == 2))) {
        while (!((token = next()).isNull())) {

          // test mot cle INITIAL_PROBABILITIES / TRANSITION_PROBABILITIES

          if (i == 0) {
            if ((type == 'o') && (read_line == 0)) {
              if (token != STAT_word[STATW_INITIAL_PROBABILITIES]) {
                status = false;
                error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_INITIAL_PROBABILITIES] , line);
              }
            }

            else {
              if (token != STAT_word[STATW_TRANSITION_PROBABILITIES]) {
                status = false;
                error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_TRANSITION_PROBABILITIES] , line);
              }
            }
          }

          i++;
        }

        if (i > 0) {
          if (((type == 'o') && (read_line == 2)) || ((type == 'e') && (read_line == 0))) {
            transition_line = in_file.tellg();
            tline = line;
          }

          if (i != 1) {
            status = false;
            error.update(STAT_parsing[STATP_FORMAT] , line);
          }
        }
      }

      else {
        cumul = 0.;

        while (!((token = next()).isNull())) {
          if (i < nb_state) {
            lstatus = locale.stringToNum(token , &proba);
            if (lstatus) {
              if ((proba < 0.) || (proba > 1. - cumul + DOUBLE_ERROR)) {
                lstatus = false;
              }

              else {
                cumul += proba;
                if ((type == 'o') && (read_line == 1)) {
                  initial[i] = proba;
                }
              }
            }

            if (!lstatus) {
              status = false;
              if ((type == 'o') && (read_line == 1)) {
                error.update(STAT_parsing[STATP_INITIAL_PROBABILITY] , line , i + 1);
              }
              else {
                error.update(STAT_parsing[STATP_TRANSITION_PROBABILITY] , line , i + 1);
              }
            }
          }

          else if ((type == 'e') || (read_line >= 3)) {
            lstatus = locale.stringToNum(token , &value);
            if (lstatus) {
              if ((value < 0) || (value >= nb_state)) {
                lstatus = false;
              }
              else if (i - nb_state < ORDER) {
                state[i - nb_state] = value;
              }
            }

            if (!lstatus) {
              status = false;
              error.update(SEQ_parsing[SEQP_STATE] , line , i + 1);
            }
          }

          i++;
        }

        if (i > 0) {
          if ((type == 'o') && (read_line == 1) && (i != nb_state)) {
            status = false;
            error.update(STAT_parsing[STATP_FORMAT] , line);
          }

          if (cumul < 1. - DOUBLE_ERROR) {
            status = false;
            error.update(STAT_parsing[STATP_PROBABILITY_SUM] , line);
          }

          if (((type == 'o') && (read_line >= 3)) || ((type == 'e') && (read_line >= 1))) {
            if (i <= nb_state) {
              status = false;
              error.update(STAT_parsing[STATP_FORMAT] , line);
            }

            else {
              order = i - nb_state;

              if (order > ORDER) {
                status = false;
                error.update(STAT_parsing[STATP_ORDER] , line);
              }

              else if (order > max_order) {
                max_order = order;
              }

              if (status) {
                for (j = 0;j < order / 2;j++) {
                  buff = state[j];
                  state[j] = state[order - j - 1];
                  state[order - j - 1] = buff;
                }

                if (read_line - (type == 'o' ? 3 : 1) == 0) {
                  for (j = 0;j < order;j++) {
                    if (state[j] != 0) {
                      status = false;
                      error.update(SEQ_parsing[SEQP_STATE] , line , nb_state + order - j);
                    }
                  }
                }

                else {

                  // verification succession des memoires

                  increase = true;
                  for (j = MIN(previous_order , order) - 1;j >= 0;j--) {
                    if (increase) {
                      if (previous_state[j] < nb_state - 1) {
                        if (state[j] != previous_state[j] + 1) {
                          status = false;
                          error.update(SEQ_parsing[SEQP_STATE] , line , nb_state + order - j);
                        }
                        increase = false;
                      }

                      else {
                        if (state[j] != 0) {
                          status = false;
                          error.update(SEQ_parsing[SEQP_STATE] , line , nb_state + order - j);
                        }
                      }
                    }

                    else {
                      if (state[j] != previous_state[j]) {
                        status = false;
                        error.update(SEQ_parsing[SEQP_STATE] , line , nb_state + order - j);
                      }
                    }
                  }

                  // cas augmentation taille de la memoire ou taille de la memoire conservee

                  if (order >= previous_order) {
                    for (j = order - 1;j >= previous_order;j--) {
                      if (state[j] != 0) {
                        status = false;
                        error.update(SEQ_parsing[SEQP_STATE] , line , nb_state + order - j);
                      }
                    }
                  }

                  // cas reduction taille de la memoire

                  else {
                    for (j = order;j < previous_order;j++) {
                      if (previous_state[j] != nb_state - 1) {
                        status = false;
                        error.update(SEQ_parsing[SEQP_STATE] , line);
                      }
                    }
                  }
                }

                // recherche derniere memoire

                for (j = 0;j < order;j++) {
                  if (state[j] != nb_state - 1) {
                    break;
                  }
                }

                if (j == order) {
                  read_line++;
                  break;
                }

                else {
                  for (j = 0;j < order;j++) {
                    previous_state[j] = state[j];
                  }
                  previous_order = order;
                }
              }
            }
          }
        }
      }

      if (i > 0) {
        read_line++;
      }
    }

    // verification nombre de memoires

    nb_terminal = read_line - (type == 'o' ? 3 : 1);
    if ((nb_state > 2) && (nb_terminal % (nb_state - 1) != 1)) {
      status = false;
      error.update(SEQ_parsing[SEQP_NB_MEMORY] , line);
    }

    if (status) {
      in_file.clear();
      in_file.seekg(transition_line);

      nb_non_terminal = (nb_terminal - 1) / (nb_state - 1);
      markov = new VariableOrderMarkov(type , nb_state , nb_non_terminal + nb_terminal , max_order);

      if (type == 'o') {
        for (i = 0;i < nb_state;i++) {
          markov->initial[i] = initial[i];
        }
      }

      // 2eme passe : lecture des probabilites de transition et des memoires

      memory = nb_non_terminal;
      line = tline;
      while (buffer.readLine(in_file , false)) {
        line++;

#       ifdef DEBUG
        cout << line << "  " << buffer << endl;
#       endif

        position = buffer.first('#');
        if (position != RW_NPOS) {
          buffer.remove(position);
        }
        i = 0;

        RWCTokenizer next(buffer);

        while (!((token = next()).isNull())) {
          if (i < nb_state) {
            locale.stringToNum(token , &proba);
            markov->transition[memory][i] = proba;
          }

          else {
            locale.stringToNum(token , &value);
            state[i - nb_state] = value;
          }

          i++;
        }

        if (i > 0) {
          markov->order[memory] = i - nb_state;
          for (j = 0;j < markov->order[memory];j++) {
            markov->state[memory][j] = state[markov->order[memory] - j - 1];
          }

          // recherche derniere memoire

          for (j = 0;j < markov->order[memory];j++) {
            if (markov->state[memory][j] != nb_state - 1) {
              break;
            }
          }

          if (j == order) {
            break;
          }
          else {
            memory++;
          }
        }
      }

      markov->build_non_terminal();

      // test etats atteignables

      logic_transition = markov->logic_transition_computation();
      status = markov->connex_component_research(error , logic_transition);

      if (status) {

        // test irreductibilite dans le cas en equilibre

        markov->Chain::component_computation(logic_transition);
        if ((type == 'e') && (markov->nb_component > 1)) {
          status = false;
          error.correction_update(STAT_parsing[STATP_CHAIN_STRUCTURE] , STAT_parsing[STATP_IRREDUCIBLE]);
        }
      }

      for (i = 0;i < nb_state;i++) {
        delete [] logic_transition[i];
      }
      delete [] logic_transition;

      if (!status) {
        delete markov;
        markov = NULL;
      }
      else {
        markov->max_order_computation();
      }
    }

    delete [] initial;
  }

  return markov;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet VariableOrderMarkov a partir d'un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              longueur des sequences.
 *
 *--------------------------------------------------------------*/

VariableOrderMarkov* variable_order_markov_ascii_read(StatError &error ,
                                                      const char *path , int length)

{
  RWCString buffer , token;
  size_t position;
  char type = 'v';
  bool status;
  register int i;
  int line;
  const VariableOrderMarkov *imarkov;
  const NonparametricProcess *observation;
  VariableOrderMarkov *markov;
  ifstream in_file(path);


  markov = NULL;
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

        // test mot cle (EQUILIBRIUM) MARKOV_CHAIN

        if (i == 0) {
          if (token == SEQ_word[SEQW_MARKOV_CHAIN]) {
            type = 'o';
          }
          else if (token == SEQ_word[SEQW_EQUILIBRIUM_MARKOV_CHAIN]) {
            type = 'e';
          }
          else {
            status = false;
            ostringstream correction_message;
            correction_message << SEQ_word[SEQW_MARKOV_CHAIN] << " or "
                               << SEQ_word[SEQW_EQUILIBRIUM_MARKOV_CHAIN];
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

      // analyse du format et lecture de la chaine de Markov d'ordre variable

      imarkov = variable_order_markov_parsing(error , in_file , line , type);

      if (imarkov) {

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

            observation = observation_parsing(error , in_file , line ,
                                              ((Chain*)imarkov)->nb_state , false);
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
          markov = new VariableOrderMarkov(imarkov , observation , length);

#         ifdef DEBUG
          imarkov->ascii_memory_tree_print(cout);
          markov->ascii_write(cout);
#         endif

        }

        delete imarkov;
        delete observation;
      }
    }
  }

  return markov;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet VariableOrderMarkov.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& VariableOrderMarkov::line_write(ostream &os) const

{
  os << nb_state << " " << STAT_word[STATW_STATES] << "   "
     << SEQ_label[SEQL_MAX_ORDER] << " " << max_order;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Sortie ASCII de l'arborescence des memoires.
 *
 *  arguments : stream, flag fichier.
 *
 *--------------------------------------------------------------*/

ostream& VariableOrderMarkov::ascii_memory_tree_print(ostream &os , bool file_flag) const

{
  register int i , j , k;
  int bnb_memory , width = column_width(nb_state);
  long old_adjust;


  old_adjust = os.setf(ios::left , ios::adjustfield);

  os << "\n";
  if (file_flag) {
    os << "# ";
  }
  os << SEQ_label[SEQL_MEMORY_TREE] << endl;

  os << "\n";
  for (i = 0;i < nb_row;i++) {
    if (!child[i]) {
      if (file_flag) {
        os << "# ";
      }

      for (j = order[i] - 1;j >= 1;j--) {
        if (state[i][j] != 0) {
          break;
        }
      }
      bnb_memory = order[i] - j;

      for (j = 0;j < order[i] - bnb_memory;j++) {
        if (state[i][j] < nb_state - 1) {
          os << "|";
        }
        else {
          os << " ";
        }

        for (k = 0;k < (j + 1) * (width + 1) + 1;k++) {
          os << " ";
        }
      }

      os << "|";
      for (j = order[i] - bnb_memory;j < order[i];j++) {
        for (k = 0;k < 3;k++) {
          os << "_";
        }
        for (k = j;k >= 0;k--) {
          os << setw(width) << state[i][k];
          if (k > 0) {
            os << " ";
          }
        }
      }

      os << endl;
    }
  }

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

# ifdef MESSAGE
  os << "\n";
  for (i = 0;i < nb_row;i++) {
    if (child[i]) {
      if (file_flag) {
        os << "# ";
      }

      for (j = max_order - 1;j >= order[i];j--) {
        os << "  ";
      }
      for (j = order[i] - 1;j >= 0;j--) {
        os << state[i][j] << " ";
      }

      for (j = 0;j < nb_state;j++) {
        os << "  ";
        for (k = max_order - 1;k >= order[child[i][j]];k--) {
          os << "  ";
        }
        for (k = order[child[i][j]] - 1;k >= 0;k--) {
          os << state[child[i][j]][k] << " ";
        }
      }

      if (memory_type[i] == COMPLETION) {
        os << "    " << SEQ_label[SEQL_COMPLETION];
      }

      os << endl;
    }
  }
# endif

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Sortie ASCII de l'arborescence correspondant a la succession des memoires
 *  de taille croissante dans le graphe des transitions.
 *
 *  arguments : stream, flag fichier.
 *
 *--------------------------------------------------------------*/

ostream& VariableOrderMarkov::ascii_transition_tree_print(ostream &os , bool file_flag) const

{
  register int i , j , k;
  int min_order , nb_root , memory , width = column_width(nb_state) , *nb_next_memory ,
      *root , *nb_leaf_memory , *nb_drawn_next_memory , *nb_drawn_leaf_memory;
  long old_adjust;


  // calcul du nombre de memoires suivantes de taille superieure

  nb_next_memory = new int[nb_row];
  root = new int[nb_row];
  nb_leaf_memory = new int[nb_row];
  nb_drawn_next_memory = new int[nb_row];
  nb_drawn_leaf_memory = new int[nb_row];

  min_order = max_order;

  for (i = 1;i < nb_row;i++) {
    nb_next_memory[i] = 0;

    if ((type == 'o') || (!child[i])) {
      for (j = 0;j < nb_state;j++) {
        if (order[next[i][j]] == order[i] + 1) {
          nb_next_memory[i]++;
        }
      }

      if (nb_next_memory[i] == 0) {
        nb_leaf_memory[i] = 1;
      }
      else {
        nb_leaf_memory[i] = 0;
      }

      if (order[i] < min_order) {
        min_order = order[i];
      }
    }
  }

  // calcul des racines

  nb_root = 0;
  for (i = 1;i < nb_row;i++) {
    if ((type == 'o') || (!child[i])) {
      for (j = 0;j < nb_memory[i];j++) {
        if (order[previous[i][j]] == order[i] - 1) {
          break;
        }
      }

      if (j == nb_memory[i]) {
        root[nb_root++] = i;
      }
    }
  }

  // calcul du nombre de memoires feuilles

  for (i = max_order;i >= min_order + 1;i--) {
    for (j = 1;j < nb_row;j++) {
      if ((order[j] == i) && (nb_leaf_memory[j] > 0) &&
          ((type == 'o') || (!child[j]))) {
        for (k = 0;k < nb_memory[j];k++) {
          if ((order[previous[j][k]] == order[j] - 1) &&
              ((type == 'o') || (!child[previous[j][k]]))) {
            nb_leaf_memory[previous[j][k]] += nb_leaf_memory[j];
            break;
          }
        }
      }
    }
  }

# ifdef DEBUG
  cout << "\n";
  for (i = 1;i < nb_row;i++) {
    if ((type == 'o') || (!child[i])) {
      for (j = max_order - 1;j >= order[i];j--) {
        cout << "  ";
      }
      for (j = order[i] - 1;j >= 0;j--) {
        cout << state[i][j] << " ";
      }

      cout << "   "  << nb_next_memory[i] << " | " << nb_leaf_memory[i];
      for (j = 0;j < nb_root;j++) {
        if (root[j] == i) {
          cout << "   root";
          break;
        }
      }
      cout << endl;
    }
  }
# endif

  for (i = 1;i < nb_row;i++) {
    if ((type == 'o') || (!child[i])) {
      nb_drawn_next_memory[i] = nb_next_memory[i];
      nb_drawn_leaf_memory[i] = nb_leaf_memory[i];
    }
  }

  old_adjust = os.setf(ios::left , ios::adjustfield);

  os << "\n";
  if (file_flag) {
    os << "# ";
  }
  os << SEQ_label[SEQL_TRANSITION_TREE] << endl;

  for (i = 0;i < nb_root;i++) {
    os << "\n";
    for (j = 0;j < nb_leaf_memory[root[i]];j++) {
      if (file_flag) {
        os << "# ";
      }

      memory = root[i];
      for (;;) {
        if (!child[memory]) {
          if (nb_drawn_leaf_memory[memory] == nb_leaf_memory[memory]) {
            if (memory != root[i]) {
              for (k = 0;k < 3;k++) {
                os << "_";
              }
            }

            for (k = order[memory] - 1;k >= 0;k--) {
              os << setw(width) << state[memory][k];
              if (k > 0) {
                os << " ";
              }
            }
          }

          else {
            if (memory != root[i]) {
              for (k = 0;k < 3;k++) {
                os << " ";
              }
            }

            for (k = 0;k < order[memory] * (width + 1) - 2;k++) {
              os << " ";
            }

            if (nb_drawn_next_memory[memory] > 0) {
              os << "|";
            }
            else {
              os << " ";
            }
          }
        }

        else {
          if (memory != root[i]) {
            for (k = 0;k < 3;k++) {
              os << " ";
            }
          }

          for (k = 0;k < order[memory] * (width + 1) - 1;k++) {
            os << " ";
          }
        }

        nb_drawn_leaf_memory[memory]--;

        if (nb_next_memory[memory] == 0) {
          os << endl;
          break;
        }

        for (k = 0;k < nb_state;k++) {
          if ((order[next[memory][k]] == order[memory] + 1) &&
              (nb_drawn_leaf_memory[next[memory][k]] > 0)) {
            if (nb_drawn_leaf_memory[next[memory][k]] == nb_leaf_memory[next[memory][k]]) {
              nb_drawn_next_memory[memory]--;
            }
            memory = next[memory][k];
            break;
          }
        }
      }
    }
  }

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  delete [] nb_next_memory;
  delete [] root;
  delete [] nb_leaf_memory;
  delete [] nb_drawn_next_memory;
  delete [] nb_drawn_leaf_memory;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Sortie ASCII des parametres de la chaine de Markov d'ordre variable.
 *
 *  arguments : stream, flag fichier.
 *
 *--------------------------------------------------------------*/

ostream& VariableOrderMarkov::ascii_print(ostream &os , bool file_flag) const

{
  register int i , j , k;
  int buff , width;
  double *stationary_probability;
  long old_adjust;


  old_adjust = os.setf(ios::left , ios::adjustfield);

  os << "\n" << nb_state << " " << STAT_word[STATW_STATES] << endl;

  // calcul des largeurs des colonnes

  width = column_width((type == 'o' ? nb_state : nb_row) , initial);

  for (i = 1;i < nb_row;i++) {
    buff = column_width(nb_state , transition[i]);
    if (buff > width) {
      width = buff;
    }
  }
  width += ASCII_SPACE;

  os << "\n";

  switch (type) {

  case 'o' : {
    os << STAT_word[STATW_INITIAL_PROBABILITIES] << endl;
    for (i = 0;i < nb_state;i++) {
      os << setw(width) << initial[i];
    }
    os << endl;
    break;
  }

  case 'e' : {

    // extraction des probabilites stationnairees correspondant aux memoires completees

    stationary_probability = new double[nb_row];

    for (i = 1;i < nb_row;i++) {
      stationary_probability[i] = initial[i];
    }
    for (i = nb_row - 1;i >= 1;i--) {
      if (memory_type[i] == COMPLETION) {
        stationary_probability[parent[i]] += stationary_probability[i];
      }
    }

    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_STATIONARY_PROBABILITIES] << endl;

    for (i = 1;i < nb_row;i++) {
      if (memory_type[i] == TERMINAL) {
        if (file_flag) {
          os << "# ";
        }
        os << setw(width) << stationary_probability[i];

        os << " ";
        for (j = max_order - 1;j >= order[i];j--) {
          os << "  ";
        }
        for (j = order[i] - 1;j >= 0;j--) {
          os << state[i][j] << " ";
        }
        os << endl;
      }
    }

    delete [] stationary_probability;
    break;
  }
  }

  os << "\n" << STAT_word[STATW_TRANSITION_PROBABILITIES] << "      ";
  if (file_flag) {
    os << "# ";
  }
  os << STAT_label[STATL_MEMORY] << endl;

  for (i = 1;i < nb_row;i++) {
    if ((memory_type[i] != TERMINAL) && (file_flag)) {
      os << "# ";
    }

    for (j = 0;j < nb_state;j++) {
      os << setw(width) << transition[i][j];
    }

    os << " ";
    for (j = max_order - 1;j >= order[i];j--) {
      os << "  ";
    }
    for (j = order[i] - 1;j >= 0;j--) {
      os << state[i][j] << " ";
    }

    if ((memory_type[i] == TERMINAL) && (file_flag)) {
      os << "# ";
    }

    switch (memory_type[i]) {

    case NON_TERMINAL : {
      os << "  " << SEQ_label[SEQL_NON_TERMINAL];
      break;
    }

    case TERMINAL : {
      os << "      " << SEQ_label[SEQL_TERMINAL];
      if (child[i]) {
        os << " (" << SEQ_label[SEQL_COMPLETED] << ")";
      }
      break;
    }

    case COMPLETION : {
      os << "    " << SEQ_label[SEQL_COMPLETION] << " ("
         << (child[i] ? SEQ_label[SEQL_NON_TERMINAL] : SEQ_label[SEQL_TERMINAL]) << ")";
      break;
    }
    }

    os << endl;
  }

  if (nb_component > 0) {
    for (i = 0;i < nb_component;i++) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }

      switch (state_type[component[i][0]]) {
      case 't' :
        os << STAT_label[STATL_TRANSIENT] << " ";
        break;
      default :
        os << STAT_label[STATL_RECURRENT] << " ";
        break;
      }
      os << STAT_label[STATL_CLASS] << ": " << STAT_label[component_nb_state[i] == 1 ? STATL_STATE : STATL_STATES];

      for (j = 0;j < component_nb_state[i];j++) {
        os << " " << component[i][j];
      }

      if (state_type[component[i][0]] == 'a') {
        os << " (" << STAT_label[STATL_ABSORBING] << " " << STAT_label[STATL_STATE] << ")";
      }
    }
    os << endl;
  }

  ascii_memory_tree_print(os , file_flag);
  ascii_transition_tree_print(os , file_flag);

  os << "\n";
  if (file_flag) {
    os << "# ";
  }
  os << SEQ_label[SEQL_MEMORY_TRANSITION_MATRIX] << endl;

  os << "\n";
  for (i = 1;i < nb_row;i++) {
    if ((type == 'o') || (!child[i])) {
      if (file_flag) {
        os << "# ";
      }

      for (j = max_order - 1;j >= order[i];j--) {
        os << "  ";
      }
      for (j = order[i] - 1;j >= 0;j--) {
        os << state[i][j] << " ";
      }

      for (j = 0;j < nb_state;j++) {
        os << "  ";
        for (k = max_order - 1;k >= order[next[i][j]];k--) {
          os << "  ";
        }
        for (k = order[next[i][j]] - 1;k >= 0;k--) {
          os << state[next[i][j]][k] << " ";
        }
      }

      switch (memory_type[i]) {

      case NON_TERMINAL : {
        os << "  " << SEQ_label[SEQL_NON_TERMINAL];
        break;
      }

      case TERMINAL : {
        os << "      " << SEQ_label[SEQL_TERMINAL];
        if (child[i]) {
          os << " (" << SEQ_label[SEQL_COMPLETED] << ")";
        }
        break;
      }

      case COMPLETION : {
        os << "    " << SEQ_label[SEQL_COMPLETION] << " ("
           << (child[i] ? SEQ_label[SEQL_NON_TERMINAL] : SEQ_label[SEQL_TERMINAL]) << ")";
        break;
      }
      }

      os << endl;
    }
  }

# ifdef DEBUG
  if (previous) {
    os << "\n";
    for (i = 1;i < nb_row;i++) {
      if ((type == 'o') || (!child[i])) {
        for (j = max_order - 1;j >= order[i];j--) {
          os << "  ";
        }
        for (j = order[i] - 1;j >= 0;j--) {
          os << state[i][j] << " ";
        }

        for (j = 0;j < nb_memory[i];j++) {
          os << "  ";
          for (k = max_order - 1;k >= order[previous[i][j]];k--) {
            os << "  ";
          }
          for (k = order[previous[i][j]] - 1;k >= 0;k--) {
            os << state[previous[i][j]][k] << " ";
          }
        }
        os << endl;
      }
    }
  }
# endif

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet VariableOrderMarkov et de la structure de donnees associee.
 *
 *  arguments : stream, pointeur sur un objet VariableOrderMarkovData,
 *              flag niveau de detail, flag fichier, flag Markov cache.
 *
 *--------------------------------------------------------------*/

ostream& VariableOrderMarkov::ascii_write(ostream &os , const VariableOrderMarkovData *seq ,
                                          bool exhaustive , bool file_flag , bool hidden) const

{
  register int i , j;
  int buff , variable , max_memory_count , *memory_count , width[2];
  double standard_normal_value , half_confidence_interval;
  long old_adjust;
  FrequencyDistribution **observation = NULL;
  SequenceCharacteristics *characteristics;


  old_adjust = os.setf(ios::left , ios::adjustfield);

  switch (hidden) {

  case false : {
    switch (type) {
    case 'o' :
      os << SEQ_word[SEQW_MARKOV_CHAIN] << endl;
      break;
    case 'e' :
      os << SEQ_word[SEQW_EQUILIBRIUM_MARKOV_CHAIN] << endl;
      break;
    }
    break;
  }

  case true : {
    switch (type) {
    case 'o' :
      os << SEQ_word[SEQW_HIDDEN_MARKOV_CHAIN] << endl;
      break;
    case 'e' :
      os << SEQ_word[SEQW_EQUILIBRIUM_HIDDEN_MARKOV_CHAIN] << endl;
      break;
    }
    break;
  }
  }

  // ecriture des parametres de la chaine de Markov d'ordre variable

  ascii_print(os , file_flag);

//  if ((nb_component == 1) && (seq) && ((!hidden) || (seq->type[0] == STATE))) {
  if ((seq) && ((!hidden) || (seq->type[0] == STATE))) {
    standard_normal_value = standard_normal_value_computation(0.025);

    if (!hidden) {
      width[0] = 0;
      for (i = 1;i < nb_row;i++) {
        if((memory_type[i] == TERMINAL) || ((type == 'o') &&
            (memory_type[i] == NON_TERMINAL))) {
          buff = column_width(nb_state , transition[i]);
          if (buff > width[0]) {
            width[0] = buff;
          }
        }
      }
      width[0]++;
    }

    memory_count = new int[nb_row];
    max_memory_count = 0;
    for (i = 1;i < nb_row;i++) {
//      if (memory_type[i] == TERMINAL) {
      if ((memory_type[i] == TERMINAL) || ((type == 'o') &&
           (memory_type[i] == NON_TERMINAL)) || (hidden)) {
        memory_count[i] = 0;
        for (j = 0;j < nb_state;j++) {
          memory_count[i] += seq->chain_data->transition[i][j];
        }

        if (memory_count[i] > max_memory_count) {
          max_memory_count = memory_count[i];
        }
      }
    }
    width[1] = column_width(max_memory_count) + ASCII_SPACE;

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    switch (hidden) {
    case false :
      os << SEQ_label[SEQL_TRANSITION_PROBABILITIY_CONFIDENCE_INTERVAL] << endl;
      break;
    case true :
      os << SEQ_label[SEQL_TRANSITION_COUNTS] << endl;
      break;
    }

    os << "\n";
    for (i = 1;i < nb_row;i++) {
//      if (memory_type[i] == TERMINAL) {
      if ((memory_type[i] == TERMINAL) || ((type == 'o') &&
           (memory_type[i] == NON_TERMINAL)) || (hidden)) {
        if (memory_count[i] > 0.) {
          if (file_flag) {
            os << "# ";
          }

          switch (hidden) {

          case false : {
            if (memory_type[i] == TERMINAL) {
              for (j = 0;j < nb_state;j++) {
                if ((transition[i][j] > 0.) && (transition[i][j] < 1.)) {
                  half_confidence_interval = standard_normal_value *
                                             sqrt(transition[i][j] * (1. - transition[i][j]) / memory_count[i]);
                  os << setw(width[0]) << MAX(transition[i][j] - half_confidence_interval , 0.)
                     << setw(width[0]) << MIN(transition[i][j] + half_confidence_interval , 1.)
                     << "| ";
                }
                else {
                  os << setw(width[0]) << " "
                     << setw(width[0]) << " "
                     << "| ";
                }
              }
            }

            else {
              for (j = 0;j < nb_state;j++) {
                os << setw(width[0]) << " "
                   << setw(width[0]) << " "
                   << "| ";
              }
            }
            break;
          }

          case true : {
            for (j = 0;j < nb_state;j++) {
              os << setw(width[1]) << seq->chain_data->transition[i][j];
            }
            os << "  ";
            break;
          }
          }

          os << setw(width[1]) << memory_count[i] << "  ";

          for (j = max_order - 1;j >= order[i];j--) {
            os << "  ";
          }
          for (j = order[i] - 1;j >= 0;j--) {
            os << state[i][j] << " ";
          }

          os << endl;
        }
      }
    }

    delete [] memory_count;
  }

  if ((seq) && (seq->type[0] == STATE)) {
    characteristics = seq->characteristics[0];
  }
  else {
    characteristics = NULL;
  }

  nonparametric_process[0]->ascii_print(os , 0 , 0 , characteristics , exhaustive , file_flag);

  if (hidden) {
    if ((type == 'o') && (nonparametric_process[0]->index_value)) {
      register int j;
      double sum , *state_marginal;


      state_marginal = new double[nb_state];

      for (i = 0;i < nb_state;i++) {
        state_marginal[i] = 0.;
      }

      for (i = 0;i < nonparametric_process[0]->length->nb_value - 1;i++) {
        for (j = 0;j < nb_state;j++) {
          state_marginal[j] += nonparametric_process[0]->index_value->point[j][i] *
                               (1. - nonparametric_process[0]->length->cumul[i]);
        }
      }

      sum = 0.;
      for (i = 0;i < nb_state;i++) {
        sum += state_marginal[i];
      }

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_STATE_PROBABILITY] << endl;
      if (file_flag) {
        os << "# ";
      }
      for (i = 0;i < nb_state;i++) {
        os << state_marginal[i] / sum  << "  ";
      }
      os << endl;

      delete state_marginal;
    }

    os << "\n" << nb_output_process << " "
       << STAT_word[nb_output_process == 1 ? STATW_OUTPUT_PROCESS : STATW_OUTPUT_PROCESSES] << endl;
  }

  // ecriture des lois associees a chaque processus d'observation

  for (i = 1;i <= nb_output_process;i++) {
    os << "\n" << STAT_word[STATW_OUTPUT_PROCESS];

    if (hidden) {
      os << " " << i;

      if (nonparametric_process[i]) {
        os << " : " << STAT_word[STATW_NONPARAMETRIC];
      }
      else {
        os << " : " << STAT_word[STATW_PARAMETRIC];
      }
    }
    os << endl;

    if (seq) {
      switch (seq->type[0]) {
      case INT_VALUE :
        variable = i - 1;
        break;
      case STATE :
        variable = i;
        break;
      }

      if (seq->observation) {
        observation = seq->observation[variable];
      }

      if (seq->characteristics[variable]) {
        characteristics = seq->characteristics[variable];
      }
      else {
        characteristics = NULL;
      }
    }

    if (nonparametric_process[i]) {
      nonparametric_process[i]->ascii_print(os , i , observation , characteristics ,
                                            exhaustive , file_flag);
    }
    else {
      parametric_process[i]->ascii_print(os , observation , exhaustive , file_flag);
    }
  }

  if (seq) {
    int nb_parameter = nb_parameter_computation(hidden ? MIN_PROBABILITY : 0.) , nb_transient_parameter;
    double information , likelihood;


    // ecriture de la loi empirique des longueurs des sequences

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    seq->hlength->ascii_characteristic_print(os , false , file_flag);

    if (exhaustive) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << "   | " << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
      seq->hlength->ascii_print(os , file_flag);
    }

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_CUMUL_LENGTH] << ": " << seq->cumul_length << endl;

    // ecriture de la quantite d'information des sequences dans le cas iid

    information = seq->iid_information_computation();

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_IID_INFORMATION] << ": " << information << " ("
       << information / seq->cumul_length << ")" << endl;

    // ecriture des vraisemblances des sequences

    switch (hidden) {

    case false : {
      likelihood = seq->likelihood;

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   ("
         << STAT_label[STATL_NORMALIZED] << ": " << likelihood / seq->cumul_length << ")" << endl;
      break;
    }

    case true : {
      if (seq->likelihood != D_INF) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << ": " << seq->likelihood << "   ("
           << STAT_label[STATL_NORMALIZED] << ": " << seq->likelihood / seq->cumul_length << ")" << endl;
      }

      likelihood = seq->hidden_likelihood;

      if (likelihood != D_INF) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_OBSERVED_SEQUENCES_LIKELIHOOD] << ": " << likelihood << "   ("
           << STAT_label[STATL_NORMALIZED] << ": " << likelihood / seq->cumul_length << ")" << endl;
      }
      break;
    }
    }

    if (likelihood != D_INF) {
      if (type == 'o') {
        nb_transient_parameter = nb_transient_parameter_computation(hidden ? MIN_PROBABILITY : 0.);

        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << nb_transient_parameter << " "
           << SEQ_label[nb_transient_parameter == 1 ? SEQL_FREE_TRANSIENT_PARAMETER : SEQL_FREE_TRANSIENT_PARAMETERS] << endl;
      }

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AIC] << "): "
         << 2 * (likelihood - nb_parameter) << endl;

      if ((type == 'o') && (nb_transient_parameter > 0)) {
        if (file_flag) {
          os << "# ";
        }
        os << nb_transient_parameter + nb_parameter << " " << STAT_label[STATL_FREE_PARAMETERS]
           << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AIC] << "): "
           << 2 * (likelihood - nb_transient_parameter - nb_parameter) << endl;
      }

      if (nb_parameter < seq->cumul_length - 1) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
           << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AICc] << "): "
           << 2 * (likelihood - (double)(nb_parameter * seq->cumul_length) /
               (double)(seq->cumul_length - nb_parameter - 1)) << endl;
      }

      if ((type == 'o') && (nb_transient_parameter > 0) &&
          (nb_transient_parameter + nb_parameter < seq->cumul_length - 1)) {
        if (file_flag) {
          os << "# ";
        }
        os << nb_transient_parameter + nb_parameter << " " << STAT_label[STATL_FREE_PARAMETERS]
           << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AICc] << "): "
           << 2 * (likelihood - (double)((nb_transient_parameter + nb_parameter) * seq->cumul_length) /
               (double)(seq->cumul_length - nb_transient_parameter - nb_parameter - 1)) << endl;
      }

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BIC] << "): "
         << 2 * likelihood - nb_parameter * log((double)seq->cumul_length) << endl;

      if ((type == 'o') && (nb_transient_parameter > 0)) {
        if (file_flag) {
          os << "# ";
        }
        os << nb_transient_parameter + nb_parameter << " " << STAT_label[STATL_FREE_PARAMETERS]
           << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BIC] << "): "
           << 2 * likelihood - (nb_transient_parameter + nb_parameter) * log((double)seq->cumul_length) << endl;
      }

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << nb_parameter + (type == 'o' ? nb_transient_parameter : 0) << " " << STAT_label[STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BICc] << "): "
         << 2 * likelihood - penalty_computation(hidden , (hidden ? MIN_PROBABILITY : 0.)) << endl;
    }
  }

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet VariableOrderMarkov.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& VariableOrderMarkov::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , markov_data , exhaustive , false);
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet VariableOrderMarkov dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool VariableOrderMarkov::ascii_write(StatError &error , const char *path ,
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
    ascii_write(out_file , markov_data , exhaustive , true);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie au format tableur des parametres de la chaine de Markov d'ordre variable.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& VariableOrderMarkov::spreadsheet_print(ostream &os) const

{
  register int i , j , k;
  double *stationary_probability;


  os << "\n" << nb_state << "\t" << STAT_word[STATW_STATES] << endl;

  switch (type) {

  case 'o' : {
    os << "\n" << STAT_word[STATW_INITIAL_PROBABILITIES] << endl;
    for (i = 0;i < nb_state;i++) {
      os << initial[i] << "\t";
    }
    os << endl;
    break;
  }

  case 'e' : {

    // extraction des probabilites stationnairees correspondant aux memoires completees

    stationary_probability = new double[nb_row];

    for (i = 1;i < nb_row;i++) {
      stationary_probability[i] = initial[i];
    }
    for (i = nb_row - 1;i >= 1;i--) {
      if (memory_type[i] == COMPLETION) {
        stationary_probability[parent[i]] += stationary_probability[i];
      }
    }

    os << "\n" << STAT_label[STATL_STATIONARY_PROBABILITIES] << endl;
    for (i = 1;i < nb_row;i++) {
      if (memory_type[i] == TERMINAL) {
        os << stationary_probability[i] << "\t\t";
        for (j = order[i] - 1;j >= 0;j--) {
          os << state[i][j] << " ";
        }
        os << endl;
      }
    }

    delete [] stationary_probability;
    break;
  }
  }

  os << "\n" << STAT_word[STATW_TRANSITION_PROBABILITIES]
     << "\t\t" << STAT_label[STATL_MEMORY] << endl;

  for (i = 1;i < nb_row;i++) {
    for (j = 0;j < nb_state;j++) {
      os << transition[i][j] << "\t";
    }

    os << "\t";
    for (j = order[i] - 1;j >= 0;j--) {
      os << state[i][j] << " ";
    }

    os << "\t";
    switch (memory_type[i]) {

    case NON_TERMINAL : {
      os << SEQ_label[SEQL_NON_TERMINAL];
      break;
    }

    case TERMINAL : {
      os << SEQ_label[SEQL_TERMINAL];
      if (child[i]) {
        os << "\t" << SEQ_label[SEQL_COMPLETED];
      }
      break;
    }

    case COMPLETION : {
      os << SEQ_label[SEQL_COMPLETION] << "\t"
         << (child[i] ? SEQ_label[SEQL_NON_TERMINAL] : SEQ_label[SEQL_TERMINAL]);
      break;
    }
    }

    os << endl;
  }

  if (nb_component > 0) {
    for (i = 0;i < nb_component;i++) {
      switch (state_type[component[i][0]]) {
      case 't' :
        os << "\n"<< STAT_label[STATL_TRANSIENT] << " ";
        break;
      default :
        os << "\n" << STAT_label[STATL_RECURRENT] << " ";
        break;
      }
      os << STAT_label[STATL_CLASS] << "\t" << STAT_label[component_nb_state[i] == 1 ? STATL_STATE : STATL_STATES];

      for (j = 0;j < component_nb_state[i];j++) {
        os << "\t" << component[i][j];
      }

      if (state_type[component[i][0]] == 'a') {
        os << " (" << STAT_label[STATL_ABSORBING] << " " << STAT_label[STATL_STATE] << ")";
      }
    }
    os << endl;
  }

  os << SEQ_label[SEQL_MEMORY_TRANSITION_MATRIX] << endl;

  os << "\n";
  for (i = 1;i < nb_row;i++) {
    if ((type == 'o') || (!child[i])) {
      for (j = order[i] - 1;j >= 0;j--) {
        os << state[i][j] << " ";
      }

      os << "\t";
      for (j = 0;j < nb_state;j++) {
        os << "\t";
        for (k = order[next[i][j]] - 1;k >= 0;k--) {
          os << state[next[i][j]][k] << " ";
        }
      }

      os << "\t";
      switch (memory_type[i]) {

      case NON_TERMINAL : {
        os << SEQ_label[SEQL_NON_TERMINAL];
        break;
      }

      case TERMINAL : {
        os << SEQ_label[SEQL_TERMINAL];
        if (child[i]) {
          os << "\t" << SEQ_label[SEQL_COMPLETED];
        }
        break;
      }

      case COMPLETION : {
        os << SEQ_label[SEQL_COMPLETION] << "\t"
           << (child[i] ? SEQ_label[SEQL_NON_TERMINAL] : SEQ_label[SEQL_TERMINAL]);
        break;
      }
      }

      os << endl;
    }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet VariableOrderMarkov et de la structure de donnees associee
 *  dans un fichier au format tableur.
 *
 *  arguments : stream, pointeur sur un objet VariableOrderMarkovData,
 *              flag Markov cache.
 *
 *--------------------------------------------------------------*/

ostream& VariableOrderMarkov::spreadsheet_write(ostream &os ,
                                                const VariableOrderMarkovData *seq ,
                                                bool hidden) const

{
  register int i;
  int variable;
  FrequencyDistribution **observation = NULL;
  SequenceCharacteristics *characteristics;


  switch (hidden) {

  case false : {
    switch (type) {
    case 'o' :
      os << SEQ_word[SEQW_MARKOV_CHAIN] << endl;
      break;
    case 'e' :
      os << SEQ_word[SEQW_EQUILIBRIUM_MARKOV_CHAIN] << endl;
      break;
    }
    break;
  }

  case true : {
    switch (type) {
    case 'o' :
      os << SEQ_word[SEQW_HIDDEN_MARKOV_CHAIN] << endl;
      break;
    case 'e' :
      os << SEQ_word[SEQW_EQUILIBRIUM_HIDDEN_MARKOV_CHAIN] << endl;
      break;
    }
    break;
  }
  }

  // ecriture des parametres de la chaine de Markov d'ordre variable

  spreadsheet_print(os);

  if ((seq) && (seq->type[0] == STATE)) {
    characteristics = seq->characteristics[0];
  }
  else {
    characteristics = NULL;
  }

  nonparametric_process[0]->spreadsheet_print(os , 0 , 0 , characteristics);

  // ecriture des lois associees a chaque processus d'observation

  if (hidden) {
    os << "\n" << nb_output_process << "\t"
       << STAT_word[nb_output_process == 1 ? STATW_OUTPUT_PROCESS : STATW_OUTPUT_PROCESSES] << endl;
  }

  for (i = 1;i <= nb_output_process;i++) {
    os << "\n" << STAT_word[STATW_OUTPUT_PROCESS];

    if (hidden) {
      os << "\t" << i;

      if (nonparametric_process[i]) {
        os << "\t" << STAT_word[STATW_NONPARAMETRIC];
      }
      else {
        os << "\t" << STAT_word[STATW_PARAMETRIC];
      }
    }
    os << endl;

    if (seq) {
      switch (seq->type[0]) {
      case INT_VALUE :
        variable = i - 1;
        break;
      case STATE :
        variable = i;
        break;
      }

      if (seq->observation) {
        observation = seq->observation[variable];
      }

      if (seq->characteristics[variable]) {
        characteristics = seq->characteristics[variable];
      }
      else {
        characteristics = NULL;
      }
    }

    if (nonparametric_process[i]) {
      nonparametric_process[i]->spreadsheet_print(os , i , observation , characteristics);
    }
    else {
      parametric_process[i]->spreadsheet_print(os , observation);
    }
  }

  if (seq) {
    int nb_parameter = nb_parameter_computation(hidden ? MIN_PROBABILITY : 0.) , nb_transient_parameter;
    double information , likelihood;


    // ecriture de la loi empirique des longueurs des sequences

    os << "\n" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
    seq->hlength->spreadsheet_characteristic_print(os);

    os << "\n\t" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
    seq->hlength->spreadsheet_print(os);

    os << "\n" << SEQ_label[SEQL_CUMUL_LENGTH] << "\t" << seq->cumul_length << endl;

    // ecriture de la quantite d'information des sequences dans le cas iid

    information = seq->iid_information_computation();

    os << "\n" << SEQ_label[SEQL_IID_INFORMATION] << "\t" << information << "\t"
       << information / seq->cumul_length << endl;

    // ecriture des vraisemblances des sequences

    switch (hidden) {

    case false : {
      likelihood = seq->likelihood;

      os << "\n" << STAT_label[STATL_LIKELIHOOD] << "\t" << likelihood << "\t"
         << STAT_label[STATL_NORMALIZED] << "\t" << likelihood / seq->cumul_length << endl;
      break;
    }

    case true : {
      if (seq->likelihood != D_INF) {
        os << "\n" << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << "\t" << seq->likelihood << "\t"
           << STAT_label[STATL_NORMALIZED] << "\t" << seq->likelihood / seq->cumul_length << endl;
      }

      likelihood = seq->hidden_likelihood;

      if (likelihood != D_INF) {
        os << "\n" << SEQ_label[SEQL_OBSERVED_SEQUENCES_LIKELIHOOD] << "\t" << likelihood << "\t"
           << STAT_label[STATL_NORMALIZED] << "\t" << likelihood / seq->cumul_length << endl;
      }
      break;
    }
    }

    if (likelihood != D_INF) {
      if (type == 'o') {
        nb_transient_parameter = nb_transient_parameter_computation(hidden ? MIN_PROBABILITY : 0.);

        os << "\n" << nb_transient_parameter << "\t"
           << SEQ_label[nb_transient_parameter == 1 ? SEQL_FREE_TRANSIENT_PARAMETER : SEQL_FREE_TRANSIENT_PARAMETERS] << endl;
      }

      os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AIC] << ")\t"
         << 2 * (likelihood - nb_parameter) << endl;
      if ((type == 'o') && (nb_transient_parameter > 0)) {
        os << nb_transient_parameter + nb_parameter << "\t" << STAT_label[STATL_FREE_PARAMETERS] << "\t"
           << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AIC] << ")\t"
           << 2 * (likelihood - nb_transient_parameter - nb_parameter) << endl;
      }

      if (nb_parameter < seq->cumul_length - 1) {
        os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
           << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AICc] << ")\t"
           << 2 * (likelihood - (double)(nb_parameter * seq->cumul_length) /
               (double)(seq->cumul_length - nb_parameter - 1)) << endl;
      }
      if ((type == 'o') && (nb_transient_parameter > 0) &&
          (nb_transient_parameter + nb_parameter < seq->cumul_length - 1)) {
        os << nb_transient_parameter + nb_parameter << "\t" << STAT_label[STATL_FREE_PARAMETERS] << "\t"
           << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AICc] << ")\t"
           << 2 * (likelihood - (double)((nb_transient_parameter + nb_parameter) * seq->cumul_length) /
               (double)(seq->cumul_length - nb_transient_parameter - nb_parameter - 1)) << endl;
      }

      os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BIC] << ")\t"
         << 2 * likelihood - nb_parameter * log((double)seq->cumul_length) << endl;
      if ((type == 'o') && (nb_transient_parameter > 0)) {
        os << nb_transient_parameter + nb_parameter << "\t" << STAT_label[STATL_FREE_PARAMETERS] << "\t"
           << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BIC] << ")\t"
           << 2 * likelihood - (nb_transient_parameter + nb_parameter) * log((double)seq->cumul_length) << endl;
      }

      os << "\n" << nb_parameter + (type == 'o' ? nb_transient_parameter : 0) << "\t" << STAT_label[STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BICc] << ")\t"
         << 2 * likelihood - penalty_computation(hidden , (hidden ? MIN_PROBABILITY : 0.)) << endl;
    }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet VariableOrderMarkov dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool VariableOrderMarkov::spreadsheet_write(StatError &error , const char *path) const

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
    spreadsheet_write(out_file , markov_data);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet VariableOrderMarkov et
 *  de la structure de donnees associee.
 *
 *  arguments : prefixe des fichiers, titre des figures,
 *              pointeur sur les sequences observees.
 *
 *--------------------------------------------------------------*/

bool VariableOrderMarkov::plot_write(const char *prefix , const char *title ,
                                     const VariableOrderMarkovData *seq) const

{
  bool status;
  register int i;
  int variable;
  FrequencyDistribution *hlength = NULL , **observation = NULL;
  SequenceCharacteristics *characteristics;


  if ((seq) && (seq->type[0] == STATE)) {
    characteristics = seq->characteristics[0];
    hlength = seq->hlength;
  }
  else {
    characteristics = NULL;
  }

  status = nonparametric_process[0]->plot_print(prefix , title , 0 , 0 , characteristics ,
                                                hlength);

  if (status) {
    if (seq) {
      hlength = seq->hlength;
    }

    for (i = 1;i <= nb_output_process;i++) {
      if (seq) {
        switch (seq->type[0]) {
        case INT_VALUE :
          variable = i - 1;
          break;
        case STATE :
          variable = i;
          break;
        }

        if (seq->observation) {
          observation = seq->observation[variable];
        }

        if (seq->characteristics[variable]) {
          characteristics = seq->characteristics[variable];
        }
        else {
          characteristics = NULL;
        }
      }

      if (nonparametric_process[i]) {
        nonparametric_process[i]->plot_print(prefix , title , i , observation ,
                                             characteristics , hlength);
      }
      else {
        parametric_process[i]->plot_print(prefix , title , i , observation);
      }
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet VariableOrderMarkov.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool VariableOrderMarkov::plot_write(StatError &error , const char *prefix ,
                                     const char *title) const

{
  bool status = plot_write(prefix , title , markov_data);

  error.init();

  if (!status) {
    error.update(STAT_error[STATR_FILE_PREFIX]);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet VariableOrderMarkov et
 *  de la structure de donnees associee.
 *
 *  arguments : pointeur sur les sequences observees.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* VariableOrderMarkov::get_plotable(const VariableOrderMarkovData *seq) const

{
  register int i , j;
  int nb_plot_set , index_length , index , variable;
  FrequencyDistribution *hlength = NULL , **observation = NULL;
  SequenceCharacteristics *characteristics;
  MultiPlotSet *plot_set;


  if ((seq) && (seq->type[0] == STATE)) {
    characteristics = seq->characteristics[0];
  }
  else {
    characteristics = NULL;
  }

  // calcul du nombre de vues

  nb_plot_set = 0;

  if ((nonparametric_process[0]->index_value) || (characteristics)) {
    nb_plot_set++;

    if (characteristics) {
      index_length = characteristics->index_value->plot_length_computation();

      if (characteristics->index_value->frequency[index_length - 1] < MAX_FREQUENCY) {
        nb_plot_set++;
      }
      nb_plot_set++;
    }
  }

  if ((nonparametric_process[0]->first_occurrence) || (characteristics)) {
    for (i = 0;i < nb_state;i++) {
      if ((nonparametric_process[0]->first_occurrence) &&
          (nonparametric_process[0]->first_occurrence[i])) {
        nb_plot_set++;
      }
      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->first_occurrence[i]->nb_element > 0)) {
        nb_plot_set++;
      }
    }
  }

  if ((nonparametric_process[0]->recurrence_time) || (characteristics)) {
    for (i = 0;i < nb_state;i++) {
      if ((nonparametric_process[0]->recurrence_time) &&
          (nonparametric_process[0]->recurrence_time[i])) {
        nb_plot_set++;
      }
      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->recurrence_time[i]->nb_element > 0)) {
        nb_plot_set++;
      }
    }
  }

  if ((nonparametric_process[0]->sojourn_time) || (characteristics)) {
    for (i = 0;i < nb_state;i++) {
      if ((nonparametric_process[0]->sojourn_time) &&
          (nonparametric_process[0]->sojourn_time[i])) {
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

      if ((characteristics) && (i < characteristics->nb_value) &&
          (characteristics->final_run[i]->nb_element > 0)) {
        nb_plot_set++;
      }
    }
  }

  if ((nonparametric_process[0]->nb_run) || (nonparametric_process[0]->nb_occurrence) ||
      ((characteristics) && (characteristics->nb_run) && (characteristics->nb_occurrence))) {
    for (i = 0;i < nb_state;i++) {
      if (nonparametric_process[0]->nb_run) {
        nb_plot_set++;
      }
      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->nb_run) && (characteristics->nb_run[i]->nb_element > 0)) {
        nb_plot_set++;
      }

      if (nonparametric_process[0]->nb_occurrence) {
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

  for (i = 1;i <= nb_output_process;i++) {
    if (seq) {
      switch (seq->type[0]) {
      case INT_VALUE :
        variable = i - 1;
        break;
      case STATE :
        variable = i;
        break;
      }

      if (seq->characteristics[variable]) {
        characteristics = seq->characteristics[variable];
      }
      else {
        characteristics = NULL;
      }
    }

    if (nonparametric_process[i]) {
      if ((nonparametric_process[i]->index_value) || (characteristics)) {
        nb_plot_set++;

        if (characteristics) {
          index_length = characteristics->index_value->plot_length_computation();

          if (characteristics->index_value->frequency[index_length - 1] < MAX_FREQUENCY) {
            nb_plot_set++;
          }
          nb_plot_set++;
        }
      }

      if ((nonparametric_process[i]->first_occurrence) || (characteristics)) {
        for (j = 0;j < nonparametric_process[i]->nb_value;j++) {
          if ((nonparametric_process[i]->first_occurrence) &&
              (nonparametric_process[i]->first_occurrence[j])) {
            nb_plot_set++;
          }
          else if ((characteristics) && (j < characteristics->nb_value) &&
                   (characteristics->first_occurrence[j]->nb_element > 0)) {
            nb_plot_set++;
          }
        }
      }

      if ((nonparametric_process[i]->recurrence_time) || (characteristics)) {
        for (j = 0;j < nonparametric_process[i]->nb_value;j++) {
          if ((nonparametric_process[i]->recurrence_time) &&
              (nonparametric_process[i]->recurrence_time[j])) {
            nb_plot_set++;
          }
          else if ((characteristics) && (i < characteristics->nb_value) &&
                   (characteristics->recurrence_time[j]->nb_element > 0)) {
            nb_plot_set++;
          }
        }
      }

      if ((nonparametric_process[i]->sojourn_time) || (characteristics)) {
        for (j = 0;j < nonparametric_process[i]->nb_value;j++) {
          if ((nonparametric_process[i]->sojourn_time) &&
              (nonparametric_process[i]->sojourn_time[j])) {
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

      if ((nonparametric_process[i]->nb_run) || (nonparametric_process[i]->nb_occurrence) ||
          ((characteristics) && (characteristics->nb_run) && (characteristics->nb_occurrence))) {
        for (j = 0;j < nonparametric_process[i]->nb_value;j++) {
          if (nonparametric_process[i]->nb_run) {
            nb_plot_set++;
          }
          else if ((characteristics) && (j < characteristics->nb_value) &&
                   (characteristics->nb_run) && (characteristics->nb_run[j]->nb_element > 0)) {
            nb_plot_set++;
          }

          if (nonparametric_process[i]->nb_occurrence) {
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

    nb_plot_set += nb_state;
  }

  plot_set = new MultiPlotSet(nb_plot_set , nb_output_process + 1);
  plot_set->border = "15 lw 0";

  if ((seq) && (seq->type[0] == STATE)) {
    characteristics = seq->characteristics[0];
    hlength = seq->hlength;
  }
  else {
    characteristics = NULL;
  }

  return plot_set;
  index = 0;
  plot_set->variable_nb_viewpoint[0] = 0;
  nonparametric_process[0]->plotable_write(*plot_set , index , 0 , 0 , characteristics ,
                                           hlength);

  if (seq) {
    hlength = seq->hlength;
  }

  for (i = 1;i <= nb_output_process;i++) {
    if (seq) {
      switch (seq->type[0]) {
      case INT_VALUE :
        variable = i - 1;
        break;
      case STATE :
        variable = i;
        break;
      }

      if (seq->observation) {
        observation = seq->observation[variable];
      }

      if (seq->characteristics[variable]) {
        characteristics = seq->characteristics[variable];
      }
      else {
        characteristics = NULL;
      }
    }

    if (nonparametric_process[i]) {
      plot_set->variable_nb_viewpoint[i] = 0;
      nonparametric_process[i]->plotable_write(*plot_set , index , i , observation ,
                                               characteristics , hlength);
    }
    else {
      parametric_process[i]->plotable_write(*plot_set , index , i , observation);
    }
  }

  return plot_set;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet VariableOrderMarkov.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* VariableOrderMarkov::get_plotable() const

{
  return get_plotable(markov_data);
}


/*--------------------------------------------------------------*
 *
 *  Calcul de l'ordre maximum des memoires.
 *
 *--------------------------------------------------------------*/

void VariableOrderMarkov::max_order_computation()

{
  register int i;


  max_order = 0;
  for (i = 0;i < nb_row;i++) {
    if ((memory_type[i] == TERMINAL) && (order[i] > max_order)) {
      max_order = order[i];
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul du nombre de parametres independants.
 *
 *  argument : probabilite minimum.
 *
 *--------------------------------------------------------------*/

int VariableOrderMarkov::nb_parameter_computation(double min_probability) const

{
  register int i , j;
  int nb_parameter = 0;


  // cas particulier ordre 0

  if (max_order == 1) {
    for (i = 0;i < nb_state - 1;i++) {
      for (j = 2;j <= nb_state;j++) {
        if (transition[j][i] != transition[1][i]) {
          break;
        }
      }

      if (j <= nb_state) {
        break;
      }
    }

    if (i == nb_state - 1) {
      for (i = 0;i < nb_state;i++) {
        if (transition[0][i] > min_probability) {
          nb_parameter++;
        }
      }

      nb_parameter--;
    }
  }

  if (nb_parameter == 0) {
    for (i = 1;i < nb_row;i++) {
      if (memory_type[i] == TERMINAL) {
//      if ((memory_type[i] == TERMINAL) || ((type == 'o') &&
//           (memory_type[i] == NON_TERMINAL))) {
        for (j = 0;j < nb_state;j++) {
          if (transition[i][j] > min_probability) {
            nb_parameter++;
          }
        }

        nb_parameter--;
      }
    }
  }

  for (i = 1;i <= nb_output_process;i++) {
    if (nonparametric_process[i]) {
      nb_parameter += nonparametric_process[i]->nb_parameter_computation(min_probability);
    }
    else {
      nb_parameter += parametric_process[i]->nb_parameter_computation();
    }
  }

  return nb_parameter;
}


/*--------------------------------------------------------------*
 *
 *  Calcul du nombre de parametres transitoires independants correspondant
 *  aux lois de transition des memoires non-terminales pour une chaine de Markov
 *  d'ordre variable ordinaire.
 *
 *  argument : probabilite minimum.
 *
 *--------------------------------------------------------------*/

int VariableOrderMarkov::nb_transient_parameter_computation(double min_probability) const

{
  register int i , j;
  int nb_parameter = 0;


  if (type == 'o') {
    for (i = 1;i < nb_row;i++) {
      if (memory_type[i] == NON_TERMINAL) {
        for (j = 0;j < nb_state;j++) {
          if (transition[i][j] > min_probability) {
            nb_parameter++;
          }
        }

        nb_parameter--;
      }
    }
  }

  return nb_parameter;
}


/*--------------------------------------------------------------*
 *
 *  Calcul d'une penalite adaptative.
 *
 *  arguments : flag Markov cache, probabilite minimum.
 *
 *--------------------------------------------------------------*/

double VariableOrderMarkov::penalty_computation(bool hidden , double min_probability) const

{
  register int i , j , k;
  int nb_parameter , sample_size;
  double sum , *state_marginal , *memory;
  double penalty = 0.;


  if (markov_data) {
    if (hidden) {
      switch (type) {

      case 'o' : {
        memory = memory_computation();

        sum = 0.;
        for (i = 1;i < nb_row;i++) {
//          if (memory_type[i] == TERMINAL) {
            sum += memory[i];
//          }
        }
        for (i = 1;i < nb_row;i++) {
//          if (memory_type[i] == TERMINAL) {
            memory[i] /= sum;
//          }
        }
        break;
      }

      case 'e' : {
        memory = new double[nb_row];
        for (i = 1;i < nb_row;i++) {
          memory[i] = initial[i];
        }
        break;
      }
      }

      state_marginal = new double[nb_state];
      for (i = 0;i < nb_state;i++) {
        state_marginal[i] = 0.;
      }
      for (i = 1;i < nb_row;i++) {
//        if (memory_type[i] == TERMINAL) {
        state_marginal[state[i][0]] += memory[i];
//        }
      }
    }

    for (i = 1;i < nb_row;i++) {
//      if (memory_type[i] == TERMINAL) {
      if ((memory_type[i] == TERMINAL) || ((type == 'o') &&
           (memory_type[i] == NON_TERMINAL))) {
        nb_parameter = 0;
        if (!hidden) {
          sample_size = 0;
        }
        for (j = 0;j < nb_state;j++) {
          if (transition[i][j] > min_probability) {
            nb_parameter++;
            if (!hidden) {
              sample_size += markov_data->chain_data->transition[i][j];
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
              penalty += nb_parameter * log(memory[i] * markov_data->cumul_length);
            }
            break;
          }
          }
        }
      }
    }

    for (i = 1;i <= nb_output_process;i++) {
      if (nonparametric_process[i]) {
        for (j = 0;j < nb_state;j++) {
          nb_parameter = 0;
          for (k = 0;k < nonparametric_process[i]->nb_value;k++) {
            if (nonparametric_process[i]->observation[j]->mass[k] > min_probability) {
              nb_parameter++;
            }
          }

          nb_parameter--;

          if (nb_parameter > 0) {
            switch (hidden) {
            case false :
              penalty += nb_parameter *
                         log((double)markov_data->marginal_distribution[0]->frequency[j]);
              break;
            case true :
              penalty += nb_parameter * log(state_marginal[j] * markov_data->cumul_length);
              break;
            }
          }
        }
      }

      else {
        for (j = 0;j < nb_state;j++) {
          nb_parameter = parametric_process[i]->observation[j]->nb_parameter_computation();

          switch (hidden) {
            case false :
            penalty += nb_parameter *
                       log((double)markov_data->marginal_distribution[0]->frequency[j]);
            break;
          case true :
            penalty += nb_parameter * log(state_marginal[j] * markov_data->cumul_length);
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
 *  Constructeur par defaut de la classe VariableOrderMarkovData.
 *
 *--------------------------------------------------------------*/

VariableOrderMarkovData::VariableOrderMarkovData()

{
  markov = NULL;
  chain_data = NULL;

  likelihood = D_INF;
  hidden_likelihood = D_INF;

  posterior_probability = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe VariableOrderMarkovData.
 *
 *  arguments : loi empirique des longueurs des sequences, nombre de variables,
 *              flag initialisation.
 *
 *--------------------------------------------------------------*/

VariableOrderMarkovData::VariableOrderMarkovData(const FrequencyDistribution &ihlength ,
                                                 int inb_variable , bool init_flag)
:MarkovianSequences(ihlength , inb_variable , init_flag)

{
  markov = NULL;
  chain_data = NULL;

  likelihood = D_INF;
  hidden_likelihood = D_INF;

  posterior_probability = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet VariableOrderMarkovData a partir
 *  d'un objet MarkovianSequences avec ajout d'une variable d'etat.
 *
 *  argument : reference sur un objet MarkovianSequences.
 *
 *--------------------------------------------------------------*/

VariableOrderMarkovData::VariableOrderMarkovData(const MarkovianSequences &seq)
:MarkovianSequences(seq , 'a' , DEFAULT)

{
  markov = NULL;
  chain_data = NULL;

  likelihood = D_INF;
  hidden_likelihood = D_INF;

  posterior_probability = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet VariableOrderMarkovData a partir
 *  d'un objet MarkovianSequences.
 *
 *  arguments : reference sur un objet MarkovianSequences, type de transformation
 *              ('c' : copie, 'a' : ajout d'une variable d'etat),
 *              ajout/suppression des lois empiriques de temps de sejour initial.
 *
 *--------------------------------------------------------------*/

VariableOrderMarkovData::VariableOrderMarkovData(const MarkovianSequences &seq ,
                                                 char transform , bool initial_run_flag)
:MarkovianSequences(seq , transform , (initial_run_flag ? ADD_INITIAL_RUN : REMOVE_INITIAL_RUN))

{
  markov = NULL;
  chain_data = NULL;

  likelihood = D_INF;
  hidden_likelihood = D_INF;

  posterior_probability = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet VariableOrderMarkovData.
 *
 *  arguments : reference sur un objet VariableOrderMarkovData,
 *              flag copie de l'objet VariableOrderMarkov.
 *
 *--------------------------------------------------------------*/

void VariableOrderMarkovData::copy(const VariableOrderMarkovData &seq ,
                                   bool model_flag)

{
  register int i;


  if ((model_flag) && (seq.markov)) {
    markov = new VariableOrderMarkov(*(seq.markov) , false);
  }
  else {
    markov = NULL;
  }

  if (seq.chain_data) {
    chain_data = new VariableOrderChainData(*(seq.chain_data));
  }
  else {
    chain_data = NULL;
  }

  likelihood = seq.likelihood;
  hidden_likelihood = seq.hidden_likelihood;

  if (seq.posterior_probability) {
    posterior_probability = new double[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      posterior_probability[i] = seq.posterior_probability[i];
    }
  }
  else {
    posterior_probability = NULL;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe VariableOrderMarkovData.
 *
 *--------------------------------------------------------------*/

VariableOrderMarkovData::~VariableOrderMarkovData()

{
  delete markov;
  delete chain_data;

  delete [] posterior_probability;
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe VariableOrderMarkovData.
 *
 *  argument : reference sur un objet VariableOrderMarkovData.
 *
 *--------------------------------------------------------------*/

VariableOrderMarkovData& VariableOrderMarkovData::operator=(const VariableOrderMarkovData &seq)

{
  if (&seq != this) {
    delete markov;
    delete chain_data;

    delete [] posterior_probability;

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

DiscreteDistributionData* VariableOrderMarkovData::extract(StatError &error , int type ,
                                                           int variable , int value) const

{
  bool status = true;
  Distribution *pdist;
  DiscreteParametric *pparam;
  FrequencyDistribution *phisto;
  DiscreteDistributionData *histo;


  histo = NULL;
  error.init();

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
                      << SEQ_error[SEQR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }

      else {
        phisto = observation[variable][value];

        if (phisto->nb_element == 0) {
          status = false;
          error.update(STAT_error[STATR_EMPTY_SAMPLE]);
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
                      << value << " " << SEQ_error[SEQR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }

      else {
        switch (type) {
        case FIRST_OCCURRENCE :
          phisto = characteristics[variable]->first_occurrence[value];
          break;
        case RECURRENCE_TIME :
          phisto = characteristics[variable]->recurrence_time[value];
          break;
        case SOJOURN_TIME :
          phisto = characteristics[variable]->sojourn_time[value];
          break;
        case FINAL_RUN :
          phisto = characteristics[variable]->final_run[value];
          break;
        case NB_RUN :
          phisto = characteristics[variable]->nb_run[value];
          break;
        case NB_OCCURRENCE :
          phisto = characteristics[variable]->nb_occurrence[value];
          break;
        }

        if (phisto->nb_element == 0) {
          status = false;
          error.update(STAT_error[STATR_EMPTY_SAMPLE]);
        }
      }
    }
  }

  if (status) {
    pdist = NULL;
    pparam = NULL;

    switch (type) {

    case OBSERVATION : {
      if (markov->nonparametric_process[variable]) {
        pdist = markov->nonparametric_process[variable]->observation[value];
      }
      else {
        pparam = markov->parametric_process[variable]->observation[value];
      }
      break;
    }

    case FIRST_OCCURRENCE : {
      pdist = markov->nonparametric_process[variable]->first_occurrence[value];
      break;
    }

    case RECURRENCE_TIME : {
      pdist = markov->nonparametric_process[variable]->recurrence_time[value];
      break;
    }

    case SOJOURN_TIME : {
      pparam = markov->nonparametric_process[variable]->sojourn_time[value];
      break;
    }

    case NB_RUN : {
      pdist = markov->nonparametric_process[variable]->nb_run[value];
      break;
    }

    case NB_OCCURRENCE : {
      pdist = markov->nonparametric_process[variable]->nb_occurrence[value];
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

VariableOrderMarkovData* VariableOrderMarkovData::remove_index_parameter(StatError &error) const

{
  VariableOrderMarkovData *seq;


  error.init();

  if (!index_parameter) {
    seq = NULL;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }
  else {
    seq = new VariableOrderMarkovData(*this , true , 'r');
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet VariableOrderMarkovData.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& VariableOrderMarkovData::ascii_write(ostream &os , bool exhaustive) const

{
  if (markov) {
    markov->ascii_write(os , this , exhaustive , false ,
                        ::test_hidden(markov->nb_output_process , markov->nonparametric_process));
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet VariableOrderMarkovData dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool VariableOrderMarkovData::ascii_write(StatError &error , const char *path ,
                                          bool exhaustive) const

{
  bool status = false;


  if (markov) {
    ofstream out_file(path);

    error.init();

    if (!out_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_NAME]);
    }

    else {
      status = true;
      markov->ascii_write(out_file , this , exhaustive , true ,
                          ::test_hidden(markov->nb_output_process , markov->nonparametric_process));
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet VariableOrderMarkovData.
 *
 *  arguments : stream, format (ligne/colonne), flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& VariableOrderMarkovData::ascii_data_write(ostream &os , char format ,
                                                   bool exhaustive) const

{
  MarkovianSequences::ascii_write(os , exhaustive , false);
  ascii_print(os , format , false , posterior_probability);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet VariableOrderMarkovData dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              format (ligne/colonne), flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool VariableOrderMarkovData::ascii_data_write(StatError &error , const char *path ,
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
    ascii_print(out_file , format , true , posterior_probability);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet VariableOrderMarkovData dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool VariableOrderMarkovData::spreadsheet_write(StatError &error , const char *path) const

{
  bool status = false;


  if (markov) {
    ofstream out_file(path);

    error.init();

    if (!out_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_NAME]);
    }

    else {
      status = true;
      markov->spreadsheet_write(out_file , this ,
                                ::test_hidden(markov->nb_output_process , markov->nonparametric_process));
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet VariableOrderMarkovData.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool VariableOrderMarkovData::plot_write(StatError &error , const char *prefix ,
                                         const char *title) const

{
  bool status = false;


  if (markov) {
    status = markov->plot_write(prefix , title , this);

    error.init();

    if (!status) {
      error.update(STAT_error[STATR_FILE_PREFIX]);
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet VariableOrderMarkovData.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* VariableOrderMarkovData::get_plotable() const

{
  MultiPlotSet *plot_set;


  if (markov) {
    plot_set = markov->get_plotable(this);
  }
  else {
    plot_set = NULL;
  }

  return plot_set;
}
