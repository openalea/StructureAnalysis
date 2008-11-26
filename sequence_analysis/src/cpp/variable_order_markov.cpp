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
 *  Constructeur par defaut de la classe Variable_order_markov.
 *
 *--------------------------------------------------------------*/

Variable_order_markov::Variable_order_markov()

{
  nb_iterator = 0;
  markov_data = 0;

  max_order = 0;

  memory_type = 0;
  order = 0;
  state = 0;

  parent = 0;
  child = 0;

  next = 0;
  nb_memory = 0;
  previous = 0;

  nb_output_process = 0;
  nonparametric_process = 0;
  parametric_process = 0;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Variable_order_markov.
 *
 *  arguments : type, nombre d'etats, nombres de memoires.
 *
 *--------------------------------------------------------------*/

Variable_order_markov::Variable_order_markov(char itype , int inb_state , int inb_row)
:Chain(itype , inb_state , inb_row , true)

{
  register int i;


  nb_iterator = 0;
  markov_data = 0;

  max_order = 0;

  memory_type = new int[nb_row];
  order = new int[nb_row];
  state = new int*[nb_row];
  parent = new int[nb_row];
  child = new int*[nb_row];

  for (i = 0;i < nb_row;i++) {
    state[i] = 0;
    child[i] = 0;
  }

  next = 0;
  nb_memory = 0;
  previous = 0;

  nb_output_process = 0;
  nonparametric_process = new Nonparametric_sequence_process*[1];
  nonparametric_process[0] = new Nonparametric_sequence_process(nb_state , nb_state);

  parametric_process = 0;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Variable_order_markov.
 *
 *  arguments : type, nombre d'etats, nombres de memoires, ordre maximum.
 *
 *--------------------------------------------------------------*/

Variable_order_markov::Variable_order_markov(char itype , int inb_state ,
                                             int inb_row , int imax_order)
:Chain(itype , inb_state , inb_row , true)

{
  register int i;


  nb_iterator = 0;
  markov_data = 0;

  max_order = imax_order;

  memory_type = new int[nb_row];
  order = new int[nb_row];
  state = new int*[nb_row];
  parent = new int[nb_row];
  child = new int*[nb_row];

  for (i = 0;i < nb_row;i++) {
    state[i] = new int[max_order];
    child[i] = 0;
  }

  next = 0;
  nb_memory = 0;
  previous = 0;

  nb_output_process = 0;
  nonparametric_process = new Nonparametric_sequence_process*[1];
  nonparametric_process[0] = new Nonparametric_sequence_process(nb_state , nb_state);

  parametric_process = 0;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Variable_order_markov d'ordre fixe.
 *
 *  arguments : type, nombre d'etats, ordre, flag initialisation,
 *              nombre de processus d'observation,
 *              nombre de valeurs observees par processus.
 *
 *--------------------------------------------------------------*/

Variable_order_markov::Variable_order_markov(char itype , int inb_state ,
                                             int iorder , bool init_flag ,
                                             int inb_output_process , int nb_value)
:Chain(itype , inb_state , (int)(pow((double)inb_state , iorder + 1) - 1) / (inb_state - 1) , init_flag)

{
  register int i , j;


  nb_iterator = 0;
  markov_data = 0;

  max_order = iorder;

  memory_type = new int[nb_row];
  order = new int[nb_row];
  state = new int*[nb_row];
  parent = new int[nb_row];
  child = new int*[nb_row];

  // racine (ordre 0)

  memory_type[0] = NON_TERMINAL;
  order[0] = 0;
  state[0] = 0;
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
      child[i] = 0;
    }
  }

  // calcul des transitions entre memoires terminales

  next = 0;
  nb_memory = 0;
  previous = 0;

  build_memory_transition();

  nb_output_process = inb_output_process;

  nonparametric_process = new Nonparametric_sequence_process*[nb_output_process + 1];
  nonparametric_process[0] = new Nonparametric_sequence_process(nb_state , nb_state , false);
  if (nb_output_process == 1) {
    nonparametric_process[1] = new Nonparametric_sequence_process(nb_state , nb_value , true);
  }

  parametric_process = 0;
}


/*--------------------------------------------------------------*
 *
 *  Completion de l'arborescence des memoires.
 *
 *  argument : reference sur un objet Variable_order_markov.
 *
 *--------------------------------------------------------------*/

void Variable_order_markov::memory_tree_completion(const Variable_order_markov &markov)

{
  bool prefix;
  register int i , j , k , m;
  int bnb_memory , border , *markov_next , *completion_next;
  Variable_order_markov *completion;


  bnb_memory = markov.nb_row;
  for (i = 1;i < markov.nb_row;i++) {
    if (markov.order[i] == markov.max_order) {
      bnb_memory--;
    }
  }

  completion = new Variable_order_markov(markov.type , markov.nb_state ,
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
    accessibility = 0;
    nb_component = 0;
    component_nb_state = 0;
    component = 0;
    state_type = 0;
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

  cumul_initial = 0;
  cumul_transition = 0;

  nb_iterator = 0;

  if (markov.markov_data) {
    markov_data = new Variable_order_markov_data(*(markov.markov_data) , false);
  }
  else {
    markov_data = 0;
  }

  max_order = markov.max_order;

  memory_type = new int[nb_row];
  order = new int[nb_row];
  state = new int*[nb_row];
  parent = new int[nb_row];
  child = new int*[nb_row];

  next = 0;
  nb_memory = 0;
  previous = 0;

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
      child[i] = 0;
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
          child[i] = 0;
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
 *  Constructeur de la classe Variable_order_markov avec completion
 *  de l'arborescence des memoires.
 *
 *  arguments : reference sur un objet Variable_order_markov,
 *              nombre de processus d'observation,
 *              nombre de valeurs observees par processus.
 *
 *--------------------------------------------------------------*/

Variable_order_markov::Variable_order_markov(const Variable_order_markov &markov ,
                                             int inb_output_process , int nb_value)

{
  register int i;
  int nb_terminal;


  memory_tree_completion(markov);

  if (type == 'e') {
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
  }

  nb_output_process = inb_output_process;

  nonparametric_process = new Nonparametric_sequence_process*[nb_output_process + 1];
  nonparametric_process[0] = new Nonparametric_sequence_process(nb_state , nb_state , false);
  if (nb_output_process == 1) {
    nonparametric_process[1] = new Nonparametric_sequence_process(nb_state , nb_value , true);
  }

  parametric_process = 0;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Variable_order_markov avec completion
 *  de l'arborescence des memoires.
 *
 *  arguments : reference sur un objet Variable_order_markov,
 *              nombre de processus d'observation,
 *              nombre de valeurs observees par processus.
 *
 *--------------------------------------------------------------*/

/* Variable_order_markov::Variable_order_markov(const Variable_order_markov &markov ,
                                             int inb_output_process , int *nb_value)

{
  register int i;
  int nb_terminal;


  memory_tree_completion(markov);

  if (type == 'e') {
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
  }

  nb_output_process = inb_output_process;

  nonparametric_process = new Nonparametric_sequence_process*[nb_output_process + 1];
  nonparametric_process[0] = new Nonparametric_sequence_process(nb_state , nb_state , false);
  parametric_process = new Parametric_process*[nb_output_process + 1];
  parametric_process[0] = 0;

  for (i = 1;i <= nb_output_process;i++) {
    if (*nb_value <= NB_OUTPUT) {
      nonparametric_process[i] = new Nonparametric_sequence_process(nb_state , *nb_value++ , true);
      parametric_process[i] = 0;
    }
    else {
      nonparametric_process[i] = 0;
      parametric_process[i] = new Parametric_process(nb_state , (int)(*nb_value++ * SAMPLE_NB_VALUE_COEFF));
    }
  }
} */


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Variable_order_markov avec completion
 *  de l'arborescence des memoires.
 *
 *  arguments : pointeur sur un objet Variable_order_markov,
 *              pointeur sur un objet Nonparametric_sequence_process,
 *              longueur des sequences.
 *
 *--------------------------------------------------------------*/

Variable_order_markov::Variable_order_markov(const Variable_order_markov *pmarkov ,
                                             const Nonparametric_process *pobservation , int length)

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
  nonparametric_process = new Nonparametric_sequence_process*[nb_output_process + 1];
  nonparametric_process[0] = new Nonparametric_sequence_process(nb_state , nb_state);
  if (pobservation) {
    nonparametric_process[1] = new Nonparametric_sequence_process(*pobservation);
  }

  parametric_process = 0;

  characteristic_computation(length , true);
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Variable_order_markov.
 *
 *  arguments : reference sur un objet Variable_order_markov.
 *              flag copie de l'objet Variable_order_markov_data.
 *
 *--------------------------------------------------------------*/

void Variable_order_markov::copy(const Variable_order_markov &markov , bool data_flag)

{
  register int i , j;


  nb_iterator = 0;

  if ((data_flag) && (markov.markov_data)) {
    markov_data = new Variable_order_markov_data(*(markov.markov_data) , false);
  }
  else {
    markov_data = 0;
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
      child[i] = 0;
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
        next[i] = 0;
      }
    }
  }
  else {
    next = 0;
  }

  if (markov.nb_memory) {
    nb_memory = new int[nb_row];
    for (i = 0;i < nb_row;i++) {
      nb_memory[i] = markov.nb_memory[i];
    }
  }
  else {
    nb_memory = 0;
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
        previous[i] = 0;
      }
    }
  }
  else {
    previous = 0;
  }

  nb_output_process = markov.nb_output_process;

  nonparametric_process = new Nonparametric_sequence_process*[nb_output_process + 1];

  if (markov.parametric_process) {
    parametric_process = new Parametric_process*[nb_output_process + 1];
    parametric_process[0] = 0;
  }
  else {
    parametric_process = 0;
  }

  nonparametric_process[0] = new Nonparametric_sequence_process(*(markov.nonparametric_process[0]));

  for (i = 1;i <= nb_output_process;i++) {
    if (markov.nonparametric_process[i]) {
      nonparametric_process[i] = new Nonparametric_sequence_process(*(markov.nonparametric_process[i]));
      if (markov.parametric_process) {
        parametric_process[i] = 0;
      }
    }
    else {
      nonparametric_process[i] = 0;
      parametric_process[i] = new Parametric_process(*(markov.parametric_process[i]));
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet Variable_order_markov.
 *
 *--------------------------------------------------------------*/

void Variable_order_markov::remove()

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
 *  Destructeur de la classe Variable_order_markov.
 *
 *--------------------------------------------------------------*/

Variable_order_markov::~Variable_order_markov()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Variable_order_markov.
 *
 *  argument : reference sur un objet Variable_order_markov.
 *
 *--------------------------------------------------------------*/

Variable_order_markov& Variable_order_markov::operator=(const Variable_order_markov &markov)

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

void Variable_order_markov::find_parent_memory(int index)

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

void Variable_order_markov::build_memory_transition()

{
  if (!next) {
    register int i , j , k;
    int bnb_memory;


#   ifdef DEBUG
    cout << "\n";
#   endif

    next = new int*[nb_row];
    next[0] = 0;

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
        next[i] = 0;
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Extraction des memoires precedentes.
 *
 *--------------------------------------------------------------*/

void Variable_order_markov::build_previous_memory()

{
  if ((next) && (!nb_memory) && (!previous)) {
    register int i , j;
    int *buffer;


    nb_memory = new int[nb_row];
    previous = new int*[nb_row];
    nb_memory[0] = 0;
    previous[0] = 0;

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
        previous[i] = 0;
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

bool Variable_order_markov::check_free_suffix() const

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

bool** Variable_order_markov::logic_transition_computation() const

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

void Variable_order_markov::component_computation()

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

void Variable_order_markov::build_non_terminal()

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
 *  arguments : reference sur un objet Format_error, type de loi,
 *              variable, etat ou observation.
 *
 *--------------------------------------------------------------*/

Parametric_model* Variable_order_markov::extract(Format_error &error , int type ,
                                                 int variable , int value) const

{
  bool status = true;
  int hvariable;
  Distribution *pdist;
  Parametric *pparam;
  Parametric_model *dist;
  Histogram *phisto;


  dist = 0;
  error.init();

  pdist = 0;
  pparam = 0;

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
    phisto = 0;

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
      dist = new Parametric_model(*pdist , phisto);
    }
    else if (pparam) {
      dist = new Parametric_model(*pparam , phisto);
    }
  }

  return dist;
}


/*--------------------------------------------------------------*
 *
 *  Extraction de la partie "donnees" d'un objet Variable_order_markov.
 *
 *  argument : reference sur un objet Format_error.
 *
 *--------------------------------------------------------------*/

Variable_order_markov_data* Variable_order_markov::extract_data(Format_error &error) const

{
  bool status = true;
  Variable_order_markov_data *seq;


  seq = 0;
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
    seq = new Variable_order_markov_data(*markov_data);
    seq->markov = new Variable_order_markov(*this , false);
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

void Variable_order_markov::threshold_application(double min_probability)

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
  accessibility = 0;

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
  next = 0;

  delete [] nb_memory;
  nb_memory = 0;

  if (previous) {
    for (i = 1;i < nb_row;i++) {
      delete [] previous[i];
    }
    delete [] previous;
  }
  previous = 0;

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

Variable_order_markov* Variable_order_markov::thresholding(double min_probability) const

{
  register int i;
  Variable_order_markov *markov;


  markov = new Variable_order_markov(*this , false);
  markov->threshold_application(min_probability);

  for (i = 1;i <= nb_output_process;i++) {
    if (nonparametric_process[i]) {
      nonparametric_process[i]->thresholding(min_probability);
    }
  }

  return markov;
}


/*--------------------------------------------------------------*
 *
 *  Analyse du format d'un objet Variable_order_markov.
 *
 *  arguments : reference sur un objet Format_error, stream,
 *              reference sur l'indice de la ligne lue, type du processus
 *              ('o' : ordinaire, 'e' : en equilibre).
 *
 *--------------------------------------------------------------*/

Variable_order_markov* variable_order_markov_parsing(Format_error &error , ifstream &in_file ,
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
  Variable_order_markov *markov;


  markov = 0;

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
      markov = new Variable_order_markov(type , nb_state , nb_non_terminal + nb_terminal , max_order);

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
        markov = 0;
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
 *  Construction d'un objet Variable_order_markov a partir d'un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              longueur des sequences.
 *
 *--------------------------------------------------------------*/

Variable_order_markov* variable_order_markov_ascii_read(Format_error &error ,
                                                        const char *path , int length)

{
  RWCString buffer , token;
  size_t position;
  char type = 'v';
  bool status;
  register int i;
  int line;
  const Variable_order_markov *imarkov;
  const Nonparametric_process *observation;
  Variable_order_markov *markov;
  ifstream in_file(path);


  markov = 0;
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

        observation = 0;

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
          markov = new Variable_order_markov(imarkov , observation , length);

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
 *  Ecriture sur une ligne d'un objet Variable_order_markov.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& Variable_order_markov::line_write(ostream &os) const

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

ostream& Variable_order_markov::ascii_memory_tree_print(ostream &os , bool file_flag) const

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

ostream& Variable_order_markov::ascii_transition_tree_print(ostream &os , bool file_flag) const

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

ostream& Variable_order_markov::ascii_print(ostream &os , bool file_flag) const

{
  register int i , j , k;
  int buff , width;
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
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_STATIONARY_PROBABILITIES] << endl;

    for (i = 1;i < nb_row;i++) {
      if (memory_type[i] == TERMINAL) {
        if (file_flag) {
          os << "# ";
        }
        os << setw(width) << initial[i];

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
 *  Ecriture d'un objet Variable_order_markov et de la structure de donnees associee.
 *
 *  arguments : stream, pointeur sur un objet Variable_order_markov_data,
 *              flag niveau de detail, flag fichier, flag Markov cache.
 *
 *--------------------------------------------------------------*/

ostream& Variable_order_markov::ascii_write(ostream &os , const Variable_order_markov_data *seq ,
                                            bool exhaustive , bool file_flag , bool hidden) const

{
  register int i , j;
  int buff , variable , max_memory_count , *memory_count , width[2];
  double standard_normal_value , half_confidence_interval;
  long old_adjust;
  Histogram **observation = 0;
  Sequence_characteristics *characteristics;


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
    characteristics = 0;
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
        characteristics = 0;
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


    // ecriture de l'histogramme des longueurs des sequences

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
    seq->hlength->ascii_characteristic_print(os , false , file_flag);

    if (exhaustive) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << "   | " << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << endl;
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
 *  Ecriture d'un objet Variable_order_markov.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Variable_order_markov::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , markov_data , exhaustive , false);
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Variable_order_markov dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Variable_order_markov::ascii_write(Format_error &error , const char *path ,
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

ostream& Variable_order_markov::spreadsheet_print(ostream &os) const

{
  register int i , j , k;


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
    os << "\n" << STAT_label[STATL_STATIONARY_PROBABILITIES] << endl;
    for (i = 1;i < nb_row;i++) {
      if (memory_type[i] == TERMINAL) {
        os << initial[i] << "\t\t";
        for (j = order[i] - 1;j >= 0;j--) {
          os << state[i][j] << " ";
        }
        os << endl;
      }
    }
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
 *  Ecriture d'un objet Variable_order_markov et de la structure de donnees associee
 *  dans un fichier au format tableur.
 *
 *  arguments : stream, pointeur sur un objet Variable_order_markov_data,
 *              flag Markov cache.
 *
 *--------------------------------------------------------------*/

ostream& Variable_order_markov::spreadsheet_write(ostream &os ,
                                                  const Variable_order_markov_data *seq ,
                                                  bool hidden) const

{
  register int i;
  int variable;
  Histogram **observation = 0;
  Sequence_characteristics *characteristics;


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
    characteristics = 0;
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
        characteristics = 0;
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


    // ecriture de l'histogramme des longueurs des sequences

    os << "\n" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << "\t";
    seq->hlength->spreadsheet_characteristic_print(os);

    os << "\n\t" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << endl;
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
 *  Ecriture d'un objet Variable_order_markov dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

bool Variable_order_markov::spreadsheet_write(Format_error &error , const char *path) const

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
 *  Sortie Gnuplot d'un objet Variable_order_markov et de la structure de donnees associee.
 *
 *  arguments : prefixe des fichiers, titre des figures,
 *              pointeur sur les sequences observees.
 *
 *--------------------------------------------------------------*/

bool Variable_order_markov::plot_write(const char *prefix , const char *title ,
                                       const Variable_order_markov_data *seq) const

{
  bool status;
  register int i;
  int variable;
  Histogram *hlength = 0 , **observation = 0;
  Sequence_characteristics *characteristics;


  if ((seq) && (seq->type[0] == STATE)) {
    characteristics = seq->characteristics[0];
    hlength = seq->hlength;
  }
  else {
    characteristics = 0;
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
          characteristics = 0;
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
 *  Sortie Gnuplot d'un objet Variable_order_markov.
 *
 *  arguments : reference sur un objet Format_error, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Variable_order_markov::plot_write(Format_error &error , const char *prefix ,
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
 *  Calcul de l'ordre maximum des memoires.
 *
 *--------------------------------------------------------------*/

void Variable_order_markov::max_order_computation()

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

int Variable_order_markov::nb_parameter_computation(double min_probability) const

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

int Variable_order_markov::nb_transient_parameter_computation(double min_probability) const

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

double Variable_order_markov::penalty_computation(bool hidden , double min_probability) const

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
              penalty += nb_parameter * log((double)markov_data->marginal[0]->frequency[j]);
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
            penalty += nb_parameter * log((double)markov_data->marginal[0]->frequency[j]);
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
 *  Constructeur par defaut de la classe Variable_order_markov_data.
 *
 *--------------------------------------------------------------*/

Variable_order_markov_data::Variable_order_markov_data()

{
  markov = 0;
  chain_data = 0;

  likelihood = D_INF;
  hidden_likelihood = D_INF;

  posterior_probability = 0;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Variable_order_markov_data.
 *
 *  arguments : histogramme des longueurs des sequences, nombre de variables,
 *              flag initialisation.
 *
 *--------------------------------------------------------------*/

Variable_order_markov_data::Variable_order_markov_data(const Histogram &ihlength , int inb_variable ,
                                                       bool init_flag)
:Markovian_sequences(ihlength , inb_variable , init_flag)

{
  markov = 0;
  chain_data = 0;

  likelihood = D_INF;
  hidden_likelihood = D_INF;

  posterior_probability = 0;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Variable_order_markov_data a partir
 *  d'un objet Markovian_sequences avec ajout d'une variable d'etat.
 *
 *  argument : reference sur un objet Markovian_sequences.
 *
 *--------------------------------------------------------------*/

Variable_order_markov_data::Variable_order_markov_data(const Markovian_sequences &seq)
:Markovian_sequences(seq , 'a' , DEFAULT)

{
  markov = 0;
  chain_data = 0;

  likelihood = D_INF;
  hidden_likelihood = D_INF;

  posterior_probability = 0;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Variable_order_markov_data a partir
 *  d'un objet Markovian_sequences.
 *
 *  arguments : reference sur un objet Markovian_sequences, type de transformation
 *              ('c' : copie, 'a' : ajout d'une variable d'etat),
 *              ajout/suppression des histogrammes de temps de sejour initial.
 *
 *--------------------------------------------------------------*/

Variable_order_markov_data::Variable_order_markov_data(const Markovian_sequences &seq ,
                                                       char transform , bool initial_run_flag)
:Markovian_sequences(seq , transform , (initial_run_flag ? ADD_INITIAL_RUN : REMOVE_INITIAL_RUN))

{
  markov = 0;
  chain_data = 0;

  likelihood = D_INF;
  hidden_likelihood = D_INF;

  posterior_probability = 0;
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Variable_order_markov_data.
 *
 *  arguments : reference sur un objet Variable_order_markov_data,
 *              flag copie de l'objet Variable_order_markov.
 *
 *--------------------------------------------------------------*/

void Variable_order_markov_data::copy(const Variable_order_markov_data &seq ,
                                      bool model_flag)

{
  register int i;


  if ((model_flag) && (seq.markov)) {
    markov = new Variable_order_markov(*(seq.markov) , false);
  }
  else {
    markov = 0;
  }

  if (seq.chain_data) {
    chain_data = new Variable_order_chain_data(*(seq.chain_data));
  }
  else {
    chain_data = 0;
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
    posterior_probability = 0;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Variable_order_markov_data.
 *
 *--------------------------------------------------------------*/

Variable_order_markov_data::~Variable_order_markov_data()

{
  delete markov;
  delete chain_data;

  delete [] posterior_probability;
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Variable_order_markov_data.
 *
 *  argument : reference sur un objet Variable_order_markov_data.
 *
 *--------------------------------------------------------------*/

Variable_order_markov_data& Variable_order_markov_data::operator=(const Variable_order_markov_data &seq)

{
  if (&seq != this) {
    delete markov;
    delete chain_data;

    delete [] posterior_probability;

    remove();
    Sequences::remove();

    Sequences::copy(seq);
    Markovian_sequences::copy(seq);
    copy(seq);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Extraction d'un histogramme.
 *
 *  arguments : reference sur un objet Format_error, type d'histogramme,
 *              variable, etat ou observation.
 *
 *--------------------------------------------------------------*/

Distribution_data* Variable_order_markov_data::extract(Format_error &error , int type ,
                                                       int variable , int value) const

{
  bool status = true;
  Distribution *pdist;
  Parametric *pparam;
  Histogram *phisto;
  Distribution_data *histo;


  histo = 0;
  error.init();

  if (type == OBSERVATION) {
    if ((variable < 2) || (variable > nb_variable)) {
      status = false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
    }

    else {
      variable--;

      if ((value < 0) || (value >= marginal[0]->nb_value)) {
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
          error.update(STAT_error[STATR_EMPTY_HISTOGRAM]);
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

      else if ((value < 0) || (value >= marginal[variable]->nb_value)) {
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
          error.update(STAT_error[STATR_EMPTY_HISTOGRAM]);
        }
      }
    }
  }

  if (status) {
    pdist = 0;
    pparam = 0;

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
      histo = new Distribution_data(*phisto , pdist);
    }
    else {
      histo = new Distribution_data(*phisto , pparam);
    }
  }

  return histo;
}


/*--------------------------------------------------------------*
 *
 *  Suppression du parametre d'index.
 *
 *  argument : reference sur un objet Format_error.
 *
 *--------------------------------------------------------------*/

Variable_order_markov_data* Variable_order_markov_data::remove_index_parameter(Format_error &error) const

{
  Variable_order_markov_data *seq;


  error.init();

  if (!index_parameter) {
    seq = 0;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }
  else {
    seq = new Variable_order_markov_data(*this , true , 'r');
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Variable_order_markov_data.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Variable_order_markov_data::ascii_write(ostream &os , bool exhaustive) const

{
  if (markov) {
    markov->ascii_write(os , this , exhaustive , false ,
                        ::test_hidden(markov->nb_output_process , markov->nonparametric_process));
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Variable_order_markov_data dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Variable_order_markov_data::ascii_write(Format_error &error , const char *path ,
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
 *  Ecriture d'un objet Variable_order_markov_data.
 *
 *  arguments : stream, format (ligne/colonne), flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Variable_order_markov_data::ascii_data_write(ostream &os , char format ,
                                                      bool exhaustive) const

{
  Markovian_sequences::ascii_write(os , exhaustive , false);
  ascii_print(os , format , false , posterior_probability);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Variable_order_markov_data dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              format (ligne/colonne), flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Variable_order_markov_data::ascii_data_write(Format_error &error , const char *path ,
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
      Markovian_sequences::ascii_write(out_file , exhaustive , true);
    }
    ascii_print(out_file , format , true , posterior_probability);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Variable_order_markov_data dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

bool Variable_order_markov_data::spreadsheet_write(Format_error &error , const char *path) const

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
 *  Sortie Gnuplot d'un objet Variable_order_markov_data.
 *
 *  arguments : reference sur un objet Format_error, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Variable_order_markov_data::plot_write(Format_error &error , const char *prefix ,
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
