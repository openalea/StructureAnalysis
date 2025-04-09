/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *
 *       Copyright 2005-2009 UMR DAP
 *
 *       File author(s): F. Chaubert-Pereira (chaubert@cirad.fr)
 *
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
#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"
#include "stat_tool/stat_tools.h"
#include "stat_tool/stat_label.h"
#include "switching_process.h"

using namespace std;


/*---------------------------------------------------------------
 *
 * Constructeur de la classe Switching_process
 *
 * arguments : nombre d'états, nombre de valeur du processus 
 *             d'observation 
 *
 *---------------------------------------------------------------*/

Switching_process::Switching_process(int inb_state, int inb_value)
{
  register int i;

  nb_state = inb_state;
  nb_value = inb_value;

  if (nb_state > 0) {
    observation = new Continuous_parametric*[nb_state];
    if (nb_value > 0) {
      for (i = 0; i < nb_state; i++) {
	observation[i] = new Continuous_parametric (0,0, 1.,0,0);
      }
    }
    else {
      for (i = 0; i < nb_state; i++) {
	observation[i] = 0;
      }
    }
  }
  else {
    observation = 0;
  }
}


/*-----------------------------------------------------------------
 *
 * Constructeur de la classe Switching_process
 *
 * arguments : nombre d'états, processus d'observation (ensemble de 
 *             modèles linéaires mixtes) 
 *
 *-----------------------------------------------------------------*/

Switching_process::Switching_process(int inb_state, Continuous_parametric **pobservation)
{
  register int i;

  nb_state = inb_state;

  nb_value = 0;
  for (i = 0; i < nb_state; i++) {
    if (pobservation[i]->nb_value > nb_value ) {
      nb_value = pobservation[i]->nb_value;
    }
  }

  observation = new Continuous_parametric*[nb_state];
  for (i = 0; i < nb_state; i++) {
    observation[i] = new Continuous_parametric(*pobservation[i]);
  }
}


/*------------------------------------------------------------------
 *
 * Copie d'un objet Switching_process
 *
 * arguments : réference sur un objet Switching_process
 *
 *------------------------------------------------------------------*/

void Switching_process::copy(const Switching_process  &process)
{
  register int i;

  nb_state = process.nb_state;
  nb_value = process.nb_value;

  observation = new Continuous_parametric*[nb_state];
  for (i = 0;i < nb_state; i++) {
    observation[i] = new Continuous_parametric(*(process.observation[i]));
  }
}


/*------------------------------------------------------------------
 *
 * Copie d'un objet Switching_process avec ajout d'un état
 *
 * arguments : reference sur un objet Switching_process, 
 *             numéro de l'état
 *
 *------------------------------------------------------------------*/

void Switching_process::add_state(const Switching_process &process, int state)
{
  register int i;

  nb_state = process.nb_state + 1;
  nb_value = process.nb_value;

  observation = new Continuous_parametric*[nb_state];
  for (i = 0; i < nb_state; i++) {
    observation[i] = new Continuous_parametric(*(process.observation[i]));
  }
  for (i = state+1; i < nb_state; i++) {
    observation[i] = new Continuous_parametric(*(process.observation[i-1]));
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par copie de la classe Switching_process.
 *
 *  arguments : reference sur un objet Nonparametric_sequence_process,
 *              type de manipulation ('c' : copy, 's' : state, 'o' :occupancy),
 *              numéro de l'état que l'on peut ajouter
 *
 *--------------------------------------------------------------*/

Switching_process::Switching_process(const Switching_process &process, char manip, int state)
{
  switch(manip) {
  case 'c' : 
    copy(process);
    break;
  case 's' :
    add_state(process, state);
    break;
  }
}


/*-----------------------------------------------------------------
 *
 * Destruction des champs d'un objet Switching_process
 *
 *-----------------------------------------------------------------*/

void Switching_process::remove()
{
  if(observation) {
    register int i;

    for( i = 0;i < nb_state;i++) {
      delete observation[i];
    }
    delete [] observation;
  }
}


/*-----------------------------------------------------------------
 *
 * Destructeur de la classe Switching_process
 *
 *-----------------------------------------------------------------*/

Switching_process::~Switching_process()
{
  remove();
}


/*-----------------------------------------------------------------
 *
 * Opérateur d'assignement de la classe Switching_process
 *
 * arguments : reference sur un objet Switching_process
 *
 *-----------------------------------------------------------------*/

Switching_process& Switching_process::operator=(const Switching_process &process)
{
  if (&process != this) {
    remove();
    copy(process);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Switching_process.
 *
 *  arguments : stream, 
 *              pointeurs sur les lois d'observation empiriques,
 *              flag niveau de detail, flag fichier
 *
 *--------------------------------------------------------------*/

ostream& Switching_process::ascii_print(ostream &os, Continuous_histo **empirical_observation,
					bool exhaustive, bool file_flag) const
{
  register int i; 

  for (i = 0; i < nb_state; i++) {
    os << "\n" << STAT_word[STATW_STATE] << " " << i << " "
       << STAT_word[STATW_OBSERVATION_DISTRIBUTION] << endl;
    observation[i]->ascii_print(os);
    //   observation[i]->ascii_characteristic_print(os , false );
    
    if (empirical_observation) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_STATE] << " " << i << " " << STAT_label[STATL_OBSERVATION] << " "
         << STAT_label[STATL_HISTOGRAM] << " - ";
      empirical_observation[i]->ascii_characteristic_print(os , false);
      
      if (exhaustive) {
	os << "\n";
	if (file_flag) {
	  os << "# ";
	}
	os << "class  " << "|  frequency  ";
      }
      if (empirical_observation) {
	os << " |  " 
	   << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_HISTOGRAM];
	os << " | " 
	   << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
	if (empirical_observation) {
	  os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM] << " ";
	}
	os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
	   << endl;
      }

      observation[i]->Continuous_distribution::ascii_print(os , true , true , true , true,
							   (empirical_observation ? empirical_observation[i] : 0));
    }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Switching_process au format tableur.
 *
 *  arguments : stream, pointeurs sur les lois d'observation 
 *              empiriques
 *
 *--------------------------------------------------------------*/

ostream& Switching_process::spreadsheet_print(ostream &os , Continuous_histo **empirical_observation) const
{
  register int i;

  for (i = 0;i < nb_state; i++) {
    os << "\n" << STAT_word[STATW_STATE] << " " << i << "\t"
       << STAT_word[STATW_OBSERVATION_DISTRIBUTION] << endl;
    observation[i]->spreadsheet_print(os);
    observation[i]->spreadsheet_characteristic_print(os);

    if (empirical_observation ) {
      os << "\n" << STAT_label[STATL_STATE] << " " << i << " "
         << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_HISTOGRAM] << "\t";
      empirical_observation[i]->spreadsheet_characteristic_print(os);
    }

    os << "\n";
    if (empirical_observation ) {
      os << "\t" 
         << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_HISTOGRAM];
      
      os << "\t" 
	 << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
      
      os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM] << " "
         << STAT_label[STATL_FUNCTION];
   
      os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
	 << STAT_label[STATL_FUNCTION] << endl;
    }

    observation[i]->Continuous_distribution::spreadsheet_print(os , true, true, true ,
								 (empirical_observation? empirical_observation[i] : 0) );
  }

  return os;
}


/*---------------------------------------------------------------
 *
 * Calcul du nombre de valeurs du processus d'observation
 *
 *---------------------------------------------------------------*/

void Switching_process::nb_value_computation()
{
  register int i;

  nb_value = 0;
  for ( i = 0; i < nb_state;i++) {
    if (observation[i]->nb_value > nb_value) {
      nb_value = observation[i]->nb_value;
    }
  }
}


/*--------------------------------------------------------------
 *
 * Calcul du nombre de parametres independants d'un processus
 * d'observation continu
 *
 *--------------------------------------------------------------*/

int Switching_process::nb_parameter_computation() const
{
  register int i;
  int nb_parameter = 0;

  for (i = 0;i < nb_state;i++) {
    nb_parameter += observation[i]->nb_parameter_computation();
  }
  
  return nb_parameter;
}


/*---------------------------------------------------------------
 *
 * Initialisation des lois d'observation par perturbation.
 *
 *--------------------------------------------------------------*/

void Switching_process::init()
{
  register int i , j;
  double noise_proba , *pdensity;

  for (i = 0;i < nb_state;i++) {
    noise_proba = SWITCHING_NOISE_PROBABILITY * nb_state / nb_value;

    pdensity = observation[i]->density;
    for (j = 0;j < i * nb_value / nb_state;j++) {
      *pdensity++ -= noise_proba / (nb_state - 1);
    }
    for (j = i * nb_value / nb_state;j < (i + 1) * nb_value / nb_state;j++) {
      *pdensity++ += noise_proba;
    }
    for (j = (i + 1) * nb_value / nb_state;j < nb_value;j++) {
      *pdensity++ -= noise_proba / (nb_state - 1);
    }
  }

  pdensity = NULL;
  delete pdensity;

}
