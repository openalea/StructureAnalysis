#include <math.h>
#include <sstream>
#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"
#include "stat_tool/stat_tools.h"
#include "stat_tool/stat_label.h"
#include "switching_process.h"

using namespace std;

//---------------------------------------------------------------
// Constructeur de la classe Switching_process
//---------------------------------------------------------------
Switching_process::Switching_process(int inb_state, int inb_value)
{
  register int i;

  nb_state = inb_state;
  nb_value = inb_value;

  if (nb_state > 0) {
    observation = new Discrete_parametric*[nb_state];
    if (nb_value > 0) {
      for (i = 0; i < nb_state; i++) {
	observation[i] = new Discrete_parametric ();
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

//-----------------------------------------------------------------
// Constructeur de la classe Switching_process
//-----------------------------------------------------------------
Switching_process::Switching_process(int inb_state, Discrete_parametric **pobservation)
{
  register int i;

  nb_state = inb_state;

  nb_value = 0;
  for (i = 0; i < nb_state; i++) {
    if (pobservation[i]->nb_value > nb_value ) {
      nb_value = pobservation[i]->nb_value;
    }
  }

  observation = new Discrete_parametric*[nb_state];
  for (i = 0; i < nb_state; i++) {
    observation[i] = new Discrete_parametric(*pobservation[i]);
  }
}

//------------------------------------------------------------------
// Copie d'un objet Switching_process
//------------------------------------------------------------------
void Switching_process::copy(const Switching_process  &process)
{
  register int i;

  nb_state = process.nb_state;
  nb_value = process.nb_value;

  observation = new Discrete_parametric*[nb_state];
  for (i = 0;i < nb_state; i++) {
    observation[i] = new Discrete_parametric(*(process.observation[i]));
  }
}

//------------------------------------------------------------------
// Copie d'un objet Switching_process avec ajout d'un état
//------------------------------------------------------------------
void Switching_process::add_state(const Switching_process &process, int state)
{
  register int i;

  nb_state = process.nb_state + 1;
  nb_value = process.nb_value;

  observation = new Discrete_parametric*[nb_state];
  for (i = 0; i < nb_state; i++) {
    observation[i] = new Discrete_parametric(*(process.observation[i]));
  }
  for (i = state+1; i < nb_state; i++) {
    observation[i] = new Discrete_parametric(*(process.observation[i-1]));
  }
}

//-----------------------------------------------------------------
// Constructeur par copie de la classe Switching_process
//-----------------------------------------------------------------
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

//-----------------------------------------------------------------
// Destruction des champs d'un objet Switching_process
//-----------------------------------------------------------------
void Switching_process::remove()
{
  if(observation) {
    register int i;

    for( i = 0;i < nb_state;i++) {
      delete observation[i];
    }
    delete [] observation;
    observation = 0;
  }
}

//-----------------------------------------------------------------
// Destructeur de la classe Switching_process
//-----------------------------------------------------------------
Switching_process::~Switching_process()
{
  remove();
}

//-----------------------------------------------------------------
// Opérateur d'assignement de la classe Switching_process
//-----------------------------------------------------------------
Switching_process& Switching_process::operator=(const Switching_process &process)
{
  if (&process != this) {
    remove();
    copy(process);
  }

  return *this;
}

//-----------------------------------------------------------------
// Ecriture d'un objet Switching_process
//-----------------------------------------------------------------
ostream& Switching_process::ascii_print(ostream &os, Histogram **empirical_observation,
					bool exhaustive, bool file_flag) const
{
  register int i; 

  for (i = 0;i < nb_state;i++) {
    os << "\n" << STAT_word[STATW_STATE] << " " << i << " "
       << STAT_word[STATW_OBSERVATION_DISTRIBUTION] << endl;
    observation[i]->ascii_print(os);
    if(observation[i]->nb_cat == 0){
      observation[i]->Parametric::ascii_print(os);
    }
    observation[i]->ascii_characteristic_print(os , false );
    
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

      observation[i]->Distribution::ascii_print(os , true , true , true,
                                                (empirical_observation ? empirical_observation[i] : 0));
    }
    os<<endl<<endl;
  }

  return os;
}

//---------------------------------------------------------------
//Ecriture d'un objet Switching_process au format tableur.
//---------------------------------------------------------------
ostream& Switching_process::spreadsheet_print(ostream &os , Histogram **empirical_observation) const
{
  register int i;

  for (i = 0;i < nb_state;i++) {
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

    observation[i]->Distribution::spreadsheet_print(os , true, true, true ,
						    (empirical_observation? empirical_observation[i] : 0) );
  }

  return os;
}

//---------------------------------------------------------------
// Calcul du nombre de valeurs du processus d'observation
//---------------------------------------------------------------
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

//--------------------------------------------------------------
// Calcul du nombre de parametres independants d'un processus
// d'observation continu
//--------------------------------------------------------------
int Switching_process::nb_parameter_computation() const
{
  register int i;
  int nb_parameter = 0;

  for (i = 0;i < nb_state;i++) {
    nb_parameter += observation[i]->nb_parameter_computation();
  }
  
  return nb_parameter;
}

//---------------------------------------------------------------
// Initialisation des lois d'observation par perturbation.
//--------------------------------------------------------------
void Switching_process::init()
{
  register int i , j;
  double noise_proba , *pmass;

  for (i = 0;i < nb_state;i++) {
    noise_proba = SWITCHING_NOISE_PROBABILITY * nb_state / nb_value;

    pmass = observation[i]->mass;
    for (j = 0;j < i * nb_value / nb_state;j++) {
      *pmass++ -= noise_proba / (nb_state - 1);
    }
    for (j = i * nb_value / nb_state;j < (i + 1) * nb_value / nb_state;j++) {
      *pmass++ += noise_proba;
    }
    for (j = (i + 1) * nb_value / nb_state;j < nb_value;j++) {
      *pmass++ -= noise_proba / (nb_state - 1);
    }
  }
}
