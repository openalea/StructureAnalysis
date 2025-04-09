#include <math.h>
#include <iostream>
#include <iomanip>
#include <limits.h>
#include "stat_tool/stat_label.h"
#include "tool/util_math.h"
#include "stat_tool/stat_tools.h"
#include "tool/config.h"
#include "stat_tool/curves.h"
#include "sequence_analysis/renewal.h"
#include "stat_tool/markovian.h"
#include "sequence_analysis/sequences.h"
#include "sequence_analysis/sequence_label.h"
#include "stat_tool/vectors.h"
#include "sequence_analysis/tops.h"
#include "switching_sequence.h"


using namespace std;

//--------------------------------------------------------------
// Constructeur de la classe Switching_sequence
//--------------------------------------------------------------
Switching_sequence::Switching_sequence()
{
  constant = 0;
  nb_covariable = 0;
  observation = 0;
  covar = 0;
  self_transition = 0;
  characteristics = 0;
}

//-------------------------------------------------------------
// Initialisation par défaut des champs de la classe Switching_sequence
//-------------------------------------------------------------
void Switching_sequence::init()
{  
  register int i;

  self_transition = 0;
  observation = 0;

  characteristics = new Sequence_characteristics*[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    characteristics[i] = 0;
  }
}

//-------------------------------------------------------------
// Constructeur de la classe Switching_sequence
//-------------------------------------------------------------
Switching_sequence::Switching_sequence(int inb_variable, int *itype, int inb_sequence, int *iidentifier, int *ilength, int inb_covariable,
				       int iconstant, bool init_flag)
  : Sequences(inb_sequence, iidentifier, ilength, inb_variable, itype, init_flag)
//:Sequences( inb_variable, itype, inb_sequence, iidentifier, ilength, init_flag)
{
  register int i, j;

  init();


  constant = iconstant;
  if (!constant){
    nb_covariable = inb_covariable + 1;
  }
  else {
    nb_covariable = inb_covariable;
  }

  covar = new double**[nb_sequence];
  for (i = 0; i < nb_sequence; i++) {
    covar[i] = new double*[length[i]];
    for (j = 0; j < length[i]; j++) {
      covar[i][j] = 0;
    }
  }

  // build_characteristic();
}

//-------------------------------------------------------------
// Constructeur de la classe Switching_sequence
//-------------------------------------------------------------
Switching_sequence:: Switching_sequence(int inb_variable, int inb_sequence, int *iidentifier, int *ilength, 
					int inb_covariable, int iconstant, bool init_flag)
 :Sequences(inb_sequence, iidentifier, ilength, inb_variable, init_flag)
//  :Sequences( inb_variable, inb_sequence, iidentifier, ilength, init_flag)
{
  register int i, j;
  init();
  // build_characteristic();  

  constant = iconstant;
  if (!constant){
    nb_covariable = inb_covariable + 1;
  }
  else {
    nb_covariable = inb_covariable;
  }
  
  covar = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    covar[i] = new double*[length[i]];
    for (j = 0;j < length[i];j++) {
      covar[i][j] = 0;
    }
  }

}

//------------------------------------------------------------
// Constructeur de la classe Switching_sequence
//------------------------------------------------------------
Switching_sequence::Switching_sequence(int inb_variable, const Histogram &ihlength, int inb_covariable, 
				       int iconstant, bool init_flag)
  : Sequences( ihlength, inb_variable, init_flag)
// :Sequences( inb_variable, ihlength, init_flag)
{
  register int i, j, k;
  double *pcovar;

  init();
  // build_characteristic();  

  constant = iconstant;
  if (!constant){
    nb_covariable = inb_covariable + 1;
  }
  else {
    nb_covariable = inb_covariable;
  }
  
  covar = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    covar[i] = new double*[length[i]];
    for (j = 0;j < length[i];j++) {
      covar[i][j] = new double[nb_covariable];

      if (init_flag) {
        pcovar = covar[i][j];
        for (k = 0;k < nb_covariable;k++) {
          *pcovar++ = 0;
        }
      }
    }
  }
}


//-------------------------------------------------------------
// Constructeur de la classe Switching_sequence
//-------------------------------------------------------------
Switching_sequence:: Switching_sequence(int inb_variable, int inb_sequence, int inb_covariable, 
					int *ilength, int ***isequence, double ***icovar, int *iidentifier,  int iconstant)
  :Sequences(inb_sequence, iidentifier, ilength, IMPLICIT_TYPE, inb_variable, INT_VALUE, isequence)
//:Sequences(inb_variable , inb_sequence , ilength , isequence)
{
  register int i, j, k;
  double *pcovar , *ccovar;

  init();

  constant = iconstant;

  if (!constant){
    nb_covariable = inb_covariable + 1;
  }
  else {
    nb_covariable = inb_covariable;
  }

  covar = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    covar[i] = new double*[length[i]];
    for (j = 0;j < length[i];j++) {
      covar[i][j] = new double[nb_covariable];
      
      pcovar = covar[i][j];
      ccovar = icovar[i][j];
      if (!constant) {
	*pcovar++ = 1;
      }
      
      if((nb_covariable > 1) || (constant)){// FAUDRA REVOIR
	for (k = 0;k < nb_covariable;k++) {
	  *pcovar++ = *ccovar++;
	  //	cout<<covar[i][j][k];
	}
      }
      //  cout<<endl;
    }
    //   cout<<endl;
   }

  //build_characteristic();
}

//-------------------------------------------------------------
// Constructeur de la classe Switching_sequence
//-------------------------------------------------------------
Switching_sequence::Switching_sequence(const Sequences &seq, int inb_covariable, double ***icovar, int iconstant)
  :Sequences(seq)
{
  register int i, j, k;
  double *pcovar , *ccovar;
  
  init();
  
  constant = iconstant;
  if (!constant){
    nb_covariable = inb_covariable + 1;
  }
  else {
    nb_covariable = inb_covariable;
  }

  covar = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    covar[i] = new double*[length[i]];
    for (j = 0;j < length[i];j++) {
      covar[i][j] = new double[nb_covariable];
      
      pcovar = covar[i][j];
      ccovar = icovar[i][j];
      if (!constant) {
	*pcovar++ = 1;
      }

      if((nb_covariable > 1) || (constant)) { //faudra revoir
	for (k = 0;k < nb_covariable;k++) {
	  *pcovar++ = *ccovar++;
	  //   	cout<<covar[i][j][k];
	}
      }
      //   cout<<endl;
    }
    //   cout<<endl;
  }
 

  build_characteristic();
}

//--------------------------------------------------------------
// Destructeur des champs de la classe Switching_sequence.
//--------------------------------------------------------------
void Switching_sequence::remove()
{
  register int i , j;

  if (self_transition) {
    for (i = 0;i < marginal[0]->nb_value;i++) {
      delete self_transition[i];
    }
    delete [] self_transition;
  }

  if (observation) {
    for (i = 1;i < nb_variable;i++) {
      for (j = 0;j < marginal[0]->nb_value;j++) {
        delete observation[i][j];
      }
      delete [] observation [i];
    }
    delete [] observation;
  }

  if (characteristics) {
    for (i = 0;i < nb_variable;i++) {
      delete characteristics[i];
    }
    delete [] characteristics;
  }

  if (covar) {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        delete [] covar[i][j];
      }
      delete [] covar[i];
    }
    delete [] covar;
  }
}

//--------------------------------------------------------------
// Destructeur de la classe Switching_sequence.
//--------------------------------------------------------------
Switching_sequence::~Switching_sequence()
{
  remove();
}

//-------------------------------------------------------------
// Remplissage de la variable covar PAS IMPLEMENTEE
//-------------------------------------------------------------
void Switching_sequence::covar_computation()
{}

//--------------------------------------------------------------
// Copie d'un objet Switching_sequence.
//--------------------------------------------------------------
void Switching_sequence::copy(const Switching_sequence &sw_seq , int param)
{
  bool initial_run_flag;
  register int i , j, k;
  double *pcovar, *ccovar;
  
  constant = sw_seq.constant;
  nb_covariable = sw_seq.nb_covariable;


  if ((sw_seq.self_transition) && (param != REVERSE)) {
    self_transition = new Self_transition*[marginal[0]->nb_value];
    for (i = 0;i < marginal[0]->nb_value;i++) {
      if (sw_seq.self_transition[i]) {
        self_transition[i] = new Self_transition(*(sw_seq.self_transition[i]));
      }
      else {
        self_transition[i] = 0;
      }
    }
  }
  
  else {
    self_transition = 0;
  }

  if (sw_seq.observation) {
    observation = new Histogram**[nb_variable];
    observation[0] = 0;

    for (i = 1;i < nb_variable;i++) {
      observation[i] = new Histogram*[marginal[0]->nb_value];
      for (j = 0;j < marginal[0]->nb_value;j++) {
        observation[i][j] = new Histogram(*(sw_seq.observation[i][j]));
      }
    }
  }

  else {
    observation = 0;
  }

  if (covar) {
    covar = new double**[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      covar[i] = new double*[sw_seq.length[i]];
      for (j = 0;j < sw_seq.length[i];j++) {
	covar[i][j] = new double[nb_covariable];
	pcovar = covar[i][j];
	
	ccovar = sw_seq.covar[i][j];
	if (sw_seq.covar[i][j]){
	  for (k = 0;k < nb_covariable;k++) {
	    *pcovar++ = *ccovar++;
	  }
	}
	else {
	  pcovar = ccovar;
	}
      }
    }
  }
  else {
    covar = 0;
  }
    
  characteristics = new Sequence_characteristics*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    if (sw_seq.characteristics[i]) {
      if ((param == ADD_INITIAL_RUN) || (param == REMOVE_INITIAL_RUN)) {
        switch (param) {
        case ADD_INITIAL_RUN :
          initial_run_flag = true;
          break;
        case REMOVE_INITIAL_RUN :
          initial_run_flag = false;
          break;
        }

        characteristics[i] = new Sequence_characteristics(*(sw_seq.characteristics[i]) , initial_run_flag);

        if (((sw_seq.characteristics[i]->initial_run) && (!initial_run_flag)) ||
            ((!(sw_seq.characteristics[i]->initial_run)) && (initial_run_flag))) {
	  build_sojourn_time_histogram(i , initial_run_flag);
        }
      }

      else if (param == REVERSE) {
        characteristics[i] = new Sequence_characteristics(*(sw_seq.characteristics[i]), 'r');

        build_index_value(i);
	build_first_occurrence_histogram(i);
	
	if (!(sw_seq.characteristics[i]->initial_run)) {
	  build_sojourn_time_histogram(i);
	}
      }

      else {
        characteristics[i] = new Sequence_characteristics(*(sw_seq.characteristics[i]));
      }
    }

    else {
      characteristics[i] = 0;
    }
  }
}

//--------------------------------------------------------------
// Copie d'un objet Switching_sequence avec ajout d'une variable.
//--------------------------------------------------------------
void Switching_sequence::add_variable(const Switching_sequence &sw_seq ,
                                       int variable , int param)
{
  bool initial_run_flag;
  register int i , j, k;
  double *pcovar, *ccovar;

  self_transition = 0;
  observation = 0;

  characteristics = new Sequence_characteristics*[nb_variable];

  nb_covariable = sw_seq.nb_covariable;
  covar = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    covar[i] = new double*[sw_seq.length[i]];
    for (j = 0;j < sw_seq.length[i];j++) {
      covar[i][j] = new double[nb_covariable];
      pcovar = covar[i][j];

      ccovar = sw_seq.covar[i][j];
      for (k = 0;k < nb_covariable;k++) {
	*pcovar++ = *ccovar++;
      }
    }
  }


  i = 0;
  for (j = 0;j < nb_variable;j++) {
    if (j != variable) {
      if (sw_seq.characteristics[i]) {
        if ((param == ADD_INITIAL_RUN) || (param == REMOVE_INITIAL_RUN)) {
          switch (param) {
          case ADD_INITIAL_RUN :
            initial_run_flag = true;
            break;
          case REMOVE_INITIAL_RUN :
            initial_run_flag = false;
            break;
          }

          characteristics[j] = new Sequence_characteristics(*(sw_seq.characteristics[i]) , initial_run_flag);

          if (((sw_seq.characteristics[i]->initial_run) && (!initial_run_flag)) ||
              ((!(sw_seq.characteristics[i]->initial_run)) && (initial_run_flag))) {
             build_sojourn_time_histogram(j , initial_run_flag);
          }
        }

        else {
          characteristics[j] = new Sequence_characteristics(*(sw_seq.characteristics[i]));
        }
      }
      
      else {
        characteristics[j] = 0;
      }
      
      i++;
    }
    
    else {
      characteristics[j] = 0;
    }
  }
}


//--------------------------------------------------------------
// Constructeur par copie de la classe Switching_sequence.
//--------------------------------------------------------------
Switching_sequence::Switching_sequence(const Switching_sequence &sw_seq , char transform ,
                                         int param1 , int param2)
{

  switch (transform) {
  case 'c' :
    Sequences::copy(sw_seq , (param1 == REVERSE ? true : false)); 
    copy(sw_seq , param1);
    break;
  case 'a' :
    Sequences::add_state_variable(sw_seq);
    add_variable(sw_seq , param1 , param2);
    break;
  default :
    Sequences::copy(sw_seq);
    copy(sw_seq);
    break;
  }
}

//--------------------------------------------------------------
// Operateur d'assignement de la classe Switching_sequence.
//--------------------------------------------------------------
Switching_sequence& Switching_sequence::operator=(const Switching_sequence &sw_seq)
{

  if (&sw_seq != this) {
    remove();  
    Sequences::remove();
    Sequences::copy(sw_seq);
    copy(sw_seq);
  }
  
  return *this;
}

//--------------------------------------------------------------
// Initialisation de la 1ere variable.
//--------------------------------------------------------------
void Switching_sequence::state_variable_init(int itype)
{
  register int i , j;

  if (itype != type[0]) {
    if (type[0] == STATE) {
      if (self_transition) {
        for (i = 0;i < marginal[0]->nb_value;i++) {
          delete self_transition[i];
        }
        delete [] self_transition;

        self_transition = 0;
      }

      if (observation) {
        for (i = 1;i < nb_variable;i++) {
          for (j = 0;j < marginal[0]->nb_value;j++) {
            delete observation[i][j];
          }
          delete [] observation[i];
        }
        delete [] observation;

        observation = 0;
      }
    }

    type[0] = itype;
  }

  for (i = 1;i < nb_variable;i++) {
    type[i] = INT_VALUE;
  }
}

//--------------------------------------------------------------
// Suppression de la 1ere variable.
//--------------------------------------------------------------
Switching_sequence* Switching_sequence::remove_variable_1() const
{
  register int i, j, k;
  int *variable , *itype;
  double *pcovar, *ccovar;
  Switching_sequence *isw_seq;


  variable = new int[nb_variable - 1];
  itype = new int[nb_variable - 1];
  for (i = 0;i < nb_variable - 1;i++) {
    variable[i] = i + 1;
    itype[i] = type[i + 1];
  }

  if(!constant){
    isw_seq = new Switching_sequence(nb_variable - 1 , itype , nb_sequence ,
				     identifier , length , nb_covariable-1, constant,  false);
  }
  else {
    isw_seq = new Switching_sequence(nb_variable - 1 , itype , nb_sequence ,
				     identifier , length , nb_covariable, constant,  false);
  }

  isw_seq->Sequences::select_variable(*this , variable);

  isw_seq->covar = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    isw_seq->covar[i] = new double*[length[i]];
    for (j = 0;j < length[i];j++) {
      isw_seq->covar[i][j] = new double[nb_covariable];
      pcovar = isw_seq->covar[i][j];

      ccovar = covar[i][j];
      for (k = 0;k < nb_covariable;k++) {
	*pcovar++ = *ccovar++;
      }
    }
  }


  for (i = 0;i < isw_seq->nb_variable;i++) {
    if (characteristics[i + 1]) {
      isw_seq->characteristics[i] = new Sequence_characteristics(*(characteristics[i + 1]));
    }
  }

  delete [] variable;
  delete [] itype;

  return isw_seq;
}

//--------------------------------------------------------------
// Calcul de la quantite d'information en faisant l'hypothese
// de variables aleatoires independantes et equidistribuees.
//--------------------------------------------------------------
double Switching_sequence::iid_information_computation() const
{
  register int i;
  double information = 0.;

  for (i = (((type[0] != STATE) || (nb_variable == 1)) ? 0 : 1);i < nb_variable;i++) {
    information += marginal[i]->information_computation();
  }

  return information;
}


//--------------------------------------------------------------
//  Extraction des probabilites de chaque valeur en fonction de l'index
// (pour une variable donnee).
//--------------------------------------------------------------
void Switching_sequence::build_index_value(int variable)
{
  register int i , j;
  int total , *frequency;

  // creation d'un objet Curves

  characteristics[variable]->index_value = new Curves(marginal[variable]->nb_value ,
                                                      max_length , true , false);
  frequency = new int[marginal[variable]->nb_value];

  // calcul des probabilites de chaque valeur en fonction de l'index

  for (i = 0;i < max_length;i++) {
    for (j = 0;j < marginal[variable]->nb_value;j++) {
      frequency[j] = 0;
    }

    for (j = 0;j < nb_sequence;j++) {
      if (i < length[j]) {
        frequency[int_sequence[j][variable][i]]++;
      }
    }

    total = 0;
    for (j = 0;j < marginal[variable]->nb_value;j++) {
      total += frequency[j];
    }
    characteristics[variable]->index_value->frequency[i] = total;
    for (j = 0;j < marginal[variable]->nb_value;j++) {
      characteristics[variable]->index_value->point[j][i] = (double)frequency[j] / (double)total;
    }
  }

  delete [] frequency;
}

//--------------------------------------------------------------
// Construction des histogrammes du temps avant la premiere occurrence
// d'une valeur (pour une variable donnee).
//--------------------------------------------------------------
void Switching_sequence::build_first_occurrence_histogram(int variable)
{
  bool *occurrence;
  register int i , j;
  int nb_value , *psequence;
  Histogram **first_occurrence;

  // creation des histogrammes

  first_occurrence = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    first_occurrence[i] = new Histogram(max_length);
  }

  // mise a jour des histogrammes

  occurrence = new bool[marginal[variable]->nb_value];

  for (i = 0;i < nb_sequence;i++) {
    nb_value = 0;
    for (j = 0;j < marginal[variable]->nb_value;j++) {
      occurrence[j] = false;
    }

    psequence = int_sequence[i][variable];
    for (j = 0;j < length[i];j++) {
      if (!occurrence[*psequence]) {
        occurrence[*psequence] = true;
        (first_occurrence[*psequence]->frequency[j])++;

        nb_value++;
        if (nb_value == marginal[variable]->nb_value) {
          break;
        }
      }

      psequence++;
    }
  }

  delete [] occurrence;

  // extraction des caracteristiques des histogrammes

  for (i = 0;i < marginal[variable]->nb_value;i++) {
    first_occurrence[i]->nb_value_computation();
    first_occurrence[i]->offset_computation();
    first_occurrence[i]->nb_element_computation();
    first_occurrence[i]->max_computation();
    first_occurrence[i]->mean_computation();
    first_occurrence[i]->variance_computation();
  }

  characteristics[variable]->first_occurrence = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->first_occurrence[i] = new Histogram(*(first_occurrence[i]));
    delete first_occurrence[i];
  }
  delete [] first_occurrence;
}

//--------------------------------------------------------------
// Construction des histogrammes du temps de retour dans une valeur
// (pour une variable donnee).
//--------------------------------------------------------------
void Switching_sequence::build_recurrence_time_histogram(int variable)
{
  register int i , j;
  int *index , *psequence;
  Histogram **recurrence_time;

  // creation des histogrammes

  recurrence_time = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    recurrence_time[i] = new Histogram(max_length);
  }

  // mise a jour des histogrammes

  index = new int[marginal[variable]->nb_value];

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < marginal[variable]->nb_value;j++) {
      index[j] = I_DEFAULT;
    }
    psequence = int_sequence[i][variable];

    for (j = 0;j < length[i];j++) {
      if (index[*psequence] != I_DEFAULT) {
        (recurrence_time[*psequence]->frequency[j - index[*psequence]])++;
      }
      index[*psequence++] = j;
    }
  }

  delete [] index;

  // extraction des caracteristiques des histogrammes

  for (i = 0;i < marginal[variable]->nb_value;i++) {
    recurrence_time[i]->nb_value_computation();
    recurrence_time[i]->offset_computation();
    recurrence_time[i]->nb_element_computation();
    recurrence_time[i]->max_computation();
    recurrence_time[i]->mean_computation();
    recurrence_time[i]->variance_computation();
  }

  characteristics[variable]->recurrence_time = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->recurrence_time[i] = new Histogram(*(recurrence_time[i]));
    delete recurrence_time[i];
  }
  delete [] recurrence_time;
}

//--------------------------------------------------------------
//  Accumulation du temps de sejour dans une valeur
// (pour une variable donnee).
//--------------------------------------------------------------
void Switching_sequence::sojourn_time_histogram_computation(int variable)
{
  register int i , j;
  int run_length , *pfrequency , *psequence;

  // initialisation des histogrammes

  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->sojourn_time[i]->offset = 1;
    characteristics[variable]->sojourn_time[i]->nb_value = characteristics[variable]->sojourn_time[i]->alloc_nb_value;

    pfrequency = characteristics[variable]->sojourn_time[i]->frequency;
    for (j = 0;j < characteristics[variable]->sojourn_time[i]->nb_value;j++) {
      *pfrequency++ = 0;
    }

    if (characteristics[variable]->initial_run) {
      characteristics[variable]->initial_run[i]->offset = 1;
      characteristics[variable]->initial_run[i]->nb_value = characteristics[variable]->initial_run[i]->alloc_nb_value;

      pfrequency = characteristics[variable]->initial_run[i]->frequency;
      for (j = 0;j < characteristics[variable]->initial_run[i]->nb_value;j++) {
        *pfrequency++ = 0;
      }
    }

    characteristics[variable]->final_run[i]->offset = 1;
    characteristics[variable]->final_run[i]->nb_value = characteristics[variable]->final_run[i]->alloc_nb_value;

    pfrequency = characteristics[variable]->final_run[i]->frequency;
    for (j = 0;j < characteristics[variable]->final_run[i]->nb_value;j++) {
      *pfrequency++ = 0;
    }
  }

  // mise a jour des histogrammes

  for (i = 0;i < nb_sequence;i++) {
    psequence = int_sequence[i][variable];
    run_length = 1;
    for (j = 1;j < length[i];j++) {
      if (*(psequence + 1) != *psequence) {
        if ((characteristics[variable]->initial_run) && (run_length == j)) {
          (characteristics[variable]->initial_run[*psequence]->frequency[run_length])++;
        }
        else {
          (characteristics[variable]->sojourn_time[*psequence]->frequency[run_length])++;
        }
        run_length = 0;
      }

      run_length++;
      psequence++;
    }

    if ((characteristics[variable]->initial_run) && (run_length == length[i])) {
      (characteristics[variable]->initial_run[*psequence]->frequency[run_length])++;
    }
    (characteristics[variable]->final_run[*psequence]->frequency[run_length])++;
  }

  // extraction des caracteristiques des histogrammes

  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->sojourn_time[i]->nb_value_computation();
    characteristics[variable]->sojourn_time[i]->offset_computation();
    characteristics[variable]->sojourn_time[i]->nb_element_computation();
    characteristics[variable]->sojourn_time[i]->max_computation();
    characteristics[variable]->sojourn_time[i]->mean_computation();
    characteristics[variable]->sojourn_time[i]->variance_computation();
  }

  if (characteristics[variable]->initial_run) {
    for (i = 0;i < marginal[variable]->nb_value;i++) {
      characteristics[variable]->initial_run[i]->nb_value_computation();
      characteristics[variable]->initial_run[i]->offset_computation();
      characteristics[variable]->initial_run[i]->nb_element_computation();
      characteristics[variable]->initial_run[i]->max_computation();
      characteristics[variable]->initial_run[i]->mean_computation();
      characteristics[variable]->initial_run[i]->variance_computation();
    }
  }

  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->final_run[i]->nb_value_computation();
    characteristics[variable]->final_run[i]->offset_computation();
    characteristics[variable]->final_run[i]->nb_element_computation();
    characteristics[variable]->final_run[i]->max_computation();
    characteristics[variable]->final_run[i]->mean_computation();
    characteristics[variable]->final_run[i]->variance_computation();
  }
}

//--------------------------------------------------------------
// Construction des histogrammes du temps de sejour dans une valeur
// (pour une variable donnee).
//--------------------------------------------------------------
void Switching_sequence::build_sojourn_time_histogram(int variable , int initial_run_flag)
{
  register int i , j;
  int run_length , *psequence;
  Histogram **sojourn_time , **initial_run , **final_run;

  // creation des histogrammes

  sojourn_time = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    sojourn_time[i] = new Histogram(max_length + 1);
  }

  if (initial_run_flag) {
    initial_run = new Histogram*[marginal[variable]->nb_value];
    for (i = 0;i < marginal[variable]->nb_value;i++) {
      initial_run[i] = new Histogram(max_length + 1);
    }
  }

  final_run = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    final_run[i] = new Histogram(max_length + 1);
  }

  // mise a jour des histogrammes

  for (i = 0;i < nb_sequence;i++) {
    psequence = int_sequence[i][variable];
    run_length = 1;
    for (j = 1;j < length[i];j++) {
      if (*(psequence + 1) != *psequence) {
        if ((initial_run_flag) && (run_length == j)) {
          (initial_run[*psequence]->frequency[run_length])++;
        }
        else {
          (sojourn_time[*psequence]->frequency[run_length])++;
        }
        run_length = 0;
      }

      run_length++;
      psequence++;
    }

    if ((initial_run_flag) && (run_length == length[i])) {
      (initial_run[*psequence]->frequency[run_length])++;
    }
    (final_run[*psequence]->frequency[run_length])++;
  }

  // extraction des caracteristiques des histogrammes

  for (i = 0;i < marginal[variable]->nb_value;i++) {
    sojourn_time[i]->nb_value_computation();
    sojourn_time[i]->offset_computation();
    sojourn_time[i]->nb_element_computation();
    sojourn_time[i]->max_computation();
    sojourn_time[i]->mean_computation();
    sojourn_time[i]->variance_computation();
  }

  if (initial_run_flag) {
    for (i = 0;i < marginal[variable]->nb_value;i++) {
      initial_run[i]->nb_value_computation();
      initial_run[i]->offset_computation();
      initial_run[i]->nb_element_computation();
      initial_run[i]->max_computation();
      initial_run[i]->mean_computation();
      initial_run[i]->variance_computation();
    }
  }

  for (i = 0;i < marginal[variable]->nb_value;i++) {
    final_run[i]->nb_value_computation();
    final_run[i]->offset_computation();
    final_run[i]->nb_element_computation();
    final_run[i]->max_computation();
    final_run[i]->mean_computation();
    final_run[i]->variance_computation();
  }

  characteristics[variable]->sojourn_time = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->sojourn_time[i] = new Histogram(*(sojourn_time[i]));
    delete sojourn_time[i];
  }
  delete [] sojourn_time;

  if (initial_run_flag) {
    characteristics[variable]->initial_run = new Histogram*[marginal[variable]->nb_value];
    for (i = 0;i < marginal[variable]->nb_value;i++) {
      characteristics[variable]->initial_run[i] = new Histogram(*(initial_run[i]));
      delete initial_run[i];
    }
    delete [] initial_run;
  }

  characteristics[variable]->final_run = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->final_run[i] = new Histogram(*(final_run[i]));
    delete final_run[i];
  }
  delete [] final_run;
}

//--------------------------------------------------------------
//  Construction des histogrammes du nombre de series d'une valeur par sequence
// (pour une variable donnee).
//--------------------------------------------------------------
void Switching_sequence::build_nb_run_histogram(int variable)
{
  register int i , j;
  int *psequence , *count;
  Histogram **nb_run;

  // creation des histogrammes

  nb_run = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    nb_run[i] = new Histogram((max_length % 2 == 0 ?
                               max_length / 2 : max_length / 2 + 1) + 1);
  }

  // mise a jour des histogrammes

  count = new int[marginal[variable]->nb_value];

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < marginal[variable]->nb_value;j++) {
      count[j] = 0;
    }

    psequence = int_sequence[i][variable];
    count[*psequence++]++;
    for (j = 1;j < length[i];j++) {
      if (*psequence != *(psequence - 1)) {
        count[*psequence]++;
      }
      psequence++;
    }

    for (j = 0;j < marginal[variable]->nb_value;j++) {
      (nb_run[j]->frequency[count[j]])++;
    }
  }

  delete [] count;

  // extraction des caracteristiques des histogrammes

  for (i = 0;i < marginal[variable]->nb_value;i++) {
    nb_run[i]->nb_value_computation();
    nb_run[i]->offset_computation();
    nb_run[i]->nb_element = nb_sequence;
    nb_run[i]->max_computation();
    nb_run[i]->mean_computation();
    nb_run[i]->variance_computation();
  }

  characteristics[variable]->nb_run = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->nb_run[i] = new Histogram(*(nb_run[i]));
    delete nb_run[i];
  }
  delete [] nb_run;
}

//--------------------------------------------------------------
// Construction des histogrammes du nombre d'occurrences
// d'une valeur par sequence (pour une variable donnee).
//--------------------------------------------------------------
void Switching_sequence::build_nb_occurrence_histogram(int variable)
{
  register int i , j;
  int *psequence , *count;
  Histogram **nb_occurrence;

  // creation des histogrammes

  nb_occurrence = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    nb_occurrence[i] = new Histogram(max_length + 1);
  }

  // mise a jour des histogrammes

  count = new int[marginal[variable]->nb_value];

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < marginal[variable]->nb_value;j++) {
      count[j] = 0;
    }

    psequence = int_sequence[i][variable];
    for (j = 0;j < length[i];j++) {
      count[*psequence++]++;
    }

    for (j = 0;j < marginal[variable]->nb_value;j++) {
      (nb_occurrence[j]->frequency[count[j]])++;
    }
  }

  delete [] count;

  // extraction des caracteristiques des histogrammes

  for (i = 0;i < marginal[variable]->nb_value;i++) {
    nb_occurrence[i]->nb_value_computation();
    nb_occurrence[i]->offset_computation();
    nb_occurrence[i]->nb_element = nb_sequence;
    nb_occurrence[i]->max_computation();
    nb_occurrence[i]->mean_computation();
    nb_occurrence[i]->variance_computation();
  }

  characteristics[variable]->nb_occurrence = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->nb_occurrence[i] = new Histogram(*(nb_occurrence[i]));
    delete nb_occurrence[i];
  }
  delete [] nb_occurrence;
}

//--------------------------------------------------------------
// Extraction des caracteristiques d'un echantillon de sequences.
// Calcul du nombre de valeurs et de l'histogramme des valeurs,
// extraction des probabilites des valeurs en fonction de l'index,
// construction des histogrammes du temps avant la premiere occurrence d'une valeur,
// des histogrammes du temps de retour dans une valeur,
// des histogrammes du temps de sejour dans une valeur,
// des histogrammes du nombre de series d'une valeur par sequence et
// des histogrammes du nombre d'occurrences d'une valeur par sequence.
//--------------------------------------------------------------
void Switching_sequence::build_characteristic(int variable , bool sojourn_time_flag ,
                                               bool initial_run_flag)
{
  register int i , j;
  bool build;

  for (i = 0;i < nb_variable;i++) {
    if ((variable == I_DEFAULT) || (i == variable)) {
      build = true;

      if (marginal[i]->nb_value > NB_OUTPUT) {
        build = false;
      }

      else {
        for (j = 0;j < marginal[i]->nb_value;j++) {
          if (marginal[i]->frequency[j] == 0) {
            build = false;
            break;
          }
        }
      }



      if (build) {
        if (sojourn_time_flag) {
          characteristics[i] = new Sequence_characteristics(marginal[i]->nb_value);
        }

        build_index_value(i);

	build_first_occurrence_histogram(i);
	build_recurrence_time_histogram(i);


        switch (sojourn_time_flag) {
        case false :
	  sojourn_time_histogram_computation(i);
          break;
        case true :
	  build_sojourn_time_histogram(i , initial_run_flag);
          break;
        }

        if (max_length <= COUNT_MAX_LENGTH) {
	  build_nb_run_histogram(i);
	  build_nb_occurrence_histogram(i);
        }
      }
    }
  }
}

//--------------------------------------------------------------
// Ecriture d'un objet Switching_sequence.
//--------------------------------------------------------------
ostream& Switching_sequence::ascii_write(ostream &os , bool exhaustive , bool comment_flag) const
{
  register int i;

  os << nb_variable << " " << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] << endl;
  os << nb_covariable << " " << "COVARIABLE(S)" << endl;
  if (!constant) {
    os << "WITH_CONSTANT" << " " << endl;
  }
  else {
    os << "WITHOUT_CONSTANT" << " " << endl;
  }

  if ((self_transition) && (exhaustive)) {
    for (i = 0;i < marginal[0]->nb_value;i++) {
      if (self_transition[i]) {
        os << "\n";
        if (comment_flag) {
          os << "# ";
        }
        os << "   | " << STAT_label[STATL_STATE] << " " << i << " - "
           << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_SELF_TRANSITION] << endl;

        self_transition[i]->ascii_print(os , comment_flag);
      }
    }
  }

  for (i = 0;i < nb_variable;i++) {
    os << "\n" << STAT_word[STATW_VARIABLE] << " " << i + 1 << " : "
       << STAT_variable_word[type[i]] << "   ";
    if (comment_flag) {
      os << "# ";
    }
    if (marginal[i])  {
      os << marginal[i]->nb_value << " ";
      
      switch (type[i]) {
      case STATE :
	os << STAT_label[marginal[i]->nb_value == 1 ? STATL_STATE : STATL_STATES] << endl;
	break;
      case INT_VALUE :
	os << STAT_label[marginal[i]->nb_value == 1 ? STATL_VALUE : STATL_VALUES] << endl;
	break;
      }

      os << "\n";
      if (comment_flag) {
	os << "# ";
      }
      os << STAT_label[type[i] == STATE ? STATL_STATE : STATL_VALUE] << " "
	 << STAT_label[STATL_HISTOGRAM] << " - ";
      marginal[i]->ascii_characteristic_print(os , false , comment_flag);
      
      if ((marginal[i]->nb_value <= ASCII_NB_VALUE) || (exhaustive)) {
	os << "\n";
	if (comment_flag) {
	  os << "# ";
	}
	os << "   | " << STAT_label[type[i] == STATE ? STATL_STATE : STATL_VALUE] << " "
	   << STAT_label[STATL_HISTOGRAM] << endl;
	marginal[i]->ascii_print(os , comment_flag);
      }
    }
 
    if (characteristics[i]) {
      characteristics[i]->ascii_print(os , type[i] , *hlength , exhaustive , comment_flag);
    }
  }

  os << "\n";
  if (comment_flag) {
    os << "# ";
  }
  os << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
  hlength->ascii_characteristic_print(os , false , comment_flag);

  if (exhaustive) {
    os << "\n";
    if (comment_flag) {
      os << "# ";
    }
    os << "   | " << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << endl;
    hlength->ascii_print(os , comment_flag);
  }

  os << "\n";
  if (comment_flag) {
    os << "# ";
  }
  os << SEQ_label[SEQL_CUMUL_LENGTH] << ": " << cumul_length << endl;

  return os;
}

//--------------------------------------------------------------
// Ecriture d'un objet Switching_sequence.
//--------------------------------------------------------------
ostream& Switching_sequence::ascii_write(ostream &os , bool exhaustive) const
{
  return ascii_write(os , exhaustive , false);
}

//--------------------------------------------------------------
// Ecriture d'un objet Switching_sequence dans un fichier.
//--------------------------------------------------------------
bool Switching_sequence::ascii_write(Format_error &error , const char *path ,
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

//--------------------------------------------------------------
// Ecriture d'un objet Switching_sequence.
//--------------------------------------------------------------
ostream& Switching_sequence::ascii_data_write(ostream &os , char format , bool exhaustive) const
{
  ascii_write(os , exhaustive , false);
  ascii_print(os , format , false);

  return os;
}

//--------------------------------------------------------------
// Ecriture d'un objet Switching_sequence dans un fichier.
//--------------------------------------------------------------
bool Switching_sequence::ascii_data_write(Format_error &error , const char *path ,
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
      ascii_write(out_file , exhaustive , true);
    }
    ascii_print(out_file , format , true);
  }

  return status;
}

//--------------------------------------------------------------
// Ecriture d'un objet Switching_sequence dans un fichier au format tableur.
//--------------------------------------------------------------
bool Switching_sequence::spreadsheet_write(Format_error &error , const char *path) const
{
  bool status;
  register int i;
  Curves *smoothed_curves;
  ofstream out_file(path);

  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    out_file << nb_variable << "\t" << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] << endl;
    out_file << nb_covariable << "\t" << "COVARIABLE(S)" << endl;
    if (!constant) {
      out_file << "WITH_CONSTANT" << "\t" << endl;
    }
    else {
      out_file << "WITHOUT_CONSTANT" << "\t" << endl;
    }
    if (self_transition) {
      for (i = 0;i < marginal[0]->nb_value;i++) {
        if (self_transition[i]) {
          out_file << "\n\t" << STAT_label[STATL_STATE] << " " << i << " - "
                   << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_SELF_TRANSITION] << endl;
          self_transition[i]->spreadsheet_print(out_file);
	  
          smoothed_curves = new Curves(*(self_transition[i]) , 's');
	  
          out_file << "\n\t" << STAT_label[STATL_STATE] << " " << i << " - "
                   << SEQ_label[SEQL_SMOOTHED] << " " << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_SELF_TRANSITION] << endl;
          smoothed_curves->spreadsheet_print(out_file);
	  
          delete smoothed_curves;
        }
      }
    }
    
    for (i = 0;i < nb_variable;i++) {
      if(marginal[i]) {
	out_file << "\n" << STAT_word[STATW_VARIABLE] << "\t" << i + 1 << "\t\t"
		 << marginal[i]->nb_value << " ";
	
	switch (type[i]) {
	case STATE :
	  out_file << STAT_label[marginal[i]->nb_value == 1 ? STATL_STATE : STATL_STATES] << endl;
	  break;
	case INT_VALUE :
	  out_file << STAT_label[marginal[i]->nb_value == 1 ? STATL_VALUE : STATL_VALUES] << endl;
	  break;
	}
	
	out_file << "\n" << STAT_label[type[i] == STATE ? STATL_STATE : STATL_VALUE] << " "
		 << STAT_label[STATL_HISTOGRAM] << "\t";
	marginal[i]->spreadsheet_characteristic_print(out_file);
	
	out_file << "\n\t" << STAT_label[type[i] == STATE ? STATL_STATE : STATL_VALUE] << " "
		 << STAT_label[STATL_HISTOGRAM] << endl;
	marginal[i]->spreadsheet_print(out_file);
      }
      if (characteristics[i]) {
	characteristics[i]->spreadsheet_print(out_file , type[i] , *hlength);
      }
    }
    
    out_file << "\n" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << "\t";
    hlength->spreadsheet_characteristic_print(out_file);
    
    out_file << "\n\t" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << endl;
    hlength->spreadsheet_print(out_file);
    
    out_file << "\n" << SEQ_label[SEQL_CUMUL_LENGTH] << "\t" << cumul_length << endl;
  }

  return status;
}


//--------------------------------------------------------------
// Accumulation des observations (pour une variable donnee).
//--------------------------------------------------------------
void Switching_sequence::observation_histogram_computation(int variable)
{
  register int i , j;
  int *pfrequency;
  int *pstate , *poutput; 

  // initialisation des histogrammes

  for (i = 0;i < marginal[0]->nb_value;i++) {
    pfrequency = observation[variable][i]->frequency;
    for (j = 0;j < marginal[variable]->nb_value;j++) {
      *pfrequency++ = 0;
    }
  }

  // mise a jour des histogrammes

  for (i = 0;i < nb_sequence;i++) {
    pstate = int_sequence[i][0];
    poutput = int_sequence[i][variable];
    for (j = 0;j < length[i];j++) {
      (observation[variable][*pstate++]->frequency[*poutput++])++;
    }
  }

  // extraction des caracteristiques des histogrammes

  for (i = 0;i < marginal[0]->nb_value;i++) {
    if (!characteristics[variable]) {
      observation[variable][i]->nb_value_computation();
    }
    observation[variable][i]->offset_computation();
    observation[variable][i]->nb_element_computation();
    observation[variable][i]->max_computation();

    if (!characteristics[variable]) {
      observation[variable][i]->mean_computation();
      observation[variable][i]->variance_computation();
    }
  }
}


//--------------------------------------------------------------
// Accumulation des observations.
//--------------------------------------------------------------
void Switching_sequence::observation_histogram_computation()
{
  register int i;

  for (i = 1;i < nb_variable;i++) {
    observation_histogram_computation(i);
  }
}

//--------------------------------------------------------------
// Construction des histogrammes correspondant aux probabilites 
// d'observation.
//--------------------------------------------------------------
void Switching_sequence::create_observation_histogram(int nb_state)
{
  if ((nb_variable > 1) && (!observation)) {
    register int i , j;

    observation = new Histogram**[nb_variable];
    observation[0] = 0;

    for (i = 1;i < nb_variable;i++) {
      observation[i] = new Histogram*[nb_state];
      for (j = 0;j < nb_state;j++) {
        observation[i][j] = new Histogram(marginal[i]->nb_value);
      }
    }
  }
}

//--------------------------------------------------------------
// Construction des histogrammes correspondant aux probabilites
// d'observation.
//--------------------------------------------------------------
void Switching_sequence::build_observation_histogram()
{
  create_observation_histogram(marginal[0]->nb_value); 
  observation_histogram_computation();
}
