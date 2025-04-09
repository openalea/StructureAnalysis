#include<math.h>
#include<iostream>
#include<iomanip>
#include<limits.h>
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
#include "switching_process.h"
#include "markov_switching.h"

using namespace std;

//--------------------------------------------------------------
// Constructeur par defaut de la classe Markov_switching_data.
//--------------------------------------------------------------
Markov_switching_data::Markov_switching_data()
{
  sw_markov = 0;
  chain_data = 0;
  likelihood = D_INF;
  hidden_likelihood = D_INF;
  posterior_probability = 0;
}

//--------------------------------------------------------------
// Constructeur de la classe Markov_switching_data. 
//--------------------------------------------------------------
Markov_switching_data::Markov_switching_data(int inb_variable , const Histogram &ihlength ,
					     int inb_covariable, int iconstant, bool init_flag)
  :Switching_sequence(inb_variable , ihlength , inb_covariable, iconstant, init_flag)
{
  sw_markov = 0;
  chain_data = 0;
  likelihood = D_INF;
  hidden_likelihood = D_INF;
  posterior_probability = 0;
}

//--------------------------------------------------------------
// Constructeur de la classe Markov_switching_data.
//--------------------------------------------------------------
Markov_switching_data::Markov_switching_data(int inb_sequence , int inb_covariable,
					     int *ilength , int ***isequence , double ***icovar,
					     int iconstant,int *iidentifier, const Markov_switching &isw_markov)
:Switching_sequence(isw_markov.nb_output_process + 1 , inb_sequence , inb_covariable, 
		    ilength , isequence, icovar, iidentifier, iconstant)
{
  register int i;
  type[0] = STATE;
  sw_markov = new Markov_switching(isw_markov, false);
  // sw_markov = 0;


  // extraction des caracteristiques des sequences simulees

  build_transition_count(*sw_markov);
  //build_observation_histogram();

  sw_markov->characteristic_computation(*this , true);

  // calcul de la vraisemblance

  likelihood = sw_markov->likelihood_computation(*this);

  posterior_probability = new double[inb_sequence];
  for(i = 0; i < inb_sequence; i++){
    posterior_probability[i] = 0.;
  }


}

//--------------------------------------------------------------
// Construction d'un objet Markov_switching_data a partir
// d'un objet Switching_sequence.
//--------------------------------------------------------------
Markov_switching_data::Markov_switching_data(const Switching_sequence &isw_seq ,
                                                       bool initial_run_flag)
:Switching_sequence(isw_seq , 'c' , (initial_run_flag ? ADD_INITIAL_RUN : REMOVE_INITIAL_RUN))
{
  sw_markov = 0;
  chain_data = 0;
  likelihood = D_INF;
  hidden_likelihood = D_INF;
  posterior_probability = 0;
}


//--------------------------------------------------------------
// Construction d'un objet Markov_switching_data a partir
// d'un objet Switching_sequence avec ajout d'une variable.
//--------------------------------------------------------------
Markov_switching_data::Markov_switching_data(const Switching_sequence &isw_seq ,
                                                       int variable)
:Switching_sequence(isw_seq , 'a' , variable , DEFAULT)
{
  sw_markov = 0;
  chain_data = 0;
  likelihood = D_INF;
  hidden_likelihood = D_INF;
  posterior_probability = 0;
}

//--------------------------------------------------------------
// Copie d'un objet Markov_switching_data.
//---------------------------------------------------------------
void Markov_switching_data::copy(const Markov_switching_data &isw_markov_data ,
                                      bool model_flag)
{
  register int i;

  if ((model_flag) && (isw_markov_data.sw_markov)) {
    sw_markov = new Markov_switching(*(isw_markov_data.sw_markov) , false);
  }
  else {
    sw_markov = 0;
  }

  if (isw_markov_data.chain_data) {
    chain_data = new Chain_data(*(isw_markov_data.chain_data));
  }
  else {
    chain_data = 0;
  }

  likelihood = isw_markov_data.likelihood;
  hidden_likelihood = isw_markov_data.hidden_likelihood;

  if (isw_markov_data.posterior_probability){
    posterior_probability = new double[nb_sequence];
    for (i =0; i < nb_sequence; i++){
      posterior_probability[i] = isw_markov_data.posterior_probability[i];
    }
  }
  else{
    posterior_probability = 0;
  }

}

//--------------------------------------------------------------
// Destructeur de la classe Markov_switching_data.
//--------------------------------------------------------------
Markov_switching_data::~Markov_switching_data()
{
  if (sw_markov){
    delete sw_markov;
  }

  if(chain_data){
    delete chain_data;
  }
 
  delete [] posterior_probability;
 
}

//--------------------------------------------------------------
// Operateur d'assignement de la classe Markov_switching_data.
//--------------------------------------------------------------
Markov_switching_data& Markov_switching_data::operator=(const Markov_switching_data &isw_markov_data)
{
  if (&isw_markov_data != this) {
    delete sw_markov;
    delete chain_data;

    remove();
    Sequences::remove();

    Sequences::copy(isw_markov_data);
    Switching_sequence::copy(isw_markov_data);
    copy(isw_markov_data);
  }

  return *this;
}


//--------------------------------------------------------------
// Extraction d'un histogramme.
//
// arguments : reference sur un objet Format_error, type d'histogramme,
//             variable, valeur (etat dans le cas des histogrammes d'observation).
//--------------------------------------------------------------
Distribution_data* Markov_switching_data::extract(Format_error &error , int type ,
                                                       int variable , int value) const
{
  bool status = true;
  Distribution *pdist;
  Parametric *pparam;
  Histogram *phisto;
  Distribution_data *histo;

  histo = 0;
  error.init();

 //  if (type == OBSERVATION) {
//     if ((variable < 2) || (variable > nb_variable)) {
//       status = false;
//       error.update(STAT_error[STATR_VARIABLE_INDEX]);
//     }

//     else {
//       variable--;

//       if ((value < 0) || (value >= marginal[0]->nb_value)) {
//         status = false;
//         ostringstream error_message;
//         error_message << STAT_label[STATL_STATE] << " " << value << " "
//                       << SEQ_error[SEQR_NOT_PRESENT];
//         error.update((error_message.str()).c_str());
//       }

//       else {
//         phisto = observation[variable][value];

//         if (phisto->nb_element == 0) {
//           status = false;
//           error.update(STAT_error[STATR_EMPTY_HISTOGRAM]);
//         }
//       }
//     }
//   }

//   else {
//     if ((variable < 1) || (variable > nb_variable)) {
//       status = false;
//       error.update(STAT_error[STATR_VARIABLE_INDEX]);
//     }

//     else {
//       variable--;

//       if (!characteristics[variable]) {
//         status = false;
//         ostringstream error_message;
//         error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
//                       << SEQ_error[SEQR_CHARACTERISTICS_NOT_COMPUTED];
//         error.update((error_message.str()).c_str());
//       }

//       else if ((value < 0) || (value >= marginal[variable]->nb_value)) {
//         status = false;
//         ostringstream error_message;
//         error_message << STAT_label[variable == 0 ? STATL_STATE : STATL_OUTPUT] << " "
//                       << value << " " << SEQ_error[SEQR_NOT_PRESENT];
//         error.update((error_message.str()).c_str());
//       }

//       else {
//         switch (type) {
//         case FIRST_OCCURRENCE :
//           phisto = characteristics[variable]->first_occurrence[value];
//           break;
//         case RECURRENCE_TIME :
//           phisto = characteristics[variable]->recurrence_time[value];
//           break;
//         case SOJOURN_TIME :
//           phisto = characteristics[variable]->sojourn_time[value];
//           break;
//         case FINAL_RUN :
//           phisto = characteristics[variable]->final_run[value];
//           break;
//         case NB_RUN :
//           phisto = characteristics[variable]->nb_run[value];
//           break;
//         case NB_OCCURRENCE :
//           phisto = characteristics[variable]->nb_occurrence[value];
//           break;
//         }

//         if (phisto->nb_element == 0) {
//           status = false;
//           error.update(STAT_error[STATR_EMPTY_HISTOGRAM]);
//         }
//       }
//     }
//   }

//   if (status) {
//     pdist = 0;
//     pparam = 0;

//     switch (type) {

//     case OBSERVATION : {
//         pparam = sw_markov->switching_process[variable]->observation[value];
//       break;
//     }

//     case FIRST_OCCURRENCE : {
//       pdist = sw_markov->switching_process[variable]->first_occurrence[value];
//       break;
//     }

//     case RECURRENCE_TIME : {
//       pdist = sw_markov->switching_process[variable]->recurrence_time[value];
//       break;
//     }

//     case SOJOURN_TIME : {
//       pparam = sw_markov->switching_process[variable]->sojourn_time[value];
//       break;
//     }

//     case NB_RUN : {
//       pdist = sw_markov->switching_process[variable]->nb_run[value];
//       break;
//     }

//     case NB_OCCURRENCE : {
//       pdist = sw_markov->switching_process[variable]->nb_occurrence[value];
//       break;
//     }
//     }

//     if (pdist) {
//       histo = new Distribution_data(*phisto , pdist);
//     }
//     else {
//       histo = new Distribution_data(*phisto , pparam);
//     }
//   }

  return histo;
}

//--------------------------------------------------------------
// Ecriture d'un objet Markov_switching_data
//--------------------------------------------------------------
ostream& Markov_switching_data::ascii_write(ostream &os , bool exhaustive) const
{
  
  Switching_sequence::ascii_write(os, exhaustive);

  if (sw_markov) {
    sw_markov->ascii_write(os , this , exhaustive ,
			   false);// A REVOIR
    // ::test_hidden(sw_markov->nb_output_process , sw_markov->sw_process));
  }

  return os;
}

//--------------------------------------------------------------
// Ecriture d'un objet Markov_switching_data dans un fichier.
//--------------------------------------------------------------
bool Markov_switching_data::ascii_write(Format_error &error , const char *path ,
                                             bool exhaustive) const
{
  bool status = false;
  
  if (sw_markov) {
    ofstream out_file(path);
    error.init();
    
    if (!out_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_NAME]);
    }
    
    else {
      status = true;
      sw_markov->ascii_write(out_file , this , exhaustive , true ); //A REVOIR
      //                          ::test_hidden(sw_markov->nb_output_process , sw_markov->sw_process));
    }
  }
  
  return status;
}

//--------------------------------------------------------------
// Ecriture d'un objet Markov_switching_data dans un fichier au format tableur.
//--------------------------------------------------------------
bool Markov_switching_data::spreadsheet_write(Format_error &error , const char *path) const
{
  bool status = false;

  if (sw_markov) {
    ofstream out_file(path);

    error.init();

    if (!out_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_NAME]);
    }

    else {
      status = true;
      sw_markov->spreadsheet_write(out_file , this ,
				   false); //A REVOIR
      //::test_hidden(sw_markov->nb_output_process , sw_markov->nonparametric_process));
    }
  }
  
  return status;
}

//--------------------------------------------------------------
// Comptage des etats initiaux et des transitions.
//--------------------------------------------------------------
void Switching_sequence::transition_count_computation(const Chain_data &chain_data ,
                                                       int order , bool begin) const
{
  register int i , j , k;
  int **ttransition , *pstate , *mstate , power[ORDER];


  for (i = 0;i < chain_data.nb_state;i++) {
    chain_data.initial[i] = 0;
  }

  for (i = 0;i < chain_data.nb_row;i++) {
    for (j = 0;j < chain_data.nb_state;j++) {
      chain_data.transition[i][j] = 0;
    }
  }

  // extraction des etats initiaux et des transitions

  i = 1;
  for (j = 0;j < order;j++) {
    power[j] = i;
    i *= chain_data.nb_state;
  }

  for (i = 0;i < nb_sequence;i++) {
    pstate = int_sequence[i][0];
    (chain_data.initial[*pstate])++;

    if (!begin) {
       pstate += order - 1;
    }
    for (j = (begin ? 1 : order);j < length[i];j++) {
      ttransition = chain_data.transition;
      mstate = pstate;
      for (k = 0;k < MIN(j - 1 , order);k++) {
        ttransition += *mstate-- * power[order - 1 - k];
      }
      for (k = j - 1;k < order;k++) {
        ttransition += *mstate * power[order - 1 - k];
      }

      (*(*ttransition + *++pstate))++;
    }
  }
}


//--------------------------------------------------------------
//Construction des comptages des etats initiaux et des transitions. 
//--------------------------------------------------------------
void Markov_switching_data::build_transition_count(const Markov_switching &isw_markov, bool begin)
{
  chain_data = new Chain_data('o' , marginal[0]->nb_value ,
                              (int)pow((double)marginal[0]->nb_value , isw_markov.order));
  transition_count_computation(*chain_data , isw_markov.order , begin);
}
