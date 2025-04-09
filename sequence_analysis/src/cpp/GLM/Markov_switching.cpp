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

//---------------------------------------------------------------
// Constructeur par defaut de la classe Markov_switching
//---------------------------------------------------------------
Markov_switching::Markov_switching()
{
  nb_iterator = 0;
  sw_markov_data = 0;
  order = 1 ;
  nb_output_process = 0;
  sw_process = 0; 

  state_subtype = 0;
  forward = 0;
  nonparametric_process = 0;
}

//--------------------------------------------------------------
// Constructeur de la classe Markov_switching
//--------------------------------------------------------------
Markov_switching::Markov_switching(int inb_state)
:Chain('o' , inb_state)
{
  nb_iterator = 0;
  sw_markov_data = 0;

  order = 1;

  state_subtype = 0;
  forward = 0;

  nb_output_process = 0;
  sw_process = new Switching_process*[1];
  sw_process[0] = new Switching_process(nb_state , nb_state);

  nonparametric_process = new Nonparametric_sequence_process*[1];
  nonparametric_process[0] = 0;
} 

//--------------------------------------------------------------
// Constructeur de la classe Markov_switching
//--------------------------------------------------------------
Markov_switching::Markov_switching(int inb_state , int iorder , int inb_output_process , int *nb_value)
:Chain('o' , inb_state , (int)pow((double)inb_state , iorder) , true)
{
  register int i;

  nb_iterator = 0;
  sw_markov_data = 0;

  state_subtype = 0;
  forward = 0;

  order = iorder;

  nb_output_process = inb_output_process;
  
  sw_process = new Switching_process*[nb_output_process + 1];
  sw_process[0] = 0;
  
  nonparametric_process = new Nonparametric_sequence_process*[nb_output_process+1];
  nonparametric_process[0] = new Nonparametric_sequence_process(nb_state, nb_state, false);

  for (i = 1;i <= nb_output_process;i++) {
    nonparametric_process[i] = 0;
    sw_process[i] = new Switching_process(nb_state , *nb_value++);
  }
}

//--------------------------------------------------------------
// Constructeur de la classe Markov_switching.
//--------------------------------------------------------------
Markov_switching::Markov_switching(const Chain *pchain , int iorder ,
				   const Nonparametric_sequence_process *poccupancy,
				   const Switching_process *pobservation , int length)
:Chain(*pchain)
{
  register int i;

  nb_iterator = 0;
  sw_markov_data = 0;

  order = iorder;

  nb_output_process = (pobservation ? 1 : 0);

  nonparametric_process = new Nonparametric_sequence_process*[nb_output_process +1];
  nonparametric_process[0] = new Nonparametric_sequence_process(*poccupancy);

  sw_process = new Switching_process*[nb_output_process + 1];
  sw_process[0] = 0;

  if (pobservation) {
    nonparametric_process[1] = 0;
    sw_process[1] = new Switching_process(*pobservation);
  }

  for(i = 0; i < nb_state; i++){
    if (transition[i][i] < 1.){
      nonparametric_process[0]->absorption[i] = 0.;
    }
    else {
      nonparametric_process[0]->absorption[i] = 1.;
    }
  }

  state_subtype = new int[nb_state];
  forward = new Forward*[nb_state];

  for(i = 0; i < nb_state; i++){
    state_subtype[i] = (nonparametric_process[0]->sojourn_time[i] ? SEMI_MARKOVIAN : MARKOVIAN);
    
    if ((state_subtype[i] == SEMI_MARKOVIAN) && (state_type[i] == 'r')){
      forward[i] = new Forward(*(nonparametric_process[0]->sojourn_time[i]));
    }
    else {
      forward[i] = 0;
    }
  }


  characteristic_computation(length , true); // pas implémentée
}

//--------------------------------------------------------------
// Construction d'une chaine de Markov_switching d'ordre superieur a 1 a partir
// d'une chaine de Markov_switching.
//--------------------------------------------------------------
Markov_switching::Markov_switching(int iorder , const Markov_switching &isw_markov)
  :Chain('o' , isw_markov.nb_state , (int)pow((double)isw_markov.nb_state , iorder),   true)
{
  register int i , j;

  accessibility = new bool*[nb_state];
  for (i = 0;i < nb_state;i++) {
    accessibility[i] = new bool[nb_state];
    for (j = 0;j < nb_state;j++) {
      accessibility[i][j] = isw_markov.accessibility[i][j];
    }
  }
  
  nb_component = isw_markov.nb_component;
  component_nb_state = new int[nb_component];
  component = new int*[nb_component];
  
  for (i = 0;i < nb_component;i++) {
    component_nb_state[i] = isw_markov.component_nb_state[i];
    component[i] = new int[component_nb_state[i]];
    for (j = 0;j < component_nb_state[i];j++) {
      component[i][j] = isw_markov.component[i][j];
    }
  }
  
  state_type = new char[nb_state];
  for (i = 0;i < nb_state;i++) {
    state_type[i] = isw_markov.state_type[i];
  }
  
  nb_iterator = 0;
  sw_markov_data = 0;

  order = iorder;
  
  nb_output_process = isw_markov.nb_output_process;
  
  sw_process = new Switching_process*[nb_output_process + 1];
  sw_process[0] = 0;
  for (i = 1;i <= nb_output_process;i++) {
    sw_process[i] = new Switching_process(*(isw_markov.sw_process[i]) , 'c' , false);
  }

  state_subtype = new int[nb_state];
  forward = new Forward*[nb_state];

  for (i = 0; i < nb_state; i++){
    state_subtype[i] = isw_markov.state_subtype[i];
    if (isw_markov.forward[i]){
      forward[i] = new Forward(*(isw_markov.forward[i]),0);
    }
    else {
      forward[i] = 0;
    }
  }

  nonparametric_process = new Nonparametric_sequence_process*[nb_output_process + 1];
  nonparametric_process[0] = new Nonparametric_sequence_process(*(isw_markov.nonparametric_process[0]), 'o', I_DEFAULT);
  for (i = 1; i <= nb_output_process; i++){
    nonparametric_process[i] = 0;
  }


}

//--------------------------------------------------------------
// Copie d'un objet Markov-switching
//--------------------------------------------------------------
void Markov_switching::copy(const Markov_switching &isw_markov , bool data_flag, int param)
{
  register int i;

  nb_iterator = 0;

  if ((data_flag) && (isw_markov.sw_markov_data)) {
    sw_markov_data = new Markov_switching_data(*(isw_markov.sw_markov_data) , false);
  }
  else {
    sw_markov_data = 0;
  }

  order = isw_markov.order;

  nb_output_process = isw_markov.nb_output_process;

  state_subtype = new int[nb_state];
  forward = new Forward*[nb_state];
  for (i = 0; i < nb_state; i++){
    state_subtype[i] = isw_markov.state_subtype[i];
    if (isw_markov.forward[i]){
      forward[i] = new Forward(*(isw_markov.forward[i]),param);
    }
    else {
      forward[i] = 0;
    }
  }

  nonparametric_process = new Nonparametric_sequence_process*[nb_output_process+1];
  sw_process = new Switching_process*[nb_output_process + 1];
  sw_process[0] = 0;

  switch(param){
  case I_DEFAULT: {
    nonparametric_process[0] = new Nonparametric_sequence_process(*(isw_markov.nonparametric_process[0]));
    break;
  }
  case 0:{
    nonparametric_process[0] = new Nonparametric_sequence_process(*(isw_markov, nonparametric_process[0]), 'o', I_DEFAULT);
    break;
  }
  default:{
    nonparametric_process[0] = new Nonparametric_sequence_process(*(isw_markov, nonparametric_process[0]), 'o', param);
    break;
  }
  }


  for (i = 1;i <= nb_output_process;i++) {
    nonparametric_process[i] = 0;
    sw_process[i] = new Switching_process(*(isw_markov.sw_process[i]) , 'c');
  }



}

//--------------------------------------------------------------
//Destruction des champs d'un objet Markov_switching
//--------------------------------------------------------------
void Markov_switching::remove()
{
  register int i;

  delete sw_markov_data;


  if (sw_process) {
    for (i = 1;i <= nb_output_process;i++) {
      delete sw_process[i];
    }
    delete [] sw_process;
  }

  delete [] state_subtype;

  if (forward){
    for (i = 0; i < nb_state; i++){
      delete forward[i];
    }
    delete [] forward;
  }

  if (nonparametric_process){
    delete nonparametric_process[0];
    delete [] nonparametric_process;
  }


}

//--------------------------------------------------------------
// Destructeur de la classe Markov_switching.
//--------------------------------------------------------------
Markov_switching::~Markov_switching()
{
  remove();
}

//--------------------------------------------------------------
// Operateur d'assignement de la classe Markov_switching.
//--------------------------------------------------------------
Markov_switching& Markov_switching::operator=(const Markov_switching &isw_markov)
{
  if (&isw_markov != this) {
    remove();
    Chain::remove();

    Chain::copy(isw_markov);
    copy(isw_markov);
  

  }
  return *this;
}

//--------------------------------------------------------------
// Extraction de la partie "donnees" d'un objet Markov_switching.
//--------------------------------------------------------------
Markov_switching_data* Markov_switching::extract_data(Format_error &error) const
{
  bool status = true;
  Markov_switching_data *psw_markov_data;

  psw_markov_data = 0;
  error.init();

  if (!sw_markov_data) {
    status = false;
    error.update(STAT_error[STATR_NO_DATA]);
  }
  else if (nb_output_process + 1 != sw_markov_data->nb_variable) {
    status = false;
    error.update(SEQ_error[SEQR_STATE_SEQUENCES]);
  }

  if (status) {
    psw_markov_data = new Markov_switching_data(*sw_markov_data);
    psw_markov_data->sw_markov = new Markov_switching(*this , false);
  }

  return psw_markov_data;
}

//--------------------------------------------------------------
// Application d'un seuil sur les parametres d'une chaine de Markov_switching.
//--------------------------------------------------------------
Markov_switching* Markov_switching::thresholding(double min_probability) const
{
  register int i;
  Markov_switching *psw_markov;

  psw_markov = new Markov_switching(*this , false); 

  psw_markov->Chain::thresholding(min_probability);
  psw_markov->component_computation();

  return psw_markov;
}

//--------------------------------------------------------------
//Ecriture sur une ligne d'un objet Markov_switching.
//--------------------------------------------------------------
ostream& Markov_switching::line_write(ostream &os) const
{
  os << nb_state << " " << STAT_word[STATW_STATES];
 
  return os;
}

//--------------------------------------------------------------
// Ecriture d'un objet Markov_switching et de la structure de donnees associee.
//--------------------------------------------------------------
ostream& Markov_switching::ascii_write(ostream &os , const Markov_switching_data *isw_markov_data , bool exhaustive ,
                             bool file_flag ) const
{
  register int i,j;
  int variable;
  Histogram **observation = 0;
  Sequence_characteristics *characteristics = 0;


  if( transition[0][0] > 0.){
    os << SEQ_word[SEQW_HIDDEN_MARKOV_CHAIN] << endl;
  }
  else {
    os << SEQ_word[SEQW_HIDDEN_SEMI_MARKOV_CHAIN] << endl;
  }

  // ecriture des parametres de la chaine de Markov
  os << "INITIAL PROBABILITIES" << endl;
  for (i = 0;i < nb_state;i++) {
    os << initial[i]<<"   ";
  }
  cout<<endl<<endl;
  
  os << "TRANSITION PROBABILITIES" << endl;
  for (i = 0;i < nb_state;i++) {
    for (j = 0; j < nb_state; j++){
      os<<transition[i][j]<<"   ";
    }
    cout<<endl;
  }
  cout<<endl;
  //ascii_print(os , file_flag );

 // ecriture des lois d'occupation des etats
  if ((isw_markov_data) && (isw_markov_data->type[0] == STATE)) {
    characteristics = isw_markov_data->characteristics[0];
  }
  nonparametric_process[0]->ascii_print(os , 0, 0, characteristics, exhaustive, file_flag, forward); 

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



  // ecriture des lois associees a chaque processus d'observation
  for (i = 1;i <= nb_output_process;i++) {
    os << "\n" << STAT_word[STATW_OUTPUT_PROCESS];
    os << " " << i;
    os << endl;

    if (isw_markov_data) { 
      switch (isw_markov_data->type[0]) {
      case INT_VALUE :
        variable = i - 1;
        break;
      case STATE :
        variable = i;
        break;
      }

      if (isw_markov_data->observation) {
        observation = isw_markov_data->observation[variable];
      }
      characteristics = isw_markov_data->characteristics[variable];
    }

    sw_process[i]->ascii_print(os , observation , exhaustive , file_flag); 
  }

  if (isw_markov_data) {
    int nb_parameter = nb_parameter_computation(MIN_PROBABILITY );
    double information , likelihood;

    // ecriture de l'histogramme des longueurs des sequences

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
    isw_markov_data->hlength->ascii_characteristic_print(os , false , file_flag);

    if (exhaustive) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << "   | " << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << endl;
      isw_markov_data->hlength->ascii_print(os , file_flag);
    }

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_CUMUL_LENGTH] << ": " << isw_markov_data->cumul_length << endl;

    // ecriture de la quantite d'information des sequences dans le cas iid

//     information = isw_markov_data->iid_information_computation();

//     os << "\n";
//     if (file_flag) {
//       os << "# ";
//     }
//     os << SEQ_label[SEQL_IID_INFORMATION] << ": " << information << " ("
//        << information / isw_markov_data->cumul_length << ")" << endl;

    // ecriture des vraisemblances des sequences

    if (isw_markov_data->likelihood != D_INF) {
      os << "\n";
      if (file_flag) {
	os << "# ";
      }
      os << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << ": " << isw_markov_data->likelihood << "   ("
           << STAT_label[STATL_NORMALIZED] << ": " << isw_markov_data->likelihood / isw_markov_data->cumul_length << ")" << endl;
    }
    
    likelihood = isw_markov_data->hidden_likelihood;

    if (likelihood != D_INF) {
      os << "\n";
      if (file_flag) {
	os << "# ";
    }
      os << "Marginal "<< SEQ_label[SEQL_OBSERVED_SEQUENCES_LIKELIHOOD] << ": " << likelihood << "   ("
	 << STAT_label[STATL_NORMALIZED] << ": " << likelihood / isw_markov_data->cumul_length << ")" << endl;
    }
    

    if ((likelihood != D_INF) || (nb_component == 1)) { // VRAi
      os << "\n";
      if (file_flag) {
	os << "# ";
      }
      os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
	 << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AIC] << "): "
	 << 2 * (likelihood - nb_parameter) << endl;
      
      if (nb_parameter < isw_markov_data->cumul_length - 1) {
	os << "\n";
	if (file_flag) {
	  os << "# ";
	}
	os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
	   << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AICc] << "): "
	   << 2 * (likelihood - (double)(nb_parameter * isw_markov_data->cumul_length) /
		   (double)(isw_markov_data->cumul_length - nb_parameter - 1)) << endl;
      }
      
      os << "\n";
      if (file_flag) {
	os << "# ";
      }
      os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
	 << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BIC] << "): "
	 << 2 * likelihood - nb_parameter * log((double)isw_markov_data->cumul_length) << endl;
    }
  }
  

  observation = NULL;
  characteristics = NULL;

  delete observation;
  delete characteristics;


  return os;
}

 
//--------------------------------------------------------------
//  Ecriture d'un objet Markov_switching
//--------------------------------------------------------------
ostream& Markov_switching::ascii_write(ostream &os , bool exhaustive) const
{
  return ascii_write(os , sw_markov_data , exhaustive , false);
}

//--------------------------------------------------------------
// Ecriture d'un objet Markov_switching dans un fichier.
//--------------------------------------------------------------
bool Markov_switching::ascii_write(Format_error &error , const char *path ,
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
    ascii_write(out_file , sw_markov_data , exhaustive , true);
  }

  return status;
}

//--------------------------------------------------------------
// Ecriture d'un objet Markov_switching et de la structure de donnees associee
// dans un fichier au format tableur.
//--------------------------------------------------------------
ostream& Markov_switching::spreadsheet_write(ostream &os , const Markov_switching_data *isw_markov_data ,
					     bool hidden) const
{
  register int i;
  int variable;
  Histogram **observation = 0;
  Sequence_characteristics *characteristics = 0;

  os << SEQ_word[SEQW_HIDDEN_MARKOV_CHAIN] << endl;
  
  // ecriture des parametres de la chaine de Markov
  
  spreadsheet_print(os);
  
 // ecriture des lois d'occupation des etats
  if ((isw_markov_data) && (isw_markov_data->type[0] == STATE)) {
    characteristics = isw_markov_data->characteristics[0];
  }
  nonparametric_process[0]->spreadsheet_print(os , 0, 0, characteristics, forward); 

 
  // ecriture des lois associees a chaque processus d'observation
  os << "\n" << nb_output_process << "\t"
     << STAT_word[nb_output_process == 1 ? STATW_OUTPUT_PROCESS : STATW_OUTPUT_PROCESSES] << endl;

  for (i = 1;i <= nb_output_process;i++) {
    os << "\n" << STAT_word[STATW_OUTPUT_PROCESS];
    os << "\t" << i;
    os << endl;

    if (isw_markov_data) {
      switch (isw_markov_data->type[0]) {
      case INT_VALUE :
        variable = i - 1;
        break;
      case STATE :
        variable = i;
        break;
      }
      
      if (isw_markov_data->observation) {
        observation = isw_markov_data->observation[variable];
      }
      characteristics = isw_markov_data->characteristics[variable];
    }
    
    sw_process[i]->spreadsheet_print(os , observation);
  }

  if (isw_markov_data) {
    int nb_parameter = nb_parameter_computation(MIN_PROBABILITY);
    double information , likelihood;
    
    // ecriture de l'histogramme des longueurs des sequences

    os << "\n" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << "\t";
    isw_markov_data->hlength->spreadsheet_characteristic_print(os);

    os << "\n\t" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << endl;
    isw_markov_data->hlength->spreadsheet_print(os);

    os << "\n" << SEQ_label[SEQL_CUMUL_LENGTH] << "\t" << isw_markov_data->cumul_length << endl;

    // ecriture de la quantite d'information des sequences dans le cas iid
    
    information = isw_markov_data->iid_information_computation();

    os << "\n" << SEQ_label[SEQL_IID_INFORMATION] << "\t" << information << "\t"
       << information / isw_markov_data->cumul_length << endl;
    
    // ecriture des vraisemblances des sequences
    
    if (isw_markov_data->likelihood != D_INF) {
      os << "\n" << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << "\t" << isw_markov_data->likelihood << "\t"
	 << STAT_label[STATL_NORMALIZED] << "\t" << isw_markov_data->likelihood / isw_markov_data->cumul_length << endl;
    }
    
    likelihood = isw_markov_data->hidden_likelihood;
    
    if (likelihood != D_INF) {
      os << "\n" << SEQ_label[SEQL_OBSERVED_SEQUENCES_LIKELIHOOD] << "\t" << likelihood << "\t"
	 << STAT_label[STATL_NORMALIZED] << "\t" << likelihood / isw_markov_data->cumul_length << endl;
    }
    
    
    if ((likelihood != D_INF) && (nb_component == 1)) {
      os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AIC] << ")\t"
         << 2 * (likelihood - nb_parameter) << endl;

      if (nb_parameter < isw_markov_data->cumul_length - 1) {
        os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
           << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AICc] << ")\t"
           << 2 * (likelihood - (double)(nb_parameter * isw_markov_data->cumul_length) /
              (double)(isw_markov_data->cumul_length - nb_parameter - 1)) << endl;
      }

      os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BIC] << ")\t"
         << 2 * likelihood - nb_parameter * log((double)isw_markov_data->cumul_length) << endl;
    }
  }
  
  observation = NULL;
  characteristics = NULL;

  delete observation;
  delete characteristics;


  return os;
}

//--------------------------------------------------------------
// Ecriture d'un objet Markov_switching dans un fichier au format tableur.
//--------------------------------------------------------------
bool Markov_switching::spreadsheet_write(Format_error &error , const char *path) const
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
    spreadsheet_write(out_file , sw_markov_data);
  }

  return status;
}

//---------------------------------------------------------------
// Calcul du nombre de paramètres indépendants d'un objet Markov_switching
//---------------------------------------------------------------
int Markov_switching::nb_parameter_computation(double min_probability) const
{
  register int i;
  int nb_parameter = Chain::nb_parameter_computation(min_probability);

  for (i = 0; i < nb_state; i++){
    if (state_subtype[i] == SEMI_MARKOVIAN){
      nb_parameter += nonparametric_process[0]->sojourn_time[i]->nb_parameter_computation();
      if(nonparametric_process[0]->sojourn_time[i]->inf_bound == 1){
	nb_parameter--;
      }
    }
  }

  for (i = 1;i <= nb_output_process;i++) {
    nb_parameter += sw_process[i]->nb_parameter_computation();
  }

  return nb_parameter;
}

//---------------------------------------------------------------
// Calcul du nombre de paramètres de transitions
// indépendants d'un objet Markov_switching
//---------------------------------------------------------------
int Markov_switching::nb_transient_parameter_computation(double min_probability) const
{
  int nb_parameter = Chain::nb_parameter_computation(min_probability);

  return nb_parameter;
}

/*--------------------------------------------------------------*
 *
 *  Calcul des probabilites de chaque etat en fonction du temps
 *  pour une semi-chaine de Markov_switching.
 *
 *--------------------------------------------------------------*/

void Markov_switching::index_state_distribution()

{
  register int i , j , k;
  double sum , *state_out , **state_in;
  Curves *index_state;
  Parametric *occupancy;


  index_state = nonparametric_process[0]->index_value;

  state_out = new double[nb_state];

  state_in = new double*[index_state->length - 1];
  for (i = 0;i < index_state->length - 1;i++) {
    state_in[i] = new double[nb_state];
  }

  for (i = 0;i < index_state->length;i++) {
    for (j = 0;j < nb_state;j++) {
      switch (state_subtype[j]) {

      // cas etat semi-markovien

      case SEMI_MARKOVIAN : {
        if (i == 0) {
          index_state->point[j][i] = initial[j];
        }
        else {
          index_state->point[j][i] = state_in[i - 1][j] - state_out[j] + index_state->point[j][i - 1];
        }

        if (i < index_state->length - 1) {
          occupancy = nonparametric_process[0]->sojourn_time[j];
          state_out[j] = 0.;
//          istate = 0.;

          for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
            if (k < i + 1) {
              state_out[j] += occupancy->mass[k] * state_in[i - k][j];
//              istate += (1. - occupancy->cumul[k - 1]) * state_in[i - k][j];
            }
            else {
              switch (type) {
              case 'o' :
                state_out[j] += occupancy->mass[k] * initial[j];
//                istate += (1. - occupancy->cumul[k - 1]) * initial[j];
                break;
              case 'e' :
                state_out[j] += forward[j]->mass[k] * initial[j];
//                istate += (1. - forward[j]->cumul[k - 1]) * initial[j];
                break;
              }
            }
          }

//          index_state->point[j][i] = istate;
        }
        break;
      }

      // cas etat markovien

      case MARKOVIAN : {
        if (i == 0) {
          index_state->point[j][i] = initial[j];
        }
        else {
          index_state->point[j][i] = state_in[i - 1][j];
        }

        if (i < index_state->length - 1) {
          state_out[j] = index_state->point[j][i];
        }
        break;
      }
      }
    }

    if (i < index_state->length - 1) {
      for (j = 0;j < nb_state;j++) {
        state_in[i][j] = 0.;
        for (k = 0;k < nb_state;k++) {
          state_in[i][j] += transition[k][j] * state_out[k];
        }
      }
    }

    // renormalisation pour tenir compte des seuils appliques sur
    // les fonctions de repartition des lois d'occupation des etats

    sum = 0.;
    for (j = 0;j < nb_state;j++) {
      sum += index_state->point[j][i];
    }

    if (sum < 1.) {
      for (j = 0;j < nb_state;j++) {
        index_state->point[j][i] /= sum;
      }
    }
  }

  delete [] state_out;

  for (i = 0;i < index_state->length - 1;i++) {
    delete [] state_in[i];
  }
  delete [] state_in;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des probabilites de chaque memoire pour une semi-chaine de Markov_switching.
 *
 *--------------------------------------------------------------*/

double* Markov_switching::memory_computation() const

{
  register int i , j , k;
  double sum , *memory , *state_out , **state_in;
  Parametric *occupancy;


  memory = new double[nb_state];
  state_out = new double[nb_state];

  switch (type) {

  case 'o' : {
    state_in = new double*[nonparametric_process[0]->length->nb_value - 3];
    for (i = 0;i < nonparametric_process[0]->length->nb_value - 3;i++) {
      state_in[i] = new double[nb_state];
    }

    for (i = 0;i < nb_state;i++) {
      memory[i] = 0.;
    }

    for (i = 0;i < nonparametric_process[0]->length->nb_value - 2;i++) {
      for (j = 0;j < nb_state;j++) {
        switch (state_subtype[j]) {

        // cas etat semi-markovien

        case SEMI_MARKOVIAN : {
          occupancy = nonparametric_process[0]->sojourn_time[j];
          state_out[j] = 0.;

          for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
            if (k < i + 1) {
              state_out[j] += occupancy->mass[k] * state_in[i - k][j];
            }
            else {
              state_out[j] += occupancy->mass[k] * initial[j];
            }
          }
          break;
        }

        // cas etat markovien

        case MARKOVIAN : {
          if (i == 0) {
            state_out[j] = initial[j];
          }
          else {
            state_out[j] = state_in[i - 1][j];
          }
          break;
        }
        }

        // accumulation des probabilites des memoires

        memory[j] += state_out[j] * (1. - nonparametric_process[0]->length->cumul[i + 1]);
      }

      if (i < nonparametric_process[0]->length->nb_value - 3) {
        for (j = 0;j < nb_state;j++) {
          state_in[i][j] = 0.;
          for (k = 0;k < nb_state;k++) {
            state_in[i][j] += transition[k][j] * state_out[k];
          }
        }
      }
    }

    for (i = 0;i < nonparametric_process[0]->length->nb_value - 3;i++) {
      delete [] state_in[i];
    }
    delete [] state_in;

    break;
  }

  case 'e' : {
    state_in = new double*[STATIONARY_PROBABILITY_LENGTH];
    for (i = 0;i < STATIONARY_PROBABILITY_LENGTH;i++) {
      state_in[i] = new double[nb_state];
    }

    i = 0;

    do {
      if (i > 0) {
        sum = 0.;
      }

      for (j = 0;j < nb_state;j++) {
        if (i > 0) {
          sum += fabs(state_in[i - 1][j] - state_out[j]);
        }

        switch (state_subtype[j]) {

        // cas etat semi-markovien

        case SEMI_MARKOVIAN : {
          occupancy = nonparametric_process[0]->sojourn_time[j];
          state_out[j] = 0.;

          for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
            if (k < i + 1) {
              state_out[j] += occupancy->mass[k] * state_in[i - k][j];
            }
            else {
              state_out[j] += forward[j]->mass[k] * initial[j];
            }
          }
          break;
        }

        // cas etat markovien

        case MARKOVIAN : {
          if (i == 0) {
            state_out[j] = initial[j];
          }
          else {
            state_out[j] = state_in[i - 1][j];
          }
          break;
        }
        }
      }

      for (j = 0;j < nb_state;j++) {
        state_in[i][j] = 0.;
        for (k = 0;k < nb_state;k++) {
          state_in[i][j] += transition[k][j] * state_out[k];
        }
      }

#     ifdef DEBUG
//      if ((i > 0) && (i % 100 == 0)) {
        cout << i << "  ";
        for (j = 0;j < nb_state;j++) {
          cout << state[j] << " ";
        }
        cout << " | " << sum / nb_state << endl;
//      }
#     endif

      i++;
    }
    while (((i == 1) || (sum / nb_state > STATIONARY_PROBABILITY_THRESHOLD)) &&
           (i < STATIONARY_PROBABILITY_LENGTH));

#   ifdef DEBUG
    cout << "\n" << SEQ_label[SEQL_LENGTH] << ": "  << i << endl;
#   endif

    for (j = 0;j < nb_state;j++) {
      memory[j] = state_in[i - 1][j];
    }

    for (i = 0;i < STATIONARY_PROBABILITY_LENGTH;i++) {
      delete [] state_in[i];
    }
    delete [] state_in;

    break;
  }
  }

  delete [] state_out;

  return memory;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la probabilite de ne pas observer un etat
 *  d'une semi-chaine de Markov_switching ordinaire.
 *
 *  arguments : etat, seuil sur la somme des probabilites de quitter un etat.
 *
 *--------------------------------------------------------------*/

void Markov_switching::state_no_occurrence_probability(int state , double increment)
{
  register int i;

  for (i = 0;i < nb_state;i++) {
    if ((i != state) && (!accessibility[i][state])) {
      break;
    }
  }

  if (i < nb_state) {
    register int j , k;
    int min_time;
    double sum , *state_out , **state_in ,
           &no_occurrence = nonparametric_process[0]->no_occurrence[state];
    Parametric *occupancy;


    state_out = new double[nb_state];

    state_in = new double*[LEAVE_LENGTH];
    state_in[0] = 0;
    for (i = 1;i < LEAVE_LENGTH;i++) {
      state_in[i] = new double[nb_state];
    }

    no_occurrence = 0.;
    for (i = 0;i < nb_state;i++) {
      if ((i != state) && (!accessibility[i][state])) {
        no_occurrence += initial[i];
      }
    }

    sum = 0.;
    for (i = 0;i < nb_state;i++) {
      if (i != state) {
        switch (state_subtype[i]) {

        case SEMI_MARKOVIAN : {
          sum += nonparametric_process[0]->sojourn_time[i]->mean;
          break;
        }

        case MARKOVIAN : {
          if (transition[i][i] < 1.) {
            sum += 1. / (1. - transition[i][i]);
          }
          break;
        }
        }
      }
    }
    min_time = (int)sum + 1;

    i = 1;

    do {

      // calcul des probabilites de quitter (semi-Markov) / d'etre dans (Markov) un etat et
      // mise a jour de la probabilite de ne pas observer l'etat selectionne

      sum = 0.;

      for (j = 0;j < nb_state;j++) {
        if ((j != state) && (accessibility[j][state])) {
          switch (state_subtype[j]) {

          // cas etat semi-markovien

          case SEMI_MARKOVIAN : {
            occupancy = nonparametric_process[0]->sojourn_time[j];
            state_out[j] = 0.;

            for (k = 1;k <= MIN(i , occupancy->nb_value - 1);k++) {
              if (k < i) {
                state_out[j] += occupancy->mass[k] * state_in[i - k][j];
              }
              else {
                state_out[j] += occupancy->mass[k] * initial[j];
              }
            }
            break;
          }

          // cas etat markovien

          case MARKOVIAN : {
            if (i == 1) {
              state_out[j] = initial[j];
            }
            else {
              state_out[j] = state_in[i - 1][j];
            }
            break;
          }
          }

          if ((transition[j][j] == 0.) || (transition[j][j] == 1.)) {
            sum += state_out[j];
          }
          else {
            sum += state_out[j] * (1. - transition[j][j]);
          }

          for (k = 0;k < nb_state;k++) {
            if ((k != state) && (!accessibility[k][state])) {
              no_occurrence += transition[j][k] * state_out[j];
            }
          }
        }
      }

      for (j = 0;j < nb_state;j++) {
        if ((j != state) && (accessibility[j][state])) {
          state_in[i][j] = 0.;
          for (k = 0;k < nb_state;k++) {
            if ((k != state) && (accessibility[k][state])) {
              state_in[i][j] += transition[k][j] * state_out[k];
            }
          }
        }
      }

      i++;
    }
    while (((sum > increment) || (i <= min_time)) && (i < LEAVE_LENGTH));

    delete [] state_out;

    for (i = 1;i < LEAVE_LENGTH;i++) {
      delete [] state_in[i];
    }
    delete [] state_in;
  }
}

/*--------------------------------------------------------------*
 *
 *  Calcul de la loi du temps avant la premiere occurrence d'un etat
 *  pour une semi-chaine de Markov_switching.
 *
 *  arguments : etat, nombre minimum de valeurs,
 *              seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

void Markov_switching::state_first_occurrence_distribution(int state , int min_nb_value ,
							   double cumul_threshold)

{
  register int i , j , k;
  double *state_out , **state_in , *pmass , *pcumul;
  Parametric *occupancy;
  Distribution *first_occurrence;


  first_occurrence = nonparametric_process[0]->first_occurrence[state];
  first_occurrence->complement = nonparametric_process[0]->no_occurrence[state];

  pmass = first_occurrence->mass;
  pcumul = first_occurrence->cumul;

  state_out = new double[nb_state];

  state_in = new double*[first_occurrence->alloc_nb_value];
  state_in[0] = 0;
  for (i = 1;i < first_occurrence->alloc_nb_value;i++) {
    state_in[i] = new double[nb_state];
  }

  *pmass = initial[state];
  *pcumul = *pmass;

  i = 1;

  while (((*pcumul < cumul_threshold - first_occurrence->complement) || (i < min_nb_value)) &&
         (i < first_occurrence->alloc_nb_value)) {

    // calcul des probabilites de quitter (semi-Markov) / d'etre dans (Markov) un etat et
    // de la valeur courante

    *++pmass = 0.;

    for (j = 0;j < nb_state;j++) {
      if (j != state) {
        switch (state_subtype[j]) {

        // cas etat semi-markovien

        case SEMI_MARKOVIAN : {
          occupancy = nonparametric_process[0]->sojourn_time[j];
          state_out[j] = 0.;

          for (k = 1;k <= MIN(i , occupancy->nb_value - 1);k++) {
            if (k < i) {
              state_out[j] += occupancy->mass[k] * state_in[i - k][j];
            }
            else {
              switch (type) {
              case 'o' :
                state_out[j] += occupancy->mass[k] * initial[j];
                break;
              case 'e' :
                state_out[j] += forward[j]->mass[k] * initial[j];
                break;
              }
            }
          }
          break;
        }

        // cas etat markovien

        case MARKOVIAN : {
          if (i == 1) {
            state_out[j] = initial[j];
          }
          else {
            state_out[j] = state_in[i - 1][j];
          }
          break;
        }
        }

        *pmass += transition[j][state] * state_out[j];
      }
    }

    for (j = 0;j < nb_state;j++) {
      if (j != state) {
        state_in[i][j] = 0.;
        for (k = 0;k < nb_state;k++) {
          if (k != state) {
            state_in[i][j] += transition[k][j] * state_out[k];
          }
        }
      }
    }

    // mise a jour de la fonction de repartition

    pcumul++;
    *pcumul = *(pcumul - 1) + *pmass;
    i++;
  }

  first_occurrence->nb_value = i;

  first_occurrence->offset_computation();
  first_occurrence->max_computation();
  first_occurrence->mean_computation();
  first_occurrence->variance_computation();

  delete [] state_out;

  for (i = 1;i < first_occurrence->alloc_nb_value;i++) {
    delete [] state_in[i];
  }
  delete [] state_in;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la probabilite de quitter un etat sans possibilite
 *  d'y revenir pour une semi-chaine de Markov_switching ordinaire.
 *
 *  arguments : etat, seuil sur la somme des probabilites de quitter un etat.
 *
 *--------------------------------------------------------------*/

void Markov_switching::state_leave_probability(int state , double increment)

{
  if (state_type[state] == 't') {
    register int i , j , k;
    int min_time;
    double sum , *state_out , **state_in , &leave = nonparametric_process[0]->leave[state];
    Parametric *occupancy;


    state_out = new double[nb_state];

    state_in = new double*[LEAVE_LENGTH];
    state_in[0] = 0;
    state_in[1] = 0;
    for (i = 2;i < LEAVE_LENGTH;i++) {
      state_in[i] = new double[nb_state];
    }

    leave = 0.;
    for (i = 0;i < nb_state;i++) {
      if ((i != state) && (!accessibility[i][state])) {
        leave += transition[state][i];
      }
    }

    sum = 0.;
    for (i = 0;i < nb_state;i++) {
      if (i != state) {
        switch (state_subtype[i]) {
        case SEMI_MARKOVIAN :
          sum += nonparametric_process[0]->sojourn_time[i]->mean;
          break;
        case MARKOVIAN :
          sum += 1. / (1. - transition[i][i]);
          break;
        }
      }
    }
    min_time = (int)sum + 1;

    i = 2;

    do {

      // calcul des probabilites de quitter (semi-Markov) / d'etre dans (Markov) un etat et
      // mise a jour de la probabilite de quitter l'etat selectionne

      sum = 0.;

      for (j = 0;j < nb_state;j++) {
        if ((j != state) && (accessibility[j][state])) {
          switch (state_subtype[j]) {

          // cas etat semi-markovien

          case SEMI_MARKOVIAN : {
            occupancy = nonparametric_process[0]->sojourn_time[j];
            state_out[j] = 0.;

            for (k = 1;k < MIN(i , occupancy->nb_value);k++) {
              if (k < i - 1) {
                state_out[j] += occupancy->mass[k] * state_in[i - k][j];
              }
              else {
                state_out[j] += occupancy->mass[k] * transition[state][j];
              }
            }
            break;
          }

          // cas etat markovien

          case MARKOVIAN : {
            if (i == 2) {
              state_out[j] = transition[state][j];
            }
            else {
              state_out[j] = state_in[i - 1][j];
            }
            break;
          }
          }

          switch (state_subtype[j]) {
          case SEMI_MARKOVIAN :
            sum += state_out[j];
            break;
          case MARKOVIAN :
            sum += state_out[j] * (1. - transition[j][j]);
            break;
          }

          for (k = 0;k < nb_state;k++) {
            if ((k != state) && (!accessibility[k][state])) {
              leave += transition[j][k] * state_out[j];
            }
          }
        }
      }

      if (transition[state][state] > 0.) {
        sum /= (1. - transition[state][state]);
      }

      for (j = 0;j < nb_state;j++) {
        if ((j != state) && (accessibility[j][state])) {
          state_in[i][j] = 0.;
          for (k = 0;k < nb_state;k++) {
            if ((k != state) && (accessibility[k][state])) {
              state_in[i][j] += transition[k][j] * state_out[k];
            }
          }
        }
      }

      i++;
    }
    while (((sum > increment) || (i <= min_time)) && (i < LEAVE_LENGTH));

    if (state_subtype[state] == SEMI_MARKOVIAN) {
      leave /= nonparametric_process[0]->sojourn_time[state]->parametric_mean_computation();
    }

    delete [] state_out;

    for (i = 2;i < LEAVE_LENGTH;i++) {
      delete [] state_in[i];
    }
    delete [] state_in;
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la loi du temps de retour dans un etat
 *  pour une semi-chaine de Markov_switching.
 *
 *  arguments : etat, nombre minimum de valeurs,
 *              seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

void Markov_switching::state_recurrence_time_distribution(int state , int min_nb_value ,
							  double cumul_threshold)

{
  register int i , j , k;
  double occupancy_mean , *state_out , **state_in , *pmass , *pcumul;
  Distribution *recurrence_time;
  Parametric *occupancy;


  recurrence_time = nonparametric_process[0]->recurrence_time[state];
  recurrence_time->complement = nonparametric_process[0]->leave[state];

  pmass = recurrence_time->mass;
  pcumul = recurrence_time->cumul;
  *pmass = 0.;
  *pcumul = 0.;

  state_out = new double[nb_state];

  state_in = new double*[recurrence_time->alloc_nb_value];
  state_in[0] = 0;
  state_in[1] = 0;
  for (i = 2;i < recurrence_time->alloc_nb_value;i++) {
    state_in[i] = new double[nb_state];
  }

  // calcul de la probabilite de la valeur 1

  switch (state_subtype[state]) {
  case SEMI_MARKOVIAN :
    occupancy_mean = nonparametric_process[0]->sojourn_time[state]->parametric_mean_computation();
    *++pmass = (occupancy_mean - 1.) / occupancy_mean;
    break;
  case MARKOVIAN :
    *++pmass = transition[state][state];
    break;
  }

  *++pcumul = *pmass;

  i = 2;

  while (((*pcumul < cumul_threshold - recurrence_time->complement) || (i < min_nb_value)) &&
         (i < recurrence_time->alloc_nb_value)) {

    // calcul des probabilites de quitter (semi-Markov) / d'etre dans (Markov) un etat et
    // de la valeur courante

    *++pmass = 0.;

    for (j = 0;j < nb_state;j++) {
      if (j != state) {
        switch (state_subtype[j]) {

        // cas etat semi-markovien

        case SEMI_MARKOVIAN : {
          occupancy = nonparametric_process[0]->sojourn_time[j];
          state_out[j] = 0.;

          for (k = 1;k < MIN(i , occupancy->nb_value);k++) {
            if (k < i - 1) {
              state_out[j] += occupancy->mass[k] * state_in[i - k][j];
            }
            else {
              state_out[j] += occupancy->mass[k] * transition[state][j];
            }
          }
          break;
        }

        // cas etat markovien

        case MARKOVIAN : {
          if (i == 2) {
            state_out[j] = transition[state][j];
          }
          else {
            state_out[j] = state_in[i - 1][j];
          }
          break;
        }
        }

        *pmass += transition[j][state] * state_out[j];
      }
    }

    for (j = 0;j < nb_state;j++) {
      if (j != state) {
        state_in[i][j] = 0.;
        for (k = 0;k < nb_state;k++) {
          if (k != state) {
            state_in[i][j] += transition[k][j] * state_out[k];
          }
        }
      }
    }

    if (state_subtype[state] == SEMI_MARKOVIAN) {
      *pmass /= occupancy_mean;
    }
    pcumul++;
    *pcumul = *(pcumul - 1) + *pmass;
    i++;
  }

  recurrence_time->nb_value = i;
  recurrence_time->nb_value_computation();

  delete [] state_out;

  for (i = 2;i < recurrence_time->alloc_nb_value;i++) {
    delete [] state_in[i];
  }
  delete [] state_in;

  if (recurrence_time->nb_value > 0) {
    recurrence_time->offset_computation();
    recurrence_time->max_computation();
    recurrence_time->mean_computation();
    recurrence_time->variance_computation();
  }

  else {
    delete nonparametric_process[0]->recurrence_time[state];
    nonparametric_process[0]->recurrence_time[state] = 0;
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul d'un melange de lois du nombre de series d'un etat ('r') ou
 *  du nombre d'occurrences d'un etat ('o') d'une semi-chaine de Markov-switching
 *  pour une distribution des longueurs de sequences donnee.
 *
 *  arguments : etat, type de forme.
 *
 *--------------------------------------------------------------*/

void Markov_switching::state_nb_pattern_mixture(int state , char pattern)

{
  register int i , j , k , m;
  int max_length , index_nb_pattern , previous_nb_pattern , increment;
  double sum , *pmass , *lmass , **state_out , *pstate_out , ***state_in;
  Distribution *pdist;
  Parametric *occupancy;


  switch (pattern) {
  case 'r' :
    pdist = nonparametric_process[0]->nb_run[state];
    break;
  case 'o' :
    pdist = nonparametric_process[0]->nb_occurrence[state];
    break;
  }

  pmass = pdist->mass;
  for (i = 0;i < pdist->nb_value;i++) {
    *pmass++ = 0.;
  }

  max_length = nonparametric_process[0]->length->nb_value - 1;

  state_out = new double*[nb_state];
  for (i = 0;i < nb_state;i++) {
    state_out[i] = new double[pattern == 'o' ? max_length : (max_length + 1) / 2 + 1];
  }

  state_in = new double**[max_length - 1];
  index_nb_pattern = 1;

  for (i = 0;i < max_length - 1;i++) {
    state_in[i] = new double*[nb_state];
    for (j = 0;j < nb_state;j++) {
      state_in[i][j] = new double[index_nb_pattern + 1];
    }
    if ((pattern == 'o') || (i % 2 == 1)) {
      index_nb_pattern++;
    }
  }

  // calcul des probabilites de quitter (semi-Markov) / d'etre dans (Markov) un etat
  // fonction du nombre de formes de l'etat selectionne

  lmass = nonparametric_process[0]->length->mass;
  index_nb_pattern = 1;

  for (i = 0;i < max_length;i++) {
    lmass++;

    for (j = 0;j < nb_state;j++) {

      // initialisation des probabilites de quitter un etat a l'instant i

      if (i < max_length - 1) {
        pstate_out = state_out[j];
        for (k = 0;k <= index_nb_pattern;k++) {
          *pstate_out++ = 0.;
        }
      }

      switch (state_subtype[j]) {

      // cas etat semi-markovien

      case SEMI_MARKOVIAN : {
        occupancy = nonparametric_process[0]->sojourn_time[j];

        for (k = (*lmass > 0. ? 1 : occupancy->offset);k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
          switch (pattern) {
          case 'r' :
            increment = 1;
            break;
          case 'o' :
            increment = k;
            break;
          }

          if (i < max_length - 1) {
            pstate_out = state_out[j];
            if (j == state) {
              pstate_out += increment;
            }
          }
          if (*lmass > 0.) {
            pmass = pdist->mass;
            if (j == state) {
              pmass += increment;
            }
          }

          if (k < i + 1) {
            switch (pattern) {

            case 'r' : {
              if ((j == state) && (k == 1) && (i % 2 == 1)) {
                previous_nb_pattern = index_nb_pattern - 1;
              }
              else {
                previous_nb_pattern = (i - k) / 2 + 1;
              }
              break;
            }

            case 'o' : {
              previous_nb_pattern = i - k + 1;
              break;
            }
            }

            if (i < max_length - 1) {
              for (m = 0;m <= previous_nb_pattern;m++) {
                *pstate_out++ += occupancy->mass[k] * state_in[i - k][j][m];
              }
            }
            if (*lmass > 0.) {
              for (m = 0;m <= previous_nb_pattern;m++) {
                *pmass++ += *lmass * (1. - occupancy->cumul[k - 1]) * state_in[i - k][j][m];
              }
            }
          }

          else {
            if (i < max_length - 1) {
              switch (type) {
              case 'o' :
                *pstate_out += occupancy->mass[k] * initial[j];
                break;
              case 'e' :
                *pstate_out += forward[j]->mass[k] * initial[j];
                break;
              }
            }

            if (*lmass > 0.) {
              switch (type) {
              case 'o' :
                *pmass += *lmass * (1. - occupancy->cumul[k - 1]) * initial[j];
                break;
              case 'e' :
                *pmass += *lmass * (1. - forward[j]->cumul[k - 1]) * initial[j];
                break;
              }
            }
          }
        }
        break;
      }

      // cas etat markovien

      case MARKOVIAN : {
        if (i < max_length - 1) {
          pstate_out = state_out[j];
          if (j == state) {
            pstate_out++;
          }
        }

        if (*lmass > 0.) {
          pmass = pdist->mass;
          if (j == state) {
            pmass++;
          }
        }

        if (i == 0) {
          *pstate_out = initial[j];
          if (*lmass > 0.) {
            *pmass += *lmass * initial[j];
          }
        }

        else {
          switch (pattern) {

          case 'r' : {
            if ((j == state) && (i % 2 == 1)) {
              previous_nb_pattern = index_nb_pattern - 1;
            }
            else {
              previous_nb_pattern = (i - 1) / 2 + 1;
            }
            break;
          }
  
          case 'o' : {
            previous_nb_pattern = i;
            break;
          }
          }

          if (i < max_length - 1) {
            for (k = 0;k <= previous_nb_pattern;k++) {
              *pstate_out++ = state_in[i - 1][j][k];
            }
          }
          if (*lmass > 0.) {
            for (k = 0;k <= previous_nb_pattern;k++) {
              *pmass++ += *lmass * state_in[i - 1][j][k];
            }
          }
        }
        break;
      }
      }
    }

    if (i < max_length - 1) {
      for (j = 0;j < nb_state;j++) {
        for (k = 0;k <= index_nb_pattern;k++) {
          state_in[i][j][k] = 0.;
          for (m = 0;m < nb_state;m++) {
            if ((pattern == 'o') || (j != state) || (j != m)) {
              state_in[i][j][k] += transition[m][j] * state_out[m][k];
            }
            else if (k < index_nb_pattern) {
              state_in[i][j][k] += transition[m][j] * state_out[m][k + 1];
            }
          }
        }
      }
    }

    if ((pattern == 'o') || (i % 2 == 1)) {
      index_nb_pattern++;
    }
  }

  // renormalisation du melange de lois du nombre de formes de l'etat
  // selectionne pour tenir compte des seuils appliques sur les fonctions
  // de repartition des lois d'occupation des etats

  pmass = pdist->mass;
  sum = 0.;
  for (i = 0;i < pdist->nb_value;i++) {
    sum += *pmass++;
  }

  if (sum < 1.) {
    pmass = pdist->mass;
    for (i = 0;i < pdist->nb_value;i++) {
      *pmass++ /= sum;
    }
  }

  pdist->nb_value_computation();
  pdist->offset_computation();
  pdist->cumul_computation();

  pdist->max_computation();
  pdist->mean_computation();
  pdist->variance_computation();

  for (i = 0;i < nb_state;i++) {
    delete [] state_out[i];
  }
  delete [] state_out;

  for (i = 0;i < max_length - 1;i++) {
    for (j = 0;j < nb_state;j++) {
      delete [] state_in[i][j];
    }
    delete [] state_in[i];
  }
  delete [] state_in;
}



/*--------------------------------------------------------------*
 *
 *  Calcul des lois caracteristiques d'un objet Semi_markov_switching.
 *
 *  arguments : longueur des sequences, flag sur le calcul des lois de comptage,
 *              indice du processus d'observation.
 *
 *--------------------------------------------------------------*/

void Markov_switching::characteristic_computation(int length , bool counting_flag , int variable)

{
  if (nb_component > 0) {
    bool computation[NB_OUTPUT_PROCESS + 1];
    register int i , j , k;
    double *memory;
    Parametric dlength(UNIFORM , length , length , D_DEFAULT , D_DEFAULT);


    memory = 0;

    // calcul des lois de type intensite et intervalle au niveau etat

    if (((variable == I_DEFAULT) || (variable == 0)) &&
        ((!(nonparametric_process[0]->length)) ||
         (dlength != *(nonparametric_process[0]->length)))) {
      computation[0] = true;
      nonparametric_process[0]->create_characteristic(dlength , false , counting_flag);

      index_state_distribution();

      for (i = 0;i < nb_state;i++) {
        if (type == 'o') {
          state_no_occurrence_probability(i);
        }
        state_first_occurrence_distribution(i);

        if (type == 'o') {
          state_leave_probability(i);
        }
        if (nonparametric_process[0]->leave[i] < 1.) {
          state_recurrence_time_distribution(i);
        }
        else {
          delete nonparametric_process[0]->recurrence_time[i];
          nonparametric_process[0]->recurrence_time[i] = 0;
        }

        if ((state_subtype[i] == MARKOVIAN) && (transition[i][i] < 1.)) {
          if (transition[i][i] > 0.) {
            nonparametric_process[0]->sojourn_time[i] = new Parametric(NEGATIVE_BINOMIAL , 1 , I_DEFAULT , 1. ,
                                                                       1. - transition[i][i] , OCCUPANCY_THRESHOLD);
            nonparametric_process[0]->sojourn_time[i]->parameter = D_DEFAULT;
            nonparametric_process[0]->sojourn_time[i]->probability = D_DEFAULT;
          }

          else {
            nonparametric_process[0]->sojourn_time[i] = new Parametric(UNIFORM , 1 , 1 , D_DEFAULT , D_DEFAULT);
            nonparametric_process[0]->sojourn_time[i]->sup_bound = I_DEFAULT;
          }

          nonparametric_process[0]->sojourn_time[i]->ident = NONPARAMETRIC;
          nonparametric_process[0]->sojourn_time[i]->inf_bound = I_DEFAULT;
        }
      }

#     ifdef MESSAGE
      if (type == 'e') {
        double sum = 0.;

        // calcul de la loi stationnaire dans le cas d'un processus en equilibre
        // avec renormalisation pour tenir compte des seuils appliques sur
        // les fonctions de repartition des lois de temps de retour dans les etats

        for (i = 0;i < nb_state;i++) {
          sum += 1. / nonparametric_process[0]->recurrence_time[i]->mean;
        }

        cout << "\n" << STAT_label[STATL_STATIONARY_PROBABILITIES] << endl;
        for (i = 0;i < nb_state;i++) {
          cout << initial[i] << " | "
               << 1. / (nonparametric_process[0]->recurrence_time[i]->mean * sum) << endl;
        }
      }
#     endif

    }

    else {
      computation[0] = false;
    }

    if (counting_flag) {
      
      // calcul des lois de comptage au niveau etat
      
      if (computation[0]) {
        for (i = 0;i < nb_state;i++) {
          state_nb_pattern_mixture(i , 'r');
          state_nb_pattern_mixture(i , 'o');
        }
      }
      
    }
  }
}

/*--------------------------------------------------------------*
 *
 *  Calcul des lois caracteristiques d'un objet Semi_markov_switching.
 *
 *  arguments : reference sur un objet Markov_switching_data,
 *              flag sur le calcul des lois de comptage,
 *              indice du processus d'observation,
 *              flag pour tenir compte des longueurs.
 *
 *--------------------------------------------------------------*/

void Markov_switching::characteristic_computation(const Markov_switching_data &seq , bool counting_flag ,
						  int variable , bool length_flag)

{
  if (nb_component > 0) {
    bool computation[NB_OUTPUT_PROCESS + 1];
    register int i , j , k;
    int seq_variable;
    double *memory;
    Distribution dlength(*(seq.hlength));


    memory = 0;

    // calcul des lois de type intensite et intervalle au niveau etat

    if (((variable == I_DEFAULT) || (variable == 0)) && ((!length_flag) ||
         ((length_flag) && ((!(nonparametric_process[0]->length)) ||
           (dlength != *(nonparametric_process[0]->length)))))) {
      computation[0] = true;
      nonparametric_process[0]->create_characteristic(dlength , false , counting_flag);

      index_state_distribution();

      for (i = 0;i < nb_state;i++) {
        if (type == 'o') {
          state_no_occurrence_probability(i);
        }
        if (seq.type[0] == STATE) {
          state_first_occurrence_distribution(i , (seq.characteristics[0] ? seq.characteristics[0]->first_occurrence[i]->nb_value : 1));
        }
        else {
          state_first_occurrence_distribution(i);
        }

        if (type == 'o') {
          state_leave_probability(i);
        }
        if (nonparametric_process[0]->leave[i] < 1.) {
          if (seq.type[0] == STATE) {
            state_recurrence_time_distribution(i , (seq.characteristics[0] ? seq.characteristics[0]->recurrence_time[i]->nb_value : 1));
          }
          else {
            state_recurrence_time_distribution(i);
          }
        }
        else {
          delete nonparametric_process[0]->recurrence_time[i];
          nonparametric_process[0]->recurrence_time[i] = 0;
        }

        if ((state_subtype[i] == MARKOVIAN) && (transition[i][i] < 1.)) {
          if (transition[i][i] > 0.) {
            nonparametric_process[0]->sojourn_time[i] = new Parametric(NEGATIVE_BINOMIAL , 1 , I_DEFAULT , 1. ,
                                                                       1. - transition[i][i] , OCCUPANCY_THRESHOLD);

            if ((seq.type[0] == STATE) && (seq.characteristics[0]) &&
                (seq.characteristics[0]->sojourn_time[i]->nb_value > nonparametric_process[0]->sojourn_time[i]->nb_value)) {
              nonparametric_process[0]->sojourn_time[i]->computation(seq.characteristics[0]->sojourn_time[i]->nb_value , OCCUPANCY_THRESHOLD);
            }
            nonparametric_process[0]->sojourn_time[i]->parameter = D_DEFAULT;
            nonparametric_process[0]->sojourn_time[i]->probability = D_DEFAULT;
          }

          else {
            nonparametric_process[0]->sojourn_time[i] = new Parametric(UNIFORM , 1 , 1 , D_DEFAULT , D_DEFAULT);
            nonparametric_process[0]->sojourn_time[i]->sup_bound = I_DEFAULT;
          }

          nonparametric_process[0]->sojourn_time[i]->ident = NONPARAMETRIC;
          nonparametric_process[0]->sojourn_time[i]->inf_bound = I_DEFAULT;
        }
      }

#     ifdef MESSAGE
      if (type == 'e') {
        double sum = 0.;

        // calcul de la loi stationnaire dans le cas d'un processus en equilibre
        // avec renormalisation pour tenir compte des seuils appliques sur
        // les fonctions de repartition des lois de temps de retour dans les etats

        for (i = 0;i < nb_state;i++) {
          sum += 1. / nonparametric_process[0]->recurrence_time[i]->mean;
        }

        cout << "\n" << STAT_label[STATL_STATIONARY_PROBABILITIES] << endl;
        for (i = 0;i < nb_state;i++) {
          cout << initial[i] << " | "
               << 1. / (nonparametric_process[0]->recurrence_time[i]->mean * sum) << endl;
        }
      }
#     endif

    }

    else {
      computation[0] = false;
    }

   
    if (counting_flag) {

      // calcul des lois de comptage au niveau etat

      if (computation[0]) {
        for (i = 0;i < nb_state;i++) {
          state_nb_pattern_mixture(i , 'r');
          state_nb_pattern_mixture(i , 'o');
        }
      }
    }
  }
}



