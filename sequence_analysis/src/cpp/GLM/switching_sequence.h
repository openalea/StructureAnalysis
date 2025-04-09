#ifndef SWITCHING_SEQUENCE_H
#define SWITCHING_SEQUENCE_H

#include <fstream>

enum{
  WITH_CONSTANT,
  WITHOUT_CONSTANT
};


class Markov_switching;

class Switching_sequence : public Sequences {//trajectoires correspondant à un processus markovien

  friend std::ostream& operator<<(std::ostream &os, const Switching_sequence &sw_seq)
  { return sw_seq.ascii_write(os); }
  
  friend class Sequence_characteristics;
  friend class Markov_switching;
  
 public:
  int constant; // presence ou absence de constante pour chaque variable
  int nb_covariable; // nombre de covariables
  Histogram ***observation; //histogrammes correspondant aux lois d'observations
  double ***covar; // Covariables (voir si double ou int) 
  Self_transition **self_transition; // probabilité de rester dans un état en fonction de l'index
  Sequence_characteristics **characteristics; // caractéristiques pour la variable d'état (dans lequel se trouve l'individu à un temps t)

  Switching_sequence();
  Switching_sequence(int inb_variable, int *itype, int inb_sequence, int *iidentifier, int *ilength,
		     int inb_covariable, int iconstant = WITHOUT_CONSTANT, bool init_flag = true);
  Switching_sequence(int inb_variable, int inb_sequence, int *iidentifier, int *ilength, int inb_covariable,
		     int iconstant = WITHOUT_CONSTANT,bool init_flag = true);
  Switching_sequence(int inb_variable, const Histogram &ihlength, int inb_covariable, int iconstant = WITHOUT_CONSTANT,
		     bool init_flag = true);
  Switching_sequence(int inb_variable, int inb_sequence, int inb_covariable, int *ilength, int ***isequence,
		     double ***icovar, int *iidentifier, int iconstant = WITHOUT_CONSTANT);
  Switching_sequence(const Sequences &seq, int inb_covariable, double ***icovar, int iconstant = WITHOUT_CONSTANT);
  Switching_sequence(const Switching_sequence &sw_seq , char transform ='c' , int param1 = DEFAULT , int param2 = DEFAULT);

  ~Switching_sequence();

  void remove();

  void covar_computation();

  void init();

  void copy(const Switching_sequence &sw_seq , int param = DEFAULT);

  void add_variable(const Switching_sequence &sw_seq , int variable , int param);

  Switching_sequence& operator=(const Switching_sequence &sw_seq);

  void state_variable_init(int itype = STATE);

  Switching_sequence* remove_variable_1() const;

  double iid_information_computation() const;


  void build_index_value(int variable);
  void build_first_occurrence_histogram(int variable);
  void build_recurrence_time_histogram(int variable);
  void sojourn_time_histogram_computation(int variable);
  void build_sojourn_time_histogram(int variable , int initial_run_flag = false);
  void build_nb_run_histogram(int variable);
  void build_nb_occurrence_histogram(int variable);
  void transition_count_computation(const Chain_data &chain_data , int order , bool begin = true) const;

  void observation_histogram_computation(int variable);
  void observation_histogram_computation();
  void create_observation_histogram(int nb_state);
  void build_observation_histogram();

  void build_characteristic(int variable = I_DEFAULT , bool sojourn_time_flag = true ,
			    bool initial_run_flag = false);

  ostream& ascii_write(ostream &os , bool exhaustive , bool comment_flag) const;

  ostream& ascii_write(ostream &os , bool exhaustive = false) const;

  bool ascii_write(Format_error &error , const char *path , bool exhaustive = false) const;

  ostream& ascii_data_write(ostream &os , char format = 'c', bool exhaustive = false) const;

  bool ascii_data_write(Format_error &error , const char *path , char format = 'c', bool exhaustive = false) const;

  bool spreadsheet_write(Format_error &error , const char *path) const;


  Markov_switching* markov_switching_estimation(Format_error &error , ostream &os ,
						const Markov_switching &isw_markov ,
						int estimator = COMPLETE_LIKELIHOOD,
						bool output_R = false,
						bool counting_flag = true , bool state_sequence = true,
						int nb_iter = I_DEFAULT, int mean_computation = COMPUTED) const;

  Markov_switching* markov_switching_estimation_multi(Format_error &error , ostream &os ,
						      const Markov_switching &isw_markov , int nb_cat,
						      int estimator = COMPLETE_LIKELIHOOD,
						      bool output_R = false,
						      bool counting_flag = true , bool state_sequence = true,
						      int nb_iter = I_DEFAULT, int mean_computation = COMPUTED) const;
  

  Markov_switching* markov_switching_estimation(Format_error &error , ostream &os ,
						int nb_state , bool left_right ,
						double **regression, double ***regmult,
						int nb_covariable, int nb_cat, int constant, int ident, int link,
						bool markov = true, bool output_R = false,
						int order = 1, bool counting_flag = true , bool state_sequence = true,
						double occupancy_mean = D_DEFAULT,
						double self_transition = D_DEFAULT, int nb_iter = I_DEFAULT,
						int mean_computation = COMPUTED) const;


  Markov_switching* markov_switching_stochastic_estimation(Format_error &error , std::ostream &os ,
							   const Markov_switching &ihsw_markov ,
							   int min_nb_state_sequence = MIN_NB_STATE_SEQUENCE ,
							   int max_nb_state_sequence = MAX_NB_STATE_SEQUENCE ,
							   bool output_R = false,
							   double parameter = NB_STATE_SEQUENCE_PARAMETER ,
							   int estimator = COMPLETE_LIKELIHOOD,
							   bool counting_flag = true ,
							   bool state_sequence = true  ,
							   int nb_iter = I_DEFAULT) const;

  Markov_switching* markov_switching_stochastic_estimation_multi(Format_error &error , std::ostream &os ,
								 const Markov_switching &ihsw_markov , int nb_cat,
								 int min_nb_state_sequence = MIN_NB_STATE_SEQUENCE ,
								 int max_nb_state_sequence = MAX_NB_STATE_SEQUENCE ,
								 bool output_R = false,
								 double parameter = NB_STATE_SEQUENCE_PARAMETER ,
								 int estimator = COMPLETE_LIKELIHOOD,
								 bool counting_flag = true ,
								 bool state_sequence = true  ,
								 int nb_iter = I_DEFAULT) const;

  Markov_switching* markov_switching_stochastic_estimation(Format_error &error , std::ostream &os ,
							   char type , int nb_state , bool left_right ,
							   double **regression, double ***regmult,
							   int nb_covariable, int nb_cat,
							   int constant, int ident, int link,
							   int order,
							   int min_nb_state_sequence = MIN_NB_STATE_SEQUENCE ,
							   int max_nb_state_sequence = MAX_NB_STATE_SEQUENCE ,
							   bool markov = true,
							   bool output_R = false,
							   double parameter = NB_STATE_SEQUENCE_PARAMETER ,
							   bool counting_flag = true ,
							   bool state_sequence = true  ,
							   double self_transition = D_DEFAULT ,
							   double occupancy_mean = D_DEFAULT,
							   int nb_iter = I_DEFAULT) const;

  //private:

  //protected:

};

#endif
