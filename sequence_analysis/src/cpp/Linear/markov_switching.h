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


#ifndef MARKOV_SWITCHING_H
#define MARKOV_SWITCHING_H


const int LEAVE_LENGTH = 10000;                        // longueur maximum pour le calcul de la probabilité de quitter un état

const double OCCUPANCY_LIKELIHODD_DIFF = 1.e-5;        // seuil pour stopper les iterations EM
const int OCCUPANCY_NB_ITER = 10000;                   // nombre maximum d'itérations EM
const int OCCUPANCY_COEFF = 10;                        // coefficient arrondi estimateur pour les lois d'occupation des états

const double MARKOV_SWITCHING_LIKELIHOOD_DIFF = 1.e-6; // seuil pour stopper les itérations EM
const int MARKOV_SWITCHING_NB_ITER = 100;              // nombre d'itérations maximums EM
const int MARKOV_SWITCHING_NB_ITER_MH = 1000;          // nombre d'itérations pour le Metropolis Hastings

const double SEMI_MARKOV_LIKELIHOOD_DIFF = 1.e-6;      // seuil pour stopper les iterations EM
const int EXPLORATION_NB_ITER = 10;                    // nombre d'iterations de la phase d'exploration
const int STOCHASTIC_EXPLORATION_NB_ITER = 5;          // nombre d'iterations de la phase d'exploration
const int SEMI_MARKOV_NB_ITER = 500;                   // nombre maximum d'iterations EM

const double MIN_SMOOTHED_PROBABILITY = 1.e-3;         // seuil sur les probabilites lissees


enum {
  SSTATE ,                                             // etat
  IN_STATE ,                                           // entree etat
  OUT_STATE                                            // sortie etat
};

enum {
  MARKOVIAN,
  SEMI_MARKOVIAN
};


/*----------------------------------------------------------------------
 *
 * Classe Markov_switching
 * Herite de STAT_interface et de Chain
 *
 *----------------------------------------------------------------------*/

class Switching_sequence;
class Markov_switching_data;

class Markov_switching: public STAT_interface, public Chain{

  friend class Switching_sequence;
  friend class Markov_switching_data;

  friend Markov_switching* markov_switching_ascii_read(Format_error &error, 
						       const char *path, 
						       int length = DEFAULT_LENGTH){}

  friend std::ostream& operator<<(std::ostream &os, const Markov_switching &isw_markov)
  { return isw_markov.ascii_write(os);}

 public:

  Markov_switching_data *sw_markov_data;                  // pointer on Markov_switching_data object
  Switching_process **sw_process;                          // output process : linear mixed model
  int *state_subtype;                                     // MARKOVIAN/SEMI_MARKOVIAN
  Forward **forward;                                      // distribution of residual time interval
  Nonparametric_sequence_process **nonparametric_process; // state process  

 
  void copy(const Markov_switching &isw_markov, bool data_flag = true, int param = I_DEFAULT);

  void remove();

  std::ostream& ascii_write(std::ostream &os, const Markov_switching_data *isw_markov_data, 
			    bool exhaustive = false, bool file_fag = false) const;
  std::ostream& spreadsheet_write(std::ostream &os, const Markov_switching_data *isw_markov_data,
				  bool hidden = false) const;

  int nb_parameter_computation(double min_probability = 0.) const;
  int nb_transient_parameter_computation(double min_probability = 0.) const;

  void index_state_distribution();
  double* memory_computation() const;
  void state_no_occurrence_probability(int state , double increment = LEAVE_INCREMENT);
  void state_first_occurrence_distribution(int state , int min_nb_value = 1 , 
					   double cumul_threshold = CUMUL_THRESHOLD);
  void state_leave_probability(int state , double increment = LEAVE_INCREMENT);
  void state_recurrence_time_distribution(int state , int min_nb_value = 1,
					  double cumul_threshold = CUMUL_THRESHOLD);
  void state_nb_pattern_mixture(int state , char pattern);


  // commun à LM et aux LMM

  double viterbi(const Markov_switching_data &isw_markov_data,
		 double *posterior_probability, bool fag = true, 
		 bool output_R = false, int index = I_DEFAULT) const ;

  double forward_backward(const Markov_switching_data &seq , int index ,
			  ostream &os , int output , char format,
			  double &max_marginal_entropy , double &entropy1) const;
    
  double forward_backward_sampling(const Markov_switching_data &seq, int index,
				   std::ostream &os, char format = 'a',
				   int nb_state_sequence = NB_STATE_SEQUENCE) const;
  

  bool state_profile_write(Format_error &error, std::ostream &os,
			   const Markov_switching_data &iseq,
			   int output, int identifier = I_DEFAULT, char format = 'a',
			   int state_sequence = GENERALIZED_VITERBI,
			   int nb_state_sequence = NB_STATE_SEQUENCE) const;

  bool state_profile_ascii_write(Format_error &error , std::ostream &os ,
				 int output, int identifier = I_DEFAULT ,
				 int state_sequence = GENERALIZED_VITERBI ,
				 int nb_state_sequence = NB_STATE_SEQUENCE) const;

  bool state_profile_write(Format_error &error , const char *path ,
			   int output, int identifier = I_DEFAULT , char format = 'a' ,
			   int state_sequence = GENERALIZED_VITERBI,
			   int nb_state_sequence = NB_STATE_SEQUENCE) const;

  bool state_profile_plot_write(Format_error &error, const char *prefix,
				const Markov_switching_data &iseq,
				int output, int identifier, const char *title = 0) const;
  
  bool state_profile_plot_write(Format_error &error , const char *prefix ,
				int output, int identifier , const char *title) const;



  double generalized_viterbi(const Markov_switching_data &seq, int index,
			     std::ostream &os, double seq_likelihood, char format,
			     int inb_state_sequence, bool output_R = false, bool posterior_proba = false) const;

  double viterbi_forward_backward(const Markov_switching_data &seq, int index,
				  std::ostream &os, int output, char format,
				  double seq_likelihood = D_INF) const;

  void log_computation();  


  Markov_switching();
  Markov_switching(int inb_state);
  Markov_switching(int inb_state, int *nb_value);
  Markov_switching(const Chain *pchain, const Nonparametric_sequence_process *poccupancy,
		   const Switching_process *pobservation , int length);
  Markov_switching(const Markov_switching &isw_markov);
  Markov_switching(const Markov_switching &isw_markov, bool data_flag = true, int param = I_DEFAULT)
      : Chain(isw_markov) { copy(isw_markov, data_flag, param);}
  virtual ~Markov_switching();

  Markov_switching& operator=(const Markov_switching &isw_markov);

  Markov_switching_data* extract_data(Format_error &error) const;
  Markov_switching* thresholding(double min_probability = MIN_PROBABILITY) const;

  std::ostream& line_write(std::ostream &os) const ;
  std::ostream& ascii_write(std::ostream &os, bool exhaustive = false) const;
  bool ascii_write(Format_error &error, const char *path, bool exhaustive = false) const;
  bool spreadsheet_write(Format_error &error, const char *prefix) const;

  bool plot_write(Format_error &error, const char *prefix, const char *title = 0) const {} //pas implémenté

  void characteristic_computation(int length, bool counting_flag, int variable = I_DEFAULT);
  void characteristic_computation(const Markov_switching_data &seq, bool counting_flag,
				  int variable = I_DEFAULT, bool length_flag = true);

  double likelihood_computation(const Switching_sequence &isw_seq, int index = I_DEFAULT) const;
  double likelihood_computation(const Markov_switching_data &isw_markov_data) const;
  double likelihood_computation(const Switching_sequence &isw_seq, double *posterior_probability,
				int index) const; 


  //simulation LMM with inter-individual heterogeneity 

  Markov_switching_data* simulation_hetero_effect(Format_error &error , const Histogram &hlength ,
						  int inb_covariable, int inb_random,
						  int iconstant, double *iresidual_variance,
						  double **iregression, double **irandom_variance,
						  double ***icovar, double **ieffect,
						  bool counting_flag , bool divergence_flag) const;

  Markov_switching_data* simulation_hetero_effect(Format_error &error , int nb_sequence ,
						  int length , int inb_covariable, int inb_random,
						  int iconstant,  double *iresidual_variance,
						  double **iregression, double **irandom_variance,
						  double ***icovar, double **ieffect, bool counting_flag) const;

  Markov_switching_data* simulation_hetero_effect(Format_error &error , int nb_sequence ,
						  const Switching_sequence &iseq ,
						  int inb_covariable, int inb_random, int iconstant,
						  double *iresidual_variance, double **iregression,
						  double **irandom_variance, double ***icovar,
						  double **ieffect, bool counting_flag) const;

  Markov_switching_data* simulation_hetero_effect(Format_error &error ,
						  const Histogram &hlength ,
						  int inb_covariable, int inb_random, int iconstant,
						  double *iresidual_variance, double **iregression,
						  double **irandom_variance, double ***icovar,
						  double **ieffect, bool counting_flag ,
						  bool divergence_flag, bool hidden) const;

  Markov_switching_data* simulation_hetero_effect(Format_error &error ,
						  int nb_sequence , int length ,
						  int inb_covariable, int inb_random, int iconstant,
						  double *iresidual_variance, double **iregression,
						  double **irandom_variance, double ***icovar,
						  double **ieffect,
						  bool counting_flag, bool hidden) const;

  Markov_switching_data* simulation_hetero_effect(Format_error &error ,
						  int nb_sequence ,
						  const Switching_sequence &iseq ,
						  int inb_covariable, int inb_random, int iconstant,
						  double *iresidual_variance, double **iregression,
						  double **irandom_variance, double ***icovar,
						  double **ieffect,
						  bool counting_flag, bool hidden) const;

  //simulation LMM with Common environment effect

  Markov_switching_data* simulation_year_effect(Format_error &error , const Histogram &hlength ,
						int inb_covariable, int inb_random,
						int iconstant, double *iresidual_variance,
						double **iregression, double **irandom_variance,
						double ***icovar, double *iyear_effect, int **iindex, int T,
						bool counting_flag , bool divergence_flag) const;

  Markov_switching_data* simulation_year_effect(Format_error &error , int nb_sequence ,
						int length , int inb_covariable, int inb_random,
						int iconstant, double *iresidual_variance,
						double **iregression, double **irandom_variance,
						double ***icovar,  double *iyear_effect, int **iindex, int T,  
						bool counting_flag) const;

  Markov_switching_data* simulation_year_effect(Format_error &error , int nb_sequence ,
						const Switching_sequence &iseq ,
						int inb_covariable, int inb_random,
						int iconstant, double *iresidual_variance,
						double **iregression, double **irandom_variance,
						double ***icovar,  double *iyear_effect, int **iindex, int T,
						bool counting_flag) const;

  Markov_switching_data* simulation_year_effect(Format_error &error ,
						const Histogram &hlength ,
						int inb_covariable, int inb_random,
						int iconstant, double *iresidual_variance,
						double **iregression, double **irandom_variance,
						double ***icovar,  double *iyear_effect, int **iindex, int T,
						bool counting_flag  , bool divergence_flag ,
						bool hidden ) const;

  Markov_switching_data* simulation_year_effect(Format_error &error ,
						int nb_sequence , int length ,
						int inb_covariable, int inb_random,
						int iconstant, double *iresidual_variance,
						double **iregression, double **irandom_variance,
						double ***icovar,  double *iyear_effect, int **iindex,int T,
						bool counting_flag , bool hidden ) const;

  Markov_switching_data* simulation_year_effect(Format_error &error ,
						int nb_sequence ,
						const Switching_sequence &iseq ,
						int inb_covariable, int inb_random,
						int iconstant, double *iresidual_variance,
						double **iregression, double **irandom_variance,
						double ***icovar,  double *iyear_effect, int **iindex,int T,
						bool counting_flag , bool hidden ) const;


// simulation LM  

  Markov_switching_data* simulation(Format_error &error , const Histogram &hlength ,
				    int inb_covariable, int iconstant,
				    double *iresidual_variance, double **iregression,
				    double ***icovar, bool counting_flag , bool divergence_flag) const;

  Markov_switching_data* simulation(Format_error &error , int nb_sequence ,
				    int length , int inb_covariable, int iconstant,
				    double *iresidual_variance, double **iregression,
				    double ***icovar, bool counting_flag) const;

  Markov_switching_data* simulation(Format_error &error , int nb_sequence ,
				    const Switching_sequence &iseq ,
				    int inb_covariable, int iconstant,
				    double *iresidual_variance, double **iregression,
				    double ***icovar,
				    bool counting_flag) const;
  
  Markov_switching_data* simulation(Format_error &error ,
				    const Histogram &hlength ,
				    int inb_covariable, int iconstant,
				    double *iresidual_variance, double **iregression,
				    double ***icovar,
				    bool counting_flag  , bool divergence_flag ,
				    bool hidden ) const;
  
  Markov_switching_data* simulation(Format_error &error ,
				    int nb_sequence , int length ,
				    int inb_covariable, int iconstant,
				    double *iresidual_variance, double **iregression,
				    double ***icovar,
				    bool counting_flag , bool hidden ) const;
  
  Markov_switching_data* simulation(Format_error &error ,
				    int nb_sequence ,
				    const Switching_sequence &iseq ,
				    int inb_covariable, int iconstant,
				    double *iresidual_variance, double **iregression,
				    double ***icovar,
				    bool counting_flag , bool hidden ) const;


  // protected :
  // private :
  
};


/*----------------------------------------------------------------------
 *
 * Classe Markov_switching_data
 * Herite de Switching_sequence
 *
 *----------------------------------------------------------------------*/

class Markov_switching_data: public Switching_sequence {

  friend class Switching_sequence;
  friend class Markov_switching;

  friend std::ostream& operator<<(std::ostream &os, const Markov_switching_data &isw_markov_data)
  {return isw_markov_data.ascii_write(os,false); }

 private:

  Markov_switching *sw_markov;             // pointer on Markov_switching object
  Chain_data *chain_data;                  // initial and transition states
  double likelihood;                       // likelihood of sequences
  double hidden_likelihood;                // likelihood of all possible sequences
  double *posterior_probability;           // a posteriori probability of most probable state sequence

  void copy(const Markov_switching_data &isw_markov_data, bool model_flag = true);

 public:

  Markov_switching_data();
  Markov_switching_data(int inb_variable, const Histogram &ihlength, int inb_covariable, 
			int inb_random, int inb_val, int iconstant, bool init_flag);
  Markov_switching_data(int inb_variable, const Histogram &ihlength, int inb_covariable, 
			int inb_random, int inb_val, int T, int iconstant, bool init_flag);
  Markov_switching_data(int inb_sequence , int inb_covariable, int inb_random, int inb_val,
			int *ilength, int **iindex, double ***ireal_sequence, double ***icovar, double **ieffect, 
			double *iyear_effect, int iconstant, int *iidentifier, const Markov_switching &isw_markov);
  Markov_switching_data(const Switching_sequence &isw_seq, bool initial_run_flag);
  Markov_switching_data(const Switching_sequence &isw_seq, int variable);
  Markov_switching_data(const Markov_switching_data &isw_markov_data, bool model_flag = true)
  : Switching_sequence(isw_markov_data){ copy(isw_markov_data, model_flag);}
  ~Markov_switching_data();

  Markov_switching_data& operator=(const Markov_switching_data &isw_markov_data);

  std::ostream& ascii_write(std::ostream &os, bool exhaustive = false) const;
  bool ascii_write(Format_error &error, const char *path, bool exhaustive = false) const;
  bool spreadsheet_write(Format_error &error, const char *path) const;
  bool plot_write(Format_error &error, const char *prefix, const char *title = 0) const {}

  void build_transition_count(const Markov_switching &isw_markov, bool begin = true);

  double random_hetero_effect_variance(int index, double *random_variance, double **ieffect, int random = 0); 
  double random_year_effect_variance(double *random_variance, double *iyear_effect, int random = 1);
  double random_hetero_effect_mean(int index, double *random_variance, double **ieffect, int random = 0); 
  double random_year_effect_mean(double *random_variance, double *iyear_effect, int random = 1);

};

#endif

