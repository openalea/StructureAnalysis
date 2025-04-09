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


#ifndef SWITCHING_SEQUENCE_H
#define SWITCHING_SEQUENCE_H

#include <fstream>

enum{
  WITH_CONSTANT,
  WITHOUT_CONSTANT
};

enum{
  SINGLE_RANDOM,
  MULTIPLE_RANDOM
};

enum{
  NORANDOM,
  HETERORANDOM,
  YEARRANDOM
};


class Markov_switching;


class Switching_sequence: public Sequences { 

  friend std::ostream& operator<<(std::ostream &os, const Switching_sequence &sw_seq)
  {  return sw_seq.ascii_write(os); }

  friend class Sequence_characteristics;
  friend class Markov_switching;

  public:
  int constant;                                  // with or without intercept (WITH_CONSTANT/WITHOUT_CONSTANT)
  int nb_covariable;                             // number of covariates
  int nb_random;                                 // number of random effects
  int nb_val;                                    // number of states (use in case of a random effect by state: MULTIPLE_RANDOM)
  Continuous_histo ***observation;               // histograms of observation process
  double ***covar;                               // covariates 
  double **effect;                               // random effects (heterogeneity case)
  double *year_effect;                           // random effects (temporal case)
  int **index;                                   // index 
  Self_transition **self_transition;             // probability to stay in a state according to index
  Sequence_characteristics **characteristics;    // state process characteristics ( state in wich individual is at time t)

  Switching_sequence();
  Switching_sequence(int inb_variable, int *itype, int inb_sequence, int *iidentifier, int *ilength,  
		     int inb_covariable, int inb_random, int inb_val, int iconstant = WITHOUT_CONSTANT, bool init_flag = true);
  Switching_sequence(int inb_variable, int inb_sequence, int *iidentifier, int *ilength, int inb_covariable,
		     int inb_random, int inb_val, int iconstant = WITHOUT_CONSTANT,bool init_flag = true);
  Switching_sequence(int inb_variable, const Histogram &ihlength, int inb_covariable, int inb_random, 
		     int inb_val, int iconstant = WITHOUT_CONSTANT, bool init_flag = true);
  Switching_sequence(int inb_variable, const Histogram &ihlength, int inb_covariable, int inb_random, 
		     int inb_val, int T, int iconstant = WITHOUT_CONSTANT, bool init_flag = true);
  Switching_sequence(int inb_variable, int inb_sequence, int inb_covariable, int inb_random, int inb_val, int *ilength, 
		     int **iindex, double ***ireal_sequence, double ***icovar, double **ieffect, 
		     double *iyear_effect, int *iidentifier, int iconstant = WITHOUT_CONSTANT);
  Switching_sequence(const Sequences &seq, int inb_covariable, int inb_random, int inb_val, int **iindex,
		     double ***icovar, double **ieffect, double *iyear_effect, int iconstant = WITHOUT_CONSTANT);
  Switching_sequence(const Switching_sequence &sw_seq, 
		     char transform ='c' , int param1 = DEFAULT , int param2 = DEFAULT);

  ~Switching_sequence();

  void remove();

  void init();

  int nb_year(); 

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


  // Estimation LM sans effet aléatoire
  Markov_switching* markov_switching_estimation(Format_error &error , ostream &os ,
						const Markov_switching &isw_markov , 
						bool VarCommune = false, int estimator = COMPLETE_LIKELIHOOD ,
						bool output_R = false,bool counting_flag = true , bool state_sequence = true,
						int nb_iter = I_DEFAULT, int mean_computation = COMPUTED) const;

  Markov_switching* markov_switching_estimation(Format_error &error , ostream &os , 
						int nb_state , bool left_right ,  
						double *residual_variance, double **regression,
						int nb_covariable, int constant, bool markov = true,
						bool VarCommune = false, bool output_R = false,
						int order = 1, bool counting_flag = true, bool state_sequence = true,
						double occupancy_mean = D_DEFAULT,
						double self_transition = D_DEFAULT, int nb_iter = I_DEFAULT,
						int mean_computation = COMPUTED) const;

  Markov_switching* markov_switching_stochastic_estimation(Format_error &error , ostream &os ,
							   const Markov_switching &ihsw_markov ,
							   int min_nb_state_sequence = MIN_NB_STATE_SEQUENCE ,
							   int max_nb_state_sequence = MAX_NB_STATE_SEQUENCE ,
							   bool VarCommune = false, bool output_R = false,
							   double parameter = NB_STATE_SEQUENCE_PARAMETER,
							   int estimator = COMPLETE_LIKELIHOOD,
							   bool counting_flag = true, bool state_sequence  = true,
							   int nb_iter = I_DEFAULT) const;
  
  Markov_switching* markov_switching_stochastic_estimation(Format_error &error , ostream &os ,
							   char type , int nb_state , bool left_right ,
							   double *residual_variance, double **regression, 
							   int nb_covariable, int constant, int order,
							   int min_nb_state_sequence = MIN_NB_STATE_SEQUENCE ,
							   int max_nb_state_sequence = MAX_NB_STATE_SEQUENCE ,
							   bool markov = true, bool VarCommune = false, bool output_R = false,
							   double parameter = NB_STATE_SEQUENCE_PARAMETER, 
							   bool counting_flag = true, bool state_sequence  = true, 
							   double occupancy_mean = D_DEFAULT,
							   double self_transition = D_DEFAULT ,
							   int nb_iter = I_DEFAULT) const;

  double* regression_parameter_no_effect(const Markov_switching *hsw_markov, int nb_sequence, 
					 int *length, double ***real_sequence, int nb_covariable, 
					 double ***covar, double ***backward, int k) const;

  double residual_variance_no_effect(const Markov_switching *hsw_markov, int nb_sequence, 
				     int *length, double ***real_sequence, int nb_covariable, 
				     double ***covar, double ***backward, double *regression, int k) const;

  double residual_variance_no_effect_VarCommune(const Markov_switching *hsw_markov, int nb_sequence, 
						int *length, double ***real_sequence, int nb_covariable, 
						double ***covar, double ***backward, double **regression) const;

  double* regression_parameter_no_effect_stochastic(const Markov_switching *hsmarkov, int nb_sequence, 
						    int *length, double ***real_sequence, int nb_covariable, 
						    double ***covar, double ****indic, int k, int nb_state_sequence) const;

  double residual_variance_no_effect_stochastic(const Markov_switching *hsmarkov, int nb_sequence, 
						int *length, double ***real_sequence, int nb_covariable, 
						double ***covar, double ****indic, double *regression, int k, 
						int nb_state_sequence) const;

  double residual_variance_no_effect_VarCommune_stochastic(const Markov_switching *hsmarkov, int nb_sequence, 
							   int *length, double ***real_sequence, int nb_covariable, 
							   double ***covar, double ****indic, double **regression,
							   int nb_state_sequence) const;

  double* standard_error_regression_parameter(const Markov_switching *hsw_markov, int nb_sequence, 
					      int *length, double ***real_sequence, int nb_covariable, 
					      double ***covar, double ***backward, int k) const;

  double standard_error_residual_variance(const Markov_switching *hsw_markov, int nb_sequence, 
					  int *length, double ***real_sequence, int nb_covariable, 
					  double ***covar, double ***backward, int k) const;

  double standard_error_residual_variance_VarCommune(const Markov_switching *hsw_markov, int nb_sequence, 
						     int *length, double ***real_sequence, int nb_covariable, 
						     double ***covar, double ***backward) const;

  // Estimation LMM avec heterogeneite inter-individuelle
  Markov_switching* markov_switching_stochastic_estimation_hetero(Format_error &error , ostream &os ,
								  const Markov_switching &ihsw_markov ,
								  int type_random, int nb_iter_conv = I_DEFAULT,
								  int min_nb_state_sequence = MIN_NB_STATE_SEQUENCE ,
								  int max_nb_state_sequence = MAX_NB_STATE_SEQUENCE ,
								  bool VarCommune = false, bool output_R = false,
								  double parameter = NB_STATE_SEQUENCE_PARAMETER, 
								  int estimator = COMPLETE_LIKELIHOOD,
								  bool counting_flag = true, bool state_sequence  = true, 
								  int nb_iter = I_DEFAULT) const;

  Markov_switching* markov_switching_stochastic_estimation_hetero(Format_error &error , ostream &os ,
								  char type , int nb_state , bool left_right ,
								  double *residual_variance, double **random_variance,
								  double **regression, int nb_covariable,
								  int nb_random, int constant, int order,
								  int type_random, int nb_iter_conv = I_DEFAULT,
								  int min_nb_state_sequence = MIN_NB_STATE_SEQUENCE ,
								  int max_nb_state_sequence = MAX_NB_STATE_SEQUENCE ,
								  bool markov = true, bool VarCommune = false,  bool output_R = false,
								  double parameter = NB_STATE_SEQUENCE_PARAMETER, 
								  bool counting_flag = true, bool state_sequence  = true,
								  double occupancy_mean = D_DEFAULT,
								  double self_transition = D_DEFAULT ,
								  int nb_iter = I_DEFAULT) const;

  double** conditional_expectation_heterogeneity_state_wise(const Markov_switching *hsw_markov, int *length, 
							     int nb_covariable, double ***covar, double **effect,
							     double ***real_sequence, double ****indic, 
							     int nb_state_sequence, int indiv) const;

  double** conditional_expectation_heterogeneity_wise(const Markov_switching *hsw_markov, int *length, 
						      int nb_covariable, double ***covar, double **effect,
						      double ***real_sequence, double ****indic, 
						      int nb_state_sequence, int indiv) const;

  double* regression_parameter_heterogeneity_stochastic(const Markov_switching *hsw_markov, int nb_sequence, 
							int *length, double ***real_sequence, int nb_covariable, 
							double ***covar, double ***random_predict, double ****indic,  
							int k, int nb_state_sequence) const;
  
  double residual_variance_heterogeneity_stochastic(const Markov_switching *hsw_markov, int nb_sequence, 
						    int *length, double ***real_sequence, int nb_covariable, 
						    double ***covar, double ***random_predict, double ****indic, 
						    double *regression, int k, int nb_state_sequence) const;

  double* random_variance_heterogeneity_stochastic(const Markov_switching *hsw_markov, int nb_sequence, 
						   int *length, double ***real_sequence, int nb_covariable, 
						   double ***covar, double ***random_predict, double ****indic,
						   double *regression, int k,  int nb_state_sequence) const;

  double* standard_error_regression_parameter_heterogeneity(const Markov_switching *hsw_markov, int nb_sequence, 
							    int *length, double ***real_sequence, int nb_covariable, 
							    double ***covar, int nb_random, double ***random_predict, 
							    double ****indic, int k,  int nb_state_sequence) const;

  double standard_error_residual_variance_heterogeneity(const Markov_switching *hsw_markov, int nb_sequence, 
							int *length, double ***real_sequence, int nb_covariable, 
							double ***covar, int nb_random, double ***random_predict, 
							double ****indic, int k,  int nb_state_sequence) const;

  double* standard_error_random_variance_heterogeneity(const Markov_switching *hsw_markov, int nb_sequence, 
						       int *length, double ***real_sequence, int nb_covariable, 
						       double ***covar, int nb_random, double ***random_predict, 
						       double ****indic, int k,  int nb_state_sequence) const;

  // Estimation LMM avec effet environnement commun
  Markov_switching* markov_switching_stochastic_estimation_MH_year(Format_error &error , std::ostream &os ,
								   const Markov_switching &ihsw_markov , 
								   int nb_iter_conv = I_DEFAULT,
								   int min_nb_random = MIN_NB_STATE_SEQUENCE ,
								   int max_nb_random = MAX_NB_STATE_SEQUENCE ,
								   double VarMH = 0.01, int estimator = COMPLETE_LIKELIHOOD ,
								   bool VarCommune = false, bool output_R = false,
								   double parameter = NB_STATE_SEQUENCE_PARAMETER ,
								   bool counting_flag = true ,
								   bool state_sequence = true  ,
								   int nb_iter = I_DEFAULT,
								   int mean_computation = COMPUTED) const;

  Markov_switching* markov_switching_stochastic_estimation_MH_year(Format_error &error , ostream &os ,
								   char type , int nb_state , bool left_right ,
								   double *residual_variance, double **random_variance,
								   double **regression, int nb_covariable,
								   int nb_random, int constant, int order,
								   int nb_iter_conv = I_DEFAULT,
								   int min_nb_random = MIN_NB_STATE_SEQUENCE ,
								   int max_nb_random = MAX_NB_STATE_SEQUENCE ,
								   double VarMH = 0.01, int estimator = COMPLETE_LIKELIHOOD, 
								   bool markov = true, bool VarCommune = false, bool output_R = false,
								   double parameter = NB_STATE_SEQUENCE_PARAMETER, 
								   bool counting_flag = true, bool state_sequence  = true, 
								   double occupancy_mean = D_DEFAULT,
								   double self_transition = D_DEFAULT ,
								   int nb_iter = I_DEFAULT,
								   int mean_computation = COMPUTED) const;  

  double* regression_parameter_temporal_stochastic(const Markov_switching *hsw_markov, int nb_sequence, int *length, 
						   double ***real_sequence, int nb_covariable, double ***covar, double *year_effect, 
						   double ***backward, int **index, int nb_random, int k) const;

  double residual_variance_temporal_stochastic(const Markov_switching *hsw_markov, int nb_sequence, int *length, double ***real_sequence, 
					       int nb_covariable, double ***covar, double *year_effect, double ***backward, int **index, 
					       double *regression, int nb_random, int k) const;
  
  double residual_variance_temporal_VarCommune_stochastic(const Markov_switching *hsw_markov, int nb_sequence, int *length, 
							  double ***real_sequence, int nb_covariable, double ***covar, double *year_effect,
							  double ***backward, int **index, double **regression, int nb_random) const;

  double* random_variance_temporal_stochastic(const Markov_switching *hsw_markov, int nb_sequence, int *length, double ***real_sequence, 
					      int nb_covariable, double ***covar, double *year_effect, double ***backward, int **index, 
					      double *regression, int nb_random, int k) const;

  double* random_variance_temporal_VarCommune_stochastic(const Markov_switching *hsw_markov, int nb_sequence, int *length, 
							 double ***real_sequence, int nb_covariable, double ***covar, double *year_effect, 
							 double ***backward, int **index, double **regression, int nb_random) const;

  //private:

  //protected:

  // modifications des attributs
  void Set_effect(double **ieffect);
  void Set_year_effect(double *iyear_effect);
};

#endif
