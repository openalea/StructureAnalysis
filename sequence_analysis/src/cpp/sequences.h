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



#ifndef SEQUENCES_H
#define SEQUENCES_H



/****************************************************************
 *
 *  Constantes :
 */


const int DEFAULT_LENGTH = 20;         // longueur par defaut d'une sequence

const int SEQUENCE_NB_VARIABLE = NB_OUTPUT_PROCESS;  // nombre maximum de variables observees

const int PLOT_NB_SEQUENCE = 200;      // nombre maximum de sequences visualisees (sortie Gnuplot)
const int PLOT_TITLE_NB_SEQUENCE = 15;  // nombre maximum de sequences identifiees (sortie Gnuplot)

enum {
  IMPLICIT_TYPE ,                      // parametre d'index implicite
  TIME ,                               // parametre d'index temps
  TIME_INTERVAL ,                      // parametre d'index intervalle de temps
  POSITION ,                           // parametre d'index position
  POSITION_INTERVAL                    // parametre d'index intervalle inter-position
};

enum {
  OBSERVED_VALUE ,
  OBSERVED_STATE ,
  THEORETICAL_STATE ,
  OBSERVED_OUTPUT ,
  THEORETICAL_OUTPUT
};

enum {
  OBSERVATION ,
  FIRST_OCCURRENCE ,
  RECURRENCE_TIME ,
  SOJOURN_TIME ,
  INITIAL_RUN ,
  FINAL_RUN ,
  NB_RUN ,
  NB_OCCURRENCE ,
  LENGTH ,
  SEQUENCE_CUMUL ,
  SEQUENCE_MEAN
};

enum {
  DEFAULT ,
  REVERSE ,
  ADD_INITIAL_RUN ,
  REMOVE_INITIAL_RUN
};

enum {
  CTM_BIC ,                            // algorithme Context Tree Maximizing/BIC
  CTM_KT ,                             // algorithme Context Tree Maximizing/Krichevsky-Trofimov
  LOCAL_BIC ,                          // algorithme d'elagage recursif/BIC 
  CONTEXT                              // algorithme Context
};

enum {
  MAXIMUM_LIKELIHOOD ,
  LAPLACE ,
  ADAPTATIVE_LAPLACE ,
  UNIFORM_SUBSET ,
  UNIFORM_CARDINALITY
};

enum {
  MULTINOMIAL_CHANGE ,
  POISSON_CHANGE ,
  ORDINAL_GAUSSIAN_CHANGE ,
  GAUSSIAN_CHANGE ,
  MEAN_CHANGE ,
  VARIANCE_CHANGE ,
  MEAN_VARIANCE_CHANGE
};

const double MAX_NB_WORD = 1.e7;       // nombre maximum de mots

const double STATIONARY_PROBABILITY_THRESHOLD = 1.e-8;  // seuil pour le calcul des probabilites
                                                        // stationnaires d'un modele en equilibre
const int STATIONARY_PROBABILITY_LENGTH = 10000;  // longueur maximum pour le calcul des probabilites
                                                  // stationnaires d'un modele en equilibre
const double LEAVE_INCREMENT = 1.e-6;  // seuil pour stopper le calcul de la probabilite
                                       // de quitter un etat/observation sans possibilite d'y revenir

const double CTM_BIC_THRESHOLD = 6.;   // seuil pour elaguer les memoires
const double CTM_KT_THRESHOLD = 12.;   // seuil pour elaguer les memoires
const double LOCAL_BIC_THRESHOLD = 10.;  // seuil pour elaguer les memoires
const double CONTEXT_THRESHOLD = 5.;  // seuil pour elaguer les memoires

const double OCCUPANCY_THRESHOLD = 0.99999;  // seuil sur la fonction de repartition
                                             // pour borner une loi d'occupation d'un etat
const double OCCUPANCY_MEAN = 10.;     // temps moyen d'occupation d'un etat

const int MIN_NB_STATE_SEQUENCE = 1;   // nombre de sequences d'etats 1ere iteration de l'algorithme MCEM
const int MAX_NB_STATE_SEQUENCE = 10;  // nombre de sequences d'etats maximum pour l'algorithme MCEM
const double NB_STATE_SEQUENCE_PARAMETER = 1.;  // parametre nombre de sequences d'etats pour l'algorithme MCEM

const int POSTERIOR_PROBABILITY_NB_SEQUENCE = 300; // nombre maximum de sequences pour la sortie des probabilites
                                                   // a posteriori des sequences d'etats les plus probables

const int NB_STATE_SEQUENCE = 10;      // nombre de sequences d'etats calculees

const int COUNT_MAX_LENGTH = 10000;    // longueur maximum des sequences pour l'extraction des comptages
const int COUNTING_MAX_LENGTH = 500;   // longueur maximum des sequences pour le calcul des lois de comptage

const int NB_SEQUENCE = 100000;        // nombre maximum de sequences simulees
const int MAX_LENGTH = 1000000;        // longueur maximum des sequences simulees
const int CUMUL_LENGTH = 1000000;      // longueur maximum cumulee des sequences simulees

enum {
  SEQUENCE ,
  TREND ,
  SUBTRACTION_RESIDUAL ,
  DIVISION_RESIDUAL ,
  STANDARDIZED_RESIDUAL
};

enum {
  ADAPTATIVE ,
  FIXED
};

enum {
  DELETION ,
  BEGIN_END_DELETION ,
  INSERTION ,
  BEGIN_END_INSERTION ,
  MATCH ,
  SUBSTITUTION ,
  TRANSPOSITION
};

enum {
  DATA ,
  GAP ,
  BEGIN_END_GAP
};

const int NB_ALIGNMENT = 1000000;      // nombre maximum d'alignements
const int DISPLAY_NB_ALIGNMENT = 30;   // nombre maximum d'alignements
                                       // pour la sortie detaillee ecran
const int FILE_NB_ALIGNMENT = 300;     // nombre maximum d'alignements
                                       // pour la sortie detaillee fichier

const double INDEL_FACTOR_1 = 0.51;    // facteur pour deduire le cout d'elision/insertion -
                                       // alignement simple
const double INDEL_FACTOR_N = 0.51;    // facteur pour deduire le cout d'elision/insertion -
                                       // alignement multiple
const double TRANSPOSITION_FACTOR = 0.;  // facteur pour deduire le cout de transposition
const double INDEL_DISTANCE = 1.0;     // cout d'elision/insertion

enum {
  CHANGE_POINT ,
  SEGMENT
};

const double ROUNDOFF_ERROR = 1.e-10;  // erreur sur une somme de doubles
const int NB_SEGMENTATION = 10;        // nombre de segmentations calculees

enum {
  APPROXIMATED ,
  EXACT
};

const double FREQUENCY_RATIO = 0.1;    // rapport des frequences pour stopper le calcul
                                       // de la fonction de correlation
// const int CORRELATION_MIN_FREQUENCY = 20;  frequence minimum pour stopper le calcul
                                           // de la fonction de correlation

const int MAX_DIFFERENCING_ORDER = 3;  // ordre maximum de differenciation
const int POINTWISE_AVERAGE_NB_SEQUENCE = 250;  // nombre maximum de sequences ecrites
                                                // pour la sortie fichier
const int ABSORBING_RUN_LENGTH = 5;    // longueur par defaut de la serie finale absorbante
const int MAX_ABSORBING_RUN_LENGTH = 20;  // longueur maximum de la serie finale absorbante



/****************************************************************
 *
 *  Definition des classes :
 */


class Nonparametric_sequence_process;
class Variable_order_markov;
class Variable_order_markov_data;
class Variable_order_markov_iterator;
class Hidden_variable_order_markov;
class Semi_markov;
class Semi_markov_data;
class Semi_markov_iterator;
class Hidden_semi_markov;
class Nonhomogeneous_markov;
class Nonhomogeneous_markov_data;
class Sequences;
class Markovian_sequences;
class Sequence_characteristics;
class Vectors;

class Switching_sequence;  // ajout par Florence Chaubert


class Nonparametric_sequence_process : public Nonparametric_process {  // processus d'observation non-parametrique
                                                                       // pour des sequences

    friend class Variable_order_markov;
    friend class Variable_order_markov_iterator;
    friend class Variable_order_markov_data;
    friend class Hidden_variable_order_markov;
    friend class Semi_markov;
    friend class Semi_markov_iterator;
    friend class Semi_markov_data;
    friend class Hidden_semi_markov;
    friend class Nonhomogeneous_markov;
    friend class Nonhomogeneous_markov_data;
    friend class Markovian_sequences;

    friend Nonparametric_sequence_process* occupancy_parsing(Format_error &error , ifstream &in_file ,
                                                             int &line , const Chain &chain ,
                                                             double cumul_threshold);
    friend bool test_hidden(int nb_output_process , Nonparametric_sequence_process **process);

private :

    Distribution *length;   // loi des longueurs des sequences
    Curves *index_value;    // probabilites des valeurs en fonction de l'index
    double *no_occurrence;  // probabilite de ne pas observer une valeur
    Distribution **first_occurrence;  // lois du temps avant la premiere occurrence d'une valeur
    double *leave;          // probabilite de quitter une valeur
    Distribution **recurrence_time;  // lois du temps de retour dans une valeur
    double *absorption;     // probabilite d'etre absorbe par une valeur
    Parametric **sojourn_time;  // lois du temps de sejour dans un valeur
    Distribution **nb_run;  // lois du nombre de series par sequence des valeurs
    Distribution **nb_occurrence;  // lois du nombre d'occurrences par sequence des valeurs

    void create_characteristic(const Distribution &ilength , bool* homogeneity ,
                               bool counting_flag = true);
    void create_characteristic(const Distribution &ilength , bool sojourn_time_flag = true ,
                               bool counting_flag = true);
    void copy(const Nonparametric_sequence_process &process , bool characteristic_flag = true);
    void init_occupancy(const Nonparametric_sequence_process &process , int occupancy_nb_value);
    void remove();

    std::ostream& ascii_print(std::ostream &os , int process , Histogram **empirical_observation ,
                              const Sequence_characteristics *characteristics ,
                              bool exhaustive , bool file_flag , Forward **forward = 0) const;
    std::ostream& spreadsheet_print(std::ostream &os , int process ,
                                    Histogram **empirical_observation = 0 ,
                                    const Sequence_characteristics *characteristics = 0 ,
                                    Forward **forward = 0) const;
    bool plot_print(const char *prefix , const char *title , int process ,
                    Histogram **empirical_observation = 0 ,
                    const Sequence_characteristics *characteristics = 0 ,
                    const Histogram *hlength = 0 , Forward **forward = 0) const;

/*    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

public :

    Nonparametric_sequence_process(int inb_state = 0 , int inb_value = 0 ,
                                   int observation_flag = false);
    Nonparametric_sequence_process(int inb_state , Parametric **occupancy);
    Nonparametric_sequence_process(const Nonparametric_process &process);
    Nonparametric_sequence_process(const Nonparametric_sequence_process &process ,
                                   char manip = 'c' , int param = true);
    ~Nonparametric_sequence_process();
    Nonparametric_sequence_process& operator=(const Nonparametric_sequence_process &process);

    // acces membres de la classe

    Distribution* get_length() const { return length; }
    Curves* get_index_value() const { return index_value; }
    double get_no_occurrence(int value) const { return no_occurrence[value]; }
    Distribution** get_first_occurrence() const { return first_occurrence; }
    Distribution* get_first_occurrence(int value) const { return first_occurrence[value]; }
    double get_leave(int value) const { return leave[value]; }
    Distribution** get_recurrence_time() const { return recurrence_time; }
    Distribution* get_recurrence_time(int value) const { return recurrence_time[value]; }
    double get_absorption(int value) const { return absorption[value]; }
    Parametric** get_sojourn_time() const { return sojourn_time; }
    Parametric* get_sojourn_time(int value) const { return sojourn_time[value]; }
    Distribution** get_nb_run() const { return nb_run; }
    Distribution* get_nb_run(int value) const { return nb_run[value]; }
    Distribution** get_nb_occurrence() const { return nb_occurrence; }
    Distribution* get_nb_occurrence(int value) const { return nb_occurrence[value]; }
};


Nonparametric_sequence_process* occupancy_parsing(Format_error &error , ifstream &in_file ,
                                                  int &line , const Chain &chain ,
                                                  double cumul_threshold = CUMUL_THRESHOLD);
bool test_hidden(int nb_output_process , Nonparametric_sequence_process **process);



// class Correlation : public STAT_interface , public Curves {
class Correlation : public STAT_interface , protected Curves {  // fonctions de correlation

    friend class Variable_order_markov;
    friend class Sequences;

    friend std::ostream& operator<<(std::ostream &os , const Correlation &correlation)
    { return correlation.ascii_write(os); }

private :

    int type;               // type de coefficient (PEARSON/SPEARMAN/KENDALL)
    int *variable_type;     // types de variables (OBSERVED/THEORETICAL STATE/OUTPUT)
    int *variable1;         // 1ere variables
    int *variable2;         // 2eme variables
    double *white_noise;    // fonction theorique pour un bruit blanc filtre

    void copy(const Correlation &correl);
    void remove();
    bool plot_print(const char *path , double *confidence_limit) const;

public :

    Correlation();
    Correlation(int itype , int max_lag , int ivariable1 , int ivariable2);
    Correlation(int inb_curve , int ilength , bool frequency_flag , int itype);
    Correlation(const Correlation &correl)
    :Curves(correl) { copy(correl); }
    virtual ~Correlation();
    Correlation& operator=(const Correlation &correl);

    Correlation* merge(Format_error &error , int nb_correl , const Correlation **icorrel) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path , bool exhaustive) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,
                    const char *title = 0) const;

/*    RWDECLARE_COLLECTABLE(Correlation);

    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

    bool white_noise_correlation(Format_error &error , int nb_point , double *filter ,
                                 int residual = true);
    bool white_noise_correlation(Format_error &error , const Distribution &dist);
    bool white_noise_correlation(Format_error &error , int order);

    // acces membres de la classe

    int get_type() const { return type; }
    int get_variable_type(int index) const { return variable_type[index]; }
    int get_variable1(int index) const { return variable1[index]; }
    int get_variable2(int index) const { return variable2[index]; }
    double get_white_noise(int lag) const { return white_noise[lag]; }
};



class Distance_matrix;
class Renewal_data;
class Tops;
class Vector_distance;
class Vectors;

class Sequences : public STAT_interface {  // sequences

    friend class Vectors;

    friend Sequences* sequences_ascii_read(Format_error &error , const char *path ,
                                           bool old_format);
    friend std::ostream& operator<<(std::ostream &os , const Sequences &seq)
    { return seq.ascii_write(os); }

protected :

    int nb_sequence;        // nombre de sequences
    int *identifier;        // identificateurs des sequences
    int max_length;         // longueur maximum des sequences
    int cumul_length;       // longueur cumulee des sequences
    int *length;            // longueurs des sequences
    Histogram *hlength;     // histogramme des longueurs des sequences
    int index_parameter_type; // type du parametre d'index (TIME/POSITION)
    Histogram *hindex_parameter;   // histogramme des parametres d'index explicites
    Histogram *index_interval;  // intervalles entre parametres d'index explicites
    int **index_parameter;  // parametres d'index explicites
    int nb_variable;        // nombre de variables
    int *type;              // type de chaque variable (INT_VALUE/REAL_VALUE/STATE/NB_INTERNODE)
    double *min_value;      // valeurs minimums
    double *max_value;      // valeurs maximums
    Histogram **marginal;   // lois marginales empiriques
    int ***int_sequence;    // sequences, variables entieres
    double ***real_sequence;  // sequences, variables reelles

    void init(int inb_sequence , int *iidentifier , int *ilength , int iindex_parameter_type ,
              int inb_variable , int *itype , bool init_flag);
    void init(int inb_sequence , int *iidentifier , int *ilength , int inb_variable ,
              bool init_flag);
    void copy(const Sequences &seq , bool reverse_flag = false);
    void add_state_variable(const Sequences &seq);
    void remove_index_parameter(const Sequences &seq);
    void remove();

    bool increasing_index_parameter_checking(Format_error &error , bool strict ,
                                             const char *pattern_label) const;
    bool increasing_sequence_checking(Format_error &error , int variable , bool strict ,
                                      const char *pattern_label , const char *variable_label) const;

    void cluster(const Sequences &seq , int variable , int step , int mode);
    void transcode(const Sequences &seq , int ivariable , int min_symbol , int max_symbol ,
                   int *symbol , bool add_flag = false);
    void cluster(const Sequences &seq , int variable , int nb_class , double *limit);
    void select_variable(const Sequences &seq , int *variable);

    bool pointwise_average_ascii_print(Format_error &error , const char *path , int *size ,
                                       bool standard_deviation , int output) const;
    bool pointwise_average_spreadsheet_print(Format_error &error , const char *path , int *size ,
                                             bool standard_deviation , int output) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive , bool comment_flag) const;
    std::ostream& ascii_print(std::ostream &os , char format , bool comment_flag ,
                              double *posterior_probability = 0 ,
                              int line_nb_character = LINE_NB_CHARACTER) const;
    bool plot_print(const char *path , int ilength) const;

    void max_length_computation();
    void cumul_length_computation();
    void build_length_histogram();

    void index_parameter_computation();
    void build_index_parameter_histogram();
    void index_interval_computation();

    void min_value_computation(int variable);
    void max_value_computation(int variable);
    void build_marginal_histogram(int variable);

    std::ostream& alignment_ascii_print(std::ostream &os , int width , int ref_index , int test_index ,
                                        const Sequences &alignment , int alignment_index) const;
    std::ostream& alignment_spreadsheet_print(std::ostream &os , int ref_index , int test_index ,
                                              const Sequences &alignment , int alignment_index) const;

    double indel_distance_computation(const Vector_distance &vector_dist ,
                                      double **rank , double **max_symbol_distance) const;
    double indel_distance_computation(const Vector_distance &vector_dist ,
                                      int index , int position , double **rank ,
                                      double **max_symbol_distance) const;
    double substitution_distance_computation(const Vector_distance &vector_dist , int ref_index ,
                                             int test_index , int ref_position , int test_position ,
                                             double **rank , const Sequences *test_seq = 0) const;
    double substitution_distance_computation(int ref_index , int test_index , int ref_position,
                                             int test_position , double substitution_distance) const;

    std::ostream& multiple_alignment_ascii_print(std::ostream &os) const;
    bool multiple_alignment_ascii_print(Format_error &error , const char *path) const;

    Sequences* multiple_alignment(const Sequences &test_seq , const Vector_distance &vector_dist ,
                                  double **rank , double **max_symbol_distance , bool begin_free ,
                                  bool end_free , int indel_cost , double indel_factor) const;

    void correlation_computation(Correlation &correl , int variable1 , int variable2 ,
                                 int normalization , bool individual_mean = false) const;

    std::ostream& profile_ascii_print(std::ostream &os , int index , int nb_segment ,
                                      double **profiles , const char *label ,
                                      double **mean = 0 , double **change_point = 0) const;
    std::ostream& profile_spreadsheet_print(std::ostream &os , int index , int nb_segment ,
                                            double **profiles , const char *label ,
                                            double **mean = 0 , double **change_point = 0) const;
    std::ostream& profile_plot_print(std::ostream &os , int index , int nb_segment ,
                                     double **profiles , double **mean = 0 ,
                                     double **change_point = 0) const;

    int nb_parameter_computation(int index , int nb_segment , int *model_type) const;
    double one_segment_likelihood(int index , int *model_type , double **rank) const;
    Sequences* segmentation_output(int *nb_segment , int *model_type , std::ostream &os ,
                                   int output = SEQUENCE , int* ichange_point = 0);
    double segmentation(int *nb_segment , int *model_type , double **rank ,
                        double *isegmentation_likelihood = 0 , int *nb_parameter = 0 ,
                        double *segment_penalty = 0);
    double forward_backward(int index , int nb_segment , int *model_type , double **rank ,
                            std::ostream *os , int output , char format ,
                            double *ilikelihood = 0 , double *ichange_point_entropy = 0 ,
                            double *isegment_entropy = 0) const;
    double forward_backward_sampling(int index , int nb_segment , int *model_type ,
                                     double **rank , std::ostream &os , char format ,
                                     int nb_segmentation) const;
    double L_segmentation(int index , int nb_segment , int *model_type , double **irank ,
                          std::ostream &os , char format , int inb_segmentation ,
                          double likelihood) const;
    double forward_backward_dynamic_programming(int index , int nb_segment , int *model_type ,
                                                double **rank , std::ostream &os , int output ,
                                                char format , double likelihood = D_INF) const;

    std::ostream& profile_ascii_print(std::ostream &os , int index , int nb_state ,
                                      double **profiles , double *begin_conditional_entropy ,
                                      double *marginal_entropy , double *begin_partial_entropy ,
                                      double *end_conditional_entropy = 0 , double *end_partial_entropy = 0) const;
    std::ostream& profile_spreadsheet_print(std::ostream &os , int index , int nb_state ,
                                            double **profiles , double *begin_conditional_entropy ,
                                            double *marginal_entropy , double *begin_partial_entropy ,
                                            double *end_conditional_entropy = 0 , double *end_partial_entropy = 0) const;
    std::ostream& profile_plot_print(std::ostream &os , int index , int nb_state ,
                                     double **profiles , double *begin_conditional_entropy ,
                                     double *marginal_entropy , double *begin_partial_entropy ,
                                     double *end_conditional_entropy = 0 , double *end_partial_entropy = 0) const;

public :

    Sequences();
    Sequences(int inb_sequence , int inb_variable);
    Sequences(int inb_sequence , int *iidentifier , int *ilength , int iindex_parameter_type ,
              int inb_variable , int *itype , bool init_flag = false)
    { init(inb_sequence , iidentifier , ilength , iindex_parameter_type , inb_variable ,
           itype , init_flag); }
    Sequences(int inb_sequence , int *iidentifier , int *ilength , int inb_variable ,
              int *itype , bool init_flag = false)
    { init(inb_sequence , iidentifier , ilength , IMPLICIT_TYPE , inb_variable ,
           itype , init_flag); }
    Sequences(int inb_sequence , int *iidentifier , int *ilength , int iindex_parameter_type ,
              int inb_variable , int itype , int ***iint_sequence);
    Sequences(int inb_sequence , int *iidentifier , int *ilength , int inb_variable ,
              double ***ireal_sequence);
    Sequences(int inb_sequence , int *iidentifier , int *ilength , int inb_variable ,
              bool init_flag = false)
    { init(inb_sequence , iidentifier , ilength , inb_variable , init_flag); }
    Sequences(const Histogram &ihlength , int inb_variable , bool init_flag = false);
    Sequences(const Renewal_data &timev);
    Sequences(const Sequences &seq , int inb_sequence , int *index);
    Sequences(const Sequences &seq , bool *segment_mean);
    Sequences(const Sequences &seq , char transform = 'c' , int param = DEFAULT);
    virtual ~Sequences();
    Sequences& operator=(const Sequences &seq);

    Distribution_data* extract(Format_error &error , int variable) const;

    Vectors* build_vectors(bool index_variable) const;
    Vectors* extract_vectors(Format_error &error , int feature_type , int variable = I_DEFAULT ,
                             int value = I_DEFAULT) const;

    Markovian_sequences* markovian_sequences(Format_error &error) const;
    Tops* tops(Format_error &error) const;

    bool check(Format_error &error , const char *pattern_label);

    Time_events* extract_time_events(Format_error &error , int variable ,
                                     int begin_date , int end_date ,
                                     int previous_date = I_DEFAULT , int next_date = I_DEFAULT) const;
    Renewal_data* extract_renewal_data(Format_error &error , int variable ,
                                       int begin_index , int end_index) const;

    Sequences* merge(Format_error &error , int nb_sample , const Sequences **iseq) const;

    Sequences* shift(Format_error &error , int variable , int shift_param) const;
    Sequences* shift(Format_error &error , int variable , double shift_param) const;
    Sequences* cluster(Format_error &error , int variable , int step ,
                       int mode = FLOOR) const;
    Sequences* transcode(Format_error &error , int variable , int *symbol) const;
    Sequences* cluster(Format_error &error , int variable , int nb_class ,
                       int *ilimit) const;
    Sequences* cluster(Format_error &error , int variable , int nb_class ,
                       double *ilimit) const;
    Sequences* scaling(Format_error &error , int variable , int scaling_coeff) const;
    Sequences* round(Format_error &error , int variable = I_DEFAULT ,
                     int mode = ROUND) const;

    Sequences* index_parameter_select(Format_error &error , std::ostream &os ,
                                      int min_index_parameter ,
                                      int max_index_parameter , bool keep) const;
    Sequences* value_select(Format_error &error , std::ostream &os , int variable ,
                            int imin_value , int imax_value , bool keep = true) const;
    Sequences* value_select(Format_error &error , std::ostream &os , int variable ,
                            double imin_value , double imax_value , bool keep = true) const;
    Sequences* select_individual(Format_error &error , int inb_sequence , int *iidentifier ,
                                 bool keep = true) const;

    Sequences* remove_index_parameter(Format_error &error) const;
    Sequences* select_variable(Format_error &error , int inb_variable , int *ivariable ,
                               bool keep = true) const;
    Sequences* merge_variable(Format_error &error , int nb_sample , const Sequences **iseq ,
                              int ref_sample = I_DEFAULT) const;

    Sequences* reverse(Format_error &error) const;
    Sequences* length_select(Format_error &error , std::ostream &os , int min_length ,
                             int imax_length , bool keep = true) const;
    Sequences* remove_run(Format_error &error , int variable , int ivalue ,
                          char position , int max_run_length = I_DEFAULT) const;
    Sequences* index_parameter_extract(Format_error &error , int min_parameter_index ,
                                       int max_parameter_index = I_DEFAULT) const;
    Sequences* segmentation_extract(Format_error &error , int variable , int nb_value ,
                                    int *ivalue , bool keep = true) const;

    Sequences* cumulate(Format_error &error , int variable = I_DEFAULT) const;
    Sequences* difference(Format_error &error , int variable = I_DEFAULT ,
                          bool first_element = false) const;
    Sequences* moving_average(Format_error &error , int nb_point , double *filter ,
                              int variable = I_DEFAULT , bool begin_end = false ,
                              int output = TREND) const;
    Sequences* moving_average(Format_error &error , const Distribution &dist ,
                              int variable = I_DEFAULT , bool begin_end = false ,
                              int output = TREND) const;

    Sequences* pointwise_average(Format_error &error , bool standard_deviation = false ,
                                 int output = SEQUENCE , const char *path = 0 ,
                                 char format = 'a') const;

    Sequences* recurrence_time_sequences(Format_error &error , int variable , int value) const;
    Sequences* sojourn_time_sequences(Format_error &error , int variable) const;

    Sequences* transform_position(Format_error &error , int step) const;

    Sequences* cross(Format_error &error) const;

    std::ostream& line_write(std::ostream &os) const;

    virtual std::ostream& ascii_data_write(std::ostream &os , char format = 'c' ,
                                           bool exhaustive = false) const;
    virtual bool ascii_data_write(Format_error &error , const char *path ,
                                  char format = 'c' , bool exhaustive = false) const;
    bool plot_data_write(Format_error &error , const char *prefix ,
                         const char *title = 0) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,
                    const char *title = 0) const;

  /*    RWDECLARE_COLLECTABLE(Sequences);

    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

    int min_index_parameter_computation() const;
    int max_index_parameter_computation(bool last_position = false) const;

    void marginal_histogram_computation(int variable);
    double mean_computation(int variable) const;
    double variance_computation(int variable , double mean) const;
    double mean_absolute_deviation_computation(int variable , double mean) const;
    double mean_absolute_difference_computation(int variable) const;
    double skewness_computation(int variable , double mean , double variance) const;
    double kurtosis_computation(int variable , double mean , double variance) const;

    Histogram* value_index_interval_computation(Format_error &error , int variable , int value) const;

    Correlation* correlation_computation(Format_error &error , int variable1 , int variable2 ,
                                         int itype = PEARSON , int max_lag = I_DEFAULT ,
                                         int normalization = EXACT , bool individual_mean = false) const;
    Correlation* partial_autocorrelation_computation(Format_error &error , int variable ,
                                                     int itype = PEARSON , int max_lag = I_DEFAULT) const;

    Distance_matrix* alignment(Format_error &error , std::ostream *os , const Vector_distance &ivector_dist ,
                               int ref_identifier = I_DEFAULT , int test_identifier = I_DEFAULT ,
                               bool begin_free = false , bool end_free = false , int indel_cost = ADAPTATIVE ,
                               double indel_factor = INDEL_FACTOR_1 , bool transposition_flag = false ,
                               double transposition_factor = TRANSPOSITION_FACTOR ,
                               const char *result_path = 0 , char result_format = 'a' ,
                               const char *alignment_path = 0 , char alignment_format = 'a') const;
    Distance_matrix* alignment(Format_error &error , std::ostream *os , int ref_identifier = I_DEFAULT ,
                               int test_identifier = I_DEFAULT , bool begin_free = false , bool end_free = false ,
                               const char *result_path = 0 , char result_format = 'a' ,
                               const char *alignment_path = 0 , char alignment_format = 'a') const;

    Sequences* multiple_alignment(Format_error &error , std::ostream &os , const Vector_distance &ivector_dist ,
                                  bool begin_free = false , bool end_free = false , int indel_cost = ADAPTATIVE ,
                                  double indel_factor = INDEL_FACTOR_N , int algorithm = AGGLOMERATIVE ,
                                  const char *path = 0) const;

    Sequences* segmentation(Format_error &error , std::ostream &os , int iidentifier ,
                            int nb_segment , int *ichange_point , int *model_type ,
                            int output = SEQUENCE) const;
    Sequences* segmentation(Format_error &error , std::ostream &os , int *nb_segment ,
                            int *model_type , int iidentifier = I_DEFAULT ,
                            int output = SEQUENCE) const;
    Sequences* segmentation(Format_error &error , std::ostream &os , int iidentifier ,
                            int max_nb_segment , int *model_type) const;

    Sequences* hierarchical_segmentation(Format_error &error , std::ostream &os , int iidentifier ,
                                         int max_nb_segment , int *model_type) const;

    Sequences* segmentation(Format_error &error , int iidentifier , int nb_segment ,
                            const Vector_distance &ivector_dist , std::ostream &os ,
                            int output = SEGMENT) const;

    bool segment_profile_write(Format_error &error , std::ostream &os , int iidentifier ,
                               int nb_segment , int *model_type , int output = SEGMENT ,
                               char format = 'a' , int segmentation = FORWARD_DYNAMIC_PROGRAMMING ,
                               int nb_segmentation = NB_SEGMENTATION) const;
    bool segment_profile_write(Format_error &error , const char *path , int iidentifier ,
                               int nb_segment , int *model_type , int output = SEGMENT ,
                               char format = 'a' , int segmentation = FORWARD_DYNAMIC_PROGRAMMING ,
                               int nb_segmentation = NB_SEGMENTATION) const;
    bool segment_profile_plot_write(Format_error &error , const char *prefix ,
                                    int iidentifier , int nb_segment , int *model_type ,
                                    int output = SEGMENT , const char *title = 0) const;

    // acces membres de la classe

    int get_nb_sequence() const { return nb_sequence; }
    int get_identifier(int iseq) const { return identifier[iseq]; }
    int get_max_length() const { return max_length; }
    int get_cumul_length() const { return cumul_length; }
    int get_length(int index_seq) const { return length[index_seq]; }
    Histogram* get_hlength() const { return hlength; }
    int get_index_parameter_type() const { return index_parameter_type; }
    Histogram* get_hindex_parameter() const { return hindex_parameter; }
    Histogram* get_index_interval() const { return index_interval; }
    int get_index_parameter(int iseq , int index) const
    { return index_parameter[iseq][index]; }
    int get_nb_variable() const { return nb_variable; }
    int get_type(int variable) const { return type[variable]; }
    double get_min_value(int variable) const { return min_value[variable]; }
    double get_max_value(int variable) const { return max_value[variable]; }
    Histogram* get_marginal(int variable) const { return marginal[variable]; }
    int get_int_sequence(int iseq , int variable , int index) const
    { return int_sequence[iseq][variable][index]; }
    double get_real_sequence(int iseq , int variable , int index) const
    { return real_sequence[iseq][variable][index]; }
};


Sequences* sequences_ascii_read(Format_error &error , const char *path ,
                                bool old_format = false);



class Sequence_characteristics {  // caracteristiques des sequences pour une variable

    friend class Nonparametric_sequence_process;
    friend class Variable_order_markov;
    friend class Variable_order_markov_data;
    friend class Semi_markov;
    friend class Semi_markov_data;
    friend class Nonhomogeneous_markov;
    friend class Nonhomogeneous_markov_data;
    friend class Markovian_sequences;

    friend class Switching_sequence;  // ajout par Florence Chaubert

private :

    int nb_value;           // nombre de valeurs a partir de 0
    Curves *index_value;    // probabilites des valeurs en fonction de l'index
    Histogram **first_occurrence;  // histogrammes du temps avant la premiere
                                    // occurrence d'une valeur
    Histogram **recurrence_time;  // histogrammes du temps de retour dans une valeur
    Histogram **sojourn_time;  // histogrammes du temps de sejour dans une valeur
                                // (series completes)
    Histogram **initial_run;  // histogrammes du temps de sejour dans une valeur
                              // (series initiales censurees a gauche)
    Histogram **final_run;  // histogrammes du temps de sejour dans une valeur
                            // (series finales censurees a droite)
    Histogram **nb_run;     // histogrammes du nombre de series par sequence des valeurs
    Histogram **nb_occurrence;  // histogrammes du nombre d'occurrences
                                // par sequence des valeurs

    void copy(const Sequence_characteristics &characteristics);
    void reverse(const Sequence_characteristics &characteristics);
    void remove();

    void create_sojourn_time_histogram(int max_length , int initial_run_flag = false);

    std::ostream& ascii_print(std::ostream &os , int type , const Histogram &hlength ,
                              bool exhaustive , bool comment_flag) const;
    std::ostream& spreadsheet_print(std::ostream &os , int type , const Histogram &hlength) const;
    bool plot_print(const char *prefix , const char *title , int variable ,
                    int nb_variable , int type , const Histogram &hlength) const;

/*    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

public :

    Sequence_characteristics(int inb_value = I_DEFAULT);
    Sequence_characteristics(const Sequence_characteristics &characteristics ,
                             bool initial_run_flag);
    Sequence_characteristics(const Sequence_characteristics &characteristics ,
                             char transform = 'c');
    ~Sequence_characteristics();
    Sequence_characteristics& operator=(const Sequence_characteristics &characteristics);

    // acces membres de la classe

    int get_nb_value() const { return nb_value; }
    Curves* get_index_value() const { return index_value; }
    Histogram* get_first_occurrence(int value) const { return first_occurrence[value]; }
    Histogram* get_recurrence_time(int value) const { return recurrence_time[value]; }
    Histogram* get_sojourn_time(int value) const { return sojourn_time[value]; }
    Histogram** get_initial_run() const { return initial_run; }
    Histogram* get_initial_run(int value) const { return initial_run[value]; }
    Histogram* get_final_run(int value) const { return final_run[value]; }
    Histogram* get_nb_run(int value) const { return nb_run[value]; }
    Histogram* get_nb_occurrence(int value) const { return nb_occurrence[value]; }
};



class Function;

class Self_transition : public Curves {  // probabilites de rester dans un etat en fonction de l'index

public :

    Self_transition(int ilength)
    :Curves(1 , ilength , true , false) {}

    Function* monomolecular_regression() const;
    Function* logistic_regression() const;
};



class Variable_order_chain_data;

class Markovian_sequences : public Sequences {  // trajectoires correspondant a
                                                // un processus markovien
    friend class Variable_order_markov;
    friend class Hidden_variable_order_markov;
    friend class Semi_markov;
    friend class Hidden_semi_markov;
    friend class Nonhomogeneous_markov;

    friend std::ostream& operator<<(std::ostream &os , const Markovian_sequences &seq)
    { return seq.ascii_write(os); }

protected :

    Self_transition **self_transition;  // probabilites de rester dans un etat
                                        // en fonction de l'index
    Histogram ***observation;  // histogrammes correspondant aux lois d'observation
    Sequence_characteristics **characteristics;  // caracteristiques pour une variable donnee

    void init();
    void copy(const Markovian_sequences &seq , int param = DEFAULT);
    void add_state_variable(const Markovian_sequences &seq , int param);
    void remove();

    Markovian_sequences* transcode(Format_error &error ,
                                   const Nonparametric_sequence_process *process) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive , bool comment_flag) const;
    bool plot_print(const char *prefix , const char *title , int variable ,
                    int nb_variable) const;

    void state_variable_init(int itype = STATE);

    void transition_count_computation(const Variable_order_chain_data &chain_data ,
                                      const Variable_order_markov &markov ,
                                      bool begin = true , bool non_terminal = false) const;
    void transition_count_computation(const Chain_data &chain_data ,
                                      const Semi_markov *smarkov = 0) const;

    void self_transition_computation(int state);
    void observation_histogram_computation(int variable);
    bool test_hidden(int variable) const;
    void build_index_value(int variable);
    void build_first_occurrence_histogram(int variable);
    void build_recurrence_time_histogram(int variable);
    void build_sojourn_time_histogram(int variable , int initial_run_flag = false);
    void build_nb_run_histogram(int variable);
    void build_nb_occurrence_histogram(int variable);

    void censored_sojourn_time_histogram_computation(Histogram **initial_run ,
                                                     Histogram **final_run ,
                                                     Histogram **single_run) const;

    std::ostream& likelihood_write(std::ostream &os , int nb_model , double **likelihood ,
                                   const char *label , bool exhaustive = false ,
                                   char algorithm = 'd') const;
    bool likelihood_write(Format_error &error , const char *path , int nb_model ,
                          double **likelihood , const char *label , char algorithm = 'd') const;

public :

    Markovian_sequences();
    Markovian_sequences(int inb_sequence , int *iidentifier , int *ilength , int iindex_parameter_type ,
                        int inb_variable , int *itype , bool init_flag = false)
    :Sequences(inb_sequence , iidentifier , ilength , iindex_parameter_type , inb_variable ,
               itype , init_flag) { init(); }
    Markovian_sequences(int inb_sequence , int *iidentifier , int *ilength ,
                        int inb_variable , bool init_flag = false)
    :Sequences(inb_sequence , iidentifier , ilength , inb_variable , init_flag) { init(); }
    Markovian_sequences(const Histogram &ihlength , int inb_variable , bool init_flag = false)
    :Sequences(ihlength , inb_variable , init_flag) { init(); }
    Markovian_sequences(const Sequences &seq);
    Markovian_sequences(const Markovian_sequences &seq , char transform = 'c' ,
                        int param = DEFAULT);
    ~Markovian_sequences();
    Markovian_sequences& operator=(const Markovian_sequences &seq);

    Distribution_data* extract(Format_error &error , int type ,
                               int variable , int value) const;

    Markovian_sequences* merge(Format_error &error , int nb_sample ,
                               const Markovian_sequences **iseq) const;

    Markovian_sequences* cluster(Format_error &error , int variable , int step ,
                                 int mode = FLOOR) const;
    Markovian_sequences* transcode(Format_error &error , int ivariable , int *symbol ,
                                   bool add_flag = false) const;
    Markovian_sequences* consecutive_values(Format_error &error , std::ostream &os ,
                                            int ivariable , bool add_flag = false) const;
    Markovian_sequences* cluster(Format_error &error , int ivariable , int nb_class ,
                                 int *ilimit , bool add_flag = false) const;
    Markovian_sequences* cluster(Format_error &error , int variable , int nb_class ,
                                 double *ilimit) const;

    Markovian_sequences* remove_index_parameter(Format_error &error) const;
    Markovian_sequences* select_variable(Format_error &error , int inb_variable ,
                                         int *ivariable , bool keep = true) const;
    Markovian_sequences* merge_variable(Format_error &error , int nb_sample ,
                                        const Markovian_sequences **iseq , int ref_sample = I_DEFAULT) const;
    Markovian_sequences* remove_variable_1() const;

    Markovian_sequences* initial_run_computation(Format_error &error) const;
    Markovian_sequences* add_absorbing_run(Format_error &error ,
                                           int sequence_length = I_DEFAULT ,
                                           int run_length = I_DEFAULT) const;

    Markovian_sequences* split(Format_error &error , int step) const;

    std::ostream& ascii_data_write(std::ostream &os , char format = 'c' ,
                                   bool exhaustive = false) const;
    bool ascii_data_write(Format_error &error , const char *path ,
                          char format = 'c' , bool exhaustive = false) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,
                    const char *title = 0) const;

    bool transition_count(Format_error &error , std::ostream &os , int max_order ,
                          bool begin = false , int estimator = MAXIMUM_LIKELIHOOD ,
                          const char *path = 0) const;
    bool word_count(Format_error &error , std::ostream &os , int variable , int word_length ,
                    int begin_state = I_DEFAULT , int end_state = I_DEFAULT ,
                    int min_frequency = 1) const;
    bool mtg_write(Format_error &error , const char *path , int *itype) const;

/*    RWDECLARE_COLLECTABLE(Markovian_sequences);

    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

    double iid_information_computation() const;

    void self_transition_computation();
    void self_transition_computation(bool *homogeneity);
    void sojourn_time_histogram_computation(int variable);
    void create_observation_histogram(int nb_state);
    void observation_histogram_computation();
    void build_observation_histogram();
    void build_characteristic(int variable = I_DEFAULT , bool sojourn_time_flag = true ,
                              bool initial_run_flag = false);

    Nonhomogeneous_markov* nonhomogeneous_markov_estimation(Format_error &error , int *ident ,
                                                            bool counting_flag = true) const;

    Variable_order_markov* variable_order_markov_estimation(Format_error &error , std::ostream &os ,
                                                            char model_type , int min_order = 0 ,
                                                            int max_order = ORDER ,
                                                            int algorithm = LOCAL_BIC ,
                                                            double threshold = LOCAL_BIC_THRESHOLD ,
                                                            int estimator = LAPLACE ,
                                                            bool global_initial_transition = true ,
                                                            bool global_sample = true ,
                                                            bool counting_flag = true) const;
    Variable_order_markov* variable_order_markov_estimation(Format_error &error ,
                                                            const Variable_order_markov &imarkov ,
                                                            bool global_initial_transition = true ,
                                                            bool counting_flag = true) const;
    Variable_order_markov* variable_order_markov_estimation(Format_error &error ,
                                                            char model_type , int order = 1 ,
                                                            bool global_initial_transition = true ,
                                                            bool counting_flag = true) const;

    Variable_order_markov* lumpability_estimation(Format_error &error , std::ostream &os , int *symbol ,
                                                  int penalty_type = BIC , int order = 1 ,
                                                  bool counting_flag = true) const;

    Semi_markov* semi_markov_estimation(Format_error &error , std::ostream &os , char model_type ,
                                        int estimator = COMPLETE_LIKELIHOOD , bool counting_flag = true ,
                                        int nb_iter = I_DEFAULT , int mean_computation = COMPUTED) const;

    Hidden_variable_order_markov* hidden_variable_order_markov_estimation(Format_error &error , std::ostream &os ,
                                                                          const Hidden_variable_order_markov &ihmarkov ,
                                                                          bool global_initial_transition = true ,
                                                                          bool counting_flag = true ,
                                                                          bool state_sequence = true ,
                                                                          int nb_iter = I_DEFAULT) const;
    Hidden_variable_order_markov* hidden_variable_order_markov_stochastic_estimation(Format_error &error , std::ostream &os ,
                                                                                     const Hidden_variable_order_markov &ihmarkov ,
                                                                                     bool global_initial_transition = true ,
                                                                                     int min_nb_state_sequence = MIN_NB_STATE_SEQUENCE ,
                                                                                     int max_nb_state_sequence = MAX_NB_STATE_SEQUENCE ,
                                                                                     double parameter = NB_STATE_SEQUENCE_PARAMETER ,
                                                                                     bool counting_flag = true ,
                                                                                     bool state_sequence = true ,
                                                                                     int nb_iter = I_DEFAULT) const;

    Hidden_semi_markov* hidden_semi_markov_estimation(Format_error &error , std::ostream &os ,
                                                      const Hidden_semi_markov &ihsmarkov ,
                                                      int estimator = COMPLETE_LIKELIHOOD ,
                                                      bool counting_flag = true ,
                                                      bool state_sequence = true ,
                                                      int nb_iter = I_DEFAULT ,
                                                      int mean_computation = COMPUTED) const;
    Hidden_semi_markov* hidden_semi_markov_estimation(Format_error &error , std::ostream &os ,
                                                      char model_type , int nb_state , bool left_right ,
                                                      int estimator = COMPLETE_LIKELIHOOD ,
                                                      bool counting_flag = true ,
                                                      bool state_sequence = true ,
                                                      double occupancy_mean = D_DEFAULT ,
                                                      int nb_iter = I_DEFAULT ,
                                                      int mean_computation = COMPUTED) const;
    Hidden_semi_markov* hidden_semi_markov_stochastic_estimation(Format_error &error , std::ostream &os ,
                                                                 const Hidden_semi_markov &ihsmarkov ,
                                                                 int min_nb_state_sequence = MIN_NB_STATE_SEQUENCE ,
                                                                 int max_nb_state_sequence = MAX_NB_STATE_SEQUENCE ,
                                                                 double parameter = NB_STATE_SEQUENCE_PARAMETER ,
                                                                 int estimator = COMPLETE_LIKELIHOOD ,
                                                                 bool counting_flag = true ,
                                                                 bool state_sequence = true ,
                                                                 int nb_iter = I_DEFAULT) const;
    Hidden_semi_markov* hidden_semi_markov_stochastic_estimation(Format_error &error , std::ostream &os ,
                                                                 char model_type , int nb_state , bool left_right ,
                                                                 int min_nb_state_sequence = MIN_NB_STATE_SEQUENCE ,
                                                                 int max_nb_state_sequence = MAX_NB_STATE_SEQUENCE ,
                                                                 double parameter = NB_STATE_SEQUENCE_PARAMETER ,
                                                                 int estimator = COMPLETE_LIKELIHOOD ,
                                                                 bool counting_flag = true ,
                                                                 bool state_sequence = true ,
                                                                 double occupancy_mean = D_DEFAULT ,
                                                                 int nb_iter = I_DEFAULT) const;

    bool lumpability_test(Format_error &error , std::ostream &os , int *symbol , int order = 1) const;

    bool comparison(Format_error &error , std::ostream &os , int nb_model ,
                    const Variable_order_markov **imarkov , const char *path = 0) const;

    bool comparison(Format_error &error , std::ostream &os , int nb_model ,
                    const Semi_markov **ismarkov , const char *path = 0) const;

    bool comparison(Format_error &error , std::ostream &os , int nb_model ,
                    const Hidden_variable_order_markov **ihmarkov ,
                    int algorithm = FORWARD , const char *path = 0) const;

    bool comparison(Format_error &error , std::ostream &os , int nb_model ,
                    const Hidden_semi_markov **ihsmarkov ,
                    int algorithm = FORWARD , const char *path = 0) const;

    // acces membres de la classe

    Curves* get_self_transition(int state) const { return self_transition[state]; }
    Histogram*** get_observation() const { return observation; }
    Histogram** get_observation(int variable) const { return observation[variable]; }
    Histogram* get_observation(int variable , int state) const
    { return observation[variable][state]; }
    Sequence_characteristics* get_characteristics(int variable) const
    { return characteristics[variable]; }
};



#endif
