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



#ifndef VECTORS_H
#define VECTORS_H



/****************************************************************
 *
 *  Constantes :
 */


const int VECTOR_NB_VARIABLE = 50;     // nombre maximum de variables
const int DISTANCE_NB_VECTOR = 2000;   // nombre maximum de vecteurs pour le calcul
                                       // d'une matrice des distances
const int CONTINGENCY_NB_VALUE = 100;  // nombre maximum de valeurs pour le calcul
                                       // d'un tableau de contingence
const int DISPLAY_CONTINGENCY_NB_VALUE = 20;  // nombre maximum de valeurs pour l'affichage
                                              // d'un tableau de contingence
const int VARIANCE_ANALYSIS_NB_VALUE = 100;  // nombre maximum de niveaux pour l'analyse de variance
const int DISPLAY_CONDITIONAL_NB_VALUE = 100;  // nombre maximum de valeurs pour l'affichage
                                               // des lois conditionnelles
const int PLOT_NB_VALUE = 30;          // seuil pour l'ecriture des frequences (sortie Gnuplot)
const double PLOT_RANGE_RATIO = 4.;    // seuil pour l'affichage a partir de 0

const int NB_SYMBOL = 50;              // nombre maximum de symboles

enum {
  ABSOLUTE_VALUE ,
  QUADRATIC
};



/****************************************************************
 *
 *  Definition des classes :
 */


class Distance_matrix;
class Regression;
class Sequences;
class TreeMatch;
class Vector_distance;
class Vectors;

class Vectors : public STAT_interface {  // vecteurs

    friend class Regression;
    friend class Sequences;

    friend Vectors* vectors_ascii_read(Format_error &error , const char *path);
    friend std::ostream& operator<<(std::ostream &os , const Vectors &vec)
    { return vec.ascii_write(os); }

protected :

    int nb_vector;          // nombre de vecteurs
    int *identifier;        // identificateurs des vecteurs
    int nb_variable;        // nombre de variables
    int *type;              // type de chaque variable (INT_VALUE/REAL_VALUE)
    double *min_value;      // valeurs minimums
    double *max_value;      // valeurs maximums
    Histogram **marginal;   // lois marginales empiriques
    double *mean;           // vecteur moyenne
    double **covariance;    // matrice de variance-covariance
    int **int_vector;       // vecteurs, variables entieres
    double **real_vector;   // vecteurs, variables reelles

    void init(int inb_vector , int *iidentifier , int inb_variable , int *itype , bool init_flag);
    void copy(const Vectors &vec);
    void remove();

    void build_real_vector(int variable = I_DEFAULT);

    void transcode(const Vectors &vec , int variable , int min_symbol , int max_symbol , int *symbol);
    void cluster(const Vectors &vec , int variable , int nb_class , double *limit);
    void select_variable(const Vectors &vec , int *variable);
    Vectors* select_variable(int explanatory_variable , int response_variable) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive , bool comment_flag) const;
    std::ostream& ascii_print(std::ostream &os , bool comment_flag) const;
    bool plot_print(const char *path , double *standard_residual = 0) const;

    void min_value_computation(int variable);
    void max_value_computation(int variable);
    int* order_computation(int variable) const;
    void build_marginal_histogram(int variable);

    void mean_computation(int variable);
    void variance_computation(int variable);
    void covariance_computation(int variable = I_DEFAULT);

    double** correlation_computation() const;

    double** spearman_rank_correlation_computation() const;
    double** kendall_rank_correlation_computation() const;

    std::ostream& rank_correlation_ascii_write(std::ostream &os , int correlation_type ,
                                               double **correlation) const;
    bool rank_correlation_ascii_write(Format_error &error , const char *path , int correlation_type ,
                                      double **correlation) const;

    double spearman_rank_single_correlation_computation() const;
    double kendall_rank_single_correlation_computation() const;

    int** joint_frequency_computation(int variable1 , int variable2) const;

    std::ostream& contingency_table_ascii_write(std::ostream &os , int variable1 , int variable2 ,
                                                int **frequency , double **deviation ,
                                                double **chi2_contribution , const Test &test ,
                                                bool file_flag = false) const;
    bool contingency_table_ascii_write(Format_error &error , const char *path , int variable1 ,
                                       int variable2 , int **frequency , double **deviation ,
                                       double **chi2_contribution , const Test &test) const;
    bool contingency_table_spreadsheet_write(Format_error &error , const char *path , int variable1 ,
                                             int variable2 , int **frequency , double **deviation ,
                                             double **chi2_contribution , const Test &test) const;

    std::ostream& variance_analysis_ascii_write(std::ostream &os , int type , const Vectors **value_vec ,
                                                bool exhaustive = false) const;
    bool variance_analysis_ascii_write(Format_error &error , const char *path , int response_type ,
                                       const Vectors **value_vec , bool exhaustive = false) const;
    bool variance_analysis_spreadsheet_write(Format_error &error , const char *path ,
                                             int response_type , const Vectors **value_vec) const;

public :

    Vectors();
    Vectors (int inb_vector , int *iidentifier , int inb_variable , int *itype ,
             bool init_flag = false)
    { init(inb_vector , iidentifier , inb_variable , itype , init_flag); }
    Vectors(int inb_vector , int *iidentifier , int inb_variable , int **iint_vector);
    Vectors(int inb_vector , int *iidentifier , int inb_variable , double **ireal_vector);
    Vectors(const Vectors &vec , int inb_vector , int *index);
    Vectors(const Vectors &vec)
    { copy(vec); }
    virtual ~Vectors();
    Vectors& operator=(const Vectors &vec);

    Distribution_data* extract(Format_error &error , int variable) const;

    bool check(Format_error &error);

    Vectors* merge(Format_error &error, int nb_sample , const Vectors **ivec) const;
    Vectors* shift(Format_error &error , int variable , int shift_param) const;
    Vectors* shift(Format_error &error , int variable , double shift_param) const;
    Vectors* transcode(Format_error &error , int variable , int *symbol) const;
    Vectors* cluster(Format_error &error , int variable , int step) const;
    Vectors* cluster(Format_error &error , int variable , int inb_value ,
                     int *ilimit) const;
    Vectors* cluster(Format_error &error , int variable , int nb_class ,
                     double *ilimit) const;
    Vectors* scaling(Format_error &error , int variable , int scaling_coeff) const;
    Vectors* round(Format_error &error , int variable = I_DEFAULT) const;

    Vectors* value_select(Format_error &error , int variable , int imin_value ,
                          int imax_value , bool keep = true) const;
    Vectors* value_select(Format_error &error , int variable , double imin_value ,
                          double imax_value , bool keep = true) const;

    Vectors* select_individual(Format_error &error , int inb_vector , int *iidentifier ,
                               bool keep = true) const;
    Vectors* select_variable(Format_error &error , int inb_variable , int *ivariable ,
                             bool keep = true) const;
    Vectors* merge_variable(Format_error &error , int nb_sample , const Vectors **ivec ,
                            int ref_sample = I_DEFAULT) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_data_write(std::ostream &os , bool exhaustive = false, 
				 bool comment_flag=false) const;
    bool ascii_data_write(Format_error &error , const char *path , bool exhaustive = false) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,
                    const char *title = 0) const;

/*    RWDECLARE_COLLECTABLE(Vectors);

    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

    double mean_absolute_deviation_computation(int variable) const;
    double mean_absolute_difference_computation(int variable) const;
    double skewness_computation(int variable) const;
    double kurtosis_computation(int variable) const;

    bool rank_correlation_computation(Format_error &error , std::ostream &os ,
                                      int correlation_type , const char *path = 0) const;

    Distance_matrix* comparison(Format_error &error , const Vector_distance &ivector_dist ,
                                bool standardization = true) const;

    bool contingency_table(Format_error &error , std::ostream &os , int variable1 , int variable2 ,
                           const char *path = 0 , char format = 'a') const;

    bool variance_analysis(Format_error &error , std::ostream &os , int class_variable ,
                           int response_variable , int response_type , const char *path = 0 ,
                           char format = 'a') const;

    Regression* linear_regression(Format_error &error , int explanatory_variable ,
                                  int response_variable) const;
    Regression* moving_average(Format_error &error , int explanatory_variable ,
                               int response_variable , int nb_point ,
                               double *filter , char algorithm = 'a') const;
    Regression* moving_average(Format_error &error , int explanatory_variable ,
                               int response_variable , const Distribution &dist ,
                               char algorithm = 'a') const;
    Regression* nearest_neighbor_smoother(Format_error &error , int explanatory_variable ,
                                          int response_variable , double span , bool weighting = true) const;

    // acces membres de la classe

    int get_nb_vector() const { return nb_vector; }
    int get_identifier(int ivec) const { return identifier[ivec]; }
    int get_nb_variable() const { return nb_variable; }
    int get_type(int variable) const { return type[variable]; }
    double get_min_value(int variable) const { return min_value[variable]; }
    double get_max_value(int variable) const { return max_value[variable]; }
    Histogram* get_marginal(int variable) const { return marginal[variable]; }
    double get_mean(int variable) const { return mean[variable]; }
    double get_covariance(int variable1, int variable2) const
    { return covariance[variable1][variable2]; }
    int get_int_vector(int ivec , int variable) const
    { return int_vector[ivec][variable]; }
    double get_real_vector(int ivec , int variable) const
    { return real_vector[ivec][variable]; }
};


Vectors* vectors_ascii_read(Format_error &error , const char *path);



class Vector_distance : public STAT_interface {  // parametres de definition
                                                 // d'une distance entre vecteurs

    friend class Vectors;
    friend class Sequences;

    friend class TreeMatch;

    friend Vector_distance* vector_distance_ascii_read(Format_error &error , const char *path);
    friend std::ostream& operator<<(std::ostream &os , const Vector_distance &param)
    { return param.ascii_write(os); }

private :

    int nb_variable;        // nombre de variables
    int distance_type;      // type de distance (ABSOLUTE_VALUE/QUADRATIC)
    int *variable_type;     // type de chaque variable (SYMBOLIC/ORDINAL/NUMERIC/CIRCULAR)
    double *weight;         // poids de chaque variable
    double *dispersion;     // quantite pour la standardisation
    int *nb_value;          // nombre de valeurs par variable
    double ***symbol_distance;  // matrice des distances entre symboles
    int *period;            // periode (variable circulaire)

    void copy(const Vector_distance &param);
    void remove();

public :

    Vector_distance();
    Vector_distance(int inb_variable , int *ivariable_type , double *iweight ,
                    int idistance_type = ABSOLUTE_VALUE);
    Vector_distance(int inb_variable , int idistance_type , int *ivariable_type ,
                    double *iweight , int *inb_value , double ***isymbol_distance ,
                    int *iperiod);
    Vector_distance(const Vector_distance &vector_dist)
    { copy(vector_dist); }
    virtual ~Vector_distance();
    Vector_distance& operator=(const Vector_distance &vector_dist);

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;

    // fonctions pour la compatibilite avec la classe STAT_interface

    bool ascii_write(Format_error &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,
                    const char *title = 0) const;

/*    RWDECLARE_COLLECTABLE(Vector_distance);

    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

    double* max_symbol_distance_computation(int variable) const;

    void dispersion_computation(int variable , const Histogram *marginal , double *rank = 0) const;

    // acces membres de la classe

    int get_nb_variable() const { return nb_variable; }
    int get_distance_type() const { return distance_type; }
    int get_variable_type(int variable) const { return variable_type[variable]; }
    double get_weight(int variable) const { return weight[variable]; }
    double get_dispersion(int variable) const { return dispersion[variable]; }
    int get_nb_value(int variable) const { return nb_value[variable]; }
    double get_symbol_distance(int variable , int symbol1 , int symbol2) const
    { return symbol_distance[variable][symbol1][symbol2]; }
    int get_period(int variable) const { return period[variable]; }
};


Vector_distance* vector_distance_ascii_read(Format_error &error , const char *path);



#endif
