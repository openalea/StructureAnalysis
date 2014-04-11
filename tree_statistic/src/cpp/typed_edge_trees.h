/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2010 CIRAD/INRIA Virtual Plants
 *
 *       File author(s): J.-B. Durand (jean-baptiste.durand@imag.fr)
 *
 *       $Source$
 *       $Id: typed_edge_trees.h 3193 2007-05-29 10:03:19Z dufourko $
 *
 *       Forum for OpenAlea developers: Openalea-devlp@lists.gforge.inria.fr
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
 *       MERCHANTABILITY or FITNESS for A PARTICULAR PURPOSE. See the
 *       GNU General Public License for more details.
 *
 *       You should have received a copy of the GNU General Public
 *       License along with this program; see the file COPYING. If not,
 *       write to the Free Software Foundation, Inc., 59
 *       Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 *  ----------------------------------------------------------------------------
 */

#ifndef TYPED_EDGE_TREES_H
#define TYPED_EDGE_TREES_H

#include "int_fl_containers.h"
#include "generic_typed_edge_tree.h"
#include "stat_tool/stat_tools.h"
#include "stat_tool/stat_label.h"
#include "stat_tool/curves.h"    // definition of Curves class for Sequences
#include "stat_tool/markovian.h" // definition of constants
#include "sequence_analysis/sequences.h" // definition of constants
#include "tree_labels.h"
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <float.h>
#include <limits.h>

extern char* label(const char*);
extern int* identifier_select(int nb_pattern , int* pattern_identifier,
                              int selected_nb_pattern, int* selected_identifier,
                              bool keep);
extern int* select_variable(int nb_variable, int selected_nb_variable,
                            int* selected_variable, bool keep);

namespace Stat_trees
{
/*! \file typed_edge_trees.h
    \brief Purpose :
              provide a specialization of the class Typed_edge_tree for
                 the Int_fl_container and
              the class "Typed_edge_trees" used for the representation of a
                 set of oriented trees
*/

/*--------------------------------------------------------------*
 *
 *  Constants :
 */
/// maximal number of simulated trees
const int NB_TREES= 100000;
/// maximal size of simulated trees
const int MAX_SIZE= 100000;
/// maximal cumulated size of simulated trees
const int CUMUL_SIZE= 100000;
/// maximal number of children for simulated trees
const int MAX_CHILDREN= 10;
/// maximal number of children for simulated trees
const int CUMUL_NB_CHILDREN= 1000000;
/// maximal tree size for the counting distribution to be computed
const int COUNTING_MAX_SIZE= 500;
/// maximal number of viewable trees (Gnuplot output)
const int PLOT_NB_TREES= 200;
/// maximal number of identified trees (Gnuplot output)
const int PLOT_TITLE_NB_TREES= 15;

/*--------------------------------------------------------------*
 *
 *  Type definitions
 */

/// nature of the characteristic quantities for trees
enum {
  FIRST_OCCURRENCE_ROOT ,
  FIRST_OCCURRENCE_LEAVES ,
  SOJOURN_SIZE ,
  NB_ZONES ,
  NB_OCCURRENCES ,
  OBSERVATION
};

/*--------------------------------------------------------------*
 *
 *  Class definitions :
 */

template<typename Generic_Int_fl_container> class Typed_edge_int_fl_tree;
typedef Typed_edge_int_fl_tree<One_int_container> Typed_edge_one_int_tree;
typedef Typed_edge_one_int_tree** pt_One_int_tree_array;
class TreeCharacteristics;
class HiddenMarkovIndOutTree;

/**
   \class Edged_typed_int_fl_tree
   \brief specialization of the Typed_edge_tree class for
          handling attributes with floating and integer types
*/

template<class Generic_Int_fl_container>
class Typed_edge_int_fl_tree : public Typed_edge_tree<Generic_Int_fl_container>
{  // specialization of the Typed_edge_tree class for
   // the Int_fl_container

protected :

    int _nb_integral;
    int _nb_float;

    void copy(const Typed_edge_int_fl_tree<Generic_Int_fl_container>& tree);

public :

    typedef Typed_edge_tree<Generic_Int_fl_container> tree_type;
    typedef typename tree_type::value value;

    typedef typename tree_type::key key;
    typedef typename tree_type::vertices_size_type vertices_size_type;
    typedef typename tree_type::children_iterator children_iterator;
    typedef typename tree_type::vertex_iterator vertex_iterator;

    typedef typename generic_visitor<tree_type>::vertex_array vertex_array;

    typedef Distribution** pt_Distribution_array;

    Typed_edge_int_fl_tree(key root= 0,
                int n= 1,
                const value& default_value= value());
    Typed_edge_int_fl_tree(int inb_integral,
                int inb_float,
                key root,
                int n);
    Typed_edge_int_fl_tree(const Typed_edge_int_fl_tree& tree);
    Typed_edge_int_fl_tree(Unlabelled_typed_edge_tree& utree, const value& default_value= value());
    Typed_edge_int_fl_tree(const Typed_edge_tree<Generic_Int_fl_container>& tree);
    Typed_edge_int_fl_tree(int inb_integral,
                int inb_float,
                Unlabelled_typed_edge_tree& utree);
    ~Typed_edge_int_fl_tree();
    Typed_edge_int_fl_tree<Generic_Int_fl_container>&
       operator=(const Typed_edge_int_fl_tree<Generic_Int_fl_container>& tree);

    key add_vertex(const value& data= value());

    Typed_edge_one_int_tree* select_int_variable(int variable);
    // selection of one single integral variable

    Generic_Int_fl_container* get_min_value(); // const
    Generic_Int_fl_container* get_max_value(); // const

    void iid_simulation(pt_Distribution_array distrib);
    // random determination of the integral labels;

    int get_nb_int() const;
    int get_nb_float() const;

    void set_structure(Unlabelled_typed_edge_tree& utree,
                       const value& default_value= value());
};

/**
   \class Typed_edge_trees
   \brief set of tree-structured data for handling attributes
          with floating and integer types
*/

template <typename Generic_Int_fl_container>
class Typed_edge_trees : public StatInterface
{  // set of tree-structured data assumed to be the realizations
   // of a common mixed multidimensional tree process
   // (nb_integral of integer type, nb_float of float type)


   // friend classes
   // friend class TreeCharacteristics;
   friend class HiddenMarkovTree;
   friend class HiddenMarkovIndOutTree;

   // template <typename Type> friend
   //    Typed_edge_trees<Type>* observed_trees_ascii_read(StatError& error,
   //                                                     const char * path);
   template <typename Type> friend
      std::ostream& operator<<(std::ostream& os, const Typed_edge_trees<Type>& otrees);

public :

   typedef Typed_edge_int_fl_tree<Generic_Int_fl_container> tree_type;
   typedef typename tree_type::value value;

   typedef typename tree_type::key key;
   typedef typename tree_type::vertices_size_type vertices_size_type;
   typedef typename tree_type::children_iterator children_iterator;
   typedef typename tree_type::vertex_iterator vertex_iterator;

   typedef TreeCharacteristics* TreeCharacteristics_array;
   typedef TreeCharacteristics** pt_TreeCharacteristics_array;
   typedef tree_type** pt_tree_type_array;
   typedef Typed_edge_trees<Generic_Int_fl_container> observed_trees;
   typedef observed_trees** pt_observed_trees_array;

   typedef int* int_array;
   typedef double* double_array;
   typedef double** double_array_2d;

protected :

   /// number of variables of integer type
   /// these variables come before that of floating type;
   int _nb_integral;
   /// number of variables of floating type;
   /// these variables come after that of integer type;
   int _nb_float;
   /// type of each of those variables
   /// (INT_VALUE / REAL_VALUE/ TIME / POSITION / NB_INTERNODE )
   /// ending with the floating variables
   int_array _type;
   /// minimal value of each variable
   Generic_Int_fl_container _min_value;
   /// maximal value of each variable
   Generic_Int_fl_container _max_value;
   /// maximal size of the trees
   unsigned int _max_size;
   /// maximal depth of the trees
   unsigned int _max_depth;
   /// number of observed trees
   int _nb_trees;
   // ? int *identifier;                // identifier of each tree

   // [quantities for which the characteristic distributions are invariant
   // - if any - to determine]
   // [frequency distributions of those quantities]
   /// frequency distribution of the tree size
   FrequencyDistribution *hsize;
   /// frequency distribution of the number of children
   FrequencyDistribution *hnb_children;

   /// observed trees
   pt_tree_type_array trees;

   /// frequency distribution of the characteristic quantities
   /// for each discrete variable
   pt_TreeCharacteristics_array characteristics;

   Typed_edge_trees(int inb_integral, int inb_float, int_array itype,
                  int inb_trees, bool init_flag= true);

   void remove();
   // deallocator
   void copy(const Typed_edge_trees<Generic_Int_fl_container>& otrees,
             bool characteristic_flag= true);
   // copy operator
   void init(int inb_integral, int inb_float, int_array itype,
             int inb_trees, bool init_flag);

   void cluster(const Typed_edge_trees& otrees, int ivariable, int step);
   void transcode(const Typed_edge_trees& otrees, int ivariable,
                  int min_symbol, int max_symbol, int_array symbol);
   void select_variable(const Typed_edge_trees& otrees, int_array variable,
                        int source_nb_integral);

    // ostream& ascii_print(ostream& os, Typed_edge_tree<double>** real_tree_sequence, char format);
//    // Typed_edge_tree<double>** real_tree_sequence : supposedly to be revised
//    bool ascii_print(StatError &error, const char* path, Typed_edge_tree<double>** real_tree_sequence, char format);
//    // idem
//    int scaling_coefficient(int variable,  Typed_edge_tree<double>** real_tree_sequence);
//    // idem
   /** Gnuplot output of \e self in the case of no other characteristic
       than the marginal frequency distribution*/
   bool plot_print(const char * prefix, const char * title,
                   int variable, int nb_variable) const;

   /** Compute min and max values of variables */
   virtual void min_max_value_computation();
   /** Compute maximal number of vertices */
   virtual void max_size_computation();
   /** Compute maximal depth */
   virtual void max_depth_computation();
   /** Compute maximal number of children */
   virtual int max_nb_children_computation(); // const
   std::ostream& ascii_write(std::ostream& os, bool exhaustive, bool comment_flag) const;
   std::ostream& ascii_print(std::ostream& os, char format, bool comment_flag,
                             int line_nb_character= LINE_NB_CHARACTER) const;
//    bool plot_print(const char *path, int ilength) const;

   // [computation of quantities for which the characteristic distributions are invariant]
   // [computation of their frequency distributions]

   // void state_variable_init(int istate_variable);

   /** Compute total number of vertices */
   virtual int cumul_size_computation() const;
   /** Compute histogram of sizes */
   virtual void build_size_frequency_distribution();
   /** Compute total number of children */
   virtual int cumul_nb_children_computation() const;
   /** Compute histogram for the number of children */
   virtual void build_nb_children_frequency_distribution();

   /** Compute statistical properties of trees for a given variable */
   virtual void build_characteristics(int variable);
   /** Compute statistical properties of trees for every variable */
   virtual void build_characteristics();
   // computation of the characteristic quantity frequency distributions for each discrete variable
   // except the observation frequency distributions

public :

   Typed_edge_trees(int inb_integral= 1, int inb_float= 0, int inb_trees= 0);
   Typed_edge_trees(const Typed_edge_trees& otrees, bool characteristic_flag= true);
   Typed_edge_trees(int inb_trees,
                    int const * itype,
                    pt_tree_type_array otrees,
                    bool frequency_distribution_flag= true);
   Typed_edge_trees(int inb_integral,
                    int inb_float,
                    const FrequencyDistribution& ihsize,
                    const FrequencyDistribution& ihnb_children,
                    bool no_child_flag= false,
                    bool init_flag= true);
   Typed_edge_trees(const Typed_edge_trees& otrees, int inb_trees, int_array index);
   virtual ~Typed_edge_trees();
   observed_trees& operator=(const observed_trees& otrees);

   DiscreteDistributionData* extract(StatError& error, int variable) const;

   /** Extract the frequency distribution of characteristic distribution from \e self*/
   DiscreteDistributionData* extract(StatError& error, int type,
                                     int variable, int value) const;
   /** Extract all vertex values from \e self*/
   Vectors* build_vectors(StatError& error) const;
   // should be a constructor of Vectors ?
   /** Extract all sequences from \e self along paths*/
   Sequences* build_sequences(StatError& error, bool all_paths= true,
                              bool auto_axis= false) const;

    // Distributions of the tree-characteristics
    // ?    Vectors* extract_vectors(StatError &error, int type, int variable = I_DEFAULT,
    //                         int value = I_DEFAULT) const;

    // ? Tops* tops(StatError &error) const;

    // ? bool check(StatError &error, const char *pattern_label);

   Typed_edge_trees* merge(StatError& error, int nb_sample,
                           const pt_observed_trees_array otrees) const;
   Typed_edge_trees* shift(StatError& error, int variable, double shift_param) const;
   Typed_edge_trees* shift(StatError& error, int variable, int shift_param) const;

   /** Cluster the values of a variable of \e self (using a given step)*/
   Typed_edge_trees* cluster(StatError& error, int variable, int step) const;
   /** Transcode the values of a variable of \e self (using an array of new values)*/
   Typed_edge_trees* transcode(StatError& error, int variable, int_array symbol) const;

   /** Cluster the values of a variable of \e self (using given classes and limits)*/
   Typed_edge_trees* cluster(StatError& error, int variable, int nb_class,
                             int_array ilimit) const;
//    Typed_edge_trees* cluster(StatError& error, int variable, int nb_class,
//                              double_array dlimit) const;
   /** Select subset of trees according to a range of values taken by a given integral variable */
   virtual Typed_edge_trees* value_select(StatError& error, int variable, int imin_value,
                                          int imax_value, bool keep= true) const;
   /** Select subset of trees according to a range of values taken by a given floating variable */
   virtual Typed_edge_trees* value_select(StatError& error, int variable, double dmin_value,
                                          double dmax_value, bool keep= true) const;
   /** Select subset of trees from their indices */
   virtual Typed_edge_trees* select_individual(StatError& error, int inb_trees,
                                               int_array iidentifier, bool keep= true) const;
   /** Select subset of variables from their indices */
   virtual Typed_edge_trees* select_variable(StatError& error, int inb_variable,
                                             int_array ivariable, bool keep= true) const;
   /** Add variables to trees */
   virtual Typed_edge_trees* merge_variable(StatError& error, int nb_sample,
                                            const pt_observed_trees_array otrees,
                                            int ref_sample= I_DEFAULT) const;
   /** Select subset of trees according to a range of sizes */
   virtual Typed_edge_trees* size_select(StatError& error, int imin_size, int imax_size,
                                         bool keep= true) const;

   /** Select subtrees using the subtree root index, the tree index
       and a flag on keeping or pruning the subtree */
   virtual Typed_edge_trees* select_subtrees(StatError& error, int iindex,
                                             int itree= I_DEFAULT, bool keep= true) const;
   // would be more interesting with true vertex and tree identifiers

//     Sequences* extract_sequences(int* tree_list) const;
//     // extraction of all possible (for instance discrete) sequences from the list of trees tree_list
//     // starting from tree root

   /** Extraction of homogeneous trees */
   Typed_edge_trees* segmentation_extract(StatError& error, int variable,
                                          int nb_value, int *ivalue,
                                          bool keep = true) const;

//     // ? Typed_edge_trees* cumulate(StatError &error, int variable = I_DEFAULT) const;
   /** First-order differentiation for trees */
   Typed_edge_trees* difference(StatError &error, int variable= I_DEFAULT) const;

   /** Type conversion for trees */
   void to_int_type(StatError &error, int variable= I_DEFAULT);
//     // Typed_edge_trees* moving_average(StatError &error, ostream& os, int nb_point, double *filter,
//     //                          int variable = I_DEFAULT, bool begin_end = false,
//     //                          int output = TREND, const char *path = 0, char format = 's') const;
//     // Typed_edge_trees* moving_average(StatError &error, ostream& os, const Distribution &dist,
//     //                         int variable = I_DEFAULT, bool begin_end = false,
//     //                         int output = TREND, const char *path = 0, char format = 's') const;
//     // somewhat tricky in the case of oriented trees : postponed or cancelled

//     Typed_edge_trees* transform_position(StatError &error, int step) const;
//     // Discretization of a variable with type POSITION

//     // Typed_edge_trees* scaling(StatError &error, int variable, int scaling_coeff) const;
//     // not most useful, unless we want to use 'moving_average' in the case of a discrete
//     // variables taking a small number of values

//     // Typed_edge_trees* cross(StatError &error) const;

    std::ostream& line_write(std::ostream& os) const;

    virtual ostream& ascii_data_write(ostream& os, char format= 'c', bool exhaustive= false) const;
    virtual bool ascii_data_write(StatError& error, const char * path,
                                  char format= 'c', bool exhaustive= false) const;
    bool plot_data_write(StatError& error, const char * prefix,
                         const char * title= NULL) const;

    std::ostream& ascii_write(std::ostream& os, bool exhaustive= false) const;
    bool ascii_write(StatError& error, const char * path,
                     bool exhaustive= false) const;
    bool spreadsheet_write(StatError& error, const char * path) const;

    /** Gnuplot output of \e self*/
    bool plot_write(StatError& error, const char * prefix,
                    const char * title= NULL) const;

    /** Graphical output of \e self*/
    MultiPlotSet* get_plotable(StatError& error,
                               int plot_type,
                               int variable) const;

    std::ostream& ascii_write_size_frequency_distribution(std::ostream& os,
                                                          bool exhaustive,
                                                          bool file_flag) const;

    std::ostream& ascii_write_nb_children_frequency_distribution(std::ostream& os,
                                                                 bool exhaustive,
                                                                 bool file_flag) const;

    double mean_computation(int variable) const;
    double variance_computation(int variable, double mean) const;
    // double mean_absolute_deviation_computation(int variable, double mean) const;
    // double mean_absolute_difference_computation(int variable) const;
    double skewness_computation(int variable, double mean, double variance) const;
    double kurtosis_computation(int variable, double mean, double variance) const;
   // computation of statistical indicators associated with each variable

   // Correlation* correlation_computation(StatError &error, int variable1, int variable2,
   //                                     int itype = PEARSON, int max_lag = I_DEFAULT,
   //                                      int normalization = APPROXIMATED) const;
   // Correlation* partial_autocorrelation_computation(StatError &error, int variable,
   //                                                 int itype = PEARSON, int max_lag = I_DEFAULT) const;
   // is this relevant for tree-structured data ? Postponed or cancelled

   double iid_information_computation() const;

   // access to the class members

   // int get_nb_variable() const;
   int get_nb_int() const;
   int get_nb_float() const;

   /** Return total number of vertices */
   virtual unsigned int get_total_size() const;

   /** Return the type of a given variable */
   int get_type(int variable) const;
   /** Return the list of types (return a pointer; object should not be modified) */
   const int* get_types_ptr() const;

   int get_min_int_value(int variable) const;
   double get_min_fl_value(int variable) const;

   int get_max_int_value(int variable) const;
   double get_max_fl_value(int variable) const;

   int get_max_size() const;
   int get_max_depth() const;

   /**Return the size of a given tree of \e self */
   virtual unsigned int get_size(int index) const;

   pt_One_int_tree_array get_identifier_trees();
   Typed_edge_one_int_tree* get_identifier_tree(int index);

   virtual DiscreteDistributionData* extract_size() const;
   virtual DiscreteDistributionData* extract_nb_children() const;
   // access to the frequency distributions of topological indicators

   /** Return the maximal value for each variable of \e self*/
   Generic_Int_fl_container get_max_value() const;
   /** Return the minimal value for each variable of \e self*/
   Generic_Int_fl_container get_min_value() const;
   /** return the number of values of a variable with finite values*/
   int get_nb_values(int variable) const;
   FrequencyDistribution* get_marginal(int variable) const;
   int get_nb_trees() const;
    // ? int get_identifier(int itree) const { return identifier[itree]; }

    // [access to the quantities for which the characteristic distributions are invariant]
    // [access to their frequency distributions]

   TreeCharacteristics_array* get_characteristics() const;
   TreeCharacteristics* get_characteristics(int variable) const;
   // access to the characteristic quantity distributions
   // for the observed variables

   // check whether given characteristic is present
   bool is_characteristic(int variable, int charac) const;

   // ? Tree_curves* get_index_value(int variable) const;
   // ? Tree_curves** get_index_value() const;
   // ? Curves* get_depth_value(int variable) const;
   // ? Curves** get_depth_value() const;
   FrequencyDistribution* get_first_occurrence_root(int variable, int value) const;
   FrequencyDistribution* get_first_occurrence_leaves(int variable, int value) const;
   FrequencyDistribution* get_sojourn_size(int variable, int value) const;
   FrequencyDistribution* get_nb_zones(int variable, int value) const;
   FrequencyDistribution* get_nb_occurrence(int variable, int value) const;

   /** Return the set of trees (a new instance is allocated) */
   pt_tree_type_array get_trees() const;
   /** Return the set of trees (return a pointer; object should not be deallocated) */
   pt_tree_type_array get_trees_ptr() const;
   /** Return a given tree (a new instance is allocated) */
   tree_type* get_tree(int itree) const;
   /** Return a given tree (return a pointer; object should not be deallocated) */
   tree_type* get_tree_ptr(int itree) const;

};


class TreeCharacteristics
{  //  set of frequency distributions for characteristic
   //  quantity distributions

   // friend classes

   // friend class Observed_trees<Int_fl_container>;
   // supposedly working on Windows compilers ?

   template <typename Generic_Int_fl_container> friend class Typed_edge_trees;
   friend class Observed_int_trees;
   friend class HiddenMarkovTree;
   friend class HiddenMarkovTreeData;
   friend class CategoricalTreeProcess;
   friend class MarkovOutTree;
   friend class MarkovOutTreeData;

public :

   typedef FrequencyDistribution** ptFrequencyDistribution_array;
   typedef Histogram** ptHistogram_array;

protected :

   /// index of the variable. 0 means a state variable
   int _variable;
   /// minimal value of the variable
   int _min_value;
   /// maximal value of the variable
   int _max_value;
   /// number of trees used for the frequency distributions
   int _nb_trees;
   /// FrequencyDistribution of the (marginal) discrete observation distribution
   FrequencyDistribution *marginal_distribution;
   /// Histogram of the (marginal) continuous observation distribution
   Histogram *marginal_histogram;
   // ?  Tree_curves *index_value;   // empirical probability of each value, as
   // a function of the index
   // ? Curves * depth_value;  empirical probability of each value, as
   // a function of the tree depth
   /// frequency distribution of the path length before the first occurrence
   /// of a given value, starting from the root
   ptFrequencyDistribution_array first_occurrence_root;
   /// frequency distribution of the path length before the first occurrence
   /// of a given value, starting from all the terminal vertices
   ptFrequencyDistribution_array  first_occurrence_leaves;
   /// frequency distribution of the sizes of the homogeneous zones
   /// for a given value of the variable
   ptFrequencyDistribution_array sojourn_size;
   /// frequency distribution of the number of homogeneous zones in the tree
   /// for a given value of the variable
   ptFrequencyDistribution_array nb_zones;
   /// frequency distribution of the number of occurrences for a given value
   // of the variable in the whole tree
   ptFrequencyDistribution_array nb_occurrences;

   // Histogram2D *children_pairs; // 2D frequency distribution of the pairs of children of vertices

   // ?  void init(bool initial_run_flag = false);
   // ?   void copy_characteristics(int variable , const Discrete_sequences &seq ,
   //                             int cvariable , bool initial_run_flag = false);
   // ?   void copy(const TreeCharacteristics &tchar , int param = DEFAULT);
   // ?   void add_variable(const TreeCharacteristics &tchar , int variable = 0);
   void remove_characteristic(ptFrequencyDistribution_array& c, int inb_values);

   void remove();
   void copy(const TreeCharacteristics& t_char);

   void init_characteristic(ptFrequencyDistribution_array& c, int inb_values, int imax_val);
   void init(int imax_size, int imax_depth);
   // allocation of the frequency distributions

   // ? void build_index_value();
   // ? void build_depth_value();

   FrequencyDistribution* get_characteristic(int value, ptFrequencyDistribution_array h) const;

   // check whether given characteristic is present
   bool is_characteristic(int charac) const;

   void build_marginal_frequency_distribution(Typed_edge_one_int_tree** otrees1, bool final_value= false);
   void build_first_occurrence_root_frequency_distribution(Typed_edge_one_int_tree** otrees1,
                                                           int imax_depth, bool final_value= false);
   void build_first_occurrence_leaves_frequency_distribution(Typed_edge_one_int_tree** otrees1,
                                                             int imax_depth, bool final_value= false);
   void build_zone_frequency_distributions(Typed_edge_one_int_tree** otrees1, int imax_size, bool final_value= false);
   void build_nb_occurrences_frequency_distribution(Typed_edge_one_int_tree** otrees1,
                                                    int imax_size, bool final_value= false);
// template <typename Generic_Int_fl_container>
// void build_children_pairs_frequency_distribution(Typed_edge_one_int_tree** otrees1, bool final_value= false);

   void ascii_write_characteristic(ptFrequencyDistribution_array c,
                                   int inb_values,
                                   std::ostream& os,
                                   std::string c1,
                                   std::string c2,
                                   bool exhaustive,
                                   bool file_flag) const;

   std::ostream& ascii_print(std::ostream& os, int type, const FrequencyDistribution& hsize,
                             bool exhaustive, bool comment_flag) const;
   std::ostream& spreadsheet_print(std::ostream& os, int type, const FrequencyDistribution& hsize) const;

  /** Gnuplot output of \e self*/
   bool plot_print(const char * prefix, const char * title, int variable,
                   int nb_variables, int type, const FrequencyDistribution& hsize) const;

  /** Graphical output of \e self*/
  MultiPlotSet* get_plotable(int plot_type,
                             int variable,
                             int nb_variables,
                             int type) const;

public :

   TreeCharacteristics();
   TreeCharacteristics(int imin_value,
                       int imax_value,
                       int imax_size,
                       int imax_depth,
                       int inb_trees,
                       Typed_edge_one_int_tree** otrees1,
                       int ivariable= I_DEFAULT,
                       bool final_value= false);

   TreeCharacteristics(const TreeCharacteristics& t_char);
   ~TreeCharacteristics();
   TreeCharacteristics&  operator=(const TreeCharacteristics &t_char);

   void  build_characteristics(int imin_value,
                               int imax_value,
                               int imax_size,
                               int imax_depth,
                               int inb_trees,
                               Typed_edge_one_int_tree** otrees1,
                               int ivariable= I_DEFAULT,
                               bool final_value= false);
   // construction of the frequency distributions

   // access to class members

   int get_variable() const;
   int get_nb_values() const;
   // ? Tree_curves* get_index_value() const
   //  { return index_value; }
   // ? Curves* get_depth_value() const
   //  { return depth_value; }
   FrequencyDistribution* get_marginal_distribution() const;
   FrequencyDistribution* get_first_occurrence_root(int value) const;
   FrequencyDistribution* get_first_occurrence_leaves(int value) const;
   FrequencyDistribution* get_sojourn_size(int value) const;
   FrequencyDistribution* get_nb_zones(int value) const;
   FrequencyDistribution* get_nb_occurrences(int value) const;

   int get_nb_value_first_occurrence_root(int value) const;
   int get_nb_value_sojourn_size(int value) const;

   std::ostream& ascii_write_marginal_distribution(std::ostream& os,
                                                   bool exhaustive,
                                                   bool file_flag) const;

   std::ostream& ascii_write_first_occurrence_root(std::ostream& os,
                                                   bool exhaustive,
                                                   bool file_flag) const;

   std::ostream& ascii_write_first_occurrence_leaves(std::ostream& os,
                                                     bool exhaustive,
                                                     bool file_flag) const;

   std::ostream& ascii_write_sojourn_size(std::ostream& os,
                                          bool exhaustive,
                                          bool file_flag) const;

   std::ostream& ascii_write_nb_zones(std::ostream& os,
                                      bool exhaustive,
                                      bool file_flag) const;

   std::ostream& ascii_write_nb_occurrences(std::ostream& os,
                                            bool exhaustive,
                                            bool file_flag) const;

};

typedef Typed_edge_int_fl_tree<Int_fl_container> Default_tree;
typedef Typed_edge_trees<Int_fl_container> Trees;
// typedef select_typed_edge_subtree<Int_fl_container> ;

#include "typed_edge_trees.hpp"
}; // end namespace

#endif
