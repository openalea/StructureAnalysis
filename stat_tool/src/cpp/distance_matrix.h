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



#ifndef DISTANCE_MATRIX_H
#define DISTANCE_MATRIX_H



/****************************************************************
 *
 *  Constantes :
 */

const int ASCII_NB_PATTERN = 10;       // nombre de formes maximum pour afficher sous forme
                                       // de resultats d'alignement
const double DISTANCE_ROUNDNESS = 1.e-12;  // arrondi sur une distance

const int GLOBAL_NB_ITER = 20;         // nombre d'iterations ou les groupes
                                       // sont recalculees globalement
const int PARTITIONING_NB_ITER_1 = 50;  // nombre maximum d'iterations
const int PARTITIONING_NB_ITER_2 = 20;  // nombre maximum d'iterations

enum {
  NEAREST_NEIGHBOR ,
  FARTHEST_NEIGHBOR ,
  AVERAGING
};

enum {
  CHILD_CLUSTER_DISTANCE ,
  DIAMETER
};



/****************************************************************
 *
 *  Definition de la classe :
 */

class Sequences;
class Clusters;
class Dendrogram;
class Distance_matrix;

class Distance_matrix : public STAT_interface {  // matrice des distances

    friend class Clusters;
    friend class Dendrogram;

    friend std::ostream& operator<<(std::ostream &os , const Distance_matrix &dist_matrix)
    { return dist_matrix.ascii_write(os); }

protected :

    int nb_row;             // nombre de lignes
    int nb_column;          // nombre de colonnes
    int *row_identifier;    // identificateurs des lignes
    int *column_identifier;  // identificateurs des colonnes
    double **distance;      // distances
    int **length;           // longueurs
    double **deletion_distance;  // distances d'elision
    int **nb_deletion;      // nombres d'elisions
    double **insertion_distance;  // distance d'insertion
    int **nb_insertion;     // nombres d'insertions
    int **nb_match;         // nombres de matchs
    double **substitution_distance;  // distances de substitution
    int **nb_substitution;  // nombres de substitutions
    double **transposition_distance;  // distances de transposition
    int **nb_transposition;  // nombres de transpositions
    int label_size;         // taille du label
    char *label;            // label

    void copy(const Distance_matrix &dist_matrix , char transform = 'c');
    void remove();

    std::ostream& property_print(double **normalized_distance , std::ostream &os ,
                                 char format) const;

    int cumul_length_computation(bool *row_flag , bool *column_flag) const;
    double cumul_distance_computation(bool *row_flag , bool *column_flag) const;
    plotable::MultiPlotSet* get_plotable(Format_error &) const;

public :

    Distance_matrix();
    Distance_matrix(int nb_pattern , const char *ilabel , int *pattern_identifier = 0);
    Distance_matrix(int nb_pattern , int irow_identifier , int icolumn_identifier ,
                    const char *ilabel , int *pattern_identifier = 0 ,
                    bool substitution_flag = true , bool transposition_flag = false);
    Distance_matrix(const Distance_matrix &dist_matrix , int inb_pattern ,
                    int *iidentifier , bool keep = true);
    Distance_matrix(const Distance_matrix &dist_matrix , int nb_cluster ,
                    const char *ilabel);
    Distance_matrix(const Distance_matrix &dist_matrix , char transform = 'c')
    { copy(dist_matrix , transform); }
    virtual ~Distance_matrix();
    Distance_matrix& operator=(const Distance_matrix &dist_matrix);

    Distance_matrix* select_individual(Format_error &error , int inb_pattern ,
                                       int *iidentifier , bool keep = true) const;
    Distance_matrix* symmetrize(Format_error &error) const;
    Distance_matrix* unnormalize(Format_error &error) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path ,
                     bool exhaustive = false) const;
    std::ostream& spreadsheet_write(std::ostream &os) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,
                    const char *title = 0) const;
    plotable::MultiPlotSet* get_plotable() const;

/*    RWDECLARE_COLLECTABLE(Distance_matrix);

    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

    bool test_symmetry() const;

    void update(int irow_identifier , int icolumn_identifier , double idistance ,
                int alignment_length , double ideletion_distance , int inb_deletion ,
                double iinsertion_distance , int inb_insertion , int inb_match ,
                double isubstitution_distance = 0. , int inb_substitution = 0 ,
                double itransposition_distance = 0. , int inb_transposition = 0);
    void update(int irow_identifier , int icolumn_identifier , double idistance , int ilength);

    Clusters* partitioning(Format_error &error , std::ostream &os , int nb_cluster ,
                           int *prototype = 0 , int initialization = 1 , int algorithm = 1) const;
    Clusters* partitioning(Format_error &error , std::ostream &os , int nb_cluster ,
                           int *cluster_nb_pattern , int **cluster_pattern) const;

    Dendrogram* agglomerative_hierarchical_clustering(int algorithm , int criterion = AVERAGING) const;
    Dendrogram* divisive_hierarchical_clustering() const;

    bool hierarchical_clustering(Format_error &error , std::ostream &os ,
                                 int algorithm = AGGLOMERATIVE , int criterion = AVERAGING ,
                                 const char *path = 0 , char format = 'a') const;

    // acces membres de la classe

    int get_nb_row() const { return nb_row; }
    int get_nb_column() const { return nb_column; }
    int get_row_identifier(int row) const { return row_identifier[row]; }
    int get_column_identifier(int column) const { return column_identifier[column]; }
    double get_distance(int row , int column) const
    { return distance[row][column]; }
    int get_length(int row , int column) const
    { return length[row][column]; }
    double get_deletion_distance(int row , int column) const
    { return deletion_distance[row][column]; }
    int get_nb_deletion(int row , int column) const
    { return nb_deletion[row][column]; }
    double get_insertion_distance(int row , int column) const
    { return insertion_distance[row][column]; }
    int get_nb_insertion(int row , int column) const
    { return nb_insertion[row][column]; }
    int get_nb_match(int row , int column) const
    { return nb_match[row][column]; }
    double get_substitution_distance(int row , int column) const
    { return substitution_distance[row][column]; }
    int get_nb_substitution(int row , int column) const
    { return nb_substitution[row][column]; }
    double get_transposition_distance(int row , int column) const
    { return transposition_distance[row][column]; }
    int get_nb_transposition(int row , int column) const
    { return nb_transposition[row][column]; }
    int get_label_size() const { return label_size; }
    char* get_label() const { return label; }
};



class Clusters : public Distance_matrix {  // resultats du clustering par partitionnement

    friend class Distance_matrix;

    friend std::ostream& operator<<(std::ostream &os , const Clusters &clusters)
    { return clusters.ascii_write(os); }

private :

    Distance_matrix *distance_matrix;  // pointeur sur un objet Distance_matrix
    int nb_pattern;         // nombre de formes
    int nb_cluster;         // nombre de groupes
    int *cluster_nb_pattern;  // effectifs des groupes
    int *assignment;        // affectations des formes
    double **pattern_distance;  // distances d'une forme a un groupe
    int **pattern_length;   // longueurs niveau forme

    void copy(const Clusters &clusters);
    void remove();

    int* pattern_sort(int cluster) const;

    int most_distant_pattern_selection(double **normalized_distance , int ipattern) const;
    int neighbor_pattern_selection(double **normalized_distance , int ipattern) const;
    int neighbor_pattern_cluster_selection(double **normalized_distance , int pattern) const;

    double max_intra_cluster_distance_computation(double **normalized_distance , int cluster) const;
    double min_inter_cluster_distance_computation(double **normalized_distance , int cluster) const;
    bool isolation_property(double **normalized_distance , int cluster , char type = 'c') const;
    double inter_cluster_distance_computation(int cluster) const;
    std::ostream& global_distance_ascii_print(std::ostream &os);

    void prototype_initialization_1();
    void prototype_initialization_2();

    void pattern_distance_update(int pattern , int old_cluster , int new_cluster);
    void algorithmic_step_1();
    void algorithmic_step_2();

public :

    Clusters();
    Clusters(const Distance_matrix &dist_matrix , int inb_cluster);
    Clusters(const Distance_matrix &dist_matrix , int inb_cluster ,
             int *icluster_nb_pattern , int **cluster_pattern);
    Clusters(const Clusters &clusters)
    :Distance_matrix(clusters) { copy(clusters); }
    ~Clusters();
    Clusters& operator=(const Clusters &clusters);

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path , bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,
                    const char *title = 0) const;

/*    RWDECLARE_COLLECTABLE(Clusters);

    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

    void cluster_nb_pattern_computation();
    void pattern_distance_computation();
    void cluster_distance_computation_1();
    void cluster_distance_computation_2();

    // acces membres de la classe

    Distance_matrix* get_distance_matrix() { return distance_matrix; }
    int get_nb_pattern() const { return nb_pattern; }
    int get_nb_cluster() const { return nb_cluster; }
    int get_cluster_nb_pattern(int cluster) const { return cluster_nb_pattern[cluster]; }
    int get_assignment(int pattern) const { return assignment[pattern]; }
    double get_pattern_distance(int pattern , int cluster) const
    { return pattern_distance[pattern][cluster]; }
    int get_pattern_length(int pattern , int cluster) const
    { return pattern_length[pattern][cluster]; }
};



class Dendrogram : public STAT_interface {      // resultats du clustering hierarchique

    friend class Distance_matrix;
    friend class Sequences;

    friend std::ostream& operator<<(std::ostream &os , const Dendrogram &dendrogram)
    { return dendrogram.ascii_write(os); }

private :

    Distance_matrix *distance_matrix;  // pointeur sur un objet Distance_matrix
    int scale;              // echelle pour representer les distances entre groupes
    int nb_cluster;         // nombre de groupes
    int *cluster_nb_pattern;  // effectifs des groupes
    int **cluster_pattern;  // compositions des groupes
    int *parent;            // noeuds parent
    int **child;            // noeuds fils
    double *child_distance;  // distances entre les deux groupes fusionnes
    double *intra_cluster_distance;  // distances intra-groupe
    double *inter_cluster_distance;  // distances inter-groupe
    double *max_intra_cluster_distance;  // diametres
    double *min_inter_cluster_distance;  // separations

    void copy(const Dendrogram &dendrogram);
    void remove();

    double* distance_ordering() const;
    double coefficient_computation(int iscale = I_DEFAULT) const;
    void tree_computation();

public :

    Dendrogram();
    Dendrogram(const Distance_matrix &dist_matrix , int iscale);
    Dendrogram(const Dendrogram &dendrogram)
    { copy(dendrogram); }
    ~Dendrogram();
    Dendrogram& operator=(const Dendrogram &dendrogram);

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path , bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,
		    const char *title = 0) const { return false; }
    // acces membres de la classe

    Distance_matrix* get_distance_matrix() { return distance_matrix; }
    int get_scale() const { return scale; }
    int get_nb_cluster() const { return nb_cluster; }
    int get_cluster_nb_pattern(int cluster) const { return cluster_nb_pattern[cluster]; }
    int get_cluster_pattern(int cluster , int index) const { return cluster_pattern[cluster][index]; }
    int get_parent(int cluster) const { return parent[cluster]; }
    int get_child(int cluster , int index) const { return child[cluster][index]; }
    double get_child_distance(int cluster) const { return child_distance[cluster]; }
    double get_intra_cluster_distance(int cluster) const { return intra_cluster_distance[cluster]; }
    double get_inter_cluster_distance(int cluster) const { return inter_cluster_distance[cluster]; }
    double get_max_intra_cluster_distance(int cluster) const { return max_intra_cluster_distance[cluster]; }
    double get_min_inter_cluster_distance(int cluster) const { return min_inter_cluster_distance[cluster]; }
};



#endif
