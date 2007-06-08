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



#ifndef REGRESSION_H
#define REGRESSION_H



/****************************************************************
 *
 *  Constantes :
 */


const int REGRESSION_NB_VECTOR = 5000;  // nombre maximum de vecteurs pour la
                                        // regression non-parametrique
const int NEIGHBORHOOD = 3;            // voisinage minimum sur les valeurs
                                       // de la variable explicative

enum {
  STAT_LINEAR ,
  STAT_LOGISTIC ,
  STAT_MONOMOLECULAR ,
  STAT_NONPARAMETRIC
};



/****************************************************************
 *
 *  Definition des classes :
 */


class Regression_kernel {  // noyau de regression

protected :

    int ident;              // identificateur de la fonction de regression
    int min_value;          // valeur minimum
    int max_value;          // valeur maximum
    double regression_df;   // degres de liberte regression
    double residual_df;     // degres de liberte residus
    int nb_parameter;       // nombre de parametres
    double *parameter;      // parametres
    double *point;          // points

    void copy(const Regression_kernel&);
    void remove();

    std::ostream& ascii_parameter_print(std::ostream &os) const;
    std::ostream& ascii_formal_print(std::ostream &os) const;
    std::ostream& ascii_print(std::ostream &os) const;
    std::ostream& spreadsheet_print(std::ostream &os) const;
    bool plot_print(const char *path) const;

/*    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

    void computation();
    double min_computation() const;
    double max_computation() const;

public :

    Regression_kernel();
    Regression_kernel(int iident , int imin_value , int imax_value);
    Regression_kernel(const Regression_kernel &regression) { copy(regression); }
    ~Regression_kernel();
    Regression_kernel& operator=(const Regression_kernel &regression);

    // acces membres de la classe

    int get_ident() const { return ident; }
    int get_min_value() const { return min_value; }
    int get_max_value() const { return max_value; }
    double get_regression_df() const { return regression_df; }
    double get_residual_df() const { return residual_df; }
    int get_nb_parameter() const { return nb_parameter; }
    double get_parameter(int index) const { return parameter[index]; }
    double get_point(int index) const { return point[index]; }
};



class Vectors;

class Regression : public STAT_interface , public Regression_kernel {  // fonction de regression

    friend class Vectors;

    friend std::ostream& operator<<(std::ostream &os , const Regression &regression)
    { return regression.ascii_write(os); }

private :

    Vectors *vectors;       // pointeur sur un objet Vectors
    int nb_vector;          // nombre de vecteurs
    double *residual;       // residus

    void copy(const Regression&);
    void remove();

    double regression_square_sum_computation() const;
    void residual_computation();
    double residual_mean_computation() const;
    double residual_variance_computation(double residual_mean) const;
    double residual_square_sum_computation() const;

public :

    Regression();
    Regression(int iident , int explanatory_variable , int response_variable , const Vectors &vec);
    Regression(const Regression &regression);
    virtual ~Regression();
    Regression& operator=(const Regression &regression);

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,
                    const char *title = 0) const;

/*    RWDECLARE_COLLECTABLE(Regression);

    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

    // acces membres de la classe

    Vectors* get_vectors() const { return vectors; }
    int get_nb_vector() const { return nb_vector; }
    double get_residual(int ivec) const { return residual[ivec]; }
};



#endif
