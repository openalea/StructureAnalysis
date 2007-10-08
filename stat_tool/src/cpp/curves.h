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



#ifndef CURVES_H
#define CURVES_H



/****************************************************************
 *
 *  Constantes :
 */


const int MAX_FREQUENCY = 50;          // effectif maximum pour lisser les courbes
const int MAX_RANGE = 2;               // demi-largeur maximum de la fenetre de lissage

const int PLOT_NB_CURVE = 12;          // nombre maximum des courbes (sortie Gnuplot)
const int PLOT_MIN_FREQUENCY = 10;     // effectif minimum pour afficher les points
                                       // des courbes (sortie Gnuplot)




/****************************************************************
 *
 *  Definition des classes :
 */


class Curves {          // famille de courbes avec effectif

/*    friend class Distribution;
    friend class Histogram;
    friend class Function;
    friend class Renewal;
    friend class Renewal_data;
    friend class Nonparametric_sequence_process;
    friend class Markov;
    friend class Variable_order_markov;
    friend class Semi_markov;
    friend class Sequences;
    friend class Sequence_characteristics;
    friend class Markovian_sequences; */

    friend std::ostream& operator<<(std::ostream& , const Curves&);

// protected :
public :

    int nb_curve;           // nombre de courbes
    int length;             // longueur des courbes
    int offset;             // abscisse des premiers points definis
    int *frequency;         // effectifs correspondant a chaque abscisse
    double **point;         // points des courbes

    void copy(const Curves&);
    void smooth(const Curves &curves , int max_frequency);
    void remove();

    std::ostream& ascii_print(std::ostream &os , bool file_flag = false ,
                              const Curves *curves = 0) const;
    std::ostream& spreadsheet_print(std::ostream &os , const Curves *curves = 0) const;
    int plot_length_computation() const;
    bool plot_print(const char *path , int ilength = I_DEFAULT ,
                    const Curves *curves_0 = 0 , const Curves *curves_1 = 0) const;
    bool plot_print_standard_residual(const char *path , double *standard_residual = 0) const;  // sequence_analysis

/*    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

    double mean_computation(int index) const;
    double total_square_sum_computation(int index , double mean) const;

// public :

    Curves();
    Curves(int inb_curve , int ilength , bool frequency_flag = false , bool init_flag = true);
    Curves(const Curves &curves , char transform = 'c' , int max_frequency = MAX_FREQUENCY);
    Curves(const Distribution &dist);
    Curves(const Histogram &histo);
    ~Curves();
    Curves& operator=(const Curves &curves);

    int max_frequency_computation() const;
    int nb_element_computation() const;

    // acces membres de la classe

/*    int get_nb_curve() const { return nb_curve; }
    int get_length() const { return length; }
    int get_offset() const { return offset; }
    int get_frequency(int ilength) const { return frequency[ilength]; }
    double get_point(int index_curve , int ilength) const
    { return point[index_curve][ilength]; } */
};



#endif
