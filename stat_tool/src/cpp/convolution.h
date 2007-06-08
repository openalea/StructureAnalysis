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



#ifndef CONVOLUTION_H
#define CONVOLUTION_H



/****************************************************************
 *
 *  Constantes :
 */


const int CONVOLUTION_NB_DISTRIBUTION = 10;  // nombre maximum de lois elementaires
const double CONVOLUTION_THRESHOLD = 0.9999;  // seuil sur la fonction de repartition
                                              // pour borner une loi

const double CONVOLUTION_INIT_PROBABILITY = 0.001;  // seuil pour l'initialisation de la probabilite
const double CONVOLUTION_LIKELIHOOD_DIFF = 1.e-5;  // seuil pour stopper les iterations EM
const int CONVOLUTION_NB_ITER = 10000;  // nombre maximum d'iterations EM
const double CONVOLUTION_DIFFERENCE_WEIGHT = 0.5;  // poids par defaut de la penalisation
                                                   // (cas des differences 1ere ou 2nde)
const double CONVOLUTION_ENTROPY_WEIGHT = 0.1;  // poids par defaut de la penalisation (cas de l'entropie)
const int CONVOLUTION_COEFF = 10;      // coefficient arrondi estimateur



/****************************************************************
 *
 *  Definition des classes :
 */


class Histogram;
class Convolution_data;

// class Convolution : public STAT_interface , public Distribution {
class Convolution : public STAT_interface , protected Distribution {  // produit de convolution
                                                                      // de lois discretes
    friend class Histogram;
    friend class Convolution_data;

    friend Convolution* convolution_building(Format_error &error , int nb_dist ,
                                             const Parametric **dist);
    friend Convolution* convolution_ascii_read(Format_error &error , const char *path ,
                                               double cumul_threshold);
    friend std::ostream& operator<<(std::ostream &os , const Convolution &convol)
    { return convol.ascii_write(os , convol.convolution_data , false , false); }

private :

    Convolution_data *convolution_data;  // pointeur sur un objet Convolution_data
    int nb_distribution;    // nombre de lois elementaires
    Parametric **distribution;  // lois elementaires

    void copy(const Convolution &convol , bool data_flag = true);
    void remove();

    std::ostream& ascii_write(std::ostream &os , const Convolution_data *convol_histo ,
                              bool exhaustive , bool file_flag) const;
    std::ostream& spreadsheet_write(std::ostream &os , const Convolution_data *convol_histo) const;
    bool plot_write(const char *prefix , const char *title ,
                    const Convolution_data *convol_histo) const;

    void expectation_step(const Histogram &histo , const Distribution **partial_convol ,
                          Reestimation<double> **reestim) const;

public :

    Convolution();
    Convolution(int nb_dist , const Parametric **pdist);
    Convolution(const Parametric &known_dist , const Parametric &unknown_dist);
    Convolution(const Convolution &convol , bool data_flag = true)
    :Distribution(convol) { copy(convol , data_flag); }
    ~Convolution();
    Convolution& operator=(const Convolution &convol);

    Parametric_model* extract(Format_error &error , int index) const;
    Convolution_data* extract_data(Format_error &error) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,
                    const char *title = 0) const;

/*    RWDECLARE_COLLECTABLE(Convolution);

    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

    void computation(int min_nb_value = 1 , double cumul_threshold = CONVOLUTION_THRESHOLD ,
                     bool *dist_flag = 0);
    Convolution_data* simulation(Format_error &error , int nb_element) const;

    // acces membres de la classe

    Convolution_data* get_convolution_data() const { return convolution_data; }
    int get_nb_distribution() const { return nb_distribution; }
    Parametric* get_distribution(int index) const { return distribution[index]; }
};


Convolution* convolution_building(Format_error &error , int nb_dist ,
                                  const Parametric **dist);
Convolution* convolution_ascii_read(Format_error &error , const char *path ,
                                    double cumul_threshold = CONVOLUTION_THRESHOLD);



// class Convolution_data : public STAT_interface , public Histogram {
class Convolution_data : public STAT_interface , protected Histogram {  // structure de donnees correspondant
                                                                        // a un produit de convolution
    friend class Histogram;
    friend class Convolution;

    friend std::ostream& operator<<(std::ostream &os , const Convolution_data &convol_histo)
    { return convol_histo.ascii_write(os , false); }

private :

    Convolution *convolution;  // pointeur sur un objet Convolution
    int nb_histogram;       // nombre d'histogrammes
    Histogram **histogram;  // histogrammes correspondant aux lois elementaires

    void copy(const Convolution_data &convol_histo , bool model_flag = true);
    void remove();

public :

    Convolution_data();
    Convolution_data(const Histogram &histo , int nb_histo);
    Convolution_data(const Convolution &convol);
    Convolution_data(const Convolution_data &convol_histo , bool model_flag = true)
    :Histogram(convol_histo) { copy(convol_histo , model_flag); }
    virtual ~Convolution_data();
    Convolution_data& operator=(const Convolution_data &convol_histo);

    Distribution_data* extract(Format_error &error , int index) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,
                    const char *title = 0) const;

/*    RWDECLARE_COLLECTABLE(Convolution_data);

    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

    // acces membres de la classe

    Convolution* get_convolution() const { return convolution; }
    int get_nb_histogram() const { return nb_histogram; }
    Histogram* get_histogram(int index) const { return histogram[index]; }
};



#endif
