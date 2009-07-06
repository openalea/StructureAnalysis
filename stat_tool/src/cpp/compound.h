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



#ifndef COMPOUND_H
#define COMPOUND_H



/****************************************************************
 *
 *  Constantes :
 */


const double COMPOUND_THRESHOLD = 0.99999;  // seuil sur la fonction de repartition
                                            // pour borner une loi

const double COMPOUND_INIT_PROBABILITY = 0.001;  // seuil pour l'initialisation de la probabilite
const double COMPOUND_LIKELIHOOD_DIFF = 1.e-5;  // seuil pour stopper les iterations EM
const int COMPOUND_NB_ITER = 10000;    // nombre maximum d'iterations EM
const double COMPOUND_DIFFERENCE_WEIGHT = 0.5;  // poids par defaut de la penalisation
                                                // (cas des differences 1ere ou 2nde)
const double COMPOUND_ENTROPY_WEIGHT = 0.1;  // poids par defaut de la penalisation (cas de l'entropie)
const int COMPOUND_COEFF = 10;         // coefficient arrondi estimateur



/****************************************************************
 *
 *  Definition des classes :
 */


class Histogram;
class Compound_data;

class Compound : public STAT_interface , public Distribution {
  //class Compound : public STAT_interface , protected Distribution {  // loi composee

    friend class Histogram;
    friend class Compound_data;

    friend Compound* compound_ascii_read(Format_error &error , const char *path ,
                                         double cumul_threshold);
    friend std::ostream& operator<<(std::ostream &os , const Compound &compound)
    { return compound.ascii_write(os , compound.compound_data , false , false); }

private :

    Compound_data *compound_data;  // pointeur sur un objet Compound_data
    Parametric *sum_distribution;  // loi de la somme
    Parametric *distribution;  // loi elementaire

    void copy(const Compound &compound , bool data_flag = true);

    std::ostream& ascii_write(std::ostream &os , const Compound_data *compound_histo ,
                              bool exhaustive , bool file_flag) const;
    std::ostream& spreadsheet_write(std::ostream &os , const Compound_data *compound_histo) const;
    bool plot_write(const char *prefix , const char *title ,
                    const Compound_data *compound_histo) const;
    MultiPlotSet* get_plotable(const Compound_data *compound_histo) const;

    void computation(Parametric **power_dist , int min_nb_value ,
                     double cumul_threshold , bool sum_flag , bool dist_flag);
    void expectation_step(const Histogram &histo , Parametric **power_dist ,
                          Reestimation<double> *sum_reestim , Reestimation<double> *reestim) const;

public :

    Compound();
    Compound(const Parametric &sum_dist , const Parametric &dist ,
             double cumul_threshold = COMPOUND_THRESHOLD);
    Compound(const Parametric &sum_dist , const Parametric &dist , char type);
    Compound(const Compound &compound , bool data_flag = true)
    :Distribution(compound) { copy(compound , data_flag); }
    ~Compound();
    Compound& operator=(const Compound &compound);

    Compound_data* extract_data(Format_error &error) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,
                    const char *title = 0) const;
    MultiPlotSet* get_plotable() const;

/*    RWDECLARE_COLLECTABLE(Compound);

    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

    void computation(int min_nb_value = 1 ,
                     double cumul_threshold = COMPOUND_THRESHOLD ,
                     bool sum_flag = true , bool dist_flag = true);
    Compound_data* simulation(Format_error &error , int nb_element) const;

    // acces membres de la classe

    Compound_data* get_compound_data() const { return compound_data; }
    Parametric* get_sum_distribution() const { return sum_distribution; }
    Parametric* get_distribution() const { return distribution; }
};


Compound* compound_ascii_read(Format_error &error , const char *path ,
                              double cumul_threshold = COMPOUND_THRESHOLD);



class Compound_data : public STAT_interface , public Histogram {
  //class Compound_data : public STAT_interface , protected Histogram {  // structure de donnees correspondant
                                                                     // a une loi composee
    friend class Histogram;
    friend class Compound;

    friend std::ostream& operator<<(std::ostream &os , const Compound_data &compound_histo)
    { return compound_histo.ascii_write(os , false); }

private :

    Compound *compound;     // pointeur sur un objet Compound
    Histogram *sum_histogram;   // histogramme correspondant a la loi de la somme
    Histogram *histogram;   // histogramme correspondant a la loi elementaire

    void copy(const Compound_data &compound_histo , bool model_flag = true);

public :

    Compound_data();
    Compound_data(const Histogram &histo , const Compound &icompound);
    Compound_data(const Compound &icompound);
    Compound_data(const Compound_data &compound_histo , bool model_flag = true)
    :Histogram(compound_histo) { copy(compound_histo , model_flag); }
    virtual ~Compound_data();
    Compound_data& operator=(const Compound_data &compound_histo);

    Distribution_data* extract(Format_error &error , char type) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,
                    const char *title = 0) const;
    MultiPlotSet* get_plotable() const;

/*    RWDECLARE_COLLECTABLE(Compound_data);

    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

    // acces membres de la classe

    Compound* get_compound() const { return compound; }
    Histogram* get_sum_histogram() const { return sum_histogram; }
    Histogram* get_histogram() const { return histogram; }
};



#endif
