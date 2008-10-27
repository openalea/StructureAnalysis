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



#ifndef MIXTURE_H
#define MIXTURE_H



/****************************************************************
 *
 *  Constantes :
 */


const int MIXTURE_NB_COMPONENT = 100;   // nombre maximum de composantes

const double MIXTURE_PARAMETER = 20.;  // parametre initial pour une loi binomiale negative
const double MIN_WEIGHT_STEP = 0.1;    // pas minimum d'initialisation des poids
const double MAX_WEIGHT_STEP = 0.5;    // pas maximum d'initialisation des poids
const int MIXTURE_COEFF = 2;           // coefficient arrondi estimateur
const double MIXTURE_LIKELIHOOD_DIFF = 1.e-5;  // seuil pour stopper les iterations EM
const int MIXTURE_NB_ITER = 500;        // nombre maximum d'iterations EM



/****************************************************************
 *
 *  Definition des classes :
 */


class Histogram;
class Mixture_data;


// melange de lois discretes parametriques

//class Mixture : public STAT_interface , protected Distribution {  
class Mixture : public STAT_interface , public Distribution {  

    friend class Histogram;
    friend class Mixture_data;

    friend Mixture* mixture_building(Format_error &error , int nb_component , double *weight ,
                                     const Parametric **component);
    friend Mixture* mixture_ascii_read(Format_error &error , const char *path ,
                                       double cumul_threshold);
    friend std::ostream& operator<<(std::ostream &os , const Mixture &mixt)
    { return mixt.ascii_write(os , mixt.mixture_data , false , false); }

private :

    Mixture_data *mixture_data;  // pointeur sur un objet Mixture_data
    int nb_component;       // nombre de composantes
    Parametric *weight;     // poids de chaque composante
    Parametric **component; // composantes

    void copy(const Mixture &mixt , bool data_flag = true);
    void remove();

    std::ostream& ascii_write(std::ostream &os , const Mixture_data *mixt_histo ,
                              bool exhaustive , bool file_flag) const;
    std::ostream& spreadsheet_write(std::ostream &os , const Mixture_data *mixt_histo) const;
    bool plot_write(const char *prefix , const char *title ,
                    const Mixture_data *mixt_histo) const;
    plotable::MultiPlotSet* get_plotable(const Mixture_data *mixt_histo) const;

    int nb_parameter_computation() const;
    double penalty_computation() const;

    void init(const Histogram &histo , bool *estimate , int min_inf_bound ,
              bool component_flag);
    void expectation_step(Mixture_data *mixt_histo , int nb_element) const;
    void variance_correction(Mixture_data *mixt_histo , bool *estimate ,
                             int min_inf_bound) const;
    bool component_order_test() const;

public :

    Mixture();
    Mixture(int inb_component , double *pweight , const Parametric **pcomponent);
    Mixture(const Mixture &mixt , bool *component_flag , int inb_value);
    Mixture(int inb_component , const Parametric **pcomponent);
    Mixture(const Mixture &mixt , bool data_flag = true)
    :Distribution(mixt) { copy(mixt , data_flag); }
    ~Mixture();
    Mixture& operator=(const Mixture &mixt);

    Parametric_model* extract(Format_error &error , int index) const;
    Mixture_data* extract_data(Format_error &error) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,
                    const char *title = 0) const;
    plotable::MultiPlotSet* get_plotable() const;

/*    RWDECLARE_COLLECTABLE(Mixture);

    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

    void computation(int min_nb_value = 1 , double cumul_threshold = CUMUL_THRESHOLD ,
                     bool component_flag = true);
    double likelihood_computation(const Mixture_data &mixt_histo) const;
    Mixture_data* simulation(Format_error &error , int nb_element) const;
  
    // acces membres de la classe

    Mixture_data* get_mixture_data() const { return mixture_data; }
    int get_nb_component() const { return nb_component; }
    Parametric* get_weight() const { return weight; }
    Parametric* get_component(int index) const { return component[index]; }
};


Mixture* mixture_building(Format_error &error , int nb_component , double *weight ,
                          const Parametric **component);
Mixture* mixture_ascii_read(Format_error &error , const char *path ,
                            double cumul_threshold = CUMUL_THRESHOLD);



// structure de donnees correspondant
//class Mixture_data : public STAT_interface , protected Histogram {  

class Mixture_data : public STAT_interface , public Histogram {  
                                                                    // a un melange
    friend class Histogram;
    friend class Mixture;

    friend std::ostream& operator<<(std::ostream &os , const Mixture_data &mixt_histo)
    { return mixt_histo.ascii_write(os , false); }

private :

    Mixture *mixture;       // pointeur sur un objet Mixture
    int nb_component;       // nombre de composantes
    Histogram *weight;      // histogramme des poids
    Histogram **component;  // histogrammes correspondant aux composantes

    void copy(const Mixture_data &mixt_histo , bool model_flag = true);
    void remove();

public :

    Mixture_data();
    Mixture_data(const Histogram &histo , int inb_component);
    Mixture_data(const Mixture &mixt);
    Mixture_data(const Mixture_data &mixt_histo , bool model_flag = true)
      :Histogram(mixt_histo) { copy(mixt_histo , model_flag); }
    virtual ~Mixture_data();
    Mixture_data& operator=(const Mixture_data &mixt_histo);

    Distribution_data* extract(Format_error &error , int index) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,
                    const char *title = 0) const;
    plotable::MultiPlotSet* get_plotable() const;

/*    RWDECLARE_COLLECTABLE(Mixture_data);

    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

    double information_computation() const;

    // acces membres de la classe

    Mixture* get_mixture() const { return mixture; }
    int get_nb_component() const { return nb_component; }
    Histogram* get_weight() const { return weight; }
    Histogram* get_component(int index) const { return component[index]; }
};



#endif
