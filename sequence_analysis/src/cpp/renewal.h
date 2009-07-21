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



#ifndef RENEWAL_H
#define RENEWAL_H



/****************************************************************
 *
 *  Constantes :
 */


const int DEFAULT_TIME = 20;           // temps d'observation par defaut
const int MAX_TIME = 500;              // temps d'observation maximum
const int PLOT_NEVENT_TIME = 10;       // nombre maximum de lois du temps avant
                                       // le n-eme evenement affichees (sortie Gnuplot)
const int PLOT_NB_TIME = 5;            // nombre maximum de lois du nombre d'evenements
                                       // affichees dans le cas de melange (sortie Gnuplot)

const double RENEWAL_THRESHOLD = 0.99999;  // seuil sur la fonction de repartition
                                           // pour borner une loi
const double RB_THRESHOLD = 2000.;     // seuil pour utiliser le calcul rapide de la loi du
                                       // nombre d'ev correspondant a une loi binomiale
const double RNB_THRESHOLD = 2000.;    // seuil pour utiliser le calcul rapide de la loi du
                                       // nombre d'ev correspondant a une loi binomiale negative

enum {
  INTER_EVENT ,
  WITHIN_OBSERVATION_PERIOD ,
  LENGTH_BIAS ,
  BACKWARD_RECURRENCE_TIME ,
  FORWARD_RECURRENCE_TIME ,
  NB_EVENT ,
  MIXTURE
};

const double MIN_NB_EVENT = 0.4;       // nombre d'evenements moyen minimum
const double MIN_INTER_EVENT = 1.;     // temps moyen minimum entre 2 evenements
const double RENEWAL_INIT_PROBABILITY = 0.001;  // seuil pour l'initialisation de la probabilite
const double RENEWAL_LIKELIHOOD_DIFF = 1.e-5;  // seuil pour stopper les iterations EM
const int RENEWAL_NB_ITER = 10000;     // nombre maximum d'iterations EM
const double RENEWAL_DIFFERENCE_WEIGHT = 0.5;  // poids par defaut de la penalisation
                                               // (cas des differences 1ere ou 2nde)
const double RENEWAL_ENTROPY_WEIGHT = 0.05;  // poids par defaut de la penalisation (cas de l'entropie)
const int RENEWAL_COEFF = 10;          // coefficient arrondi estimateur

const int NB_COMPLETE_INTERVAL = 3;    // nombre minimum d'intervalles de temps complets
const double MEAN_COEFF = 2.;          // coefficient sur la moyenne pour compenser le biais par la longueur
const double MAX_VALUE_COEFF = 10.;    // coefficient pour deduire la valeur maximum de la loi inter-evenement

const int RENEWAL_NB_ELEMENT = 1000000;  // taille maximum de l'echantillon pour la simulation



/****************************************************************
 *
 *  Definition des classes :
 */


class Length_bias : public Parametric {  // loi biaisee par la longueur

/*    friend class Renewal;
    friend class Renewal_data; */

public :

    Length_bias(int inb_value = 0 , int iident = NONPARAMETRIC ,
                int iinf_bound = I_DEFAULT , int isup_bound = I_DEFAULT ,
                double iparameter = D_DEFAULT, double iprobability = D_DEFAULT)
    :Parametric(inb_value , iident , iinf_bound , isup_bound , iparameter , iprobability) {}
    Length_bias(const Parametric &inter_event)
    :Parametric(inter_event) { computation(inter_event); }
    Length_bias(const Length_bias &length_bias)
    :Parametric((Parametric&)length_bias) {}

    void computation(const Parametric&);
};



class Backward : public Parametric {  // loi de l'intervalle de temps apres le dernier evenement

/*    friend class Renewal;
    friend class Renewal_data; */

public :

    Backward(int inb_value = 0 , int iident = NONPARAMETRIC ,
             int iinf_bound = I_DEFAULT , int isup_bound = I_DEFAULT ,
             double iparameter = D_DEFAULT, double iprobability = D_DEFAULT)
    :Parametric(inb_value , iident , iinf_bound , isup_bound , iparameter , iprobability) {}
    Backward(const Backward &dist , int ialloc_nb_value = I_DEFAULT)
    :Parametric(dist , 'c' , ialloc_nb_value) {}

    void computation(const Parametric &inter_event , const Distribution &time);
};



class Nb_event : public Parametric {  // loi du nombre d'evenements

//    friend class Renewal;

// private :
public :

    char type;              // 'o' : ordinaire, 'e' : en equilibre
    int time;               // temps d'observation

/*    RWspace binaryStoreSize(int ialloc_nb_value = I_DEFAULT) const;
    void restoreGuts(RWvistream &is);
    void restoreGuts(RWFile &file);
    void saveGuts(RWvostream &os , int ialloc_nb_value = I_DEFAULT) const;
    void saveGuts(RWFile &file , int ialloc_nb_value = I_DEFAULT) const; */

    void binomial_computation();
    void negative_binomial_computation();

// public :

    Nb_event(char itype = 'v' , int itime = 0 , int inb_value = 0 , int iident = NONPARAMETRIC ,
             int iinf_bound = I_DEFAULT , int isup_bound = I_DEFAULT ,
             double iparameter = D_DEFAULT , double iprobability = D_DEFAULT);
    Nb_event(char itype , int itime , Parametric &inter_event);
    Nb_event(const Nb_event &nb_event , int ialloc_nb_value = I_DEFAULT);

    void ordinary_computation(Parametric &inter_event);
    void computation(Parametric &inter_event);

    // acces membres de la classe

/*    char get_type() const { return type; }
    int get_time() const { return time; } */
};



class Renewal_iterator;
class Time_events;
class Renewal_data;

class Renewal : public STAT_interface {  // processus de renouvellement

    friend class Renewal_iterator;
    friend class Time_events;
    friend class Renewal_data;

    friend Renewal* renewal_building(Format_error &error , const Parametric &inter_event ,
                                     char type, int time);
    friend Renewal* renewal_ascii_read(Format_error& error , const char *path ,
                                       char type, int time,
                                       double cumul_threshold);
    friend std::ostream& operator<<(std::ostream &os , const Renewal &renew)
    { return renew.ascii_write(os , renew.renewal_data , false , false); }

private :

    int nb_iterator;        // nombre d'iterateurs pointant sur l'objet
    Renewal_data *renewal_data;  // pointeur sur un objet Renewal_data
    char type;              // 'o' : ordinaire, 'e' : en equilibre
    int nb_event_max;       // borne max sur le nombre d'evenements
    Distribution *time;     // loi du temps d'observation
    Parametric *inter_event;  // loi inter-evenement
    Length_bias *length_bias;  // loi biaisee par la longueur
    Backward *backward;     // loi de l'intervalle de temps apres le dernier evenement
    Forward *forward;       // loi de l'intervalle de temps residuel
    Parametric **nevent_time;  // lois du temps avant le n-eme evenement
    Nb_event **nb_event;    // lois du nombre d'evenements
                            // pour un temps d'observation donne
    Distribution *mixture;  // melange de lois du nombre d'evenements
    Curves *index_event;    // probabilites de non-evenement/evenement fonction du temps

    void init(int inf_bound , int sup_bound , double parameter , double probability);
    void init(int ident , int inf_bound , int sup_bound , double parameter , double probability);
    void copy(const Renewal &renew , bool data_flag = true);
    void remove();
    void type_init(int itype);

    std::ostream& ascii_write(std::ostream &os , const Renewal_data *timev ,
                              bool exhaustive , bool file_flag) const;
    std::ostream& spreadsheet_write(std::ostream &os , const Renewal_data *timev) const;
    bool plot_write(const char *prefix , const char *title ,
                    const Renewal_data *timev) const;
    MultiPlotSet* get_plotable(const Renewal_data *timev) const;

    void index_event_computation();

    void expectation_step(const Time_events &timev , Reestimation<double> *reestim) const;
    void expectation_step(const Time_events &timev , Reestimation<double> *inter_event_reestim ,
                          Reestimation<double> *length_bias_reestim , int estimator ,
                          bool combination = false , int mean_computation = COMPUTED) const;

public :

    Renewal();
    Renewal(char itype , const Histogram &htime , const Parametric &iinter_event);
    Renewal(char itype , const Distribution &itime , const Parametric &iinter_event);
    Renewal(const Renewal_data &irenewal_data , const Parametric &iinter_event);
    Renewal(const Renewal &renew , bool data_flag = true)
    { copy(renew , data_flag); }
    virtual ~Renewal();
    void conditional_delete();
    Renewal& operator=(const Renewal &renew);

    Parametric_model* extract(Format_error &error , int dist_type , int itime = I_DEFAULT) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,
                    const char *title = 0) const;
    MultiPlotSet* get_plotable() const;

/*    RWDECLARE_COLLECTABLE(Renewal);

    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

    void computation(bool inter_event_flag = true , char itype = 'v' ,
                     const Distribution *dtime = 0);

    double likelihood_computation(const Time_events &timev) const;

    Renewal_data* simulation(Format_error &error , char itype ,
                             const Histogram &ihtime) const;
    Renewal_data* simulation(Format_error &error , char itype ,
                             int nb_element , int itime) const;
    Renewal_data* simulation(Format_error &error , char itype ,
                             int nb_element , const Time_events &itimev) const;

    // acces membres de la classe

    int get_nb_iterator() const { return nb_iterator; }
    Renewal_data* get_renewal_data() const { return renewal_data; }
    char get_type() const { return type; }
    Distribution* get_time() const { return time; }
    Parametric* get_inter_event() const { return inter_event; }
    Length_bias* get_length_bias() const { return length_bias; }
    Backward* get_backward() const { return backward; }
    Forward* get_forward() const { return forward; }
    Parametric* get_nevent_time(int inb_event) const { return nevent_time[inb_event]; }
    Nb_event* get_nb_event(int itime) const { return nb_event[itime]; }
    Distribution* get_mixture() const { return mixture; }
    Curves* get_index_event() const { return index_event; }
};


Renewal* renewal_building(Format_error &error , const Parametric &inter_event ,
                          char type = 'e' , int time = DEFAULT_TIME);
Renewal* renewal_ascii_read(Format_error& error , const char *path ,
                            char type = 'e' , int time = DEFAULT_TIME ,
                            double cumul_threshold = RENEWAL_THRESHOLD);



class Renewal_iterator {  // iterateur processus de renouvellement

private :

    Renewal *renewal;       // pointeur sur un objet Renewal
    int interval;           // intervalle de temps
    int counter;            // compteur
    int length;             // longueur de la sequence
    int *sequence;          // sequence

    void copy(const Renewal_iterator &iterator);

public :

    Renewal_iterator(Renewal *irenewal , int ilength = 1);
    Renewal_iterator(const Renewal_iterator &iterator)
    { copy(iterator); }
    ~Renewal_iterator();
    Renewal_iterator& operator=(const Renewal_iterator &iterator);

    void simulation(int ilength = 1 , char type = 'v');

    // acces membres de la classe

    Renewal* get_renewal() const { return renewal; }
    int get_interval() const { return interval; }
    int get_counter() const { return counter; }
    int get_length() const { return length; }
    int get_sequence(int index) const { return sequence[index]; }
};



class Time_events : public STAT_interface {  // echantillons {temps, nombre d'evenements, effectif}

    friend class Histogram;
    friend class Renewal;

    friend Time_events* time_events_ascii_read(Format_error &error , const char *path);
    friend Time_events* old_time_events_ascii_read(Format_error &error , const char *path);
    friend std::ostream& operator<<(std::ostream &os , const Time_events &timev)
    { return timev.ascii_write(os , true); }

protected :

    int nb_element;         // effectif total
    int nb_class;           // nombre de classes
    int *time;              // temps d'observation
    int *nb_event;          // nombre d'evenements
    int *frequency;         // effectif de chacune des classes
                            // {temps, nombre d'evenements}
    Histogram *htime;       // histogramme du temps d'observation
    Histogram **hnb_event;  // histogrammes du nombre d'evenements
                            // pour un temps d'observation donne
    Histogram *mixture;     // histogramme du nombre d'evenements

    void build_histogram();
    void build_sample();
    void build(int inb_element , int *itime , int *inb_event);
    void copy(const Time_events&);
    void merge(int nb_sample , const Time_events **ptimev);
    void remove();

    std::ostream& ascii_write(std::ostream &os , bool exhaustive , char type) const;
    std::ostream& ascii_file_write(std::ostream &os , bool exhaustive , char type = 'v') const;
    std::ostream& spreadsheet_write(std::ostream &os , char type = 'v') const;

    void nb_element_computation();
    double min_inter_event_computation() const;

public :

    Time_events(int inb_class = 0);
    Time_events(int inb_element , int *itime , int *inb_event)
    { build(inb_element , itime , inb_event); }
    Time_events(int nb_sample , const Time_events **ptimev) { merge(nb_sample , ptimev); }
    Time_events(const Time_events &timev) { copy(timev); }
    virtual ~Time_events();
    Time_events& operator=(const Time_events &timev);

    Distribution_data* extract(Format_error &error , int histo_type , int itime = I_DEFAULT) const;

    Time_events* time_scaling(Format_error &error , int scaling_coeff) const;
    Time_events* time_select(Format_error &error , int min_time ,
                             int max_time) const;
    Time_events* nb_event_select(Format_error &error , int min_nb_event ,
                                 int max_nb_event) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = true) const;
    bool ascii_write(Format_error &error , const char *path ,
                     bool exhaustive = true) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,
                    const char *title = 0) const;
    MultiPlotSet* get_plotable() const;

/*    RWDECLARE_COLLECTABLE(Time_events);

    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

    double information_computation() const;

    Renewal* estimation(Format_error &error , std::ostream &os , char type ,
                        const Parametric &iinter_event , int estimator = LIKELIHOOD ,
                        int nb_iter = I_DEFAULT , int equilibrium_estimator = COMPLETE_LIKELIHOOD ,
                        int mean_computation = COMPUTED , double weight = D_DEFAULT ,
                        int penalty_type = SECOND_DIFFERENCE , int outside = ZERO) const;
    Renewal* estimation(Format_error &error , std::ostream &os , char type , int estimator = LIKELIHOOD ,
                        int nb_iter = I_DEFAULT , int equilibrium_estimator = COMPLETE_LIKELIHOOD ,
                        int mean_computation = COMPUTED , double weight = D_DEFAULT ,
                        int penalty_type = SECOND_DIFFERENCE , int outside = ZERO) const;

    // acces membres de la classe

    int get_nb_element() const { return nb_element; }
    int get_nb_class() const { return nb_class; }
    Histogram* get_htime() const { return htime; }
    Histogram* get_hnb_event(int itime) const { return hnb_event[itime]; }
    Histogram* get_mixture() const { return mixture; }
};


Time_events* time_events_ascii_read(Format_error &error , const char *path);
Time_events* old_time_events_ascii_read(Format_error &error , const char *path);



class Renewal_data : public Time_events {  // donnees correspondant a
                                           // un processus de renouvellement
    friend class Renewal;
    friend class Sequences;

    friend std::ostream& operator<<(std::ostream &os , Renewal_data &timev)
    { return timev.ascii_write(os , false); }

private :

    Renewal *renewal;       // pointeur sur un objet Renewal
    char type;              // 'o' : ordinaire, 'e' : en equilibre
    int *length;            // longueurs des sequences
    int **sequence;         // sequences
    Histogram *inter_event; // intervalles de temps entre 2 evenements
    Histogram *within;      // intervalles de temps entre 2 evenements a l'interieur
                            // de la periode d'observation
    Histogram *length_bias;  // intervalles de temps entre 2 evenements
                             // recouvrant une date d'observation
    Histogram *backward;    // intervalles de temps apres le dernier evenement
    Histogram *forward;     // intervalles de temps residuel
    Curves *index_event;    // probabilites de non-evenement/evenement fonction du temps

    void copy(const Renewal_data &timev , bool model_flag = true);
    void remove();

    std::ostream& ascii_write(std::ostream &os , bool exhaustive , bool file_flag) const;
    std::ostream& spreadsheet_write(std::ostream &os) const;

    void build_index_event(int offset = 1);

public :

    Renewal_data();
    Renewal_data(int nb_element , int itime);
    Renewal_data(int itype , const Renewal &renew);
    Renewal_data(const Time_events &timev , int itype);
    Renewal_data(int nb_sample , const Renewal_data **itimev);
    Renewal_data(const Renewal_data &timev , bool model_flag = true)
    :Time_events(timev) { copy(timev , model_flag); }
    ~Renewal_data();
    Renewal_data& operator=(const Renewal_data&);

    Renewal_data* merge(Format_error &error , int nb_sample , const Renewal_data **itimev) const;
    Distribution_data* extract(Format_error &error , int histo_type , int itime = I_DEFAULT) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,
                    const char *title = 0) const;
    MultiPlotSet* get_plotable() const;

/*    RWDECLARE_COLLECTABLE(Renewal_data);

    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

    Renewal* estimation(Format_error &error , std::ostream &os , const Parametric &iinter_event ,
                        int estimator = LIKELIHOOD , int nb_iter = I_DEFAULT ,
                        int mean_computation = COMPUTED , double weight = D_DEFAULT ,
                        int penalty_type = SECOND_DIFFERENCE , int outside = ZERO) const;
    Renewal* estimation(Format_error &error , std::ostream &os , int estimator = LIKELIHOOD ,
                        int nb_iter = I_DEFAULT , int mean_computation = COMPUTED ,
                        double weight = D_DEFAULT , int penalty_type = SECOND_DIFFERENCE ,
                        int outside = ZERO) const;

    // acces membres de la classe

    Renewal* get_renewal() const { return renewal; }
    char get_type() const { return type; }
    int get_length(int index_seq) const { return length[index_seq]; }
    int get_sequence(int index_seq , int index) const
    { return sequence[index_seq][index]; }
    Histogram* get_inter_event() const { return inter_event; }
    Histogram* get_within() const { return within; }
    Histogram* get_length_bias() const { return length_bias; }
    Histogram* get_backward() const { return backward; }
    Histogram* get_forward() const { return forward; }
    Curves* get_index_event() const { return index_event; }
};



#endif
