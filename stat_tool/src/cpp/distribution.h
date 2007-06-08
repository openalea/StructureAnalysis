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



#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H



// class Parametric_model : public STAT_interface , public Parametric {
class Parametric_model : public STAT_interface , protected Parametric {  // loi de probabilite parametrique

    friend class Distribution;  // Hack pour Windows

    friend class Histogram;
    friend class Distribution_data;

    friend Parametric_model* parametric_ascii_read(Format_error &error , const char *path ,
                                                   double cumul_threshold);
    friend std::ostream& operator<<(std::ostream &os , const Parametric_model &dist)
    { return dist.ascii_write(os , dist.histogram , false , false); }

private :

    Distribution_data *histogram;  // pointeur sur un objet Distribution_data

    std::ostream& ascii_write(std::ostream &os , const Distribution_data *histo ,
                              bool exhaustive , bool file_flag) const;
    std::ostream& spreadsheet_write(std::ostream &os , const Distribution_data *histo) const;
    bool plot_write(Format_error &error , const char *prefix , const char *title ,
                    const Distribution_data *histo) const;

public :

    Parametric_model(int inb_value = 0 , int iident = NONPARAMETRIC ,
                     int iinf_bound = I_DEFAULT , int isup_bound = I_DEFAULT ,
                     double iparameter = D_DEFAULT, double iprobability = D_DEFAULT)
    :Parametric(inb_value , iident , iinf_bound , isup_bound , iparameter , iprobability)
    { histogram = 0; }
    Parametric_model(int iident , int iinf_bound , int isup_bound , double iparameter ,
                     double iprobability , double cumul_threshold = CUMUL_THRESHOLD)
    :Parametric(iident , iinf_bound , isup_bound , iparameter , iprobability , cumul_threshold)
    { histogram = 0; }
    Parametric_model(const Histogram &histo);
    Parametric_model(const Distribution &dist)
    :Parametric(dist) { histogram = 0; }
    Parametric_model(const Parametric &dist)
    :Parametric(dist) { histogram = 0; }
    Parametric_model(const Distribution &dist , const Histogram *histo);
    Parametric_model(const Parametric &dist , const Histogram *histo);
    Parametric_model(const Parametric_model &dist , bool data_flag = true);
    ~Parametric_model();
    Parametric_model& operator=(const Parametric_model &dist);

    Distribution_data* extract_data(Format_error &error) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,
                    const char *title = 0) const;

/*    RWDECLARE_COLLECTABLE(Parametric_model);

    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

    Distribution_data* simulation(Format_error &error , int nb_element) const;

    Distribution_data* get_histogram() const { return histogram; }
};


Parametric_model* parametric_ascii_read(Format_error &error , const char *path ,
                                        double cumul_threshold = CUMUL_THRESHOLD);



// class Distribution_data : public STAT_interface , protected Histogram {
class Distribution_data : public STAT_interface , public Histogram {  // histogramme

    friend class Parametric_model;

    friend std::ostream& operator<<(std::ostream &os , const Distribution_data &histo)
    { return histo.ascii_write(os); }

private :

    Parametric_model *distribution;  // pointeur sur un objet Parametric_model

    std::ostream& ascii_write(std::ostream &os , bool exhaustive , bool file_flag) const;

public :

    Distribution_data(int inb_value = 0)
    :Histogram(inb_value) { distribution = 0; }
    Distribution_data(const Distribution &dist)
    :Histogram(dist) { distribution = 0; }
    Distribution_data(const Histogram &histo)
    :Histogram(histo) { distribution = 0; }
    Distribution_data(int inb_element , int *pelement)
    :Histogram(inb_element , pelement) { distribution = 0; }
    Distribution_data(const Histogram &histo , char transform , int param)
    :Histogram(histo , transform , param) { distribution = 0; }
    Distribution_data(int nb_histo , const Histogram **phisto)
    :Histogram(nb_histo , phisto) { distribution = 0; }
    Distribution_data(const Histogram &histo , const Distribution *dist);
    Distribution_data(const Histogram &histo , const Parametric *dist);
    Distribution_data(const Distribution_data &histo , bool model_flag = true);
    virtual ~Distribution_data();
    Distribution_data& operator=(const Distribution_data &histo);

    Parametric_model* extract_model(Format_error &error) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,
                    const char *title = 0) const;

/*    RWDECLARE_COLLECTABLE(Distribution_data);

    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

    Parametric* get_distribution() const { return distribution; }
};



#endif
