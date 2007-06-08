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



#include <sstream>
#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
// #include <rw/vstream.h>
// #include <rw/rwfile.h>
#include "stat_tools.h"
#include "distribution.h"
#include "compound.h"
#include "stat_label.h"

using namespace std;


extern char* label(const char *file_name);



/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe Compound.
 *
 *--------------------------------------------------------------*/

Compound::Compound()

{
  compound_data = 0;
  sum_distribution = 0;
  distribution = 0;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Compound.
 *
 *  arguments : references sur la loi de la somme et sur la loi elementaire.
 *              seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

Compound::Compound(const Parametric &sum_dist , const Parametric &dist ,
                   double cumul_threshold)

{
  compound_data = 0;

  sum_distribution = new Parametric(sum_dist , 'n');
  distribution = new Parametric(dist , 'n');

  Distribution::init((sum_distribution->nb_value - 1) * (distribution->nb_value - 1) + 1);

  computation(1 , cumul_threshold , false , false);
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Compound.
 *
 *  arguments : references sur la loi de la somme et sur la loi elementaire.
 *              type de la loi inconnue.
 *
 *--------------------------------------------------------------*/

Compound::Compound(const Parametric &sum_dist , const Parametric &dist , char type)

{
  compound_data = 0;

  switch (type) {

  case 's' : {
    sum_distribution = new Parametric(sum_dist , 'c' , (int)(sum_dist.nb_value * NB_VALUE_COEFF));
    if ((dist.ident == POISSON) || (dist.ident == NEGATIVE_BINOMIAL)) {
      distribution = new Parametric(dist.ident , dist.inf_bound , dist.sup_bound ,
                                    dist.parameter , dist.probability , COMPOUND_THRESHOLD);
    }
    else {
      distribution = new Parametric(dist , 'n');
    }
    break;
  }

  case 'e' : {
    sum_distribution = new Parametric(sum_dist , 'n');
    distribution = new Parametric(dist , 'c' , (int)(dist.nb_value * NB_VALUE_COEFF));
    break;
  }
  }

  Distribution::init((sum_distribution->alloc_nb_value - 1) * (distribution->alloc_nb_value - 1) + 1);
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Compound.
 *
 *  arguments : reference sur un objet Compound,
 *              flag copie de l'objet Compound_data.
 *
 *--------------------------------------------------------------*/

void Compound::copy(const Compound &compound , bool data_flag)

{
  if ((data_flag) && (compound.compound_data)) {
    compound_data = new Compound_data(*(compound.compound_data) , false);
  }
  else {
    compound_data = 0;
  }

  sum_distribution = new Parametric(*(compound.sum_distribution));
  distribution = new Parametric(*(compound.distribution));
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Compound.
 *
 *--------------------------------------------------------------*/

Compound::~Compound()

{
  delete compound_data;

  delete sum_distribution;
  delete distribution;
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Compound.
 *
 *  argument : reference sur un objet Compound.
 *
 *--------------------------------------------------------------*/

Compound& Compound::operator=(const Compound &compound)

{
  if (&compound != this) {
    delete compound_data;

    delete sum_distribution;
    delete distribution;

    delete [] mass;
    delete [] cumul;

    Distribution::copy(compound);
    copy(compound);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Extraction de la partie "donnees" d'un objet Compound.
 *
 *  argument : reference sur un objet Format_error.
 *
 *--------------------------------------------------------------*/

Compound_data* Compound::extract_data(Format_error &error) const

{
  Compound_data *compound_histo;


  error.init();

  if (!compound_data) {
    compound_histo = 0;
    error.update(STAT_error[STATR_NO_DATA]);
  }

  else {
    compound_histo = new Compound_data(*compound_data);
    compound_histo->compound = new Compound(*this , false);
  }

  return compound_histo;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Compound a partir d'un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

Compound* compound_ascii_read(Format_error &error , const char *path ,
                              double cumul_threshold)

{
  RWCString buffer , token;

  size_t position;
  bool status;
  register int i;
  int line , read_line;
  Parametric *sum_dist , *dist;
  Compound *compound;
  ifstream in_file(path);


  compound = 0;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    line = 0;
    read_line = 0;

    sum_dist = 0;
    dist = 0;

    while (buffer.readLine(in_file , false)) {
      line++;

#     ifdef DEBUG
      cout << line << "  " << buffer << endl;
#     endif

      position = buffer.first('#');
      if (position != RW_NPOS) {
        buffer.remove(position);
      }
      i = 0;

      RWCTokenizer next(buffer);

      while (!((token = next()).isNull())) {

        // test mots cles COMPOUND_DISTRIBUTION / SUM_DISTRIBUTION /
        // ELEMENTARY_DISTRIBUTION

        if (i == 0) {
          switch (read_line) {

          case 0 : {
            if (token != STAT_word[STATW_COMPOUND]) {
              status = false;
              error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_COMPOUND] , line);
            }
            break;
          }

          case 1 : {
            if (token != STAT_word[STATW_SUM]) {
              status = false;
              error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_SUM] , line);
            }
            break;
          }

          case 2 : {
            if (token != STAT_word[STATW_ELEMENTARY]) {
              status = false;
              error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_ELEMENTARY] , line);
            }
            break;
          }
          }
        }

        i++;
      }

      if (i > 0) {
        if (i != 1) {
          status = false;
          error.update(STAT_parsing[STATP_FORMAT] , line);
        }

        switch (read_line) {

        case 1 : {
          sum_dist = parametric_parsing(error , in_file , line ,
                                        NEGATIVE_BINOMIAL , CUMUL_THRESHOLD);
          if (!sum_dist) {
            status = false;
          }
          break;
        }

        case 2 : {
          dist = parametric_parsing(error , in_file , line ,
                                    NEGATIVE_BINOMIAL , cumul_threshold);
          if (!dist) {
            status = false;
          }
          break;
        }
        }

        read_line++;
        if (read_line == 3) {
          break;
        }
      }
    }

    if (read_line != 3) {
      status = false;
      error.update(STAT_parsing[STATP_FORMAT] , line);
    }

    while (buffer.readLine(in_file , false)) {
      line++;

#     ifdef DEBUG
      cout << line << "  " << buffer << endl;
#     endif

      position = buffer.first('#');
      if (position != RW_NPOS) {
        buffer.remove(position);
      }

      if (!(buffer.isNull())) {
        status = false;
        error.update(STAT_parsing[STATP_FORMAT] , line);
      }
    }

    if (status) {
      compound = new Compound(*sum_dist , *dist , cumul_threshold);
    }

    delete sum_dist;
    delete dist;
  }

  return compound;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet Compound.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& Compound::line_write(ostream &os) const

{
  os << STAT_label[STATL_MEAN] << ": " << mean << "   "
     << STAT_label[STATL_VARIANCE] << ": " << variance;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'une loi composee et de la structure de donnees associee
 *  dans un fichier.
 *
 *  arguments : stream, pointeur sur un objet Compound_data,
 *              flag niveau de detail, flag fichier.
 *
 *--------------------------------------------------------------*/

ostream& Compound::ascii_write(ostream &os , const Compound_data *compound_histo ,
                               bool exhaustive , bool file_flag) const

{
  os << STAT_word[STATW_COMPOUND] << endl;
  ascii_characteristic_print(os , exhaustive , file_flag);

  if (compound_histo) {
    double likelihood , information;
    Test test(CHI2);


    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_HISTOGRAM] << " - ";
    compound_histo->ascii_characteristic_print(os , exhaustive , file_flag);

    likelihood = likelihood_computation(*compound_histo);
    information = compound_histo->information_computation();

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   ("
       << STAT_label[STATL_NORMALIZED] << ": " << likelihood / compound_histo->nb_element << ")" << endl;

    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_MAX_LIKELIHOOD] << ": " << information << "   ("
       << STAT_label[STATL_INFORMATION] << ": " << information / compound_histo->nb_element << endl;

    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_DEVIANCE] << ": " << 2 * (information - likelihood) << endl;

    chi2_fit(*compound_histo , test);
    os << "\n";
    test.ascii_print(os , file_flag);
  }

  if (exhaustive) {
    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << "  ";
    if (compound_histo) {
      os << " | " << STAT_label[STATL_HISTOGRAM];
    }
    os << " | " << STAT_label[STATL_COMPOUND] << " " << STAT_label[STATL_DISTRIBUTION];
    if (compound_histo) {
      os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM]
         << " " << STAT_label[STATL_FUNCTION];
    }
    os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_COMPOUND]
       << " " << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION] << endl;
    ascii_print(os , file_flag , true , false , compound_histo);
  }

  os << "\n" << STAT_word[STATW_SUM] << endl;
  sum_distribution->ascii_print(os);
  sum_distribution->ascii_characteristic_print(os , exhaustive , file_flag);

  if (compound_histo) {
    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_SUM] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
    compound_histo->sum_histogram->ascii_characteristic_print(os , exhaustive , file_flag);
  }

  if (exhaustive) {
    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << "  ";
    if (compound_histo) {
      os << " | " << STAT_label[STATL_SUM] << " " << STAT_label[STATL_HISTOGRAM];
    }
    os << " | " << STAT_label[STATL_SUM] << " " << STAT_label[STATL_DISTRIBUTION];
    if (compound_histo) {
      os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_SUM]
         << " " << STAT_label[STATL_HISTOGRAM] << " " << STAT_label[STATL_FUNCTION];
    }
    os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_SUM]
       << " " << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION] << endl;

    sum_distribution->Distribution::ascii_print(os , file_flag , true , false ,
                                                (compound_histo ? compound_histo->sum_histogram : 0));
  }

  os << "\n" << STAT_word[STATW_ELEMENTARY] << endl;
  distribution->ascii_print(os);
  distribution->ascii_characteristic_print(os , exhaustive , file_flag);

  if (compound_histo) {
    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_ELEMENTARY] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
    compound_histo->histogram->ascii_characteristic_print(os , exhaustive , file_flag);
  }

  if (exhaustive) {
    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << "  ";
    if ((compound_histo) && (compound_histo->histogram->nb_element > 0)) {
      os << " | " << STAT_label[STATL_ELEMENTARY] << " " << STAT_label[STATL_HISTOGRAM];
    }
    os << " | " << STAT_label[STATL_ELEMENTARY] << " " << STAT_label[STATL_DISTRIBUTION];
    if ((compound_histo) && (compound_histo->histogram->nb_element > 0)) {
      os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_ELEMENTARY]
         << " " << STAT_label[STATL_HISTOGRAM] << " " << STAT_label[STATL_FUNCTION];
    }
    os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_ELEMENTARY]
       << " " << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION] << endl;

    distribution->Distribution::ascii_print(os , file_flag , true , false ,
                                            (((compound_histo) && (compound_histo->histogram->nb_element > 0)) ?
                                             compound_histo->histogram : 0));
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Compound.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Compound::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , compound_data , exhaustive , false);
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Compound dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Compound::ascii_write(Format_error &error , const char *path ,
                           bool exhaustive) const

{
  bool status;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    ascii_write(out_file , compound_data , exhaustive , true);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'une loi composee  et de la structure de donnees associee
 *  dans un fichier au format tableur.
 *
 *  arguments : stream, pointeur sur une famille d'histogrammes.
 *
 *--------------------------------------------------------------*/

ostream& Compound::spreadsheet_write(ostream &os , const Compound_data *compound_histo) const

{
  os << STAT_word[STATW_COMPOUND] << endl;
  spreadsheet_characteristic_print(os , true);

  if (compound_histo) {
    double likelihood , information;
    Test test(CHI2);


    os << "\n" << STAT_label[STATL_HISTOGRAM] << "\t";
    compound_histo->spreadsheet_characteristic_print(os , true);

    likelihood = likelihood_computation(*compound_histo);
    information = compound_histo->information_computation();

    os << "\n" << STAT_label[STATL_LIKELIHOOD] << "\t" << likelihood << "\t"
       << STAT_label[STATL_NORMALIZED] << "\t" << likelihood / compound_histo->nb_element << endl;
    os << STAT_label[STATL_MAX_LIKELIHOOD] << "\t" << information << "\t"
       << STAT_label[STATL_INFORMATION] << "\t" << information / compound_histo->nb_element << endl;
    os << STAT_label[STATL_DEVIANCE] << "\t" << 2 * (information - likelihood) << endl;

    chi2_fit(*compound_histo , test);
    os << "\n";
    test.spreadsheet_print(os);
  }

  os << "\n";
  if (compound_histo) {
    os << "\t" << STAT_label[STATL_HISTOGRAM];
  }
  os << "\t" << STAT_label[STATL_COMPOUND] << " " << STAT_label[STATL_DISTRIBUTION];
  if (compound_histo) {
    os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM]
       << " " << STAT_label[STATL_FUNCTION];
  }
  os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_COMPOUND]
     << " " << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION] << endl;

  spreadsheet_print(os , true , false , false , compound_histo);

  os << "\n" << STAT_word[STATW_SUM] << endl;
  sum_distribution->spreadsheet_print(os);
  sum_distribution->spreadsheet_characteristic_print(os , true);

  if (compound_histo) {
    os << "\n" << STAT_label[STATL_SUM] << " " << STAT_label[STATL_HISTOGRAM] << "\t";
    compound_histo->sum_histogram->spreadsheet_characteristic_print(os , true);
  }

  os << "\n";
  if (compound_histo) {
    os << "\t" << STAT_label[STATL_SUM] << " " << STAT_label[STATL_HISTOGRAM];
  }
  os << "\t" << STAT_label[STATL_SUM] << " " << STAT_label[STATL_DISTRIBUTION];
  if (compound_histo) {
    os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_SUM]
       << " " << STAT_label[STATL_HISTOGRAM] << " " << STAT_label[STATL_FUNCTION];
  }
  os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_SUM]
     << " " << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION] << endl;

  sum_distribution->Distribution::spreadsheet_print(os , true , false , false ,
                                                    (compound_histo ? compound_histo->sum_histogram : 0));

  os << "\n" << STAT_word[STATW_ELEMENTARY] << endl;
  distribution->spreadsheet_print(os);
  distribution->spreadsheet_characteristic_print(os , true);

  if (compound_histo) {
    os << "\n" << STAT_label[STATL_ELEMENTARY] << " " << STAT_label[STATL_HISTOGRAM] << "\t";
    compound_histo->histogram->spreadsheet_characteristic_print(os , true);
  }

  os << "\n";
  if ((compound_histo) && (compound_histo->histogram->nb_element > 0)) {
    os << "\t" << STAT_label[STATL_ELEMENTARY] << " " << STAT_label[STATL_HISTOGRAM];
  }
  os << "\t" << STAT_label[STATL_ELEMENTARY] << " " << STAT_label[STATL_DISTRIBUTION];
  if ((compound_histo) && (compound_histo->histogram->nb_element > 0)) {
    os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_ELEMENTARY]
       << " " << STAT_label[STATL_HISTOGRAM] << " " << STAT_label[STATL_FUNCTION];
  }
  os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_ELEMENTARY]
     << " " << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION] << endl;

  distribution->Distribution::spreadsheet_print(os , true , false , false ,
                                                (((compound_histo) && (compound_histo->histogram->nb_element > 0)) ?
                                                 compound_histo->histogram : 0));

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Compound dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

bool Compound::spreadsheet_write(Format_error &error , const char *path) const

{
  bool status;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    spreadsheet_write(out_file , compound_data);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'une loi composee et de la structure de donnees associee.
 *
 *  arguments : prefixe des fichiers, titre des figures,
 *              pointeur sur une famille d'histogrammes.
 *
 *--------------------------------------------------------------*/

bool Compound::plot_write(const char *prefix , const char *title ,
                          const Compound_data *compound_histo) const

{
  bool status;
  register int i;
  int nb_histo = 0 , index_dist[3];
  double scale[3];
  const Distribution *pdist[3];
  const Histogram *phisto[3];
  ostringstream data_file_name;


  // ecriture du fichier de donnees

  data_file_name << prefix << ".dat";

  pdist[0] = this;
  pdist[1] = sum_distribution;
  pdist[2] = distribution;

  if (compound_histo) {
    phisto[nb_histo] = compound_histo;
    index_dist[nb_histo++] = 0;
    scale[0] = compound_histo->nb_element;

    phisto[nb_histo] = compound_histo->sum_histogram;
    index_dist[nb_histo++] = 1;
    scale[1] = compound_histo->sum_histogram->nb_element;

    if (compound_histo->histogram->nb_element > 0) {
      phisto[nb_histo] = compound_histo->histogram;
      index_dist[nb_histo++] = 2;
      scale[2] = compound_histo->histogram->nb_element;
    }
    else {
      scale[2] = 1.;
    }

    status = ::plot_print((data_file_name.str()).c_str() , 3 , pdist , scale , 0 ,
                          nb_histo , phisto , index_dist);
  }

  else {
    scale[0] = 1.;
    scale[1] = 1.;
    scale[2] = 1.;

    status = ::plot_print((data_file_name.str()).c_str() , 3 , pdist , scale , 0 , 0 , 0 , 0);
  }

  // ecriture du fichier de commandes et du fichier d'impression

  if (status) {
    for (i = 0;i < 2;i++) {
      ostringstream file_name[2];

      switch (i) {
      case 0 :
        file_name[0] << prefix << ".plot";
        break;
      case 1 :
        file_name[0] << prefix << ".print";
        break;
      }

      ofstream out_file((file_name[0].str()).c_str());

      if (i == 1) {
        out_file << "set terminal postscript" << endl;
        file_name[1] << label(prefix) << ".ps";
        out_file << "set output \"" << file_name[1].str() << "\"\n\n";
      }

      out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
               << "set title";
      if (title) {
        out_file << " \"" << title << "\"";
      }
      out_file << "\n\n";

      if (nb_value - 1 < TIC_THRESHOLD) {
        out_file << "set xtics 0,1" << endl;
      }

      if (compound_histo) {
        out_file << "plot [0:" << nb_value - 1 << "] [0:"
                 << (int)(MAX(compound_histo->max , max * compound_histo->nb_element) * YSCALE) + 1
                 << "] \"" << label((data_file_name.str()).c_str()) << "\" using 1 title \""
                 << STAT_label[STATL_HISTOGRAM] << "\" with impulses,\\" << endl;
        out_file << "\"" << label((data_file_name.str()).c_str()) << "\" using " << nb_histo + 1 << " title \""
                 << STAT_label[STATL_COMPOUND] << " " << STAT_label[STATL_DISTRIBUTION]
                 << "\" with linespoints" << endl;
      }

      else {
        out_file << "plot [0:" << nb_value - 1 << "] [0:"
                 << MIN(max * YSCALE , 1.) << "] \""
                 << label((data_file_name.str()).c_str()) << "\" using 1 title \""
                 << STAT_label[STATL_COMPOUND] << " " << STAT_label[STATL_DISTRIBUTION]
                 << "\" with linespoints" << endl;
      }

      if (nb_value - 1 < TIC_THRESHOLD) {
        out_file << "set xtics autofreq" << endl;
      }

      if (i == 0) {
        out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
      }
      out_file << endl;

      if (sum_distribution->nb_value - 1 < TIC_THRESHOLD) {
        out_file << "set xtics 0,1" << endl;
      }

      if (compound_histo) {
        out_file << "plot [0:" << sum_distribution->nb_value - 1 << "] [0:"
                 << (int)(MAX(compound_histo->sum_histogram->max ,
                              sum_distribution->max * compound_histo->sum_histogram->nb_element) * YSCALE) + 1
                 << "] \"" << label((data_file_name.str()).c_str()) << "\" using 2 title \""
                 << STAT_label[STATL_SUM] << " " << STAT_label[STATL_HISTOGRAM]
                 << "\" with impulses,\\" << endl;
        out_file << "\"" << label((data_file_name.str()).c_str()) << "\" using " << nb_histo + 2 << " title \""
                 << STAT_label[STATL_SUM] << " " << STAT_label[STATL_DISTRIBUTION];
        sum_distribution->plot_title_print(out_file);
        out_file << "\" with linespoints" << endl;
      }

      else {
        out_file << "plot [0:" << sum_distribution->nb_value - 1 << "] [0:"
                 << MIN(sum_distribution->max * YSCALE , 1.) << "] \""
                 << label((data_file_name.str()).c_str()) << "\" using 2 title \""
                 << STAT_label[STATL_SUM] << " " << STAT_label[STATL_DISTRIBUTION];
        sum_distribution->plot_title_print(out_file);
        out_file << "\" with linespoints" << endl;
      }

      if (sum_distribution->nb_value - 1 < TIC_THRESHOLD) {
        out_file << "set xtics autofreq" << endl;
      }

      if (i == 0) {
        out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
      }
      out_file << endl;

      if (distribution->nb_value - 1 < TIC_THRESHOLD) {
        out_file << "set xtics 0,1" << endl;
      }

      if ((compound_histo) && (compound_histo->histogram->nb_element > 0)) {
        out_file << "plot [0:" << distribution->nb_value - 1 << "] [0:"
                 << (int)(MAX(compound_histo->histogram->max ,
                              distribution->max * compound_histo->histogram->nb_element) * YSCALE) + 1
                 << "] \"" << label((data_file_name.str()).c_str()) << "\" using 3 title \""
                 << STAT_label[STATL_ELEMENTARY] << " " << STAT_label[STATL_HISTOGRAM]
                 << "\" with impulses,\\" << endl;
        out_file << "\"" << label((data_file_name.str()).c_str()) << "\" using " << nb_histo + 3 << " title \""
                 << STAT_label[STATL_ELEMENTARY] << " " << STAT_label[STATL_DISTRIBUTION];
        distribution->plot_title_print(out_file);
        out_file << "\" with linespoints" << endl;
      }

      else {
        out_file << "plot [0:" << distribution->nb_value - 1 << "] [0:"
                 << MIN(distribution->max * YSCALE , 1.) << "] \""
                 << label((data_file_name.str()).c_str()) << "\" using " << nb_histo + 3 << " title \""
                 << STAT_label[STATL_ELEMENTARY] << " " << STAT_label[STATL_DISTRIBUTION];
        distribution->plot_title_print(out_file);
        out_file << "\" with linespoints" << endl;
      }

      if (distribution->nb_value - 1 < TIC_THRESHOLD) {
        out_file << "set xtics autofreq" << endl;
      }

      if (i == 1) {
        out_file << "\nset terminal x11" << endl;
      }

      out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet Compound.
 *
 *  arguments : reference sur un objet Format_error, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Compound::plot_write(Format_error &error , const char *prefix ,
                          const char *title) const

{
  bool status = plot_write(prefix , title , compound_data);

  error.init();

  if (!status) {
    error.update(STAT_error[STATR_FILE_PREFIX]);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Fonctions pour la persistance.
 *
 *--------------------------------------------------------------*/

/* RWDEFINE_COLLECTABLE(Compound , STATI_COMPOUND);


RWspace Compound::binaryStoreSize() const

{
  RWspace size = Distribution::binaryStoreSize() +
                 sum_distribution->binaryStoreSize() + distribution->binaryStoreSize();

  if (compound_data) {
    size += compound_data->recursiveStoreSize();
  }

  return size;
}


void Compound::restoreGuts(RWvistream &is)

{
  delete compound_data;

  delete sum_distribution;
  delete distribution;

  Distribution::restoreGuts(is);

  sum_distribution = new Parametric();
  sum_distribution->restoreGuts(is);
  distribution = new Parametric();
  distribution->restoreGuts(is);

  is >> compound_data;
  if (compound_data == RWnilCollectable) {
    compound_data = 0;
  }
}


void Compound::restoreGuts(RWFile &file)

{
  delete compound_data;

  delete sum_distribution;
  delete distribution;

  Distribution::restoreGuts(file);

  sum_distribution = new Parametric();
  sum_distribution->restoreGuts(file);
  distribution = new Parametric();
  distribution->restoreGuts(file);

  file >> compound_data;
  if (compound_data == RWnilCollectable) {
    compound_data = 0;
  }
}


void Compound::saveGuts(RWvostream &os) const

{
  Distribution::saveGuts(os);

  sum_distribution->saveGuts(os);
  distribution->saveGuts(os);

  if (compound_data) {
    os << compound_data;
  }
  else {
    os << RWnilCollectable;
  }
}


void Compound::saveGuts(RWFile &file) const

{
  Distribution::saveGuts(file);

  sum_distribution->saveGuts(file);
  distribution->saveGuts(file);

  if (compound_data) {
    file << compound_data;
  }
  else {
    file << RWnilCollectable;
  }
} */


/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe Compound_data.
 *
 *--------------------------------------------------------------*/

Compound_data::Compound_data()

{
  compound = 0;
  sum_histogram = 0;
  histogram = 0;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Compound_data.
 *
 *  arguments : reference sur un objet Histogram et
 *              sur un objet Compound.
 *
 *--------------------------------------------------------------*/

Compound_data::Compound_data(const Histogram &histo , const Compound &icompound)
:Histogram(histo)

{
  compound = 0;

  sum_histogram = new Histogram(icompound.sum_distribution->alloc_nb_value);
  histogram = new Histogram(icompound.distribution->alloc_nb_value);
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Compound_data.
 *
 *  argument : reference sur un objet Compound.
 *
 *--------------------------------------------------------------*/

Compound_data::Compound_data(const Compound &icompound)
:Histogram(icompound)

{
  compound = 0;

  sum_histogram = new Histogram(*(icompound.sum_distribution));
  histogram = new Histogram(*(icompound.distribution));
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Compound_data.
 *
 *  arguments : reference sur un objet Compound_data,
 *              flag copie de l'objet Compound.
 *
 *--------------------------------------------------------------*/

void Compound_data::copy(const Compound_data &compound_histo , bool model_flag)

{
  if ((model_flag) && (compound_histo.compound)) {
    compound = new Compound(*(compound_histo.compound) , false);
  }
  else {
    compound = 0;
  }

  sum_histogram = new Histogram(*(compound_histo.sum_histogram));
  histogram = new Histogram(*(compound_histo.histogram));
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Compound_data.
 *
 *--------------------------------------------------------------*/

Compound_data::~Compound_data()

{
  delete compound;

  delete sum_histogram;
  delete histogram;
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Compound_data.
 *
 *  argument : reference sur un objet Compound_data.
 *
 *--------------------------------------------------------------*/

Compound_data& Compound_data::operator=(const Compound_data &compound_histo)

{
  if (&compound_histo != this) {
    delete compound;

    delete sum_histogram;
    delete histogram;

    delete [] frequency;

    Histogram::copy(compound_histo);
    copy(compound_histo);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Extraction de l'histogramme correspondant a la loi de la somme ou
 *  de l'histogramme correspondant a la loi elementaire.
 *
 *  arguments : reference sur un objet Format_error, type de l'histogramme.
 *
 *--------------------------------------------------------------*/

Distribution_data* Compound_data::extract(Format_error &error , char type) const

{
  Distribution_data *phisto;


  error.init();

  switch (type) {

  case 's' : {
    phisto = new Distribution_data(*sum_histogram , (compound ? compound->sum_distribution : 0));
    break;
  }

  case 'e' : {
    if (histogram->nb_element == 0) {
      phisto = 0;
      error.update(STAT_error[STATR_EMPTY_HISTOGRAM]);
    }

    else {
      phisto = new Distribution_data(*histogram , (compound ? compound->distribution : 0));
    }
    break;
  }

  case 'c' : {
    phisto = new Distribution_data(*this , compound);
    break;
  }
  }

  return phisto;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet Compound_data.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& Compound_data::line_write(ostream &os) const

{
  os << STAT_label[STATL_SAMPLE_SIZE] << ": " << nb_element << "   "
     << STAT_label[STATL_MEAN] << ": " << mean << "   "
     << STAT_label[STATL_VARIANCE] << ": " << variance;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Compound_data.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Compound_data::ascii_write(ostream &os , bool exhaustive) const

{
  if (compound) {
    compound->ascii_write(os , this , exhaustive , false);
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Compound_data dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Compound_data::ascii_write(Format_error &error , const char *path ,
                                bool exhaustive) const

{
  bool status = false;


  if (compound) {
    ofstream out_file(path);

    error.init();

    if (!out_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_NAME]);
    }

    else {
      status = true;
      compound->ascii_write(out_file , this , exhaustive , true);
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Compound_data dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

bool Compound_data::spreadsheet_write(Format_error &error , const char *path) const

{
  bool status = false;


  if (compound) {
    ofstream out_file(path);

    error.init();

    if (!out_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_NAME]);
    }

    else {
      status = true;
      compound->spreadsheet_write(out_file , this);
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet Compound_data.
 *
 *  arguments : reference sur un objet Format_error, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Compound_data::plot_write(Format_error &error , const char *prefix ,
                               const char *title) const

{
  bool status = false;


  if (compound) {
    status = compound->plot_write(prefix , title , this);

    error.init();

    if (!status) {
      error.update(STAT_error[STATR_FILE_PREFIX]);
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Fonctions pour la persistance.
 *
 *--------------------------------------------------------------*/

/* RWDEFINE_COLLECTABLE(Compound_data , STATI_COMPOUND_DATA);


RWspace Compound_data::binaryStoreSize() const

{
  RWspace size = Histogram::binaryStoreSize() +
                 sum_histogram->binaryStoreSize() + histogram->binaryStoreSize();

  if (compound) {
    size += compound->recursiveStoreSize();
  }

  return size;
}


void Compound_data::restoreGuts(RWvistream &is)

{
  delete compound;

  delete sum_histogram;
  delete histogram;

  Histogram::restoreGuts(is);

  sum_histogram = new Histogram();
  sum_histogram->restoreGuts(is);
  histogram = new Histogram();
  histogram->restoreGuts(is);

  is >> compound;
  if (compound == RWnilCollectable) {
    compound = 0;
  }
}


void Compound_data::restoreGuts(RWFile &file)

{
  delete compound;

  delete sum_histogram;
  delete histogram;

  Histogram::restoreGuts(file);

  sum_histogram = new Histogram();
  sum_histogram->restoreGuts(file);
  histogram = new Histogram();
  histogram->restoreGuts(file);

  file >> compound;
  if (compound == RWnilCollectable) {
    compound = 0;
  }
}


void Compound_data::saveGuts(RWvostream &os) const

{
  Histogram::saveGuts(os);

  sum_histogram->saveGuts(os);
  histogram->saveGuts(os);

  if (compound) {
    os << compound;
  }
  else {
    os << RWnilCollectable;
  }
}


void Compound_data::saveGuts(RWFile &file) const

{
  Histogram::saveGuts(file);

  sum_histogram->saveGuts(file);
  histogram->saveGuts(file);

  if (compound) {
    file << compound;
  }
  else {
    file << RWnilCollectable;
  }
} */
