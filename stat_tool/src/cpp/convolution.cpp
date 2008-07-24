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



#include <math.h>
#include <sstream>
#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"
// #include <rw/vstream.h>
// #include <rw/rwfile.h>
#include "stat_tools.h"
#include "distribution.h"
#include "convolution.h"
#include "stat_label.h"

using namespace std;


extern char* label(const char *file_name);



/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe Convolution.
 *
 *--------------------------------------------------------------*/

Convolution::Convolution()

{
  convolution_data = 0;
  nb_distribution = 0;
  distribution = 0;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Convolution.
 *
 *  arguments : nombre de lois elementaires, pointeurs sur les lois elementaires.
 *
 *--------------------------------------------------------------*/

Convolution::Convolution(int nb_dist , const Parametric **pdist)

{
  register int i;
  int cnb_value = 1;


  convolution_data = 0;
  nb_distribution = nb_dist;

  distribution = new Parametric*[nb_distribution];
  for (i = 0;i < nb_distribution;i++) {
    distribution[i] = new Parametric(*pdist[i] , 'n');
    cnb_value += distribution[i]->nb_value - 1;
  }

  Distribution::init(cnb_value);
  computation(1 , CONVOLUTION_THRESHOLD);
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Convolution.
 *
 *  arguments : references sur la loi connue et sur la loi inconnue.
 *
 *--------------------------------------------------------------*/

Convolution::Convolution(const Parametric &known_dist , const Parametric &unknown_dist)

{
  convolution_data = 0;

  nb_distribution = 2;
  distribution = new Parametric*[nb_distribution];

  if ((known_dist.ident == POISSON) || (known_dist.ident == NEGATIVE_BINOMIAL)) {
    distribution[0] = new Parametric(known_dist.ident , known_dist.inf_bound , known_dist.sup_bound ,
                                     known_dist.parameter , known_dist.probability , CONVOLUTION_THRESHOLD);
  }
  else {
    distribution[0] = new Parametric(known_dist , 'n');
  }
  distribution[1] = new Parametric(unknown_dist , 'c' , (int)(unknown_dist.nb_value * NB_VALUE_COEFF));

  Distribution::init(distribution[0]->alloc_nb_value + distribution[1]->alloc_nb_value - 1);
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Convolution.
 *
 *  arguments : reference sur un objet Convolution,
 *              flag copie de l'objet Convolution_data.
 *
 *--------------------------------------------------------------*/

void Convolution::copy(const Convolution &convol , bool data_flag)

{
  register int i;


  if ((data_flag) && (convol.convolution_data)) {
    convolution_data = new Convolution_data(*(convol.convolution_data) , false);
  }
  else {
    convolution_data = 0;
  }

  nb_distribution = convol.nb_distribution;

  distribution = new Parametric*[nb_distribution];
  for (i = 0;i < nb_distribution;i++) {
    distribution[i] = new Parametric(*(convol.distribution[i]));
  }
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet Convolution.
 *
 *--------------------------------------------------------------*/

void Convolution::remove()

{
  delete convolution_data;

  if (distribution) {
    register int i;

    for (i = 0;i < nb_distribution;i++) {
      delete distribution[i];
    }
    delete [] distribution;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Convolution.
 *
 *--------------------------------------------------------------*/

Convolution::~Convolution()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Convolution.
 *
 *  argument : reference sur un objet Convolution.
 *
 *--------------------------------------------------------------*/

Convolution& Convolution::operator=(const Convolution &convol)

{
  if (&convol != this) {
    remove();
    delete [] mass;
    delete [] cumul;

    Distribution::copy(convol);
    copy(convol);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Extraction d'une loi elementaire.
 *
 *  arguments : reference sur un objet Format_error, indice de la loi.
 *
 *--------------------------------------------------------------*/

Parametric_model* Convolution::extract(Format_error &error , int index) const

{
  Parametric_model *pdist;


  if ((index < 1) || (index > nb_distribution)) {
    pdist = 0;
    error.update(STAT_error[STATR_DISTRIBUTION_INDEX]);
  }

  else {
    index--;
    pdist = new Parametric_model(*distribution[index] ,
                                 (convolution_data ? convolution_data->histogram[index] : 0));
  }

  return pdist;
}


/*--------------------------------------------------------------*
 *
 *  Extraction de la partie "donnees" d'un objet Convolution.
 *
 *  argument : reference sur un objet Format_error.
 *
 *--------------------------------------------------------------*/

Convolution_data* Convolution::extract_data(Format_error &error) const

{
  Convolution_data *convol_histo;


  error.init();

  if (!convolution_data) {
    convol_histo = 0;
    error.update(STAT_error[STATR_NO_DATA]);
  }

  else {
    convol_histo = new Convolution_data(*convolution_data);
    convol_histo->convolution = new Convolution(*this , false);
  }

  return convol_histo;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Convolution a partir de lois elementaires.
 *
 *  arguments : reference sur un objet Format_error, nombre de lois elementaires,
 *              pointeurs sur les lois elementaires.
 *
 *--------------------------------------------------------------*/

Convolution* convolution_building(Format_error &error , int nb_dist , const Parametric **dist)

{
  Convolution *convol;


  error.init();

  if ((nb_dist < 2) || (nb_dist > CONVOLUTION_NB_DISTRIBUTION)) {
    convol = 0;
    error.update(STAT_parsing[STATP_NB_DISTRIBUTION]);
  }

  else {
    convol = new Convolution(nb_dist , dist);
  }

  return convol;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Convolution a partir d'un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

Convolution* convolution_ascii_read(Format_error &error , const char *path ,
                                    double cumul_threshold)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  bool status , lstatus;
  register int i , j;
  int line;
  long index , nb_dist;
  const Parametric **dist;
  Convolution *convol;
  ifstream in_file(path);


  convol = 0;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    line = 0;
    nb_dist = 0;

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
        switch (i) {

        // test mot cle CONVOLUTION

        case 0 : {
          if (token != STAT_word[STATW_CONVOLUTION]) {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_CONVOLUTION] , line , i + 1);
          }
          break;
        }

        // test nombre de lois

        case 1 : {
          lstatus = locale.stringToNum(token , &nb_dist);
          if ((lstatus) && ((nb_dist < 2) || (nb_dist > CONVOLUTION_NB_DISTRIBUTION))) {
            lstatus = false;
          }

          if (!lstatus) {
            status = false;
            error.update(STAT_parsing[STATP_NB_DISTRIBUTION] , line , i + 1);
          }
          break;
        }

        // test mot cle DISTRIBUTIONS

        case 2 : {
          if (token != STAT_word[STATW_DISTRIBUTIONS]) {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_DISTRIBUTIONS] , line , i + 1);
          }
          break;
        }
        }

        i++;
      }

      if (i > 0) {
        if (i != 3) {
          status = false;
          error.update(STAT_parsing[STATP_FORMAT] , line);
        }

        break;
      }
    }

    if (nb_dist == 0) {
      status = false;
      error.update(STAT_parsing[STATP_FORMAT] , line);
    }

    if (status) {
      dist = new const Parametric*[nb_dist];
      for (i = 0;i < nb_dist;i++) {
        dist[i] = 0;
      }

      for (i = 0;i < nb_dist;i++) {
        while (buffer.readLine(in_file , false)) {
          line++;

#         ifdef DEBUG
          cout << line << "  " << buffer << endl;
#         endif

          position = buffer.first('#');
          if (position != RW_NPOS) {
            buffer.remove(position);
          }
          j = 0;

          RWCTokenizer next(buffer);

          while (!((token = next()).isNull())) {
            switch (j) {

            // test mot cle DISTRIBUTION

            case 0 : {
              if (token != STAT_word[STATW_DISTRIBUTION]) {
                status = false;
                error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_DISTRIBUTION] , line , j + 1);
              }
              break;
            }

            // test indice de la loi

            case 1 : {
              lstatus = locale.stringToNum(token , &index);
              if ((lstatus) && (index != i + 1)) {
                lstatus = false;
              }

              if (!lstatus) {
                status = false;
                error.correction_update(STAT_parsing[STATP_DISTRIBUTION_INDEX] , i + 1 , line , j + 1);
              }
              break;
            }
            }

            j++;
          }

          if (j > 0) {
            if (j != 2) {
              status = false;
              error.update(STAT_parsing[STATP_FORMAT] , line);
            }

            dist[i] = parametric_parsing(error , in_file , line ,
                                         NEGATIVE_BINOMIAL , cumul_threshold);
            break;
          }
        }

        if (!dist[i]) {
          status = false;
          error.update(STAT_parsing[STATP_FORMAT] , line);
        }
      }

      while (buffer.readLine(in_file , false)) {
        line++;

#       ifdef DEBUG
        cout << line << "  " << buffer << endl;
#       endif

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
        convol = new Convolution(nb_dist , dist);
      }

      for (i = 0;i < nb_dist;i++) {
        delete dist[i];
      }
      delete [] dist;
    }
  }

  return convol;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet Convolution.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& Convolution::line_write(ostream &os) const

{
  os << nb_distribution << " " << STAT_word[STATW_DISTRIBUTIONS] << "   "
     << STAT_label[STATL_MEAN] << ": " << mean << "   "
     << STAT_label[STATL_VARIANCE] << ": " << variance;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un produit de convolution et de la structure
 *  de donnees associee dans un fichier.
 *
 *  arguments : stream, pointeur sur un objet Convolution_data,
 *              flag niveau de detail, flag fichier.
 *
 *--------------------------------------------------------------*/

ostream& Convolution::ascii_write(ostream &os , const Convolution_data *convol_histo ,
                                  bool exhaustive , bool file_flag) const

{
  register int i;
  double scale[CONVOLUTION_NB_DISTRIBUTION];
  const Distribution *pdist[CONVOLUTION_NB_DISTRIBUTION];


  os << STAT_word[STATW_CONVOLUTION] << " " << nb_distribution << " " << STAT_word[STATW_DISTRIBUTIONS] << endl;
  ascii_characteristic_print(os , exhaustive , file_flag);

  if (convol_histo) {
    double likelihood , information;
    Test test(CHI2);


    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_HISTOGRAM] << " - ";
    convol_histo->ascii_characteristic_print(os , exhaustive , file_flag);

    likelihood = likelihood_computation(*convol_histo);
    information = convol_histo->information_computation();

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   ("
       << STAT_label[STATL_NORMALIZED] << ": " << likelihood / convol_histo->nb_element << ")" << endl;

    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_MAX_LIKELIHOOD] << ": " << information << "   ("
       << STAT_label[STATL_INFORMATION] << ": " << information / convol_histo->nb_element << ")" << endl;

    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_DEVIANCE] << ": " << 2 * (information - likelihood) << endl;

    chi2_fit(*convol_histo , test);
    os << "\n";
    test.ascii_print(os , file_flag);

    if (exhaustive) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << "   | " << STAT_label[STATL_HISTOGRAM] << " | " << STAT_label[STATL_CONVOLUTION]
         << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM] << " "
         << STAT_label[STATL_FUNCTION] << " | " << STAT_label[STATL_CUMULATIVE] << " "
         << STAT_label[STATL_CONVOLUTION] << " " << STAT_label[STATL_FUNCTION] << endl;

      ascii_print(os , file_flag , true , false , convol_histo);
    }

    for (i = 0;i < nb_distribution;i++) {
      os << "\n" << STAT_word[STATW_DISTRIBUTION] << " " << i + 1 << endl;
      distribution[i]->ascii_print(os);
      distribution[i]->ascii_characteristic_print(os , exhaustive , file_flag);
      if (exhaustive) {
        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_VARIANCE_MEAN_RATIO] << ": "
           << distribution[i]->variance / distribution[i]->mean << "   "
           << STAT_label[STATL_VARIATION_COEFF] << ": "
           << sqrt(distribution[i]->variance) / distribution[i]->mean << endl;
      }

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_HISTOGRAM] << " " << i + 1 << " - ";
      convol_histo->histogram[i]->ascii_characteristic_print(os , exhaustive , file_flag);

      if (exhaustive) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << "   | " << STAT_label[STATL_HISTOGRAM] << " " << i + 1
           << " | " << STAT_label[STATL_DISTRIBUTION] << " " << i + 1
           << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM] << " " << i + 1
           << " " << STAT_label[STATL_FUNCTION] << " | " << STAT_label[STATL_CUMULATIVE] << " "
           << STAT_label[STATL_DISTRIBUTION] << " " << i + 1 << " " << STAT_label[STATL_FUNCTION] << endl;

        distribution[i]->Distribution::ascii_print(os , file_flag , true , false , convol_histo->histogram[i]);
      }
    }
  }

  else {
    for (i = 0;i < nb_distribution;i++) {
      os << "\n" << STAT_word[STATW_DISTRIBUTION] << " " << i + 1 << endl;
      distribution[i]->ascii_print(os);
      distribution[i]->ascii_characteristic_print(os , exhaustive , file_flag);
    }
  }

  if (exhaustive) {
    for (i = 0;i < nb_distribution;i++) {
      pdist[i] = distribution[i];
      scale[i] = 1.;
    }

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << "   | " << STAT_label[STATL_CONVOLUTION];
    for (i = 0;i < nb_distribution;i++) {
      os << " | " << STAT_label[STATL_DISTRIBUTION] << " " << i + 1;
    }
    os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_CONVOLUTION] << " "
       << STAT_label[STATL_FUNCTION] << endl;

    ascii_print(os , nb_distribution , pdist , scale , file_flag , true);
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Convolution.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Convolution::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , convolution_data , exhaustive , false);
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Convolution dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Convolution::ascii_write(Format_error &error , const char *path ,
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
    ascii_write(out_file , convolution_data , exhaustive , true);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un produit de convolution et de la structure
 *  de donnees associee dans un fichier au format tableur.
 *
 *  arguments : stream, pointeur sur un objet Convolution_data.
 *
 *--------------------------------------------------------------*/

ostream& Convolution::spreadsheet_write(ostream &os , const Convolution_data *convol_histo) const

{
  register int i;
  double scale[CONVOLUTION_NB_DISTRIBUTION];
  const Distribution *pdist[CONVOLUTION_NB_DISTRIBUTION];


  os << STAT_word[STATW_CONVOLUTION] << "\t" << nb_distribution << "\t" << STAT_word[STATW_DISTRIBUTIONS] << endl;
  spreadsheet_characteristic_print(os , true);

  if (convol_histo) {
    double likelihood , information;
    Test test(CHI2);


    os << "\n" << STAT_label[STATL_HISTOGRAM] << "\t";
    convol_histo->spreadsheet_characteristic_print(os , true);

    likelihood = likelihood_computation(*convol_histo);
    information = convol_histo->information_computation();

    os << "\n" << STAT_label[STATL_LIKELIHOOD] << "\t" << likelihood << "\t"
       << STAT_label[STATL_NORMALIZED] << "\t" << likelihood / convol_histo->nb_element << endl;
    os << STAT_label[STATL_MAX_LIKELIHOOD] << "\t" << information << "\t"
       << STAT_label[STATL_INFORMATION] << "\t" << information / convol_histo->nb_element << endl;
    os << STAT_label[STATL_DEVIANCE] << "\t" << 2 * (information - likelihood) << endl;

    chi2_fit(*convol_histo , test);
    os << "\n";
    test.spreadsheet_print(os);

    os << "\n\t" << STAT_label[STATL_HISTOGRAM] << "\t" << STAT_label[STATL_CONVOLUTION]
       << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM] << " "
       << STAT_label[STATL_FUNCTION] << "\t" << STAT_label[STATL_CUMULATIVE] << " "
       << STAT_label[STATL_CONVOLUTION] << " " << STAT_label[STATL_FUNCTION] << endl;

    spreadsheet_print(os , true , false , false , convol_histo);

    for (i = 0;i < nb_distribution;i++) {
      os << "\n" << STAT_word[STATW_DISTRIBUTION] << "\t" << i + 1 << endl;
      distribution[i]->spreadsheet_print(os);
      distribution[i]->spreadsheet_characteristic_print(os , true);
      os << STAT_label[STATL_VARIANCE_MEAN_RATIO] << "\t"
         << distribution[i]->variance / distribution[i]->mean << "\t"
         << STAT_label[STATL_VARIATION_COEFF] << "\t"
         << sqrt(distribution[i]->variance) / distribution[i]->mean << endl;

      os << "\n" << STAT_label[STATL_HISTOGRAM] << " " << i + 1 << "\t";
      convol_histo->histogram[i]->spreadsheet_characteristic_print(os , true);

      os << "\n\t" << STAT_label[STATL_HISTOGRAM] << " " << i + 1
         << "\t" << STAT_label[STATL_DISTRIBUTION] << " " << i + 1
         << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM] << " "
         << i + 1 << " " << STAT_label[STATL_FUNCTION] << "\t" << STAT_label[STATL_CUMULATIVE] << " "
         << STAT_label[STATL_DISTRIBUTION] << " " << i + 1 << " " << STAT_label[STATL_FUNCTION] << endl;

      distribution[i]->Distribution::spreadsheet_print(os , true , false , false , convol_histo->histogram[i]);
    }
  }

  else {
    for (i = 0;i < nb_distribution;i++) {
      os << "\n" << STAT_word[STATW_DISTRIBUTION] << "\t" << i + 1 << endl;
      distribution[i]->spreadsheet_print(os);
      distribution[i]->spreadsheet_characteristic_print(os , true);
    }
  }

  for (i = 0;i < nb_distribution;i++) {
    pdist[i] = distribution[i];
    scale[i] = 1.;
  }

  os << "\n\t" << STAT_label[STATL_CONVOLUTION];
  for (i = 0;i < nb_distribution;i++) {
    os << "\t" << STAT_label[STATL_DISTRIBUTION] << " " << i + 1;
  }
  os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_CONVOLUTION] << " "
     << STAT_label[STATL_FUNCTION] << endl;

  spreadsheet_print(os , nb_distribution , pdist , scale , true);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Convolution dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

bool Convolution::spreadsheet_write(Format_error &error , const char *path) const

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
    spreadsheet_write(out_file , convolution_data);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un produit de convolution et
 *  de la structure de donnees associee.
 *
 *  arguments : prefixe des fichiers, titre des figures,
 *              pointeur sur un objet Convolution_data.
 *
 *--------------------------------------------------------------*/

bool Convolution::plot_write(const char *prefix , const char *title ,
                             const Convolution_data *convol_histo) const

{
  bool status;
  register int i , j;
  int index_dist[CONVOLUTION_NB_DISTRIBUTION + 1];
  double plot_max , scale[CONVOLUTION_NB_DISTRIBUTION + 2];
  const Distribution *pdist[CONVOLUTION_NB_DISTRIBUTION + 2];
  const Histogram *phisto[CONVOLUTION_NB_DISTRIBUTION + 1];
  ostringstream data_file_name[CONVOLUTION_NB_DISTRIBUTION + 1];


  // ecriture des fichiers de donnees

  data_file_name[0] << prefix << ".dat";

  if (convol_histo) {
    pdist[0] = this;
    scale[0] = 1.;

    pdist[1] = this;
    phisto[0] = convol_histo;
    index_dist[0] = 1;
    scale[1] = convol_histo->nb_element;

    for (i = 0;i < nb_distribution;i++) {
      pdist[i + 2] = distribution[i];
      phisto[i + 1] = convol_histo->histogram[i];
      index_dist[i + 1] = i + 2;
      scale[i + 2] = convol_histo->histogram[i]->nb_element;
    }

    status = ::plot_print((data_file_name[0].str()).c_str() , nb_distribution + 2 , pdist ,
                          scale , 0 , nb_distribution + 1 , phisto , index_dist);
  }

  else {
    status = plot_print((data_file_name[0].str()).c_str());
  }

  if (status) {
    for (i = 0;i < nb_distribution;i++) {
      data_file_name[i + 1] << prefix << i + 1 << ".dat";
      distribution[i]->plot_print((data_file_name[i + 1].str()).c_str());
    }

    // ecriture du fichier de commandes et du fichier d'impression

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

      plot_max = distribution[0]->max;
      for (j = 1;j < nb_distribution;j++) {
        if (distribution[j]->max > plot_max) {
          plot_max = distribution[j]->max;
        }
      }

      out_file << "plot [0:" << nb_value - 1 << "] [0:"
               << MIN(plot_max * YSCALE , 1.) << "] \""
               << label((data_file_name[0].str()).c_str()) << "\"";
      if (convol_histo) {
        out_file << " using " << nb_distribution + 2;
      }
      out_file << " title \"" << STAT_label[STATL_CONVOLUTION] << "\" with linespoints";

      for (j = 0;j < nb_distribution;j++) {
        out_file << ",\\" << endl;
        out_file << "\"" << label((data_file_name[j + 1].str()).c_str()) << "\" title \""
                 << STAT_label[STATL_DISTRIBUTION] << " " << j + 1;
        distribution[j]->plot_title_print(out_file);
        out_file << "\" with linespoints";
      }
      out_file << endl;

      if (convol_histo) {
        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;

        out_file << "plot [0:" << nb_value - 1 << "] [0:"
                 << (int)(MAX(convol_histo->max , max * convol_histo->nb_element) * YSCALE) + 1
                 << "] \"" << label((data_file_name[0].str()).c_str()) << "\" using 1 title \""
                 << STAT_label[STATL_HISTOGRAM] << "\" with impulses,\\" << endl;
        out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << nb_distribution + 3
                 << " title \"" << STAT_label[STATL_CONVOLUTION] << "\" with linespoints" << endl;
      }

      if (nb_value - 1 < TIC_THRESHOLD) {
        out_file << "set xtics autofreq" << endl;
      }

      if (convol_histo) {
        for (j = 0;j < nb_distribution;j++) {
          if (i == 0) {
            out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
          out_file << endl;

          if (distribution[j]->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }

          out_file << "plot [0:" << distribution[j]->nb_value - 1 << "] [0:"
                   << (int)(MAX(convol_histo->histogram[j]->max ,
                                distribution[j]->max * convol_histo->histogram[j]->nb_element) * YSCALE) + 1
                   << "] \"" << label((data_file_name[0].str()).c_str()) << "\" using " << j + 2
                   << " title \"" << STAT_label[STATL_HISTOGRAM] << " " << j + 1
                   << "\" with impulses,\\" << endl;
          out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << nb_distribution + j + 4
                   << " title \"" << STAT_label[STATL_DISTRIBUTION] << " " << j + 1;
          distribution[j]->plot_title_print(out_file);
          out_file << "\" with linespoints" << endl;

          if (distribution[j]->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }
        }
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
 *  Sortie Gnuplot d'un objet Convolution.
 *
 *  arguments : reference sur un objet Format_error, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Convolution::plot_write(Format_error &error , const char *prefix ,
                             const char *title) const

{
  bool status = plot_write(prefix , title , convolution_data);

  error.init();

  if (!status) {
    error.update(STAT_error[STATR_FILE_PREFIX]);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un produit de convolution et
 *  de la structure de donnees associee.
 *
 *  argument : pointeur sur un objet Convolution_data.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* Convolution::get_plotable(const Convolution_data *convol_histo) const

{
  register int i;
  double plot_max;
  std::stringstream ss;


  // nombre de fenetres: nb_distribution + 2 si ajustement

  MultiPlotSet *plotset = new MultiPlotSet(convol_histo ? nb_distribution + 2 : 1);
  MultiPlotSet &set = *plotset;

  set.title = "Convolution fit";
  set.border = "15 lw 0";

  // 1ere vue : produit de convolution

  if (nb_value - 1 < TIC_THRESHOLD) {
    set[0].xtics = 1;
  }

  plot_max = distribution[0]->max;
  for (i = 1;i < nb_distribution;i++) {
    if (distribution[i]->max > plot_max) {
      plot_max = distribution[i]->max;
    }
  }

  set[0].xrange = Range(0, nb_value - 1);
  set[0].yrange = Range(0, MIN(plot_max * YSCALE , 1.));

  // definition du nombre de SinglePlot 

  set[0].resize(nb_distribution + 1);

  set[0][0].legend = STAT_label[STATL_CONVOLUTION];
  set[0][0].style = "linespoints";

  plotable_mass_write(set[0][0]);

  for (i = 0;i < nb_distribution;i++) {
    ss.str("");
    ss << STAT_label[STATL_DISTRIBUTION] << " " << i + 1;
    set[0][i + 1].legend = ss.str();

    set[0][i + 1].style = "linespoints";

    distribution[i]->plotable_mass_write(set[0][i + 1]);
  }

  if (convol_histo) {

    // 2eme vue : produit de convolution ajuste

    if (nb_value - 1 < TIC_THRESHOLD) {
      set[1].xtics = 1;
    }

    set[1].xrange = Range(0, nb_value - 1);
    set[1].yrange = Range(0, (MAX(convol_histo->max , max * convol_histo->nb_element)
                              * YSCALE) + 1);

    set[1].resize(2);

    ss.str("");
    ss << STAT_label[STATL_HISTOGRAM];
    set[1][0].legend = ss.str();
    set[1][0].style = "impulses";

    ss.str("");
    ss << STAT_label[STATL_CONVOLUTION];
    set[1][1].legend = ss.str();
    set[1][1].style = "linespoints";

    convol_histo->plotable_frequency_write(set[1][0]);
    plotable_mass_write(set[1][1] , convol_histo->nb_element);

    for (i = 0;i < nb_distribution;i++) {

      // vues suivantes : distributions ajustees

      if (distribution[i]->nb_value - 1 < TIC_THRESHOLD) {
        set[i + 2].xtics = 1;
      }

      set[i + 2].xrange = Range(0 , distribution[i]->nb_value - 1);
      set[i + 2].yrange = Range(0 , (MAX(convol_histo->histogram[i]->max ,
                                         distribution[i]->max * convol_histo->histogram[i]->nb_element)
                                     * YSCALE) + 1);

      set[i + 2].resize(2);

      ss.str("");
      ss << STAT_label[STATL_HISTOGRAM] << " " << i + 1;
      set[i + 2][0].legend = ss.str();
      set[i + 2][0].style = "impulses";

      ss.str("");
      ss << STAT_label[STATL_DISTRIBUTION] << " " << i + 1;
      set[i + 2][1].legend = ss.str();
      set[i + 2][1].style = "linespoints";

      convol_histo->histogram[i]->plotable_frequency_write(set[i + 2][0]);
      distribution[i]->plotable_mass_write(set[i + 2][1] , convol_histo->histogram[i]->nb_element);
    }
  }

  return plotset;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet Convolution.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* Convolution::get_plotable() const

{
  MultiPlotSet *set;


  set = get_plotable(convolution_data);

  return set;
}


/*--------------------------------------------------------------*
 *
 *  Fonctions pour la persistance.
 *
 *--------------------------------------------------------------*/

/* RWDEFINE_COLLECTABLE(Convolution , STATI_CONVOLUTION);


RWspace Convolution::binaryStoreSize() const

{
  register int i;
  RWspace size;


  size = Distribution::binaryStoreSize() + sizeof(nb_distribution);
  for (i = 0;i < nb_distribution;i++) {
    size += distribution[i]->binaryStoreSize();
  }

  if (convolution_data) {
    size += convolution_data->recursiveStoreSize();
  }

  return size;
}


void Convolution::restoreGuts(RWvistream &is)

{
  register int i;


  remove();

  Distribution::restoreGuts(is);

  is >> nb_distribution;

  distribution = new Parametric*[nb_distribution];
  for (i = 0;i < nb_distribution;i++) {
    distribution[i] = new Parametric();
    distribution[i]->restoreGuts(is);
  }

  is >> convolution_data;
  if (convolution_data == RWnilCollectable) {
    convolution_data = 0;
  }
}


void Convolution::restoreGuts(RWFile &file)

{
  register int i;


  remove();

  Distribution::restoreGuts(file);

  file.Read(nb_distribution);

  distribution = new Parametric*[nb_distribution];
  for (i = 0;i < nb_distribution;i++) {
    distribution[i] = new Parametric();
    distribution[i]->restoreGuts(file);
  }

  file >> convolution_data;
  if (convolution_data == RWnilCollectable) {
    convolution_data = 0;
  }
}


void Convolution::saveGuts(RWvostream &os) const

{
  register int i;


  Distribution::saveGuts(os);

  os << nb_distribution;

  for (i = 0;i < nb_distribution;i++) {
    distribution[i]->saveGuts(os);
  }

  if (convolution_data) {
    os << convolution_data;
  }
  else {
    os << RWnilCollectable;
  }
}


void Convolution::saveGuts(RWFile &file) const

{
  register int i;


  Distribution::saveGuts(file);

  file.Write(nb_distribution);

  for (i = 0;i < nb_distribution;i++) {
    distribution[i]->saveGuts(file);
  }

  if (convolution_data) {
    file << convolution_data;
  }
  else {
    file << RWnilCollectable;
  }
} */


/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe Convolution_data.
 *
 *--------------------------------------------------------------*/

Convolution_data::Convolution_data()

{
  convolution = 0;
  nb_histogram = 0;
  histogram = 0;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Convolution_data.
 *
 *  arguments : reference sur un objet Histogram,
 *              nombre d'histogrammes.
 *
 *--------------------------------------------------------------*/

Convolution_data::Convolution_data(const Histogram &histo , int nb_histo)
:Histogram(histo)

{
  register int i;


  convolution = 0;
  nb_histogram = nb_histo;

  histogram = new Histogram*[nb_histogram];
  for (i = 0;i < nb_histogram;i++) {
    histogram[i] = new Histogram(nb_value);
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Convolution_data.
 *
 *  argument : reference sur un objet Convolution.
 *
 *--------------------------------------------------------------*/

Convolution_data::Convolution_data(const Convolution &convol)
:Histogram(convol)

{
  register int i;


  convolution = 0;
  nb_histogram = convol.nb_distribution;

  histogram = new Histogram*[nb_histogram];
  for (i = 0;i < nb_histogram;i++) {
    histogram[i] = new Histogram(*(convol.distribution[i]));
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Convolution_data.
 *
 *  arguments : reference sur un objet Convolution_data,
 *              flag copie de l'objet Convolution.
 *
 *--------------------------------------------------------------*/

void Convolution_data::copy(const Convolution_data &convol_histo , bool model_flag)

{
  register int i;


  if ((model_flag) && (convol_histo.convolution)) {
    convolution = new Convolution(*(convol_histo.convolution) , false);
  }
  else {
    convolution = 0;
  }

  nb_histogram = convol_histo.nb_histogram;

  histogram = new Histogram*[nb_histogram];
  for (i = 0;i < nb_histogram;i++) {
    histogram[i] = new Histogram(*(convol_histo.histogram[i]));
  }
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet Convolution_data.
 *
 *--------------------------------------------------------------*/

void Convolution_data::remove()

{
  delete convolution;

  if (histogram) {
    register int i;

    for (i = 0;i < nb_histogram;i++) {
      delete histogram[i];
    }
    delete [] histogram;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Convolution_data.
 *
 *--------------------------------------------------------------*/

Convolution_data::~Convolution_data()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Convolution_data.
 *
 *  argument : reference sur un objet Convolution_data.
 *
 *--------------------------------------------------------------*/

Convolution_data& Convolution_data::operator=(const Convolution_data &convol_histo)

{
  if (&convol_histo != this) {
    remove();
    delete [] frequency;

    Histogram::copy(convol_histo);
    copy(convol_histo);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Extraction d'un sous-histogramme.
 *
 *  arguments : reference sur un objet Format_error, indice de l'histogramme.
 *
 *--------------------------------------------------------------*/

Distribution_data* Convolution_data::extract(Format_error &error , int index) const

{
  Distribution_data *phisto;


  error.init();

  if ((index < 1) || (index > nb_histogram)) {
    phisto = 0;
    error.update(STAT_error[STATR_HISTOGRAM_INDEX]);
  }

  else {
    index--;
    phisto = new Distribution_data(*histogram[index] ,
                                   (convolution ? convolution->distribution[index] : 0));
  }

  return phisto;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet Convolution_data.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& Convolution_data::line_write(ostream &os) const

{
  os << STAT_label[STATL_SAMPLE_SIZE] << ": " << nb_element << "   "
     << STAT_label[STATL_MEAN] << ": " << mean << "   "
     << STAT_label[STATL_VARIANCE] << ": " << variance;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Convolution_data.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Convolution_data::ascii_write(ostream &os , bool exhaustive) const

{
  if (convolution) {
    convolution->ascii_write(os , this , exhaustive , false);
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Convolution_data dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Convolution_data::ascii_write(Format_error &error , const char *path ,
                                   bool exhaustive) const

{
  bool status = false;


  if (convolution) {
    ofstream out_file(path);

    error.init();

    if (!out_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_NAME]);
    }

    else {
      status = true;
      convolution->ascii_write(out_file , this , exhaustive , true);
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Convolution_data dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

bool Convolution_data::spreadsheet_write(Format_error &error , const char *path) const

{
  bool status = false;


  if (convolution) {
    ofstream out_file(path);

    error.init();

    if (!out_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_NAME]);
    }

    else {
      status = true;
      convolution->spreadsheet_write(out_file , this);
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet Convolution_data.
 *
 *  arguments : reference sur un objet Format_error, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Convolution_data::plot_write(Format_error &error , const char *prefix ,
                                  const char *title) const

{
  bool status = false;


  if (convolution) {
    status = convolution->plot_write(prefix , title , this);

    error.init();

    if (!status) {
      error.update(STAT_error[STATR_FILE_PREFIX]);
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet Convolution_data.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* Convolution_data::get_plotable() const

{
  MultiPlotSet *set;


  if (convolution) {
    set = convolution->get_plotable(this);
  }
  else {
    set = 0;
  }

  return set;
}


/*--------------------------------------------------------------*
 *
 *  Fonctions pour la persistance.
 *
 *--------------------------------------------------------------*/

/* RWDEFINE_COLLECTABLE(Convolution_data , STATI_CONVOLUTION_DATA);


RWspace Convolution_data::binaryStoreSize() const

{
  register int i;
  RWspace size;


  size = Histogram::binaryStoreSize() + sizeof(nb_histogram);
  for (i = 0;i < nb_histogram;i++) {
    size += histogram[i]->binaryStoreSize();
  }

  if (convolution) {
    size += convolution->recursiveStoreSize();
  }

  return size;
}


void Convolution_data::restoreGuts(RWvistream &is)

{
  register int i;


  remove();

  Histogram::restoreGuts(is);

  is >> nb_histogram;

  histogram = new Histogram*[nb_histogram];
  for (i = 0;i < nb_histogram;i++) {
    histogram[i] = new Histogram();
    histogram[i]->restoreGuts(is);
  }

  is >> convolution;
  if (convolution == RWnilCollectable) {
    convolution = 0;
  }
}


void Convolution_data::restoreGuts(RWFile &file)

{
  register int i;


  remove();

  Histogram::restoreGuts(file);

  file.Read(nb_histogram);

  histogram = new Histogram*[nb_histogram];
  for (i = 0;i < nb_histogram;i++) {
    histogram[i] = new Histogram();
    histogram[i]->restoreGuts(file);
  }

  file >> convolution;
  if (convolution == RWnilCollectable) {
    convolution = 0;
  }
}


void Convolution_data::saveGuts(RWvostream &os) const

{
  register int i;


  Histogram::saveGuts(os);

  os << nb_histogram;

  for (i = 0;i < nb_histogram;i++) {
    histogram[i]->saveGuts(os);
  }

  if (convolution) {
    os << convolution;
  }
  else {
    os << RWnilCollectable;
  }
}


void Convolution_data::saveGuts(RWFile &file) const

{
  register int i;


  Histogram::saveGuts(file);

  file.Write(nb_histogram);

  for (i = 0;i < nb_histogram;i++) {
    histogram[i]->saveGuts(file);
  }

  if (convolution) {
    file << convolution;
  }
  else {
    file << RWnilCollectable;
  }
} */
