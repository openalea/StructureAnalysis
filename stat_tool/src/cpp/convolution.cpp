/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2015 CIRAD/INRA/Inria Virtual Plants
 *
 *       File author(s): Yann Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id: convolution.cpp 17984 2015-04-23 06:42:25Z guedon $
 *
 *       Forum for V-Plants developers:
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

#include "stat_tools.h"
#include "distribution.h"
#include "convolution.h"
#include "stat_label.h"

using namespace std;


namespace stat_tool {



/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe Convolution.
 *
 *--------------------------------------------------------------*/

Convolution::Convolution()

{
  convolution_data = NULL;
  nb_distribution = 0;
  distribution = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Convolution.
 *
 *  arguments : nombre de lois elementaires, pointeurs sur les lois elementaires.
 *
 *--------------------------------------------------------------*/

Convolution::Convolution(int nb_dist , const DiscreteParametric **pdist)

{
  register int i;
  int cnb_value = 1;


  convolution_data = NULL;
  nb_distribution = nb_dist;

  distribution = new DiscreteParametric*[nb_distribution];
  for (i = 0;i < nb_distribution;i++) {
    distribution[i] = new DiscreteParametric(*pdist[i] , 'n');
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

Convolution::Convolution(const DiscreteParametric &known_dist ,
                         const DiscreteParametric &unknown_dist)

{
  convolution_data = NULL;

  nb_distribution = 2;
  distribution = new DiscreteParametric*[nb_distribution];

  if ((known_dist.ident == POISSON) || (known_dist.ident == NEGATIVE_BINOMIAL)) {
    distribution[0] = new DiscreteParametric(known_dist.ident , known_dist.inf_bound ,
                                             known_dist.sup_bound , known_dist.parameter ,
                                             known_dist.probability , CONVOLUTION_THRESHOLD);
  }
  else {
    distribution[0] = new DiscreteParametric(known_dist , 'n');
  }
  distribution[1] = new DiscreteParametric(unknown_dist , 'c' , (int)(unknown_dist.nb_value * NB_VALUE_COEFF));

  Distribution::init(distribution[0]->alloc_nb_value + distribution[1]->alloc_nb_value - 1);
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Convolution.
 *
 *  arguments : reference sur un objet Convolution,
 *              flag copie de l'objet ConvolutionData.
 *
 *--------------------------------------------------------------*/

void Convolution::copy(const Convolution &convol , bool data_flag)

{
  register int i;


  if ((data_flag) && (convol.convolution_data)) {
    convolution_data = new ConvolutionData(*(convol.convolution_data) , false);
  }
  else {
    convolution_data = NULL;
  }

  nb_distribution = convol.nb_distribution;

  distribution = new DiscreteParametric*[nb_distribution];
  for (i = 0;i < nb_distribution;i++) {
    distribution[i] = new DiscreteParametric(*(convol.distribution[i]));
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
 *  arguments : reference sur un objet StatError, indice de la loi.
 *
 *--------------------------------------------------------------*/

DiscreteParametricModel* Convolution::extract(StatError &error , int index) const

{
  DiscreteParametricModel *pdist;


  if ((index < 1) || (index > nb_distribution)) {
    pdist = NULL;
    error.update(STAT_error[STATR_DISTRIBUTION_INDEX]);
  }

  else {
    index--;
    pdist = new DiscreteParametricModel(*distribution[index] ,
                                        (convolution_data ? convolution_data->frequency_distribution[index] : NULL));
  }

  return pdist;
}


/*--------------------------------------------------------------*
 *
 *  Extraction de la partie "donnees" d'un objet Convolution.
 *
 *  argument : reference sur un objet StatError.
 *
 *--------------------------------------------------------------*/

ConvolutionData* Convolution::extract_data(StatError &error) const

{
  ConvolutionData *convol_histo;


  error.init();

  if (!convolution_data) {
    convol_histo = NULL;
    error.update(STAT_error[STATR_NO_DATA]);
  }

  else {
    convol_histo = new ConvolutionData(*convolution_data);
    convol_histo->convolution = new Convolution(*this , false);
  }

  return convol_histo;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Convolution a partir de lois elementaires.
 *
 *  arguments : reference sur un objet StatError, nombre de lois elementaires,
 *              pointeurs sur les lois elementaires.
 *
 *--------------------------------------------------------------*/

Convolution* convolution_building(StatError &error , int nb_dist ,
                                  const DiscreteParametric **dist)

{
  Convolution *convol;


  error.init();

  if ((nb_dist < 2) || (nb_dist > CONVOLUTION_NB_DISTRIBUTION)) {
    convol = NULL;
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
 *  arguments : reference sur un objet StatError, path,
 *              seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

Convolution* convolution_ascii_read(StatError &error , const char *path ,
                                    double cumul_threshold)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  bool status , lstatus;
  register int i , j;
  int line;
  long index , nb_dist;
  const DiscreteParametric **dist;
  Convolution *convol;
  ifstream in_file(path);


  convol = NULL;
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
      dist = new const DiscreteParametric*[nb_dist];
      for (i = 0;i < nb_dist;i++) {
        dist[i] = NULL;
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

            dist[i] = discrete_parametric_parsing(error , in_file , line ,
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
 *  arguments : stream, pointeur sur un objet ConvolutionData,
 *              flag niveau de detail, flag fichier.
 *
 *--------------------------------------------------------------*/

ostream& Convolution::ascii_write(ostream &os , const ConvolutionData *convol_histo ,
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
    os << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
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
      os << "   | " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " | " << STAT_label[STATL_CONVOLUTION]
         << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
         << " " << STAT_label[STATL_FUNCTION] << " | " << STAT_label[STATL_CUMULATIVE]
         << " " << STAT_label[STATL_CONVOLUTION] << " " << STAT_label[STATL_FUNCTION] << endl;

      ascii_print(os , file_flag , true , false , convol_histo);
    }

    for (i = 0;i < nb_distribution;i++) {
      os << "\n" << STAT_word[STATW_DISTRIBUTION] << " " << i + 1 << endl;
      distribution[i]->ascii_print(os);
      distribution[i]->ascii_parametric_characteristic_print(os , exhaustive , file_flag);
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
      os << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i + 1 << " - ";
      convol_histo->frequency_distribution[i]->ascii_characteristic_print(os , exhaustive , file_flag);

      if (exhaustive) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << "   | " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i + 1
           << " | " << STAT_label[STATL_DISTRIBUTION] << " " << i + 1
           << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
           << " " << i + 1 << " " << STAT_label[STATL_FUNCTION] << " | " << STAT_label[STATL_CUMULATIVE]
           << " " << STAT_label[STATL_DISTRIBUTION] << " " << i + 1 << " " << STAT_label[STATL_FUNCTION] << endl;

        distribution[i]->Distribution::ascii_print(os , file_flag , true , false ,
                                                   convol_histo->frequency_distribution[i]);
      }
    }
  }

  else {
    for (i = 0;i < nb_distribution;i++) {
      os << "\n" << STAT_word[STATW_DISTRIBUTION] << " " << i + 1 << endl;
      distribution[i]->ascii_print(os);
      distribution[i]->ascii_parametric_characteristic_print(os , exhaustive , file_flag);
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
    os << "  ";
    for (i = 0;i < nb_distribution;i++) {
      os << " | " << STAT_label[STATL_DISTRIBUTION] << " " << i + 1;
    }
    os << " | " << STAT_label[STATL_CONVOLUTION] << " | " << STAT_label[STATL_CUMULATIVE]
       << " " << STAT_label[STATL_CONVOLUTION] << " " << STAT_label[STATL_FUNCTION] << endl;

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
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Convolution::ascii_write(StatError &error , const char *path ,
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
 *  arguments : stream, pointeur sur un objet ConvolutionData.
 *
 *--------------------------------------------------------------*/

ostream& Convolution::spreadsheet_write(ostream &os , const ConvolutionData *convol_histo) const

{
  register int i;
  double scale[CONVOLUTION_NB_DISTRIBUTION];
  const Distribution *pdist[CONVOLUTION_NB_DISTRIBUTION];


  os << STAT_word[STATW_CONVOLUTION] << "\t" << nb_distribution << "\t" << STAT_word[STATW_DISTRIBUTIONS] << endl;
  spreadsheet_characteristic_print(os , true);

  if (convol_histo) {
    double likelihood , information;
    Test test(CHI2);


    os << "\n" << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
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

    os << "\n\t" << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t" << STAT_label[STATL_CONVOLUTION]
       << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
       << " " << STAT_label[STATL_FUNCTION] << "\t" << STAT_label[STATL_CUMULATIVE]
       << " " << STAT_label[STATL_CONVOLUTION] << " " << STAT_label[STATL_FUNCTION] << endl;

    spreadsheet_print(os , true , false , false , convol_histo);

    for (i = 0;i < nb_distribution;i++) {
      os << "\n" << STAT_word[STATW_DISTRIBUTION] << "\t" << i + 1 << endl;
      distribution[i]->spreadsheet_print(os);
      distribution[i]->spreadsheet_parametric_characteristic_print(os , true);
      os << STAT_label[STATL_VARIANCE_MEAN_RATIO] << "\t"
         << distribution[i]->variance / distribution[i]->mean << "\t"
         << STAT_label[STATL_VARIATION_COEFF] << "\t"
         << sqrt(distribution[i]->variance) / distribution[i]->mean << endl;

      os << "\n" << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i + 1 << "\t";
      convol_histo->frequency_distribution[i]->spreadsheet_characteristic_print(os , true);

      os << "\n\t" << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i + 1
         << "\t" << STAT_label[STATL_DISTRIBUTION] << " " << i + 1
         << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
         << " " << i + 1 << " " << STAT_label[STATL_FUNCTION] << "\t" << STAT_label[STATL_CUMULATIVE]
         << " " << STAT_label[STATL_DISTRIBUTION] << " " << i + 1 << " " << STAT_label[STATL_FUNCTION] << endl;

      distribution[i]->Distribution::spreadsheet_print(os , true , false , false ,
                                                       convol_histo->frequency_distribution[i]);
    }
  }

  else {
    for (i = 0;i < nb_distribution;i++) {
      os << "\n" << STAT_word[STATW_DISTRIBUTION] << "\t" << i + 1 << endl;
      distribution[i]->spreadsheet_print(os);
      distribution[i]->spreadsheet_parametric_characteristic_print(os , true);
    }
  }

  for (i = 0;i < nb_distribution;i++) {
    pdist[i] = distribution[i];
    scale[i] = 1.;
  }

  os << "\n";
  for (i = 0;i < nb_distribution;i++) {
    os << "\t" << STAT_label[STATL_DISTRIBUTION] << " " << i + 1;
  }
  os << "\t" << STAT_label[STATL_CONVOLUTION] << "\t" << STAT_label[STATL_CUMULATIVE]
     << " " << STAT_label[STATL_CONVOLUTION] << " " << STAT_label[STATL_FUNCTION] << endl;

  spreadsheet_print(os , nb_distribution , pdist , scale , true);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Convolution dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool Convolution::spreadsheet_write(StatError &error , const char *path) const

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
 *              pointeur sur un objet ConvolutionData.
 *
 *--------------------------------------------------------------*/

bool Convolution::plot_write(const char *prefix , const char *title ,
                             const ConvolutionData *convol_histo) const

{
  bool status;
  register int i , j;
  double plot_max , scale[CONVOLUTION_NB_DISTRIBUTION + 2];
  const Distribution *pdist[CONVOLUTION_NB_DISTRIBUTION + 2];
  const FrequencyDistribution *phisto[CONVOLUTION_NB_DISTRIBUTION + 1];
  ostringstream data_file_name[CONVOLUTION_NB_DISTRIBUTION + 1];


  // ecriture des fichiers de donnees

  data_file_name[0] << prefix << ".dat";

  if (convol_histo) {
    pdist[0] = this;
    scale[0] = 1.;

    pdist[1] = this;
    phisto[0] = convol_histo;
    scale[1] = convol_histo->nb_element;

    for (i = 0;i < nb_distribution;i++) {
      pdist[i + 2] = distribution[i];
      phisto[i + 1] = convol_histo->frequency_distribution[i];
      scale[i + 2] = convol_histo->frequency_distribution[i]->nb_element;
    }

    status = stat_tool::plot_print((data_file_name[0].str()).c_str() , nb_distribution + 2 , pdist ,
                                   scale , NULL , nb_distribution + 1 , phisto);
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
               << MIN(plot_max * YSCALE , 1.) << "] ";
      for (j = 0;j < nb_distribution;j++) {
        out_file << "\"" << label((data_file_name[j + 1].str()).c_str()) << "\" title \""
                 << STAT_label[STATL_DISTRIBUTION] << " " << j + 1;
        distribution[j]->plot_title_print(out_file);
        out_file << "\" with linespoints,\\" << endl;
      }

      out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\"";
      if (convol_histo) {
        out_file << " using " << nb_distribution + 2;
      }
      out_file << " title \"" << STAT_label[STATL_CONVOLUTION] << "\" with linespoints" << endl;

      if (convol_histo) {
        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;

        out_file << "plot [0:" << nb_value - 1 << "] [0:"
                 << (int)(MAX(convol_histo->max , max * convol_histo->nb_element) * YSCALE) + 1
                 << "] \"" << label((data_file_name[0].str()).c_str()) << "\" using 1 title \""
                 << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses,\\" << endl;
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
                   << (int)(MAX(convol_histo->frequency_distribution[j]->max ,
                                distribution[j]->max * convol_histo->frequency_distribution[j]->nb_element) * YSCALE) + 1
                   << "] \"" << label((data_file_name[0].str()).c_str()) << "\" using " << j + 2
                   << " title \"" << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << j + 1
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
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Convolution::plot_write(StatError &error , const char *prefix ,
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
 *  argument : pointeur sur un objet ConvolutionData.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* Convolution::get_plotable(const ConvolutionData *convol_histo) const

{
  register int i;
  int xmax;
  double ymax;
  ostringstream title , legend;


  MultiPlotSet *plot_set = new MultiPlotSet(convol_histo ? nb_distribution + 3 : 1);
  MultiPlotSet &plot = *plot_set;

  title.str("");
  title << STAT_label[STATL_CONVOLUTION];
  if (convol_histo) {
    title << " " << STAT_label[STATL_FIT];
  }
  plot.title = title.str();

  plot.border = "15 lw 0";

  // 1ere vue : produit de convolution

  xmax = nb_value - 1;
  if ((cumul[xmax] > 1. - DOUBLE_ERROR) &&
      (mass[xmax] > PLOT_MASS_THRESHOLD)) {
    xmax++;
  }
  plot[0].xrange = Range(0 , xmax);

  ymax = distribution[0]->max;
  for (i = 1;i < nb_distribution;i++) {
    if (distribution[i]->max > ymax) {
      ymax = distribution[i]->max;
    }
  }
  plot[0].yrange = Range(0. , MIN(ymax * YSCALE , 1.));

  if (nb_value - 1 < TIC_THRESHOLD) {
    plot[0].xtics = 1;
  }

  plot[0].resize(nb_distribution + 1);

  for (i = 0;i < nb_distribution;i++) {
    legend.str("");
    legend << STAT_label[STATL_DISTRIBUTION] << " " << i + 1;
    distribution[i]->plot_title_print(legend);
    plot[0][i].legend = legend.str();

    plot[0][i].style = "linespoints";

    distribution[i]->plotable_mass_write(plot[0][i]);
  }

  plot[0][nb_distribution].legend = STAT_label[STATL_CONVOLUTION];

  plot[0][nb_distribution].style = "linespoints";

  plotable_mass_write(plot[0][nb_distribution]);

  if (convol_histo) {

    // 2eme vue : produit de convolution ajuste

    plot[1].xrange = Range(0 , xmax);
    plot[1].yrange = Range(0. , ceil(MAX(convol_histo->max ,
                                     max * convol_histo->nb_element) * YSCALE));

    if (nb_value - 1 < TIC_THRESHOLD) {
      plot[1].xtics = 1;
    }

    plot[1].resize(2);

    plot[1][0].legend = STAT_label[STATL_FREQUENCY_DISTRIBUTION];

    plot[1][0].style = "impulses";

    convol_histo->plotable_frequency_write(plot[1][0]);

    plot[1][1].legend = STAT_label[STATL_CONVOLUTION];

    plot[1][1].style = "linespoints";

    plotable_mass_write(plot[1][1] , convol_histo->nb_element);

    // 3eme vue : fonctions de repartition

    plot[2].xrange = Range(0 , xmax);
    plot[2].yrange = Range(0. , 1.);

    if (nb_value - 1 < TIC_THRESHOLD) {
      plot[2].xtics = 1;
    }

    plot[2].resize(2);

    legend.str("");
    legend << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
           << " " << STAT_label[STATL_FUNCTION];
    plot[2][0].legend = legend.str();

    plot[2][0].style = "linespoints";

    convol_histo->plotable_cumul_write(plot[2][0]);

    legend.str("");
    legend << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_CONVOLUTION]
           << " " << STAT_label[STATL_FUNCTION];
    plot[2][1].legend = legend.str();

    plot[2][1].style = "linespoints";

    plotable_cumul_write(plot[2][1]);

    for (i = 0;i < nb_distribution;i++) {

      // vues suivantes : distributions ajustees

      xmax = distribution[i]->nb_value - 1;
      if ((distribution[i]->cumul[xmax] > 1. - DOUBLE_ERROR) &&
          (distribution[i]->mass[xmax] > PLOT_MASS_THRESHOLD)) {
        xmax++;
      }
      plot[i + 3].xrange = Range(0 , xmax);

      plot[i + 3].yrange = Range(0. , ceil(MAX(convol_histo->frequency_distribution[i]->max ,
                                           distribution[i]->max * convol_histo->frequency_distribution[i]->nb_element)
                                           * YSCALE));

      if (distribution[i]->nb_value - 1 < TIC_THRESHOLD) {
        plot[i + 3].xtics = 1;
      }

      plot[i + 3].resize(2);

      legend.str("");
      legend << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i + 1;
      plot[i + 3][0].legend = legend.str();

      plot[i + 3][0].style = "impulses";

      convol_histo->frequency_distribution[i]->plotable_frequency_write(plot[i + 3][0]);

      legend.str("");
      legend << STAT_label[STATL_DISTRIBUTION] << " " << i + 1;
      distribution[i]->plot_title_print(legend);
      plot[i + 3][1].legend = legend.str();

      plot[i + 3][1].style = "linespoints";

      distribution[i]->plotable_mass_write(plot[i + 3][1] , convol_histo->frequency_distribution[i]->nb_element);
    }
  }

  return plot_set;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet Convolution.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* Convolution::get_plotable() const

{
  return get_plotable(convolution_data);
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe ConvolutionData.
 *
 *--------------------------------------------------------------*/

ConvolutionData::ConvolutionData()

{
  convolution = NULL;
  nb_distribution = 0;
  frequency_distribution = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe ConvolutionData.
 *
 *  arguments : reference sur un objet FrequencyDistribution,
 *              nombre de lois empiriques.
 *
 *--------------------------------------------------------------*/

ConvolutionData::ConvolutionData(const FrequencyDistribution &histo , int nb_dist)
:FrequencyDistribution(histo)

{
  register int i;


  convolution = NULL;
  nb_distribution = nb_dist;

  frequency_distribution = new FrequencyDistribution*[nb_distribution];
  for (i = 0;i < nb_distribution;i++) {
    frequency_distribution[i] = new FrequencyDistribution(nb_value);
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe ConvolutionData.
 *
 *  argument : reference sur un objet Convolution.
 *
 *--------------------------------------------------------------*/

ConvolutionData::ConvolutionData(const Convolution &convol)
:FrequencyDistribution(convol)

{
  register int i;


  convolution = NULL;
  nb_distribution = convol.nb_distribution;

  frequency_distribution = new FrequencyDistribution*[nb_distribution];
  for (i = 0;i < nb_distribution;i++) {
    frequency_distribution[i] = new FrequencyDistribution(*(convol.distribution[i]));
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet ConvolutionData.
 *
 *  arguments : reference sur un objet ConvolutionData,
 *              flag copie de l'objet Convolution.
 *
 *--------------------------------------------------------------*/

void ConvolutionData::copy(const ConvolutionData &convol_histo , bool model_flag)

{
  register int i;


  if ((model_flag) && (convol_histo.convolution)) {
    convolution = new Convolution(*(convol_histo.convolution) , false);
  }
  else {
    convolution = NULL;
  }

  nb_distribution = convol_histo.nb_distribution;

  frequency_distribution = new FrequencyDistribution*[nb_distribution];
  for (i = 0;i < nb_distribution;i++) {
    frequency_distribution[i] = new FrequencyDistribution(*(convol_histo.frequency_distribution[i]));
  }
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet ConvolutionData.
 *
 *--------------------------------------------------------------*/

void ConvolutionData::remove()

{
  delete convolution;

  if (frequency_distribution) {
    register int i;

    for (i = 0;i < nb_distribution;i++) {
      delete frequency_distribution[i];
    }
    delete [] frequency_distribution;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe ConvolutionData.
 *
 *--------------------------------------------------------------*/

ConvolutionData::~ConvolutionData()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe ConvolutionData.
 *
 *  argument : reference sur un objet ConvolutionData.
 *
 *--------------------------------------------------------------*/

ConvolutionData& ConvolutionData::operator=(const ConvolutionData &convol_histo)

{
  if (&convol_histo != this) {
    remove();
    delete [] frequency;

    FrequencyDistribution::copy(convol_histo);
    copy(convol_histo);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Extraction de lois empiriques elementaires.
 *
 *  arguments : reference sur un objet StatError, indice de la loi empirique.
 *
 *--------------------------------------------------------------*/

DiscreteDistributionData* ConvolutionData::extract(StatError &error , int index) const

{
  DiscreteDistributionData *phisto;


  error.init();

  if ((index < 1) || (index > nb_distribution)) {
    phisto = NULL;
    error.update(STAT_error[STATR_FREQUENCY_DISTRIBUTION_INDEX]);
  }

  else {
    index--;
    phisto = new DiscreteDistributionData(*frequency_distribution[index] ,
                                          (convolution ? convolution->distribution[index] : NULL));
  }

  return phisto;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet ConvolutionData.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& ConvolutionData::line_write(ostream &os) const

{
  os << STAT_label[STATL_SAMPLE_SIZE] << ": " << nb_element << "   "
     << STAT_label[STATL_MEAN] << ": " << mean << "   "
     << STAT_label[STATL_VARIANCE] << ": " << variance;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet ConvolutionData.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& ConvolutionData::ascii_write(ostream &os , bool exhaustive) const

{
  if (convolution) {
    convolution->ascii_write(os , this , exhaustive , false);
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet ConvolutionData dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool ConvolutionData::ascii_write(StatError &error , const char *path ,
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
 *  Ecriture d'un objet ConvolutionData dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool ConvolutionData::spreadsheet_write(StatError &error , const char *path) const

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
 *  Sortie Gnuplot d'un objet ConvolutionData.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool ConvolutionData::plot_write(StatError &error , const char *prefix ,
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
 *  Sortie graphique d'un objet ConvolutionData.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* ConvolutionData::get_plotable() const

{
  MultiPlotSet *plot_set;


  if (convolution) {
    plot_set = convolution->get_plotable(this);
  }
  else {
    plot_set = NULL;
  }

  return plot_set;
}


};  // namespace stat_tool
