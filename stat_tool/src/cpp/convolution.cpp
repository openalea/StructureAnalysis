/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2017 CIRAD/INRA/Inria Virtual Plants
 *
 *       File author(s): Yann Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id$
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
#include <string>
#include <vector>

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "convolution.h"
#include "stat_label.h"

using namespace std;
using namespace boost;


namespace stat_tool {



/*--------------------------------------------------------------*/
/**
 *  \brief Default constructor of the Convolution class.
 */
/*--------------------------------------------------------------*/

Convolution::Convolution()

{
  convolution_data = NULL;
  nb_distribution = 0;
  distribution = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the Convolution class.
 *
 *  \param[in] nb_dist number of distributions,
 *  \param[in] pdist   pointer on the elementary distributions.
 */
/*--------------------------------------------------------------*/

Convolution::Convolution(int nb_dist , const DiscreteParametric **pdist)

{
  int i;
  int cnb_value = 1;


  convolution_data = NULL;
  nb_distribution = nb_dist;

  distribution = new DiscreteParametric*[nb_distribution];
  for (i = 0;i < nb_distribution;i++) {
    distribution[i] = new DiscreteParametric(*pdist[i] , NORMALIZATION);
    cnb_value += distribution[i]->nb_value - 1;
  }

  Distribution::init(cnb_value);
  computation(1 , CONVOLUTION_THRESHOLD);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the Convolution class.
 *
 *  \param[in] nb_dist number of distributions,
 *  \param[in] pdist   pointer on the elementary distributions.
 */
/*--------------------------------------------------------------*/

Convolution::Convolution(int nb_dist , const vector<DiscreteParametric> idist)

{
  int i;
  int cnb_value = 1;


  convolution_data = NULL;
  nb_distribution = nb_dist;

  distribution = new DiscreteParametric*[nb_distribution];
  for (i = 0;i < nb_distribution;i++) {
    distribution[i] = new DiscreteParametric(idist[i] , NORMALIZATION);
    cnb_value += distribution[i]->nb_value - 1;
  }

  Distribution::init(cnb_value);
  computation(1 , CONVOLUTION_THRESHOLD);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the Convolution class.
 *
 *  \param[in] known_dist   reference on the known distribution,
 *  \param[in] unknown_dist reference on the unknown distribution.
 */
/*--------------------------------------------------------------*/

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
    distribution[0] = new DiscreteParametric(known_dist , NORMALIZATION);
  }
  distribution[1] = new DiscreteParametric(unknown_dist , DISTRIBUTION_COPY ,
                                           (int)(unknown_dist.nb_value * NB_VALUE_COEFF));

  Distribution::init(distribution[0]->alloc_nb_value + distribution[1]->alloc_nb_value - 1);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a Convolution object.
 *
 *  \param[in] convol    reference on a Convolution object,
 *  \param[in] data_flag flag copy of the included ConvolutionData object.
 */
/*--------------------------------------------------------------*/

void Convolution::copy(const Convolution &convol , bool data_flag)

{
  int i;


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


/*--------------------------------------------------------------*/
/**
 *  \brief Destruction of the data members of a Convolution object.
 */
/*--------------------------------------------------------------*/

void Convolution::remove()

{
  delete convolution_data;

  if (distribution) {
    int i;

    for (i = 0;i < nb_distribution;i++) {
      delete distribution[i];
    }
    delete [] distribution;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the Convolution class.
 */
/*--------------------------------------------------------------*/

Convolution::~Convolution()

{
  remove();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the Convolution class.
 *
 *  \param[in] convol reference on a Convolution object.
 *
 *  \return           Convolution object.
 */
/*--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of an elementary distribution.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] index distribution index.
 *
 *  \return          DiscreteParametricModel object.
 */
/*--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of the data part of a Convolution object.
 *
 *  \param[in] error reference on a StatError object.
 *
 *  \return          ConvolutionData object.
 */
/*--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a Convolution object on the basis of
 *         elementary distributions.
 *
 *  \param[in] error   reference on a StatError object,
 *  \param[in] nb_dist number of distributions,
 *  \param[in] dist    pointer on the distributions.
 *
 *  \return            Convolution object.
 */
/*--------------------------------------------------------------*/

Convolution* Convolution::building(StatError &error , int nb_dist , const DiscreteParametric **dist)

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


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a Convolution object on the basis of
 *         elementary distributions.
 *
 *  \param[in] error   reference on a StatError object,
 *  \param[in] nb_dist number of distributions,
 *  \param[in] dist    pointer on the distributions.
 *
 *  \return            Convolution object.
 */
/*--------------------------------------------------------------*/

Convolution* Convolution::building(StatError &error , int nb_dist , const vector<DiscreteParametric> dist)

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


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a Convolution object from a file.
 *
 *  \param[in]  error           reference on a StatError object,
 *  \param[in]  path            file path,
 *  \param[in]  cumul_threshold threshold on the cumulative distribution function.
 *
 *  \return                     Convolution object.
 */
/*--------------------------------------------------------------*/

Convolution* Convolution::ascii_read(StatError &error , const string path , double cumul_threshold)

{
  string buffer;
  size_t position;
  typedef tokenizer<char_separator<char>> tokenizer;
  char_separator<char> separator(" \t");
  bool status , lstatus;
  int i , j;
  int line , nb_dist , index;
  const DiscreteParametric **dist;
  Convolution *convol;
  ifstream in_file(path.c_str());


  convol = NULL;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    line = 0;
    nb_dist = 0;

    while (getline(in_file , buffer)) {
      line++;

#     ifdef DEBUG
      cout << line << "  " << buffer << endl;
#     endif

      position = buffer.find('#');
      if (position != string::npos) {
        buffer.erase(position);
      }
      i = 0;

      tokenizer tok_buffer(buffer , separator);

      for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
        switch (i) {

        // test CONVOLUTION keyword

        case 0 : {
          if (*token != STAT_word[STATW_CONVOLUTION]) {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEYWORD] , STAT_word[STATW_CONVOLUTION] , line , i + 1);
          }
          break;
        }

        // test number of distributions

        case 1 : {
          lstatus = true;

/*          try {
            nb_dist = stoi(*token);   in C++ 11
          }
          catch(invalid_argument &arg) {
            lstatus = false;
          } */
          nb_dist = atoi(token->c_str());

          if ((lstatus) && ((nb_dist < 2) || (nb_dist > CONVOLUTION_NB_DISTRIBUTION))) {
            lstatus = false;
          }

          if (!lstatus) {
            status = false;
            error.update(STAT_parsing[STATP_NB_DISTRIBUTION] , line , i + 1);
          }
          break;
        }

        // test DISTRIBUTIONS keyword

        case 2 : {
          if (*token != STAT_word[STATW_DISTRIBUTIONS]) {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEYWORD] , STAT_word[STATW_DISTRIBUTIONS] , line , i + 1);
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
        while (getline(in_file , buffer)) {
          line++;

#         ifdef DEBUG
          cout << line << "  " << buffer << endl;
#         endif

          position = buffer.find('#');
          if (position != string::npos) {
            buffer.erase(position);
          }
          j = 0;

          tokenizer tok_buffer(buffer , separator);

          for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
            switch (j) {

            // test DISTRIBUTION keyword

            case 0 : {
              if (*token != STAT_word[STATW_DISTRIBUTION]) {
                status = false;
                error.correction_update(STAT_parsing[STATP_KEYWORD] , STAT_word[STATW_DISTRIBUTION] , line , j + 1);
              }
              break;
            }

            // test distribution index

            case 1 : {
              lstatus = true;

/*              try {
                index = stoi(*token);   in C++ 11
              }
              catch(invalid_argument &arg) {
                lstatus = false;
              } */
              index = atoi(token->c_str());

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

            dist[i] = DiscreteParametric::parsing(error , in_file , line ,
                                                  NEGATIVE_BINOMIAL , cumul_threshold);
            break;
          }
        }

        if (!dist[i]) {
          status = false;
          error.update(STAT_parsing[STATP_FORMAT] , line);
        }
      }

      while (getline(in_file , buffer)) {
        line++;

#       ifdef DEBUG
        cout << line << "  " << buffer << endl;
#       endif

        position = buffer.find('#');
        if (position != string::npos) {
          buffer.erase(position);
        }
        if (!(trim_right_copy_if(buffer , is_any_of(" \t")).empty())) {
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


/*--------------------------------------------------------------*/
/**
 *  \brief Writing on a single line of a Convolution object.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& Convolution::line_write(ostream &os) const

{
  os << nb_distribution << " " << STAT_word[STATW_DISTRIBUTIONS] << "   "
     << STAT_label[STATL_MEAN] << ": " << mean << "   "
     << STAT_label[STATL_VARIANCE] << ": " << variance;

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a convolution of distributions and the associated
 *         frequency distributions in a file.
 *
 *  \param[in,out] os           stream,
 *  \param[in]     convol_histo pointer on a ConvolutionData object,
 *  \param[in]     exhaustive   flag detail level,
 *  \param[in]     file_flag    flag file output.
 */
/*--------------------------------------------------------------*/

ostream& Convolution::ascii_write(ostream &os , const ConvolutionData *convol_histo ,
                                  bool exhaustive , bool file_flag) const

{
  int i;
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


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a Convolution object.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     exhaustive flag detail level.
 */
/*--------------------------------------------------------------*/

ostream& Convolution::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , convolution_data , exhaustive , false);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a Convolution object in a file.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] path       file path,
 *  \param[in] exhaustive flag detail level.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool Convolution::ascii_write(StatError &error , const string path ,
                              bool exhaustive) const

{
  bool status;
  ofstream out_file(path.c_str());


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


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a convolution of distributions and the associated
 *         frequency distributions in a file at the spreadsheet format.
 *
 *  \param[in,out] os           stream,
 *  \param[in]     convol_histo pointer on a ConvolutionData object.
 */
/*--------------------------------------------------------------*/

ostream& Convolution::spreadsheet_write(ostream &os , const ConvolutionData *convol_histo) const

{
  int i;
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


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a Convolution object in a file at the spreadsheet format.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 *
 *  \return          error status.
 */
/*--------------------------------------------------------------*/

bool Convolution::spreadsheet_write(StatError &error , const string path) const

{
  bool status;
  ofstream out_file(path.c_str());


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


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a convolution of distributions and the associated
 *         frequency distributions using Gnuplot.
 *
 *  \param[in] prefix       file prefix,
 *  \param[in] title        figure title,
 *  \param[in] convol_histo pointer on a ConvolutionData object.
 *
 *  \return                 error status.
 */
/*--------------------------------------------------------------*/

bool Convolution::plot_write(const char *prefix , const char *title ,
                             const ConvolutionData *convol_histo) const

{
  bool status;
  int i , j;
  double plot_max , scale[CONVOLUTION_NB_DISTRIBUTION + 2];
  const Distribution *pdist[CONVOLUTION_NB_DISTRIBUTION + 2];
  const FrequencyDistribution *phisto[CONVOLUTION_NB_DISTRIBUTION + 1];
  ostringstream data_file_name[CONVOLUTION_NB_DISTRIBUTION + 1];


  // writing of the data files

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

    // writing of the script files

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


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a Convolution object using Gnuplot.
 *
 *  \param[in] error  reference on a StatError object,
 *  \param[in] prefix file prefix,
 *  \param[in] title  figure title.
 *
 *  \return           error status.
 */
/*--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a convolution of distributions and the associated
 *         frequency distributions.
 *
 *  \param[in] convol_histo pointer on a ConvolutionData object.
 *
 *  \return                 MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* Convolution::get_plotable(const ConvolutionData *convol_histo) const

{
  int i;
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

  // convolution of distributions

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

    // fit of convolution

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

    // cumulative distribution functions

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

      // fit of elementary distributions

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


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a Convolution object.
 *
 *  \return MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* Convolution::get_plotable() const

{
  return get_plotable(convolution_data);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Default constructor of the ConvolutionData class.
 */
/*--------------------------------------------------------------*/

ConvolutionData::ConvolutionData()

{
  convolution = NULL;
  nb_distribution = 0;
  frequency_distribution = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the ConvolutionData class.
 *
 *  \param[in] histo   reference on a FrequencyDistribution object,
 *  \param[in] nb_dist number of frequency distributions.
 */
/*--------------------------------------------------------------*/

ConvolutionData::ConvolutionData(const FrequencyDistribution &histo , int nb_dist)
:FrequencyDistribution(histo)

{
  int i;


  convolution = NULL;
  nb_distribution = nb_dist;

  frequency_distribution = new FrequencyDistribution*[nb_distribution];
  for (i = 0;i < nb_distribution;i++) {
    frequency_distribution[i] = new FrequencyDistribution(nb_value);
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the ConvolutionData class.
 *
 *  \param[in] convol reference on a Convolution object.
 */
/*--------------------------------------------------------------*/

ConvolutionData::ConvolutionData(const Convolution &convol)
:FrequencyDistribution(convol)

{
  int i;


  convolution = NULL;
  nb_distribution = convol.nb_distribution;

  frequency_distribution = new FrequencyDistribution*[nb_distribution];
  for (i = 0;i < nb_distribution;i++) {
    frequency_distribution[i] = new FrequencyDistribution(*(convol.distribution[i]));
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a ConvolutionData object.
 *
 *  \param[in] convol_histo reference on a ConvolutionData object,
 *  \param[in] model_flag   flag copy of the included Convolution object.
 */
/*--------------------------------------------------------------*/

void ConvolutionData::copy(const ConvolutionData &convol_histo , bool model_flag)

{
  int i;


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


/*--------------------------------------------------------------*/
/**
 *  \brief Destruction of the data members of a ConvolutionData object.
 */
/*--------------------------------------------------------------*/

void ConvolutionData::remove()

{
  delete convolution;

  if (frequency_distribution) {
    int i;

    for (i = 0;i < nb_distribution;i++) {
      delete frequency_distribution[i];
    }
    delete [] frequency_distribution;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the ConvolutionData class.
 */
/*--------------------------------------------------------------*/

ConvolutionData::~ConvolutionData()

{
  remove();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the ConvolutionData class.
 *
 *  \param[in] convol_histo reference on a ConvolutionData object.
 *
 *  \return                 ConvolutionData object.
 */
/*--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of an elementary frequency distribution.
 *
 *  \param[in] error  reference on a StatError object,
 *  \param[in] index  frequency distribution index.
 *
 *  \return           DiscreteDistributionData object.
 */
/*--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*/
/**
 *  \brief Writing on a single line of a ConvolutionData object.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& ConvolutionData::line_write(ostream &os) const

{
  os << STAT_label[STATL_SAMPLE_SIZE] << ": " << nb_element << "   "
     << STAT_label[STATL_MEAN] << ": " << mean << "   "
     << STAT_label[STATL_VARIANCE] << ": " << variance;

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a ConvolutionData object.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     exhaustive flag detail level.
 */
/*--------------------------------------------------------------*/

ostream& ConvolutionData::ascii_write(ostream &os , bool exhaustive) const

{
  if (convolution) {
    convolution->ascii_write(os , this , exhaustive , false);
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a ConvolutionData object in a file.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] path       file path,
 *  \param[in] exhaustive flag detail level.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool ConvolutionData::ascii_write(StatError &error , const string path ,
                                  bool exhaustive) const

{
  bool status = false;


  if (convolution) {
    ofstream out_file(path.c_str());

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


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a ConvolutionData object in a file at the spreadsheet format.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 *
 *  \return          error status.
 */
/*--------------------------------------------------------------*/

bool ConvolutionData::spreadsheet_write(StatError &error , const string path) const

{
  bool status = false;


  if (convolution) {
    ofstream out_file(path.c_str());

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


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a ConvolutionData object using Gnuplot.
 *
 *  \param[in] error  reference on a StatError object,
 *  \param[in] prefix file prefix,
 *  \param[in] title  figure title.
 *
 *  \return           error status.
 */
/*--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a ConvolutionData object.
 *
 *  \return MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

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
