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
 *       $Id: discrete_mixture.cpp 17990 2015-04-23 06:45:46Z guedon $
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



#include <sstream>
#include <iomanip>

#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"
#include "tool/config.h"

#include "stat_tools.h"
#include "distribution.h"
#include "discrete_mixture.h"
#include "stat_label.h"

using namespace std;


namespace stat_tool {



/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe DiscreteMixture.
 *
 *--------------------------------------------------------------*/

DiscreteMixture::DiscreteMixture()

{
  mixture_data = NULL;
  nb_component = 0;
  weight = NULL;
  component = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe DiscreteMixture.
 *
 *  arguments : nombre de composantes, poids, pointeurs sur les composantes.
 *
 *--------------------------------------------------------------*/

DiscreteMixture::DiscreteMixture(int inb_component , double *pweight ,
                                 const DiscreteParametric **pcomponent)

{
  register int i;


  mixture_data = NULL;
  nb_component = inb_component;

  weight = new DiscreteParametric(nb_component);
  for (i = 0;i < nb_component;i++) {
    weight->mass[i] = *pweight++;
  }
  weight->cumul_computation();
  weight->max_computation();

  nb_value = 0;
  component = new DiscreteParametric*[nb_component];

  for (i = 0;i < nb_component;i++) {
    component[i] = new DiscreteParametric(*pcomponent[i] , 'n');
    if (component[i]->nb_value > nb_value) {
      nb_value = component[i]->nb_value;
    }
  }

  Distribution::init(nb_value);

  computation(1 , CUMUL_THRESHOLD , false);
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe DiscreteMixture.
 *
 *  arguments : reference sur un objet DiscreteMixture, flag sur les composantes a copier,
 *              nombre de valeurs.
 *
 *--------------------------------------------------------------*/

DiscreteMixture::DiscreteMixture(const DiscreteMixture &mixt , bool *component_flag , int inb_value)

{
  register int i;


  mixture_data = NULL;
  nb_component = mixt.nb_component;

  weight = new DiscreteParametric(nb_component);

  nb_value = inb_value;
  component = new DiscreteParametric*[nb_component];

  for (i = 0;i < nb_component;i++) {
    if (component_flag[i]) {
      component[i] = new DiscreteParametric(inb_value , mixt.component[i]->ident);
    }
    else {
      component[i] = new DiscreteParametric(*(mixt.component[i]));
      if (component[i]->nb_value > nb_value) {
        nb_value = component[i]->nb_value;
      }
    }
  }

  Distribution::init(nb_value);
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe DiscreteMixture.
 *
 *  arguments : nombre de composantes, pointeurs sur les composantes.
 *
 *--------------------------------------------------------------*/

DiscreteMixture::DiscreteMixture(int inb_component , const DiscreteParametric **pcomponent)

{
  register int i;


  mixture_data = NULL;
  nb_component = inb_component;

  weight = NULL;

  component = new DiscreteParametric*[nb_component];
  for (i = 0;i < nb_component;i++) {
    component[i] = new DiscreteParametric(*pcomponent[i] , 'n');
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet DiscreteMixture.
 *
 *  arguments : reference sur un objet DiscreteMixture,
 *              flag copie de l'objet DiscreteMixtureData.
 *
 *--------------------------------------------------------------*/

void DiscreteMixture::copy(const DiscreteMixture &mixt , bool data_flag)

{
  register int i;


  if ((data_flag) && (mixt.mixture_data)) {
    mixture_data = new DiscreteMixtureData(*(mixt.mixture_data) , false);
  }
  else {
    mixture_data = NULL;
  }

  nb_component = mixt.nb_component;

  weight = new DiscreteParametric(*(mixt.weight));

  component = new DiscreteParametric*[nb_component];
  for (i = 0;i < nb_component;i++) {
    component[i] = new DiscreteParametric(*(mixt.component[i]));
  }
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet DiscreteMixture.
 *
 *--------------------------------------------------------------*/

void DiscreteMixture::remove()

{
  delete mixture_data;

  delete weight;

  if (component) {
    register int i;

    for (i = 0;i < nb_component;i++) {
      delete component[i];
    }
    delete [] component;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe DiscreteMixture.
 *
 *--------------------------------------------------------------*/

DiscreteMixture::~DiscreteMixture()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe DiscreteMixture.
 *
 *  argument : reference sur un objet DiscreteMixture.
 *
 *--------------------------------------------------------------*/

DiscreteMixture& DiscreteMixture::operator=(const DiscreteMixture &mixt)

{
  if (&mixt != this) {
    remove();
    delete [] mass;
    delete [] cumul;

    Distribution::copy(mixt);
    copy(mixt);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Extraction d'une composante.
 *
 *  arguments : reference sur un objet StatError, indice de la composante.
 *
 *--------------------------------------------------------------*/

DiscreteParametricModel* DiscreteMixture::extract(StatError &error , int index) const

{
  DiscreteParametricModel *pcomponent;


  if ((index < 1) || (index > nb_component)) {
    pcomponent = NULL;
    error.update(STAT_error[STATR_DISTRIBUTION_INDEX]);
  }

  else {
    index--;
    pcomponent = new DiscreteParametricModel(*component[index] ,
                                             (mixture_data ? mixture_data->component[index] : NULL));
  }

  return pcomponent;
}


/*--------------------------------------------------------------*
 *
 *  Extraction de la partie "donnees" d'un objet DiscreteMixture.
 *
 *  argument : reference sur un objet StatError.
 *
 *--------------------------------------------------------------*/

DiscreteMixtureData* DiscreteMixture::extract_data(StatError &error) const

{
  DiscreteMixtureData *mixt_histo;


  error.init();

  if (!mixture_data) {
    mixt_histo = NULL;
    error.update(STAT_error[STATR_NO_DATA]);
  }

  else {
    mixt_histo = new DiscreteMixtureData(*mixture_data);
    mixt_histo->mixture = new DiscreteMixture(*this , false);
  }

  return mixt_histo;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet DiscreteMixture a partir de poids et de composantes.
 *
 *  arguments : reference sur un objet StatError, nombre de composantes,
 *              poids, pointeurs sur les composantes.
 *
 *--------------------------------------------------------------*/

DiscreteMixture* discrete_mixture_building(StatError &error , int nb_component , double *weight ,
                                           const DiscreteParametric **component)

{
  bool status;
  register int i;
  double cumul;
  DiscreteMixture *mixt;


  mixt = NULL;
  error.init();

  if ((nb_component < 2) || (nb_component > DISCRETE_MIXTURE_NB_COMPONENT)) {
    error.update(STAT_parsing[STATP_NB_DISTRIBUTION]);
  }

  else {
    status = true;
    cumul = 0.;

    for (i = 0;i < nb_component;i++) {
      if ((weight[i] <= 0.) || (weight[i] > 1. - cumul + DOUBLE_ERROR)) {
        status = false;
        error.update(STAT_parsing[STATP_WEIGHT_VALUE]);
      }
      else {
        cumul += weight[i];
      }
    }

    if (status) {
      mixt = new DiscreteMixture(nb_component , weight , component);
    }
  }

  return mixt;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet DiscreteMixture a partir d'un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

DiscreteMixture* discrete_mixture_ascii_read(StatError &error , const char *path ,
                                             double cumul_threshold)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  bool status , lstatus;
  register int i , j;
  int line;
  long index , nb_component;
  double cumul , weight[DISCRETE_MIXTURE_NB_COMPONENT];
  const DiscreteParametric **component;
  DiscreteMixture *mixt;
  ifstream in_file(path);


  mixt = NULL;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    line = 0;
    nb_component = 0;

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

        // test mot cle MIXTURE

        case 0 : {
          if (token != STAT_word[STATW_MIXTURE]) {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_MIXTURE] , line , i + 1);
          }
          break;
        }

        // test nombre de composantes

        case 1 : {
          lstatus = locale.stringToNum(token , &nb_component);
          if ((lstatus) && ((nb_component < 2) || (nb_component > DISCRETE_MIXTURE_NB_COMPONENT))) {
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

    if (nb_component == 0) {
      status = false;
      error.update(STAT_parsing[STATP_FORMAT] , line);
    }

    if (status) {
      component = new const DiscreteParametric*[nb_component];
      for (i = 0;i < nb_component;i++) {
        component[i] = NULL;
      }

      cumul = 0.;

      for (i = 0;i < nb_component;i++) {
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

            // test indice de la composante

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

            // test nom du parametre (WEIGHT)

            case 2 : {
              if (token != STAT_word[STATW_WEIGHT]) {
                status = false;
                error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_WEIGHT] , line , j + 1);
              }
              break;
            }

            // test separateur

            case 3 : {
              if (token != ":") {
                status = false;
                error.update(STAT_parsing[STATP_SEPARATOR] , line , j + 1);
              }
              break;
            }

            // test valeur du poids

            case 4 : {
              lstatus = locale.stringToNum(token , weight + i);
              if (lstatus) {
                if ((weight[i] <= 0.) || (weight[i] > 1. - cumul + DOUBLE_ERROR)) {
                  lstatus = false;
                }
                else {
                  cumul += weight[i];
                }
              }

              if (!lstatus) {
                status = false;
                error.update(STAT_parsing[STATP_WEIGHT_VALUE] , line , j + 1);
              }
              break;
            }
            }

            j++;
          }

          if (j > 0) {
            if (j != 5) {
              status = false;
              error.update(STAT_parsing[STATP_FORMAT] , line);
            }

            component[i] = discrete_parametric_parsing(error , in_file , line ,
                                                       NEGATIVE_BINOMIAL , cumul_threshold);
            break;
          }
        }

        if (!component[i]) {
          status = false;
          error.update(STAT_parsing[STATP_FORMAT] , line);
        }
      }

      if (cumul < CUMUL_THRESHOLD) {
        status = false;
        error.update(STAT_parsing[STATP_PROBABILITY_SUM]);
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
        mixt = new DiscreteMixture(nb_component , weight , component);
      }

      for (i = 0;i < nb_component;i++) {
        delete component[i];
      }
      delete [] component;
    }
  }

  return mixt;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet DiscreteMixture.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& DiscreteMixture::line_write(ostream &os) const

{
  os << nb_component << " " << STAT_word[STATW_DISTRIBUTIONS] << "   "
     << STAT_label[STATL_MEAN] << ": " << mean << "   "
     << STAT_label[STATL_VARIANCE] << ": " << variance;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un melange de lois discretes et de la structure de donnees associee.
 *
 *  arguments : stream, pointeur sur un objet DiscreteMixtureData,
 *              flag niveau de detail, flag fichier.
 *
 *--------------------------------------------------------------*/

ostream& DiscreteMixture::ascii_write(ostream &os , const DiscreteMixtureData *mixt_histo ,
                                      bool exhaustive , bool file_flag) const

{
  register int i , j;
  int bnb_parameter , buff , width;
  double **sup_norm_dist , scale[DISCRETE_MIXTURE_NB_COMPONENT];
  const Distribution *pcomponent[DISCRETE_MIXTURE_NB_COMPONENT];
  long old_adjust;


  os << STAT_word[STATW_MIXTURE] << " " << nb_component << " " << STAT_word[STATW_DISTRIBUTIONS] << endl;
  ascii_characteristic_print(os , false , file_flag);

  if (mixt_histo) {
    double likelihood , information;
    Test test(CHI2);


    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    mixt_histo->ascii_characteristic_print(os , false , file_flag);

    likelihood = Distribution::likelihood_computation(*mixt_histo);
    information = mixt_histo->FrequencyDistribution::information_computation();

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   ("
       << STAT_label[STATL_NORMALIZED] << ": " << likelihood / mixt_histo->nb_element << ")" << endl;

    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_MAX_LIKELIHOOD] << ": " << information << "   ("
       << STAT_label[STATL_INFORMATION] << ": " << information / mixt_histo->nb_element << ")" << endl;

    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_DEVIANCE] << ": " << 2 * (information - likelihood) << endl;

    bnb_parameter = nb_parameter_computation();

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << bnb_parameter << " " << STAT_label[STATL_FREE_PARAMETERS] << "   2 * "
       << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " (" << STAT_criterion_word[AIC] << "): "
       << 2 * (likelihood - bnb_parameter) << endl;

    if (0 < bnb_parameter < mixt_histo->nb_element - 1) {
      if (file_flag) {
        os << "# ";
      }
      os << bnb_parameter << " " << STAT_label[STATL_FREE_PARAMETERS] << "   2 * "
         << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " (" << STAT_criterion_word[AICc] << "): "
         << 2 * (likelihood - (double)(bnb_parameter * mixt_histo->nb_element) /
            (double)(mixt_histo->nb_element - bnb_parameter - 1)) << endl;
    }

    if (file_flag) {
      os << "# ";
    }
    os << bnb_parameter << " " << STAT_label[STATL_FREE_PARAMETERS] << "   2 * "
       << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " (" << STAT_criterion_word[BIC] << "): "
       << 2 * likelihood - bnb_parameter * log((double)mixt_histo->nb_element) << endl;

    if (file_flag) {
      os << "# ";
    }
    os << bnb_parameter << " " << STAT_label[STATL_FREE_PARAMETERS] << "   2 * "
       << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " (" << STAT_criterion_word[BICc] << "): "
       << 2 * likelihood - penalty_computation() << endl;

    likelihood = likelihood_computation(*mixt_histo);
    information = mixt_histo->information_computation();

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_CLASSIFICATION_LIKELIHOOD] << ": " << likelihood << "   ("
       << STAT_label[STATL_NORMALIZED] << ": " << likelihood / mixt_histo->nb_element << ")" << endl;

    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_MAX_CLASSIFICATION_LIKELIHOOD] << ": " << information << "   ("
       << STAT_label[STATL_INFORMATION] << ": " << information / mixt_histo->nb_element << ")" << endl;

    chi2_fit(*mixt_histo , test);
    os << "\n";
    test.ascii_print(os , file_flag);
  }

  else {
    for (i = 0;i < nb_component;i++) {
      os << "\n" << STAT_word[STATW_DISTRIBUTION] << " " << i + 1 << " "
         << STAT_word[STATW_WEIGHT] << " : "  << weight->mass[i] << endl;

      component[i]->ascii_print(os);
      component[i]->ascii_parametric_characteristic_print(os , exhaustive , file_flag);
    }
  }

  if (exhaustive) {
    for (i = 0;i < nb_component;i++) {
      pcomponent[i] = component[i];

      if (mixt_histo) {
        scale[i] = weight->mass[i] * mixt_histo->nb_element;
      }
      else {
        scale[i] = weight->mass[i];
      }
    }

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << "  ";
    if (mixt_histo) {
      os << " | " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    }
    for (i = 0;i < nb_component;i++) {
      os << " | " << STAT_label[STATL_DISTRIBUTION] << " " << i + 1;
    }
    os << " | " << STAT_label[STATL_MIXTURE];
    if (mixt_histo) {
      os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
         << " " << STAT_label[STATL_FUNCTION];
    }
    os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_MIXTURE]
       << " " << STAT_label[STATL_FUNCTION] << endl;

    ascii_print(os , nb_component , pcomponent , scale , file_flag , true , mixt_histo);
  }

  if (mixt_histo) {
    if (exhaustive) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
      mixt_histo->weight->ascii_characteristic_print(os , false , file_flag);

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << "   | " << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
         << " | " << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_DISTRIBUTION] << endl;

      weight->Distribution::ascii_print(os , file_flag , false , false , mixt_histo->weight);
    }

    for (i = 0;i < nb_component;i++) {
      os << "\n" << STAT_word[STATW_DISTRIBUTION] << " " << i + 1 << "  "
         << STAT_word[STATW_WEIGHT] << " : "  << weight->mass[i] << endl;
      component[i]->ascii_print(os);
      component[i]->ascii_parametric_characteristic_print(os , exhaustive , file_flag);

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i + 1 << " - ";
      mixt_histo->component[i]->ascii_characteristic_print(os , exhaustive , file_flag);

      if ((exhaustive) && (mixt_histo->component[i]->nb_element > 0)) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << "   | " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i + 1
           << " | " << STAT_label[STATL_DISTRIBUTION] << " " << i + 1
           << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
           << " " << i + 1 << " " << STAT_label[STATL_FUNCTION]
           << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION]
           << " " << i + 1 << " " << STAT_label[STATL_FUNCTION] << endl;

        component[i]->Distribution::ascii_print(os , file_flag , true , false ,
                                                mixt_histo->component[i]);
      }
    }
  }

  old_adjust = os.setf(ios::left , ios::adjustfield);

  sup_norm_dist = new double*[nb_component];
  for (i = 0;i < nb_component;i++) {
    sup_norm_dist[i] = new double[nb_component];
  }

  for (i = 0;i < nb_component;i++) {
    sup_norm_dist[i][i] = 0.;
    for (j = i + 1;j < nb_component;j++) {
      sup_norm_dist[i][j] = component[i]->sup_norm_distance_computation(*(component[j]));
      sup_norm_dist[j][i] = sup_norm_dist[i][j];
    }
  }

  width = column_width(nb_component , sup_norm_dist[0]);
  for (i = 1;i < nb_component;i++) {
    buff = column_width(nb_component , sup_norm_dist[i]);
    if (buff > width) {
      width = buff;
    }
  }
  width += ASCII_SPACE;

  os << "\n";
  if (file_flag) {
    os << "# ";
  }
  os << STAT_label[STATL_COMPONENT_DISTANCE] << endl;

  for (i = 0;i < nb_component;i++) {
    if (file_flag) {
      os << "# ";
    }
    for (j = 0;j < nb_component;j++) {
      if (j != i) {
        os << setw(width) << sup_norm_dist[i][j];
      }
      else {
        os << setw(width) << " ";
      }
    }
    os << endl;
  }

  for (i = 0;i < nb_component;i++) {
    delete [] sup_norm_dist[i];
  }
  delete [] sup_norm_dist;

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet DiscreteMixture.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& DiscreteMixture::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , mixture_data , exhaustive , false);
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet DiscreteMixture dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool DiscreteMixture::ascii_write(StatError &error , const char *path ,
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
    ascii_write(out_file , mixture_data , exhaustive , true);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un melange de lois discretes et de la structure de donnees associee
 *  dans un fichier au format tableur.
 *
 *  arguments : stream, pointeur sur un objet DiscreteMixtureData.
 *
 *--------------------------------------------------------------*/

ostream& DiscreteMixture::spreadsheet_write(ostream &os , const DiscreteMixtureData *mixt_histo) const

{
  register int i , j;
  int bnb_parameter;
  double **sup_norm_dist , scale[DISCRETE_MIXTURE_NB_COMPONENT];
  const Distribution *pcomponent[DISCRETE_MIXTURE_NB_COMPONENT];


  os << STAT_word[STATW_MIXTURE] << "\t" << nb_component << "\t" << STAT_word[STATW_DISTRIBUTIONS] << endl;
  spreadsheet_characteristic_print(os);

  if (mixt_histo) {
    double likelihood , information;
    Test test(CHI2);


    os << "\n" << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
    mixt_histo->spreadsheet_characteristic_print(os);

    likelihood = Distribution::likelihood_computation(*mixt_histo);
    information = mixt_histo->FrequencyDistribution::information_computation();

    os << "\n" << STAT_label[STATL_LIKELIHOOD] << "\t" << likelihood << "\t"
       << STAT_label[STATL_NORMALIZED] << "\t" << likelihood / mixt_histo->nb_element << endl;
    os << STAT_label[STATL_MAX_LIKELIHOOD] << "\t" << information << "\t"
       << STAT_label[STATL_INFORMATION] << "\t" << information / mixt_histo->nb_element << endl;
    os << STAT_label[STATL_DEVIANCE] << "\t" << 2 * (information - likelihood) << endl;

    bnb_parameter = nb_parameter_computation();

    os << "\n" << bnb_parameter << "\t" << STAT_label[STATL_FREE_PARAMETERS] << "\t"
       << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " (" << STAT_criterion_word[AIC] << ")\t"
       << 2 * (likelihood - bnb_parameter) << endl;

    if (0 < bnb_parameter < mixt_histo->nb_element - 1) {
      os << bnb_parameter << "\t" << STAT_label[STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " (" << STAT_criterion_word[AICc] << ")\t"
         << 2 * (likelihood - (double)(bnb_parameter * mixt_histo->nb_element) /
            (double)(mixt_histo->nb_element - bnb_parameter - 1)) << endl;
    }

    os << bnb_parameter << "\t" << STAT_label[STATL_FREE_PARAMETERS] << "\t"
       << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " (" << STAT_criterion_word[BIC] << ")\t"
       << 2 * likelihood - bnb_parameter * log((double)mixt_histo->nb_element) << endl;

    os << bnb_parameter << "\t" << STAT_label[STATL_FREE_PARAMETERS] << "\t"
       << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " (" << STAT_criterion_word[BICc] << ")\t"
       << 2 * likelihood - penalty_computation() << endl;

    likelihood = likelihood_computation(*mixt_histo);
    information = mixt_histo->information_computation();

    os << "\n" << STAT_label[STATL_CLASSIFICATION_LIKELIHOOD] << "\t" << likelihood << "\t"
       << STAT_label[STATL_NORMALIZED] << "\t" << likelihood / mixt_histo->nb_element << endl;
    os << STAT_label[STATL_MAX_CLASSIFICATION_LIKELIHOOD] << "\t" << information << "\t"
       << STAT_label[STATL_INFORMATION] << "\t" << information / mixt_histo->nb_element << endl;

    chi2_fit(*mixt_histo , test);
    os << "\n";
    test.spreadsheet_print(os);
  }

  else {
    for (i = 0;i < nb_component;i++) {
      os << "\n" << STAT_word[STATW_DISTRIBUTION] << "\t" << i + 1 << "\t"
         << STAT_word[STATW_WEIGHT] << "\t"  << weight->mass[i] << endl;
      component[i]->spreadsheet_print(os);
      component[i]->spreadsheet_parametric_characteristic_print(os , true);
    }
  }

  for (i = 0;i < nb_component;i++) {
    pcomponent[i] = component[i];

    if (mixt_histo) {
      scale[i] = weight->mass[i] * mixt_histo->nb_element;
    }
    else {
      scale[i] = weight->mass[i];
    }
  }

  os << "\n";
  if (mixt_histo) {
    os << "\t" << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
  }
  for (i = 0;i < nb_component;i++) {
    os << "\t" << STAT_label[STATL_DISTRIBUTION] << " " << i + 1;
  }
  os << "\t" << STAT_label[STATL_MIXTURE];
  if (mixt_histo) {
    os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
       << " " << STAT_label[STATL_FUNCTION];
  }
  os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_MIXTURE]
     << " " << STAT_label[STATL_FUNCTION] << endl;

  spreadsheet_print(os , nb_component , pcomponent , scale , true , mixt_histo);

  if (mixt_histo) {
    os << "\n" << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
    mixt_histo->weight->spreadsheet_characteristic_print(os);

    os << "\n\t" << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
       << "\t" << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_DISTRIBUTION] << endl;

    weight->Distribution::spreadsheet_print(os , false , false , false , mixt_histo->weight);

    for (i = 0;i < nb_component;i++) {
      os << "\n" << STAT_word[STATW_DISTRIBUTION] << "\t" << i + 1 << "\t"
         << STAT_word[STATW_WEIGHT] << "\t"  << weight->mass[i] << endl;
      component[i]->spreadsheet_print(os);
      component[i]->spreadsheet_parametric_characteristic_print(os , true);

      os << "\n" << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i + 1 << "\t";
      mixt_histo->component[i]->spreadsheet_characteristic_print(os , true);

      if (mixt_histo->component[i]->nb_element > 0) {
        os << "\n\t" << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i + 1
           << "\t" << STAT_label[STATL_DISTRIBUTION] << " " << i + 1
           << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
           << " " << i + 1 << " " << STAT_label[STATL_FUNCTION]
           << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION]
           << " " << i + 1 << " " << STAT_label[STATL_FUNCTION] << endl;

        component[i]->Distribution::spreadsheet_print(os , true , false , false ,
                                                      mixt_histo->component[i]);
      }
    }
  }

  sup_norm_dist = new double*[nb_component];
  for (i = 0;i < nb_component;i++) {
    sup_norm_dist[i] = new double[nb_component];
  }

  for (i = 0;i < nb_component;i++) {
    for (j = i + 1;j < nb_component;j++) {
      sup_norm_dist[i][j] = component[i]->sup_norm_distance_computation(*(component[j]));
      sup_norm_dist[j][i] = sup_norm_dist[i][j];
    }
  }

  os << "\n" << STAT_label[STATL_COMPONENT_DISTANCE] << endl;

  for (i = 0;i < nb_component;i++) {
    for (j = 0;j < nb_component;j++) {
      if (j != i) {
        os << sup_norm_dist[i][j];
      }
      os << "\t";
    }
    os << endl;
  }

  for (i = 0;i < nb_component;i++) {
    delete [] sup_norm_dist[i];
  }
  delete [] sup_norm_dist;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet DiscreteMixture dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool DiscreteMixture::spreadsheet_write(StatError &error , const char *path) const

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
    spreadsheet_write(out_file , mixture_data);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un melange de lois discretes et de la structure de donnees associee.
 *
 *  arguments : prefixe des fichiers, titre des figures,
 *              pointeur sur un objet DiscreteMixtureData.
 *
 *--------------------------------------------------------------*/

bool DiscreteMixture::plot_write(const char *prefix , const char *title ,
                                 const DiscreteMixtureData *mixt_histo) const

{
  bool status;
  register int i , j , k;
  int nb_histo = 0;
  double scale[DISCRETE_MIXTURE_NB_COMPONENT + 2];
  const Distribution *pdist[DISCRETE_MIXTURE_NB_COMPONENT + 2];
  const FrequencyDistribution *phisto[DISCRETE_MIXTURE_NB_COMPONENT + 2];
  ostringstream data_file_name[DISCRETE_MIXTURE_NB_COMPONENT + 1];


  // ecriture des fichiers de donnees

  data_file_name[0] << prefix << ".dat";

  if (mixt_histo) {
    pdist[0] = this;
    phisto[nb_histo++] = mixt_histo;
    scale[0] = mixt_histo->nb_element;

    pdist[1] = weight;
    phisto[nb_histo++] = mixt_histo->weight;
    scale[1] = mixt_histo->weight->nb_element;

    for (i = 0;i < nb_component;i++) {
      pdist[i + 2] = component[i];
      if (mixt_histo->component[i]->nb_element > 0) {
        phisto[nb_histo++] = mixt_histo->component[i];
        scale[i + 2] = mixt_histo->component[i]->nb_element;
      }
      else {
        scale[i + 2] = 1.;
      }
    }

    status = stat_tool::plot_print((data_file_name[0].str()).c_str() , nb_component + 2 , pdist ,
                                   scale , NULL , nb_histo , phisto);
  }

  else {
    status = plot_print((data_file_name[0].str()).c_str());
  }

  if (status) {
    for (i = 0;i < nb_component;i++) {
      data_file_name[i + 1] << prefix << i + 1 << ".dat";

      pdist[0] = component[i];

      if (mixt_histo) {
        scale[0] = weight->mass[i] * mixt_histo->nb_element;
      }
      else {
        scale[0] = weight->mass[i];
      }

      stat_tool::plot_print((data_file_name[i + 1].str()).c_str() , 1 , pdist ,
                            scale , NULL , 0 , NULL);
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

      if (mixt_histo) {
        out_file << "plot [0:" << nb_value - 1 << "] [0:"
                 << (int)(MAX(mixt_histo->max , max * mixt_histo->nb_element) * YSCALE) + 1
                 << "] \"" << label((data_file_name[0].str()).c_str()) << "\" using 1 title \""
                 << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses,\\" << endl;
      }

      else {
        out_file << "plot [0:" << nb_value - 1 << "] [0:"
                 << MIN(max * YSCALE , 1.) << "] ";
      }

      for (j = 0;j < nb_component;j++) {
        out_file << "\"" << label((data_file_name[j + 1].str()).c_str()) << "\" title \""
                 << STAT_label[STATL_DISTRIBUTION] << " " << j + 1;
        component[j]->plot_title_print(out_file);
        out_file << "\" with linespoints,\\" << endl;
      }

      out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\"";
      if (mixt_histo) {
        out_file << " using " << nb_histo + 1;
      }
      out_file << " title \"" << STAT_label[STATL_MIXTURE] << "\" with linespoints" << endl;

      if (nb_value - 1 < TIC_THRESHOLD) {
        out_file << "set xtics autofreq" << endl;
      }

      if (mixt_histo) {
        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;

        if (weight->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }

        out_file << "plot [0:" << weight->nb_value - 1 << "] [0:"
                 << (int)(MAX(mixt_histo->weight->max ,
                              weight->max * mixt_histo->weight->nb_element) * YSCALE) + 1
                 << "] \"" << label((data_file_name[0].str()).c_str()) << "\" using " << 2
                 << " title \"" << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                 << "\" with impulses,\\" << endl;
        out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << nb_histo + 2
                 << " title \"" << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_DISTRIBUTION]
                 << "\" with linespoints" << endl;

        if (weight->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }

        j = 3;
        for (k = 0;k < nb_component;k++) {
          if (mixt_histo->component[k]->nb_element > 0) {
            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;

            if (component[k]->nb_value - 1 < TIC_THRESHOLD) {
              out_file << "set xtics 0,1" << endl;
            }

            out_file << "plot [0:" << component[k]->nb_value - 1 << "] [0:"
                     << (int)(MAX(mixt_histo->component[k]->max ,
                                  component[k]->max * mixt_histo->component[k]->nb_element) * YSCALE) + 1
                     << "] \"" << label((data_file_name[0].str()).c_str()) << "\" using " << j++
                     << " title \"" << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << k + 1
                     << "\" with impulses,\\" << endl;
            out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << nb_histo + k + 3
                     << " title \"" << STAT_label[STATL_DISTRIBUTION] << " " << k + 1;
            component[k]->plot_title_print(out_file);
            out_file << "\" with linespoints" << endl;

            if (component[k]->nb_value - 1 < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }
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
 *  Sortie Gnuplot d'un objet DiscreteMixture.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool DiscreteMixture::plot_write(StatError &error , const char *prefix ,
                                 const char *title) const

{
  bool status = plot_write(prefix , title , mixture_data);

  error.init();

  if (!status) {
    error.update(STAT_error[STATR_FILE_PREFIX]);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un melange de lois discretes et de la structure de donnees associee.
 *
 *  argument : pointeur sur un objet DiscreteMixtureData.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* DiscreteMixture::get_plotable(const DiscreteMixtureData *mixt_histo) const

{
  register int i , j;
  int nb_plot_set , xmax;
  double scale;
  ostringstream title , legend;
  MultiPlotSet *plot_set;


  // calcul du nombre de graphiques

  nb_plot_set = 1;
  if (mixt_histo) {
    nb_plot_set += 2;
    for (i = 0;i < nb_component;i++) {
      if (mixt_histo->component[i]->nb_element > 0) {
        nb_plot_set++;
      }
    }
  }

  plot_set = new MultiPlotSet(nb_plot_set);
  MultiPlotSet &plot = *plot_set;

  title.str("");
  title << STAT_label[STATL_MIXTURE];
  if (mixt_histo) {
    title << " " << STAT_label[STATL_FIT];
  }
  plot.title = title.str();

  plot.border = "15 lw 0";

  // 1ere vue : melange ajuste

  xmax = nb_value - 1;
  if ((cumul[xmax] > 1. - DOUBLE_ERROR) &&
      (mass[xmax] > PLOT_MASS_THRESHOLD)) {
    xmax++;
  }
  plot[0].xrange = Range(0 , xmax);

  if (nb_value - 1 < TIC_THRESHOLD) {
    plot[0].xtics = 1;
  }

  if (mixt_histo) {
    plot[0].yrange = Range(0. , ceil(MAX(mixt_histo->max , max * mixt_histo->nb_element)
                                     * YSCALE));

    plot[0].resize(nb_component + 2);

    plot[0][0].legend = STAT_label[STATL_FREQUENCY_DISTRIBUTION];

    plot[0][0].style = "impulses";

    mixt_histo->plotable_frequency_write(plot[0][0]);

    scale = mixt_histo->nb_element;
    i = 1;
  }

  else {
    plot[0].yrange = Range(0. , MIN(max * YSCALE , 1.));

    plot[0].resize(nb_component + 1);

    scale = 1.;
    i = 0;
  }

  for (j = 0;j < nb_component;j++) {
    legend.str("");
    legend << STAT_label[STATL_DISTRIBUTION] << " " << j + 1;
    component[j]->plot_title_print(legend);
    plot[0][i + j].legend = legend.str();

    plot[0][i + j].style = "linespoints";

    if (mixt_histo) {
      component[j]->plotable_mass_write(plot[0][i + j] , weight->mass[j] * mixt_histo->nb_element);
    }
    else {
      component[j]->plotable_mass_write(plot[0][i + j] , weight->mass[j]);
    }
  }

  plot[0][i + nb_component].legend = STAT_label[STATL_MIXTURE];

  plot[0][i + nb_component].style = "linespoints";

  plotable_mass_write(plot[0][i + nb_component] , scale);

  if (mixt_histo) {

    // 2eme vue : fonctions de repartition

    plot[1].xrange = Range(0 , xmax);
    plot[1].yrange = Range(0. , 1.);

    if (nb_value - 1 < TIC_THRESHOLD) {
      plot[1].xtics = 1;
    }

    plot[1].resize(2);

    legend.str("");
    legend << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
           << " " << STAT_label[STATL_FUNCTION];
    plot[1][0].legend = legend.str();

    plot[1][0].style = "linespoints";

    mixt_histo->plotable_cumul_write(plot[1][0]);

    legend.str("");
    legend << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_MIXTURE]
           << " " << STAT_label[STATL_FUNCTION];
    plot[1][1].legend = legend.str();

    plot[1][1].style = "linespoints";

    plotable_cumul_write(plot[1][1]);

    // 3eme vue : poids

    plot[2].xrange = Range(0 , weight->nb_value - 1);
    plot[2].yrange = Range(0 , ceil(MAX(mixt_histo->weight->max ,
                                    weight->max * mixt_histo->weight->nb_element)
                                    * YSCALE));

    if (weight->nb_value - 1 < TIC_THRESHOLD) {
      plot[2].xtics = 1;
    }

    plot[2].resize(2);

    legend.str("");
    legend << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    plot[2][0].legend = legend.str();

    plot[2][0].style = "impulses";

    mixt_histo->weight->plotable_frequency_write(plot[2][0]);

    legend.str("");
    legend << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_DISTRIBUTION];
    plot[2][1].legend = legend.str();

    plot[2][1].style = "linespoints";

    weight->plotable_mass_write(plot[2][1] , mixt_histo->weight->nb_element);

    i = 3;
    for (j = 0;j < nb_component;j++) {
      if (mixt_histo->component[j]->nb_element > 0) {

        // vues suivantes : composantes ajustees

        xmax = component[j]->nb_value - 1;
        if ((component[j]->cumul[xmax] > 1. - DOUBLE_ERROR) &&
            (component[j]->mass[xmax] > PLOT_MASS_THRESHOLD)) {
          xmax++;
        }
        plot[i].xrange = Range(0 , xmax);

        plot[i].yrange = Range(0 , ceil(MAX(mixt_histo->component[j]->max ,
                                        component[j]->max * mixt_histo->component[j]->nb_element)
                                        * YSCALE));

        if (component[j]->nb_value - 1 < TIC_THRESHOLD) {
          plot[i].xtics = 1;
        }

        plot[i].resize(2);

        legend.str("");
        legend << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << j + 1;
        plot[i][0].legend = legend.str();

        plot[i][0].style = "impulses";

        mixt_histo->component[j]->plotable_frequency_write(plot[i][0]);

        legend.str("");
        legend << STAT_label[STATL_DISTRIBUTION] << " " << j + 1;
        plot[i][1].legend = legend.str();

        plot[i][1].style = "linespoints";

        component[j]->plotable_mass_write(plot[i][1] , mixt_histo->component[j]->nb_element);
  
        i++;
      }
    }
  }

  return plot_set;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet DiscreteMixture.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* DiscreteMixture::get_plotable() const

{
  return get_plotable(mixture_data);
}


/*--------------------------------------------------------------*
 *
 *  Calcul du nombre de parametres independants.
 *
 *--------------------------------------------------------------*/

int DiscreteMixture::nb_parameter_computation() const

{
  register int i;
  int bnb_parameter = nb_component - 1;


  for (i = 0;i < nb_component;i++) {
    bnb_parameter += component[i]->nb_parameter_computation();
  }

  return bnb_parameter;
}


/*--------------------------------------------------------------*
 *
 *  Calcul d'une penalite adaptative.
 *
 *--------------------------------------------------------------*/

double DiscreteMixture::penalty_computation() const

{
  register int i;
  double penalty = 0.;


  if (mixture_data) {
    penalty += (nb_component - 1) * log((double)mixture_data->nb_element);
    for (i = 0;i < nb_component;i++) {
      penalty += component[i]->nb_parameter_computation() *
                 log(weight->mass[i] * mixture_data->nb_element);
    }
  }

  return penalty;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe DiscreteMixtureData.
 *
 *--------------------------------------------------------------*/

DiscreteMixtureData::DiscreteMixtureData()

{
  mixture = NULL;
  nb_component = 0;
  weight = NULL;
  component = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe DiscreteMixtureData.
 *
 *  arguments : reference sur un objet FrequencyDistribution,
 *              nombre de composantes.
 *
 *--------------------------------------------------------------*/

DiscreteMixtureData::DiscreteMixtureData(const FrequencyDistribution &histo , int inb_component)
:FrequencyDistribution(histo)

{
  register int i;


  mixture = NULL;
  nb_component = inb_component;

  weight = new FrequencyDistribution(nb_component);

  component = new FrequencyDistribution*[nb_component];
  for (i = 0;i < nb_component;i++) {
    component[i] = new FrequencyDistribution(nb_value);
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe DiscreteMixtureData.
 *
 *  arguments : reference sur un objet FrequencyDistribution,
 *              pointeur sur un objet DiscreteMixture.
 *
 *--------------------------------------------------------------*/

DiscreteMixtureData::DiscreteMixtureData(const FrequencyDistribution &histo , const DiscreteMixture *pmixture)
:FrequencyDistribution(histo)

{
  register int i;


  mixture = new DiscreteMixture(*pmixture , false); 
  nb_component = pmixture->nb_component;

  weight = new FrequencyDistribution(nb_component);

  component = new FrequencyDistribution*[nb_component];
  for (i = 0;i < nb_component;i++) {
    component[i] = new FrequencyDistribution(nb_value);
  }

//  mixture->mixture_data = new DiscreteMixtureData(*this, false);   introduit par J.-B. mais pas coherent avec le design
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe DiscreteMixtureData.
 *
 *  argument : reference sur un objet DiscreteMixture.
 *
 *--------------------------------------------------------------*/

DiscreteMixtureData::DiscreteMixtureData(const DiscreteMixture &mixt)
:FrequencyDistribution(mixt)

{
  register int i;


  mixture = NULL;
  nb_component = mixt.nb_component;

  weight = new FrequencyDistribution(*(mixt.weight));

  component = new FrequencyDistribution*[nb_component];
  for (i = 0;i < mixt.nb_component;i++) {
    component[i] = new FrequencyDistribution(*(mixt.component[i]));
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet DiscreteMixtureData.
 *
 *  arguments : reference sur un objet DiscreteMixtureData,
 *              flag copie de l'objet DiscreteMixture.
 *
 *--------------------------------------------------------------*/

void DiscreteMixtureData::copy(const DiscreteMixtureData &mixt_histo , bool model_flag)

{
  register int i;


  if ((model_flag) && (mixt_histo.mixture)) {
    mixture = new DiscreteMixture(*(mixt_histo.mixture) , false);
  }
  else {
    mixture = NULL;
  }

  nb_component = mixt_histo.nb_component;

  weight = new FrequencyDistribution(*(mixt_histo.weight));

  component = new FrequencyDistribution*[nb_component];
  for (i = 0;i < nb_component;i++) {
    component[i] = new FrequencyDistribution(*(mixt_histo.component[i]));
  }
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet DiscreteMixtureData.
 *
 *--------------------------------------------------------------*/

void DiscreteMixtureData::remove()

{
  delete mixture;

  delete weight;

  if (component) {
    register int i;

    for (i = 0;i < nb_component;i++) {
      delete component[i];
    }
    delete [] component;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe DiscreteMixtureData.
 *
 *--------------------------------------------------------------*/

DiscreteMixtureData::~DiscreteMixtureData()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe DiscreteMixtureData.
 *
 *  argument : reference sur un objet DiscreteMixtureData.
 *
 *--------------------------------------------------------------*/

DiscreteMixtureData& DiscreteMixtureData::operator=(const DiscreteMixtureData &mixt_histo)

{
  if (&mixt_histo != this) {
    remove();
    delete [] frequency;

    FrequencyDistribution::copy(mixt_histo);
    copy(mixt_histo);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Extraction d'une composante empirique.
 *
 *  arguments : reference sur un objet StatError, indice de la composante.
 *
 *--------------------------------------------------------------*/

DiscreteDistributionData* DiscreteMixtureData::extract(StatError &error , int index) const

{
  bool status = true;
  DiscreteDistributionData *pcomponent;


  pcomponent = NULL;
  error.init();

  if ((index < 1) || (index > nb_component)) {
    status = false;
    error.update(STAT_error[STATR_FREQUENCY_DISTRIBUTION_INDEX]);
  }
  else {
    index--;
    if (component[index]->nb_element == 0) {
      status = false;
      error.update(STAT_error[STATR_EMPTY_SAMPLE]);
    }
  }

  if (status) {
    pcomponent = new DiscreteDistributionData(*component[index] ,
                                              (mixture ? mixture->component[index] : NULL));
  }

  return pcomponent;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet DiscreteMixtureData.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& DiscreteMixtureData::line_write(ostream &os) const

{
  os << STAT_label[STATL_SAMPLE_SIZE] << ": " << nb_element << "   "
     << STAT_label[STATL_MEAN] << ": " << mean << "   "
     << STAT_label[STATL_VARIANCE] << ": " << variance;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet DiscreteMixtureData.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& DiscreteMixtureData::ascii_write(ostream &os , bool exhaustive) const

{
  if (mixture) {
    mixture->ascii_write(os , this , exhaustive , false);
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet DiscreteMixtureData dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool DiscreteMixtureData::ascii_write(StatError &error , const char *path ,
                                      bool exhaustive) const

{
  bool status = false;


  if (mixture) {
    ofstream out_file(path);

    error.init();

    if (!out_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_NAME]);
    }

    else {
      status = true;
      mixture->ascii_write(out_file , this , exhaustive , true);
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet DiscreteMixtureData dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool DiscreteMixtureData::spreadsheet_write(StatError &error , const char *path) const

{
  bool status = false;


  if (mixture) {
    ofstream out_file(path);

    error.init();

    if (!out_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_NAME]);
    }

    else {
      status = true;
      mixture->spreadsheet_write(out_file , this);
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet DiscreteMixtureData.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool DiscreteMixtureData::plot_write(StatError &error , const char *prefix ,
                                     const char *title) const

{
  bool status = false;


  if (mixture) {
    status = mixture->plot_write(prefix , title , this);

    error.init();

    if (!status) {
      error.update(STAT_error[STATR_FILE_PREFIX]);
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet DiscreteMixtureData.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* DiscreteMixtureData::get_plotable() const

{
  MultiPlotSet *plot_set;


  if (mixture) {
    plot_set = mixture->get_plotable(this);
  }
  else {
    plot_set = NULL;
  }

  return plot_set;
}


};  // namespace stat_tool
