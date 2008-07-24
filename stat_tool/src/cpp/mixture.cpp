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
#include "tool/rw_locale.h"
// #include <rw/vstream.h>
// #include <rw/rwfile.h>
#include "stat_tools.h"
#include "distribution.h"
#include "mixture.h"
#include "stat_label.h"

using namespace std;


extern char* label(const char *file_name);



/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe Mixture.
 *
 *--------------------------------------------------------------*/

Mixture::Mixture()

{
  mixture_data = 0;
  nb_component = 0;
  weight = 0;
  component = 0;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Mixture.
 *
 *  arguments : nombre de composantes, poids, pointeurs sur les composantes.
 *
 *--------------------------------------------------------------*/

Mixture::Mixture(int inb_component , double *pweight , const Parametric **pcomponent)

{
  register int i;


  mixture_data = 0;
  nb_component = inb_component;

  weight = new Parametric(nb_component);
  for (i = 0;i < nb_component;i++) {
    weight->mass[i] = *pweight++;
  }
  weight->cumul_computation();
  weight->max_computation();

  nb_value = 0;
  component = new Parametric*[nb_component];

  for (i = 0;i < nb_component;i++) {
    component[i] = new Parametric(*pcomponent[i] , 'n');
    if (component[i]->nb_value > nb_value) {
      nb_value = component[i]->nb_value;
    }
  }

  Distribution::init(nb_value);

  computation(1 , CUMUL_THRESHOLD , false);
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Mixture.
 *
 *  arguments : reference sur un objet Mixture, flag sur les composantes a copier,
 *              nombre de valeurs.
 *
 *--------------------------------------------------------------*/

Mixture::Mixture(const Mixture &mixt , bool *component_flag , int inb_value)

{
  register int i;


  mixture_data = 0;
  nb_component = mixt.nb_component;

  weight = new Parametric(nb_component);

  nb_value = inb_value;
  component = new Parametric*[nb_component];

  for (i = 0;i < nb_component;i++) {
    if (component_flag[i]) {
      component[i] = new Parametric(inb_value , mixt.component[i]->ident);
    }
    else {
      component[i] = new Parametric(*(mixt.component[i]));
      if (component[i]->nb_value > nb_value) {
        nb_value = component[i]->nb_value;
      }
    }
  }

  Distribution::init(nb_value);
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Mixture.
 *
 *  arguments : nombre de composantes, pointeurs sur les composantes.
 *
 *--------------------------------------------------------------*/

Mixture::Mixture(int inb_component , const Parametric **pcomponent)

{
  register int i;


  mixture_data = 0;
  nb_component = inb_component;

  weight = 0;

  component = new Parametric*[nb_component];
  for (i = 0;i < nb_component;i++) {
    component[i] = new Parametric(*pcomponent[i] , 'n');
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Mixture.
 *
 *  arguments : reference sur un objet Mixture,
 *              flag copie de l'objet Mixture_data.
 *
 *--------------------------------------------------------------*/

void Mixture::copy(const Mixture &mixt , bool data_flag)

{
  register int i;


  if ((data_flag) && (mixt.mixture_data)) {
    mixture_data = new Mixture_data(*(mixt.mixture_data) , false);
  }
  else {
    mixture_data = 0;
  }

  nb_component = mixt.nb_component;

  weight = new Parametric(*(mixt.weight));

  component = new Parametric*[nb_component];
  for (i = 0;i < nb_component;i++) {
    component[i] = new Parametric(*(mixt.component[i]));
  }
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet Mixture.
 *
 *--------------------------------------------------------------*/

void Mixture::remove()

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
 *  Destructeur de la classe Mixture.
 *
 *--------------------------------------------------------------*/

Mixture::~Mixture()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Mixture.
 *
 *  argument : reference sur un objet Mixture.
 *
 *--------------------------------------------------------------*/

Mixture& Mixture::operator=(const Mixture &mixt)

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
 *  arguments : reference sur un objet Format_error, indice de la composante.
 *
 *--------------------------------------------------------------*/

Parametric_model* Mixture::extract(Format_error &error , int index) const

{
  Parametric_model *pcomponent;


  if ((index < 1) || (index > nb_component)) {
    pcomponent = 0;
    error.update(STAT_error[STATR_DISTRIBUTION_INDEX]);
  }

  else {
    index--;
    pcomponent = new Parametric_model(*component[index] ,
                                      (mixture_data ? mixture_data->component[index] : 0));
  }

  return pcomponent;
}


/*--------------------------------------------------------------*
 *
 *  Extraction de la partie "donnees" d'un objet Mixture.
 *
 *  argument : reference sur un objet Format_error.
 *
 *--------------------------------------------------------------*/

Mixture_data* Mixture::extract_data(Format_error &error) const

{
  Mixture_data *mixt_histo;


  error.init();

  if (!mixture_data) {
    mixt_histo = 0;
    error.update(STAT_error[STATR_NO_DATA]);
  }

  else {
    mixt_histo = new Mixture_data(*mixture_data);
    mixt_histo->mixture = new Mixture(*this , false);
  }

  return mixt_histo;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Mixture a partir de poids et de composantes.
 *
 *  arguments : reference sur un objet Format_error, nombre de composantes,
 *              poids, pointeurs sur les composantes.
 *
 *--------------------------------------------------------------*/

Mixture* mixture_building(Format_error &error , int nb_component , double *weight ,
                          const Parametric **component)

{
  bool status;
  register int i;
  double cumul;
  Mixture *mixt;


  mixt = 0;
  error.init();

  if ((nb_component < 2) || (nb_component > MIXTURE_NB_COMPONENT)) {
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
      mixt = new Mixture(nb_component , weight , component);
    }
  }

  return mixt;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Mixture a partir d'un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

Mixture* mixture_ascii_read(Format_error &error , const char *path ,
                            double cumul_threshold)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  bool status , lstatus;
  register int i , j;
  int line;
  long index , nb_component;
  double cumul , weight[MIXTURE_NB_COMPONENT];
  const Parametric **component;
  Mixture *mixt;
  ifstream in_file(path);


  mixt = 0;
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
          if ((lstatus) && ((nb_component < 2) || (nb_component > MIXTURE_NB_COMPONENT))) {
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
      component = new const Parametric*[nb_component];
      for (i = 0;i < nb_component;i++) {
        component[i] = 0;
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

            component[i] = parametric_parsing(error , in_file , line ,
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
        mixt = new Mixture(nb_component , weight , component);
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
 *  Ecriture sur une ligne d'un objet Mixture.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& Mixture::line_write(ostream &os) const

{
  os << nb_component << " " << STAT_word[STATW_DISTRIBUTIONS] << "   "
     << STAT_label[STATL_MEAN] << ": " << mean << "   "
     << STAT_label[STATL_VARIANCE] << ": " << variance;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un melange et de la structure de donnees associee.
 *
 *  arguments : stream, pointeur sur un objet Mixture_data,
 *              flag niveau de detail, flag fichier.
 *
 *--------------------------------------------------------------*/

ostream& Mixture::ascii_write(ostream &os , const Mixture_data *mixt_histo ,
                              bool exhaustive , bool file_flag) const

{
  register int i;
  int bnb_parameter;
  double scale[MIXTURE_NB_COMPONENT];
  const Distribution *pcomponent[MIXTURE_NB_COMPONENT];


  os << STAT_word[STATW_MIXTURE] << " " << nb_component << " " << STAT_word[STATW_DISTRIBUTIONS] << endl;
  ascii_characteristic_print(os , false , file_flag);

  if (mixt_histo) {
    double likelihood , information;
    Test test(CHI2);


    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_HISTOGRAM] << " - ";
    mixt_histo->ascii_characteristic_print(os , false , file_flag);

    likelihood = Distribution::likelihood_computation(*mixt_histo);
    information = mixt_histo->Histogram::information_computation();

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
      component[i]->ascii_characteristic_print(os , exhaustive , file_flag);
    }
  }

  if (exhaustive) {
    for (i = 0;i < nb_component;i++) {
      pcomponent[i] = component[i];

      if (mixt_histo) {
        scale[i] = mixt_histo->nb_element * weight->mass[i];
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
      os << " | " << STAT_label[STATL_HISTOGRAM];
    }
    os << " | " << STAT_label[STATL_MIXTURE];
    for (i = 0;i < nb_component;i++) {
      os << " | " << STAT_label[STATL_DISTRIBUTION] << " " << i + 1;
    }
    if (mixt_histo) {
      os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM] << " "
         << STAT_label[STATL_FUNCTION];
    }
    os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_MIXTURE] << " "
       << STAT_label[STATL_FUNCTION] << endl;

    ascii_print(os , nb_component , pcomponent , scale , file_flag , true , mixt_histo);
  }

  if (mixt_histo) {
    if (exhaustive) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
      mixt_histo->weight->ascii_characteristic_print(os , false , file_flag);

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << "   | " << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_HISTOGRAM] << " "
         << " | " << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_DISTRIBUTION] << endl;

      weight->Distribution::ascii_print(os , file_flag , false , false , mixt_histo->weight);
    }

    for (i = 0;i < nb_component;i++) {
      os << "\n" << STAT_word[STATW_DISTRIBUTION] << " " << i + 1 << "  "
         << STAT_word[STATW_WEIGHT] << " : "  << weight->mass[i] << endl;
      component[i]->ascii_print(os);
      component[i]->ascii_characteristic_print(os , exhaustive , file_flag);

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_HISTOGRAM] << " " << i + 1 << " - ";
      mixt_histo->component[i]->ascii_characteristic_print(os , exhaustive , file_flag);

      if ((exhaustive) && (mixt_histo->component[i]->nb_element > 0)) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << "   | " << STAT_label[STATL_HISTOGRAM] << " " << i + 1
           << " | " << STAT_label[STATL_DISTRIBUTION] << " " << i + 1
           << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM]
           << " " << i + 1 << " " << STAT_label[STATL_FUNCTION]
           << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION]
           << " " << i + 1 << " " << STAT_label[STATL_FUNCTION] << endl;

        component[i]->Distribution::ascii_print(os , file_flag , true , false , mixt_histo->component[i]);
      }
    }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Mixture.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Mixture::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , mixture_data , exhaustive , false);
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Mixture dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Mixture::ascii_write(Format_error &error , const char *path ,
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
 *  Ecriture d'un melange et de la structure de donnees associee
 *  dans un fichier au format tableur.
 *
 *  arguments : stream, pointeur sur un objet Mixture_data.
 *
 *--------------------------------------------------------------*/

ostream& Mixture::spreadsheet_write(ostream &os , const Mixture_data *mixt_histo) const

{
  register int i;
  int bnb_parameter;
  double scale[MIXTURE_NB_COMPONENT];
  const Distribution *pcomponent[MIXTURE_NB_COMPONENT];


  os << STAT_word[STATW_MIXTURE] << "\t" << nb_component << "\t" << STAT_word[STATW_DISTRIBUTIONS] << endl;
  spreadsheet_characteristic_print(os);

  if (mixt_histo) {
    double likelihood , information;
    Test test(CHI2);


    os << "\n" << STAT_label[STATL_HISTOGRAM] << "\t";
    mixt_histo->spreadsheet_characteristic_print(os);

    likelihood = Distribution::likelihood_computation(*mixt_histo);
    information = mixt_histo->Histogram::information_computation();

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
      component[i]->spreadsheet_characteristic_print(os , true);
    }
  }

  for (i = 0;i < nb_component;i++) {
    pcomponent[i] = component[i];

    if (mixt_histo) {
      scale[i] = mixt_histo->nb_element * weight->mass[i];
    }
    else {
      scale[i] = weight->mass[i];
    }
  }

  os << "\n";
  if (mixt_histo) {
    os << "\t" << STAT_label[STATL_HISTOGRAM];
  }
  os << "\t" << STAT_label[STATL_MIXTURE];
  for (i = 0;i < nb_component;i++) {
    os << "\t" << STAT_label[STATL_DISTRIBUTION] << " " << i + 1;
  }
  if (mixt_histo) {
    os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM] << " "
       << STAT_label[STATL_FUNCTION];
  }
  os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_MIXTURE] << " "
     << STAT_label[STATL_FUNCTION] << endl;

  spreadsheet_print(os , nb_component , pcomponent , scale , true , mixt_histo);

  if (mixt_histo) {
    os << "\n" << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_HISTOGRAM] << "\t";
    mixt_histo->weight->spreadsheet_characteristic_print(os);

    os << "\n\t" << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_HISTOGRAM] << " "
       << "\t" << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_DISTRIBUTION] << endl;

    weight->Distribution::spreadsheet_print(os , false , false , false , mixt_histo->weight);

    for (i = 0;i < nb_component;i++) {
      os << "\n" << STAT_word[STATW_DISTRIBUTION] << "\t" << i + 1 << "\t"
         << STAT_word[STATW_WEIGHT] << "\t"  << weight->mass[i] << endl;
      component[i]->spreadsheet_print(os);
      component[i]->spreadsheet_characteristic_print(os , true);

      os << "\n" << STAT_label[STATL_HISTOGRAM] << " " << i + 1 << "\t";
      mixt_histo->component[i]->spreadsheet_characteristic_print(os , true);

      if (mixt_histo->component[i]->nb_element > 0) {
        os << "\n\t" << STAT_label[STATL_HISTOGRAM] << " " << i + 1
           << "\t" << STAT_label[STATL_DISTRIBUTION] << " " << i + 1
           << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM]
           << " " << i + 1 << " " << STAT_label[STATL_FUNCTION]
           << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION]
           << " " << i + 1 << " " << STAT_label[STATL_FUNCTION] << endl;

        component[i]->Distribution::spreadsheet_print(os , true , false , false , mixt_histo->component[i]);
      }
    }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Mixture dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

bool Mixture::spreadsheet_write(Format_error &error , const char *path) const

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
 *  Sortie Gnuplot d'un melange et de la structure de donnees associee.
 *
 *  arguments : prefixe des fichiers, titre des figures,
 *              pointeur sur un objet Mixture_data.
 *
 *--------------------------------------------------------------*/

bool Mixture::plot_write(const char *prefix , const char *title ,
                         const Mixture_data *mixt_histo) const

{
  bool status;
  register int i , j , k;
  int nb_histo = 0 , index_dist[MIXTURE_NB_COMPONENT + 2];
  double scale[MIXTURE_NB_COMPONENT + 2];
  const Distribution *pdist[MIXTURE_NB_COMPONENT + 2];
  const Histogram *phisto[MIXTURE_NB_COMPONENT + 2];
  ostringstream data_file_name[MIXTURE_NB_COMPONENT + 1];


  // ecriture des fichiers de donnees

  data_file_name[0] << prefix << ".dat";

  if (mixt_histo) {
    pdist[0] = this;
    phisto[nb_histo] = mixt_histo;
    index_dist[nb_histo++] = 0;
    scale[0] = mixt_histo->nb_element;

    pdist[1] = weight;
    phisto[nb_histo] = mixt_histo->weight;
    index_dist[nb_histo++] = 1;
    scale[1] = mixt_histo->weight->nb_element;

    for (i = 0;i < nb_component;i++) {
      pdist[i + 2] = component[i];
      if (mixt_histo->component[i]->nb_element > 0) {
        phisto[nb_histo] = mixt_histo->component[i];
        index_dist[nb_histo++] = i + 2;
        scale[i + 2] = mixt_histo->component[i]->nb_element;
      }
      else {
        scale[i + 2] = 1.;
      }
    }

    status = ::plot_print((data_file_name[0].str()).c_str() , nb_component + 2 , pdist ,
                          scale , 0 , nb_histo , phisto , index_dist);
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

      ::plot_print((data_file_name[i + 1].str()).c_str() , 1 , pdist , scale , 0 , 0 , 0 , 0);
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
                 << STAT_label[STATL_HISTOGRAM] << "\" with impulses,\\" << endl;
        out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << nb_histo + 1
                 << " title \"" << STAT_label[STATL_MIXTURE] << "\" with linespoints";
      }

      else {
        out_file << "plot [0:" << nb_value - 1 << "] [0:"
                 << MIN(max * YSCALE , 1.) << "] \""
                 << label((data_file_name[0].str()).c_str()) << "\" title \""
                 << STAT_label[STATL_MIXTURE] << "\" with linespoints";
      }

      for (j = 0;j < nb_component;j++) {
        out_file << ",\\" << endl;
        out_file << "\"" << label((data_file_name[j + 1].str()).c_str()) << "\" title \""
                 << STAT_label[STATL_DISTRIBUTION] << " " << j + 1;
        component[j]->plot_title_print(out_file);
        out_file << "\" with linespoints";
      }
      out_file << endl;

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
                 << " title \"" << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_HISTOGRAM]
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
                     << " title \"" << STAT_label[STATL_HISTOGRAM] << " " << k + 1
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
 *  Sortie Gnuplot d'un objet Mixture.
 *
 *  arguments : reference sur un objet Format_error, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Mixture::plot_write(Format_error &error , const char *prefix ,
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
 *  Sortie graphique d'un melange et de la structure de donnees associee.
 *
 *  argument : pointeur sur un objet Mixture_data.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* Mixture::get_plotable(const Mixture_data *mixt_histo) const

{
  register int i , j;
  std::stringstream ss;


  // nombre de fenetres: nb_component + 2 si ajustement

  MultiPlotSet *plotset = new MultiPlotSet(mixt_histo ? nb_component + 2 : 1);
  MultiPlotSet &set = *plotset;

  set.title = "Mixture fit";
  set.border = "15 lw 0";

  // 1ere vue : melange ajuste

  if (nb_value - 1 < TIC_THRESHOLD) {
    set[0].xtics = 1;
  }

  set[0].xrange = Range(0, nb_value - 1);

  if (mixt_histo) {
    set[0].yrange = Range(0, (MAX(mixt_histo->max , max * mixt_histo->nb_element)
                              * YSCALE) + 1);
  }
  else {
    set[0].yrange = Range(0, MIN(max * YSCALE , 1.));
  }

  // definition du nombre de SinglePlot 

  i = 0;

  if (mixt_histo) {
    set[0].resize(nb_component + 2);
    set[0][i].legend = STAT_label[STATL_HISTOGRAM];
    set[0][i].style = "impulses";

    mixt_histo->plotable_frequency_write(set[0][i++]);
  }

  else {
   set[0].resize(nb_component + 1);
  }

  set[0][i].legend = STAT_label[STATL_MIXTURE];
  set[0][i].style = "linespoints";

  plotable_mass_write(set[0][i++] , (mixt_histo ? mixt_histo->nb_element : 1));
  
  for (j = 0;j < nb_component;j++) {
    ss.str("");
    ss << STAT_label[STATL_DISTRIBUTION] << " " << j + 1;
    set[0][i + j].legend = ss.str();

    set[0][i + j].style = "linespoints";

    if (mixt_histo) {
      component[j]->plotable_mass_write(set[0][i + j] , weight->mass[j] * mixt_histo->nb_element);
    }
    else {
      component[j]->plotable_mass_write(set[0][i + j] , weight->mass[j]);
    }
  }

  if (mixt_histo) {

    // 2eme vue : poids

    if (weight->nb_value - 1 < TIC_THRESHOLD) {
      set[1].xtics = 1;
    }

    set[1].xrange = Range(0, weight->nb_value - 1);
    set[1].yrange = Range(0, (MAX(mixt_histo->weight->max ,
                                  weight->max * mixt_histo->weight->nb_element) 
                              * YSCALE) + 1);

    set[1].resize(2);

    ss.str("");
    ss << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_HISTOGRAM];
    set[1][0].legend = ss.str();
    set[1][0].style = "impulses";

    ss.str("");
    ss << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_DISTRIBUTION];
    set[1][1].legend = ss.str();
    set[1][1].style = "linespoints";

    mixt_histo->weight->plotable_frequency_write(set[1][0]);
    weight->plotable_mass_write(set[1][1] , mixt_histo->weight->nb_element);

    i = 2;
    for (j = 0;j < nb_component;j++) {
      if (mixt_histo->component[j]->nb_element > 0) {

        // vues suivantes : composantes ajustees

        if (component[j]->nb_value - 1 < TIC_THRESHOLD) {
          set[i].xtics = 1;
        }

        set[i].xrange = Range(0 , component[j]->nb_value - 1);
        set[i].yrange = Range(0 , (MAX(mixt_histo->component[j]->max ,
                                       component[j]->max * mixt_histo->component[j]->nb_element)
                                   * YSCALE) + 1);

        set[i].resize(2);

        ss.str("");
        ss << STAT_label[STATL_HISTOGRAM] << " " << j + 1;
        set[i][0].legend = ss.str();
        set[i][0].style = "impulses";

        ss.str("");
        ss << STAT_label[STATL_DISTRIBUTION] << " " << j + 1;
        set[i][1].legend = ss.str();
        set[i][1].style = "linespoints";

        mixt_histo->component[j]->plotable_frequency_write(set[i][0]);
        component[j]->plotable_mass_write(set[i][1] , mixt_histo->component[j]->nb_element);
  
        i++;
      }
    }
  }

  return plotset;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet Mixture.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* Mixture::get_plotable() const

{
  MultiPlotSet *set;


  set = get_plotable(mixture_data);

  return set;
}


/*--------------------------------------------------------------*
 *
 *  Calcul du nombre de parametres independants.
 *
 *--------------------------------------------------------------*/

int Mixture::nb_parameter_computation() const

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

double Mixture::penalty_computation() const

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
 *  Fonctions pour la persistance.
 *
 *--------------------------------------------------------------*/

/* RWDEFINE_COLLECTABLE(Mixture , STATI_MIXTURE);


RWspace Mixture::binaryStoreSize() const

{
  register int i;
  RWspace size;


  size = Distribution::binaryStoreSize() + sizeof(nb_component) + weight->binaryStoreSize();
  for (i = 0;i < nb_component;i++) {
    size += component[i]->binaryStoreSize();
  }

  if (mixture_data) {
    size += mixture_data->recursiveStoreSize();
  }

  return size;
}


void Mixture::restoreGuts(RWvistream &is)

{
  register int i;


  remove();

  Distribution::restoreGuts(is);

  is >> nb_component;

  weight = new Parametric();
  weight->restoreGuts(is);

  component = new Parametric*[nb_component];
  for (i = 0;i < nb_component;i++) {
    component[i] = new Parametric();
    component[i]->restoreGuts(is);
  }

  is >> mixture_data;
  if (mixture_data == RWnilCollectable) {
    mixture_data = 0;
  }
}


void Mixture::restoreGuts(RWFile &file)

{
  register int i;


  remove();

  Distribution::restoreGuts(file);

  file.Read(nb_component);

  weight = new Parametric();
  weight->restoreGuts(file);

  component = new Parametric*[nb_component];
  for (i = 0;i < nb_component;i++) {
    component[i] = new Parametric();
    component[i]->restoreGuts(file);
  }

  file >> mixture_data;
  if (mixture_data == RWnilCollectable) {
    mixture_data = 0;
  }
}


void Mixture::saveGuts(RWvostream &os) const

{
  register int i;


  Distribution::saveGuts(os);

  os << nb_component;
  weight->saveGuts(os);

  for (i = 0;i < nb_component;i++) {
    component[i]->saveGuts(os);
  }

  if (mixture_data) {
    os << mixture_data;
  }
  else {
    os << RWnilCollectable;
  }
}


void Mixture::saveGuts(RWFile &file) const

{
  register int i;


  Distribution::saveGuts(file);

  file.Write(nb_component);
  weight->saveGuts(file);

  for (i = 0;i < nb_component;i++) {
    component[i]->saveGuts(file);
  }

  if (mixture_data) {
    file << mixture_data;
  }
  else {
    file << RWnilCollectable;
  }
} */




// ###################### Mixture_data #########################################



/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe Mixture_data.
 *
 *--------------------------------------------------------------*/

Mixture_data::Mixture_data()

{
  mixture = 0;
  nb_component = 0;
  weight = 0;
  component = 0;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Mixture_data.
 *
 *  arguments : reference sur un objet Histogram,
 *              nombre d'histogrammes.
 *
 *--------------------------------------------------------------*/

Mixture_data::Mixture_data(const Histogram &histo , int inb_component)
:Histogram(histo)

{
  register int i;


  mixture = 0;
  nb_component = inb_component;

  weight = new Histogram(nb_component);

  component = new Histogram*[nb_component];
  for (i = 0;i < nb_component;i++) {
    component[i] = new Histogram(nb_value);
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Mixture_data.
 *
 *  argument : reference sur un objet Mixture.
 *
 *--------------------------------------------------------------*/

Mixture_data::Mixture_data(const Mixture &mixt)
:Histogram(mixt)

{
  register int i;


  mixture = 0;
  nb_component = mixt.nb_component;

  weight = new Histogram(*(mixt.weight));

  component = new Histogram*[nb_component];
  for (i = 0;i < mixt.nb_component;i++) {
    component[i] = new Histogram(*(mixt.component[i]));
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Mixture_data.
 *
 *  arguments : reference sur un objet Mixture_data,
 *              flag copie de l'objet Mixture.
 *
 *--------------------------------------------------------------*/

void Mixture_data::copy(const Mixture_data &mixt_histo , bool model_flag)

{
  register int i;


  if ((model_flag) && (mixt_histo.mixture)) {
    mixture = new Mixture(*(mixt_histo.mixture) , false);
  }
  else {
    mixture = 0;
  }

  nb_component = mixt_histo.nb_component;

  weight = new Histogram(*(mixt_histo.weight));

  component = new Histogram*[nb_component];
  for (i = 0;i < nb_component;i++) {
    component[i] = new Histogram(*(mixt_histo.component[i]));
  }
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet Mixture_data.
 *
 *--------------------------------------------------------------*/

void Mixture_data::remove()

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
 *  Destructeur de la classe Mixture_data.
 *
 *--------------------------------------------------------------*/

Mixture_data::~Mixture_data()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Mixture_data.
 *
 *  argument : reference sur un objet Mixture_data.
 *
 *--------------------------------------------------------------*/

Mixture_data& Mixture_data::operator=(const Mixture_data &mixt_histo)

{
  if (&mixt_histo != this) {
    remove();
    delete [] frequency;

    Histogram::copy(mixt_histo);
    copy(mixt_histo);
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

Distribution_data* Mixture_data::extract(Format_error &error , int index) const

{
  bool status = true;
  Distribution_data *pcomponent;


  pcomponent = 0;
  error.init();

  if ((index < 1) || (index > nb_component)) {
    status = false;
    error.update(STAT_error[STATR_HISTOGRAM_INDEX]);
  }
  else {
    index--;
    if (component[index]->nb_element == 0) {
      status = false;
      error.update(STAT_error[STATR_EMPTY_HISTOGRAM]);
    }
  }

  if (status) {
    pcomponent = new Distribution_data(*component[index] ,
                                       (mixture ? mixture->component[index] : 0));
  }

  return pcomponent;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet Mixture_data.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& Mixture_data::line_write(ostream &os) const

{
  os << STAT_label[STATL_SAMPLE_SIZE] << ": " << nb_element << "   "
     << STAT_label[STATL_MEAN] << ": " << mean << "   "
     << STAT_label[STATL_VARIANCE] << ": " << variance;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Mixture_data.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Mixture_data::ascii_write(ostream &os , bool exhaustive) const

{
  if (mixture) {
    mixture->ascii_write(os , this , exhaustive , false);
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Mixture_data dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Mixture_data::ascii_write(Format_error &error , const char *path ,
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
 *  Ecriture d'un objet Mixture_data dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

bool Mixture_data::spreadsheet_write(Format_error &error , const char *path) const

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
 *  Sortie Gnuplot d'un objet Mixture_data.
 *
 *  arguments : reference sur un objet Format_error, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Mixture_data::plot_write(Format_error &error , const char *prefix ,
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
 *  Sortie graphique d'un objet Mixture_data.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* Mixture_data::get_plotable() const

{
  MultiPlotSet *set;


  if (mixture) {
    set = mixture->get_plotable(this);
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

/* RWDEFINE_COLLECTABLE(Mixture_data , STATI_MIXTURE_DATA);


RWspace Mixture_data::binaryStoreSize() const

{
  register int i;
  RWspace size;


  size = Histogram::binaryStoreSize() + sizeof(nb_component) + weight->binaryStoreSize();
  for (i = 0;i < nb_component;i++) {
    size += component[i]->binaryStoreSize();
  }

  if (mixture) {
    size += mixture->recursiveStoreSize();
  }

  return size;
}


void Mixture_data::restoreGuts(RWvistream &is)

{
  register int i;


  remove();

  Histogram::restoreGuts(is);

  is >> nb_component;

  weight = new Histogram();
  weight->restoreGuts(is);

  component = new Histogram*[nb_component];
  for (i = 0;i < nb_component;i++) {
    component[i] = new Histogram();
    component[i]->restoreGuts(is);
  }

  is >> mixture;
  if (mixture == RWnilCollectable) {
    mixture = 0;
  }
}


void Mixture_data::restoreGuts(RWFile &file)

{
  register int i;


  remove();

  Histogram::restoreGuts(file);

  file.Read(nb_component);

  weight = new Histogram();
  weight->restoreGuts(file);

  component = new Histogram*[nb_component];
  for (i = 0;i < nb_component;i++) {
    component[i] = new Histogram();
    component[i]->restoreGuts(file);
  }

  file >> mixture;
  if (mixture == RWnilCollectable) {
    mixture = 0;
  }
}


void Mixture_data::saveGuts(RWvostream &os) const

{
  register int i;


  Histogram::saveGuts(os);

  os << nb_component;
  weight->saveGuts(os);

  for (i = 0;i < nb_component;i++) {
    component[i]->saveGuts(os);
  }

  if (mixture) {
    os << mixture;
  }
  else {
    os << RWnilCollectable;
  }
}


void Mixture_data::saveGuts(RWFile &file) const

{
  register int i;


  Histogram::saveGuts(file);

  file.Write(nb_component);
  weight->saveGuts(file);

  for (i = 0;i < nb_component;i++) {
    component[i]->saveGuts(file);
  }

  if (mixture) {
    file << mixture;
  }
  else {
    file << RWnilCollectable;
  }
} */
