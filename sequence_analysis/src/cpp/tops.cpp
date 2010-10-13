/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2010 CIRAD/INRIA Virtual Plants
 *
 *       File author(s): Y. Guedon (yann.guedon@cirad.fr)
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



#include <sstream>

#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/distribution.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/stat_label.h"

#include "sequences.h"
#include "tops.h"
#include "sequence_label.h"

using namespace std;


extern char* label(const char *file_name);
extern int* identifier_select(int nb_pattern , int *pattern_identifier , int selected_nb_pattern ,
                              int *selected_identifier , bool keep);



/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe TopParameters.
 *
 *  arguments : probabilite axe porteur, probabilite axes portes,
 *              rapport de rythme, position maximum.
 *
 *--------------------------------------------------------------*/

TopParameters::TopParameters(double iprobability , double iaxillary_probability ,
                             double irhythm_ratio , int imax_position)

{
  tops = NULL;

  probability = iprobability;
  axillary_probability = iaxillary_probability;
  rhythm_ratio = irhythm_ratio;

  max_position = 0;

  if (imax_position == 0) {
    axillary_nb_internode = NULL;
  }
  else {
    axillary_nb_internode_computation(imax_position);
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet TopParameters.
 *
 *  arguments : reference sur un objet TopParameters,
 *              flag copie de l'objet Tops.
 *
 *--------------------------------------------------------------*/

void TopParameters::copy(const TopParameters &parameters , bool data_flag)

{
  register int i;


  if ((data_flag) && (parameters.tops)) {
    tops = new Tops(*(parameters.tops) , false);
  }
  else {
    tops = NULL;
  }

  probability = parameters.probability;
  axillary_probability = parameters.axillary_probability;
  rhythm_ratio = parameters.rhythm_ratio;

  max_position = parameters.max_position;

  axillary_nb_internode = new Distribution*[max_position + 1];
  axillary_nb_internode[0] = NULL;
  for (i = 1;i <= max_position;i++) {
    axillary_nb_internode[i] = new Distribution(*(parameters.axillary_nb_internode[i]));
  }
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet TopParameters.
 *
 *--------------------------------------------------------------*/

void TopParameters::remove()

{
  delete tops;

  if (axillary_nb_internode) {
    register int i;

    for (i = 1;i <= max_position;i++) {
      delete axillary_nb_internode[i];
    }
    delete [] axillary_nb_internode;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe TopParameters.
 *
 *--------------------------------------------------------------*/

TopParameters::~TopParameters()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Extraction de la loi du nombre d'entrenoeuds des axes portes
 *  pour une position donnee.
 *
 *  arguments : reference sur un objet StatError, position.
 *
 *--------------------------------------------------------------*/

DiscreteParametricModel* TopParameters::extract(StatError &error , int position) const

{
  DiscreteParametricModel *dist;


  error.init();

  if ((position < 1) || (position > max_position)) {
    dist = NULL;
    error.update(SEQ_error[SEQR_POSITION]);
  }

  else {
    dist = new DiscreteParametricModel(*axillary_nb_internode[position] ,
                                       (tops ? tops->axillary_nb_internode[position] : NULL));
  }

  return dist;
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe TopParameters.
 *
 *  argument : reference sur un objet TopParameters.
 *
 *--------------------------------------------------------------*/

TopParameters& TopParameters::operator=(const TopParameters &parameters)

{
  if (&parameters != this) {
    remove();
    copy(parameters);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet TopParameters a partir d'un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              position maximum.
 *
 *--------------------------------------------------------------*/

TopParameters* top_parameters_ascii_read(StatError &error , const char *path ,
                                         int max_position)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  bool status , lstatus;
  register int i;
  int line , read_line;
  double fvalue , probability , axillary_probability , rhythm_ratio;
  TopParameters *parameters;
  ifstream in_file(path);


  parameters = NULL;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    line = 0;

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

        // test mot cle TOP_PARAMETERS

        if (i == 0) {
          if (token != SEQ_word[SEQW_TOP_PARAMETERS]) {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEY_WORD] , SEQ_word[SEQW_TOP_PARAMETERS] , line);
          }
        }

        i++;
      }

      if (i > 0) {
        if (i != 1) {
          status = false;
          error.update(STAT_parsing[STATP_FORMAT] , line);
        }
        break;
      }
    }

    read_line = 0;

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

        // test nom du parametre

        case 0 : {
          switch (read_line) {

          case 0 : {
            if (token != STAT_word[STATW_PROBABILITY]) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_PROBABILITY] , line , i + 1);
            }
            break;
          }

          case 1 : {
            if (token != SEQ_word[SEQW_AXILLARY_PROBABILITY]) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , SEQ_word[SEQW_AXILLARY_PROBABILITY] , line , i + 1);
            }
            break;
          }

          case 2 : {
            if (token != SEQ_word[SEQW_RHYTHM_RATIO]) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , SEQ_word[SEQW_RHYTHM_RATIO] , line , i + 1);
            }
            break;
          }
          }

          break;
        }

        // test separateur

        case 1 : {
          if (token != ":") {
            status = false;
            error.update(STAT_parsing[STATP_SEPARATOR] , line , i + 1);
          }
          break;
        }

        // test valeur du parametre

        case 2 : {
          lstatus = locale.stringToNum(token , &fvalue);
          if (lstatus) {
            if (read_line <= 1) {
              if ((fvalue < TOP_MIN_PROBABILITY) || (fvalue > 1.)) {
                lstatus = false;
              }

              else {
                switch (read_line) {
                case 0 :
                  probability = fvalue;
                  break;
                case 1 :
                  axillary_probability = fvalue;
                  break;
                }
              }
            }

            else {
              if ((fvalue < MIN_RHYTHM_RATIO) || (fvalue > 1. / MIN_RHYTHM_RATIO)) {
                lstatus = false;
              }
              else {
                rhythm_ratio = fvalue;
              }
            }
          }

          if (!lstatus) {
            status = false;
            error.update(STAT_parsing[STATP_PARAMETER_VALUE] , line , i + 1);
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

        read_line++;
        if (read_line == 3) {
          break;
        }
      }
    }

    if (read_line < 3) {
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

    if ((max_position < 1) || (max_position > MAX_POSITION)) {
      status = false;
      error.update(SEQ_error[SEQR_MAX_POSITION]);
    }

    if (status) {
      parameters = new TopParameters(probability , axillary_probability ,
                                     rhythm_ratio , max_position);
    }
  }

  return parameters;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet TopParameters.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& TopParameters::line_write(ostream &os) const

{
  os << STAT_word[STATW_PROBABILITY] << " : " << probability << "   "
     << SEQ_word[SEQW_AXILLARY_PROBABILITY] << " : " << axillary_probability << "   "
     << SEQ_word[SEQW_RHYTHM_RATIO] << " : " << rhythm_ratio;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'une cime et de ces parametres dans un fichier.
 *
 *  arguments : stream, pointeur sur un objet Top,
 *              flag niveau de detail, flag fichier.
 *
 *--------------------------------------------------------------*/

ostream& TopParameters::ascii_write(ostream &os , const Tops *itops ,
                                    bool exhaustive , bool file_flag) const

{
  register int i;
  double *scale;
  const Distribution **pdist;


  os << SEQ_word[SEQW_TOP_PARAMETERS] << endl;

  // ecriture des parametres

  os << "\n" << STAT_word[STATW_PROBABILITY] << " : " << probability << endl;
  os << SEQ_word[SEQW_AXILLARY_PROBABILITY] << " : " << axillary_probability << endl;
  os << SEQ_word[SEQW_RHYTHM_RATIO] << " : " << rhythm_ratio << endl;

  if (itops) {

    // ecriture de la loi empirique du nombre d'entrenoeuds de l'axe porteur

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    itops->nb_internode->ascii_characteristic_print(os , exhaustive , file_flag);

    if (exhaustive) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << "   | " << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
      itops->nb_internode->ascii_print(os , file_flag);
    }
  }

  // ecriture des lois et des lois empiriques du nombre d'entrenoeuds des axes portes

  for (i = 1;i <= max_position;i++) {
    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_DISTRIBUTION] << " " << i << endl;
    axillary_nb_internode[i]->ascii_characteristic_print(os , exhaustive , file_flag);

    if ((itops) && (itops->axillary_nb_internode[i])) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i << " - ";
      itops->axillary_nb_internode[i]->ascii_characteristic_print(os , exhaustive , file_flag);

      if (exhaustive) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << "   | " << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i
           << " | " << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_DISTRIBUTION] << " " << i
           << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i
           << " " << STAT_label[STATL_FUNCTION] << " | " << STAT_label[STATL_CUMULATIVE] << " "
           << STAT_label[STATL_DISTRIBUTION] << " " << i << " " << STAT_label[STATL_FUNCTION] << endl;

        axillary_nb_internode[i]->ascii_print(os , file_flag , true , false ,
                                              itops->axillary_nb_internode[i]);
      }
    }
  }

  if (exhaustive) {
    pdist = new const Distribution*[max_position - 1];
    scale = new double[max_position - 1];
    for (i = 2;i <= max_position;i++) {
      pdist[i - 2] = axillary_nb_internode[i];
      scale[i - 2] = 1.;
    }

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << "  ";
    for (i = 1;i <= max_position;i++) {
      os << " | " << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_DISTRIBUTION] << " " << i;
    }
    os << endl;

    axillary_nb_internode[1]->ascii_print(os , max_position - 1 , pdist ,
                                          scale , file_flag , false);
    delete [] pdist;
    delete [] scale;
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet TopParameters.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& TopParameters::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , tops , exhaustive , false);
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet TopParameters dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool TopParameters::ascii_write(StatError &error , const char *path ,
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
    ascii_write(out_file , tops , exhaustive , true);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'une cime et de ces parametres dans un fichier au format tableur.
 *
 *  arguments : stream, pointeur sur un objet Tops.
 *
 *--------------------------------------------------------------*/

ostream& TopParameters::spreadsheet_write(ostream &os , const Tops *itops) const

{
  register int i;
  double *scale;
  const Distribution **pdist;


  os << SEQ_word[SEQW_TOP_PARAMETERS] << endl;

  // ecriture des parametres

  os << "\n" << STAT_word[STATW_PROBABILITY] << "\t" << probability << endl;
  os << SEQ_word[SEQW_AXILLARY_PROBABILITY] << "\t" << axillary_probability << endl;
  os << SEQ_word[SEQW_RHYTHM_RATIO] << "\t" << rhythm_ratio << endl;

  if (itops) {

    // ecriture de la loi empirique du nombre d'entrenoeuds de l'axe porteur

    os << "\n" << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
    itops->nb_internode->spreadsheet_characteristic_print(os , true);

    os << "\n\t" << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
    itops->nb_internode->spreadsheet_print(os);
  }

  // ecriture des lois et des lois empiriques du nombre d'entrenoeuds des axes portes

  for (i = 1;i <= max_position;i++) {
    os << "\n" << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_DISTRIBUTION] << " " << i << endl;
    axillary_nb_internode[i]->spreadsheet_characteristic_print(os , true);

    if ((itops) && (itops->axillary_nb_internode[i])) {
      os << "\n" << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i << "\t";
      itops->axillary_nb_internode[i]->spreadsheet_characteristic_print(os , true);

      os << "\n\t" << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i
         << "\t" << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_DISTRIBUTION] << " " << i
         << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i
         << " " << STAT_label[STATL_FUNCTION] << "\t" << STAT_label[STATL_CUMULATIVE] << " "
         << STAT_label[STATL_DISTRIBUTION] << " " << i << " " << STAT_label[STATL_FUNCTION] << endl;

      axillary_nb_internode[i]->spreadsheet_print(os , true , false , false ,
                                                  itops->axillary_nb_internode[i]);
    }
  }

  pdist = new const Distribution*[max_position - 1];
  scale = new double[max_position - 1];
  for (i = 2;i <= max_position;i++) {
    pdist[i - 2] = axillary_nb_internode[i];
    scale[i - 2] = 1.;
  }

  os << "\n";
  for (i = 1;i <= max_position;i++) {
    os << "\t" << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_DISTRIBUTION] << " " << i;
  }
  os << endl;

  axillary_nb_internode[1]->spreadsheet_print(os , max_position - 1 , pdist ,
                                              scale , false);
  delete [] pdist;
  delete [] scale;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet TopParameters dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool TopParameters::spreadsheet_write(StatError &error , const char *path) const

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
    spreadsheet_write(out_file , tops);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'une cime et de ces parametres.
 *
 *  arguments : prefixe des fichiers, titre des figures,
 *              pointeur sur un objet Tops.
 *
 *--------------------------------------------------------------*/

bool TopParameters::plot_write(const char *prefix , const char *title ,
                               const Tops *itops) const

{
  bool status;
  register int i , j , k;
  int nb_dist , nb_histo;
  double *scale;
  const Distribution **pdist;
  const FrequencyDistribution **phisto;
  ostringstream *data_file_name;


  // ecriture des fichiers de donnees

  nb_dist = MIN(max_position , PLOT_NB_AXILLARY);
  data_file_name = new ostringstream[nb_dist + 1];

  data_file_name[0] << prefix << 0 << ".dat";
  status = axillary_nb_internode[1]->plot_print((data_file_name[0].str()).c_str());

  if (status) {
    for (i = 2;i <= nb_dist;i++) {
      data_file_name[i - 1] << prefix << i - 1 << ".dat";
      axillary_nb_internode[i]->plot_print((data_file_name[i - 1].str()).c_str());
    }

    if (itops) {
      data_file_name[nb_dist] << prefix << nb_dist << ".dat";

      nb_histo = 1;
      for (i = 1;i <= max_position;i++) {
        if (itops->axillary_nb_internode[i]) {
          nb_histo++;
        }
      }

      phisto = new const FrequencyDistribution*[nb_histo];
      pdist = new const Distribution*[nb_histo - 1];
      scale = new double[nb_histo - 1];

      nb_histo = 0;
      phisto[nb_histo++] = itops->nb_internode;

      for (i = 1;i <= max_position;i++) {
        if (itops->axillary_nb_internode[i]) {
          phisto[nb_histo] = itops->axillary_nb_internode[i];
          pdist[nb_histo - 1] = axillary_nb_internode[i];
          scale[nb_histo++ - 1] = itops->axillary_nb_internode[i]->nb_element;
        }
      }

      plot_print((data_file_name[nb_dist].str()).c_str() , nb_histo - 1 , pdist ,
                 scale , 0 , nb_histo , phisto);
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

      if (axillary_nb_internode[nb_dist]->nb_value - 1 < TIC_THRESHOLD) {
        out_file << "set xtics 0,1" << endl;
      }

      out_file << "plot [0:" << axillary_nb_internode[nb_dist]->nb_value - 1 << "] [0:"
               << MIN(1. , axillary_nb_internode[1]->max * YSCALE) << "] ";
      for (j = 1;j <= nb_dist;j++) {
        out_file << "\"" << label((data_file_name[j - 1].str()).c_str()) << "\" title \""
                 << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_DISTRIBUTION]
                 << " "  << j << "\" with linespoints";
        if (j < nb_dist) {
          out_file << ",\\";
        }
        out_file << endl;
      }

      if (axillary_nb_internode[nb_dist]->nb_value - 1 < TIC_THRESHOLD) {
        out_file << "set xtics autofreq" << endl;
      }

      if (itops) {
        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;

        if (itops->nb_internode->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }
        if ((int)(itops->nb_internode->max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics 0,1" << endl;
        }

        out_file << "plot [0:" << itops->nb_internode->nb_value - 1 << "] [0:"
                 << (int)(itops->nb_internode->max * YSCALE) + 1 << "] \""
                 << label((data_file_name[nb_dist].str()).c_str()) << "\" using 1 title \""
                 << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                 << "\" with impulses" << endl;

        if (itops->nb_internode->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }
        if ((int)(itops->nb_internode->max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics autofreq" << endl;
        }

        j = 2;
        for (k = 1;k <= max_position;k++) {
          if (itops->axillary_nb_internode[k]) {
            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;

            if (axillary_nb_internode[k]->nb_value - 1 < TIC_THRESHOLD) {
              out_file << "set xtics 0,1" << endl;
            }

            out_file << "plot [0:" << axillary_nb_internode[k]->nb_value - 1 << "] [0:"
                     << (int)(MAX(itops->axillary_nb_internode[k]->max ,
                                  axillary_nb_internode[k]->max * itops->axillary_nb_internode[k]->nb_element) * YSCALE) + 1
                     << "] \"" << label((data_file_name[nb_dist].str()).c_str()) << "\" using " << j
                     << " title \"" << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                     << " " << k << "\" with impulses,\\" << endl;
            out_file << "\"" << label((data_file_name[nb_dist].str()).c_str()) << "\" using " << nb_histo - 1 + j++
                     << " title \"" << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_DISTRIBUTION]
                     << " " << k << "\" with linespoints" << endl;

            if (axillary_nb_internode[k]->nb_value - 1 < TIC_THRESHOLD) {
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

    if (itops) {
      delete [] phisto;
      delete [] pdist;
      delete [] scale;
    }
  }

  delete [] data_file_name;

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet TopParameters.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool TopParameters::plot_write(StatError &error , const char *prefix ,
                               const char *title) const

{
  bool status = plot_write(prefix , title , tops);

  error.init();

  if (!status) {
    error.update(STAT_error[STATR_FILE_PREFIX]);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'une cime et de ces parametres.
 *
 *  argument : pointeur sur un objet Tops.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* TopParameters::get_plotable(const Tops *itops) const

{
  register int i , j;
  int nb_dist , nb_histo;
  ostringstream legend;
  MultiPlotSet *plot_set;


  nb_dist = MIN(max_position , PLOT_NB_AXILLARY);

  if (itops) {
    nb_histo = 0;
    for (i = 1;i <= max_position;i++) {
      if (itops->axillary_nb_internode[i]) {
        nb_histo++;
      }
    }

    plot_set = new MultiPlotSet(nb_histo + 2);
  }

  else {
    plot_set = new MultiPlotSet(1);
  }

  MultiPlotSet &plot = *plot_set;

  plot.border = "15 lw 0";

  // 1ere vue : lois du nombre d'entrenoeuds des axes portes

  plot[0].xrange = Range(0 , axillary_nb_internode[nb_dist]->nb_value - 1);
  plot[0].yrange = Range(0. , MIN(1. , axillary_nb_internode[1]->max * YSCALE));

  if (axillary_nb_internode[nb_dist]->nb_value - 1 < TIC_THRESHOLD) {
    plot[0].xtics = 1;
  }

  plot[0].resize(nb_dist);

  for (i = 1;i <= nb_dist;i++) {
    legend.str("");
    legend << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_DISTRIBUTION]
           << " " << i;
    plot[0][i - 1].legend = legend.str();

    plot[0][i - 1].style = "linespoints";

    axillary_nb_internode[i]->plotable_mass_write(plot[0][i - 1]);
  }

  if (itops) {

    // vue suivante : loi empirique du nombre d'entrenoeuds de l'axe porteur

    plot[1].xrange = Range(0 , itops->nb_internode->nb_value - 1);
    plot[1].yrange = Range(0 , ceil(itops->nb_internode->max * YSCALE));

    if (itops->nb_internode->nb_value - 1 < TIC_THRESHOLD) {
      plot[1].xtics = 1;
    }
    if (ceil(itops->nb_internode->max * YSCALE) < TIC_THRESHOLD) {
      plot[1].ytics = 1;
    }

    plot[1].resize(1);

    legend.str("");
    legend << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    plot[1][0].legend = legend.str();

    plot[1][0].style = "impulses";

    itops->nb_internode->plotable_frequency_write(plot[1][0]);

    // vues suivantes : ajustement des lois du nombre d'entrenoeuds des axes portes

    i = 2;
    for (j = 1;j <= max_position;j++) {
      if (itops->axillary_nb_internode[j]) {
        plot[i].xrange = Range(0 , axillary_nb_internode[j]->nb_value - 1);
        plot[i].yrange = Range(0 , ceil(MAX(itops->axillary_nb_internode[j]->max ,
                                            axillary_nb_internode[j]->max * itops->axillary_nb_internode[j]->nb_element) * YSCALE));

        if (axillary_nb_internode[j]->nb_value - 1 < TIC_THRESHOLD) {
          plot[i].xtics = 1;
        }

        plot[i].resize(2);

        legend.str("");
        legend << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
               << " " << j;
        plot[i][0].legend = legend.str();

        plot[i][0].style = "impulses";

        itops->axillary_nb_internode[j]->plotable_frequency_write(plot[i][0]);

        legend.str("");
        legend << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_DISTRIBUTION]
               << " " << j;
        plot[i][1].legend = legend.str();

        plot[i][1].style = "linespoints";

        axillary_nb_internode[j]->plotable_mass_write(plot[i][1] , itops->axillary_nb_internode[j]->nb_element);
        i++;
      }
    }
  }

  return plot_set;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet TopParameters.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* TopParameters::get_plotable() const

{
  return get_plotable(tops);
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe Tops.
 *
 *--------------------------------------------------------------*/

Tops::Tops()

{
  top_parameters = NULL;

  nb_internode = NULL;
  max_position = 0;
  axillary_nb_internode = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Tops.
 *
 *  arguments : nombre de cimes, identificateurs des cimes,
 *              nombres de positions des axes porteurs, flag initialisation.
 *
 *--------------------------------------------------------------*/

Tops::Tops(int nb_top , int *iidentifier , int *nb_position , bool init_flag)

{
  int itype[1] = {NB_INTERNODE};


  init(nb_top , iidentifier , nb_position , NULL , POSITION ,
       1 , itype , false , init_flag);

  top_parameters = NULL;

  nb_internode = NULL;
  max_position = 0;
  axillary_nb_internode = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Tops a partir d'un objet Sequences.
 *
 *  argument : reference sur un objet Sequences.
 *
 *--------------------------------------------------------------*/

Tops::Tops(const Sequences &seq)
:Sequences(seq)

{
  top_parameters = NULL;

  max_position = 0;
  build_nb_internode_frequency_distribution();
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Tops.
 *
 *  arguments : reference sur un objet Tops, nombre de cimes,
 *              indices des cimes selectionnees.
 *
 *--------------------------------------------------------------*/

Tops::Tops(const Tops &tops , int inb_sequence , int *index)
:Sequences(tops , inb_sequence , index)

{
  top_parameters = NULL;

  max_position = 0;
  build_nb_internode_frequency_distribution();
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Tops.
 *
 *  arguments : reference sur un objet Tops,
 *              flag copie de l'objet TopParameters.
 *
 *--------------------------------------------------------------*/

void Tops::copy(const Tops &tops , bool model_flag)

{
  register int i;


  if ((model_flag) && (tops.top_parameters)) {
    top_parameters = new TopParameters(*(tops.top_parameters) , false);
  }
  else {
    top_parameters = NULL;
  }

  nb_internode = new FrequencyDistribution(*(tops.nb_internode));

  max_position = tops.max_position;
  axillary_nb_internode = new FrequencyDistribution*[max_position + 1];
  axillary_nb_internode[0] = NULL;

  for (i = 1;i <= max_position;i++) {
    if (tops.axillary_nb_internode[i]) {
      axillary_nb_internode[i] = new FrequencyDistribution(*(tops.axillary_nb_internode[i]));
    }
    else {
      axillary_nb_internode[i] = NULL;
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par copie de la classe Tops.
 *
 *  arguments : reference sur un objet Tops,
 *              flag copie de l'objet TopParameters et flag inversion.
 *
 *--------------------------------------------------------------*/

Tops::Tops(const Tops &tops , bool model_flag , bool reverse_flag)

{
  Sequences::copy(tops , reverse_flag);

  switch (reverse_flag) {

  case false : {
    copy(tops , model_flag);
    break;
  }

  case true : {
    top_parameters = NULL;

    max_position = 0;
    build_nb_internode_frequency_distribution();
    break;
  }
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par fusion de la classe Tops.
 *
 *  arguments : nombre d'objets Tops, pointeur sur les objets Tops.
 *
 *--------------------------------------------------------------*/

Tops::Tops(int nb_sample , const Tops **ptops)

{
  register int i , j , k , m , n;
  int nb_histo , *pposition , *cposition , *plength , *clength;
  const FrequencyDistribution **phisto;


  phisto = new const FrequencyDistribution*[nb_sample];

  // calcul du nombre et de la longueur des cimes

  nb_sequence = 0;
  for (i = 0;i < nb_sample;i++) {
    nb_sequence += ptops[i]->nb_sequence;
  }

  identifier = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    identifier[i] = i + 1;
  }

  length = new int[nb_sequence];

  max_length = 0;
  cumul_length = 0;

  i = 0;
  for (j = 0;j < nb_sample;j++) {
    for (k = 0;k < ptops[j]->nb_sequence;k++) {
      length[i++] = ptops[j]->length[k];
    }

    if (ptops[j]->max_length > max_length) {
      max_length = ptops[j]->max_length;
    }
    cumul_length += ptops[j]->cumul_length;
    phisto[j] = ptops[j]->hlength;
  }

  hlength = new FrequencyDistribution(nb_sample , phisto);

  index_parameter_type = POSITION;

  for (i = 0;i < nb_sample;i++) {
    phisto[i] = ptops[i]->hindex_parameter;
  }
  hindex_parameter = new FrequencyDistribution(nb_sample , phisto);

  index_interval = NULL;

  nb_variable = 1;

  type = new int[nb_variable];
  type[0] = NB_INTERNODE;

  min_value = new double[nb_variable];
  max_value = new double[nb_variable];

  min_value[0] = ptops[0]->min_value[0];
  max_value[0] = ptops[0]->max_value[0];
  for (i = 1;i < nb_sample;i++) {
    if (ptops[i]->min_value[0] < min_value[0]) {
      min_value[0] = ptops[i]->min_value[0];
    }
    if (ptops[i]->max_value[0] > max_value[0]) {
      max_value[0] = ptops[i]->max_value[0];
    }
  }

  marginal_distribution = new FrequencyDistribution*[nb_variable];

  for (i = 0;i < nb_sample;i++) {
    phisto[i] = ptops[i]->marginal_distribution[0];
  }
  marginal_distribution[0] = new FrequencyDistribution(nb_sample , phisto);

  // copie des sequences

  i = 0;
  index_parameter = new int*[nb_sequence];
  for (j = 0;j < nb_sample;j++) {
    for (k = 0;k < ptops[j]->nb_sequence;k++) {
      index_parameter[i] = new int[length[i] + 1];

      pposition = index_parameter[i];
      cposition = ptops[j]->index_parameter[k];
      for (m = 0;m <= length[i];m++) {
        *pposition++ = *cposition++;
      }
      i++;
    }
  }

  i = 0;
  int_sequence = new int**[nb_sequence];
  for (j = 0;j < nb_sample;j++) {
    for (k = 0;k < ptops[j]->nb_sequence;k++) {
      int_sequence[i] = new int*[nb_variable];
      int_sequence[i][0] = new int[length[i]];

      plength = int_sequence[i][0];
      clength = ptops[j]->int_sequence[k][0];
      for (m = 0;m < length[i];m++) {
        *plength++ = *clength++;
      }
      i++;
    }
  }

  top_parameters = NULL;

  // calcul de la position maximum et fusion des lois empiriques
  // du nombre d'entrenoeuds axe porteur/axes portes

  max_position = 0;
  for (i = 0;i < nb_sample;i++) {
    if (ptops[i]->max_position > max_position) {
      max_position = ptops[i]->max_position;
    }
    phisto[i] = ptops[i]->nb_internode;
  }

  nb_internode = new FrequencyDistribution(nb_sample , phisto);

  axillary_nb_internode = new FrequencyDistribution*[max_position + 1];
  axillary_nb_internode[0] = NULL;

  for (i = 1;i <= max_position;i++) {
    nb_histo = 0;
    for (j = 0;j < nb_sample;j++) {
      if ((i <= ptops[j]->max_position) && (ptops[j]->axillary_nb_internode[i])) {
        phisto[nb_histo++] = ptops[j]->axillary_nb_internode[i];
      }
    }

    if (nb_histo > 0) {
      axillary_nb_internode[i] = new FrequencyDistribution(nb_histo , phisto);
    }
    else {
      axillary_nb_internode[i] = NULL;
    }
  }

  delete [] phisto;
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet Tops.
 *
 *--------------------------------------------------------------*/

void Tops::remove()

{
  delete top_parameters;

  delete nb_internode;

  if (axillary_nb_internode) {
    register int i;

    for (i = 1;i <= max_position;i++) {
      delete axillary_nb_internode[i];
    }
    delete [] axillary_nb_internode;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Tops.
 *
 *--------------------------------------------------------------*/

Tops::~Tops()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Tops.
 *
 *  argument : reference sur un objet Tops.
 *
 *--------------------------------------------------------------*/

Tops& Tops::operator=(const Tops &tops)

{
  if (&tops != this) {
    remove();
    Sequences::remove();

    Sequences::copy(tops);
    copy(tops);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Extraction de la loi empirique du nombre d'entrenoeuds des axes portes
 *  pour une position donnee.
 *
 *  arguments : reference sur un objet StatError, position.
 *
 *--------------------------------------------------------------*/

DiscreteDistributionData* Tops::extract(StatError &error , int position) const

{
  bool status = true;
  DiscreteDistributionData *histo;


  histo = NULL;
  error.init();

  if ((position < 1) || (position > max_position)) {
    status = false;
    error.update(SEQ_error[SEQR_POSITION]);
  }
  else if (!axillary_nb_internode[position]) {
    status = false;
    error.update(STAT_error[STATR_EMPTY_SAMPLE]);
  }

  if (status) {
    histo = new DiscreteDistributionData(*axillary_nb_internode[position] ,
                                         (top_parameters ? top_parameters->axillary_nb_internode[position] : NULL));
  }

  return histo;
}


/*--------------------------------------------------------------*
 *
 *  Suppression des premiers entrenoeuds de l'axe porteur.
 *
 *  arguments : reference sur un objet StatError, nombre d'entrenoeuds.
 *
 *--------------------------------------------------------------*/

Tops* Tops::shift(StatError &error , int inb_internode) const

{
  bool status = true;
  register int i , j;
  int *pposition , *cposition , *plength , *clength;
  Tops *tops;


  tops = NULL;
  error.init();

  if (inb_internode < 1) {
    status = false;
    error.update(SEQ_error[SEQR_NB_INTERNODE]);
  }
  else {
    for (i = 0;i < nb_sequence;i++) {
      if (index_parameter[i][0] <= inb_internode) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_TOP] << " " << i + 1 << ": "
                      << SEQ_error[SEQR_NB_INTERNODE];
        error.update((error_message.str()).c_str());
      }
    }
  }

  if (status) {

    // copie des cimes

    tops = new Tops(nb_sequence , identifier , length , false);

    for (i = 0;i < tops->nb_sequence;i++) {
      pposition = tops->index_parameter[i];
      cposition = index_parameter[i];
      for (j = 0;j <= tops->length[i];j++) {
        *pposition++ = *cposition++ - inb_internode;
      }

      plength = tops->int_sequence[i][0];
      clength = int_sequence[i][0];
      for (j = 0;j < tops->length[i];j++) {
        *plength++ = *clength++;
      }
    }

    tops->min_value[0] = min_value[0];
    tops->max_value[0] = max_value[0];
    tops->marginal_distribution[0] = new FrequencyDistribution(*marginal_distribution[0]);

    tops->build_nb_internode_frequency_distribution();
  }

  return tops;
}


/*--------------------------------------------------------------*
 *
 *  Selection de cimes par l'identificateur.
 *
 *  arguments : reference sur un objet StatError, nombre de cimes,
 *              identificateur des cimes, flag pour conserver ou rejeter
 *              les cimes selectionnees.
 *
 *--------------------------------------------------------------*/

Tops* Tops::select_individual(StatError &error , int inb_sequence ,
                              int *iidentifier , bool keep) const

{
  bool status = true , *selected_top;
  register int i , j;
  int max_identifier , *index;
  Tops *tops;


  tops = NULL;
  error.init();

  if ((inb_sequence < 1) || (inb_sequence > (keep ? nb_sequence : nb_sequence - 1))) {
    status = false;
    error.update(SEQ_error[SEQR_NB_TOP]);
  }

  else {
    max_identifier = 1;
    for (i = 0;i < inb_sequence;i++) {
      if (iidentifier[i] > max_identifier) {
        max_identifier = iidentifier[i];
      }
    }

    selected_top = new bool[max_identifier + 1];
    for (i = 0;i <= max_identifier;i++) {
      selected_top[i] = false;
    }

    for (i = 0;i < inb_sequence;i++) {
      for (j = 0;j < nb_sequence;j++) {
        if (iidentifier[i] == identifier[j]) {
          break;
        }
      }

      if (j == nb_sequence) {
        status = false;
        ostringstream error_message;
        error_message << iidentifier[i] << ": " << SEQ_error[SEQR_TOP_IDENTIFIER];
        error.update((error_message.str()).c_str());
      }

      else if (selected_top[iidentifier[i]]) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_TOP] << " " << iidentifier[i] << " "
                      << STAT_error[STATR_ALREADY_SELECTED];
        error.update((error_message.str()).c_str());
      }
      else {
        selected_top[iidentifier[i]] = true;
      }
    }

    delete [] selected_top;
  }

  if (status) {
    index = identifier_select(nb_sequence , identifier , inb_sequence , iidentifier , keep);

    tops = new Tops(*this , (keep ? inb_sequence : nb_sequence - inb_sequence) , index);

    delete [] index;
  }

  return tops;
}


/*--------------------------------------------------------------*
 *
 *  Inversion du sens de parcours de cimes.
 *
 *  argument : reference sur un objet StatError.
 *
 *--------------------------------------------------------------*/

Tops* Tops::reverse(StatError &error) const

{
  bool status = true;
  register int i;
  Tops *tops;


  tops = NULL;
  error.init();

  for (i = 0;i < nb_sequence;i++) {
    if (index_parameter[i][length[i]] == index_parameter[i][length[i] - 1]) {
      status = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_TOP] << " " << i + 1 << ": "
                    << SEQ_error[SEQR_MAIN_AXE_NB_INTERNODE];
      error.update((error_message.str()).c_str());
    }
  }

  if (status) {
    tops = new Tops(*this , false , true);
  }

  return tops;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Tops a partir d'un fichier.
 *
 *  arguments : reference sur un objet StatError, path, flag format.
 *
 *--------------------------------------------------------------*/

Tops* tops_ascii_read(StatError &error , const char *path , bool old_format)

{
  Sequences *seq;
  Tops *tops;


  tops = NULL;

  seq = sequences_ascii_read(error , path , old_format);

  if (seq) {
    tops = seq->tops(error);
    delete seq;
  }

  return tops;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet Tops.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& Tops::line_write(ostream &os) const

{
  os << nb_sequence << " " << SEQ_label[SEQL_TOPS] << " - "
     << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "   "
     << STAT_label[STATL_MEAN] << ": " << nb_internode->mean << "   "
     << STAT_label[STATL_VARIANCE] << ": " << nb_internode->variance;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Tops.
 *
 *  arguments : stream, flag niveau de detail, flag commentaire.
 *
 *--------------------------------------------------------------*/

ostream& Tops::ascii_write(ostream &os , bool exhaustive , bool comment_flag) const

{
  register int i;


  // ecriture de la loi empirique du nombre d'entrenoeuds de l'axe porteur

  if (comment_flag) {
    os << "# ";
  }
  os << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
  nb_internode->ascii_characteristic_print(os , exhaustive , comment_flag);

  if (exhaustive) {
    os << "\n";
    if (comment_flag) {
      os << "# ";
    }
    os << "   | " << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
    nb_internode->ascii_print(os , comment_flag);
  }

  // ecriture des lois empiriques du nombre d'entrenoeuds des axes portes

  for (i = 1;i <= max_position;i++) {
    if (axillary_nb_internode[i]) {
      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i << " - ";
      axillary_nb_internode[i]->ascii_characteristic_print(os , exhaustive , comment_flag);

      if (exhaustive) {
        os << "\n";
        if (comment_flag) {
          os << "# ";
        }
        os << "   | " << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i << endl;

        axillary_nb_internode[i]->ascii_print(os , comment_flag);
      }
    }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Tops.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Tops::ascii_write(ostream &os , bool exhaustive) const

{
  if (top_parameters) {
    top_parameters->ascii_write(os , this , exhaustive , false);
  }
  else {
    ascii_write(os , exhaustive , false);
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Tops dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Tops::ascii_write(StatError &error , const char *path , bool exhaustive) const

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
    if (top_parameters) {
      top_parameters->ascii_write(out_file , this , exhaustive , true);
    }
    else {
      ascii_write(out_file , exhaustive , false);
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Tops.
 *
 *  arguments : stream, format (lignes/colonnes), flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Tops::ascii_data_write(ostream &os , char format , bool exhaustive) const

{
  register int i;


  os << SEQ_word[SEQW_INDEX_PARAMETER] << " : "
     << SEQ_index_parameter_word[index_parameter_type];

  os << "\n" << nb_variable << " " << STAT_word[STATW_VARIABLE] << endl;

  os << "\n" << STAT_word[STATW_VARIABLE] << " " << 1 << " : "
     << STAT_variable_word[type[0]] << endl;
  os << "\n";

  ascii_write(os , exhaustive , false);
  ascii_print(os , format , false);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Tops dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              format (lignes/colonnes), flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Tops::ascii_data_write(StatError &error , const char *path ,
                            char format , bool exhaustive) const

{
  bool status = false;
  register int i;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    out_file << nb_variable << " " << STAT_word[STATW_VARIABLES] << endl;

    out_file << "\n" << SEQ_word[SEQW_INDEX_PARAMETER] << " : "
             << SEQ_index_parameter_word[index_parameter_type]
             << STAT_word[STATW_VARIABLE] << " " << 1 << " : "
             << STAT_variable_word[type[0]] << endl;
    out_file << "\n";

    ascii_write(out_file , exhaustive , true);
    ascii_print(out_file , format , true);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Tops dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool Tops::spreadsheet_write(StatError &error , const char *path) const

{
  bool status;
  register int i;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    if (top_parameters) {
      top_parameters->spreadsheet_write(out_file , this);
    }

    else {

      // ecriture de la loi empirique du nombre d'entrenoeuds de l'axe porteur

      out_file << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
      nb_internode->spreadsheet_characteristic_print(out_file , true);

      out_file << "\n\t" << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
      nb_internode->spreadsheet_print(out_file);

      // ecriture des lois et des lois empiriques du nombre d'entrenoeuds des axes portes

      for (i = 1;i <= max_position;i++) {
        if (axillary_nb_internode[i]) {
          out_file << "\n" << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i << "\t";
          axillary_nb_internode[i]->spreadsheet_characteristic_print(out_file , true);

          out_file << "\n\t" << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i << endl;
          axillary_nb_internode[i]->spreadsheet_print(out_file);
        }
      }
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet Tops.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Tops::plot_write(StatError &error , const char *prefix ,
                      const char *title) const

{
  bool status;


  error.init();

  if (top_parameters) {
    status = top_parameters->plot_write(prefix , title , this);
  }

  else {
    register int i , j , k;
    int nb_histo;
    const FrequencyDistribution **phisto;
    ostringstream data_file_name;


    // ecriture du fichier de donnees

    data_file_name << prefix << ".dat";

    nb_histo = 0;
    for (i = 1;i <= max_position;i++) {
      if (axillary_nb_internode[i]) {
        nb_histo++;
      }
    }

    phisto = new const FrequencyDistribution*[nb_histo];

    nb_histo = 0;
    for (i = 1;i <= max_position;i++) {
      if (axillary_nb_internode[i]) {
        phisto[nb_histo++] = axillary_nb_internode[i];
      }
    }

    status = nb_internode->plot_print((data_file_name.str()).c_str() , nb_histo , phisto);

    if (status) {

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

        if (nb_internode->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }
        if ((int)(nb_internode->max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics 0,1" << endl;
        }

        out_file << "plot [0:" << nb_internode->nb_value - 1 << "] [0:"
                 << (int)(nb_internode->max * YSCALE) + 1 << "] \""
                 << label((data_file_name.str()).c_str()) << "\" using 1 title \""
                 << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                 << "\" with impulses" << endl;

        if (nb_internode->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }
        if ((int)(nb_internode->max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics autofreq" << endl;
        }

        j = 2;
        for (k = 1;k <= max_position;k++) {
          if (axillary_nb_internode[k]) {
            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;

            if (axillary_nb_internode[k]->nb_value - 1 < TIC_THRESHOLD) {
              out_file << "set xtics 0,1" << endl;
            }
            if ((int)(axillary_nb_internode[k]->max * YSCALE) + 1 < TIC_THRESHOLD) {
              out_file << "set ytics 0,1" << endl;
            }

            out_file << "plot [0:" << axillary_nb_internode[k]->nb_value - 1 << "] [0:"
                     << (int)(axillary_nb_internode[k]->max * YSCALE) + 1 << "] \""
                     << label((data_file_name.str()).c_str()) << "\" using " << j++
                     << " title \"" << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                     << " " << k << "\" with impulses" << endl;

            if (axillary_nb_internode[k]->nb_value - 1 < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }
            if ((int)(axillary_nb_internode[k]->max * YSCALE) + 1 < TIC_THRESHOLD) {
              out_file << "set ytics autofreq" << endl;
            }
          }
        }

        if (i == 1) {
          out_file << "\nset terminal x11" << endl;
        }

        out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
      }
    }

    delete [] phisto;
  }

  if (!status) {
    error.update(STAT_error[STATR_FILE_PREFIX]);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet Tops.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* Tops::get_plotable() const

{
  MultiPlotSet *plot_set;


  if (top_parameters) {
    plot_set = top_parameters->get_plotable(this);
  }

  else {
    register int i , j;
    int nb_histo;
    const FrequencyDistribution *phisto[2] , **merged_histo;
    ostringstream legend;


    nb_histo = 0;
    for (i = 1;i <= max_position;i++) {
      if (axillary_nb_internode[i]) {
        nb_histo++;
      }
    }

    plot_set = new MultiPlotSet(nb_histo + 2);
    MultiPlotSet &plot = *plot_set;

    plot.border = "15 lw 0";

    // 1ere vue : loi empirique du nombre d'entrenoeuds de l'axe porteur

    plot[0].xrange = Range(0 , nb_internode->nb_value - 1);
    plot[0].yrange = Range(0 , ceil(nb_internode->max * YSCALE));

    if (nb_internode->nb_value - 1 < TIC_THRESHOLD) {
      plot[0].xtics = 1;
    }
    if (ceil(nb_internode->max * YSCALE) < TIC_THRESHOLD) {
      plot[0].ytics = 1;
    }

    plot[0].resize(1);

    legend.str("");
    legend << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    plot[0][0].legend = legend.str();

    plot[0][0].style = "impulses";

    nb_internode->plotable_frequency_write(plot[0][0]);

    // 2eme vue : lois empiriques du nombre d'entrenoeuds des axes portes superposees

    merged_histo = new const FrequencyDistribution*[nb_histo];

    i = nb_histo - 1;
    for (j = max_position;j >= 1;j--) {
      if (axillary_nb_internode[j]) {
        if (i == nb_histo - 1) {
          merged_histo[i] = new FrequencyDistribution(*axillary_nb_internode[j]);
        }

        else {
          phisto[0] = merged_histo[i + 1];
          phisto[1] = axillary_nb_internode[j];
          merged_histo[i] = new FrequencyDistribution(2 , phisto);
        }

        i--;
      }
    }

    plot[1].xrange = Range(0 , merged_histo[0]->nb_value - 1);
    plot[1].yrange = Range(0 , ceil(merged_histo[0]->max * YSCALE));

    if (merged_histo[0]->nb_value - 1 < TIC_THRESHOLD) {
      plot[1].xtics = 1;
    }
    if (ceil(merged_histo[0]->max * YSCALE) < TIC_THRESHOLD) {
      plot[1].ytics = 1;
    }

    plot[1].resize(nb_histo);

    i = 0;
    for (j = 1;j <= max_position;j++) {
      if (axillary_nb_internode[j]) {
        legend.str("");
        legend << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
               << " " << j;
        plot[1][i].legend = legend.str();

        plot[1][i].style = "impulses";

        merged_histo[i]->plotable_frequency_write(plot[1][i]);
        i++;
      }
    }

    for (i = 0;i < nb_histo;i++) {
      delete merged_histo[i];
    }
    delete [] merged_histo;

    // vues suivantes : lois empiriques du nombre d'entrenoeuds des axes portes

    i = 2;
    for (j = 1;j <= max_position;j++) {
      if (axillary_nb_internode[j]) {
        plot[i].xrange = Range(0 , axillary_nb_internode[j]->nb_value - 1);
        plot[i].yrange = Range(0 , ceil(axillary_nb_internode[j]->max * YSCALE));

        if (axillary_nb_internode[j]->nb_value - 1 < TIC_THRESHOLD) {
          plot[i].xtics = 1;
        }
        if (ceil(axillary_nb_internode[j]->max * YSCALE) < TIC_THRESHOLD) {
          plot[i].ytics = 1;
        }

        plot[i].resize(1);

        legend.str("");
        legend << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
               << " " << j;
        plot[i][0].legend = legend.str();

        plot[i][0].style = "impulses";

        axillary_nb_internode[j]->plotable_frequency_write(plot[i][0]);
        i++;
      }
    }
  }

  return plot_set;
}
