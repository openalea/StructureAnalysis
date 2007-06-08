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
 *  Constructeur de la classe Top_parameters.
 *
 *  arguments : probabilite axe porteur, probabilite axes portes,
 *              rapport de rythme, position maximum.
 *
 *--------------------------------------------------------------*/

Top_parameters::Top_parameters(double iprobability , double iaxillary_probability ,
                               double irhythm_ratio , int imax_position)

{
  tops = 0;

  probability = iprobability;
  axillary_probability = iaxillary_probability;
  rhythm_ratio = irhythm_ratio;

  max_position = 0;

  if (imax_position == 0) {
    axillary_nb_internode = 0;
  }
  else {
    axillary_nb_internode_computation(imax_position);
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Top_parameters.
 *
 *  arguments : reference sur un objet Top_parameters,
 *              flag copie de l'objet Tops.
 *
 *--------------------------------------------------------------*/

void Top_parameters::copy(const Top_parameters &parameters , bool data_flag)

{
  register int i;


  if ((data_flag) && (parameters.tops)) {
    tops = new Tops(*(parameters.tops) , false);
  }
  else {
    tops = 0;
  }

  probability = parameters.probability;
  axillary_probability = parameters.axillary_probability;
  rhythm_ratio = parameters.rhythm_ratio;

  max_position = parameters.max_position;

  axillary_nb_internode = new Distribution*[max_position + 1];
  axillary_nb_internode[0] = 0;
  for (i = 1;i <= max_position;i++) {
    axillary_nb_internode[i] = new Distribution(*(parameters.axillary_nb_internode[i]));
  }
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet Top_parameters.
 *
 *--------------------------------------------------------------*/

void Top_parameters::remove()

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
 *  Destructeur de la classe Top_parameters.
 *
 *--------------------------------------------------------------*/

Top_parameters::~Top_parameters()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Extraction de la loi du nombre d'entrenoeuds des axes portes
 *  pour une position donnee.
 *
 *  arguments : reference sur un objet Format_error, position.
 *
 *--------------------------------------------------------------*/

Parametric_model* Top_parameters::extract(Format_error &error , int position) const

{
  Parametric_model *dist;


  error.init();

  if ((position < 1) || (position > max_position)) {
    dist = 0;
    error.update(SEQ_error[SEQR_POSITION]);
  }

  else {
    dist = new Parametric_model(*axillary_nb_internode[position] ,
                                (tops ? tops->axillary_nb_internode[position] : 0));
  }

  return dist;
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Top_parameters.
 *
 *  argument : reference sur un objet Top_parameters.
 *
 *--------------------------------------------------------------*/

Top_parameters& Top_parameters::operator=(const Top_parameters &parameters)

{
  if (&parameters != this) {
    remove();
    copy(parameters);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Top_parameters a partir d'un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              position maximum.
 *
 *--------------------------------------------------------------*/

Top_parameters* top_parameters_ascii_read(Format_error &error , const char *path ,
                                          int max_position)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  bool status , lstatus;
  register int i;
  int line , read_line;
  double fvalue , probability , axillary_probability , rhythm_ratio;
  Top_parameters *parameters;
  ifstream in_file(path);


  parameters = 0;
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
      parameters = new Top_parameters(probability , axillary_probability , rhythm_ratio , max_position);
    }
  }

  return parameters;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet Top_parameters.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& Top_parameters::line_write(ostream &os) const

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

ostream& Top_parameters::ascii_write(ostream &os , const Tops *itops ,
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

    // ecriture de l'histogramme du nombre d'entrenoeuds de l'axe porteur

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
    itops->nb_internode->ascii_characteristic_print(os , exhaustive , file_flag);

    if (exhaustive) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << "   | " << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_HISTOGRAM] << endl;
      itops->nb_internode->ascii_print(os , file_flag);
    }
  }

  // ecriture des lois et des histogrammes du nombre d'entrenoeuds des axes portes

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
      os << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_HISTOGRAM] << " " << i << " - ";
      itops->axillary_nb_internode[i]->ascii_characteristic_print(os , exhaustive , file_flag);

      if (exhaustive) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << "   | " << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_HISTOGRAM] << " " << i
           << " | " << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_DISTRIBUTION] << " " << i
           << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM] << " " << i
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
 *  Ecriture d'un objet Top_parameters.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Top_parameters::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , tops , exhaustive , false);
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Top_parameters dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Top_parameters::ascii_write(Format_error &error , const char *path ,
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

ostream& Top_parameters::spreadsheet_write(ostream &os , const Tops *itops) const

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

    // ecriture de l'histogramme du nombre d'entrenoeuds de l'axe porteur

    os << "\n" << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_HISTOGRAM] << "\t";
    itops->nb_internode->spreadsheet_characteristic_print(os , true);

    os << "\n\t" << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_HISTOGRAM] << endl;
    itops->nb_internode->spreadsheet_print(os);
  }

  // ecriture des lois et des histogrammes du nombre d'entrenoeuds des axes portes

  for (i = 1;i <= max_position;i++) {
    os << "\n" << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_DISTRIBUTION] << " " << i << endl;
    axillary_nb_internode[i]->spreadsheet_characteristic_print(os , true);

    if ((itops) && (itops->axillary_nb_internode[i])) {
      os << "\n" << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_HISTOGRAM] << " " << i << "\t";
      itops->axillary_nb_internode[i]->spreadsheet_characteristic_print(os , true);

      os << "\n\t" << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_HISTOGRAM] << " " << i
         << "\t" << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_DISTRIBUTION] << " " << i
         << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM] << " " << i
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
 *  Ecriture d'un objet Top_parameters dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

bool Top_parameters::spreadsheet_write(Format_error &error , const char *path) const

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

bool Top_parameters::plot_write(const char *prefix , const char *title ,
                                const Tops *itops) const

{
  bool status;
  register int i , j , k;
  int nb_histo , nb_dist , *index_dist;
  double *scale;
  const Distribution **pdist;
  const Histogram **phisto;
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

      phisto = new const Histogram*[nb_histo];
      index_dist = new int[nb_histo];
      pdist = new const Distribution*[nb_histo - 1];
      scale = new double[nb_histo - 1];

      nb_histo = 0;
      phisto[nb_histo] = itops->nb_internode;
      index_dist[nb_histo++] = I_DEFAULT;

      for (i = 1;i <= max_position;i++) {
        if (itops->axillary_nb_internode[i]) {
          phisto[nb_histo] = itops->axillary_nb_internode[i];
          index_dist[nb_histo] = nb_histo - 1;
          pdist[nb_histo - 1] = axillary_nb_internode[i];
          scale[nb_histo++ - 1] = itops->axillary_nb_internode[i]->nb_element;
        }
      }

      plot_print((data_file_name[nb_dist].str()).c_str() , nb_histo - 1 , pdist ,
                 scale , 0 , nb_histo , phisto , index_dist);
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
                 << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_HISTOGRAM]
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
                     << " title \"" << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_HISTOGRAM]
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
      delete [] index_dist;
      delete [] pdist;
      delete [] scale;
    }
  }

  delete [] data_file_name;

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet Top_parameters.
 *
 *  arguments : reference sur un objet Format_error, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Top_parameters::plot_write(Format_error &error , const char *prefix ,
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
 *  Fonctions pour la persistance.
 *
 *--------------------------------------------------------------*/

/* RWDEFINE_COLLECTABLE(Top_parameters , STATI_TOP_PARAMETERS);


RWspace Top_parameters::binaryStoreSize() const

{
  register int i;
  RWspace size;


  size = sizeof(probability) + sizeof(axillary_probability) + sizeof(rhythm_ratio) + sizeof(max_position);

  for (i = 1;i <= max_position;i++) {
    size += axillary_nb_internode[i]->binaryStoreSize();
  }

  if (tops) {
    size += tops->recursiveStoreSize();
  }

  return size;
}


void Top_parameters::restoreGuts(RWvistream &is)

{
  register int i;


  remove();

  is >> probability >> axillary_probability >> rhythm_ratio >> max_position;

  axillary_nb_internode = new Distribution*[max_position + 1];
  axillary_nb_internode[0] = 0;
  for (i = 1;i <= max_position;i++) {
    axillary_nb_internode[i] = new Distribution();
    axillary_nb_internode[i]->restoreGuts(is);
  }

  is >> tops;
  if (tops == RWnilCollectable) {
    tops = 0;
  }
}


void Top_parameters::restoreGuts(RWFile &file)

{
  register int i;


  remove();

  file.Read(probability);
  file.Read(axillary_probability);
  file.Read(rhythm_ratio);
  file.Read(max_position);

  axillary_nb_internode = new Distribution*[max_position + 1];
  axillary_nb_internode[0] = 0;
  for (i = 1;i <= max_position;i++) {
    axillary_nb_internode[i] = new Distribution();
    axillary_nb_internode[i]->restoreGuts(file);
  }

  file >> tops;
  if (tops == RWnilCollectable) {
    tops = 0;
  }
}


void Top_parameters::saveGuts(RWvostream &os) const

{
  register int i;


  os << probability << axillary_probability << rhythm_ratio << max_position;

  for (i = 1;i <= max_position;i++) {
    axillary_nb_internode[i]->saveGuts(os);
  }

  if (tops) {
    os << tops;
  }
  else {
    os << RWnilCollectable;
  }
}


void Top_parameters::saveGuts(RWFile &file) const

{
  register int i;


  file.Write(probability);
  file.Write(axillary_probability);
  file.Write(rhythm_ratio);
  file.Write(max_position);

  for (i = 1;i <= max_position;i++) {
    axillary_nb_internode[i]->saveGuts(file);
  }

  if (tops) {
    file << tops;
  }
  else {
    file << RWnilCollectable;
  }
} */


/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe Tops.
 *
 *--------------------------------------------------------------*/

Tops::Tops()

{
  top_parameters = 0;

  nb_internode = 0;
  max_position = 0;
  axillary_nb_internode = 0;
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
  int itype[2] = {POSITION , NB_INTERNODE};


  init(2 , itype , nb_top , iidentifier , nb_position , init_flag);

  top_parameters = 0;

  nb_internode = 0;
  max_position = 0;
  axillary_nb_internode = 0;
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
  top_parameters = 0;

  max_position = 0;
  build_nb_internode_histogram();
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
  top_parameters = 0;

  max_position = 0;
  build_nb_internode_histogram();
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Tops.
 *
 *  arguments : reference sur un objet Tops,
 *              flag copie de l'objet Top_parameters.
 *
 *--------------------------------------------------------------*/

void Tops::copy(const Tops &tops , bool model_flag)

{
  register int i;


  if ((model_flag) && (tops.top_parameters)) {
    top_parameters = new Top_parameters(*(tops.top_parameters) , false);
  }
  else {
    top_parameters = 0;
  }

  nb_internode = new Histogram(*(tops.nb_internode));

  max_position = tops.max_position;
  axillary_nb_internode = new Histogram*[max_position + 1];
  axillary_nb_internode[0] = 0;

  for (i = 1;i <= max_position;i++) {
    if (tops.axillary_nb_internode[i]) {
      axillary_nb_internode[i] = new Histogram(*(tops.axillary_nb_internode[i]));
    }
    else {
      axillary_nb_internode[i] = 0;
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par copie de la classe Tops.
 *
 *  arguments : reference sur un objet Tops,
 *              flag copie de l'objet Top_parameters et flag inversion.
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
    top_parameters = 0;

    max_position = 0;
    build_nb_internode_histogram();
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
  int nb_histo , blength , *psequence , *csequence;
  const Histogram **phisto;


  phisto = new const Histogram*[nb_sample];

  nb_variable = 2;

  type = new int[nb_variable];
  type[0] = POSITION;
  type[1] = NB_INTERNODE;

  min_value = new int[nb_variable];
  max_value = new int[nb_variable];
  min_value[0] = 0;
  max_value[0] = 0;

  min_value[1] = ptops[0]->min_value[1];
  max_value[1] = ptops[0]->max_value[1];
  for (i = 1;i < nb_sample;i++) {
    if (ptops[i]->min_value[1] < min_value[1]) {
      min_value[1] = ptops[i]->min_value[1];
    }
    if (ptops[i]->max_value[1] > max_value[1]) {
      max_value[1] = ptops[i]->max_value[1];
    }
  }

  marginal = new Histogram*[nb_variable];
  marginal[0] = 0;

  for (i = 0;i < nb_sample;i++) {
    phisto[i] = ptops[i]->marginal[1];
  }
  marginal[1] = new Histogram(nb_sample , phisto);

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

  hlength = new Histogram(nb_sample , phisto);

  index_interval = 0;

  // copie des sequences

  i = 0;
  sequence = new int**[nb_sequence];
  for (j = 0;j < nb_sample;j++) {
    for (k = 0;k < ptops[j]->nb_sequence;k++) {
      sequence[i] = new int*[nb_variable];
      for (m = 0;m < nb_variable;m++) {
        blength = (type[m] == POSITION ? length[i] + 1 :  length[i]);
        sequence[i][m] = new int[blength];

        psequence = sequence[i][m];
        csequence = ptops[j]->sequence[k][m];
        for (n = 0;n < blength;n++) {
          *psequence++ = *csequence++;
        }
      }
      i++;
    }
  }

  top_parameters = 0;

  // calcul de la position maximum et fusion des histogrammes
  // du nombre d'entrenoeuds axe porteur/axes portes

  max_position = 0;
  for (i = 0;i < nb_sample;i++) {
    if (ptops[i]->max_position > max_position) {
      max_position = ptops[i]->max_position;
    }
    phisto[i] = ptops[i]->nb_internode;
  }

  nb_internode = new Histogram(nb_sample , phisto);

  axillary_nb_internode = new Histogram*[max_position + 1];
  axillary_nb_internode[0] = 0;

  for (i = 1;i <= max_position;i++) {
    nb_histo = 0;
    for (j = 0;j < nb_sample;j++) {
      if ((i <= ptops[j]->max_position) && (ptops[j]->axillary_nb_internode[i])) {
        phisto[nb_histo++] = ptops[j]->axillary_nb_internode[i];
      }
    }

    if (nb_histo > 0) {
      axillary_nb_internode[i] = new Histogram(nb_histo , phisto);
    }
    else {
      axillary_nb_internode[i] = 0;
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
 *  Extraction de l'histogramme du nombre d'entrenoeuds des axes portes
 *  pour une position donnee.
 *
 *  arguments : reference sur un objet Format_error, position.
 *
 *--------------------------------------------------------------*/

Distribution_data* Tops::extract(Format_error &error , int position) const

{
  bool status = true;
  Distribution_data *histo;


  histo = 0;
  error.init();

  if ((position < 1) || (position > max_position)) {
    status = false;
    error.update(SEQ_error[SEQR_POSITION]);
  }
  else if (!axillary_nb_internode[position]) {
    status = false;
    error.update(STAT_error[STATR_EMPTY_HISTOGRAM]);
  }

  if (status) {
    histo = new Distribution_data(*axillary_nb_internode[position] ,
                                  (top_parameters ? top_parameters->axillary_nb_internode[position] : 0));
  }

  return histo;
}


/*--------------------------------------------------------------*
 *
 *  Suppression des premiers entrenoeuds de l'axe porteur.
 *
 *  arguments : reference sur un objet Format_error, nombre d'entrenoeuds.
 *
 *--------------------------------------------------------------*/

Tops* Tops::shift(Format_error &error , int inb_internode) const

{
  bool status = true;
  register int i , j;
  int *pposition , *cposition , *plength , *clength;
  Tops *tops;


  tops = 0;
  error.init();

  if (inb_internode < 1) {
    status = false;
    error.update(SEQ_error[SEQR_NB_INTERNODE]);
  }
  else {
    for (i = 0;i < nb_sequence;i++) {
      if (sequence[i][0][0] <= inb_internode) {
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
      pposition = tops->sequence[i][0];
      cposition = sequence[i][0];
      for (j = 0;j <= tops->length[i];j++) {
        *pposition++ = *cposition++ - inb_internode;
      }

      plength = tops->sequence[i][1];
      clength = sequence[i][1];
      for (j = 0;j < tops->length[i];j++) {
        *plength++ = *clength++;
      }
    }

    tops->min_value[1] = min_value[1];
    tops->max_value[1] = max_value[1];
    tops->marginal[1] = new Histogram(*(marginal[1]));

    tops->build_nb_internode_histogram();
  }

  return tops;
}


/*--------------------------------------------------------------*
 *
 *  Selection de cimes par l'identificateur.
 *
 *  arguments : reference sur un objet Format_error, nombre de cimes,
 *              identificateur des cimes, flag pour conserver ou rejeter
 *              les cimes selectionnees.
 *
 *--------------------------------------------------------------*/

Tops* Tops::select_individual(Format_error &error , int inb_sequence ,
                              int *iidentifier , bool keep) const

{
  bool status = true , *selected_top;
  register int i , j;
  int max_identifier , *index;
  Tops *tops;


  tops = 0;
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
 *  argument : reference sur un objet Format_error.
 *
 *--------------------------------------------------------------*/

Tops* Tops::reverse(Format_error &error) const

{
  bool status = true;
  register int i;
  Tops *tops;


  tops = 0;
  error.init();

  for (i = 0;i < nb_sequence;i++) {
    if (sequence[i][0][length[i]] == sequence[i][0][length[i] - 1]) {
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
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

Tops* tops_ascii_read(Format_error &error , const char *path)

{
  Sequences *seq;
  Tops *tops;


  tops = 0;

  seq = sequences_ascii_read(error , path);

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
     << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_HISTOGRAM] << "   "
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


  // ecriture de l'histogramme du nombre d'entrenoeuds de l'axe porteur

  if (comment_flag) {
    os << "# ";
  }
  os << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
  nb_internode->ascii_characteristic_print(os , exhaustive , comment_flag);

  if (exhaustive) {
    os << "\n";
    if (comment_flag) {
      os << "# ";
    }
    os << "   | " << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_HISTOGRAM] << endl;
    nb_internode->ascii_print(os , comment_flag);
  }

  // ecriture des histogrammes du nombre d'entrenoeuds des axes portes

  for (i = 1;i <= max_position;i++) {
    if (axillary_nb_internode[i]) {
      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_HISTOGRAM] << " " << i << " - ";
      axillary_nb_internode[i]->ascii_characteristic_print(os , exhaustive , comment_flag);

      if (exhaustive) {
        os << "\n";
        if (comment_flag) {
          os << "# ";
        }
        os << "   | " << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_HISTOGRAM] << " " << i << endl;

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
 *  arguments : reference sur un objet Format_error, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Tops::ascii_write(Format_error &error , const char *path , bool exhaustive) const

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


  os << nb_variable << " " << STAT_word[STATW_VARIABLES] << endl;

  os << "\n";
  for (i = 0;i < nb_variable;i++) {
    os << STAT_word[STATW_VARIABLE] << " " << i + 1 << " : "
       << STAT_sequence_word[type[i]] << endl;
  }
  os << endl;

  ascii_write(os , exhaustive , false);
  ascii_print(os , format , false);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Tops dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              format (lignes/colonnes), flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Tops::ascii_data_write(Format_error &error , const char *path ,
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

    out_file << "\n";
    for (i = 0;i < nb_variable;i++) {
      out_file << STAT_word[STATW_VARIABLE] << " " << i + 1 << " : "
               << STAT_sequence_word[type[i]] << endl;
    }
    out_file << endl;

    ascii_write(out_file , exhaustive , true);
    ascii_print(out_file , format , true);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Tops dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

bool Tops::spreadsheet_write(Format_error &error , const char *path) const

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

      // ecriture de l'histogramme du nombre d'entrenoeuds de l'axe porteur

      out_file << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_HISTOGRAM] << "\t";
      nb_internode->spreadsheet_characteristic_print(out_file , true);

      out_file << "\n\t" << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_HISTOGRAM] << endl;
      nb_internode->spreadsheet_print(out_file);

      // ecriture des lois et des histogrammes du nombre d'entrenoeuds des axes portes

      for (i = 1;i <= max_position;i++) {
        if (axillary_nb_internode[i]) {
          out_file << "\n" << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_HISTOGRAM] << " " << i << "\t";
          axillary_nb_internode[i]->spreadsheet_characteristic_print(out_file , true);

          out_file << "\n\t" << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_HISTOGRAM] << " " << i << endl;
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
 *  arguments : reference sur un objet Format_error, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Tops::plot_write(Format_error &error , const char *prefix ,
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
    const Histogram **phisto;
    ostringstream data_file_name;


    // ecriture du fichier de donnees

    data_file_name << prefix << ".dat";

    nb_histo = 0;
    for (i = 1;i <= max_position;i++) {
      if (axillary_nb_internode[i]) {
        nb_histo++;
      }
    }

    phisto = new const Histogram*[nb_histo];

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
                 << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_HISTOGRAM]
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
                     << " title \"" << SEQ_label[SEQL_NB_INTERNODE] << " " << STAT_label[STATL_HISTOGRAM]
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
 *  Fonctions pour la persistance.
 *
 *--------------------------------------------------------------*/

/* RWDEFINE_COLLECTABLE(Tops , STATI_TOPS);


RWspace Tops::binaryStoreSize() const

{
  register int i;
  RWspace size;


  size = Sequences::binaryStoreSize() + nb_internode->binaryStoreSize() +
         sizeof(max_position);

  for (i = 1;i <= max_position;i++) {
    size += sizeof(true);
    if (axillary_nb_internode[i]) {
      size += axillary_nb_internode[i]->binaryStoreSize();
    }
  }

  if (top_parameters) {
    size += top_parameters->recursiveStoreSize();
  }

  return size;
}


void Tops::restoreGuts(RWvistream &is)

{
  bool status;
  register int i;


  remove();

  Sequences::restoreGuts(is);

  nb_internode = new Histogram();
  nb_internode->restoreGuts(is);

  is >> max_position;

  axillary_nb_internode = new Histogram*[max_position + 1];
  axillary_nb_internode[0] = 0;
  for (i = 1;i <= max_position;i++) {
    is >> status;
    if (status) {
      axillary_nb_internode[i] = new Histogram();
      axillary_nb_internode[i]->restoreGuts(is);
    }
    else {
      axillary_nb_internode[i] = 0;
    }
  }

  is >> top_parameters;
  if (top_parameters == RWnilCollectable) {
    top_parameters = 0;
  }
}


void Tops::restoreGuts(RWFile &file)

{
  bool status;
  register int i;


  remove();

  Sequences::restoreGuts(file);

  nb_internode = new Histogram();
  nb_internode->restoreGuts(file);

  file.Read(max_position);

  axillary_nb_internode = new Histogram*[max_position + 1];
  axillary_nb_internode[0] = 0;
  for (i = 1;i <= max_position;i++) {
    file.Read(status);
    if (status) {
      axillary_nb_internode[i] = new Histogram();
      axillary_nb_internode[i]->restoreGuts(file);
    }
    else {
      axillary_nb_internode[i] = 0;
    }
  }

  file >> top_parameters;
  if (top_parameters == RWnilCollectable) {
    top_parameters = 0;
  }
}


void Tops::saveGuts(RWvostream &os) const

{
  register int i;


  Sequences::saveGuts(os);

  nb_internode->saveGuts(os);
  os << max_position;

  for (i = 1;i <= max_position;i++) {
    if (axillary_nb_internode[i]) {
      os << true;
      axillary_nb_internode[i]->saveGuts(os);
    }
    else {
      os << false;
    }
  }

  if (top_parameters) {
    os << top_parameters;
  }
  else {
    os << RWnilCollectable;
  }
}


void Tops::saveGuts(RWFile &file) const

{
  register int i;


  Sequences::saveGuts(file);

  nb_internode->saveGuts(file);
  file.Write(max_position);

  for (i = 1;i <= max_position;i++) {
    if (axillary_nb_internode[i]) {
      file.Write(true);
      axillary_nb_internode[i]->saveGuts(file);
    }
    else {
      file.Write(false);
    }
  }

  if (top_parameters) {
    file << top_parameters;
  }
  else {
    file << RWnilCollectable;
  }
} */
