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
 *       $Id: compound.cpp 17978 2015-04-23 06:36:40Z guedon $
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

#include "stat_tools.h"
#include "distribution.h"
#include "compound.h"
#include "stat_label.h"

using namespace std;


namespace stat_tool {



/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe Compound.
 *
 *--------------------------------------------------------------*/

Compound::Compound()

{
  compound_data = NULL;
  sum_distribution = NULL;
  distribution = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Compound.
 *
 *  arguments : references sur la loi de la somme et sur la loi elementaire.
 *              seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

Compound::Compound(const DiscreteParametric &sum_dist , const DiscreteParametric &dist ,
                   double cumul_threshold)

{
  compound_data = NULL;

  sum_distribution = new DiscreteParametric(sum_dist , 'n');
  distribution = new DiscreteParametric(dist , 'n');

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

Compound::Compound(const DiscreteParametric &sum_dist , const DiscreteParametric &dist , char type)

{
  compound_data = NULL;

  switch (type) {

  case 's' : {
    sum_distribution = new DiscreteParametric(sum_dist , 'c' , (int)(sum_dist.nb_value * NB_VALUE_COEFF));
    if ((dist.ident == POISSON) || (dist.ident == NEGATIVE_BINOMIAL)) {
      distribution = new DiscreteParametric(dist.ident , dist.inf_bound , dist.sup_bound ,
                                            dist.parameter , dist.probability , COMPOUND_THRESHOLD);
    }
    else {
      distribution = new DiscreteParametric(dist , 'n');
    }
    break;
  }

  case 'e' : {
    sum_distribution = new DiscreteParametric(sum_dist , 'n');
    distribution = new DiscreteParametric(dist , 'c' , (int)(dist.nb_value * NB_VALUE_COEFF));
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
 *              flag copie de l'objet CompoundData.
 *
 *--------------------------------------------------------------*/

void Compound::copy(const Compound &compound , bool data_flag)

{
  if ((data_flag) && (compound.compound_data)) {
    compound_data = new CompoundData(*(compound.compound_data) , false);
  }
  else {
    compound_data = NULL;
  }

  sum_distribution = new DiscreteParametric(*(compound.sum_distribution));
  distribution = new DiscreteParametric(*(compound.distribution));
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
 *  argument : reference sur un objet StatError.
 *
 *--------------------------------------------------------------*/

CompoundData* Compound::extract_data(StatError &error) const

{
  CompoundData *compound_histo;


  error.init();

  if (!compound_data) {
    compound_histo = NULL;
    error.update(STAT_error[STATR_NO_DATA]);
  }

  else {
    compound_histo = new CompoundData(*compound_data);
    compound_histo->compound = new Compound(*this , false);
  }

  return compound_histo;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Compound a partir d'un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

Compound* compound_ascii_read(StatError &error , const char *path ,
                              double cumul_threshold)

{
  RWCString buffer , token;

  size_t position;
  bool status;
  register int i;
  int line , read_line;
  DiscreteParametric *sum_dist , *dist;
  Compound *compound;
  ifstream in_file(path);


  compound = NULL;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    line = 0;
    read_line = 0;

    sum_dist = NULL;
    dist = NULL;

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
          sum_dist = discrete_parametric_parsing(error , in_file , line ,
                                                 NEGATIVE_BINOMIAL , CUMUL_THRESHOLD);
          if (!sum_dist) {
            status = false;
          }
          break;
        }

        case 2 : {
          dist = discrete_parametric_parsing(error , in_file , line ,
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
 *  arguments : stream, pointeur sur un objet CompoundData,
 *              flag niveau de detail, flag fichier.
 *
 *--------------------------------------------------------------*/

ostream& Compound::ascii_write(ostream &os , const CompoundData *compound_histo ,
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
    os << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
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
      os << " | " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    }
    os << " | " << STAT_label[STATL_COMPOUND] << " " << STAT_label[STATL_DISTRIBUTION];
    if (compound_histo) {
      os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
         << " " << STAT_label[STATL_FUNCTION];
    }
    os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_COMPOUND]
       << " " << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION] << endl;
    ascii_print(os , file_flag , true , false , compound_histo);
  }

  os << "\n" << STAT_word[STATW_SUM] << endl;
  sum_distribution->ascii_print(os);
  sum_distribution->ascii_parametric_characteristic_print(os , exhaustive , file_flag);

  if (compound_histo) {
    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_SUM] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    compound_histo->sum_frequency_distribution->ascii_characteristic_print(os , exhaustive , file_flag);
  }

  if (exhaustive) {
    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << "  ";
    if (compound_histo) {
      os << " | " << STAT_label[STATL_SUM] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    }
    os << " | " << STAT_label[STATL_SUM] << " " << STAT_label[STATL_DISTRIBUTION];
    if (compound_histo) {
      os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_SUM]
         << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION];
    }
    os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_SUM]
       << " " << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION] << endl;

    sum_distribution->Distribution::ascii_print(os , file_flag , true , false ,
                                                (compound_histo ? compound_histo->sum_frequency_distribution : NULL));
  }

  os << "\n" << STAT_word[STATW_ELEMENTARY] << endl;
  distribution->ascii_print(os);
  distribution->ascii_parametric_characteristic_print(os , exhaustive , file_flag);

  if (compound_histo) {
    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_ELEMENTARY] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    compound_histo->frequency_distribution->ascii_characteristic_print(os , exhaustive , file_flag);
  }

  if (exhaustive) {
    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << "  ";
    if ((compound_histo) && (compound_histo->frequency_distribution->nb_element > 0)) {
      os << " | " << STAT_label[STATL_ELEMENTARY] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    }
    os << " | " << STAT_label[STATL_ELEMENTARY] << " " << STAT_label[STATL_DISTRIBUTION];
    if ((compound_histo) && (compound_histo->frequency_distribution->nb_element > 0)) {
      os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_ELEMENTARY]
         << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION];
    }
    os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_ELEMENTARY]
       << " " << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION] << endl;

    distribution->Distribution::ascii_print(os , file_flag , true , false ,
                                            (((compound_histo) && (compound_histo->frequency_distribution->nb_element > 0)) ?
                                             compound_histo->frequency_distribution : NULL));
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
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Compound::ascii_write(StatError &error , const char *path ,
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
 *  arguments : stream, pointeur sur une famille de lois empiriques.
 *
 *--------------------------------------------------------------*/

ostream& Compound::spreadsheet_write(ostream &os , const CompoundData *compound_histo) const

{
  os << STAT_word[STATW_COMPOUND] << endl;
  spreadsheet_characteristic_print(os , true);

  if (compound_histo) {
    double likelihood , information;
    Test test(CHI2);


    os << "\n" << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
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
    os << "\t" << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
  }
  os << "\t" << STAT_label[STATL_COMPOUND] << " " << STAT_label[STATL_DISTRIBUTION];
  if (compound_histo) {
    os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
       << " " << STAT_label[STATL_FUNCTION];
  }
  os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_COMPOUND]
     << " " << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION] << endl;

  spreadsheet_print(os , true , false , false , compound_histo);

  os << "\n" << STAT_word[STATW_SUM] << endl;
  sum_distribution->spreadsheet_print(os);
  sum_distribution->spreadsheet_parametric_characteristic_print(os , true);

  if (compound_histo) {
    os << "\n" << STAT_label[STATL_SUM] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
    compound_histo->sum_frequency_distribution->spreadsheet_characteristic_print(os , true);
  }

  os << "\n";
  if (compound_histo) {
    os << "\t" << STAT_label[STATL_SUM] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
  }
  os << "\t" << STAT_label[STATL_SUM] << " " << STAT_label[STATL_DISTRIBUTION];
  if (compound_histo) {
    os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_SUM]
       << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION];
  }
  os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_SUM]
     << " " << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION] << endl;

  sum_distribution->Distribution::spreadsheet_print(os , true , false , false ,
                                                    (compound_histo ? compound_histo->sum_frequency_distribution : NULL));

  os << "\n" << STAT_word[STATW_ELEMENTARY] << endl;
  distribution->spreadsheet_print(os);
  distribution->spreadsheet_parametric_characteristic_print(os , true);

  if (compound_histo) {
    os << "\n" << STAT_label[STATL_ELEMENTARY] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
    compound_histo->frequency_distribution->spreadsheet_characteristic_print(os , true);
  }

  os << "\n";
  if ((compound_histo) && (compound_histo->frequency_distribution->nb_element > 0)) {
    os << "\t" << STAT_label[STATL_ELEMENTARY] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
  }
  os << "\t" << STAT_label[STATL_ELEMENTARY] << " " << STAT_label[STATL_DISTRIBUTION];
  if ((compound_histo) && (compound_histo->frequency_distribution->nb_element > 0)) {
    os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_ELEMENTARY]
       << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION];
  }
  os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_ELEMENTARY]
     << " " << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION] << endl;

  distribution->Distribution::spreadsheet_print(os , true , false , false ,
                                                (((compound_histo) && (compound_histo->frequency_distribution->nb_element > 0)) ?
                                                 compound_histo->frequency_distribution : NULL));

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Compound dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool Compound::spreadsheet_write(StatError &error , const char *path) const

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
 *              pointeur sur une famille de lois empiriques.
 *
 *--------------------------------------------------------------*/

bool Compound::plot_write(const char *prefix , const char *title ,
                          const CompoundData *compound_histo) const

{
  bool status;
  register int i;
  int nb_histo = 0;
  double scale[3];
  const Distribution *pdist[3];
  const FrequencyDistribution *phisto[3];
  ostringstream data_file_name;


  // ecriture du fichier de donnees

  data_file_name << prefix << ".dat";

  pdist[0] = this;
  pdist[1] = sum_distribution;
  pdist[2] = distribution;

  if (compound_histo) {
    phisto[nb_histo++] = compound_histo;
    scale[0] = compound_histo->nb_element;

    phisto[nb_histo++] = compound_histo->sum_frequency_distribution;
    scale[1] = compound_histo->sum_frequency_distribution->nb_element;

    if (compound_histo->frequency_distribution->nb_element > 0) {
      phisto[nb_histo++] = compound_histo->frequency_distribution;
      scale[2] = compound_histo->frequency_distribution->nb_element;
    }
    else {
      scale[2] = 1.;
    }

    status = stat_tool::plot_print((data_file_name.str()).c_str() , 3 , pdist , scale ,
                                   NULL , nb_histo , phisto);
  }

  else {
    scale[0] = 1.;
    scale[1] = 1.;
    scale[2] = 1.;

    status = stat_tool::plot_print((data_file_name.str()).c_str() , 3 , pdist , scale ,
                                   NULL , 0 , NULL);
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
                 << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses,\\" << endl;
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
                 << (int)(MAX(compound_histo->sum_frequency_distribution->max ,
                              sum_distribution->max * compound_histo->sum_frequency_distribution->nb_element) * YSCALE) + 1
                 << "] \"" << label((data_file_name.str()).c_str()) << "\" using 2 title \""
                 << STAT_label[STATL_SUM] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
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

      if ((compound_histo) && (compound_histo->frequency_distribution->nb_element > 0)) {
        out_file << "plot [0:" << distribution->nb_value - 1 << "] [0:"
                 << (int)(MAX(compound_histo->frequency_distribution->max ,
                              distribution->max * compound_histo->frequency_distribution->nb_element) * YSCALE) + 1
                 << "] \"" << label((data_file_name.str()).c_str()) << "\" using 3 title \""
                 << STAT_label[STATL_ELEMENTARY] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
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
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Compound::plot_write(StatError &error , const char *prefix ,
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
 *  Sortie graphique d'une loi composee et de la structure de donnees associee.
 *
 *  argument : pointeur sur un objet CompoundData.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* Compound::get_plotable(const CompoundData *compound_histo) const

{
  register int i , j;
  int xmax;
  double scale = 1.;
  ostringstream title , legend;


  MultiPlotSet *plot_set = new MultiPlotSet(compound_histo ? 4 : 3);
  MultiPlotSet &plot = *plot_set;

  title.str("");
  title << STAT_label[STATL_COMPOUND] << " " << STAT_label[STATL_DISTRIBUTION];
  if (compound_histo) {
    title << " " << STAT_label[STATL_FIT];
  }
  plot.title = title.str();

  plot.border = "15 lw 0";

  // 1ere vue : loi composee

  xmax = nb_value - 1;
  if ((cumul[xmax] > 1. - DOUBLE_ERROR) &&
      (mass[xmax] > PLOT_MASS_THRESHOLD)) {
    xmax++;
  }
  plot[0].xrange = Range(0 , xmax);

  if (nb_value - 1 < TIC_THRESHOLD) {
    plot[0].xtics = 1;
  }

  if (compound_histo) {
    plot[0].yrange = Range(0. , ceil(MAX(compound_histo->max ,
                                         max * compound_histo->nb_element) * YSCALE));

    plot[0].resize(2);

    plot[0][0].legend = STAT_label[STATL_FREQUENCY_DISTRIBUTION];

    plot[0][0].style = "impulses";

    compound_histo->plotable_frequency_write(plot[0][0]);

    scale = compound_histo->nb_element;
    i = 1;
  }

  else {
    plot[0].yrange = Range(0. , MIN(max * YSCALE , 1.));

    plot[0].resize(1);
    i = 0;
  }

  legend.str("");
  legend << STAT_label[STATL_COMPOUND] << " " << STAT_label[STATL_DISTRIBUTION];
  plot[0][i].legend = legend.str();

  plot[0][i].style = "linespoints";

  plotable_mass_write(plot[0][i] , scale);

  if (compound_histo) {

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

    compound_histo->plotable_cumul_write(plot[1][0]);

    legend.str("");
    legend << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_COMPOUND] << " "
           << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION];
    plot[1][1].legend = legend.str();

    plot[1][1].style = "linespoints";

    plotable_cumul_write(plot[1][1]);

    i = 2;
  }

  else {
    i = 1;
  }

  // 2eme vue : loi de la somme

  xmax = sum_distribution->nb_value - 1;
  if ((sum_distribution->cumul[xmax] > 1. - DOUBLE_ERROR) &&
      (sum_distribution->mass[xmax] > PLOT_MASS_THRESHOLD)) {
    xmax++;
  }
  plot[i].xrange = Range(0 , xmax);

  if (sum_distribution->nb_value - 1 < TIC_THRESHOLD) {
    plot[i].xtics = 1;
  }

  if (compound_histo) {
    plot[i].yrange = Range(0. , ceil(MAX(compound_histo->sum_frequency_distribution->max ,
                                         sum_distribution->max * compound_histo->sum_frequency_distribution->nb_element) * YSCALE));

    plot[i].resize(2);

    legend.str("");
    legend << STAT_label[STATL_SUM] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    plot[i][0].legend = legend.str();

    plot[i][0].style = "impulses";

    compound_histo->sum_frequency_distribution->plotable_frequency_write(plot[i][0]);

    scale = compound_histo->sum_frequency_distribution->nb_element;
    j = 1;
  }

  else {
    plot[i].yrange = Range(0. , MIN(sum_distribution->max * YSCALE , 1.));

    plot[i].resize(1);
    j = 0;
  }

  legend.str("");
  legend << STAT_label[STATL_SUM] << " " << STAT_label[STATL_DISTRIBUTION];
  sum_distribution->plot_title_print(legend);
  plot[i][j].legend = legend.str();

  plot[i][j].style = "linespoints";

  sum_distribution->plotable_mass_write(plot[i][j] , scale);

  i++;

  // 3eme vue : loi elementaire

  xmax = distribution->nb_value - 1;
  if ((distribution->cumul[xmax] > 1. - DOUBLE_ERROR) &&
      (distribution->mass[xmax] > PLOT_MASS_THRESHOLD)) {
    xmax++;
  }
  plot[i].xrange = Range(0 , xmax);

  if (distribution->nb_value - 1 < TIC_THRESHOLD) {
    plot[i].xtics = 1;
  }

  if ((compound_histo) && (compound_histo->frequency_distribution->nb_element > 0)) {
    plot[i].yrange = Range(0. , ceil(MAX(compound_histo->frequency_distribution->max ,
                                         distribution->max * compound_histo->frequency_distribution->nb_element) * YSCALE));

    plot[i].resize(2);
    legend.str("");
    legend << STAT_label[STATL_ELEMENTARY] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    plot[i][0].legend = legend.str();

    plot[i][0].style = "impulses";

    compound_histo->frequency_distribution->plotable_frequency_write(plot[i][0]);

    scale = compound_histo->frequency_distribution->nb_element;
    j = 1;
  }

  else {
    plot[i].yrange = Range(0. , MIN(distribution->max * YSCALE , 1.));

    plot[i].resize(1);

    scale = 1.;
    j = 0;
  }

  legend.str("");
  legend << STAT_label[STATL_ELEMENTARY] << " " << STAT_label[STATL_DISTRIBUTION];
  distribution->plot_title_print(legend);
  plot[i][j].legend = legend.str();

  plot[i][j].style = "linespoints";

  distribution->plotable_mass_write(plot[i][j] , scale);

  return plot_set;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet Compound.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* Compound::get_plotable() const

{
  return get_plotable(compound_data);
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe CompoundData.
 *
 *--------------------------------------------------------------*/

CompoundData::CompoundData()

{
  compound = NULL;
  sum_frequency_distribution = NULL;
  frequency_distribution = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe CompoundData.
 *
 *  arguments : reference sur un objet FrequencyDistribution et
 *              sur un objet Compound.
 *
 *--------------------------------------------------------------*/

CompoundData::CompoundData(const FrequencyDistribution &histo , const Compound &icompound)
:FrequencyDistribution(histo)

{
  compound = NULL;

  sum_frequency_distribution = new FrequencyDistribution(icompound.sum_distribution->alloc_nb_value);
  frequency_distribution = new FrequencyDistribution(icompound.distribution->alloc_nb_value);
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe CompoundData.
 *
 *  argument : reference sur un objet Compound.
 *
 *--------------------------------------------------------------*/

CompoundData::CompoundData(const Compound &icompound)
:FrequencyDistribution(icompound)

{
  compound = NULL;

  sum_frequency_distribution = new FrequencyDistribution(*(icompound.sum_distribution));
  frequency_distribution = new FrequencyDistribution(*(icompound.distribution));
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet CompoundData.
 *
 *  arguments : reference sur un objet CompoundData,
 *              flag copie de l'objet Compound.
 *
 *--------------------------------------------------------------*/

void CompoundData::copy(const CompoundData &compound_histo , bool model_flag)

{
  if ((model_flag) && (compound_histo.compound)) {
    compound = new Compound(*(compound_histo.compound) , false);
  }
  else {
    compound = NULL;
  }

  sum_frequency_distribution = new FrequencyDistribution(*(compound_histo.sum_frequency_distribution));
  frequency_distribution = new FrequencyDistribution(*(compound_histo.frequency_distribution));
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe CompoundData.
 *
 *--------------------------------------------------------------*/

CompoundData::~CompoundData()

{
  delete compound;

  delete sum_frequency_distribution;
  delete frequency_distribution;
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe CompoundData.
 *
 *  argument : reference sur un objet CompoundData.
 *
 *--------------------------------------------------------------*/

CompoundData& CompoundData::operator=(const CompoundData &compound_histo)

{
  if (&compound_histo != this) {
    delete compound;

    delete sum_frequency_distribution;
    delete frequency_distribution;

    delete [] frequency;

    FrequencyDistribution::copy(compound_histo);
    copy(compound_histo);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Extraction de la loi empirique de la somme ou
 *  de la loi empirique elementaire.
 *
 *  arguments : reference sur un objet StatError, type de loi empirique.
 *
 *--------------------------------------------------------------*/

DiscreteDistributionData* CompoundData::extract(StatError &error , char type) const

{
  DiscreteDistributionData *phisto;


  error.init();

  switch (type) {

  case 's' : {
    phisto = new DiscreteDistributionData(*sum_frequency_distribution ,
                                          (compound ? compound->sum_distribution : NULL));
    break;
  }

  case 'e' : {
    if (frequency_distribution->nb_element == 0) {
      phisto = NULL;
      error.update(STAT_error[STATR_EMPTY_SAMPLE]);
    }

    else {
      phisto = new DiscreteDistributionData(*frequency_distribution ,
                                            (compound ? compound->distribution : NULL));
    }
    break;
  }

  case 'c' : {
    phisto = new DiscreteDistributionData(*this , compound);
    break;
  }
  }

  return phisto;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet CompoundData.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& CompoundData::line_write(ostream &os) const

{
  os << STAT_label[STATL_SAMPLE_SIZE] << ": " << nb_element << "   "
     << STAT_label[STATL_MEAN] << ": " << mean << "   "
     << STAT_label[STATL_VARIANCE] << ": " << variance;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet CompoundData.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& CompoundData::ascii_write(ostream &os , bool exhaustive) const

{
  if (compound) {
    compound->ascii_write(os , this , exhaustive , false);
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet CompoundData dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool CompoundData::ascii_write(StatError &error , const char *path ,
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
 *  Ecriture d'un objet CompoundData dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool CompoundData::spreadsheet_write(StatError &error , const char *path) const

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
 *  Sortie Gnuplot d'un objet CompoundData.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool CompoundData::plot_write(StatError &error , const char *prefix ,
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
 *  Sortie graphique d'un objet CompoundData.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* CompoundData::get_plotable() const

{
  MultiPlotSet *plot_set;


  if (compound) {
    plot_set = compound->get_plotable(this);
  }
  else {
    plot_set = NULL;
  }

  return plot_set;
}


};  // namespace stat_tool
