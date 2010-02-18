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



#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"

#include "stat_tools.h"
#include "markovian.h"
#include "stat_label.h"

using namespace std;


extern char* label(const char *file_name);



/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe DiscreteParametricProcess.
 *
 *  argument : nombre d'etats, nombre de valeurs.
 *
 *--------------------------------------------------------------*/

DiscreteParametricProcess::DiscreteParametricProcess(int inb_state , int inb_value)

{
  register int i;


  nb_state = inb_state;
  nb_value = inb_value;

  if (nb_state > 0) {
    observation = new DiscreteParametric*[nb_state];
    if (nb_value > 0) {
      for (i = 0;i < nb_state;i++) {
        observation[i] = new DiscreteParametric(UNIFORM, 0, nb_value , D_DEFAULT , D_DEFAULT);
      }
    }

    else {
      for (i = 0;i < nb_state;i++) {
        observation[i] = NULL;
      }
    }
  }

  else {
    observation = NULL;
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe DiscreteParametricProcess.
 *
 *  arguments : nombre d'etats, lois d'observation.
 *
 *--------------------------------------------------------------*/

DiscreteParametricProcess::DiscreteParametricProcess(int inb_state , DiscreteParametric **pobservation)

{
  register int i;


  nb_state = inb_state;

  nb_value = 0;
  for (i = 0;i < nb_state;i++) {
    if (pobservation[i]->nb_value > nb_value) {
      nb_value = pobservation[i]->nb_value;
    }
  }

  observation = new DiscreteParametric*[nb_state];
  for (i = 0;i < nb_state;i++) {
    observation[i] = new DiscreteParametric(*pobservation[i] , 'c' , nb_value);
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet DiscreteParametricProcess.
 *
 *  argument : reference sur un objet DiscreteParametricProcess.
 *
 *--------------------------------------------------------------*/

void DiscreteParametricProcess::copy(const DiscreteParametricProcess &process)

{
  register int i;


  nb_state = process.nb_state;
  nb_value = process.nb_value;

  observation = new DiscreteParametric*[nb_state];
  for (i = 0;i < nb_state;i++) {
    observation[i] = new DiscreteParametric(*(process.observation[i]) , 'c' , nb_value);
  }
}


/*--------------------------------------------------------------*
 *
 *  Application d'une permutation des etats. La validite de la 
 *  permutation doit etre verifiee par la procedure appelante.
 *
 *--------------------------------------------------------------*/

void DiscreteParametricProcess::state_permutation(int *perm) const

{
  register int i;
  DiscreteParametric **pobservation= new DiscreteParametric*[nb_state];


  for (i= 0;i < nb_state;i++) {
    pobservation[perm[i]]= observation[i];
  }
  for (i= 0;i < nb_state;i++) {
    observation[i]= pobservation[i];
  }
  delete [] pobservation;
  pobservation= NULL;
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet DiscreteParametricProcess avec ajout d'un etat.
 *
 *  arguments : reference sur un objet DiscreteParametricProcess, etat de reference.
 *
 *--------------------------------------------------------------*/

void DiscreteParametricProcess::add_state(const DiscreteParametricProcess &process , int state)

{
  register int i;


  nb_state = process.nb_state + 1;
  nb_value = process.nb_value;

  observation = new DiscreteParametric*[nb_state];
  for (i = 0;i <= state;i++) {
    observation[i] = new DiscreteParametric(*(process.observation[i]) , 'c' , nb_value);
  }
  for (i = state + 1;i < nb_state;i++) {
    observation[i] = new DiscreteParametric(*(process.observation[i - 1]) , 'c' , nb_value);
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par copie de la classe DiscreteParametricProcess.
 *
 *  arguments : reference sur un objet DiscreteParametricProcess,
 *              type de manipulation ('c' : copy, 's' : state), etat de reference.
 *
 *--------------------------------------------------------------*/

DiscreteParametricProcess::DiscreteParametricProcess(const DiscreteParametricProcess &process ,
                                                     char manip , int state)

{
  switch (manip) {
  case 'c' :
    copy(process);
    break;
  case 's' :
    add_state(process , state);
    break;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet DiscreteParametricProcess.
 *
 *--------------------------------------------------------------*/

void DiscreteParametricProcess::remove()

{
  if (observation) {
    register int i;


    for (i = 0;i < nb_state;i++) {
      delete observation[i];
    }
    delete [] observation;

    observation = NULL;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe DiscreteParametricProcess.
 *
 *--------------------------------------------------------------*/

DiscreteParametricProcess::~DiscreteParametricProcess()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe DiscreteParametricProcess.
 *
 *  argument : reference sur un objet DiscreteParametricProcess.
 *
 *--------------------------------------------------------------*/

DiscreteParametricProcess& DiscreteParametricProcess::operator=(const DiscreteParametricProcess &process)

{
  if (&process != this) {
    remove();
    copy(process);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Analyse du format des lois d'observation.
 *
 *  arguments : reference sur un objet StatError, stream,
 *              reference sur l'indice de la ligne lue, nombre d'etats,
 *              seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

DiscreteParametricProcess* observation_parsing(StatError &error , ifstream &in_file ,
                                               int &line , int nb_state , double cumul_threshold)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  bool status = true , lstatus;
  register int i , j;
  long index;
  DiscreteParametric **dist;
  DiscreteParametricProcess *process;


  process = NULL;

  dist = new DiscreteParametric*[nb_state];
  for (i = 0;i < nb_state;i++) {
    dist[i] = NULL;
  }

  for (i = 0;i < nb_state;i++) {
    while (buffer.readLine(in_file , false)) {
      line++;

#     ifdef DEBUG
      cout << line << "  " << buffer << endl;
#     endif

      position = buffer.first('#');
      if (position != RW_NPOS) {
        buffer.remove(position);
      }
      j = 0;

      RWCTokenizer next(buffer);

      while (!((token = next()).isNull())) {
        switch (j) {

        // test mot cle STATE

        case 0 : {
          if (token != STAT_word[STATW_STATE]) {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_STATE] , line , j + 1);
          }
          break;
        }

        // test indice de l'etat

        case 1 : {
          lstatus = locale.stringToNum(token , &index);
          if ((lstatus) && (index != i)) {
            lstatus = false;
          }

          if (!lstatus) {
            status = false;
            error.correction_update(STAT_parsing[STATP_STATE_INDEX] , i , line , j + 1);
          }
          break;
        }

        // test mot cle OBSERVATION_DISTRIBUTION

        case 2 : {
          if (token != STAT_word[STATW_OBSERVATION_DISTRIBUTION]) {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_OBSERVATION_DISTRIBUTION] , line , j + 1);
          }
          break;
        }
        }

        j++;
      }

      if (j > 0) {
        if (j != 3) {
          status = false;
          error.update(STAT_parsing[STATP_FORMAT] , line);
        }

        dist[i] = discrete_parametric_parsing(error , in_file , line ,
                                              UNIFORM , cumul_threshold);
        if (!dist[i]) {
          status = false;
        }

        break;
      }
    }
  }

  if (status) {
    process = new DiscreteParametricProcess(nb_state , dist);
  }

  for (i = 0;i < nb_state;i++) {
    delete dist[i];
  }
  delete [] dist;

  return process;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet DiscreteParametricProcess.
 *
 *  arguments : stream, pointeurs sur les lois d'observation empiriques,
 *              flag niveau de detail, flag fichier.
 *
 *--------------------------------------------------------------*/

ostream& DiscreteParametricProcess::ascii_print(ostream &os , FrequencyDistribution **empirical_observation ,
                                                bool exhaustive , bool file_flag) const

{
  register int i;


  for (i = 0;i < nb_state;i++) {
    os << "\n" << STAT_word[STATW_STATE] << " " << i << " "
       << STAT_word[STATW_OBSERVATION_DISTRIBUTION] << endl;
    observation[i]->ascii_print(os);
    observation[i]->ascii_characteristic_print(os , false , file_flag);

    if (empirical_observation) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_STATE] << " " << i << " " << STAT_label[STATL_OBSERVATION] << " "
         << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
      empirical_observation[i]->ascii_characteristic_print(os , false , file_flag);
    }

    if (exhaustive) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << "  ";
      if (empirical_observation) {
        os << " | " << STAT_label[STATL_STATE] << " " << i << " "
           << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      }
      os << " | " << STAT_label[STATL_STATE] << " " << i << " "
         << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
      if (empirical_observation) {
        os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
           << STAT_label[STATL_FUNCTION];
      }
      os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
         << STAT_label[STATL_FUNCTION] << endl;

      observation[i]->Distribution::ascii_print(os , file_flag , true , false ,
                                                (empirical_observation ? empirical_observation[i] : 0));
    }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet DiscreteParametricProcess au format tableur.
 *
 *  arguments : stream, pointeurs sur les lois d'observation empiriques.
 *
 *--------------------------------------------------------------*/

ostream& DiscreteParametricProcess::spreadsheet_print(ostream &os , FrequencyDistribution **empirical_observation) const

{
  register int i;


  for (i = 0;i < nb_state;i++) {
    os << "\n" << STAT_word[STATW_STATE] << " " << i << "\t"
       << STAT_word[STATW_OBSERVATION_DISTRIBUTION] << endl;
    observation[i]->spreadsheet_print(os);
    observation[i]->spreadsheet_characteristic_print(os);

    if (empirical_observation) {
      os << "\n" << STAT_label[STATL_STATE] << " " << i << " "
         << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
      empirical_observation[i]->spreadsheet_characteristic_print(os);
    }

    os << "\n";
    if (empirical_observation) {
      os << "\t" << STAT_label[STATL_STATE] << " " << i << " "
         << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    }
    os << "\t" << STAT_label[STATL_STATE] << " " << i << " "
       << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
    if (empirical_observation) {
      os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
         << STAT_label[STATL_FUNCTION];
    }
    os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
       << STAT_label[STATL_FUNCTION] << endl;

    observation[i]->Distribution::spreadsheet_print(os , true , false , false ,
                                                    (empirical_observation ? empirical_observation[i] : 0));
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet DiscreteParametricProcess.
 *
 *  arguments : prefixe des fichiers, titre des figures,
 *              indice du processus d'observation,
 *              pointeurs sur les lois d'observation empiriques.
 *
 *--------------------------------------------------------------*/

bool DiscreteParametricProcess::plot_print(const char *prefix , const char *title , int process ,
                                           FrequencyDistribution **empirical_observation) const

{
  bool status;
  register int i , j;
  int *dist_nb_value , *index_dist;
  double *scale;
  const Distribution **pdist;
  const FrequencyDistribution **phisto;
  ostringstream data_file_name;


  // ecriture du fichier de donnees

  data_file_name << prefix << process << ".dat";

  pdist = new const Distribution*[nb_state];
  dist_nb_value = new int[nb_state];
  scale = new double[nb_state];

  if (empirical_observation) {
    phisto = new const FrequencyDistribution*[nb_state];
    index_dist = new int[nb_state];
  }
  else {
    phisto = NULL;
    index_dist = NULL;
  }

  for (i = 0;i < nb_state;i++) {
    pdist[i] = observation[i];
    dist_nb_value[i] = observation[i]->nb_value;

    if ((empirical_observation) && (empirical_observation[i]->nb_element > 0)) {
      phisto[i] = empirical_observation[i];
      index_dist[i] = i;
      scale[i] = phisto[i]->nb_element;
    }
    else {
      scale[i] = 1.;
    }
  }

  status = ::plot_print((data_file_name.str()).c_str() , nb_state , pdist , scale ,
                        dist_nb_value , (empirical_observation ? nb_state : 0) ,
                        phisto , index_dist);

  if (status) {

    // ecriture des fichiers de commandes et des fichiers d'impression

    for (i = 0;i < 2;i++) {
      ostringstream file_name[2];

      switch (i) {
      case 0 :
        file_name[0] << prefix << process << ".plot";
        break;
      case 1 :
        file_name[0] << prefix << process << ".print";
        break;
      }

      ofstream out_file((file_name[0].str()).c_str());

      if (i == 1) {
        out_file << "set terminal postscript" << endl;
        file_name[1] << label(prefix) << process << ".ps";
        out_file << "set output \"" << file_name[1].str() << "\"\n\n";
      }

      out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
               << "set title";
      if (title) {
        out_file << " \"" << title << " - " << STAT_label[STATL_OUTPUT_PROCESS]
                 << " " << process << "\"";
      }
      out_file << "\n\n";

      for (j = 0;j < nb_state;j++) {
        if (observation[j]->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }

        if ((empirical_observation) && (empirical_observation[j]->nb_element > 0)) {
          out_file << "plot [0:" << observation[j]->nb_value - 1 << "] [0:"
                   << (int)(MAX(empirical_observation[j]->max , observation[j]->max * scale[j]) * YSCALE) + 1
                   << "] \"" << label((data_file_name.str()).c_str()) << "\" using " << j + 1
                   << " title \"" << STAT_label[STATL_STATE] << " " << j << " "
                   << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                   << "\" with impulses,\\" << endl;
          out_file << "\"" << label((data_file_name.str()).c_str()) << "\" using " << nb_state + j + 1
                   << " title \"" << STAT_label[STATL_STATE] << " " << j << " "
                   << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
          observation[j]->plot_title_print(out_file);
          out_file << "\" with linespoints" << endl;
        }

        else {
          out_file << "plot [0:" << observation[j]->nb_value - 1 << "] [0:"
                   << MIN(observation[j]->max * YSCALE , 1.) << "] \""
                   << label((data_file_name.str()).c_str()) << "\" using " << j + 1
                   << " title \"" << STAT_label[STATL_STATE] << " " << j << " "
                   << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
          observation[j]->plot_title_print(out_file);
          out_file << "\" with linespoints" << endl;
        }

        if (observation[j]->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }

        if ((i == 0) && (j < nb_state - 1)) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;
      }

      if (i == 1) {
        out_file << "\nset terminal x11" << endl;
      }

      out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
    }
  }

  delete [] pdist;
  delete [] dist_nb_value;
  delete [] scale;
  if (empirical_observation) {
    delete [] phisto;
    delete [] index_dist;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet DiscreteParametricProcess.
 *
 *  arguments : reference sur un objet MultiPlotSet, indice du MultiPlot,
 *              indice du processus d'observation,
 *              pointeurs sur les lois d'observation empiriques.
 *
 *--------------------------------------------------------------*/

void DiscreteParametricProcess::plotable_write(MultiPlotSet &plot , int &index , int process ,
                                               FrequencyDistribution **empirical_observation) const

{
  register int i , j;
  double scale;
  ostringstream title , legend;


  plot.variable_nb_viewpoint[process] = 1;

  for (i = 0;i < nb_state;i++) {

    // vue : ajustement loi d'observation

    plot.variable[index] = process;
//    plot.viewpoint[index] = OBSERVATION;

    title.str("");
    title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
    plot[index].title = title.str();

    plot[index].xrange = Range(0 , observation[i]->nb_value - 1);
    if (observation[i]->nb_value - 1 < TIC_THRESHOLD) {
      plot[index].xtics = 1;
    }

    if ((empirical_observation) && (empirical_observation[i]->nb_element > 0)) {
      scale = empirical_observation[i]->nb_element;
      plot[index].yrange = Range(0 , ceil(MAX(empirical_observation[i]->max ,
                                              observation[i]->max * scale) * YSCALE));

      plot[index].resize(2);

      legend.str("");
      legend << STAT_label[STATL_STATE] << " " << i << " "
             << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      plot[index][0].legend = legend.str();

      plot[index][0].style = "impulses";

      empirical_observation[i]->plotable_frequency_write(plot[index][0]);
      j = 1;
    }

    else {
      scale = 1.;
      plot[index].yrange = Range(0 , MIN(observation[i]->max * YSCALE , 1.));

      plot[index].resize(1);
      j = 0;
    }

    legend.str("");
    legend << STAT_label[STATL_STATE] << " " << i << " "
           << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
    plot[index][j].legend = legend.str();

    plot[index][j].style = "linespoints";

    observation[i]->plotable_mass_write(plot[index][j] , scale);
    index++;
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul du nombre de valeurs du processus d'observation.
 *
 *--------------------------------------------------------------*/

void DiscreteParametricProcess::nb_value_computation()

{
  register int i;


  nb_value = 0;
  for (i = 0;i < nb_state;i++) {
    if (observation[i]->nb_value > nb_value) {
      nb_value = observation[i]->nb_value;
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul du nombre de parametres independants d'un processus
 *  d'observation parametrique.
 *
 *--------------------------------------------------------------*/

int DiscreteParametricProcess::nb_parameter_computation() const

{
  register int i;
  int nb_parameter = 0;


  for (i = 0;i < nb_state;i++) {
    nb_parameter += observation[i]->nb_parameter_computation();
  }

  return nb_parameter;
}


/*--------------------------------------------------------------*
 *
 *  initialisation des lois d'observation par perturbation.
 *
 *--------------------------------------------------------------*/

void DiscreteParametricProcess::init()

{
  register int i , j;
  double noise_proba , *pmass;


  for (i = 0;i < nb_state;i++) {
    noise_proba = NOISE_PROBABILITY * nb_state / nb_value;

    pmass = observation[i]->mass;
    for (j = 0;j < i * nb_value / nb_state;j++) {
      *pmass++ -= noise_proba / (nb_state - 1);
    }
    for (j = i * nb_value / nb_state;j < (i + 1) * nb_value / nb_state;j++) {
      *pmass++ += noise_proba;
    }
    for (j = (i + 1) * nb_value / nb_state;j < nb_value;j++) {
      *pmass++ -= noise_proba / (nb_state - 1);
    }
  }
}
