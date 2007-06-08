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



#include <string.h>
#include <sstream>
// #include <rw/rwfile.h>
#include "stat_tools.h"
#include "stat_label.h"

using namespace std;



/*--------------------------------------------------------------*
 *
 *  Lecture d'un objet du module STAT dans un fichier au format binaire.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

/* RWCollectable* STAT_binary_read(Format_error &error , const char *path)

{
  RWCollectable *obj;
  RWFile in_file(path , "r");


  obj = 0;
  error.init();

  if (!(in_file.isValid())) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    in_file >> obj;
  }

  return obj;
} */


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet STAT_interface dans un fichier au format binaire.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

/* bool STAT_interface::binary_write(Format_error &error , const char *path) const

{
  bool status;
  RWFile out_file(path , "w");


  error.init();

  if (!(out_file.isValid())) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    out_file << *this;
  }

  return status;
} */


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Format_error.
 *
 *  argument : nombre maximum d'erreurs.
 *
 *--------------------------------------------------------------*/

Format_error::Format_error(int imax_nb_error)

{
  register int i;


  nb_error = 0;
  max_nb_error = imax_nb_error;

  line = new int[max_nb_error];
  column = new int[max_nb_error];

  label = new char*[max_nb_error];
  for (i = 0;i < max_nb_error;i++) {
    label[i] = new char[ERROR_LENGTH];
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Format_error.
 *
 *--------------------------------------------------------------*/

Format_error::~Format_error()

{
  register int i;


  delete [] line;
  delete [] column;

  if (label) {
    for(i = 0;i < max_nb_error;i++) {
      delete [] label[i];
    }
    delete [] label;
  }
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Format_error.
 *
 *  arguments : stream, type des erreurs (ERROR / WARNING).
 *
 *--------------------------------------------------------------*/

ostream& Format_error::ascii_write(ostream &os , int type) const

{
  register int i;


  if (nb_error > 0) {
    for (i = 0;i < nb_error;i++) {
      switch (type) {
      case ERROR :
        os << "*** ERROR : ";
        break;
      case WARNING :
        os << "*** WARNING : ";
        break;
      }

      if (line[i] > 0) {
        os << STAT_label[STATL_LINE] << " " << line[i];
        if (column[i] > 0) {
          os << ": " << STAT_label[STATL_WORD] << " " << column[i];
        }
        os << ": ";
      }

      os << label[i] << endl;
    }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Mise a jour d'un objet Format_error.
 *
 *  arguments : label, ligne, colonne.
 *
 *--------------------------------------------------------------*/

void Format_error::update(const char *ilabel , int iline , int icolumn)

{
  if (nb_error < max_nb_error) {
    line[nb_error] = iline;
    column[nb_error] = icolumn;
    strcpy(label[nb_error] , ilabel);
    nb_error++;
  }
}


/*--------------------------------------------------------------*
 *
 *  Mise a jour d'un objet Format_error.
 *
 *  arguments : label, correction, ligne, colonne.
 *
 *--------------------------------------------------------------*/

void Format_error::correction_update(const char *ilabel , const char *correction , int iline , int icolumn)

{
  if (nb_error < max_nb_error) {
    ostringstream blabel;


    line[nb_error] = iline;
    column[nb_error] = icolumn;

    blabel << ilabel << ": should be " << correction;
    strcpy(label[nb_error] , (blabel.str()).c_str());
    nb_error++;
  }
}


/*--------------------------------------------------------------*
 *
 *  Mise a jour d'un objet Format_error.
 *
 *  arguments : label, correction, ligne, colonne.
 *
 *--------------------------------------------------------------*/

void Format_error::correction_update(const char *ilabel , int correction , int iline , int icolumn)

{
  if (nb_error < max_nb_error) {
    ostringstream blabel;


    line[nb_error] = iline;
    column[nb_error] = icolumn;

    blabel << ilabel << ": should be " << correction;
    strcpy(label[nb_error] , (blabel.str()).c_str());
    nb_error++;
  }
}
