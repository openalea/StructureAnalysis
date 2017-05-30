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



#include <string.h>
#include <string>
#include <sstream>

#include "stat_tools.h"
#include "stat_label.h"

using namespace std;


namespace stat_tool {



/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the StatError class.
 *
 *  \param[in] imax_nb_error maximum number of errors.
 */
/*--------------------------------------------------------------*/

StatError::StatError(int imax_nb_error)

{
  int i;


  nb_error = 0;
  max_nb_error = imax_nb_error;

  line = new int[max_nb_error];
  column = new int[max_nb_error];

  label = new char*[max_nb_error];
  for (i = 0;i < max_nb_error;i++) {
    label[i] = new char[ERROR_LENGTH];
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the StatError class.
 */
/*--------------------------------------------------------------*/

StatError::~StatError()

{
  int i;


  delete [] line;
  delete [] column;

  if (label) {
    for(i = 0;i < max_nb_error;i++) {
      delete [] label[i];
    }
    delete [] label;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a StatError object.
 *
 *  \param[in,out] os   stream,
 *  \param[in]     type type (ERROR/WARNING).
 */
/*--------------------------------------------------------------*/

ostream& StatError::ascii_write(ostream &os , error_type type) const

{
  int i;


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


/*--------------------------------------------------------------*/
/**
 *  \brief Update of a StatError object.
 *
 *  \param[in] ilabel  label,
 *  \param[in] iline   line,
 *  \param[in] icolumn column.
 */
/*--------------------------------------------------------------*/

void StatError::update(const char *ilabel , int iline , int icolumn)

{
  if (nb_error < max_nb_error) {
    line[nb_error] = iline;
    column[nb_error] = icolumn;
    strcpy(label[nb_error] , ilabel);
    nb_error++;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Update of a StatError object.
 *
 *  \param[in] ilabel     label,
 *  \param[in] correction correction,
 *  \param[in] iline      line,
 *  \param[in] icolumn    column.
 */
/*--------------------------------------------------------------*/

void StatError::correction_update(const char *ilabel , const char *correction ,
                                  int iline , int icolumn)

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


/*--------------------------------------------------------------*/
/**
 *  \brief Update of a StatError object.
 *
 *  \param[in] ilabel     label,
 *  \param[in] correction correction,
 *  \param[in] iline      line,
 *  \param[in] icolumn    column.
 */
/*--------------------------------------------------------------*/

void StatError::correction_update(const char *ilabel , int correction ,
                                  int iline , int icolumn)

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


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a StructureAnalysis object (for display on the console).
 *
 *  \param[in] exhaustive flag detail level,
 *
 *  \return    string.
 */
/*--------------------------------------------------------------*/

string StatInterface::ascii_write(bool exhaustive) const

{
  ostringstream oss;
      

  ascii_write(oss , exhaustive);

  return oss.str();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a StructureAnalysis object in a file.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] path       file path,
 *  \param[in] exhaustive flag detail level.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

/* bool StatInterface::ascii_write(StatError &error , const string path ,
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
    ascii_write(out_file , exhaustive);
  }

  return status;
} */


};  // namespace stat_tool
