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



#ifndef CURVES_H
#define CURVES_H


#include "stat_tools.h"


namespace stat_tool {


/****************************************************************
 *
 *  Constants
 */


  const int MAX_FREQUENCY = 50;          // maximum frequency for smoothing the curves
  const int MAX_RANGE = 2;               // maximum half-width of the smoothing window

  const int PLOT_NB_CURVE = 12;          // maximum number of curves (Gnuplot output)
  const int PLOT_MIN_FREQUENCY = 10;     // minimum frequency for plotting curve points (Gnuplot output)

  enum curve_transformation {
    CURVE_COPY ,
    SMOOTHING
  };



/****************************************************************
 *
 *  Class definition
 */


  /// \brief Family of curves with frequencies

  class Curves {

    friend std::ostream& operator<<(std::ostream& , const Curves&);

  public :

    int nb_curve;           ///< number of curves
    int length;             ///< curve length
    int offset;             ///< index of the beginning of the curves
    int *index_parameter;   ///< explicit index parameters
    int *frequency;         ///< frequency for each index
    double **point;         ///< curves

    void copy(const Curves&);
    void smooth(const Curves &curves , int max_frequency);
    void remove();

    Curves();
    Curves(int inb_curve , int ilength , bool frequency_flag = false ,
           bool index_parameter_flag = false , bool init_flag = true);
    Curves(const Curves &curves , curve_transformation transform = CURVE_COPY ,
           int max_frequency = MAX_FREQUENCY);
    Curves(const Distribution &dist);
    Curves(const FrequencyDistribution &histo);
    ~Curves();
    Curves& operator=(const Curves &curves);

    std::ostream& ascii_print(std::ostream &os , bool file_flag = false ,
                              const Curves *curves = NULL) const;
    std::ostream& spreadsheet_print(std::ostream &os , const Curves *curves = NULL) const;
    int plot_length_computation() const;
    bool plot_print(const char *path , int ilength = I_DEFAULT ,
                    const Curves *curves_0 = NULL , const Curves *curves_1 = NULL) const;
    bool plot_print_standard_residual(const char *path , double *standard_residual = NULL) const;  // sequence_analysis
    void plotable_write(int index , SinglePlot &plot) const;
    void plotable_write(MultiPlot &plot) const;
    void plotable_frequency_write(SinglePlot &plot) const;

    int max_frequency_computation() const;
    int nb_element_computation() const;
    double mean_computation(int index) const;
    double total_square_sum_computation(int index , double mean) const;
  };


};  // namespace stat_tool



#endif
