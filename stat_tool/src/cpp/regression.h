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



#ifndef REGRESSION_H
#define REGRESSION_H


#include "vectors.h"


namespace stat_tool {



/****************************************************************
 *
 *  Constants
 */


  const int REGRESSION_NB_VECTOR = 10000;  // maximum number of individuals for the nonparametric regression
  const int NEIGHBORHOOD = 3;            // minimum neighborhood on the values of the explanatory variable

  enum parametric_function {
    LINEAR_FUNCTION ,
    LOGISTIC ,
    MONOMOLECULAR ,
    NONPARAMETRIC_FUNCTION ,
    CONSTANT_FUNCTION
  };



/****************************************************************
 *
 *  Class definition
 */


  /// \brief Regression kernel class

  class RegressionKernel {

  public :

    parametric_function ident;  ///< identifier of the regression function
    int min_value;          ///< minimum value
    int max_value;          ///< maximum value
    double regression_df;   ///< degrees of freedom regression
    double residual_df;     ///< degrees of freedom residuals
    int nb_parameter;       ///< number of parameters
    double *parameter;      ///< parameters
//    double step;             step for representing the regression function
    double *point;          ///< points

    void copy(const RegressionKernel&);
    void remove();

    RegressionKernel();
//    RegressionKernel(parametric_function iident , double imin_value , double imax_value , double istep = 1);
    RegressionKernel(parametric_function iident , int imin_value , int imax_value);
    RegressionKernel(const RegressionKernel &regression) { copy(regression); }
    ~RegressionKernel();
    RegressionKernel& operator=(const RegressionKernel &regression);

    std::ostream& ascii_parameter_print(std::ostream &os) const;
    std::ostream& ascii_formal_print(std::ostream &os) const;
    std::ostream& ascii_print(std::ostream &os) const;
    std::ostream& spreadsheet_print(std::ostream &os) const;
    bool plot_print(const char *path) const;
    void plotable_write(SinglePlot &plot) const;

    void computation();
    double min_computation() const;
    double max_computation() const;
  };


  /// \brief Regression function

  class Regression : public StatInterface , public RegressionKernel {

    friend class Vectors;

    friend std::ostream& operator<<(std::ostream &os , const Regression &regression)
    { return regression.ascii_write(os); }

  private :

    Vectors *vectors;       ///< pointer on a Vectors object
    int nb_vector;          ///< number of individuals
    double *residual;       ///< residuals

    void copy(const Regression&);
    void remove();

    double regression_square_sum_computation() const;
    void residual_computation();
    double residual_mean_computation() const;
    double residual_variance_computation(double residual_mean) const;
    double residual_square_sum_computation() const;
    
  public :

    Regression();
    Regression(parametric_function iident , int explanatory_variable ,
               int response_variable , const Vectors &vec);
    Regression(const Regression &regression);
    ~Regression();
    Regression& operator=(const Regression &regression);

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const std::string path , bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const std::string path) const;
    bool plot_write(StatError &error , const char *prefix , const char *title = NULL) const;
    MultiPlotSet* get_plotable() const;

    // class member access

    Vectors* get_vectors() const { return vectors; }
    int get_nb_vector() const { return nb_vector; }
    double get_residual(int ivec) const { return residual[ivec]; }
  };


};  // namespace stat_tool



#endif
