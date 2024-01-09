/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2015 CIRAD/INRA/Inria Virtual Plants
 *
 *       File author(s): Samuel Dufour-Kowalski, Yann Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id: plotable.h 16076 2014-03-17 14:47:49Z guedon $
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



#ifndef PLOTABLE_H
#define PLOTABLE_H


#include <list>
#include <vector>
#include <utility>
#include <iterator>
#include <string>

using std::string;


namespace stat_tool



{
  typedef std::pair<float, float> PlotPoint;
  typedef std::pair<float, float> Range;
  typedef std::pair<PlotPoint, string> Label;


  // Class SinglePlot: list of (X, Y) pair with properties

  class SinglePlot
  {
    typedef std::list<PlotPoint> Points;
    typedef std::vector<Label> Labels;

  public:

    Points data;
    Labels data_and_text;
    string legend;
    string style;

    // Combination with
    // impulse, linepoints or line
    // linestyle : -- : -. -
    // marker : + , o . s v x > < ^

    string color;
    // b  : blue,
    // g  : green,
    // r  : red,
    // c  : cyan,
    // m  : magenta,
    // y  : yellow,
    // k  : black,
    // w  : white.

    bool label;
    // if True, then we are in the label mode and the variable label is used instead of data.
    
  public:

    SinglePlot() {};

    void add_point(float x, float y);
    void add_point(const PlotPoint& p);
  
    // related to the label mode
    void add_text(float x, float y, const string& p);   
    float get_x(int);
    float get_y(int);
    std::string get_label(int);
    int get_size() { return data_and_text.size(); }
 
    int size()
    { return data.size(); };   

    Points::const_iterator begin()
    { return data.begin(); };

    Points::const_iterator end()
    { return data.end(); };
  };


  // Class MultiPlot: set of SinglePlot

  class MultiPlot
  {
    std::vector<SinglePlot> plots;
    
  public:

    string title;

    float xtics; 
    float ytics;

    Range xrange;
    Range yrange;

    string xlabel;
    string ylabel;

    bool grid;

    int group; // Window id (default is 0)

    MultiPlot(int size = 1)
    :plots(size), xtics(0), ytics(0), 
    xrange(0, 0), yrange(0, 0), 
    grid(false), group(0)
    {};

    SinglePlot& operator[](int index)
    { return plots[index]; };

    void resize(int newsize)
    { plots.resize(newsize); };

    int size()
    { return plots.size(); };

    std::vector<SinglePlot>::const_iterator begin()
    { return plots.begin(); };

    std::vector<SinglePlot>::const_iterator end()
    { return plots.end(); };
  };


  // Class MultiPlotSet: list of MultiPlot

  class MultiPlotSet
  {
    std::vector<MultiPlot> multiplots;

  public:

    string title;
    string border;
    int nb_variable;
    std::vector<int> variable_nb_viewpoint;
    std::vector<int> variable;
    std::vector<int> viewpoint;

    MultiPlotSet(int size = 1)
    :multiplots(size)
    {};

    MultiPlotSet(int size, int inb_variable);

    MultiPlot& operator[](int index)
    { return multiplots[index]; };

    int size()
    { return multiplots.size(); };

    std::vector<MultiPlot>::const_iterator begin()
    { return multiplots.begin(); };

    std::vector<MultiPlot>::const_iterator end()
    { return multiplots.end(); };
  };


};  // namespace stat_tool



#endif
