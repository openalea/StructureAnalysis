/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       StructureAnalysis: Identifying patterns in plant architecture and development
 *
 *       Copyright 1995-2018 CIRAD AGAP
 *
 *       File author(s): Samuel Dufour-Kowalski, Yann Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id: plotable.cpp 16835 2014-06-12 06:30:51Z guedon $
 *
 *       Forum for StructureAnalysis developers:
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


#include <iostream>
#include <cstdlib>

#include "stat_tools.h"
#include "plotable.h"

using namespace std;


namespace stat_tool {


//////////////////////// SinglePlot /////////////////////////////////////////////

SinglePlot::SinglePlot()
{
  legend = "";
  style = "";
}

SinglePlot::SinglePlot(const SinglePlot& plot)
{
    data = plot.data;
    data_and_text = plot.data_and_text;
    legend = plot.legend;
    style = plot.style;
}

void SinglePlot::add_point(float x, float y)
{
  // data.push_back(std::pair<float, float>(x,y));
    add_text(x, y, "");
}

void SinglePlot::add_point(const PlotPoint& p)
{
  data.push_back(p);
}

void SinglePlot::add_text(float x, float y, const string& text)
{
  PlotPoint pt;
  pt.first = x;
  pt.second = y;
  data_and_text.push_back(std::pair<PlotPoint, string>(pt,text));
}

float SinglePlot::get_x(int i)
{
  if (i >= data_and_text.size())
  {
    cerr << "index larger than size of the data array" << std::endl;
    return 0;
  }

  Label pt = data_and_text[i];
  PlotPoint pp;
  pp = pt.first;

  return pp.first;
}

float SinglePlot::get_y(int i)
{
  if (i >= data_and_text.size())
  {
    cerr << "index larger than size of the data array" << std::endl;
    return 0;
  }

  Label pt = data_and_text[i];
  PlotPoint pp;
  pp = pt.first;

  return pp.second;
}

std::string SinglePlot::get_label(int i)
{
  if (i >= data_and_text.size())
  {
    cerr << "index larger than size of the data array" << std::endl;
    return 0;
  }

  Label pt = data_and_text[i];
  return pt.second;
}

int SinglePlot::get_size() { return data_and_text.size(); }
 
int SinglePlot::size()
{ 
  return data_and_text.size(); 
}   

//////////////////////// MultiPlot /////////////////////////////////////////////
MultiPlot::MultiPlot(int size)
{
  plots = std::vector<SinglePlot>(size);
  xtics = 0.;
  ytics = 0.;
  xrange = Range(0., 0.);
  yrange = Range(0., 0.);
  grid = false;
  group = 0;
}

MultiPlot::MultiPlot(const MultiPlot& multiplot)
{
	plots = multiplot.plots;
    title = multiplot.title;

    xtics = multiplot.xtics;
    ytics = multiplot.ytics;

    xrange = multiplot.xrange;
    yrange = multiplot.yrange;

    xlabel = multiplot.xlabel;
    ylabel = multiplot.ylabel;

    grid = multiplot.grid;

    group = multiplot.group;
}

MultiPlot::MultiPlot(const MultiPlot *multiplot)
:MultiPlot(*multiplot)
{}

SinglePlot& MultiPlot::operator[](int index)
  { return plots[index]; };

// Resize the vector
void MultiPlot::resize(int newsize)
{
    plots.resize(newsize);
}

// Return the number of plots
int MultiPlot::size()
{
    return static_cast<int>(plots.size());
}

// Begin iterator
std::vector<SinglePlot>::const_iterator MultiPlot::begin()
{
    return plots.begin();
}

// End iterator
std::vector<SinglePlot>::const_iterator MultiPlot::end()
{
    return plots.end();
}

//////////////////////// MultiPlotSet /////////////////////////////////////////////

/*
MultiPlotSet::MultiPlotSet(int size, int inb_variable)
:multiplots(size), variable_nb_viewpoint(inb_variable) , variable(size), viewpoint(size)
{
  nb_variable = inb_variable;
}
*/

};  // namespace stat_tool
