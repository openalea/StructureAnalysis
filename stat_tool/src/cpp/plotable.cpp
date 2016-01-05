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
 *       $Id: plotable.cpp 16835 2014-06-12 06:30:51Z guedon $
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


#include <iostream>
#include <stdlib.h>

#include "plotable.h"

using namespace std;


namespace stat_tool {


//////////////////////// SinglePlot /////////////////////////////////////////////

void SinglePlot::add_point(float x, float y)
{
  data.push_back(std::pair<float, float>(x,y));
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


//////////////////////// MultiPlotSet /////////////////////////////////////////////

/*
MultiPlotSet::MultiPlotSet(int size, int inb_variable)
:multiplots(size), variable_nb_viewpoint(inb_variable) , variable(size), viewpoint(size)
{
  nb_variable = inb_variable;
}
*/

};  // namespace stat_tool
