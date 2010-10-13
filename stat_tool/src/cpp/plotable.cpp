#include <iostream>
#include <stdlib.h>

#include "plotable.h"

using namespace plotable;
using namespace std;



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


MultiPlotSet::MultiPlotSet(int size, int inb_variable)
:multiplots(size), variable_nb_viewpoint(inb_variable) , variable(size), viewpoint(size)
{
  nb_variable = inb_variable;
}
