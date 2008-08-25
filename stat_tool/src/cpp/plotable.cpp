#include "plotable.h"

using namespace plotable;
using namespace std;

//////////////////////// SinglePlot /////////////////////////////////////////////


void SinglePlot::add_point(float x, float y)
{
  data.push_back(std::pair<float, float>(x,y));
}

void SinglePlot::add_point(const PlotPoint&  p)
{
  data.push_back(p);
}


