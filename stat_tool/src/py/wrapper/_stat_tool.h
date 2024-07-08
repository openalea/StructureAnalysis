#ifndef AUTOWIG__STAT_TOOL
#define AUTOWIG__STAT_TOOL

#include <boost/python.hpp>
#include <type_traits>
#include <stat_tool/markovian.h>
#include <stat_tool/stat_tools.h>
#include <stat_tool/distance_matrix.h>
#include <stat_tool/stat_label.h>
#include <stat_tool/regression.h>
#include <stat_tool/compound.h>
#include <stat_tool/vectors.h>
#include <stat_tool/distribution.h>
#include <stat_tool/mixture.h>
#include <stat_tool/plotable.h>
#include <stat_tool/curves.h>
#include <stat_tool/discrete_mixture.h>
#include <stat_tool/chain_reestimation.h>
#include <stat_tool/convolution.h>
#include <stat_tool/config.h>
#include <stat_tool/reestimation.h>

namespace autowig
{
     template<class T> struct Held {
        typedef T* Type;
        static bool const is_class = false;
    };
}

#endif