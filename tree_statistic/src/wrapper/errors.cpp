// Includes ====================================================================

#include "stat_tool/stat_tools.h"

#include "wrapper_util.h"

#include <boost/python.hpp>

// Using =======================================================================
using namespace boost::python;
using namespace tree_statistic::wrap_util;

// Declarations ================================================================

namespace  {

// Module ======================================================================
BOOST_PYTHON_MODULE(_errors)
{

    // Error initialisation
    // Import the stat_tool module
    static object StatTreeErrorClass;
    object stat_tool = import("vplants.stat_tool");
    object StatError = stat_tool.attr("StatError");
    char error_name[] = "ctrees.StatTreeError";

    StatTreeErrorClass = object(handle<>(PyErr_NewException(error_name,
                                                            StatError.ptr(), NULL)));
    scope().attr("StatTreeError") = StatTreeErrorClass;

    }
}