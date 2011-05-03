// Includes ====================================================================

#include "stat_tool/stat_tools.h"

#include "errors.h"

#include <boost/python.hpp>

// Using =======================================================================
using namespace boost::python;
using namespace tree_statistic;

// Declarations ================================================================

namespace  {

    // Module ======================================================================
    BOOST_PYTHON_MODULE(_errors)
    {

        // Error initialisation
        // Import the stat_tool module
        object stat_tool = import("openalea.stat_tool");
        object StatError = stat_tool.attr("StatError");
        char error_name[] = "_errors.StatTreeError";

        StatTreeError = object(handle<>(PyErr_NewException(error_name,
                                                           StatError.ptr(), NULL)));
        scope().attr("StatTreeError") = StatTreeError;

    }
}
