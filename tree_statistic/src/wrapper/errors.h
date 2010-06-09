// WRAPPER COMMON CLASSES


#ifndef __WRAPPER_UTIL
#define __WRAPPER_UTIL

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>

#include "stat_tool/stat_tools.h"

using namespace boost::python;

namespace tree_statistic
{
     // New error Type
     static object StatTreeError;

     inline void throw_stat_tree_error(const char* error_message)
     {
       PyErr_SetString(StatTreeError.ptr(), error_message);
       throw_error_already_set();
     };

     // Throw exception StatTreeError with a message
     inline void throw_stat_tree_error(ostringstream& error_message)
     {
       PyErr_SetString(StatTreeError.ptr(), (error_message.str()).c_str());
       throw_error_already_set();
     };

    /*
     void throw_stat_tree_error(object StatTreeError, const char* error_message)
     {
       PyErr_SetString(StatTreeError.ptr(), error_message);
       boost::python::throw_error_already_set();
     };

     // Throw exception StatTreeError with a message
     void throw_stat_tree_error(object StatTreeError, ostringstream& error_message)
     {
       PyErr_SetString(StatTreeError.ptr(), (error_message.str()).c_str());
       boost::python::throw_error_already_set();
     };
    */
     // Throw python predefined exception with a message
     void throw_python_error(PyObject* exception_class,
                             ostringstream& error_message)
     {
       PyErr_SetString(exception_class, (error_message.str()).c_str());
       throw_error_already_set();
     };
}


#endif
