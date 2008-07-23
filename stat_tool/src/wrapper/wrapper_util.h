// WRAPPER COMMON CLASSES


#ifndef __WRAPPER_UTIL
#define __WRAPPER_UTIL

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>

#include "stat_tool/stat_tools.h"


namespace stat_tool
{
  namespace wrap_util
  {

    template<int num, int id> 
      struct UniqueInt
      {
	int v;
	enum { value=num };
	
      UniqueInt(int _v) : v(_v) { }
	operator int() const { return v; }
      };
    
    

    inline void throw_error(Format_error &error)
    {
      ostringstream error_message;
      error_message << error;
      PyErr_SetString(PyExc_Exception, (error_message.str()).c_str());
      boost::python::throw_error_already_set();
    };

    
    inline void throw_error(const char* error_message)
    {
      PyErr_SetString(PyExc_Exception, error_message);
      boost::python::throw_error_already_set();
    };

  };
};
 
#endif
