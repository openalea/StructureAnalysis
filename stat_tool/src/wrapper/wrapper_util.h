// WRAPPER COMMON CLASSES


#ifndef __WRAPPER_UTIL
#define __WRAPPER_UTIL

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>

#include "stat_tool/stat_tools.h"

using namespace boost::python;


namespace stat_tool
{
  namespace wrap_util
  {
    // This new error is instantiated in the file stat_tool_wrap
    // during the initialisation of the module.
    static object StatErrorClass;


    template<int num, int id>
      struct UniqueInt
      {
        int v;
        enum { value=num };

        UniqueInt(int _v) : v(_v) { }
        operator int() const { return v; }
      };



    inline void throw_stat_error(StatError &error)
    {
      ostringstream error_message;
      error_message << error;
      PyErr_SetString(StatErrorClass.ptr(), (error_message.str()).c_str());
      boost::python::throw_error_already_set();
    };

    inline void throw_stat_error(const char* error_message)
    {
      PyErr_SetString(StatErrorClass.ptr(), error_message);
      boost::python::throw_error_already_set();
    };


    inline void throw_error(StatError &error)
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



    // Redefine auto_ptr for array (use delete[] instead)
    template<class T>
      class auto_ptr_array
      {
      private:
        T* ap;    // refers to the actual owned object (if any)
      public:
        typedef T element_type;

        // constructor
        explicit auto_ptr_array (T* ptr = 0) throw() : ap(ptr) { }

        // copy constructors (with implicit conversion)
        // - note: nonconstant parameter
        auto_ptr_array (auto_ptr_array& rhs) throw() : ap(rhs.release()) { }

        template<class Y>
        auto_ptr_array (auto_ptr_array<Y>& rhs) throw() : ap(rhs.release()) { }

        // assignments (with implicit conversion)
        // - note: nonconstant parameter
        auto_ptr_array& operator= (auto_ptr_array& rhs) throw()
        {
            reset(rhs.release());
            return *this;
        }
        template<class Y>
        auto_ptr_array& operator= (auto_ptr_array<Y>& rhs) throw()
        {
            reset(rhs.release());
            return *this;
        }

        // destructor
        ~auto_ptr_array() throw()
        {
            delete[] ap;
        }

        // value access
        T* get() const throw()
        {
            return ap;
        }

    T& operator[](int index) throw()
    {
      return ap[index];
    }

        T& operator*() const throw()
        {
            return *ap;
        }
        T* operator->() const throw()
        {
            return ap;
        }

        // release ownership
        T* release() throw()
        {
            T* tmp(ap);
            ap = 0;
            return tmp;
        }

        // reset value
        void reset (T* ptr=0) throw()
        {
            if (ap != ptr)
            {
                delete ap;
                ap = ptr;
            }
        }
      };



  };
};

#endif
