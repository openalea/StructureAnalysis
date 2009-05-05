// WRAPPER COMMON CLASSES
// copy of stat_tool/wrap_util.h

#ifndef __WRAPPER_UTIL_SEQUENCE_ANALYSIS
#define __WRAPPER_UTIL_SEQUENCE_ANALYSIS

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>

#include "stat_tool/stat_tools.h"


namespace sequence_analysis
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


   /* template<class T>
    list_2_array(boost::python::list data)
    {
      int size = len(data);
      sequence_analysis::wrap_util::auto_ptr_array<const T *>
             sequens(new const T*[size]);

      for (int i = 0; i < size; i++)
        sequens[i] = boost::python::extract<T *>(data[i]);

      return sequens;
    }
*/
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
