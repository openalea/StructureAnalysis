#include "_stat_tool.h"


namespace autowig
{
    PyObject* error_167c53cdfe3c52b182c9f8fb3ce1bf67 = NULL;

    void translate_167c53cdfe3c52b182c9f8fb3ce1bf67(class ::std::system_error const & error) { PyErr_SetString(error_167c53cdfe3c52b182c9f8fb3ce1bf67, error.what()); };
}



void wrapper_167c53cdfe3c52b182c9f8fb3ce1bf67()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    std::string name_167c53cdfe3c52b182c9f8fb3ce1bf67 = boost::python::extract< std::string >(boost::python::scope().attr("__name__"));
    name_167c53cdfe3c52b182c9f8fb3ce1bf67 = name_167c53cdfe3c52b182c9f8fb3ce1bf67 + "." + "SystemError";
    autowig::error_167c53cdfe3c52b182c9f8fb3ce1bf67 = PyErr_NewException(strdup(name_167c53cdfe3c52b182c9f8fb3ce1bf67.c_str()), PyExc_RuntimeError, NULL);
    boost::python::scope().attr("SystemError") = boost::python::object(boost::python::handle<>(boost::python::borrowed(autowig::error_167c53cdfe3c52b182c9f8fb3ce1bf67)));
    boost::python::register_exception_translator< class ::std::system_error >(&autowig::translate_167c53cdfe3c52b182c9f8fb3ce1bf67);

}