#include "_stat_tool.h"


namespace autowig
{
    PyObject* error_a636d351527b51298739e1ccb4313967 = NULL;

    void translate_a636d351527b51298739e1ccb4313967(class ::std::runtime_error const & error) { PyErr_SetString(error_a636d351527b51298739e1ccb4313967, error.what()); };
}



void wrapper_a636d351527b51298739e1ccb4313967()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    std::string name_a636d351527b51298739e1ccb4313967 = boost::python::extract< std::string >(boost::python::scope().attr("__name__"));
    name_a636d351527b51298739e1ccb4313967 = name_a636d351527b51298739e1ccb4313967 + "." + "RuntimeError";
    autowig::error_a636d351527b51298739e1ccb4313967 = PyErr_NewException(strdup(name_a636d351527b51298739e1ccb4313967.c_str()), PyExc_RuntimeError, NULL);
    boost::python::scope().attr("RuntimeError") = boost::python::object(boost::python::handle<>(boost::python::borrowed(autowig::error_a636d351527b51298739e1ccb4313967)));
    boost::python::register_exception_translator< class ::std::runtime_error >(&autowig::translate_a636d351527b51298739e1ccb4313967);

}