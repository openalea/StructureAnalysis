#include "_stat_tool.h"


namespace autowig
{
    PyObject* error_b6a067f70ca259909a6411bfb14cfdca = NULL;

    void translate_b6a067f70ca259909a6411bfb14cfdca(class ::std::ios_base::failure const & error) { PyErr_SetString(error_b6a067f70ca259909a6411bfb14cfdca, error.what()); };
}



void wrapper_b6a067f70ca259909a6411bfb14cfdca()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    std::string name_5647113ef4105dfab0588ffcaf6c479b = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + "._ios_base");
    boost::python::object module_5647113ef4105dfab0588ffcaf6c479b(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_5647113ef4105dfab0588ffcaf6c479b.c_str()))));
    boost::python::scope().attr("_ios_base") = module_5647113ef4105dfab0588ffcaf6c479b;
    boost::python::scope scope_5647113ef4105dfab0588ffcaf6c479b = module_5647113ef4105dfab0588ffcaf6c479b;
    std::string name_b6a067f70ca259909a6411bfb14cfdca = boost::python::extract< std::string >(boost::python::scope().attr("__name__"));
    name_b6a067f70ca259909a6411bfb14cfdca = name_b6a067f70ca259909a6411bfb14cfdca + "." + "Failure";
    autowig::error_b6a067f70ca259909a6411bfb14cfdca = PyErr_NewException(strdup(name_b6a067f70ca259909a6411bfb14cfdca.c_str()), PyExc_RuntimeError, NULL);
    boost::python::scope().attr("Failure") = boost::python::object(boost::python::handle<>(boost::python::borrowed(autowig::error_b6a067f70ca259909a6411bfb14cfdca)));
    boost::python::register_exception_translator< class ::std::ios_base::failure >(&autowig::translate_b6a067f70ca259909a6411bfb14cfdca);

}