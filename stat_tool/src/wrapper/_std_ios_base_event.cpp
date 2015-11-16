#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_ios_base_event()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;        std::string _ios_base_5647113ef4105dfab0588ffcaf6c479b_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + "._ios_base");
        boost::python::object _ios_base_5647113ef4105dfab0588ffcaf6c479b_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(_ios_base_5647113ef4105dfab0588ffcaf6c479b_name.c_str()))));
        boost::python::scope().attr("_ios_base") = _ios_base_5647113ef4105dfab0588ffcaf6c479b_module;
        boost::python::scope _ios_base_5647113ef4105dfab0588ffcaf6c479b_scope = _ios_base_5647113ef4105dfab0588ffcaf6c479b_module;
        boost::python::enum_< enum ::std::ios_base::event >("event");
}