#include "_stat_tool.h"


void wrapper_376f9a4e43f35cfb980b57231ed359d0()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    std::string name_5647113ef4105dfab0588ffcaf6c479b = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + "._ios_base");
    boost::python::object module_5647113ef4105dfab0588ffcaf6c479b(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_5647113ef4105dfab0588ffcaf6c479b.c_str()))));
    boost::python::scope().attr("_ios_base") = module_5647113ef4105dfab0588ffcaf6c479b;
    boost::python::scope scope_5647113ef4105dfab0588ffcaf6c479b = module_5647113ef4105dfab0588ffcaf6c479b;
    boost::python::enum_< enum ::std::ios_base::seekdir > enum_376f9a4e43f35cfb980b57231ed359d0("seekdir");
    enum_376f9a4e43f35cfb980b57231ed359d0.value("BEG", ::std::ios_base::beg);
    enum_376f9a4e43f35cfb980b57231ed359d0.value("CUR", ::std::ios_base::cur);
    enum_376f9a4e43f35cfb980b57231ed359d0.value("END", ::std::ios_base::end);

}