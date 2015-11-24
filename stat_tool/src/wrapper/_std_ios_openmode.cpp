#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_ios_openmode()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::enum_< enum ::std::_Ios_Openmode >("ios_openmode")
            .value("S__APP", ::std::_Ios_Openmode::_S_app)
            .value("S__ATE", ::std::_Ios_Openmode::_S_ate)
            .value("S__BIN", ::std::_Ios_Openmode::_S_bin)
            .value("S__IN", ::std::_Ios_Openmode::_S_in)
            .value("S__OUT", ::std::_Ios_Openmode::_S_out)
            .value("S__TRUNC", ::std::_Ios_Openmode::_S_trunc)
            .value("S__IOS_OPENMODE_END", ::std::_Ios_Openmode::_S_ios_openmode_end);
}