#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_locale_impl()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;        std::string _locale_2f3439617e035c41b1282a03e900ef19_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + "._locale");
        boost::python::object _locale_2f3439617e035c41b1282a03e900ef19_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(_locale_2f3439617e035c41b1282a03e900ef19_name.c_str()))));
        boost::python::scope().attr("_locale") = _locale_2f3439617e035c41b1282a03e900ef19_module;
        boost::python::scope _locale_2f3439617e035c41b1282a03e900ef19_scope = _locale_2f3439617e035c41b1282a03e900ef19_module;
        boost::python::class_< class ::std::locale::_Impl, std::shared_ptr< class ::std::locale::_Impl >, boost::noncopyable >("Impl", boost::python::no_init);
}