#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_basic_ostream_sentry_5a580b025e1e50d782c5b88bc514ef45()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;        std::string __basic_ostream_d2368b7a8f79566f82141f26afccb03c_d2368b7a8f79566f82141f26afccb03c_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".__basic_ostream_d2368b7a8f79566f82141f26afccb03c");
        boost::python::object __basic_ostream_d2368b7a8f79566f82141f26afccb03c_d2368b7a8f79566f82141f26afccb03c_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(__basic_ostream_d2368b7a8f79566f82141f26afccb03c_d2368b7a8f79566f82141f26afccb03c_name.c_str()))));
        boost::python::scope().attr("__basic_ostream_d2368b7a8f79566f82141f26afccb03c") = __basic_ostream_d2368b7a8f79566f82141f26afccb03c_d2368b7a8f79566f82141f26afccb03c_module;
        boost::python::scope __basic_ostream_d2368b7a8f79566f82141f26afccb03c_d2368b7a8f79566f82141f26afccb03c_scope = __basic_ostream_d2368b7a8f79566f82141f26afccb03c_d2368b7a8f79566f82141f26afccb03c_module;
        boost::python::class_< class ::std::basic_ostream<char, std::char_traits<char> >::sentry, std::shared_ptr< class ::std::basic_ostream<char, std::char_traits<char> >::sentry > >("Sentry", boost::python::no_init)
            .def(boost::python::init< class ::std::basic_ostream<char, std::char_traits<char> > & >());
}