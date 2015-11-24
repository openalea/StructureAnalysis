#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_basic_istream_sentry_c702204f91205b45a91cfb9aa50bbeac()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;        std::string __basic_istream_4c9e5eab62d358b9b9e7849871764a79_4c9e5eab62d358b9b9e7849871764a79_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".__basic_istream_4c9e5eab62d358b9b9e7849871764a79");
        boost::python::object __basic_istream_4c9e5eab62d358b9b9e7849871764a79_4c9e5eab62d358b9b9e7849871764a79_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(__basic_istream_4c9e5eab62d358b9b9e7849871764a79_4c9e5eab62d358b9b9e7849871764a79_name.c_str()))));
        boost::python::scope().attr("__basic_istream_4c9e5eab62d358b9b9e7849871764a79") = __basic_istream_4c9e5eab62d358b9b9e7849871764a79_4c9e5eab62d358b9b9e7849871764a79_module;
        boost::python::scope __basic_istream_4c9e5eab62d358b9b9e7849871764a79_4c9e5eab62d358b9b9e7849871764a79_scope = __basic_istream_4c9e5eab62d358b9b9e7849871764a79_4c9e5eab62d358b9b9e7849871764a79_module;
        boost::python::class_< class ::std::basic_istream<char, std::char_traits<char> >::sentry, std::shared_ptr< class ::std::basic_istream<char, std::char_traits<char> >::sentry > >("Sentry", boost::python::no_init)
            .def(boost::python::init< class ::std::basic_istream<char, std::char_traits<char> > &, bool >());
}