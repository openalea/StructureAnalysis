#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_iterator_4d8f54b7d04b551ba2228a04d19d1819()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::class_< struct ::std::iterator<std::random_access_iterator_tag, std::pair<std::pair<float, float>, std::basic_string<char> >, long, const std::pair<std::pair<float, float>, std::basic_string<char> > *, const std::pair<std::pair<float, float>, std::basic_string<char> > &>, std::shared_ptr< struct ::std::iterator<std::random_access_iterator_tag, std::pair<std::pair<float, float>, std::basic_string<char> >, long, const std::pair<std::pair<float, float>, std::basic_string<char> > *, const std::pair<std::pair<float, float>, std::basic_string<char> > &> > >("_Iterator_4d8f54b7d04b551ba2228a04d19d1819", boost::python::no_init);
}