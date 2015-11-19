#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_iterator_9c0826bb654f5325bfb0de5de24ffca4()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::class_< struct ::std::iterator<std::random_access_iterator_tag, stat_tool::process_distribution, long, stat_tool::process_distribution *, stat_tool::process_distribution &>, std::shared_ptr< struct ::std::iterator<std::random_access_iterator_tag, stat_tool::process_distribution, long, stat_tool::process_distribution *, stat_tool::process_distribution &> > >("_Iterator_9c0826bb654f5325bfb0de5de24ffca4", boost::python::no_init);
}