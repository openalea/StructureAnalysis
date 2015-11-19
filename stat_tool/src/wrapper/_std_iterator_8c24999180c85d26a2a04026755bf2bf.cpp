#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_iterator_8c24999180c85d26a2a04026755bf2bf()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::class_< struct ::std::iterator<std::random_access_iterator_tag, stat_tool::MultiPlot, long, stat_tool::MultiPlot *, stat_tool::MultiPlot &>, std::shared_ptr< struct ::std::iterator<std::random_access_iterator_tag, stat_tool::MultiPlot, long, stat_tool::MultiPlot *, stat_tool::MultiPlot &> > >("_Iterator_8c24999180c85d26a2a04026755bf2bf", boost::python::no_init);
}