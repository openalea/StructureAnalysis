#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_iterator_3c1edee0a8ac5394aef746b3bee6ccc1()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::class_< struct ::std::iterator<std::random_access_iterator_tag, stat_tool::SinglePlot, long, stat_tool::SinglePlot *, stat_tool::SinglePlot &>, std::shared_ptr< struct ::std::iterator<std::random_access_iterator_tag, stat_tool::SinglePlot, long, stat_tool::SinglePlot *, stat_tool::SinglePlot &> > >("_Iterator_3c1edee0a8ac5394aef746b3bee6ccc1", boost::python::no_init);
}