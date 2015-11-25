#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_iterator_00810bf83c9a57e992c8045efcf3e5f0()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::class_< struct ::std::iterator<std::random_access_iterator_tag, stat_tool::discrete_parametric, long, const stat_tool::discrete_parametric *, const stat_tool::discrete_parametric &>, std::shared_ptr< struct ::std::iterator<std::random_access_iterator_tag, stat_tool::discrete_parametric, long, const stat_tool::discrete_parametric *, const stat_tool::discrete_parametric &> > >("_Iterator_00810bf83c9a57e992c8045efcf3e5f0", boost::python::no_init);
}