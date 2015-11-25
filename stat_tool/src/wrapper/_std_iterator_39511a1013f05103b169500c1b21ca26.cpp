#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/discrete_mixture.h>
#include <stat_tool/stat_tools.h>
#include <stat_tool/convolution.h>

void _std_iterator_39511a1013f05103b169500c1b21ca26()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::class_< struct ::std::iterator<std::random_access_iterator_tag, stat_tool::Reestimation<double> *, long, stat_tool::Reestimation<double> *const *, stat_tool::Reestimation<double> *const &>, std::shared_ptr< struct ::std::iterator<std::random_access_iterator_tag, stat_tool::Reestimation<double> *, long, stat_tool::Reestimation<double> *const *, stat_tool::Reestimation<double> *const &> > >("_Iterator_39511a1013f05103b169500c1b21ca26", boost::python::no_init);
}