#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_iterator_7bc2f4b4aa385af9a1b2250c3e69f597()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::class_< struct ::std::iterator<std::random_access_iterator_tag, bool, long, bool *, bool &>, std::shared_ptr< struct ::std::iterator<std::random_access_iterator_tag, bool, long, bool *, bool &> > >("_Iterator_7bc2f4b4aa385af9a1b2250c3e69f597", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< struct ::std::iterator<std::random_access_iterator_tag, bool, long, bool *, bool &> const & >());
}