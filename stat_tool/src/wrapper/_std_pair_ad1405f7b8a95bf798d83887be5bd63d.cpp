#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_pair_ad1405f7b8a95bf798d83887be5bd63d()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        void (::std::pair<std::pair<float, float>, std::basic_string<char> >::*method_pointer_271d2b72a2525fe885f92796f9dfde71)(struct ::std::pair<std::pair<float, float>, std::basic_string<char> > &) = &::std::pair<std::pair<float, float>, std::basic_string<char> >::swap;
        boost::python::class_< struct ::std::pair<std::pair<float, float>, std::basic_string<char> >, std::shared_ptr< struct ::std::pair<std::pair<float, float>, std::basic_string<char> > > >("_Pair_ad1405f7b8a95bf798d83887be5bd63d", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< struct ::std::pair<float, float> const &, class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const & >())
            .def(boost::python::init< struct ::std::pair<std::pair<float, float>, std::basic_string<char> > const & >())
            .def("swap", method_pointer_271d2b72a2525fe885f92796f9dfde71)
            .def_readwrite("first", &::std::pair<std::pair<float, float>, std::basic_string<char> >::first)
            .def_readwrite("second", &::std::pair<std::pair<float, float>, std::basic_string<char> >::second);
}