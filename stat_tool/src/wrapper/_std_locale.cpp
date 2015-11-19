#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_locale()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > (::std::locale::*method_pointer_8991966518a054baa9f3c4bf41f9c27d)() const = &::std::locale::name;
        bool (::std::locale::*method_pointer_d6f72d3b43955c1bbf1e4531b69d14a9)(class ::std::locale const &) const = &::std::locale::operator==;
        bool (::std::locale::*method_pointer_4bb8912bb35551caaf4fba37f7b53426)(class ::std::locale const &) const = &::std::locale::operator!=;
        class ::std::locale (*method_pointer_5bf034e284795fd18ef049d2d37994ea)(class ::std::locale const &) = ::std::locale::global;
        class ::std::locale const & (*method_pointer_37beb60d04255aa09f819d3e2545c8ef)() = ::std::locale::classic;
        boost::python::class_< class ::std::locale, std::shared_ptr< class ::std::locale > >("Locale", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< class ::std::locale const & >())
            .def(boost::python::init< class ::std::locale const &, class ::std::locale const &, int >())
            .def("name", method_pointer_8991966518a054baa9f3c4bf41f9c27d)
            .def("__eq__", method_pointer_d6f72d3b43955c1bbf1e4531b69d14a9)
            .def("__neq__", method_pointer_4bb8912bb35551caaf4fba37f7b53426)
            .def("global", method_pointer_5bf034e284795fd18ef049d2d37994ea)
            .def("classic", method_pointer_37beb60d04255aa09f819d3e2545c8ef, boost::python::return_value_policy< boost::python::return_by_value >())
            .staticmethod("global")
            .staticmethod("classic");
}