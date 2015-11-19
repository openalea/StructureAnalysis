#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_ostreambuf_iterator_8ed59c366b1d5196ad34b13b32ead60f()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        class ::std::ostreambuf_iterator<char, std::char_traits<char> > & (::std::ostreambuf_iterator<char, std::char_traits<char> >::*method_pointer_89787c9a48f05f71868e8147122a3834)() = &::std::ostreambuf_iterator<char, std::char_traits<char> >::operator*;
        class ::std::ostreambuf_iterator<char, std::char_traits<char> > & (::std::ostreambuf_iterator<char, std::char_traits<char> >::*method_pointer_b4832a06446959ed8962b8f4795f2bed)(int) = &::std::ostreambuf_iterator<char, std::char_traits<char> >::operator++;
        bool (::std::ostreambuf_iterator<char, std::char_traits<char> >::*method_pointer_5e66be1634685ecbb19a1289c5b77655)() const = &::std::ostreambuf_iterator<char, std::char_traits<char> >::failed;
        boost::python::class_< class ::std::ostreambuf_iterator<char, std::char_traits<char> >, std::shared_ptr< class ::std::ostreambuf_iterator<char, std::char_traits<char> > > >("_OstreambufIterator_8ed59c366b1d5196ad34b13b32ead60f", boost::python::no_init)
            .def(boost::python::init< class ::std::basic_ostream<char, std::char_traits<char> > & >())
            .def(boost::python::init< class ::std::basic_streambuf<char, std::char_traits<char> > * >())
            .def("__mul__", method_pointer_89787c9a48f05f71868e8147122a3834, boost::python::return_internal_reference<>())
            .def("__next__", method_pointer_b4832a06446959ed8962b8f4795f2bed, boost::python::return_internal_reference<>())
            .def("failed", method_pointer_5e66be1634685ecbb19a1289c5b77655);
}