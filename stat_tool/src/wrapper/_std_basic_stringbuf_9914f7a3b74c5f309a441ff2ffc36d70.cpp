#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_basic_stringbuf_9914f7a3b74c5f309a441ff2ffc36d70()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > (::std::basic_stringbuf<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_8e4a0081b508568c908a5ef90032657c)() const = &::std::basic_stringbuf<char, std::char_traits<char>, std::allocator<char> >::str;
        void (::std::basic_stringbuf<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_ef86af5d59cb5c85819ed5d64bd2e75d)(class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &) = &::std::basic_stringbuf<char, std::char_traits<char>, std::allocator<char> >::str;
        boost::python::class_< class ::std::basic_stringbuf<char, std::char_traits<char>, std::allocator<char> >, std::shared_ptr< class ::std::basic_stringbuf<char, std::char_traits<char>, std::allocator<char> > >, boost::python::bases< class ::std::basic_streambuf<char, std::char_traits<char> > >, boost::noncopyable >("_BasicStringbuf_9914f7a3b74c5f309a441ff2ffc36d70", boost::python::no_init)
            .def(boost::python::init< enum ::std::_Ios_Openmode >())
            .def(boost::python::init< class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &, enum ::std::_Ios_Openmode >())
            .def("str", method_pointer_8e4a0081b508568c908a5ef90032657c)
            .def("str", method_pointer_ef86af5d59cb5c85819ed5d64bd2e75d);
        boost::python::implicitly_convertible< std::shared_ptr< class ::std::basic_stringbuf<char, std::char_traits<char>, std::allocator<char> > >, std::shared_ptr< class ::std::basic_streambuf<char, std::char_traits<char> > > >();
}