#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_basic_ostringstream_316571d6f3fa564aab8b7a00819385a8()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        class ::std::basic_stringbuf<char, std::char_traits<char>, std::allocator<char> > * (::std::basic_ostringstream<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_ef83b25edf5c57379d5208fbe3099d69)() const = &::std::basic_ostringstream<char, std::char_traits<char>, std::allocator<char> >::rdbuf;
        class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > (::std::basic_ostringstream<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_3efe3244c58352f5ab681cf61abfd520)() const = &::std::basic_ostringstream<char, std::char_traits<char>, std::allocator<char> >::str;
        void (::std::basic_ostringstream<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_d9948c3716cf572fb52a495d9ce8dbcc)(class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &) = &::std::basic_ostringstream<char, std::char_traits<char>, std::allocator<char> >::str;
        boost::python::class_< class ::std::basic_ostringstream<char, std::char_traits<char>, std::allocator<char> >, std::shared_ptr< class ::std::basic_ostringstream<char, std::char_traits<char>, std::allocator<char> > >, boost::python::bases< class ::std::basic_ostream<char, std::char_traits<char> >, class ::std::basic_ios<char, std::char_traits<char> > >, boost::noncopyable >("_BasicOstringstream_316571d6f3fa564aab8b7a00819385a8", boost::python::no_init)
            .def(boost::python::init< enum ::std::_Ios_Openmode >())
            .def(boost::python::init< class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &, enum ::std::_Ios_Openmode >())
            .def("rdbuf", method_pointer_ef83b25edf5c57379d5208fbe3099d69, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("str", method_pointer_3efe3244c58352f5ab681cf61abfd520)
            .def("str", method_pointer_d9948c3716cf572fb52a495d9ce8dbcc);
        boost::python::implicitly_convertible< std::shared_ptr< class ::std::basic_ostringstream<char, std::char_traits<char>, std::allocator<char> > >, std::shared_ptr< class ::std::basic_ostream<char, std::char_traits<char> > > >();
        boost::python::implicitly_convertible< std::shared_ptr< class ::std::basic_ostringstream<char, std::char_traits<char>, std::allocator<char> > >, std::shared_ptr< class ::std::basic_ios<char, std::char_traits<char> > > >();
}