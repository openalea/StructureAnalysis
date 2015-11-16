#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_basic_ostream_d2368b7a8f79566f82141f26afccb03c()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::std::basic_ostream<char, std::char_traits<char> >::*method_pointer_f270a6a07f125dd0956ea3e8a64ac78b)(long) = &::std::basic_ostream<char, std::char_traits<char> >::operator<<;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::std::basic_ostream<char, std::char_traits<char> >::*method_pointer_5e3f3df97dd85f19b7075a67b863b56e)(unsigned long) = &::std::basic_ostream<char, std::char_traits<char> >::operator<<;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::std::basic_ostream<char, std::char_traits<char> >::*method_pointer_8a5fce9ed87d5ea1bbd38da46d1c5239)(bool) = &::std::basic_ostream<char, std::char_traits<char> >::operator<<;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::std::basic_ostream<char, std::char_traits<char> >::*method_pointer_d62bd81172d650f0a25fa4b1f5074dd4)(short) = &::std::basic_ostream<char, std::char_traits<char> >::operator<<;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::std::basic_ostream<char, std::char_traits<char> >::*method_pointer_23eba9fc0b8259c182b4296b601e65e2)(unsigned short) = &::std::basic_ostream<char, std::char_traits<char> >::operator<<;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::std::basic_ostream<char, std::char_traits<char> >::*method_pointer_e8ebedf7785c5b47adbd1985807571c9)(int) = &::std::basic_ostream<char, std::char_traits<char> >::operator<<;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::std::basic_ostream<char, std::char_traits<char> >::*method_pointer_af0be78a8c115608b9f0badb035215ca)(unsigned int) = &::std::basic_ostream<char, std::char_traits<char> >::operator<<;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::std::basic_ostream<char, std::char_traits<char> >::*method_pointer_47589b94c6935d8b81004f78968b4f10)(long long) = &::std::basic_ostream<char, std::char_traits<char> >::operator<<;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::std::basic_ostream<char, std::char_traits<char> >::*method_pointer_973e0b33ce185d31b1a5605ce1d90819)(unsigned long long) = &::std::basic_ostream<char, std::char_traits<char> >::operator<<;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::std::basic_ostream<char, std::char_traits<char> >::*method_pointer_11892bd195e15433bc6d9a33711830d8)(double) = &::std::basic_ostream<char, std::char_traits<char> >::operator<<;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::std::basic_ostream<char, std::char_traits<char> >::*method_pointer_795f8dffe56450809fcd9bac9022c0cd)(float) = &::std::basic_ostream<char, std::char_traits<char> >::operator<<;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::std::basic_ostream<char, std::char_traits<char> >::*method_pointer_3ff2f6b5a4ce57efa84e2b04fb3e593e)(long double) = &::std::basic_ostream<char, std::char_traits<char> >::operator<<;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::std::basic_ostream<char, std::char_traits<char> >::*method_pointer_d47bf010575d5b0a892080be15714273)(class ::std::basic_streambuf<char, std::char_traits<char> > *) = &::std::basic_ostream<char, std::char_traits<char> >::operator<<;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::std::basic_ostream<char, std::char_traits<char> >::*method_pointer_0a0356821c655592b2df94b2d3c1c194)(char) = &::std::basic_ostream<char, std::char_traits<char> >::put;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::std::basic_ostream<char, std::char_traits<char> >::*method_pointer_104ea009de895af0ba17c1b0e5120691)() = &::std::basic_ostream<char, std::char_traits<char> >::flush;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::std::basic_ostream<char, std::char_traits<char> >::*method_pointer_9fd3bb45734e566ab847be918fb76f69)(long, enum ::std::_Ios_Seekdir) = &::std::basic_ostream<char, std::char_traits<char> >::seekp;
        boost::python::class_< class ::std::basic_ostream<char, std::char_traits<char> >, std::shared_ptr< class ::std::basic_ostream<char, std::char_traits<char> > >, boost::python::bases< class ::std::basic_ios<char, std::char_traits<char> > >, boost::noncopyable >("_BasicOstream_d2368b7a8f79566f82141f26afccb03c", boost::python::no_init)
            .def(boost::python::init< class ::std::basic_streambuf<char, std::char_traits<char> > * >())
            .def("__lshift__", method_pointer_f270a6a07f125dd0956ea3e8a64ac78b, boost::python::return_internal_reference<>())
            .def("__lshift__", method_pointer_5e3f3df97dd85f19b7075a67b863b56e, boost::python::return_internal_reference<>())
            .def("__lshift__", method_pointer_8a5fce9ed87d5ea1bbd38da46d1c5239, boost::python::return_internal_reference<>())
            .def("__lshift__", method_pointer_d62bd81172d650f0a25fa4b1f5074dd4, boost::python::return_internal_reference<>())
            .def("__lshift__", method_pointer_23eba9fc0b8259c182b4296b601e65e2, boost::python::return_internal_reference<>())
            .def("__lshift__", method_pointer_e8ebedf7785c5b47adbd1985807571c9, boost::python::return_internal_reference<>())
            .def("__lshift__", method_pointer_af0be78a8c115608b9f0badb035215ca, boost::python::return_internal_reference<>())
            .def("__lshift__", method_pointer_47589b94c6935d8b81004f78968b4f10, boost::python::return_internal_reference<>())
            .def("__lshift__", method_pointer_973e0b33ce185d31b1a5605ce1d90819, boost::python::return_internal_reference<>())
            .def("__lshift__", method_pointer_11892bd195e15433bc6d9a33711830d8, boost::python::return_internal_reference<>())
            .def("__lshift__", method_pointer_795f8dffe56450809fcd9bac9022c0cd, boost::python::return_internal_reference<>())
            .def("__lshift__", method_pointer_3ff2f6b5a4ce57efa84e2b04fb3e593e, boost::python::return_internal_reference<>())
            .def("__lshift__", method_pointer_d47bf010575d5b0a892080be15714273, boost::python::return_internal_reference<>())
            .def("put", method_pointer_0a0356821c655592b2df94b2d3c1c194, boost::python::return_internal_reference<>())
            .def("flush", method_pointer_104ea009de895af0ba17c1b0e5120691, boost::python::return_internal_reference<>())
            .def("seekp", method_pointer_9fd3bb45734e566ab847be918fb76f69, boost::python::return_internal_reference<>());
        boost::python::implicitly_convertible< std::shared_ptr< class ::std::basic_ostream<char, std::char_traits<char> > >, std::shared_ptr< class ::std::basic_ios<char, std::char_traits<char> > > >();
}