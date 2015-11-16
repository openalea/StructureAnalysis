#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_basic_ios_dc0bed07c88f59b7b3e4f11a1fde9779()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        bool (::std::basic_ios<char, std::char_traits<char> >::*method_pointer_6bc1f908c7dc571eb767e24056cf6e24)() const = &::std::basic_ios<char, std::char_traits<char> >::operator!;
        enum ::std::_Ios_Iostate (::std::basic_ios<char, std::char_traits<char> >::*method_pointer_782e5997f99f535d92a745d45b71d60e)() const = &::std::basic_ios<char, std::char_traits<char> >::rdstate;
        void (::std::basic_ios<char, std::char_traits<char> >::*method_pointer_8e8d726cb074552d9af16b48432d41bf)(enum ::std::_Ios_Iostate) = &::std::basic_ios<char, std::char_traits<char> >::clear;
        void (::std::basic_ios<char, std::char_traits<char> >::*method_pointer_38e2cdbc9116522595a702004fcfdfe2)(enum ::std::_Ios_Iostate) = &::std::basic_ios<char, std::char_traits<char> >::setstate;
        void (::std::basic_ios<char, std::char_traits<char> >::*method_pointer_f9e2edab94cc585aa855888bf7e328bf)(enum ::std::_Ios_Iostate) = &::std::basic_ios<char, std::char_traits<char> >::_M_setstate;
        bool (::std::basic_ios<char, std::char_traits<char> >::*method_pointer_46730aa08c215c3c99d5ee8cec6a7ce9)() const = &::std::basic_ios<char, std::char_traits<char> >::good;
        bool (::std::basic_ios<char, std::char_traits<char> >::*method_pointer_0629626d20335aeda86c9bf2e222bcef)() const = &::std::basic_ios<char, std::char_traits<char> >::eof;
        bool (::std::basic_ios<char, std::char_traits<char> >::*method_pointer_8218cb5d736b5e7a949b37b05febc396)() const = &::std::basic_ios<char, std::char_traits<char> >::fail;
        bool (::std::basic_ios<char, std::char_traits<char> >::*method_pointer_d0c771645fe25f909f28aef2ac138e1b)() const = &::std::basic_ios<char, std::char_traits<char> >::bad;
        enum ::std::_Ios_Iostate (::std::basic_ios<char, std::char_traits<char> >::*method_pointer_e96abfee93565ef6995bafdf7e6eee08)() const = &::std::basic_ios<char, std::char_traits<char> >::exceptions;
        void (::std::basic_ios<char, std::char_traits<char> >::*method_pointer_c62d5646f5115297b1e4c0b6005bc69d)(enum ::std::_Ios_Iostate) = &::std::basic_ios<char, std::char_traits<char> >::exceptions;
        class ::std::basic_ostream<char, std::char_traits<char> > * (::std::basic_ios<char, std::char_traits<char> >::*method_pointer_ad29ae6111e6577b8f29eec8d14dfb37)() const = &::std::basic_ios<char, std::char_traits<char> >::tie;
        class ::std::basic_ostream<char, std::char_traits<char> > * (::std::basic_ios<char, std::char_traits<char> >::*method_pointer_d3f3e45fd5c25f8384ef8013ce0eb55b)(class ::std::basic_ostream<char, std::char_traits<char> > *) = &::std::basic_ios<char, std::char_traits<char> >::tie;
        class ::std::basic_streambuf<char, std::char_traits<char> > * (::std::basic_ios<char, std::char_traits<char> >::*method_pointer_8361471080435ad3b93ef0771c3ccd6c)() const = &::std::basic_ios<char, std::char_traits<char> >::rdbuf;
        class ::std::basic_streambuf<char, std::char_traits<char> > * (::std::basic_ios<char, std::char_traits<char> >::*method_pointer_106e378f6a6857478d7c21de6e9ff6c8)(class ::std::basic_streambuf<char, std::char_traits<char> > *) = &::std::basic_ios<char, std::char_traits<char> >::rdbuf;
        class ::std::basic_ios<char, std::char_traits<char> > & (::std::basic_ios<char, std::char_traits<char> >::*method_pointer_e47206a5f39a578289cdfa18dbbc9cfc)(class ::std::basic_ios<char, std::char_traits<char> > const &) = &::std::basic_ios<char, std::char_traits<char> >::copyfmt;
        char (::std::basic_ios<char, std::char_traits<char> >::*method_pointer_826cf7aa5a9353e9ae959fa000ab9eac)() const = &::std::basic_ios<char, std::char_traits<char> >::fill;
        char (::std::basic_ios<char, std::char_traits<char> >::*method_pointer_b10f7f3ff7975c99b912af0fc02fbcc9)(char) = &::std::basic_ios<char, std::char_traits<char> >::fill;
        class ::std::locale (::std::basic_ios<char, std::char_traits<char> >::*method_pointer_43bc9d1d3a1c54618a9ea674e0480598)(class ::std::locale const &) = &::std::basic_ios<char, std::char_traits<char> >::imbue;
        char (::std::basic_ios<char, std::char_traits<char> >::*method_pointer_3d0fcfa0f93b544b843fc23cda6849d8)(char, char) const = &::std::basic_ios<char, std::char_traits<char> >::narrow;
        char (::std::basic_ios<char, std::char_traits<char> >::*method_pointer_89da866d5a03589c823dfc23399dbc33)(char) const = &::std::basic_ios<char, std::char_traits<char> >::widen;
        boost::python::class_< class ::std::basic_ios<char, std::char_traits<char> >, std::shared_ptr< class ::std::basic_ios<char, std::char_traits<char> > >, boost::python::bases< class ::std::ios_base >, boost::noncopyable >("_BasicIos_dc0bed07c88f59b7b3e4f11a1fde9779", boost::python::no_init)
            .def(boost::python::init< class ::std::basic_streambuf<char, std::char_traits<char> > * >())
            .def("__not__", method_pointer_6bc1f908c7dc571eb767e24056cf6e24)
            .def("rdstate", method_pointer_782e5997f99f535d92a745d45b71d60e)
            .def("clear", method_pointer_8e8d726cb074552d9af16b48432d41bf)
            .def("setstate", method_pointer_38e2cdbc9116522595a702004fcfdfe2)
            .def("m__setstate", method_pointer_f9e2edab94cc585aa855888bf7e328bf)
            .def("good", method_pointer_46730aa08c215c3c99d5ee8cec6a7ce9)
            .def("eof", method_pointer_0629626d20335aeda86c9bf2e222bcef)
            .def("fail", method_pointer_8218cb5d736b5e7a949b37b05febc396)
            .def("bad", method_pointer_d0c771645fe25f909f28aef2ac138e1b)
            .def("exceptions", method_pointer_e96abfee93565ef6995bafdf7e6eee08)
            .def("exceptions", method_pointer_c62d5646f5115297b1e4c0b6005bc69d)
            .def("tie", method_pointer_ad29ae6111e6577b8f29eec8d14dfb37, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("tie", method_pointer_d3f3e45fd5c25f8384ef8013ce0eb55b, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("rdbuf", method_pointer_8361471080435ad3b93ef0771c3ccd6c, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("rdbuf", method_pointer_106e378f6a6857478d7c21de6e9ff6c8, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("copyfmt", method_pointer_e47206a5f39a578289cdfa18dbbc9cfc, boost::python::return_internal_reference<>())
            .def("fill", method_pointer_826cf7aa5a9353e9ae959fa000ab9eac)
            .def("fill", method_pointer_b10f7f3ff7975c99b912af0fc02fbcc9)
            .def("imbue", method_pointer_43bc9d1d3a1c54618a9ea674e0480598)
            .def("narrow", method_pointer_3d0fcfa0f93b544b843fc23cda6849d8)
            .def("widen", method_pointer_89da866d5a03589c823dfc23399dbc33);
        boost::python::implicitly_convertible< std::shared_ptr< class ::std::basic_ios<char, std::char_traits<char> > >, std::shared_ptr< class ::std::ios_base > >();
}