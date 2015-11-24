#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_codecvt_abstract_base_2d0c5972479455e0a41e25a3a417b6a7()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        int (::std::__codecvt_abstract_base<char, char, __mbstate_t>::*method_pointer_34905e5735a651eb88f9b0b8a1ca59a3)() const = &::std::__codecvt_abstract_base<char, char, __mbstate_t>::encoding;
        bool (::std::__codecvt_abstract_base<char, char, __mbstate_t>::*method_pointer_888d0ed92ca85c32aad80d404bad51fa)() const = &::std::__codecvt_abstract_base<char, char, __mbstate_t>::always_noconv;
        int (::std::__codecvt_abstract_base<char, char, __mbstate_t>::*method_pointer_b1b78f9ddec95eb6a5fe04599e39bbef)() const = &::std::__codecvt_abstract_base<char, char, __mbstate_t>::max_length;
        boost::python::class_< class ::std::__codecvt_abstract_base<char, char, __mbstate_t>, std::shared_ptr< class ::std::__codecvt_abstract_base<char, char, __mbstate_t> >, boost::noncopyable >("_CodecvtAbstractBase_2d0c5972479455e0a41e25a3a417b6a7", boost::python::no_init)
            .def("encoding", method_pointer_34905e5735a651eb88f9b0b8a1ca59a3)
            .def("always_noconv", method_pointer_888d0ed92ca85c32aad80d404bad51fa)
            .def("max_length", method_pointer_b1b78f9ddec95eb6a5fe04599e39bbef);
}