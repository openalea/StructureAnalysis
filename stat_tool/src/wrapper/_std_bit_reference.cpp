#include <boost/python.hpp>
#include <stat_tool/reestimation.h>

void _std_bit_reference()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        bool (::std::_Bit_reference::*method_pointer_54df9b4d13085d5584c79d44ed761dfa)(struct ::std::_Bit_reference const &) const = &::std::_Bit_reference::operator==;
        bool (::std::_Bit_reference::*method_pointer_8eee2a34be3e5fc08b75d6dfb689721c)(struct ::std::_Bit_reference const &) const = &::std::_Bit_reference::operator<;
        void (::std::_Bit_reference::*method_pointer_16768fd13062518d94fe3e033f0686b2)() = &::std::_Bit_reference::flip;
        boost::python::class_< struct ::std::_Bit_reference, std::shared_ptr< struct ::std::_Bit_reference > >("BitReference", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< struct ::std::_Bit_reference const & >())
            .def("__eq__", method_pointer_54df9b4d13085d5584c79d44ed761dfa)
            .def("__lt__", method_pointer_8eee2a34be3e5fc08b75d6dfb689721c)
            .def("flip", method_pointer_16768fd13062518d94fe3e033f0686b2)
            .def_readwrite("m__mask", &::std::_Bit_reference::_M_mask);
}