#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_bit_iterator()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        struct ::std::_Bit_iterator (::std::_Bit_iterator::*method_pointer_bea69dc8dde1504c81faf15494ad6175)() const = &::std::_Bit_iterator::_M_const_cast;
        struct ::std::_Bit_reference (::std::_Bit_iterator::*method_pointer_4f9014171b115c4e823586be6be0d51d)() const = &::std::_Bit_iterator::operator*;
        struct ::std::_Bit_iterator (::std::_Bit_iterator::*method_pointer_4e6ed95df1a45765a5931e4a0b487a35)(int) = &::std::_Bit_iterator::operator++;
        struct ::std::_Bit_iterator (::std::_Bit_iterator::*method_pointer_e614999933d051c4a8914467445ca301)(int) = &::std::_Bit_iterator::operator--;
        struct ::std::_Bit_iterator & (::std::_Bit_iterator::*method_pointer_1723b78182035cd7977cdd7c9d673c72)(long) = &::std::_Bit_iterator::operator+=;
        struct ::std::_Bit_iterator & (::std::_Bit_iterator::*method_pointer_ffbbcf2a236a50e7b027e093539ddc36)(long) = &::std::_Bit_iterator::operator-=;
        struct ::std::_Bit_iterator (::std::_Bit_iterator::*method_pointer_13f989495d3b543d85086dad2222a8d5)(long) const = &::std::_Bit_iterator::operator+;
        struct ::std::_Bit_iterator (::std::_Bit_iterator::*method_pointer_8739c349f58d5d048364809bb1a15361)(long) const = &::std::_Bit_iterator::operator-;
        struct ::std::_Bit_reference (::std::_Bit_iterator::*method_pointer_7db0bbb1094b55849f7dfd9d25c0d5f2)(long) const = &::std::_Bit_iterator::operator[];
        boost::python::class_< struct ::std::_Bit_iterator, std::shared_ptr< struct ::std::_Bit_iterator >, boost::python::bases< struct ::std::_Bit_iterator_base > >("BitIterator", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< struct ::std::_Bit_iterator const & >())
            .def("m__const_cast", method_pointer_bea69dc8dde1504c81faf15494ad6175)
            .def("__mul__", method_pointer_4f9014171b115c4e823586be6be0d51d)
            .def("__next__", method_pointer_4e6ed95df1a45765a5931e4a0b487a35)
            .def("__prev__", method_pointer_e614999933d051c4a8914467445ca301)
            .def("__iadd__", method_pointer_1723b78182035cd7977cdd7c9d673c72, boost::python::return_internal_reference<>())
            .def("__isub__", method_pointer_ffbbcf2a236a50e7b027e093539ddc36, boost::python::return_internal_reference<>())
            .def("__add__", method_pointer_13f989495d3b543d85086dad2222a8d5)
            .def("__sub__", method_pointer_8739c349f58d5d048364809bb1a15361)
            .def("__getitem__", method_pointer_7db0bbb1094b55849f7dfd9d25c0d5f2);
        boost::python::implicitly_convertible< std::shared_ptr< struct ::std::_Bit_iterator >, std::shared_ptr< struct ::std::_Bit_iterator_base > >();
}