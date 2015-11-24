#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_bit_const_iterator()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        struct ::std::_Bit_iterator (::std::_Bit_const_iterator::*method_pointer_cbbdeb3cd5d05855beab78cc498c3b0a)() const = &::std::_Bit_const_iterator::_M_const_cast;
        bool (::std::_Bit_const_iterator::*method_pointer_d2b3981a24a45330a021fa73eb0e1f19)() const = &::std::_Bit_const_iterator::operator*;
        struct ::std::_Bit_const_iterator (::std::_Bit_const_iterator::*method_pointer_350a82e399c153f0a0fea357deb735ef)(int) = &::std::_Bit_const_iterator::operator++;
        struct ::std::_Bit_const_iterator (::std::_Bit_const_iterator::*method_pointer_8c5c23a6a4465aa89e2b0525f83acb8e)(int) = &::std::_Bit_const_iterator::operator--;
        struct ::std::_Bit_const_iterator & (::std::_Bit_const_iterator::*method_pointer_cde3abc4e00951e99230923b09b31af1)(long) = &::std::_Bit_const_iterator::operator+=;
        struct ::std::_Bit_const_iterator & (::std::_Bit_const_iterator::*method_pointer_9d0ea07b94c35f23abc6edbac2db20bc)(long) = &::std::_Bit_const_iterator::operator-=;
        struct ::std::_Bit_const_iterator (::std::_Bit_const_iterator::*method_pointer_20a4881e804d503faac1b4ff89f72f13)(long) const = &::std::_Bit_const_iterator::operator+;
        struct ::std::_Bit_const_iterator (::std::_Bit_const_iterator::*method_pointer_7dd921d7fd4d578da2a5c6d6ac866399)(long) const = &::std::_Bit_const_iterator::operator-;
        bool (::std::_Bit_const_iterator::*method_pointer_60a4016bf6f95387bd596d9eb2870d10)(long) const = &::std::_Bit_const_iterator::operator[];
        boost::python::class_< struct ::std::_Bit_const_iterator, std::shared_ptr< struct ::std::_Bit_const_iterator >, boost::python::bases< struct ::std::_Bit_iterator_base > >("BitConstIterator", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< struct ::std::_Bit_iterator const & >())
            .def(boost::python::init< struct ::std::_Bit_const_iterator const & >())
            .def("m__const_cast", method_pointer_cbbdeb3cd5d05855beab78cc498c3b0a)
            .def("__mul__", method_pointer_d2b3981a24a45330a021fa73eb0e1f19)
            .def("__next__", method_pointer_350a82e399c153f0a0fea357deb735ef)
            .def("__prev__", method_pointer_8c5c23a6a4465aa89e2b0525f83acb8e)
            .def("__iadd__", method_pointer_cde3abc4e00951e99230923b09b31af1, boost::python::return_internal_reference<>())
            .def("__isub__", method_pointer_9d0ea07b94c35f23abc6edbac2db20bc, boost::python::return_internal_reference<>())
            .def("__add__", method_pointer_20a4881e804d503faac1b4ff89f72f13)
            .def("__sub__", method_pointer_7dd921d7fd4d578da2a5c6d6ac866399)
            .def("__getitem__", method_pointer_60a4016bf6f95387bd596d9eb2870d10);
        boost::python::implicitly_convertible< std::shared_ptr< struct ::std::_Bit_const_iterator >, std::shared_ptr< struct ::std::_Bit_iterator_base > >();
}