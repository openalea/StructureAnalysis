#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_bit_iterator_base()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        void (::std::_Bit_iterator_base::*method_pointer_763caa4cbc2c549eb648f796106fc3c0)() = &::std::_Bit_iterator_base::_M_bump_up;
        void (::std::_Bit_iterator_base::*method_pointer_ce8b9ea492465c7cbea27bc9fbe1f650)() = &::std::_Bit_iterator_base::_M_bump_down;
        void (::std::_Bit_iterator_base::*method_pointer_c8144c67388e5522b62390ac8db8af06)(long) = &::std::_Bit_iterator_base::_M_incr;
        bool (::std::_Bit_iterator_base::*method_pointer_506d69f349cf5ee3858527e0a29a5f58)(struct ::std::_Bit_iterator_base const &) const = &::std::_Bit_iterator_base::operator==;
        bool (::std::_Bit_iterator_base::*method_pointer_46cd91b4661c5e1bbbca310cedad82d0)(struct ::std::_Bit_iterator_base const &) const = &::std::_Bit_iterator_base::operator<;
        bool (::std::_Bit_iterator_base::*method_pointer_50380632726e5c24b4d4fd7cda5c45e5)(struct ::std::_Bit_iterator_base const &) const = &::std::_Bit_iterator_base::operator!=;
        bool (::std::_Bit_iterator_base::*method_pointer_76d5227d870d5229b5f09e3c4bf9f2fc)(struct ::std::_Bit_iterator_base const &) const = &::std::_Bit_iterator_base::operator>;
        bool (::std::_Bit_iterator_base::*method_pointer_82b6e0a2012059e0a61a2f9327f23dac)(struct ::std::_Bit_iterator_base const &) const = &::std::_Bit_iterator_base::operator<=;
        bool (::std::_Bit_iterator_base::*method_pointer_ec4d726afe2e539587ac7669e7687bd7)(struct ::std::_Bit_iterator_base const &) const = &::std::_Bit_iterator_base::operator>=;
        boost::python::class_< struct ::std::_Bit_iterator_base, std::shared_ptr< struct ::std::_Bit_iterator_base >, boost::python::bases< struct ::std::iterator<std::random_access_iterator_tag, bool, long, bool *, bool &> > >("BitIteratorBase", boost::python::no_init)
            .def(boost::python::init< struct ::std::_Bit_iterator_base const & >())
            .def("m__bump_up", method_pointer_763caa4cbc2c549eb648f796106fc3c0)
            .def("m__bump_down", method_pointer_ce8b9ea492465c7cbea27bc9fbe1f650)
            .def("m__incr", method_pointer_c8144c67388e5522b62390ac8db8af06)
            .def("__eq__", method_pointer_506d69f349cf5ee3858527e0a29a5f58)
            .def("__lt__", method_pointer_46cd91b4661c5e1bbbca310cedad82d0)
            .def("__neq__", method_pointer_50380632726e5c24b4d4fd7cda5c45e5)
            .def("__gt__", method_pointer_76d5227d870d5229b5f09e3c4bf9f2fc)
            .def("__le__", method_pointer_82b6e0a2012059e0a61a2f9327f23dac)
            .def("__ge__", method_pointer_ec4d726afe2e539587ac7669e7687bd7)
            .def_readwrite("m__offset", &::std::_Bit_iterator_base::_M_offset);
        boost::python::implicitly_convertible< std::shared_ptr< struct ::std::_Bit_iterator_base >, std::shared_ptr< struct ::std::iterator<std::random_access_iterator_tag, bool, long, bool *, bool &> > >();
}