#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_list_iterator_71ce44b5b2f65ccf9e56437976179a0b()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        struct ::std::_List_iterator<std::pair<float, float> > (::std::_List_iterator<std::pair<float, float> >::*method_pointer_a1eb6e62fb9a5fe1b4877226184abecd)() const = &::std::_List_iterator<std::pair<float, float> >::_M_const_cast;
        struct ::std::pair<float, float> & (::std::_List_iterator<std::pair<float, float> >::*method_pointer_4422a1a21892528799732125bf8af93d)() const = &::std::_List_iterator<std::pair<float, float> >::operator*;
        struct ::std::_List_iterator<std::pair<float, float> > (::std::_List_iterator<std::pair<float, float> >::*method_pointer_29c0c1cba411588981e57b9ad7c7db94)(int) = &::std::_List_iterator<std::pair<float, float> >::operator++;
        struct ::std::_List_iterator<std::pair<float, float> > (::std::_List_iterator<std::pair<float, float> >::*method_pointer_0f865a4853f1538a92063aced195db4a)(int) = &::std::_List_iterator<std::pair<float, float> >::operator--;
        bool (::std::_List_iterator<std::pair<float, float> >::*method_pointer_9eb4283d26605868828a3279a2aae131)(struct ::std::_List_iterator<std::pair<float, float> > const &) const = &::std::_List_iterator<std::pair<float, float> >::operator==;
        bool (::std::_List_iterator<std::pair<float, float> >::*method_pointer_6c4c25aad4f353168df180dce0e5b03e)(struct ::std::_List_iterator<std::pair<float, float> > const &) const = &::std::_List_iterator<std::pair<float, float> >::operator!=;
        boost::python::class_< struct ::std::_List_iterator<std::pair<float, float> >, std::shared_ptr< struct ::std::_List_iterator<std::pair<float, float> > > >("_ListIterator_71ce44b5b2f65ccf9e56437976179a0b", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< struct ::std::__detail::_List_node_base * >())
            .def("m__const_cast", method_pointer_a1eb6e62fb9a5fe1b4877226184abecd)
            .def("__mul__", method_pointer_4422a1a21892528799732125bf8af93d, boost::python::return_internal_reference<>())
            .def("__next__", method_pointer_29c0c1cba411588981e57b9ad7c7db94)
            .def("__prev__", method_pointer_0f865a4853f1538a92063aced195db4a)
            .def("__eq__", method_pointer_9eb4283d26605868828a3279a2aae131)
            .def("__neq__", method_pointer_6c4c25aad4f353168df180dce0e5b03e)
            .def_readwrite("m__node", &::std::_List_iterator<std::pair<float, float> >::_M_node);
}