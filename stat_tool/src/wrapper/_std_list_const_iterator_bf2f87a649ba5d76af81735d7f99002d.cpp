#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_list_const_iterator_bf2f87a649ba5d76af81735d7f99002d()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        struct ::std::_List_iterator<std::pair<float, float> > (::std::_List_const_iterator<std::pair<float, float> >::*method_pointer_b265422706e450a78e55d1668584871d)() const = &::std::_List_const_iterator<std::pair<float, float> >::_M_const_cast;
        struct ::std::pair<float, float> const & (::std::_List_const_iterator<std::pair<float, float> >::*method_pointer_f11d25967321576f9fa1a1aaccfc6cf8)() const = &::std::_List_const_iterator<std::pair<float, float> >::operator*;
        struct ::std::_List_const_iterator<std::pair<float, float> > (::std::_List_const_iterator<std::pair<float, float> >::*method_pointer_6fdf00c00589547e9234241b9da275d5)(int) = &::std::_List_const_iterator<std::pair<float, float> >::operator++;
        struct ::std::_List_const_iterator<std::pair<float, float> > (::std::_List_const_iterator<std::pair<float, float> >::*method_pointer_4dc156d920d258838e3bef44b0ba9f39)(int) = &::std::_List_const_iterator<std::pair<float, float> >::operator--;
        bool (::std::_List_const_iterator<std::pair<float, float> >::*method_pointer_e8e3013fc0565b989e698f792812ee35)(struct ::std::_List_const_iterator<std::pair<float, float> > const &) const = &::std::_List_const_iterator<std::pair<float, float> >::operator==;
        bool (::std::_List_const_iterator<std::pair<float, float> >::*method_pointer_32e1a59cd6a05f31a02a748305ab1549)(struct ::std::_List_const_iterator<std::pair<float, float> > const &) const = &::std::_List_const_iterator<std::pair<float, float> >::operator!=;
        boost::python::class_< struct ::std::_List_const_iterator<std::pair<float, float> >, std::shared_ptr< struct ::std::_List_const_iterator<std::pair<float, float> > > >("_ListConstIterator_bf2f87a649ba5d76af81735d7f99002d", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< struct ::std::__detail::_List_node_base const * >())
            .def(boost::python::init< struct ::std::_List_iterator<std::pair<float, float> > const & >())
            .def(boost::python::init< struct ::std::_List_const_iterator<std::pair<float, float> > const & >())
            .def("m__const_cast", method_pointer_b265422706e450a78e55d1668584871d)
            .def("__mul__", method_pointer_f11d25967321576f9fa1a1aaccfc6cf8, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("__next__", method_pointer_6fdf00c00589547e9234241b9da275d5)
            .def("__prev__", method_pointer_4dc156d920d258838e3bef44b0ba9f39)
            .def("__eq__", method_pointer_e8e3013fc0565b989e698f792812ee35)
            .def("__neq__", method_pointer_32e1a59cd6a05f31a02a748305ab1549)
            .def_readwrite("m__node", &::std::_List_const_iterator<std::pair<float, float> >::_M_node);
}