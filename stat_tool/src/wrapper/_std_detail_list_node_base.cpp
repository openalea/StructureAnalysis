#include <boost/python.hpp>
#include <stat_tool/plotable.h>

void _std_detail_list_node_base()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;        std::string detail_b537fe2d10595157ad939f6518c2cc0e_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".detail");
        boost::python::object detail_b537fe2d10595157ad939f6518c2cc0e_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(detail_b537fe2d10595157ad939f6518c2cc0e_name.c_str()))));
        boost::python::scope().attr("detail") = detail_b537fe2d10595157ad939f6518c2cc0e_module;
        boost::python::scope detail_b537fe2d10595157ad939f6518c2cc0e_scope = detail_b537fe2d10595157ad939f6518c2cc0e_module;
        void (*method_pointer_1427d93ccd8c59288d94c659c4ad3cef)(struct ::std::__detail::_List_node_base &, struct ::std::__detail::_List_node_base &) = ::std::__detail::_List_node_base::swap;
        void (::std::__detail::_List_node_base::*method_pointer_c3032c1171355925a70881b475d8886c)(struct ::std::__detail::_List_node_base * const, struct ::std::__detail::_List_node_base * const) = &::std::__detail::_List_node_base::_M_transfer;
        void (::std::__detail::_List_node_base::*method_pointer_92939f6b79b05c38b93d4ed105f03ee1)() = &::std::__detail::_List_node_base::_M_reverse;
        void (::std::__detail::_List_node_base::*method_pointer_4bfe7cada9db5a1ea9d9e363b8fc24e8)(struct ::std::__detail::_List_node_base * const) = &::std::__detail::_List_node_base::_M_hook;
        void (::std::__detail::_List_node_base::*method_pointer_4d50dce7c1e15b48af1cbc6ea84554a7)() = &::std::__detail::_List_node_base::_M_unhook;
        boost::python::class_< struct ::std::__detail::_List_node_base, std::shared_ptr< struct ::std::__detail::_List_node_base > >("ListNodeBase", boost::python::no_init)
            .def("swap", method_pointer_1427d93ccd8c59288d94c659c4ad3cef)
            .def("m__transfer", method_pointer_c3032c1171355925a70881b475d8886c)
            .def("m__reverse", method_pointer_92939f6b79b05c38b93d4ed105f03ee1)
            .def("m__hook", method_pointer_4bfe7cada9db5a1ea9d9e363b8fc24e8)
            .def("m__unhook", method_pointer_4d50dce7c1e15b48af1cbc6ea84554a7)
            .staticmethod("swap")
            .def_readwrite("m__next", &::std::__detail::_List_node_base::_M_next)
            .def_readwrite("m__prev", &::std::__detail::_List_node_base::_M_prev);
}