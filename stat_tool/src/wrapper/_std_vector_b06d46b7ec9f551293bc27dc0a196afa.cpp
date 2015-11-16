#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_vector_b06d46b7ec9f551293bc27dc0a196afa()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        void (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_208fbea201e25b8884bf6c5e64906edf)(unsigned long, enum ::stat_tool::process_distribution const &) = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::assign;
        unsigned long (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_2976530ccb5657b6b6c0a54e0f5c9bf8)() const = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::size;
        unsigned long (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_1aacab760834510ea8d4153fc713ba92)() const = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::max_size;
        void (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_70ee84578153504299edef13598e052d)(unsigned long) = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::resize;
        void (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_bd273cd8cf6052bab101f81c2e8e2af8)(unsigned long, enum ::stat_tool::process_distribution const &) = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::resize;
        void (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_7e1a87abab0e55769d56b92f87de9364)() = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::shrink_to_fit;
        unsigned long (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_c4a00f13846e5f839d297a023ef32936)() const = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::capacity;
        bool (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_ca8611fc69aa580da7d1319be7a1b2af)() const = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::empty;
        void (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_a51bd8eeb45656289869d33118b75250)(unsigned long) = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::reserve;
        enum ::stat_tool::process_distribution const & (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_fdeddab37b135a37b33c6270959703f7)(unsigned long) const = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::operator[];
        enum ::stat_tool::process_distribution const & (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_563d298b97c95b388324b19f0f74ab68)(unsigned long) const = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::at;
        enum ::stat_tool::process_distribution const & (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_73ae975684b25694bbeee2d13496bdba)() const = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::front;
        enum ::stat_tool::process_distribution const & (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_6acc6f10ad2e54568819075f729606fa)() const = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::back;
        void (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_dd51d9a7dfbd547ea2f1392263c449cb)(enum ::stat_tool::process_distribution const &) = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::push_back;
        void (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_2d363ed2fc5950f68a03d934d57851cc)() = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::pop_back;
        void (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_f617e6b977e651a6ae8985dd17f99bda)(class ::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> > &) = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::swap;
        void (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_823a942d3fe95fb9bee873042664f0b9)() = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::clear;
        boost::python::class_< class ::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >, std::shared_ptr< class ::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> > > >("_Vector_b06d46b7ec9f551293bc27dc0a196afa", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< class ::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> > const & >())
            .def("assign", method_pointer_208fbea201e25b8884bf6c5e64906edf)
            .def("size", method_pointer_2976530ccb5657b6b6c0a54e0f5c9bf8)
            .def("max_size", method_pointer_1aacab760834510ea8d4153fc713ba92)
            .def("resize", method_pointer_70ee84578153504299edef13598e052d)
            .def("resize", method_pointer_bd273cd8cf6052bab101f81c2e8e2af8)
            .def("shrink_to_fit", method_pointer_7e1a87abab0e55769d56b92f87de9364)
            .def("capacity", method_pointer_c4a00f13846e5f839d297a023ef32936)
            .def("empty", method_pointer_ca8611fc69aa580da7d1319be7a1b2af)
            .def("reserve", method_pointer_a51bd8eeb45656289869d33118b75250)
            .def("__getitem__", method_pointer_fdeddab37b135a37b33c6270959703f7, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("at", method_pointer_563d298b97c95b388324b19f0f74ab68, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("front", method_pointer_73ae975684b25694bbeee2d13496bdba, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("back", method_pointer_6acc6f10ad2e54568819075f729606fa, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("push_back", method_pointer_dd51d9a7dfbd547ea2f1392263c449cb)
            .def("pop_back", method_pointer_2d363ed2fc5950f68a03d934d57851cc)
            .def("swap", method_pointer_f617e6b977e651a6ae8985dd17f99bda)
            .def("clear", method_pointer_823a942d3fe95fb9bee873042664f0b9);
}