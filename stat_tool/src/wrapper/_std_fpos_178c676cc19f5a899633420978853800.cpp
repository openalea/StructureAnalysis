#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_fpos_178c676cc19f5a899633420978853800()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        void (::std::fpos<__mbstate_t>::*method_pointer_ed1a6c4e05c5553e96148e06f3ee768b)(::__mbstate_t) = &::std::fpos<__mbstate_t>::state;
        ::__mbstate_t (::std::fpos<__mbstate_t>::*method_pointer_35027a8c13575f9b9099346815865029)() const = &::std::fpos<__mbstate_t>::state;
        class ::std::fpos<__mbstate_t> & (::std::fpos<__mbstate_t>::*method_pointer_adbf6586a17c5fb4bae583a30f9d708b)(long) = &::std::fpos<__mbstate_t>::operator+=;
        class ::std::fpos<__mbstate_t> & (::std::fpos<__mbstate_t>::*method_pointer_db9583d4a8435d0eac696439cb98b7d7)(long) = &::std::fpos<__mbstate_t>::operator-=;
        class ::std::fpos<__mbstate_t> (::std::fpos<__mbstate_t>::*method_pointer_0b022b170a795f31bce964c4da680e9b)(long) const = &::std::fpos<__mbstate_t>::operator+;
        class ::std::fpos<__mbstate_t> (::std::fpos<__mbstate_t>::*method_pointer_311975ce3d2456fdbd6117a44b8271ea)(long) const = &::std::fpos<__mbstate_t>::operator-;
        long (::std::fpos<__mbstate_t>::*method_pointer_b6240d6e3b375a1eb7cb4f24fb85d5d3)(class ::std::fpos<__mbstate_t> const &) const = &::std::fpos<__mbstate_t>::operator-;
        boost::python::class_< class ::std::fpos<__mbstate_t>, std::shared_ptr< class ::std::fpos<__mbstate_t> > >("_Fpos_178c676cc19f5a899633420978853800", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< long >())
            .def("state", method_pointer_ed1a6c4e05c5553e96148e06f3ee768b)
            .def("state", method_pointer_35027a8c13575f9b9099346815865029)
            .def("__iadd__", method_pointer_adbf6586a17c5fb4bae583a30f9d708b, boost::python::return_internal_reference<>())
            .def("__isub__", method_pointer_db9583d4a8435d0eac696439cb98b7d7, boost::python::return_internal_reference<>())
            .def("__add__", method_pointer_0b022b170a795f31bce964c4da680e9b)
            .def("__sub__", method_pointer_311975ce3d2456fdbd6117a44b8271ea)
            .def("__sub__", method_pointer_b6240d6e3b375a1eb7cb4f24fb85d5d3);
}