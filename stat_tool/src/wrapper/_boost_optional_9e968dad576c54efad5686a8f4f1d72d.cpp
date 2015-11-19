#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _boost_optional_9e968dad576c54efad5686a8f4f1d72d()
{
        std::string boost_21ee8db290f35815a57c7bf74dca851d_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".boost");
        boost::python::object boost_21ee8db290f35815a57c7bf74dca851d_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(boost_21ee8db290f35815a57c7bf74dca851d_name.c_str()))));
        boost::python::scope().attr("boost") = boost_21ee8db290f35815a57c7bf74dca851d_module;
        boost::python::scope boost_21ee8db290f35815a57c7bf74dca851d_scope = boost_21ee8db290f35815a57c7bf74dca851d_module;
        void (::boost::optional<std::locale>::*method_pointer_ba1c7d30735154b1bad91b233fb7da7f)(class ::boost::optional<std::locale> &) = &::boost::optional<std::locale>::swap;
        class ::std::locale const & (::boost::optional<std::locale>::*method_pointer_f8862d9d880b5382bc799e1594f37b72)() const = &::boost::optional<std::locale>::get;
        class ::std::locale & (::boost::optional<std::locale>::*method_pointer_7077796a6e1c551c9c9ae72dc13b7fc0)() = &::boost::optional<std::locale>::get;
        class ::std::locale const & (::boost::optional<std::locale>::*method_pointer_d4913a65315553fd803e4251f2171c23)(class ::std::locale const &) const = &::boost::optional<std::locale>::get_value_or;
        class ::std::locale & (::boost::optional<std::locale>::*method_pointer_1e6abb6bf7425f9eb7b0e017e47ce633)(class ::std::locale &) = &::boost::optional<std::locale>::get_value_or;
        class ::std::locale const & (::boost::optional<std::locale>::*method_pointer_4668087f73d653f69bd6054f07fa8c8a)() const = &::boost::optional<std::locale>::operator*;
        class ::std::locale & (::boost::optional<std::locale>::*method_pointer_d44e8e8f04a150a293fe6e4830288da6)() = &::boost::optional<std::locale>::operator*;
        bool (::boost::optional<std::locale>::*method_pointer_ae060ed350415fa8bbaacef25706e994)() const = &::boost::optional<std::locale>::operator!;
        boost::python::class_< class ::boost::optional<std::locale>, std::shared_ptr< class ::boost::optional<std::locale> >, boost::python::bases< class ::boost::optional_detail::optional_base<std::locale> > >("_Optional_9e968dad576c54efad5686a8f4f1d72d", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< class ::std::locale const & >())
            .def(boost::python::init< bool, class ::std::locale const & >())
            .def(boost::python::init< class ::boost::optional<std::locale> const & >())
            .def("swap", method_pointer_ba1c7d30735154b1bad91b233fb7da7f)
            .def("get", method_pointer_f8862d9d880b5382bc799e1594f37b72, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("get", method_pointer_7077796a6e1c551c9c9ae72dc13b7fc0, boost::python::return_internal_reference<>())
            .def("get_value_or", method_pointer_d4913a65315553fd803e4251f2171c23, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("get_value_or", method_pointer_1e6abb6bf7425f9eb7b0e017e47ce633, boost::python::return_internal_reference<>())
            .def("__mul__", method_pointer_4668087f73d653f69bd6054f07fa8c8a, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("__mul__", method_pointer_d44e8e8f04a150a293fe6e4830288da6, boost::python::return_internal_reference<>())
            .def("__not__", method_pointer_ae060ed350415fa8bbaacef25706e994);
        boost::python::implicitly_convertible< std::shared_ptr< class ::boost::optional<std::locale> >, std::shared_ptr< class ::boost::optional_detail::optional_base<std::locale> > >();
}