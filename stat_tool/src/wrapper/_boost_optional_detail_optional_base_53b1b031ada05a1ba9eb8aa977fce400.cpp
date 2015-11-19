#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _boost_optional_detail_optional_base_53b1b031ada05a1ba9eb8aa977fce400()
{
        std::string boost_21ee8db290f35815a57c7bf74dca851d_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".boost");
        boost::python::object boost_21ee8db290f35815a57c7bf74dca851d_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(boost_21ee8db290f35815a57c7bf74dca851d_name.c_str()))));
        boost::python::scope().attr("boost") = boost_21ee8db290f35815a57c7bf74dca851d_module;
        boost::python::scope boost_21ee8db290f35815a57c7bf74dca851d_scope = boost_21ee8db290f35815a57c7bf74dca851d_module;        std::string optional_detail_8b9e69d15fe354e087a1a4561a039920_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".optional_detail");
        boost::python::object optional_detail_8b9e69d15fe354e087a1a4561a039920_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(optional_detail_8b9e69d15fe354e087a1a4561a039920_name.c_str()))));
        boost::python::scope().attr("optional_detail") = optional_detail_8b9e69d15fe354e087a1a4561a039920_module;
        boost::python::scope optional_detail_8b9e69d15fe354e087a1a4561a039920_scope = optional_detail_8b9e69d15fe354e087a1a4561a039920_module;
        void (::boost::optional_detail::optional_base<std::locale>::*method_pointer_f26b17ec7142528090efcd0f90542109)() = &::boost::optional_detail::optional_base<std::locale>::reset;
        void (::boost::optional_detail::optional_base<std::locale>::*method_pointer_2d8acf53ac78502a8a08cde477d0a6ff)(class ::std::locale const &) = &::boost::optional_detail::optional_base<std::locale>::reset;
        class ::std::locale const * (::boost::optional_detail::optional_base<std::locale>::*method_pointer_49af2e396a1553a4a2ed124e51a61fe3)() const = &::boost::optional_detail::optional_base<std::locale>::get_ptr;
        class ::std::locale * (::boost::optional_detail::optional_base<std::locale>::*method_pointer_1906131205235bb390b13c960502869d)() = &::boost::optional_detail::optional_base<std::locale>::get_ptr;
        bool (::boost::optional_detail::optional_base<std::locale>::*method_pointer_d0fdd913cb905149a55766798d19d61e)() const = &::boost::optional_detail::optional_base<std::locale>::is_initialized;
        boost::python::class_< class ::boost::optional_detail::optional_base<std::locale>, std::shared_ptr< class ::boost::optional_detail::optional_base<std::locale> >, boost::noncopyable >("_OptionalBase_53b1b031ada05a1ba9eb8aa977fce400", boost::python::no_init)
            .def("reset", method_pointer_f26b17ec7142528090efcd0f90542109)
            .def("reset", method_pointer_2d8acf53ac78502a8a08cde477d0a6ff)
            .def("get_ptr", method_pointer_49af2e396a1553a4a2ed124e51a61fe3, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("get_ptr", method_pointer_1906131205235bb390b13c960502869d, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("is_initialized", method_pointer_d0fdd913cb905149a55766798d19d61e);
}