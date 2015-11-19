#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _boost_io_detail_format_item_f16c2d027d035688b9a271d47e6c36f7()
{
        std::string boost_21ee8db290f35815a57c7bf74dca851d_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".boost");
        boost::python::object boost_21ee8db290f35815a57c7bf74dca851d_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(boost_21ee8db290f35815a57c7bf74dca851d_name.c_str()))));
        boost::python::scope().attr("boost") = boost_21ee8db290f35815a57c7bf74dca851d_module;
        boost::python::scope boost_21ee8db290f35815a57c7bf74dca851d_scope = boost_21ee8db290f35815a57c7bf74dca851d_module;        std::string io_c231905b0f2c5e7bb2a34adfcdef6d49_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".io");
        boost::python::object io_c231905b0f2c5e7bb2a34adfcdef6d49_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(io_c231905b0f2c5e7bb2a34adfcdef6d49_name.c_str()))));
        boost::python::scope().attr("io") = io_c231905b0f2c5e7bb2a34adfcdef6d49_module;
        boost::python::scope io_c231905b0f2c5e7bb2a34adfcdef6d49_scope = io_c231905b0f2c5e7bb2a34adfcdef6d49_module;        std::string detail_899ac85328d357009cfe07d243e97c3d_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".detail");
        boost::python::object detail_899ac85328d357009cfe07d243e97c3d_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(detail_899ac85328d357009cfe07d243e97c3d_name.c_str()))));
        boost::python::scope().attr("detail") = detail_899ac85328d357009cfe07d243e97c3d_module;
        boost::python::scope detail_899ac85328d357009cfe07d243e97c3d_scope = detail_899ac85328d357009cfe07d243e97c3d_module;
        void (::boost::io::detail::format_item<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_cf7dc14d93915affa6cb67796ffcd2e9)(char) = &::boost::io::detail::format_item<char, std::char_traits<char>, std::allocator<char> >::reset;
        void (::boost::io::detail::format_item<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_820e855cf53a540eb4d83a9d6411bb23)() = &::boost::io::detail::format_item<char, std::char_traits<char>, std::allocator<char> >::compute_states;
        long (*method_pointer_142f5705975f53fdacda52a1934e2c52)() = ::boost::io::detail::format_item<char, std::char_traits<char>, std::allocator<char> >::max_streamsize;
        boost::python::class_< struct ::boost::io::detail::format_item<char, std::char_traits<char>, std::allocator<char> >, std::shared_ptr< struct ::boost::io::detail::format_item<char, std::char_traits<char>, std::allocator<char> > > >("_FormatItem_f16c2d027d035688b9a271d47e6c36f7", boost::python::no_init)
            .def(boost::python::init< char >())
            .def("reset", method_pointer_cf7dc14d93915affa6cb67796ffcd2e9)
            .def("compute_states", method_pointer_820e855cf53a540eb4d83a9d6411bb23)
            .def("max_streamsize", method_pointer_142f5705975f53fdacda52a1934e2c52)
            .staticmethod("max_streamsize")
            .def_readwrite("arg_n_", &::boost::io::detail::format_item<char, std::char_traits<char>, std::allocator<char> >::argN_)
            .def_readwrite("res_", &::boost::io::detail::format_item<char, std::char_traits<char>, std::allocator<char> >::res_)
            .def_readwrite("appendix_", &::boost::io::detail::format_item<char, std::char_traits<char>, std::allocator<char> >::appendix_)
            .def_readwrite("fmtstate_", &::boost::io::detail::format_item<char, std::char_traits<char>, std::allocator<char> >::fmtstate_)
            .def_readwrite("truncate_", &::boost::io::detail::format_item<char, std::char_traits<char>, std::allocator<char> >::truncate_)
            .def_readwrite("pad_scheme_", &::boost::io::detail::format_item<char, std::char_traits<char>, std::allocator<char> >::pad_scheme_);
}