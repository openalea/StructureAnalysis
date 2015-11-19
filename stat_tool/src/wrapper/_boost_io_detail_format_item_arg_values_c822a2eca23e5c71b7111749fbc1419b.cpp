#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _boost_io_detail_format_item_arg_values_c822a2eca23e5c71b7111749fbc1419b()
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
        boost::python::scope detail_899ac85328d357009cfe07d243e97c3d_scope = detail_899ac85328d357009cfe07d243e97c3d_module;        std::string __format_item_f16c2d027d035688b9a271d47e6c36f7_f16c2d027d035688b9a271d47e6c36f7_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".__format_item_f16c2d027d035688b9a271d47e6c36f7");
        boost::python::object __format_item_f16c2d027d035688b9a271d47e6c36f7_f16c2d027d035688b9a271d47e6c36f7_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(__format_item_f16c2d027d035688b9a271d47e6c36f7_f16c2d027d035688b9a271d47e6c36f7_name.c_str()))));
        boost::python::scope().attr("__format_item_f16c2d027d035688b9a271d47e6c36f7") = __format_item_f16c2d027d035688b9a271d47e6c36f7_f16c2d027d035688b9a271d47e6c36f7_module;
        boost::python::scope __format_item_f16c2d027d035688b9a271d47e6c36f7_f16c2d027d035688b9a271d47e6c36f7_scope = __format_item_f16c2d027d035688b9a271d47e6c36f7_f16c2d027d035688b9a271d47e6c36f7_module;
        boost::python::enum_< enum ::boost::io::detail::format_item<char, std::char_traits<char>, std::allocator<char> >::arg_values >("arg_values");
}