#include <boost/python.hpp>
#include <stat_tool/char_ptr_to_string.h>
#include <stat_tool/stat_tools.h>

void _stat_tool_plot_print_3ac1e584c77a5bfdb16843bbe1851fed()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        bool (*function_pointer_bc2f9e35db9254d6a6ba9577e1b838ea)(class ::stat_tool::RegressionKernel const &, class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &) = ::stat_tool::plot_print_3ac1e584c77a5bfdb16843bbe1851fed;
        boost::python::def("plot_print_3ac1e584c77a5bfdb16843bbe1851fed", function_pointer_bc2f9e35db9254d6a6ba9577e1b838ea);
}