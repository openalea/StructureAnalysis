#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/char_ptr_to_string.h>
#include <stat_tool/stat_tools.h>

void _stat_tool_contingency_table_53063479d90a5bb5bad885693e2e5f52()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        bool (*function_pointer_ebdcbba919fd5e7a9af8bccab4b78e8e)(class ::stat_tool::Vectors const &, class ::stat_tool::StatError &, class ::std::basic_ostream<char, std::char_traits<char> > &, int, int, class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &, enum ::stat_tool::output_format) = ::stat_tool::contingency_table_53063479d90a5bb5bad885693e2e5f52;
        boost::python::def("contingency_table_53063479d90a5bb5bad885693e2e5f52", function_pointer_ebdcbba919fd5e7a9af8bccab4b78e8e);
}