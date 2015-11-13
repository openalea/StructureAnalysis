#include <boost/python.hpp>
#include <stat_tool/char_ptr_to_string.h>
#include <stat_tool/stat_tools.h>

void _stat_tool_spreadsheet_write_5245e9312629559f86ec2be0f532879d()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        bool (*function_pointer_ad2067bda1d851a9a1d76e9b0cf9b4de)(class ::stat_tool::StatInterface const &, class ::stat_tool::StatError &, class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &) = ::stat_tool::spreadsheet_write_5245e9312629559f86ec2be0f532879d;
        boost::python::def("spreadsheet_write_5245e9312629559f86ec2be0f532879d", function_pointer_ad2067bda1d851a9a1d76e9b0cf9b4de);
}