#include <boost/python.hpp>
#include <stat_tool/char_ptr_to_string.h>
#include <stat_tool/stat_tools.h>

void _stat_tool_spreadsheet_write_606ac800ccf550858b6d60ec81c6e17e()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        bool (*function_pointer_2f7c63cf636e537394dee48973eb14d0)(class ::stat_tool::Dendrogram const &, class ::stat_tool::StatError &, class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &) = ::stat_tool::spreadsheet_write_606ac800ccf550858b6d60ec81c6e17e;
        boost::python::def("spreadsheet_write_606ac800ccf550858b6d60ec81c6e17e", function_pointer_2f7c63cf636e537394dee48973eb14d0);
}