#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/char_ptr_to_string.h>
#include <stat_tool/stat_tools.h>

void _stat_tool_spreadsheet_write_b486a5f6069e5ee1b978538a2f8fee74()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        bool (*function_pointer_67c3d10f782a5fd59713d89c7d5acf43)(class ::stat_tool::MixtureData const &, class ::stat_tool::StatError &, class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &) = ::stat_tool::spreadsheet_write_b486a5f6069e5ee1b978538a2f8fee74;
        boost::python::def("spreadsheet_write_b486a5f6069e5ee1b978538a2f8fee74", function_pointer_67c3d10f782a5fd59713d89c7d5acf43);
}