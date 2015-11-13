#include <boost/python.hpp>
#include <stat_tool/char_ptr_to_string.h>
#include <stat_tool/stat_tools.h>

void _stat_tool_get_label_565170fe479d5248a4b31c8011783af4()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > (*function_pointer_1d95ab919a9f562a9c519de5ce95e260)(class ::stat_tool::DistanceMatrix const &) = ::stat_tool::get_label_565170fe479d5248a4b31c8011783af4;
        boost::python::def("get_label_565170fe479d5248a4b31c8011783af4", function_pointer_1d95ab919a9f562a9c519de5ce95e260);
}