#include <boost/python.hpp>
#include <stat_tool/char_ptr_to_string.h>
#include <stat_tool/stat_tools.h>

void _stat_tool_ascii_write_0d934e54b58f523aac0af7ab3df12737()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        bool (*function_pointer_c12a17ef188d5ba981379c4e4f7d68c7)(class ::stat_tool::DistanceMatrix const &, class ::stat_tool::StatError &, class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &, bool) = ::stat_tool::ascii_write_0d934e54b58f523aac0af7ab3df12737;
        boost::python::def("ascii_write_0d934e54b58f523aac0af7ab3df12737", function_pointer_c12a17ef188d5ba981379c4e4f7d68c7);
}