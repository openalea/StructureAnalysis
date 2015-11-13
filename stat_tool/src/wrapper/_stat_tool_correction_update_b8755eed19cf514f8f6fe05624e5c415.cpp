#include <boost/python.hpp>
#include <stat_tool/char_ptr_to_string.h>
#include <stat_tool/stat_tools.h>

void _stat_tool_correction_update_b8755eed19cf514f8f6fe05624e5c415()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        void (*function_pointer_acc87474f0f055068bf3a346072d44c8)(class ::stat_tool::StatError &, class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &, class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &, int, int) = ::stat_tool::correction_update_b8755eed19cf514f8f6fe05624e5c415;
        boost::python::def("correction_update_b8755eed19cf514f8f6fe05624e5c415", function_pointer_acc87474f0f055068bf3a346072d44c8);
}