#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/char_ptr_to_string.h>
#include <stat_tool/stat_tools.h>

void _stat_tool_plot_print_eb16d8e35f49592fb427e566f41924a5()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        bool (*function_pointer_62794368bbb55ae5b51cef49671310b0)(class ::stat_tool::Distribution const &, class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &, class ::stat_tool::FrequencyDistribution const *) = ::stat_tool::plot_print_eb16d8e35f49592fb427e566f41924a5;
        boost::python::def("plot_print_eb16d8e35f49592fb427e566f41924a5", function_pointer_62794368bbb55ae5b51cef49671310b0);
}