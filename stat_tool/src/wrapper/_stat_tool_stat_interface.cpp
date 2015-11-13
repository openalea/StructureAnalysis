#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _stat_tool_stat_interface()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::StatInterface::*method_pointer_3142c767cdb557b1a7f57cb7cdcb0d1d)(class ::std::basic_ostream<char, std::char_traits<char> > &) const = &::stat_tool::StatInterface::line_write;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::StatInterface::*method_pointer_46dda22d835a5df7afe6ca38c125d95d)(class ::std::basic_ostream<char, std::char_traits<char> > &, bool) const = &::stat_tool::StatInterface::ascii_write;
        class ::stat_tool::MultiPlotSet * (::stat_tool::StatInterface::*method_pointer_dcd5f5dbeaed52b1a46267335ba3f07b)() const = &::stat_tool::StatInterface::get_plotable;
        boost::python::class_< class ::stat_tool::StatInterface, std::shared_ptr< class ::stat_tool::StatInterface >, boost::noncopyable >("StatInterface", boost::python::no_init)            .def("line_write", method_pointer_3142c767cdb557b1a7f57cb7cdcb0d1d, boost::python::return_internal_reference<>())            .def("ascii_write", method_pointer_46dda22d835a5df7afe6ca38c125d95d, boost::python::return_internal_reference<>())            .def("get_plotable", method_pointer_dcd5f5dbeaed52b1a46267335ba3f07b, boost::python::return_value_policy< boost::python::reference_existing_object >());
}