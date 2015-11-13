#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _stat_tool_stat_error()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::StatError::*method_pointer_9a50c46e29225e6482a864d8c627313f)(class ::std::basic_ostream<char, std::char_traits<char> > &, enum ::stat_tool::error_type) const = &::stat_tool::StatError::ascii_write;
        void (::stat_tool::StatError::*method_pointer_6aca1e4328215b179690c13e4babbc77)() = &::stat_tool::StatError::init;
        int (::stat_tool::StatError::*method_pointer_428a45e45b335e5ca2005b2914d07101)() const = &::stat_tool::StatError::get_nb_error;
        int (::stat_tool::StatError::*method_pointer_9007bd748c2c5e138d772eb3fd7ec146)() const = &::stat_tool::StatError::get_max_nb_error;
        boost::python::class_< class ::stat_tool::StatError, std::shared_ptr< class ::stat_tool::StatError > >("StatError", boost::python::no_init)
            .def(boost::python::init< int >())            .def("ascii_write", method_pointer_9a50c46e29225e6482a864d8c627313f, boost::python::return_internal_reference<>())            .def("init", method_pointer_6aca1e4328215b179690c13e4babbc77)            .def("get_nb_error", method_pointer_428a45e45b335e5ca2005b2914d07101)            .def("get_max_nb_error", method_pointer_9007bd748c2c5e138d772eb3fd7ec146);
}