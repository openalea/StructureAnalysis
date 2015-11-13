#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>
#include <stat_tool/markovian.h>

void _stat_tool_chain_data()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        int (::stat_tool::ChainData::*method_pointer_4f13feefcea2566cb8fdaa172ed5ab53)() const = &::stat_tool::ChainData::nb_parameter_computation;
        void (::stat_tool::ChainData::*method_pointer_576a1605ed7f5cf090c80db6f04828bf)(class ::stat_tool::Chain &) const = &::stat_tool::ChainData::estimation;
        boost::python::class_< class ::stat_tool::ChainData, std::shared_ptr< class ::stat_tool::ChainData >, boost::python::bases< class ::stat_tool::ChainReestimation<int> > >("ChainData", boost::python::no_init)
            .def(boost::python::init< enum ::stat_tool::process_type, int, int, bool >())
            .def(boost::python::init< class ::stat_tool::ChainData const & >())            .def("nb_parameter_computation", method_pointer_4f13feefcea2566cb8fdaa172ed5ab53)            .def("estimation", method_pointer_576a1605ed7f5cf090c80db6f04828bf);
        boost::python::implicitly_convertible< std::shared_ptr< class ::stat_tool::ChainData >, std::shared_ptr< class ::stat_tool::ChainReestimation<int> > >();
}