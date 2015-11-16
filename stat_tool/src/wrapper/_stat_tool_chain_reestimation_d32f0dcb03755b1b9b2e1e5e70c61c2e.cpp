#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>
#include <stat_tool/chain_reestimation.h>

void _stat_tool_chain_reestimation_d32f0dcb03755b1b9b2e1e5e70c61c2e()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        void (::stat_tool::ChainReestimation<int>::*method_pointer_6075593975375673a3b0c7e3b22cff24)() = &::stat_tool::ChainReestimation<int>::init;
        void (::stat_tool::ChainReestimation<int>::*method_pointer_8747980cd3bd5a45ac85c9a3644bd0ef)(class ::stat_tool::ChainReestimation<int> const &) = &::stat_tool::ChainReestimation<int>::copy;
        void (::stat_tool::ChainReestimation<int>::*method_pointer_c27152443a63579581036b6f2788d6c6)() = &::stat_tool::ChainReestimation<int>::remove;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::ChainReestimation<int>::*method_pointer_c86942b6943e5e4489c08e90caa4d02e)(class ::std::basic_ostream<char, std::char_traits<char> > &) const = &::stat_tool::ChainReestimation<int>::print;
        boost::python::class_< class ::stat_tool::ChainReestimation<int>, std::shared_ptr< class ::stat_tool::ChainReestimation<int> > >("_ChainReestimation_d32f0dcb03755b1b9b2e1e5e70c61c2e", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< enum ::stat_tool::process_type, int, int, bool >())
            .def(boost::python::init< class ::stat_tool::ChainReestimation<int> const & >())
            .def("init", method_pointer_6075593975375673a3b0c7e3b22cff24)
            .def("copy", method_pointer_8747980cd3bd5a45ac85c9a3644bd0ef)
            .def("remove", method_pointer_c27152443a63579581036b6f2788d6c6)
            .def("print", method_pointer_c86942b6943e5e4489c08e90caa4d02e, boost::python::return_internal_reference<>())
            .def_readwrite("type", &::stat_tool::ChainReestimation<int>::type)
            .def_readwrite("nb_state", &::stat_tool::ChainReestimation<int>::nb_state)
            .def_readwrite("nb_row", &::stat_tool::ChainReestimation<int>::nb_row);
}