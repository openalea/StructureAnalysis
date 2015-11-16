#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/stat_tools.h>
#include <stat_tool/convolution.h>
#include <stat_tool/discrete_mixture.h>

void _stat_tool_convolution()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        class ::stat_tool::DiscreteParametricModel * (::stat_tool::Convolution::*method_pointer_877bbe1644c65ee6ab38b0d65b2525e0)(class ::stat_tool::StatError &, int) const = &::stat_tool::Convolution::extract;
        class ::stat_tool::ConvolutionData * (::stat_tool::Convolution::*method_pointer_22c6701b93ab560b966fd9b4b7acc860)(class ::stat_tool::StatError &) const = &::stat_tool::Convolution::extract_data;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::Convolution::*method_pointer_36fe1cf3107b51dbba30766ac4b96f7a)(class ::std::basic_ostream<char, std::char_traits<char> > &) const = &::stat_tool::Convolution::line_write;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::Convolution::*method_pointer_2e8a8bf8a886594b89b7f47cbf9eee09)(class ::std::basic_ostream<char, std::char_traits<char> > &, bool) const = &::stat_tool::Convolution::ascii_write;
        class ::stat_tool::MultiPlotSet * (::stat_tool::Convolution::*method_pointer_86dcef2fa3c258658f5aa5f2f7bd634c)() const = &::stat_tool::Convolution::get_plotable;
        class ::stat_tool::ConvolutionData * (::stat_tool::Convolution::*method_pointer_30d362d92bec5625af53132bae339399)(class ::stat_tool::StatError &, int) const = &::stat_tool::Convolution::simulation;
        class ::stat_tool::ConvolutionData * (::stat_tool::Convolution::*method_pointer_7ef5bd984eb15ea2a7279ab780fc6836)() const = &::stat_tool::Convolution::get_convolution_data;
        int (::stat_tool::Convolution::*method_pointer_b8c7b4200a91527cb38ad8f217089c7d)() const = &::stat_tool::Convolution::get_nb_distribution;
        class ::stat_tool::DiscreteParametric * (::stat_tool::Convolution::*method_pointer_35714d9d7c605e439c910b173beb37a6)(int) const = &::stat_tool::Convolution::get_distribution;
        boost::python::class_< class ::stat_tool::Convolution, std::shared_ptr< class ::stat_tool::Convolution >, boost::python::bases< class ::stat_tool::StatInterface, class ::stat_tool::Distribution > >("Convolution", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< class ::stat_tool::DiscreteParametric const &, class ::stat_tool::DiscreteParametric const & >())
            .def(boost::python::init< class ::stat_tool::Convolution const &, bool >())
            .def("extract", method_pointer_877bbe1644c65ee6ab38b0d65b2525e0, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("extract_data", method_pointer_22c6701b93ab560b966fd9b4b7acc860, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("line_write", method_pointer_36fe1cf3107b51dbba30766ac4b96f7a, boost::python::return_internal_reference<>())
            .def("ascii_write", method_pointer_2e8a8bf8a886594b89b7f47cbf9eee09, boost::python::return_internal_reference<>())
            .def("get_plotable", method_pointer_86dcef2fa3c258658f5aa5f2f7bd634c, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("simulation", method_pointer_30d362d92bec5625af53132bae339399, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("get_convolution_data", method_pointer_7ef5bd984eb15ea2a7279ab780fc6836, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("get_nb_distribution", method_pointer_b8c7b4200a91527cb38ad8f217089c7d)
            .def("get_distribution", method_pointer_35714d9d7c605e439c910b173beb37a6, boost::python::return_value_policy< boost::python::reference_existing_object >());
        boost::python::implicitly_convertible< std::shared_ptr< class ::stat_tool::Convolution >, std::shared_ptr< class ::stat_tool::StatInterface > >();
        boost::python::implicitly_convertible< std::shared_ptr< class ::stat_tool::Convolution >, std::shared_ptr< class ::stat_tool::Distribution > >();
}