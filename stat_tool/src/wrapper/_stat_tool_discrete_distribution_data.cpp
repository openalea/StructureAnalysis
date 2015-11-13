#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/stat_tools.h>
#include <stat_tool/convolution.h>
#include <stat_tool/discrete_mixture.h>

void _stat_tool_discrete_distribution_data()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        class ::stat_tool::DiscreteParametricModel * (::stat_tool::DiscreteDistributionData::*method_pointer_09f4b4034c225520bfc427bc9432f05a)(class ::stat_tool::StatError &) const = &::stat_tool::DiscreteDistributionData::extract_model;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::DiscreteDistributionData::*method_pointer_215ec6ce42c4555896b7dba2d4cbb279)(class ::std::basic_ostream<char, std::char_traits<char> > &) const = &::stat_tool::DiscreteDistributionData::line_write;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::DiscreteDistributionData::*method_pointer_3c58cd08063e56a28914d1f1fe842bb0)(class ::std::basic_ostream<char, std::char_traits<char> > &, bool) const = &::stat_tool::DiscreteDistributionData::ascii_write;
        class ::stat_tool::MultiPlotSet * (::stat_tool::DiscreteDistributionData::*method_pointer_3f4f4826dfb5525f8f9d15225dfdf431)() const = &::stat_tool::DiscreteDistributionData::get_plotable;
        class ::stat_tool::DiscreteParametric * (::stat_tool::DiscreteDistributionData::*method_pointer_3dc769f178ff5cacaf81003605ec49dc)() const = &::stat_tool::DiscreteDistributionData::get_distribution;
        boost::python::class_< class ::stat_tool::DiscreteDistributionData, std::shared_ptr< class ::stat_tool::DiscreteDistributionData >, boost::python::bases< class ::stat_tool::StatInterface, class ::stat_tool::FrequencyDistribution > >("DiscreteDistributionData", boost::python::no_init)
            .def(boost::python::init< int >())
            .def(boost::python::init< class ::stat_tool::Distribution const & >())
            .def(boost::python::init< class ::stat_tool::FrequencyDistribution const & >())
            .def(boost::python::init< class ::stat_tool::FrequencyDistribution const &, enum ::stat_tool::frequency_distribution_transformation, int, enum ::stat_tool::rounding >())
            .def(boost::python::init< class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::Distribution const * >())
            .def(boost::python::init< class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::DiscreteParametric const * >())
            .def(boost::python::init< class ::stat_tool::DiscreteDistributionData const &, bool >())            .def("extract_model", method_pointer_09f4b4034c225520bfc427bc9432f05a, boost::python::return_value_policy< boost::python::reference_existing_object >())            .def("line_write", method_pointer_215ec6ce42c4555896b7dba2d4cbb279, boost::python::return_internal_reference<>())            .def("ascii_write", method_pointer_3c58cd08063e56a28914d1f1fe842bb0, boost::python::return_internal_reference<>())            .def("get_plotable", method_pointer_3f4f4826dfb5525f8f9d15225dfdf431, boost::python::return_value_policy< boost::python::reference_existing_object >())            .def("get_distribution", method_pointer_3dc769f178ff5cacaf81003605ec49dc, boost::python::return_value_policy< boost::python::reference_existing_object >());
        boost::python::implicitly_convertible< std::shared_ptr< class ::stat_tool::DiscreteDistributionData >, std::shared_ptr< class ::stat_tool::StatInterface > >();
        boost::python::implicitly_convertible< std::shared_ptr< class ::stat_tool::DiscreteDistributionData >, std::shared_ptr< class ::stat_tool::FrequencyDistribution > >();
}