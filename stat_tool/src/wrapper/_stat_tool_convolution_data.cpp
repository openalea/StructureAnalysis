#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/discrete_mixture.h>
#include <stat_tool/stat_tools.h>
#include <stat_tool/convolution.h>

void _stat_tool_convolution_data()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        class ::stat_tool::DiscreteDistributionData * (::stat_tool::ConvolutionData::*method_pointer_9095edf2ab535d02b6e5dedd623f48d2)(class ::stat_tool::StatError &, int) const = &::stat_tool::ConvolutionData::extract;
        class ::stat_tool::Convolution * (::stat_tool::ConvolutionData::*method_pointer_3d736d0b87aa50ea9a7b6f951a72aee0)() const = &::stat_tool::ConvolutionData::get_convolution;
        int (::stat_tool::ConvolutionData::*method_pointer_88f77df9758959e5af42130b5e24340a)() const = &::stat_tool::ConvolutionData::get_nb_distribution;
        class ::stat_tool::FrequencyDistribution * (::stat_tool::ConvolutionData::*method_pointer_407c3928cc115eb7b55d1ae35f3c8cc9)(int) const = &::stat_tool::ConvolutionData::get_frequency_distribution;
        boost::python::class_< class ::stat_tool::ConvolutionData, std::shared_ptr< class ::stat_tool::ConvolutionData >, boost::python::bases< class ::stat_tool::StatInterface, class ::stat_tool::FrequencyDistribution > >("ConvolutionData", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< class ::stat_tool::FrequencyDistribution const &, int >())
            .def(boost::python::init< class ::stat_tool::Convolution const & >())
            .def(boost::python::init< class ::stat_tool::ConvolutionData const &, bool >())
            .def("extract", method_pointer_9095edf2ab535d02b6e5dedd623f48d2, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("get_convolution", method_pointer_3d736d0b87aa50ea9a7b6f951a72aee0, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("get_nb_distribution", method_pointer_88f77df9758959e5af42130b5e24340a)
            .def("get_frequency_distribution", method_pointer_407c3928cc115eb7b55d1ae35f3c8cc9, boost::python::return_value_policy< boost::python::reference_existing_object >());
        boost::python::implicitly_convertible< std::shared_ptr< class ::stat_tool::ConvolutionData >, std::shared_ptr< class ::stat_tool::StatInterface > >();
        boost::python::implicitly_convertible< std::shared_ptr< class ::stat_tool::ConvolutionData >, std::shared_ptr< class ::stat_tool::FrequencyDistribution > >();
}