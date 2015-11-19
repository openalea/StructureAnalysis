#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/stat_tools.h>
#include <stat_tool/convolution.h>
#include <stat_tool/discrete_mixture.h>

void _stat_tool_compound_data()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        class ::stat_tool::DiscreteDistributionData * (::stat_tool::CompoundData::*method_pointer_f030d5b9a4fd569c81737fcc1a40e44d)(class ::stat_tool::StatError &, enum ::stat_tool::compound_distribution) const = &::stat_tool::CompoundData::extract;
        class ::stat_tool::Compound * (::stat_tool::CompoundData::*method_pointer_b4247691758657039c75f6e8a5005fc5)() const = &::stat_tool::CompoundData::get_compound;
        class ::stat_tool::FrequencyDistribution * (::stat_tool::CompoundData::*method_pointer_71f25a4d9a3d5dd48e19ef196515a723)() const = &::stat_tool::CompoundData::get_sum_frequency_distribution;
        class ::stat_tool::FrequencyDistribution * (::stat_tool::CompoundData::*method_pointer_7eb4fdc1dfe05a79952798f599446ac0)() const = &::stat_tool::CompoundData::get_frequency_distribution;
        boost::python::class_< class ::stat_tool::CompoundData, std::shared_ptr< class ::stat_tool::CompoundData >, boost::python::bases< class ::stat_tool::StatInterface, class ::stat_tool::FrequencyDistribution > >("CompoundData", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::Compound const & >())
            .def(boost::python::init< class ::stat_tool::Compound const & >())
            .def(boost::python::init< class ::stat_tool::CompoundData const &, bool >())
            .def("extract", method_pointer_f030d5b9a4fd569c81737fcc1a40e44d, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("get_compound", method_pointer_b4247691758657039c75f6e8a5005fc5, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("get_sum_frequency_distribution", method_pointer_71f25a4d9a3d5dd48e19ef196515a723, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("get_frequency_distribution", method_pointer_7eb4fdc1dfe05a79952798f599446ac0, boost::python::return_value_policy< boost::python::reference_existing_object >());
        boost::python::implicitly_convertible< std::shared_ptr< class ::stat_tool::CompoundData >, std::shared_ptr< class ::stat_tool::StatInterface > >();
        boost::python::implicitly_convertible< std::shared_ptr< class ::stat_tool::CompoundData >, std::shared_ptr< class ::stat_tool::FrequencyDistribution > >();
}