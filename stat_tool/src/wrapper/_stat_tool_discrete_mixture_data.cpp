#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/stat_tools.h>
#include <stat_tool/convolution.h>
#include <stat_tool/discrete_mixture.h>

void _stat_tool_discrete_mixture_data()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        class ::stat_tool::DiscreteDistributionData * (::stat_tool::DiscreteMixtureData::*method_pointer_7d5bf00793c45f8ca8afc677d8825fd7)(class ::stat_tool::StatError &, int) const = &::stat_tool::DiscreteMixtureData::extract;
        double (::stat_tool::DiscreteMixtureData::*method_pointer_4bb3b235714c596f88fa6dbda61c3583)() const = &::stat_tool::DiscreteMixtureData::information_computation;
        class ::stat_tool::DiscreteMixture * (::stat_tool::DiscreteMixtureData::*method_pointer_87c6ed1036fe58e8a09e4a6b3dea1ac9)() const = &::stat_tool::DiscreteMixtureData::get_mixture;
        int (::stat_tool::DiscreteMixtureData::*method_pointer_ccc8e0ad33b95412bd8532b25f0424ff)() const = &::stat_tool::DiscreteMixtureData::get_nb_component;
        class ::stat_tool::FrequencyDistribution * (::stat_tool::DiscreteMixtureData::*method_pointer_709feca86c75592499578eefef265ea6)() const = &::stat_tool::DiscreteMixtureData::get_weight;
        class ::stat_tool::FrequencyDistribution * (::stat_tool::DiscreteMixtureData::*method_pointer_bc26366684395f8c85296fb6178556d0)(int) const = &::stat_tool::DiscreteMixtureData::get_component;
        boost::python::class_< class ::stat_tool::DiscreteMixtureData, std::shared_ptr< class ::stat_tool::DiscreteMixtureData >, boost::python::bases< class ::stat_tool::StatInterface, class ::stat_tool::FrequencyDistribution > >("DiscreteMixtureData", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< class ::stat_tool::FrequencyDistribution const &, int >())
            .def(boost::python::init< class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::DiscreteMixture const * >())
            .def(boost::python::init< class ::stat_tool::DiscreteMixture const & >())
            .def(boost::python::init< class ::stat_tool::DiscreteMixtureData const &, bool >())
            .def("extract", method_pointer_7d5bf00793c45f8ca8afc677d8825fd7, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("information_computation", method_pointer_4bb3b235714c596f88fa6dbda61c3583)
            .def("get_mixture", method_pointer_87c6ed1036fe58e8a09e4a6b3dea1ac9, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("get_nb_component", method_pointer_ccc8e0ad33b95412bd8532b25f0424ff)
            .def("get_weight", method_pointer_709feca86c75592499578eefef265ea6, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("get_component", method_pointer_bc26366684395f8c85296fb6178556d0, boost::python::return_value_policy< boost::python::reference_existing_object >());
        boost::python::implicitly_convertible< std::shared_ptr< class ::stat_tool::DiscreteMixtureData >, std::shared_ptr< class ::stat_tool::StatInterface > >();
        boost::python::implicitly_convertible< std::shared_ptr< class ::stat_tool::DiscreteMixtureData >, std::shared_ptr< class ::stat_tool::FrequencyDistribution > >();
}