#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/multivariate_mixture.h>
#include <stat_tool/regression.h>
#include <stat_tool/convolution.h>
#include <stat_tool/stat_tools.h>
#include <stat_tool/mixture.h>

void _stat_tool_multivariate_mixture_data()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        class ::stat_tool::DiscreteDistributionData * (::stat_tool::MultivariateMixtureData::*method_pointer_e8cc8ae3b7fe5d778fb576c2ac3e1f86)(class ::stat_tool::StatError &, int, int) const = &::stat_tool::MultivariateMixtureData::extract;
        class ::stat_tool::DiscreteDistributionData * (::stat_tool::MultivariateMixtureData::*method_pointer_a06334f0426158cda4050e6a9c721b28)(class ::stat_tool::StatError &, int) const = &::stat_tool::MultivariateMixtureData::extract_marginal;
        double (::stat_tool::MultivariateMixtureData::*method_pointer_1dd7cda444d35238a9a2409de29f093f)() const = &::stat_tool::MultivariateMixtureData::information_computation;
        class ::stat_tool::MultivariateMixture * (::stat_tool::MultivariateMixtureData::*method_pointer_768dbcc427815947a358a6288177ea30)() const = &::stat_tool::MultivariateMixtureData::get_mixture;
        int (::stat_tool::MultivariateMixtureData::*method_pointer_60d122f2172454d990bab9f59ac9d2f7)() const = &::stat_tool::MultivariateMixtureData::get_nb_component;
        class ::stat_tool::FrequencyDistribution * (::stat_tool::MultivariateMixtureData::*method_pointer_628611c055045c15aa3facd568cdccdd)() const = &::stat_tool::MultivariateMixtureData::get_weight;
        class ::stat_tool::FrequencyDistribution * (::stat_tool::MultivariateMixtureData::*method_pointer_9d73537d6e6655ca831f71b0021af566)(int, int) const = &::stat_tool::MultivariateMixtureData::get_component;
        boost::python::class_< class ::stat_tool::MultivariateMixtureData, std::shared_ptr< class ::stat_tool::MultivariateMixtureData >, boost::python::bases< class ::stat_tool::Vectors > >("MultivariateMixtureData", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< class ::stat_tool::Vectors const &, int >())
            .def(boost::python::init< class ::stat_tool::MultivariateMixture const & >())
            .def(boost::python::init< class ::stat_tool::MultivariateMixtureData const &, bool >())
            .def("extract", method_pointer_e8cc8ae3b7fe5d778fb576c2ac3e1f86, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("extract_marginal", method_pointer_a06334f0426158cda4050e6a9c721b28, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("information_computation", method_pointer_1dd7cda444d35238a9a2409de29f093f)
            .def("get_mixture", method_pointer_768dbcc427815947a358a6288177ea30, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("get_nb_component", method_pointer_60d122f2172454d990bab9f59ac9d2f7)
            .def("get_weight", method_pointer_628611c055045c15aa3facd568cdccdd, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("get_component", method_pointer_9d73537d6e6655ca831f71b0021af566, boost::python::return_value_policy< boost::python::reference_existing_object >());
        boost::python::implicitly_convertible< std::shared_ptr< class ::stat_tool::MultivariateMixtureData >, std::shared_ptr< class ::stat_tool::Vectors > >();
}