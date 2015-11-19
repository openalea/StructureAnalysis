#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/stat_tools.h>
#include <stat_tool/convolution.h>
#include <stat_tool/discrete_mixture.h>

void _stat_tool_discrete_parametric_model()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        class ::stat_tool::DiscreteDistributionData * (::stat_tool::DiscreteParametricModel::*method_pointer_32615e4dbd38500b914c96f21f414665)(class ::stat_tool::StatError &) const = &::stat_tool::DiscreteParametricModel::extract_data;
        class ::stat_tool::DiscreteDistributionData * (::stat_tool::DiscreteParametricModel::*method_pointer_10571acabe925b6aaff7cc61716a1496)(class ::stat_tool::StatError &, int) const = &::stat_tool::DiscreteParametricModel::simulation;
        class ::stat_tool::DiscreteDistributionData * (::stat_tool::DiscreteParametricModel::*method_pointer_4c807bb6b90358228851305f7e6bc996)() const = &::stat_tool::DiscreteParametricModel::get_frequency_distribution;
        boost::python::class_< class ::stat_tool::DiscreteParametricModel, std::shared_ptr< class ::stat_tool::DiscreteParametricModel >, boost::python::bases< class ::stat_tool::StatInterface, class ::stat_tool::DiscreteParametric > >("DiscreteParametricModel", boost::python::no_init)
            .def(boost::python::init< int, enum ::stat_tool::discrete_parametric, int, int, double, double >())
            .def(boost::python::init< enum ::stat_tool::discrete_parametric, int, int, double, double, double >())
            .def(boost::python::init< class ::stat_tool::FrequencyDistribution const & >())
            .def(boost::python::init< class ::stat_tool::Distribution const & >())
            .def(boost::python::init< class ::stat_tool::DiscreteParametric const & >())
            .def(boost::python::init< class ::stat_tool::Distribution const &, class ::stat_tool::FrequencyDistribution const * >())
            .def(boost::python::init< class ::stat_tool::DiscreteParametric const &, class ::stat_tool::FrequencyDistribution const * >())
            .def(boost::python::init< class ::stat_tool::DiscreteParametricModel const &, bool >())
            .def("extract_data", method_pointer_32615e4dbd38500b914c96f21f414665, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("simulation", method_pointer_10571acabe925b6aaff7cc61716a1496, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("get_frequency_distribution", method_pointer_4c807bb6b90358228851305f7e6bc996, boost::python::return_value_policy< boost::python::reference_existing_object >());
        boost::python::implicitly_convertible< std::shared_ptr< class ::stat_tool::DiscreteParametricModel >, std::shared_ptr< class ::stat_tool::StatInterface > >();
        boost::python::implicitly_convertible< std::shared_ptr< class ::stat_tool::DiscreteParametricModel >, std::shared_ptr< class ::stat_tool::DiscreteParametric > >();
}