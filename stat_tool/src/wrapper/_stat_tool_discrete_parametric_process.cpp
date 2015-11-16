#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/stat_tools.h>
#include <stat_tool/convolution.h>
#include <stat_tool/discrete_mixture.h>
#include <stat_tool/markovian.h>

void _stat_tool_discrete_parametric_process()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        void (::stat_tool::DiscreteParametricProcess::*method_pointer_0659c17596ec518882570be258267c93)(class ::stat_tool::DiscreteParametricProcess const &) = &::stat_tool::DiscreteParametricProcess::copy;
        void (::stat_tool::DiscreteParametricProcess::*method_pointer_850a40312e435b39ad131196d1a87f09)() = &::stat_tool::DiscreteParametricProcess::remove;
        void (::stat_tool::DiscreteParametricProcess::*method_pointer_1f5be9cbc73c53a5ad2264314197ec16)() = &::stat_tool::DiscreteParametricProcess::nb_value_computation;
        int (::stat_tool::DiscreteParametricProcess::*method_pointer_bf5091e7332b5c1f80763c568d769400)() const = &::stat_tool::DiscreteParametricProcess::nb_parameter_computation;
        double (::stat_tool::DiscreteParametricProcess::*method_pointer_464f8d86b1d95924b46b13dd966b6773)(class ::stat_tool::Distribution *) const = &::stat_tool::DiscreteParametricProcess::mean_computation;
        double (::stat_tool::DiscreteParametricProcess::*method_pointer_2fd07342c4735ca5ae18a3164e60c237)(class ::stat_tool::Distribution *, double) const = &::stat_tool::DiscreteParametricProcess::variance_computation;
        class ::stat_tool::Distribution * (::stat_tool::DiscreteParametricProcess::*method_pointer_0a3b99321bbe5f59907fd46f2775d1c3)(class ::stat_tool::Distribution *) = &::stat_tool::DiscreteParametricProcess::mixture_computation;
        void (::stat_tool::DiscreteParametricProcess::*method_pointer_ef05876447e35b788dd014ff886c147b)() = &::stat_tool::DiscreteParametricProcess::init;
        boost::python::class_< class ::stat_tool::DiscreteParametricProcess, std::shared_ptr< class ::stat_tool::DiscreteParametricProcess > >("DiscreteParametricProcess", boost::python::no_init)
            .def(boost::python::init< int, int >())
            .def(boost::python::init< class ::stat_tool::DiscreteParametricProcess const & >())
            .def("copy", method_pointer_0659c17596ec518882570be258267c93)
            .def("remove", method_pointer_850a40312e435b39ad131196d1a87f09)
            .def("nb_value_computation", method_pointer_1f5be9cbc73c53a5ad2264314197ec16)
            .def("nb_parameter_computation", method_pointer_bf5091e7332b5c1f80763c568d769400)
            .def("mean_computation", method_pointer_464f8d86b1d95924b46b13dd966b6773)
            .def("variance_computation", method_pointer_2fd07342c4735ca5ae18a3164e60c237)
            .def("mixture_computation", method_pointer_0a3b99321bbe5f59907fd46f2775d1c3, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("init", method_pointer_ef05876447e35b788dd014ff886c147b)
            .def_readwrite("nb_state", &::stat_tool::DiscreteParametricProcess::nb_state)
            .def_readwrite("nb_value", &::stat_tool::DiscreteParametricProcess::nb_value)
            .def_readwrite("weight", &::stat_tool::DiscreteParametricProcess::weight)
            .def_readwrite("mixture", &::stat_tool::DiscreteParametricProcess::mixture)
            .def_readwrite("restoration_weight", &::stat_tool::DiscreteParametricProcess::restoration_weight)
            .def_readwrite("restoration_mixture", &::stat_tool::DiscreteParametricProcess::restoration_mixture);
}