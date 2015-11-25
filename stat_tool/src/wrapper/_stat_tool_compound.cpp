#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/discrete_mixture.h>
#include <stat_tool/stat_tools.h>
#include <stat_tool/convolution.h>

void _stat_tool_compound()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        class ::stat_tool::CompoundData * (::stat_tool::Compound::*method_pointer_e2e33880d5e45709a3165f902811fef7)(class ::stat_tool::StatError &) const = &::stat_tool::Compound::extract_data;
        void (::stat_tool::Compound::*method_pointer_40b9d2106f415ac190bab7e1b8715f95)(int, double, bool, bool) = &::stat_tool::Compound::computation;
        class ::stat_tool::CompoundData * (::stat_tool::Compound::*method_pointer_b0550c167c9d5c97a464129d6e723e37)(class ::stat_tool::StatError &, int) const = &::stat_tool::Compound::simulation;
        class ::stat_tool::CompoundData * (::stat_tool::Compound::*method_pointer_b25b77fef5d259b79275ea462cfa8cea)() const = &::stat_tool::Compound::get_compound_data;
        class ::stat_tool::DiscreteParametric * (::stat_tool::Compound::*method_pointer_9c6954c5fb705e09ad736fb70b9ed568)() const = &::stat_tool::Compound::get_sum_distribution;
        class ::stat_tool::DiscreteParametric * (::stat_tool::Compound::*method_pointer_aef6ebd1e1e252a491bf74f73217bdf4)() const = &::stat_tool::Compound::get_distribution;
        boost::python::class_< class ::stat_tool::Compound, std::shared_ptr< class ::stat_tool::Compound >, boost::python::bases< class ::stat_tool::StatInterface, class ::stat_tool::Distribution > >("Compound", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< class ::stat_tool::DiscreteParametric const &, class ::stat_tool::DiscreteParametric const &, double >())
            .def(boost::python::init< class ::stat_tool::DiscreteParametric const &, class ::stat_tool::DiscreteParametric const &, enum ::stat_tool::compound_distribution >())
            .def(boost::python::init< class ::stat_tool::Compound const &, bool >())
            .def("extract_data", method_pointer_e2e33880d5e45709a3165f902811fef7, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("computation", method_pointer_40b9d2106f415ac190bab7e1b8715f95)
            .def("simulation", method_pointer_b0550c167c9d5c97a464129d6e723e37, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("get_compound_data", method_pointer_b25b77fef5d259b79275ea462cfa8cea, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("get_sum_distribution", method_pointer_9c6954c5fb705e09ad736fb70b9ed568, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("get_distribution", method_pointer_aef6ebd1e1e252a491bf74f73217bdf4, boost::python::return_value_policy< boost::python::reference_existing_object >());
        boost::python::implicitly_convertible< std::shared_ptr< class ::stat_tool::Compound >, std::shared_ptr< class ::stat_tool::StatInterface > >();
        boost::python::implicitly_convertible< std::shared_ptr< class ::stat_tool::Compound >, std::shared_ptr< class ::stat_tool::Distribution > >();
}