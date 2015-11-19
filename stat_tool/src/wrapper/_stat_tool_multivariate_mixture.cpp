#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/multivariate_mixture.h>
#include <stat_tool/regression.h>
#include <stat_tool/convolution.h>
#include <stat_tool/stat_tools.h>
#include <stat_tool/mixture.h>

void _stat_tool_multivariate_mixture()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        class ::stat_tool::DiscreteParametricModel * (::stat_tool::MultivariateMixture::*method_pointer_1963b3f8041454a0a7a97dd725b93e86)(class ::stat_tool::StatError &, int, int) const = &::stat_tool::MultivariateMixture::extract_parametric_model;
        class ::stat_tool::Distribution * (::stat_tool::MultivariateMixture::*method_pointer_7b7a6d59de805fdca5e66b8347871266)(class ::stat_tool::StatError &, int, int) const = &::stat_tool::MultivariateMixture::extract_categorical_model;
        class ::stat_tool::Distribution * (::stat_tool::MultivariateMixture::*method_pointer_0bc6a3bbe55f51d5b2ee84e3f9a09791)(class ::stat_tool::StatError &, int) const = &::stat_tool::MultivariateMixture::extract_distribution;
        class ::stat_tool::MultivariateMixtureData * (::stat_tool::MultivariateMixture::*method_pointer_f538f32fc03a5eb3a55cf97b6ac6b92a)(class ::stat_tool::StatError &) const = &::stat_tool::MultivariateMixture::extract_data;
        double (::stat_tool::MultivariateMixture::*method_pointer_71f1fbb2602759d4b0bb16fc5694a80f)(class ::stat_tool::Vectors const &, bool) const = &::stat_tool::MultivariateMixture::likelihood_computation;
        class ::stat_tool::MultivariateMixtureData * (::stat_tool::MultivariateMixture::*method_pointer_93428291e97650838a9c17d6894d26e5)(class ::stat_tool::StatError &, int) const = &::stat_tool::MultivariateMixture::simulation;
        class ::stat_tool::MultivariateMixtureData * (::stat_tool::MultivariateMixture::*method_pointer_5236af2599fd5fe187872fa02b640655)(class ::stat_tool::StatError &, class ::stat_tool::Vectors const &, bool) const = &::stat_tool::MultivariateMixture::cluster;
        bool (::stat_tool::MultivariateMixture::*method_pointer_a7d34676826857ca81186bc07c31304a)(int) const = &::stat_tool::MultivariateMixture::is_parametric;
        class ::stat_tool::MultivariateMixtureData * (::stat_tool::MultivariateMixture::*method_pointer_6d0df484ac0c541799eac5910fe6b979)() const = &::stat_tool::MultivariateMixture::get_mixture_data;
        int (::stat_tool::MultivariateMixture::*method_pointer_f02035b89a8b553eb304584488ccf307)() const = &::stat_tool::MultivariateMixture::get_nb_component;
        int (::stat_tool::MultivariateMixture::*method_pointer_cbdbed965dd358c89762acc255133aa8)() const = &::stat_tool::MultivariateMixture::get_nb_variable;
        class ::stat_tool::DiscreteParametric * (::stat_tool::MultivariateMixture::*method_pointer_7aed7d0db9935a4ea2a1c374faa720da)() const = &::stat_tool::MultivariateMixture::get_weight;
        class ::stat_tool::DiscreteParametricProcess * (::stat_tool::MultivariateMixture::*method_pointer_0d97c07e61595fd68e800b5d26ad4c51)(int) const = &::stat_tool::MultivariateMixture::get_parametric_process;
        class ::stat_tool::CategoricalProcess * (::stat_tool::MultivariateMixture::*method_pointer_780e2a118e8e56ca9fbf57b628e8d199)(int) const = &::stat_tool::MultivariateMixture::get_categorical_process;
        class ::stat_tool::DiscreteParametric * (::stat_tool::MultivariateMixture::*method_pointer_f5c98428d2b85f62842b6f386dd00cc1)(int, int) const = &::stat_tool::MultivariateMixture::get_parametric_component;
        class ::stat_tool::Distribution * (::stat_tool::MultivariateMixture::*method_pointer_9273e04671165096b543ef784b8f49dc)(int, int) const = &::stat_tool::MultivariateMixture::get_categorical_component;
        boost::python::class_< class ::stat_tool::MultivariateMixture, std::shared_ptr< class ::stat_tool::MultivariateMixture >, boost::python::bases< class ::stat_tool::StatInterface > >("MultivariateMixture", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< class ::stat_tool::MultivariateMixture const &, bool >())
            .def("extract_parametric_model", method_pointer_1963b3f8041454a0a7a97dd725b93e86, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("extract_categorical_model", method_pointer_7b7a6d59de805fdca5e66b8347871266, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("extract_distribution", method_pointer_0bc6a3bbe55f51d5b2ee84e3f9a09791, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("extract_data", method_pointer_f538f32fc03a5eb3a55cf97b6ac6b92a, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("likelihood_computation", method_pointer_71f1fbb2602759d4b0bb16fc5694a80f)
            .def("simulation", method_pointer_93428291e97650838a9c17d6894d26e5, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("cluster", method_pointer_5236af2599fd5fe187872fa02b640655, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("is_parametric", method_pointer_a7d34676826857ca81186bc07c31304a)
            .def("get_mixture_data", method_pointer_6d0df484ac0c541799eac5910fe6b979, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("get_nb_component", method_pointer_f02035b89a8b553eb304584488ccf307)
            .def("get_nb_variable", method_pointer_cbdbed965dd358c89762acc255133aa8)
            .def("get_weight", method_pointer_7aed7d0db9935a4ea2a1c374faa720da, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("get_parametric_process", method_pointer_0d97c07e61595fd68e800b5d26ad4c51, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("get_categorical_process", method_pointer_780e2a118e8e56ca9fbf57b628e8d199, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("get_parametric_component", method_pointer_f5c98428d2b85f62842b6f386dd00cc1, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("get_categorical_component", method_pointer_9273e04671165096b543ef784b8f49dc, boost::python::return_value_policy< boost::python::reference_existing_object >());
        boost::python::implicitly_convertible< std::shared_ptr< class ::stat_tool::MultivariateMixture >, std::shared_ptr< class ::stat_tool::StatInterface > >();
}