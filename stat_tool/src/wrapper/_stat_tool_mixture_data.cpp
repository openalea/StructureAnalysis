#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/multivariate_mixture.h>
#include <stat_tool/regression.h>
#include <stat_tool/convolution.h>
#include <stat_tool/stat_tools.h>
#include <stat_tool/mixture.h>

void _stat_tool_mixture_data()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        class ::stat_tool::DiscreteDistributionData * (::stat_tool::MixtureData::*method_pointer_22b0a6e9b8985815a70d8c14b0968354)(class ::stat_tool::StatError &, int, int) const = &::stat_tool::MixtureData::extract;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::MixtureData::*method_pointer_9c739f6c664a5fa5988c97abf15afc0d)(class ::std::basic_ostream<char, std::char_traits<char> > &, bool) const = &::stat_tool::MixtureData::ascii_data_write;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::MixtureData::*method_pointer_d6d234225317561587b622f00be4464d)(class ::std::basic_ostream<char, std::char_traits<char> > &, bool) const = &::stat_tool::MixtureData::ascii_write;
        class ::stat_tool::MultiPlotSet * (::stat_tool::MixtureData::*method_pointer_3c47dfb85330571f819725339cebc09e)() const = &::stat_tool::MixtureData::get_plotable;
        void (::stat_tool::MixtureData::*method_pointer_65b10f92a726577cbf6564411c83c993)(enum ::stat_tool::variable_nature) = &::stat_tool::MixtureData::state_variable_init;
        double (::stat_tool::MixtureData::*method_pointer_33ca36d613b75e3ebc8bbd312abb4ab5)() const = &::stat_tool::MixtureData::classification_information_computation;
        double (::stat_tool::MixtureData::*method_pointer_500bb41fb5f55d4e8a20e25372760385)() const = &::stat_tool::MixtureData::information_computation;
        void (::stat_tool::MixtureData::*method_pointer_e3d6ff48072e511caff022cc29ffe844)(int) = &::stat_tool::MixtureData::build_observation_frequency_distribution;
        void (::stat_tool::MixtureData::*method_pointer_acc63112985659a790213c952e1d2bdf)(int, int, double) = &::stat_tool::MixtureData::build_observation_histogram;
        void (::stat_tool::MixtureData::*method_pointer_d2c29ecee265561ca35772e3fb810044)(int) = &::stat_tool::MixtureData::build_observation_histogram;
        bool (::stat_tool::MixtureData::*method_pointer_2841172517a05fbc8539aeddde1ea10c)(class ::stat_tool::StatError &, int, double, double) = &::stat_tool::MixtureData::select_step;
        class ::stat_tool::Mixture * (::stat_tool::MixtureData::*method_pointer_7f8d731d35af5419bd666194a383fb26)() const = &::stat_tool::MixtureData::get_mixture;
        class ::stat_tool::FrequencyDistribution * (::stat_tool::MixtureData::*method_pointer_a03eb68b67645232b08139d24f6872bc)(int, int) const = &::stat_tool::MixtureData::get_observation_distribution;
        class ::stat_tool::Histogram * (::stat_tool::MixtureData::*method_pointer_7a2d10c5ce5952b68902366170c32ad0)(int, int) const = &::stat_tool::MixtureData::get_observation_histogram;
        double (::stat_tool::MixtureData::*method_pointer_fb0825b198655748bbb63293623f2952)() const = &::stat_tool::MixtureData::get_likelihood;
        double (::stat_tool::MixtureData::*method_pointer_297cd32e02845d91946c609b7dabe7a0)() const = &::stat_tool::MixtureData::get_restoration_likelihood;
        double (::stat_tool::MixtureData::*method_pointer_d994a13faeee5eecbd9fd6c88cba2457)() const = &::stat_tool::MixtureData::get_sample_entropy;
        double (::stat_tool::MixtureData::*method_pointer_ba471cb1da16567d8a412133323a293d)(int) const = &::stat_tool::MixtureData::get_posterior_probability;
        double (::stat_tool::MixtureData::*method_pointer_c2cf941a4bdb513dac8d83e9192ce107)(int) const = &::stat_tool::MixtureData::get_entropy;
        boost::python::class_< class ::stat_tool::MixtureData, std::shared_ptr< class ::stat_tool::MixtureData >, boost::python::bases< class ::stat_tool::Vectors > >("MixtureData", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< class ::stat_tool::Vectors const &, enum ::stat_tool::vector_transformation >())
            .def(boost::python::init< class ::stat_tool::MixtureData const &, bool, enum ::stat_tool::vector_transformation >())
            .def("extract", method_pointer_22b0a6e9b8985815a70d8c14b0968354, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("ascii_data_write", method_pointer_9c739f6c664a5fa5988c97abf15afc0d, boost::python::return_internal_reference<>())
            .def("ascii_write", method_pointer_d6d234225317561587b622f00be4464d, boost::python::return_internal_reference<>())
            .def("get_plotable", method_pointer_3c47dfb85330571f819725339cebc09e, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("state_variable_init", method_pointer_65b10f92a726577cbf6564411c83c993)
            .def("classification_information_computation", method_pointer_33ca36d613b75e3ebc8bbd312abb4ab5)
            .def("information_computation", method_pointer_500bb41fb5f55d4e8a20e25372760385)
            .def("build_observation_frequency_distribution", method_pointer_e3d6ff48072e511caff022cc29ffe844)
            .def("build_observation_histogram", method_pointer_acc63112985659a790213c952e1d2bdf)
            .def("build_observation_histogram", method_pointer_d2c29ecee265561ca35772e3fb810044)
            .def("select_step", method_pointer_2841172517a05fbc8539aeddde1ea10c)
            .def("get_mixture", method_pointer_7f8d731d35af5419bd666194a383fb26, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("get_observation_distribution", method_pointer_a03eb68b67645232b08139d24f6872bc, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("get_observation_histogram", method_pointer_7a2d10c5ce5952b68902366170c32ad0, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("get_likelihood", method_pointer_fb0825b198655748bbb63293623f2952)
            .def("get_restoration_likelihood", method_pointer_297cd32e02845d91946c609b7dabe7a0)
            .def("get_sample_entropy", method_pointer_d994a13faeee5eecbd9fd6c88cba2457)
            .def("get_posterior_probability", method_pointer_ba471cb1da16567d8a412133323a293d)
            .def("get_entropy", method_pointer_c2cf941a4bdb513dac8d83e9192ce107);
        boost::python::implicitly_convertible< std::shared_ptr< class ::stat_tool::MixtureData >, std::shared_ptr< class ::stat_tool::Vectors > >();
}