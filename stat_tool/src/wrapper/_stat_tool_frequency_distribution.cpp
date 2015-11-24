#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/stat_tools.h>
#include <stat_tool/convolution.h>
#include <stat_tool/discrete_mixture.h>

void _stat_tool_frequency_distribution()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        bool (::stat_tool::FrequencyDistribution::*method_pointer_1e33eb1d4cac5053a5297cdd62f11f86)(class ::stat_tool::FrequencyDistribution const &) const = &::stat_tool::FrequencyDistribution::operator==;
        bool (::stat_tool::FrequencyDistribution::*method_pointer_27baba88c70154d1a16dd9895e36ee28)(class ::stat_tool::FrequencyDistribution const &) const = &::stat_tool::FrequencyDistribution::operator!=;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::FrequencyDistribution::*method_pointer_fe56e7a238fd57bab9b0d879dad60138)(class ::std::basic_ostream<char, std::char_traits<char> > &, int, bool) const = &::stat_tool::FrequencyDistribution::ascii_print;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::FrequencyDistribution::*method_pointer_47742b5d162a572196d7afae32ba3de5)(class ::std::basic_ostream<char, std::char_traits<char> > &, bool, bool) const = &::stat_tool::FrequencyDistribution::ascii_write;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::FrequencyDistribution::*method_pointer_77c38099f1665e368624a3b1a0fdeaf0)(class ::std::basic_ostream<char, std::char_traits<char> > &, bool) const = &::stat_tool::FrequencyDistribution::spreadsheet_characteristic_print;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::FrequencyDistribution::*method_pointer_b738b069f80e5ddfb00f3c74384e9988)(class ::std::basic_ostream<char, std::char_traits<char> > &) const = &::stat_tool::FrequencyDistribution::spreadsheet_circular_characteristic_print;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::FrequencyDistribution::*method_pointer_cdfaa12eae3154d7a3b55265b7e2e800)(class ::std::basic_ostream<char, std::char_traits<char> > &, bool, bool) const = &::stat_tool::FrequencyDistribution::spreadsheet_print;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::FrequencyDistribution::*method_pointer_222a2135eb8157e19c940ec1d11d60f5)(class ::std::basic_ostream<char, std::char_traits<char> > &) const = &::stat_tool::FrequencyDistribution::plot_title_print;
        void (::stat_tool::FrequencyDistribution::*method_pointer_9c4ca4152b24578185b165f5eada378f)(class ::stat_tool::SinglePlot &) const = &::stat_tool::FrequencyDistribution::plotable_frequency_write;
        void (::stat_tool::FrequencyDistribution::*method_pointer_cdf19052baa0516db9e6c59e9ad197cd)(class ::stat_tool::SinglePlot &) const = &::stat_tool::FrequencyDistribution::plotable_mass_write;
        void (::stat_tool::FrequencyDistribution::*method_pointer_8e060e16db8b56009d6d4f42547088bd)(class ::stat_tool::SinglePlot &) const = &::stat_tool::FrequencyDistribution::plotable_survivor_write;
        double (::stat_tool::FrequencyDistribution::*method_pointer_0b8608a35c3550a99ac8a0958d802fd4)() const = &::stat_tool::FrequencyDistribution::concentration_computation;
        void (::stat_tool::FrequencyDistribution::*method_pointer_9fe9a16abbbb56babbba6a548924f162)(class ::stat_tool::Reestimation<double> const *, int) = &::stat_tool::FrequencyDistribution::update;
        class ::stat_tool::FrequencyDistribution * (::stat_tool::FrequencyDistribution::*method_pointer_5882a600801d5cc5ac78f111a990f6b4)(int) const = &::stat_tool::FrequencyDistribution::frequency_scale;
        int (::stat_tool::FrequencyDistribution::*method_pointer_2b8a299165215e44ad3ec5cb8e6c85b0)() const = &::stat_tool::FrequencyDistribution::min_interval_computation;
        class ::stat_tool::DiscreteParametric * (::stat_tool::FrequencyDistribution::*method_pointer_cada63185bcd53a9992002c7b742c69f)(int, int, bool, double) const = &::stat_tool::FrequencyDistribution::parametric_estimation;
        double (::stat_tool::FrequencyDistribution::*method_pointer_b0a0fd68ce7952998ac6c6e0c971980f)(class ::stat_tool::ContinuousParametric const &, int) const = &::stat_tool::FrequencyDistribution::likelihood_computation;
        void (::stat_tool::FrequencyDistribution::*method_pointer_01abb08d00865de783d5f214e2fb750c)(class ::stat_tool::FrequencyDistribution const &, int) = &::stat_tool::FrequencyDistribution::shift;
        void (::stat_tool::FrequencyDistribution::*method_pointer_2010daaf8405536eb6244ab6f6ae7e30)(class ::stat_tool::FrequencyDistribution const &, int, enum ::stat_tool::rounding) = &::stat_tool::FrequencyDistribution::cluster;
        class ::stat_tool::DiscreteDistributionData * (::stat_tool::FrequencyDistribution::*method_pointer_258945c32f1c5ddfafe44231fd2ce510)(class ::stat_tool::StatError &, int) const = &::stat_tool::FrequencyDistribution::shift;
        class ::stat_tool::DiscreteDistributionData * (::stat_tool::FrequencyDistribution::*method_pointer_8d4e878e460e5fe7ad1180daa4d9591e)(class ::stat_tool::StatError &, int, enum ::stat_tool::rounding) const = &::stat_tool::FrequencyDistribution::cluster;
        class ::stat_tool::DiscreteDistributionData * (::stat_tool::FrequencyDistribution::*method_pointer_2cc4c46f861e5de0a617e1d4376dbc8f)(class ::stat_tool::StatError &, double, class ::std::basic_ostream<char, std::char_traits<char> > &) const = &::stat_tool::FrequencyDistribution::cluster;
        class ::stat_tool::DiscreteDistributionData * (::stat_tool::FrequencyDistribution::*method_pointer_e924ee1e64a45a2ebf02d8b85778ec21)(class ::stat_tool::StatError &, int, int, bool) const = &::stat_tool::FrequencyDistribution::value_select;
        class ::stat_tool::MultiPlotSet * (::stat_tool::FrequencyDistribution::*method_pointer_22f6f58186ea54b6adfad7525b30c70c)() const = &::stat_tool::FrequencyDistribution::get_plotable;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::FrequencyDistribution::*method_pointer_f37ed0d0b42b5aa894702eb85129cbcb)(class ::std::basic_ostream<char, std::char_traits<char> > &) const = &::stat_tool::FrequencyDistribution::survival_ascii_write;
        class ::stat_tool::MultiPlotSet * (::stat_tool::FrequencyDistribution::*method_pointer_bed5046815d5536db150188d88618460)(class ::stat_tool::StatError &) const = &::stat_tool::FrequencyDistribution::survival_get_plotable;
        void (::stat_tool::FrequencyDistribution::*method_pointer_5e6d53b27504567f93fafd3d6edf27d4)(class ::std::basic_ostream<char, std::char_traits<char> > &, class ::stat_tool::FrequencyDistribution const &) const = &::stat_tool::FrequencyDistribution::F_comparison;
        void (::stat_tool::FrequencyDistribution::*method_pointer_7ef560611e0e57208c1aa16c12990c7c)(class ::std::basic_ostream<char, std::char_traits<char> > &, class ::stat_tool::FrequencyDistribution const &) const = &::stat_tool::FrequencyDistribution::t_comparison;
        bool (::stat_tool::FrequencyDistribution::*method_pointer_896194f1dcca5cc9aa7d08a6826cd90f)(class ::stat_tool::StatError &, class ::std::basic_ostream<char, std::char_traits<char> > &, class ::stat_tool::FrequencyDistribution const &) const = &::stat_tool::FrequencyDistribution::wilcoxon_mann_whitney_comparison;
        class ::stat_tool::DiscreteParametricModel * (::stat_tool::FrequencyDistribution::*method_pointer_0e6e66d328825e58ac42478d950343a9)(class ::stat_tool::StatError &, class ::stat_tool::DiscreteParametric const &) const = &::stat_tool::FrequencyDistribution::fit;
        class ::stat_tool::DiscreteParametricModel * (::stat_tool::FrequencyDistribution::*method_pointer_4a6b4adefa695ef284f6ff789e2957a7)(class ::stat_tool::StatError &, int, int, bool, double) const = &::stat_tool::FrequencyDistribution::parametric_estimation;
        class ::stat_tool::DiscreteParametricModel * (::stat_tool::FrequencyDistribution::*method_pointer_93be0d5f583551e48f38113575b8cdcb)(class ::stat_tool::StatError &, int, bool, double) const = &::stat_tool::FrequencyDistribution::type_parametric_estimation;
        class ::stat_tool::DiscreteMixture * (::stat_tool::FrequencyDistribution::*method_pointer_6aebbbfa54915d9983501cd6ca26dbf0)(class ::stat_tool::StatError &, class ::stat_tool::DiscreteMixture const &, class ::std::vector<bool, std::allocator<bool> > &, int, bool, bool, double) const = &::stat_tool::FrequencyDistribution::discrete_mixture_estimation;
        class ::stat_tool::DiscreteMixture * (::stat_tool::FrequencyDistribution::*method_pointer_3c7877e75df35d599ee7ea6ab878440e)(class ::stat_tool::StatError &, class ::stat_tool::DiscreteMixture const &, int, bool, bool, double) const = &::stat_tool::FrequencyDistribution::discrete_mixture_estimation;
        class ::stat_tool::DiscreteMixture * (::stat_tool::FrequencyDistribution::*method_pointer_c045947962f951478cb56d4ed8d9e77d)(class ::stat_tool::StatError &, int, class ::std::vector<int, std::allocator<int> > &, int, bool, bool, double) const = &::stat_tool::FrequencyDistribution::discrete_mixture_estimation;
        class ::stat_tool::DiscreteMixture * (::stat_tool::FrequencyDistribution::*method_pointer_d2266288a16c5cd88230776b8a347eb1)(class ::stat_tool::StatError &, class ::std::basic_ostream<char, std::char_traits<char> > &, int, int, class ::std::vector<int, std::allocator<int> > &, int, bool, bool, enum ::stat_tool::model_selection_criterion, double) const = &::stat_tool::FrequencyDistribution::discrete_mixture_estimation;
        class ::stat_tool::Convolution * (::stat_tool::FrequencyDistribution::*method_pointer_9d4ed33554365e1aab874caf0433a719)(class ::stat_tool::StatError &, class ::std::basic_ostream<char, std::char_traits<char> > &, class ::stat_tool::DiscreteParametric const &, class ::stat_tool::DiscreteParametric const &, enum ::stat_tool::estimation_criterion, int, double, enum ::stat_tool::penalty_type, enum ::stat_tool::side_effect) const = &::stat_tool::FrequencyDistribution::convolution_estimation;
        class ::stat_tool::Convolution * (::stat_tool::FrequencyDistribution::*method_pointer_089ef3e1597659f2ab6fdd581f19f25d)(class ::stat_tool::StatError &, class ::std::basic_ostream<char, std::char_traits<char> > &, class ::stat_tool::DiscreteParametric const &, int, enum ::stat_tool::estimation_criterion, int, double, enum ::stat_tool::penalty_type, enum ::stat_tool::side_effect) const = &::stat_tool::FrequencyDistribution::convolution_estimation;
        class ::stat_tool::Compound * (::stat_tool::FrequencyDistribution::*method_pointer_1b8ec7ac94b65332bc5431f9f8f7b611)(class ::stat_tool::StatError &, class ::std::basic_ostream<char, std::char_traits<char> > &, class ::stat_tool::DiscreteParametric const &, class ::stat_tool::DiscreteParametric const &, enum ::stat_tool::compound_distribution, enum ::stat_tool::estimation_criterion, int, double, enum ::stat_tool::penalty_type, enum ::stat_tool::side_effect) const = &::stat_tool::FrequencyDistribution::compound_estimation;
        class ::stat_tool::Compound * (::stat_tool::FrequencyDistribution::*method_pointer_b3e7b182e0115b8c8aab8609b120a386)(class ::stat_tool::StatError &, class ::std::basic_ostream<char, std::char_traits<char> > &, class ::stat_tool::DiscreteParametric const &, enum ::stat_tool::compound_distribution, int, enum ::stat_tool::estimation_criterion, int, double, enum ::stat_tool::penalty_type, enum ::stat_tool::side_effect) const = &::stat_tool::FrequencyDistribution::compound_estimation;
        class ::stat_tool::DiscreteParametricModel * (::stat_tool::FrequencyDistribution::*method_pointer_8f9161c281e15f15bb8e27c32d0b2368)(class ::stat_tool::StatError &, class ::std::basic_ostream<char, std::char_traits<char> > &, class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const *, class ::stat_tool::DiscreteParametric const &, enum ::stat_tool::estimation_criterion, int, enum ::stat_tool::duration_distribution_mean_estimator, double, enum ::stat_tool::penalty_type, enum ::stat_tool::side_effect, double) const = &::stat_tool::FrequencyDistribution::estimation;
        class ::stat_tool::DiscreteParametricModel * (::stat_tool::FrequencyDistribution::*method_pointer_8142dce6ba7150a09399f5593b839551)(class ::stat_tool::StatError &, class ::std::basic_ostream<char, std::char_traits<char> > &, class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const *, enum ::stat_tool::estimation_criterion, int, enum ::stat_tool::duration_distribution_mean_estimator, double, enum ::stat_tool::penalty_type, enum ::stat_tool::side_effect) const = &::stat_tool::FrequencyDistribution::estimation;
        boost::python::class_< class ::stat_tool::FrequencyDistribution, std::shared_ptr< class ::stat_tool::FrequencyDistribution >, boost::python::bases< class ::stat_tool::Reestimation<int> > >("FrequencyDistribution", boost::python::no_init)
            .def(boost::python::init< int >())
            .def(boost::python::init< class ::stat_tool::Distribution const & >())
            .def(boost::python::init< class ::std::vector<stat_tool::Reestimation<int> *, std::allocator<stat_tool::Reestimation<int> *> > const & >())
            .def(boost::python::init< class ::stat_tool::FrequencyDistribution const &, enum ::stat_tool::frequency_distribution_transformation, int, enum ::stat_tool::rounding >())
            .def(boost::python::init< class ::stat_tool::FrequencyDistribution const & >())
            .def("__eq__", method_pointer_1e33eb1d4cac5053a5297cdd62f11f86)
            .def("__neq__", method_pointer_27baba88c70154d1a16dd9895e36ee28)
            .def("ascii_print", method_pointer_fe56e7a238fd57bab9b0d879dad60138, boost::python::return_internal_reference<>())
            .def("ascii_write", method_pointer_47742b5d162a572196d7afae32ba3de5, boost::python::return_internal_reference<>())
            .def("spreadsheet_characteristic_print", method_pointer_77c38099f1665e368624a3b1a0fdeaf0, boost::python::return_internal_reference<>())
            .def("spreadsheet_circular_characteristic_print", method_pointer_b738b069f80e5ddfb00f3c74384e9988, boost::python::return_internal_reference<>())
            .def("spreadsheet_print", method_pointer_cdfaa12eae3154d7a3b55265b7e2e800, boost::python::return_internal_reference<>())
            .def("plot_title_print", method_pointer_222a2135eb8157e19c940ec1d11d60f5, boost::python::return_internal_reference<>())
            .def("plotable_frequency_write", method_pointer_9c4ca4152b24578185b165f5eada378f)
            .def("plotable_mass_write", method_pointer_cdf19052baa0516db9e6c59e9ad197cd)
            .def("plotable_survivor_write", method_pointer_8e060e16db8b56009d6d4f42547088bd)
            .def("concentration_computation", method_pointer_0b8608a35c3550a99ac8a0958d802fd4)
            .def("update", method_pointer_9fe9a16abbbb56babbba6a548924f162)
            .def("frequency_scale", method_pointer_5882a600801d5cc5ac78f111a990f6b4, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("min_interval_computation", method_pointer_2b8a299165215e44ad3ec5cb8e6c85b0)
            .def("parametric_estimation", method_pointer_cada63185bcd53a9992002c7b742c69f, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("likelihood_computation", method_pointer_b0a0fd68ce7952998ac6c6e0c971980f)
            .def("shift", method_pointer_01abb08d00865de783d5f214e2fb750c)
            .def("cluster", method_pointer_2010daaf8405536eb6244ab6f6ae7e30)
            .def("shift", method_pointer_258945c32f1c5ddfafe44231fd2ce510, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("cluster", method_pointer_8d4e878e460e5fe7ad1180daa4d9591e, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("cluster", method_pointer_2cc4c46f861e5de0a617e1d4376dbc8f, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("value_select", method_pointer_e924ee1e64a45a2ebf02d8b85778ec21, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("get_plotable", method_pointer_22f6f58186ea54b6adfad7525b30c70c, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("survival_ascii_write", method_pointer_f37ed0d0b42b5aa894702eb85129cbcb, boost::python::return_internal_reference<>())
            .def("survival_get_plotable", method_pointer_bed5046815d5536db150188d88618460, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("f__comparison", method_pointer_5e6d53b27504567f93fafd3d6edf27d4)
            .def("t_comparison", method_pointer_7ef560611e0e57208c1aa16c12990c7c)
            .def("wilcoxon_mann_whitney_comparison", method_pointer_896194f1dcca5cc9aa7d08a6826cd90f)
            .def("fit", method_pointer_0e6e66d328825e58ac42478d950343a9, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("parametric_estimation", method_pointer_4a6b4adefa695ef284f6ff789e2957a7, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("type_parametric_estimation", method_pointer_93be0d5f583551e48f38113575b8cdcb, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("discrete_mixture_estimation", method_pointer_6aebbbfa54915d9983501cd6ca26dbf0, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("discrete_mixture_estimation", method_pointer_3c7877e75df35d599ee7ea6ab878440e, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("discrete_mixture_estimation", method_pointer_c045947962f951478cb56d4ed8d9e77d, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("discrete_mixture_estimation", method_pointer_d2266288a16c5cd88230776b8a347eb1, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("convolution_estimation", method_pointer_9d4ed33554365e1aab874caf0433a719, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("convolution_estimation", method_pointer_089ef3e1597659f2ab6fdd581f19f25d, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("compound_estimation", method_pointer_1b8ec7ac94b65332bc5431f9f8f7b611, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("compound_estimation", method_pointer_b3e7b182e0115b8c8aab8609b120a386, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("estimation", method_pointer_8f9161c281e15f15bb8e27c32d0b2368, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("estimation", method_pointer_8142dce6ba7150a09399f5593b839551, boost::python::return_value_policy< boost::python::reference_existing_object >());
        boost::python::implicitly_convertible< std::shared_ptr< class ::stat_tool::FrequencyDistribution >, std::shared_ptr< class ::stat_tool::Reestimation<int> > >();
}