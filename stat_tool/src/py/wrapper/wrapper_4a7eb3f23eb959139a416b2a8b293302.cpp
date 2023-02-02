#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::Reestimation< double > const volatile * get_pointer<class ::stat_tool::Reestimation< double > const volatile >(class ::stat_tool::Reestimation< double > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_4a7eb3f23eb959139a416b2a8b293302()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    void  (::stat_tool::Reestimation< double >::*method_pointer_a0a43540e1f75b71a89bc38dfb8496c3)(int ) = &::stat_tool::Reestimation< double >::init;
    void  (::stat_tool::Reestimation< double >::*method_pointer_0484efbfc4e851b285075ca7089d2e42)(class ::stat_tool::Reestimation< double > const &) = &::stat_tool::Reestimation< double >::copy;
    void  (::stat_tool::Reestimation< double >::*method_pointer_5fef560562df5f659de16de2e7d7f22b)() = &::stat_tool::Reestimation< double >::nb_value_computation;
    void  (::stat_tool::Reestimation< double >::*method_pointer_3e9109db93eb51119033323b7d78541c)() = &::stat_tool::Reestimation< double >::offset_computation;
    void  (::stat_tool::Reestimation< double >::*method_pointer_1696a44557aa5f8b95b2ff465c5698d7)() = &::stat_tool::Reestimation< double >::nb_element_computation;
    void  (::stat_tool::Reestimation< double >::*method_pointer_1fb3b445c362510b96bf568cbc25e465)() = &::stat_tool::Reestimation< double >::max_computation;
    double  (::stat_tool::Reestimation< double >::*method_pointer_dadcd7943bf65dd390b8c9d88d928bf9)() const = &::stat_tool::Reestimation< double >::mode_computation;
    void  (::stat_tool::Reestimation< double >::*method_pointer_7d7a4d483e9f5c2395fc5e7b6a45c82e)() = &::stat_tool::Reestimation< double >::mean_computation;
    double  (::stat_tool::Reestimation< double >::*method_pointer_09921e18c9a05dda97887e4f4909d062)(double ) const = &::stat_tool::Reestimation< double >::quantile_computation;
    void  (::stat_tool::Reestimation< double >::*method_pointer_72025524ab915fe69fd174cb9aed56d1)(bool ) = &::stat_tool::Reestimation< double >::variance_computation;
    double  (::stat_tool::Reestimation< double >::*method_pointer_929601c120585c9ca554c76aac0f7f13)(double ) const = &::stat_tool::Reestimation< double >::mean_absolute_deviation_computation;
    double  (::stat_tool::Reestimation< double >::*method_pointer_b2dd3129882d560f99a7dd0e3dcbd2ba)() const = &::stat_tool::Reestimation< double >::log_geometric_mean_computation;
    double  (::stat_tool::Reestimation< double >::*method_pointer_120abcbade6d5b3890be43312644a98c)() const = &::stat_tool::Reestimation< double >::skewness_computation;
    double  (::stat_tool::Reestimation< double >::*method_pointer_645f645aa9a5599ca39ec3543ab955c1)() const = &::stat_tool::Reestimation< double >::kurtosis_computation;
    double  (::stat_tool::Reestimation< double >::*method_pointer_e9b78065ca7f50d6bc08511a48e0c89c)() const = &::stat_tool::Reestimation< double >::information_computation;
    double  (::stat_tool::Reestimation< double >::*method_pointer_635bb10b5e43521c8c74f55556422bfb)(class ::stat_tool::Distribution const &) const = &::stat_tool::Reestimation< double >::likelihood_computation;
    void  (::stat_tool::Reestimation< double >::*method_pointer_afefaec20b875ec2b938b4e278d79dac)(class ::stat_tool::Distribution *) const = &::stat_tool::Reestimation< double >::distribution_estimation;
    double  (::stat_tool::Reestimation< double >::*method_pointer_30208c54cc5851f5bf835d4cf2786e5b)(class ::stat_tool::DiscreteParametric *, int , bool ) const = &::stat_tool::Reestimation< double >::binomial_estimation;
    double  (::stat_tool::Reestimation< double >::*method_pointer_3cfc934e930d5db2ac8dfe94fcbb05e4)(class ::stat_tool::DiscreteParametric *, int , bool , double ) const = &::stat_tool::Reestimation< double >::poisson_estimation;
    double  (::stat_tool::Reestimation< double >::*method_pointer_c615f9e5c9835224951a4564b77c873c)(class ::stat_tool::DiscreteParametric *, int , bool , double ) const = &::stat_tool::Reestimation< double >::negative_binomial_estimation;
    double  (::stat_tool::Reestimation< double >::*method_pointer_613761bddc3d560aaf95c1e0932e0811)(class ::stat_tool::DiscreteParametric *, int , bool , double ) const = &::stat_tool::Reestimation< double >::poisson_geometric_estimation;
    double  (::stat_tool::Reestimation< double >::*method_pointer_8223667131165dd8960618140cecaa26)(class ::stat_tool::DiscreteParametric *, int , bool , double , bool ) const = &::stat_tool::Reestimation< double >::parametric_estimation;
    double  (::stat_tool::Reestimation< double >::*method_pointer_769343cc51b456a095b21d33ade773e3)(class ::stat_tool::DiscreteParametric *, int , bool , double , bool ) const = &::stat_tool::Reestimation< double >::type_parametric_estimation;
    class ::stat_tool::DiscreteParametric * (::stat_tool::Reestimation< double >::*method_pointer_f18ae085a8655e00a6f93c8fd09a97fc)(int , bool , double ) const = &::stat_tool::Reestimation< double >::type_parametric_estimation;
    void  (::stat_tool::Reestimation< double >::*method_pointer_f450909b3083577b81440b484af234b2)(class ::stat_tool::Reestimation< double > const *, double ) = &::stat_tool::Reestimation< double >::equilibrium_process_combination;
    void  (::stat_tool::Reestimation< double >::*method_pointer_eb722d07c130523f8f60fdffc218c2d8)(class ::stat_tool::Reestimation< double > const *, class ::stat_tool::Distribution *, double ) const = &::stat_tool::Reestimation< double >::equilibrium_process_estimation;
    void  (::stat_tool::Reestimation< double >::*method_pointer_d3979f2affdd51c492ebee7114f6f2ef)(class ::stat_tool::ContinuousParametric *, int ) const = &::stat_tool::Reestimation< double >::gamma_estimation;
    void  (::stat_tool::Reestimation< double >::*method_pointer_ffbe1ec3da41552c96458e473744603f)(class ::stat_tool::ContinuousParametric *, int ) const = &::stat_tool::Reestimation< double >::zero_inflated_gamma_estimation;
    void  (::stat_tool::Reestimation< double >::*method_pointer_b373ca23612d53eaa9f3a77500783036)(class ::stat_tool::ContinuousParametric *) const = &::stat_tool::Reestimation< double >::inverse_gaussian_estimation;
    boost::python::class_< class ::stat_tool::Reestimation< double >, autowig::Held< class ::stat_tool::Reestimation< double > >::Type > class_4a7eb3f23eb959139a416b2a8b293302("_Reestimation_4a7eb3f23eb959139a416b2a8b293302", "", boost::python::no_init);
    class_4a7eb3f23eb959139a416b2a8b293302.def(boost::python::init< int  >(""));
    class_4a7eb3f23eb959139a416b2a8b293302.def(boost::python::init< class ::stat_tool::Reestimation< double > const & >(""));
    class_4a7eb3f23eb959139a416b2a8b293302.def("init", method_pointer_a0a43540e1f75b71a89bc38dfb8496c3, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def("copy", method_pointer_0484efbfc4e851b285075ca7089d2e42, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def("nb_value_computation", method_pointer_5fef560562df5f659de16de2e7d7f22b, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def("offset_computation", method_pointer_3e9109db93eb51119033323b7d78541c, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def("nb_element_computation", method_pointer_1696a44557aa5f8b95b2ff465c5698d7, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def("max_computation", method_pointer_1fb3b445c362510b96bf568cbc25e465, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def("mode_computation", method_pointer_dadcd7943bf65dd390b8c9d88d928bf9, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def("mean_computation", method_pointer_7d7a4d483e9f5c2395fc5e7b6a45c82e, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def("quantile_computation", method_pointer_09921e18c9a05dda97887e4f4909d062, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def("variance_computation", method_pointer_72025524ab915fe69fd174cb9aed56d1, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def("mean_absolute_deviation_computation", method_pointer_929601c120585c9ca554c76aac0f7f13, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def("log_geometric_mean_computation", method_pointer_b2dd3129882d560f99a7dd0e3dcbd2ba, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def("skewness_computation", method_pointer_120abcbade6d5b3890be43312644a98c, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def("kurtosis_computation", method_pointer_645f645aa9a5599ca39ec3543ab955c1, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def("information_computation", method_pointer_e9b78065ca7f50d6bc08511a48e0c89c, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def("likelihood_computation", method_pointer_635bb10b5e43521c8c74f55556422bfb, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def("distribution_estimation", method_pointer_afefaec20b875ec2b938b4e278d79dac, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def("binomial_estimation", method_pointer_30208c54cc5851f5bf835d4cf2786e5b, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def("poisson_estimation", method_pointer_3cfc934e930d5db2ac8dfe94fcbb05e4, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def("negative_binomial_estimation", method_pointer_c615f9e5c9835224951a4564b77c873c, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def("poisson_geometric_estimation", method_pointer_613761bddc3d560aaf95c1e0932e0811, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def("parametric_estimation", method_pointer_8223667131165dd8960618140cecaa26, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def("type_parametric_estimation", method_pointer_769343cc51b456a095b21d33ade773e3, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def("type_parametric_estimation", method_pointer_f18ae085a8655e00a6f93c8fd09a97fc, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_4a7eb3f23eb959139a416b2a8b293302.def("equilibrium_process_combination", method_pointer_f450909b3083577b81440b484af234b2, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def("equilibrium_process_estimation", method_pointer_eb722d07c130523f8f60fdffc218c2d8, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def("gamma_estimation", method_pointer_d3979f2affdd51c492ebee7114f6f2ef, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def("zero_inflated_gamma_estimation", method_pointer_ffbe1ec3da41552c96458e473744603f, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def("inverse_gaussian_estimation", method_pointer_b373ca23612d53eaa9f3a77500783036, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def_readwrite("nb_value", &::stat_tool::Reestimation< double >::nb_value, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def_readwrite("alloc_nb_value", &::stat_tool::Reestimation< double >::alloc_nb_value, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def_readwrite("offset", &::stat_tool::Reestimation< double >::offset, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def_readwrite("nb_element", &::stat_tool::Reestimation< double >::nb_element, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def_readwrite("max", &::stat_tool::Reestimation< double >::max, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def_readwrite("mean", &::stat_tool::Reestimation< double >::mean, "");
    class_4a7eb3f23eb959139a416b2a8b293302.def_readwrite("variance", &::stat_tool::Reestimation< double >::variance, "");

}