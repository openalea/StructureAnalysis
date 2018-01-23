#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::Reestimation< int > const volatile * get_pointer<class ::stat_tool::Reestimation< int > const volatile >(class ::stat_tool::Reestimation< int > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_710d7ee5573c5d7f8f9127a08b4f3dfd()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    void  (::stat_tool::Reestimation< int >::*method_pointer_09f5f68263f3506d9d082e14dcc2a6e9)(int ) = &::stat_tool::Reestimation< int >::init;
    void  (::stat_tool::Reestimation< int >::*method_pointer_0e9d7d4cc2ae53f5bab80c00665ca2d4)(class ::stat_tool::Reestimation< int > const &) = &::stat_tool::Reestimation< int >::copy;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::Reestimation< int >::*method_pointer_a4aad21b6790571c93675176c29392af)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &, bool , bool ) const = &::stat_tool::Reestimation< int >::ascii_characteristic_print;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::Reestimation< int >::*method_pointer_39af36c4a73953f3a40f16210221bc8d)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &, bool ) const = &::stat_tool::Reestimation< int >::ascii_circular_characteristic_print;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::Reestimation< int >::*method_pointer_9eac2967c2a65b95a6ec2efce075edd3)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &) const = &::stat_tool::Reestimation< int >::print;
    void  (::stat_tool::Reestimation< int >::*method_pointer_08c79701de3f58609887e2fbf225b2b5)() = &::stat_tool::Reestimation< int >::nb_value_computation;
    void  (::stat_tool::Reestimation< int >::*method_pointer_fae34c98fafd5884bb3fff56394048be)() = &::stat_tool::Reestimation< int >::offset_computation;
    void  (::stat_tool::Reestimation< int >::*method_pointer_4f0ea42a7b0353c2815b56f8a6cb2145)() = &::stat_tool::Reestimation< int >::nb_element_computation;
    void  (::stat_tool::Reestimation< int >::*method_pointer_f7e8187bd6f251f991a20a81eafc5da1)() = &::stat_tool::Reestimation< int >::max_computation;
    double  (::stat_tool::Reestimation< int >::*method_pointer_4b825a60b4e75353b3f1b5d02cbc0db5)() const = &::stat_tool::Reestimation< int >::mode_computation;
    void  (::stat_tool::Reestimation< int >::*method_pointer_403634937e735ba29c79f285687fa465)() = &::stat_tool::Reestimation< int >::mean_computation;
    double  (::stat_tool::Reestimation< int >::*method_pointer_da53bfee85e758bc937d0b3cc033b19c)(double ) const = &::stat_tool::Reestimation< int >::quantile_computation;
    void  (::stat_tool::Reestimation< int >::*method_pointer_aec14b5616e452118d6e1cca66fa7784)(bool ) = &::stat_tool::Reestimation< int >::variance_computation;
    double  (::stat_tool::Reestimation< int >::*method_pointer_60831f20d2f35b52bf21cc8eb00e9b45)(double ) const = &::stat_tool::Reestimation< int >::mean_absolute_deviation_computation;
    double  (::stat_tool::Reestimation< int >::*method_pointer_c8307a9f72c658449a5b00afae2e395b)() const = &::stat_tool::Reestimation< int >::log_geometric_mean_computation;
    double  (::stat_tool::Reestimation< int >::*method_pointer_ed267b607496559aa65188c0fc8f2383)() const = &::stat_tool::Reestimation< int >::skewness_computation;
    double  (::stat_tool::Reestimation< int >::*method_pointer_83cd780e01845844b1785a28fc783e4d)() const = &::stat_tool::Reestimation< int >::kurtosis_computation;
    double  (::stat_tool::Reestimation< int >::*method_pointer_e966d769588a5879b0837a21328c359e)() const = &::stat_tool::Reestimation< int >::information_computation;
    double  (::stat_tool::Reestimation< int >::*method_pointer_da2fe2021d4356b2a8f8b8d469c66d88)(class ::stat_tool::Distribution const &) const = &::stat_tool::Reestimation< int >::likelihood_computation;
    void  (::stat_tool::Reestimation< int >::*method_pointer_fa2b802f19cc55a78f2efdbc4e1531ec)(class ::stat_tool::Distribution *) const = &::stat_tool::Reestimation< int >::distribution_estimation;
    double  (::stat_tool::Reestimation< int >::*method_pointer_10fcb109ba1350e080f108b31aed30b2)(class ::stat_tool::DiscreteParametric *, int , bool ) const = &::stat_tool::Reestimation< int >::binomial_estimation;
    double  (::stat_tool::Reestimation< int >::*method_pointer_6decafe5f5115415a4a49f92995a09c2)(class ::stat_tool::DiscreteParametric *, int , bool , double ) const = &::stat_tool::Reestimation< int >::poisson_estimation;
    double  (::stat_tool::Reestimation< int >::*method_pointer_c33901e482005cf98d4f2c93a931f092)(class ::stat_tool::DiscreteParametric *, int , bool , double ) const = &::stat_tool::Reestimation< int >::negative_binomial_estimation;
    double  (::stat_tool::Reestimation< int >::*method_pointer_e059770a2ead5e73a90f19d59199a4e7)(class ::stat_tool::DiscreteParametric *, int , bool , double ) const = &::stat_tool::Reestimation< int >::poisson_geometric_estimation;
    double  (::stat_tool::Reestimation< int >::*method_pointer_55738da3ea3a534d882b1830a3725082)(class ::stat_tool::DiscreteParametric *, int , bool , double , bool ) const = &::stat_tool::Reestimation< int >::parametric_estimation;
    double  (::stat_tool::Reestimation< int >::*method_pointer_d0a58313593258579c53cdb1d134ef24)(class ::stat_tool::DiscreteParametric *, int , bool , double , bool ) const = &::stat_tool::Reestimation< int >::type_parametric_estimation;
    class ::stat_tool::DiscreteParametric * (::stat_tool::Reestimation< int >::*method_pointer_40c64aef492b5f2f908b1fe767f51ffb)(int , bool , double ) const = &::stat_tool::Reestimation< int >::type_parametric_estimation;
    void  (::stat_tool::Reestimation< int >::*method_pointer_969025eb71ba524086575bbcbc0e90ae)(class ::stat_tool::Reestimation< int > const *, double ) = &::stat_tool::Reestimation< int >::equilibrium_process_combination;
    void  (::stat_tool::Reestimation< int >::*method_pointer_b164d588de2f5515ae36b350e22788c6)(class ::stat_tool::Reestimation< int > const *, class ::stat_tool::Distribution *, double ) const = &::stat_tool::Reestimation< int >::equilibrium_process_estimation;
    void  (::stat_tool::Reestimation< int >::*method_pointer_b7f0f4bf43515f50b41c59907193dea1)(class ::stat_tool::ContinuousParametric *, int ) const = &::stat_tool::Reestimation< int >::gamma_estimation;
    void  (::stat_tool::Reestimation< int >::*method_pointer_636052c7c86a5d409a449b97991cb254)(class ::stat_tool::ContinuousParametric *, int ) const = &::stat_tool::Reestimation< int >::zero_inflated_gamma_estimation;
    void  (::stat_tool::Reestimation< int >::*method_pointer_7a6825284cb7585d897e8ae2847483a9)(class ::stat_tool::ContinuousParametric *) const = &::stat_tool::Reestimation< int >::inverse_gaussian_estimation;
    boost::python::class_< class ::stat_tool::Reestimation< int >, autowig::Held< class ::stat_tool::Reestimation< int > >::Type > class_710d7ee5573c5d7f8f9127a08b4f3dfd("_Reestimation_710d7ee5573c5d7f8f9127a08b4f3dfd", "", boost::python::no_init);
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def(boost::python::init< int  >(""));
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def(boost::python::init< class ::stat_tool::Reestimation< int > const & >(""));
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("init", method_pointer_09f5f68263f3506d9d082e14dcc2a6e9, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("copy", method_pointer_0e9d7d4cc2ae53f5bab80c00665ca2d4, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("ascii_characteristic_print", method_pointer_a4aad21b6790571c93675176c29392af, boost::python::return_internal_reference<>(), "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("ascii_circular_characteristic_print", method_pointer_39af36c4a73953f3a40f16210221bc8d, boost::python::return_internal_reference<>(), "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("print", method_pointer_9eac2967c2a65b95a6ec2efce075edd3, boost::python::return_internal_reference<>(), "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("nb_value_computation", method_pointer_08c79701de3f58609887e2fbf225b2b5, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("offset_computation", method_pointer_fae34c98fafd5884bb3fff56394048be, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("nb_element_computation", method_pointer_4f0ea42a7b0353c2815b56f8a6cb2145, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("max_computation", method_pointer_f7e8187bd6f251f991a20a81eafc5da1, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("mode_computation", method_pointer_4b825a60b4e75353b3f1b5d02cbc0db5, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("mean_computation", method_pointer_403634937e735ba29c79f285687fa465, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("quantile_computation", method_pointer_da53bfee85e758bc937d0b3cc033b19c, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("variance_computation", method_pointer_aec14b5616e452118d6e1cca66fa7784, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("mean_absolute_deviation_computation", method_pointer_60831f20d2f35b52bf21cc8eb00e9b45, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("log_geometric_mean_computation", method_pointer_c8307a9f72c658449a5b00afae2e395b, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("skewness_computation", method_pointer_ed267b607496559aa65188c0fc8f2383, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("kurtosis_computation", method_pointer_83cd780e01845844b1785a28fc783e4d, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("information_computation", method_pointer_e966d769588a5879b0837a21328c359e, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("likelihood_computation", method_pointer_da2fe2021d4356b2a8f8b8d469c66d88, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("distribution_estimation", method_pointer_fa2b802f19cc55a78f2efdbc4e1531ec, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("binomial_estimation", method_pointer_10fcb109ba1350e080f108b31aed30b2, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("poisson_estimation", method_pointer_6decafe5f5115415a4a49f92995a09c2, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("negative_binomial_estimation", method_pointer_c33901e482005cf98d4f2c93a931f092, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("poisson_geometric_estimation", method_pointer_e059770a2ead5e73a90f19d59199a4e7, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("parametric_estimation", method_pointer_55738da3ea3a534d882b1830a3725082, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("type_parametric_estimation", method_pointer_d0a58313593258579c53cdb1d134ef24, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("type_parametric_estimation", method_pointer_40c64aef492b5f2f908b1fe767f51ffb, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("equilibrium_process_combination", method_pointer_969025eb71ba524086575bbcbc0e90ae, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("equilibrium_process_estimation", method_pointer_b164d588de2f5515ae36b350e22788c6, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("gamma_estimation", method_pointer_b7f0f4bf43515f50b41c59907193dea1, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("zero_inflated_gamma_estimation", method_pointer_636052c7c86a5d409a449b97991cb254, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("inverse_gaussian_estimation", method_pointer_7a6825284cb7585d897e8ae2847483a9, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def_readwrite("nb_value", &::stat_tool::Reestimation< int >::nb_value, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def_readwrite("alloc_nb_value", &::stat_tool::Reestimation< int >::alloc_nb_value, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def_readwrite("offset", &::stat_tool::Reestimation< int >::offset, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def_readwrite("nb_element", &::stat_tool::Reestimation< int >::nb_element, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def_readwrite("max", &::stat_tool::Reestimation< int >::max, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def_readwrite("mean", &::stat_tool::Reestimation< int >::mean, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def_readwrite("variance", &::stat_tool::Reestimation< int >::variance, "");

}