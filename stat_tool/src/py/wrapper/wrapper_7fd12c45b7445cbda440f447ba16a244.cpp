#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::Mixture const volatile * get_pointer<class ::stat_tool::Mixture const volatile >(class ::stat_tool::Mixture const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_7fd12c45b7445cbda440f447ba16a244()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    class ::stat_tool::DiscreteParametricModel * (::stat_tool::Mixture::*method_pointer_c8e5ae61c41751f6b367998e3a1c9f5b)(class ::stat_tool::StatError &, int , int ) const = &::stat_tool::Mixture::extract;
    class ::stat_tool::MixtureData * (::stat_tool::Mixture::*method_pointer_ef68cb2853d15909b42f4b8a9f14ec5a)(class ::stat_tool::StatError &) const = &::stat_tool::Mixture::extract_data;
    class ::stat_tool::Mixture * (::stat_tool::Mixture::*method_pointer_496e0d85f6695277bc677b0002f4ded1)(double ) const = &::stat_tool::Mixture::thresholding;
    class ::stat_tool::Mixture * (*method_pointer_bfda23f7a14a5062b8040d1e9e4b8128)(class ::stat_tool::StatError &, class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > const, double ) = ::stat_tool::Mixture::ascii_read;
    double  (::stat_tool::Mixture::*method_pointer_a9a77ac4837d5614af40e9748c503f80)(class ::stat_tool::MixtureData const &) const = &::stat_tool::Mixture::classification_likelihood_computation;
    double  (::stat_tool::Mixture::*method_pointer_73eb4eb85305564ea2fc9b91dc32941b)(class ::stat_tool::Vectors const &, int ) const = &::stat_tool::Mixture::likelihood_computation;
    class ::stat_tool::MixtureData * (::stat_tool::Mixture::*method_pointer_cac2788f9500518aa5edca6b080bd7f5)(class ::stat_tool::StatError &, int ) const = &::stat_tool::Mixture::simulation;
    class ::stat_tool::MixtureData * (::stat_tool::Mixture::*method_pointer_2ee28a917b47549abf14ad4d66d2435f)() const = &::stat_tool::Mixture::get_mixture_data;
    int  (::stat_tool::Mixture::*method_pointer_0cfa8cca28b151a29cb99a793d4598a1)() const = &::stat_tool::Mixture::get_nb_component;
    class ::stat_tool::DiscreteParametric * (::stat_tool::Mixture::*method_pointer_cceb5d6e7b9755c88bde614966c07952)() const = &::stat_tool::Mixture::get_weight;
    int  (::stat_tool::Mixture::*method_pointer_fd6e36af37df590491e51920f8a5e59c)() const = &::stat_tool::Mixture::get_nb_output_process;
    class ::stat_tool::CategoricalProcess * (::stat_tool::Mixture::*method_pointer_4a3affa993e55db88afe1cac8f8303d6)(int ) const = &::stat_tool::Mixture::get_categorical_process;
    class ::stat_tool::DiscreteParametricProcess * (::stat_tool::Mixture::*method_pointer_f6bef98d8a0652e2bcd0dca1597eb4b2)(int ) const = &::stat_tool::Mixture::get_discrete_parametric_process;
    class ::stat_tool::ContinuousParametricProcess * (::stat_tool::Mixture::*method_pointer_4b463eddb9cc5b5fb92002c7d34bffc5)(int ) const = &::stat_tool::Mixture::get_continuous_parametric_process;
    boost::python::class_< class ::stat_tool::Mixture, autowig::Held< class ::stat_tool::Mixture >::Type, boost::python::bases< class ::stat_tool::StatInterface > > class_7fd12c45b7445cbda440f447ba16a244("Mixture", "Multivariate mixture of distributions\n\n", boost::python::no_init);
    class_7fd12c45b7445cbda440f447ba16a244.def(boost::python::init<  >(""));
    class_7fd12c45b7445cbda440f447ba16a244.def(boost::python::init< int , double , double , double , bool  >(""));
    class_7fd12c45b7445cbda440f447ba16a244.def(boost::python::init< int , int , double , double , bool , enum ::stat_tool::tying_rule  >(""));
    class_7fd12c45b7445cbda440f447ba16a244.def(boost::python::init< class ::stat_tool::Mixture const &, bool  >(""));
    class_7fd12c45b7445cbda440f447ba16a244.def("extract", method_pointer_c8e5ae61c41751f6b367998e3a1c9f5b, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_7fd12c45b7445cbda440f447ba16a244.def("extract_data", method_pointer_ef68cb2853d15909b42f4b8a9f14ec5a, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_7fd12c45b7445cbda440f447ba16a244.def("thresholding", method_pointer_496e0d85f6695277bc677b0002f4ded1, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_7fd12c45b7445cbda440f447ba16a244.def("ascii_read", method_pointer_bfda23f7a14a5062b8040d1e9e4b8128, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_7fd12c45b7445cbda440f447ba16a244.def("classification_likelihood_computation", method_pointer_a9a77ac4837d5614af40e9748c503f80, "");
    class_7fd12c45b7445cbda440f447ba16a244.def("likelihood_computation", method_pointer_73eb4eb85305564ea2fc9b91dc32941b, "");
    class_7fd12c45b7445cbda440f447ba16a244.def("simulation", method_pointer_cac2788f9500518aa5edca6b080bd7f5, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_7fd12c45b7445cbda440f447ba16a244.def("get_mixture_data", method_pointer_2ee28a917b47549abf14ad4d66d2435f, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_7fd12c45b7445cbda440f447ba16a244.def("get_nb_component", method_pointer_0cfa8cca28b151a29cb99a793d4598a1, "");
    class_7fd12c45b7445cbda440f447ba16a244.def("get_weight", method_pointer_cceb5d6e7b9755c88bde614966c07952, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_7fd12c45b7445cbda440f447ba16a244.def("get_nb_output_process", method_pointer_fd6e36af37df590491e51920f8a5e59c, "");
    class_7fd12c45b7445cbda440f447ba16a244.def("get_categorical_process", method_pointer_4a3affa993e55db88afe1cac8f8303d6, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_7fd12c45b7445cbda440f447ba16a244.def("get_discrete_parametric_process", method_pointer_f6bef98d8a0652e2bcd0dca1597eb4b2, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_7fd12c45b7445cbda440f447ba16a244.def("get_continuous_parametric_process", method_pointer_4b463eddb9cc5b5fb92002c7d34bffc5, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_7fd12c45b7445cbda440f447ba16a244.staticmethod("ascii_read");

    if(autowig::Held< class ::stat_tool::Mixture >::is_class)
    {
        boost::python::implicitly_convertible< autowig::Held< class ::stat_tool::Mixture >::Type, autowig::Held< class ::stat_tool::StatInterface >::Type >();
    }

}