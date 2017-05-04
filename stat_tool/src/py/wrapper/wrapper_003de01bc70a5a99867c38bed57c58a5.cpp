#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::DiscreteParametricProcess const volatile * get_pointer<class ::stat_tool::DiscreteParametricProcess const volatile >(class ::stat_tool::DiscreteParametricProcess const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_003de01bc70a5a99867c38bed57c58a5()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    void  (::stat_tool::DiscreteParametricProcess::*method_pointer_0659c17596ec518882570be258267c93)(class ::stat_tool::DiscreteParametricProcess const &) = &::stat_tool::DiscreteParametricProcess::copy;
    void  (::stat_tool::DiscreteParametricProcess::*method_pointer_850a40312e435b39ad131196d1a87f09)() = &::stat_tool::DiscreteParametricProcess::remove;
    void  (::stat_tool::DiscreteParametricProcess::*method_pointer_1f5be9cbc73c53a5ad2264314197ec16)() = &::stat_tool::DiscreteParametricProcess::nb_value_computation;
    int  (::stat_tool::DiscreteParametricProcess::*method_pointer_bf5091e7332b5c1f80763c568d769400)() const = &::stat_tool::DiscreteParametricProcess::nb_parameter_computation;
    double  (::stat_tool::DiscreteParametricProcess::*method_pointer_464f8d86b1d95924b46b13dd966b6773)(class ::stat_tool::Distribution *) const = &::stat_tool::DiscreteParametricProcess::mean_computation;
    double  (::stat_tool::DiscreteParametricProcess::*method_pointer_2fd07342c4735ca5ae18a3164e60c237)(class ::stat_tool::Distribution *, double ) const = &::stat_tool::DiscreteParametricProcess::variance_computation;
    class ::stat_tool::Distribution * (::stat_tool::DiscreteParametricProcess::*method_pointer_0a3b99321bbe5f59907fd46f2775d1c3)(class ::stat_tool::Distribution *) = &::stat_tool::DiscreteParametricProcess::mixture_computation;
    void  (::stat_tool::DiscreteParametricProcess::*method_pointer_ef05876447e35b788dd014ff886c147b)() = &::stat_tool::DiscreteParametricProcess::init;
    boost::python::class_< class ::stat_tool::DiscreteParametricProcess, autowig::Held< class ::stat_tool::DiscreteParametricProcess >::Type > class_003de01bc70a5a99867c38bed57c58a5("DiscreteParametricProcess", "Discrete parametric observation process\n\n", boost::python::no_init);
    class_003de01bc70a5a99867c38bed57c58a5.def(boost::python::init< int , int  >(""));
    class_003de01bc70a5a99867c38bed57c58a5.def(boost::python::init< class ::stat_tool::DiscreteParametricProcess const & >(""));
    class_003de01bc70a5a99867c38bed57c58a5.def("copy", method_pointer_0659c17596ec518882570be258267c93, "");
    class_003de01bc70a5a99867c38bed57c58a5.def("remove", method_pointer_850a40312e435b39ad131196d1a87f09, "");
    class_003de01bc70a5a99867c38bed57c58a5.def("nb_value_computation", method_pointer_1f5be9cbc73c53a5ad2264314197ec16, "");
    class_003de01bc70a5a99867c38bed57c58a5.def("nb_parameter_computation", method_pointer_bf5091e7332b5c1f80763c568d769400, "");
    class_003de01bc70a5a99867c38bed57c58a5.def("mean_computation", method_pointer_464f8d86b1d95924b46b13dd966b6773, "");
    class_003de01bc70a5a99867c38bed57c58a5.def("variance_computation", method_pointer_2fd07342c4735ca5ae18a3164e60c237, "");
    class_003de01bc70a5a99867c38bed57c58a5.def("mixture_computation", method_pointer_0a3b99321bbe5f59907fd46f2775d1c3, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_003de01bc70a5a99867c38bed57c58a5.def("init", method_pointer_ef05876447e35b788dd014ff886c147b, "");
    class_003de01bc70a5a99867c38bed57c58a5.def_readwrite("nb_state", &::stat_tool::DiscreteParametricProcess::nb_state, "");
    class_003de01bc70a5a99867c38bed57c58a5.def_readwrite("nb_value", &::stat_tool::DiscreteParametricProcess::nb_value, "");

}