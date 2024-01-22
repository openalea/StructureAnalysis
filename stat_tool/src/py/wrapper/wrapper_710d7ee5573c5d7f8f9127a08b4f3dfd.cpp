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
    double  (::stat_tool::Reestimation< int >::*method_pointer_da2fe2021d4356b2a8f8b8d469c66d88)(class ::stat_tool::Distribution const &) const = &::stat_tool::Reestimation< int >::likelihood_computation;
    void  (::stat_tool::Reestimation< int >::*method_pointer_09f5f68263f3506d9d082e14dcc2a6e9)(int ) = &::stat_tool::Reestimation< int >::init;
    void  (::stat_tool::Reestimation< int >::*method_pointer_fae34c98fafd5884bb3fff56394048be)() = &::stat_tool::Reestimation< int >::offset_computation;
    void  (::stat_tool::Reestimation< int >::*method_pointer_f7e8187bd6f251f991a20a81eafc5da1)() = &::stat_tool::Reestimation< int >::max_computation;
    void  (::stat_tool::Reestimation< int >::*method_pointer_403634937e735ba29c79f285687fa465)() = &::stat_tool::Reestimation< int >::mean_computation;
    void  (::stat_tool::Reestimation< int >::*method_pointer_aec14b5616e452118d6e1cca66fa7784)(bool ) = &::stat_tool::Reestimation< int >::variance_computation;
    void  (::stat_tool::Reestimation< int >::*method_pointer_0e9d7d4cc2ae53f5bab80c00665ca2d4)(class ::stat_tool::Reestimation< int > const &) = &::stat_tool::Reestimation< int >::copy;
    boost::python::class_< class ::stat_tool::Reestimation< int >, autowig::Held< class ::stat_tool::Reestimation< int > >::Type > class_710d7ee5573c5d7f8f9127a08b4f3dfd("_Reestimation_710d7ee5573c5d7f8f9127a08b4f3dfd", "", boost::python::no_init);
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def(boost::python::init< int  >(""));
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def(boost::python::init< class ::stat_tool::Reestimation< int > const & >(""));
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("likelihood_computation", method_pointer_da2fe2021d4356b2a8f8b8d469c66d88, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("init", method_pointer_09f5f68263f3506d9d082e14dcc2a6e9, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("offset_computation", method_pointer_fae34c98fafd5884bb3fff56394048be, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("max_computation", method_pointer_f7e8187bd6f251f991a20a81eafc5da1, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("mean_computation", method_pointer_403634937e735ba29c79f285687fa465, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("variance_computation", method_pointer_aec14b5616e452118d6e1cca66fa7784, "");
    class_710d7ee5573c5d7f8f9127a08b4f3dfd.def("copy", method_pointer_0e9d7d4cc2ae53f5bab80c00665ca2d4, "");

}