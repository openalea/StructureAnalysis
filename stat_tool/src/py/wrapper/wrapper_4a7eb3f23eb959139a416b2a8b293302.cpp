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
    double  (::stat_tool::Reestimation< double >::*method_pointer_635bb10b5e43521c8c74f55556422bfb)(class ::stat_tool::Distribution const &) const = &::stat_tool::Reestimation< double >::likelihood_computation;
    boost::python::class_< class ::stat_tool::Reestimation< double >, autowig::Held< class ::stat_tool::Reestimation< double > >::Type > class_4a7eb3f23eb959139a416b2a8b293302("_Reestimation_4a7eb3f23eb959139a416b2a8b293302", "", boost::python::no_init);
    class_4a7eb3f23eb959139a416b2a8b293302.def("likelihood_computation", method_pointer_635bb10b5e43521c8c74f55556422bfb, "");

}