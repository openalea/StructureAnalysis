#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::Forward const volatile * get_pointer<class ::stat_tool::Forward const volatile >(class ::stat_tool::Forward const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_e4b407b97f515e8bbcdf5a6f0e8aa41c()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    void  (::stat_tool::Forward::*method_pointer_bfc33d5dcf2453368e8c72f45ac32592)(class ::stat_tool::DiscreteParametric const &) = &::stat_tool::Forward::computation;
    boost::python::class_< class ::stat_tool::Forward, autowig::Held< class ::stat_tool::Forward >::Type, boost::python::bases< class ::stat_tool::DiscreteParametric > > class_e4b407b97f515e8bbcdf5a6f0e8aa41c("Forward", "Forward recurrence or sojourn time distribution\n\n", boost::python::no_init);
    class_e4b407b97f515e8bbcdf5a6f0e8aa41c.def(boost::python::init< int , enum ::stat_tool::discrete_parametric , int , int , double , double  >(""));
    class_e4b407b97f515e8bbcdf5a6f0e8aa41c.def(boost::python::init< class ::stat_tool::DiscreteParametric const &, int  >(""));
    class_e4b407b97f515e8bbcdf5a6f0e8aa41c.def(boost::python::init< class ::stat_tool::Forward const &, int  >(""));
    class_e4b407b97f515e8bbcdf5a6f0e8aa41c.def("computation", method_pointer_bfc33d5dcf2453368e8c72f45ac32592, "");

    if(autowig::Held< class ::stat_tool::Forward >::is_class)
    {
        boost::python::implicitly_convertible< autowig::Held< class ::stat_tool::Forward >::Type, autowig::Held< class ::stat_tool::DiscreteParametric >::Type >();
    }

}