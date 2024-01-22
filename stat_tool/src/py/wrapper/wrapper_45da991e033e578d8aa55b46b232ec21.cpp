#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::locale::facet const volatile * get_pointer<class ::std::locale::facet const volatile >(class ::std::locale::facet const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_45da991e033e578d8aa55b46b232ec21()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    std::string name_2f3439617e035c41b1282a03e900ef19 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + "._locale");
    boost::python::object module_2f3439617e035c41b1282a03e900ef19(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_2f3439617e035c41b1282a03e900ef19.c_str()))));
    boost::python::scope().attr("_locale") = module_2f3439617e035c41b1282a03e900ef19;
    boost::python::scope scope_2f3439617e035c41b1282a03e900ef19 = module_2f3439617e035c41b1282a03e900ef19;
    boost::python::class_< class ::std::locale::facet, autowig::Held< class ::std::locale::facet >::Type, boost::noncopyable > class_45da991e033e578d8aa55b46b232ec21("Facet", "", boost::python::no_init);

}