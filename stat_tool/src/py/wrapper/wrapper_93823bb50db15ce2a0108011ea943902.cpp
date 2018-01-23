#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> struct ::std::integral_constant< bool, true > const volatile * get_pointer<struct ::std::integral_constant< bool, true > const volatile >(struct ::std::integral_constant< bool, true > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_93823bb50db15ce2a0108011ea943902()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::class_< struct ::std::integral_constant< bool, true >, autowig::Held< struct ::std::integral_constant< bool, true > >::Type > class_93823bb50db15ce2a0108011ea943902("_IntegralConstant_93823bb50db15ce2a0108011ea943902", "", boost::python::no_init);
    class_93823bb50db15ce2a0108011ea943902.def(boost::python::init<  >(""));
    class_93823bb50db15ce2a0108011ea943902.def(boost::python::init< struct ::std::integral_constant< bool, true > const & >(""));
    class_93823bb50db15ce2a0108011ea943902.def_readonly("value", ::std::integral_constant< bool, true >::value, "");

}