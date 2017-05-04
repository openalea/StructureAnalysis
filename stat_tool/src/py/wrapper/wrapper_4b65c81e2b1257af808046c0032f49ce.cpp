#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::initializer_list< double > const volatile * get_pointer<class ::std::initializer_list< double > const volatile >(class ::std::initializer_list< double > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_4b65c81e2b1257af808046c0032f49ce()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    ::std::initializer_list< double >::size_type  (::std::initializer_list< double >::*method_pointer_e4601aa8fff55f44991a28baedce1b50)() const = &::std::initializer_list< double >::size;
    boost::python::class_< class ::std::initializer_list< double >, autowig::Held< class ::std::initializer_list< double > >::Type > class_4b65c81e2b1257af808046c0032f49ce("_InitializerList_4b65c81e2b1257af808046c0032f49ce", "", boost::python::no_init);
    class_4b65c81e2b1257af808046c0032f49ce.def("__len__", method_pointer_e4601aa8fff55f44991a28baedce1b50, "");

}