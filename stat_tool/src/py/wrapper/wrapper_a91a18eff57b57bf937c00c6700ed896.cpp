#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> struct ::std::input_iterator_tag const volatile * get_pointer<struct ::std::input_iterator_tag const volatile >(struct ::std::input_iterator_tag const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_a91a18eff57b57bf937c00c6700ed896()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::class_< struct ::std::input_iterator_tag, autowig::Held< struct ::std::input_iterator_tag >::Type > class_a91a18eff57b57bf937c00c6700ed896("InputIteratorTag", "", boost::python::no_init);
    class_a91a18eff57b57bf937c00c6700ed896.def(boost::python::init<  >(""));
    class_a91a18eff57b57bf937c00c6700ed896.def(boost::python::init< struct ::std::input_iterator_tag const & >(""));

}