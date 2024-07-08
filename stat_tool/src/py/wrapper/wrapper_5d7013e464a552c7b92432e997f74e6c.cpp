#include "_stat_tool.h"



namespace autowig
{

    void method_decorator_af406258214b5e0d8168d234bb932acb(class ::stat_tool::MultiPlot & instance, int  param_in_0, const class ::stat_tool::SinglePlot & param_out) { instance.operator[](param_in_0) = param_out; }
}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::MultiPlot const volatile * get_pointer<class ::stat_tool::MultiPlot const volatile >(class ::stat_tool::MultiPlot const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_5d7013e464a552c7b92432e997f74e6c()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    class ::stat_tool::SinglePlot & (::stat_tool::MultiPlot::*method_pointer_af406258214b5e0d8168d234bb932acb)(int ) = &::stat_tool::MultiPlot::operator[];
    void  (::stat_tool::MultiPlot::*method_pointer_95aa6653df98514e8e5b7e3f8fbf1b4b)(int ) = &::stat_tool::MultiPlot::resize;
    int  (::stat_tool::MultiPlot::*method_pointer_eab7c2749bae5da3903dfe580bd77831)() = &::stat_tool::MultiPlot::size;
    boost::python::class_< class ::stat_tool::MultiPlot, autowig::Held< class ::stat_tool::MultiPlot >::Type > class_5d7013e464a552c7b92432e997f74e6c("MultiPlot", "", boost::python::no_init);
    class_5d7013e464a552c7b92432e997f74e6c.def(boost::python::init< int  >(""));
    class_5d7013e464a552c7b92432e997f74e6c.def("__getitem__", method_pointer_af406258214b5e0d8168d234bb932acb, boost::python::return_internal_reference<>(), "");
    class_5d7013e464a552c7b92432e997f74e6c.def("__getitem__", autowig::method_decorator_af406258214b5e0d8168d234bb932acb);
    class_5d7013e464a552c7b92432e997f74e6c.def("resize", method_pointer_95aa6653df98514e8e5b7e3f8fbf1b4b, "");
    class_5d7013e464a552c7b92432e997f74e6c.def("__len__", method_pointer_eab7c2749bae5da3903dfe580bd77831, "");
    class_5d7013e464a552c7b92432e997f74e6c.def_readwrite("title", &::stat_tool::MultiPlot::title, "");
    class_5d7013e464a552c7b92432e997f74e6c.def_readwrite("xtics", &::stat_tool::MultiPlot::xtics, "");
    class_5d7013e464a552c7b92432e997f74e6c.def_readwrite("ytics", &::stat_tool::MultiPlot::ytics, "");
    class_5d7013e464a552c7b92432e997f74e6c.def_readwrite("xrange", &::stat_tool::MultiPlot::xrange, "");
    class_5d7013e464a552c7b92432e997f74e6c.def_readwrite("yrange", &::stat_tool::MultiPlot::yrange, "");
    class_5d7013e464a552c7b92432e997f74e6c.def_readwrite("xlabel", &::stat_tool::MultiPlot::xlabel, "");
    class_5d7013e464a552c7b92432e997f74e6c.def_readwrite("ylabel", &::stat_tool::MultiPlot::ylabel, "");
    class_5d7013e464a552c7b92432e997f74e6c.def_readwrite("grid", &::stat_tool::MultiPlot::grid, "");
    class_5d7013e464a552c7b92432e997f74e6c.def_readwrite("group", &::stat_tool::MultiPlot::group, "");

}