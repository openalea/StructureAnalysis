#include "_stat_tool.h"



namespace autowig
{

    void method_decorator_7365756e35f3565592c28fbb79abbc51(class ::stat_tool::TemplateMultiPlotSet< enum ::stat_tool::process_distribution > & instance, int  param_in_0, const class ::stat_tool::MultiPlot & param_out) { instance.operator[](param_in_0) = param_out; }
}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::TemplateMultiPlotSet< enum ::stat_tool::process_distribution > const volatile * get_pointer<class ::stat_tool::TemplateMultiPlotSet< enum ::stat_tool::process_distribution > const volatile >(class ::stat_tool::TemplateMultiPlotSet< enum ::stat_tool::process_distribution > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_4f445857b49f50029570f15c5ecf9324()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    class ::stat_tool::MultiPlot & (::stat_tool::TemplateMultiPlotSet< enum ::stat_tool::process_distribution >::*method_pointer_7365756e35f3565592c28fbb79abbc51)(int ) = &::stat_tool::TemplateMultiPlotSet< enum ::stat_tool::process_distribution >::operator[];
    int  (::stat_tool::TemplateMultiPlotSet< enum ::stat_tool::process_distribution >::*method_pointer_ab6753b51ae25240aa7951c0d1db07d8)() = &::stat_tool::TemplateMultiPlotSet< enum ::stat_tool::process_distribution >::size;
    boost::python::class_< class ::stat_tool::TemplateMultiPlotSet< enum ::stat_tool::process_distribution >, autowig::Held< class ::stat_tool::TemplateMultiPlotSet< enum ::stat_tool::process_distribution > >::Type > class_4f445857b49f50029570f15c5ecf9324("_TemplateMultiPlotSet_4f445857b49f50029570f15c5ecf9324", "", boost::python::no_init);
    class_4f445857b49f50029570f15c5ecf9324.def("__getitem__", method_pointer_7365756e35f3565592c28fbb79abbc51, boost::python::return_internal_reference<>(), "");
    class_4f445857b49f50029570f15c5ecf9324.def("__getitem__", autowig::method_decorator_7365756e35f3565592c28fbb79abbc51);
    class_4f445857b49f50029570f15c5ecf9324.def("__len__", method_pointer_ab6753b51ae25240aa7951c0d1db07d8, "");

}