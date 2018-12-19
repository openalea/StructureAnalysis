#include "_stat_tool.h"



namespace autowig
{

    void method_decorator_76e1e28f2dcd5a39991614335fac9e25(class ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > > & instance, ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::size_type  param_in_0, const ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::reference  param_out) { instance.at(param_in_0) = param_out; }
    void method_decorator_132e20a77b3255c48587a899b67b8872(class ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > > & instance, const ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::reference  param_out) { instance.front() = param_out; }
    void method_decorator_08e5ff75ace55fab8537144a41fe71d7(class ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > > & instance, const ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::reference  param_out) { instance.back() = param_out; }
}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > > const volatile * get_pointer<class ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > > const volatile >(class ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_f0840635c9e2594cb061b2ce1bc3514a()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    void  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_7cc449c669f75214a6afb35c69cdffef)(::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::size_type , ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::const_reference ) = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::assign;
    ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::size_type  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_c67ba0db26e351fe86cb9916cbf32950)() const = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::size;
    ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::size_type  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_ace27348c09f580a9cdd0f580ee38597)() const = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::capacity;
    bool  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_58677708c8cf52ddb7d5489daa8f5a5c)() const = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::empty;
    ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::size_type  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_eaba44dc1e065f75898d5488c271c7d8)() const = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::max_size;
    void  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_39d2842b2f745185acffea75cc7e9f2a)(::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::size_type ) = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::reserve;
    ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::reference  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_76e1e28f2dcd5a39991614335fac9e25)(::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::size_type ) = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::at;
    ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::const_reference  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_7a76521fedf65cd7acbb0e6d4bd73b51)(::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::size_type ) const = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::at;
    ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::reference  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_132e20a77b3255c48587a899b67b8872)() = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::front;
    ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::const_reference  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_37db95842c715205a0eb68138abfd925)() const = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::front;
    ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::reference  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_08e5ff75ace55fab8537144a41fe71d7)() = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::back;
    ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::const_reference  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_c8631a20cd355aafbe65477431797ecc)() const = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::back;
    ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::value_type * (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_55c22c66377c5476ba312498f21898ef)() = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::data;
    ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::value_type const * (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_96fe6d4c32f95e03a47c09ea5ef59f88)() const = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::data;
    void  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_7e62a1f8caae5074b47acf0ab87114b7)(::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::const_reference ) = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::push_back;
    void  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_880c093d3a54565cab3ef0feb35c4fc8)() = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::pop_back;
    void  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_bf9613453d42516693d7e8c99c31a881)() = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::clear;
    void  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_f3ec1ac6858d538dab23a22ff57a9afc)(class ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > > &) = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::swap;
    bool  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_dddfa9916cbd5e8aaf2a6dd103cfafdd)() const = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::__invariants;
    void  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_d144e11524295ac9bcf15f7be190b347)(class ::stat_tool::FrequencyDistribution *, class ::stat_tool::FrequencyDistribution *) = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::assign;
    void  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_599948a55230599697ac719eadeb6137)(class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::FrequencyDistribution * > > , class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::FrequencyDistribution * > > ) = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::assign;
    void  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_0ad2af814e325367bb15f54f999ee32c)(class ::stat_tool::FrequencyDistribution const *, class ::stat_tool::FrequencyDistribution const *) = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::assign;
    boost::python::class_< class ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >, autowig::Held< class ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > > >::Type > class_f0840635c9e2594cb061b2ce1bc3514a("_Vector_f0840635c9e2594cb061b2ce1bc3514a", "", boost::python::no_init);
    class_f0840635c9e2594cb061b2ce1bc3514a.def("assign", method_pointer_7cc449c669f75214a6afb35c69cdffef, "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("__len__", method_pointer_c67ba0db26e351fe86cb9916cbf32950, "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("capacity", method_pointer_ace27348c09f580a9cdd0f580ee38597, "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("empty", method_pointer_58677708c8cf52ddb7d5489daa8f5a5c, "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("max_size", method_pointer_eaba44dc1e065f75898d5488c271c7d8, "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("reserve", method_pointer_39d2842b2f745185acffea75cc7e9f2a, "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("at", method_pointer_76e1e28f2dcd5a39991614335fac9e25, boost::python::return_internal_reference<>(), "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("at", autowig::method_decorator_76e1e28f2dcd5a39991614335fac9e25);
    class_f0840635c9e2594cb061b2ce1bc3514a.def("at", method_pointer_7a76521fedf65cd7acbb0e6d4bd73b51, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("front", method_pointer_132e20a77b3255c48587a899b67b8872, boost::python::return_internal_reference<>(), "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("front", autowig::method_decorator_132e20a77b3255c48587a899b67b8872);
    class_f0840635c9e2594cb061b2ce1bc3514a.def("front", method_pointer_37db95842c715205a0eb68138abfd925, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("back", method_pointer_08e5ff75ace55fab8537144a41fe71d7, boost::python::return_internal_reference<>(), "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("back", autowig::method_decorator_08e5ff75ace55fab8537144a41fe71d7);
    class_f0840635c9e2594cb061b2ce1bc3514a.def("back", method_pointer_c8631a20cd355aafbe65477431797ecc, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("data", method_pointer_55c22c66377c5476ba312498f21898ef, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("data", method_pointer_96fe6d4c32f95e03a47c09ea5ef59f88, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("push_back", method_pointer_7e62a1f8caae5074b47acf0ab87114b7, "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("pop_back", method_pointer_880c093d3a54565cab3ef0feb35c4fc8, "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("clear", method_pointer_bf9613453d42516693d7e8c99c31a881, "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("swap", method_pointer_f3ec1ac6858d538dab23a22ff57a9afc, "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("invariants", method_pointer_dddfa9916cbd5e8aaf2a6dd103cfafdd, "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("assign", method_pointer_d144e11524295ac9bcf15f7be190b347, "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("assign", method_pointer_599948a55230599697ac719eadeb6137, "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("assign", method_pointer_0ad2af814e325367bb15f54f999ee32c, "");

    struct vector_f0840635c9e2594cb061b2ce1bc3514a_from_python
    {
        vector_f0840635c9e2594cb061b2ce1bc3514a_from_python()
        {
            boost::python::converter::registry::push_back(
                &convertible,
                &construct,
                boost::python::type_id< class ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > > >());
        }

        static void* convertible(PyObject* obj_ptr)
        { return obj_ptr; }

        static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
        {
            boost::python::handle<> obj_iter(PyObject_GetIter(obj_ptr));
            void* storage = ((boost::python::converter::rvalue_from_python_storage< class ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > > >*)data)->storage.bytes;
            new (storage) class ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >();
            data->convertible = storage;
            class ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >& result = *((class ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >*)storage);
            unsigned int i = 0;
            for(;; i++)
            {
                boost::python::handle<> py_elem_hdl(boost::python::allow_null(PyIter_Next(obj_iter.get())));
                if(PyErr_Occurred())
                { boost::python::throw_error_already_set(); }
                if(!py_elem_hdl.get())
                { break; }
                boost::python::object py_elem_obj(py_elem_hdl);
                result.push_back(boost::python::extract< class ::stat_tool::FrequencyDistribution  >(py_elem_obj));
            }
        }
    };

    vector_f0840635c9e2594cb061b2ce1bc3514a_from_python();
}