#include "_stat_tool.h"



namespace autowig
{

    void method_decorator_4409fa6cc9db5f3fb94f9e8f5668d282(class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > > & instance, ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::size_type  param_in_0, const ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::reference  param_out) { instance.at(param_in_0) = param_out; }
    void method_decorator_9a3908795a5354a7a134f96cc3d4d40e(class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > > & instance, const ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::reference  param_out) { instance.front() = param_out; }
    void method_decorator_600215fa341d52bc8d5f3cee9b6e2d14(class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > > & instance, const ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::reference  param_out) { instance.back() = param_out; }
}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > > const volatile * get_pointer<class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > > const volatile >(class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_c5260d1a7cbb5ef3bc32b582df09a801()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    void  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_cbf4a8cdf73e546ba473bcd57daec19b)(::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::size_type , ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::const_reference ) = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::assign;
    ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::size_type  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_858e418f2f7951d2980b51c18987dd86)() const = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::size;
    ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::size_type  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_5833e4057c0b5e00b51873cb64c31736)() const = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::capacity;
    bool  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_8d18a26b56ad5ab888cbbb440a055a60)() const = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::empty;
    ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::size_type  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_e496883c9e3c539e9d6a09e64f63e56d)() const = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::max_size;
    void  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_176ce459777b586ba4a5f49851a0b4bd)(::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::size_type ) = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::reserve;
    ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::reference  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_4409fa6cc9db5f3fb94f9e8f5668d282)(::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::size_type ) = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::at;
    ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::const_reference  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_055d6a4e563653b3951db04999b25161)(::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::size_type ) const = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::at;
    ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::reference  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_9a3908795a5354a7a134f96cc3d4d40e)() = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::front;
    ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::const_reference  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_acf0fc6d8906510a8ed0c34da83e8071)() const = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::front;
    ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::reference  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_600215fa341d52bc8d5f3cee9b6e2d14)() = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::back;
    ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::const_reference  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_83166a375e8d58b29414279856a8ad40)() const = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::back;
    ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::value_type * (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_3cb2aec1a95655058ad1cd48bd0f8bc4)() = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::data;
    ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::value_type const * (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_33155e1a3cc356cc876abeead3163cfb)() const = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::data;
    void  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_8112bd9e318a54caaae1c775dc6f1324)(::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::const_reference ) = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::push_back;
    void  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_b9195bbce5e95056878cdad532101c58)() = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::pop_back;
    void  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_ac62c78786695f46b99b59a0b81e4ec9)() = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::clear;
    void  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_6b50cf3b00b5598ab313843dd6550b90)(class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > > &) = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::swap;
    bool  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_8e36a7a8b8b457ed925567ec57abaf7e)() const = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::__invariants;
    void  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_01fbf78b8de351a2a9056041f3da49b7)(class ::stat_tool::DiscreteParametric *, class ::stat_tool::DiscreteParametric *) = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::assign;
    void  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_4d6745df83155e5aa65ea12e418af070)(class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::DiscreteParametric * > > , class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::DiscreteParametric * > > ) = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::assign;
    void  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_da48e521c5df5024a6e1c48e9fe81ce8)(class ::stat_tool::DiscreteParametric const *, class ::stat_tool::DiscreteParametric const *) = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::assign;
    boost::python::class_< class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >, autowig::Held< class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > > >::Type > class_c5260d1a7cbb5ef3bc32b582df09a801("_Vector_c5260d1a7cbb5ef3bc32b582df09a801", "", boost::python::no_init);
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("assign", method_pointer_cbf4a8cdf73e546ba473bcd57daec19b, "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("__len__", method_pointer_858e418f2f7951d2980b51c18987dd86, "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("capacity", method_pointer_5833e4057c0b5e00b51873cb64c31736, "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("empty", method_pointer_8d18a26b56ad5ab888cbbb440a055a60, "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("max_size", method_pointer_e496883c9e3c539e9d6a09e64f63e56d, "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("reserve", method_pointer_176ce459777b586ba4a5f49851a0b4bd, "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("at", method_pointer_4409fa6cc9db5f3fb94f9e8f5668d282, boost::python::return_internal_reference<>(), "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("at", autowig::method_decorator_4409fa6cc9db5f3fb94f9e8f5668d282);
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("at", method_pointer_055d6a4e563653b3951db04999b25161, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("front", method_pointer_9a3908795a5354a7a134f96cc3d4d40e, boost::python::return_internal_reference<>(), "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("front", autowig::method_decorator_9a3908795a5354a7a134f96cc3d4d40e);
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("front", method_pointer_acf0fc6d8906510a8ed0c34da83e8071, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("back", method_pointer_600215fa341d52bc8d5f3cee9b6e2d14, boost::python::return_internal_reference<>(), "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("back", autowig::method_decorator_600215fa341d52bc8d5f3cee9b6e2d14);
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("back", method_pointer_83166a375e8d58b29414279856a8ad40, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("data", method_pointer_3cb2aec1a95655058ad1cd48bd0f8bc4, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("data", method_pointer_33155e1a3cc356cc876abeead3163cfb, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("push_back", method_pointer_8112bd9e318a54caaae1c775dc6f1324, "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("pop_back", method_pointer_b9195bbce5e95056878cdad532101c58, "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("clear", method_pointer_ac62c78786695f46b99b59a0b81e4ec9, "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("swap", method_pointer_6b50cf3b00b5598ab313843dd6550b90, "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("invariants", method_pointer_8e36a7a8b8b457ed925567ec57abaf7e, "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("assign", method_pointer_01fbf78b8de351a2a9056041f3da49b7, "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("assign", method_pointer_4d6745df83155e5aa65ea12e418af070, "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("assign", method_pointer_da48e521c5df5024a6e1c48e9fe81ce8, "");

    struct vector_c5260d1a7cbb5ef3bc32b582df09a801_from_python
    {
        vector_c5260d1a7cbb5ef3bc32b582df09a801_from_python()
        {
            boost::python::converter::registry::push_back(
                &convertible,
                &construct,
                boost::python::type_id< class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > > >());
        }

        static void* convertible(PyObject* obj_ptr)
        { return obj_ptr; }

        static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
        {
            boost::python::handle<> obj_iter(PyObject_GetIter(obj_ptr));
            void* storage = ((boost::python::converter::rvalue_from_python_storage< class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > > >*)data)->storage.bytes;
            new (storage) class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >();
            data->convertible = storage;
            class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >& result = *((class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >*)storage);
            unsigned int i = 0;
            for(;; i++)
            {
                boost::python::handle<> py_elem_hdl(boost::python::allow_null(PyIter_Next(obj_iter.get())));
                if(PyErr_Occurred())
                { boost::python::throw_error_already_set(); }
                if(!py_elem_hdl.get())
                { break; }
                boost::python::object py_elem_obj(py_elem_hdl);
                result.push_back(boost::python::extract< class ::stat_tool::DiscreteParametric  >(py_elem_obj));
            }
        }
    };

    vector_c5260d1a7cbb5ef3bc32b582df09a801_from_python();
}