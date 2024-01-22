#include "_stat_tool.h"



namespace autowig
{

    void method_decorator_94cd69d22a215710bbc9d5d47d9ef407(class ::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > > & instance, ::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::size_type  param_in_0, const ::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::reference  param_out) { instance.at(param_in_0) = param_out; }
    void method_decorator_9fd4c67b61e95859bd4ff0d340f4e561(class ::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > > & instance, const ::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::reference  param_out) { instance.front() = param_out; }
    void method_decorator_6529949466095ba4b601d75b91975843(class ::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > > & instance, const ::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::reference  param_out) { instance.back() = param_out; }
}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > > const volatile * get_pointer<class ::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > > const volatile >(class ::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_52efeeab642a50c286578b88496f1362()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    void  (::std::vector< ::std::vector< int, ::std::allocator< int > >, ::std::allocator< ::std::vector< int, ::std::allocator< int > > > >::*method_pointer_e3635631042c5f3093b45f077afd13f5)(::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::size_type , ::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::const_reference ) = &::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::assign;
    ::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::size_type  (::std::vector< ::std::vector< int, ::std::allocator< int > >, ::std::allocator< ::std::vector< int, ::std::allocator< int > > > >::*method_pointer_8871c22805c3594b8073fb7381c3fd23)() const = &::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::size;
    ::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::size_type  (::std::vector< ::std::vector< int, ::std::allocator< int > >, ::std::allocator< ::std::vector< int, ::std::allocator< int > > > >::*method_pointer_523864aa8a165ceeb15b2517c80ab992)() const = &::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::capacity;
    bool  (::std::vector< ::std::vector< int, ::std::allocator< int > >, ::std::allocator< ::std::vector< int, ::std::allocator< int > > > >::*method_pointer_11d324798c245f008e2ee109a78bc5d7)() const = &::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::empty;
    ::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::size_type  (::std::vector< ::std::vector< int, ::std::allocator< int > >, ::std::allocator< ::std::vector< int, ::std::allocator< int > > > >::*method_pointer_2a0f4e69cdce5c4182d8067495d859cf)() const = &::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::max_size;
    void  (::std::vector< ::std::vector< int, ::std::allocator< int > >, ::std::allocator< ::std::vector< int, ::std::allocator< int > > > >::*method_pointer_a09d8c9dff1956cc897b82aa563569cc)(::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::size_type ) = &::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::reserve;
    ::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::reference  (::std::vector< ::std::vector< int, ::std::allocator< int > >, ::std::allocator< ::std::vector< int, ::std::allocator< int > > > >::*method_pointer_94cd69d22a215710bbc9d5d47d9ef407)(::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::size_type ) = &::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::at;
    ::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::const_reference  (::std::vector< ::std::vector< int, ::std::allocator< int > >, ::std::allocator< ::std::vector< int, ::std::allocator< int > > > >::*method_pointer_6c3a26f257ab537ea221dde6ef510e2e)(::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::size_type ) const = &::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::at;
    ::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::reference  (::std::vector< ::std::vector< int, ::std::allocator< int > >, ::std::allocator< ::std::vector< int, ::std::allocator< int > > > >::*method_pointer_9fd4c67b61e95859bd4ff0d340f4e561)() = &::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::front;
    ::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::const_reference  (::std::vector< ::std::vector< int, ::std::allocator< int > >, ::std::allocator< ::std::vector< int, ::std::allocator< int > > > >::*method_pointer_773d6b903bbc5b4f8ec6bb8f0d20a934)() const = &::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::front;
    ::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::reference  (::std::vector< ::std::vector< int, ::std::allocator< int > >, ::std::allocator< ::std::vector< int, ::std::allocator< int > > > >::*method_pointer_6529949466095ba4b601d75b91975843)() = &::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::back;
    ::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::const_reference  (::std::vector< ::std::vector< int, ::std::allocator< int > >, ::std::allocator< ::std::vector< int, ::std::allocator< int > > > >::*method_pointer_0cebfc5bb1d55c16970a0dc48da20824)() const = &::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::back;
    ::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::value_type * (::std::vector< ::std::vector< int, ::std::allocator< int > >, ::std::allocator< ::std::vector< int, ::std::allocator< int > > > >::*method_pointer_b7a50f4c0c205e7080f3f325c35c4aa1)() = &::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::data;
    ::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::value_type const * (::std::vector< ::std::vector< int, ::std::allocator< int > >, ::std::allocator< ::std::vector< int, ::std::allocator< int > > > >::*method_pointer_05eb8c2eb7935cabae8ee45bb2373556)() const = &::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::data;
    void  (::std::vector< ::std::vector< int, ::std::allocator< int > >, ::std::allocator< ::std::vector< int, ::std::allocator< int > > > >::*method_pointer_71c8b779e5ad54efa6625312113a5a60)(::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::const_reference ) = &::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::push_back;
    void  (::std::vector< ::std::vector< int, ::std::allocator< int > >, ::std::allocator< ::std::vector< int, ::std::allocator< int > > > >::*method_pointer_f6a58a17fb735efe9168d4eba95615b2)() = &::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::pop_back;
    void  (::std::vector< ::std::vector< int, ::std::allocator< int > >, ::std::allocator< ::std::vector< int, ::std::allocator< int > > > >::*method_pointer_3a71f1b847b859b7828e1165e97db42f)() = &::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::clear;
    void  (::std::vector< ::std::vector< int, ::std::allocator< int > >, ::std::allocator< ::std::vector< int, ::std::allocator< int > > > >::*method_pointer_c5e4c77a0a7c598f935c54f4c17c8fe4)(class ::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > > &) = &::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::swap;
    bool  (::std::vector< ::std::vector< int, ::std::allocator< int > >, ::std::allocator< ::std::vector< int, ::std::allocator< int > > > >::*method_pointer_5583f22d57975390ba49bfd1dcda4e00)() const = &::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::__invariants;
    void  (::std::vector< ::std::vector< int, ::std::allocator< int > >, ::std::allocator< ::std::vector< int, ::std::allocator< int > > > >::*method_pointer_9370d8cc1ac85bbd945ca6bcd7872e66)(class ::std::vector< int, class ::std::allocator< int > > *, class ::std::vector< int, class ::std::allocator< int > > *) = &::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::assign;
    void  (::std::vector< ::std::vector< int, ::std::allocator< int > >, ::std::allocator< ::std::vector< int, ::std::allocator< int > > > >::*method_pointer_cda97ac3dbe957e3afa498c964c9916f)(class ::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< int, class ::std::allocator< int > > * > > , class ::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< int, class ::std::allocator< int > > * > > ) = &::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::assign;
    void  (::std::vector< ::std::vector< int, ::std::allocator< int > >, ::std::allocator< ::std::vector< int, ::std::allocator< int > > > >::*method_pointer_526ebbafea3559bea02845838a7136bb)(class ::std::vector< int, class ::std::allocator< int > > const *, class ::std::vector< int, class ::std::allocator< int > > const *) = &::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >::assign;
    boost::python::class_< class ::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >, autowig::Held< class ::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > > >::Type > class_52efeeab642a50c286578b88496f1362("_Vector_52efeeab642a50c286578b88496f1362", "", boost::python::no_init);
    class_52efeeab642a50c286578b88496f1362.def("assign", method_pointer_e3635631042c5f3093b45f077afd13f5, "");
    class_52efeeab642a50c286578b88496f1362.def("__len__", method_pointer_8871c22805c3594b8073fb7381c3fd23, "");
    class_52efeeab642a50c286578b88496f1362.def("capacity", method_pointer_523864aa8a165ceeb15b2517c80ab992, "");
    class_52efeeab642a50c286578b88496f1362.def("empty", method_pointer_11d324798c245f008e2ee109a78bc5d7, "");
    class_52efeeab642a50c286578b88496f1362.def("max_size", method_pointer_2a0f4e69cdce5c4182d8067495d859cf, "");
    class_52efeeab642a50c286578b88496f1362.def("reserve", method_pointer_a09d8c9dff1956cc897b82aa563569cc, "");
    class_52efeeab642a50c286578b88496f1362.def("at", method_pointer_94cd69d22a215710bbc9d5d47d9ef407, boost::python::return_internal_reference<>(), "");
    class_52efeeab642a50c286578b88496f1362.def("at", autowig::method_decorator_94cd69d22a215710bbc9d5d47d9ef407);
    class_52efeeab642a50c286578b88496f1362.def("at", method_pointer_6c3a26f257ab537ea221dde6ef510e2e, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_52efeeab642a50c286578b88496f1362.def("front", method_pointer_9fd4c67b61e95859bd4ff0d340f4e561, boost::python::return_internal_reference<>(), "");
    class_52efeeab642a50c286578b88496f1362.def("front", autowig::method_decorator_9fd4c67b61e95859bd4ff0d340f4e561);
    class_52efeeab642a50c286578b88496f1362.def("front", method_pointer_773d6b903bbc5b4f8ec6bb8f0d20a934, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_52efeeab642a50c286578b88496f1362.def("back", method_pointer_6529949466095ba4b601d75b91975843, boost::python::return_internal_reference<>(), "");
    class_52efeeab642a50c286578b88496f1362.def("back", autowig::method_decorator_6529949466095ba4b601d75b91975843);
    class_52efeeab642a50c286578b88496f1362.def("back", method_pointer_0cebfc5bb1d55c16970a0dc48da20824, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_52efeeab642a50c286578b88496f1362.def("data", method_pointer_b7a50f4c0c205e7080f3f325c35c4aa1, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_52efeeab642a50c286578b88496f1362.def("data", method_pointer_05eb8c2eb7935cabae8ee45bb2373556, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_52efeeab642a50c286578b88496f1362.def("push_back", method_pointer_71c8b779e5ad54efa6625312113a5a60, "");
    class_52efeeab642a50c286578b88496f1362.def("pop_back", method_pointer_f6a58a17fb735efe9168d4eba95615b2, "");
    class_52efeeab642a50c286578b88496f1362.def("clear", method_pointer_3a71f1b847b859b7828e1165e97db42f, "");
    class_52efeeab642a50c286578b88496f1362.def("swap", method_pointer_c5e4c77a0a7c598f935c54f4c17c8fe4, "");
    class_52efeeab642a50c286578b88496f1362.def("invariants", method_pointer_5583f22d57975390ba49bfd1dcda4e00, "");
    class_52efeeab642a50c286578b88496f1362.def("assign", method_pointer_9370d8cc1ac85bbd945ca6bcd7872e66, "");
    class_52efeeab642a50c286578b88496f1362.def("assign", method_pointer_cda97ac3dbe957e3afa498c964c9916f, "");
    class_52efeeab642a50c286578b88496f1362.def("assign", method_pointer_526ebbafea3559bea02845838a7136bb, "");

    struct vector_52efeeab642a50c286578b88496f1362_from_python
    {
        vector_52efeeab642a50c286578b88496f1362_from_python()
        {
            boost::python::converter::registry::push_back(
                &convertible,
                &construct,
                boost::python::type_id< class ::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > > >());
        }

        static void* convertible(PyObject* obj_ptr)
        { return obj_ptr; }

        static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
        {
            boost::python::handle<> obj_iter(PyObject_GetIter(obj_ptr));
            void* storage = ((boost::python::converter::rvalue_from_python_storage< class ::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > > >*)data)->storage.bytes;
            new (storage) class ::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >();
            data->convertible = storage;
            class ::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >& result = *((class ::std::vector< class ::std::vector< int, class ::std::allocator< int > >, class ::std::allocator< class ::std::vector< int, class ::std::allocator< int > > > >*)storage);
            unsigned int i = 0;
            for(;; i++)
            {
                boost::python::handle<> py_elem_hdl(boost::python::allow_null(PyIter_Next(obj_iter.get())));
                if(PyErr_Occurred())
                { boost::python::throw_error_already_set(); }
                if(!py_elem_hdl.get())
                { break; }
                boost::python::object py_elem_obj(py_elem_hdl);
                result.push_back(boost::python::extract< class ::std::vector< int, class ::std::allocator< int > >  >(py_elem_obj));
            }
        }
    };

    vector_52efeeab642a50c286578b88496f1362_from_python();
}