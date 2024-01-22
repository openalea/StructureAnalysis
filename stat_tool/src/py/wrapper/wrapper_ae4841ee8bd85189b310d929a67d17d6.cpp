#include "_stat_tool.h"



namespace autowig
{

    void method_decorator_6495e8031e6a586397e5a5b434c10381(class ::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > > & instance, ::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::size_type  param_in_0, const ::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::reference  param_out) { instance.at(param_in_0) = param_out; }
    void method_decorator_00b8f41fa0f45c55bd15b1778a4cd00e(class ::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > > & instance, const ::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::reference  param_out) { instance.front() = param_out; }
    void method_decorator_478e344e518b5381a3c45e535f065b59(class ::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > > & instance, const ::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::reference  param_out) { instance.back() = param_out; }
}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > > const volatile * get_pointer<class ::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > > const volatile >(class ::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_ae4841ee8bd85189b310d929a67d17d6()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    void  (::std::vector< ::std::vector< double, ::std::allocator< double > >, ::std::allocator< ::std::vector< double, ::std::allocator< double > > > >::*method_pointer_5219a7f9dd6e5bc182199bcf34e23a1b)(::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::size_type , ::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::const_reference ) = &::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::assign;
    ::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::size_type  (::std::vector< ::std::vector< double, ::std::allocator< double > >, ::std::allocator< ::std::vector< double, ::std::allocator< double > > > >::*method_pointer_e82624f3a0c850ea86b6188e72924462)() const = &::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::size;
    ::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::size_type  (::std::vector< ::std::vector< double, ::std::allocator< double > >, ::std::allocator< ::std::vector< double, ::std::allocator< double > > > >::*method_pointer_81bb1441552c50c79761ba6be07f5bfa)() const = &::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::capacity;
    bool  (::std::vector< ::std::vector< double, ::std::allocator< double > >, ::std::allocator< ::std::vector< double, ::std::allocator< double > > > >::*method_pointer_22a8bc2066a05cc096f4b0784d96a924)() const = &::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::empty;
    ::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::size_type  (::std::vector< ::std::vector< double, ::std::allocator< double > >, ::std::allocator< ::std::vector< double, ::std::allocator< double > > > >::*method_pointer_e227b09d898750208e1a33e79fe5c359)() const = &::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::max_size;
    void  (::std::vector< ::std::vector< double, ::std::allocator< double > >, ::std::allocator< ::std::vector< double, ::std::allocator< double > > > >::*method_pointer_83ef44bec57455448b7e487a168437f8)(::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::size_type ) = &::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::reserve;
    ::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::reference  (::std::vector< ::std::vector< double, ::std::allocator< double > >, ::std::allocator< ::std::vector< double, ::std::allocator< double > > > >::*method_pointer_6495e8031e6a586397e5a5b434c10381)(::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::size_type ) = &::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::at;
    ::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::const_reference  (::std::vector< ::std::vector< double, ::std::allocator< double > >, ::std::allocator< ::std::vector< double, ::std::allocator< double > > > >::*method_pointer_47fe2297249154a5a0b6bdf86d6fc848)(::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::size_type ) const = &::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::at;
    ::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::reference  (::std::vector< ::std::vector< double, ::std::allocator< double > >, ::std::allocator< ::std::vector< double, ::std::allocator< double > > > >::*method_pointer_00b8f41fa0f45c55bd15b1778a4cd00e)() = &::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::front;
    ::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::const_reference  (::std::vector< ::std::vector< double, ::std::allocator< double > >, ::std::allocator< ::std::vector< double, ::std::allocator< double > > > >::*method_pointer_025c4ba3f3fa53bab0dbcada0ca359f1)() const = &::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::front;
    ::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::reference  (::std::vector< ::std::vector< double, ::std::allocator< double > >, ::std::allocator< ::std::vector< double, ::std::allocator< double > > > >::*method_pointer_478e344e518b5381a3c45e535f065b59)() = &::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::back;
    ::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::const_reference  (::std::vector< ::std::vector< double, ::std::allocator< double > >, ::std::allocator< ::std::vector< double, ::std::allocator< double > > > >::*method_pointer_8402a2ee7ae85f688e50b48a8d2633cd)() const = &::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::back;
    ::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::value_type * (::std::vector< ::std::vector< double, ::std::allocator< double > >, ::std::allocator< ::std::vector< double, ::std::allocator< double > > > >::*method_pointer_0e308e3b5cee5630a41cdd5760b3ea50)() = &::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::data;
    ::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::value_type const * (::std::vector< ::std::vector< double, ::std::allocator< double > >, ::std::allocator< ::std::vector< double, ::std::allocator< double > > > >::*method_pointer_5ea394f283175a28b533129752257a3e)() const = &::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::data;
    void  (::std::vector< ::std::vector< double, ::std::allocator< double > >, ::std::allocator< ::std::vector< double, ::std::allocator< double > > > >::*method_pointer_3c2aed4cc2e7526eae27d6fcf142e32e)(::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::const_reference ) = &::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::push_back;
    void  (::std::vector< ::std::vector< double, ::std::allocator< double > >, ::std::allocator< ::std::vector< double, ::std::allocator< double > > > >::*method_pointer_bf0adbdf1bb8501e842960a0d96161f2)() = &::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::pop_back;
    void  (::std::vector< ::std::vector< double, ::std::allocator< double > >, ::std::allocator< ::std::vector< double, ::std::allocator< double > > > >::*method_pointer_d491e2e8cad55edc816efcce40dfdaba)() = &::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::clear;
    void  (::std::vector< ::std::vector< double, ::std::allocator< double > >, ::std::allocator< ::std::vector< double, ::std::allocator< double > > > >::*method_pointer_a6d568faf5c45dc4ac25ca463e37e9f3)(class ::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > > &) = &::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::swap;
    bool  (::std::vector< ::std::vector< double, ::std::allocator< double > >, ::std::allocator< ::std::vector< double, ::std::allocator< double > > > >::*method_pointer_ff6af1af31db560babcbdd102f79a354)() const = &::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::__invariants;
    void  (::std::vector< ::std::vector< double, ::std::allocator< double > >, ::std::allocator< ::std::vector< double, ::std::allocator< double > > > >::*method_pointer_539a166ed5dd55639d7eb9eaeec18ff2)(class ::std::vector< double, class ::std::allocator< double > > *, class ::std::vector< double, class ::std::allocator< double > > *) = &::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::assign;
    void  (::std::vector< ::std::vector< double, ::std::allocator< double > >, ::std::allocator< ::std::vector< double, ::std::allocator< double > > > >::*method_pointer_c1a9da50f7405c088097559461c1a892)(class ::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< double, class ::std::allocator< double > > * > > , class ::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< double, class ::std::allocator< double > > * > > ) = &::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::assign;
    void  (::std::vector< ::std::vector< double, ::std::allocator< double > >, ::std::allocator< ::std::vector< double, ::std::allocator< double > > > >::*method_pointer_0c396dc8a69d519995ee0ebe3d23830b)(class ::std::vector< double, class ::std::allocator< double > > const *, class ::std::vector< double, class ::std::allocator< double > > const *) = &::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >::assign;
    boost::python::class_< class ::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >, autowig::Held< class ::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > > >::Type > class_ae4841ee8bd85189b310d929a67d17d6("_Vector_ae4841ee8bd85189b310d929a67d17d6", "", boost::python::no_init);
    class_ae4841ee8bd85189b310d929a67d17d6.def("assign", method_pointer_5219a7f9dd6e5bc182199bcf34e23a1b, "");
    class_ae4841ee8bd85189b310d929a67d17d6.def("__len__", method_pointer_e82624f3a0c850ea86b6188e72924462, "");
    class_ae4841ee8bd85189b310d929a67d17d6.def("capacity", method_pointer_81bb1441552c50c79761ba6be07f5bfa, "");
    class_ae4841ee8bd85189b310d929a67d17d6.def("empty", method_pointer_22a8bc2066a05cc096f4b0784d96a924, "");
    class_ae4841ee8bd85189b310d929a67d17d6.def("max_size", method_pointer_e227b09d898750208e1a33e79fe5c359, "");
    class_ae4841ee8bd85189b310d929a67d17d6.def("reserve", method_pointer_83ef44bec57455448b7e487a168437f8, "");
    class_ae4841ee8bd85189b310d929a67d17d6.def("at", method_pointer_6495e8031e6a586397e5a5b434c10381, boost::python::return_internal_reference<>(), "");
    class_ae4841ee8bd85189b310d929a67d17d6.def("at", autowig::method_decorator_6495e8031e6a586397e5a5b434c10381);
    class_ae4841ee8bd85189b310d929a67d17d6.def("at", method_pointer_47fe2297249154a5a0b6bdf86d6fc848, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_ae4841ee8bd85189b310d929a67d17d6.def("front", method_pointer_00b8f41fa0f45c55bd15b1778a4cd00e, boost::python::return_internal_reference<>(), "");
    class_ae4841ee8bd85189b310d929a67d17d6.def("front", autowig::method_decorator_00b8f41fa0f45c55bd15b1778a4cd00e);
    class_ae4841ee8bd85189b310d929a67d17d6.def("front", method_pointer_025c4ba3f3fa53bab0dbcada0ca359f1, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_ae4841ee8bd85189b310d929a67d17d6.def("back", method_pointer_478e344e518b5381a3c45e535f065b59, boost::python::return_internal_reference<>(), "");
    class_ae4841ee8bd85189b310d929a67d17d6.def("back", autowig::method_decorator_478e344e518b5381a3c45e535f065b59);
    class_ae4841ee8bd85189b310d929a67d17d6.def("back", method_pointer_8402a2ee7ae85f688e50b48a8d2633cd, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_ae4841ee8bd85189b310d929a67d17d6.def("data", method_pointer_0e308e3b5cee5630a41cdd5760b3ea50, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_ae4841ee8bd85189b310d929a67d17d6.def("data", method_pointer_5ea394f283175a28b533129752257a3e, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_ae4841ee8bd85189b310d929a67d17d6.def("push_back", method_pointer_3c2aed4cc2e7526eae27d6fcf142e32e, "");
    class_ae4841ee8bd85189b310d929a67d17d6.def("pop_back", method_pointer_bf0adbdf1bb8501e842960a0d96161f2, "");
    class_ae4841ee8bd85189b310d929a67d17d6.def("clear", method_pointer_d491e2e8cad55edc816efcce40dfdaba, "");
    class_ae4841ee8bd85189b310d929a67d17d6.def("swap", method_pointer_a6d568faf5c45dc4ac25ca463e37e9f3, "");
    class_ae4841ee8bd85189b310d929a67d17d6.def("invariants", method_pointer_ff6af1af31db560babcbdd102f79a354, "");
    class_ae4841ee8bd85189b310d929a67d17d6.def("assign", method_pointer_539a166ed5dd55639d7eb9eaeec18ff2, "");
    class_ae4841ee8bd85189b310d929a67d17d6.def("assign", method_pointer_c1a9da50f7405c088097559461c1a892, "");
    class_ae4841ee8bd85189b310d929a67d17d6.def("assign", method_pointer_0c396dc8a69d519995ee0ebe3d23830b, "");

    struct vector_ae4841ee8bd85189b310d929a67d17d6_from_python
    {
        vector_ae4841ee8bd85189b310d929a67d17d6_from_python()
        {
            boost::python::converter::registry::push_back(
                &convertible,
                &construct,
                boost::python::type_id< class ::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > > >());
        }

        static void* convertible(PyObject* obj_ptr)
        { return obj_ptr; }

        static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
        {
            boost::python::handle<> obj_iter(PyObject_GetIter(obj_ptr));
            void* storage = ((boost::python::converter::rvalue_from_python_storage< class ::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > > >*)data)->storage.bytes;
            new (storage) class ::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >();
            data->convertible = storage;
            class ::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >& result = *((class ::std::vector< class ::std::vector< double, class ::std::allocator< double > >, class ::std::allocator< class ::std::vector< double, class ::std::allocator< double > > > >*)storage);
            unsigned int i = 0;
            for(;; i++)
            {
                boost::python::handle<> py_elem_hdl(boost::python::allow_null(PyIter_Next(obj_iter.get())));
                if(PyErr_Occurred())
                { boost::python::throw_error_already_set(); }
                if(!py_elem_hdl.get())
                { break; }
                boost::python::object py_elem_obj(py_elem_hdl);
                result.push_back(boost::python::extract< class ::std::vector< double, class ::std::allocator< double > >  >(py_elem_obj));
            }
        }
    };

    vector_ae4841ee8bd85189b310d929a67d17d6_from_python();
}