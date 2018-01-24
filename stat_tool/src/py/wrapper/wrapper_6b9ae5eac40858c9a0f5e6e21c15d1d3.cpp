#include "_stat_tool.h"



namespace autowig
{

    void method_decorator_7debf7c14b9b59bda0df7817656d79e8(class ::std::vector< int, class ::std::allocator< int > > & instance, ::std::vector< int, class ::std::allocator< int > >::size_type  param_in_0, int param_out) { instance.operator[](param_in_0) = param_out; }
    void method_decorator_bb1e0852f2ca56c094260a03787426c7(class ::std::vector< int, class ::std::allocator< int > > & instance, ::std::vector< int, class ::std::allocator< int > >::size_type  param_in_0, int param_out) { instance.at(param_in_0) = param_out; }
    void method_decorator_7ec1ac72b0b05f3a9707175bcd5da0bd(class ::std::vector< int, class ::std::allocator< int > > & instance, int param_out) { instance.front() = param_out; }
    void method_decorator_ed1cf37568ed54cbbd326e6ccbe5f27d(class ::std::vector< int, class ::std::allocator< int > > & instance, int param_out) { instance.back() = param_out; }
}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::vector< int, class ::std::allocator< int > > const volatile * get_pointer<class ::std::vector< int, class ::std::allocator< int > > const volatile >(class ::std::vector< int, class ::std::allocator< int > > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_6b9ae5eac40858c9a0f5e6e21c15d1d3()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    void  (::std::vector< int, ::std::allocator< int > >::*method_pointer_3ee60599950b5555a32f72572e8ff771)(::std::vector< int, class ::std::allocator< int > >::size_type , ::std::vector< int, class ::std::allocator< int > >::value_type const &) = &::std::vector< int, class ::std::allocator< int > >::assign;
    ::std::vector< int, class ::std::allocator< int > >::size_type  (::std::vector< int, ::std::allocator< int > >::*method_pointer_2f0bd94041965427ab114d1ec9369eb1)() const = &::std::vector< int, class ::std::allocator< int > >::size;
    ::std::vector< int, class ::std::allocator< int > >::size_type  (::std::vector< int, ::std::allocator< int > >::*method_pointer_03cb2a43c5ae5df48ecc631a008fa511)() const = &::std::vector< int, class ::std::allocator< int > >::max_size;
    void  (::std::vector< int, ::std::allocator< int > >::*method_pointer_b661720aed4355f39b433f04f40c652d)(::std::vector< int, class ::std::allocator< int > >::size_type ) = &::std::vector< int, class ::std::allocator< int > >::resize;
    void  (::std::vector< int, ::std::allocator< int > >::*method_pointer_656ab36255e95c4eacd64f33eba6a02e)(::std::vector< int, class ::std::allocator< int > >::size_type , ::std::vector< int, class ::std::allocator< int > >::value_type const &) = &::std::vector< int, class ::std::allocator< int > >::resize;
    void  (::std::vector< int, ::std::allocator< int > >::*method_pointer_f9cfe6149ce85d2a9c11d29a2ff6ef88)() = &::std::vector< int, class ::std::allocator< int > >::shrink_to_fit;
    ::std::vector< int, class ::std::allocator< int > >::size_type  (::std::vector< int, ::std::allocator< int > >::*method_pointer_2d96cb90afc35aaaa142783706900e63)() const = &::std::vector< int, class ::std::allocator< int > >::capacity;
    bool  (::std::vector< int, ::std::allocator< int > >::*method_pointer_829beec6ac39542092370174938c108d)() const = &::std::vector< int, class ::std::allocator< int > >::empty;
    void  (::std::vector< int, ::std::allocator< int > >::*method_pointer_bb2b15e55a165e4590a962713b38756e)(::std::vector< int, class ::std::allocator< int > >::size_type ) = &::std::vector< int, class ::std::allocator< int > >::reserve;
    ::std::vector< int, class ::std::allocator< int > >::reference  (::std::vector< int, ::std::allocator< int > >::*method_pointer_7debf7c14b9b59bda0df7817656d79e8)(::std::vector< int, class ::std::allocator< int > >::size_type ) = &::std::vector< int, class ::std::allocator< int > >::operator[];
    ::std::vector< int, class ::std::allocator< int > >::const_reference  (::std::vector< int, ::std::allocator< int > >::*method_pointer_be70b700bb335ba8847833cc620fa92e)(::std::vector< int, class ::std::allocator< int > >::size_type ) const = &::std::vector< int, class ::std::allocator< int > >::operator[];
    ::std::vector< int, class ::std::allocator< int > >::reference  (::std::vector< int, ::std::allocator< int > >::*method_pointer_bb1e0852f2ca56c094260a03787426c7)(::std::vector< int, class ::std::allocator< int > >::size_type ) = &::std::vector< int, class ::std::allocator< int > >::at;
    ::std::vector< int, class ::std::allocator< int > >::const_reference  (::std::vector< int, ::std::allocator< int > >::*method_pointer_a36af7a241c15d6887cc6c239cd0d230)(::std::vector< int, class ::std::allocator< int > >::size_type ) const = &::std::vector< int, class ::std::allocator< int > >::at;
    ::std::vector< int, class ::std::allocator< int > >::reference  (::std::vector< int, ::std::allocator< int > >::*method_pointer_7ec1ac72b0b05f3a9707175bcd5da0bd)() = &::std::vector< int, class ::std::allocator< int > >::front;
    ::std::vector< int, class ::std::allocator< int > >::const_reference  (::std::vector< int, ::std::allocator< int > >::*method_pointer_b7cadb076a605b51b2601b9b3480c6b5)() const = &::std::vector< int, class ::std::allocator< int > >::front;
    ::std::vector< int, class ::std::allocator< int > >::reference  (::std::vector< int, ::std::allocator< int > >::*method_pointer_ed1cf37568ed54cbbd326e6ccbe5f27d)() = &::std::vector< int, class ::std::allocator< int > >::back;
    ::std::vector< int, class ::std::allocator< int > >::const_reference  (::std::vector< int, ::std::allocator< int > >::*method_pointer_aaa6ab4cb09b56e3adee1ae72ce60d90)() const = &::std::vector< int, class ::std::allocator< int > >::back;
    void  (::std::vector< int, ::std::allocator< int > >::*method_pointer_341df5e2719858f39df095cac9121eaf)(::std::vector< int, class ::std::allocator< int > >::value_type const &) = &::std::vector< int, class ::std::allocator< int > >::push_back;
    void  (::std::vector< int, ::std::allocator< int > >::*method_pointer_bbaecaa6c9535f04a1ffda1223792c23)() = &::std::vector< int, class ::std::allocator< int > >::pop_back;
    void  (::std::vector< int, ::std::allocator< int > >::*method_pointer_433012ce9fb655529590bcbd5d85150b)(class ::std::vector< int, class ::std::allocator< int > > &) = &::std::vector< int, class ::std::allocator< int > >::swap;
    void  (::std::vector< int, ::std::allocator< int > >::*method_pointer_201a5d8f6cc15fd1b83d45af764f3905)() = &::std::vector< int, class ::std::allocator< int > >::clear;
    boost::python::class_< class ::std::vector< int, class ::std::allocator< int > >, autowig::Held< class ::std::vector< int, class ::std::allocator< int > > >::Type > class_6b9ae5eac40858c9a0f5e6e21c15d1d3("_Vector_6b9ae5eac40858c9a0f5e6e21c15d1d3", "", boost::python::no_init);
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def(boost::python::init<  >(""));
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def(boost::python::init< ::std::vector< int, class ::std::allocator< int > >::allocator_type const & >(""));
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def(boost::python::init< ::std::vector< int, class ::std::allocator< int > >::size_type , ::std::vector< int, class ::std::allocator< int > >::allocator_type const & >(""));
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def(boost::python::init< ::std::vector< int, class ::std::allocator< int > >::size_type , ::std::vector< int, class ::std::allocator< int > >::value_type const &, ::std::vector< int, class ::std::allocator< int > >::allocator_type const & >(""));
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def(boost::python::init< class ::std::vector< int, class ::std::allocator< int > > const & >(""));
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def(boost::python::init< class ::std::vector< int, class ::std::allocator< int > > const &, ::std::vector< int, class ::std::allocator< int > >::allocator_type const & >(""));
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("assign", method_pointer_3ee60599950b5555a32f72572e8ff771, "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("__len__", method_pointer_2f0bd94041965427ab114d1ec9369eb1, "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("max_size", method_pointer_03cb2a43c5ae5df48ecc631a008fa511, "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("resize", method_pointer_b661720aed4355f39b433f04f40c652d, "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("resize", method_pointer_656ab36255e95c4eacd64f33eba6a02e, "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("shrink_to_fit", method_pointer_f9cfe6149ce85d2a9c11d29a2ff6ef88, "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("capacity", method_pointer_2d96cb90afc35aaaa142783706900e63, "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("empty", method_pointer_829beec6ac39542092370174938c108d, "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("reserve", method_pointer_bb2b15e55a165e4590a962713b38756e, "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("__getitem__", method_pointer_7debf7c14b9b59bda0df7817656d79e8, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("__getitem__", autowig::method_decorator_7debf7c14b9b59bda0df7817656d79e8);
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("__getitem__", method_pointer_be70b700bb335ba8847833cc620fa92e, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("at", method_pointer_bb1e0852f2ca56c094260a03787426c7, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("at", autowig::method_decorator_bb1e0852f2ca56c094260a03787426c7);
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("at", method_pointer_a36af7a241c15d6887cc6c239cd0d230, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("front", method_pointer_7ec1ac72b0b05f3a9707175bcd5da0bd, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("front", autowig::method_decorator_7ec1ac72b0b05f3a9707175bcd5da0bd);
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("front", method_pointer_b7cadb076a605b51b2601b9b3480c6b5, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("back", method_pointer_ed1cf37568ed54cbbd326e6ccbe5f27d, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("back", autowig::method_decorator_ed1cf37568ed54cbbd326e6ccbe5f27d);
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("back", method_pointer_aaa6ab4cb09b56e3adee1ae72ce60d90, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("push_back", method_pointer_341df5e2719858f39df095cac9121eaf, "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("pop_back", method_pointer_bbaecaa6c9535f04a1ffda1223792c23, "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("swap", method_pointer_433012ce9fb655529590bcbd5d85150b, "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("clear", method_pointer_201a5d8f6cc15fd1b83d45af764f3905, "");

    struct vector_6b9ae5eac40858c9a0f5e6e21c15d1d3_from_python
    {
        vector_6b9ae5eac40858c9a0f5e6e21c15d1d3_from_python()
        {
            boost::python::converter::registry::push_back(
                &convertible,
                &construct,
                boost::python::type_id< class ::std::vector< int, class ::std::allocator< int > > >());
        }

        static void* convertible(PyObject* obj_ptr)
        { return obj_ptr; }

        static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
        {
            boost::python::handle<> obj_iter(PyObject_GetIter(obj_ptr));
            void* storage = ((boost::python::converter::rvalue_from_python_storage< class ::std::vector< int, class ::std::allocator< int > > >*)data)->storage.bytes;
            new (storage) class ::std::vector< int, class ::std::allocator< int > >();
            data->convertible = storage;
            class ::std::vector< int, class ::std::allocator< int > >& result = *((class ::std::vector< int, class ::std::allocator< int > >*)storage);
            unsigned int i = 0;
            for(;; i++)
            {
                boost::python::handle<> py_elem_hdl(boost::python::allow_null(PyIter_Next(obj_iter.get())));
                if(PyErr_Occurred())
                { boost::python::throw_error_already_set(); }
                if(!py_elem_hdl.get())
                { break; }
                boost::python::object py_elem_obj(py_elem_hdl);
                result.push_back(boost::python::extract< int  >(py_elem_obj));
            }
        }
    };

    vector_6b9ae5eac40858c9a0f5e6e21c15d1d3_from_python();
}