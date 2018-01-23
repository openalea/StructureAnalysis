#include "_stat_tool.h"



namespace autowig
{

    void method_decorator_045401984933540cb0e3a863f5977fdb(class ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > > & instance, ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::size_type  param_in_0, const ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::reference  param_out) { instance.operator[](param_in_0) = param_out; }
    void method_decorator_61c66d0617825b53853bf6b6eb563ce6(class ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > > & instance, ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::size_type  param_in_0, const ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::reference  param_out) { instance.at(param_in_0) = param_out; }
    void method_decorator_7018524989f25691a29daef5c6f34123(class ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > > & instance, const ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::reference  param_out) { instance.front() = param_out; }
    void method_decorator_38961475510c5bc688c8584ceb5c8447(class ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > > & instance, const ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::reference  param_out) { instance.back() = param_out; }
}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > > const volatile * get_pointer<class ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > > const volatile >(class ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_d3fbcd4a393754ca8d7be71823564225()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    void  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_b35900984de1535ba9c3272e7fbfe936)(::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::size_type , ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::value_type const &) = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::assign;
    void  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_9f0240edc2ef5942870cfaefdd5a8a0f)(class ::std::initializer_list< enum ::stat_tool::discrete_parametric > ) = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::assign;
    ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::size_type  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_8891cae2099e5949b7a6b1040345147a)() const = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::size;
    ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::size_type  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_954ec1290d81584ead36bee8220bc98c)() const = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::max_size;
    void  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_28c69a59dff55ff8b348fd6b53761520)(::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::size_type ) = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::resize;
    void  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_853de5268a015969bfca351361740099)(::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::size_type , ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::value_type const &) = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::resize;
    void  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_da0c7ebe9eb05808a0085de787316881)() = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::shrink_to_fit;
    ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::size_type  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_197be9db4592550ca0ae1739b11715a0)() const = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::capacity;
    bool  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_ec51b8c1baf05bd0ab3564545c93240c)() const = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::empty;
    void  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_dc1256f36b6a5314b6263f85191380cd)(::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::size_type ) = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::reserve;
    ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::reference  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_045401984933540cb0e3a863f5977fdb)(::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::size_type ) = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::operator[];
    ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::const_reference  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_c0a9f33260545d51a9e8b657badb8a5f)(::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::size_type ) const = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::operator[];
    ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::reference  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_61c66d0617825b53853bf6b6eb563ce6)(::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::size_type ) = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::at;
    ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::const_reference  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_4fda175be5e75e40bcacb6cb407d8040)(::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::size_type ) const = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::at;
    ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::reference  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_7018524989f25691a29daef5c6f34123)() = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::front;
    ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::const_reference  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_e3df47bc7ad3536f86e94af015548669)() const = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::front;
    ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::reference  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_38961475510c5bc688c8584ceb5c8447)() = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::back;
    ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::const_reference  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_93a2d6f1674354db9aed90f7d24e71e3)() const = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::back;
    void  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_ec5d679e21de5aa691d6c64505b382b0)(::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::value_type const &) = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::push_back;
    void  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_b4bf56b9553b5d53834d8335dea11816)() = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::pop_back;
    void  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_c80628635bb05c60bdcdc1047af036be)(class ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > > &) = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::swap;
    void  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_bdd2666fcceb50ea943bfb48822362eb)() = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::clear;
    void  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_1d8980dc2fa05fadbd88bad0f19397d9)(enum ::stat_tool::discrete_parametric const *, enum ::stat_tool::discrete_parametric const *) = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::assign;
    boost::python::class_< class ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >, autowig::Held< class ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > > >::Type > class_d3fbcd4a393754ca8d7be71823564225("_Vector_d3fbcd4a393754ca8d7be71823564225", "", boost::python::no_init);
    class_d3fbcd4a393754ca8d7be71823564225.def(boost::python::init<  >(""));
    class_d3fbcd4a393754ca8d7be71823564225.def(boost::python::init< ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::allocator_type const & >(""));
    class_d3fbcd4a393754ca8d7be71823564225.def(boost::python::init< ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::size_type , ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::allocator_type const & >(""));
    class_d3fbcd4a393754ca8d7be71823564225.def(boost::python::init< ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::size_type , ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::value_type const &, ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::allocator_type const & >(""));
    class_d3fbcd4a393754ca8d7be71823564225.def(boost::python::init< class ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > > const & >(""));
    class_d3fbcd4a393754ca8d7be71823564225.def(boost::python::init< class ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > > const &, ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::allocator_type const & >(""));
    class_d3fbcd4a393754ca8d7be71823564225.def(boost::python::init< class ::std::initializer_list< enum ::stat_tool::discrete_parametric > , ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::allocator_type const & >(""));
    class_d3fbcd4a393754ca8d7be71823564225.def("assign", method_pointer_b35900984de1535ba9c3272e7fbfe936, "");
    class_d3fbcd4a393754ca8d7be71823564225.def("assign", method_pointer_9f0240edc2ef5942870cfaefdd5a8a0f, "");
    class_d3fbcd4a393754ca8d7be71823564225.def("__len__", method_pointer_8891cae2099e5949b7a6b1040345147a, "");
    class_d3fbcd4a393754ca8d7be71823564225.def("max_size", method_pointer_954ec1290d81584ead36bee8220bc98c, "");
    class_d3fbcd4a393754ca8d7be71823564225.def("resize", method_pointer_28c69a59dff55ff8b348fd6b53761520, "");
    class_d3fbcd4a393754ca8d7be71823564225.def("resize", method_pointer_853de5268a015969bfca351361740099, "");
    class_d3fbcd4a393754ca8d7be71823564225.def("shrink_to_fit", method_pointer_da0c7ebe9eb05808a0085de787316881, "");
    class_d3fbcd4a393754ca8d7be71823564225.def("capacity", method_pointer_197be9db4592550ca0ae1739b11715a0, "");
    class_d3fbcd4a393754ca8d7be71823564225.def("empty", method_pointer_ec51b8c1baf05bd0ab3564545c93240c, "");
    class_d3fbcd4a393754ca8d7be71823564225.def("reserve", method_pointer_dc1256f36b6a5314b6263f85191380cd, "");
    class_d3fbcd4a393754ca8d7be71823564225.def("__getitem__", method_pointer_045401984933540cb0e3a863f5977fdb, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_d3fbcd4a393754ca8d7be71823564225.def("__getitem__", autowig::method_decorator_045401984933540cb0e3a863f5977fdb);
    class_d3fbcd4a393754ca8d7be71823564225.def("__getitem__", method_pointer_c0a9f33260545d51a9e8b657badb8a5f, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_d3fbcd4a393754ca8d7be71823564225.def("at", method_pointer_61c66d0617825b53853bf6b6eb563ce6, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_d3fbcd4a393754ca8d7be71823564225.def("at", autowig::method_decorator_61c66d0617825b53853bf6b6eb563ce6);
    class_d3fbcd4a393754ca8d7be71823564225.def("at", method_pointer_4fda175be5e75e40bcacb6cb407d8040, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_d3fbcd4a393754ca8d7be71823564225.def("front", method_pointer_7018524989f25691a29daef5c6f34123, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_d3fbcd4a393754ca8d7be71823564225.def("front", autowig::method_decorator_7018524989f25691a29daef5c6f34123);
    class_d3fbcd4a393754ca8d7be71823564225.def("front", method_pointer_e3df47bc7ad3536f86e94af015548669, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_d3fbcd4a393754ca8d7be71823564225.def("back", method_pointer_38961475510c5bc688c8584ceb5c8447, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_d3fbcd4a393754ca8d7be71823564225.def("back", autowig::method_decorator_38961475510c5bc688c8584ceb5c8447);
    class_d3fbcd4a393754ca8d7be71823564225.def("back", method_pointer_93a2d6f1674354db9aed90f7d24e71e3, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_d3fbcd4a393754ca8d7be71823564225.def("push_back", method_pointer_ec5d679e21de5aa691d6c64505b382b0, "");
    class_d3fbcd4a393754ca8d7be71823564225.def("pop_back", method_pointer_b4bf56b9553b5d53834d8335dea11816, "");
    class_d3fbcd4a393754ca8d7be71823564225.def("swap", method_pointer_c80628635bb05c60bdcdc1047af036be, "");
    class_d3fbcd4a393754ca8d7be71823564225.def("clear", method_pointer_bdd2666fcceb50ea943bfb48822362eb, "");
    class_d3fbcd4a393754ca8d7be71823564225.def("assign", method_pointer_1d8980dc2fa05fadbd88bad0f19397d9, "");

    struct vector_d3fbcd4a393754ca8d7be71823564225_from_python
    {
        vector_d3fbcd4a393754ca8d7be71823564225_from_python()
        {
            boost::python::converter::registry::push_back(
                &convertible,
                &construct,
                boost::python::type_id< class ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > > >());
        }

        static void* convertible(PyObject* obj_ptr)
        { return obj_ptr; }

        static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
        {
            boost::python::handle<> obj_iter(PyObject_GetIter(obj_ptr));
            void* storage = ((boost::python::converter::rvalue_from_python_storage< class ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > > >*)data)->storage.bytes;
            new (storage) class ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >();
            data->convertible = storage;
            class ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >& result = *((class ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >*)storage);
            unsigned int i = 0;
            for(;; i++)
            {
                boost::python::handle<> py_elem_hdl(boost::python::allow_null(PyIter_Next(obj_iter.get())));
                if(PyErr_Occurred())
                { boost::python::throw_error_already_set(); }
                if(!py_elem_hdl.get())
                { break; }
                boost::python::object py_elem_obj(py_elem_hdl);
                result.push_back(boost::python::extract< enum ::stat_tool::discrete_parametric  >(py_elem_obj));
            }
        }
    };

    vector_d3fbcd4a393754ca8d7be71823564225_from_python();
}