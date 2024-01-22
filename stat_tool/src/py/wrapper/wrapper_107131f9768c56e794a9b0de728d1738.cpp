#include "_stat_tool.h"



namespace autowig
{

    void method_decorator_807b8927890750b8a86c6ad1f6b386c2(class ::std::vector< double, class ::std::allocator< double > > & instance, ::std::vector< double, class ::std::allocator< double > >::size_type  param_in_0, double param_out) { instance.at(param_in_0) = param_out; }
    void method_decorator_1339251591c5515a98e83a0050706699(class ::std::vector< double, class ::std::allocator< double > > & instance, double param_out) { instance.front() = param_out; }
    void method_decorator_9d019f4c78a1536db0954a15a9df10be(class ::std::vector< double, class ::std::allocator< double > > & instance, double param_out) { instance.back() = param_out; }
}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::vector< double, class ::std::allocator< double > > const volatile * get_pointer<class ::std::vector< double, class ::std::allocator< double > > const volatile >(class ::std::vector< double, class ::std::allocator< double > > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_107131f9768c56e794a9b0de728d1738()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    void  (::std::vector< double, ::std::allocator< double > >::*method_pointer_bd3b0a01510b53d6a23ffb47de810871)(::std::vector< double, class ::std::allocator< double > >::size_type , ::std::vector< double, class ::std::allocator< double > >::const_reference ) = &::std::vector< double, class ::std::allocator< double > >::assign;
    ::std::vector< double, class ::std::allocator< double > >::size_type  (::std::vector< double, ::std::allocator< double > >::*method_pointer_532fefd4df105fd5b03b5056bfbc96d2)() const = &::std::vector< double, class ::std::allocator< double > >::size;
    ::std::vector< double, class ::std::allocator< double > >::size_type  (::std::vector< double, ::std::allocator< double > >::*method_pointer_840d124381ff523a8a307316e85df988)() const = &::std::vector< double, class ::std::allocator< double > >::capacity;
    bool  (::std::vector< double, ::std::allocator< double > >::*method_pointer_194feced398858e6b0520c20d9b56f05)() const = &::std::vector< double, class ::std::allocator< double > >::empty;
    ::std::vector< double, class ::std::allocator< double > >::size_type  (::std::vector< double, ::std::allocator< double > >::*method_pointer_9a5efa6fac8555b9bbc3931a8e59029e)() const = &::std::vector< double, class ::std::allocator< double > >::max_size;
    void  (::std::vector< double, ::std::allocator< double > >::*method_pointer_57072bfde93f5e9ebb4b00b299e11b85)(::std::vector< double, class ::std::allocator< double > >::size_type ) = &::std::vector< double, class ::std::allocator< double > >::reserve;
    ::std::vector< double, class ::std::allocator< double > >::reference  (::std::vector< double, ::std::allocator< double > >::*method_pointer_807b8927890750b8a86c6ad1f6b386c2)(::std::vector< double, class ::std::allocator< double > >::size_type ) = &::std::vector< double, class ::std::allocator< double > >::at;
    ::std::vector< double, class ::std::allocator< double > >::const_reference  (::std::vector< double, ::std::allocator< double > >::*method_pointer_9c305f05d46d542b8424268d38c945d0)(::std::vector< double, class ::std::allocator< double > >::size_type ) const = &::std::vector< double, class ::std::allocator< double > >::at;
    ::std::vector< double, class ::std::allocator< double > >::reference  (::std::vector< double, ::std::allocator< double > >::*method_pointer_1339251591c5515a98e83a0050706699)() = &::std::vector< double, class ::std::allocator< double > >::front;
    ::std::vector< double, class ::std::allocator< double > >::const_reference  (::std::vector< double, ::std::allocator< double > >::*method_pointer_b7286d66b5185ce3a82408c1c08d0321)() const = &::std::vector< double, class ::std::allocator< double > >::front;
    ::std::vector< double, class ::std::allocator< double > >::reference  (::std::vector< double, ::std::allocator< double > >::*method_pointer_9d019f4c78a1536db0954a15a9df10be)() = &::std::vector< double, class ::std::allocator< double > >::back;
    ::std::vector< double, class ::std::allocator< double > >::const_reference  (::std::vector< double, ::std::allocator< double > >::*method_pointer_cfc57d90d8f95ead97edfdb47ab2158f)() const = &::std::vector< double, class ::std::allocator< double > >::back;
    void  (::std::vector< double, ::std::allocator< double > >::*method_pointer_dbfb61dd68505d29bb323f7f7721e7b8)(::std::vector< double, class ::std::allocator< double > >::const_reference ) = &::std::vector< double, class ::std::allocator< double > >::push_back;
    void  (::std::vector< double, ::std::allocator< double > >::*method_pointer_67d448aeb7eb55bcbb1a622756f299cf)() = &::std::vector< double, class ::std::allocator< double > >::pop_back;
    void  (::std::vector< double, ::std::allocator< double > >::*method_pointer_198c30f4c28f58918a7a6006b9664857)() = &::std::vector< double, class ::std::allocator< double > >::clear;
    void  (::std::vector< double, ::std::allocator< double > >::*method_pointer_e961430f7bb95b32a0e5f04c4005d642)(class ::std::vector< double, class ::std::allocator< double > > &) = &::std::vector< double, class ::std::allocator< double > >::swap;
    bool  (::std::vector< double, ::std::allocator< double > >::*method_pointer_8217412730d054eb8e62f4e4052e2174)() const = &::std::vector< double, class ::std::allocator< double > >::__invariants;
    void  (::std::vector< double, ::std::allocator< double > >::*method_pointer_c7692583b02e5a13abd18a32600e5a96)(class ::std::move_iterator< class ::std::__wrap_iter< double * > > , class ::std::move_iterator< class ::std::__wrap_iter< double * > > ) = &::std::vector< double, class ::std::allocator< double > >::assign;
    boost::python::class_< class ::std::vector< double, class ::std::allocator< double > >, autowig::Held< class ::std::vector< double, class ::std::allocator< double > > >::Type > class_107131f9768c56e794a9b0de728d1738("_Vector_107131f9768c56e794a9b0de728d1738", "", boost::python::no_init);
    class_107131f9768c56e794a9b0de728d1738.def("assign", method_pointer_bd3b0a01510b53d6a23ffb47de810871, "");
    class_107131f9768c56e794a9b0de728d1738.def("__len__", method_pointer_532fefd4df105fd5b03b5056bfbc96d2, "");
    class_107131f9768c56e794a9b0de728d1738.def("capacity", method_pointer_840d124381ff523a8a307316e85df988, "");
    class_107131f9768c56e794a9b0de728d1738.def("empty", method_pointer_194feced398858e6b0520c20d9b56f05, "");
    class_107131f9768c56e794a9b0de728d1738.def("max_size", method_pointer_9a5efa6fac8555b9bbc3931a8e59029e, "");
    class_107131f9768c56e794a9b0de728d1738.def("reserve", method_pointer_57072bfde93f5e9ebb4b00b299e11b85, "");
    class_107131f9768c56e794a9b0de728d1738.def("at", method_pointer_807b8927890750b8a86c6ad1f6b386c2, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_107131f9768c56e794a9b0de728d1738.def("at", autowig::method_decorator_807b8927890750b8a86c6ad1f6b386c2);
    class_107131f9768c56e794a9b0de728d1738.def("at", method_pointer_9c305f05d46d542b8424268d38c945d0, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_107131f9768c56e794a9b0de728d1738.def("front", method_pointer_1339251591c5515a98e83a0050706699, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_107131f9768c56e794a9b0de728d1738.def("front", autowig::method_decorator_1339251591c5515a98e83a0050706699);
    class_107131f9768c56e794a9b0de728d1738.def("front", method_pointer_b7286d66b5185ce3a82408c1c08d0321, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_107131f9768c56e794a9b0de728d1738.def("back", method_pointer_9d019f4c78a1536db0954a15a9df10be, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_107131f9768c56e794a9b0de728d1738.def("back", autowig::method_decorator_9d019f4c78a1536db0954a15a9df10be);
    class_107131f9768c56e794a9b0de728d1738.def("back", method_pointer_cfc57d90d8f95ead97edfdb47ab2158f, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_107131f9768c56e794a9b0de728d1738.def("push_back", method_pointer_dbfb61dd68505d29bb323f7f7721e7b8, "");
    class_107131f9768c56e794a9b0de728d1738.def("pop_back", method_pointer_67d448aeb7eb55bcbb1a622756f299cf, "");
    class_107131f9768c56e794a9b0de728d1738.def("clear", method_pointer_198c30f4c28f58918a7a6006b9664857, "");
    class_107131f9768c56e794a9b0de728d1738.def("swap", method_pointer_e961430f7bb95b32a0e5f04c4005d642, "");
    class_107131f9768c56e794a9b0de728d1738.def("invariants", method_pointer_8217412730d054eb8e62f4e4052e2174, "");
    class_107131f9768c56e794a9b0de728d1738.def("assign", method_pointer_c7692583b02e5a13abd18a32600e5a96, "");

    struct vector_107131f9768c56e794a9b0de728d1738_from_python
    {
        vector_107131f9768c56e794a9b0de728d1738_from_python()
        {
            boost::python::converter::registry::push_back(
                &convertible,
                &construct,
                boost::python::type_id< class ::std::vector< double, class ::std::allocator< double > > >());
        }

        static void* convertible(PyObject* obj_ptr)
        { return obj_ptr; }

        static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
        {
            boost::python::handle<> obj_iter(PyObject_GetIter(obj_ptr));
            void* storage = ((boost::python::converter::rvalue_from_python_storage< class ::std::vector< double, class ::std::allocator< double > > >*)data)->storage.bytes;
            new (storage) class ::std::vector< double, class ::std::allocator< double > >();
            data->convertible = storage;
            class ::std::vector< double, class ::std::allocator< double > >& result = *((class ::std::vector< double, class ::std::allocator< double > >*)storage);
            unsigned int i = 0;
            for(;; i++)
            {
                boost::python::handle<> py_elem_hdl(boost::python::allow_null(PyIter_Next(obj_iter.get())));
                if(PyErr_Occurred())
                { boost::python::throw_error_already_set(); }
                if(!py_elem_hdl.get())
                { break; }
                boost::python::object py_elem_obj(py_elem_hdl);
                result.push_back(boost::python::extract< double  >(py_elem_obj));
            }
        }
    };

    vector_107131f9768c56e794a9b0de728d1738_from_python();
}