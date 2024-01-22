#include "_stat_tool.h"



namespace autowig
{

    void method_decorator_69c0c592e3fe59d8bdebd1d4198a02a8(class ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > > & instance, ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::size_type  param_in_0, const ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::reference  param_out) { instance.at(param_in_0) = param_out; }
    void method_decorator_c912a0ba414255519660d70e1e31d76a(class ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > > & instance, const ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::reference  param_out) { instance.front() = param_out; }
    void method_decorator_76c22ca1a9ff58faad1916d1a482e144(class ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > > & instance, const ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::reference  param_out) { instance.back() = param_out; }
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
    void  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_f1d6d85989cb57f78b4956da867c8852)(::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::size_type , ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::const_reference ) = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::assign;
    ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::size_type  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_9676a9c6810650889b9086336b17e5bb)() const = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::size;
    ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::size_type  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_44bbf6e8c648552186300f65bc147e07)() const = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::capacity;
    bool  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_11bbfe1b53f357f991113f0bf98ff25f)() const = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::empty;
    ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::size_type  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_700bb88ac01f57cb84abb0a3c9fc5bcf)() const = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::max_size;
    void  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_1b0de924887053a697e9f6616245b3c8)(::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::size_type ) = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::reserve;
    ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::reference  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_69c0c592e3fe59d8bdebd1d4198a02a8)(::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::size_type ) = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::at;
    ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::const_reference  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_00cb4a69c5f45dc6b1afc06b5c2d06a7)(::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::size_type ) const = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::at;
    ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::reference  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_c912a0ba414255519660d70e1e31d76a)() = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::front;
    ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::const_reference  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_e960959cc83e58f4bc6de8091821add3)() const = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::front;
    ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::reference  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_76c22ca1a9ff58faad1916d1a482e144)() = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::back;
    ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::const_reference  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_b9b4aeed075c5d1db57e8991446e2e8b)() const = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::back;
    void  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_b964eb1d6273562a9a3105a6e2de988e)(::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::const_reference ) = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::push_back;
    void  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_6629667764ae5533bad150ebfdc42db4)() = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::pop_back;
    void  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_763528fe19935424b071be4f79cf279a)() = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::clear;
    void  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_05db417bf20a54259824709ca71411fb)(class ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > > &) = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::swap;
    bool  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_825c4ed16b58583ba3b038844596ade1)() const = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::__invariants;
    void  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_e8453f36a48a58c5b1584366019fb5ce)(enum ::stat_tool::discrete_parametric *, enum ::stat_tool::discrete_parametric *) = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::assign;
    void  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_7923177fd23152ba842a7fef14a360fb)(class ::std::move_iterator< class ::std::__wrap_iter< enum ::stat_tool::discrete_parametric * > > , class ::std::move_iterator< class ::std::__wrap_iter< enum ::stat_tool::discrete_parametric * > > ) = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::assign;
    void  (::std::vector< enum ::stat_tool::discrete_parametric, ::std::allocator< enum ::stat_tool::discrete_parametric > >::*method_pointer_671136fe8e9853fcb08c4c166ee3be80)(enum ::stat_tool::discrete_parametric const *, enum ::stat_tool::discrete_parametric const *) = &::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >::assign;
    boost::python::class_< class ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > >, autowig::Held< class ::std::vector< enum ::stat_tool::discrete_parametric, class ::std::allocator< enum ::stat_tool::discrete_parametric > > >::Type > class_d3fbcd4a393754ca8d7be71823564225("_Vector_d3fbcd4a393754ca8d7be71823564225", "", boost::python::no_init);
    class_d3fbcd4a393754ca8d7be71823564225.def("assign", method_pointer_f1d6d85989cb57f78b4956da867c8852, "");
    class_d3fbcd4a393754ca8d7be71823564225.def("__len__", method_pointer_9676a9c6810650889b9086336b17e5bb, "");
    class_d3fbcd4a393754ca8d7be71823564225.def("capacity", method_pointer_44bbf6e8c648552186300f65bc147e07, "");
    class_d3fbcd4a393754ca8d7be71823564225.def("empty", method_pointer_11bbfe1b53f357f991113f0bf98ff25f, "");
    class_d3fbcd4a393754ca8d7be71823564225.def("max_size", method_pointer_700bb88ac01f57cb84abb0a3c9fc5bcf, "");
    class_d3fbcd4a393754ca8d7be71823564225.def("reserve", method_pointer_1b0de924887053a697e9f6616245b3c8, "");
    class_d3fbcd4a393754ca8d7be71823564225.def("at", method_pointer_69c0c592e3fe59d8bdebd1d4198a02a8, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_d3fbcd4a393754ca8d7be71823564225.def("at", autowig::method_decorator_69c0c592e3fe59d8bdebd1d4198a02a8);
    class_d3fbcd4a393754ca8d7be71823564225.def("at", method_pointer_00cb4a69c5f45dc6b1afc06b5c2d06a7, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_d3fbcd4a393754ca8d7be71823564225.def("front", method_pointer_c912a0ba414255519660d70e1e31d76a, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_d3fbcd4a393754ca8d7be71823564225.def("front", autowig::method_decorator_c912a0ba414255519660d70e1e31d76a);
    class_d3fbcd4a393754ca8d7be71823564225.def("front", method_pointer_e960959cc83e58f4bc6de8091821add3, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_d3fbcd4a393754ca8d7be71823564225.def("back", method_pointer_76c22ca1a9ff58faad1916d1a482e144, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_d3fbcd4a393754ca8d7be71823564225.def("back", autowig::method_decorator_76c22ca1a9ff58faad1916d1a482e144);
    class_d3fbcd4a393754ca8d7be71823564225.def("back", method_pointer_b9b4aeed075c5d1db57e8991446e2e8b, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_d3fbcd4a393754ca8d7be71823564225.def("push_back", method_pointer_b964eb1d6273562a9a3105a6e2de988e, "");
    class_d3fbcd4a393754ca8d7be71823564225.def("pop_back", method_pointer_6629667764ae5533bad150ebfdc42db4, "");
    class_d3fbcd4a393754ca8d7be71823564225.def("clear", method_pointer_763528fe19935424b071be4f79cf279a, "");
    class_d3fbcd4a393754ca8d7be71823564225.def("swap", method_pointer_05db417bf20a54259824709ca71411fb, "");
    class_d3fbcd4a393754ca8d7be71823564225.def("invariants", method_pointer_825c4ed16b58583ba3b038844596ade1, "");
    class_d3fbcd4a393754ca8d7be71823564225.def("assign", method_pointer_e8453f36a48a58c5b1584366019fb5ce, "");
    class_d3fbcd4a393754ca8d7be71823564225.def("assign", method_pointer_7923177fd23152ba842a7fef14a360fb, "");
    class_d3fbcd4a393754ca8d7be71823564225.def("assign", method_pointer_671136fe8e9853fcb08c4c166ee3be80, "");

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