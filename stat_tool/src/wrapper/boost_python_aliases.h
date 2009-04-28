 



#define ARGS python::args


 

// calls METHOD_NAME of the INPUT_CLASS
// :param input_class: a class instance
// :param var1: a variable 
// :returns: OUTPUT_TYPE
#define WRAP_METHOD1(INPUT_TYPE, METHOD_NAME, OUTPUT_TYPE, VARTYPE1) \
static OUTPUT_TYPE* METHOD_NAME(const Convolution& input_class, VARTYPE1 var1) \
{\
    Format_error error; \
    OUTPUT_TYPE* ret = NULL; \
    ret = input_class.METHOD_NAME(error, var1); \
    if(!ret) stat_tool::wrap_util::throw_error(error); \
    return ret; \
}
// same but no arguments
#define WRAP_METHOD0(INPUT_TYPE, METHOD_NAME, OUTPUT_TYPE) \
static OUTPUT_TYPE* METHOD_NAME(const Convolution& input_class) \
{\
    Format_error error; \
    OUTPUT_TYPE* ret = NULL; \
    ret = input_class.METHOD_NAME(error); \
    if(!ret) stat_tool::wrap_util::throw_error(error); \
    return ret; \
}


// same but  cast the result !!!!! histo->getdata is too specifi
#define  WRAP_METHOD0_CAST(Convolution, METHOD_NAME, Parametric_model); \
  static Parametric_model* METHODS_NAME(const INPUT_CLASS& convol)\
  {\
    Parametric_model* ret;\
    INPUT_CLASS_data* convol_histo = NULL;\
    convol_histo = convol.METHOD_NAME_data();\
    ret = new OUTPUT_CLASS(convol,(convol_histo ? convol_histo->METHOD_NAME() : NULL));\
    return ret;\
  }





// Specific methods


#define WRAP_METHOD_FILE_ASCII_WRITE(INPUT_CLASS) \
  static void file_ascii_write(const INPUT_CLASS& m, const char* path, bool exhaustive)\
  {\
     bool result = true;\
     Format_error error;\
     result = m.ascii_write(error, path, exhaustive);\
     if (!result)\
        stat_tool::wrap_util::throw_error(error);\
   }\


// boost python declarations -------------------------------------------
//
// quick alias to avoid writting the return_value_policy...
#define DEF_RETURN_VALUE(NAME, REFERENCE, ARGS, DOCSTRING) \
    .def(NAME, REFERENCE, return_value_policy< manage_new_object >(), ARGS, DOCSTRING)

#define DEF_RETURN_VALUE_NO_ARGS(NAME, REFERENCE, DOCSTRING) \
    .def(NAME, REFERENCE, return_value_policy< manage_new_object >(), DOCSTRING)
     




// def("__len__", &Class::method, "docstring")
#define DEF_LEN(CLASS, FUNCTION_NAME) \
    .def("__len__", &CLASS::FUNCTION_NAME, "Return the size of the Class instance") 

#define DEF_STR() .def(self_ns::str(self)) // __str__

