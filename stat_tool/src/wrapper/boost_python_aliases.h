
// two internal macros used by other macros
#define METHOD_HEADER(OUTPUT_TYPE) \
	{\
		StatError error; \
		OUTPUT_TYPE* ret = NULL;\

#define ERROR \
  stat_tool::wrap_util::throw_error(error); \

#define METHOD_FOOTER\
    return ret; \
	}\


// wrapping of simple methods with variable number of arguments
// :param INPUT_TYPE: self on a class
// :param METHOD_NAME: the name of the method to be called on INPUT_TYPE class
// :param OUTPUT_TYPE: type of the class that will be returned
// :param var1: optional variables
// no argument case
#define WRAP_METHOD0(INPUT_TYPE, METHOD_NAME, OUTPUT_TYPE) \
static OUTPUT_TYPE* METHOD_NAME(const INPUT_TYPE& input_class) \
    METHOD_HEADER(OUTPUT_TYPE) \
    ret = input_class.METHOD_NAME(error); \
    if(!ret) ERROR;\
    METHOD_FOOTER \
// 1 variable
#define WRAP_METHOD1(INPUT_TYPE, METHOD_NAME, OUTPUT_TYPE, VARTYPE1) \
static OUTPUT_TYPE* METHOD_NAME(const INPUT_TYPE& input_class, VARTYPE1 var1) \
    METHOD_HEADER(OUTPUT_TYPE)\
    ret = input_class.METHOD_NAME(error, var1); \
    if(!ret) ERROR;\
    METHOD_FOOTER \
// 2 variables
#define WRAP_METHOD2(INPUT_TYPE, METHOD_NAME, OUTPUT_TYPE, VARTYPE1, VARTYPE2) \
static OUTPUT_TYPE* METHOD_NAME(const INPUT_TYPE& input_class, VARTYPE1 var1, VARTYPE2 var2) \
    METHOD_HEADER(OUTPUT_TYPE)\
    ret = input_class.METHOD_NAME(error, var1, var2); \
    if(!ret) ERROR;\
    METHOD_FOOTER \
// 3 variables
#define WRAP_METHOD3(INPUT_TYPE, METHOD_NAME, OUTPUT_TYPE, VARTYPE1,VARTYPE2,VARTYPE3) \
static OUTPUT_TYPE* METHOD_NAME(const INPUT_TYPE& input_class, VARTYPE1 var1, VARTYPE2 var2, VARTYPE3 var3) \
    METHOD_HEADER(OUTPUT_TYPE)\
    ret = input_class.METHOD_NAME(error, var1, var2, var3); \
    if(!ret) ERROR;\
    METHOD_FOOTER \



// Specific methods
// ---------------------- file_ascii_write (ascii_write) method ----------------
#define WRAP_METHOD_FILE_ASCII_WRITE(INPUT_CLASS) \
  static void file_ascii_write(const INPUT_CLASS& m, const char* path, bool exhaustive)\
  {\
     bool result = true;\
     StatError error;\
     result = m.ascii_write(error, path, exhaustive);\
     if (!result) ERROR;\
   }

// ---------------------- ascii_write method -----------------------------------
#define WRAP_METHOD_ASCII_WRITE(INPUT_CLASS) \
static std::string ascii_write(const INPUT_CLASS& input, bool exhaustive)\
  {\
    std::stringstream s;\
    std::string res;\
    input.ascii_write(s, exhaustive);\
    res = s.str();\
    return res;\
  }

// ---------------------- plot_write method ------------------------------------
#define WRAP_METHOD_PLOT_WRITE(INPUT_CLASS) \
  static void plot_write(const INPUT_CLASS& input, const std::string& prefix, const std::string& title)\
  {\
	  StatError error;\
      if(! input.plot_write(error, prefix.c_str(), title.c_str()))\
         ERROR;\
  }

// ---------------------- spreadshhet_write method -----------------------------
#define WRAP_METHOD_SPREADSHEET_WRITE(INPUT_CLASS) \
  static void spreadsheet_write(const INPUT_CLASS& input, const std::string& filename)\
  {\
    StatError error;\
    if(! input.spreadsheet_write(error, filename.c_str()))\
       ERROR;\
  }


// ---------------------- survival_get_plotable --------------------------------
#define WRAP_METHOD_SURVIVAL_GET_PLOTABLE(INPUT_CLASS) \
static MultiPlotSet* survival_get_plotable(const INPUT_CLASS& p) \
  { \
    StatError error; \
    MultiPlotSet* ret = p.survival_get_plotable(error); \
    if (!ret) ERROR;\
    return ret;\
  }

// ---------------------- survival_plot_write ---------------------------------
#define WRAP_METHOD_SURVIVAL_PLOT_WRITE(INPUT_CLASS) \
  static void survival_plot_write(const INPUT_CLASS& p, const std::string& prefix, const std::string& title) \
  {\
    StatError error;\
    if(!p.survival_plot_write(error, prefix.c_str(), title.c_str()))\
     ERROR;\
  }


// ---------------------- survival_spreadsheet_write ---------------------------------------
#define WRAP_METHOD_SURVIVAL_SPREADSHEET_WRITE(INPUT_CLASS) \
 static void survival_spreadsheet_write(const INPUT_CLASS& p, const std::string& filename) \
  {\
    StatError error;\
    if(!p.survival_spreadsheet_write(error, filename.c_str()))\
      ERROR;\
  }

// ---------------------- survival_ascii_write ---------------------------------------
#define WRAP_METHOD_SURVIVAL_ASCII_WRITE(INPUT_CLASS) \
static std::string survival_ascii_write(const INPUT_CLASS& p) \
  {\
    std::stringstream s; \
    std::string res; \
    p.survival_ascii_write(s);\
    res = s.str();\
    return res;\
  }






// INPUT PARAMETER VALIDATION -------------



#define CHECK(VALUE, MIN, MAX) \
	if (VALUE< MIN || VALUE>=MAX)\
     {\
    	 PyErr_SetString(PyExc_TypeError,(error_message.str()).c_str());\
    	 throw_error_already_set();\
     }\




// boost python declarations -------------------------------------------
//
// quick alias to avoid writting the return_value_policy...


#define DEF_RETURN_VALUE(NAME, REFERENCE, ARGS, DOCSTRING) \
    .def(NAME, REFERENCE, return_value_policy< manage_new_object >(), ARGS, DOCSTRING)

#define DEF_RETURN_VALUE_NO_ARGS(NAME, REFERENCE, DOCSTRING) \
    .def(NAME, REFERENCE, return_value_policy< manage_new_object >(), DOCSTRING)

#define DEF_INIT_MAKE_CONSTRUCTOR(WRAPPED_FUNCTION, DOCSTRING) \
	.def("__init__", make_constructor(WRAPPED_FUNCTION), DOCSTRING)



