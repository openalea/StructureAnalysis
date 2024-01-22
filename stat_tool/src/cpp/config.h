#if defined WIN32 || defined _WIN32 || defined __CYGWIN__
  #ifdef STAT_TOOL_MAKEDLL
    #ifdef __GNUC__
      #define STAT_TOOL_API __attribute__ ((dllexport))
    #else
      #define STAT_TOOL_API __declspec(dllexport)
    #endif
  #else
    #ifdef __GNUC__
      #define STAT_TOOL_API __attribute__ ((dllimport))
    #else
      #define STAT_TOOL_API __declspec(dllimport)
    #endif
  #endif
#else
  #if __GNUC__ >= 4
    #define STAT_TOOL_API __attribute__ ((visibility ("default")))
  #else
    #define STAT_TOOL_API
  #endif
#endif