/// \file

#ifndef asilva_defs_h
#define asilva_defs_h

// "For this is how God loved the world:
// he gave his only Son, so that everyone
// who believes in him may not perish
// but may have eternal life."
//
// John 3:16

///\mainpage
/// **aRibeiroCore**
///
/// Math/General Utils lib.
///

#include <aRibeiroCore/buildFlags.h>

#include <aRibeiroCore/constants.h>
#include <stdlib.h> // NULL definition

//#pragma warning( disable : <number> )

#ifdef _WIN32
    #pragma warning(disable:4996)
    #pragma warning(disable:4244)
    #pragma warning(disable:4309)
    #pragma warning(disable:4018)
#endif

#ifdef _WIN32
    #define WIN32_LEAN_AND_MEAN
#include <WinSock2.h>
#include <WS2tcpip.h>
    #include <Windows.h>
    #include <wchar.h>
    #ifndef swprintf
        #define swprintf _snwprintf
    #endif
#endif

//#include <stdint.h>
#include <stdint.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include <stdarg.h>
#include <stdexcept>


#if defined(OS_TARGET_linux) || defined(OS_TARGET_mac)

    #include <sys/errno.h>

#endif

namespace aRibeiro {
    extern void (*OnAbortBeforeExit)();
    extern void (*OnAbortFNC)(const char* file, int line, const char* format, ...);
    void DefaultAbortFNC(const char* file, int line, const char* format, ...);
}


/// \brief Exit application if boolean expression is true.
///
/// Tests if a boolean expression is true or false.
///
/// If the expression is true then it calls `exit(-1);` from stdlib.
///
/// Example:
///
/// \code
///    ARIBEIRO_ABORT(determinant == 0, "trying to invert a singular matrix\n");
/// \endcode
///
/// \author Alessandro Ribeiro
///
#define ARIBEIRO_ABORT(bool_exp, ...) \
    if (bool_exp) {\
        aRibeiro::OnAbortFNC(__FILE__, __LINE__, __VA_ARGS__); \
        exit(-1); \
    }
        // fprintf(stderr, "[%s:%i]\n", __FILE__, __LINE__);\
        // fprintf(stderr, __VA_ARGS__);\

#define ARIBEIRO_ABORT_LF(file,line,bool_exp, ...) \
    if (bool_exp) {\
        aRibeiro::OnAbortFNC(file, line, __VA_ARGS__); \
        exit(-1); \
    }
        // fprintf(stderr, "[%s:%i]\n", file, line);\
        // fprintf(stderr, __VA_ARGS__);\

#include <aRibeiroCore/SSE2.h>
#include <aRibeiroCore/MethodPointer.h>

namespace aRibeiro {
    /// \class aribeiro_OnDataMethodPtrType
    /// \brief Callback pattern with data and size parameters
    ///
    /// Definition of a standard data/size method pointer pattern.
    ///
    /// Example of use with functions:
    ///
    /// \code
    ///    void callbackFunction(const uint8_t *data, size_t s){
    ///        ...
    ///    }
    ///
    ///    aribeiro_OnDataMethodPtrType OnData;
    ///
    ///    OnData = &callbackFunction;
    ///
    ///    uint8_t *data;
    ///    size_t size;
    ///
    ///    ...
    ///
    ///    if (OnData != NULL)
    ///        OnData(data,size);
    /// \endcode
    ///
    /// Example of use with method:
    ///
    /// \code
    ///    class ExampleClass {
    ///    public:
    ///        void callbackFunction(const uint8_t *data, size_t s){
    ///            ...
    ///        }
    ///    };
    ///
    ///    ExampleClass obj;
    ///
    ///    aribeiro_OnDataMethodPtrType OnData;
    ///
    ///    OnData = aribeiro_OnDataMethodPtrType( &obj, &ExampleClass::callbackFunction );
    ///
    ///    uint8_t *data;
    ///    size_t size;
    ///
    ///    ...
    ///
    ///    if (OnData != NULL)
    ///        OnData(data,size);
    /// \endcode
    ///
    /// \author Alessandro Ribeiro
    ///
    DefineMethodPointer(aribeiro_OnDataMethodPtrType, void, const uint8_t *data, size_t s) VoidMethodCall(data,s)


    ARIBEIRO_INLINE
    static void* malloc_aligned(size_t size, size_t N = 32) {
        size_t complete_16bytes = (N - size % N) % N;
    #if defined(ARIBEIRO_SSE2)
        return (void*)_mm_malloc(size + complete_16bytes, N);
    #else
        return (void*)aligned_alloc(N, size + complete_16bytes);
    #endif
    }

    /*
    template <typename __T__>
    static void free_aligned(__T__ &buffer) {
    #if defined(ARIBEIRO_SSE2)
        _mm_free(buffer);
    #else
        free(buffer);
    #endif
        buffer = NULL;
    }
    */

    ARIBEIRO_INLINE
    static void free_aligned(void* buffer) {
    #if defined(ARIBEIRO_SSE2)
        _mm_free(buffer);
    #else
        free(buffer);
    #endif
    }


    ARIBEIRO_INLINE
    static int array_index_of(const uint8_t* input, int start_input, int input_size, const uint8_t* pattern, int pattern_size) {
        int test_limit = input_size - pattern_size;
        for (int i = start_input; i <= test_limit; i++) {
            if (memcmp(&input[i], pattern, pattern_size) == 0)
                return i;
        }
        return input_size;
    }
}

#endif
