/// \file

#ifndef asilva_defs_h
#define asilva_defs_h

///\mainpage
/// **OpenGL Starter**
///
/// C++ game, math and multimedia framework.<br />
/// It is written to work with OpenGL,<br />
/// but has several other implementations as<br />
/// SIMD vector math, threads, timers, video, etc...
///
/// Compatible with MacOS X, Linux, Windows and<br />
/// Rasberry Pi (Raspbian) .
///
/// <a href="https://github.com/A-Ribeiro/OpenGLStarter" target=�_blank� >github.com/A-Ribeiro/OpenGLStarter</a>
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
        fprintf(stderr, "[%s:%i]\n", __FILE__, __LINE__);\
        fprintf(stderr, __VA_ARGS__);\
        exit(-1); \
    }

#define ARIBEIRO_ABORT_LF(file,line,bool_exp, ...) \
    if (bool_exp) {\
        fprintf(stderr, "[%s:%i]\n", file, line);\
        fprintf(stderr, __VA_ARGS__);\
        exit(-1); \
    }

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


    static void free_aligned(void* buffer) {
    #if defined(ARIBEIRO_SSE2)
        _mm_free(buffer);
    #else
        free(buffer);
    #endif
    }
}

#endif
