#include "common.h"

namespace aRibeiro {

    bool __OnAbortBeforeExit_called = false;
    void (*OnAbortBeforeExit)() = NULL;
    void (*OnAbortFNC)(const char* file, int line, const char* format, ...) = DefaultAbortFNC;

    void DefaultAbortFNC(const char* file, int line, const char* format, ...) {
        
        va_list args;

        // std::vector<char> char_buffer;
        // va_start(args, format);
        // char_buffer.resize(vsnprintf(NULL, 0, format, args) + 1);
        // va_end(args);

        // va_start(args, format);
        // int len = vsnprintf(&char_buffer[0], char_buffer.size(), format, args);
        // va_end(args);

        fprintf(stderr, "[%s:%i]\n", file, line);
        va_start(args, format);
        vfprintf(stderr, format, args);
        va_end(args);
        fprintf(stderr, "\n");

        if (OnAbortBeforeExit != NULL && !__OnAbortBeforeExit_called) {
            __OnAbortBeforeExit_called = true;
            OnAbortBeforeExit();
        }
        exit(-1);
    }
}
