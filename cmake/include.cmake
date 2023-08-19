if(NOT OpenGLStarter_Integration)

    if (NOT LIBS_REPOSITORY_URL)
        set(LIBS_REPOSITORY_URL "https://github.com/A-Ribeiro/public_libs/raw/main")
    endif()

    list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake-modules)

    include(cmake/00_DetectOSAndArchitecture.cmake)
    include(cmake/01_tools.cmake)
    include(cmake/02_copy_headers_to_include_directory.cmake)
    include(cmake/03_ide_setup.cmake)
    include(cmake/detect_neon.cmake)
    include(cmake/detect_openmp.cmake)
    include(cmake/detect_rpi.cmake)
    include(cmake/detect_sse2.cmake)
    include(cmake/aribeiro_options.cmake)

    set( OpenGLStarter_Integration ON )

endif()
