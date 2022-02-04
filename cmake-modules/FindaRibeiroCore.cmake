
find_path(ARIBEIROCORE_INCLUDE_DIR aRibeiroCore/aRibeiroCore.h)
find_library(ARIBEIROCORE_LIBRARIES NAMES aRibeiroCore)

if(ARIBEIROCORE_INCLUDE_DIR AND ARIBEIROCORE_LIBRARIES)
	set(ARIBEIROCORE_FOUND ON)
endif()

if(ARIBEIROCORE_FOUND)
	if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
		MESSAGE(STATUS "Found aRibeiroCore include:  ${ARIBEIROCORE_INCLUDE_DIR}/aRibeiroCore/aRibeiroCore.h")
		MESSAGE(STATUS "Found aRibeiroCore library: ${ARIBEIROCORE_LIBRARIES}")
	endif()
else()
	if(${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED)
		MESSAGE(FATAL_ERROR "Could NOT find aRibeiroCore development files")
	endif()
endif()
