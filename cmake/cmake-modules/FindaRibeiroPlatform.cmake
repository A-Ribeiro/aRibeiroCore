# "For this is how God loved the world:
# he gave his only Son, so that everyone
# who believes in him may not perish
# but may have eternal life."
#
# John 3:16

if(ARIBEIROPLATFORM_INCLUDE_DIR AND ARIBEIROPLATFORM_LIBRARIES)
	unset(ARIBEIROPLATFORM_INCLUDE_DIR)
	unset(ARIBEIROPLATFORM_LIBRARIES)
endif()

find_path(ARIBEIROPLATFORM_INCLUDE_DIR aRibeiroPlatform/aRibeiroPlatform.h)
find_library(ARIBEIROPLATFORM_LIBRARIES NAMES aRibeiroPlatform)

if(ARIBEIROPLATFORM_INCLUDE_DIR AND ARIBEIROPLATFORM_LIBRARIES)
	set(ARIBEIROPLATFORM_FOUND ON)
endif()

if(ARIBEIROPLATFORM_FOUND)

	find_package(aRibeiroCore REQUIRED QUIET)

	#list(APPEND ARIBEIROPLATFORM_LIBRARIES pthread)
	#list(APPEND ARIBEIROPLATFORM_LIBRARIES rt)
	list(APPEND ARIBEIROPLATFORM_LIBRARIES ${ARIBEIROCORE_LIBRARIES})

	if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
		MESSAGE(STATUS "Found ${CMAKE_FIND_PACKAGE_NAME} include:  ${ARIBEIROPLATFORM_INCLUDE_DIR}/aRibeiroPlatform/aRibeiroPlatform.h")
		MESSAGE(STATUS "Found ${CMAKE_FIND_PACKAGE_NAME} library: ${ARIBEIROPLATFORM_LIBRARIES}")
	endif()
else()
	if(${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED)
		MESSAGE(FATAL_ERROR "Could NOT find ${CMAKE_FIND_PACKAGE_NAME} development files")
	endif()
endif()