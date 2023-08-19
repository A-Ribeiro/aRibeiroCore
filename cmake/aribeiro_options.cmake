############################################################################
# Compile with fast rsqrt
############################################################################

set( ARIBEIRO_RSQRT_MODE STDLIB CACHE STRING "Use the rsqrt special implementation." )
set_property(CACHE ARIBEIRO_RSQRT_MODE PROPERTY STRINGS STDLIB SSE2 CARMACK)

if(ARIBEIRO_RSQRT_MODE STREQUAL SSE2)
    add_definitions(-DARIBEIRO_RSQRT_SSE2)
    message(STATUS "*********************************")
    message(STATUS "* RSQRT - USING SSE INSTRUCTION *")
    message(STATUS "*********************************")
elseif(ARIBEIRO_RSQRT_MODE STREQUAL CARMACK)
    add_definitions(-DARIBEIRO_RSQRT_CARMACK)
    message(STATUS "***********************************")
    message(STATUS "* RSQRT - USING CARMACK ALGORITHM *")
    message(STATUS "***********************************")
else()
    message(STATUS "************************")
    message(STATUS "* RSQRT - USING STDLIB *")
    message(STATUS "************************")
endif()
