cmake_minimum_required(VERSION 3.15)
project(wcurvint_project LANGUAGES CXX C)

add_library(wcurvint SHARED)

target_include_directories(wcurvint
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/thirdpart/wsolver/include
)

file(GLOB_RECURSE WCURVINT_SRC_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.c
)
target_sources(wcurvint PRIVATE ${WCURVINT_SRC_FILES})

if(WIN32)
    set(CMAKE_CONFIGURATION_TYPES "Release" CACHE STRING "Supported build configurations" FORCE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)

    target_compile_definitions(wcurvint
        PRIVATE
            _WINDOWS
            WCURVINT_EXPORTS
    )
    
    target_link_directories(wcurvint
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/thirdpart/wsolver/lib
    )
    target_link_libraries(wcurvint
        PRIVATE
            wsolver.lib
    )
endif()

set_target_properties(wcurvint PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib
)
set_target_properties(wcurvint PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_BINARY_DIR}/bin
    ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_BINARY_DIR}/lib
)