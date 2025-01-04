# we need cmake >= 3.11 to support fetch_content that downloads dependencies
cmake_minimum_required(VERSION 3.11)

# set name of the project
project(SCP)

set(CMAKE_CXX_STANDARD 17)

# data folder path
set(SCP_DATA_FOLDER
        "\"${CMAKE_CURRENT_LIST_DIR}/data\""
        CACHE STRING ""
        )


# -----------------------------------------------------------------------------
# macro and functions
# -----------------------------------------------------------------------------

# create static library named LIB_NAME
function(
    create_SCP_library #
    LIB_NAME #
    SOURCE #
    DEPENDENCY #
    INCLUDE_DIRS #
    LINK_LIBS #
    COMPILE_DEFINITIONS #
)
    add_library(${LIB_NAME} STATIC ${SOURCE})
    add_library(scp::${LIB_NAME} ALIAS ${LIB_NAME})
    add_dependencies(${LIB_NAME} ${DEPENDENCY})
    target_include_directories(${LIB_NAME} ${INCLUDE_DIRS})
    target_link_libraries(${LIB_NAME} ${LINK_LIBS})

    if(COMPILE_DEFINITIONS)
        # scope is PUBLIC so these definitions are propagated
        target_compile_definitions(${LIB_NAME} ${COMPILE_DEFINITIONS})
    endif()

    # For solution explorer in visual studios
    set_property(TARGET ${LIB_NAME} PROPERTY FOLDER "libs")
endfunction()

# create executable named APP_NAME
function(
    create_SCP_app #
    APP_NAME #
    SOURCE #
    DEPENDENCY #
    INCLUDE_DIRS #
    LINK_LIBS #
    COMPILE_DEFINITIONS #
)
    add_executable(${APP_NAME} ${SOURCE})
    add_dependencies(${APP_NAME} ${DEPENDENCY})
    target_include_directories(${APP_NAME} ${INCLUDE_DIRS})
    target_link_libraries(${APP_NAME} ${LINK_LIBS})

    if(COMPILE_DEFINITIONS)
        target_compile_definitions(${APP_NAME} ${COMPILE_DEFINITIONS})
    endif()

    # For solution explorer in visual studios
    set_property(TARGET ${APP_NAME} PROPERTY FOLDER "Apps")
endfunction()

# -----------------------------------------------------------------------------
# code
# -----------------------------------------------------------------------------

# our code
add_subdirectory(code)

# thirdparty code
add_subdirectory(ext)


set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER "thirdparty")

