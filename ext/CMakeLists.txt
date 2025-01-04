# FetchContent command is available with cmake >= 3.11
include(FetchContent)

# -----------------------------------------------------------------------------
# helper function to download a 'fetch content declare'
# -----------------------------------------------------------------------------

macro(fetch what)
    FetchContent_GetProperties("${what}")
    if(NOT ${${what}_POPULATED})
        message(STATUS "fetching ${what} ...")
        FetchContent_Populate(${what})
    endif()
    mark_as_advanced(${${what}_SOURCE_DIR})
endmacro()

# -----------------------------------------------------------------------------
# dependencies
# -----------------------------------------------------------------------------

# glad
set(OpenGL_GL_PREFERENCE GLVND)
add_subdirectory(glad)

# -----------------------------------------------------------------------------
# glfw
FetchContent_Declare(
    glfw
    URL https://github.com/glfw/glfw/archive/3.3.2.zip
    URL_HASH MD5=69d0aa90a9278c0dd8be975d1ae4aaf2
    UPDATE_COMMAND "" CONFIGURE_COMMAND "" BUILD_COMMAND "" INSTALL_COMMAND ""
)
fetch(glfw)

# build options
set(GLFW_BUILD_EXAMPLES
    OFF
    CACHE BOOL "GLFW lib only"
)
set(GLFW_BUILD_TESTS
    OFF
    CACHE BOOL "GLFW lib only"
)
set(GLFW_BUILD_DOCS
    OFF
    CACHE BOOL "GLFW lib only"
)
set(GLFW_INSTALL
    OFF
    CACHE BOOL "GLFW lib only"
)

add_subdirectory(${glfw_SOURCE_DIR} glfw)
set_property(TARGET glfw PROPERTY FOLDER "thirdparty")

# -----------------------------------------------------------------------------
# glm
FetchContent_Declare(
    glm #
    URL https://github.com/g-truc/glm/releases/download/0.9.9.7/glm-0.9.9.7.zip #
    URL_HASH MD5=da4740b5b41a1e2ccb2955af74e7c298 #
)
fetch(glm)
set(GLM_TEST_ENABLE
    OFF
    CACHE BOOL ""
)
add_subdirectory(${glm_SOURCE_DIR} glm)
target_include_directories(glm INTERFACE ${glm_SOURCE_DIR}/glm>)


# -----------------------------------------------------------------------------
# imgui
#
# imgui for graphical user interface. since imgui is not a cmake project, we
# need to add target by our own. we use imgui implementation for glfw and
# opengl3

FetchContent_Declare(
    imgui #
    GIT_REPOSITORY https://github.com/ocornut/imgui.git #
    GIT_TAG v1.78 #
)
fetch(imgui)

file(
    GLOB
    IMGUI_SRC #
    ${imgui_SOURCE_DIR}/*.cpp #
    ${imgui_SOURCE_DIR}/examples/imgui_impl_glfw.cpp #
    ${imgui_SOURCE_DIR}/examples/imgui_impl_opengl3.cpp #
)

add_library(imgui STATIC ${IMGUI_SRC})

if(NOT WIN32)
    # Add this line, if GLAD should be used (is needed for fresh Ubuntu)
    target_compile_definitions(imgui PUBLIC -DIMGUI_IMPL_OPENGL_LOADER_GLAD)
endif(NOT WIN32)

set_target_properties(imgui PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories(
    imgui PUBLIC ${imgui_SOURCE_DIR} #
                 ${imgui_SOURCE_DIR}/examples #
)

target_link_libraries(
    imgui PUBLIC glfw #
                 glad #
)

set(SCP_IMGUI_FONT_FOLDER
    "\"${imgui_SOURCE_DIR}/misc/fonts\""
    CACHE STRING ""
)
set_property(TARGET imgui PROPERTY FOLDER "thirdparty")

# -----------------------------------------------------------------------------
# stb_image
add_subdirectory(stb_image)

add_subdirectory(tinyobjloader)
set_property(TARGET tinyobjloader PROPERTY FOLDER "thirdparty")

# -----------------------------------------------------------------------------
# portable file dialog
add_library(portable_file_dialogs INTERFACE)
target_include_directories(
    portable_file_dialogs
    INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/portable-file-dialog
)

# -----------------------------------------------------------------------------
# dependencies for SCP
# -----------------------------------------------------------------------------

# Eigen library
#
# download as zip, so it's a bit smaller than downloading the whole repo
FetchContent_Declare(
    eigen #
    URL https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.zip #
    URL_HASH MD5=888aab45512cc0c734b3e8f60280daba #
)
fetch(eigen)
add_library(eigen INTERFACE)
add_library(Eigen3::Eigen ALIAS eigen)
target_include_directories(eigen INTERFACE ${eigen_SOURCE_DIR})


