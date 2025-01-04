cmake_minimum_required(VERSION 3.5)

project(gui)

set(SCP_TARGET_NAME ${PROJECT_NAME})

set(SCP_SHADER_FOLDER
    "\"${CMAKE_CURRENT_SOURCE_DIR}/shaders\""
    CACHE STRING ""
)

file(
    GLOB
    SCP_SOURCES #
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp" #
    "${CMAKE_CURRENT_SOURCE_DIR}/imgui_widgets/*.cpp" #
    "${CMAKE_CURRENT_SOURCE_DIR}/shaders/*.vert" #
    "${CMAKE_CURRENT_SOURCE_DIR}/shaders/*.frag" #
)

# dependencies
list(
    APPEND
    SCP_TARGET_DEPENDENCIES #
    "glm" #
    "imgui" #
    "stb_image" #
    "portable_file_dialogs" #
    "tinyobjloader" #
    "eigen" #
    "scp::utils" #
)

# target include dirs
list(APPEND SCP_TARGET_INCLUDE_DIRS #
     PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include" #
     PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/imgui_widgets/include" #
)

# target link libs
list(
    APPEND
    SCP_TARGET_LINK_LIBS #
    PUBLIC
    "glm" #
    PUBLIC
    "imgui" #
    PUBLIC
    "stb_image" #
    PUBLIC
    "portable_file_dialogs" #
    PUBLIC
    "tinyobjloader" #
    PUBLIC
    "eigen" #
    PUBLIC
    "scp::utils" #
)

# compile definitions
list(
    APPEND
    SCP_COMPILE_DEFINITIONS #
    PUBLIC
    "SCP_DATA_FOLDER=${SCP_DATA_FOLDER}" #
    PUBLIC
    "SCP_SHADER_FOLDER=${SCP_SHADER_FOLDER}" #
    PUBLIC
    "SCP_IMGUI_FONT_FOLDER=${SCP_IMGUI_FONT_FOLDER}" #
    PUBLIC
    "IMGUI_IMPL_OPENGL_LOADER_GLAD" #
)

# we use filesystem which is "experimental" for some compilers for those we need
# to add fs lib to link libraries
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # gcc < 7: supports fs as experimentals
    # gcc = 8: supports fs but needs to link stdc++fs
    # gcc > 9: supports fs as std library
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
        list(APPEND SCP_TARGET_LINK_LIBS PUBLIC "stdc++fs")
    endif()
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    # msvc
    # TODO:
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    # llvm clang
    # ~~~
    # llvm clang needs to link stdc++fs
    list(APPEND SCP_TARGET_LINK_LIBS PUBLIC "stdc++fs")
endif()

# create target
create_scp_library(
    ${SCP_TARGET_NAME}
    "${SCP_SOURCES}" #
    "${SCP_TARGET_DEPENDENCIES}" #
    "${SCP_TARGET_INCLUDE_DIRS}" #
    "${SCP_TARGET_LINK_LIBS}" #
    "${SCP_COMPILE_DEFINITIONS}" #
)

