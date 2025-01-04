cmake_minimum_required(VERSION 3.5)

project(RBSim)

set(CMAKE_CXX_STANDARD 11)

set(SCP_TARGET_NAME ${PROJECT_NAME})

file(GLOB SCP_SOURCES #
     "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp" #
     "${CMAKE_CURRENT_SOURCE_DIR}/src/lcp/*.cpp" #
)

# dependencies
list(APPEND SCP_TARGET_DEPENDENCIES #
     "eigen" #
     "scp::utils" #
     "scp::gui" #
     "scp::optimization" #
)

# target include dirs
list(APPEND SCP_TARGET_INCLUDE_DIRS #
     PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include" #
)

# target link libs
list(APPEND SCP_TARGET_LINK_LIBS #
    PUBLIC "eigen" #
    PUBLIC "scp::utils" #
    PUBLIC "scp::gui" #
    PUBLIC "scp::optimization" #
)

# compile definitions
list(APPEND SCP_COMPILE_DEFINITIONS #
     PUBLIC "SCP_DATA_FOLDER=${SCP_DATA_FOLDER}" #
)

# create target
create_SCP_library(
    ${SCP_TARGET_NAME}
    "${SCP_SOURCES}" #
    "${SCP_TARGET_DEPENDENCIES}" #
    "${SCP_TARGET_INCLUDE_DIRS}" #
    "${SCP_TARGET_LINK_LIBS}" #
    "${SCP_COMPILE_DEFINITIONS}" #
)

