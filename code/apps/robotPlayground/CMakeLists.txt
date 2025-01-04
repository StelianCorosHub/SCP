cmake_minimum_required(VERSION 3.5)

project(app_robotPlayground)

file(GLOB SCP_SOURCES #
     "*.h" #
     "*.cpp" #
)

list(
    APPEND
    SCP_TARGET_DEPENDENCIES #
    "scp::utils" #
    "scp::gui" #
    "scp::RBSim" #
    "scp::control" #
)

list(APPEND SCP_TARGET_INCLUDE_DIRS #
     PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

list(APPEND SCP_TARGET_LINK_LIBS #
     PUBLIC "scp::utils" #
     PUBLIC "scp::gui" #
     PUBLIC "scp::RBSim" #
     PUBLIC "scp::control"
)

list(APPEND SCP_COMPILE_DEFINITIONS #
     PUBLIC "SCP_DATA_FOLDER=${SCP_DATA_FOLDER}" #
)

create_scp_app(
    ${PROJECT_NAME}
    "${SCP_SOURCES}" #
    "${SCP_TARGET_DEPENDENCIES}" #
    "${SCP_TARGET_INCLUDE_DIRS}" #
    "${SCP_TARGET_LINK_LIBS}" #
    "${SCP_COMPILE_DEFINITIONS}"
)
