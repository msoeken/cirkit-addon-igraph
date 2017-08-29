find_package(PkgConfig)
pkg_check_modules(PC_LIBIGRAPH QUIET igraph)

add_cirkit_library(
  NAME cirkit_igraph
  AUTO_DIRS src
  USE
    cirkit_classical
    cirkit_cli
    ${PC_LIBIGRAPH_LIBRARY_DIRS}/libigraph${CMAKE_SHARED_LIBRARY_SUFFIX}
  INCLUDE
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src
    PUBLIC ${PC_LIBIGRAPH_INCLUDE_DIRS}
  DEFINE
    PUBLIC ADDON_IGRAPH
  COMMANDS
    classical/cli/commands/igraph_commands.hpp
)
