set(COMMANDS_DIR "commands")

set(OPENGE_COMMANDS_SOURCES 
  ${COMMANDS_DIR}/commands.cpp
  ${COMMANDS_DIR}/command_count.cpp
  ${COMMANDS_DIR}/command_dedup.cpp
  ${COMMANDS_DIR}/command_help.cpp
  ${COMMANDS_DIR}/command_localrealign.cpp
  ${COMMANDS_DIR}/command_mergesort.cpp
  ${COMMANDS_DIR}/command_stats.cpp
  ${COMMANDS_DIR}/command_version.cpp
  ${COMMANDS_DIR}/command_view.cpp
  PARENT_SCOPE)

