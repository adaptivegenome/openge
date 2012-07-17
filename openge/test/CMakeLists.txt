set( OPENGE_TEST_DATA ${PROJECT_SOURCE_DIR}/openge/test/data/)

## Test file formats
# SAM
# add_test(NAME oge_sam_pipe COMMAND bash -c "cat ${OPENGE_TEST_DATA}simple.sam | ${CMAKE_CURRENT_BINARY_DIR}/openge count")
# add_test(NAME oge_sam_console_redirect COMMAND openge count < ${OPENGE_TEST_DATA}simple.sam)
add_test(NAME oge_sam_input COMMAND openge count ${OPENGE_TEST_DATA}simple.sam)

# BAM
add_test(NAME oge_bam_input COMMAND openge count ${OPENGE_TEST_DATA}simple.bam)

## Test count command
add_test(NAME oge_count COMMAND openge count ${OPENGE_TEST_DATA}simple.bam)

## Test coverage command
add_test(NAME oge_coverage COMMAND openge coverage ${OPENGE_TEST_DATA}simple.bam)

## Test dedup command
add_test(NAME oge_dedup COMMAND openge dedup ${OPENGE_TEST_DATA}simple.bam -o /dev/null)

## Test help command
add_test(NAME oge_help_count COMMAND openge help count)
add_test(NAME oge_help_coverage COMMAND openge help coverage)
add_test(NAME oge_help_dedup COMMAND openge help dedup)
add_test(NAME oge_help_mergesort COMMAND openge help mergesort)
add_test(NAME oge_help_stats COMMAND openge help stats)
add_test(NAME oge_help_version COMMAND openge help version)
add_test(NAME oge_help_view COMMAND openge help view)

## Test mergesort command
add_test(NAME oge_mergesort COMMAND openge mergesort ${OPENGE_TEST_DATA}simple.bam -o /dev/null)

## Test stats command
add_test(NAME oge_stats COMMAND openge stats ${OPENGE_TEST_DATA}simple.bam)
add_test(NAME oge_stats_inserts COMMAND openge stats -I ${OPENGE_TEST_DATA}simple.bam)
add_test(NAME oge_stats_length COMMAND openge stats -L ${OPENGE_TEST_DATA}simple.bam)

## Test version command
add_test(NAME oge_version COMMAND openge version)

## Test view command
add_test(NAME oge_view COMMAND openge view ${OPENGE_TEST_DATA}simple.bam -o /dev/null)
add_test(NAME oge_view_length COMMAND openge view ${OPENGE_TEST_DATA}simple.bam -o /dev/null -n 1) #TODO- check length
