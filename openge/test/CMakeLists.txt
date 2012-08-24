set( OPENGE_TEST_DATA ${PROJECT_SOURCE_DIR}/openge/test/data)
set( OPENGE_TEST_SCRIPT ${PROJECT_SOURCE_DIR}/openge/test)

## Test file formats
# SAM
# TODO: test console redirection as well
add_test(NAME oge_sam_input COMMAND openge count ${OPENGE_TEST_DATA}/simple.sam)
add_test(NAME oge_sam_input_pipe COMMAND bash ${OPENGE_TEST_SCRIPT}/pipe_assist.sh ${OPENGE_TEST_DATA}/simple.sam "${EXECUTABLE_OUTPUT_PATH}/openge count")
add_test(NAME oge_sam_input_multi COMMAND openge count ${OPENGE_TEST_DATA}/simple.sam ${OPENGE_TEST_DATA}/simple.sam)
add_test(NAME oge_sam_input_diff_sam COMMAND bash ${OPENGE_TEST_SCRIPT}/compare_sam.sh "${EXECUTABLE_OUTPUT_PATH}/openge view ${OPENGE_TEST_DATA}/simple.sam" "${OPENGE_TEST_DATA}/simple.sam" )
add_test(NAME oge_sam_input_diff_bam COMMAND bash ${OPENGE_TEST_SCRIPT}/compare_sam.sh "${EXECUTABLE_OUTPUT_PATH}/openge view ${OPENGE_TEST_DATA}/simple.bam" "${OPENGE_TEST_DATA}/simple.sam" )

# BAM
add_test(NAME oge_bam_input COMMAND openge count ${OPENGE_TEST_DATA}/simple.bam)
add_test(NAME oge_bam_input_pipe COMMAND bash ${OPENGE_TEST_SCRIPT}/pipe_assist.sh ${OPENGE_TEST_DATA}/simple.bam "${EXECUTABLE_OUTPUT_PATH}/openge count")
add_test(NAME oge_bam_input_multi COMMAND openge count ${OPENGE_TEST_DATA}/simple.bam ${OPENGE_TEST_DATA}/simple.bam)
add_test(NAME oge_bam_truncated COMMAND openge count ${OPENGE_TEST_DATA}/208.truncated.bam)
add_test(NAME oge_bam_samtools_truncated COMMAND samtools view ${OPENGE_TEST_DATA}/208.truncated.bam)
set_tests_properties(oge_bam_truncated oge_bam_samtools_truncated PROPERTIES WILL_FAIL true)

## Test count command
add_test(NAME oge_count COMMAND openge count ${OPENGE_TEST_DATA}/simple.bam)

## Test coverage command
add_test(NAME oge_coverage COMMAND openge coverage ${OPENGE_TEST_DATA}/simple.bam)

## Test dedup command
add_test(NAME oge_dedup COMMAND openge dedup ${OPENGE_TEST_DATA}/simple.bam -o /dev/null)

## Test help command
add_test(NAME oge_help_count COMMAND openge help count)
add_test(NAME oge_help_coverage COMMAND openge help coverage)
add_test(NAME oge_help_dedup COMMAND openge help dedup)
add_test(NAME oge_help_mergesort COMMAND openge help mergesort)
add_test(NAME oge_help_stats COMMAND openge help stats)
add_test(NAME oge_help_version COMMAND openge help version)
add_test(NAME oge_help_view COMMAND openge help view)

## Test mergesort command
add_test(NAME oge_mergesort COMMAND openge mergesort ${OPENGE_TEST_DATA}/simple.bam -o /dev/null)

## Test stats command
add_test(NAME oge_stats COMMAND openge stats ${OPENGE_TEST_DATA}/simple.bam)
add_test(NAME oge_stats_inserts COMMAND openge stats -I ${OPENGE_TEST_DATA}/simple.bam)
add_test(NAME oge_stats_length COMMAND openge stats -L ${OPENGE_TEST_DATA}/simple.bam)

## Test version command
add_test(NAME oge_version COMMAND openge version)

## Test view command
add_test(NAME oge_view COMMAND openge view ${OPENGE_TEST_DATA}/simple.bam -o /dev/null)
add_test(NAME oge_view_length COMMAND openge view ${OPENGE_TEST_DATA}/simple.bam -o /dev/null -n 1) #TODO- check length
