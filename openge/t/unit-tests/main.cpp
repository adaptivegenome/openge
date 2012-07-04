#include <gtest/gtest.h>

/**
 * @brief Main entry point for test-suite
 * @param argc number of input args
 * @param argv args
 * @return 0 on success, non-zero on failure
 */
int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
