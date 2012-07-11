#include <gtest/gtest.h>
#include <util/gatk/AlignmentBlock.h>

/**
 * @test AlignmentBlock test
 * @brief Tests the class 'AlignmentBlock'
 */
TEST(AlignmentBlock, test) {
    AlignmentBlock ab(1, 2, 3);
    EXPECT_EQ(1, ab.getReadStart());
    EXPECT_EQ(2, ab.getReferenceStart());
    EXPECT_EQ(3, ab.getLength());
}
