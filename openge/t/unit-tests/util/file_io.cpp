#include <gtest/gtest.h>
#include <util/file_io.h>



/**
 * @test file_io detectFileFormatFromFilename
 * @brief Tests the function detectFileFormatFromFilename
 */
TEST(file_io, detectFileFormatFromFilename) {
    // lower case
    {
        EXPECT_EQ(FORMAT_BAM, detectFileFormatFromFilename("my/file.bam"));
        EXPECT_EQ(FORMAT_SAM, detectFileFormatFromFilename("my/file.sam"));
        EXPECT_EQ(FORMAT_CRAM, detectFileFormatFromFilename("my/file.cram"));
        EXPECT_EQ(FORMAT_FASTQ, detectFileFormatFromFilename("my/file.fastq"));
        EXPECT_EQ(FORMAT_UNKNOWN, detectFileFormatFromFilename("my/text.file"));
        EXPECT_EQ(FORMAT_UNKNOWN, detectFileFormatFromFilename("ab"));
    }
    // upper case
    {
        EXPECT_EQ(FORMAT_BAM, detectFileFormatFromFilename("my/file.BAM"));
        EXPECT_EQ(FORMAT_SAM, detectFileFormatFromFilename("my/file.SAM"));
        EXPECT_EQ(FORMAT_CRAM, detectFileFormatFromFilename("my/file.CRAM"));
        EXPECT_EQ(FORMAT_FASTQ, detectFileFormatFromFilename("my/file.FASTQ"));
        EXPECT_EQ(FORMAT_UNKNOWN, detectFileFormatFromFilename("my/text.FILE"));
        EXPECT_EQ(FORMAT_UNKNOWN, detectFileFormatFromFilename("AB"));
    }
    // mixed case
    {
        EXPECT_EQ(FORMAT_BAM, detectFileFormatFromFilename("my/file.BaM"));
        EXPECT_EQ(FORMAT_SAM, detectFileFormatFromFilename("my/file.SAm"));
        EXPECT_EQ(FORMAT_CRAM, detectFileFormatFromFilename("my/file.CraM"));
        EXPECT_EQ(FORMAT_FASTQ, detectFileFormatFromFilename("my/file.fASTq"));
        EXPECT_EQ(FORMAT_UNKNOWN, detectFileFormatFromFilename("my/text.FIlE"));
        EXPECT_EQ(FORMAT_UNKNOWN, detectFileFormatFromFilename("aB"));
    }
}
