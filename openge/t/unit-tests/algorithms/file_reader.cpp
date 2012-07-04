#include <gtest/gtest.h>

#include <algorithms/file_reader.h>
#include <vector>

using namespace std;



/**
 * @test FileReader deduceFileFormat
 * @brief Tests the class method deduceFileFormat
 */
TEST(FileReader, deduceFileFormat) {
    // file not exists
    {
        FileReader absent;
        absent.addFile("i/dont/exist");
        EXPECT_EQ(FileReader::FORMAT_UNKNOWN, absent.deduceFileFormat());
    }
    // unknown format
    {
        FileReader unknown;
        unknown.addFile("openge/t/test-cases/algorithms/sample.unknown");
        EXPECT_EQ(FileReader::FORMAT_UNKNOWN, unknown.deduceFileFormat());
    }
    // SAM
    {
        FileReader sam;
        sam.addFile("openge/t/test-cases/algorithms/sample.sam");
        EXPECT_EQ(FileReader::FORMAT_SAM, sam.deduceFileFormat());
    }
    // BAM
    {
        FileReader bam;
        bam.addFile("openge/t/test-cases/algorithms/sample.bam");
        EXPECT_EQ(FileReader::FORMAT_BAM, bam.deduceFileFormat());
    }
}
