#include <gtest/gtest.h>
#include <util/picard_structures.h>
#include <vector>
using namespace std;


/**
 * @test ReadEnds test
 * @brief Tests the class 'ReadEnds'
 */
TEST(ReadEnds, test) {
    // test compare function
    {
        ReadEnds a, b;
        EXPECT_EQ(0, ReadEnds::compare(a, b));
        EXPECT_FALSE(a < b);
    }
    {
        ReadEnds a, b;
        b.libraryId = 1;
        EXPECT_EQ(-2, ReadEnds::compare(a, b));
        EXPECT_TRUE(a < b);
    }
    {
        ReadEnds a, b;
        b.read1Sequence = 1;
        EXPECT_EQ(-2, ReadEnds::compare(a, b));
        EXPECT_TRUE(a < b);
    }
    {
        ReadEnds a, b;
        b.read1Coordinate = 1;
        EXPECT_EQ(-2, ReadEnds::compare(a, b));
        EXPECT_TRUE(a < b);
    }
    {
        ReadEnds a, b;
        b.orientation = RE_F;
        EXPECT_EQ(-1, ReadEnds::compare(a, b));
        EXPECT_TRUE(a < b);
    }
    {
        ReadEnds a, b;
        b.read2Sequence = 1;
        EXPECT_EQ(-2, ReadEnds::compare(a, b));
        EXPECT_TRUE(a < b);
    }
    {
        ReadEnds a, b;
        b.read2Coordinate = 1;
        EXPECT_EQ(-2, ReadEnds::compare(a, b));
        EXPECT_TRUE(a < b);
    }
    {
        ReadEnds a, b;
        b.read1IndexInFile = 1;
        EXPECT_EQ(-2, ReadEnds::compare(a, b));
        EXPECT_TRUE(a < b);
    }
    {
        ReadEnds a, b;
        b.read2IndexInFile = 1;
        EXPECT_EQ(-2, ReadEnds::compare(a, b));
        EXPECT_TRUE(a < b);
    }
    {
        ReadEnds a, b;
        a.read2IndexInFile = 1;
        EXPECT_EQ(2, ReadEnds::compare(a, b));
        EXPECT_FALSE(a < b);
    }
    {
        ReadEnds a, b;
        a.read1IndexInFile = 2;
        a.read2IndexInFile = 1;
        EXPECT_EQ(3, ReadEnds::compare(a, b));
        EXPECT_FALSE(a < b);
    }
    // test isPaired
    {
        ReadEnds b;
        EXPECT_FALSE(b.isPaired());
        b.read2Sequence = 1;
        EXPECT_TRUE(b.isPaired());
    }
}


/**
 * @test ReadEndsMap test
 * @brief Tests the class 'ReadEndsMap'
 */
TEST(ReadEndsMap, test) {
    ReadEndsMap rem;
    EXPECT_EQ((size_t)0, rem.size());
    vector<ReadEnds*> list = rem.allReadEnds();
    EXPECT_EQ((size_t)0, list.size());
    ReadEnds* a1 = new ReadEnds;
    ReadEnds* a2 = new ReadEnds;
    ReadEnds* a3 = new ReadEnds;
    rem.put(0, "a1", a1);
    EXPECT_EQ((size_t)1, rem.size());
    rem.put(0, "a2", a2);
    EXPECT_EQ((size_t)2, rem.size());
    rem.put(0, "a3", a3);
    EXPECT_EQ((size_t)3, rem.size());
    rem.remove(0, "a3");
    EXPECT_EQ((size_t)2, rem.size());
    rem.remove(0, "a2");
    EXPECT_EQ((size_t)1, rem.size());
    rem.remove(0, "a1");
    EXPECT_EQ((size_t)0, rem.size());
    rem.put(0, "a2", a2);
    list = rem.allReadEnds();
    EXPECT_EQ((size_t)1, list.size());
    EXPECT_EQ(list[0], a2);
    delete a3;
    delete a2;
    delete a1;
}
