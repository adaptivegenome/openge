#ifndef PICARD_STRUCTURES_H
#define PICARD_STRUCTURES_H

#include <map>
#include <vector>
#include <iostream>

typedef enum
{
    RE_NONE, RE_F, RE_R, RE_FF, RE_RR, RE_FR, RE_RF
} readends_orientation_t;

class ReadEnds{
public:
    
    short libraryId;
    short score;
    readends_orientation_t orientation;
    int read1Sequence;
    int read1Coordinate;
    long read1IndexInFile;
    int read2Sequence;
    int read2Coordinate;
    long read2IndexInFile;
    
    ReadEnds()
    : libraryId(-1)
    , score(-1)
    , orientation(RE_NONE)
    , read1Sequence(-1)
    , read1Coordinate(-1)
    , read1IndexInFile(-1)
    , read2Sequence(-1)
    , read2Coordinate(-1)
    , read2IndexInFile(-1)
    {}
    
    bool isPaired() { return read2Sequence != -1; }
    
    static int compare(const ReadEnds & lhs, const ReadEnds & rhs) {
        int retval = 0;
        if (retval == 0) retval = lhs.libraryId - rhs.libraryId;
        if (retval == 0) retval = lhs.read1Sequence - rhs.read1Sequence;
        if (retval == 0) retval = lhs.read1Coordinate - rhs.read1Coordinate;
        if (retval == 0) retval = lhs.orientation - rhs.orientation;
        if (retval == 0) retval = lhs.read2Sequence   - rhs.read2Sequence;
        if (retval == 0) retval = lhs.read2Coordinate - rhs.read2Coordinate;
        if (retval == 0) retval = (int) (lhs.read1IndexInFile - rhs.read1IndexInFile);
        if (retval == 0) retval = (int) (lhs.read2IndexInFile - rhs.read2IndexInFile);
        
        return retval;
    }
    
    bool operator<(const ReadEnds & a) const
    {
        return 0 > compare(*this, a);
    }
};
std::ostream& operator<< (std::ostream& out, const ReadEnds & re );

struct compareReadEnds {
    bool operator ()(const ReadEnds *lhs, const ReadEnds *rhs) { return *lhs < *rhs; }
    bool operator ()(const ReadEnds & lhs, const ReadEnds & rhs) { return lhs < rhs; }
};

class ReadEndsMap
{
protected:
    std::map<std::string, ReadEnds *> m;
public:
    void put(int index, std::string key, ReadEnds * val) {
        m[key] = val;
    }
    
    ReadEnds * remove(int index, std::string key)
    {
        ReadEnds * ret = m[key];
        m.erase(key);
        return ret;
    }
    
    size_t size() { return m.size();};
    
    std::vector<ReadEnds *> allReadEnds()
    {
        std::vector<ReadEnds *> ret;
        ret.reserve(m.size());
        for(std::map<std::string, ReadEnds *>::iterator i = m.begin(); i != m.end(); i++)
            ret.push_back(i->second);
        return ret;
    }
    
};

#endif