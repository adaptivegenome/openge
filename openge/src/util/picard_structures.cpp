#include "picard_structures.h"


std::ostream& operator<< (std::ostream& out, const ReadEnds & re )
{
    out << "ReadEnds (LID " << re.libraryId << ")" << std::endl;
    out << " Seq: " << re.read1Sequence << "/" << re.read2Sequence << std::endl;
    out << " Coord: " << re.read1Coordinate << "/" << re.read2Coordinate << std::endl;
    out << " Orientation: " << re.orientation << std::endl;
    out << " Score: " << re.score << std::endl;
    
    return out;
}