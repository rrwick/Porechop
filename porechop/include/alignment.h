
#ifndef ALIGNMENT_H
#define ALIGNMENT_H


#include <string>
#include <seqan/basic.h>
#include <seqan/align.h>

using namespace seqan;


class ScoredAlignment {
public:
    ScoredAlignment(Align<Dna5String, ArrayGaps> & alignment,
                    int readLength, int adapterLength, int score);
    std::string getString();

    int m_readLength;
    int m_adapterLength;
    int m_readStartPos;
    int m_readEndPos;
    int m_adapterStartPos;
    int m_adapterEndPos;
    int m_rawScore;
    double m_alignedRegionPercentIdentity;
    double m_fullAdapterPercentIdentity;
};

#endif // ALIGNMENT_H
