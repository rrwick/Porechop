
#include "alignment.h"

#include <iostream>

ScoredAlignment::ScoredAlignment(Align<Dna5String, ArrayGaps> & alignment,
                                 int readLength, int adapterLength, int score):
    m_readLength(readLength), m_adapterLength(adapterLength),
    m_readStartPos(-1), m_adapterStartPos(-1), m_rawScore(score)
{
    // Extract the alignment sequences into C++ strings for constant time random access.
    std::ostringstream stream1;
    stream1 << row(alignment, 0);
    std::string readAlignment =  stream1.str();
    std::ostringstream stream2;
    stream2 << row(alignment, 1);
    std::string adapterAlignment =  stream2.str();

    int alignmentLength = std::max(readAlignment.size(), adapterAlignment.size());
    if (alignmentLength == 0)
        return;

    int alignmentStartPos = -1;
    int alignmentEndPos = -1;

    // We consider the alignment to have started when we've encountered a base in both
    // sequences (though not necessarily at the same time).
    bool readStarted = false;
    bool adapterStarted = false;
    for (int i = 0; i < alignmentLength; ++i) {
        if (readAlignment[i] != '-')
            readStarted = true;
        if (adapterAlignment[i] != '-')
            adapterStarted = true;
        if (readStarted && adapterStarted) {
            alignmentStartPos = i;
            break;
        }
    }

    // We use the same logic to see when the alignment has ended.
    bool readEnded = false;
    bool adapterEnded = false;
    for (int i = alignmentLength - 1; i >= 0; --i) {
        if (readAlignment[i] != '-')
            readEnded = true;
        if (adapterAlignment[i] != '-')
            adapterEnded = true;
        if (readEnded && adapterEnded) {
            alignmentEndPos = i;
            break;
        }
    }

    if (alignmentStartPos == -1 || alignmentEndPos == -1)
        return;

    int readBases = 0, adapterBases = 0;
    for (int i = 0; i < alignmentLength; ++i) {
        char base1 = readAlignment[i];
        char base2 = adapterAlignment[i];

        if (i == alignmentStartPos) {
            m_readStartPos = readBases;
            m_adapterStartPos = adapterBases;
        }
        if (i == alignmentEndPos) {
            m_readEndPos = readBases;
            m_adapterEndPos = adapterBases;
        }

        if (base1 != '-')
            ++readBases;
        if (base2 != '-')
            ++adapterBases;
    }
}

std::string ScoredAlignment::getString() {
    return std::to_string(m_readStartPos) + "," + 
           std::to_string(m_readEndPos) + "," + 
           std::to_string(m_adapterStartPos) + "," + 
           std::to_string(m_adapterEndPos) + "," + 
           std::to_string(m_rawScore);
}
