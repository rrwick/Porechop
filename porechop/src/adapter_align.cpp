#include "adapter_align.h"

#include <seqan/align.h>
#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>
#include <utility>


char * adapterAlignment(char * readSeq, char * adapterSeq,
                        int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore) {
    std::string output;

    Dna5String sequenceH = readSeq;
    Dna5String sequenceV = adapterSeq;
    std::string readName = "";
    std::string refName = "";

    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), sequenceH);
    assignSource(row(alignment, 1), sequenceV);
    Score<int, Simple> scoringScheme(matchScore, mismatchScore, gapExtensionScore, gapOpenScore);

    AlignConfig<true, true, true, true> alignConfig;
    int score = globalAlignment(alignment, scoringScheme, alignConfig);

    ScoredAlignment scoredAlignment(alignment, length(readSeq), length(adapterSeq), score);
    return cppStringToCString(scoredAlignment.getString());
}


void freeCString(char * p) {
    free(p);
}


char * cppStringToCString(std::string cpp_string) {
    char * c_string = (char*)malloc(sizeof(char) * (cpp_string.size() + 1));
    std::copy(cpp_string.begin(), cpp_string.end(), c_string);
    c_string[cpp_string.size()] = '\0';
    return c_string;
}
