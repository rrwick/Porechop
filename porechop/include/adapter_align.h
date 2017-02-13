#ifndef ADAPTER_ALIGN_H
#define ADAPTER_ALIGN_H

#include <seqan/sequence.h>
#include <string>
#include <vector>
#include "alignment.h"

using namespace seqan;

// Functions that are called by the Python script must have C linkage, not C++ linkage.
extern "C" {
    char * adapterAlignment(char * readSeq, char * adapterSeq,
                            int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore);
    void freeCString(char * p);
}

char * cppStringToCString(std::string cpp_string);


#endif // ADAPTER_ALIGN_H
