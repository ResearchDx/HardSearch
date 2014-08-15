#ifndef LIB_INCLUDES
#define LIB_INCLUDES
#include <cstdlib>
#include "ReadSequence.h"
#endif

struct BlastAlignment
{
	float bitScore = 0.0;
	int matchStart = NULL;
	int matchEnd = NULL;
	int gaps = NULL;
	int length = NULL;
	long alignedStart = NULL;
	long alignedEnd = NULL;
	int bitFlag = NULL;
	std::string chromosome;
	std::string alignedSequence;
};


