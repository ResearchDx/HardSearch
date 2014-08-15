/*
 * AlignmentNode.h
 *
 *  Created on: May 27, 2014
 *      Author: galaxy
 */
#ifndef LIB_INCLUDES
#define LIB_INCLUDES
#include <string>
#include <cstring>
#include <vector>
#include "AlignmentStructs.h"
#include "Utils.h"
#endif


#ifndef ALIGNMENTNODE_H_
#define ALIGNMENTNODE_H_

struct BlastAlignment;


class AlignmentNode {
public:

	AlignmentNode(BlastAlignment blastAlignment);
	AlignmentNode();

	virtual ~AlignmentNode();

	float bit_score;
	int unsigned match_start;
	int unsigned match_end;
	int unsigned gaps;
	int unsigned length;
	long unsigned aligned_start;
	long unsigned aligned_end;
	enum {
		INVERTED = 0x0001,
		ALIGNED_LEFT = 0x0002,
		ALIGNED_RIGHT = 0x0004,
		GAPPED_ALIGNMENT = 0x0008,
		BLAST_PAIRED = 0x0010
	};
	int bitflag;
	std::string chromosome;
	std::string aligned_sequence;
	AlignmentNode *mate_alignment;

	AlignmentNode MateAlignment();
	void SetAlignmentMate(AlignmentNode *alignmentNode);

	// Return true if this alignment occurs on the forward strand
	// related to the reference
	bool IsForwardStrand();

	// Return true if this alignment occurs to the right of the
	// mate alignment
	bool IsRightOfMate();

	// Returns the orientation of the current alignment in string
	// format
	std::string Orientation();

	// Returns the orientation of the mate alignment if there is
	// one in string format
	std::string MateOrientation();

	std::string GetAlignmentCoordinates();



	// Check the secondary alignment against the current alignment
	// to see if it spans a breakpoing or if there is/are gaps between
	// the two
	bool CheckAlignmentMate();

};

#endif /* ALIGNMENTNODE_H_ */
