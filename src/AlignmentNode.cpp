/*
 * AlignmentNode.cpp
 *
 *  Created on: May 27, 2014
 *      Author: galaxy
 */

#include "AlignmentNode.h"

AlignmentNode::AlignmentNode() {
	this->bit_score = float();
	this->match_start = int();
	this->match_end = int();
	this->gaps = int();
	this->length = int();
	this->aligned_start = long();
	this->aligned_end = long();
	this->aligned_sequence = std::string();
	this->mate_alignment = NULL;
	this->bitflag = 0;
	this->chromosome = std::string();
}

AlignmentNode::AlignmentNode(BlastAlignment blastAlignment) {
	this->bit_score = blastAlignment.bitScore;
	this->match_start = blastAlignment.matchStart;
	this->match_end = blastAlignment.matchEnd;
	this->gaps = blastAlignment.gaps;
	this->length = blastAlignment.length;
	this->aligned_start = blastAlignment.alignedStart;
	this->aligned_end = blastAlignment.alignedEnd;
	this->aligned_sequence = blastAlignment.alignedSequence;
	this->mate_alignment = NULL;
	this->bitflag = blastAlignment.bitFlag;
	this->chromosome = blastAlignment.chromosome;

}

void AlignmentNode::SetAlignmentMate(AlignmentNode *alignment)
{
	this->mate_alignment = alignment;
	this->bitflag |= AlignmentNode::BLAST_PAIRED;

	if(!this->CheckAlignmentMate())
		this->bitflag |= AlignmentNode::GAPPED_ALIGNMENT;
}

bool AlignmentNode::IsForwardStrand()
{
	return !(this->bitflag & AlignmentNode::INVERTED);
}

bool AlignmentNode::IsRightOfMate()
{
	return this->match_start > this->mate_alignment->match_start;
}

std::string AlignmentNode::Orientation()
{
	if(this->IsForwardStrand())
		return "forward";
	else
		return "reverse";
}

std::string AlignmentNode::MateOrientation()
{
	if(this->mate_alignment == NULL)
		return NULL;
	else
		return this->mate_alignment->Orientation();
}

std::string AlignmentNode::GetAlignmentCoordinates()
{
	std::string coordinates;
	coordinates = this->chromosome + ": " + HardSearch::NumToString<long unsigned>(this->aligned_start);
	coordinates += " - " + this->chromosome + ": " + HardSearch::NumToString<long unsigned>(this->aligned_end);
	return coordinates;
}

bool AlignmentNode::CheckAlignmentMate()
{
	/* if the first alignment begins at the start of the sequence
	 * check if the mate alignment begins at the base following
	 * the end of the first alignment
	 * return true if condition applies
	 */
	if(this->match_start < 25)
	{
		if(this->mate_alignment->match_start == 1)
			return true;
		else
			return false;
	}
	/* if the first alignment begins in the middle of the sequence
	 * check if the mate alignment ends at the base preceding the
	 * start of the first alignment
	 * return true if condition applies
	 */
	else if(this->match_start >= 25)
	{
		if(this->match_start - 1 == this->mate_alignment->match_end)
			return true;
		else
			return false;
	}
	return false;
}

AlignmentNode::~AlignmentNode() {
	// TODO Auto-generated destructor stub
}



