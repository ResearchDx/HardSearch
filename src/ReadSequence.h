#ifndef LIB_INCLUDES
#define LIB_INCLUDES
#include <string>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cctype>
#include <stdio.h>


#include "AlignmentNode.h"
#include "BreakPointNode.h"
#include "SoftClipRead.h"
#endif

#ifndef READSEQUENCE_H_
#define READSEQUENCE_H_

#include <stdint.h>


class AlignmentNode;
class BreakPointNode;
//class SoftClipRead;

struct BreakPoint
{
	std::string leftchromosome;
	long unsigned leftcoordinate;
	std::string rightchromosome;
	long unsigned rightcoordinate;
	BreakPointNode *bpcA;
	BreakPointNode *bpcB;
	int unsigned counter;
	bool found;

};


class ReadSequence
{
	public:

	enum {
		READ_PAIRED = 0x0001,
		READ_DISCORDANT = 0x0001,			// discordant read
		READ_SOFTCLIP = 0x0002,				// softclip read
		READ_MAPPED_IN_PAIR = 0x0002,
		READ_REALIGNED_SOFTCLIP = 0x0004,
		READ_UNMAPPED = 0x0004,				// unmapped read
		MATE_UNMAPPED = 0x0008,
		READ_REVERSE_STRAND = 0x0010,
		MATE_REVERSE_STRAND = 0x0020,
		FIRST_IN_PAIR = 0x0040,
		SECOND_IN_PAIR = 0x0080,
		NOT_PRIMARY_ALIGNMENT = 0x0100,
		FAILED_PV_QC = 0x200,
		PCR_DUP = 0x400
	};

		ReadSequence();
		ReadSequence(std::string readID);
		void SetReadValues(std::vector<std::string> read_data);
		void SetAlignmentValues(AlignmentNode blastAlignment);

		void SetReadValues(	std::string read_id,
							uint32_t bitflag,
							int32_t tid,
							std::string chromosome,
							int32_t start_coordinate,
							std::string cigar,
							uint32_t quality_score,
							std::string sequence,
							int32_t isize);

		// Set the breakpoint of this read
		void SetBreakPoint(BreakPointNode *bpc);
		// Set the Mate sequence of this read
		void SetReadMate(ReadSequence *ptrReadSequence);

		// Parse the cigar of this read, split up and store the results
		void FindCigarList();
		bool HasSoftClips();
		bool HasDeletions();


		bool SpansBreakPoint(BreakPoint bp);

		// Get the bitflag as an int for comparison
		uint32_t GetBitFlag();

		// Return true if the read is unmapped from the bitflag
		bool IsUnmapped();


		// SAM Alignment information
		// Missing quality score field
		// Missing additional info section fields
		std::string readID;
		uint32_t bitflag;
		int32_t tid;
		std::string chromosome;
		int32_t start_pos;
		std::string cigar;
		std::string sequence;
		uint32_t qual_score;
		std::vector<std::string> cigarList;
		int32_t insert_size;

		int32_t read_type;
		bool discordant_read;

		ReadSequence *mate;
		BreakPoint breakPoint;

		// TODO: possible to capture all the alignments using a tree structure
		std::vector<AlignmentNode> blastAlignments;

		// Determine if a breakpoint exists within this read sequence and if it does
		// find and set the genomic coordinates of the breakpoint
		void FindBreakpoint();

		void SetDiscordantRead();

		bool IsMateOf(ReadSequence rs);

		int* GetDeletionStart();
		int* GetDeletionEnd();


		bool operator ==(const ReadSequence other) const
		{
			return this->readID == other.readID;
		}

		bool operator !=(const ReadSequence other) const
		{
			return !(this->readID == other.readID);
		}

		void PrintAlignmentInSAMFormat();

	private:
		int32_t _deletion_start;
		int32_t _deletion_end;
		void _InitMateAndBreakPoint();
		void _FindDeletionCoordinates();



};

#endif
