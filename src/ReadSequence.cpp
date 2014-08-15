#include "ReadSequence.h"

ReadSequence::ReadSequence(std::string inputreadID)
{
	this->readID = inputreadID;
	_InitMateAndBreakPoint();
}

ReadSequence::ReadSequence()
{
	this->readID = std::string();
	_InitMateAndBreakPoint();
}

void ReadSequence::SetReadValues(std::vector<std::string> read_data)
{

	if(this->readID =="")
		this->readID = read_data[0];
	this->bitflag = atoi(read_data[1].c_str());
	this->chromosome = read_data[2];
	this->start_pos = atoi(read_data[3].c_str());
	this->cigar = read_data[5];
	this->insert_size = atoi(read_data[8].c_str());
	this->qual_score = atoi(read_data[10].c_str());
	std::string sequence;
	sequence = read_data[9];
	sequence.erase(std::remove_if(sequence.begin(), sequence.end(), isspace), sequence.end());
	this->sequence = sequence;


	_InitMateAndBreakPoint();
}

void ReadSequence::SetReadValues(	std::string read_id,
									uint32_t bitflag,
									int32_t tid,
									std::string chromosome,
									int32_t start_coordinate,
									std::string cigar,
									uint32_t quality_score,
									std::string sequence,
									int32_t isize)
{
	if(this->readID =="")
		this->readID = read_id;
	this->bitflag = bitflag;
	this->tid = tid;
	this->chromosome = chromosome;
	this->start_pos = start_coordinate;
	this->cigar = cigar;
	this->qual_score = quality_score;
	this->sequence = sequence;
	this->insert_size = isize;

	_InitMateAndBreakPoint();

}

void ReadSequence::SetAlignmentValues(AlignmentNode blastAlignment)
{
	this->blastAlignments.push_back(blastAlignment);
}

void ReadSequence::FindBreakpoint()
{
	long unsigned leftBreakPoint = long();
	long unsigned rightBreakPoint = long();
	std::string leftchr = std::string();
	std::string rightchr = std::string();


	if (((this->bitflag & READ_UNMAPPED) && !(this->blastAlignments.empty())) ||
			((this->cigar.find("S") != std::string::npos) && !(this->blastAlignments.empty())))
	{

		if (this->blastAlignments[0].bitflag & AlignmentNode::ALIGNED_LEFT)
		{
			leftchr = this->blastAlignments[0].chromosome;
			leftBreakPoint = this->blastAlignments[0].aligned_end;
			if(this->blastAlignments[0].mate_alignment != NULL)
			{
				rightBreakPoint = this->blastAlignments[0].mate_alignment->aligned_start;
				rightchr = this->blastAlignments[0].mate_alignment->chromosome;
			}
		} else if (this->blastAlignments[0].bitflag & AlignmentNode::ALIGNED_RIGHT)
		{
			leftchr = this->blastAlignments[0].chromosome;
			leftBreakPoint = this->blastAlignments[0].aligned_start;
			if(this->blastAlignments[0].mate_alignment != NULL)
			{
				rightBreakPoint = this->blastAlignments[0].mate_alignment->aligned_end;
				rightchr = this->blastAlignments[0].mate_alignment->chromosome;
			}
		}

		// if the read does span a breakpoint
		// ie. maps to two different locations then report it

			if(leftchr != std::string() && rightchr != std::string())
			{
				if(this->blastAlignments[0].CheckAlignmentMate())
				{
					this->breakPoint.leftchromosome = leftchr;
					this->breakPoint.leftcoordinate = leftBreakPoint;
					this->breakPoint.rightchromosome = rightchr;
					this->breakPoint.rightcoordinate = rightBreakPoint;
					this->breakPoint.found = true;
				}

				// these alignments don't begin directly after the detected breakpoint
				// ignore them
				/*
				else
				{
					printf("The read: %s\n", this->readID.c_str());
					printf("First Alignment From %d To %d\n", this->blastAlignments[0].matchStart, this->blastAlignments[0].matchEnd);
					printf("Coordinates %s: %lu - %lu\n", this->blastAlignments[0].chromosome.c_str(), this->blastAlignments[0].alignedStart,this->blastAlignments[0].alignedEnd);
					printf("Mate Alignment From %d To %d\n", this->blastAlignments[0].mateAlignment->matchStart, this->blastAlignments[0].mateAlignment->matchEnd);
					printf("Coordinates %s: %lu - %lu\n", this->blastAlignments[0].mateAlignment->chromosome.c_str(), this->blastAlignments[0].mateAlignment->alignedStart,this->blastAlignments[0].mateAlignment->alignedEnd);
				}
				*/
		}
#ifdef VERBOSE_FINDBREAKPOINT
		printf("For Read ID: %s \n The breakpoint is between %s:%d - %s:%d \n", this->readID.c_str(), leftchr.c_str(), leftBreakPoint, rightchr.c_str(), rightBreakPoint);
#endif
	}
}


void ReadSequence::SetDiscordantRead()
{
	this->discordant_read = true;
}

bool ReadSequence::IsMateOf(ReadSequence rs)
{
	if(this->readID != rs.readID)
		return false;
	if((this->bitflag & ReadSequence::FIRST_IN_PAIR) &&
			(rs.bitflag & ReadSequence::SECOND_IN_PAIR))
		return true;
	else if ((this->bitflag & ReadSequence::SECOND_IN_PAIR) &&
			(rs.bitflag & ReadSequence::FIRST_IN_PAIR))
		return true;
	return false;
}

int* ReadSequence::GetDeletionStart()
{
	if(this->_deletion_start == -1)
		this->_FindDeletionCoordinates();
	if(this->_deletion_start == -2)
		return NULL;
	return &this->_deletion_start;
}

int* ReadSequence::GetDeletionEnd()
{
	if(this->_deletion_start == -1)
		this->_FindDeletionCoordinates();
	if(this->_deletion_start == -2)
		return NULL;
	return &this->_deletion_end;
}

void ReadSequence::_FindDeletionCoordinates()
{
	int32_t start = this->start_pos;
	int32_t base_counter = 0;
	bool first = true;

	this->FindCigarList();

	if(!this->HasDeletions())
	{
		this->_deletion_start = -2;
		return;
	}


	std::vector<std::string>::iterator c_iter;
	for ( c_iter =  this->cigarList.begin();
			c_iter != this->cigarList.end();
			c_iter++)
	{
		if (SoftClipRead::IsNumber(*c_iter))
			base_counter = atoi((*c_iter).c_str());
		else
		{
			// if softclip occurs in the middle
			// add to the length of the read else
			// ignore it
			if (first && (*c_iter) == "S")
			{
				first = false;
			}
			else
			{
				std::string types = "MIHSD";
				std::size_t found = types.find((*c_iter));
				if (found != std::string::npos && found != 4)
				{
					start += base_counter;
				}
				else
				{
					this->_deletion_start = start;
					this->_deletion_end = start + base_counter - 1;
				}
			}
		}
	}
}

void ReadSequence::PrintAlignmentInSAMFormat()
{
	printf("%s\t%d\t%s\t%d\t\%s\n", this->readID.c_str(),
									this->bitflag,
									this->chromosome.c_str(),
									this->start_pos,
									this->sequence.c_str());
}

void ReadSequence::FindCigarList()
{
	std::string value;
	std::string type;

	this->cigarList = *(new std::vector<std::string>);

	// loop through the CIGAR string and extract the alignment
	// information
	for(std::string::iterator index = this->cigar.begin(); index != this->cigar.end(); index++)
	{
		// if the character at the current position in the CIGAR
		// is a digit append it to the string
		if(std::isdigit((*index)))
			value += (*index);

		// else, the character at the current position in the CIGAR
		// represents the alignment of the preceding bases
		else
		{
			type = (*index);
			this->cigarList.push_back(value);
			this->cigarList.push_back(type);
			value = "";
		}
	}

}


bool ReadSequence::HasSoftClips()
{
	int softclip_count = std::count(this->cigar.begin(), this->cigar.end(), 'S');
	if (softclip_count > 0)
		return true;
	return false;

}

bool ReadSequence::HasDeletions()
{
	int deletion_count = std::count(this->cigar.begin(), this->cigar.end(), 'D');
	if (deletion_count > 0)
		return true;
	return false;

}

bool ReadSequence::SpansBreakPoint(BreakPoint bp)
{
	if (this->breakPoint.bpcA == bp.bpcA && this->breakPoint.bpcB == bp.bpcB)
		return true;
	else if (this->breakPoint.bpcB == bp.bpcA && this->breakPoint.bpcA == bp.bpcB)
		return true;
	else
		return false;

}

uint32_t ReadSequence::GetBitFlag()
{
	return this->bitflag;
}

bool ReadSequence::IsUnmapped()
{
	int bitflag = this->GetBitFlag();
	return bitflag & ReadSequence::READ_UNMAPPED;
}



void ReadSequence::SetReadMate(ReadSequence *ptrReadSequence)
{
	(*ptrReadSequence).mate = (this);
	this->mate = ptrReadSequence;
}


// Set the BreakPoints of the read
// Sets A first and then B
void ReadSequence::SetBreakPoint(BreakPointNode *bpc)
{
	if(this->breakPoint.bpcA == NULL)
		this->breakPoint.bpcA = bpc;
	else
		this->breakPoint.bpcB = bpc;
}

void ReadSequence::_InitMateAndBreakPoint()
{
	this->discordant_read = false;
	this->read_type = 0;
	this->mate = NULL;
	this->_deletion_start = -1;

	this->blastAlignments = std::vector<AlignmentNode>();
	this->breakPoint.bpcA = NULL;
	this->breakPoint.bpcB = NULL;
	this->breakPoint.leftchromosome = std::string();
	this->breakPoint.rightchromosome = std::string();
	this->breakPoint.leftcoordinate = long();
	this->breakPoint.rightcoordinate = long();
	this->breakPoint.found = false;
}


