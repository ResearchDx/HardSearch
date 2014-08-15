/*
 * SoftClipRead.cpp
 *
 *  Created on: Jun 20, 2014
 *      Author: galaxy
 */

#include "SoftClipRead.h"

SoftClipRead::SoftClipRead() {
	this->softclip_start = 0;
	this->softclip_end = 0;
	this->softclip_flag = 0;
	this->softclip_size = 0;
}

SoftClipRead::SoftClipRead(ReadSequence rs) {
	this->softclip_start = 0;
	this->softclip_end = 0;
	this->softclip_flag = 0;
	this->softclip_size = 0;

	this->SetReadValues(rs.readID,
						rs.bitflag,
						rs.tid,
						rs.chromosome,
						rs.start_pos,
						rs.cigar,
						rs.qual_score,
						rs.sequence,
						rs.insert_size);
	//TODO: put this in better location...
	if(rs.mate != NULL)
		this->SetReadMate(rs.mate);
	this->FindCigarList();
	this->_DetermineSoftClips();
	this->_FindDirection();
	this->_FindPosition();
	this->_FindSoftClipSequence();
}

bool SoftClipRead::IsLeftClip()
{
	return this->softclip_flag & SoftClipRead::LEFT_CLIP;
}

bool SoftClipRead::IsRightClip()
{
	return this->softclip_flag & SoftClipRead::RIGHT_CLIP;
}

bool SoftClipRead::HasMultipleSoftClips()
{
	return this->softclip_flag & SoftClipRead::MULTIPLE_CLIPS;
}

bool SoftClipRead::SizeIsGreaterThan(int unsigned size)
{
	return this->softclip_size > size;
}

void SoftClipRead::_FindDirection()
{

	if(!(this->softclip_flag & SoftClipRead::MULTIPLE_CLIPS))
	{
		if(this->cigarList[1] == "S")
			this->softclip_flag |= SoftClipRead::LEFT_CLIP;
		else if(this->cigarList[this->cigarList.size() - 1] == "S")
			this->softclip_flag |= SoftClipRead::RIGHT_CLIP;
		else
			this->softclip_flag |= SoftClipRead::MIDDLE_CLIP;
	}
}

void SoftClipRead::_FindPosition()
{
	if(this->IsLeftClip())
	{
		this->softclip_start = this->start_pos - this->softclip_size;
		this->softclip_end = this->start_pos - 1;
	}
	else if(this->IsRightClip())
	{
		this->softclip_start = this->start_pos + strlen(this->sequence.c_str()) - this->softclip_size;
		this->softclip_end = this->start_pos + strlen(this->sequence.c_str()) - 1;
	}
}



void SoftClipRead::_DetermineSoftClips()
{
	int unsigned softclip_count = 0;
	int unsigned value = 0;

	for(std::vector<std::string>::iterator index = this->cigarList.begin(); index != this->cigarList.end(); index++)
	{
		if(IsNumber((*index)))
		{
			value = atoi((*index).c_str());
		}
		else
		{
			if((*index) == "S")
			{
				softclip_count += 1;
			// Set the softclip size
				if (value > this->softclip_size)
					this->softclip_size = value;\
			}
		}
	}
	if(softclip_count > 1)
	{
		this->softclip_flag |= SoftClipRead::MULTIPLE_CLIPS;
	}


}


void SoftClipRead::_FindSoftClipSequence()
{
	if (!(this->softclip_flag & SoftClipRead::MULTIPLE_CLIPS))
	{
		if (this->softclip_flag & SoftClipRead::LEFT_CLIP)
		{
			this->softclip_sequence = this->sequence.substr(0, this->softclip_size);
		}
		else if (this->softclip_flag & SoftClipRead::RIGHT_CLIP)
		{
			int start = strlen(this->sequence.c_str()) - this->softclip_size;
			this->softclip_sequence = this->sequence.substr(start, this->softclip_size);
		}
	}
}

std::string SoftClipRead::GetSoftClipBedFormat()
{
	std::string output_format;
	std::string direction;
	if(this->softclip_flag & SoftClipRead::LEFT_CLIP)
		direction = "Left";
	else if(this->softclip_flag & SoftClipRead::RIGHT_CLIP)
		direction = "Right";

	output_format = this->chromosome + "\t" + SoftClipRead::NumToString<long>(this->softclip_start) +
			"\t" + SoftClipRead::NumToString<long>(this->softclip_start + 1) +
			"\t" + this->softclip_sequence + "|" + direction + "|" + SoftClipRead::NumToString<long>(this->tid) +  "\n";

	return output_format;
}

void SoftClipRead::FixSoftClip(std::string reference)
{
	// Length of the softclip sequence
	int softclip_length = strlen(this->softclip_sequence.c_str());
	// Length of the softclip sequence(?)
	int softclip_size = atoi(this->cigarList[this->cigarList.size() - 2].c_str());
	int index = 1;
	// Fix Left Clipped Reads
	if((this->softclip_flag & SoftClipRead::LEFT_CLIP) &&
			(atoi(this->cigarList[0].c_str()) > 5))
	{
		std::string scBase = this->softclip_sequence.substr(softclip_length - 1,1);
		std::string refBase = reference.substr(softclip_length - 1, 1);
		while (refBase == scBase)
		{
			index += 1;
			this->softclip_end -= 1;
			this->softclip_size -= 1;
			this->start_pos = this->start_pos - 1;
			//this->softclipSequence = this->softclipSequence.substr(0,softclip_length - index + 1);
			this->softclip_sequence = this->softclip_sequence.substr(0, this->softclip_size);
			scBase = this->softclip_sequence.substr(softclip_length - index,1);
			refBase = reference.substr(softclip_length - index, 1);
			if (index >= softclip_length)
				break;
		}
	}
	// Fix Right Clipped Reads
	// Everything should be fixed...

	else if ((this->softclip_flag & SoftClipRead::RIGHT_CLIP) &&
			(softclip_size > 5))
	{
		std::string scBase = this->softclip_sequence.substr(0,1);
		std::string refBase = reference.substr(0,1);
		while (refBase == scBase)
		{
			// why wasn't the index incremented here?
			index += 1;
			this->softclip_start += 1;
			this->softclip_size -= 1;
			this->softclip_sequence = this->softclip_sequence.substr(1, this->softclip_size);
			scBase = this->softclip_sequence.substr(0,1);
			refBase = reference.substr(index, 1);

			if (index >= softclip_length)
				break;

		}
	}
}


std::string SoftClipRead::GetReference(std::string sReferenceGenome)
{
	std::string cmd = "samtools faidx " + sReferenceGenome + " " + this->chromosome + ":" + NumToString<long>(this->softclip_start) + "-" + NumToString<long>(this->softclip_end);
	std::string reference;

	FILE *ls = popen(cmd.c_str(), "r");
	char buf[100];

	while(fgets(buf, sizeof(buf), ls) != 0)
	{
		if(buf[0] != '>')
			reference += std::string(buf, strlen(buf));
		memset(&buf[0],0,sizeof(buf));
	}
	pclose(ls);

	reference.erase(std::remove(reference.begin(), reference.end(), '\n'), reference.end());
	return reference;
}

bool SoftClipRead::IsEndClip()
{
	if(this->softclip_flag & SoftClipRead::LEFT_CLIP)
		return true;
	else if(this->softclip_flag & SoftClipRead::RIGHT_CLIP)
		return true;
	else
		return false;
}

bool SoftClipRead::SortSoftClippedReads(SoftClipRead scRead1, SoftClipRead scRead2)
{
	std::string chromosome1 = scRead1.chromosome.substr(3,std::strlen(scRead1.chromosome.c_str()) - 2);
	std::string chromosome2 = scRead2.chromosome.substr(3,std::strlen(scRead2.chromosome.c_str()) - 2);


	if(IsNumber(chromosome1) && IsNumber(chromosome2))
	{
		if (atoi(chromosome1.c_str()) == atoi(chromosome2.c_str()))
			return scRead1.softclip_start > scRead2.softclip_start;
		return atoi(chromosome1.c_str()) > atoi(chromosome2.c_str());
		return false;
	}
	else
	{
		if(chromosome1 == chromosome2)
			return scRead1.softclip_start > scRead2.softclip_start;
		return chromosome1 > chromosome2;
		return false;
	}
}


bool SoftClipRead::IsNumber(std::string query)
{

	std::string test;

	for(std::string::iterator queryIter = query.begin(); queryIter != query.end(); queryIter++)
	{
		if(std::isalpha((*queryIter)))
			return false;
	}
	return true;
}


SoftClipRead::~SoftClipRead() {
	// TODO Auto-generated destructor stub
}

