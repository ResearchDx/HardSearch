/*
 * SoftClipRead.h
 *
 *  Created on: Jun 20, 2014
 *      Author: galaxy
 */
#include "ReadSequence.h"
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdint.h>

#ifndef SOFTCLIPREAD_H_
#define SOFTCLIPREAD_H_

class SoftClipRead : public ReadSequence {
public:

	template<typename T>
	static std::string NumToString(T number)
	{
		std::string sResult;
		std::ostringstream tempStream;
		tempStream << std::fixed << std::setprecision(0) << number;
		sResult = tempStream.str();
		return sResult;
	}


	SoftClipRead();
	SoftClipRead(ReadSequence rs);
	virtual ~SoftClipRead();

	enum {
		LEFT_CLIP = 0x0001,
		RIGHT_CLIP = 0x0002,
		MIDDLE_CLIP = 0x0004,
		MULTIPLE_CLIPS = 0x0008
	};

	//std::vector<std::string> cigarList;
	std::string softclip_sequence;
	int softclip_flag;
	long unsigned softclip_start;
	long unsigned softclip_end;
	uint16_t softclip_size;


	// Returns the bed format representation of the softclip
	// Example: CHR	START	START+1	SEQUENCE|DIRECTION
	std::string GetSoftClipBedFormat();

	void FixSoftClip(std::string reference);


	std::string GetReference(std::string sReferenceGenome);

	// check to see if the softclip occurs on the left
	// or right end of the read, if there are multiple softclip
	// events then return false. Only interested in softclips
	// that occurs at the end
	bool IsEndClip();
	bool IsLeftClip();
	bool IsRightClip();
	bool SizeIsGreaterThan(int unsigned size);
	bool HasMultipleSoftClips();

	// pass into std::sort to sort softclip reads by
	// chromosome and then start location
	static bool SortSoftClippedReads(SoftClipRead scRead1, SoftClipRead scRead2);
	static bool IsNumber(std::string query);

private:
	// sets the direction of the softclip
	void _FindDirection();

	// sets the bitflag for multiple softclips if present
	void _DetermineSoftClips();

	// sets the start and end coordinate of the softclip
	void _FindPosition();

	// sets the softclip sequence based on the start and end coordinates
	// of the softclip
	void _FindSoftClipSequence();



};

#endif /* SOFTCLIPREAD_H_ */
