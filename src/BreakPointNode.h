/*
 * BreakPointNode.h
 *
 *  Created on: May 28, 2014
 *      Author: galaxy
 */
#ifndef LIB_INCLUDES
#define LIB_INCLUDES
#include <string>
#include <cstring>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include "AlignmentStructs.h"
#include "ReadSequence.h"
#endif


#ifndef BREAKPOINTNODE_H_
#define BREAKPOINTNODE_H_

class ReadSequence;

class BreakPointNode {
public:

	long unsigned coordinate;
	long coordinateTotal;
	std::string chromosome;
	std::vector<ReadSequence *> readsSupporting;
	int count;
	int dRange;
	bool inTargetRegion;


	BreakPointNode();
	BreakPointNode(ReadSequence &readSequence, bool direction);
	bool InGroup(std::string readChromosome, long readCoordinate);
	void AddSupportingRead(ReadSequence &readSequence);

	bool operator ==(const BreakPointNode &b);

	virtual ~BreakPointNode();


};


#endif /* BREAKPOINTNODE_H_ */
