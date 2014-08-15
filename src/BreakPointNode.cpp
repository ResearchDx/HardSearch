/*
 * BreakPointNode.cpp
 *
 *  Created on: May 28, 2014
 *      Author: galaxy
 */

#include "BreakPointNode.h"


BreakPointNode::BreakPointNode() {
	this->coordinate = 0.0;
	this->coordinateTotal = 0;
	this->chromosome = std::string();
	this->count = int();
	this->dRange = int();
}

BreakPointNode::BreakPointNode(ReadSequence &readSequence, bool direction = false)
{
	if(direction)
	{
		this->coordinate = readSequence.breakPoint.rightcoordinate;
		this->coordinateTotal = readSequence.breakPoint.rightcoordinate;
		this->chromosome = readSequence.breakPoint.rightchromosome;
	}
	else
	{
	this->coordinate = readSequence.breakPoint.leftcoordinate;
	this->coordinateTotal = readSequence.breakPoint.leftcoordinate;
	this->chromosome = readSequence.breakPoint.leftchromosome;
	}
	this->count = 1;
	this->dRange = int();
	this->AddSupportingRead(readSequence);
}

bool BreakPointNode::InGroup(std::string readChromosome, long readCoordinate)
{
	// if inside...
	// recalculate mean
	if(this->chromosome == readChromosome) {
		if(abs(readCoordinate - this->coordinate) < 1)
		{

			this->coordinateTotal += readCoordinate;
			this->count += 1;
			this->coordinate = this->coordinateTotal/this->count;


			//printf("breakpoint: %s - %lu, coordinate total: %lu, read count: %d\n", this->chromosome.c_str(), this->coordinate, this->coordinateTotal, this->count);

			return true;
		}
	}
	return false;
}

void BreakPointNode::AddSupportingRead(ReadSequence &readSequence)
{
	ReadSequence *pReadSequence;
	pReadSequence = &readSequence;
	this->readsSupporting.push_back(pReadSequence);
}

bool BreakPointNode::operator ==(const BreakPointNode &b)
{
	if(this->chromosome == b.chromosome && this->coordinate == b.coordinate) return true;
	return false;

}


BreakPointNode::~BreakPointNode() {
	// TODO Auto-generated destructor stub
}
