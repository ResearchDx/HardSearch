#include <string.h>
#include <iostream>
#include <fstream>

#include <vector>
#include <stdio.h>
#include <stdlib.h>


#ifndef BED_UTILS
#define BED_UTILS

struct GenomicRegion
{
	std::string chromosome;
	long unsigned start;
	long unsigned end;
};

class BedUtils
{
public:

	std::vector<GenomicRegion> bed_regions;
	
	BedUtils();

	// Load the regions into a list of GenomicRegions
	void LoadBedFile(std::string file_path);
	void PrintBedRegions();

	// Given the components of a GenomicRegion determine if 
	// the query intersects with the regions defined in the bed file
	bool InBedRegion(std::string chromosome, long unsigned start_coordinate, long unsigned end_coordinate);

	// Given a query and reference GenomicRegion determine if 
	// the query intersects with the reference
	bool IsIntersecting(GenomicRegion query, GenomicRegion reference);

	static std::vector<std::string>SplitByDelimiter(char *input, char *delimiter);

};

#endif
