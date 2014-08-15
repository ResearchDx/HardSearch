#include "BedUtils.h"

BedUtils::BedUtils()
{
}

void BedUtils::LoadBedFile(std::string file_path)
{
	std::ifstream bed_file;
	bed_file.open(file_path.c_str());

	std::string raw_data;
	std::vector<std::string> data;
	while (std::getline(bed_file, raw_data))
	{
		GenomicRegion gr;
		data = SplitByDelimiter((char*) raw_data.c_str(), "\t");
		gr.chromosome = data[0];
		gr.start = atol(data[1].c_str());
		gr.end = atol(data[2].c_str());	
		this->bed_regions.push_back(gr);
	}
}

bool BedUtils::InBedRegion(std::string chromosome, long unsigned start_coordinate, long unsigned end_coordinate)
{
	// convert the parameters into a genomic region
	GenomicRegion query;
	query.chromosome = chromosome;
	query.start = start_coordinate;
	query.end = end_coordinate;

	for (std::vector<GenomicRegion>::iterator region_iterator = this->bed_regions.begin(); region_iterator != this->bed_regions.end(); region_iterator++)
	{
		if(IsIntersecting(query, *region_iterator))
			return true;
	}
	return false;
}

bool BedUtils::IsIntersecting(GenomicRegion query, GenomicRegion reference)
{
	if (query.chromosome == reference.chromosome)
	{
		// if the query starts before the reference and ends in the middle
		// of the reference
		if (query.start <= reference.start && query.end >= reference.start)
			return true;

		// if the query starts in the middle of the reference and ends after
		// the reference
		if (query.end >= reference.end && query.start <= reference.end)
			return true;

		// if the query starts in the middle of the reference and ends in the
		// middle of the reference
		if (query.start >= reference.start && query.end <= reference.end)
			return true;
	}
	return false;
}

void BedUtils::PrintBedRegions()
{
	for (std::vector<GenomicRegion>::iterator region_iterator = this->bed_regions.begin(); region_iterator != this->bed_regions.end(); region_iterator++)
	{
		printf("%s\n", region_iterator->chromosome.c_str());
	}
}

std::vector<std::string> BedUtils::SplitByDelimiter(char *input, char *delimiter)
{
	std::vector<std::string> output;
	char* input_array = strtok(input, delimiter);
	while (input_array)
	{
		output.push_back(input_array);
		input_array = strtok(NULL, delimiter);
	}
	return output;
}

