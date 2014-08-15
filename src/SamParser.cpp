/*
 * SamParser.cpp
 *
 *  Created on: Jun 20, 2014
 *      Author: galaxy
 */

#include "SamParser.h"

SamParser::SamParser(	float min_standard_deviation,
						uint32_t min_softclip_size,
						uint32_t min_mapping_quality) {

	_InitDefaultValues();
	SetParameters(	min_standard_deviation,
					min_softclip_size,
					min_mapping_quality);
	}

SamParser::SamParser()
{
	_InitDefaultValues();
}

void SamParser::SetParameters(	float min_standard_deviation,
								uint32_t min_softclip_size,
								uint32_t min_mapping_quality)
{
	this->_min_standard_deviation = min_standard_deviation;
	this->_min_mapping_quality = min_mapping_quality;
	this->_min_softclip_size = min_softclip_size;
}

void SamParser::Parse(std::string sam_file_path, std::string reference_file_path)
{

	this->_reference_file_path = reference_file_path;
	this->_sam_file_path = sam_file_path;

	samfile_t *fp_in = NULL;
	bam_header_t *bh = NULL;
	bam1_t *b = NULL;


	_GetAlignmentStatistics();

	fp_in = samopen(sam_file_path.c_str(), "rb", 0);

	if(NULL == fp_in)
	{
		printf("Could not open file\n");
	}

	b = bam_init1();
	bh = fp_in->header;

	float upper_limit = this->mean
			+ (this->_min_standard_deviation * this->deviation);
	float lower_limit = this->mean
			- (this->_min_standard_deviation * this->deviation);

	while(samread(fp_in, b) > 0)
	{

		ReadSequence rs = _GetReadFromAlignment(b, bh);

		ReadSequence *rs_ptr = new ReadSequence(_GetReadFromAlignment(b, bh));

		if(rs.qual_score > this->_min_mapping_quality)
		{
			if(rs.HasSoftClips())
			{
				this->target_reads_ptr[rs.readID].push_back(rs);
				this->target_reads_ptr2.Add(rs_ptr);
			}
			if(rs.HasDeletions())
			{
				this->deletion_reads_ptr.push_back(rs_ptr);
				this->target_reads_ptr2.Add(rs_ptr);
			}

			// if it is not a softclip or unmapped check to see if it is discordant
			if(rs.insert_size != 0)
			{

				if ((rs.insert_size >= (upper_limit)) || (rs.insert_size <= lower_limit))
				{
					rs_ptr->SetDiscordantRead();
					bool read_is_mate = false;
					std::map<std::string,int> discordant_hash;

					ReadSequence *prev_rs = this->discordant_reads_ptr.Get(rs_ptr);

					if(prev_rs == NULL)
					{
						this->discordant_reads_ptr.Add(rs_ptr);
					}
					else
					{
						if(prev_rs->IsMateOf(*rs_ptr))
							prev_rs->SetReadMate(rs_ptr);
					}
				}
			}
		}

		bam_destroy1(b);
		b = bam_init1();
	}


	//TODO: Replace this by iterating over the discordant_reads_ptr
	// in the future for xml output
	for(unsigned int x = 0; x < this->discordant_reads_ptr.size; x++ )
	{
		if(this->discordant_reads_ptr.GetMap()[x]->size() > 0)
		{
			for(std::vector< HardSearch::HashEntry* >::iterator h_iter = this->discordant_reads_ptr.GetMap()[x]->begin(); h_iter != this->discordant_reads_ptr.GetMap()[x]->end(); h_iter++)
			{
				this->discordant_reads.push_back(*(*h_iter)->GetValue());
			}
		}
	}

	bam_destroy1(b);
	samclose(fp_in);


	printf("Number of Discordant Reads: %lu\n", this->discordant_reads.size());
	printf("Retrieving Read IDs\n");

	_GetMateReadID(sam_file_path);
	_ParseReads();

	printf("the number of target reads found %lu\n",
			this->target_reads_ptr.size());

}

void SamParser::ParseForSoftClips(	ReadSequence rs,
									pthread_mutex_t *mut,
									uint64_t &left_counter,
									uint64_t &right_counter,
									int min_softclip_size,
									std::string reference_fp,
									std::vector<SoftClipRead *> *softclip_reads_ptr)
{
	if(rs.HasSoftClips())
	{

		SoftClipRead softclip = SoftClipRead(rs);
		SoftClipRead *scr_ptr = new SoftClipRead(rs);

		if((softclip.IsLeftClip()) || (softclip.IsRightClip()))
		{
//			FixSoftClip(softclip);
			softclip.FixSoftClip(softclip.GetReference(reference_fp));

			if(softclip.SizeIsGreaterThan(min_softclip_size))
			{
				if (softclip.IsLeftClip())
				{
					pthread_mutex_lock(mut);
					left_counter += 1;
					pthread_mutex_unlock(mut);
				}

				if (softclip.IsRightClip())
				{
					pthread_mutex_lock(mut);
					right_counter += 1;
					pthread_mutex_unlock(mut);
				}

				softclip.read_type &= ReadSequence::READ_SOFTCLIP;

				//this might actually go away...
				pthread_mutex_lock(mut);
				softclip_reads_ptr->push_back(scr_ptr);
			//	this->softclip_reads.push_back(softclip);
				pthread_mutex_unlock(mut);
			}

		}

	}

}

void SamParser::ParseForDeletions(ReadSequence rs, pthread_mutex_t *mut)
{
	if(rs.HasDeletions())
	{
		ReadSequence *deletion_read_ptr = new ReadSequence(rs);
		pthread_mutex_lock(mut);
		this->deletion_reads_ptr.push_back(deletion_read_ptr);
		pthread_mutex_unlock(mut);
	}
}


void SamParser::ParseForUnmappedReads(ReadSequence rs)
{
	if(rs.bitflag & ReadSequence::READ_UNMAPPED)
	{
		ReadSequence *unmapped_read_ptr = new ReadSequence(rs);
		this->unmapped_reads_ptr.push_back(unmapped_read_ptr);
	}
}

void SamParser::PrintNumberOfSoftClips()
{
	printf("Total Softclipped reads: %lu\n", this->softclip_reads_ptr.size());
	printf("Left softclipped reads: %lu\n", this->left_softclip_reads);
	printf("Right softclipped reads: %lu\n", this->right_softclip_reads);
}

void SamParser::PrintAlignmentStatistics()
{
	printf(	"the lower bounds is: %f"
			"the upper bounds is: %f\n", this->_lower_limit, this->_upper_limit);
}

void SamParser::FixSoftClip(SoftClipRead rs)
{
	rs.FixSoftClip(rs.GetReference(this->_reference_file_path));
}

void SamParser::_GetMateReadID(std::string sam_file_path, HardSearch::HashMap reads)
{
	samfile_t *fp_in = NULL;
	bam1_t *b = NULL;
	bam_header_t *bh = NULL;

	fp_in = samopen(this->_sam_file_path.c_str(), "rb", 0);

	if(NULL == fp_in)
	{
		printf("Could not open file\n");
		return;
	}

	bh = fp_in->header;
	b = bam_init1();

	int counter = 0;

	while(samread(fp_in, b) > 0)
	{
		ReadSequence *rs = new ReadSequence(_GetReadFromAlignment(b, bh));

		ReadSequence *prev_rs = reads.Get(rs);

		if(prev_rs != NULL)
		{
			if(prev_rs->IsMateOf(*rs))
			{
				counter++;
				prev_rs->SetReadMate(rs);
				printf("Set Mate for %d of %d\r", counter, reads.item_count );
			}
		}

		bam_destroy1(b);
		b = bam_init1();
	}
	bam_destroy1(b);
	samclose(fp_in);

	for(unsigned int x = 0; x < this->unmapped_reads_ptr2.size; x++ )
	{
		if(this->discordant_reads_ptr.GetMap()[x]->size() > 0)
		{
			for(std::vector< HardSearch::HashEntry* >::iterator h_iter = this->unmapped_reads_ptr2.GetMap()[x]->begin(); h_iter != this->unmapped_reads_ptr2.GetMap()[x]->end(); h_iter++)
			{
				this->unmapped_reads_ptr.push_back((*h_iter)->GetValue());
			}
		}
	}

}

void SamParser::_GetMateReadID(std::string sam_file_path)
{
	samfile_t *fp_in = NULL;
	bam1_t *b = NULL;
	bam_header_t *bh = NULL;

	fp_in = samopen(sam_file_path.c_str(), "rb", 0);

	if(NULL == fp_in)
	{
		printf("Could not open file\n");
		return;
	}

	bh = fp_in->header;
	b = bam_init1();

	while(samread(fp_in, b) > 0)
	{
		ReadSequence rs = _GetReadFromAlignment(b, bh);

		ReadSequence *rs_ptr = new ReadSequence(rs);

		ReadSequence *prev_rs = this->target_reads_ptr2.Get(rs_ptr);
		if(prev_rs != NULL)
		{
			if(prev_rs->IsMateOf(rs))
			{
				prev_rs->SetReadMate(rs_ptr);
			}
		}

		bam_destroy1(b);
		b = bam_init1();
	}
	bam_destroy1(b);
	samclose(fp_in);
}

ReadSequence SamParser::_GetReadFromAlignment( const bam1_t *b, bam_header_t *bh)
{
	std::string read_id = GetReadID(b);
	uint32_t bitflag = GetBitflag(b);
	std::string chromosome = GetChromosome(b, bh);
	uint32_t start_coordinate = GetStartCoordinate(b);
	std::string cigar = GetCigar(b);
	uint32_t quality_score = GetQualityScore(b);
	std::string sequence = GetReadSequence(b);
	uint32_t isize = GetInsertSize(b);
	int32_t tid = b->core.tid;

	ReadSequence rs;
	rs.SetReadValues(	read_id,
						bitflag,
						tid,
						chromosome,
						start_coordinate,
						cigar,
						quality_score,
						sequence,
						isize);
	return rs;
}

void SamParser::_InitDefaultValues()
{
	this->_sam_file_path = "";
	this->total_softclip_reads = 0;
	this->left_softclip_reads = 0;
	this->right_softclip_reads = 0;
	this->mean = 0;
	this->deviation = 0;
	this->_lower_limit = 0;
	this->_upper_limit = 0;
	this->_min_standard_deviation = 3.0;
	this->_min_softclip_size = 5;
	this->_min_mapping_quality = 20;

}

void *SamParser::_Parse(void *arguments)
{
	struct job_list *params;
	params = (struct job_list *) arguments;

	pthread_mutex_lock(params->mut);
	int unsigned total_jobs = params->jobs.size() - 1;
	pthread_mutex_unlock(params->mut);

	while(total_jobs != 0)
	{

		pthread_mutex_lock(params->mut);
		int unsigned cur_job = params->completed_jobs;
		params->completed_jobs++;
		pthread_mutex_unlock(params->mut);

		if(cur_job > total_jobs)
			break;


		struct parse_args *job_params;
		job_params = params->jobs[cur_job];

		printf("Fixing SoftClip Read %d of %d\r", job_params->counter, total_jobs);



		ReadSequence rs = *job_params->first_read;

		ParseForSoftClips(	rs, params->mut,
							*params->left_softclip_reads,
							*params->right_softclip_reads,
							job_params->min_softclip_size,
							job_params->reference_fp,
							job_params->softclip_reads_ptr);

		if(rs.mate != NULL)
		{
			ParseForSoftClips(	*job_params->second_read,
								params->mut,
								*params->left_softclip_reads,
								*params->right_softclip_reads,
								job_params->min_softclip_size,
								job_params->reference_fp,
								job_params->softclip_reads_ptr);
		}



	}

	return 0;
}

void SamParser::_ParseReads()
{
	int y = 0;
	int unsigned max_threads = 36;
	pthread_t threads[max_threads];
	// initialize thread output ptr vector

	std::vector< struct parse_args *> jobs;
	// initialize thread status array
	std::vector<bool> thread_status;
	for(int unsigned tc = 0; tc < max_threads; tc++)
	{
		thread_status.push_back(true);
	}

	pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;

	// This loops through the hashmap and sets up all of the jobs
	// in the queue
	for(unsigned int x = 0; x < this->target_reads_ptr2.size; x++ )
	{
		if(this->target_reads_ptr2.GetMap()[x]->size() > 0)
		{
			std::vector< HardSearch::HashEntry* >::iterator h_iter;
			std::vector< HardSearch::HashEntry* >::iterator iter_start;
			std::vector< HardSearch::HashEntry* >::iterator iter_end;

			iter_start = this->target_reads_ptr2.GetMap()[x]->begin();
			iter_end = this->target_reads_ptr2.GetMap()[x]->end();

			for(h_iter = iter_start; h_iter != iter_end; h_iter++)
			{

				//prepare and add job to queue;

				y += 1;
				ReadSequence *first_read = (*h_iter)->GetValue();
				ReadSequence *second_read = (*h_iter)->GetValue()->mate;


				// Sample1.bam file.... is the information still relavent...

				struct parse_args *args = new parse_args();
				args->mut = &mut;
				args->first_read = first_read;
				args->second_read = second_read;
				args->reference_fp = this->_reference_file_path;
				args->softclip_reads_ptr = &this->softclip_reads_ptr;
				args->min_softclip_size = this->_min_softclip_size;
				args->counter = y;

				jobs.push_back(args);

			}
		}
	}

	struct job_list *jl = new job_list();
	jl->completed_jobs = 0;
	jl->jobs = jobs;
	jl->mut = &mut;
	jl->left_softclip_reads = &this->left_softclip_reads;
	jl->right_softclip_reads = &this->right_softclip_reads;

	for(int unsigned t = 0; t < max_threads; t++)
	{
		pthread_create(&threads[t], NULL, &_Parse, jl);
	}
	for(int unsigned t = 0; t< max_threads; t++)
	{
		pthread_join(threads[t], NULL);
	}

	delete jl;
}


bool SamParser::OutputSoftClipReads(std::string output_file_name)
{
	std::ofstream output_file;
	output_file.open(output_file_name.c_str());

	if(output_file)
	{
		for(std::vector<SoftClipRead*>::iterator scReadIter = this->softclip_reads_ptr.begin(); scReadIter != this->softclip_reads_ptr.end(); scReadIter++)
		{
			if((*scReadIter)->IsEndClip() && (*scReadIter)->SizeIsGreaterThan(this->_min_softclip_size))
			{
				std::string softclip = (*scReadIter)->GetSoftClipBedFormat();
				output_file.write(softclip.c_str(), strlen(softclip.c_str()));
			}
		}

		for(std::vector<ReadSequence>::iterator drpReadIter = this->discordant_reads.begin(); drpReadIter != this->discordant_reads.end(); drpReadIter++)
		{
			if(drpReadIter->discordant_read)
			{
				std::string softclip = drpReadIter->chromosome +
						"\t" + HardSearch::NumToString<long>(drpReadIter->start_pos) +
						"\t" + HardSearch::NumToString<long>(drpReadIter->start_pos + 1) +
						"\t" + "unknown|left" + "|" + HardSearch::NumToString<int32_t>(drpReadIter->tid)  + "\n";
				output_file.write(softclip.c_str(), strlen(softclip.c_str()));
			}
		}
		output_file.flush();
		return true;
		output_file.close();
	}
	return false;
}

std::string SamParser::GetReadID(const bam1_t *b)
{
	char *name = bam1_qname(b);
	return std::string(name, strlen(name));
}

uint32_t SamParser::GetBitflag(const bam1_t *b)
{
	uint32_t bitflag = b->core.flag;
	return bitflag;
}

std::string SamParser::GetChromosome(const bam1_t *b, bam_header_t *bh)
{
	if(b->core.tid >= 0)
		return std::string(bh->target_name[b->core.tid], strlen(bh->target_name[b->core.tid]));
	else
		return std::string("1/aefhimopqrstv[]*",1);
}

int32_t SamParser::GetStartCoordinate(const bam1_t *b)
{
	return b->core.pos + 1;
}

std::string SamParser::GetCigar(const bam1_t *b)
{
	std::string cigar_string;
	uint32_t *cigar = bam1_cigar(b);

	for(uint32_t n=0; n<(b->core.n_cigar);n++)
	{
		int32_t cop = cigar[n] & BAM_CIGAR_MASK;
		int32_t cl = cigar[n] >> BAM_CIGAR_SHIFT;

		cigar_string += HardSearch::NumToString<uint32_t>(cl);

		switch (cop)
		{
		case BAM_CMATCH:
			cigar_string += "M";
			break;
		case BAM_CSOFT_CLIP:
			cigar_string += "S";
			break;
		case BAM_CINS:
			cigar_string += "I";
			break;
		case BAM_CDEL:
			cigar_string += "D";
			break;
		case BAM_CREF_SKIP:
			cigar_string += "N";
			break;
		case BAM_CHARD_CLIP:
			cigar_string += "H";
			break;
		case BAM_CPAD:
			cigar_string += "=";
			break;
		case BAM_CDIFF:
			cigar_string += "X";
			break;
		default:
			cigar_string += "?";
			break;
		}
	}

	if (b->core.n_cigar == 0)
		cigar_string += "*"; // if the read is unmapped and no cigar is present

	return cigar_string;
}

uint32_t SamParser::GetQualityScore(const bam1_t *b)
{
	return b->core.qual;
}

std::string SamParser::GetReadSequence(const bam1_t *b)
{
	std::string sequence;

	char *qseq = (char *) malloc(b->core.l_qseq+1);
	uint8_t *s = bam1_seq(b);
	for(int n=0; n<(b->core.l_qseq); n++)
	{
		char v = bam1_seqi(s,n);
		qseq[n] = bam_nt16_rev_table[v];
	}

	sequence = std::string(qseq, b->core.l_qseq);
	delete qseq;
	return sequence;
}

uint32_t SamParser::GetInsertSize(const bam1_t *b)
{
	return std::abs(b->core.isize);
}

void SamParser::_GetAlignmentStatistics()
{
	samfile_t *fp_in = NULL;
	bam1_t *b = NULL;

	uint64_t	total_size = 0,
				sum_sq = 0;
	uint32_t	reads = 0,
				insert_size = 0;
	fp_in = samopen(this->_sam_file_path.c_str(), "rb", 0);

	if(NULL == fp_in)
	{
		printf("Could not open file\n");
		return;
	}

	b = bam_init1();

	while(samread(fp_in, b) > 0)
	{

		insert_size = std::abs(b->core.isize);

		if(insert_size > 10000)
		{
			continue;
		}

		total_size += insert_size;
		sum_sq += insert_size * insert_size;
		reads++;
		bam_destroy1(b);
		b = bam_init1();
	}

	this->mean = total_size / reads;
	this->deviation = std::sqrt(sum_sq/reads - std::pow(this->mean, 2));

	this->_upper_limit = this->mean
			+ (this->_min_standard_deviation * this->deviation);
	this->_lower_limit = this->mean
			- (this->_min_standard_deviation * this->deviation);

	printf("size: %lu / reads: %d  +/- deviation %d\n", total_size, reads, this->deviation);

	bam_destroy1(b);
	samclose(fp_in);

}

void SamParser::RetrieveReadsForRegion(	int tid,
										int beg,
										int end,
										std::string direction)
{
	bam_fetch_f funct = SamParser::_UnmappedMatesCallBack;

	samfile_t *fp_in = NULL;
	fp_in = samopen(this->_sam_file_path.c_str(), "rb", 0);

	struct args *data = new args();
	data->bh = fp_in->header;
	data->direction = direction;
	data->quality = this->_min_mapping_quality;
//	data->unmapped_reads_ptr = &this->unmapped_reads_ptr;
//	data->unmapped_reads_set = &this->unmapped_reads_set;
	data->unmapped_reads_ptr2 = &this->unmapped_reads_ptr2;

	BGZF *bam_file = NULL;
	bam_index_t *bam_index;

	bam_file = bam_open(this->_sam_file_path.c_str(), "rb");

	bam_index = bam_index_load(this->_sam_file_path.c_str());

	bam_fetch(bam_file, bam_index, tid, beg, end, data, (*funct));

	samclose(fp_in);
	bam_close(bam_file);

}


int SamParser::_UnmappedMatesCallBack(const bam1_t *b, void *data)
{
	struct args *param;
	param = (struct args *) data;
	bam_header_t *bh = param->bh;
	ReadSequence *rs = new ReadSequence(SamParser::_GetReadFromAlignment(b, bh));

	if(rs->qual_score >= param->quality)
	{
		if(!(rs->bitflag & ReadSequence::FAILED_PV_QC) &&
				!(rs->bitflag & ReadSequence::PCR_DUP) &&
				(rs->bitflag & ReadSequence::READ_PAIRED))
		{
			if(param->direction == "left")
			{
				if((rs->bitflag & ReadSequence::MATE_UNMAPPED) &&
						(rs->bitflag & ReadSequence::READ_REVERSE_STRAND))
				{
						param->unmapped_reads_ptr2->Add(rs);
				}
			}
			else
			{
				if((rs->bitflag & ReadSequence::MATE_UNMAPPED) &&
						!(rs->bitflag & ReadSequence::READ_REVERSE_STRAND))
				{
					param->unmapped_reads_ptr2->Add(rs);
				}
			}
		}
	}


	return 0;
}


SamParser::~SamParser() {
	// TODO Auto-generated destructor stub
	for(std::vector<ReadSequence *>::iterator it = this->deletion_reads_ptr.begin(); it != this->deletion_reads_ptr.end(); it++)
	{
		delete *it;
	}

	for(std::vector<SoftClipRead *>::iterator it = this->softclip_reads_ptr.begin(); it != this->softclip_reads_ptr.end(); it++)
	{
		delete *it;
	}

}
