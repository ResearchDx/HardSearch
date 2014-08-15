#include "main.h"


//#define DEBUG
//#define VERBOSE
#define CMD_PARSER
//#define VERBOSE_FINDBREAKPOINT
int unsigned UNMATCHEDBASES = 25; // minimum number of unaligned bases remaining after initial blast
int unsigned MIN_SOFTCLIPS = 2; // minimum number of softclips per target region
int unsigned MIN_SOFTCLIP_SIZE = 5; // minimum size of softclip to include
int unsigned MIN_MAPPING_QUALITY = 20; // minimum read quality of aligment
int unsigned THREADS = 1; // number of threads to run
float MIN_STANDARD_DEVIATION = 2.5;

typedef rapidxml::xml_node<>* XmlNode;




int main(int argc, char *argv[]) {

#ifdef CMD_PARSER
	if(SizeOfParameterList(argv) < 4)
	{
		DisplayHelp();
		return 0;
	}
	char cCurrentPath[FILENAME_MAX];
	getcwd(cCurrentPath, sizeof(cCurrentPath));
	boost::filesystem::path cwd(cCurrentPath);

	boost::program_options::options_description desc("Allowed Options");
	_AddProgramParameters(desc);

	boost::program_options::variables_map vm;
	try {
		boost::program_options::store(
				boost::program_options::parse_command_line(argc, argv, desc),
				vm);
		if (vm.count("help")) {
			DisplayHelp();
			return 0;
		}
	} catch (boost::program_options::error &e) {
		printf("ERROR: %s \n", e.what());
		return 1;
	}

	std::string path_to_sam;
	std::string path_to_reference;
	std::string output_directory;
	std::string path_to_bed;
	std::string path_to_blast;
	std::string path_to_blast_reference;
	std::string path_to_blasthelper;


	std::string exec_path = get_selfpath();
	path_to_sam = CheckFilePath(vm["input"], cwd, "input");
	path_to_reference = CheckFilePath(vm["reference"], cwd, "reference");
	path_to_bed = CheckFilePath(vm["bed"], cwd, "bed");
	path_to_blast = CheckFilePath(vm["blast"], cwd, "blast");
	path_to_blast_reference = CheckFilePath(vm["blast_reference"], cwd, "blast_reference", true);
	path_to_blasthelper = GetDirectory((char *)exec_path.c_str()) + "/BlastHelper";


	if(path_to_sam == "")
	{
		printf("\tMissing argument -i Input bam file\n");
		DisplayHelp();
		return 0;
	}
	if(path_to_reference == "")
	{
		printf("\tMissing argument -r Path to reference file\n");
		DisplayHelp();
		return 0;
	}
	if(path_to_bed == "")
	{
		printf("\tMissing argument -b Path to bed file\n");
		DisplayHelp();
		return 0;
	}
	if(path_to_blast == "")
	{
		printf("\tMissing argument -l Path to blastn\n");
		DisplayHelp();
		return 0;
	}
	if(path_to_blast_reference == "")
	{
		printf("\tMissing argument --blast_reference Path to blast reference\n");
		DisplayHelp();
		return 0;
	}

	_ParseProgramParameters(vm);



	output_directory = cwd.c_str();
#endif

	std::time_t start_time;
	start_time = std::time(&start_time);

	// retrieve softclip reads from the BAM alignment file
	// the same way that softsearch does
	// TODO: identify discordant read pairs from the bam file
	printf("Retrieving Soft Clip Reads from %s\n", path_to_sam.c_str());
	SamParser *sp = new SamParser;
	sp->SetParameters(MIN_STANDARD_DEVIATION, MIN_SOFTCLIP_SIZE,
			MIN_MAPPING_QUALITY);

	sp->Parse(path_to_sam, path_to_reference);

	if (!sp->OutputSoftClipReads(
			output_directory + "/hardsearch.direction.bed")) {
		printf("Error Writing to file %s \n",
				(output_directory + "/hardsearch.direction.bed").c_str());
		return 0;
	}

	MergeSoftClips(output_directory);

	std::vector<std::string> sam_read_IDs;

	// get the sequences with unmapped mates from the regions with soft clipped reads
	// above the minimum

	//put target_regions here
	std::vector<std::string> target_regions;

	// retrieve target regions from the /output/directory/hardsearch.fixed.bed file
	// if the file does not exists throw fileinput exception
	// TODO: catch fileinput exception when the file does not exist
	// in the event that no softclip events were found
	target_regions = GetTargetRegions(
			output_directory + "/hardsearch.fixed.bed", MIN_SOFTCLIPS);

	if (!_RetrieveUnmappedMateReads(path_to_sam, output_directory,
			target_regions, MIN_MAPPING_QUALITY, sp)) {
		return 0;
	}

	sp->_GetMateReadID(sp->GetSamFilePath(), sp->unmapped_reads_ptr2);

	printf("Number of Unmapped Reads Found: %lu\n",
			sp->unmapped_reads_ptr.size());

	// Initialize the Output Variables outside of scope
	HardSearch::RapidXMLWrapper results_xml;
	rapidxml::xml_node<>* node_events = results_xml.AddNode("Events");
	(*results_xml.doc).append_node(node_events);

	printf("Gathering Reads IDs for reads with Unmapped Mates\n");
	if (sp->unmapped_reads_ptr.size() != 0
				|| sp->softclip_reads_ptr.size() != 0)
	{

		std::vector<ReadSequence> readList;
		std::vector<ReadSequence *>::iterator rs_iter;
		std::vector<SoftClipRead *>::iterator sc_iter;

		for (rs_iter = sp->unmapped_reads_ptr.begin();
				rs_iter != sp->unmapped_reads_ptr.end(); rs_iter++)
		{
			ReadSequence rs = ReadSequence(**rs_iter);
			readList.push_back(rs);
		}

		for (sc_iter = sp->softclip_reads_ptr.begin();
				sc_iter != sp->softclip_reads_ptr.end(); sc_iter++)
		{
			ReadSequence rs = ReadSequence(**sc_iter);
			readList.push_back(rs);

		}

		printf("Mapping Reads with BLAST\n");
		// Blast Initial Half of Read

		pthread_t threads[THREADS];
		struct job_queue *jl = new struct job_queue(); // data structure that contains the job
											// queue with all the jobs to be done

		std::vector<struct arg_struct*> j;// job queue containing information about
											// each job
		pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
		for (int unsigned i = 0; i < readList.size(); i++) {

			std::string blast_xml_string;
			// if the read is unmapped from the bitflag set by the aligner
			// or the read contains a softclip

			ReadSequence rs = readList[i];
			std::string read_sequence;
			if ((rs.IsUnmapped()) || (rs.HasSoftClips())) {
				read_sequence = ConstructFastaRead(readList[i].readID,
						readList[i].sequence);
			} else if ((rs.mate->IsUnmapped() || (rs.mate->HasSoftClips()))) {
				read_sequence = ConstructFastaRead(readList[i].mate->readID,
						readList[i].mate->sequence);
			} else {
				continue;
			}

			struct arg_struct *args = new arg_struct();

			args->read_sequence = read_sequence;
			args->read_list = &readList;
			args->mate_mode = false;			args->mut = &mut;
			args->counter = i;

			j.push_back(args);
		}

		jl->jobs = &j;
		jl->completed_jobs = 0;
		jl->mut = &mut;
		jl->path_to_blasthelper = path_to_blasthelper;
		jl->path_to_blast = path_to_blast;
		jl->path_to_blast_reference = path_to_blast_reference;



		for(int unsigned t = 0; t < THREADS; t++)
		{
			pthread_create(&threads[t], NULL, &_BlastReads, jl);
		}
		for(int unsigned t = 0; t < THREADS; t++)
		{
			pthread_join(threads[t], NULL);
		}

		struct job_queue *jl2 = new struct job_queue();
		std::vector<struct arg_struct*> j2;

		printf("Re-Mapping Unaligned Remainder of Reads with BLAST\n");
		//Blast Remaining Half of Read

		for (int unsigned i = 0; i < readList.size(); i++) {
			// if the readis unmapped and the bitflag is 4?
			std::string blast_xml_string;

			ReadSequence rs = readList[i];

			std::string read_sequence;
			if ((rs.IsUnmapped()) || (rs.HasSoftClips())) {
				read_sequence = _ConstructFastaFileMates(readList[i]);
			} else if ((rs.mate->IsUnmapped() || (rs.mate->HasSoftClips()))) {
				read_sequence = _ConstructFastaFileMates(*readList[i].mate);
			} else {
				continue;
			}

			if (strlen(read_sequence.c_str()) != 0) {
				struct arg_struct *args = new arg_struct();

				args->read_sequence = read_sequence;
				args->read_list = &readList;
				args->mate_mode = true;
				args->mut = &mut;
				args->counter = i;

				j2.push_back(args);
			}
		}

		jl2->jobs = &j2;
		jl2->completed_jobs = 0;
		jl2->mut = &mut;
		jl2->path_to_blasthelper = path_to_blasthelper;
		jl2->path_to_blast = path_to_blast;
		jl2->path_to_blast_reference = path_to_blast_reference;


		for(int unsigned t = 0; t < THREADS; t++)
		{
			pthread_create(&threads[t], NULL, &_BlastReads, jl2);
		}
		for(int unsigned t = 0; t < THREADS; t++)
		{
			pthread_join(threads[t], NULL);
		}

		rapidxml::xml_node<>* node_softclip_reads_event = results_xml.AddNode("Softclip_Reads");
		node_events->append_node(node_softclip_reads_event);


		// output softclip reads realigned by blast

		for (std::vector<ReadSequence>::iterator sc_iter =
				readList.begin();
				sc_iter != readList.end(); sc_iter++)
		{

			if(sc_iter->blastAlignments.size() !=0 || sc_iter->mate->blastAlignments.size() != 0)
			{
				XmlNode node_softclip_read = results_xml.AddNode("Softclip_Read_Pair");
				XmlNode n_sc_1;
				XmlNode n_sc_2;

				if(sc_iter->blastAlignments.size() == 0)
				{
					n_sc_1 = results_xml.AddReadNode(*sc_iter);
					if((*sc_iter).mate != NULL)
					{
						n_sc_2 = results_xml.AddReadNode(*(*sc_iter).mate);
						results_xml.AddReadAlignmentToNode((*sc_iter).mate->blastAlignments[0], n_sc_2);
					}

				}
				else
				{
					n_sc_1 = results_xml.AddReadNode(*sc_iter->mate);
					n_sc_2 = results_xml.AddReadNode((*sc_iter));
					results_xml.AddReadAlignmentToNode((*sc_iter).blastAlignments[0], n_sc_2);
				}
				node_softclip_read->append_node(n_sc_1);
				node_softclip_read->append_node(n_sc_2);
				node_softclip_reads_event->append_node(node_softclip_read);
			}
		}

		// output deletions

		XmlNode node_deletions = results_xml.AddNode("Deletions");
		node_events->append_node(node_deletions);

		for (std::vector<ReadSequence *>::iterator rs_iter =
				sp->deletion_reads_ptr.begin();
				rs_iter != sp->deletion_reads_ptr.end(); rs_iter++)
		{
			XmlNode node_deletion_pair = results_xml.AddNode("Deletion_Read_Pair");
			node_deletions->append_node(node_deletion_pair);

			XmlNode n_dr_1 = results_xml.AddReadNode(**rs_iter);
			node_deletion_pair->append_node(n_dr_1);



			if((**rs_iter).GetDeletionStart() != NULL)
			{
				int32_t del_start = *((**rs_iter).GetDeletionStart());
				int32_t del_end = *((**rs_iter).GetDeletionEnd());

				XmlNode deletion_coordinates = results_xml.AddNode("Deletion_Region");
				std::string deletion_region = (**rs_iter).chromosome + ":"
						+ HardSearch::NumToString<int32_t>(del_start) + "-"
						+ HardSearch::NumToString<int32_t>(del_end);

				HardSearch::RapidXMLWrapper::SetNodeValue(deletion_coordinates, deletion_region);
				n_dr_1->append_node(deletion_coordinates);
			}


			if((**rs_iter).mate != NULL)
			{
				XmlNode n_dr_2 = results_xml.AddReadNode(*(**rs_iter).mate);
				node_deletion_pair->append_node(n_dr_2);

				if((**rs_iter).mate->GetDeletionStart() != NULL)
				{
					int32_t del_start = *((**rs_iter).mate->GetDeletionStart());
					int32_t del_end = *((**rs_iter).mate->GetDeletionEnd());

					XmlNode deletion_coordinates = results_xml.AddNode("Deletion_Region");
					std::string deletion_region = (**rs_iter).chromosome + ":"
							+ HardSearch::NumToString<int32_t>(del_start) + "-"
							+ HardSearch::NumToString<int32_t>(del_end);
					HardSearch::RapidXMLWrapper::SetNodeValue(deletion_coordinates, deletion_region);
					n_dr_2->append_node(deletion_coordinates);
				}
			}
		}
	}

	rapidxml::xml_node<>* node_discordant_reads = results_xml.AddNode(
			"Discordant_Reads");
	node_events->append_node(node_discordant_reads);
	rapidxml::xml_attribute<>* discordant_mean = results_xml.AddAttribute("mean", HardSearch::NumToString<int unsigned>(sp->mean));
	rapidxml::xml_attribute<>* discordant_dev = results_xml.AddAttribute("dev", HardSearch::NumToString<int unsigned>(sp->deviation));
	rapidxml::xml_attribute<>* discordant_distance = results_xml.AddAttribute("num_dev", HardSearch::NumToString<float>(MIN_STANDARD_DEVIATION, 2));

	node_discordant_reads->append_attribute(discordant_mean);
	node_discordant_reads->append_attribute(discordant_dev);
	node_discordant_reads->append_attribute(discordant_distance);


	for (std::vector<ReadSequence>::iterator sc_iter =
			sp->discordant_reads.begin(); sc_iter != sp->discordant_reads.end();
			sc_iter++) {


		XmlNode node_discordant_pair = results_xml.AddNode("Discordant_Read_Pair");
		node_discordant_reads->append_node(node_discordant_pair);

		XmlNode n_dc_1 = results_xml.AddReadNode(*sc_iter);
		node_discordant_pair->append_node(n_dc_1);

		if((*sc_iter).mate != NULL)
		{
			XmlNode n_dc_2 = results_xml.AddReadNode(*(*sc_iter).mate);
			node_discordant_pair->append_node(n_dc_2);
		}
	}
	std::string outputFileName = output_directory + "/results.xml";
	results_xml.OutputXml(outputFileName);

	std::time_t end_time;
	std::time(&end_time);
	printf("\nTime\t%lu\ts\n", end_time - start_time);
	return 0;
}


std::vector<std::string> GetTargetRegions(std::string input_file,
		int min_soft_clipped_count) {
	/* Read from the softclip fixed.bed file and grab the target regions
	 * to investigate
	 */

	std::ifstream soft_clipped_reads(input_file.c_str());
	if (!soft_clipped_reads) {
		throw std::ios_base::failure(input_file + " not found");
	}

	std::string line;
	std::vector<std::string> data;
	std::vector<std::string> reads;
	std::vector<std::string> target_region;
	while (std::getline(soft_clipped_reads, line)) {
		data = SplitByDelimiter((char *) line.c_str(), "\t");
		reads = SplitByDelimiter((char *) data[3].c_str(), ";");
		// count the number of right and left reads
		// for each event...
		int count = 0;
		std::string read_direction;
		std::string chromosome;
		std::string chr_tid;
		std::vector<std::string> chr_data = SplitByDelimiter(
				(char *) data[0].c_str(), ":");
		chromosome = chr_data[0];
		chr_tid = chr_data[1];

		for (std::vector<std::string>::iterator read_iter = reads.begin();
				read_iter != reads.end(); read_iter++) {
//			if( strlen(SplitByDelimiter((char *) (*read_iter).c_str(), "|")[0].c_str()) > MIN_SOFTCLIP_SIZE )
//			{
			count += 1;
			if (GetSoftClippedDirection(*read_iter) == 1)
				read_direction = "left";
			else
				read_direction = "right";
//			}
		}

		if (count >= min_soft_clipped_count) {
			int unsigned start;
			int unsigned end;
			if (read_direction == "left") {
				start = atoi(data[1].c_str());
				end = start + 300;
			} else {
				end = atoi(data[1].c_str());
				start = end - 300;
				if (start < 0)
					start = 0;
			}
			target_region.push_back(
					chromosome + "\t" + NumToString<int>(start) + "\t"
							+ NumToString<int>(end) + "\t" + read_direction
							+ "\t" + chr_tid);
		}
	}
	soft_clipped_reads.close();
	return target_region;
}

int GetSoftClippedDirection(std::string soft_clipped_read) {
	int bitflag = 0;
	std::vector<std::string> read_data = SplitByDelimiter(
			(char *) soft_clipped_read.c_str(), "|");
	if (read_data.size() != 2) {
		bitflag = 0;
		return bitflag;
	}
	std::string direction = read_data[1];
	if (read_data[1] == "Left")
		bitflag = 1;
	return bitflag;
}

void *_BlastReads(void *arguments) {
	struct job_queue *args;
	args = (struct job_queue *) arguments;
	struct arg_struct *params;
	uint32_t *jobs_completed;
	uint32_t current_job;
	uint32_t total_jobs;
	std::string cmd;

	pid_t e_pid;
	int pipefd[2];
	FILE* output;

	char line[256];
	int status;

	pthread_mutex_lock(args->mut);
	jobs_completed = &args->completed_jobs;
	total_jobs = args->jobs->size();
	pthread_mutex_unlock(args->mut);

	while(*jobs_completed < args->jobs->size())
	{
		pthread_mutex_lock(args->mut);
		current_job = *jobs_completed;
		++*jobs_completed;
		pthread_mutex_unlock(args->mut);

		printf("Started Blast job %u of %u\r", current_job, total_jobs);
		params = (struct arg_struct*) (*args->jobs)[current_job];

		std::string *xml_string;


		if(pipe(pipefd) == -1)
		{
			perror("Pipe Failed");
		}
		e_pid = vfork();
		if(e_pid < 0 )
		{
			perror("Fork Failed");
			continue;
		}

		else if(e_pid == (pid_t) 0)
		{
			// Child
			close(pipefd[0]);
			dup2(pipefd[1], STDOUT_FILENO);
			close(pipefd[1]);

			execle(	(*args).path_to_blasthelper.c_str(), "BlastHelper",
					(params->read_sequence).c_str(), (*args).path_to_blast.c_str(),
					(*args).path_to_blast_reference.c_str(),
					NULL, NULL);
			printf("BlastHelper process did not execute properly\n");
			_exit(1);
		}

		else
		{
			close(pipefd[1]);
			output = fdopen(pipefd[0], "r");

			std::string stream_string;

			xml_string = new std::string();

			while(fgets(line, sizeof(line), output))
			{
				stream_string = std::string(line, strlen(line));
				*xml_string += stream_string;
			}
			pclose(output);

		}

		waitpid(e_pid, NULL, 0);
		ReadSequence *rs = &(*params->read_list)[params->counter];
		std::vector<ReadSequence> *rl = &(*params->read_list);
		bool *mm = &(params->mate_mode);
		ParseBlastXML(*xml_string, *rl, *rs, *mm);

	}
	return 0;
}


// Pass by reference works with this, just need to include the iterator for the
// vector to know which one to edit...
void ParseNode(rapidxml::xml_node<> *iterationNode, ReadSequence &reads,
		bool MateMode) {
	// Check for presence of Hits to avoid segmentation faults
	if (iterationNode->first_node("Iteration_hits")->first_node("Hit") != 0) {

		rapidxml::xml_node<> *iterationHits = iterationNode->first_node(
				"Iteration_hits");

		BlastAlignment bestAlignment;
		for (rapidxml::xml_node<> *hit = iterationHits->first_node("Hit"); hit;
				hit = hit->next_sibling()) {
			std::string chromosome = hit->first_node("Hit_def")->value();
			// get the hits for each chromosome
			// and iterate over them to retrieve the highest score
			// TODO: possible store the information in a struct before merging
			// with the read sequence

			for (rapidxml::xml_node<> *hsp =
					hit->first_node("Hit_hsps")->first_node("Hsp"); hsp; hsp =
					hsp->next_sibling("Hsp")) {
				// Store the alignment with the highest bit score and keep
				// TODO fix this for loop to use the last_node -> value and first_node -> value
				if (atoi(hsp->first_node("Hsp_bit-score")->value())
						> bestAlignment.bitScore) {
					bestAlignment.bitScore = atof(
							hsp->first_node("Hsp_bit-score")->value());
					bestAlignment.alignedStart = atol(
							hsp->first_node("Hsp_hit-from")->value());
					bestAlignment.alignedEnd = atol(
							hsp->first_node("Hsp_hit-to")->value());
					bestAlignment.matchStart = atoi(
							hsp->first_node("Hsp_query-from")->value());
					bestAlignment.matchEnd = atoi(
							hsp->first_node("Hsp_query-to")->value());
					bestAlignment.gaps = atoi(
							hsp->first_node("Hsp_gaps")->value());
					bestAlignment.length = atoi(
							hsp->first_node("Hsp_align-len")->value());
					bestAlignment.alignedSequence =
							hsp->first_node("Hsp_qseq")->value();
					bestAlignment.chromosome = chromosome;

					// Set bitFlag
					if (bestAlignment.alignedEnd < bestAlignment.alignedStart)
						bestAlignment.bitFlag |= AlignmentNode::INVERTED;
					if (!MateMode)
					{
						bestAlignment.bitFlag |=
								(atoi(hsp->first_node("Hsp_query-from")->value())
										<= UNMATCHEDBASES) ?
										AlignmentNode::ALIGNED_LEFT :
										AlignmentNode::ALIGNED_RIGHT;
					}

				}

				if (hsp == hit->first_node("Hit_hsps")->last_node("Hsp"))
					break;
			}

			// if the hit is the last one then break out of the loop
			if (hit == iterationHits->last_node("Hit"))
				break;
		}

		// Save alignment information
		// if this is the first pass than MateMode should be false
		if (!MateMode) {
			AlignmentNode saveAlignment = AlignmentNode(bestAlignment);
			reads.SetAlignmentValues(saveAlignment);
		} else {
			if(reads.blastAlignments[0].bitflag & AlignmentNode::ALIGNED_LEFT)
				bestAlignment.bitFlag |= AlignmentNode::ALIGNED_RIGHT;
			else
				bestAlignment.bitFlag |= AlignmentNode::ALIGNED_LEFT;
			AlignmentNode *pSaveAlignment;
			pSaveAlignment = new AlignmentNode(bestAlignment);
			reads.blastAlignments[0].SetAlignmentMate(pSaveAlignment);
		}
		// Should have the alignment with the highest score now
		// class method that sets the read information

	}
}

static std::vector<std::string> SplitByDelimiter(char *input, char *delimiter) {
	std::vector<std::string> output;
	char* input_array = strtok(input, delimiter);
	while (input_array) {
		output.push_back(input_array);
		input_array = strtok(NULL, delimiter);
	}
	return output;
}

std::string _ConstructFastaFileMates(ReadSequence rs) {
	std::string fastaLine;
	for (int unsigned alignment = 0; alignment < rs.blastAlignments.size();
			alignment++) {
		std::string subString = std::string();

		// Get the substring of the sequence if over the
		// required number of base pairs
		// if first alignment on the left side of read
		// the mate will be from the base after the match end to the end
		if (std::strlen(rs.sequence.c_str())
				- rs.blastAlignments[alignment].match_end >= UNMATCHEDBASES) {
			subString = rs.sequence.substr(
					rs.blastAlignments[alignment].match_end,
					(std::strlen(rs.sequence.c_str())
							- rs.blastAlignments[alignment].length));
		}
		// if first alignment on the right side of read
		// the mate will be from the start to the base prior to the match start
		else if (rs.blastAlignments[alignment].match_start >= UNMATCHEDBASES) {

			subString = rs.sequence.substr(0,
					rs.blastAlignments[alignment].match_start - 1);

		}

		if (std::strlen(subString.c_str()) != 0) {

			fastaLine = ConstructFastaRead(rs.readID, subString).c_str();

			return fastaLine;
		}
	}
	fastaLine = "";
	return fastaLine;

}

std::string ConstructFastaRead(std::string readName, std::string readSequence) {
	std::string fasta;
	fasta = ">" + readName + "\\n" + readSequence;
	return fasta;
}

std::string GetDirectory(char *filepath) {
	std::vector<std::string> filePathTokens = SplitByDelimiter(filepath, "/");
	std::string result;

	for (int unsigned i = 0; i + 1 < filePathTokens.size(); i++) {
		result += "/";
		result += filePathTokens[i];
	}
	return result;
}

void ParseBlastXML(std::string xml_string, std::vector<ReadSequence> &readList,
		ReadSequence &rs, bool mateMode) {

	std::vector<char> xmlCopy(xml_string.begin(), xml_string.end());
	xmlCopy.push_back('\0');

	rapidxml::xml_document<> doc;
	doc.parse<0>(&xmlCopy[0]);

	std::string nodeName;
	rapidxml::xml_node<char> *node = doc.first_node()->first_node(
			"BlastOutput_iterations");

	// Parse the xml for the nodes of interests
	for (rapidxml::xml_node<> *child = node->first_node(); child;
			child = child->next_sibling()) {
		if (rs.readID == child->first_node("Iteration_query-def")->value()) {
			if (rs.IsUnmapped() || rs.HasSoftClips())
				ParseNode(child, rs, mateMode);
			else
				ParseNode(child, *(rs.mate), mateMode);

		}

	}
}

// compare the chromosomes between two read sequences for sorting
// this compares the chromosomes of the blast alignment
bool CompareChr(ReadSequence read1, ReadSequence read2) {
	if (read1.blastAlignments.empty()) {
		return false;
	} else if (read2.blastAlignments.empty()) {
		return true;
	}
	std::string chromosome1 = read1.blastAlignments[0].chromosome.substr(3,
			std::strlen(read1.blastAlignments[0].chromosome.c_str()) - 2);
	std::string chromosome2 = read2.blastAlignments[0].chromosome.substr(3,
			std::strlen(read1.blastAlignments[0].chromosome.c_str()) - 2);

	if (atoi(chromosome1.c_str()) == atoi(chromosome2.c_str()))
		return read1.breakPoint.leftcoordinate > read2.breakPoint.leftcoordinate;
	return atoi(chromosome1.c_str()) > atoi(chromosome2.c_str());
}

void GetAdditionalTargetRegions(std::ifstream &target_region_file,
		std::vector<std::string> &target_regions) {
	std::string line;
	while (std::getline(target_region_file, line)) {
		target_regions.push_back(line);
	}
}

std::string get_selfpath() {
	char buff[FILENAME_MAX];
	ssize_t len = readlink("/proc/self/exe", buff, sizeof(buff) - 1);
	if (len != -1) {
		buff[len] = '\0';
		return std::string(buff, strlen(buff));
	} else {
		return std::string("\0");
	}
}

std::string _MergeSoftClipRegion(std::vector<std::string> softclips,
		std::vector<std::string> softclip_data, std::string _tid,
		std::string direction) {
	std::string softclip_region = softclip_data[0] + ":" + _tid + "\t"
			+ softclip_data[1] + "\t" + softclip_data[2] + "\t";
	std::string softclip_reads;

	std::vector<std::string>::iterator finalIter;
	for (std::vector<std::string>::iterator scIter = softclips.begin();
			scIter != softclips.end(); scIter++) {
		finalIter = softclips.end();
		--finalIter;

		if (scIter != finalIter) {
			softclip_region.append((*scIter).c_str());
			softclip_region.append("|" + direction + ";");
		} else {
			softclip_region.append((*scIter).c_str());
			softclip_region.append("|" + direction + "\n");
		}
	}

	return softclip_region;
}

std::string CheckFilePath(boost::program_options::variable_value val,
		FS::path cwd, std::string param, bool blast_database)
{
	std::string path;
	FS::path input;
	path = "";

	try
	{
		if(blast_database)
		{
			input = FS::path(GetParam<std::string>(val));
			path = std::string(input.c_str(), strlen(input.c_str()));
			return path;
		}
		else
		{
			input = FS::path(GetParam<std::string>(val));
		}
	}
	catch (boost::bad_any_cast &e)
	{
		return path;
	}

	if (input.is_relative()) {
		FS::path full_path(cwd / input);


		if (!FS::exists(full_path)) {
			printf("File %s not found...\n", full_path.c_str());
			return path;
		}
		else if (!FS::is_directory(full_path))
		{
			path = std::string(full_path.c_str(), strlen(full_path.c_str()));
		}
	}
	else
	{
		path = std::string(input.c_str(), strlen(input.c_str()));

			if (FS::exists(path.c_str()) && !FS::is_directory(path.c_str()))

			{
				return path;
			}
			else
			{
				path = "";
				return path;
			}
	}

	return path;
}

void DisplayHelp()
{
	printf("\nSWTools help\n");
	printf("usage:\t./SWTools -i file.bam -r reference.fa -b regions.bed \n");
	printf("\t\t-l /full/path/to/blastn --blast_reference path/to/blast/reference \n\n");
	printf("\t-i  input\t\tPath to BAM file\n");
	printf("\t-r  reference\t\tPath to Alignment Reference (.fasta)\n");
	printf("\t-b  bed\t\t\tPath to the Bed file to restrict search regions\n");
	printf("\t-l  blast\t\tPath to blasn executable to run\n");
	printf("\t--blast_reference\tPath to the reference to use with blast\n");
	printf("\t-t  threads\t\t[1] - Number of threads to use\n");
	printf("\t--min_softclips\t\t[2] - Minimum number of softclips required to count region\n");
	printf("\t--min_unmatched_bases\t[25] - Number of unmatched bases required to rerun BLASTN\n");
	printf("\t--min_softclip_size\t[5] - Minimum length of softclip\n");
	printf("\t--min_mapping_quality\t[20] - Minimum mapping quality to consider read\n");
	printf("\t--min_deviation\t\t[2.5] - Minimum standard deviations from mean to consider discordant\n");
	printf("\n");
}

void _AddProgramParameters(boost::program_options::options_description &desc)
{
	desc.add_options()("help,h", "Print help message")
			("input,i", boost::program_options::value<std::string>(), "Input BAM file")
			("bed,b", boost::program_options::value<std::string>(),"Bed region file")
			("reference,r",	boost::program_options::value<std::string>(),"Reference fasta file")
			("blast,l", boost::program_options::value<std::string>(),"Path to blast program")
			("blast_reference", boost::program_options::value<std::string>(), "Path to blast reference")
			("threads,t", boost::program_options::value<int unsigned>(), "[1] - Number of threads to use")
			("min_softclips", boost::program_options::value<int unsigned>(), "[2] - Minimum number of softclips to consider region")
			("min_unmatched_bases", boost::program_options::value<int unsigned>(), "[25] - Minimum number of bases to run BLASTN on")
			("min_softclip_size", boost::program_options::value<int unsigned>(), "[5] - Minimum lenght of softclip to keep")
			("min_mapping_quality", boost::program_options::value<int unsigned>(), "[20] - Minimum mapping quality to consider read")
			("min_deviation", boost::program_options::value<float>(), "[2.5] - Minimum deviation from the mean to consider discordant");

}

void _ParseProgramParameters(boost::program_options::variables_map vm)
{
	UNMATCHEDBASES = (vm["min_unmatched_bases"].empty()) ? UNMATCHEDBASES : GetParam<int unsigned>(vm["min_unmatched_bases"]);
	MIN_SOFTCLIPS = (vm["min_softclips"].empty()) ? MIN_SOFTCLIPS : GetParam<int unsigned>(vm["min_softclips"]);
	MIN_SOFTCLIP_SIZE = (vm["min_softclip_size"].empty()) ? MIN_SOFTCLIP_SIZE : GetParam<int unsigned>(vm["min_softclip_size"]);
	MIN_MAPPING_QUALITY = (vm["min_mapping_quality"].empty()) ? MIN_MAPPING_QUALITY : GetParam<int unsigned>(vm["min_mapping_quality"]);
	THREADS = (vm["threads"].empty()) ? THREADS : GetParam<int unsigned>(vm["threads"]);
	MIN_STANDARD_DEVIATION = (vm["min_deviation"].empty()) ? MIN_STANDARD_DEVIATION : GetParam<float>(vm["min_deviation"]);
}

int SizeOfParameterList(char** argv)
{
	int counter = 0;
	while(*(++argv))
	{
		counter++;
	}
	return counter;
}


void MergeSoftClips(std::string output_directory) {

	std::string cmd;

	// split the direction.bed file by delimiter and check validate the number of fields
	cmd =
			"sortBed -i " + output_directory
					+ "/hardsearch.direction.bed | mergeBed -nms -i stdin | egrep \";|,\" | awk '{OFS=\"\t\"}(NF==4)'";

	FILE *ls = popen(cmd.c_str(), "r");
	int character;
	character = fgetc(ls);

	std::string line;
	std::vector<std::string> softclip_data;
	std::vector<std::string> softclip_in;
	std::ofstream ofs_softclip_fixed;
	ofs_softclip_fixed.open(
			(output_directory + "/hardsearch.fixed.bed").c_str());

	// Get each line from the stream containing the mergeBed result of the
	// soft clipped reads and output them to a file containing the coordinates
	// and clipped sequence
	while (character != EOF) {
		if ((char) character == '\n') {
			softclip_data = SplitByDelimiter((char *) line.c_str(), "\t");
			softclip_in = SplitByDelimiter((char *) softclip_data[3].c_str(),
					";");
			std::vector<std::string> rightClips;
			std::vector<std::string> leftClips;
			std::string _tid;
			for (std::vector<std::string>::iterator scDataIter =
					softclip_in.begin(); scDataIter != softclip_in.end();
					scDataIter++) {

				std::vector<std::string> softclip_info = SplitByDelimiter(
						(char *) (*scDataIter).c_str(), "|");
				_tid = softclip_info[2];

				if (softclip_info.size() < 3)
					continue;
				std::string softclip_direction = softclip_info[1];
				if (softclip_direction == "Right") {
					rightClips.push_back((*scDataIter));
				} else {
					leftClips.push_back((*scDataIter));
				}
			}

			std::vector<std::string>::iterator finalIter;
			if (rightClips.size() != 0) {

				// construct string in the following format
				// CHR	POS_START	POS_END	SOFTCLIP_SEQUENCES|SOFTCLIP_DIRECTION
				std::string softclipOut = _MergeSoftClipRegion(rightClips,
						softclip_data, _tid, "Right");
				ofs_softclip_fixed.write(softclipOut.c_str(),
						strlen(softclipOut.c_str()));

			}

			if (leftClips.size() != 0) {
				std::string softclipOut = _MergeSoftClipRegion(leftClips,
						softclip_data, _tid, "Left");
				ofs_softclip_fixed.write(softclipOut.c_str(),
						strlen(softclipOut.c_str()));

			}

			line = "";
			character = fgetc(ls);
			continue;
		}
		line += character;
		character = fgetc(ls);
	}
	pclose(ls);
	ofs_softclip_fixed.flush();
	ofs_softclip_fixed.close();
}

/* Given the path to the BAM file, the output directory, and the target regions
 * identified by the soft clip filtering this function will output a .sam file containing
 * all of the unmapped read in the target regions
 */
bool _RetrieveUnmappedMateReads(std::string path_to_sam,
		std::string output_directory, std::vector<std::string> target_regions,
		int unsigned mapping_quality, SamParser *sp) {

	for (std::vector<std::string>::iterator target_iter =
			target_regions.begin(); target_iter != target_regions.end();
			target_iter++) {

		//printf("%s\n", target_iter->c_str());

		std::vector<std::string> target_data = SplitByDelimiter(
				(char *) target_iter->c_str(), "\t");

		int tid = atoi(target_data[4].c_str());
		int beg = atoi(target_data[1].c_str());
		int end = atoi(target_data[2].c_str());
		std::string direction = target_data[3];

		sp->RetrieveReadsForRegion(tid, beg, end, direction);

	}

	return true;
}
