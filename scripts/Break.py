import collections
import Read
import jellyfish
import bs4

class Base:
    def __init__(self):
        self._A = 0
        self._T = 0
        self._G = 0
        self._C = 0
        self._N = 0
        self._INDEL = 0
        self.count = 0
        self.coordinate = None
        self.chromosome = None

    def __iadd__(self, other):
        if other.upper() == 'A':
            self._A += 1
        elif other.upper() == 'T':
            self._T += 1
        elif other.upper() =='G':
            self._G += 1
        elif other.upper() == 'C':
            self._C += 1
        elif other.upper() == 'N':
            self._N += 1

        self.count += 1
        return self

    def SetCoordinates(self, reference_name, coord):
        if self.chromosome == None:
            self.chromosome = reference_name
        if self.coordinate == None:
            self.coordinate = coord

        else:
            if (self.chromosome != reference_name):
                raise Exception("ReferenceNameError", "new value {} is not equivalent to current value {}".format(reference_name, self.chromosome))
            if (self.coordinate != coord):
                raise Exception("CoordinateError", "new value {:d} is not equivalent to current value {:d}".format(coord, self.coordinate))



    def GetConsensusBase(self):
        ''' @ abstract  Get the base with the highest count at a particular position
            @ return    String containing the base with the greatest count at that position
        '''
        bases = []
        bases.append(['A',self._A])
        bases.append(['T',self._T])
        bases.append(['G',self._G])
        bases.append(['C',self._C])
        bases.append(['N',self._N])

        bases.sort(key=lambda x:x[1])
        #if bases[-1][1] == 0:
        #    return '-'
        if bases[-1][1] >= 20:
            return bases[-1][0]
        else:
            return '-'

class ConsensusSequence:
    def __init__(self):
        self.left_chromosome = ""
        self.left_start = -1
        self.left_end = -1
        self.left_sequence = ""
        self.left_orientation = ""


        self.right_chromosome = None
        self.right_start = None
        self.right_end = None
        self.right_sequence = None
        self.right_orientation = None


    def AddRead(self, read, orientation):
        chromosome, coordinates = read.split(":")
        start, end = [int(x) for x in coordinates.split("-")]

        if self.left_chromosome == "":
            self.left_chromosome = chromosome

        if self.left_start == -1:
            self.left_start = start

        if self.left_end == -1:
            self.left_end = end

        if self.left_orientation == "":
            self.left_orientation = orientation

        if orientation == "forward":
            if self.left_start > start:
                self.left_start = start
            if self.left_end < end:
                self.left_end = end
        elif orientation == "reverse":
            if self.left_start < start:
                self.left_start = start
            if self.left_end > end:
                self.left_end = end


class Break:
    def __init__(self):
        '''
            @param  _min_sep_distance     distance of the end of the read to the
                                          coordinates of the breakpoint to consider
        '''
        self.left_coordinate = None
        self.left_chromosome = None
        self.right_coordinate = None
        self.right_chromosome = None
        self.left_breakpoints = set()
        self.right_breakpoints = set()

        self.events = collections.defaultdict(int)


        self._min_sep_distance = 5
        self._read_length = 155
        self._read_dict = collections.defaultdict(int)
        self.reads = []
        self.count = 1

        # initialze class instance variables before calling the following functions
        # some functions depend on certain instance variables being set...
        self._InitConsensusSequenceContainers()

    def __eq__(self, other):
        if(self.left_chromosome == other.left_chromosome and
           self.right_chromosome == other.right_chromosome):
            if(abs(self.left_coordinate - other.left_coordinate) < self._min_sep_distance and
               abs(self.right_coordinate - other.right_coordinate) < self._min_sep_distance):
                return True
        if(self.left_chromosome == other.right_chromosome and
           self.right_chromosome == other.left_chromosome):
            if(abs(self.left_coordinate - other.right_coordinate) < self._min_sep_distance and
               abs(self.right_coordinate - other.left_coordinate) < self._min_sep_distance):
                return True
        return False

    def AddRead(self, read):
        self.left_coordinate = read.left_coordinate
        self.right_coordinate = read.right_coordinate
        self.left_chromosome = read.left_chromosome
        self.right_chromosome = read.right_chromosome
        self.AddSupportingRead(read)
#        self.reads.append(read)

    def AddSupportingRead(self,read):
        '''
            @abstract   Add a supporting read to the list of reads spanning this breakpoint equivalent
        '''
        if not read.is_supporting_read:
            if self._read_dict[read.read.read_id.string] == 0:
                self.reads.append(read)
                self.count += 1
                self._read_dict[read.read.read_id.string] += 1

    def _DetermineEventFromReads(self):
        for x in self.reads:
            #if x.discordant:
            #    continue
            if x._HasDeletions():
                self.events["deletion"] += 1
            elif x.right_chromosome != None and x.left_chromosome != x.right_chromosome:
                self.events["translocation"] += 1
            elif x.read.flag_1 != None:

                if(int(x.read.findAll('flag')[0].string) & Read.Read.READ_REVERSE_STRAND):
                    if(int(x.read.flag_1.string) & Read.Read.ALIGNED_LEFT):
                        if(x.read.orientation_1.string == 'forward'):
                            if(x.read.orientation_2 != None and x.read.orientation_2.string == 'reverse'):
                                self.events["inversion"] += 1
                        elif(x.read.orientation_1.string == 'reverse'):
                            self.events["inversion"] += 1

                    elif(int(x.read.flag_1.string) & Read.Read.ALIGNED_RIGHT):
                        if(x.read.orientation_1.string == 'reverse'):
                            if(x.read.orientation_2 != None and x.read.orientation_2.string == 'forward'):
                                self.events["inversion"] += 1
                        elif(x.read.orientation_2.string == 'reverse'):
                            self.events["inversion"] += 1

                else:
                    if(int(x.read.flag_1.string) & Read.Read.ALIGNED_LEFT):
                        if(x.read.orientation_1 == 'reverse'):
                            if(x.read.orientation_2 != None and x.read.orientation_2.string == 'forward'):
                                self.events["inversion"] += 1
                        elif(x.read.orientation_1.string == 'forward'):
                            self.events["inversion"] += 1
                    elif(int(x.read.flag_1.string) & Read.Read.ALIGNED_RIGHT):
                        if(x.read.orientation_1.string == 'forward'):
                            if(x.read.orientation_2 != None and x.read.orientation_2.string == 'reverse'):
                                self.events["inversion"] += 1
                        elif(x.read.orientation_2.string == 'forward'):
                            self.events["inversion"] += 1

            elif x.discordant:
                self.events["discordant"] += 1

    def _InitConsensusSequenceContainers(self):
        ''' @abstarct   set the list container for the breakpoint
        '''
        self.left_consensus_sequence = ConsensusSequence()
        self.left_sequence = []
        self.right_sequence = []
        read_length = self._read_length + (2 * self._min_sep_distance)
        for x in range( read_length):
            self.left_sequence.append(Base())
            self.right_sequence.append(Base())


    def _FindConsensus(self):
        '''
            @TODO   Refactoring required, currently the left and right sequences are parsed
                    in a non intuitive matter, and are also appended as such...
                    This is currently not readable...

            @discussion     this function iterates through the supporting reads in the read
                            list and then counts the nucletide base in the left and right
                            side of the realigned segment. All the counts are added to a list
                            that keeps track of the position and base count.
        '''
        for read in self.reads:
            offset = 0
            c_index = 0

#            l_start, l_end = [int(x) for x in read.match_coordinates_1.split("-")]

            if read._HasDeletions():
                continue

            if(read.left_chromosome == self.left_chromosome and
               abs(read.left_coordinate - self.left_coordinate) < self._min_sep_distance):
                offset = read.left_coordinate - self.left_coordinate
                self.left_breakpoints.add(read.left_coordinate)

                if(int(read.read.flag_1.string) & Read.Read.ALIGNED_RIGHT):
                    for i in range(self._min_sep_distance - offset, len(read.left_sequence) + self._min_sep_distance - offset):
                        char = read.left_sequence[c_index]
                        self.left_sequence[i] += char
                        self.left_sequence[i].SetCoordinates(read.left_chromosome, read.left_coordinate - c_index)
                        c_index += 1
                else:
                    for i in range(self._min_sep_distance - offset, len(read.left_sequence) + self._min_sep_distance - offset):
                        char = read.left_sequence[len(read.left_sequence) - 1 - c_index]
                        self.left_sequence[i] += char
                        self.left_sequence[i].SetCoordinates(read.left_chromosome, read.left_coordinate - c_index)
                        c_index += 1

            elif(read.left_chromosome == self.right_chromosome and
               abs(read.left_coordinate - self.right_coordinate) < self._min_sep_distance):
                offset = read.left_coordinate - self.right_coordinate
                self.right_breakpoints.add(read.left_coordinate)
                # TODO this should only work for 'left' reads that are ALIGNED_RIGHT...
                # fix that...
                if(int(read.read.flag_1.string) & Read.Read.ALIGNED_RIGHT):
                    for i in range(self._min_sep_distance - offset, len(read.left_sequence) + self._min_sep_distance - offset):
                        char = read.left_sequence[c_index]
                        self.right_sequence[i] += char
                        self.right_sequence[i].SetCoordinates(read.left_chromosome, read.left_coordinate - c_index)
                        c_index += 1
                else:
                    for i in range(self._min_sep_distance - offset, len(read.left_sequence) + self._min_sep_distance - offset):
                        char = read.left_sequence[len(read.left_sequence) - 1 - c_index]
                        self.right_sequence[i] += char
                        self.right_sequence[i].SetCoordinates(read.left_chromosome, read.left_coordinate - c_index)
                        c_index += 1

            c_index = 0
            offset = 0

            if(read.right_chromosome == self.right_chromosome and
                abs(read.right_coordinate - self.right_coordinate) < self._min_sep_distance):
                offset = read.right_coordinate - self.right_coordinate
                self.right_breakpoints.add(read.right_coordinate)

                if(int(read.read.flag_1.string) & Read.Read.ALIGNED_RIGHT):
                    for i in range(self._min_sep_distance - offset, len(read.right_sequence) + self._min_sep_distance - offset):
                        char = read.right_sequence[len(read.right_sequence) - 1 - c_index]
                        self.right_sequence[i] += char
                        self.right_sequence[i].SetCoordinates(read.right_chromosome, read.right_coordinate - c_index)
                        c_index += 1
                else:
                    for i in range(self._min_sep_distance - offset, len(read.right_sequence) + self._min_sep_distance - offset):
                        char = read.right_sequence[c_index]
                        self.right_sequence[i] += char
                        self.right_sequence[i].SetCoordinates(read.right_chromosome, read.right_coordinate - c_index)
                        c_index += 1

            elif(read.right_chromosome == self.left_chromosome and
                abs(read.right_coordinate - self.left_coordinate) < self._min_sep_distance):
                offset = read.right_coordinate - self.left_coordinate
                self.left_breakpoints.add(read.right_coordinate)

                if(int(read.read.flag_1.string) & Read.Read.ALIGNED_RIGHT):
                    for i in range(self._min_sep_distance - offset, len(read.right_sequence) + self._min_sep_distance - offset):
                        char = read.right_sequence[len(read.right_sequence) - 1 - c_index]
                        self.left_sequence[i] += char
                        self.left_sequence[i].SetCoordinates(read.right_chromosome, read.right_coordinate - c_index)
                        c_index += 1
                else:
                    for i in range(self._min_sep_distance - offset, len(read.right_sequence) + self._min_sep_distance - offset):
                        char = read.right_sequence[len(read.right_sequence) - 1 - c_index]
                        self.left_sequence[i] += char
                        self.left_sequence[i].SetCoordinates(read.right_chromosome, read.right_coordinate - c_index)
                        c_index += 1


            # this is always going to look at the left read but the left read will not always be the same
#            for i in range(self._min_sep_distance - offset, len(read.left_sequence) + self._min_sep_distance):
#                char = read.left_sequence[len(read.left_sequence) - 1 - c_index]
#                self.left_sequence[i] += char
#                c_index += 1

            self.left_consensus_sequence.AddRead(read.read.start_coordinate_1.string, read.read.orientation_1.string)





    def _FlipSequence(self, sequence):
        ''' TODO: only the left sequence is flipped / flipped right as well

            Flip the consensus sequence so that it is no longer in the opposite
            orientation
        '''
        flipped_sequence = ""

        index = 0
        for x in range(len(sequence)):
            flipped_sequence += sequence[-x - 1]

        return flipped_sequence


    def _GetComparisonSequences(self, reads, chr_dict, debug=False):
        '''
            @abstract           Adds supporting reads with levenshtein_distance > 0.95 
                                to breakpoint consensus
            @param  reads       List containing reads with single Blast mapping
            @param  chr_dict    Dictionary containing all the chromosomes with left_breakpoints
            @param  debug       Boolean used to print output

        '''
        found = False
        if chr_dict[self.left_chromosome] == 0 and chr_dict[self.right_chromosome] == 0:
            return
        for read in reads:
            sequence = ""
            start = 0
            end = 0
            left = False

            if read.left_chromosome != self.left_chromosome and read.left_chromosome != self.right_chromosome:
                if found:
                    break
                continue
            else:
                found = True
                if(not abs(read.left_coordinate - self.left_coordinate) < 5 or
                   not abs(read.left_coordinate - self.left_coordinate) < 5 ):
                    continue
#                TODO:  this works for the current use case (i.e when the breakpoint is set such that
#                       the smaller one is left 29448093 versus 42493955) Might not work if it is reversed...
#                       because the left sequence is always reversed currently... will there be an instance
#                       where this is not true...?
                if(abs(read.left_coordinate - self.left_coordinate) < 5):
                    left = True

                for y in range(len(read.left_sequence)):
                    sequence += self._GetBaseAtPosition(read.left_chromosome, read.left_coordinate - y)

                if(left):
                    sequence = self._FlipSequence(sequence)


                start, end = self._GetReferenceIndex(sequence)
                sequence = sequence if (sequence.strip("-")== "") else sequence.strip("-")

                sequence_r = ""

                for y in range(len(read.right_sequence)):
                    if left:
                        sequence_r += self._GetBaseAtPosition(self.right_chromosome, self.right_coordinate - y)
                    else:
                        sequence_r += self._GetBaseAtPosition(self.left_chromosome, self.left_coordinate - y)

                if (sequence_r + sequence).strip("-") == "":
                    continue

                if not left:
                    sequence_r = self._FlipSequence(sequence_r)
                    sequence = sequence_r + sequence
                    lev_dist = jellyfish.levenshtein_distance(sequence, read.right_sequence + read.left_sequence[start:end] )
                    if debug:
                        print("Comparing levenshtein_distance for:\n" + sequence + "\n" + read.right_sequence + read.left_sequence[start:end])

                else:
                    sequence = sequence + sequence_r
                    lev_dist = jellyfish.levenshtein_distance(sequence, read.left_sequence[start:end] + read.right_sequence)
                    if debug:
                        print("Comparing levenshtein_distance for:\n" + sequence + "\n" + read.left_sequence[start:end] + read.right_sequence)

                rel_lev = 1 - (lev_dist/max(len(sequence), len(read.left_sequence)))

                if rel_lev > 0.95:
                    if not read.is_supporting_read:
                        self.AddSupportingRead(read)
                        read.is_supporting_read = True


    def _GetReferenceIndex(self, sequence):
        start = 0
        end = len(sequence) - 1
        for x in sequence:
            if x == "-":
                start += 1
            else:
                break

        for x in range(len(sequence)):
            if sequence[len(sequence) - 1 - x] == "-":
                end -= 1
        end += start
        return (start, end)

    def _GetMajorityBaseStart(self, chromosome, coordinates , sequence):
        coordinate = None
        count = None
        result_chr = None

        for x in sequence:
            if x.chromosome == None:
                continue
            for y in coordinates:
                if chromosome == x.chromosome:
                    if y == x.coordinate:


                        if count == None:
                            count = x.count
                            coordinate = y
                        if x.count > count:
                            count = x.count
                            coordinate = y
                        break
        return coordinate



    def _GetBaseAtPosition(self, chromosome, coordinate):
        for x in self.left_sequence:
            if x.chromosome == None:
                continue
            if chromosome == x.chromosome:
                if coordinate == x.coordinate:
                    return x.GetConsensusBase()
            else:
                break

        for x in self.right_sequence:
            if x.chromosome == None:
                continue
            if chromosome == x.chromosome:
                if coordinate == x.coordinate:
                    return x.GetConsensusBase()
        return "-"



    def GetLeftSequence(self):
        result = ""
        for x in self.left_sequence:
            result += x.GetConsensusBase()
        return result

    def GetRightSequence(self):
        result = ""
        for x in self.right_sequence:
            result += x.GetConsensusBase()
        return result

    def PrintBreakStatistics(self):
        statistics = "There are {count} softclip reads spanning this breakpoint: ".format(count = self.count)
        statistics += "{chr_left}:{coord_left} - {chr_right}:{coord_right}".format( chr_left = self.left_chromosome,
                                                                                    coord_left = self.left_coordinate,
                                                                                    chr_right = self.right_chromosome,
                                                                                    coord_right = self.right_coordinate)