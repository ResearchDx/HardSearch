import bs4
#from enum import Enum

class ReadSequence:
    def __init__(self, read):
        pass


class Read:

    INVERTED = 0x0001
    ALIGNED_LEFT = 0x0002
    ALIGNED_RIGHT = 0x0004
    GAPPED_ALIGNMENT = 0x0008
    BLAST_PAIRED = 0x0010
    READ_REVERSE_STRAND = 0x0010
    
    
    def __init__(self):
        self.read = None
        self.left_coordinate = None
        self.right_coordinate = None
        self.left_chromosome = None
        self.right_chromosome = None
        self.left_sequence = ""
        self.right_sequence = ""
        self.cigar_list = []
        self.discordant = False
        self.is_supporting_read = False

    def SetRead(self, read):
        self.read = read

        if self._HasDeletions():
            self._ParseForDeletionCoordinates()
        # else if discordant
        # ignore for now...
        elif self._HasBlastAlignments():
            self._ParseForBreakPointCoordinate()
            self._SetBreakPointSequences()
        #parse read for breakpoint coordinate

    def _HasDeletions(self):
        if self.read.deletion_region != None:
            if self.read.deletion_region.string != None:
                return True
        return False

    def _HasBlastAlignments(self):
        if self.read.alignments != None:
            return True
        return False

    def _ParseForDeletionCoordinates(self):
        # what if theres a deletion on each read
        # then we will have a problem

        left, right = [x for x in self.read.deletion_region.string.split("-")]
        self.left_chromosome = left.split(":")[0]
        self.left_coordinate = int(left.split(":")[1])

        self.right_chromosome = self.left_chromosome
        self.right_coordinate = int(right)





    def _ParseForBreakPointCoordinate(self):
        chromosome_1, coordinate_1 = self.read.start_coordinate_1.string.split(":")
        a_start_1, a_end_1 = [int(x) for x in self.read.match_coordinates_1.string.split("-")]
        m_start_1, m_end_1 = [int(x) for x in coordinate_1.split("-")]

        coordinate_1 = None
        coordinate_2 = None

        if(int(self.read.flag_1.string) & self.ALIGNED_LEFT):
            coordinate_1 = m_end_1
        else:
            coordinate_1 = m_start_1

        if(int(self.read.flag_1.string) & self.BLAST_PAIRED):
            chromosome_2, coordinate_2 = self.read.start_coordinate_2.string.split(":")
            a_start_2, a_end_2 = [int(x) for x in self.read.match_coordinates_2.string.split("-")]
            m_start_2, m_end_2 = [int(x) for x in coordinate_2.split("-")]

            if(int(self.read.flag_2.string) & self.ALIGNED_LEFT):
                coordinate_2 = m_end_2
            else:
                coordinate_2 = m_start_2

        if(coordinate_2 == None):
            self.left_coordinate = coordinate_1

            self.left_chromosome = chromosome_1
        else:
            self.left_coordinate = coordinate_1
            self.left_chromosome = chromosome_1
            self.right_coordinate = coordinate_2
            self.right_chromosome = chromosome_2

    def _SetBreakPointSequences(self):
        a_start_1, a_end_1 = [int(x) for x in self.read.match_coordinates_1.string.split("-")]

        a_start_1 -= 1

        read_seq = self.read.findAll("read_sequence")[1].string
        
        if(int(self.read.flag_1.string) & self.ALIGNED_LEFT):
            self.left_sequence = read_seq[a_start_1:a_end_1]
        else:
            self.left_sequence = read_seq[a_start_1:a_end_1]
            
        if(int(self.read.flag_1.string) & self.BLAST_PAIRED):
            a_start_2, a_end_2 = [int(x) for x in self.read.match_coordinates_2.string.split("-")]
            a_start_2 -= 1

            if(int(self.read.flag_1.string) & self.ALIGNED_LEFT):
                remainder = read_seq[a_end_1:]
            else:
                remainder = read_seq[:a_start_1]
            self.right_sequence = remainder[a_start_2:a_end_2]

        # set both a right and left sequence
        else:
            if(int(self.read.flag_1.string) & self.ALIGNED_LEFT):
                self.right_sequence = read_seq[a_end_1:]
            else:
                self.right_sequence = read_seq[:a_start_1]