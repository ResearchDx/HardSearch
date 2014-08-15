import bs4
import Break
import Read

    raw_xml = """<softclip_read_pair>
<read>
<read_id>ALK_EML4:FUSION-4810</read_id>
<flag>99</flag>
<reference_name>chr2</reference_name>
<start_coordinate>29447680</start_coordinate>
<mapping_quality>29</mapping_quality>
<cigar>150M</cigar>
<insert_size>414</insert_size>
<read_sequence>TACCGGCAGATCCCTTTGCCTGCAGGGGCCTGGCCTGCGAGGGCTCTCAAGAGCCTTTCCCTCTGCCCTTTTCAAGCCTCTGCCCATCTGTCCTGGGCATGTCTCTGCCAGCAGTAAGAGCTGGTTGGGAACACACTGAGTTCTCTGTGA</read_sequence>
<quality></quality>
</read>
<read>
<read_id>ALK_EML4:FUSION-4810</read_id>
<flag>147</flag>
<reference_name>chr2</reference_name>
<start_coordinate>29448016</start_coordinate>
<mapping_quality>29</mapping_quality>
<cigar>78M72S</cigar>
<insert_size>414</insert_size>
<read_sequence>TAGGGGTGGGGGCGAGCTTTCACCATCGTGATGGACACTGAAGGAGCTCCCCACCCCCTGATCAGCCAGGAGGATGCAGTTACTATATTTAGAGTAGAAAATCACTCAGAGACATTGTCATTCAGCCTAGGTCATAAGGCTATTATGCAT</read_sequence>
<quality></quality>
<alignments>
<start_coordinate_1>chr2:29448016-29448090</start_coordinate_1>
<orientation_1>forward</orientation_1>
<match_coordinates_1>1-75</match_coordinates_1>
<bitscore_1>145</bitscore_1>
<length_1>75</length_1>
<flag_1>26</flag_1>
<start_coordinate_2>chr2:42493955-42493883</start_coordinate_2>
<orientation_2>reverse</orientation_2>
<match_coordinates_2>3-75</match_coordinates_2>
<bitscore_2>135</bitscore_2>
<length_2>73</length_2>
<flag_2>5</flag_2>
</alignments>
</read>
</softclip_read_pair>"""

    bp = Read.Read()
    xml = bs4.BeautifulSoup(raw_xml, "html.parser")
    bp.SetRead(xml)

    b = Break()
    b.AddRead(bp)

    b._FindConsensus()
    assert(b.GetLeftSequence() == "-----TAGGAGGACCGACTAGTCCCCCACCCCTCGAGGAAGTCACAGGTAGTGCTACCACTTTCGAGCGGGGGTGGGGAT-------------------------------------------------------------------------------------")
    assert(b.GetRightSequence() == "-----AGTTACTATATTTAGAGTAGAAAATCACTCAGAGACATTGTCATTCAGCCTAGGTCATAAGGCTATTATGCAT---------------------------------------------------------------------------------------")
