import Read
import bs4


if __name__ == "__main__":

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

    raw_xml2 = """<softclip_read_pair>
<read>
<read_id>ALK_EML4:FUSION-3906</read_id>
<flag>163</flag>
<reference_name>chr2</reference_name>
<start_coordinate>29447631</start_coordinate>
<mapping_quality>29</mapping_quality>
<cigar>150M</cigar>
<insert_size>461</insert_size>
<read_sequence>TGAGGGATGGCACCATATGGGGACACAGTGTGTGCTGCCATCTTCCTTCTCCCGGCAGATCCCTTTGCCTGCAGGGGCCTGGCCTGCGAGGGCTCTCAAGAGCCTTTCCCTCTGCCCTTTTCAAGCCTCTGCCCATCTGTCCTGGGCATG</read_sequence>
<quality></quality>
</read>
<read>
<read_id>ALK_EML4:FUSION-3906</read_id>
<flag>83</flag>
<reference_name>chr2</reference_name>
<start_coordinate>29448019</start_coordinate>
<mapping_quality>29</mapping_quality>
<cigar>73M77S</cigar>
<insert_size>461</insert_size>
<read_sequence>GGGTGGGGGCGAGCTTTCACCATCGTGATGGACACTGAAGGAGCTCCCCACCCCCTGATCAGCCAGGAGGATATAGTTACTATATTTAGAGTAGAAAACCACTCAGAGACATTGTCATTCAGCCTAGGTCATAAGGCTATAATGCATAAG</read_sequence>
<quality></quality>
<alignments>
<start_coordinate_1>chr2:42493958-42493880</start_coordinate_1>
<orientation_1>reverse</orientation_1>
<match_coordinates_1>72-150</match_coordinates_1>
<bitscore_1>147</bitscore_1>
<length_1>79</length_1>
<flag_1>21</flag_1>
<start_coordinate_2>chr2:29448019-29448089</start_coordinate_2>
<orientation_2>forward</orientation_2>
<match_coordinates_2>1-71</match_coordinates_2>
<bitscore_2>137</bitscore_2>
<length_2>71</length_2>
<flag_2>2</flag_2>
</alignments>
</read>
</softclip_read_pair>"""

    bp = Read.Read()
    xml = bs4.BeautifulSoup(raw_xml, "html.parser")
    bp.SetRead(xml)

    bp._ParseForBreakPointCoordinate()
    bp._SetBreakPointSequences()

    assert(bp.left_coordinate == 29448090)
    assert(bp.right_coordinate == 42493955)
    assert(bp.left_sequence == "TAGGGGTGGGGGCGAGCTTTCACCATCGTGATGGACACTGAAGGAGCTCCCCACCCCCTGATCAGCCAGGAGGAT")
    assert(bp.right_sequence == "AGTTACTATATTTAGAGTAGAAAATCACTCAGAGACATTGTCATTCAGCCTAGGTCATAAGGCTATTATGCAT")

    bp = Read.Read()
    xml = bs4.BeautifulSoup(raw_xml2, "html.parser")
    bp.SetRead(xml)

    bp._ParseForBreakPointCoordinate()
    bp._SetBreakPointSequences()

    assert(bp.left_coordinate == 42493958)
    assert(bp.right_coordinate == 29448089)
    assert(bp.left_sequence == "TATAGTTACTATATTTAGAGTAGAAAACCACTCAGAGACATTGTCATTCAGCCTAGGTCATAAGGCTATAATGCATAAG")
    assert(bp.right_sequence == "GGGTGGGGGCGAGCTTTCACCATCGTGATGGACACTGAAGGAGCTCCCCACCCCCTGATCAGCCAGGAGGA")



    a = '''     <Softclip_Read_Pair>
            <Read>
                <Read_ID>M01801:26:000000000-A43C6:1:2110:8249:22370</Read_ID>
                <Flag>137</Flag>
                <Reference_Name>chr2</Reference_Name>
                <Start_Coordinate>29447927</Start_Coordinate>
                <Mapping_Quality>37</Mapping_Quality>
                <Cigar>151M</Cigar>
                <Insert_Size>0</Insert_Size>
                <Read_Sequence>TTGGGGAAGAGTGGGCTAGTGCATTACATAGGGTGGGAGCCAAACAGGAGCTGCGCCGGTGGAAGCATGTGGGAGCTAGAAGTGACGTCTAGGGGTGGGGGCGAGCTTTCACCATCGTGATGGACACTGAAGGAGCTCCCCACCCCCTGAT</Read_Sequence>
                <Quality/>
            </Read>
            <Read>
                <Read_ID>M01801:26:000000000-A43C6:1:2110:8249:22370</Read_ID>
                <Flag>69</Flag>
                <Reference_Name>chr2</Reference_Name>
                <Start_Coordinate>29447927</Start_Coordinate>
                <Mapping_Quality>0</Mapping_Quality>
                <Cigar>*</Cigar>
                <Insert_Size>0</Insert_Size>
                <Read_Sequence>ATGGGGAATAACTGTGTTTAGTGGTTGGTCCTGGCTCTTATGCATAATAGCCTTATGACCTAAGCTGAATGACAATGTCTCTGAGTGGTTTTCTACTCTAAATATAGTAACTAGTATCCTCCTGGCTGATCAGGGGGTGGGGAGCTCCTTC</Read_Sequence>
                <Quality/>
                <Alignments>
                    <Start_Coordinate_1>chr2:42493844-42493956</Start_Coordinate_1>
                    <Orientation_1>forward</Orientation_1>
                    <Match_Coordinates_1>1-113</Match_Coordinates_1>
                    <Bitscore_1>212</Bitscore_1>
                    <Length_1>113</Length_1>
                    <Flag_1>18</Flag_1>
                    <Start_Coordinate_2>chr2:29448092-29448055</Start_Coordinate_2>
                    <Orientation_2>reverse</Orientation_2>
                    <Match_Coordinates_2>1-38</Match_Coordinates_2>
                    <Bitscore_2>74</Bitscore_2>
                    <Length_2>38</Length_2>
                    <Flag_2>5</Flag_2>
                </Alignments>
            </Read>
        </Softclip_Read_Pair>'''
    b = Read.Read()
    xml_b = bs4.BeautifulSoup(a, 'html.parser')
    b.SetRead(xml_b)

    assert(b.left_coordinate == 42493956)
    assert(b.right_coordinate == 29448092)


    c = '''     <Softclip_Read_Pair>
            <Read>
                <Read_ID>KMT2A_ABL1_TRANS-8832</Read_ID>
                <Flag>147</Flag>
                <Reference_Name>chr9</Reference_Name>
                <Start_Coordinate>133730414</Start_Coordinate>
                <Mapping_Quality>29</Mapping_Quality>
                <Cigar>150M</Cigar>
                <Insert_Size>326</Insert_Size>
                <Read_Sequence>GAGGTCCATCTCGCTGAGATACGAAGGGAGGGTGTACCATTACAGGATCAACACTGCTTCTGATGGCAAGGTAGGGGACCCTTGGCAGGGGGCGCTGATGGGCCCAGGGCAGGGGAACCAGAGGTCCTGCTGTCGGATTGATAAATTATT</Read_Sequence>
                <Quality/>
            </Read>
            <Read>
                <Read_ID>KMT2A_ABL1_TRANS-8832</Read_ID>
                <Flag>99</Flag>
                <Reference_Name>chr9</Reference_Name>
                <Start_Coordinate>133730238</Start_Coordinate>
                <Mapping_Quality>29</Mapping_Quality>
                <Cigar>108S42M</Cigar>
                <Insert_Size>326</Insert_Size>
                <Read_Sequence>TTCCCCTTTGCAAATAGAGTCAACATCTCCCACAGAACCAATTTCAGCCTCTGAAAATCCAGGAGATGGTCCAGTGGCCCAACCAAGCCCCAATAATACCTCATGCCTGCCCAAACCAAAAATGGCCAAGGCTGGGTCCCAAGCAACTAC</Read_Sequence>
                <Quality/>
                <Alignments>
                    <Start_Coordinate_1>chr11:118374251-118374357</Start_Coordinate_1>
                    <Orientation_1>forward</Orientation_1>
                    <Match_Coordinates_1>1-107</Match_Coordinates_1>
                    <Bitscore_1>206</Bitscore_1>
                    <Length_1>107</Length_1>
                    <Flag_1>26</Flag_1>
                    <Start_Coordinate_2>chr9:133730238-133730279</Start_Coordinate_2>
                    <Orientation_2>forward</Orientation_2>
                    <Match_Coordinates_2>2-43</Match_Coordinates_2>
                    <Bitscore_2>81</Bitscore_2>
                    <Length_2>42</Length_2>
                    <Flag_2>4</Flag_2>
                </Alignments>
            </Read>
        </Softclip_Read_Pair>
        '''
    xml_c = bs4.BeautifulSoup(c, 'html.parser')
    d = Read.Read()
    d.SetRead(xml_c)

    assert(d.left_coordinate == 118374357)
    assert(d.right_coordinate == 133730238)
    print("assertions passed")