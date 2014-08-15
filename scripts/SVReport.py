from __future__ import print_function
import bs4
import random
import Read
import collections
import Break
import LoadBed
import argparse
import os
import sys

def read_xml_file(xml_fp):
    raw_xml = ""
    
    try:
        f = open(xml_fp, 'r')
        line = f.readline()
        while(line):
            raw_xml += line
            line = f.readline()
    except IOError:
        raw_xml = ""
        return raw_xml
    else:
        f.close()
    return raw_xml

if __name__== "__main__":
    parser = argparse.ArgumentParser(description="Run SV Detection Script")
    parser.add_argument('-i', required=True, help="Reads in XML format")
    parser.add_argument('-b', required=True, help="Bed containing target region")
    parser.add_argument('--softclip_reads', type=int, help="Supporting softclip reads needed to investigate event")
    parser.add_argument('-c', type=int, help="Total read count required to report event")
    args = parser.parse_args()

    xml_file_path = args.i
    bed_file_path = args.b
    t_count = 20 if (args.c is None) else args.c
    sc_reads = 2 if (args.softclip_reads is None) else args.softclip_reads

    if not os.path.exists(xml_file_path):
        xml_file_path = os.path.join(os.getcwd(), xml_file_path)
    if not os.path.exists(bed_file_path):
        bed_file_path = os.path.join(os.getcwd(), bed_file_path)

    raw_xml = read_xml_file(xml_file_path)
    if (raw_xml == ""):
        print("XML file empty")
        sys.exit()
   

    bed_dict = LoadBed.LoadBedFile(bed_file_path)

    xml = bs4.BeautifulSoup(raw_xml, "html.parser")

    discordant_reads = xml.findAll("discordant_read_pair")
    softclip_reads = xml.findAll("softclip_read_pair")
    deletion_reads = xml.findAll("deletion_read_pair")

    read_list = []

    for x in softclip_reads:
        a = Read.Read()
        a.SetRead(x)
        read_list.append(a)

    for x in deletion_reads:
        a = Read.Read()
        a.SetRead(x)
        read_list.append(a)

    breakpoint_list= []
    breakpoint_dict = collections.defaultdict(int)

    read_list.sort(key=lambda x:x.left_coordinate)

    for x in read_list:
        entry = x.read.reference_name.string + ":" + str(x.left_coordinate)
        breakpoint_dict[entry] += 1

    break_list = []
    single_list = []

    # Determine if the reads support a breakpoint
    # Add a new breakpoint to the list if it does not support one
    for x in read_list:
        if(not x.left_coordinate is None and not x.right_coordinate is None):
            b = Break.Break()
            b.AddRead(x)
            match = False
            for y in break_list:
                if(y == b):
                    y.AddSupportingRead(x)
                    x.is_supporting_read = True
                    match = True
                    break
            if not match:
                if bed_dict != None:
                    if LoadBed.IsIntersectingBed((b.left_chromosome.lower(), b.left_coordinate, b.left_coordinate + 1), bed_dict) or LoadBed.IsIntersectingBed((b.right_chromosome.lower(), b.right_coordinate, b.right_coordinate + 1), bed_dict):
                        break_list.append(b)
                else:
                    break_list.append(b)

        if(x.left_coordinate is None or x.right_coordinate is None):
            single_list.append(x)


    # sort single list

    single_list.sort(key = lambda x:x.left_chromosome)
    chr_dict = collections.defaultdict(int)
    for x in single_list:
        chr_dict[x.left_chromosome] += 1

    # for each breakpoint in the list
    # check to see if any softclip / unmapped reads with
    # a single Blast mapping are related by levenshtein distance
    for x in break_list:
        x._FindConsensus()
        x._GetComparisonSequences(single_list, chr_dict)


    # can only parse for discordant when all the breakpoints have been found
    # unless we group discordant reads up as well...

    left_chromosome = None
    left_coordinate = None
    right_chromosome = None
    right_coordinate = None

    deviation = int(xml.discordant_reads['mean']) + (int(xml.discordant_reads['dev']) * int(xml.discordant_reads['num_dev']))
    for y in discordant_reads:
        if len(y.findAll("read")) != 2:
            continue
        if int(y.insert_size.string) < 200:
            continue
        left_chromosome = y.findAll("reference_name")[0].string
        left_coordinate = int(y.findAll("start_coordinate")[0].string)

        right_chromosome = y.findAll("reference_name")[1].string
        right_coordinate = int(y.findAll("start_coordinate")[1].string)

        last_high_l = None
        last_high_r = None
        difference_l = None
        difference_r = None
        high_index = None
        difference = None
        cur_index = 0
        diff = None
        lh = None
        for x in break_list:
            if left_chromosome == x.left_chromosome and right_chromosome == x.right_chromosome:

                difference_l = abs(left_coordinate - x.left_coordinate)
                difference_r = abs(right_coordinate - x.right_coordinate)
                if difference_l > deviation or difference_r > deviation:
                    cur_index += 1
                    continue

                diff = difference_r + difference_l
                if last_high_l == None:
                    last_high_l = difference_l
                    last_high_r = difference_r
                    lh = diff
                    high_index = cur_index

                elif diff < lh:
                    last_high_l = difference_l
                    last_high_r = difference_r
                    lh = diff
                    high_index = cur_index

            if left_chromosome == x.right_chromosome and right_chromosome == x.left_chromosome:
                difference_l = abs(left_coordinate - x.right_coordinate)
                difference_r = abs(right_coordinate - x.left_coordinate)
                if difference_l > deviation or difference_r > deviation:
                    cur_index += 1  
                    continue


                diff = difference_r + difference_l
                if last_high_l == None:
                    last_high_l = difference_l
                    last_high_r = difference_r
                    lh = diff
                    high_index = cur_index
                elif diff < lh:
                    last_high_l = difference_l
                    last_high_r = difference_r
                    lh = diff
                    high_index = cur_index
            cur_index += 1

        if high_index != None:
            a = Read.Read()
            a.SetRead(y)
            a.discordant = True
            #        print("left: {} right: {} index {}".format(str(last_high_l), str(last_high_r), str(high_index)))
            break_list[high_index].AddSupportingRead(a)
            a.is_supporting_read = True


    index = 0

    print("Event ID", "Reference Name", "Coordinate", "Reference Name", "Coordinate", "Event Type", "Reads",sep="\t")

    for x in break_list:
        # narrows down the range of reads to 1
        # majority of false positives will not have high read count...
        if(len(x.reads) > t_count):
            event_reads = ""
            greatest = 0
            greatest_event = ""
            x._DetermineEventFromReads()
            for key in x.events.keys():
                event_reads += "{}:{};".format(key[:3], str(x.events[key]))
                if key == "discordant":
                    continue
                if x.events[key] > greatest:
                    greatest = x.events[key]
                    greatest_event = key
            event_reads += "{}:{}".format("total", str(len(x.reads)))
            print(str(index), x.left_chromosome, str(x.left_coordinate),
                x.right_chromosome, str(x.right_coordinate), greatest_event, event_reads, sep="\t")
        index += 1
