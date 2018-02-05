#!/usr/bin/env python

"""
This script makes microscopy probes, and the primers to amplify those probes.
Basically, it takes 1kb on a 500bp sliding window, and blasts it against your
assembly. If there are no other blast hits, it then takes that region +- 1kb
and designs primers to amplify it.
"""

import sys
import subprocess
import re
from Bio.Seq import Seq

if len(sys.argv) != 7:
    print "Usage: create.microscopy.probes.py <fasta file> <blast db> \
           <lowest Tm> <highest Tm> <Mininum number of mismatches> \
           <Maximum occurrences of primer in assembly (default = 1)>\n"
    quit()

low_tm = int(sys.argv[3])
hi_tm = int(sys.argv[4])
min_mismatches = int(sys.argv[5])
max_occurrences = int(sys.argv[6])

# Takes the target sequences and reads them into a dictionary
seqs = open(sys.argv[1], "r")
seqs_dict = {}
for line in seqs:
    if line.startswith(">"):
        header = line.strip()
        seqs_dict[header] = ""
    else:
        seqs_dict[header] += line.strip()
seqs.close()

# Takes the sequences and chops them up into 1kb pieces
probes_dict = {}
for item in seqs_dict:
    for x in range(0, len(seqs_dict[item]), 500):
        probe = seqs_dict[item][x:x + 1000]
        header = "%s_%s" % (item, x + 1)
        probes_dict[header] = probe

temp_file = open("probes.temp.fasta", "w")
for item in probes_dict:
    temp_file.write("%s\n%s\n" % (item, probes_dict[item]))

# Blasts the putative probes against the assembly
blastn = "blastn -task blastn -query probes.temp.fasta -db %s -outfmt 5 \
          -evalue 1 -num_threads 10 -soft_masking false -qcov_hsp_perc 10 \
          -dust no -out ./probes.blastn.xml" % (sys.argv[2])
parse = "parse.xml.py ./probes.blastn.xml"
subprocess.call(blastn, stdin=None, stdout=None, stderr=None, shell = True)
subprocess.call(parse, stdin=None, stdout=None, stderr=None, shell = True)

# Determines how many BLAST hits each probe had, and writes those with only
# one hit (to itself) to a file
file = open("probes.blastn.xml.parsed", "r")
out = open("good.probes.fasta", "w")
targets = {}
for line in file:
    if line.startswith("Query_"):
        if len(targets) != 1:
            targets = {}
            regex = re.match("Query.*?: (.*?)\n", line)
            query = regex.group(1)
        else:
            out.write(">%s\n%s\n" % (query, probes_dict[">" + query]))
            regex = re.match("Query.*?: (.*?)\n", line)
            query = regex.group(1)
            targets = {}
    elif line.startswith("Target"):
        regex = re.match("(Target.*?):.*?\n", line)
        targets[regex.group(1)] = ""
else:
    if len(targets) == 1:
        out.write(">%s\n%s\n" % (query, probes_dict[">" + query]))

rm = "rm probes.temp.fasta"
subprocess.call(rm, stdin=None, stdout=None, stderr=None, shell = True)
out.close()

probes = open("good.probes.fasta", "r")
targets = {}
for line in probes:
    if line.startswith(">"):
        header = line.strip()
        targets[header] = ""
    else:
        targets[header] += line.strip()
probes.close()

# Takes the good probes and 1kb up/downstream and writes it to a file to check
# for primers
seqs = open(sys.argv[1], "r")
out = open("probes.for.primers", "w")
for line in seqs:
    if line.startswith(">"):
        pass
    else:
        for item in targets:
            probe = "(\w{1000})" + str(targets[item]) + "(\w{1000})"
            regex = re.search(probe, line)
            if regex:
                out.write("%s\n%sN%s\n" % (item, regex.group(1), regex.group(2)))
seqs.close()
out.close()

potential_primers = open("probes.for.primers", "r")
seqs_dict = {}
primers_dict = {}
for line in potential_primers:
    if line.startswith(">"):
        header = line.strip()
        seqs_dict[header] = ""
        primers_dict[header] = []
    else:
        seqs_dict[header] += line.strip()
        for letter in line.strip():
            primers_dict[header].append(letter)


def check_clamp(primer, reverse):
    """Function to test if the putative primer has a GC clamp"""

    if reverse == True:
        clamp = primer[-2:]
    else:
        clamp = primer[0:2]
    if clamp == "GG" or clamp == "GC" or clamp == "CG" or clamp == "CC":
        return True


def check_Tm(GC, primer):
    """Function to calculate Tm and determine if it is acceptable"""

    Tm = 64.9 + ((41 * (GC - 16.4)) / (len(primer)))
    if Tm >= low_tm and Tm <= hi_tm:
        return Tm, True
    else:
        return None, False


def check_targets(primer):
    """Function to determine how many times the putative primer occurs in the
       full assembly"""

    temp_blastn = open("temp_blastn", "r")
    targets = 0
    for line in temp_blastn:
        if line.startswith(" Identities ="):
            regex = re.match(" Identities = ([0-9]*)/([0-9]*) .*?\n", line)
            if len(primer) - int(regex.group(1)) <= min_mismatches:
                targets += 1
    temp_blastn.close()
    return targets


# Code (basically) taken from create_pcr_primers.py. It looks for all possible
# PCR primers that could amplify the good probe regions
out = open("probes.good.primers", "w")
for item in seqs_dict:
    forward = seqs_dict[item][:1000].upper()
    reverse = seqs_dict[item][-1000:].upper()
    for x in range(18, 28):
        for y in range(0, 1000 - x):
            f_primer = forward[y:y + x]
            if check_clamp(f_primer, False) == True:
                GC = f_primer.count("G") + f_primer.count("C")
                Tm, good_tm = check_Tm(GC, f_primer)
                if good_tm == True:
                    temp_query = open("temp_query", "w")
                    temp_query.write(">primer\n")
                    temp_query.write(f_primer)
                    temp_query.close()
                    blastn = "blastn -query temp_query -db %s -word_size 7 -out\
                              temp_blastn" % sys.argv[7]
                    subprocess.call(blastn, stdin=None, stdout=None,
                                            stderr=None, shell = True)
                    targets = check_targets(f_primer)
                    if targets <= max_occurrences:
                        out.write("%s_%sF Tm: %s GC: %s Occurrences in\
                                  assembly: %s\n%s\n" % (item, y, Tm,
                                  (float(GC) / y) * 100, targets, f_primer))

            r_primer = reverse[y:y + x]
            if check_clamp(r_primer, True) == True:
                GC = r_primer.count("G") + r_primer.count("C")
                Tm, good_tm = check_Tm(GC, r_primer)
                if good_tm == True:
                    temp_query = open("temp_query", "w")
                    temp_query.write(">primer\n")
                    temp_query.write(r_primer)
                    temp_query.close()
                    blastn = "blastn -query temp_query -db %s -word_size 7\
                    -out temp_blastn" % sys.argv[7]
                    subprocess.call(blastn, stdin=None, stdout=None,
                                            stderr=None, shell = True)
                    targets = check_targets(r_primer)
                    if targets <= max_occurrences:
                        r_primer_rc = Seq(r_primer)
                        r_primer_rc = r_primer_rc.reverse_complement()
                        out.write("%s_%sR Tm: %s GC: %s Occurrences in\
                                  assembly: %s\n%s\n" % (item, y - 1000,
                                  Tm, (float(GC) / x) * 100, targets, r_primer_rc))

# Some more cleanup
rm = "rm temp_query temp_blastn"
subprocess.call(rm, stdin=None, stdout=None, stderr=None, shell = True)
