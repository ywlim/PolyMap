#!/usr/bin/python3

import sys

# this script gets de novo CDR3 amino acid sequences and links it to cell barcode/UMI in Read 2.

print("usage: denovo_cdr3_cellbc_UMI.py <fastq with cellbc/umi> <CDR3 file> <output file>")

read2_infile = open(sys.argv[1], "r")  # fastq with cellbc/UMI. could be Read1 or Read2 depending on sequencing strategy
outfile = open(sys.argv[3],"w+") # output file name
outfile.write("cdr3\tcell_barcode\tUMI") # output file header

# specify maximum number of low-quality nucleotides in seq read
# we used 5 here to maximize number of usable reads
maxee=5

# Parse cell bc/umi fastq files into a dictionary, {read header: [cell barcode, UMI]} 
header_bc_dict = {}
while True:
	line = read2_infile.readline().strip()
	if line == "":
		break
	header = line.split(" ")[0]
	header = header.replace("@","")
	line = read2_infile.readline().strip()
	cell_bc = line[0:12]
	umi = line[12:20]
	header_bc_dict[header] = [cell_bc, umi]
	read2_infile.readline()
	read2_infile.readline()
read2_infile.close()

# parse join file to get CDR3
join_infile = open(sys.argv[2],"r")
for line in join_infile:
	line = line.strip()
	if "cdr3" in line: # only parse line if CDR3 was found in the line
		line = line.split("\t")
		header = line[0]
		# get ee for forward and rev reads
		for i in range(len(line)):
			if "eef:" in line[i]:
				eef = line[i]
			if "eer:" in line[i]:
				eer = line[i]
		eef = float(eef.split("=")[-1].strip(";"))
		eer = float(eer.split("=")[-1].strip(";"))
		# only parse line if error is not >= 2
		if eef < maxee and eer < maxee:
		#print(line)
			for item in line: # get CDR3 amino acid sequence
				if "cdr3=" in item:
					#print(item)
					cdr3 = item.split(";")[0].split("=")[1]
					#print(cdr3)
			# print corresponding cell barcode/UMI info
			outfile.write("\n" + cdr3 + "\t" + header_bc_dict[header][0] + "\t" + header_bc_dict[header][1])
