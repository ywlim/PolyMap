#!/usr/bin/python3

import sys
from scipy.spatial import distance

# this scripts looks for barcodes in Read 1 to match to the provided barcode/antigen variant file. 2 mismatches are allowed
# links barcodes from Read 1 to the UMI/cell barcode in Read 2

print("usage: dropseq_custombc_cellbc_UMI.py <path to R1.fastq with ag bc> <path to R2.fastq with umi/cell bc> <path to tab-delimited barcode/variant file> <output file name>")

read1_infile = open(sys.argv[1],"r")
read2_infile = open(sys.argv[2], "r")
bc_infile = open(sys.argv[3],"r")
outfile = open(sys.argv[4],"w+")
outfile.write("ag_variant\tcell_barcode\tUMI")

# read custom barcode file. dictionary {antigen barcode: corresponding antigen}
# this script is set up for 20bp. must change indices in read_ag_bc (line 39) if ag bc is different len
ag_bc_dict = {}
for line in bc_infile:
	line = line.strip().split("\t")
	ag_bc_dict[line[1]] = line[0]
ag_bcs = list(ag_bc_dict.keys())
ag_bcs_split = []
for bc in ag_bcs:
	ag_bcs_split.append(list(bc))
bc_infile.close()

# Parse Read1/Read2 fastq files and merge
ag_in_read = False # default no antigen in read
count = 0
while True:
	# Read 1 fastq - get antigen barcode and corresponding antigen using dict barcode:antigen
	line = read1_infile.readline().strip() # header
	if line == "": # stop if EOF
		break
	header = line[1:]
	line = read1_infile.readline().strip() # read
	read_ag_bc = list(line[56:76]) # 20bp
	# record ag variant in hamming distance <= 2 to the known ag barcode list
	for bc in ag_bcs_split:
		hamming_dist = distance.hamming(read_ag_bc, bc)*len(read_ag_bc)
		if hamming_dist <= 2:
			ag_in_read = True
			ag_variant = ag_bc_dict["".join(bc)]
	line = read1_infile.readline() # +
	line = read1_infile.readline() # qual
	
	# Read 2 fastq - get cell barcode and UMI barcode
	line = read2_infile.readline().strip() # header
	line = read2_infile.readline().strip() # read
	cell_bc = line[0:12]
	umi = line[12:]
	read2_infile.readline() # +
	read2_infile.readline() # qual
	
	if ag_in_read: # if antigen is found in read, write read header, antigen, cell barcode, and UMI to output file
		outfile.write("\n" + ag_variant + "\t" + cell_bc + "\t" + umi)
	else:
		read_ag_bc = "".join(read_ag_bc)
	count += 1
	ag_in_read = False

print(count)
	
outfile.close()
read1_infile.close()
read2_infile.close()
