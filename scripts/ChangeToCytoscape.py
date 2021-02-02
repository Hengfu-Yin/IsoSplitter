#!/usr/bin/env python3

import os
import sys

fname='GeneCluster.txt'
if not os.path.isfile(fname):
	print("error: The file "+fname+" is not found\n")
	sys.exit()

with open('Cytoscape.txt','a') as sobj:
	sobj.write('#seq1'+'\t'+'#seq2'+'\n')
with open(fname,'r') as obj:
	for line in obj:
		if line:
			List=line.split()
			for i in range(1,len(List)):
				with open('Cytoscape.txt','a') as sobj:
					sobj.write(List[0]+'\t'+List[i]+'\n')
		