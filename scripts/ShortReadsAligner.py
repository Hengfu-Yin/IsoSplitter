#!/usr/bin/env python3

import re
import multiprocessing
import os
import sys
import subprocess
import argparse
import distutils.spawn
import time
import math

#List all sequence labels and the position of 'breaking points'
def ListName(file,dirName):
	#Open the file ('Breakpoint_out.txt') produced by ShortReadsAligner
	with open(file,"r") as obj:
		r = obj.readline()
		#Convert the format of 'Breakpoint_out.txt' to that of 's1.txt'(seqenceName breakponit1 breakpoint2 ...)
		while r:
			flag = re.match("seq",r)
			if flag:
				#s = r.split()[2].replace("_","/")
				s = r.split()[2]
				r = obj.readline()
				rr = r.split()
				s = s+" "+rr[1][1:-1]
				i = 3
				while i < len(rr):
					s = s+" "+rr[i][:-1]
					i = i + 2
				with open(dirName+"/s1.txt","a") as sobj:
					sobj.write(s+"\n")
			r = obj.readline()

#Extract sequences that include splitting sites for short reads alignment using bowtie2
def FindFaToAlign(file,dirName):
	#Save 's1.txt' to dictTemp ({seqenceName:[breakpoint1 breakpoint2 ...]})
	with open(dirName+"/s1.txt","r") as obj:
		dictTemp = dict()
		for i in obj:
			listT = i.split()
			dictTemp[listT[0]] = listT
	#Open the input file 'longReadsFile'
	with open(file,"r") as obj:
		r = obj.readline()
		while not re.match('>', r):
			r = obj.readline()
			if not r:
				break
		#For each sequence, 135bp is taken before and after breakpoints.
		while r:
			strR = r.split(maxsplit=1)[0][1:].strip()
			flag = 0
			s = ""
			if strR in dictTemp:
				flag = 1
				listBreak = dictTemp[strR]
				listTemp = dictTemp[strR][1:]	
				del dictTemp[strR]
			r = obj.readline()
			while not re.match('>', r):
				if flag == 1:
					s = s + r.strip()			
				r = obj.readline()
				if not r:
					break
			if flag == 1:
				#Processing of the first breakpoint
				Dict = {}
				listTemp = [int(i) for i in listTemp]
				listTemp.sort()
				if listTemp[0] - 135 <= 0:
					start = 0
					num = listTemp[0]			
				else:
					start = listTemp[0] - 135
					num = 135
				end = listTemp[0] + 135		
				if end >= len(s):
					end = -1
				listBreak[1] = num
				Dict[num] = listTemp[0]
				if end == -1:
					ss = s[start:]
				else:
					ss = s[start:end]
				#If there are other breakpoints, discuss them in categories.
				if len(listTemp) > 1:
					i = 1
					ss = ""
					while(i < len(listTemp)):
						if end != -1:
							if listTemp[i] - listTemp[i-1] <= 270 and i != len(listTemp) - 1:
								end = listTemp[i] + 135
								if end >= len(s):
									end = -1
								num = listBreak[i] + listTemp[i] - listTemp[i-1]
							elif listTemp[i] - listTemp[i-1] <= 270 and i == len(listTemp) - 1:
								end = listTemp[i] + 135
								if end >= len(s):
									ss = ss + s[start:]
								else:
									ss = ss + s[start:end]
								num = listBreak[i] + listTemp[i] - listTemp[i-1]
							elif listTemp[i] - listTemp[i-1] > 270 and i != len(listTemp) - 1:
								ss = ss + s[start:end]
								start = listTemp[i] - 135
								end = listTemp[i] + 135
								if end >= len(s):
									end = -1
								num = listBreak[i] + 270
							else:
								ss = ss + s[start:end]
								start = listTemp[i] - 135
								end = listTemp[i] + 135
								if end >= len(s):
									ss = ss + s[start:]
								else:
									ss = ss + s[start:end]
								num = listBreak[i] + 270				
						elif end == -1 and i != len(listTemp) - 1:
							num = listBreak[i] + listTemp[i] - listTemp[i-1]
						else:
							ss = ss + s[start:]
							num = listBreak[i] + listTemp[i] - listTemp[i-1]
						listBreak[i+1] = num
						Dict[num] = listTemp[i]
						i = i + 1	
				#According to the interception, the position of breakpoints are changed to a new position.
				#Write to 's3.txt', format is 'seqenceName breakpoint1 breakpoint2 ...'					
				with open(dirName+"/s3.txt","a") as subobj:
					sr = listBreak[0]
					for i in listBreak[1:]:
						sr = sr+" "+str(i)
					subobj.write(sr+"\n")	
				#Splice the intercepted fragments into one line and write them to file 's2.txt'
				#format is 'first line: >seqence name  second line: seqence'			
				with open(dirName+"/s2.txt","a") as suobj:
					suobj.write(">"+listBreak[0]+"\n")
					suobj.write(ss+'\n')
				#Correspondence between the old and new positions of breakpoints in each sequence
				#format is 'seqenceName {new : old}'
				with open(dirName+"/s4.txt","a") as ssobj:
					ssobj.write(listBreak[0]+"  "+str(Dict)+"\n")

#Divide a file (s2.txt) into n parts
def SplitFileByfa(num,dirName):
	#Count the rows of 's2.txt'
	count = subprocess.Popen("wc -l "+dirName+"/s2.txt",stdout=subprocess.PIPE,shell=True)
	hang = int(count.stdout.read().decode().split()[0]) / 2
	with open(dirName+"/s2.txt","r") as obj:
		r = obj.readline()
		k = 1
		j = 1
		while r:
			s = r
			r = obj.readline()
			if not r:
				break
			s = s + r
			#Divide the file into equal parts
			if j < num and k <= hang / num :
				k = k + 1
			elif j < num and k > hang / num :
				k = 2
				j = j + 1
			'''
			#Divide the file into proportions
			tnum = math.ceil(num / 4)
			if j <= tnum and k <= hang / (tnum * 2):
				k = k + 1
			elif tnum < j <= tnum * 2 and k <= hang / num:
				k = k + 1
			elif tnum * 2 < j < num and k <= hang / (num * 2):
				k = k + 1
			elif (j <= tnum and k > hang / (tnum * 2)) or (tnum < j <= tnum * 2 and k > hang / num) or (tnum * 2 < j < num and k > hang / (num * 2)):
				k = 2
				j = j + 1
			'''
			#The each file is named 's2[*].txt'
			name = dirName+"/s2"+str(j)+".txt"
			with open(name,"a") as sobj:
				sobj.write(s)
			r = obj.readline()

#Pick data produced by bowtie2: (>=15)M(>=15)S, write to 'ShortReadsMapped.sam'
def FindMatched(r,length,dirName):
	listR = r.split()
	num = re.search("M",listR[5]).start()
	numberM = int(listR[5][:num])
	if len(listR[5])-1==num:
		for i in length:
			if int(listR[3])<=int(i)<=int(listR[3])+numberM:
				with open(dirName+"/ACMMapped.sam","a") as sobj:
					sobj.write(r)
				break
	elif numberM >= 15 and int(listR[5][num+1:-1]) >= 15:	
		for i in length:
			if int(listR[3]) + numberM == int(i) :
				with open(dirName+"/ShortReadsMapped.sam","a") as sobj:
					sobj.write(r)
				break		

#Pick data produced by bowtie2: be satisfied with 's3.txt'
def CountSites(dirName):
	#Save 's3.txt' to dictTemp ({seqenceName:[breakpoint1 breakpoint2 ...]})
	with open(dirName+"/s3.txt","r") as obj:
		dictTemp = dict()
		for i in obj:
			listT = i.split()
			dictTemp[listT[0]] = listT[1:]
	#Open the merged file 's1.sam'
	with open(dirName+"/s1.sam","r") as obj:
		r = obj.readline()	
		while r:
			listR1 = r.split()
			if listR1[2] in dictTemp:
				length = dictTemp[listR1[2]]	
				del dictTemp[listR1[2]]									
			FindMatched(r,length,dirName)
			r = obj.readline()
			if not r:
				break
			listR2 = r.split()
			while listR1[2] == listR2[2]:
				FindMatched(r,length,dirName)
				r = obj.readline()
				listR2 = r.split()
				if not r:
					break

#Count breakpoints on matched reference sequences and their occurrence times
def AllMatchedPoint(dirName):
	if not os.path.exists(dirName+'/ShortReadsMapped.sam'):
		print("There is no meet the condition of split site !")
		sys.exit()
	with open(dirName+"/ShortReadsMapped.sam",'r') as file_obj :
		r = file_obj.readline()
		Dict = {}
		while r:
			Dict.clear()
			first = r
			list1 = r.split("\t")
			key = int(list1[3]) + int(list1[5][0:re.search("M",list1[5]).start()])
			Dict[key] = 1
			r = file_obj.readline()
			if not r:
				with open(dirName+"/s5.txt","a") as obj:
					obj.write(first.split("\t")[2]+"\t"+str(Dict)+"\n")
				break
			list2 = r.split("\t")
			while list1[2] == list2[2]:
				key = int(list2[3]) + int(list2[5][0:re.search("M",list2[5]).start()])
				if  key in Dict:
					Dict[key] = Dict[key] + 1
				else:
					Dict[key] = 1
				r = file_obj.readline()
				list2 = r.split("\t")
				if not r:
					break
			#Write the statistical results to file 's5.txt'
			#format is 'seqenceName {breakpoint : occurrence number}'
			with open(dirName+"/s5.txt","a") as obj:
				obj.write(first.split("\t")[2]+"\t"+str(Dict)+"\n")


def ACM_AllMatchedPoint(dirName):
	# if not os.path.exists(dirName+'/ACMMapped.sam'):
	# 	print("There is no meet the condition of ACM split site !")
	# 	sys.exit()
	with open(dirName+"/s3.txt","r") as obj:
		dictTemp = dict()
		for i in obj:
			listT = i.split()
			dictTemp[listT[0]] = listT[1:]
	with open(dirName+"/ACMMapped.sam",'r') as file_obj :
		r = file_obj.readline()
		Dict = {}
		while r:
			Dict.clear()
			first = r
			list1 = r.split("\t")
			length = dictTemp[list1[2]]
			for i in length:
				if int(list1[3])<=int(i)<=(int(list1[3]) + int(list1[5][0:re.search("M",list1[5]).start()])):
			#key = int(list1[3]) + int(list1[5][0:re.search("M",list1[5]).start()])
					Dict[int(i)] = 1
					break
			r = file_obj.readline()
			if not r:
				with open(dirName+"/s6.txt","a") as obj:
					obj.write(first.split("\t")[2]+"\t"+str(Dict)+"\n")
				break
			list2 = r.split("\t")
			while list1[2] == list2[2]:
				#key = int(list2[3]) + int(list2[5][0:re.search("M",list2[5]).start()])
				for i in length:
					if int(list2[3])<=int(i)<=(int(list2[3]) + int(list1[5][0:re.search("M",list1[5]).start()])) :
						if  int(i) in Dict:
							Dict[int(i)] = Dict[int(i)] + 1
						else:
							Dict[int(i)] = 1
						break
				r = file_obj.readline()
				list2 = r.split("\t")
				if not r:
					break
			#Write the statistical results to file 's5.txt'
			#format is 'seqenceName {breakpoint : occurrence number}'
			with open(dirName+"/s6.txt","a") as obj:
				obj.write(first.split("\t")[2]+"\t"+str(Dict)+"\n")

#Change the modified breakpoint position into the original predicted value
def ChangeBreakPoint(dirName,inputFile,outputFile):
	#Save 's4.txt' to dictTemp ({seqenceName:'{new1:old1}...'})
	with open(dirName+"/s4.txt","r") as obj:
		dictTemp = dict()
		for i in obj:
			listT = i.split(maxsplit=1)
			dictTemp[listT[0]] = listT[1]
	with open(dirName+"/"+inputFile,"r") as obj:
		r = obj.readline()
		while r:
			flag = 0
			dictW = {}
			strList = r.split(maxsplit=1)
			if strList[0] in dictTemp:
				dictR = strList[1]
				dictR = eval(dictR)
				dictI = dictTemp[strList[0]]
				dictI = eval(dictI)
				for k in dictR.keys():
					if k in dictI:
						dictW[dictI[k]] = dictR[k]
				flag = 1
			if flag == 1:
				#format is 'seqenceName {breakpoint : occurrence number}'
				with open(dirName+"/"+outputFile,"a") as suobj:
					if dictW:
						suobj.write(strList[0]+"\t"+str(dictW)+"\n")
			r = obj.readline()

def ACM_count(dirName):
	with open(dirName+"/s7.txt",'r') as obj:
		for line in obj:
			if line:
				count=line.count(':')
				list=line.split()
				sum=0
				for i in range(1,len(list)):
					if i%2==0:
						sum=sum+int(list[i][:-1])
				with open(dirName+"/Average_counts_per.txt",'a') as sobj:
					sobj.write(list[0]+'\t'+str(round(sum/count,3))+'\n')

#check file whether exit
def CheckFileIsExit(exec_dir,fileName):
	if not os.path.exists(exec_dir+r'/'+fileName):
		print("error:The file "+ fileName +" is not exit!")
		sys.exit()

#check the input file whether exit
def CheckSeqfile(parser,fname):
	if not os.path.isfile(fname):
		parser.error("file "+fname+" is not found\n")
	else:
		return fname

#Remove temporary files produced by bowtie2-build,and the files (s2[*].txt) splited with -n
def DeleteTempFile(num,dirName):
    for i in range(1,num+1):
        subprocess.run("rm -rf "+dirName+"/index"+str(i)+".1.bt2",shell=True)
        subprocess.run("rm -rf "+dirName+"/index"+str(i)+".2.bt2",shell=True)
        subprocess.run("rm -rf "+dirName+"/index"+str(i)+".3.bt2",shell=True)
        subprocess.run("rm -rf "+dirName+"/index"+str(i)+".4.bt2",shell=True)
        subprocess.run("rm -rf "+dirName+"/index"+str(i)+".rev.1.bt2",shell=True)
        subprocess.run("rm -rf "+dirName+"/index"+str(i)+".rev.2.bt2",shell=True) 
        subprocess.run("rm -rf "+dirName+"/s2"+str(i)+".txt",shell=True)

#Check your environment whether installed bowtie2-build and bowtie2, and them is in PATH.
def check_Bowtie2():
    ret1 = distutils.spawn.find_executable("bowtie2-build")
    ret2 = distutils.spawn.find_executable("bowtie2")
    if not ret1:
    	print("error:The bowtie2-build is required but not installed or not in the PATH!")
    	sys.exit()
    if not ret2:
    	print("error:The bowtie2 is required but not installed or not in the PATH!")
    	sys.exit()

#Define command line parameters		
def CommandLineArgs(maxThread):
	useage="ShortReadsAligner [options] longReadsFile shortReadsFile BreakPointFile"
	des="This program is short reads mapping.The total number of CPUs you can use is "+maxThread+",so the option of -t multiply by -n should less than or equal to "+maxThread
	parser = argparse.ArgumentParser(usage=useage,description=des,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('longReadsFile',help="the file is your input of long-read transcriptome sequences in the first step",type=lambda x: CheckSeqfile(parser,x))
	parser.add_argument('shortReadsFile',help="second generation transcriptome sequence",type=lambda x: CheckSeqfile(parser,x))
	parser.add_argument('BreakPointFile',help="the file is output of first step",type=lambda x: CheckSeqfile(parser,x))
	parser.add_argument('-v','--version', action='version', version='ShortReadsAligner.py 1.0')
	parser.add_argument('-t',dest='thread',metavar='<integer>',help='the number of threads',type=int,default=4)
	parser.add_argument('-n',dest='fileNum',metavar='<integer>',help='the number of split files',type=int,default=3)
	parser.add_argument('-q','--fastq',action='store_true',help='If this switch is present, the shortReadsFile must be in fastq format(default is in fasta format)', default=False)
	args = parser.parse_args()
	return args

#Call bowtie2 with the shortReads file in fasta format
def Bowtie2ByFa(i,thread,rawData,dirName):
	subprocess.run("bowtie2 --local --no-unal --no-hd --quiet -a -p "+thread+" -x "+dirName+"/index"+str(i)+" -f -U "+rawData+" | grep -Ew '([0-9]{2,3}M[0-9]{2,3}S|[0-9]{2,3}M)' | grep --color=auto 'NM:i:0' > "+dirName+"/s"+str(i)+".sam ",shell=True)

#Call bowtie2 with the shortReads file in fastq format
def Bowtie2ByFq(i,thread,rawData,dirName):
	subprocess.run("bowtie2 --local --no-unal --no-hd --quiet -a -p "+thread+" -x "+dirName+"/index"+str(i)+" -U "+rawData+" | grep -Ew '([0-9]{2,3}M[0-9]{2,3}S|[0-9]{2,3}M)' | grep --color=auto 'NM:i:0' > "+dirName+"/s"+str(i)+".sam ",shell=True)

#The main process of script
def Main():
	#check dependence bowtie2 and bowtie2-build.
	check_Bowtie2()
	#command-line argument parsing.
	maxThread = multiprocessing.cpu_count()
	args = CommandLineArgs(str(maxThread))
	if args.thread < 1 or args.thread > maxThread:
		print("error:The thread's number should between 1 and "+str(maxThread))
		sys.exit()
	else:
		thread = str(args.thread)
	if args.fileNum < 1 or args.fileNum > maxThread:
		print("error:The number of split files should between 1 and "+str(maxThread))
		sys.exit()
	else:
		num = args.fileNum
	rawData = args.shortReadsFile
	poolnum = num*args.thread
	if poolnum > maxThread:
		print('error:The number of thread or split files is too big, please keep \"t*n<='+str(maxThread)+'\"!')
	#Find the directory where the script is running and create the output directory (Sout[time]) under that directory
	exec_dir = os.path.dirname(os.path.realpath(__file__))
	dirName = "Sout"+str(int(time.time()))
	#subprocess.run("mkdir "+dirName,shell=True)
	if os.path.isdir(exec_dir+"/"+dirName):
		print("error:The output directory"+dirName+"is exit, please delete it !")
		sys.exit()
	os.mkdir(dirName)
	#List the sequence names of all breakpoints to be aligned
	ListName(args.BreakPointFile,dirName)
	#Sequence names are used to extract the parts of the genome to be aligned
	FindFaToAlign(args.longReadsFile,dirName)
	#Divide files into multiple files for alignment
	SplitFileByfa(num,dirName)
	#Call the core program bowtie2-build
	for i in range(1,num+1):
		name = dirName+"/s2"+str(i)+".txt "
		index = dirName+"/index"+str(i)
		subprocess.run("bowtie2-build -q "+name+index,shell=True)
	#Call the core program bowtie2
	pool = multiprocessing.Pool(poolnum)
	if args.fastq == False:
		for i in range(1,num + 1):
			pool.apply_async(Bowtie2ByFa,args=(i,thread,rawData,dirName,))
	else:
		for i in range(1,num+1):
			pool.apply_async(Bowtie2ByFq,args=(i,thread,rawData,dirName,))
	pool.close()
	pool.join()	
	#Remove temporary files
	DeleteTempFile(num,dirName)
	#Merge the result files generated by each bowtie2 into 's1.sam'
	for i in range(2,num+1):
		subprocess.run("cat "+dirName+"/s"+str(i)+".sam >> "+dirName+"/s1.sam",shell=True)
		subprocess.run("rm -rf "+dirName+"/s"+str(i)+".sam",shell=True)
	subprocess.run("sort -t '\t' -k 3 "+dirName+"/s1.sam -o "+dirName+"/s1.sam",shell=True)
	#Screen qualified records from 's1.sam'
	CountSites(dirName)
	#Count the qualified breakpoints
	AllMatchedPoint(dirName)
	ChangeBreakPoint(dirName,"s5.txt","JunctionReadsCount.txt")
	if not os.path.exists(dirName+'/ACMMapped.sam'):
		print("There is no meet the condition of ACM split site !")
		subprocess.run("rm -rf "+dirName+"/s1.sam",shell=True)
		for i in range(1,6):
			subprocess.run("rm -rf "+dirName+"/s"+str(i)+".txt",shell=True)
	else:
		ACM_AllMatchedPoint(dirName)
		ChangeBreakPoint(dirName,"s6.txt","s7.txt")
		ACM_count(dirName)
		#Remove the temporary output file 's1.sam'
		subprocess.run("rm -rf "+dirName+"/s1.sam",shell=True)
		#Remove the temporary file 's[*].txt'
		for i in range(1,8):
			subprocess.run("rm -rf "+dirName+"/s"+str(i)+".txt",shell=True)

if __name__ == '__main__':
	Main()
