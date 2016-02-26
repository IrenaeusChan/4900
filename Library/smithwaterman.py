"""
Irenaeus Chan
1/27/2015

Sequence Alignment Algorithm
Smith-Waterman Algorithm

Adapted from: http://puriney.github.io/2013/08/22/smith-waterman-algorithm/
"""

import math
import string
from phipsistore import PhiPsiStore
from atom import Atom

def readSeq (file_name):
	seq = []
	with open(file_name, "r") as stream:
		for line in stream:
			aminoAcid = line.split()
			seq.append(PhiPsiStore(int(aminoAcid[0]), aminoAcid[1], float(aminoAcid[2]), float(aminoAcid[3]), float(aminoAcid[4])))
	return seq

def printMatrix(matrix, mainSeq, targetSeq):
	print '\t'+('\t'.join(map(str,list(targetSeq))))
	i = 0
	for line in matrix:
		print mainSeq[i]+"\t"+('\t'.join(map(str,line)))
		i +=1

def align(mSeq, tSeq, mSeq_filename, tSeq_filename):
	matrix = []
	path = []
	mainSeq = ""
	targetSeq = ""

	#Convert Main and Target Sequence into Strings
	for everyAminoAcid in mSeq:
		mainSeq += everyAminoAcid.getLetter()
	for everyAminoAcid in tSeq:
		targetSeq += everyAminoAcid.getLetter()

	row = len(mainSeq)
	col = len(targetSeq)

	mainSeq = "^"+mainSeq
	targetSeq = "^"+targetSeq
	mSeq.insert(0,"^")
	tSeq.insert(0,"^")

	for i in range(row+1):
		matrix.append([0]*(col+1))		#Creates Columns for every Row
		path.append([" "]*(col+1))		#?? Not Sure

	insertDeleteValue = -1
	matchValue = 2

	for i in range(1,row+1):
		for j in range(1,col+1):
			#Penalty Map
			fromLeft = matrix[i][j-1] + insertDeleteValue
			fromTop = matrix[i-1][j] + insertDeleteValue
			if mainSeq[i]==targetSeq[j]:
				fromDiagonal = matrix[i-1][j-1] + matchValue
			else:
				fromDiagonal = matrix[i-1][j-1] + insertDeleteValue

			matrix[i][j]= max(fromLeft,fromTop,fromDiagonal)
			
			#Path Map
			if matrix[i][j]==fromLeft:
				path[i][j]="-"
			elif matrix[i][j]==fromTop:
				path[i][j] = "|"
			elif matrix[i][j] == fromDiagonal:
				path[i][j] = "M"
			else:
				pass

			if matrix[i][j]<0:
				matrix[i][j]=0

	#Trace Back
	row = len(matrix)-1
	col = len(matrix[0])-1

	while row>=0:
		maxPenaltyValue = max(matrix[row])		#Largest Number in the Row
		while col>=0:
			if (matrix[row][col] == maxPenaltyValue):
				rowPair = row
				colPair = col
				outputMain=[]
				outputTarget=[]
				while (rowPair != 0 and colPair != 0):
					if path[rowPair][colPair]=="M":
						outputMain.append(mSeq[rowPair])
						outputTarget.append(tSeq[colPair])
						rowPair -= 1
						colPair -= 1
					elif path[rowPair][colPair]=="-":
						outputMain.append("-")
						outputTarget.append(tSeq[colPair])
						colPair -= 1
					elif path[rowPair][colPair]=="|":
						outputMain.append(mSeq[rowPair])
						outputTarget.append("-")
						rowPair -= 1
					elif path[rowPair][colPair]==" ":
						if rowPair > 0:
							outputMain.append(mSeq[rowPair])
							rowPair -=1
						if colPair > 0:
							outputTarget.append(tSeq[colPair])
							colPair -= 1

				main = list(reversed(outputMain))
				target = list(reversed(outputTarget))
			col -=1
		row -=1

	with open('{0}'.format(mSeq_filename), 'w') as output:
		for everyAA in main:
			output.write(str(everyAA)+"\n")
	with open('{0}'.format(tSeq_filename), 'w') as output:
		for everyAA in target:
			output.write(str(everyAA)+"\n")