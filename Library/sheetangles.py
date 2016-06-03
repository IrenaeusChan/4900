"""
Irenaeus Chan
5/11/2016

Sheet Angles

Checklist:
- Evaluate a systematic approach to geometrically organizing the information based on sequences
"""
import math
import vector

def printSheetBackbone(sheetBackbone):
	with open('testing.txt', "w") as temp:
		temp.write("x,y,z\n")
		for what in sheetBackbone:
			temp.write(str(what[0])+","+str(what[1])+","+str(what[2])+"\n")

#Code to convert the SheetList into a single List of Coordinates
def organizeSheet (sheet):
	sheetBackbone = []
	#position = []					Not sure if we need these yet...
	#seqres = []
	for AA in sheet.amino_acids:
		for atom in AA.backbone:
			sheetBackbone.append([atom.x, atom.y, atom.z])
		#position += [AA.position, AA.position, AA.position]
		#seqres += [AA.seqres, AA.seqres, AA.seqres]
	return sheetBackbone#, position, seqres

def angleCalculation(atomList):
	for i in range(0,len(atomList)-1):
		print vector.dihedralAngle(atomList[i], atomList[i+1])

def evaluateAngles (sheet, filename, atomNumber):
	topAtoms, bottomAtoms = [], []
	allAtoms = []
	sheetBackbone = organizeSheet(sheet)
	#printSheetBackbone(sheetBackbone)
	regressionVector, regressionPoint = vector.orthogonalDistanceRegression(sheetBackbone)
	print regressionVector
	print regressionPoint
	for pos, atom in enumerate(sheetBackbone):
		orthogonalLine = vector.orthogonalLineCalculation(regressionVector, atom)
		orthogonalVector = vector.orthogonalVectorCalculation(orthogonalLine, atom, regressionVector, regressionPoint)
		#allAtoms.append(orthogonalVector)
		#if pos % 2 == 0:			#if the number is even
		#	topAtoms.append(orthogonalVector)
		#else:
		#	bottomAtoms.append(orthogonalVector)

	#angleCalculation(topAtoms)
	#print "SPACE"
	#angleCalculation(bottomAtoms)
	#print "SPACE"
	#angleCalculation(allAtoms)

	filePrep = filename.split('.')