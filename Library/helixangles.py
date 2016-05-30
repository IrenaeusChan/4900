"""
Irenaeus Chan
5/11/2016

Helix Model Check

Checklist:
- Evaluate a systematic approach to geometrically organizing the information based on sequences
"""

import string
import sys
from protein import buildProtein
import vector
import extra
import aminoacid
import atom

#Code to store turn (Large number of small segments)
def organizeStructure (helix):
	helixBackbone = []
	position = []
	seqres = []
	for AA in helix.amino_acids:
		helixBackbone += AA.backbone
		position += [AA.position, AA.position, AA.position]
		seqres += [AA.seqres, AA.seqres, AA.seqres]
	return helixBackbone, position, seqres
#1 Turn = 13 Atoms
#turn large list of 1 helix into smaller list of 13, run through list...

#Code to evaluate angles from each atom to it's next respective atom
def evaluateAngles (secondaryStructure, filename, atomNumber):
	structureBackbone, position, seqres = organizeStructure(secondaryStructure)
	filePrep = filename.split('.')
	for i in range(atomNumber,len(structureBackbone)+1):
		#We want to initialize these variables for EACH segment of atoms
		totalx, totaly, totalz = 0, 0 ,0
		Nvector2, CAvector2, Cvector2 = 0, 0, 0
		secondNpos, secondCApos, secondCpos = 0, 0, 0

		#Calculates the Average Center of the Atom segments
		for atom in structureBackbone[i-atomNumber:i]:
		 	totalx += atom.x
		 	totaly += atom.y
		 	totalz += atom.z
		totalx /= atomNumber
		totaly /= atomNumber
		totalz /= atomNumber
		center = [totalx, totaly, totalz]

		with open('{0}N.txt'.format(filePrep[0]), "a") as N, open('{0}CA.txt'.format(filePrep[0]), "a") as Ca, open('{0}C.txt'.format(filePrep[0]), "a") as C:
			#For each atom in the 13 Atom Turn, we want to calculate the Angles
			#There is a bug where if the list of atoms is less than a single "Turn" it won't do the angle calculations
			for count, atom in enumerate(structureBackbone[i-atomNumber:i]):
				if (atom.atom == "N"):
					try:
						#The first vector should be the last N that was found, compared to...
						Nvector1 = Nvector2
						firstNpos = secondNpos
						# the next N that will be in the sequence
						Nvector2 = vector.vectorCalculation(center, [atom.x, atom.y, atom.z])
						secondNpos = position[i-atomNumber:i][count]
						seqresN = seqres[i-atomNumber:i][count]
						#As long asiceber the we have at least found TWO Nitrogens, we can calculate their angles
						if Nvector1 != 0:
							N.write("N " + str(vector.dihedralAngle(Nvector1, Nvector2)) + " " + 
								str(firstNpos) + " " + str(secondNpos) + " " + str(seqresN) +"\n")
					except IndexError:
						print "Error Detected in PDB File: N in Residue {0}".format(secondNpos)
				elif (atom.atom == "CA"):
					try:
						#The first vector should be the last Ca that was found, compared to...
						CAvector1 = CAvector2
						firstCApos = secondCApos
						# the next N that will be in the sequence
						CAvector2 = vector.vectorCalculation(center, [atom.x, atom.y, atom.z])
						secondCApos = position[i-atomNumber:i][count]
						seqresCA = seqres[i-atomNumber:i][count]
						#As long as the we have at least found TWO alpha-Carbons, we can calculate their angles
						if CAvector1 != 0:
							Ca.write("Ca " + str(vector.dihedralAngle(CAvector1, CAvector2)) + " " + 
								str(firstCApos) + " " + str(secondCApos) + " " + str(seqresCA) + "\n")
					except IndexError:
						print "Error Detected in PDB File: Ca in Residue {0}".format(secondCApos)
				elif (atom.atom == "C"):
					try:
						#The first vector should be the last C that was found, compared to...
						Cvector1 = Cvector2
						firstCpos = secondCpos
						# the next N that will be in the sequence
						Cvector2 = vector.vectorCalculation(center, [atom.x, atom.y, atom.z])
						secondCpos = position[i-atomNumber:i][count]
						seqresC = seqres[i-atomNumber:i][count]
						#As long as the we have at least found TWO Carbons, we can calculate their angles
						if Cvector1 != 0:
							C.write("C " + str(vector.dihedralAngle(Cvector1, Cvector2)) + " " + 
								str(firstCpos) + " " + str(secondCpos) + " " + str(seqresC) + "\n")
					except IndexError:
						print "Error Detected in PDB File: C in Residue {0}".format(secondCpos)

def transitionAngle(atomFile, infoFile, secondaryType, helixOrsheet):
	atomList = []
	outPrep = atomFile.split(".")
	with open(atomFile, "r") as atoms, open(infoFile, "r") as stream, open("{0}parsed.txt".format(outPrep[0]), "a") as out:
		for line in atoms:
			atomList.append(line)
		for line in stream:
			#Reads from the file the information regarding the HELIX structures
			if helixOrsheet == "helix":
				if (line[0:5] == "HELIX"):
					start = int(line[21:25])
					stop = int(line[33:37])
					seqres = str(line[19:20])
					helixType = int(line[39:40])
					if helixType != secondaryType:
						continue
					for a in atomList:
						atom = a.split()		#[0]-ATOM [1]-ANGLE [2]-POS1 [3]-POS2 [4]-SEQRES
						if (int(atom[2]) == (start-1) and int(atom[3]) == start and atom[4] == seqres):
							out.write(atom[0] + " " + atom[1] + " " + atom[2] + " " + atom[3] + "\n")
						if (int(atom[2]) == stop and int(atom[3]) == (stop+1) and atom[4] == seqres):
							out.write(atom[0] + " " + atom[1] + " " + atom[2] + " " + atom[3] + "\n")
			elif helixOrsheet == "sheet":
				if (line[0:5] == "SHEET"):
					start = int(line[22:26])
					stop = int(line[33:37])
					seqres = str(line[21:22])
					sheetType = int(line[38:40])
					if sheetType != secondaryType:
						continue
					for a in atomList:
						atom = a.split()		#[0]-ATOM [1]-ANGLE [2]-POS1 [3]-POS2 [4]-SEQRES
						if (int(atom[2]) == (start-1) and int(atom[3]) == start and atom[4] == seqres):
							out.write(atom[0] + " " + atom[1] + " " + atom[2] + " " + atom[3] + "\n")
						if (int(atom[2]) == stop and int(atom[3]) == (stop+1) and atom[4] == seqres):
							out.write(atom[0] + " " + atom[1] + " " + atom[2] + " " + atom[3] + "\n")

	return 0

#Code to compare relative angles to each segment (Variation in curvature)
def relativeAngle():
	return 0