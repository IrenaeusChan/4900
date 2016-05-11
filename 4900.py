import sys
import string
import glob
import os
sys.path.append(os.path.realpath("Library"))
from phipsi import calculatePhiPsi
import protein
import smithwaterman
from newCoord import newCalculate
import extra
from helixangles import evaluateAngles

MAX_VALUE = 100
path = os.getcwd()

#Prepares the protein for calculations
def prepare(file_name):
	p = protein.buildProtein(file_name)
	cenx, ceny, cenz = protein.weightedAverage(p)
	center = [cenx, ceny, cenz]
	return p, center

#Function for computing all the PDB Files in a Directory
def computeAll(file_name):
	only_100 = 0
	for filename in glob.glob(os.path.join(path, '*.pdb')):
		if (only_100 == MAX_VALUE):
			break
		else:
			print "{0} files completed...".format(only_100)
			p, center = prepare(filename)
			calculatePhiPsi(p, center, filename)
			protein.relativeToCenter(p, center)
			only_100+=1

#Checks to see if the file exists and is formatted properly
def check(file1, file2):
	if os.path.isfile(file1) == True and file1.endswith(".txt"):
		mainSeq = smithwaterman.readSeq(file1)
	else:
		if os.path.isfile(file1[:-4]) == True:
			p, center = prepare(file1[:-4])
			calculatePhiPsi(p, center, file1[:-4])
			mainSeq = smithwaterman.readSeq(file1)
		else:
			print "\nERROR: The FIRST File you provided is incorrect"
			mainSeq = False
	if os.path.isfile(file2) == True  and file2.endswith(".txt"):
		targetSeq = smithwaterman.readSeq(file1)
	else:
		if os.path.isfile(file2[:-4]) == True:
			p, center = prepare(file2[:-4])
			calculatePhiPsi(p, center, file2[:-4])
			targetSeq = smithwaterman.readSeq(file2)
		else:
			print "\nERROR: The SECOND File you provided is incorrect"
			targetSeq = False
	return mainSeq, targetSeq

#Checks normal formatting
def format(filename):
	if os.path.isfile(filename) == True and filename.endswith(".pdb"):
		return True
	else:
		return False

if __name__ == '__main__':
	if len(sys.argv) < 2:
		print "\nERROR: No file was provided"
		sys.exit(1)
	#Does a single file Ramachandran Plot
	elif format(sys.argv[1]):
		print "\nComputing Using File: {0}".format(sys.argv[1])
		p, center = prepare(sys.argv[1])
		calculatePhiPsi(p, center, sys.argv[1])
	#Multiple file Ramachandran Plot
	elif (sys.argv[1] == "all" and len(sys.argv) > 2):
		path += sys.argv[2]
		print "\nComputing All Files."
		if (os.path.isdir(path) == True):
			computeAll(sys.argv[2])
		else:
			print "\nERROR: This folder does not exist"
			print "ERROR: Please make sure directory format is correct - \Directory"
    #New Geometric Calculations
	elif (sys.argv[1] == "new" and len(sys.argv) > 2 and format(sys.argv[2])):
		p = protein.buildProtein(sys.argv[2])
		coilList = extra.buildCoil(sys.argv[2], p)
		if (len(coilList)<1):
			print "\nThere were no COILS in this protein"
		else:
			for coil in coilList:
				newCalculate(coil, sys.argv[2])
	#Utilizes Smith-Waterman Algorithm to Locally Align two proteins
	elif (sys.argv[1] == "align" and len(sys.argv) > 3):
		print "\nAligning your two files."
		mainSeq, targetSeq = check(sys.argv[2], sys.argv[3])
		if (mainSeq == False or targetSeq == False):
			print "\nERROR: Please make sure the file type is correct - filename.pdb.txt"
		else:
			smithwaterman.align(mainSeq, targetSeq, sys.argv[2], sys.argv[3])
	#If the user wants to calculate the PhiPsi angles for ONLY the Helix or SHEETS	
	elif (len(sys.argv) > 2 and (sys.argv[1] == "Helix" or sys.argv[1] == "helix") and format(sys.argv[2])):
		#We need to read the PDB file once again to get more information about positions
		p, center = prepare(sys.argv[2])
		helixList = extra.buildHelix(sys.argv[2])
		if (len(helixList)<1):
			print "\nThere were no HELICES in this protein"
		else:
			print "\nComputing with only the HELICES in this Protein"
			for helix in helixList:
				calculatePhiPsi(helix, center, sys.argv[2])
	#Same for HELIX
	elif (len(sys.argv) > 2 and (sys.argv[1] == "Sheets" or sys.argv[1] == "sheets") and format(sys.argv[2])):
		p, center = prepare(sys.argv[2])
		sheetList = extra.buildSheet(sys.argv[2], p)
		if (len(sheetList)<1):
			print "\nThere were no SHEETS in this protein"
		else:
			print "\nComputing with only the SHEETS in this Protein"
			for sheet in sheetList:
				calculatePhiPsi(sheet, center, sys.argv[2])
	elif (len(sys.argv) > 2 and (sys.argv[1] == "coil" or sys.argv[1] == "Coil") and format(sys.argv[2])):
		p, center = prepare(sys.argv[2])
		coilList = extra.buildCoil(sys.argv[2], p)
		#Run with CoilList. Does not generate BEFORE and AFTER yet. But should.
		if (len(coilList)<1):
			print "\nThere were no COILS in this protein"
		else:
			print "\nComputing with only the COILS in this Protein"
			for coil in coilList:
				calculatePhiPsi(coil, center, sys.argv[2])
	elif (len(sys.argv) > 2 and (sys.argv[1] == "Angles" or sys.argv[1] == "angles") and format(sys.argv[2])):
		helixList = extra.buildHelix(sys.argv[2], protein.buildProtein(sys.argv[2]))
		for helix in helixList:
			if helix.helixType == 1:
				evaluateAngles(helix, sys.argv[2])
				sys.exit(1)
	elif (len(sys.argv) > 2 and (sys.argv[1] == "Center" or sys.argv[1] == "center") and format(sys.argv[2])):
		p, center = prepare(sys.argv[2])
		print "\nProtein's Center (X, Y, Z): " + str(protein.weightedAverage(p))
	elif (len(sys.argv) > 2 and (sys.argv[1] == "Distance" or sys.argv[1] == "distance") and format(sys.argv[2])):
		p, center = prepare(sys.argv[2])
		protein.relativeToCenter(p, center)
	else:
		print "\nERROR: File type is incorrect or does not exist"
		sys.exit(1)

	print "\nOther Options Available:"
	print "\tCenter\t\t- Calculates Center of Protein"
	print "\tDistance\t- Amino Acid Residue from Center"
	print "\tHelix\t\t- Applies Ramachandran Plots to JUST Helix"
	print "\tSheets\t\t- Applies Ramachandran Plots to JUST Sheets"
	print "\tCoil\t\t- Applies Ramachandran Plots to JUST Coils"
	print ""
	print "If you would like to calculate entire directories of PDB Files:"
	print "\tall folder"
	print "To align two PDB files that have existing PhiPsi Calculations already done"
	print "using Smith-Waterman Algorithm:"
	print "\talign file_name1 file_name2"
	print ""
	print "In order to use the New Geometric Analysis on a single PDB file (In progress):"
	print "\t new file_name"