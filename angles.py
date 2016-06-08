import sys
import glob
import os
sys.path.append(os.path.realpath("Library"))
from protein import buildProtein
from extra import buildHelix, buildSheet
import helixangles
import sheetangles
import time
import vector

MAX_VALUE = 100 		#Used for cases where there are... like... 1 million files
path = os.getcwd()

#Checks normal formatting
def format(filename):
	if os.path.isfile(filename) == True and filename.endswith(".pdb"):
		return True
	else:
		return False

#This sets the amount of atoms that will be checked PER turn
def setSecondaryType(secondaryType):
	if secondaryType == 1:
		return 11		
		#Alpha Helix is 3.6 Amino Acids with 13 Atoms involved in the turn
		# not counting the Hydrogen and the Oxygen that form the hydrogen bond, it is 11
		# backbone atoms that make up a single turn
	elif secondaryType == 5:
		return 8 		
		#Pi Helix is 3 Amino Acids with 10 Atoms involved in the turn (3/10)
		# again, not counting the Hydrogen and the Oxygen that form the hydrogen bond, it is
		# 8 backbone atoms that make up a single turn
	else:
		return 1

def helixOrSheet(filename, secondaryType, helixOrsheet):
	if (helixOrsheet == "helix"):
		helixList = buildHelix(filename, buildProtein(filename))
		atomNumber = setSecondaryType(secondaryType)
		for helix in helixList:
			if helix.helixType == secondaryType:
				helixangles.evaluateAngles(helix, filename, atomNumber)
	elif (helixOrsheet == "sheet"):
		sheetList = buildSheet(filename, buildProtein(filename))
		for sheet in sheetList:
			if sheet.sheetType == secondaryType:
				#vector.orthogonalVectorCalculation([1,3,2], [-6,-3,2], [-1,-5,2])
				sheetangles.evaluateAngles(sheet, sys.argv[2])
				#This 10 is temporary until I figure out how to deal with the angle

#Function for computing all the PDB Files in a Directory
def computeAll(secondaryType, helixOrsheet):
	only_100 = 0
	for filename in glob.glob(os.path.join(path, '*.pdb')):
		if (only_100 == MAX_VALUE):
			break
		else:
			only_100+=1
			helixOrSheet(filename, secondaryType, helixOrsheet)
			print "{0} files completed...".format(only_100)

def transFunction(filename, secondaryType, helixOrsheet):
	p = buildProtein(filename)
	atomNumber = setSecondaryType(secondaryType)
	helixangles.evaluateAngles(p, filename, atomNumber)		
	#This one has to be the whole protein because we're doing TRANSITIONS
	#The type of Helix does matter in this case because the angles are calculated differently
	helixangles.transitionAngle("{0}N.txt".format(filename.split('.')[0]), filename, secondaryType, helixOrsheet)
	helixangles.transitionAngle("{0}CA.txt".format(filename.split('.')[0]), filename, secondaryType, helixOrsheet)
	helixangles.transitionAngle("{0}C.txt".format(filename.split('.')[0]), filename, secondaryType, helixOrsheet)

	#From PDB
	#---------------------HELIX-----------------------------------------
	# 1. Right Handed Alpha		6. Left-Handed Alpha
	# 2. Right Handed Omega		7. Left-Handed Omega
	# 3. Right-Handed Pi 		8. Left-Handed Gamma
	# 4. Right-Handed Gamma		9. 2/7 ribbon/Helix
	# 5. Right-Handed 3/10 		10. Polyproline *Proline breaks Helices
	#--------------------------------------------------------------------
	#
	#----------------------SHEET-----------------------------------------
	# 1. Parallel Sheet 		-1. Anti-Parallel Sheet
	# 0. The First Strand in a Sheet

if __name__ == '__main__':
	if (len(sys.argv) > 3 and (int(sys.argv[1]) in range(1,11)) and (sys.argv[3] == "helix" or sys.argv[3] == "sheet")):		# and format(sys.argv[2])
		path += sys.argv[2]
		if (os.path.isdir(path) == True):
			print "\nComputing All Files."
			computeAll(int(sys.argv[1]), sys.argv[3])
		elif format(sys.argv[2]):
			helixOrSheet(sys.argv[2], int(sys.argv[1]), sys.argv[3])
		else:
			print "\nERROR: File type is incorrect or does not exist"
	#This one only calculates the transition from the Helices to Coils
	elif (len(sys.argv) > 3 and (int(sys.argv[1]) in range(1,11)) and (sys.argv[3] == "trans")):		# and format(sys.argv[2])
		path += sys.argv[2]
		if (os.path.isdir(path) == True):
			print "\nComputing All files"
			only_100 = 0
			for filename in glob.glob(os.path.join(path, "*.pdb")):
				only_100+=1
				transFunction(filename, int(sys.argv[1]), "helix")
				print "{0} files completed...".format(only_100)
		elif format(sys.argv[2]):
			transFunction(sys.argv[2], int(sys.argv[1]), "helix")
		else:
			print "\nERROR: File type is incorrect or does not exist"
	else:
		print "\nERROR: 2 Arguments Required SecondaryType File HelixOrSheetOrTrans"