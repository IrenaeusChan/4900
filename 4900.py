import sys
import string
import glob
import os
sys.path.append(os.path.realpath("Library"))
import phipsi
import protein
import smithwaterman
import newCoord

path = os.getcwd()

def calculate(file_name):
	p = protein.buildProtein(file_name)
	cenx, ceny, cenz = protein.weightedAverage(p)
	center = [cenx, ceny, cenz]
	phipsi.calculatePhiPsi(p, center, file_name)
	return p, center

if __name__ == '__main__':
	only_100 = 0
	if len(sys.argv) < 2:
		print "\nERROR: No file was provided"
		sys.exit(1)
	elif (sys.argv[1].endswith(".pdb")):
		print "\nComputing Using File: {0}".format(sys.argv[1])
		p, center = calculate(sys.argv[1])
	elif (sys.argv[1] == "all" and len(sys.argv) > 2):
		print "\nComputing All Files."
		path += sys.argv[2]
		for filename in glob.glob(os.path.join(path, '*.pdb')):
			if (only_100 == 100):
				break
			else:
				print "{0} files completed...".format(only_100)
				p, center = calculate(filename)
				protein.relativeToCenter(p, center)
				only_100+=1
	elif (sys.argv[1] == "new" and len(sys.argv) > 2): 	#New Geometric Calculations
														#Working Here
		p = protein.buildProtein(sys.argv[2])				
		newCoord.calculateCoordinates(p, sys.argv[2])
	elif (sys.argv[1] == "align" and len(sys.argv) > 3):
		print "\nAligning your two files."
		if os.path.isfile(sys.argv[2]) == True:
			mainSeq = smithwaterman.readSeq(sys.argv[2])
		else:
			p, center = calculate(sys.argv[2][:-4])
			mainSeq = smithwaterman.readSeq(sys.argv[2])
		if os.path.isfile(sys.argv[3]) == True:
			targetSeq = smithwaterman.readSeq(sys.argv[3])
		else:
			p, center = calculate(sys.argv[3][:-4])
			targetSeq = smithwaterman.readSeq(sys.argv[3])
		smithwaterman.align(mainSeq, targetSeq, sys.argv[2], sys.argv[3])
	else:
		print "\nERROR: File type is incorrect."
		print "\nERROR: If ALL selected, set path directory."
		sys.exit(1)
		
	"""
	#If the user wants to calculate the PhiPsi angles for ONLY the Helix or SHEETS	
	#Temperarily unavailable
	if (len(sys.argv) > 2 and (sys.argv[2] == "Helix" or sys.argv[2] == "helix")):
		#We need to read the PDB file once again to get more information about positions
		helixList, sheetList = readFileExtra(sys.argv[1], protein)
		if (len(helixList)<1):
			print "\nThere were no HELICES in this protein"
		for helix in helixList:
			calculatePhiPsi(helix, center)
	#Same for HELIX
	#Temperarily unavailable
	elif (len(sys.argv) > 2 and (sys.argv[2] == "Sheets" or sys.argv[2] == "sheets")):
		helixList, sheetList = readFileExtra(sys.argv[1], protein)
		if (len(sheetList)<1):
			print "\nThere were no SHEETS in this protein"
		for sheet in sheetList:
			calculatePhiPsi(sheet, center)
	"""

	if (len(sys.argv) > 2 and (sys.argv[2] == "Center" or sys.argv[2] == "center")):
		print "Protein's Center (X, Y, Z): " + str(protein.weightedAverage(p))
	elif (len(sys.argv) > 2 and (sys.argv[2] == "Distance" or sys.argv[2] == "distance")):
		protein.relativeToCenter(p, center)
	else:
		print "Other Options Available:"
		print "\tCenter"
		print "\tDistance"
		print ""
		print "If you would like to calculate entire files of PDB Files:"
		print "\tall file_path"
		print "To align two PDB files that have existing PhiPsi Calculations already done:"
		print "\talign file_name1 file_name2"

