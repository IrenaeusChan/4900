import sys
import string
import glob
import os
sys.path.append(os.path.realpath("Library"))
import phipsi
import protein
import smithwaterman
import newCoord
import extra

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
	#If the user wants to calculate the PhiPsi angles for ONLY the Helix or SHEETS	
	elif (len(sys.argv) > 2 and (sys.argv[1] == "Helix" or sys.argv[1] == "helix")):
		#We need to read the PDB file once again to get more information about positions
		p = protein.buildProtein(sys.argv[2])
		cenx, ceny, cenz = protein.weightedAverage(p)
		center = [cenx, ceny, cenz]
		helixList, sheetList, coilList = extra.readFileExtra(sys.argv[2], p)
		if (len(helixList)<1):
			print "\nThere were no HELICES in this protein"
		for helix in helixList:
			phipsi.calculatePhiPsi(helix, center, sys.argv[2])
	#Same for HELIX
	elif (len(sys.argv) > 2 and (sys.argv[1] == "Sheets" or sys.argv[1] == "sheets")):
		p = protein.buildProtein(sys.argv[2])
		cenx, ceny, cenz = protein.weightedAverage(p)
		center = [cenx, ceny, cenz]
		helixList, sheetList, coilList = extra.readFileExtra(sys.argv[2], p)
		if (len(sheetList)<1):
			print "\nThere were no SHEETS in this protein"
		for sheet in sheetList:
			phipsi.calculatePhiPsi(sheet, center, sys.argv[2])
	elif (sys.argv[1] == "coil" and len(sys.argv) > 2):
		p = protein.buildProtein(sys.argv[2])
		helixList, sheetList, coilList = extra.readFileExtra(sys.argv[2], p)
		#Run with CoilList. Does not generate BEFORE and AFTER yet. But should.
		if (len(coilList)<1):
			print "\nThere were no COILS in this protein"
		for coil in coilList:
			newCoord.calculate(coil, sys.argv[2])
	else:
		print "\nERROR: File type is incorrect."
		print "\nERROR: If ALL selected, set path directory."
		sys.exit(1)

	if (len(sys.argv) > 2 and (sys.argv[2] == "Center" or sys.argv[2] == "center")):
		print "Protein's Center (X, Y, Z): " + str(protein.weightedAverage(p))
	elif (len(sys.argv) > 2 and (sys.argv[2] == "Distance" or sys.argv[2] == "distance")):
		protein.relativeToCenter(p, center)
	else:
		print "Other Options Available:"
		print "\tCenter"
		print "\tDistance"
		print "\tNew Geometric"
		print "\tHelix"
		print "\tSheet"
		print ""
		print "If you would like to calculate entire files of PDB Files:"
		print "\tall file_path"
		print "To align two PDB files that have existing PhiPsi Calculations already done:"
		print "\talign file_name1 file_name2"

