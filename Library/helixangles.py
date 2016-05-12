"""
Irenaeus Chan
5/11/2016

Helix Model Check
"""

import sys
from protein import buildProtein
import vector
import extra
import aminoacid
import atom

#Code to store turn (Large number of small segments)
def organizeHelix (helix):
	helixBackbone = []
	position = []
	for AA in helix.amino_acids:
		helixBackbone += AA.backbone
		position += [AA.position, AA.position, AA.position]
	return helixBackbone, position
#1 Turn = 13 Atoms
#turn large list of 1 helix into smaller list of 13, run through list...

#Code to evaluate angles from each atom to it's next respective atom
def evaluateAngles (helix, filename):
	helixBackbone, position = organizeHelix(helix)
	for i in range(13,len(helixBackbone)+1):
		#We want to initialize these variables for EACH 13 Atom Section
		totalx, totaly, totalz = 0, 0 ,0
		Nvector2, CAvector2, Cvector2 = 0, 0, 0
		secondNpos, secondCApos, secondCpos = 0, 0, 0	

		#Calculates the Average Center of the 13 Atoms	
		for atom in helixBackbone[i-13:i]:
		 	totalx += atom.x
		 	totaly += atom.y
		 	totalz += atom.z
		totalx /= 13
		totaly /= 13
		totalz /= 13
		center = [totalx, totaly, totalz]

		
		#For each atom in the 13 Atom Turn, we want to calculate the Angles
		for count, atom in enumerate(helixBackbone[i-13:i]):
			if (atom.atom == "N"):
				#output.write(str(count))
				#The first vector should be the last N that was found, compared to...
				Nvector1 = Nvector2
				firstNpos = secondNpos
				# the next N that will be in the sequence
				Nvector2 = vector.vectorCalculation(center, [atom.x, atom.y, atom.z])
				secondNpos = position[i-13:i][count]
				#As long as the we have at least found TWO Nitrogens, we can calculate their angles
				if Nvector1 != 0:
					with open('N.txt', "a") as output:
						output.write("N " + str(vector.dihedralAngle(Nvector1, Nvector2)) + " " + 
							str(firstNpos) + " " + str(secondNpos) + "\n")
					#print "N " + str(vector.dihedralAngle(Nvector1, Nvector2))
			elif (atom.atom == "CA"):
				#output.write(str(count))
				#The first vector should be the last Ca that was found, compared to...
				CAvector1 = CAvector2
				firstCApos = secondCApos
				# the next N that will be in the sequence
				CAvector2 = vector.vectorCalculation(center, [atom.x, atom.y, atom.z])
				secondCApos = position[i-13:i][count]
				#As long as the we have at least found TWO alpha-Carbons, we can calculate their angles
				if CAvector1 != 0:
					with open('CA.txt', "a") as output:
						output.write("Ca " + str(vector.dihedralAngle(CAvector1, CAvector2)) + " " + 
							str(firstCApos) + " " + str(secondCApos) + "\n")
					#print "Ca " + str(vector.dihedralAngle(CAvector1, CAvector2))
			elif (atom.atom == "C"):
				#output.write(str(count))
				#The first vector should be the last C that was found, compared to...
				Cvector1 = Cvector2
				firstCpos = secondCpos
				# the next N that will be in the sequence
				Cvector2 = vector.vectorCalculation(center, [atom.x, atom.y, atom.z])
				secondCpos = position[i-13:i][count]
				#As long as the we have at least found TWO Carbons, we can calculate their angles
				if Cvector1 != 0:
					with open('C.txt', "a") as output:
						output.write("C " + str(vector.dihedralAngle(Cvector1, Cvector2)) + " " + 
							str(firstCpos) + " " + str(secondCpos) + "\n")
					#print "C " + str(vector.dihedralAngle(Cvector1, Cvector2))

#Code to compare relative angles to each segment (Variation in curvature)
def relativeAngle():
	return 0