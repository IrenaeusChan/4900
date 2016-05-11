"""
Irenaeus Chan
5/11/2016

Helix Model Check
"""

from protein import buildProtein
import vector
import extra
import aminoacid
import atom

#Code to store turn (Large number of small segments)
def organizeHelix (helix):
	helixBackbone = []
	for AA in helix.amino_acids:
		helixBackbone += AA.backbone
	return helixBackbone
#1 Turn = 13 Atoms
#turn large list of 1 helix into smaller list of 13, run through list...

#Code to evaluate angles from each atom to it's next respective atom
def evaluateAngles (helix, filename):
	helixBackbone = organizeHelix(helix)
	totalx, totaly, totalz = 0, 0 ,0
	Nvector2, CAvector2, Cvector2 = 0, 0, 0
	for i in range(13,len(helixBackbone)+1):
		for atom in helixBackbone[i-13:i]:
		 	totalx += atom.x
		 	totaly += atom.y
		 	totalz += atom.z
		totalx /= 13
		totaly /= 13
		totalz /= 13
		center = [totalx, totaly, totalz]

		with open('{0}.txt'.format(filename), "a") as output:
			for count, atom in enumerate(helixBackbone[i-13:i]):
				if atom.atom == "N":
					output.write(str(count))
					Nvector1 = Nvector2
					Nvector2 = vector.vectorCalculation(center, [atom.x, atom.y, atom.z])
					if Nvector1 != 0:
						output.write("\nN " + str(vector.dihedralAngle(Nvector1, Nvector2)))
						#print "N " + str(vector.dihedralAngle(Nvector1, Nvector2))
				elif atom.atom == "CA":
					output.write(str(count))
					CAvector1 = CAvector2
					CAvector2 = vector.vectorCalculation(center, [atom.x, atom.y, atom.z])
					if CAvector1 != 0:
						output.write("\nCa " + str(vector.dihedralAngle(CAvector1, CAvector2)))
						#print "Ca " + str(vector.dihedralAngle(CAvector1, CAvector2))
				elif atom.atom == "C":
					output.write(str(count))
					Cvector1 = Cvector2
					Cvector2 = vector.vectorCalculation(center, [atom.x, atom.y, atom.z])
					if Cvector1 != 0:
						output.write("\nC " + str(vector.dihedralAngle(Cvector1, Cvector2)))
						#print "C " + str(vector.dihedralAngle(Cvector1, Cvector2))
		return 0
#Code to compare relative angles to each segment (Variation in curvature)
def relativeAngle ():
	return 0