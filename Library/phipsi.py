"""
Irenaeus Chan
11/27/2015

Ramachandran Plot Generator
"""

import sys
from itertools import tee, islice, chain, izip
import vector

def previousAndNext(some_iterable):
#http://stackoverflow.com/questions/1011938/python-previous-and-next-values-inside-a-loop
	prevs, items, nexts = tee(some_iterable, 3)
	prevs = chain([None], prevs)
	nexts = chain(islice(nexts, 1, None), [None])
	return izip(prevs, items, nexts)

def calculatePhiPsi(protein, center, filename):
	write = 'w'
	if (len(sys.argv) > 2 and (sys.argv[1] == "Helix" or sys.argv[1] == "helix")):
		write = 'a'
	elif (len(sys.argv) > 2 and (sys.argv[1] == "Sheets" or sys.argv[1] == "sheets")):
		write = 'a'
	elif (sys.argv[1] == "all" and len(sys.argv) > 2):
		write = 'a'

	with open('{0}.txt'.format(filename), write) as output:
		#Sets an iterator to examine the previous, current, and next values
		for prev, AA, nxt in previousAndNext(protein.amino_acids):
			#Due to how PhiPsi angles are calculated we can't calculate the beginning and end of residues
			# The first argument ensures it's not the first in the sequence, the second ensures it's not the beginning of the residue
			if (prev is None or prev.seqres != AA.seqres):
				C = [AA.backbone[2].x, AA.backbone[2].y, AA.backbone[2].z]
				continue
			#This checks if it's the end of the ENTIRE sequence
			elif (nxt is None):
				continue
			#This checks whether or not it is at the end of the residue
			elif (AA.seqres !=  nxt.seqres):
				continue

			aa = AA.amino_acid
			pos = AA.position
			
			N = [AA.backbone[0].x, AA.backbone[0].y, AA.backbone[0].z]
			Ca = [AA.backbone[1].x, AA.backbone[1].y, AA.backbone[1].z]
			vectorCN = vector.vectorCalculation(C, N)
			vectorNCa = vector.vectorCalculation(N, Ca)
			normalVector1 = vector.crossProduct(vectorCN, vectorNCa)

			C = [AA.backbone[2].x, AA.backbone[2].y, AA.backbone[2].z]
			vectorCaC = vector.vectorCalculation(Ca, C)
			normalVector2 = vector.crossProduct(vectorNCa, vectorCaC)

			phi = vector.dihedralAngle(normalVector1, normalVector2)
			#The cross product vectors are both normal to the axis vectorNCa (central vector),
			# so the angle between them is the dihedral angle that we are looking for.  
			# However, since "angle" only returns values between 0 and pi, we need to make
			# sure we get the right sign relative to the rotation axis
			if vector.dotProduct(vector.crossProduct(normalVector1, normalVector2), vectorNCa) < 0:
				phi = -phi

			normalVector1 = vector.crossProduct(vectorNCa,vectorCaC)
			N = [nxt.backbone[0].x, nxt.backbone[0].y, nxt.backbone[0].z]
			vectorCN = vector.vectorCalculation(C, N)
			normalVector2 = vector.crossProduct(vectorCaC, vectorCN)

			psi = vector.dihedralAngle(normalVector1, normalVector2)
			#The cross product vectors are both normal to the axis vectorNCa (central vector),
			# so the angle between them is the dihedral angle that we are looking for.  
			# However, since "angle" only returns values between 0 and pi, we need to make
			# sure we get the right sign relative to the rotation axis
			if vector.dotProduct(vector.crossProduct(normalVector1, normalVector2), vectorCaC) < 0:
				psi = -psi

			aminoacid = [AA.avgx, AA.avgy, AA.avgz]
			d = vector.vectorMagnitude(vector.vectorCalculation(center, aminoacid))
			#Writes the Phi, Psi, and Distances for the specific Amino Acid
			output.write(str(pos) + ' ' + aa + ' ' + str(phi) + ' ' + str(psi) + ' ' + str(d) + '\n')