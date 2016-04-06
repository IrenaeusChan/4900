"""
Irenaeus Chan
3/23/2016

New Coordinate Generation System
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

def calculateCoordinates(protein, filename):
	write = 'w'
	if (sys.argv[1] == "all" and len(sys.argv) > 2):
		write = 'a'

	with open('{0}.txt'.format(filename), write) as output:
		#Sets an iterator to examine the previous, current, and next values
		for prev, AA, nxt in previousAndNext(protein.amino_acids):
			#We can't calculate the beginning of residues because we need the "Wing" to exist
			# The first argument ensures it's not the first in the sequence, the second ensures it's not the beginning of the residue
			if (prev is None or prev.seqres != AA.seqres):
				#[0] = N, [1] = Ca, [2] = C
				center = [AA.backbone[1].x, AA.backbone[1].y, AA.backbone[1].z]
				front = [AA.backbone[2].x, AA.backbone[2].y, AA.backbone[2].z]
				wing = [AA.backbone[0].x, AA.backbone[0].y, AA.backbone[0].z]

				#We want to develop the initial coordinate system using the first AminoAcid
				vectorX = vector.vectorCalculation(center, front)
				vectorZ = vector.vectorCalculation(center, wing)
				#When Z is positive, the cross product has to be done differently where Z then X
				# if Z is negative, the cross product has to be done with X then Z
				if vectorZ[2] < 0:
					vectorY = vector.crossProduct(vectorX, vectorZ)
				else:
					vectorY = vector.crossProduct(vectorZ, vectorX)
				
				#Vector Z may not be perpindicular to both Vector X and Vector Y
				vectorZ = vector.crossProduct(vectorX, vectorY)
				#This should guaruntee that vector Z is now perpindicular to BOTH X and Y and not just Y
				continue
			#This checks if it's the end of the ENTIRE sequence and whether or not it is at the end of the residue
			elif (nxt is None or AA.seqres !=  nxt.seqres):
				continue
			else:
				for i in range(0,3):
					#This is because, we calculate PhiPsi as NCa, but this Vector is calculated by doing CaN
					normalVector1 = vector.crossProduct(vector.vectorCalculation(wing,center),vectorX)

					C = front
					N = [AA.backbone[i].x, AA.backbone[i].y, AA.backbone[i].z]
					#if i == 2:
					#	Ca = [nxt.backbone[0].x, nxt.backbone[0].y, nxt.backbone[0].z]
					#else:
					#	Ca = [AA.backbone[i+1].x, AA.backbone[i+1].y, AA.backbone[i+1].z]
					vectorCN = vector.vectorCalculation(C, N)
					#vectorNCa = vector.vectorCalculation(N, Ca)
					#^This is supposed to be vectorX
					normalVector2 = vector.crossProduct(vectorCN, vectorX)

					#Same algorithm as calculating PhiPsi Angles
					angle = vector.dihedralAngle(normalVector1, normalVector2)
					if vector.dotProduct(vector.crossProduct(normalVector2, normalVector1), vectorX) < 0:
						angle = -angle

					#Just to make it easier for myself to remember which is what
					origin = center
					wing = center
					center = front
					front = [AA.backbone[i].x, AA.backbone[i].y, AA.backbone[i].z]

					#I want center to move to origin and I want front to move to new vector relative to origin
					#new = [front[0] - (center[0] - origin[0]), front[1] - (center[1] - origin[1]), front[2] - (center[2] - origin[2])]

					#The "scalar projection" of the VectorPoint onto Vector(X, Y, or Z) IS the LENGTH of (x, y, or z)
					vectorPoint = vector.vectorCalculation(origin, center)
					x1 = (vector.dotProduct(vectorX, vectorPoint))/vector.vectorMagnitude(vectorX)
					y1 = (vector.dotProduct(vectorY, vectorPoint))/vector.vectorMagnitude(vectorY)
					z1 = (vector.dotProduct(vectorZ, vectorPoint))/vector.vectorMagnitude(vectorZ)

					vectorPoint = vector.vectorCalculation(origin, front)
					x2 = (vector.dotProduct(vectorX, vectorPoint))/vector.vectorMagnitude(vectorX)
					y2 = (vector.dotProduct(vectorY, vectorPoint))/vector.vectorMagnitude(vectorY)
					z2 = (vector.dotProduct(vectorZ, vectorPoint))/vector.vectorMagnitude(vectorZ)

					newCoordVector = vector.vectorCalculation([x1,y1,z1], [x2,y2,z2])

					output.write(str(AA.backbone[i].atom) + ' ' + str(newCoordVector[0]) + ' ' + str(newCoordVector[1]) + ' ' + str(newCoordVector[2]) + ' ' + str(angle) + '\n')

					#Now we change the coordinate system to the new line
					vectorX = vector.vectorCalculation(center, front)
					vectorZ = vector.vectorCalculation(center, wing)
					#When Z is positive, the cross product has to be done differently where Z then X
					# if Z is negative, the cross product has to be done with X then Z
					if vectorZ[2] < 0:
						vectorY = vector.crossProduct(vectorX, vectorZ)
					else:
						vectorY = vector.crossProduct(vectorZ, vectorX)
					
					#Vector Z may not be perpindicular to both Vector X and Vector Y
					vectorZ = vector.crossProduct(vectorX, vectorY)
					#This should guaruntee that vector Z is now perpindicular to BOTH X and Y and not just Y