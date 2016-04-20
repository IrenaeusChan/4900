"""
Irenaeus Chan
3/23/2016

New Coordinate Generation System
"""

import sys
from itertools import tee, islice, chain, izip
import vector

SYMBOL = {'GLY':'G', 'ALA':'A', 'SER':'S', 'THR':'T', 'CYS':'C', 'VAL':'V', 'LEU':'L', 'ILE':'I', 'MET':'M', 'PRO':'P', 'PHE':'F', 'TYR':'Y', 'TRP':'W', 'ASP':'D', 'GLU':'E', 'ASN':'N', 'GLN':'Q', 'HIS':'H', 'LYS':'K', 'ARG':'R'}

def formatSeq(sequence):
	seq = ""
	for AA in sequence:
		seq += SYMBOL[AA]
	return seq

def previousAndNext(some_iterable):
#http://stackoverflow.com/questions/1011938/python-previous-and-next-values-inside-a-loop
	prevs, items, nexts = tee(some_iterable, 3)
	prevs = chain([None], prevs)
	nexts = chain(islice(nexts, 1, None), [None])
	return izip(prevs, items, nexts)

def calculateCoordinates(protein, filename):
	write = 'w'

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
			else:
				for i in range(0,3):
					#This is because, we calculate PhiPsi as NCa, but this Vector is calculated by doing CaN
					normalVector1 = vector.crossProduct(vector.vectorCalculation(wing,center),vectorX)
					C = front
					N = [AA.backbone[i].x, AA.backbone[i].y, AA.backbone[i].z]
					vectorCN = vector.vectorCalculation(C, N)
					normalVector2 = vector.crossProduct(vectorX, vectorCN)

					#Same algorithm as calculating PhiPsi Angles
					angle1 = vector.dihedralAngle(normalVector1, normalVector2)
					if vector.dotProduct(vector.crossProduct(normalVector1, normalVector2), vectorX) < 0:
						angle1 = -angle1
					"""
					if nxt is None:
						angle2 = 0
					else:
					"""
					if i == 2:
						if nxt is not None:
							Ca = [nxt.backbone[0].x, nxt.backbone[0].y, nxt.backbone[0].z]
						else:
							Ca = None
					else:
						Ca = [AA.backbone[i+1].x, AA.backbone[i+1].y, AA.backbone[i+1].z]
					if Ca != None:
						vectorNCa = vector.vectorCalculation(N, Ca)
						normalVector3 = vector.crossProduct(vectorCN, vectorNCa)

						#Same algorithm as calculating PhiPsi Angles
						angle2 = vector.dihedralAngle(normalVector1, normalVector3)
						if vector.dotProduct(vector.crossProduct(normalVector1, normalVector3), vectorX) < 0:
							angle2 = -angle2
					else:
						angle2 = "NULL"

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

					output.write(str(AA.position) + ' ' + str(AA.backbone[i].atom) + ' ' + str(newCoordVector[0]) + ' ' + str(newCoordVector[1]) + ' ' + str(newCoordVector[2]) + ' ' + str(angle1) + ' ' + str(angle2) + '\n')

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

def calculate(protein, filename):
	write = 'a'
	with open('NEW_{0}.txt'.format(filename), write) as output:
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
				#This is because, we calculate PhiPsi as NCa, but this Vector is calculated by doing CaN
				normalVector1 = vector.crossProduct(vector.vectorCalculation(wing,center),vectorX)
				origin = center
				seq = [AA.amino_acid]
				continue
			else:
				seq.append(AA.amino_acid)
				for i in range(0,3):
					C = front
					N = [AA.backbone[i].x, AA.backbone[i].y, AA.backbone[i].z]
					vectorCN = vector.vectorCalculation(C, N)
					normalVector2 = vector.crossProduct(vectorX, vectorCN)

					#Same algorithm as calculating PhiPsi Angles
					angle1 = vector.dihedralAngle(normalVector1, normalVector2)
					if vector.dotProduct(vector.crossProduct(normalVector1, normalVector2), vectorX) < 0:
						angle1 = -angle1

					if i == 2:
						if nxt is not None:
							Ca = [nxt.backbone[0].x, nxt.backbone[0].y, nxt.backbone[0].z]
						else:
							Ca = None
					else:
						Ca = [AA.backbone[i+1].x, AA.backbone[i+1].y, AA.backbone[i+1].z]
					
					if Ca != None:
						vectorNCa = vector.vectorCalculation(N, Ca)
						normalVector3 = vector.crossProduct(vectorCN, vectorNCa)

						#Same algorithm as calculating PhiPsi Angles
						angle2 = vector.dihedralAngle(normalVector1, normalVector3)
						if vector.dotProduct(vector.crossProduct(normalVector1, normalVector3), vectorX) < 0:
							angle2 = -angle2
					else:
						angle2 = "NULL"

					#Just to make it easier for myself to remember which is what
					wing = center
					center = front
					front = [AA.backbone[i].x, AA.backbone[i].y, AA.backbone[i].z]

					#I want center to move to origin and I want front to move to new vector relative to origin
					#new = [front[0] - (center[0] - origin[0]), front[1] - (center[1] - origin[1]), front[2] - (center[2] - origin[2])]

					#The "scalar projection" of the VectorPoint onto Vector(X, Y, or Z) IS the LENGTH of (x, y, or z)
					vectorPoint = vector.vectorCalculation(origin, center)
					magnitude = vector.vectorMagnitude(vectorPoint)
					x1 = (vector.dotProduct(vectorX, vectorPoint))/vector.vectorMagnitude(vectorX)
					y1 = (vector.dotProduct(vectorY, vectorPoint))/vector.vectorMagnitude(vectorY)
					z1 = (vector.dotProduct(vectorZ, vectorPoint))/vector.vectorMagnitude(vectorZ)

					vectorPoint = vector.vectorCalculation(origin, front)
					x2 = (vector.dotProduct(vectorX, vectorPoint))/vector.vectorMagnitude(vectorX)
					y2 = (vector.dotProduct(vectorY, vectorPoint))/vector.vectorMagnitude(vectorY)
					z2 = (vector.dotProduct(vectorZ, vectorPoint))/vector.vectorMagnitude(vectorZ)

					newCoordVector = vector.vectorCalculation([x1,y1,z1], [x2,y2,z2])

					if nxt is None:
						sequence = formatSeq(seq)
						output.write(str(AA.position) + ' ' + str(AA.backbone[i].atom) + ' ' + str(magnitude) + ' ' + 
							str(newCoordVector[0]) + ' ' + str(newCoordVector[1]) + ' ' + str(newCoordVector[2]) + ' ' + 
							str(angle1) + ' ' + str(angle2) + ' ' + str(sequence) + '\n')