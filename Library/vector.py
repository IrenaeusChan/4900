"""
Irenaeus Chan
11/27/2015

Math Calculations Required for Protein Analysis
"""

import sys
import math

def vectorCalculation(coord1, coord2):
#Converts the coordinates for the atoms into a vector between the two atoms representing the "bond"
	return [coord1[0]-coord2[0], coord1[1]-coord2[1], coord1[2]-coord2[2]]

def crossProduct(vector1, vector2):
	x = vector1[1]*vector2[2] - vector1[2]*vector2[1]
	y = vector1[2]*vector2[0] - vector1[0]*vector2[2]
	z = vector1[0]*vector2[1] - vector1[1]*vector2[0]
	return [x,y,z]

def dotProduct(vector1, vector2):
	return vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2]

def vectorMagnitude(vector):
	return ((vector[0])**2 + (vector[1])**2 + (vector[2])**2)**0.5

def dihedralAngle(normalVector1, normalVector2):
	return math.degrees(math.acos(dotProduct(normalVector1, normalVector2)/(vectorMagnitude(normalVector1) * vectorMagnitude(normalVector2))))