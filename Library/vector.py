"""
Irenaeus Chan
11/27/2015

Math Calculations Required for Protein Analysis
"""

import numpy as np
from decimal import *

def vectorCalculation(coord1, coord2):
#Converts the coordinates for the atoms into a vector between the two atoms representing the "bond"
	return [coord2[0]-coord1[0], coord2[1]-coord1[1], coord2[2]-coord1[2]]

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

def orthogonalDistanceRegression(listOfCoord):
#Original code by @dwf from http://stackoverflow.com/questions/2298390/fitting-a-line-in-3d
#Edited and changed accordingly by Irenaeus Chan
	data = np.array(listOfCoord)			#As long as the list of coordinates is [[x, y, z], [x, y, z], [x, y, z]]
	mean = data.mean(axis=0)				#Calculate the mean of the points, i.e. the "center" of the cloud
	uu, dd, vv = np.linalg.svd(data-mean)	#Do an SVD on the mean-centered data

	#According to @dwf
	#vv[0] contains the first principal component, i.e. the direction vector of the "best fit" line
	#Therefore in order to generate two points on the line, we will take the vector and extrapolate
	# a set amount of distance ahead and behind the vector (vector projection)
	# ---------------------------------------------------------------------------------------------
	#This means np.mgrid creates a matrix that contains the distance away from the central cloud
	# we reorganize that matrix into a single vertical array using the np.newaxis function on the
	# whole array
	pointsOnLine = vv[0] * np.mgrid[-50:50:2j][:, np.newaxis]	
	#The -10:10 accounts for a rough spread of how far the start and stop atoms are from each 
	# other, the distance was eyeballed using R statistical software and is arbitrary
	pointsOnLine += mean
	pointsOnLine = pointsOnLine.tolist()
	regressionVector = vectorCalculation(pointsOnLine[0], pointsOnLine[1])
	#regressionVector = vectorCalculation([decimal.Decimal(pointsOnLine[0][0]), decimal.Decimal(pointsOnLine[0][1]), decimal.Decimal(pointsOnLine[0][2])], 
	#	[decimal.Decimal(pointsOnLine[1][0]), decimal.Decimal(pointsOnLine[1][1]), decimal.Decimal(pointsOnLine[1][2])])
	regressionVector = (round(Decimal(regressionVector[0]), 3), round(Decimal(regressionVector[1]), 3), round(Decimal(regressionVector[2]), 3))
	pointsOnLine[0] = [round(Decimal(pointsOnLine[0][0]), 3), round(Decimal(pointsOnLine[0][1]), 3), round(Decimal(pointsOnLine[0][2]), 3)]
	return regressionVector, pointsOnLine[0]

def orthogonalLineCalculation(vector, point):
	#Vector V defined as (x, y, z) and point P defined as (a, b, c).
	#Any point along vector V can be represented by (xk, yk, zk)
	#Let us assume that (xk-a, xk-b, xk-c) is a NEW vector, L from the point, P to the vector, V
	# that is also orthogonal to vector V
	#To find the values of the orthogonal vector, L, the dot product between the orthogonal
	# vector, L and vector V must be equivalent to 0
	#Therefore, (xk-a, yk-b, zk-c) . (x, y, z) = 0 would yeild the appropriate orthogonal vector, L
	#As such, k = (ax + by + cz)/(x^2 + y^2 + z^2)
	point = [round(Decimal(point[0]),3), round(Decimal(point[1]),3), round(Decimal(point[2]),3)]
	k = ((point[0]*vector[0]) + (point[1]*vector[1]) + (point[2]*vector[2]))/((vector[0]**2) + (vector[1]**2) + (vector[2]**2))
	orthogonalLine = (point[0]-vector[0]*k, point[1]-vector[1]*k, point[2]-vector[2]*k)
	orthogonalLine = (round(Decimal(orthogonalLine[0]), 3), round(Decimal(orthogonalLine[1]), 3), round(Decimal(orthogonalLine[2]), 3))
	#Since (xk-a, xk-b, xk-c) is a vector from point, P to the vector, V, we want to look for
	# the vector that goes from the point on vector, V to point, P, meaning (a-xk, b-xk, c-xk)
	return orthogonalLine

def orthogonalVectorCalculation(Vorth, Patom, Vregr, Pregr):
	print "HERE"
	print Vorth
	print Patom
	print Vregr
	print Pregr

	l = ((Vregr[1]*Patom[0]) - (Vregr[1]*Pregr[0]) - (Vregr[0]*Patom[1]) + (Vregr[0]*Pregr[1]))/((Vregr[0]*Vorth[1]) - (Vregr[1]*Vorth[0]))
	m = (Patom[0] + (l*Vorth[0]) - Pregr[0])/(Vregr[0])
	print "ANSWERS"
	print l
	print m

	print "COORDINATES"
	print (Patom[2]+(l*Vorth[2]))
	print (Pregr[2]+(m*Vregr[2]))

	if ((Patom[2]+(l*Vorth[2])) == (Pregr[2]+(m*Vregr[2]))):
		pointOfIntersection = (Patom[0]+(l*Vorth[0]), Patom[1]+(l*Vorth[1]), Patom[2]+(l*Vorth[2]))
		print "SUCCESS"
		print pointOfIntersection
		return vectorCalculation(pointOfIntersection, Patom)
	else:
		print "There was no intersection..."
	print "END"

#R Code
#data<-read.csv(choose.files(), header=T)
#s<-scatterplot3d(data$x, data$y, data$z)
#p1<-s$xyz.convert(data$xpoints[1], data$ypoints[1], data$zpoints[1])
#p2<-s$xyz.convert(data$xpoints[2], data$ypoints[2], data$zpoints[2])
#segments(p1$x, p1$y, p2$x, p2$y, lwd=2, col=2)