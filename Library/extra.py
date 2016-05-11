"""
Irenaeus Chan
11/27/2015

Helix and Sheets
"""

from atom import Atom
from aminoacid import AminoAcid

class Helix(object):
	def __init__(self, start, stop, seqres, helixType, amino_acids):
		"""Creates a new Helix made up of Amino Acids

		Arguments:
			start: The starting position of the Helix
			stop: The ending position of the Helix
			seqres: Which chain does the Helix belong to
			amino_acids: All the Amino Acids that build up the Helix

		Exceptions:
			ValuError: If given invalid start, stop, seqres, helixType, or amino_acids
		"""

		if isinstance(start, int):
			self.start = start
		else: 
			raise ValueError('Invalid Start Position {0}'.format(start))

		if isinstance(stop, int):
			self.stop = stop
		else:
			raise ValueError('Invalid Stop Position {0}'.format(stop))

		if isinstance(seqres, basestring):
			self.seqres = seqres
		else:
			raise ValueError('Invalid SEQRES {0}'.format(seqres))

		if isinstance(helixType, int):
			self.helixType = helixType
		else:
			raise ValueError('Invalid Type {0}'.format(helixType))

		self.amino_acids = amino_acids

	def __hash__(self):
		return hash(self.__repr__())

	def __eq__(self,other):
		return self.amino_acids == other.amino_acids

	def __ne__(self, other):
		return not self.__eq__(other)

	def __repr__(self):
		helix_sequence = ""
		for AA in self.amino_acids:
			helix_sequence += "{0}\n".format(AA)
		return helix_sequence

class Sheet(object):
	def __init__(self, start, stop, seqres, sheetType, amino_acids):
		"""Creates a new Sheet made up of Amino Acids

		Arguments:
			start: The starting position of the Sheet
			stop: The ending position of the Sheet
			seqres: Which chain does the Sheet belong to
			amino_acids: All the Amino Acids that build up the Sheet

		Exceptions:
			ValuError: If given invalid start, stop, seqres, sheetType, or amino_acids
		"""
		if isinstance(start, int):
			self.start = start
		else: 
			raise ValueError('Invalid Start Position {0}'.format(start))

		if isinstance(stop, int):
			self.stop = stop
		else:
			raise ValueError('Invalid Stop Position {0}'.format(stop))

		if isinstance(seqres, basestring):
			self.seqres = seqres
		else:
			raise ValueError('Invalid SEQRES {0}'.format(seqres))

		if isinstance(sheetType, int):
			self.sheetType = sheetType
		else:
			raise ValueError('Invalid Type {0}'.format(sheetType))
			
		self.amino_acids = amino_acids

	def __hash__(self):
		return hash(self.__repr__())

	def __eq__(self,other):
		return self.amino_acids == other.amino_acids

	def __ne__(self, other):
		return not self.__eq__(other)

	def __repr__(self):
		sheet_sequence = ""
		for AA in self.amino_acids:
			sheet_sequence += "{0}\n".format(AA)
		return sheet_sequence

class Coil(object):
	def __init__(self, start, stop, amino_acids):
		"""Creates a new Coil made up of Amino Acids

		Arguments:
			start: The starting position of the Coil
			stop: The ending position of the Coil
			seqres: Which chain does the Coil belong to
			amino_acids: All the Amino Acids that build up the Coil

		Exceptions:
			ValuError: If given invalid start, stop, seqres, or amino_acids
		"""

		if isinstance(start, int):
			self.start = start
		else: 
			raise ValueError('Invalid Start Position {0}'.format(start))

		if isinstance(stop, int):
			self.stop = stop
		else:
			raise ValueError('Invalid Stop Position {0}'.format(stop))

		self.amino_acids = amino_acids

	def __hash__(self):
		return hash(self.__repr__())

	def __eq__(self,other):
		return self.amino_acids == other.amino_acids

	def __ne__(self, other):
		return not self.__eq__(other)

	def __repr__(self):
		coil_sequence = ""
		for AA in self.amino_acids:
			coil_sequence += "{0}\n".format(AA)
		return coil_sequence

def buildHelix(file_name, protein):
	sequence = []
	helixList = []

	with open(file_name, "r") as stream:
		for line in stream:
			#Reads from the file the information regarding the HELIX structures
			if (line[0:5] == "HELIX"):
				start = int(line[21:25])
				stop = int(line[33:37])
				seqres = str(line[19:20])
				helixType = int(line[39:40])

				#Looks for the position of the Amino Acid using the already parsed sequence and copies the sequence
				for AA in protein.amino_acids:
					if (AA.position >= start and AA.position <= stop):
						sequence.append(AA)
					if (AA.position == stop):
						break

				#Appends the HELIX sequence to a list of other sequences
				helixList.append(Helix(start, stop, seqres, helixType, sequence))
				sequence = []
	return helixList

def buildSheet(file_name, protein):
	sequence = []
	sheetList = []

	with open(file_name, "r") as stream:
		for line in stream:
			#Reads from the file the information regarding the SHEET structures
			if (line[0:5] == "SHEET"):
				start = int(line[22:26])
				stop = int(line[33:37])
				seqres = str(line[21:22])
				sheetType = int(line[38:40])
				#There's a problem here where the beta turns are not "coils"
				# but are a part of the sheet...

				#Looks for the position of the Amino Acid using the already parsed sequence and copies the sequence
				for AA in protein.amino_acids:
					if (AA.position >= start and AA.position <= stop):
						sequence.append(AA)
					if (AA.position == stop):
						break

				#Appends the SHEET sequence to a list of other sequences
				sheetList.append(Sheet(start, stop, seqres, sheetType, sequence))
	return sheetList

def buildCoil(file_name, protein):
	sequence = []
	coilList = []
	startList = [0]
	stopList = []

	with open(file_name, "r") as stream:
		for line in stream:
			#Reads from the file the information regarding the HELIX structures
			if (line[0:5] == "HELIX"):
				start = int(line[21:25])
				stop = int(line[33:37])
				stopList.append(start)
				startList.append(stop)
			
			#Reads from the file the information regarding the SHEET structures
			if (line[0:5] == "SHEET"):
				start = int(line[22:26])
				stop = int(line[33:37])
				stopList.append(start)
				startList.append(stop)
				#There's a problem here where the beta turns are not "coils"

		stopList.append(protein.length)

		for start, stop in zip(startList, stopList):
			for AA in protein.amino_acids:
				if (AA.position > start and AA.position < stop):
					sequence.append(AA)
			coilList.append(Coil(start, stop, sequence))
			sequence = []
	return coilList