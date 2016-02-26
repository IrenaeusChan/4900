"""
Irenaeus Chan
1/27/2015

Amino Acid Phi Psi Store Class
"""

AMINO_ACIDS = {'GLY', 'ALA', 'SER', 'THR', 'CYS', 'VAL', 'LEU', 'ILE', 'MET', 'PRO', 'PHE', 'TYR', 'TRP', 'ASP', 'GLU', 'ASN', 'GLN', 'HIS', 'LYS', 'ARG'}
AA_SINGLE = {'GLY':'G', 'ALA':'A', 'SER':'S', 'THR':'T', 'CYS':'C', 'VAL':'V', 'LEU':'L', 'ILE':'I', 'MET':'M', 'PRO':'P', 'PHE':'F', 'TYR':'Y', 'TRP':'W', 'ASP':'D', 'GLU':'E', 'ASN':'N', 'GLN':'Q', 'HIS':'H', 'LYS':'K', 'ARG':'R'}

class PhiPsiStore(object):
	"""A configuration for a single Amino Acid and its Phi Psi Angles"""

	def __init__(self, pos, amino_acid, phi, psi, distance):
		"""Creates a new Amino Acid

		Arguments:
			position: The residual position
			amino_acid: The specific Amino Acid
			phi: The calculated phi angle
			psi: The calculated psi angle
			distance: The distance the Amino Acid is from the Center

		Exceptions:
			ValuError: If given invalid amino_acid, phi, psi, or distance
		"""

		if isinstance (pos, int):
			self.pos = pos
		else:
			raise ValueError('Invalid Residual Position {0}'.format(pos)) 

		if amino_acid in AMINO_ACIDS:
			self.amino_acid = amino_acid
		else:
			raise ValueError('Invalid Amino Acid {0}'.format(amino_acid))

		if isinstance(phi, float):
			self.phi = phi
		else:
			raise ValueError('Invalid Position {0}'.format(phi))

		if isinstance(psi, float):
			self.psi = psi
		else:
			raise ValueError('Invalid Position {0}'.format(psi))

		if isinstance(distance, float):
			self.distance = distance
		else:
			raise ValueError('Invalid Position {0}'.format(distance))

	def __hash__(self):
		return hash(self.__repr__())

	def __eq__(self, other):
		return self.__dict__ == other.__dict__

	def __ne__(self, other):
		return not self.__eq__(other)

	def __repr__(self):
		return "{0} {1} {2} {3} {4}".format(self.pos, self.amino_acid, self.phi, self.psi, self.distance)

	def getLetter(self):
		return AA_SINGLE[self.amino_acid]