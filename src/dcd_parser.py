from struct import *
from rmf_parser import parse_all_rmfs
import numpy as np
import os
import glob

INT_SIZE = 4

class DCDReader():
	def __init__(self, fl):
		self.fl = fl
		self.reader = open(self.fl, 'rb')

	def read_remarks(self):
		_ = self.reader.read(2*INT_SIZE)
		n_frames = unpack("i", self.reader.read(INT_SIZE))[0]
		# print("Number of models in dcd file: {}".format(n_frames))
		_ = self.reader.read(64*INT_SIZE)
		n_atoms = unpack("i", self.reader.read(INT_SIZE))[0]
		# print("Number of particles in dcd file: {}".format(n_atoms))
		_ = self.reader.read(1*INT_SIZE)
		frame_size = unpack("i", self.reader.read(INT_SIZE))[0]
		return n_frames, n_atoms, frame_size

	def return_coords(self):
		self.n_frames, self.n_atoms, self.frame_size = self.read_remarks()
		coords = []
		for frame in range(self.n_frames):
			x = list(unpack("%df" % self.n_atoms, self.reader.read(self.n_atoms*INT_SIZE)))
			_ = self.reader.read(2*INT_SIZE)
			y = list(unpack("%df" % self.n_atoms, self.reader.read(self.n_atoms*INT_SIZE)))
			_ = self.reader.read(2*INT_SIZE)
			z = list(unpack("%df" % self.n_atoms, self.reader.read(self.n_atoms*INT_SIZE)))
			_ = self.reader.read(2*INT_SIZE)
			coords.append(np.vstack((x,y,z)).T)
		return np.reshape(coords, (self.n_frames, self.n_atoms, 3))


def parse_all_dcds(folder, resolution, subunit=None, selection=None):
	try:
		dcd_file = glob.glob(os.path.join(folder, '*.dcd'))[0]
	except:
		raise ValueError('No DCD file found')
	reader = DCDReader(dcd_file)
	coords = reader.return_coords()
	try:
		rmf_coords, mass, radii, ps_names = parse_all_rmfs(folder, resolution=resolution, subunit=subunit, selection=selection)
	except:
		raise ValueError('No rmf file found')
	if rmf_coords.shape[1] != coords.shape[1]:
		raise ValueError('Number of particles in rmf file {} does not match \
			number of particles in dcd file {}'.format(rmf_coords.shape[1], coords.shape[1]))
	return coords, mass, radii, ps_names

