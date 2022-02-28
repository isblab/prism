import os
import pytest
import numpy as np
import sys

def test_density_calc():
	sys.path.append('./src/')
	from src.sparse_grid import SparseGrid
	from src.bead_density import BeadDensity
	from src.utils import _get_bounding_box
	arr = np.load('tests/data/input.npz')
	coords = arr['arr_0']
	grid = SparseGrid(voxel_size=4)
	grid.create_grid(coords)
	bead_density = BeadDensity(coords.shape[0], grid=grid, voxel_size=4)
	rand_id = np.random.choice(coords.shape[1], 1, replace=False)
	k1, k2 = _get_bounding_box(coords[:,rand_id[0],:])
	bead_density.construct_kernel(k1,k2)
	assert bead_density.return_density_opt(coords[:,rand_id[0],:], 3, 5, 50)[0].size != 0
	assert bead_density.return_density_opt(coords[:,rand_id[0],:], 3, 5, 50)[1].size != 0


