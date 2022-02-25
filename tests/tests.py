import os
import pytest
import numpy as np

from src.main import main_density_calc

def check_density_calc():
	from src.sparse_grid import SparseGrid
	from src.bead_density import BeadDensity
	from src.utils import _get_bounding_box
	arr = np.load('tests/data/input.npz')
	coords = arr['arr_0']
	grid = SparseGrid(voxel_size=4)
	grid.create_grid(coords)
	bead_density = BeadDensity(coords.shape[0], grid=grid, voxel_size=4)
	rand_id = np.random.choice(coords.shape[0], 1, replace=False)
	k1, k2 = _get_bounding_box(coords[:,rand_id,:])
	bead_density.construct_kernel(k1,k2)
	assert bead_density.return_density_opt(coords[:,rand_id,:], 3, 5, 50)[0].size != 0
	assert bead_density.return_density_opt(coords[:,rand_id,:], 3, 5, 50)[1].size != 0

