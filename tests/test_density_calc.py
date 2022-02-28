import os
import pytest
import numpy as np
import pytest
import sys
import pandas as pd
sys.path.append('./src/')


@pytest.fixture(scope='session')
def output_dir(tmp_path_factory):
    return tmp_path_factory.mktemp("output_dir")


def test_density_calc():
	from src.sparse_grid import SparseGrid
	from src.bead_density import BeadDensity
	from src.utils import _get_bounding_box
	from src.patch_computer import calc_bead_spread
	arr = np.load('tests/data/prism_input.npz')
	coords = arr['arr_0']
	grid = SparseGrid(voxel_size=4)
	grid.create_grid(coords)
	bead_density = BeadDensity(coords.shape[0], grid=grid, voxel_size=4)
	rand_id = np.random.choice(coords.shape[1], 1, replace=False)
	k1, k2 = _get_bounding_box(coords[:,rand_id[0],:])
	bead_density.construct_kernel(k1,k2)
	densities = bead_density.return_density_opt(coords[:,rand_id[0],:], 3, 5, 50)
	assert densities[0].shape != 0
	bead_spread = calc_bead_spread(densities, grid)
	assert bead_spread >= 0


def test_patch_calc(output_dir):
	from patch_computer import get_patches, annotate_patches
	with open('tests/data/bead_spreads.txt') as f:
		bs = f.read().splitlines()
	bs = [float(i) for i in bs]
	patches = get_patches(bs, 2, np.random.randn(100, len(bs), 3), [3]*len(bs))
	annotated_patches = annotate_patches(patches, 2, ['test']*len(bs) , len(bs))
	high_prec, low_prec = patches[:2], patches[2+1:]	
	annot_df = pd.DataFrame(np.array(annotated_patches), columns = ['Bead', 'Bead Name', 'Type', 'Class', 'Patch'])
	annot_df['Patch'] = pd.to_numeric(annot_df["Patch"])
	annot_df.sort_values(['Patch'], ascending=[True])
	annot_df.to_csv(str(output_dir) + '/annotations.txt', index=None)
	assert os.path.exists(os.path.join( str(output_dir), 'annotations.txt'))









