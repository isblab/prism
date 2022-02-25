import numpy as np
import pandas as pd
import os
from multiprocessing import Pool
from functools import partial
from sparse_grid import SparseGrid
from bead_density import BeadDensity
from patch_computer import calc_bead_spread, get_patches, annotate_patches
from utils import _get_bounding_box
import argparse


def main_density_calc(i, grid, voxel_size, n_breaks):
  bead_density = BeadDensity(coords.shape[0], grid=grid, voxel_size=voxel_size)
  k1, k2 = _get_bounding_box(coords[:,i,:])
  bead_density.construct_kernel(k1,k2)
  return bead_density.return_density_opt(coords[:,i,:], radius[i], mass[i], n_breaks)

if __name__ == '__main__':
	parser = argparse.ArgumentParser("PrISM")
	parser.add_argument("--input", "-i", help="Input npz file", required=True, type=str)
	parser.add_argument("--voxel_size", "-v", help="Voxel size for density calculations", default=4, type=int)
	parser.add_argument("--return_spread", "-rs", help="Return the spread bead_spread", action='store_true', default = True)
	parser.add_argument("--output", "-o", help="Output directory", required = True, type=str)
	parser.add_argument("--classes", "-cl", help="Number of classes", default = 2, type=int)
	parser.add_argument("--cores", "-co", help="Number of cores to use", default = 16, type=int)
	parser.add_argument("--models", "-m", help="Percentage of total models to use", default = 1, type=float)
	parser.add_argument("--n_breaks", "-n", help="Number of breaks to use for cDist calculation", default = 50, type=int)
	args = parser.parse_args()

	arr = np.load(args.input)
	coords = arr['arr_0']
	mass = arr['arr_1']
	radius = arr['arr_2']
	ps_names = arr['arr_3']
	del arr
	models = round(args.models*coords.shape[0])
	if args.models != 1:
		selected_models = np.random.choice(coords.shape[0], models, replace=False)
		coords = coords[selected_models]
	print("Number of Models = {}".format(coords.shape[0]))
	print("Number of Beads = {}".format(coords.shape[1]))

	grid = SparseGrid(voxel_size=args.voxel_size)
	grid.create_grid(coords)
	grid.pad_grid(0)

	with Pool(args.cores) as p:
	    densities = p.map( partial(main_density_calc, grid=grid, voxel_size=args.voxel_size, n_breaks=args.n_breaks), range(0, coords.shape[1] ))

	print('Density calculation done')

	with Pool(args.cores) as p:
		bead_spread = p.map( partial(calc_bead_spread, grid=grid), densities) 

	print('Bead Spread calculation done')

	if not os.path.exists(args.output):
		os.makedirs(args.output)

	if args.return_spread == 1:
		with open(args.output + "/bead_spreads_cl" + str(args.classes) + ".txt", "w") as fl:
			for bs in bead_spread:
				fl.write(str(bs))
				fl.write("\n")

	patches = get_patches(bead_spread, args.classes, coords, radius)
	annotated_patches = annotate_patches(patches, args.classes, ps_names, coords.shape[1])
	high_prec, low_prec = patches[:args.classes], patches[args.classes+1:]
	print('Patch computing done')

	annot_df = pd.DataFrame(np.array(annotated_patches), columns = ['Bead', 'Bead Name', 'Type', 'Class', 'Patch'])
	annot_df['Patch'] = pd.to_numeric(annot_df["Patch"])
	annot_df.sort_values(['Patch'], ascending=[True])
	annot_df.to_csv(args.output + '/annotations_cl' + str(args.classes) + '.txt', index=None)

	with open(args.output + "/low_prec_cl" + str(args.classes) + ".txt", "w") as fl:
		lev = 1
		fl.write("Level" + "\t" + "Bead Indices" + "\t" + "Bead Names")
		fl.write("\n")
		for level in low_prec:
			for l in level:
				fl.write(str(lev))
				fl.write("\t")
				fl.write(",".join(str(item) for item in l))
				fl.write("\t")
				fl.write(",".join(ps_names[name] for name in l))
				fl.write("\n")
			lev=lev+1

	with open(args.output + "/high_prec_cl" + str(args.classes) + ".txt", "w") as fl:
		lev = 1
		fl.write("Level" + "\t" + "Bead Indices" + "\t" + "Bead Names")
		fl.write("\n")
		for level in high_prec:
			for l in level:
				fl.write(str(lev))
				fl.write("\t")
				fl.write(",".join(str(item) for item in l))
				fl.write("\t")
				fl.write(",".join(ps_names[name] for name in l))
				fl.write("\n")
			lev=lev+1









