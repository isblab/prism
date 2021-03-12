import numpy as np
import sys, os, glob
from joblib import Parallel, delayed

import IMP
import IMP.atom
import IMP.rmf
import RMF
from tqdm import tqdm
import argparse
import warnings


def _get_number_of_beads(rmf_file, resolution, molecule):
    m = IMP.Model()
    inf = RMF.open_rmf_file_read_only(rmf_file)
    h = IMP.rmf.create_hierarchies(inf, m)[0]
    IMP.rmf.load_frame(inf, 0)

    # s0 = IMP.atom.Selection(h, resolution=resolution, molecule=molecule)
    s0 = IMP.atom.Selection(h, resolution=resolution)

    return (len(s0.get_selected_particles()))


# TODO
# Get number of residues.

def get_coordinates(path, output_base_path, output_path, resolution,
                    molecule):  # path to directory containing RMF3 files

    num_models = len(glob.glob("%s/*.rmf3" % (path)))
    num_beads = _get_number_of_beads(glob.glob("%s/*.rmf3" % (path))[0], resolution, molecule=molecule)
    with open(os.path.join(output_base_path, 'meta_info.txt'), 'w') as f:
        f.write('Number of Models: {} \n Number of bead in each model: {}'.format(num_models, num_beads))

    # for mod_id,str_file in tqdm(enumerate(sorted(glob.glob("%s/*.rmf3" % (path)),key=lambda x:int(x.split('/')[-1].split('.')[0])))):
    def get_coordinates(mod_id, str_file):
        conform = np.empty([num_beads, 3])
        radii = None
        m = IMP.Model()

        # print str_file
        # print(str_file, mod_id)

        inf = RMF.open_rmf_file_read_only(str_file)
        h = IMP.rmf.create_hierarchies(inf, m)[0]
        IMP.rmf.load_frame(inf, 0)

        m.update()

        # s0 = IMP.atom.Selection(h, resolution=resolution, molecule=molecule)
        s0 = IMP.atom.Selection(h, resolution=resolution)

        for i, leaf in enumerate(s0.get_selected_particles()):
            # print leaf, make sure multiple beads are not present in multiscale systems
            p = IMP.core.XYZR(leaf)
            pxyz = p.get_coordinates()

            conform[i][0] = pxyz[0]
            conform[i][1] = pxyz[1]
            conform[i][2] = pxyz[2]
            # print leaf,conform[mod_id][i]
            radii = p.get_radius()

        np.savez(os.path.join(output_path, "{}.npz".format(mod_id)), conform, radii)
        return conform, radii

    path = glob.glob("%s/*.rmf3" % (path))
    sorted_paths = sorted(path, key=lambda x: int(x.split('/')[-1].split('.')[0]))
    # print(sorted_paths[449])
    # sys.exit(0)
    Parallel(n_jobs=-1)(
        delayed(get_coordinates)(mod_id, str_file) for mod_id, str_file in tqdm(enumerate(sorted_paths)))

    conform = np.empty([num_models, num_beads, 3])
    radii = np.empty([num_models])
    # This step is an additional overhead because while using joblib, the same array cannot be manipulated.
    # Adding it as a separate step will have less time overhead as compared to running everything in a for-loop.
    output_models = glob.glob(output_path + "/*")
    sorted_paths = sorted(output_models, key=lambda x: int(x.split('/')[-1].split('.')[0]))
    for index, model_path in enumerate(sorted_paths):
        model = np.load(model_path)
        conform[index] = model['arr_0']
        radii[index] = model['arr_1']
    out_path = os.path.join(output_base_path, 'all_models_coordinates.npz')
    np.savez(out_path, conform, radii)
    return out_path


if __name__ == '__main__':
    warnings.warn("This module is deprected.")
    parser = argparse.ArgumentParser(description='Generate distance matrices')
    parser.add_argument('--input',
                        help="The path to rmf3 or files",
                        required=True)
    parser.add_argument('--type',
                        default='rmf',
                        help="rmf or pdb input.",
                        required=True)
    parser.add_argument('--output_dir',
                        help="Path to save the model_cooridnate files",
                        required=True)
    parser.add_argument('--resolution',
                        help="The resolution to sample the beads",
                        default=1,
                        required=False)
    parser.add_argument('--molecule',
                        help="Molecule that needs to be sampled",
                        default="B",
                        required=False)
    # TODO
    # Take PDB or RMF as input

    args = parser.parse_args()
    output_base_path = args.output_dir
    output_path = os.path.join(output_base_path, 'models_coordinate_npz')
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    resolution = int(args.resolution)

    get_coordinates(args.input, output_base_path, output_path, resolution, args.molecule)