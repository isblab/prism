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


def _get_number_of_beads(input_type,input_file, resolution, subunit):
    m = IMP.Model()

    if input_type =="rmf":
        inf = RMF.open_rmf_file_read_only(input_file)
        h = IMP.rmf.create_hierarchies(inf, m)[0]
        IMP.rmf.load_frame(inf, 0)
        m.update()

        if subunit:
            s0 = IMP.atom.Selection(h, resolution=resolution, molecule=subunit)

        else:
            s0 = IMP.atom.Selection(h, resolution=resolution)


    elif input_type == "pdb":

        h = IMP.atom.read_pdb(input_file, m, IMP.atom.CAlphaPDBSelector())

        s0 = IMP.atom.Selection(h)

    return (len(s0.get_selected_particles()))



def get_coordinates(input_type, path, output_base_path, output_path, resolution,
                    subunit):  # path to directory containing RMF/PDB files

    if input_type=="rmf":
        input_suffix = "rmf3"
    elif input_type=="pdb":
        input_suffix = "pdb"

    num_models = len(glob.glob("{}/*.".format(path) + input_suffix))

    num_beads = _get_number_of_beads(input_type, glob.glob("{}/*.".format(path)+input_suffix)[0], resolution, subunit)

    with open(os.path.join(output_base_path, 'meta_info.txt'), 'w') as f:
        f.write('Number of Models: {} \n Number of bead in each model: {}'.format(num_models, num_beads))

    # for mod_id,str_file in tqdm(enumerate(sorted(glob.glob("%s/*.rmf3" % (path)),key=lambda x:int(x.split('/')[-1].split('.')[0])))):
    def get_coordinates(input_type,mod_id, str_file):
        conform = np.empty([num_beads, 3])
        radii = None
        m = IMP.Model()

        if input_type =="rmf":
            inf = RMF.open_rmf_file_read_only(str_file)
            h = IMP.rmf.create_hierarchies(inf, m)[0]
            IMP.rmf.load_frame(inf, 0)
            m.update()

            if subunit:
                s0 = IMP.atom.Selection(h, resolution=resolution, molecule=subunit)

            else:
                s0 = IMP.atom.Selection(h, resolution=resolution)

        elif input_type == "pdb":

            h = IMP.atom.read_pdb(str_file, m, IMP.atom.CAlphaPDBSelector())

            s0 = IMP.atom.Selection(h)

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

    path = glob.glob("{}/*.".format(path) + input_suffix)
    sorted_paths = sorted(path, key=lambda x: int(x.split('/')[-1].split('.')[0]))

    Parallel(n_jobs=-1)(
        delayed(get_coordinates)(input_type, mod_id, str_file) for mod_id, str_file in tqdm(enumerate(sorted_paths)))

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
