import argparse
import os

import IMP.sampcon.precision_rmsd

from omegaconf import OmegaConf

from .generate_input.generate_distance_vectors import run as generate_dist_vectors
from .generate_input.get_model_coordinates import get_coordinates
from .train import start_training

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Train PrISM on given ensemble of integrative structure models.')

    parser.add_argument('--input',
                        help="The input file (if npz input) or directory containing input files (if rmf/pdb input)",
                        required=True)
    parser.add_argument('--output_dir',
                        help="Path to save the intermediate and final output files",
                        required=True)
    parser.add_argument('--type',
                        default='npz',
                        help="rmf, pdb or npz input.")
    parser.add_argument('--resolution',
                        help="The resolution to sample the beads",
                        default=1,
                        required=False)
    parser.add_argument('--subunit',
                        help="Subunit (protein) for which precision is to be calculated (applicable only on rmf or pdb input types). Only one subunit can be specified. By default precision is calculated on all subunits.",
                        default=None,
                        required=False)
    parser.add_argument(
        '--selection', '-sn', dest="selection",
        help='file containing dictionary'
        'of selected subunits and residues'
        'for RMSD and clustering calculation'
        "each entry in the dictionary takes the form"
        "'selection name': [(residue_start, residue_end, protein name)]. See example.",
        default=None)
    parser.add_argument('--config',
                        help="Config file containing details of training parameters",
                        default='test_config.yml')
    parser.add_argument('--skip_input_generation',
                        help="If the input has already been generated and you would like to train with different hyperparamters set this flag to 1",
                        default=0)
    parser.add_argument('--gpu',
                        help="Set to 1, to use GPU.",
                        type=int, default=0)
    args = parser.parse_args()

    # Get model coordinates.
    output_base_path = args.output_dir
    output_path = os.path.join(output_base_path, 'models_coordinate_npz')
    os.makedirs(output_path, exist_ok=True)
    type = args.type

    if args.skip_input_generation != 1:
        npz_path = None

        if args.selection:
            rmsd_custom_ranges = IMP.sampcon.precision_rmsd.parse_custom_ranges(args.selection)

        else:
            rmsd_custom_ranges = None

        if type == 'rmf':
            resolution = int(args.resolution)

            print("Getting the bead coordinates from RMF files with following parameters")
            print("Input: {}, Output: {}".format(args.input, output_base_path))
            print("Resolution: {}, Selected single subunit: {}, Selected group of subunits: {} (if both selections are 'None', by default all subunits are selected).".format(
                resolution, args.subunit, args.selection))

            npz_path = get_coordinates(type, args.input, output_base_path, output_path, resolution, args.subunit,rmsd_custom_ranges)
            type = 'npz'

        elif type == 'pdb':
            resolution = int(args.resolution)
            print("Getting the CA  coordinates from PDB files with following parameters")
            print("Input: {}, Output: {}".format(args.input, output_base_path))
            npz_path = get_coordinates(type, args.input, output_base_path, output_path, resolution, args.subunit,rmsd_custom_ranges)
            type = 'npz'

        # Generate bead_wise distance.
        if type == 'npz':
            if not npz_path:
                npz_path = args.input
            print("Generating bead-wise distances")
            beadwise_distances_path = generate_dist_vectors(npz_path, args.output_dir)
    else:
        beadwise_distances_path = args.input

    conf = OmegaConf.load(args.config)

    # Add input and output directory details
    conf['dataset'] = {}
    conf['dataset']['input_dir'] = beadwise_distances_path
    conf['dataset']['output_dir'] = args.output_dir
    conf['use_gpu'] = True if args.gpu == 1 else False

    print("Starting training")
    start_training(conf, args.gpu)
