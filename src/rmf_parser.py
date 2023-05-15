import numpy as np
import sys, os, glob
from functools import partial
import IMP
import IMP.atom
import IMP.rmf
import IMP.sampcon.precision_rmsd
import IMP.sampcon.rmsd_calculation
import RMF
import argparse
from multiprocessing import Pool


def get_selected_particles(m, input_file, frame_index, input_type, resolution, subunit=None,selection=None):
    s0 = None
    if input_type =="rmf":
        inf = RMF.open_rmf_file_read_only( input_file )
        h = IMP.rmf.create_hierarchies( inf, m )[0]
        IMP.rmf.load_frame( inf, frame_index )
        m.update()
        resolution = float( resolution )
        if subunit:
            s0 = IMP.atom.Selection( h, resolution=resolution, molecule=subunit )
        elif selection:
            s0 = IMP.sampcon.rmsd_calculation.parse_rmsd_selection( h, selection, resolution )
        else:
            s0 = IMP.atom.Selection( h, resolution=resolution )
        del inf
        selected_particles = []
        [selected_particles.append( p )  for p in s0.get_selected_particles() if "gaussian" not in str(p)]

    elif input_type == "pdb":
        h = IMP.atom.read_pdb(input_file, m, IMP.atom.CAlphaPDBSelector())
        s0 = IMP.atom.Selection(h)
        selected_particles = s0.get_selected_particles()
    return selected_particles

def _get_number_of_beads(input_file,input_type, resolution, subunit, selection):
    m = IMP.Model()
    s0 = get_selected_particles(m,input_file, input_type, resolution, subunit, selection)
    if s0:
        return (len(s0))
    return 0


def get_bead_name(p, input_type):
    ''' Input: particle
    Output: bead name in the format moleculename_copynumber_startresidue_endresidue
    '''
    if input_type=="rmf":
        mol_name = IMP.atom.get_molecule_name(IMP.atom.Hierarchy(p))
        copy_number=IMP.atom.get_copy_index(IMP.atom.Hierarchy(p))
        if IMP.atom.Fragment.get_is_setup(p):
            residues_in_bead = IMP.atom.Fragment(p).get_residue_indexes()
            bead_name = mol_name+":"+str(copy_number)+":"+str(min(residues_in_bead))+":"+str(max(residues_in_bead))
        else:
            residue_in_bead = str(IMP.atom.Residue(p).get_index())
            bead_name = mol_name+":"+str(copy_number)+":"+residue_in_bead+":"+residue_in_bead
        return bead_name

    elif input_type =="pdb":
        mol_name = IMP.atom.get_molecule_name(IMP.atom.Hierarchy(p))
        residue_number = IMP.atom.Residue(IMP.atom.Hierarchy(p).get_parent()).get_index()
        bead_name = mol_name+":"+str(residue_number)
        return bead_name



def get_coordinates(str_file, frame_index, input_type, resolution, subunit, selection):
    m = IMP.Model()
    s0 = get_selected_particles(m, str_file, 0, input_type,resolution, subunit, selection)
    conform = np.empty([ len(s0), 3])
    for i, leaf in enumerate(s0):
        p = IMP.core.XYZR(leaf)
        conform[i] = p.get_coordinates()  
    return conform

def get_attributes(str_file, input_type, resolution, subunit, selection):
    m = IMP.Model()
    s0 = get_selected_particles(m, str_file, 0, input_type, resolution, subunit, selection)
    radii = np.empty([ len(s0) ])
    mass = np.empty([len(s0)])
    bead_names = [0]*len(s0)
    for i, leaf in enumerate(s0):
        p = IMP.core.XYZR(leaf)
        radii[i] = p.get_radius()
        mass[i] = IMP.atom.Mass(leaf).get_mass()
        bead_names[i] = get_bead_name(leaf,input_type)
    return mass, radii, bead_names


def parse_all_rmfs(path, resolution, subunit, selection):
    # Read the selected particles.
    selection = parse_custom_ranges( selection )
    files_path = glob.glob(os.path.join(path, "*.rmf3" ))
    # with Pool(16) as p:
    #     coords = p.map(partial(get_coordinates, input_type="rmf", resolution=resolution, subunit=subunit, selection=selection), files_path)
    print('Files Detected: ', files_path)
    with Pool(16) as p:
        if len(files_path) > 1:
            coordinates = p.map(partial(get_coordinates, frame_index=0, input_type="rmf", resolution=resolution, subunit=subunit, selection=selection), files_path)
        if len(files_path) == 1:  # load all frames if only a single rmf found
            inf = RMF.open_rmf_file_read_only(files_path[0])
            all_frames = inf.get_number_of_frames()
            iterobj = p.imap(partial(get_coordinates, files_path[0], input_type='rmf', resolution=resolution, subunit=subunit, selection=selection), range(all_frames), chunksize=20)
            coordinates = []
            for coords in tqdm.tqdm(iterobj, total=all_frames):
                coordinates.append( coords )
    mass, radii, bead_names = get_attributes(files_path[0], input_type="rmf", resolution=resolution, subunit=subunit, selection=selection)
    return np.array(coordinates), np.array(mass), np.array(radii), np.array(bead_names)


def main(input_type, path, output_base_path, resolution, subunit, selection):  # path to directory containing RMF/PDB files
    if input_type=="rmf":
        input_suffix = "rmf3"
    elif input_type=="pdb":
        input_suffix = "pdb"
    files_path = glob.glob(os.path.join(path, "*." + input_suffix ))
    with Pool(16) as p:
        coords = p.map(partial(get_coordinates, input_type=input_type, resolution=resolution, subunit=subunit, selection=selection), files_path)
    mass, radii, bead_names = get_attributes(files_path[0], input_type=input_type, resolution=resolution, subunit=subunit, selection=selection)
    out_path = os.path.join(output_base_path, 'prism.npz')
    np.savez(out_path, np.array(coords), np.array(mass), np.array(radii), np.array(bead_names))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate PrISM inputs from rmf or pdb files')
    parser.add_argument('--input',
                        help="The path to rmf3 or files",
                        required=True)
    parser.add_argument('--type',
                        default='rmf',
                        help="rmf or pdb input.",
                        required=True)
    parser.add_argument('--output_dir',
                        help="Path to save the prism input files.",
                        required=True)
    parser.add_argument('--resolution',
                        help="The resolution at which to sample the beads.",
                        default=30,
                        required=False)
    parser.add_argument('--subunit',
                        help="Subunit that needs to be sampled.",
                        default="B",
                        required=False)
    parser.add_argument('--selection', 
                        help='File containing dictionary of selected subunits and residues.', 
                        default=None, 
                        required=False)

    args = parser.parse_args()
    output_base_path = args.output_dir
    if not os.path.exists(output_base_path):
        os.makedirs(output_base_path)

    main(args.type, args.input, output_base_path, args.resolution, args.subunit, args.selection)


