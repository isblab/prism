from __future__ import print_function
from IMP import ArgumentParser
import os
import IMP
import RMF
import IMP.rmf
import IMP.sampcon.precision_rmsd
import IMP.sampcon.rmsd_calculation

from generate_input.get_model_coordinates import get_selected_particles

__doc__ = "Get an RMF where beads of the cluster representative model are colored \
based on their precision as reported by PRISM."

###########################################################

def parse_args():
    parser = ArgumentParser(
            description="Color the regions of a representative model based on the precision as output from PRISM")
    parser.add_argument('--input', '-i', dest="input",
            help='representative model in RMF or PDB format', default="cluster.0/cluster_center_model.rmf3",required=True)
    parser.add_argument('--subunit', '-su', dest="subunit",
            help='annotate variation in precision over this subunit only', default=None)
    parser.add_argument('--resolution', '-r', dest="resolution", type=int,
            help='bead size (residues per bead) for annotating precision', default=1)
    parser.add_argument(
        '--selection', '-sn', dest="selection",
        help='file containing dictionary'
        'of selected subunits and residues'
        'for RMSD and clustering calculation'
        "each entry in the dictionary takes the form"
        "'selection name': [(residue_start, residue_end, protein name)",
        default=None)
    parser.add_argument('--precision_file','-pf',dest="precision_file",required=True,type=str,
            help='location of output from the autoencoder; one precision per line; required argument')
    parser.add_argument('--output', '-o', dest="output",
            help='precision-colored model in RMF format. Visualize using Chimera.', default="precision_colored_cluster_center_model.rmf3")

    return parser.parse_args()

#Improvement This function is duplicated in IMP.sampcon.rmsd_calculation's get_rmfs_coordinates_one_rmf
# Consider moving to common function

def get_bead_name(p, input_type):
    ''' Input: particle
    Output: bead name in the format moleculename_copynumber_startresidue_endresidue
    '''

    if input_type=="rmf":

        mol_name = IMP.atom.get_molecule_name(IMP.atom.Hierarchy(p))

        copy_number=IMP.atom.get_copy_index(IMP.atom.Hierarchy(p))

        if IMP.atom.Fragment.get_is_setup(p):
            residues_in_bead = IMP.atom.Fragment(p).get_residue_indexes()

            bead_name = mol_name+":"+str(copy_number)+":"+str(min(residues_in_bead))+"-"+str(max(residues_in_bead))

        else:
            residue_in_bead = str(IMP.atom.Residue(p).get_index())

            bead_name = mol_name+":"+str(copy_number)+":"+residue_in_bead+"-"+residue_in_bead

        return bead_name

    elif input_type =="pdb":

        mol_name = IMP.atom.get_molecule_name(IMP.atom.Hierarchy(p))

        residue_number = IMP.atom.Residue(IMP.atom.Hierarchy(p).get_parent()).get_index()

        bead_name = mol_name+":"+str(residue_number)

        return bead_name

def main():
    args = parse_args()

    # Open representative model

    if not os.path.exists(args.input):
        print("Cluster representative file not found %s",args.input)
        exit(1)

    if args.selection:
        rmsd_custom_ranges = IMP.sampcon.precision_rmsd.parse_custom_ranges(args.selection)
    else:
        rmsd_custom_ranges = None
    # create model and hierarchy for input file
    m = IMP.Model()
    if args.input.lower().endswith("rmf") or args.input.lower().endswith("rmf3"):

        input_type = "rmf"

    elif args.input.lower().endswith("pdb"):

        input_type = "pdb"

    else:
        print("Input file is not in RMF or PDB format.")
        exit(1)

    s0 = get_selected_particles(m,args.input,input_type,args.resolution,args.subunit,rmsd_custom_ranges)

    particles = s0.get_selected_particles()

    # Now get precisions from file
    pf = open(args.precision_file,'r')
    precisions = [float(pr.strip()) for pr in pf.readlines()]
    pf.close()

    if len(precisions)!=len(particles):
        print("Number of precision values is not equal to the number of selected particles.Check that the same selection arguments (e.g. subunit/resolution) are used in the sampcon code that generated the superposed particles, and this code")
        exit(1)

    # Create a new RMF with beads of same XYZR as selected particles, but colored according to precision
    # This is a reduced version of the representative cluster center model
    m_new = IMP.Model()

    # Create new hierarchy to add new beads to
    p_root = m_new.add_particle("System")
    h_root = IMP.atom.Hierarchy.setup_particle(m_new,p_root)

    # Simultaneously output precisions with bead names to a text file for user to see
    precisions_out_file = open('bead_precisions.txt','w')

    # default protein
    prev_prot ="DUMMY.0"


    for i,leaf in enumerate(particles):

        #ASSUMPTION Assuming single state models for now
        # One can make this multi-state by adding state name in the bead name

        bead_name = get_bead_name(leaf,input_type)

        ## see if we need to create a new protein
        prot_base_name = bead_name.split(':')[0]
        copy_number = bead_name.split(':')[1]
        curr_prot = prot_base_name+"."+copy_number

        if curr_prot != prev_prot:
            # Create a new hierarchy particle for the protein
            p_curr_prot = m_new.add_particle(curr_prot)

            # Add to the new model's hierarchy
            h_curr_prot = IMP.atom.Hierarchy.setup_particle(m_new,p_curr_prot)
            h_root.add_child(h_curr_prot)

            prev_prot = curr_prot

        # Create a new particle
        start_res = bead_name.split(':')[2]
        end_res = bead_name.split(':')[3]

        if end_res ==start_res:
            particle_name = start_res
        else:
            particle_name = start_res+"-"+end_res
            
        p_new = m_new.add_particle(particle_name)

        # Decorate it with the same XYZR (sphere) as original particle
        xyzr_new = IMP.core.XYZR.setup_particle(m_new,p_new,IMP.core.XYZR(leaf).get_sphere())

        # Color particle based on its precision
        c_new  = IMP.display.Colored.setup_particle(m_new,p_new,IMP.display.get_hot_color(precisions[i]))
        # other options include get_gray_color, get_gnuplot_color

        # Add particle to the new model's hierarchy
        h_new = IMP.atom.Hierarchy.setup_particle(m_new,p_new)
        h_curr_prot.add_child(h_new)

        # Also, output bead name and precision to a text file so user can see the values
        print(bead_name,"%.4f" % precisions[i], file=precisions_out_file)

    precisions_out_file.close()

    # Now create a new RMF file with the new model
    rmf_new = RMF.create_rmf_file(args.output)

    IMP.rmf.add_hierarchy(rmf_new,h_root)
    IMP.rmf.save_frame(rmf_new)

    del rmf_new

if __name__ == '__main__':
    main()
