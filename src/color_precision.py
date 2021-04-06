from __future__ import print_function
from IMP import ArgumentParser
import os
import IMP
import RMF
import IMP.rmf

__doc__ = "Color the beads of the cluster representative model based on the precision."

###########################################################

def parse_args():
    parser = ArgumentParser(
            description="Color the beads of a cluster center model based on the bead precision as output from PRISM")
    parser.add_argument('--path', '-p', dest="path",
            help='path to the good-scoring models', default="./")
    parser.add_argument('--clusternum', '-c', dest="cluster_number",type=int,
            help='cluster index of model to be colored; required argument', default=0)
    parser.add_argument('--subunit', '-su', dest="subunit",
            help='calculate RMSD/sampling and cluster precision/densities '
                 'etc over this subunit only', default=None)
    parser.add_argument('--resolution', '-r', dest="resolution", type=int,
            help='bead size (residues per bead) for RMSD and precision calculations; '
                 'use biggest bead size for faster calculations', default=1)
    parser.add_argument('--precision_file','-pf',dest="precision_file",required=True,type=str,
            help='location of file containing precision value per bead; one precision per line; required argument')

    return parser.parse_args()

def get_bead_name(p):
    ''' Input: particle
    Output: bead name in the format moleculename_copynumber_startresidue_endresidue
    '''

    mol_name = IMP.atom.get_molecule_name(IMP.atom.Hierarchy(p))

    copy_number=IMP.atom.get_copy_index(IMP.atom.Hierarchy(p))

    if IMP.atom.Fragment.get_is_setup(p):
        residues_in_bead = IMP.atom.Fragment(p).get_residue_indexes()

        bead_name = mol_name+":"+str(copy_number)+":"+str(min(residues_in_bead))+"-"+str(max(residues_in_bead))

    else:
        residue_in_bead = str(IMP.atom.Residue(p).get_index())

        bead_name = mol_name+":"+str(copy_number)+":"+residue_in_bead+"-"+residue_in_bead

    return bead_name

def main():
    args = parse_args()

    # Open representative model
    cluster_centroid = os.path.join(args.path,"cluster.%s/cluster_center_model.rmf3" %str(args.cluster_number))

    if not os.path.exists(cluster_centroid):
        print("Cluster centroid file not found for cluster %d.",args.cluster_number)
        exit(1)

    rmf_fh = RMF.open_rmf_file_read_only(cluster_centroid)

    # Build hierarchy from the RMF file
    m = IMP.Model()
    h = IMP.rmf.create_hierarchies(rmf_fh, m)[0]
    IMP.rmf.load_frame(rmf_fh, 0)
    m.update()

    # Select particles to color by precision
    if args.subunit:
        s0 = IMP.atom.Selection(h, resolution=args.resolution, molecule=args.subunit)
    else:
        s0 = IMP.atom.Selection(h, resolution=args.resolution)

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
    precisions_out_file = open(os.path.join(args.path,"cluster.%s/bead_precisions.txt"  %str(args.cluster_number)),'w')

    for i,leaf in enumerate(particles):

        #ASSUMPTION Assuming single state models for now
        # One can make this multi-state by adding state name in the bead name

        bead_name = get_bead_name(leaf)

        # Create a new particle
        p_new = m_new.add_particle(bead_name)

        # Decorate it with the same XYZR (sphere) as original particle
        xyzr_new = IMP.core.XYZR.setup_particle(m_new,p_new,IMP.core.XYZR(leaf).get_sphere())

        # Color particle based on its precision
        c_new  = IMP.display.Colored.setup_particle(m_new,p_new,IMP.display.get_hot_color(precisions[i]))
        # other options include get_gray_color, get_gnuplot_color

        # Add particle to the new model's hierarchy
        h_new = IMP.atom.Hierarchy.setup_particle(m_new,p_new)
        h_root.add_child(h_new)

        # Also, output bead name and precision to a text file so user can see the values
        print(bead_name,precisions[i], file=precisions_out_file)

    precisions_out_file.close()

    # Now create a new RMF file with the new model
    rmf_new = RMF.create_rmf_file(os.path.join(args.path,"cluster.%s/precision_colored_cluster_center_model.rmf3" %str(args.cluster_number)))

    IMP.rmf.add_hierarchy(rmf_new,h_root)
    IMP.rmf.save_frame(rmf_new)

    del rmf_fh
    del rmf_new

if __name__ == '__main__':
    main()
