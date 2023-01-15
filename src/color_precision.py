from IMP import ArgumentParser
import os
import IMP
import RMF
import IMP.rmf
import IMP.sampcon.precision_rmsd
import IMP.sampcon.rmsd_calculation
import numpy as np
from rmf_parser import *

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
    parser.add_argument('--annotations_file','-a',dest="annot",required=True,type=str,
            help='location of annotation outputs')
    parser.add_argument('--output', '-o', dest="output",
            help='precision-colored model in RMF format. Visualize using Chimera.', default="precision_colored_cluster_center_model.rmf3")

    return parser.parse_args()


def colour_rmf( tp, m_new, p_new ):
    low_cols = [(1,0,0), (1,0.5, 0.5), (1,0.75,0.75)] #red 
    high_cols = [(0,0.39,0), ( 0.20,0.80,0.20), (0.56,0.93,0.56)] #green

    if type(tp) == int:
        c_new  = IMP.display.Colored.setup_particle(m_new,p_new,IMP.display.Color(1, 1, 1))
    elif tp[0] == 'low':
        c_ind = int(tp[1]) - 1
        c_rgb = low_cols[c_ind]
        c_new  = IMP.display.Colored.setup_particle(m_new,p_new,IMP.display.Color(c_rgb[0], c_rgb[1], c_rgb[2]))
    elif tp[0] == 'high':
        c_ind = int(tp[1]) - 1
        c_rgb = high_cols[c_ind]
        c_new  = IMP.display.Colored.setup_particle(m_new,p_new,IMP.display.Color(c_rgb[0], c_rgb[1], c_rgb[2]))
    else:
        c_new  = IMP.display.Colored.setup_particle(m_new,p_new,IMP.display.Color(1, 1, 1))

    return c_new, m_new, p_new


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

    m1 = IMP.Model()
    sTotal = get_selected_particles(m1,args.input,input_type, args.resolution,None,None)
    
    s0_particles = s0
    sTotal_particles = sTotal
    s0_beads = [get_bead_name(i, input_type) for i in s0_particles]
    sTotal_beads = [get_bead_name(i, input_type) for i in sTotal_particles]
    

    ind_dict = dict( (k,0) for k in sTotal_beads)
    pf = open(args.annot,'r')
    precisions = [pr.split(',')[2:4] for pr in pf.readlines()[1:]]
    pf.close()

    if len(precisions)!=len(s0_particles):
        print("Number of precision values: " + str(len(precisions)) + " is not equal to the number of selected particles " + str(len(s0_particles)) +".Check that the same selection arguments (e.g. subunit/resolution) are used in the sampcon code that generated the superposed particles, and this code")
        exit(1)
    low_cols = [(1,0,0), (1,0.5, 0.5), (1,0.75,0.75)] #red 
    high_cols = [(0,0.39,0), ( 0.20,0.80,0.20), (0.56,0.93,0.56)] #green

    
    for i,prec in enumerate(precisions):
        ind_dict[s0_beads[i]] = prec

    # Create a new RMF with beads of same XYZR as selected particles, but colored according to precision
    # This is a reduced version of the representative cluster center model
    m_new = IMP.Model()

    # Create new hierarchy to add new beads to
    p_root = m_new.add_particle("System")
    h_root = IMP.atom.Hierarchy.setup_particle(m_new,p_root)

    # Simultaneously output precisions with bead names to a text file for user to see


    # default protein
    prev_prot ="DUMMY.0"

    
    for i,leaf in enumerate(sTotal_particles):
        #ASSUMPTION Assuming single state models for now
        # One can make this multi-state by adding state name in the bead name
        # Select only if the particle is not a gaussian.          
        bead_name = get_bead_name(leaf,input_type)

        #print(bead_name, precisions[i])
        ## see if we need to create a new protein
        prot_base_name = bead_name.split(':')[0]

        if input_type == "rmf":
            copy_number = bead_name.split(':')[1]
            curr_prot = prot_base_name+"."+copy_number
            start_res = bead_name.split(':')[2]
            end_res = bead_name.split(':')[3]
            if start_res ==end_res:
                res_range =start_res
            else:
                res_range = start_res+"-"+end_res
        elif input_type =="pdb":
            curr_prot = prot_base_name
            res_range = bead_name.split(':')[1]
        if curr_prot != prev_prot:
            # Create a new hierarchy particle for the protein
            p_curr_prot = m_new.add_particle(curr_prot)
            # Add to the new model's hierarchy
            h_curr_prot = IMP.atom.Hierarchy.setup_particle(m_new,p_curr_prot)
            h_root.add_child(h_curr_prot)
            prev_prot = curr_prot

        # create new particle
        p_new = m_new.add_particle(res_range)
        
        # Decorate with arbitrary mass
        mass_new = IMP.atom.Mass.setup_particle(m_new,p_new,1.0)

        # Decorate it with the same XYZR (sphere) as original particle
        xyzr_new = IMP.core.XYZR.setup_particle(m_new,p_new,IMP.core.XYZR(leaf).get_sphere())

        # Color particle based on its precision

        #IMP.display.color(r,g,b)   (255,0,0)
        tp = ind_dict[bead_name]
        c_new, m_new, p_new = colour_rmf( tp, m_new, p_new )
        # if type(tp) == int:
        #     c_new  = IMP.display.Colored.setup_particle(m_new,p_new,IMP.display.Color(1, 1, 1))
        # elif tp[0] == 'low':
        #     c_ind = int(tp[1]) - 1
        #     c_rgb = low_cols[c_ind]
        #     c_new  = IMP.display.Colored.setup_particle(m_new,p_new,IMP.display.Color(c_rgb[0], c_rgb[1], c_rgb[2]))
        # elif tp[0] == 'high':
        #     c_ind = int(tp[1]) - 1
        #     c_rgb = high_cols[c_ind]
        #     c_new  = IMP.display.Colored.setup_particle(m_new,p_new,IMP.display.Color(c_rgb[0], c_rgb[1], c_rgb[2]))
        # else:
        #     c_new  = IMP.display.Colored.setup_particle(m_new,p_new,IMP.display.Color(1, 1, 1))
        
        # other options include get_gray_color, get_gnuplot_color

        # Add particle to the new model's hierarchy
        h_new = IMP.atom.Hierarchy.setup_particle(m_new,p_new)
        h_curr_prot.add_child(h_new)

        # Also, output bead name and precision to a text file so user can see the values

    # Now create a new RMF file with the new model
    rmf_new = RMF.create_rmf_file(args.output)
    IMP.rmf.add_hierarchy(rmf_new,h_root)
    IMP.rmf.save_frame(rmf_new)
    del rmf_new

if __name__ == '__main__':
    main()
