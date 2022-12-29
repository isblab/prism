import numpy as np
import os
import argparse

from ihm_parser import *


def ihm_set_bfactor( modelFile ):
	with open( modelFile ) as fh: # "../example/PDBDEV/PDBDEV_00000001/PDBDEV_00000001.cif"
		# s is an object of class ihm.System.
	    s, = ihm.reader.read(fh)
	_id = []
	_id = [atom.seq_id for model in s._all_starting_models() for atom in model.get_atoms() if atom.atom_id == "CA"]

	i = 0
	for model in s._all_starting_models():
		for atom in model.get_atoms():
			if _id[i] == atom.seq_id:
				atom.biso = b
			elif _id[i] != atom.seq_id:
				i+=1
				b +=1
				atom.biso = b

	with open("modified.cif", "w") as fh:
		ihm.dumper.write( fh, [s] )

def get_chain_length_ihm( system ):
	# Obtain the no. of residues in the system.
	# s --> an object of class ihm.System.
	res_no = 0
	for model in s._all_starting_models():
		for atom in model.get_atoms():
			if atom.atom_id == "CA":
				res_no += 1
	return res_no


def get_patch_coloured_model( annotations ):
	with open( annotations,'r') as pf:
		precisions = [pr.split(',')[2:4] for pr in pf.readlines()[1:]]

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





def main():
	parser = argparse.ArgumentParser("Coloring pdb/cif files")
	parser.add_argument('--input', '-i',
			help='representative model in PDB/CIF format',required=True)
	parser.add_argument("--input_type", "-t", help="Type of input: pdb/mmcif", required=True, type=str)
	parser.add_argument('--precision_file','-pf',required=True,type=str,
			help='location of bead spread outputs')
	parser.add_argument('--output', '-o',
			help='precision-colored model in PDB format. Visualize using Chimera.', default="precision_colored_cluster_center_model.pdb")
	args = parser.parse_args()

	if not os.path.exists(args.input):
		print("Cluster representative file not found %s",args.input)
		exit(1)

	struct = get_structure(args.input, args.input_type)
	pf = open(args.precision_file,'r')
	precisions = [float(pr) for pr in pf.readlines()]
	pf.close()
	precisions = scale(np.array(precisions))
	precisions = [round(p,2) for p in precisions]
	if len(precisions) != get_chain_length( struct ):
		print("Number of precision values: {} is not equal to the number of particles {}".format(len(precisions),get_chain_length(struct) ))
	
	count = 0
	for chain in struct:
		for residue in chain.get_residues():
			if residue.id[0] == ' ':
				for atom in residue:
					atom.set_bfactor(precisions[count])
			count += 1
	if args.input_type == "cif":
		save_structure(struct, args.output, cif = True)
	else:
		save_structure(struct, args.output, cif = False)


if __name__ == '__main__':
	main()

