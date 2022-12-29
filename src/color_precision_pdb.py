from Bio.PDB import PDBIO, PDBParser
from Bio.PDB import MMCIFIO, MMCIFParser
import os
import numpy as np
import argparse

def get_structure(fl, cif = True):
	if cif:
		parser=MMCIFParser() 
		modelStructure=parser.get_structure('complexStruct', fl)
	else:
		parser=PDBParser() 
		modelStructure=parser.get_structure('complexStruct', fl)
	return modelStructure[0]

def get_chain_length(model):
	res_no = 0
	for chain in model:
		for r in chain.get_residues():
			if r.id[0] == ' ':
				res_no +=1
	return res_no


def scale(v):
	return (v - v.min()) / (v.max() - v.min())

def save_structure(struct, output, cif = False):
	if cif:
		io = MMCIFIO()
	else:		
		io = PDBIO()
	
	io.set_structure(struct)
	io.save(output)

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

	if args.input_type == "cif":
		struct = get_structure(args.input, cif = True)
	else:
		struct = get_structure(args.input, cif = False)

	pf = open(args.precision_file,'r')
	precisions = [float(pr) for pr in pf.readlines()]
	pf.close()
	precisions = scale(np.array(precisions))
	precisions = [round(p,2) for p in precisions]
	if len(precisions) != get_chain_length(struct):
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

