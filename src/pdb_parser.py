from Bio.PDB import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
import numpy as np
import glob
import os

MASS = {'ALA': 71.03, 'ARG': 156.10, 'ASN':114.04, 'ASP': 115.02, 'CYS': 103.00, 'GLU': 129.04, 'GLN': 128.05, 'GLY': 57.02,
'HIS': 137.05, 'ILE': 113.08, 'LEU': 113.08, 'LYS':128.9, 'MET': 131.04, 'PHE': 147.06, 'PRO':97.05, 'SER': 87.03, 'THR': 101.04,
'TRP': 186.07, 'TYR': 163.06, 'VAL': 99.06, 'TPO': 101.04, 'SEP': 87.03}

RADIUS = {'ALA': 3.20, 'ARG': 5.60, 'ASN': 4.04, 'ASP': 4.04,'CYS': 3.65, 'GLN': 4.64, 'GLU': 4.63, 'GLY':1.72, 'HIS':4.73, 'ILE':3.94, 
'LEU': 4.24, 'LYS': 5.02, 'MET': 4.47, 'PHE': 4.99, 'PRO': 3.61, 'SER': 3.39, 'THR': 3.56, 'TRP': 5.38, 'TYR': 5.36, 'VAL':3.55, 
'TPO': 3.56, 'SEP': 3.39}

def get_chain_length(model):
	res_no = 0
	for chain in model:
		for r in chain.get_residues():
			if r.id[0] == ' ':
				res_no +=1
	return res_no
    
def get_features( modelFile, cif = False ):
	
	if cif:
		parser = MMCIFParser()
		models=parser.get_structure('cif',modelFile)
	else:
		parser=PDBParser()
		models=parser.get_structure('pdb',modelFile)
	
	coords = []
	mass = []
	radius = []
	res_ids = []
	for i, model in enumerate(models):
		for chain in model:
			for residue in chain.get_residues():
				if residue.id[0] == ' ':
					atom = residue['CA']
					coords.append(atom.get_coord())
					if i == 0:
						mass.append(MASS[residue.get_resname()])
						radius.append(RADIUS[residue.get_resname()])
						t_id = residue.full_id[2] + '_' + str(residue.id[1])
						res_ids.append(t_id)
	chain_len = get_chain_length(models[0])
	coords = np.array(coords).reshape(-1, chain_len, 3)
	mass = np.array(mass).reshape(chain_len, 1)
	radius = np.array(radius).reshape(chain_len, 1)
	return [coords, mass, radius, res_ids]


def parse_all_struct(folder, _type = "pdb" ):
	if _type == "cif":
		files = glob.glob(os.path.join(folder, '*.cif'))

		if len(files) == 0:
			print('No CIF file found')
		all_feats = [get_features(i, cif = True) for i in files]
		coords = np.vstack([i[0] for i in all_feats])
	
	elif _type == "pdb":
		files = glob.glob(os.path.join(folder, '*.pdb'))

		if len(files) == 0:
			print('No PDB file found')
		all_feats = [get_features(i, cif = False) for i in files]
		coords = np.vstack([i[0] for i in all_feats])

	mass = all_feats[0][1]
	radius = all_feats[0][2]
	res_ids = all_feats[0][3]
	del all_feats
	return coords, mass, radius, res_ids