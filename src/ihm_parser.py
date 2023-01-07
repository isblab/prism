import numpy as np
from math import pi
import os
import wget
import requests
import json
from collections import defaultdict

import ihm
import ihm.reader

import IMP
import IMP.atom
import RMF
from dcd_parser import *
from main import *
from pdb_parser import *


def check_software( mmcif ):
	# Check if the entry is IMP or non-IMP
	# Very naive check for IMP
	IMP_entry = False

	for software in mmcif.software:
		try:
			print( "Software used = ", software.name)
			if 'IMP' in software.name:
				IMP_entry = True
				print('IMP entry:', IMP_entry)
				break
		except TypeError:
			continue
	return IMP_entry


def check_dcd( mmcif ):
	# Check if the entry contains a dcd file or not.
	found_dcd = False
	ensembles = {}
	for i, e in enumerate(mmcif.ensembles):
		try:
			print(f'Ens: {i} | Models_deposited: {e.num_models_deposited} | Models_total: {e.num_models}')
		except AttributeError:
			print(f'Ens: {i} | Missing model_group')
		except TypeError:
			print(f'Ens: {i} | Models_deposited: N/A | Models_total: {e.num_models}')

	    # Sometimes we don't have info about where ensemble is stored
		fl = e.file
		if fl is None:
			print(f'Missing data for en: {i} {e.name}')
		else:
			found_dcd = True
	return found_dcd


def download_dcd( url, fl ):
	# Download the dcd file using the url.
	# url = "/".join(e.file.repo.url.split("/")[:-1])
	print( "Downloading dcd file..." )
	response = wget.download( url, fl )	


# Loop over all the model groups to get the required info.
def get_hierarchy_from_model( model ):
	"""return nested dictionaries of following format [chain_id][res_id][atom_name]"""
	infinite_defaultdict = lambda: defaultdict(infinite_defaultdict)
	root = infinite_defaultdict()
	count = 0

	for a in model.get_atoms():
		root[a.asym_unit.id][a.seq_id][a.atom_id] = a

	for s in model.get_spheres():
		# Consider only by-residue spheres
		seq_ids = list(set(s.seq_id_range))
		# if len(seq_ids) != 1:
		# 	continue

		seq_id = seq_ids[0]

		# if root[s.asym_unit.id][seq_id]['CA'] is None:
		count +=1
		root[s.asym_unit.id][seq_id] = s
	
	return(root)


def parse_models( system, model_group ):
	"""parse individiual models in the mmcif file"""
	# parse all models
	models = {}
	# models counter
	gimg = 0
	gim = 0
	# iterate over state groups
	for istg, stg in enumerate(system.state_groups):
		# iterate over states in state groups
		for ist, st in enumerate(stg):
			# iterate over model groups in states
			for img, mg in enumerate(st):
				gimg += 1
				if mg == model_group:
					# iterate over models in model groups
					for im, m in enumerate(mg):
						gim += 1
						m_ = get_hierarchy_from_model(m)
						models[gim] = m_
				else:
					continue
	print( "mg --> ", gimg, "\tm --> ", gim)
	return(models)


def get_all_attributes( models, IMP_ ):
	# Get the x, y, z coords and radius for all models.
	coords, radius, mass, ps_names = [], [], [], []
	avg_mass = 110
	volume = IMP.atom.get_volume_from_mass( avg_mass )
	avg_radius = 0.8 * (3.0 / 4.0 / pi * volume) ** (1.0 / 3.0)

	# iterate over all models in the model group.
	num_models = 0
	for m in models.keys():
		num_models += 1
		num_spheres = 0
		# iterate over all asym_units in a model.
		for asym_unit in models[m].keys():
			# iterate over all residues in an asym_unit.
			for res in models[m][asym_unit].keys():
				num_spheres += 1
				if IMP_:
					obj = models[m][asym_unit][res]
				else:
					obj = models[m][asym_unit][res]["CA"]
				# Check if it is a sphere or atom object.
				if isinstance(obj, ihm.model.Sphere) or isinstance(obj, ihm.model.Atom):
					coords.append( [obj.x, obj.y, obj.z])
					mass.append( avg_mass )
					ps_names.append( f"{m}_{asym_unit}_{res}" )
					try:
						radius.append( obj.radius )
					# non-IMP entries do not have a radius.
					except AttributeError:
						radius.append( avg_radius )						
						
	coords = np.array( coords ).reshape( num_models, num_spheres, 3 )
	radius = np.array( radius ).reshape( num_models*num_spheres, 1 )
	mass = np.array( mass ).reshape( num_models*num_spheres, 1 )

	return [coords, radius, mass, ps_names]


def get_patch_coloured_model( annotations, coords, radius ):
	with open( f"./{annotations}/annotations_cl2.txt",'r') as pf:
		file = pf.readlines()[1:]
	precisions = [pr.split(',')[2:4] for pr in file]
	particles = [pr.split(',')[1] for pr in file]

	low_cols = [(1,0,0), (1,0.5, 0.5), (1,0.75,0.75)] #red 
	high_cols = [(0,0.39,0), ( 0.20,0.80,0.20), (0.56,0.93,0.56)] #green

	ind_dict = {}

	for i,prec in enumerate( precisions ):
		ind_dict[particles[i]] = prec

	# Create a new RMF with beads of same XYZR as selected particles, but colored according to precision
	# This is a reduced version of the representative cluster center model
	m_new = IMP.Model()

	# Create new hierarchy to add new beads to
	p_root = m_new.add_particle("System")
	h_root = IMP.atom.Hierarchy.setup_particle( m_new, p_root )
    
	for i,bead_name in enumerate( particles ):
		# bead_name = model_asym_id_residue.
		res = bead_name.split( "_" )[2]

		# create new particle
		# Decorate with arbitrary mass, coordinates, radius.
		p_new = IMP.Particle( m_new, res )
		x, y, z = coords[0][i]
		p_new = IMP.atom.Mass.setup_particle( m_new, p_new, 110.0 )
		p_new = IMP.core.XYZR.setup_particle( p_new, IMP.algebra.Sphere3D( IMP.algebra.Vector3D(x, y, z), 2.7 ) )

		# Color particle based on its precision
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

		# Add particle to the new model's hierarchy
		h_new = IMP.atom.Hierarchy.setup_particle( m_new, p_new )
		h_root.add_child(h_new)

	# Now create a new RMF file with the new model
	rmf_new = RMF.create_rmf_file( f"./{annotations}/patch_coloured_model.rmf3" )
	IMP.rmf.add_hierarchy( rmf_new, h_root )
	IMP.rmf.save_frame( rmf_new )
	del rmf_new


def parse_all_models( args ):
	try:
	    with open( args.input, encoding='utf8') as fh:
	        mmcif, = ihm.reader.read(fh, model_class=ihm.model.Model)
	except UnicodeDecodeError:
	    with open( args.input, encoding='ascii', errors='ignore') as fh:
	        mmcif, = ihm.reader.read(fh, model_class=ihm.model.Model)

	IMP_ = check_software( mmcif )

	if IMP_:		
		for istg, stg in enumerate( mmcif.state_groups ):
			# iterate over states in state groups
			for ist, st in enumerate( stg ):
				# iterate over model groups in states
				for img, mg in enumerate( st ):
					found_dcd = check_dcd( mmcif )
					
					if found_dcd:	
						for i, e in enumerate(mmcif.ensembles):
							if e.model_group == mg:
								model_group = parse_models( mmcif, mg )
								
								url = e.file.repo.url
								fl = f"Ensemble{img}_{i}.dcd"
								if not os.path.isfile( fl ):
									download_dcd( url, fl )
								coords = DCDReader( fl ).return_coords()
								_, mass, radius, ps_names = get_all_attributes( model_group, IMP_ )

								if coords.shape[1] != radius.shape[0]:
									raise ValueError('Number of particles in cif file does not match \
										number of particles in dcd file {}'.format(coords.shape[0], radius.shape[0]))

								# Run PrISM t get bead precision.
								run_prism( coords, mass, radius, ps_names, args, f"output_{ist}_{img}" )

								# Get a patch coloured bead model.
					else:
						print( "IMP used for modelling. \n No dcd file found." )

	else:
		for istg, stg in enumerate( mmcif.state_groups ):
			# iterate over states in state groups
			for ist, st in enumerate( stg ):
				# iterate over model groups in states
				for img, mg in enumerate( st ):
					model_groups = parse_models( mmcif, mg )
					coords, mass, radius, ps_names = get_all_attributes( model_groups, IMP_ )
					# Run PrISM t get bead precision.
					run_prism( coords, mass, radius, ps_names, args, f"output_{ist}_{img}" )
					print( "Creating patch coloured model..." )
					get_patch_coloured_model( f"output_{ist}_{img}", coords, radius )













# exit()

# # Iterate over state groups
# for istg, stg in enumerate(mmcif.state_groups):
#     print(f'State group: {istg}')
#     # Iterate over states
#     for ist, st in enumerate(stg):
#         print(' ' * 20, f'State: {ist}')
#         # Iterate over model groups
#         for img, mg in enumerate(st):
#             # Basic info about model group
#             print(' ' * 40, f'Model group: {img} Size: {len(mg)}')
            
#             # Find corresponding ensemble
#             found_ensemble = False
#             for i, e in enumerate(mmcif.ensembles):
#                 if e.model_group == mg:
#                     found_ensemble = True
#                     break
            
#             # Detailed info about ensembles.
#             if found_ensemble:
#                 try:
#                     print(f'Associated ensemble: {i} | Models_deposited: {e.num_models_deposited} | Models_total: {e.num_models}')
#                 except AttributeError:
#                     print(f'Ens: {i} | Missing model_group')
#                 except TypeError:
#                     print(f'Ens: {i} | Models_deposited: N/A | Models_total: {e.num_models}')
#             else:
#                 print('No associated ensembles')            




# out = False
# for k in ["25"]:   # , "26", "37", "38", "41"
# 	with open( f"../example/PDBDEV/PDBDEV_000000{k}/PDBDEV_000000{k}.cif" ) as fh:
# 		# s is an object of class ihm.System.
# 	    s, = ihm.reader.read(fh)

# 	for i,e in enumerate( s.ensembles ):
# 		try:
# 			print( e.file.repo.url )
# 		except:
# 			print( f"No dcd file link present for {k}..." )
# 			out = True
# 			break

# 	if out:
# 		continue


# print( dir(e), "\n" )
# print(dir(e.model_group))
# # print( dir(e.file), "\n" )
# url = "/".join(e.file.repo.url.split("/")[:-1])
# response = requests.get( url )
# data = response.text 
# print("%s contains %d models at %s"% (e, e.num_models, e.file.repo.url))
# response = wget.download( e.file.repo.url, f"Ensemble{i}.dcd" )


# def ihm_cif_parser( modelFile ):
# 	# Parse an ihm cif file to extract the coords, radius, mass for all residues in all models.
# 	# modelFile --> cif file path.

# 	with open( modelFile ) as fh: # '../example/PDBDEV/PDBDEV_00000001/PDBDEV_00000001.cif'
# 		# s is an object of class ihm.System.
# 	    s, = ihm.reader.read(fh)

# 	coords, mass, radius, res_ids = [], [], [], []
# 	for model in s._all_starting_models():
# 		seq = model.asym_unit.entity.sequence
# 		for atom in model.get_atoms():
# 			if atom.atom_id == "CA":
# 				res_ids.append( seq[atom.seq_id-1].code + "_" + str( atom.seq_id ) )
# 				coords.append( [atom.x, atom.y, atom.z] )
# 				mass.append( MASS[seq[atom.seq_id-1].id] )
# 				radius.append( RADIUS[seq[atom.seq_id-1].id] )

# 	chain_len = len( coords )
# 	coords = np.array(coords).reshape(-1, chain_len, 3)
# 	radius = np.array( radius ).reshape( chain_len, 1 )
# 	mass = np.array( mass ).reshape( chain_len, 1 )
# 	return [coords, mass, radius, res_ids]


# def parse_all_ihm_cif(folder, _type = "pdb" ):
# 	if _type == "cif":
# 		files = glob.glob(os.path.join(folder, '*.cif'))

# 		if len(files) == 0:
# 			print('No CIF file found')
# 		all_feats = [get_features(i, cif = True) for i in files]
# 		coords = np.vstack([i[0] for i in all_feats])

# 	elif _type == "ihm":
# 		files = glob.glob(os.path.join(folder, '*.cif'))

# 		if len(files) == 0:
# 			print('No CIF file found')
# 		all_feats = [ihm_cif_parser( i ) for i in files]
# 		coords = np.vstack([i[0] for i in all_feats])
	
# 	elif _type == "pdb":
# 		files = glob.glob(os.path.join(folder, '*.pdb'))

# 		if len(files) == 0:
# 			print('No PDB file found')
# 		all_feats = [get_features(i, cif = False) for i in files]
# 		coords = np.vstack([i[0] for i in all_feats])

# 	mass = all_feats[0][1]
# 	radius = all_feats[0][2]
# 	res_ids = all_feats[0][3]
# 	del all_feats
# 	return coords, mass, radius, res_ids

# # parse_all_ihm_cif(  "../example/PDBDEV/PDBDEV_00000001/", _type = "ihm" )


# def ihm_set_bfactor( modelFile ):
# 	with open( modelFile ) as fh: # "../example/PDBDEV/PDBDEV_00000001/PDBDEV_00000001.cif"
# 		# s is an object of class ihm.System.
# 	    s, = ihm.reader.read(fh)
# 	_id = []
# 	_id = [atom.seq_id for model in s._all_starting_models() for atom in model.get_atoms() if atom.atom_id == "CA"]

# 	i = 0
# 	for model in s._all_starting_models():
# 		for atom in model.get_atoms():
# 			if _id[i] == atom.seq_id:
# 				atom.biso = b
# 			elif _id[i] != atom.seq_id:
# 				i+=1
# 				b +=1
# 				atom.biso = b

# 	with open("modified.cif", "w") as fh:
# 		ihm.dumper.write( fh, [s] )

# def get_chain_length_ihm( s ):
# 	# Obtain the no. of residues in the system.
# 	# s --> an object of class ihm.System.
# 	res_no = 0
# 	for model in s._all_starting_models():
# 		for atom in model.get_atoms():
# 			if atom.atom_id == "CA":
# 				res_no += 1
# 	return res_no
