import numpy as np
# from math import pi
import os
import wget
import requests
import json
from collections import defaultdict

import ihm
import ihm.reader
import ihm.dumper

import IMP
# import IMP.atom
import RMF
from dcd_parser import *
from main import *
# from pdb_parser import *
from color_precision import colour_rmf
from color_precision_pdb import scale


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
			exit()
		else:
			found_dcd = True
	return found_dcd


def download_dcd( url, fl ):
	# Download the dcd file using the url.
	# url = "/".join(e.file.repo.url.split("/")[:-1])
	print( "Downloading dcd file..." )
	response = wget.download( url, fl )
	print( "\n" )


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
	avg_radius = 2.73

	# iterate over all models in the model group.
	num_models = 0
	for m in models.keys():
		num_models += 1
		num_spheres = 0
		a = 0
		# iterate over all asym_units in a model.
		for asym_unit in models[m].keys():
			# iterate over all residues in an asym_unit.
			for res in models[m][asym_unit].keys():				
				if IMP_:
					obj = models[m][asym_unit][res]
				else:
					obj = models[m][asym_unit][res]["CA"]
				# Check if it is a sphere or atom object.
				if isinstance(obj, ihm.model.Sphere) or isinstance(obj, ihm.model.Atom):
					num_spheres += 1
					coords.append( [obj.x, obj.y, obj.z])
					try:
						n_res = obj.seq_id_range
						bead_mass = 110 if n_res[0] == n_res[1] else ( n_res[1] - n_res[0] + 1 )*110
						mass.append( bead_mass )
						radius.append( obj.radius )
						ps_names.append( f"{m}_{asym_unit}_{n_res[0]}-{n_res[1]}" )
					# non-IMP entries do not have a radius.
					except AttributeError:
						radius.append( avg_radius )
						mass.append( avg_mass )
						ps_names.append( f"{m}_{asym_unit}_{res}" )
						
	coords = np.array( coords ).reshape( num_models, num_spheres, 3 )
	radius = np.array( radius ).reshape( num_models*num_spheres, 1 )
	mass = np.array( mass ).reshape( num_models*num_spheres, 1 )
	return [coords, radius, mass, ps_names]


def get_patch_coloured_rmf( annotations, coords, mass, radius ):
	with open( f"./{annotations}/annotations_cl2.txt",'r') as pf:
		file = pf.readlines()[1:]
	precisions = [pr.split(',')[2:4] for pr in file]
	particles = [pr.split(',')[1] for pr in file]
	
	ind_dict = {}
	for i,prec in enumerate( precisions ):
		ind_dict[particles[i]] = prec

	# Create a new RMF with beads of same XYZR as selected particles, but colored according to precision
	m_new = IMP.Model()

	# Create new hierarchy to add new beads to
	p_root = m_new.add_particle("System")
	h_root = IMP.atom.Hierarchy.setup_particle( m_new, p_root )
	for i,bead_name in enumerate( particles ):
		res = bead_name.split( "_" )[2]

		# create new particle
		# Decorate with arbitrary mass, coordinates, radius.
		p_new = IMP.Particle( m_new, res )
		x, y, z = coords[0][i]
		r = radius[i][0]		
		p_new = IMP.atom.Mass.setup_particle( m_new, p_new, mass[i][0] )
		p_new = IMP.core.XYZR.setup_particle( p_new, IMP.algebra.Sphere3D( IMP.algebra.Vector3D(x, y, z), r ) )

		# Color particle based on its precision
		tp = ind_dict[bead_name]
		c_new, m_new, p_new = colour_rmf( tp, m_new, p_new )
		
		# Add particle to the new model's hierarchy
		h_new = IMP.atom.Hierarchy.setup_particle( m_new, p_new )
		h_root.add_child(h_new)
	# Now create a new RMF file with the new model
	rmf_new = RMF.create_rmf_file( f"./{annotations}/patch_coloured_model.rmf3" )
	IMP.rmf.add_hierarchy( rmf_new, h_root )
	IMP.rmf.save_frame( rmf_new )
	del rmf_new


def ihm_set_bfactor( precision_file, mmcif, model_group ):
	# Set the bfactor as the precision values.
	with open(f"./{precision_file}/bead_spreads_cl2.txt",'r') as pf:
		precision = [float(pr) for pr in pf.readlines()]
	precision = scale( np.array( precision ) )
	precision = [round( p,2 ) for p in precision]

	for istg, stg in enumerate(mmcif.state_groups):
		for ist, st in enumerate(stg):
			for img, mg in enumerate(st):
				if mg == model_group:
					for im, m in enumerate(mg):
						_id = [atom.seq_id for atom in m.get_atoms() if atom.atom_id == "CA"]
						i = 0
						for atom in m.get_atoms():
							if _id[i] == atom.seq_id:
								atom.biso = precision[i]
							elif _id[i] != atom.seq_id:
								i = i+1 if i != len( precision )-1 else i
								atom.biso = precision[i]
								
								
	with open( f"./{precision_file}/patch_coloured_model.cif", "w" ) as fh:
		ihm.dumper.write( fh, [mmcif] )


def parse_ihm_models( args ):
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
					
					# Run Prism for IMP entries only if dcd file is present.
					# Can be changed later.
					if found_dcd:	
						for i, e in enumerate(mmcif.ensembles):
							if e.model_group == mg:
								model_group = parse_models( mmcif, mg )
								
								url = e.file.repo.url
								fl = f"Ensemble{img}_{i}.dcd"
								if not os.path.isfile( fl ):
									download_dcd( url, fl )
								coords = DCDReader( fl ).return_coords()
								_, radius, mass, ps_names = get_all_attributes( model_group, IMP_ )

								if coords.shape[1] != radius.shape[0]:
									raise ValueError('Number of particles in cif file does not match \
										number of particles in dcd file {}'.format(coords.shape[0], radius.shape[0]))
									exit()
								# Run PrISM to get bead precision.
								run_prism( coords, mass, radius, ps_names, args, f"output_{ist}_{img}" )
								# Get a patch coloured bead model.
								print( "Creating patch coloured model..." )
								get_patch_coloured_rmf( f"output_{ist}_{img}", coords, mass, radius )
					else:
						print( "IMP used for modelling. \n No dcd file found." )
						exit()
	else:
		for istg, stg in enumerate( mmcif.state_groups ):
			# iterate over states in state groups
			for ist, st in enumerate( stg ):
				# iterate over model groups in states
				for img, mg in enumerate( st ):
					model_groups = parse_models( mmcif, mg )
					coords, radius, mass, ps_names = get_all_attributes( model_groups, IMP_ )
					# Run PrISM to get bead precision.
					run_prism( coords, mass, radius, ps_names, args, f"output_{ist}_{img}" )
					# Get a patch coloured ihm cif file.
					print( "Creating patch coloured model..." )
					# get_patch_coloured_rmf( f"output_{ist}_{img}", coords, mass, radius )
					ihm_set_bfactor( f"output_{ist}_{img}", mmcif, mg )
