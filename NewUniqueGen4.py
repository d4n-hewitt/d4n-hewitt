from pymatgen.symmetry.analyzer import SpacegroupAnalyzer, PointGroupAnalyzer
from pymatgen.core import Lattice, Structure
from bsym.interface.pymatgen import unique_structure_substitutions
from pymatgen.io.cif import CifParser
import numpy as np
import os
import sys
import math
import glob
import re
sys.path.insert(0,"/Users/dan/Documents/Python_calls/")
from Reader import atomic_data, utility_functions


a_number_dict={'H':1,'He':2,'Li':3,'Be':4,'B':5,'C':6,'N':7,'O':8,'F':9,'Ne':10,
'Na':11,'Mg':12,'Al':13,'Si':14,'P':15,'S':16,'Cl':17,'Ar':18,'K':19,'Ca':20,
'Sc':21,'Ti':22,'V':23,'Cr':24,'Mn':25,'Fe':26,'Co':27,'Ni':28,'Cu':29,'Zn':30,
'Ga':31,'Ge':32,'As':33,'Se':34,'Br':35,'Kr':36,'Rb':37,'Sr':38,'Y':39,'Zr':40,
'Nb':41,'Mo':42,'Tc':43,'Ru':44,'Rh':45,'Pd':46,'Ag':47,'Cd':48,'In':49,'Sn':50,
'Sb':51,'Te':52,'I':53,'Xe':54,'Cs':55,'Ba':56,'La':57,'Ce':58,'Pr':59,'Nd':60,
'Pm':61,'Sm':62,'Eu':63,'Gd':64,'Tb':65,'Dy':66,'Ho':67,'Er':68,'Tm':69,'Yb':70,
'Lu':71,'Hf':72,'Ta':73,'W':74,'Re':75,'Os':76,'Ir':77,'Pt':78,'Au':79,'Hg':80,
'Tl':81,'Pb':82,'Bi':83,'Po':84,'At':85,'Rn':86,'Fr':87,'Ra':88,'Ac':89,'Th':90,
'Pa':91,'U':92,'Np':93,'Pu':94,'Am':95,'Cm':96,'Bk':97,'Cf':98,'Es':99,'Fm':100,
'Md':101,'No':102,'Lr':103,'Rf':104,'Db':105,'Sg':106,'Bh':107,'Hs':108,
'Mt':109,'Ds':110,'Rg':111,'Cn':112,'Nh':113,'Fl':114,'Mc':115,'Lv':116,
'Ts':117,'Og':118}
inv_a_number = dict([[v,k] for k,v in a_number_dict.items()])


file=sys.argv[1]
Topology_code=file.split('.')[0]

Al_num=int(input('Number of Al defects\n'))


def _main(file):

	#####Sort data so that it can be used to generate a pymatgen structures object#####
	structure_data, cell_params=utility_functions.reader(file)
	sym_group = get_sym_group(file)
	labels, coords=pymatgen_formatting(structure_data)
	matrix=Params_to_vector(cell_params)

	##Generate pymatgen structure object##
	mystructure = Structure.from_spacegroup(sym_group, Lattice(matrix), labels, coords, tol=0.1)
	##structure generated##
	Si_num=get_Si_number(mystructure)

	##Generate the symmetry unique structures using bsym#
	unique_structures=unique_structure_substitutions(
		mystructure,'Si',
		{'Al':Al_num,'Si':Si_num-Al_num},
		atol=0.01)

	
	Structure_count=0
	##Reformat the bsym structure object into my own class so that it's easier to use##
	for x in unique_structures:
		Structure_count+=1
		unique_reformatted=reformat_data_to_my_classes(x)
		
		Metrical=utility_functions.MetricTensor(cell_params)
		##Generate the perdiodic replicas of atoms near cell boundaries##
		ghost_atoms=utility_functions.PBC(unique_reformatted, cell_params)


		##Calculate all bond distances##
		distances = utility_functions.distance_compiler(unique_reformatted, Metrical)
		distances.sort(key=lambda x:x[4])



		##Get a list of the indexes of the aluminiums in the current structure##
		Al_list=[]
		for j in range(len(unique_reformatted)):
			if unique_reformatted[j].elem == 13 and unique_reformatted[j].img==False:
				Al_list.append(unique_reformatted[j].i)

		all_sodium_positions=[]
		##Loop through the list of Al indexes, so that we only look at 1 Al at a time##
		for index in Al_list:
			

			##Calculate the four shortest bond distances to Al, identifying the four bonded oxygens##
			Al_O_distances=[]
			for i in range(len(distances)):
				if distances[i][2]==13 and distances[i][4] != 0 and distances[i][3]==8 and distances[i][4] < 2 and distances[i][0] == index:
					Al_O_distances.append(distances[i])
					print(distances[i])
			Al_O_distances.sort(key=lambda x:x[1])
			Al_O_neighbours=Al_O_distances
			#print(len(Al_O_neighbours))

			##Determine which Si is bonded to the O bonded to the Al##
			Si_O_distances=[]
			for i in range(len(distances)):
				for ii in range(len(Al_O_neighbours)):
					if distances[i][1]==Al_O_neighbours[ii][1] and distances[i][4] < 2 and distances[i][4] != 0 and distances [i][0] != Al_O_neighbours[ii][0]:
						Si_O_distances.append(distances[i])
						print(distances[i])
			Si_O_distances.sort(key=lambda x:x[1])
			Si_O_neighbours=Si_O_distances

			#print(Si_O_neighbours)

			current_sodium_positions=generate_sodium_positions(Al_O_neighbours, Si_O_neighbours, unique_reformatted, Metrical, cell_params)
			#print(all_sodium_positions)

			all_sodium_positions.append(current_sodium_positions)
		#print(all_sodium_positions)
		##Calculate the length of the longest list of Na positions. Max should normally be 9 per T site, assuming one 180 degree angle##
		longest_list=longest_list_of_sodium_posisitions(all_sodium_positions, Al_num)

		#print(longest_list)

		sodim_coord_index_combinations=[[]]
		for i in range(Al_num):
			sodim_coord_index_combinations=appender(sodim_coord_index_combinations, longest_list)
		print(sodim_coord_index_combinations)

		length_of_combinations=[]
		for i in range(len(all_sodium_positions)):
			length_of_combinations.append(len(all_sodium_positions[i]))

		garbage=[]
		for i in range(len(length_of_combinations)):
			for element in sodim_coord_index_combinations:
				if element[i]>=length_of_combinations[i]:
					garbage.append(element)

		sodim_coord_index_combinations=[i for i in sodim_coord_index_combinations if i not in garbage]
		print(sodim_coord_index_combinations)



		na_number=0
		print(all_sodium_positions)
		
		for i in sodim_coord_index_combinations:
			sodium_coordinate_data=[]
			for j in range(len(i)):
				print(all_sodium_positions[j][i[j]])
				sodium_coordinate_data.append(atomic_data(
											0,
											False,
											11,
											[all_sodium_positions[j][i[j]][0], all_sodium_positions[j][i[j]][1], all_sodium_positions[j][i[j]][2]]
											)
									)

			na_number+=1
			_makecif(cell_params, unique_reformatted, na_number, Structure_count, sodium_coordinate_data, Topology_code)
			write_inp(cell_params,unique_reformatted, na_number, Structure_count, sodium_coordinate_data, Topology_code)
			del sodium_coordinate_data
		






		print('Next structure...')

	return unique_structures
	












def write_inp(cell_params,unique_reformatted, na_number, Structure_count, sodium_coordinate_data, Topology_code):


	outfile=str(Topology_code)+str(Structure_count)+'na'+str(na_number)+'.inp'
	print('writing '+str(outfile))
	with open(outfile,'w') as file:
		file.write(' &GLOBAL\n')
		file.write('	PRINT_LEVEL LOW\n')
		file.write('PROJECT {}{}na{}\n'.format(Topology_code,Structure_count, na_number))
		file.write('  RUN_TYPE CELL_OPT\n')
		file.write(' &END GLOBAL\n')

		file.write(' &MOTION\n')
		file.write('  &GEO_OPT\n')
		file.write('   TYPE MINIMIZATION\n')
		file.write('   OPTIMIZER LBFGS\n')
		file.write('   MAX_ITER 3000\n')
		file.write('   MAX_DR 2.9999999999999997E-04\n')
		file.write('   MAX_FORCE 4.5000000000000003E-05\n')
		file.write('   RMS_DR 1.4999999999999999E-04\n')
		file.write('   RMS_FORCE 3.0000000000000001E-05\n')
		file.write('  &END GEO_OPT\n')

		file.write('  &CELL_OPT\n')
		file.write('   OPTIMIZER BFGS\n')
		file.write('   MAX_ITER 1000\n')
		file.write('   MAX_DR 3.0000000000000001E-03\n')	
		file.write('   MAX_FORCE 4.4999999999999999E-04\n')	
		file.write('   RMS_DR 1.5000000000000000E-03\n')			
		file.write('   RMS_FORCE 2.9999999999999997E-04\n')
		file.write('   STEP_START_VAL 1\n')
		file.write('   TYPE DIRECT_CELL_OPT\n')
		file.write('   EXTERNAL_PRESSURE 1.0000000000000000E+00\n')
		file.write('   PRESSURE_TOLERANCE 1.0000000000000000E+01\n')
		file.write('  &END CELL_OPT\n')
		file.write(' &END MOTION\n')
					
		file.write(' &FORCE_EVAL\n')
		file.write('  METHOD QS\n')
		file.write('  STRESS_TENSOR ANALYTICAL\n')
		file.write('  &DFT\n')
		file.write('   BASIS_SET_FILE_NAME ./GTH_BASIS_SETS\n')
		file.write('   POTENTIAL_FILE_NAME ./POTENTIAL\n')
		file.write('   CHARGE 0\n')
		file.write('   &SCF\n')
		file.write('	MAX_SCF 30\n')
		file.write('	EPS_SCF 9.9999999999999995E-08\n')
		file.write('	SCF_GUESS ATOMIC\n')
		file.write('	&OT T\n')
		file.write('	 MINIMIZER BROYDEN\n')
		file.write('	 PRECONDITIONER FULL_KINETIC\n')
		file.write('	&END OT\n')
		file.write('	&OUTER_SCF T\n')
		file.write('	 EPS_SCF 9.9999999999999995E-08\n')
		file.write('	&END OUTER_SCF\n')
		file.write('   &END SCF\n')
		file.write('   &QS\n')
		file.write('	EPS_DEFAULT 9.9999999999999998E-13\n')
		file.write('   &END QS\n')
			
		file.write('   &MGRID\n')
		file.write('   CUTOFF 650\n')
		file.write('   &END MGRID\n')
		file.write('   &XC\n')
		file.write('	DENSITY_CUTOFF 1.0000000000000000E-10\n')
		file.write('	GRADIENT_CUTOFF 1.0000000000000000E-10\n')
		file.write('	TAU_CUTOFF 1.0000000000000000E-10\n')
		file.write('	&XC_FUNCTIONAL NO_SHORTCUT\n')
		file.write('	 &PBE T\n')
		file.write('	  PARAMETRIZATION ORIG\n')
		file.write('	 &END PBE\n')
		file.write('	&END XC_FUNCTIONAL\n')
		file.write('   &END XC\n')
		file.write('   &POISSON\n')
		file.write('	POISSON_SOLVER PERIODIC\n')
		file.write('	PERIODIC XYZ\n')
		file.write('   &END POISSON\n')
		file.write('  &END DFT\n')

		file.write('  &SUBSYS\n')
		file.write('  &CELL\n')
		file.write('MULTIPLE_UNIT_CELL 1 1 1\n')
		file.write('ABC {} {} {}\n'.format(cell_params[0],cell_params[1],cell_params[2]))
		file.write('ALPHA_BETA_GAMMA {} {} {}\n'.format(cell_params[3],cell_params[4],cell_params[5]))
		file.write('&END CELL\n')
		file.write('&COORD\n')
		file.write(' SCALED\n')
		for i in range(len(unique_reformatted)):
			if unique_reformatted[i].img==False:
				file.write('{} {} {} {}\n'.format(inv_a_number[unique_reformatted[i].elem], unique_reformatted[i].coord[0], unique_reformatted[i].coord[1], unique_reformatted[i].coord[2]))
		for i in range(len(sodium_coordinate_data)):
			file.write('{} {} {} {}\n'.format(inv_a_number[sodium_coordinate_data[i].elem], sodium_coordinate_data[i].coord[0], sodium_coordinate_data[i].coord[1], sodium_coordinate_data[i].coord[2]))
		file.write('&END COORD\n')
		file.write('	&KIND Si\n')
		file.write('	  BASIS_SET TZV2P-GTH\n')
		file.write('	  POTENTIAL GTH-PBE-q4\n')
		file.write('	&END KIND\n')
		file.write('	&KIND O\n')
		file.write('	  BASIS_SET TZV2P-GTH\n')
		file.write('	  POTENTIAL GTH-PBE-q6\n')
		file.write('	&END KIND\n')
		file.write('	&KIND Al\n')
		file.write('	  BASIS_SET TZV2P-GTH\n')
		file.write('	  POTENTIAL GTH-PBE-q3\n')
		file.write('	&END KIND\n')
		file.write('	&KIND H\n')
		file.write('	  BASIS_SET TZV2P-GTH\n')
		file.write('	  POTENTIAL GTH-PBE-q1\n')
		file.write('	&END KIND\n')
		file.write('	&KIND Na\n')
		file.write('	  BASIS_SET TZV2P-GTH\n')
		file.write('	  POTENTIAL GTH-PBE-q9\n')
		file.write('	&END KIND\n')
	
		file.write('	 &TOPOLOGY\n')
		file.write('	  MULTIPLE_UNIT_CELL 1 1 1\n')
		file.write('	 &END TOPOLOGY\n')
		file.write('   &END SUBSYS\n')
		file.write(' &END FORCE_EVAL\n')





def _makecif(cell_params,unique_reformatted, na_number, Structure_count, sodium_coordinate_data, Topology_code):


		outfile=str(Topology_code)+str(Structure_count)+'na'+str(na_number)+'.cif'
		with open(outfile,'w') as file:

			#Pre-position block.
			file.write('data_image0\n')
			file.write('_symmetry_space_group_name_H-M    "P 1"\n')
			file.write('_symmetry_int_tables_number       1\n')
			file.write('\n')
			file.write('loop_\n')
			file.write(' _symmetry_equiv_pos_as_xyz\n')
			file.write(" 'x, y, z'\n")
			file.write('\n')

			#Cell parameter block.
			file.write('_cell_length_a\t{}\n'.format(cell_params[0]))
			file.write('_cell_length_b\t{}\n'.format(cell_params[1]))
			file.write('_cell_length_c\t{}\n'.format(cell_params[2]))
			file.write('_cell_angle_alpha\t{}\n'.format(cell_params[3]))
			file.write('_cell_angle_beta\t{}\n'.format(cell_params[4]))
			file.write('_cell_angle_gamma\t{}\n'.format(cell_params[5]))
			file.write('\n')

			#Position block.
			file.write('loop_\n')
			#file.write(' _atom_site_label\n')
			file.write(' _atom_site_type_symbol\n')
			file.write(' _atom_site_fract_x\n')
			file.write(' _atom_site_fract_y\n')
			file.write(' _atom_site_fract_z\n')
		

			for i in range(len(unique_reformatted)):
				if unique_reformatted[i].img==False:
					file.write('{} {} {} {}\n'.format(inv_a_number[unique_reformatted[i].elem], unique_reformatted[i].coord[0], unique_reformatted[i].coord[1], unique_reformatted[i].coord[2]))
			for i in range(len(sodium_coordinate_data)):
				file.write('{} {} {} {}\n'.format(inv_a_number[sodium_coordinate_data[i].elem], sodium_coordinate_data[i].coord[0], sodium_coordinate_data[i].coord[1], sodium_coordinate_data[i].coord[2]))


def longest_list_of_sodium_posisitions(all_sodium_positions, Al_num):
	longest_list = 0
	for i in range(Al_num):
		if len(all_sodium_positions[i]) > longest_list:
			longest_list = len(all_sodium_positions[i])
	return longest_list


def appender(vectors, largest):
   newvectors=[]
   for j in range(len(vectors)):
	   for i in range(largest):
		   newvector=vectors[j]+[i]
		   newvectors.append(newvector)
   return newvectors



def generate_sodium_positions(Al_O_neighbours, Si_O_neighbours, unique_reformatted, Metrical, cell_params):

	oxygen_data=[]
	silicon_data=[]
	aluminium_data=[]
	for i in range(len(Al_O_neighbours)):
		oxygen_data.append(unique_reformatted[Al_O_neighbours[i][1]])
		silicon_data.append(unique_reformatted[Si_O_neighbours[i][0]])
		aluminium_data.append(unique_reformatted[Al_O_neighbours[i][0]])
		#print(unique_reformatted[Al_O_neighbours[i][0]].i)
		#print(unique_reformatted[Si_O_neighbours[i][0]].i)
		#print(unique_reformatted[Al_O_neighbours[i][1]].i)

	all_sodium_positions=[]
	#count=0
	for i in range(len(oxygen_data)):
		#print(oxygen_data[i].coord)
		theta=angle(silicon_data[i].coord, oxygen_data[i].coord, aluminium_data[i].coord, Metrical)
		if 180 - theta < 0.1:
			print('this one')
			#all_sodium_positions.append([0, 0, 0])

			Frac2Cartarray = Fract2Cart(cell_params)

			a_cart = np.dot(Frac2Cartarray, silicon_data[i].coord)
			b_cart = np.dot(Frac2Cartarray, aluminium_data[i].coord)


			mean_vec = a_cart - b_cart
			mean_vec=np.array(mean_vec)

			mean_vec = list(map(lambda x: x, mean_vec))
			print(mean_vec)

			vec_alpha=np.array(mean_vec)
			vec_alpha=vec_alpha/np.linalg.norm(vec_alpha)
			print(vec_alpha)


			contains_zero=False
			for ii in range(len(vec_alpha)):
				if vec_alpha[ii] == 0:
					index_to_change=ii
					print(ii)
					contains_zero=True
				#else:
					#index_to_change=100

			if contains_zero==False:
				vec_beta=[(vec_alpha[1]-vec_alpha[2])/vec_alpha[0], (vec_alpha[2]-vec_alpha[0])/vec_alpha[1], (vec_alpha[0]-vec_alpha[1])/vec_alpha[2]]
			elif contains_zero==True and index_to_change==0:
				vec_beta=[0, (vec_alpha[1]-vec_alpha[2])/vec_alpha[1], (vec_alpha[2]-vec_alpha[1])/vec_alpha[2]]
			elif contains_zero==True and index_to_change==1:
				vec_beta=[(vec_alpha[0]-vec_alpha[2])/vec_alpha[0], 0, (vec_alpha[2]-vec_alpha[0])/vec_alpha[2]]
			elif contains_zero==True and index_to_change==2:
				vec_beta=[(vec_alpha[0]-vec_alpha[1])/vec_alpha[0], (vec_alpha[1]-vec_alpha[0])/vec_alpha[1], 0]


			
			vec_beta=vec_beta/np.linalg.norm(vec_beta)
			print(vec_beta)
			vec_gamma=np.cross(vec_alpha,vec_beta)

			print(np.dot(vec_alpha, vec_beta))

			new_axes=[vec_beta,vec_gamma,vec_alpha]

			sod_pos=[]
			for ii in range(0, 6):
				r = 2.5
				theta = 2 * np.pi / 6 * ii

				x_coord = r * np.cos(theta)
				y_coord = r * np.sin(theta)
				sod_coord=[x_coord, y_coord, 0]

				oriented_pos = np.transpose(new_axes)@sod_coord


				Cart2Frac = np.linalg.inv(Frac2Cartarray)
				frac_ori_pos = Cart2Frac@oriented_pos

				
				all_sodium_positions.append(oxygen_data[i].coord + frac_ori_pos)

				#all_sodium_positions.append(sod_pos)

		else:
			#count+=1
			#print(count)
			Frac2Cartarray = Fract2Cart(cell_params)

			SiO_bond=oxygen_data[i].coord - silicon_data[i].coord
			AlO_bond=oxygen_data[i].coord - aluminium_data[i].coord

			SiO_bond_cart=Frac2Cartarray@SiO_bond
			AlO_bond_cart=Frac2Cartarray@AlO_bond

			print('{} {}'.format(SiO_bond_cart, AlO_bond_cart))

			SiO_mod=np.linalg.norm(SiO_bond_cart)
			AlO_mod=np.linalg.norm(AlO_bond_cart)
			
			print('{} {}'.format(SiO_mod, AlO_mod))

			SiO_normalised=SiO_bond/SiO_mod
			AlO_normalised=AlO_bond/AlO_mod

			new_bond=[SiO_normalised + AlO_normalised]
			#new_bond.append(SiO_bond*AlO_mod + AlO_bond*SiO_mod)
			new_bond=np.array(new_bond[0])
			#print(new_bond)
			magnitude=new_bond@Metrical@new_bond
			#print(magnitude)
			scaler=2.5/(np.sqrt(magnitude))
			sodium_position=(scaler*new_bond)+oxygen_data[i].coord
			#print(sodium_position)
			all_sodium_positions.append(sodium_position)

	return all_sodium_positions



def angle(a,b,c,Metric_Tensor): # angle between two vector calculator
	# b is midpoint and angle calc between ba and bc
	ba=a-b
	bc=c-b
	mod_ba=np.sqrt(ba@Metric_Tensor@ba)
	mod_bc=np.sqrt(bc@Metric_Tensor@bc)
	theta=np.arccos((ba@Metric_Tensor@bc)/(mod_ba*mod_bc))
	theta=theta*(180/math.pi)
	if theta < 90:
		theta=180-theta
	return theta


def reformat_data_to_my_classes(x):
	unique_reformatted=[]
	for j in range(len(x.species)):
			unique_reformatted.append(atomic_data(
												j,
												False,
												a_number_dict[str(x.species[j])],
												[x.frac_coords[j][0], x.frac_coords[j][1], x.frac_coords[j][2]]
												)
										)
	return unique_reformatted

def get_Si_number(mystructure):
	Si_num=0
	for i in mystructure.species:
		if str(i) == 'Si':
			Si_num+=1
	return Si_num


def pymatgen_formatting(structure_data):
	labels=[]
	coords=[]
	for i in range(len(structure_data)):
		labels.append(a_number_dict[structure_data[i].elem])
		coords.append(structure_data[i].coord)
	return labels, coords


def get_sym_group(file):
	with open(file) as f:
		for line in f:
			if line.find('_symmetry_space_group_name_H-M')==0:
				sym_group=line.split("'")[1]
	return sym_group


def Fract2Cart(cell_params):
	a = float(cell_params[0])
	b = float(cell_params[1])
	c = float(cell_params[2])
	alpha = np.radians(float(cell_params[3]))
	beta = np.radians(float(cell_params[4]))
	gamma = np.radians(float(cell_params[5]))
	omega = a*b*c*np.sqrt(1-np.power(np.cos(alpha),2)-np.power(np.cos(beta),2)
	-np.power(np.cos(gamma),2)+
	2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))


	Frac2Cartarray = np.asarray([[a, b*np.cos(gamma), c*np.cos(beta)],
					[0, b*np.sin(gamma), c*((np.cos(alpha)-(np.cos(beta)*np.cos(gamma)))/np.sin(gamma))],
					[0, 0, omega/(a*b*np.sin(gamma))]])

	Cart2Frac = np.linalg.inv(Frac2Cartarray)

	return Frac2Cartarray

def Params_to_vector(cell_params):
	Frac2Cartarray = Fract2Cart(cell_params)
	matrix = np.transpose(Frac2Cartarray)
	return matrix





labels = _main(file)

#print(labels)









