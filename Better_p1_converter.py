import sys
import numpy as np
import os
import glob
import re
from progress.bar import IncrementalBar

##########################################################
## Program written by Daniel Hewitt and Martin Hutereau ##
######## Any questions feel free to email us at ##########
######## zccadhe@ucl.ac.uk or uccamhu@ucl.ac.uk ##########
######## Members of the Ben Slater group at UCL ##########
##########################################################

class atomic_data:
	def __init__(self,index,image,element,coordinates):
		self.i=index
		self.img=image
		self.elem=element
		self.coord=np.asarray(coordinates,dtype=np.float64)


#Dictionary for converting element labels to atomic numbers
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


#Dictionary for converting plus or minus coordinates into a matrix transformation
sym_dict={'x':np.array([1,0,0],np.float64),
'-x':np.array([-1,0,0],np.float64),
'y':np.array([0,1,0],np.float64),
'-y':np.array([0,-1,0],np.float64),
'z':np.array([0,0,1],np.float64),
'-z':np.array([0,0,-1],np.float64),
'X':np.array([1,0,0],np.float64),
'-X':np.array([-1,0,0],np.float64),
'Y':np.array([0,1,0],np.float64),
'-Y':np.array([0,-1,0],np.float64),
'Z':np.array([0,0,1],np.float64),
'-Z':np.array([0,0,-1],np.float64),
}


#Read all cif files from current dir
files=glob.glob('*.cif')

#Reads in cell params, atom coords, and symmetry operations from cif file
#Atom coordinations stored in atomic data class
def reader(file):
	#Some cif files put errors on their cell params, so remove brackets to avoid this.
	flag=re.compile('[(][0-9]{1,3}[)]')

	cell_params=[]
	site_label=[]
	atom_group=False
	atom_index=0
	structure_data=[]
	symmy = False
	symmetry_ops=[]


	with open(file, 'r') as file:
		for line in file:


			if len(line.split())==0 or 'loop_' in line:
				symmy = False

			if symmy == True:
				symmetry_ops.append(line.strip("',\n\"").split(','))
			if '_space_group_symop_operation_xyz' in line:
				symmy = True

			if line.find('_cell_length')==0 or line.find('_cell_angle_')==0:
				if len(cell_params)<6:
					cell_params.append(float(flag.split(line.split()[1])[0]))


			#Append site lables to allow correct looping when grabbing coords. 
			#Some cif files have different length lines for e.g. charge or spin on atoms.
			if '_atom_site_' in line:
				site_label.append(line.split()[0])
				atom_group=True



			if atom_group==True:
				if len(line.split())>3:
					atom_index+=1
					atom_label=line.split()[site_label.index('_atom_site_type_symbol')]
					x_coord=line.split()[site_label.index('_atom_site_fract_x')]
					y_coord=line.split()[site_label.index('_atom_site_fract_y')]
					z_coord=line.split()[site_label.index('_atom_site_fract_z')]

					structure_data.append(atomic_data(
													atom_index,
													False,
													atom_label,
													[x_coord, y_coord, z_coord]
													)
											)
									
				elif atom_index>=1 and len(line.split())<3:
					atom_group=False

	#Convert all coordinates to be between 0 and 1 in fractional space, to avoid double counting of atom positions
	for tt in range(0,len(structure_data)):
		for s in range(3):
			while structure_data[tt].coord[s]<0.0:
				structure_data[tt].coord[s]=structure_data[tt].coord[s]+1.0
			while structure_data[tt].coord[0]<0.0:
				structure_data[tt].coord[s]=structure_data[tt].coord[s]-1.0        


	#print(symmetry_ops)
	return structure_data, cell_params, symmetry_ops




#Generate all the symmetry operations from the cif file in a usable form.
def Operations(symmetry_ops):

	numbers=[0,1,2,3,4,5,6,7,8,9]

	flag = re.compile('[-]{0,1}[xyzXYZ]{1,1}')
	numflag = re.compile('[-]{0,1}[0-9]{1,1}/{0,1}[0-9]{0,1}')
	all_sym_ops=[]

	for line in symmetry_ops:

		operations=[0,0,0]
		operations2=[[0,0,0],[0,0,0],[0,0,0]]


		for i in range(0,3):

			new_strings=flag.findall(line[i])
			#print(new_strings)
			num_strings=numflag.findall(line[i])
			#print(num_strings)

			for j in new_strings:

				operations2[i]+=sym_dict[j]

			if len(num_strings)==1 and len(num_strings[0])==3:

				changer=float(num_strings[0][0]) / float(num_strings[0][2])
				operations[i]=changer

		all_sym_ops.append([operations, operations2])
		#print(all_sym_ops)

	return(all_sym_ops)


#Apply all symmetry operations to all atom coordinates.
def AtomLooper(atom_coord, all_sym_ops):

	all_coords=[]

	for coord in atom_coord:
		for i in range(len(all_sym_ops)):
			#print(coord.coord)
			new_coord = [a_number_dict[coord.elem], all_sym_ops[i][1]@coord.coord[0:3]]

			#print(new_coord)

			new_coord[1] = new_coord[1] + all_sym_ops[i][0]
			new_coord = [new_coord[0], new_coord[1][0], new_coord[1][1], new_coord[1][2]]


			all_coords.append(new_coord)
			#print(all_coords)

	#Again, change all coordinates to be withint he range 0 and 1.
	for i in range(len(all_coords)):
		for j in range(1, 4):
			if all_coords[i][j] < 0:
				all_coords[i][j]=all_coords[i][j] + 1
			if all_coords[i][j] > 1:
				all_coords[i][j]=all_coords[i][j] - 1





	return(all_coords)


#Metric tensor allows change between cartesian and fractional space.
#Here it is used to calculate cartesian distances from fractional coordinates.
def MetricTensor(cell_params):
	
	a=float(cell_params[0])
	b=float(cell_params[1])
	c=float(cell_params[2])
	alpha=np.radians(float(cell_params[3]))
	beta=np.radians(float(cell_params[4]))
	gamma=np.radians(float(cell_params[5]))

	Metrical=[[a*a, a*b*np.cos(gamma), a*c*np.cos(beta)],
			[b*a*np.cos(gamma), b*b, b*c*np.cos(alpha)],
			[c*a*np.cos(beta), c*b*np.cos(alpha), c*c]]

	return Metrical



#Calculate distances between atoms.
def distance(p,q,Metrical):
	pq = p-q
	distance = np.sqrt(pq@Metrical@np.transpose(pq))
	return distance




#Remove all duplicate atoms within a given tolerance. Here, we use 0.1 angstroms.
def RemoveDuplicates(all_coords, Metrical):

	remove_these=[]
	unique_coords=np.asarray(all_coords, dtype=np.float64)
	unique_coords=np.unique(unique_coords, axis=0)

	for j in range(len(unique_coords)):
		for i in range(j+1, len(unique_coords)):
			length = distance(unique_coords[i][1:4], unique_coords[j][1:4], Metrical)
			#Remove coordinates below the 0.1 angstr cut off
			if length < 0.1:
				remove_these.append(i)

	#print(remove_these)
	#print(unique_coords)


	unique_coords = [unique_coords[g] for g in range(0, len(unique_coords)) if g not in remove_these]



   # print(unique_coords)
	
	return unique_coords


#Make cif file in p1 symmetry
def _makecif(outfile,cell_params,atoms):

	with open(outfile,'w') as file:


		#Pre-position block.
		file.write('\ndata_image0\n'\
		+'_symmetry_space_group_name_H-M    "P 1"\n'\
		+'_symmetry_int_tables_number       1\n'\
		+'\n'\
		+'loop_\n'\
		+' _symmetry_equiv_pos_as_xyz\n'\
		+" 'x, y, z'\n"\
		+'\n')

		#Cell parameter block.
		file.write('_cell_length_a\t{}\n'.format(cell_params[0])\
		+'_cell_length_b\t{}\n'.format(cell_params[1])\
		+'_cell_length_c\t{}\n'.format(cell_params[2])\
		+'_cell_angle_alpha\t{}\n'.format(cell_params[3])\
		+'_cell_angle_beta\t{}\n'.format(cell_params[4])\
		+'_cell_angle_gamma\t{}\n'.format(cell_params[5])\
		+'\n')

		#Position block.
		file.write('loop_\n'\
		+' _atom_site_type_symbol\n'\
		+' _atom_site_fract_x\n'\
		+' _atom_site_fract_y\n'\
		+' _atom_site_fract_z\n')

		for g in atoms:
			file.write(
				' {} {:.4f} {:.4f} {:.4f}\n'.format(
					inv_a_number[g[0]],g[1],g[2],g[3]
					)
				)




def wrapper(filename):
	atom_coord, cell_params, symmetry_ops = reader(filename)

	Metrical = MetricTensor(cell_params)

	all_sym_ops = Operations(symmetry_ops)

	all_coords = AtomLooper(atom_coord, all_sym_ops)
	#print('AtomLooper')

	unique_coords = RemoveDuplicates(all_coords, Metrical)

	#print('RemoveDuplicates')

	outfile = filename.replace('.cif', '_p1.cif')

	_makecif(outfile, cell_params, unique_coords)
	#print('_makecif')

#Progress bar to show how far through the file list you are.
#Useful to see if no progress is being made !
bar = IncrementalBar('Progress', max=len(files), suffix = '%(elapsed)d seconds', width = 130)
#Loop over all files unless structure is in p1 symmetry already
for filename in files:
	#print(filename)
	try: 
		wrapper(filename)
		bar.next()
	except ValueError as e:
		print(e)
		continue

bar.finish()


























