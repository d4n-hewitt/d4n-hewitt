import glob
import numpy as np
import multiprocessing as mp

class RASPA_data:
	def __init__(self, framework, guest_molecules, L_mol_kg, L_mol_kg_err, E_Ads, L_mil_g, L_mil_g_err, time, finished):
		self.framework_name=framework
		self.guests=guest_molecules
		self.L_mol_kg=np.asarray(L_mol_kg, dtype=np.float64)
		self.L_mol_kg_err=np.asarray(L_mol_kg_err, dtype=np.float64)
		self.E_Ads=np.asarray(E_Ads, dtype=np.float64)
		self.L_mil_g=np.asarray(L_mil_g, dtype=np.float64)
		self.L_mil_g_err=np.asarray(L_mil_g_err, dtype=np.float64)
		self.time=time
		self.fin=finished

def Get_data(filename):

	framework=filename.split('/')[1]
	with open(filename, 'r') as file:



		Molecule_Definitions_region=False
		Loading_region=False
		Time_block=False
		guest_molecules=[]
		L_mol_kg=[]
		L_mol_kg_error=[]
		L_mil_g=[]
		L_mil_g_error=[]
		finished=False
		Total_time=0
		Enthalpy_region=False
		Enth_data_region=False
		E_Ads=[]


		for line in file:

			if 'Simulation finished,' in line:
				finished=True

			if 'THE SYSTEM HAS A NET CHARGE' in line or 'INAPPROPRIATE NUMBER OF UNIT CELLS USED' in line:
				finished=False


			if len(line.split())==0:
				continue

##################################################################################

			if 'MoleculeDefinitions' in line:
				Molecule_Definitions_region=True

			if Molecule_Definitions_region==True:
				if '(Adsorbate molecule)' in line:
					guest_molecules.append(line.split()[2].strip('[').strip(']'))

			if 'Framework Status' in line:
				Molecule_Definitions_region=False

##################################################################################

			if 'Number of molecules:' in line:
				Loading_region=True

			if Loading_region==True:
				if 'Average loading absolute [mol/kg framework]' in line:
					L_mol_kg.append(line.split()[5])
					L_mol_kg_error.append(line.split()[7])

				if 'Average loading absolute [milligram/gram framework]' in line:
					L_mil_g.append(line.split()[5])
					L_mil_g_error.append(line.split()[7])


##################################################################################


			if 'Total CPU timings:' in line:
				Time_block=True
			if Time_block==True:
				if 'total time:' in line:
					Total_time=line.split()[2]

			if 'Production run CPU timings of the MC moves:' in line:
				Time_block=False

##################################################################################


			if 'Enthalpy of adsorption:' in line:
				Enthalpy_region=True
			if Enthalpy_region==True:
				if 'Component' in line:
					Enth_data_region=True

			if Enthalpy_region and Enth_data_region == True:
				#print(line)
				if len(line.split()) ==4 and 'Block[ ' not in line and 'Total' not in line:
					#print(line)
					E_Ads.append(line.split()[0])
					#print(E_Ads)
			if 'Total enthalpy of adsorption from components and measured mol-fraction' in line:
				Enth_data_region = False
			if 'derivative of the chemical potential with respect to density (constant T,V):' in line:
				Enthalpy_region = False



##################################################################################


		structure_data=RASPA_data(
									framework,
									guest_molecules,
									L_mol_kg,
									L_mol_kg_error,
									E_Ads,
									L_mil_g,
									L_mil_g_error,
									Total_time,
									finished)

		return structure_data


def write_headers(structure_data):

	with open('Data_enth.txt', 'w') as file:

		#Write out the headers for the csv

		file.write('framework_name,')
		for i in range(len(structure_data.guests)):
			file.write('{}_L_mol_kg,'.format(structure_data.guests[i]))
			file.write('{}_L_mol_kg_err,'.format(structure_data.guests[i]))
			
		for i in range(len(structure_data.guests)):
			file.write('{}_E_Ads,'.format(structure_data.guests[i]))

		for i in range(len(structure_data.guests)):
			file.write('{}_L_mil_g,'.format(structure_data.guests[i]))
			file.write('{}_L_mil_g_err,'.format(structure_data.guests[i]))
			
		file.write('CPU_time\n')



def write_csv(structure_data):

	with open('Data_enth.txt', 'a') as file:

		file.write('{},'.format(structure_data.framework_name))

		for i in range(len(structure_data.guests)):
			file.write('{},'.format(structure_data.L_mol_kg[i]))
			file.write('{},'.format(structure_data.L_mol_kg_err[i]))

		for i in range(len(structure_data.guests)):
			file.write('{},'.format(structure_data.E_Ads[i]))
			

		for i in range(len(structure_data.guests)):
			file.write('{},'.format(structure_data.L_mil_g[i]))
			file.write('{},'.format(structure_data.L_mil_g_err[i]))
			
		file.write('{}\n'.format(structure_data.time))




def wrapper(filename):
		structure_data = Get_data(filename)
		if structure_data.fin==True:
			write_csv(structure_data)
	

def _multiprocess(processes,filenames):
	pool=mp.Pool(processes=processes)
	results=[pool.apply_async(wrapper,args=(ww,)) for ww in filenames]
	results=[p.get() for p in results]
	return results

if __name__=="__main__":
	files=glob.glob('*/*/Output/System_0/*.data')
	structure_data=Get_data(files[0])
	write_headers(structure_data)
	filenames = files
	_multiprocess(20, filenames)




