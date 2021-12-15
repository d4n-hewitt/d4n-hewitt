import numpy as np
import os
import time
import shutil




def main():

	start=time.perf_counter()

	path, Adsorbates, simulation_files, molecule_files, zeolite_names = Init_lists()
	for zeolite in zeolite_names:
		zeolite_name=zeolite.split('.')[0]
		Make_directory(zeolite_name)
		for ads in Adsorbates:
			current_dir=os.path.join(path, zeolite_name, ads)
			Make_directory(current_dir)

			Move_files(simulation_files, current_dir)
			Move_files(molecule_files, current_dir)
			shutil.copy(os.path.join(path, 'Zeolites', zeolite), current_dir)



			cell_exp = UnitCellExpander(os.path.join(current_dir, zeolite))
			Write_input(current_dir, zeolite_name, cell_exp, ads)
			Write_runner(current_dir, zeolite_name, ads)
			Write_subber(current_dir)

	end=time.perf_counter()
	print(f'os time = {end - start}')




def Init_lists():

	path=os.getcwd()

	Adsorbates=[x.split('_')[0] for x in os.listdir(os.path.join(path,'molecule_files'))]
	#print(Adsorbates)

	simulation_files=[os.path.join(path, 'simulation_files', x) for x in os.listdir(os.path.join(path,'simulation_files'))]
	molecule_files=[os.path.join(path, 'molecule_files', x) for x in os.listdir(os.path.join(path,'molecule_files'))]
	#zeolites=[os.path.join(path, 'Zeolites', x) for x in os.listdir(os.path.join(path,'Zeolites'))]
	zeolite_names=os.listdir(os.path.join(path,'Zeolites'))
	return path, Adsorbates, simulation_files, molecule_files, zeolite_names






def Make_directory(name):
	os.mkdir(name)

def Move_files(files_to_move, current_dir):
	[shutil.copy(i, current_dir) for i in files_to_move]



def Write_input(current_dir, zeolite, supercell_list, adsorbate):
	with open(f'{current_dir}/simulation.input', 'w') as file:
		file.write(f'''SimulationType MonteCarlo
NumberOfCycles 1000000000000
NumberOfInitializationCycles 1000
PrintEvery 10000

Framework 0
FrameworkName {zeolite}
UnitCells {supercell_list[0]} {supercell_list[1]} {supercell_list[2]}
ExternalTemperature 523.0
ExternalPressure 1.5E6
RemoveAtomNumberCodeFromLabel yes

WriteFreeEnergyProfileEvery 5000

ContinueAfterCrash no
WriteBinaryRestartFileEvery 1000

ChargeMethod Ewald
CutOff 14
EwaldPrecision 1e-6
CutOffVDW 14.0

Component 0 MoleculeName              {adsorbate}_xylene
            StartingBead              0
            MoleculeDefinition        TraPPE
            ComputeFreeEnergyProfile  yes
            TranslationProbability    1.0
            RotationProbability       1.0
            ReinsertionProbability    1.0
            CreateNumberOfMolecules   0''')


def Write_runner(current_dir, zeolite_name, ads):
	with open(f'{current_dir}/run_me.sh', 'w') as file:
		file.write(f'''#!/bin/bash -l

# 1. Force bash as the executing shell.
#$ -S /bin/bash

# 2. Request one hour of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=48:00:00

# 3. Request 1 gigabyte of RAM per core
#$ -l mem=1G

# 4. Set the name of the job.
#$ -N J{zeolite_name.split('_')[0]}_{ads}

# 6. Select the MPI parallel environment and 12 processors.
#$ -pe smp 1

#$ -wd {current_dir}
# 7. Set the working directory to somewhere in your scratch space.
# Replace <your userid> with your UCL userid.
~/RASPA_INSTALL/bin/simulate simulation.input''')






def UnitCellExpander(filename):
        cell_params=[]
        cell_block=False
        rad_cut=29
        cell='_cell_'
        with open(filename,'r') as f:
                for line in f:
                        if len(line.split())==2:
                                cell_block=True
                        if cell_block==True:
                                if line.find(cell)==0:
                                        if len(cell_params)<6:
                                                cell_params.append(line.split()[1])

        a=float(cell_params[0])
        b=float(cell_params[1])
        c=float(cell_params[2])
        alpha=np.radians(float(cell_params[3]))
        beta=np.radians(float(cell_params[4]))
        gamma=np.radians(float(cell_params[5]))

        cell_l=cell_params[0]

        by=np.sin(gamma)
        cell_w=b*by

        cy=(np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma)
        cz=np.sqrt(np.sin(beta)**2-cy**2)
        cell_h=c*cz
        ######new stuff#########
        ######new stuff#########
        ######new stuff#########
        x_dim=int(rad_cut/float(cell_l))+(rad_cut%float(cell_l)>0)
        y_dim=int(rad_cut/float(cell_w))+(rad_cut%float(cell_w)>0)
        z_dim=int(rad_cut/float(cell_h))+(rad_cut%float(cell_h)>0)

        cell_exp = [x_dim,y_dim,z_dim]

        return cell_exp





def Write_subber(current_dir):
	with open('qsubber.sh', 'a') as file:
		file.write(f'qsub {current_dir}/run_me.sh\n')






if __name__ == '__main__':
	main()





