import re
import numpy as np
import sys
import glob
import pickle
import pandas as pd
import numpy as np
import plotly.graph_objects as go



dic={'ortho': None, 'meta': None, 'para': None}


class atomic_data:
	def __init__(self,index,image,element,coordinates):
		self.i=index
		self.img=image
		self.elem=element
		self.coord=coordinates

framework_file='Framework_final.pdb'
files=glob.glob('Movie_*_component_*.pdb')
precision=10
min_probability=3


def reader(file, precision):

	isomer=file.split('_')[7]
	structure=file.split('_')[1]
	with open(file, 'r') as file:

		atom_index=0
		structure_data=[]

		for line in file:

			if 'CRYST1' in line:
				cell_params=[int( float(line.split()[1]) * precision) + 1, 
							int( float(line.split()[2]) * precision) + 1, 
							int( float(line.split()[3]) * precision) + 1, 
							line.split()[4], 
							line.split()[5], 
							line.split()[6]]

			if 'ATOM' in line:
				atom_index+=1
				atom_label=line.split()[2]
				x_coord=int( float(line.split()[4]) * precision)
				y_coord=int( float(line.split()[5]) * precision)
				z_coord=int( float(line.split()[6]) * precision)

				structure_data.append(atomic_data(
												atom_index,
												False,
												atom_label,
												[x_coord, y_coord, z_coord]
												)
										)

	return cell_params, structure_data, isomer, structure



def Generate_grid(cell_params):
	grid=np.zeros((cell_params[0], cell_params[1], cell_params[2]))
	return grid

def Populate_grid(structure_data, grid, min_probability):

	for atom in structure_data:
		try:
			grid[atom.coord[0], atom.coord[1], atom.coord[2]] += 1
		except IndexError:
			continue

	max_value=np.max(grid)
	grid[grid < min_probability] = 0

	return grid




#def Convert_grid(grid):
#	coordinates = []

#	coordinates=[index for index in grid if grid[grid > 0]]

#	return coordinates


def Generate_framework(framework_file):

	coords=[]
	with open('Framework_final.pdb') as file:
		for line in file:
		
			if 'CRYST1' in line:
				cell_params=[float(line.split()[1]), 
						float(line.split()[2]), 
						float(line.split()[3])]

		
			if 'ATOM' in line:
				coords.append([line.split()[2], 
					float(line.split()[4]), 
					float(line.split()[5]), 
					float(line.split()[6])])

	cell_params=[int(cell_params[0] *10  ) +1, 
				int(cell_params[1] *10   )+1, 
				int(cell_params[2]  *10  )+1]

	coords= [[x[0], int(x[1] *10   ), 
			int(x[2]  *10  ), 
			int(x[3]  *10  )] for x in coords]

	df=pd.DataFrame(coords, columns=['atom', 'x', 'y', 'z'])
	Odf=df[df['atom']=='O']
	Sidf=df[df['atom']=='Si']

	return Odf, Sidf


def Sort_isomer_data(grid, isomer):
	
	o_coords=np.argwhere(dic['ortho'] > 0)
	m_coords=np.argwhere(dic['meta'] > 0)
	p_coords=np.argwhere(dic['para'] > 0)

	return o_coords, m_coords, p_coords


def write_figure(Odf, Sidf, o_coords, m_coords, p_coords, structure):
	fig=go.Figure()

	fig.layout.scene.camera.projection.type = "orthographic"
	fig.update_layout(legend_itemsizing="constant")

	fig.add_trace(go.Scatter3d(
		x=Odf['x'],
		y=Odf['y'],
		z=Odf['z'],
		name='Oxygen',
		mode='markers',
		marker=dict(
			size=4,
			color='red',
			opacity=0.3))
		)

	fig.add_trace(go.Scatter3d(
		x=Sidf['x'],
		y=Sidf['y'],
		z=Sidf['z'],
		name='Silicon',
		mode='markers',
		marker=dict(
			size=6,
			color='blue',
			opacity=0.3))
		)

	fig.add_trace(go.Scatter3d(
		x=m_coords[:,0],
		y=m_coords[:,1],
		z=m_coords[:,2],
		name='Meta xylene',
		mode='markers',
		marker=dict(
			size=1,
			color='green',
			opacity=0.3))
		)
	fig.add_trace(go.Scatter3d(
		x=p_coords[:,0],
		y=p_coords[:,1],
		z=p_coords[:,2],
		name='Para xylene',
		mode='markers',
		marker=dict(
			size=1,
			color='pink',
			opacity=0.3))
		)
	fig.add_trace(go.Scatter3d(
		x=o_coords[:,0],
		y=o_coords[:,1],
		z=o_coords[:,2],
		name='Ortho xylene',
		mode='markers',
		marker=dict(
			size=1,
			color='gold',
			opacity=0.3))
		)

	fig.update_scenes(xaxis_visible=False, yaxis_visible=False,zaxis_visible=False )

	fig.update_layout(
		updatemenus=[
			dict(
				type="buttons",
				buttons=[
					dict(label="Reset",
						 method="update",
						 args=[{"visible": [True, True, True, True, True]},
							   {"trace": "trace", 
								"trace": "trace", 
								"trace": "trace", 
								"trace": "trace", 
								"trace": "trace"}]),
					dict(label="Remove Framework",
						 method="update",
						 args=[{"visible": [False, False, True, True, True]},
							   {"trace": "trace", 
								"trace": "trace", 
								"trace": "trace", 
								"trace": "trace", 
								"trace": "trace"}]),
					dict(label="Remove Meta",
						 method="update",
						 args=[{"visible": [True, True, False, True, None]},
							   {"trace": "trace", 
								"trace": "trace", 
								"trace": "trace", 
								"trace": "trace", 
								"trace": "trace"}]),
					dict(label="Remove Para",
						 method="update",
						 args=[{"visible": [True, True, True, False, True]},
							   {"trace": "trace", 
								"trace": "trace", 
								"trace": "trace", 
								"trace": "trace", 
								"trace": "trace"}]),
					dict(label="Remove Ortho",
						 method="update",
						 args=[{"visible": [True, True, True, True, False]},
							   {"trace": "trace", 
								"trace": "trace", 
								"trace": "trace", 
								"trace": "trace", 
								"trace": "trace"}])
				],
			)
		]
	)

	fig.write_html(f'{structure}.html')





def write_grid(file, grid):
	print(file)

	structure_name=file.split('_')[1]
	component=file.split('_')[7]
	with open(f'{structure_name}_{component}_grid.pkl', 'wb') as file:
		pickle.dump(grid, file)



Odf, Sidf = Generate_framework(framework_file)
for file in files:
	cell_params, structure_data, isomer, structure = reader(file, precision)
	grid=Generate_grid(cell_params)
	grid=Populate_grid(structure_data, grid, min_probability)
	dic[isomer]=grid
o_coords, m_coords, p_coords = Sort_isomer_data(grid, isomer)
write_figure(Odf, Sidf, o_coords, m_coords, p_coords, structure)







