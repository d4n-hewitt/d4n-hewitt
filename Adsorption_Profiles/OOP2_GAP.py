import sys
import glob
import pandas as pd
import numpy as np
import plotly.graph_objects as go


framework_file='Framework_final.pdb'
files=glob.glob('Movie_*_component_*.pdb')
precision=10
min_probability=4


class atomic_data:
	def __init__(self,element,coordinates):
		self.elem=element
		self.coord=coordinates

	def to_dict(self):
		return {
			'elem' : self.elem,
			'x' : self.coord[0],
			'y' : self.coord[1],
			'z' : self.coord[2]
		}

class reader:
	def __init__(self):
		self.box_params=None
		self.cell_params=None
		self.isomer=None
		self.zeolite=None
		self.structure_data=[]
		self.original_sd=[]

	def extras(self, file):
		self.isomer=file.split('_')[7]
		self.zeolite=file.split('_')[1]

	def read_pdb(self, file, precision):
		with open(file, 'r') as file:
			self.structure_data=[]
			for line in file:
				if 'CRYST1' in line:
					self.box_params=[int(float(line.split()[x])) * precision + 1 for x in range(1,4)]
					self.cell_params=[float(line.split()[x]) for x in range(1,7)]

				if 'ATOM' in line:
					atom_label=line.split()[2]
					coords=[int(float(line.split()[x]) * precision) for x in range(4,7)]
					orig_coords=[float(line.split()[x]) for x in range(4,7)]

					self.structure_data.append(atomic_data(
													atom_label,
													coords
													)
											)
					self.original_sd.append(atomic_data(
													atom_label,
													orig_coords)
											)


class data_grid:
	def __init__(self, box_params, structure_data, min_probability, isomer):
		self.grid=[]
		self.box_params=box_params
		self.Generate_grid()
		self.Populate_grid(structure_data)
		self.Reduce_grid(min_probability)
		self.isomer=isomer

	def Generate_grid(self):
		self.grid=np.zeros((self.box_params[0], self.box_params[1], self.box_params[2]))

	def Populate_grid(self, structure_data):
		for atom in structure_data:
			try:
				self.grid[atom.coord[0], atom.coord[1], atom.coord[2]] += 1
			except IndexError:
				continue

	def Reduce_grid(self, min_probability):
		self.grid[self.grid < min_probability] = 0
		self.grid[self.grid >= min_probability] = 1
		self.grid.astype(dtype=np.int8)

	def Extract_points(self):
		grid_points = np.argwhere(self.grid == 1)
		grid_points.astype(dtype=np.int16)
		df = pd.DataFrame(grid_points, columns=['x', 'y', 'z'])
		ads = [self.isomer for x in grid_points]
		df['adsorbate'] = ads
		return df


class utilities:
	def __init__(self, cell_params, structure_data):
		self.cell_params=cell_params
		self.structure_data=structure_data

	def Fract2Cart(self):
		a, b, c = self.cell_params[0:3]
		alpha, beta, gamma = self.cell_params[3:6]

		omega = a*b*c*np.sqrt(1-np.power(np.cos(alpha),2)-np.power(np.cos(beta),2)
		-np.power(np.cos(gamma),2)+
		2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))

		Frac2Cartarray = np.asarray([[a, b*np.cos(gamma), c*np.cos(beta)],
						[0, b*np.sin(gamma), c*((np.cos(alpha)-(np.cos(beta)*np.cos(gamma)))/np.sin(gamma))],
						[0, 0, omega/(a*b*np.sin(gamma))]])
		return Frac2Cartarray
	
	def Convert_to_fractional(self):
		Frac2Cartarray = self.Fract2Cart()
		Cart2Frac = np.linalg.inv(Frac2Cartarray)
		frac_coords=[]
		for atom in self.structure_data:
			t_coords=np.dot(Cart2Frac, atom.coord[:])
			frac_coords.append(atomic_data(
										atom.elem,
										t_coords))
		return frac_coords

	def Supersize(self):
		frac_coords = self.Convert_to_fractional()
		super_frac=[]
		for atom in frac_coords:
			for x in range(0,3):
				for y in range(0,3):
					for z in range(0,3):
						super_frac.append(atomic_data(
													atom.elem,
													[atom.coord[0]+x-1,
													atom.coord[1]+y-1,
													atom.coord[2]+z-1]))
		return super_frac

	def Convert_to_cartesian(self):
		super_frac = self.Supersize()
		Frac2Cartarray = self.Fract2Cart()
		super_cart=[]
		for atom in super_frac:
			t_coords=np.dot(Frac2Cartarray, atom.coord[:])
			super_cart.append(atomic_data(
											atom.elem,
											t_coords))
		return super_cart

	def Cubify_cell(self):
		super_cart = self.Convert_to_cartesian()
		cubic_cell=[]
		for atom in super_cart:
			if 0 < atom.coord[0] < self.cell_params[0]:
				if 0 < atom.coord[1] < self.cell_params[1]:
					if 0 < atom.coord[2] < self.cell_params[2]:
						cubic_cell.append(atomic_data(
													atom.elem,
													atom.coord))
		return cubic_cell

	def Cell_to_grid(self):
		cubic_cell = self.Cubify_cell()
		grid_cell=[]
		for atom in cubic_cell:
			grid_cell.append(atomic_data(
										atom.elem,
										[int(atom.coord[0] * 10),
										int(atom.coord[1] * 10),
										int(atom.coord[2] * 10)]))
		return grid_cell

	def Object_to_df(self):
		grid_cell = self.Cell_to_grid()
		df = pd.DataFrame.from_records([s.to_dict() for s in grid_cell])
		df=df.astype({'x':'int16', 'y':'int16', 'z':'int16'})
		return df


class Figure:
	def __init__(self):
		self.fig=go.Figure()
		self.button_dict=[]

	def Update_settings(self):
		self.fig.layout.scene.camera.projection.type = 'orthographic'
		self.fig.update_layout(legend_itemsizing='constant')
		self.fig.update_scenes(xaxis_visible=False, yaxis_visible=False,zaxis_visible=False )
		self.fig.update_layout(
				updatemenus=[
						dict(type = 'buttons',
							buttons = self.button_dict)])

	def Create_basic_buttons(self, Elements):
		True_array=[True for x in range(len(Elements))]
		self.button_dict.append(dict(label='Reset',
								method='update',
								args=[{'visible': True_array}]))
		True_array[0], True_array[1] = False, False
		self.button_dict.append(dict(label='Remove Framework',
								method='update',
								args=[{'visible': True_array}]))

	def Create_ads_buttons(self, Elements):
		True_array=[True for x in range(len(Elements))]
		for idx, element in enumerate(Elements):
			RemEl_array = True_array.copy()
			RemEl_array[idx] = False
			args = [{'visible': RemEl_array}]
			self.button_dict.append(dict(label=f'Remove {element}',
								method='update',
								args=args))

	def Add_Framework(self, data):
		Element=str(data['elem'].iloc[0])
		dic={'O': 'red',
			'Si': 'blue'}
		col=dic[Element]
		self.fig.add_trace(go.Scatter3d(
			x=data['x'],
			y=data['y'],
			z=data['z'],
			name=Element,
			mode='markers',
			marker=dict(
				size=4,
				color=col,
				opacity=0.3)))

	def Add_Adsorbate(self, data, colour):
		Adsorbate = data['adsorbate'].iloc[0]
		self.fig.add_trace(go.Scatter3d(
			x=data['x'],
			y=data['y'],
			z=data['z'],
			name=Adsorbate,
			mode='markers',
			marker=dict(
				size=1,
				color=colour,
				opacity=0.5)))


class Wrapper:
	def __init__(self, framework_file, precision, files, min_probability):
		self.f_dfs=None
		self.ads_dfs=[]
		self.Elements=[]
		self.read=reader()
		self.Generate_framework(framework_file, precision)
		self.Generate_adsorbate(files, precision, min_probability)
		self.Generate_Elements_list()
		

	def Generate_framework(self, framework_file, precision):
		self.read.read_pdb(framework_file, precision)
		utility=utilities(self.read.cell_params, self.read.original_sd)
		framework_df = utility.Object_to_df()
		O_df = framework_df[framework_df['elem']=='O']
		Si_df = framework_df[framework_df['elem']=='Si']
		self.f_dfs = [O_df, Si_df]

	def Generate_adsorbate(self, files, precision, min_probability):
		ads_dfs=[]
		for file in files:
			self.read.extras(file)
			self.read.read_pdb(file, precision)
			Grid = data_grid(self.read.box_params, self.read.structure_data, min_probability, self.read.isomer)
			adsorbate_df = Grid.Extract_points()
			self.ads_dfs.append(adsorbate_df)

	def Generate_Elements_list(self):
		self.Elements.append([x['elem'].iloc[0] for x in self.f_dfs])
		self.Elements.append([x['adsorbate'].iloc[0] for x in self.ads_dfs])
		self.Elements = [item for sublist in self.Elements for item in sublist]

	def Generate_figure(self):
		Fig=Figure()
		for df in self.f_dfs:
			Fig.Add_Framework(df)
		for idx, df in enumerate(self.ads_dfs):
			Adsorbate = df['adsorbate'].iloc[0]
			colour = input(f'Choose colour for representation of {Adsorbate}:\t')
			Fig.Add_Adsorbate(df, colour)

		Fig.Create_basic_buttons(self.Elements)
		Fig.Create_ads_buttons(self.Elements)
		Fig.Update_settings()
		fig=Fig.fig
		zeolite = self.read.zeolite

		fig.write_html(f'{zeolite}.html')

if __name__ == '__main__':
	wrapper=Wrapper(framework_file, precision, files, min_probability)
	wrapper.Generate_figure()







