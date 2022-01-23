Hey, I'm Dan and I'm currently a PhD student at UCL working on modelling nanoporous materials for their adsorption properties. I use a mix of classical and quantum simulations in order to model zeolites for their potential use in numerous industrial applications.  
This GitHub is a repository for some of the scripts I've written to automate my work. Check out my Adsorption Profile generator for simple visulation of adsorption behaviour in zeolites!

# Generate Adsorption Profiles

Generate Adsorption Profiles (GAP) is a tool I've made to allow visualisation of adsorbates in zeolite structures using the output pdb files from RASPA. The **OOP2_GAP.py** script available in the **Adsorption_Profiles** directory can be run in the Movies/System_0 directory of any RASPA calculation where guest molecule positions are printed to the Movie files.  

![alt text](https://github.com/d4n-hewitt/d4n-hewitt/blob/main/Adsorption_Profiles/Example_image.png?raw=true)  

An example of the '.pdb' files created from a RASPA output are provided for the three xylene isomers adsorbing in a hypothetical zeolite framework during a competitive adsorption simulation. Provided also is an example of the html file produced by running the script in this directory, which can be viewed in a browser allowing interactive viewing of the adsorption profile in 3D. Buttons have been provided to allow the removal of the Framework or adsorbates in order to more easily view the adsorption behaviour.

# dcTST Setup  
This directory contains an example of how to set up the first stages of a dcTST calculation on a high-throughput scale for many different zeolites. The calculation essentially generates a free energy profile for the adsorbate in the system, however a similar workflow could be used to produce input files for any other type of RASPA simulation on a high-throughput scale.  

# Utility Scripts  
This directory contains a few useful utility scripts I've written to help my work.  
The p1 converter script allows rapid transformation of any CIF file with symmetry transformations present into the p1 symmmetry group - incredibly useful if you're using the Deem or IZA-SC databases of zeolites and wish to see the full coordinates of a unit cell.  
The GCMC extractor tool can be used at the end of a series of RASPA calculations in order to extract useful data such as the loading and errors in the loading, and writes the data into a '.csv' format. This allows the file to be easily read into a Pandas dataframe for analysis.  
The Unique Crystal generator script produces unique aluminosilicate structures given a base silicate structure. This will see use in future publications...



