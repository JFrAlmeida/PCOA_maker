## How to Use PCOA_maker

This script makes PCOAs based on raw counts files where each row is features (COG number, PFAM number, etc.) and columns are 
samples (Counts of each feature per genome)

It can take any number of counts files, as long as the dataset is the same

It requires a file controling the groups you want to represent and colour in the PCOA

The PCOA generated has the legend inside the figure, couldnt figure out an easy way to make it work, I recommend using 
Inkscape or similar to then work the legend however you see fit


### Installation:
Open the terminal, and navigate to where you want your PCOA maker to be in, then type in the temrinal:
```
	git clone https://github.com/JFrAlmeida/PCOA_maker.git
	cd PCOA_maker

	conda env create -f environment.yml
	conda activate pcoa_maker
```
### Settings:
!! Check the settings inside PCOA_Maker.py !!

allows:
- pcoa group filename
- distance metric control
- image format
- image dpi
- option to include names of samples in the PCOA

### What do I prepare for the PCOA

A goups file, following the same format as that in the PCOA_maker/groups/ folder  
Your counts files, place them in PCOA_maker/counts_files

### Usage:
An example counts file and groups file is already in your folders, just run the following once, and check your Outputs:
```
python PCOA_Maker.py
```




