# 4900


For the CIS 4900 Research Project Course



File Navigation

1GZX BLAST                      Contains the 30 proteins identified by BLAST that was used

Library	                        Contains the object files and class files used for 4900.py

PDB Files	                Six proteins used to test and create the dihedral angle script

4900.py			        The main program

Magnitudes.xlsx		        The ANOVA results as well as post-hoc Tukey HSD on magnitude

Ramachandran Plots.xlsx	  	The result from the six proteins and rough work involving 4900.py



All code is primarily written in python and uses objects and libraries created by myself. 
There are some functionalities such as analyzing helices and sheets that are currently not updated, but the majority of the code remains functional.



To generate the phi and psi angles for Ramachandran plots, simply use the command: python 4900.py valid_pdb_file.pdb and the code should create 
a name_of_valid_pdb_file.txt that outputs something similar to:



Residue Sequence Number | Amino Acid | Phi | Psi | Distance from Center

2 LEU -86.9199524082 121.663292073 21.3849995069

3 SER -76.2210055178 171.134822774 24.4887606323

4 PRO -59.6996415098 -40.6502638334 27.8152472225

5 ALA -67.4581643516 -37.6316116506 26.7488192941

6 ASP -62.5102906578 -42.4955902714 23.0363722521

…


It is possible to calculate the dihedral angles for multiple files at once. In order to do this, simply use the command: 
python 4900.py all folder_name_with_pdb_files¬ and the program will calculate the dihedral angles for every single PDB file located within the 
folder and create the same output.txt as the one seen above. Note that the program will not create individual dihedral angles for every single PDB file, 
but will combine every single file into one large text file. 

It is possible to set a set amount of files for the program to calculate. 


On line 31, there is a variable called: “only_100” which is currently set at 100, by changing the numerical value of this variable, the program will 
calculate the number of files indicated by this numerical value whether it be 1000 or 10.



In addition to generating the dihedral angles for the creation of a Ramachandran plot, the program can also align two sequences either given the 
original PDB file or the produced output.txt that was created containing the dihedral angles. To align the two sequences using the Smith-Waterman algorithm, 
simply use the command: python 4900.py align name_of_1st_file.pdb.txt name_of_2nd_file.pdb.txt. This will produce an output that is exactly the same as 
the output shown above, however there will be dashes wherever the algorithm determines there is an insertion or deletion event that has occurred. Similar to:



Residue Sequence Number | Amino Acid | Phi | Psi | Distance from Center

8 THR -76.76965996 -43.2136435723 21.7505985422
9 ASN -60.5833946488 -48.1678156214 18.0376890921

10 VAL -56.819949344 -46.8358462565 17.0042297568

11 LYS -57.0985494695 -47.6622868463 21.5732900989

-

-

-

12 GLY -52.322497643 -54.8021946285 19.6792601157

13 VAL -76.0858751628 -24.2011010745 16.4231124152

...



There are other functionalities of the program as well including the ability to calculate the x, y, and z coordinates of the “center” of the protein 
based on molecular mass and the distance that each amino acid within the protein is from that calculated center. The command for both of these functions 
are: python 4900.py PDB_file.pdb center and python 4900.py PDB_file.pdb distance respectively. The output for the center functionality will be a single 
sentence stating the x, y, and z coordinates of the proteins center. The output for the distance functionality will be a text file called, distance.txt 
that will contain the amino acid name and their relative distance away from the center of the protein.