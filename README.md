# molecular_feature_generation
A collection of codes to generate novel steric features for homogeneous catalysis to serve as descriptors for both activity and selectivity. More details on the development of these descriptors can be found at https://doi.org/10.26434/chemrxiv-2024-l2jgc
These molecular features are initially developed for hydrogen atom abstraction reactions, but can be further extended to other solution-phase reactions. 


## solid angle
Solid angle is quantified by the angle formed at any vertex atom by any three atoms. The script can be run with:

```python solid_angle.py *.xlsx atom1 atom2 atom3 atom4```

Here, atom1 is the index of the designated vertex atom, and atom2-4 are the indexes of the three atoms forming this solid angle. *.xlsx is the input file containing corrdinates of the molecule.

![Alt text](https://github.com/tanwarmayank524/molecular_feature_generation/blob/main/solid_angle/solid_angle.png)

## ground state intersection volume
This descriptor quantifies the sterics at any particular atom in a molecule due to the rest of the atoms by taking the pair-wise intersection between atoms and summing them up. Atoms here are assumed as soft van der Waals spheres and their radius corresponding to their van der Waals radii.
The higher the intersection volume at a particular atom, the higher the sterics. The script can be run with:

```python intersection_volume.py *.xlsx atom1```

Here, atom1 is the index of the atom whose sterics are being quantified. *.xlsx is the input file containing corrdinates of the molecule.

![Alt text](https://github.com/tanwarmayank524/molecular_feature_generation/blob/main/ground_state_intersection_volume/Intersection_Volume.png)

## transition state intersection volume
This descriptor quantifies the transition state sterics between two molecules. The higher the value, the higher the sterics. The script can be run with:

```python intersection_volume.py *.xlsx```

Here, *.xlsx is the input file containing corrdinates of both the molecules.

![Alt text](https://github.com/tanwarmayank524/molecular_feature_generation/blob/main/transition_state_intersection_volume/Intersection_Volume.png)




# Authors
Mayank Tanwar
email: tanwa008@umn.edu
GitHub: tanwarmayank524

# Citation
Publications relevant to the code: https://doi.org/10.26434/chemrxiv-2024-l2jgc
# Acknowledgements
NSF Center for Synthetic Organic Electrochemistry (https://cci.utah.edu/)

 
