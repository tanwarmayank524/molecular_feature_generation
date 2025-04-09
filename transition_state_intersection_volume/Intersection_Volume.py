import pandas as pd
import math
import sys

#read input file and create lists of atoms and their cartesian coordinates
input_file = sys.argv[1]
dataset = pd.read_excel(input_file)
dataset = dataset.fillna(0)

#substrate 1
Atom1 = (dataset.iloc[:,[0]]).values.tolist()
X1 = (dataset.iloc[:,[1]]).values.tolist()
Y1 = (dataset.iloc[:,[2]]).values.tolist()
Z1 = (dataset.iloc[:,[3]]).values.tolist()

#substrate 2
Atom2 = (dataset.iloc[:,[4]]).values.tolist()
X2 = (dataset.iloc[:,[5]]).values.tolist()
Y2 = (dataset.iloc[:,[6]]).values.tolist()
Z2 = (dataset.iloc[:,[7]]).values.tolist()

#element dictionary containing vanderwaal radius
element_dict = {
    "H" : 1.20,
    "C" : 1.70,
    "N" : 1.55,
    "O" : 1.52,
    "F" : 1.47,
    "Si" : 2.10, 
    "P" : 1.80,
    "S" : 1.80,
    "Cl" : 1.75,
}

#define scaling factors
scaling = [1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
#initialize 
distance = [[0 for j in range(len(Atom2))] for i in range(len(Atom1))]
volume_scaling = [[0 for f in range(len(scaling))] for i in range(len(Atom1))]#this gives intersection volume wrt each atom for every scaling factor
intersection_volume = [[0 for f in range(len(scaling))] for i in range(1)]#this gives total intersection volume for all atoms for every scaling factor
pi = math.pi

#calculate intersection volume
for f in range(len(scaling)):
    for i in range(len(Atom1)):
        V = 0.
        atom1 = Atom1[i][0]
        if atom1 == 0:
            continue
        r1 = element_dict.get(atom1)
        r1 = r1*scaling[f]
        for j in range(len(Atom2)):
            atom2 = Atom2[j][0]
            if atom2 == 0:
                continue
            r2 = element_dict.get(atom2)
            r2 = r2*scaling[f]
            dist = math.dist((X1[i][0], Y1[i][0], Z1[i][0]), (X2[j][0], Y2[j][0], Z2[j][0]))
            distance[i][j] = dist
            if r1 + r2 <= dist:
                V = V
            elif dist == 0:
                V = V
            elif abs(r1-r2) >= dist:
                if r1 < r2:
                    V = 4*pi*r1*r1*r1/3 + V
                else:
                    V = 4*pi*r2*r2*r2/3 + V
            else:
                V = pi*((r1+r2-dist)**2)*(dist**2+2*dist*(r1+r2)+6*r1*r2-3*(r1**2+r2**2))/(12*dist) + V

        volume_scaling[i][f] = V
        intersection_volume[0][f] += V     

#write intersection volume
df1 = pd.DataFrame(intersection_volume)
df2 = pd.DataFrame(volume_scaling)
with pd.ExcelWriter('Volume.xlsx', engine='openpyxl') as writer:
    df1.to_excel(writer, sheet_name='Sheet1', index=True)
    df2.to_excel(writer, sheet_name='Sheet2', index=True)