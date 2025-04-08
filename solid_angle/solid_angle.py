import pandas as pd
import numpy as np
import sys

#read input file
input_file = sys.argv[1]
dataset = pd.read_excel(input_file)
df = dataset.to_numpy()
df = df.T

#define atoms and atomic coordinates
Atom = df[0]
X = df[1]
Y = df[2]
Z = df[3]

#specify the four atoms for solid angle calculation
P0 = np.array([X[int(sys.argv[2]) - 1], Y[int(sys.argv[2]) - 1], Z[int(sys.argv[2]) - 1]])#out of plane atom
P1 = np.array([X[int(sys.argv[3]) - 1], Y[int(sys.argv[3]) - 1], Z[int(sys.argv[3]) - 1]])#first in-plane atom
P2 = np.array([X[int(sys.argv[4]) - 1], Y[int(sys.argv[4]) - 1], Z[int(sys.argv[4]) - 1]])#second in-plane atom
P3 = np.array([X[int(sys.argv[5]) - 1], Y[int(sys.argv[5]) - 1], Z[int(sys.argv[5]) - 1]])#third in-plane atom

#define solid angle calculation
def calculate_solid_angle(point, plane_points):

    vectors = [point - p for p in plane_points]

    cross_product = np.cross(vectors[0], vectors[1])

    cross_product_mag = np.linalg.norm(cross_product)

    dot_product = np.dot(vectors[0], vectors[1])

    solid_angle = 2 * np.arctan2(cross_product_mag, dot_product)

    return solid_angle

#calculate solid angle
point = P0
plane_points = P1, P2, P3

solid_angle_rad = calculate_solid_angle(point, plane_points)
solid_angle_deg = np.degrees(solid_angle_rad)

#print solid angle
statement = f"Solid angle from point {point} to plane defined by {plane_points} is {solid_angle_rad:.4f} radians ({solid_angle_deg:.4f} degrees)."
np.savetxt("output.txt", [statement], fmt='%s')