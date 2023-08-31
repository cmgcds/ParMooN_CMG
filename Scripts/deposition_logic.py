#!/usr/bin/env python3

# Generates deposition fraction
# Usage: python deposition_fraction.py <latest particle position csv> <save/show> <particle diameter>
# save/show -> optional, default is show
# particle diameter -> optional, default is 4.3; 0 -> don't add from paper

# %%
import sys
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
from cycler import cycler
from tqdm import tqdm

tqdm.pandas()

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Computer Modern"]
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['legend.fontsize'] = 20
plt.rcParams['ytick.labelsize'] = 20
plt.rcParams['axes.prop_cycle'] = cycler(color=['darkblue', '#d62728', '#2ca02c', '#ff7f0e', '#bcbd22', '#8c564b', '#17becf', '#9467bd', '#e377c2', '#7f7f7f'])

def sameSideOfPlane(point1, point2, point3, refPoint, point):
    x1, y1, z1 = point1
    x2, y2, z2 = point2
    x3, y3, z3 = point3
    x, y, z = point
    rX, rY, rZ = refPoint
    a = (y2 - y1)*(z3 - z1) - (z2 - z1)*(y3 - y1)
    b = (z2 - z1)*(x3 - x1) - (x2 - x1)*(z3 - z1)
    c = (x2 - x1)*(y3 - y1) - (y2 - y1)*(x3 - x1)
    d = -a*x1 - b*y1 - c*z1
    val = a*x + b*y + c*z + d
    refVal = a*rX + b*rY + c*rZ + d
    return (val > 0 and refVal > 0) or (val < 0 and refVal < 0)

def inSphere(point, seg):
    x, y, z = point
    center, radius = seg
    cX, cY, cZ = center
    return (x - cX)**2 + (y - cY)**2 + (z - cZ)**2 < radius**2


seg3 = ([-0.00916789, 0.107428, -0.191593], 0.010)
seg4 = ([0.0118495, 0.109477, -0.204223], 0.019)
seg5 = ([0.0316239, 0.113915, -0.214249], 0.0066)
seg6a = ([0.0358753, 0.116738, -0.224577], 0.008)
seg6b = ([0.0333019, 0.126099, -0.224295], 0.0065)
seg6c = ([0.0408924, 0.119268, -0.231801], 0.007)
seg7a = ([0.0445801, 0.116936, -0.21015], 0.010)
seg7b = ([0.0500298, 0.115545, -0.214051], 0.006)
seg8 = ([-0.0185997, 0.107727, -0.202463], 0.008)
seg9a = ([-0.0347457, 0.112414, -0.19602], 0.0127)
seg9b = ([-0.0428033, 0.121789, -0.194985], 0.00625265)
seg10 = ([-0.0280901, 0.105663, -0.219729], 0.012)
seg11 = ([-0.0469778, 0.0961967, -0.223908], 0.0136637)
seg12 = ([-0.0357696, 0.113489, -0.233305], 0.009)
seg13 = ([0.0614592, 0.133665, -0.179065], 0.031)
seg14 = ([0.0257088, 0.166302, -0.228401], 0.0388246)
seg15a = ([0.048863, 0.125997, -0.25725], 0.025)
seg15b = ([0.0578417, 0.13668, -0.278319], 0.033)
seg16a = ([0.0678228, 0.12313, -0.218891], 0.0187)
seg16b = ([0.0889196, 0.132538, -0.227953], 0.028)
seg17a = ([-0.0670284, 0.0713526, -0.234222], 0.0212804)
seg17b = ([-0.0871651, 0.0550605, -0.247243], 0.03)
seg18 = ([-0.0396406, 0.102922, -0.272769], 0.0361585)
seg19 = ([-0.0723341, 0.101592, -0.173525], 0.0381007)
seg20 = ([-0.0772369, 0.110594, -0.230988], 0.0245135)
seg21a = ([-0.0448059, 0.150873, -0.187616], 0.0276968)
seg21b = ([-0.0537398, 0.164288, -0.188217], 0.0359187)
seg22a = ([-0.0407379, 0.138604, -0.23724], 0.0200263)
seg22b = ([-0.0450925, 0.156375, -0.245063], 0.0263275)

class LeftLung:
    def __init__(self, point):
        self.point = point
        pass

    def seg34(self):
        val = sameSideOfPlane([0.00396, 0.1099, -0.191382],
                              [0.001735, 0.1133, -0.2011],
                              [0.001649, 0.10254, -0.2000],
                              [-0.007, 0.1, -0.193],
                              self.point)
        return 3 if val else 4

    def seg45(self):
        val = sameSideOfPlane([0.0257677,0.114642,-0.211365],
                              [0.0264381,0.107334,-0.211896],
                              [0.0242488,0.112508,-0.215532],
                              [0.0164216, 0.108844, -0.207705],
                              self.point)
        return 4 if val else 5

    def seg56(self):
        val = sameSideOfPlane([0.0343477, 0.111998, -0.21962],
                              [0.03176, 0.114367, -0.221207],
                              [0.036208, 0.115116, -0.217883],
                              [0.0332826, 0.117274, -0.214308],
                              self.point)
        return 5 if val else 6

    def seg57(self):
        val = sameSideOfPlane([0.0387031, 0.111458, -0.211436],
                              [0.0408733, 0.114347, -0.215259],
                              [0.038741, 0.117294, -0.212645],
                              [0.0332826, 0.117274, -0.214308],
                              self.point)
        return 5 if val else 7

    def seg614(self):
        val = sameSideOfPlane([0.0326649, 0.125194, -0.224243],
                              [0.0316696, 0.124109, -0.22695],
                              [0.0339185, 0.125326, -0.226154],
                              [0.367313, 0.113588, -0.224294],
                              self.point)
        return 6 if val else 14

    def seg615(self):
        val = sameSideOfPlane([0.0403648, 0.115793, -0.231076],
                              [0.0426591, 0.118383, -0.22913],
                              [0.0392978, 0.120591, -0.231464],
                              [0.367313, 0.113588, -0.224294],
                              self.point)
        return 6 if val else 15

    def seg713(self):
        val = sameSideOfPlane([0.0451224, 0.119878, -0.203707],
                              [0.0477829, 0.121177, -0.206496],
                              [0.0482063, 0.119171, -0.206787],
                              [0.0434574, 0.112952, -0.210408],
                              self.point)
        return 7 if val else 13

    def seg716(self):
        val = sameSideOfPlane([0.0519382, 0.114605, -0.212265],
                              [0.515506, 0.116596, -0.213827],
                              [0.0511742, 0.116374, -0.215987],
                              [0.0434574, 0.112952, -0.210408],
                              self.point)
        return 7 if val else 16

    def seg1323(self):
        val = sameSideOfPlane([0.059, 0.121, -0.178],
                              [0.050, 0.126, -0.172],
                              [0.055, 0.149, -0.187],
                              [0.054, 0.119, -0.198],
                              self.point)
        return 13 if val else 23

    def seg1424(self):
        val = sameSideOfPlane([0.036, 0.168, -0.221],
                              [0.046, 0.166, -0.235],
                              [0.044, 0.162, -0.247],
                              [0.029, 0.133, -0.235],
                              self.point)
        return 14 if val else 24

    def seg1525(self):
        val = sameSideOfPlane([0.062, 0.108, -0.276],
                              [0.041, 0.119, -0.284],
                              [0.067, 0.155, -0.266],
                              [0.045, 0.121, -0.246],
                              self.point)
        return 15 if val else 25

    def seg1626(self):
        val = sameSideOfPlane([0.089, 0.108, -0.226],
                              [0.078, 0.118, -0.242],
                              [0.074, 0.141, -0.217],
                              [0.068, 0.111, -0.219],
                              self.point)
        return 16 if val else 26

class RightLung:
    def __init__(self, point):
        self.point = point
        pass

    def seg38(self):
        val = sameSideOfPlane([-0.0155714, 0.101882, -0.1944],
                              [-0.0116865, 0.109269, -0.20228],
                              [-0.0151642, 0.114347, -0.196612],
                              [-0.009, 0.102, -0.192],
                              self.point)
        return 3 if val else 8

    def seg89(self):
        val = sameSideOfPlane([-0.0262068, 0.111116, -0.196478],
                              [-0.0283781, 0.108044, -0.20268],
                              [-0.0279246, 0.11303, -0.20205],
                              [-0.014958, 0.103675, -0.202574],
                              self.point)
        return 8 if val else 9

    def seg810(self):
        val = sameSideOfPlane([-0.0284162, 0.107697, -0.21279],
                              [-0.0226891, 0.108775, -0.214457],
                              [-0.0236845, 0.106937, -0.216261],
                              [-0.014958, 0.103675, -0.202574],
                              self.point)
        return 8 if val else 10

    def seg919(self):
        val = sameSideOfPlane([-0.0451173, 0.113083, -0.191222],
                              [-0.0453529, 0.114524, -0.192305],
                              [-0.0457442, 0.112277, -0.195197],
                              [-0.0335852, 0.110919, -0.195308],
                              self.point)
        return 9 if val else 19

    def seg921(self):
        val = sameSideOfPlane([-0.0422019, 0.120862, -0.193467],
                              [-0.0427378, 0.120182, -0.195737],
                              [-0.0401117, 0.121948, -0.195337],
                              [-0.0335852, 0.110919, -0.195308],
                              self.point)
        return 9 if val else 21

    def seg1011(self):
        val = sameSideOfPlane([-0.0415394, 0.0974918, -0.227024],
                              [-0.0416943, 0.100341, -0.224749],
                              [-0.0421243, 0.102724, -0.228819],
                              [-0.0323595, 0.103872, -0.227816],
                              self.point)
        return 10 if val else 11

    def seg1012(self):
        val = sameSideOfPlane([-0.036299, 0.108656, -0.230138],
                              [-0.0326029, 0.106679, -0.232364],
                              [-0.033017, 0.111906, -0.229614],
                              [-0.0323595, 0.103872, -0.227816],
                              self.point)
        return 10 if val else 12

    def seg1117(self):
        val = sameSideOfPlane([-0.0525026, 0.08941, -0.222551],
                              [-0.0518987, 0.0889323, -0.226285],
                              [-0.0535424, 0.0903605, -0.225974],
                              [-0.0497583, 0.0984493, -0.225411],
                              self.point)
        return 11 if val else 17

    def seg1120(self):
        val = sameSideOfPlane([-0.0562036, 0.10034, -0.225532],
                              [-0.0561035, 0.0996118, -0.227661],
                              [-0.055674, 0.103075, -0.227932],
                              [-0.0497583, 0.0984493, -0.225411],
                              self.point)
        return 11 if val else 20

    def seg1218(self):
        val = sameSideOfPlane([-0.0380397, 0.107117, -0.237818],
                              [-0.0388701, 0.110965, -0.238563],
                              [-0.0342797, 0.108929, -0.23908],
                              [-0.0311773, 0.110894, -0.232258],
                              self.point)
        return 12 if val else 18

    def seg1222(self):
        val = sameSideOfPlane([-0.0351501, 0.118375, -0.230849],
                              [-0.0373866, 0.117661, -0.234033],
                              [-0.0329818, 0.117217, -0.234101],
                              [-0.0311773, 0.110894, -0.232258],
                              self.point)
        return 12 if val else 22

    def seg1727(self):
        val = sameSideOfPlane([-0.076, 0.055, -0.232],
                              [-0.080, 0.057, -0.229],
                              [-0.070, 0.058, -0.245],
                              [-0.069, 0.073, -0.233],
                              self.point)
        return 17 if val else 27

    def seg1828(self):
        val = sameSideOfPlane([-0.047, 0.085, -0.270],
                              [-0.062, 0.094, -0.267],
                              [-0.031, 0.088, -0.272],
                              [-0.038, 0.092, -0.261],
                              self.point)
        return 18 if val else 28

    def seg1929(self):
        val = sameSideOfPlane([-0.078, 0.084, -0.195],
                              [-0.066, 0.091, -0.173],
                              [-0.057, 0.100, -0.154],
                              [-0.057, 0.112, -0.171],
                              self.point)
        return 19 if val else 29

    def seg2030(self):
        val = sameSideOfPlane([-0.080, 0.109, -0.222],
                              [-0.083, 0.095, -0.235],
                              [-0.080, 0.098, -0.244],
                              [-0.072, 0.097, -0.232],
                              self.point)
        return 20 if val else 30

    def seg2131(self):
        val = sameSideOfPlane([-0.056, 0.154, -0.165],
                              [-0.076, 0.146, -0.191],
                              [-0.065, 0.157, -0.219],
                              [-0.073, 0.142, -0.201],
                              self.point)
        return 21 if val else 31

    def seg2232(self):
        val = sameSideOfPlane([-0.055, 0.146, -0.225],
                              [-0.061, 0.143, -0.235],
                              [-0.049, 0.139, -0.255],
                              [-0.050, 0.144, -0.234],
                              self.point)
        return 22 if val else 32


moreSeg = False
def categorise(row):  
    if row['z'] >= -0.06:
        return 1
    if row['z'] >= -0.187 and row['x'] >= -0.0165 and row['x'] < 0.020:
        return 2

    point = [row['x'], row['y'], row['z']]
    ll = LeftLung(point)
    rl = RightLung(point)

    # Left Lung
    if inSphere(point, seg3):
        return 3
    if inSphere(point, seg4):
        return 4
    if inSphere(point, seg5):
        return 5

    if inSphere(point, seg6a) or inSphere(point, seg6b) or inSphere(point, seg6c):
        return 6
    if inSphere(point, seg14):
        # return 14 if ll.seg1424() == 14 else 24
        return 14
    if inSphere(point, seg15a) or inSphere(point, seg15b):
        # return 15 if ll.seg1525() == 15 else 25
        return 15

    if inSphere(point, seg7a) or inSphere(point, seg7b):
        return 7
    if inSphere(point, seg13):
        # return 13 if ll.seg1323() == 13 else 23
        return 13
    if inSphere(point, seg16a) or inSphere(point, seg16b):
        # return 16 if ll.seg1626() == 16 else 26
        return 16

    # Right Lung
    if inSphere(point, seg8):
        return 8
    if inSphere(point, seg9a) or inSphere(point, seg9b):
        return 9
    if inSphere(point, seg19):
        # return 19 if rl.seg1929() == 19 else 29
        return 19
    if inSphere(point, seg21a) or inSphere(point, seg21b):
        # return 21 if rl.seg2131() == 21 else 31
        return 21

    if inSphere(point, seg10):
        return 10
    if inSphere(point, seg11):
        return 11
    if inSphere(point, seg17a) or inSphere(point, seg17b):
        # return 17 if rl.seg1727() == 17 else 27
        return 17
    if inSphere(point, seg20):
        # return 20 if rl.seg2030() == 20 else 30
        return 20

    if inSphere(point, seg12):
        return 12
    if inSphere(point, seg18):
        # return 18 if rl.seg1828() == 18 else 28
        return 18
    if inSphere(point, seg22a) or inSphere(point, seg22b):
        # return 22 if rl.seg2232() == 22 else 32
        return 22

    return -1


argc = len(sys.argv)

if(argc < 3):
    print("Error: Not enough arguments")
    print("Sample : python deposition_fraction.py <latest particle position csv> <particle diameter>")
    exit(0)

filename = sys.argv[1]
particle = sys.argv[2] 

#available particle sizes for monodisperse
available_particles = ["2.5", "4.3", "8", "10"]
if(str(particle) not in available_particles):
    print("Error: Invalid particle size")
    print("Available particle sizes: ", available_particles)
    exit(0);
 
 
# Read the csv file. 
df = pd.read_csv(filename)
df['section'] = df.apply(lambda row: categorise(row), axis=1)



# Get only the particles, which are deposited and not marked as errors
df_grouped = pd.DataFrame(df[(df['error'] == 0) & (df['deposition'] == 1)])


print("[INFO] Particles without Error : ", df_grouped.shape[0])
effective_particles = df_grouped.shape[0]

escaped_particles = df_grouped[df_grouped['escaped'] == 1].shape[0]
print("[INFO] Particles escaped : ", escaped_particles)

deposited_particles = effective_particles - escaped_particles
print("[INFO] Particles deposited : ", deposited_particles)

deposition_fraction = ( (effective_particles - escaped_particles) / effective_particles ) 
print("[INFO] Deposition fraction : ", deposition_fraction)

# get only the deposited particles and ignore the escaped particles
df_grouped = df_grouped[(df_grouped['deposition'] == 1) & (df_grouped['escaped'] == 0) ]

particles = df_grouped.shape[0] #Archive code for backward compatibility

# get only section column
df_grouped = df_grouped[['section']]

df_grouped['count'] = 1
df_grouped = df_grouped.groupby('section').sum('count').reset_index()

#Remove all the particles which are tagged as -1
df_grouped = df_grouped[df_grouped['section'] != -1]

# print the number of particles, which are tagged as -1
print("[INFO] Particles which are not zoned : ", df_grouped[df_grouped['section'] == -1].shape[0])


# Sum count column
final_total = df_grouped['count'].sum()

df_grouped['Deposition fraction'] = (df_grouped['count'] / final_total)    * 100 * deposition_fraction



print("Sum of Deposition fraction : ", df_grouped['Deposition fraction'].sum())

# print("Total particles: ", df.shape[0])
# print("Total deposited: ", particles)
# print("Deposition percentage: ", particles / df.shape[0] * 100)
# print("Total escaped: ", df[df['escaped'] == 1].shape[0])
# print("Escape percentage: ", df[df['escaped'] == 1].shape[0] / df.shape[0] * 100)
# print("Total stagnant: ", df[df['deposition'] == 0].shape[0])
# print("Stagnant percentage: ", df[df['deposition'] == 0].shape[0] / df.shape[0] * 100)
# print("Total error: ", df[df['error'] == 1].shape[0])
# print("Error percentage: ", df[df['error'] == 1].shape[0] / df.shape[0] * 100)
print(df_grouped)

def les1(size):
    if size == "2.5":
        return [1.383, 0.162, 0.281, 0.077, 0.355, 0.206, 0.180, 0.330, 0.285, 0.246, 0.623, 0.289, 0.090, 0.121, 0.400, 0.272, 0.492, 0.361, 0.355, 0.642, 0.928, 0.412]
    if size == "4.3":
        return [3.113, 0.454, 0.749, 0.295, 1.695, 1.347, 1.259, 0.972, 1.280, 1.079, 3.840, 2.167, 0.257, 0.651, 1.829, 1.291, 2.654, 1.686, 4.334, 4.634, 2.537, 4.704]
    if size == "8":
        return [22.495, 4.569, 5.376, 3.664, 4.706, 2.848, 1.602, 6.231, 3.502, 5.142, 7.768, 4.992, 0.167, 0.082, 6.044, 1.361, 1.554, 6.900, 4.991, 2.385, 1.381, 1.123]
    if size == "10":
        return [42.115, 8.121, 6.752, 4.776, 2.232, 1.488, 0.446, 5.784, 1.613, 4.180, 4.272, 2.492, 0.023, 0.005, 2.264, 0.297, 0.192, 3.767, 1.131, 0.272, 0.156, 0.076]

def les2(size):
    if size == "2.5":
        return [12.478, 3.354, 1.092, 0.478, 0.901, 1.060, 0.766, 0.744, 0.862, 0.801, 1.404, 1.247, 0.888, 0.813, 4.064, 2.122, 2.571, 2.810, 3.354, 2.154, 2.285, 2.353]
    if size == "4.3":
        return [17.977, 5.771, 1.797, 0.849, 1.376, 1.841, 1.507, 1.368, 1.747, 1.660, 3.122, 2.834, 0.857, 0.733, 7.614, 3.547, 3.960, 6.137, 7.182, 3.707, 4.274, 4.637]
    if size == "8":
        return [37.706, 10.597, 2.978, 2.726, 2.645, 1.750, 0.999, 4.371, 2.316, 4.183, 6.048, 3.402, 0.218, 0.004, 4.917, 1.321, 1.600, 6.412, 4.115, 1.227, 0.998, 0.998]
    if size == "10":
        return [48.105, 9.008, 3.854, 3.205, 1.567, 0.812, 0.295, 5.217, 1.166, 4.057, 4.336, 2.567, 0.013, 0.001, 2.073, 0.306, 0.291, 3.603, 0.893, 0.144, 0.123, 0.026]

def rans1(size):
    if size == "2.5":
        return [18.460, 6.385, 1.178, 1.161, 1.015, 1.061, 0.601, 1.589, 1.127, 1.389, 1.873, 1.368, 0.510, 0.517, 2.526, 1.030, 1.791, 2.310, 2.310, 1.818, 1.368, 1.738]
    if size == "4.3":
        return [20.087, 7.173, 1.506, 1.534, 1.328, 1.765, 0.885, 1.831, 1.409, 1.928, 3.059, 2.123, 0.548, 0.403, 2.970, 1.285, 2.227, 3.234, 3.293, 2.219, 1.540, 2.630]
    if size == "8":
        return [25.332, 7.774, 3.308, 3.983, 2.982, 2.299, 1.053, 4.103, 1.996, 4.356, 4.387, 5.204, 0.311, 0.165, 5.321, 1.093, 1.259, 7.005, 3.117, 1.287, 0.843, 1.952]
    if size == "10":
        return [34.614, 7.092, 4.738, 5.420, 2.275, 1.960, 0.558, 4.528, 1.178, 4.331, 2.642, 5.257, 0.099, 0.104, 3.358, 0.378, 0.302, 5.927, 1.214, 0.272, 0.168, 0.502]

def rans3(size):
    if size == "2.5":
        return [1.475, 0.123, 0.103, 0.049, 0.109, 0.073, 0.039, 0.176, 0.112, 0.114, 0.128, 0.071, 0.054, 0.053, 0.211, 0.132, 0.163, 0.136, 0.190, 0.284, 0.231, 0.204]
    if size == "4.3":
        return [2.992, 0.503, 0.462, 0.307, 0.878, 0.869, 0.967, 0.722, 1.146, 1.117, 1.646, 1.030, 0.186, 0.397, 1.670, 1.425, 1.452, 1.133, 2.054, 2.123, 1.866, 1.845]
    if size == "8":
        return [7.271, 1.177, 1.317, 0.971, 1.720, 1.562, 1.356, 1.366, 1.882, 1.516, 4.519, 4.104, 0.635, 0.721, 3.616, 2.512, 2.727, 2.299, 3.725, 3.187, 2.768, 2.830]
    if size == "10":
        return [14.112, 2.275, 2.346, 2.144, 3.958, 2.763, 1.989, 1.765, 2.208, 1.791, 6.580, 5.420, 0.887, 0.887, 4.017, 3.069, 2.526, 2.311, 3.459, 2.979, 2.490, 2.310]
    
def paperDf(size):
    size = str(size)
    return les1(size), les2(size), rans1(size), rans3(size)


# 0.5 polydisperse paper data
poly = {}
poly["median"] = [2872.984833, 992.0367886, 329.1251811, 299.0160116, 265.222956, 218.9161379, 210.3374801, 737.9986299, 303.8357614, 326.5042877, 514.9980155, 287.2984833, 157.7308583, 359.3813664, 540.3048785, 365.1741273, 395.5689752, 395.5689752, 575.9922854, 644.2104683, 408.4238653, 386.1940297]
poly["mean"] = [2120.258871, 599.4842503, 301.416253, 218.9161379, 220.6734069, 183.6068538, 177.827941, 549.0138914, 250.7873091, 273.8419634, 377.0412694, 210.3374801, 104.0785208, 121.1527659, 194.1755055, 213.7278471, 164.1639442, 164.1639442, 506.8285968, 359.3813664, 235.2489955, 197.3053625]
poly["max"] = [3622.661685, 2103.374801, 438.8954829, 389.2940606, 992.0367886, 2004.856687, 278.2559402, 2468.090675, 398.7442601, 380.0678299, 793.059105, 415.0071275, 859.0684664, 540.3048785, 1281.264782, 467.8847485, 692.2736124, 901.2828581, 726.2917502, 901.2828581, 618.9658189, 649.3816316]

# Generate Plot
plots = list(range(1, 24, 2))
plt.figure(figsize=(6.4,4.8),dpi=300)
if (particle != 0):
    if particle == "poly":
        plt.plot(list(range(1, 23)), poly["median"], marker="x", label="Polydisperse(median)")
        plt.scatter(list(range(1, 23)), poly["min"], marker="x", label="Polydisperse(min)")
        plt.scatter(list(range(1, 23)), poly["max"], marker="x", label="Polydisperse(max)")
    else:
        les1Values, les2Values, rans1Values, rans3Values = paperDf(particle)
        # print(les1Values)
        plt.plot(list(range(1, 23)), les1Values, marker="o", label="LES1",linestyle='dashed')
        plt.plot(list(range(1, 23)), les2Values, marker="s", label="LES2",linestyle='dashed')
        # plt.plot(list(range(1, 23)), rans1Values, marker="v", label="RANS1")
        # plt.plot(list(range(1, 23)), rans3Values, marker="1", label="RANS3")
plt.plot(df_grouped['section'].values, df_grouped['Deposition fraction'].values, marker="^", label="VMS")
plt.xlabel("Segments")
plt.ylabel("Deposition fraction " +"(\%)")
plt.yscale("log")
plt.title("Deposition fraction - Particle size " + str(particle) + " $\mu m$")
plt.legend(loc="lower left", ncol=2, fontsize=10)
plt.xticks(list(range(1, 24, 2)))
plt.yticks([0.001, 0.01, 0.1, 1, 10, 100])
# plt.grid(True, which="both", ls="-", alpha=0.5)
# plt.grid(which='major', axis='both', color='k', alpha=0.6)
# plt.grid(which='minor', axis='both', color='k', alpha=0.2)
plt.grid(alpha=0.5)
plt.tight_layout()
plt.savefig(f'deposition_fraction_{particle}.png', dpi=300)
plt.show()
# %%
