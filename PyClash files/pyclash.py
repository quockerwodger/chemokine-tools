# -*- coding: utf-8 -*-

from Bio.PDB import *
import numpy as np
import math
import warnings
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
warnings.simplefilter('ignore')


###### Lookup table for element radii ######
vdw_radii = {
        "H": 1.200,
        "HE": 1.400,
        "C": 1.700,
        "N": 1.550,
        "O": 1.520,
        "F": 1.470,
        "NA": 2.270,
        "MG": 1.730,
        "P": 1.800,
        "S": 1.800,
        "CL": 1.750,
        "K": 2.750,
        "CA": 2.310,
        "NI": 1.630,
        "CU": 1.400,
        "ZN": 1.390,
        "SE": 1.900,
        "BR": 1.850,
        "CD": 1.580,
        "I": 1.980,
        "HG": 1.550,
    }

###### Functions for steric clash calculation ######

def compute_sphere(val, res):
        '''
        This function will create a sphere for each atom based on a given resolution
        '''
        z_ang = math.radians(180/res) #Angle of layers along 'z' axis
        r = val #radius
        arc = (math.pi*r)/res #Distance separating each point on a layer
        
        d = 0 #angle between Z axis and plane
        t = 0 #angle between Y axis & point
        
        points = 0 # no. points to be generated
        for value in range(2,(res)): #calculates the number of points for the whole sphere
            theta = math.radians(value*(180/res))
            x = (2*(math.pi)*(math.sin(theta))*r)//arc
            points += int(x)

        coords = np.zeros(((int(points)+2),3), dtype=np.float32) #create empty array of right size
        
        #generates the points
        for k in range(int(points)+2): 
            if k == 0: #The top point of the sphere
                coords[k,0]= 0
                coords[k,1]=0
                coords[k,2]=r
                d+=z_ang
            elif k < (int(points)+2): #Calculating the middle layers
                z = math.cos(d)*r
                coords[k,2]=z
                w = math.sin(d)*r
                n_points = (2*math.pi*w)//arc
                a = (math.radians(360))/n_points #angle for a given layer
                coords[k,0]=math.sin(t)*w
                coords[k,1]=math.cos(t)*w
                t += a
                if t >= (math.radians(360)):
                    d += z_ang
                    if (math.degrees(d))%4 == 2:
                        t = 0 
                    else:
                        t = a/2
                #print(d)
            if k == (int(points)+2-1): #Bottom of sphere
                coords[k,0]=0
                coords[k,1]=0
                coords[k,2]=-r
        return coords

def compute_a_value(target_coord, target_radii, invading_sphere, res, invading_coords):
    '''
    This function will compute the height of the overlapping sphere cap (a-value) for the invading sphere 
    i.e. it takes the centre of one sphere (target) and calculates the size of the cap of an intersecting sphere using the points of that sphere (invading_sphere)
    '''
    volume_vars = []
    threshold = 1 / res #The threshold is set depending on the resolution, i.e. fewer points means you'll need a larger window to ensure there are points within the threshold

    for index, (centre, radius) in enumerate(zip(target_coord, target_radii)):
        target_distance = int(radius)*3 #This ensures that only atoms that are actually close to each other are measured, speeds up computation
        for i, (sphere, origin) in enumerate(zip(invading_sphere, invading_coords)):
            if(abs(math.dist(centre, origin)))<target_distance: #Select only spheres within range
                point_distances = np.array([math.dist(point, centre) for point in sphere]) #distance from each point in invading sphere to target centre
                within_threshold = (np.abs(point_distances - radius) <= threshold) & (point_distances <= radius) #filters out those points that are sitting on the outer edge of the target radius
                key_points = sphere[within_threshold]
                distances = point_distances[within_threshold]
                passed_threshold = distances - radius
                if len(distances) > 1: #This will filter out the few edge cases where only one point from the invading sphere is actually within range
                    key_triple = list(zip(distances, key_points, passed_threshold))
                    sorted_triple = sorted(key_triple, key=lambda x: x[0], reverse=False)
                    sorted_distance, sorted_key_point, sorted_threshold = zip(*sorted_triple)
                    upper_point = sorted_key_point[0]
                    all_a_values = [math.dist(upper_point, point) for point in sorted_key_point[1:]]
                    a_value = max(all_a_values) / 2
                    holding_vars = np.array([i, a_value, index])
                    volume_vars.append(holding_vars)

    return volume_vars

def res_a_value(target_coord, target_radii, invading_sphere, res):
    '''
    A modified version of compute_a_value that specifically calculates the intersects within residues to allow the residue volume to be calculated
    '''
    volume_vars = []
    for index, (centre, radius) in enumerate(zip(target_coord, target_radii)):
        for i, sphere in enumerate(invading_sphere):
            if i != index: #This is the key difference. compute_a_value works for comparing between chains. Because these are within the same chain, this test ensures that there is no same-same comparison
                distance = []
                key_point = []
                threshold = 1/res
                passed_threshold = []
                for point in sphere:
                    if (abs(math.dist(point, centre)-radius) <= threshold) and (math.dist(point, centre) <= radius): #This tests for being within the radius & being close to the radius. It has to be less than the radius as points outside the target radius will break the formula
                        hold_distance = (math.dist(point, centre))
                        distance.append(hold_distance)
                        key_point.append(point)
                        hold_threshold = (math.dist(point, centre)-radius)
                        passed_threshold.append(hold_threshold)
                if len(distance)>1:
                    key_triple = list(zip(distance, key_point, passed_threshold))
                    sorted_triple = sorted(key_triple, key = lambda x: x[0], reverse = False)
                    sorted_distance, sorted_key_point, sorted_threshold = zip(*sorted_triple)
                    upper_point = sorted_key_point[0]
                    all_a_values = []
                    for point in sorted_key_point[1:]:
                        holding_a = math.dist(upper_point, point)
                        all_a_values.append(holding_a)
                        a_value = max(all_a_values)/2
                        holding_vars = np.array([i, a_value, index])
                    volume_vars.append(holding_vars)

    return(volume_vars)

def calc_cap_vol(a, r):
    '''
    Calculates the volume of the spherical cap
    '''
    A = (r**2 - a**2)**0.5
    vol = ((math.pi)/3)*((r-A)**2)*(2*r+A)
    return(vol)

def calc_sphere_vol(rad):
    '''
    Calculates volume of sphere
    '''
    vol = (4/3)*math.pi*rad**3
    return(vol)

def res_vol(rad, coord, res):
    '''
    Calculates the sum volume of each residue
    '''
    #First calculates the volume of each atom
    vol_atom = calc_sphere_vol(rad)
    #Then generates the spheres
    spheres = []
    for index, item in enumerate(rad):
        sphere = compute_sphere(item, res)
        spheres.append(sphere)
    #Moves the spheres to the correct coordinates so that the atom centre is the sphere centre
    for index, item in enumerate(spheres):
        for key, value in enumerate(coord):
            if index == key:
                for a, b in enumerate(item):
                    item[a]= np.add(b, value)
    #Computing the a-value for each overlap
    overlaps = res_a_value(coord, rad, spheres, res)
    cap_vol = []
    for item in overlaps:
        for i, radius in enumerate(rad):
            if item[0]==i:
                holding_cap_vol = (calc_cap_vol(item[1], radius))
                cap_vol.append([holding_cap_vol, item[0], item[2]])
    #Calculating total intersect volume:
    intersect_vol = []
    for i, A in enumerate(cap_vol):
        for j, B in enumerate(cap_vol):
            if A[1] == B[2] and A[2] == B[1]:
                vol_sum = A[0]+B[0]
                intersect_vol.append([vol_sum, A[1], B[1]])
    #Removing duplicate intersections (i.e. overlap of A:B is the same as B:A)
    pairs_to_remove = set()
    matching_pairs = set()
    for i, pair1 in enumerate(intersect_vol):
        for j, pair2 in enumerate(intersect_vol):
            if i != j and pair1[1] == pair2[2] and pair1[2] == pair2[1]:
                if tuple(pair2) not in matching_pairs:
                    matching_pairs.add(tuple(pair1))
                    pairs_to_remove.add(tuple(pair2))
    intersect_vol = [pair for pair in intersect_vol if tuple(pair) not in pairs_to_remove]
    #Calculating the full residue volume by first adding the atom volumes, then subtracting the intersect volumes
    full_vol = np.sum(vol_atom)
    all_intersect = np.sum(intersect_vol, axis = 0)
    final_vol = full_vol - all_intersect[0]
    return(final_vol)

def percent_clash(intersect_vol, chain_atoms, chain_res, res_vol):
    '''
    Calculates the proportion of the residue that clashes somewhere (i.e. sum of all clashes divided by the total volume)
    '''
    vol_list = []
    #Add residue name to each atom
    for value in intersect_vol:
        for a, atom in enumerate(chain_atoms):        
            if value[1] == a:
                residue = [value[0], atom.get_parent()]
        for p, ref in enumerate(chain_res):
            if residue[1] == ref:
                residue.append(p)
        vol_list.append(residue)
    #Sum together intersect volumes for each residue
    sums_dict = {}
    for item in vol_list:
        value, key1, key2 = item
        key = (key1, key2)
        if key in sums_dict:
            sums_dict[key] += value
        else:
            sums_dict[key] = value
    summed_array = [[value, key[0], key[1]] for key, value in sums_dict.items()]
    final_proportion = []
    #Divide the total intersect volume by the volume of the relevant residue
    for total in summed_array:
        for ref in res_vol:
            if total[2] == ref[1]:
                proportion = total[0]/ref[0]
                final_proportion.append([proportion, total[1]])
    return(final_proportion)

###### Functions for charge clash calculation ######
def sum_forces(forces):
    '''
    Will sum together the input values that have the same residue/atom
    '''
    sums_dict = {}
    for item in forces:
        value, key = item
        if key in sums_dict:
            sums_dict[key] += value
        else:
            sums_dict[key] = value
    summed_array = [[value, key] for key, value in sums_dict.items()]
    return summed_array

def calc_forces(a_coords, a_charges, b_coords, b_charges):
   '''
   Calculate the forces acting on each atom in nN. 
   If repulsive, force will be +ve. If attractive, force will be -ve
   Sums the forces together as if all acting in same direction. 
       Made this assumption as majority of forces will be miniscule relative to the main ones
   '''
   net_a = []
   net_b = []
   for index, (a_loc, a_char) in enumerate(zip(a_coords, a_charges)):
       for i, (b_loc, b_char) in enumerate(zip(b_coords, b_charges)):
           if abs(math.dist(a_loc, b_loc))<15: #Takes only those atoms within 15A of each other
               dist = (math.dist(a_loc, b_loc))/(10**10)
               e_charge = 1.602*(10**-19) #Charge of an electron in coulombs, needed to convert charge from pqr file (stored as elementary charge)
               force = ((8.987*10**9)*(((a_char*e_charge)*(b_char*e_charge))/(dist**2)))*(10**9)
               net_a.append([force, index])
               net_b.append([force, i])
   a_net_sum = sum_forces(net_a)
   b_net_sum = sum_forces(net_b)
   return a_net_sum, b_net_sum

def convert_index_res(input_array, atom_list):
    '''
    Replaces atom_index with residue name to allow forces on each residue to be calculated
    '''
    output_list = []
    for value in input_array:
        for a, atom in enumerate(atom_list):        
            if value[1] == a:
                residue = [value[0], atom.get_parent()]
        for p, ref in enumerate(atom_list):
            if residue[1] == ref:
                residue.append(p)
        output_list.append(residue)
    return output_list

###### Missing interface functions ######
def get_interface(structure):
    chains = list(structure.get_chains())
    chain_a_res = list(chains[0].get_residues())
    chain_b_res = list(chains[1].get_residues())
    a_interface_list = []
    b_interface_list = []
    for residue_A in chain_a_res:
        for residue_B in chain_b_res:
            if residue_A.resname == "GLY":
                atom_A = "CA"
            else:
                atom_A = "CB"
            if residue_B.resname == "GLY":
                atom_B = "CA"
            else:
                atom_B = "CB"	
            distance = residue_A[atom_A] - residue_B[atom_B]
            if distance < 8.0:
				#print(f"A:{residue_A.resname}{residue_A.get_id()[1]},B:{residue_B.resname}{residue_B.get_id()[1]},{atom_B}, Distance:{distance}")
                a_interface_list.append([residue_A])
                b_interface_list.append([residue_B])
    a_short = []
    [a_short.append(x) for x in a_interface_list if x not in a_short]
    b_short = []
    [b_short.append(x) for x in b_interface_list if x not in b_short]
    return a_short, b_short

###### Input values ######
#Resolution
res = int(input("Provide an even number as the resolution: "))

#Structure names
clash_input = input("Provide the clash pdb: ")
clash_pqr_input = input("Provide the clash pqr: ")
template_input = input("Provide the template pdb: ")

#Load in structures
pdb_parser = PDBParser()
pqr_parser = PDBParser(is_pqr = True)
clash_pdb = pdb_parser.get_structure("Clash", clash_input)
clash_pqr = pqr_parser.get_structure("Clash_pqr", clash_pqr_input)
template_pdb = pdb_parser.get_structure("Template", template_input)

###### Intersect calculations ######

#Look up lists
clash_pdb_chains = list(clash_pdb.get_chains())

chain_a_atoms = list(clash_pdb_chains[0].get_atoms())
chain_b_atoms = list(clash_pdb_chains[1].get_atoms())

a_atom_coords = np.array([a.coord for a in chain_a_atoms], dtype = np.float64)
a_atom_radii = np.array([vdw_radii[a.element] for a in chain_a_atoms], dtype = np.float64)

b_atom_coords = np.array([a.coord for a in chain_b_atoms], dtype = np.float64)
b_atom_radii = np.array([vdw_radii[a.element] for a in chain_b_atoms], dtype = np.float64)

chain_a_res_vol = []
chain_b_res_vol = []
chain_a_res = list(clash_pdb_chains[0].get_residues())
chain_b_res = list(clash_pdb_chains[1].get_residues())

#Creating list of residue volumes
for i, residue in enumerate(chain_a_res):
    residue_atoms = list(residue.get_atoms())
    atom_coords = np.array([a.coord for a in residue_atoms], dtype = np.float64)
    atom_radii = np.array([vdw_radii[a.element] for a in residue_atoms], dtype = np.float64)
    residue_volume = res_vol(atom_radii, atom_coords, res)
    chain_a_res_vol.append([residue_volume, i])

for i, residue in enumerate(chain_b_res):
    residue_atoms = list(residue.get_atoms())
    atom_coords = np.array([a.coord for a in residue_atoms], dtype = np.float64)
    atom_radii = np.array([vdw_radii[a.element] for a in residue_atoms], dtype = np.float64)
    residue_volume = res_vol(atom_radii, atom_coords, res)
    chain_b_res_vol.append([residue_volume, i])

# Generating spheres for input atoms
spheres_a = []
spheres_b = []
for index, item in enumerate(a_atom_radii):
    sphere = compute_sphere(item, res)
    spheres_a.append(sphere)

for index, item in enumerate(b_atom_radii):
    sphere = compute_sphere(item, res)
    spheres_b.append(sphere)

spheres_b = np.array(spheres_b)

spheres_a = np.array(spheres_a)

#Moving spheres to correct coordinates

for index, item in enumerate(spheres_a):
    for key, value in enumerate(a_atom_coords):
        if index == key:
            for a, b in enumerate(item):
                item[a]= np.add(b, value)


for index, item in enumerate(spheres_b):
    for key, value in enumerate(b_atom_coords):
        if index == key:
            for a, b in enumerate(item):
                item[a]= np.add(b, value)

#Calculate the 'a' values for any chain A atoms that overlap with chain B atoms (and vice versa)
sphere_a_overlaps = compute_a_value(b_atom_coords, b_atom_radii, spheres_a, res, a_atom_coords)

sphere_b_overlaps = compute_a_value(a_atom_coords, a_atom_radii, spheres_b, res, b_atom_coords)

#Calculate cap volumes
a_cap_vol = []
for item in sphere_a_overlaps:
    for i, radius in enumerate(a_atom_radii):
        if item[0]==i:
            holding_cap_vol = (calc_cap_vol(item[1], radius))
            a_cap_vol.append([holding_cap_vol, item[0], item[2]])

b_cap_vol = []
for item in sphere_b_overlaps:
    for i, radius in enumerate(b_atom_radii):
        if item[0]==i:
            holding_cap_vol = (calc_cap_vol(item[1], radius))
            b_cap_vol.append([holding_cap_vol, item[0], item[2]])
            

#Getting arrays in same order
a_cap_order = []
b_cap_order = []
for i, A in enumerate(a_cap_vol):
    for j, B in enumerate(b_cap_vol):
        if A[1] == B[2] and A[2] == B[1]:
            a_cap_order.append(A)
            b_cap_order.append(B)

#Calculating intersect volumes
a_intersect_vol = []
b_intersect_vol = []
for i, A in enumerate(a_cap_order):
    for j, B in enumerate(b_cap_order):
        if i == j:
            vol_sum = A[0]+B[0]
            a_intersect_vol.append([vol_sum, A[1], A[2]])
            b_intersect_vol.append([vol_sum, B[1], B[2]])

#Calculating final percentages

a_clash_score = percent_clash(a_intersect_vol, chain_a_atoms, chain_a_res, chain_a_res_vol)
b_clash_score = percent_clash(b_intersect_vol, chain_b_atoms, chain_b_res, chain_b_res_vol)

#Storing clash_score as b-factor. Sets all b-factors to 0, then overwrites the relevant values with the clash percentage

all_res = list(clash_pdb.get_residues())
for residue in all_res:
    for atom in residue:
        atom.set_bfactor(0)
    for val in a_clash_score:
        if val[1] == residue:
            for atom in residue:
                atom.set_bfactor(val[0]*100)
    for val in b_clash_score:
        if val[1] == residue:
            for atom in residue:
                atom.set_bfactor(val[0]*100)

#Saving results

f = open('clash_scores.txt', 'w')
b_factors = []
for residue in all_res:
    for a, atom in enumerate(residue):
        if a == 0:
            b_factors.append(str(atom.get_bfactor()))

f.write(','.join(b_factors))
f.close()

io = PDBIO()
io.set_structure(clash_pdb)
io.save("steric_clashes.pdb")
print("Clashes saved")

###### End of intersect volume calculation ######

###### Charge clash calculation ######

#Look up lists
pqr_chains = list(clash_pqr.get_chains())

pqr_a_atoms = list(pqr_chains[0].get_atoms())
pqr_b_atoms = list(pqr_chains[1].get_atoms())

pqr_a_coords = np.array([a.coord for a in pqr_a_atoms], dtype = np.float64)
pqr_b_coords = np.array([a.coord for a in pqr_b_atoms], dtype = np.float64)

pqr_a_charges = np.array([a.get_charge() for a in pqr_a_atoms], dtype = np.float64)
pqr_b_charges = np.array([a.get_charge() for a in pqr_b_atoms], dtype = np.float64)

#Calculate the forces between chains

all_forces = calc_forces(pqr_a_coords, pqr_a_charges, pqr_b_coords, pqr_b_charges)

a_net = all_forces[0]
b_net = all_forces[1]

#Converting atom index to residue then summing forces by residue
a_net_res = convert_index_res(a_net, pqr_a_atoms)
b_net_res = convert_index_res(b_net, pqr_b_atoms)

a_sum = sum_forces(a_net_res)
b_sum = sum_forces(b_net_res)

#Storing clash_score as b-factor. Sets all b-factors to 0, then overwrites the relevant values with the clash percentage

all_res = list(clash_pdb.get_residues())
for residue in all_res:
    for atom in residue:
        atom.set_bfactor(0)
    for val in a_sum:
        if val[1] == residue:
            for atom in residue:
                atom.set_bfactor(val[0])
    for val in b_sum:
        if val[1] == residue:
            for atom in residue:
                atom.set_bfactor(val[0])

#Saving results
f = open('charge_forces.txt', 'w')
b_factors = []
for residue in all_res:
    for a, atom in enumerate(residue):
        if a == 0:
            b_factors.append(str(atom.get_bfactor()))

f.write(','.join(b_factors))
f.close()

io.set_structure(clash_pdb)
io.save("net_forces.pdb")
print("Charges saved")

###### End of charge clash calculation ######

###### Calculating missing interface contacts ######

#Finding interface residues
template_interfaces = get_interface(template_pdb)

template_ckbp = template_interfaces[0]

clash_interfaces = get_interface(clash_pdb)

clash_ckbp = clash_interfaces[0]
clash_chem = clash_interfaces[1]

#Updating b-factors with scores:
    #If in both template & clash, score = -2.0
    #If in clash only, score = -1.0
    #If in template but not clash, score = 1.0
    #Otherwise, score = 0.0

all_res = list(clash_pdb.get_residues())
for residue in all_res:
    for atom in residue:
        atom.set_bfactor(0)
    for val in template_ckbp:
        if val[0] == residue:
            for atom in residue:
                atom.set_bfactor(1)
    for val in clash_chem:
        if val[0] == residue:
            for atom in residue:
                atom.set_bfactor(-1)
    for val in clash_ckbp:
        if val[0] == residue:
            for atom in residue:
                atom.set_bfactor(-2)

#Saving file

io.set_structure(clash_pdb)
io.save("missing_contacts.pdb")
print("Missing contacts saved")