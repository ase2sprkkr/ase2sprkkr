import sys
import numpy as np
import scipy.constants
import pymatgen.core
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from copy import deepcopy
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial.transform import Rotation as R
from random import seed
from random import randint
from itertools import permutations
from datetime import datetime
import warnings
from matplotlib import rc
from matplotlib.backend_bases import MouseEvent
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import Axes3D, proj3d

print(' Using Matplotlib version {}'.format(matplotlib.__version__))


class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        super().draw(renderer)

    def do_3d_projection(self, renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        return np.min(zs)

rc('font',**{'family':'sans-serif', 'size':12})
#rc('text', usetex=True)


def on_pick(event):
    if isinstance(event, MouseEvent):
        return  # We want pick_event, not mouse press

    ind = event.ind[0]  # Index of picked point
    picked_point = highsymmpts[ind]
    #print(f"Clicked on point: {picked_point}")
    #annot.set_text(f"Picked: {picked_point}")
    fig_kkr.canvas.draw_idle()
    
    # Toggle: remove if last in path
    if selected_path and selected_path[-1] == ind:
        #print(f"Removed point {ind} from path")
        selected_path.pop()
    else:
        #print(f"Added point {ind} to path")
        selected_path.append(ind)

    # Update annotation
    #annot.set_text(f"Picked: {picked_point}")

    # Update highlight marker. Should this be (1) the last clicked point or (2) the last point on path? I think the latter...
    # Option (1)
    #highlight.set_data([picked_point[0]], [picked_point[1]])
    #highlight.set_3d_properties([picked_point[2]])
    #highlight.set_visible(True)
    # Option (2)
    if selected_path:
      highlight.set_data([highsymmpts[selected_path[-1]][0]], [highsymmpts[selected_path[-1]][1]])
      highlight.set_3d_properties([highsymmpts[selected_path[-1]][2]])
      highlight.set_visible(True)
    else:
#      highlight.set_data([highsymmpts[selected_path[-1]][0]], [highsymmpts[selected_path[-1]][1]])
#      highlight.set_3d_properties([highsymmpts[selected_path[-1]][2]])
      highlight.set_visible(False)

    # Update the path line
    update_path_line()

    fig_kkr.canvas.draw_idle()



def update_path_line():
  if selected_path:
    path_coords = highsymmpts[selected_path]
    path_line.set_data(path_coords[:, 0], path_coords[:, 1])
    path_line.set_3d_properties(path_coords[:, 2])
  else:
    path_line.set_data([], [])
    path_line.set_3d_properties([])



"""
Tool to analyze the symmetry of structure described in SPRKKR potential file and constructing a k-point path for band structure calculations. This version is a GUI based approach that allows the user to choose the k-path also for non-primitive unit cells.
As always, please make sure you understand what you are doing and use at your own risk.
The author does not take responsibility if your results are wrong :)
"""

def dot_norm(a,b):
    return np.dot(a,b)/np.sqrt(np.dot(a,a)*np.dot(b,b))




__author__ = 'AP'

a0 = scipy.constants.value('Bohr radius')*1e10

potfile = sys.argv[1]

ALAT = 0.0
cell = np.zeros((3,3))
rec_cell = np.zeros((3,3))


###########################################
#                                         #
# Read unit cell data from potential file #
#                                         #
###########################################

with open(potfile) as f:
    line = f.readline()
    Nlinemax = 10000
    Nline = 0
    while('NQ    ' not in line or Nline > Nlinemax):
        line = f.readline()
        Nline += 1
        if(Nline>Nlinemax):
            print('Could not find NQ.')
            sys.exit(1)
    NQ = int(line.split()[1])
    line = f.readline()
    NT = int(line.split()[1])
    line = f.readline()
    NM = int(line.split()[1])
    # Initialize arrays
    species_Z = np.zeros(NT)
    species_labels = []
    positions = np.zeros((NQ,3))
  
    
    Nline = 0
    while('ALAT  ' not in line or Nline > Nlinemax):
        line = f.readline()
        Nline += 1
        if(Nline>Nlinemax):
            print('Could not find ALAT.')
            sys.exit(1)
    ALAT = float(line.split()[1])
    for i in range(3):
        line = f.readline()
        cell[i,0] = float(line.split()[1])
        cell[i,1] = float(line.split()[2])
        cell[i,2] = float(line.split()[3])
    
    Nline = 0
    while('BASSCALE' not in line or Nline > Nlinemax):
        line = f.readline()
        Nline += 1
        if(Nline>Nlinemax):
            print('Could not find site positions.')
            sys.exit(1)
    line = f.readline()
    for i in range(NQ):
        line = f.readline()
        positions[i,0] = float(line.split()[1])
        positions[i,1] = float(line.split()[2])
        positions[i,2] = float(line.split()[3])

    IREFQ = np.zeros(NQ)
    Nline = 0
    while('OCCUPATION' not in line or Nline > Nlinemax):
        line = f.readline()
        Nline += 1
        if(Nline>Nlinemax):
            print('Could not find site occupations array.')
            sys.exit(1)
    line = f.readline()
    for i in range(NQ):
        line = f.readline()
        IREFQ[i] = int(line.split()[1])


    Nline = 0
    while('TYPES' not in line or Nline > Nlinemax):
        line = f.readline()
        Nline += 1
        if(Nline>Nlinemax):
            print('Could not find site information about atom types.')
            sys.exit(1)
    line = f.readline()
    for i in range(NT):
        line = f.readline()
        species_labels.append(line.split()[1])
        species_Z[i] = int(line.split()[2])


##############################
#                            #
# Convert to pymatgen format #
#                            #
##############################

# Create a pymatgen object. Pymatgen does not recognize empty cells so I change them to H.

pmg_cell = ALAT*a0*cell
pmg_coords = ALAT*a0*positions

if any(species_Z==0):
    print('\nNOTE: Empty cells are replaced by H!\n')
    pmg_Z = np.where(species_Z==0, species_Z+1, species_Z)
else:
    pmg_Z = species_Z

# Make sure that the species list has the same length as coordinates list
pmg_Z_full = np.zeros(np.shape(IREFQ))
for i in range(len(IREFQ)):
    pmg_Z_full[i] = pmg_Z[int(IREFQ[i]-1)]

pmg_structure = pymatgen.core.Structure(cell*2*np.pi, pmg_Z_full, positions*2*np.pi, coords_are_cartesian=True)
print(pymatgen.core.Structure(pmg_cell, pmg_Z_full, pmg_coords, coords_are_cartesian=True))
 

#sga = SpacegroupAnalyzer(pmg_structure, symprec=0.001)
#kpath_structure = sga.get_primitive_standard_structure(international_monoclinic=False,keep_site_properties=True)

#symmetries = sga.get_symmetry_operations()
#print('Detected symmetries:')
#for symmetry in symmetries:
#    print(symmetry)
    
#symmetry_dataset = sga.get_symmetry_dataset()
#print(symmetry_dataset)


#sg_info = kpath_structure.get_space_group_info()

#print('\nThe structure is {}'.format(sga.get_lattice_type()))
#print('Space group #{} ({})\n'.format(sg_info[1],sg_info[0]))

# Some problematic structures cause warnings from pymatgen and the results are nonsense. Therefore, stop execution if a warning is issued.
#with warnings.catch_warnings():
#    warnings.simplefilter('error')
#    kpath = HighSymmKpath(kpath_structure,symprec=0.1,path_type='setyawan_curtarolo',atol=1e-5,has_magmoms=True)


"""
################################################################
#                                                              #
# Checking the matching between database and SPRKKR structures #
#                                                              #
################################################################

tol = 1.e-4
supercell_found=False
if (np.abs(pmg_structure.lattice.volume-kpath_structure.lattice.volume) > tol):
    print('Unit cell is not primitive. Trying to construct k-path for corresponding primitive unit cell.')
    vol_ratio = pmg_structure.lattice.volume/kpath_structure.lattice.volume
    print('Volume ratio (sc/pc) is {:.3f}.'.format(pmg_structure.lattice.volume/kpath_structure.lattice.volume))
    #primitive_structure = sga.get_primitive_standard_structure()
    primitive_structure = sga.find_primitive() # This is a primitive version of the original structure, however the orientation can be different.
    #primitive_structure = pmg_structure.get_primitive_structure()
    
    print(primitive_structure.lattice.matrix)
    
    
    print(pmg_structure)
    print('Transformation matrix from primitive unit cell to supercell:') 
#    print(primitive_structure.lattice.matrix)
#    print(pmg_structure.lattice.matrix)
    # sc_matrix is transformation from primitive to supercell, so that SC = sc_matrix × PC
    # TODO: NOTE: This is not the correct supercell matrix, because the lattice vectors can be chosen differently. Or is it? Needs verification.
    sc_matrix = np.matmul(pmg_structure.lattice.matrix,primitive_structure.lattice.inv_matrix)
#    print(np.round(sc_matrix).astype(int))
    print(sc_matrix)
    print(np.linalg.det(sc_matrix))
    det_T = int(vol_ratio)
    
    SC_cell = pmg_structure.lattice.matrix
    PC_cell = primitive_structure.lattice.matrix
    
    # Construct supercell matrix
    a11 = np.linalg.norm(SC_cell[0,:])/np.linalg.norm(PC_cell[0,:])
    a22 = np.linalg.norm(SC_cell[1,:])/np.linalg.norm(PC_cell[1,:])
    a33 = np.linalg.norm(SC_cell[2,:])/np.linalg.norm(PC_cell[2,:])
    #a11 = np.linalg.norm(SC_cell[:,0])/np.linalg.norm(PC_cell[:,0])
    #a22 = np.linalg.norm(SC_cell[:,1])/np.linalg.norm(PC_cell[:,1])
    #a33 = np.linalg.norm(SC_cell[:,2])/np.linalg.norm(PC_cell[:,2])

    print(np.linalg.det(SC_cell)/np.linalg.det(PC_cell))
    
    sc_matrix = np.array([[a11,0,0],[0,a22,0,],[0,0,a33]])
    print(SC_cell)
    print(PC_cell)
    
    # From now on, work with the primitive structure. 
    # TODO: Need to make sure the primitive cell orientation corresponds to the supercell orientation!!!
    supercell_structure = deepcopy(pmg_structure)
    pmg_structure = deepcopy(primitive_structure)
    supercell_found = True
    #sys.exit(1)
else:
    print('Unit cell is primitive. Continue...')


find_transform = False
tol_mat = 1.e-4
if (not np.allclose(pmg_structure.lattice.matrix,kpath_structure.lattice.matrix,atol=tol_mat)):
    print('Reciprocal lattice vectors are not the same. Need to find transformation.')
    find_transform = True
else:
    print('Input structure has same reciprocal lattice vectors as the database structure.')

"""

BZ_pmg   =   pmg_structure.lattice.get_brillouin_zone()
#BZ_kpath = kpath_structure.lattice.get_brillouin_zone()


kx_BZ_pmg   = []
#kx_BZ_kpath = []
ky_BZ_pmg   = []
#ky_BZ_kpath = []
kz_BZ_pmg   = []
#kz_BZ_kpath = []


for face in BZ_pmg:
    for vertex in face: 
        kx_BZ_pmg.append(vertex[0])
        ky_BZ_pmg.append(vertex[1])
        kz_BZ_pmg.append(vertex[2])

#for face in BZ_kpath:
#    for vertex in face: 
#        kx_BZ_kpath.append(vertex[0])
#        ky_BZ_kpath.append(vertex[1])
#        kz_BZ_kpath.append(vertex[2])


vert_BZ_pmg = np.array([np.array(kx_BZ_pmg), np.array(ky_BZ_pmg), np.array(kz_BZ_pmg)]).T
#vert_BZ_kpath = np.array([np.array(kx_BZ_kpath), np.array(ky_BZ_kpath), np.array(kz_BZ_kpath)]).T

#c, c_i = np.unique(vert_BZ_kpath,axis=0,return_index=True)
d, d_i = np.unique(vert_BZ_pmg,axis=0,return_index=True)

#kpath_BZ_unique_verts = vert_BZ_kpath[c_i]
pmg_BZ_unique_verts = vert_BZ_pmg[d_i]

"""
# Go through all rows in kpath_BZ_unique_verts and compare each to rows in pmg_BZ_unique_verts. 
# Default: No rotation
min_rot = R.identity()
BZ_match = False

if(find_transform):
    matches = []
    for vert_kpath in kpath_BZ_unique_verts:
        for vert_pmg in pmg_BZ_unique_verts:
            if np.isclose(vert_kpath,vert_pmg,atol=1e-8).all(): 
               print('vertices match')
               matches.append(1)


    if len(matches) == len(kpath_BZ_unique_verts):
        print('Brillouin zones have the same orientation.')

        # In this case the reciprocal and real space lattice vectors can still be connected with 
        # a simple transform matrix (with integer elements).  
        M = np.matmul(kpath.structure.lattice.matrix,np.linalg.inv(pmg_structure.lattice.matrix))
        print('The transformation matrix connecting the real space lattice vectors is:')
        print(np.round(M,6))

    else: 
        Niter = 0
        Nitermax = 100
        while (not BZ_match and Niter < Nitermax):
            Niter+=1
        # In this case there is a rotation between the BZs, and we need to find it. The problem is that align_vectors() works pairwise, 
        # i.e. it tries to match the first row of a with the first row of b rather than finding a solution for any maches. 
        # Therefore the vertices should be ordered in so that the corresponding points are on same rows. That is difficult to achieve, 
        # because there can be Nvert > 20 individual vertices and if we consider all possible orders that we should try, 
        # the number of permutations is Nvert! (=factorial(Nvert)). 
        # Let's take 3 'smartly' chosen values from the first and try permutations of 3 points from the second one. 
        # As an example, if there are 24 vertices, you can choose the 3 points in 24!/(24-3)! = 12144 ways. 
        # When we have a match based on the rssd of align_vectors, we need to check the matching of all unique vertices. 
        # The point would be to find the rotation matrix that minimizes rssd (optimally it is zero).
        # Hopefully the matrix will then apply to the rest of the vertices as well. 
        # If not, we choose another 3 vertices to try until a match is found or max number of trials is reached. 
            Nvert = 3
            max_iter_rssd = 1000
            rssd_tol = 1e-12 
            min_rssd = 1e4
            iteration = 0
            while (min_rssd > rssd_tol or iteration > max_iter_rssd):
                iteration+=1
                indices = []
                j = 0
                # Try to be smart choosing the vertices and make sure that they are located in different parts of the BZ:
                # 1 with y, z coordinates positive and x negative
                # 1 with x, z positive and y negative
                # 1 with x, y positive and z negative
                # Additionally check that dot(a,b)/sqrt(dot(a,a)*dot(b,b)) is not close to -1 (antiparallel).
                # This can surely be done better, but the point is to avoid choosing e.g. antiparallel vectors.
                # After this I added the condition to try until a match is found, so this step is not really necessary anymore. 

                while j< Nvert:
                    r = randint(0,len(kpath_BZ_unique_verts)-1)
#                    if r not in indices:
#                        if j==0 and kpath_BZ_unique_verts[r][0] < 0.0 and kpath_BZ_unique_verts[r][1] > 0.0 and kpath_BZ_unique_verts[r][2] > 0.0:
#                            indices.append(r)
#                            j+=1
#                        elif j==1 and kpath_BZ_unique_verts[r][0] > 0.0 and kpath_BZ_unique_verts[r][1] < 0.0 and kpath_BZ_unique_verts[r][2] > 0.0 and dot_norm(kpath_BZ_unique_verts[indices[0]],kpath_BZ_unique_verts[r]) > -0.9:
#                            indices.append(r)
#                            j+=1
#                        elif j==2 and kpath_BZ_unique_verts[r][0] > 0.0 and kpath_BZ_unique_verts[r][1] > 0.0 and kpath_BZ_unique_verts[r][2] < 0.0 and dot_norm(kpath_BZ_unique_verts[indices[0]],kpath_BZ_unique_verts[r]) > -0.9 and dot_norm(kpath_BZ_unique_verts[indices[1]],kpath_BZ_unique_verts[r]) > -0.9:
#                            indices.append(r)
#                            j+=1
# Remove this if above is uncommented
                    indices.append(r)
                    j+=1



                       
                print('Chosen indices: ', indices)
                A = kpath_BZ_unique_verts[indices,:]

                print('Chosen vectors:')
                print(A)
    
                # Loop over possible permutations of Nvert numbers 
                print('Looping over permutations of indices')
                perm = permutations(range(len(pmg_BZ_unique_verts)),Nvert)
                B = np.zeros((len(indices),3))
                for p in perm:
                    i = 0
                    for index in p:
                        B[i,:] = pmg_BZ_unique_verts[index,:]
                        i+=1

                    rot, rssd, sens = R.align_vectors(A, B, return_sensitivity=True)
                    if (rssd < min_rssd):
                        print('Better match found with rssd {:.4E} for indices {}'.format(np.round(rssd,6), p))
                        min_rssd = rssd
                        min_rot = rot


            # Apply min_rot to the full list of vertices and see if the BZs match
            print('Best rotation found:')
            print(np.round(min_rot.as_matrix(),6))
            print('Testing it for all BZ vertices.')


            matches = []
            for vert_kpath in kpath_BZ_unique_verts:
                for vert_pmg in pmg_BZ_unique_verts:
                    if np.isclose(vert_kpath,min_rot.apply(vert_pmg),atol=1e-10).all(): 
                        matches.append(1)

            print(len(matches), 'vertices match')
            if len(matches) == len(kpath_BZ_unique_verts):
                print('Brillouin zones have the same orientation.')
                BZ_match = True
            else: 
                print('There is still mismatch between BZ orientations.', len(kpath_BZ_unique_verts)-len(matches), 'vertices do not match.')
                print('Cannot continue with this choice.')
                print('Trying another set of vectors.')
                BZ_match = False
            

        M = np.matmul(kpath.structure.lattice.reciprocal_lattice.matrix,np.linalg.inv(min_rot.apply(pmg_structure.lattice.reciprocal_lattice.matrix)))
        print('The transformation matrix connecting the reciprocal space lattice vectors is:')
        print(np.round(M,6))



# Now we have the transformation min_rot for rotating the BZ to the 'standard' orientation
# and for transforming the reciprocal lattice vectors (M).
# We need to apply them to the standardized k-path to have them in the KKR basis.

"""

#############################
#                           #
# Visualize Brillouin zones #
#                           #
#############################

# Matplotlib axes trickery copied online to have the same scale along x, y, and z.
def set_axes_equal(ax: plt.Axes):
    #Set 3D plot axes to equal scale.
    #
    #Make axes of 3D plot have equal scale so that spheres appear as
    #spheres and cubes as cubes.  Required since `ax.axis('equal')`
    #and `ax.set_aspect('equal')` don't work on 3D.
    limits = np.array([
        ax.get_xlim3d(),
        ax.get_ylim3d(),
        ax.get_zlim3d(),
    ])
    origin = np.mean(limits, axis=1)
    radius = 0.5 * np.max(np.abs(limits[:, 1] - limits[:, 0]))
    _set_axes_radius(ax, origin, radius)

def _set_axes_radius(ax, origin, radius):
    x, y, z = origin
    ax.set_xlim3d([x - radius, x + radius])
    ax.set_ylim3d([y - radius, y + radius])
    ax.set_zlim3d([z - radius, z + radius])

# End of trickery copied online



def set_equal_lims(ax):
    # Get current axis limits
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    zlim = ax.get_zlim()

    # Get the max range and center
    all_lims = np.array([xlim, ylim, zlim])
    centers = np.mean(all_lims, axis=1)
    max_range = 2*np.max(all_lims[:,1] - all_lims[:,0])

    # Set all axes to the same range around the center
    for ctr, setlim in zip(centers, [ax.set_xlim, ax.set_ylim, ax.set_zlim]):
        setlim(ctr - max_range/2, ctr + max_range/2)
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_zlim(-1, 1)
    ax.set_box_aspect([1, 1, 1])


#fig_std = plt.figure(figsize=(6,6))
#ax_std = fig_std.add_subplot(111, projection='3d')
#ax_std.set_box_aspect([1.0, 1.0, 1.0])

fig_kkr = plt.figure(figsize=(8,6))
ax_kkr = fig_kkr.add_subplot(111, projection='3d')

ax_kkr.set_axis_off()
ax_kkr.set_position([0, 0, 1, 1])
#ax_kkr.dist = 20
set_equal_lims(ax_kkr)



# Draw KKR structure BZ
BZ = pmg_structure.lattice.get_brillouin_zone()
for face in BZ:
    x = np.zeros(len(face))
    y = np.zeros(len(face))
    z = np.zeros(len(face))
    i=0
    for vertex in face: 
        x[i]=vertex[0]
        y[i]=vertex[1]
        z[i]=vertex[2]
        i+=1

    pc = Poly3DCollection([list(zip(x,y,z))],facecolor=[0.,0.,1.,.05],edgecolor=[0.,0.,0.,.8],linewidth=1)
    ax_kkr.add_collection3d(pc)

# Draw reciprocal space vectors of KKR structure
draw_reciprocal_basis = True
draw_k_basis = True

rec_basis = pmg_structure.lattice.reciprocal_lattice.matrix
if draw_reciprocal_basis:
  arrow = Arrow3D([0,rec_basis[0,0]],[0,rec_basis[0,1]],[0,rec_basis[0,2]], mutation_scale=20, lw=.5, arrowstyle="-|>", color="blue")
  ax_kkr.add_artist(arrow)
  ax_kkr.text(rec_basis[0,0] + 0.02, rec_basis[0,1] + 0.02, rec_basis[0,2] + 0.02, r'$b_1$', fontsize=10)
  
  arrow = Arrow3D([0,rec_basis[1,0]],[0,rec_basis[1,1]],[0,rec_basis[1,2]], mutation_scale=20, lw=.5, arrowstyle="-|>", color="blue")
  ax_kkr.add_artist(arrow)
  ax_kkr.text(rec_basis[1,0] + 0.02, rec_basis[1,1] + 0.02, rec_basis[1,2] + 0.02, r'$b_2$', fontsize=10)
  
  arrow = Arrow3D([0,rec_basis[2,0]],[0,rec_basis[2,1]],[0,rec_basis[2,2]], mutation_scale=20, lw=.5, arrowstyle="-|>", color="blue")
  ax_kkr.add_artist(arrow)
  ax_kkr.text(rec_basis[2,0] + 0.02, rec_basis[2,1] + 0.02, rec_basis[2,2] + 0.02, r'$b_3$', fontsize=10)

# Draw k-space unit vectors
if draw_k_basis:
  length=0.4
  arrow = Arrow3D([0,length],[0,0],[0,0], mutation_scale=20, lw=1, arrowstyle="-|>", color="red")
  ax_kkr.add_artist(arrow)
  ax_kkr.text(length + 0.02, 0 + 0.02, 0 + 0.02, r'$k_x$', fontsize=10)
  
  arrow = Arrow3D([0,0],[0,length],[0,0], mutation_scale=20, lw=1, arrowstyle="-|>", color="green")
  ax_kkr.add_artist(arrow)
  ax_kkr.text(0 + 0.02, length + 0.02, 0 + 0.02, r'$k_y$', fontsize=10)
  
  arrow = Arrow3D([0,0],[0,0],[0,length], mutation_scale=20, lw=1, arrowstyle="-|>", color="blue")
  ax_kkr.add_artist(arrow)
  ax_kkr.text(0 + 0.02, 0 + 0.02, length + 0.02, r'$k_z$', fontsize=10)





# Collect possibly interesting high-symmetry points of the BZ that is in correct orientation. 
# These are:
# 1. Gamma point
# 2. Vertices of the BZ
# 3. Centers of BZ faces
# 4. Midpoints of BZ edges.
# Make sure that each point is included only once. 

#print(BZ)

highsymmpnt_size = 15

sz_gamma = 20
sz_vertex = 20
sz_face = 20
sz_edge = 20

c_gamma = [0,0,0]
c_vertex = [0,0,0]
c_face = [0,0,0]
c_edge = [0,0,0]

sz = []
colors = []

# 1. Gamma point
#ax_kkr.scatter(0,0,0,s=highsymmpnt_size,color=highsymmvertex_color)
highsymmpts = np.array([[0,0,0]])
sz.append(sz_gamma)
colors.append(c_gamma)

# 2. Vertices of BZ
#ax_kkr.scatter(pmg_BZ_unique_verts[:,0],pmg_BZ_unique_verts[:,1],pmg_BZ_unique_verts[:,2],s=highsymmpnt_size,color=highsymmvertex_color)
highsymmpts = np.append(highsymmpts,pmg_BZ_unique_verts,axis=0)
for i in range(len(pmg_BZ_unique_verts)):
  sz.append(sz_vertex)
  colors.append(c_vertex)

# 3. Centers of BZ faces
for face in BZ:
  x = 0
  y = 0
  z = 0
  nvert = 0
  for vertex in face: 
    x+=vertex[0]
    y+=vertex[1]
    z+=vertex[2]
    nvert+=1
  x=x/nvert
  y=y/nvert
  z=z/nvert
  #ax_kkr.scatter(x,y,z,s=highsymmpnt_size,color=highsymmface_color)
  highsymmpts = np.append(highsymmpts,np.array([[x,y,z]]),axis=0)
  sz.append(sz_face)
  colors.append(c_face)
  
# 4. Midpoints of BZ edges.
# Need to make sure every point is added only once because the midpoints belong to several faces.
edge_midpoints_list = []
for face in BZ:
  coords = np.zeros((len(face)+1,3))
  i=0
  for vertex in face: 
    coords[i,0]=vertex[0]
    coords[i,1]=vertex[1]
    coords[i,2]=vertex[2]
    i+=1
  # Add the first coordinate as last to have every edge covered
  coords[i,0]=face[0][0]
  coords[i,1]=face[0][1]
  coords[i,2]=face[0][2]
  
  # Go through each list of cooridnates and add midpoints of edges to the list
  for i in range(len(coords[:,0])-1):
    x_mid = (coords[i,0]+coords[i+1,0])/2.
    y_mid = (coords[i,1]+coords[i+1,1])/2.
    z_mid = (coords[i,2]+coords[i+1,2])/2.
    edge_midpoints_list.append([x_mid,y_mid,z_mid])
    


edge_midpoints_list = np.array(edge_midpoints_list)

d, d_i = np.unique(edge_midpoints_list,axis=0,return_index=True)

edge_midpoints = edge_midpoints_list[d_i]
  
highsymmpts = np.append(highsymmpts,edge_midpoints,axis=0)

for i in range(len(edge_midpoints)):
  sz.append(sz_edge)
  colors.append(c_edge)
sz = np.array(sz)
#print(sz)
#print(highsymmpts)

#ax_kkr.scatter(edge_midpoints[:,0],edge_midpoints[:,1],edge_midpoints[:,2],s=highsymmpnt_size,color=highsymmedge_color)


# Scatter plot highsymmpts

ax_kkr.scatter(highsymmpts[:,0],highsymmpts[:,1],highsymmpts[:,2],s=sz,color=colors,picker=True)


#kpath_original = deepcopy(kpath)

# Add high-symmetry k-points to the BZ.
#for kpt in kpath.kpath['kpoints']:
#    k_kkr = kpath.kpath['kpoints'][kpt]
#    # Convert from the direct coordinates to (kx,ky,kz)
#    kpath.kpath['kpoints'][kpt] = np.matmul(k_kkr,kpath_structure.lattice.reciprocal_lattice.matrix)
#    k_kkr = kpath.kpath['kpoints'][kpt]
#    ax_std.scatter(k_kkr[0],k_kkr[1],k_kkr[2],color='red')
#    ax_std.text(k_kkr[0]+0.05,k_kkr[1],k_kkr[2], '$'+kpt+'$', horizontalalignment='left',verticalalignment='center',color='black')

"""
# Draw the path
i = 1
for segment in kpath.kpath['path']:
    kx = []
    ky = []
    kz = []
    for highsymmpoint in segment:
        kx.append(kpath.kpath['kpoints'][highsymmpoint][0])
        ky.append(kpath.kpath['kpoints'][highsymmpoint][1])
        kz.append(kpath.kpath['kpoints'][highsymmpoint][2])
    ax_std.plot(kx,ky,kz,color='red',lw=2)
"""

"""
# Transform and add high-symmetry k-points to the BZ of the KKR orientation. 
# In the previous step, they were already transformed to the kpath basis, now
# it is only needed to rotate them.
kpath_kkr = deepcopy(kpath)
for kpt in kpath_kkr.kpath['kpoints']:
    k_kkr = kpath_kkr.kpath['kpoints'][kpt]
    
    kpath_kkr.kpath['kpoints'][kpt] = min_rot.apply(kpath_kkr.kpath['kpoints'][kpt],inverse=True)

    k_kkr = kpath_kkr.kpath['kpoints'][kpt]
    ax_kkr.scatter(k_kkr[0],k_kkr[1],k_kkr[2],color='red')
    ax_kkr.text(k_kkr[0]+0.05,k_kkr[1],k_kkr[2], '$'+kpt+'$', horizontalalignment='left',verticalalignment='center',color='black')


# Draw the path
i = 1
for segment in kpath_kkr.kpath['path']:
    kx = []
    ky = []
    kz = []
    for highsymmpoint in segment:
        kx.append(kpath_kkr.kpath['kpoints'][highsymmpoint][0])
        ky.append(kpath_kkr.kpath['kpoints'][highsymmpoint][1])
        kz.append(kpath_kkr.kpath['kpoints'][highsymmpoint][2])
    ax_kkr.plot(kx,ky,kz,color='red',lw=2)
"""







"""
if supercell_found:
    BZ_sc = supercell_structure.lattice.get_brillouin_zone()
    for face in BZ_sc:
        x = np.zeros(len(face))
        y = np.zeros(len(face))
        z = np.zeros(len(face))
        i=0
        for vertex in face:
            x[i]=vertex[0]
            y[i]=vertex[1]
            z[i]=vertex[2]
            i+=1

        pc = Poly3DCollection([list(zip(x,y,z))],facecolor=[0.,0.,1.,.1],edgecolor=[0.,0.,0.,.8],linewidth=1)
        ax_kkr.add_collection3d(pc)

    # Draw reciprocal space vectors of KKR structure
    rec_basis = pmg_structure.lattice.reciprocal_lattice.matrix
    ax_kkr.plot([0,rec_basis[0,0]],[0,rec_basis[0,1]],[0,rec_basis[0,2]],lw=2,color='blue')
    ax_kkr.plot([0,rec_basis[1,0]],[0,rec_basis[1,1]],[0,rec_basis[1,2]],lw=2,color='green')
    ax_kkr.plot([0,rec_basis[2,0]],[0,rec_basis[2,1]],[0,rec_basis[2,2]],lw=2,color='black')

    kpath_kkr = deepcopy(kpath)
    for kpt in kpath_kkr.kpath['kpoints']:
        k_kkr = kpath_kkr.kpath['kpoints'][kpt]

        #kpath_kkr.kpath['kpoints'][kpt] = min_rot.apply(kpath_kkr.kpath['kpoints'][kpt],inverse=True)
        # Convert from the direct coordinates to (kx,ky,kz)
        kpath_kkr.kpath['kpoints'][kpt] = np.matmul(k_kkr,np.linalg.inv(sc_matrix))

        #kpath_kkr.kpath['kpoints'][kpt] = np.matmul(np.linalg.inv(sc_matrix),k_kkr)

        k_kkr = kpath_kkr.kpath['kpoints'][kpt]
        ax_kkr.scatter(k_kkr[0],k_kkr[1],k_kkr[2],color=[0,1,0])
        ax_kkr.text(k_kkr[0]+0.05,k_kkr[1],k_kkr[2], '$'+kpt+'$', horizontalalignment='left',verticalalignment='center',color='black')


    # Draw the path
    i = 1
    for segment in kpath_kkr.kpath['path']:
        kx = []
        ky = []
        kz = []
        for highsymmpoint in segment:
            kx.append(kpath_kkr.kpath['kpoints'][highsymmpoint][0])
            ky.append(kpath_kkr.kpath['kpoints'][highsymmpoint][1])
            kz.append(kpath_kkr.kpath['kpoints'][highsymmpoint][2])
        ax_kkr.plot(kx,ky,kz,color=[0,1,0],lw=2)


"""














"""
# Draw standard structure BZ
BZ = kpath_structure.lattice.get_brillouin_zone()
for face in BZ:
    x = np.zeros(len(face))
    y = np.zeros(len(face))
    z = np.zeros(len(face))
    i=0
    for vertex in face:
        x[i]=vertex[0]
        y[i]=vertex[1]
        z[i]=vertex[2]
        i+=1

    pc = Poly3DCollection([list(zip(x,y,z))],facecolor=[0.,1.,0.,.1],edgecolor=[0.,0.,0.,.8],linewidth=1)
    ax_std.add_collection3d(pc)

# Draw reciprocal space vectors of standard structure
rec_basis = kpath_structure.lattice.reciprocal_lattice.matrix


ax_std.plot([0,rec_basis[0,0]],[0,rec_basis[0,1]],[0,rec_basis[0,2]],lw=1.5,color='blue')
ax_std.plot([0,rec_basis[1,0]],[0,rec_basis[1,1]],[0,rec_basis[1,2]],lw=1.5,color='green')
ax_std.plot([0,rec_basis[2,0]],[0,rec_basis[2,1]],[0,rec_basis[2,2]],lw=1.5,color='black')
"""




"""
# Convert the k-point coordinates into 1/Angstrom and print
print('\nCoordinates of high-symmetry k-points in 1/Angstrom:\n')
for kpt in kpath_kkr.kpath['kpoints']:
    factor = (2*np.pi)/(ALAT*a0)
    print('{}:\t ({:9.6f},{:9.6f},{:9.6f} )'.format(kpt,factor*kpath_kkr.kpath['kpoints'][kpt][0],factor*kpath_kkr.kpath['kpoints'][kpt][1],factor*kpath_kkr.kpath['kpoints'][kpt][2]))
"""





"""
# Convert the k-point coordinates into SPRKKR format and print
print('\nCoordinates of high-symmetry k-points in 2*pi/ALAT for BLOCHSF calculations:\n')
for kpt in kpath_kkr.kpath['kpoints']:
    print('{}:\t {{{:9.6f},{:9.6f},{:9.6f} }}'.format(kpt,kpath_kkr.kpath['kpoints'][kpt][0],kpath_kkr.kpath['kpoints'][kpt][1],kpath_kkr.kpath['kpoints'][kpt][2]))





# Construct and print a string representing the high-symmetry k-path
highsymmstring=''

i = 1
for segment in kpath_kkr.kpath['path']:
    j = 1
    for highsymmpoint in segment:
        highsymmstring+=highsymmpoint
        if j < len(segment): 
            highsymmstring+='-'
            j+=1

    if i < len(kpath_kkr.kpath['path']): 
        highsymmstring+=', '
    i+=1

print('\nHigh-symmetry k-point path:')
print(highsymmstring+'\n')





################
#              #
# Print k-path #
#              # 
################

j = 0
k = 1
print('')

for segment in kpath_kkr.kpath['path']:
    for i in range(0,len(segment)-1):
        KA = kpath_kkr.kpath['path'][j][i]
        KE = kpath_kkr.kpath['path'][j][i+1]
        str_KA = 'KA'+str(k)+'='
        str_KE = 'KE'+str(k)+'='
        #print('    {}{{{:9.6f},{:9.6f},{:9.6f} }}  {}{{{:9.6f},{:9.6f},{:9.6f} }}'.format(str_KA,kpath_kkr.kpath['kpoints'][KA][0],kpath_kkr.kpath['kpoints'][KA][1],kpath_kkr.kpath['kpoints'][KA][2],str_KE,kpath_kkr.kpath['kpoints'][KE][0],kpath_kkr.kpath['kpoints'][KE][1],kpath_kkr.kpath['kpoints'][KE][2]))
        print('    {}{{{},{},{}}} {}{{{},{},{}}}'.format(str_KA,np.round(kpath_kkr.kpath['kpoints'][KA][0],5),np.round(kpath_kkr.kpath['kpoints'][KA][1],5),np.round(kpath_kkr.kpath['kpoints'][KA][2],5),str_KE,np.round(kpath_kkr.kpath['kpoints'][KE][0],5),np.round(kpath_kkr.kpath['kpoints'][KE][1],5),np.round(kpath_kkr.kpath['kpoints'][KE][2],5)))
        k+=1
    j+=1
"""






#ax_std.set_title('BZ of standard primitive cell')

#ax_kkr.set_title('BZ of structure from potential file')


#set_axes_equal(ax_std)
#ax_std.set_proj_type('ortho')

#set_axes_equal(ax_kkr)
ax_kkr.set_proj_type('ortho')

# Marker for highlighting the selected point
highlight, = ax_kkr.plot([0], [0], [0], 'o', color='red', markersize=8, visible=False)

# Line to show selected path
path_line, = ax_kkr.plot([], [], [], '-', color='red', linewidth=3)
path_scatter = ax_kkr.scatter([], [], [], s=50, color='red')
# Store path as a list of indices
selected_path = []

fig_kkr.canvas.mpl_connect('pick_event', on_pick)

plt.tight_layout()
plt.show()

#print(highsymmpts[selected_path])


################
#              #
# Print k-path #
#              # 
################


kpath_length = 0.0
# Find length of k-path in 1/angstrom
if len(highsymmpts[selected_path]) > 1:
  for i in range(len(highsymmpts[selected_path])-1):
    kpath_length += np.sqrt((highsymmpts[selected_path][i+1][0] - highsymmpts[selected_path][i][0])**2 + (highsymmpts[selected_path][i+1][1] - highsymmpts[selected_path][i][1])**2 + (highsymmpts[selected_path][i+1][2] - highsymmpts[selected_path][i][2])**2)
    
kpath_length *= (2.0*np.pi)/(ALAT*a0) # In 1/angstrom

K_resolution = 0.01  # I set this as default resolution of k-sampling, in 1/angstrom. It is just a guiding value and only affects the NK printed out.

NK = int(np.ceil(kpath_length/K_resolution))
print('')
print('K-path length is {:.3f} 1/A.'.format(kpath_length))
print('The suggested {} k-points along the path corresponds to k-resolution of {} 1/Angstrom.'.format(NK,K_resolution))

NKDIR = len(highsymmpts[selected_path])-1


j = 0
k = 1
print('')

if len(highsymmpts[selected_path]) > 1:
  print('    TASK BSF')
  print('    NKDIR={:2d}  NK={:4d}'.format(NKDIR,NK))
  for i in range(len(highsymmpts[selected_path])-1):
    KA = highsymmpts[selected_path][i]
    KE = highsymmpts[selected_path][i+1]
    str_KA = 'KA'+str(k)+'='
    str_KE = 'KE'+str(k)+'='
    print('    {}{{{},{},{}}} {}{{{},{},{}}}'.format(str_KA,np.round(KA[0],5),np.round(KA[1],5),np.round(KA[2],5),str_KE,np.round(KE[0],5),np.round(KE[1],5),np.round(KE[2],5)))
    k+=1
else:
  print('No k-path selected')






