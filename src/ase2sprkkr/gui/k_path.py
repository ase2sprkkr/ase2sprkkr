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
from ase import Atoms

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


"""
Tool to analyze the symmetry of structure described in SPRKKR potential file and constructing a k-point path for band structure calculations. This version is a GUI based approach that allows the user to choose the k-path also for non-primitive unit cells.
As always, please make sure you understand what you are doing and use at your own risk.
The author does not take responsibility if your results are wrong :)
"""

def dot_norm(a,b):
    return np.dot(a,b)/np.sqrt(np.dot(a,a)*np.dot(b,b))


def k_path_gui(atoms:Atoms, verbose=False):

    def on_pick(event):
        if isinstance(event, MouseEvent):
            return  # We want pick_event, not mouse press

        ind = event.ind[0]  # Index of picked point
        picked_point = highsymmpts[ind]
        fig_kkr.canvas.draw_idle()

        # Toggle: remove if last in path
        if selected_path and selected_path[-1] == ind:
            #print(f"Removed point {ind} from path")
            selected_path.pop()
        else:
            #print(f"Added point {ind} to path")
            selected_path.append(ind)

        if selected_path:
          highlight.set_data([highsymmpts[selected_path[-1]][0]], [highsymmpts[selected_path[-1]][1]])
          highlight.set_3d_properties([highsymmpts[selected_path[-1]][2]])
          highlight.set_visible(True)
        else:
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



    # Create a pymatgen object. Pymatgen does not recognize empty cells so I change them to H.

    species_Z = atoms.get_atomic_numbers()
    cell = atoms.cell
    positions = atoms.positions
    ALAT = cell.get_bravais_lattice().a
    cell = cell.copy() / ALAT

    pmg_Z_full = atoms.get_atomic_numbers()
    pmg_Z_full[pmg_Z_full == 0] = 1

    #pymatgen does not support vacuum
    pmg_structure = pymatgen.core.Structure(cell*2*np.pi, pmg_Z_full, positions*2*np.pi, coords_are_cartesian=True)
    if verbose:
        a0 = 1e10
        pmg_cell = a0*cell
        pmg_coords = a0*positions
        print(pymatgen.core.Structure(pmg_cell, pmg_Z_full, pmg_coords, coords_are_cartesian=True))

    BZ_pmg   =   pmg_structure.lattice.get_brillouin_zone()

    kx_BZ_pmg   = []
    ky_BZ_pmg   = []
    kz_BZ_pmg   = []

    for face in BZ_pmg:
        for vertex in face:
            kx_BZ_pmg.append(vertex[0])
            ky_BZ_pmg.append(vertex[1])
            kz_BZ_pmg.append(vertex[2])

    vert_BZ_pmg = np.array([np.array(kx_BZ_pmg), np.array(ky_BZ_pmg), np.array(kz_BZ_pmg)]).T
    d, d_i = np.unique(vert_BZ_pmg,axis=0,return_index=True)
    pmg_BZ_unique_verts = vert_BZ_pmg[d_i]

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


    fig_kkr = plt.figure(figsize=(8,6))
    ax_kkr = fig_kkr.add_subplot(111, projection='3d')

    ax_kkr.set_axis_off()
    ax_kkr.set_position([0, 0, 1, 1])
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

    # Scatter plot highsymmpts
    ax_kkr.scatter(highsymmpts[:,0],highsymmpts[:,1],highsymmpts[:,2],s=sz,color=colors,picker=True)
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

    kpath_length *= (2.0*np.pi)/ALAT # In 1/angstrom

    K_resolution = 0.01  # I set this as default resolution of k-sampling, in 1/angstrom. It is just a guiding value and only affects the NK printed out.

    if len(selected_path):

        NK = int(np.ceil(kpath_length/K_resolution))
        if verbose:
            print('')
            print('K-path length is {:.3f} 1/A.'.format(kpath_length))
            print('The suggested {} k-points along the path corresponds to k-resolution of {} 1/Angstrom.'.format(NK,K_resolution))

        path_points = highsymmpts[selected_path]
        return {
            'NKDIR': len(path_points)-1,
            'NK': NK,
            'KA' : path_points[:-1],
            'KE' : path_points[1:],
        }
    else:
        return None


def k_path_to_string(k_path:dict):
    out = []
    if k_path:
      out.append('    TASK BSF')
      out.append('    NKDIR={:2d}  NK={:4d}'.format(k_path['NKDIR'],k_path['NK']))
      k=0
      for KA, KE in zip(k_path['KA'], k_path['KE']):
        k+=1
        str_KA = 'KA'+str(k)+'='
        str_KE = 'KE'+str(k)+'='
        out.append('    {}{{{},{},{}}} {}{{{},{},{}}}'.format(str_KA,np.round(KA[0],5),np.round(KA[1],5),np.round(KA[2],5),str_KE,np.round(KE[0],5),np.round(KE[1],5),np.round(KE[2],5)))
    return '\n'.join(out)
