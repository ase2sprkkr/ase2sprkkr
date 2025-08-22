#!/usr/bin/env python3

import numpy as np
import glob
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm as cm
import argparse
import os
import sys

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

class PotentialFileReader:
    def __init__(self, file_name):
        self.file_name = file_name
        self.compound = ""
        self.num_nq = 0
        self.num_nt = 0
        self.iq = np.array([])
        self.pos = np.array([])
        self.icl = np.array([])
        self.itoq = np.array([])
        self.conc = np.array([])
        self.data_labels = np.array([])
        self.labels = []
        self.trim_labels = []
        self.spin_mom = np.array([])

    def file_length(self, fname):
        with open(fname) as f:
            for i, l in enumerate(f):
                pass
        return i + 1

    def read_potential_file(self):
        if not os.path.exists(self.file_name):
            raise FileNotFoundError(f"Potential file not found: {self.file_name}")
            
        file_size = self.file_length(self.file_name)
        pot_file = open(self.file_name)    
        
        # First pass: get basic dimensions
        for ln in range(file_size):
            line = pot_file.readline() 
            data = line.split()
            if not data:
                continue
                
            if len(data) > 1 and data[0] == 'NQ':
                self.num_nq = int(data[1])
            elif len(data) > 1 and data[0] == 'NT':
                self.num_nt = int(data[1])
        
        # Reset and read all data
        pot_file.seek(0)
        for ln in range(file_size):
            line = pot_file.readline() 
            data = line.split()
            if not data:
                continue
                
            if len(data) > 1 and data[0] == 'SYSTEM':
                self.compound = self.process_system_info(data)
            elif len(data) > 1 and data[0] == 'NQ':
                # Already got this, just skip
                pass
            elif len(data) > 1 and data[1] == 'QBAS(X)':
                self.iq, self.pos = self.process_base_vec(pot_file, data, self.num_nq, ln)
            elif len(data) > 1 and data[1] == 'IREFQ':
                self.icl, self.itoq, self.conc = self.process_irefq_info(pot_file, data, self.num_nq)
            elif len(data) > 1 and data[1] == 'TXT_T':
                self.data_labels, self.labels, self.trim_labels = self.process_txt_t_info(pot_file, data, self.num_nt)
            elif len(data) > 1 and data[0] == 'MOMENTS':
                self.spin_mom = self.process_moments_info(pot_file, data, ln, file_size)
                
        pot_file.close()
        
        return (self.compound, self.num_nq, self.num_nt, self.iq, self.pos, 
                self.icl, self.itoq, self.conc, self.data_labels, self.labels, 
                self.trim_labels, self.spin_mom)
    

    def process_system_info(self, data):
        compound = data[1]
        print(f"*************Creating input files of {compound} for UppASD simulations*************")
        return compound
        
    def process_base_vec(self, pot_file, data, num_nq, ln):
        iq = np.zeros(num_nq, dtype=np.int64)
        pos = np.zeros((num_nq, 3), dtype=np.float64)
        
        for i in range(num_nq):
            line = pot_file.readline()            
            base_data = line.split()
            if not base_data or len(base_data) < 4:
                continue
            iq[i] = int(base_data[0])
            pos[i, 0] = float(base_data[1])
            pos[i, 1] = float(base_data[2])
            pos[i, 2] = float(base_data[3])
        return iq, pos

    def process_irefq_info(self, pot_file, data, num_nq):
        icl = np.zeros(num_nq, dtype=np.int64)
        itoq = np.zeros(num_nq, dtype=np.int64)
        conc = np.zeros(num_nq, dtype=np.float64)
        
        for i in range(num_nq):
            line = pot_file.readline()
            base_data = line.split()
            if len(base_data) < 6:
                continue
            icl[i] = int(base_data[1])
            itoq[i] = int(base_data[4])
            conc[i] = float(base_data[5])
        return icl, itoq, conc

    def process_txt_t_info(self, pot_file, data, num_nt):
        data_labels = np.empty(num_nt, dtype='object')
        labels = []
        trim_labels = []
        
        for i in range(num_nt):
            line = pot_file.readline()
            txt_t_data = line.split()
            if len(txt_t_data) < 2:
                continue
            data_labels[i] = str(txt_t_data[1])
            labels.append(f'${data_labels[i]}$')
            ind = str(data_labels[i]).find('_')
            if ind > 0:
                trim_labels.append(str(data_labels[i])[:ind])
            else:
                trim_labels.append(str(data_labels[i]))
        return data_labels, labels, trim_labels
        
    def process_moments_info(self, pot_file, data, ln, file_size):
        spin_mom = []
        reading_moments = False
        
        # Continue from current position
        for _ in range(ln, file_size):
            line = pot_file.readline()
            if not line:
                break
                
            data = line.split()
            if not data:
                continue
                
            if data[0] == 'MOMENTS':
                reading_moments = True
                continue
                
            if reading_moments and data[0] == 'TYPE':
                # Read the next line which contains the moment data
                line = pot_file.readline()
                data = line.split()
                if data and len(data) >= 3:
                    try:
                        spin_mom.append(float(data[2]))
                    except ValueError:
                        pass
        
        return np.array(spin_mom)


class JXCProcessor:
    def __init__(self, maptype=1):
        self.maptype = maptype
        self.jfile_format = "{:5.0f}  {:5.0f}  {: 4.10f}  {: 4.10f}  {: 4.10f}  {: 4.8f}  {:4.10f}\n"
        self.dmfile_format = "{:5.0f}  {:5.0f}  {: 4.10f}  {: 4.10f}  {: 4.10f}  {: 4.8f}  {: 4.8f}  {: 4.8f}  {:4.10f}\n"
        
        # Data attributes
        self.itype = np.array([])
        self.isite = np.array([])
        self.jtype = np.array([])
        self.jsite = np.array([])
        self.bond_type1 = np.array([])
        self.bond_type2 = np.array([])
        self.mod_ij = np.array([])
        self.Jij_meV = np.array([])
        self.Dij_x = np.array([])
        self.Dij_y = np.array([])
        self.Dij_z = np.array([])
            
    def file_length(self, fname):
        with open(fname) as f:
            for i, l in enumerate(f):
                pass
        return i + 1
            
    def process_JXC_file(self, file_path):
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"JXC file not found: {file_path}")
            
        JXCFile = open(file_path)
        file_size = self.file_length(file_path)
        
        for lines in range(file_size):
            line = JXCFile.readline()
            data = line.split()
            if len(data) > 1 and data[0] == 'IT':
                data = pd.read_csv(JXCFile, skiprows=0, header=None, sep='\s+').values
                non_zero_ind = np.where(data[:, 10] > 0)
                
                self.bond_type1 = np.zeros([len(non_zero_ind[0]), 3], dtype=np.float64)
                self.bond_type1[:, 0] = np.asarray(data[non_zero_ind[0], 7], dtype=np.float64)
                self.bond_type1[:, 1] = np.asarray(data[non_zero_ind[0], 8], dtype=np.float64)
                self.bond_type1[:, 2] = np.asarray(data[non_zero_ind[0], 9], dtype=np.float64)
                
                self.bond_type2 = np.zeros([len(non_zero_ind[0]), 3], dtype=np.int64)
                self.bond_type2[:, 0] = np.asarray(data[non_zero_ind[0], 4], dtype=np.int64)
                self.bond_type2[:, 1] = np.asarray(data[non_zero_ind[0], 5], dtype=np.int64)
                self.bond_type2[:, 2] = np.asarray(data[non_zero_ind[0], 6], dtype=np.int64)
                
                self.itype = np.asarray(data[non_zero_ind[0], 0], dtype=np.int64)
                self.isite = np.asarray(data[non_zero_ind[0], 1], dtype=np.int64)
                self.jtype = np.asarray(data[non_zero_ind[0], 2], dtype=np.int64)
                self.jsite = np.asarray(data[non_zero_ind[0], 3], dtype=np.int64)
                self.mod_ij = np.asarray(data[non_zero_ind[0], 10], dtype=np.float64)
                self.Jij_meV = np.asarray(data[non_zero_ind[0], 11], dtype=np.float64)
        
        JXCFile.close()

    def process_DMI_file(self, file_path):
        """Process Dzyaloshinskii-Moriya interaction file"""
        if not os.path.exists(file_path):
            print(f"No DMI file found: {file_path}")
            return False
            
        DMIFile = open(file_path)
        file_size = self.file_length(file_path)
        
        for lines in range(file_size):
            line = DMIFile.readline()
            data = line.split()
            if len(data) > 1 and data[0] == 'IT':
                data = pd.read_csv(DMIFile, skiprows=0, header=None, sep='\s+').values
                non_zero_ind = np.where(data[:, 10] > 0)
                
                self.Dij_x = np.asarray(data[non_zero_ind[0], 11], dtype=np.float64)
                self.Dij_y = np.asarray(data[non_zero_ind[0], 12], dtype=np.float64)
                self.Dij_z = np.asarray(data[non_zero_ind[0], 13], dtype=np.float64)
                DMIFile.close()
                return True
        
        DMIFile.close()
        return False

    def write_jfile(self, jfile_name, exclude_list, trim_labels):
        if len(self.itype) == 0:
            print("No JXC data to write")
            return False
            
        jfile = open(jfile_name, 'w')
        for ii in range(len(self.itype)):
            if (trim_labels[self.itype[ii] - 1] not in exclude_list) and (
                    trim_labels[self.jtype[ii] - 1] not in exclude_list):
                if self.maptype == 1:
                    jfile.write(self.jfile_format.format(self.isite[ii], self.jsite[ii], 
                                                         self.bond_type1[ii, 0], self.bond_type1[ii, 1],
                                                         self.bond_type1[ii, 2], 
                                                         self.Jij_meV[ii], self.mod_ij[ii]))
                else:
                    jfile.write(self.jfile_format.format(self.isite[ii], self.jsite[ii], 
                                                         self.bond_type2[ii, 0], self.bond_type2[ii, 1],
                                                         self.bond_type2[ii, 2], 
                                                         self.Jij_meV[ii], self.mod_ij[ii]))
        jfile.close()
        return True

    def write_dmfile(self, dmfile_name, exclude_list, trim_labels):
        if len(self.Dij_x) == 0:
            print("No DMI data to write")
            return False
            
        dmfile = open(dmfile_name, 'w')
        for ii in range(len(self.itype)):
            if (trim_labels[self.itype[ii] - 1] not in exclude_list) and (
                    trim_labels[self.jtype[ii] - 1] not in exclude_list):
                if self.maptype == 1:
                    dmfile.write(self.dmfile_format.format(self.isite[ii], self.jsite[ii],
                                                           self.bond_type1[ii, 0], self.bond_type1[ii, 1],
                                                           self.bond_type1[ii, 2],
                                                           self.Dij_x[ii], self.Dij_y[ii], self.Dij_z[ii],
                                                           self.mod_ij[ii]))
                else:
                    dmfile.write(self.dmfile_format.format(self.isite[ii], self.jsite[ii],
                                                           self.bond_type2[ii, 0], self.bond_type2[ii, 1],
                                                           self.bond_type2[ii, 2],
                                                           self.Dij_x[ii], self.Dij_y[ii], self.Dij_z[ii],
                                                           self.mod_ij[ii]))
        dmfile.close()
        return True


class Plotter:
    def __init__(self, font_size=28):
        self.font_size = font_size
        
    def plot_exchange_interactions(self, processor, num_nt, trim_labels, spin_mom, exclude_list, 
                                 exchange_radius=4.0, tol=0.01, output_dir="."):
        """Plot exchange interactions as function of distance"""
        if len(processor.Jij_meV) == 0:
            print("No exchange data to plot")
            return False
            
        # Apply distance cutoff
        ind_cut = np.where(processor.mod_ij <= exchange_radius)
        mod_ij = processor.mod_ij[ind_cut[0]]
        Jij_meV = processor.Jij_meV[ind_cut[0]] / 13.6  # Convert to mRy
        itype = processor.itype[ind_cut[0]] - 1  # Convert to 0-based indexing
        jtype = processor.jtype[ind_cut[0]] - 1
        bond_type1 = processor.bond_type1[ind_cut[0], :]
        
        colors = cm.Paired(np.linspace(0, 1, 2 * num_nt + 2))
        
        plots_created = 0
        for ii in range(num_nt):
            if trim_labels[ii] not in exclude_list:
                i_ind = np.where(itype == ii)
                
                if len(i_ind[0]) == 0:
                    continue
                    
                fig = plt.figure()
                counter = 0
                
                for jj in range(num_nt):
                    j_ind = np.where(jtype[i_ind[0]] == jj)
                    
                    if len(j_ind[0]) == 0:
                        continue
                        
                    if len(spin_mom) > ii and len(spin_mom) > jj:
                        sign_i = np.sign(spin_mom[ii])
                        sign_j = np.sign(spin_mom[jj])
                        sign_ij = sign_i * sign_j
                    else:
                        sign_ij = 1.0
                        
                    if trim_labels[jj] not in exclude_list:
                        counter += 1
                        plt.plot(mod_ij[i_ind[0][j_ind[0]]], sign_ij * Jij_meV[i_ind[0][j_ind[0]]],
                                alpha=0.75, lw=4, c=colors[counter], 
                                label=f'${trim_labels[ii]}$-${trim_labels[jj]}$')
                        plt.scatter(mod_ij[i_ind[0][j_ind[0]]], sign_ij * Jij_meV[i_ind[0][j_ind[0]]],
                                  color=colors[counter], alpha=0.75, s=300, lw=1.00, edgecolor='black')
                
                # Plotting options
                if len(Jij_meV[i_ind[0]]) > 0:
                    extremum = max(abs(np.max(Jij_meV[i_ind[0]])), abs(np.min(Jij_meV[i_ind[0]]))) * 1.10
                    
                    plt.legend(fontsize=self.font_size, loc='upper right')
                    plt.ylabel(r'$J_{ij}$ [mRy]', fontsize=self.font_size)
                    plt.xlabel(r'$r_{ij}/a_{lat}$', fontsize=self.font_size)
                    
                    ax = plt.gca()
                    ax.set_facecolor((1, 1, 1))
                    ax.tick_params(axis='x', colors='black', labelsize=self.font_size, width=2)
                    ax.tick_params(axis='y', colors='black', labelsize=self.font_size, width=2)
                    ax.set_ylim(-extremum, extremum)
                    plt.axhline(0, color='black', linestyle='--')
                    
                    for axis in ['top', 'bottom', 'left', 'right']:
                        ax.spines[axis].set_linewidth(3)
                    
                    plt.grid(False)
                    fig.set_size_inches(18.5, 10.5)
                    output_path = os.path.join(output_dir, f'Jij_{trim_labels[ii]}.pdf')
                    plt.savefig(output_path, transparent=False, dpi=300, bbox_inches='tight')
                    plt.close(fig)
                    plots_created += 1
                    print(f"Created plot: {output_path}")
        
        return plots_created > 0


def write_pos_file(trim_labels, exclude_list, num_nt, labels, num_nq, icl, iq, pos, conc, output_file='posfile.dat'):
    posfile_format = "{:5.0f} {:5.0f} {: 4.10f} {: 4.10f} {: 4.10f}\n"
    
    posfile = open(output_file, 'w')
    atoms_written = 0
    
    for i in range(num_nq):
        if trim_labels[icl[i] - 1] not in exclude_list:
            posfile.write(posfile_format.format(iq[i], icl[i], pos[i, 0], pos[i, 1], pos[i, 2]))
            atoms_written += 1
            
    posfile.close()
    print(f"Written {atoms_written} atoms to {output_file}")
    return atoms_written > 0


def write_mom_file(trim_labels, exclude_list, num_nt, spin_mom, icl, iq, output_file='momfile.dat'):
    if len(spin_mom) == 0:
        print("Warning: spin_mom array is empty. No MOMENTS information to write.")
        return False
        
    momfile_format = "{:5.0f} {:5.0f} {: 4.10f} {: 4.10f} {: 4.10f} {: 4.10f}\n"
    
    momfile = open(output_file, 'w')
    moments_written = 0
    
    for i in range(num_nt):
        if i < len(trim_labels) and trim_labels[i] not in exclude_list:
            ind = np.where(icl == (i + 1))
            for j in range(len(ind[0])):
                momfile.write(momfile_format.format(iq[ind[0][j]], 1, spin_mom[i], 0, 0, 1))
                moments_written += 1
                
    momfile.close()
    print(f"Written {moments_written} moments to {output_file}")
    return moments_written > 0


def parse_arguments():
    parser = argparse.ArgumentParser(description='Convert SPRKKR output to UppASD input files and plots the computed Jij')
    
    # File arguments
    parser.add_argument('-p', '--pot-file', type=str, default=None, 
                       help='Input potential file (.pot_new)')
    parser.add_argument('-j', '--jxc-file', type=str, default=None,
                       help='Input exchange interaction file (_XCPLTEN_Jij.dat)')
    parser.add_argument('-d', '--dmi-file', type=str, default=None,
                       help='Input DMI file (_DMIVEC_Dij.dat)')
    
    # Output control
    parser.add_argument('-o', '--output-dir', type=str, default='.',
                       help='Output directory for files')
    parser.add_argument('-n', '--no-write', action='store_true',
                       help='Skip writing output files')
    
    # Plotting options
    parser.add_argument('--plot', action='store_true',
                       help='Generate exchange interaction plots')
    parser.add_argument('--no-plot', action='store_true',
                       help='Skip generating plots')
    parser.add_argument('-r', '--exchange-radius', type=float, default=4.0,
                       help='Maximum distance for exchange interactions')
    
    # Processing options
    parser.add_argument('-m', '--maptype', type=int, default=2, choices=[1, 2],
                       help='Mapping type for coordinates (1: Cartesian, 2: Lattice)')
    parser.add_argument('-e', '--exclude', type=str, nargs='+', default=['Vc'],
                       help='Element types to exclude')
    parser.add_argument('-f', '--font-size', type=int, default=28,
                       help='Font size for plots')
    
    return parser.parse_args()


def find_files_by_pattern(args):
    """Find files by pattern if not explicitly provided"""
    if args.pot_file is None:
        pot_files = glob.glob("*.pot_new")
        if pot_files:
            args.pot_file = pot_files[0]
    
    if args.jxc_file is None:
        jxc_files = glob.glob("*_XCPLTEN_Jij.dat")
        if jxc_files:
            args.jxc_file = jxc_files[0]
    
    if args.dmi_file is None:
        dmi_files = glob.glob("*_DMIVEC_Dij.dat")
        if dmi_files:
            args.dmi_file = dmi_files[0]
    
    return args


def main():
    args = parse_arguments()
    args = find_files_by_pattern(args)
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Read potential file
    if args.pot_file is None:
        print("Error: No potential file found!")
        sys.exit(1)
    
    try:
        reader = PotentialFileReader(args.pot_file)
        compound, num_nq, num_nt, iq, pos, icl, itoq, conc, data_labels, labels, trim_labels, spin_mom = reader.read_potential_file()
        
        print(f"Compound: {compound}")
        print(f"Number of atoms: {num_nq}")
        print(f"Number of types: {num_nt}")
        print(f"Spin moments: {spin_mom}")
        
        # Write position and moment files unless disabled
        if not args.no_write:
            pos_file = os.path.join(args.output_dir, 'posfile.dat')
            mom_file = os.path.join(args.output_dir, 'momfile.dat')
            
            write_pos_file(trim_labels, args.exclude, num_nt, labels, num_nq, icl, iq, pos, conc, pos_file)
            write_mom_file(trim_labels, args.exclude, num_nt, spin_mom, icl, iq, mom_file)
        
        # Process exchange interactions if file provided
        if args.jxc_file:
            processor = JXCProcessor(maptype=args.maptype)
            processor.process_JXC_file(args.jxc_file)
            
            # Write exchange file unless disabled
            if not args.no_write:
                j_file = os.path.join(args.output_dir, 'jfile.dat')
                processor.write_jfile(j_file, args.exclude, trim_labels)
            
            # Process DMI if available
            if args.dmi_file:
                has_dmi = processor.process_DMI_file(args.dmi_file)
                if has_dmi and not args.no_write:
                    dm_file = os.path.join(args.output_dir, 'dmfile.dat')
                    processor.write_dmfile(dm_file, args.exclude, trim_labels)
            
            # Plot exchange interactions if requested
            should_plot = args.plot or (not args.no_plot and args.plot is None)
            if should_plot:
                plotter = Plotter(font_size=args.font_size)
                plotter.plot_exchange_interactions(processor, num_nt, trim_labels, spin_mom, 
                                                 args.exclude, args.exchange_radius, args.output_dir)
        else:
            print("No JXC file found - skipping exchange processing")
            
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
