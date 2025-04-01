import sys
import numpy as np
import re
import matplotlib.pyplot as plt

class calculate_E_formation:
    def __init__(self, filein):
        self.filein = filein
        self.fileout = filein + ".out"
        self.debug = True#False

        self.VBM = 0.0
        self.band_gap = 0.0
        self.host_supercell_energy = 0.0
        self.host_type = 0  #1 - mono-atomic, 2-binary, 3-ternary crystal, etc.
        self.host_atom = []
        self.vacancies = []
        self.impurities = []
        self.charge_start = 0
        self.charge_end = 0
        self.chem_potentials = {} #dictionary type
        self.supercell_energy_list = []
        self.correction_terms_list = []

        # Read the input data
        self.read_data()

        self.n_charges = int(self.charge_end - self.charge_start + 1)
        self.supercell_energy = np.array(self.supercell_energy_list).reshape(self.n_charges, 2)
        if self.debug == True:
            print(f"self.supercell_energy: {self.supercell_energy}")
        self.correction_terms = np.array(self.correction_terms_list).reshape(self.n_charges, 3)
        if self.debug == True:
            print(f"self.correction_terms: {self.correction_terms}")

        # Write the formation energy for each charge state
        self.energy_step = 0.0001 # eV
        self.energy_part = 0.0  # eV
        self.charge_dependent_part = 0.0 # eV
        self.fermi_level = []   #x-axis Fermi-level
        for i in range(0,int(self.band_gap/self.energy_step)):
            self.fermi_level.append(i * self.energy_step)
        self.all_data = np.array([])  # will stack actual the data array in self.write_data().
        self.formation_energy_data = [] #the lowest energy value for a given fermi level

        self.write_data()

        self.CTL_list = []   #charge transition levels
        self.CTL_list = self.find_CTL() #read the self.formation_energy_data point-by-point, see where the slope changes.
        self.CTL_count = int(len(self.CTL_list)/2)
        self.CTL_array = np.array(self.CTL_list).reshape(self.CTL_count,2)



        # Plot the formation energy
        self.fig = self.plot_data()

    def read_data(self):
        with open(self.filein, 'r') as data_file:
            for line in data_file:
                if "&VBM" in line:
                    tmp_line =  data_file.readline()
                    self.VBM = float(tmp_line)
                    if self.debug == True:
                        print(f"VBM={self.VBM} is read.")

                elif "&band_gap" in line:
                    tmp_line = data_file.readline()
                    self.band_gap = float(tmp_line)
                    if self.debug == True:
                        print(f"band gap={self.band_gap} is read.")

                elif "&Host_type" in line:
                    tmp_line = data_file.readline()
                    re_findall = re.findall(r"[\w]+", tmp_line)
                    self.host_type = len(re_findall)
                    for i in range(0, len(re_findall)):
                        self.host_atom.append(re_findall[i])
                    if self.debug == True:
                        print(f"host_type={self.host_type} is read.")
                    if self.debug == True:
                        print(f"host_atom={self.host_atom} is read.")

                elif  "&Vacancies" in line:
                    tmp_line = data_file.readline()
                    re_findall = re.findall(r"[\w]+", tmp_line)
                    for i in range(0, len(re_findall)):
                        self.vacancies.append(re_findall[i])
                    if self.debug == True:
                        print(f"vacancies={self.vacancies} is read.")

                elif "&Impurities" in line:
                    tmp_line = data_file.readline()
                    re_findall = re.findall(r"[\w]+", tmp_line)
                    for i in range(0, len(re_findall)):
                        self.impurities.append(re_findall[i])
                    if self.debug == True:
                        print(f"impurities={self.impurities} is read.")

                elif "&Chemical_potentials" in line:
                    chk = True
                    while chk == True:
                        tmp_line = data_file.readline()
                        tmp_line_split = tmp_line.split()
                        if len(tmp_line_split) > 1:
                            self.chem_potentials[tmp_line_split[0]] = float(tmp_line_split[1])
                            if self.debug == True:
                                print(f"chem_potentials={self.chem_potentials} is read.")
                        elif len(tmp_line_split) == 0:  # hit the end of the chem_potential data.
                            chk = False
                        elif "&" in tmp_line:
                            print( "Warning: there should be an empty line between every section. \n")
                            chk = False
                        else:
                            continue

                elif "&Host_supercell_energy" in line:
                    tmp_line = data_file.readline()
                    self.host_supercell_energy = float(tmp_line)
                    if self.debug == True:
                        print(f"host_supercell_energy={self.host_supercell_energy} is read.")

                elif "&Charge_state_range" in line:
                    tmp_line = data_file.readline()
                    self.charge_start = int(tmp_line.split()[0])
                    self.charge_end = int(tmp_line.split()[1])
                    if self.debug == True:
                        print(f"charge_start and charge_end is read.")

                elif "&Defective_supercell_energy" in line:
                    for i in range(self.charge_start, self.charge_end+1):
                        tmp_line = data_file.readline()
                        if i == int(tmp_line.split()[0]):
                            self.supercell_energy_list.append(int(tmp_line.split()[0])) #charge state
                            if self.debug == True:
                                print(f"CS:supercell_energy_list={self.supercell_energy_list} is read.")
                            self.supercell_energy_list.append(float(tmp_line.split()[1])) #supercell E
                            if self.debug == True:
                                print(f"SE:host_supercell_energy={self.supercell_energy_list} is read.")
                        else:
                            sys.exit("check the charge state range and defective supercell energy\n")
                            
                elif "&Correction_terms" in line:
                    for i in range(self.charge_start, self.charge_end+1):
                        tmp_line = data_file.readline()
                        self.correction_terms_list.append(int(i))   #charge state
                        if self.debug == True:
                            print(f"CS:correction_terms_list={self.correction_terms_list} is read.")
                        self.correction_terms_list.append(float(tmp_line.split()[0])) #short-range
                        if self.debug == True:
                            print(f"SR:correction_terms_list={self.correction_terms_list} is read.")
                        self.correction_terms_list.append(float(tmp_line.split()[-1])) #Makov-Payne or any E_corr
                        if self.debug == True:
                            print(f"LR:correction_terms_list={self.correction_terms_list} is read.")

                else:
                    continue

    def write_data(self):
        Ry = 13.605692  #eV
        for i in range(0, self.n_charges):
            q = self.charge_start + i   # e.g. -2, -1, 0, 1, 2
            print ("\n" + "charge_state =" + str(q) + "\n")

            if q == self.correction_terms[i,0]: # Should contain the charge state, q.
                E_correction = self.correction_terms[i,2] # Makov-Payne or other E_corr
                print ("E_correction = " + str(E_correction) + "\n")
            else:
                sys.exit("Check the self.correction_terms np array. The first column should contain the charge state.\n")

            if q == self.correction_terms[i,0]: # Should contain the charge state, q.
                DeltaV_pot_align = self.correction_terms[i,1] # potential_alignment
                print ("DeltaV_pot_align = " + str(DeltaV_pot_align) + "\n")
            else:
                sys.exit("Check the self.correction_terms np array. The first column should contain the charge state.\n")

            if q == self.supercell_energy[i,0]:
                Energy_part = self.supercell_energy[i,1] *Ry + E_correction - 1*self.host_supercell_energy*Ry + q*DeltaV_pot_align
                print ("Energy_part = self.supercell_energy + E_correction - self.host_supercell_energy + q*DeltaV_pot_align = " + str(Energy_part) + "\n")
            else:
                sys.exit("Check the self.supercell_energy np array.\n")

            formation_energy =  Energy_part + q* self.VBM
            print ("Formation energy = Energy_part + q* self.VBM = " + str(formation_energy) + "\n")
            #add the chemical potential parts.

            if self.vacancies:
                for c in self.vacancies:
                    formation_energy = formation_energy - (-1)* self.chem_potentials[c]
                    print ("chem_potential for " + str(c) + " (vacancy)  has been added. \n")
                    print ("Energy_part + q* VBM + chem_pot (so far) = " + str(formation_energy) + "\n")
            else:
                print('no vacancies are present\n')

            if self.impurities:
                for c in self.impurities:
                    formation_energy = formation_energy - (+1)* self.chem_potentials[c]
                    print ("chem_potential for " + str(c) + " (impurity)  has been added. \n")
                    print ("Energy_part + q* VBM + chem_pot (so far) = " + str(formation_energy) + "\n")
            else:
                print('no foreign impurities are present \n')

            data_list =[]

            for j in range(0,int(self.band_gap/self.energy_step)):
                data_list.append(formation_energy + q* self.fermi_level[j])

            if len(self.all_data) == 0:
                if self.debug == True:
                    print(f"self.all_data=0, = {self.all_data}")
                self.all_data = np.append(self.all_data, np.array(data_list))
                if self.debug == True:
                    print(f"self.all_data=0a, = {self.all_data}")
            else:
                if self.debug == True:
                    print(f"self.all_data not 0, = {self.all_data}")
                self.all_data = np.vstack((self.all_data, np.array(data_list)))
                if self.debug == True:
                    print(f"self.all_data not 0a, = {self.all_data}")

        for i in range(0,int(self.band_gap/self.energy_step)):
            tmp_array = self.all_data[:,i]
            minimum_E = np.sort(tmp_array)[0]
            self.formation_energy_data.append(minimum_E)

        with open(self.fileout, 'w') as fout:
            for i in range(0,int(self.band_gap/self.energy_step)):
                fout.write('{0:<15.7f} {1:<15.7f} \n'.format(self.fermi_level[i], self.formation_energy_data[i]))

    def find_CTL(self):
        #find the charge transition levels where the slope changes
        CTL = []
        slope_to_compare = self.charge_end
        slope_threshold = 0.0001
        delta_x = self.energy_step

        # fout_str = str(self.filein) + "_CTL"
        # with open(fout_str, 'w') as fout:
        #     fout.write("charge transition levels for " + str(self.filein) + "\n")

        for i in range(1, int(self.band_gap/self.energy_step)):
            delta_y = self.formation_energy_data[i]-self.formation_energy_data[i-1]
            slope = delta_y/delta_x
            #print("slope = " + str(slope) + "\n")
            slope_diff = slope_to_compare - slope
            if abs(slope_diff) < slope_threshold:
                continue
            else:
                print("CTL from " + str(slope_to_compare) + " to " + str (slope_to_compare - 1) + " = " + \
                    str((self.fermi_level[i-1] + self.fermi_level[i])/2.0) + "\n")
                CTL_point = self.fermi_level[i-1] + delta_x/2.0
                CTL.append(CTL_point)
                CTL_energy = self.formation_energy_data[i-1]+slope*delta_x/2.0
                CTL.append(CTL_energy)
                # with open(fout_str, 'a') as fout:
                #     fout.write('{0:<3d} to  {1:<3d} = {2:<10.5f} \n'.format(slope_to_compare,\
                #         slope_to_compare-1, CTL_point))
                slope_to_compare = slope_to_compare - 1
        return CTL

    def plot_data(self):
        x_data = self.fermi_level
        x_min = self.fermi_level[0]
        x_max = self.fermi_level[int(self.band_gap/self.energy_step)-1]+self.energy_step

        fig = plt.figure(figsize=(10,8))
        axes = fig.add_subplot(1,1,1)
        axes.tick_params(axis='both', which='major', labelsize=15)
        axes.tick_params(axis='both', which='minor', labelsize=12)
        plt.xlim (x_min, x_max)
        #plt.ylim(-2,3)

        #plot all data
        for i in range(0, self.n_charges):
            q = self.charge_start + i   # e.g. -2, -1, 0, 1, 2
            y_data = self.all_data[i,:]
            label_str = fr'q={q}' #$Cu_{{Zn}}$ with
            axes.plot(x_data, y_data, '-', label=label_str)

        #plot the lowest-energy data on top
        y_data = self.formation_energy_data
        axes.plot(x_data,y_data, color = "k"  , linewidth = 3, linestyle = '-')

        axes.legend(loc=0)

        axes.set_xlabel('Fermi level (eV)', fontsize = 18)
        axes.set_ylabel('Defect formation energy (eV)', fontsize = 18)
        #axes.set_title('Co-vacancy in ZnS, relaxed', fontsize = 18)

        x_CTL_data = self.CTL_array[:,0]
        y_CTL_data = self.CTL_array[:,1]
        axes.scatter(x_CTL_data, y_CTL_data, color = 'k', marker = 'o', s = 50)

        #fig.savefig(str(self.fileout)+'.png')
        plt.legend(prop={'size': 18})
        # plt.show(fig)
        return fig
