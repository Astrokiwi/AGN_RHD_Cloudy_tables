import numpy as np
import glob
import parameter_ranges as pr
from astropy import constants
from astropy import units
import sys
import os
import pandas as pd

import cProfile, pstats, io
from pstats import SortKey

def logn_to_rho(logn):
    """logn_to_rho(logn)

    simple conversion from hydrogen density (log base 10 of cm**-2) to mass density (g cm**-3, not logged)

    Assumes that n=1 cm**-3 means rho=2.35E-24 g cm**-3, which is more or less typical for an ISM mix
    """
    return 2.35e-24*10**logn


def ryd2microns(x):
    return (1./(constants.Ryd*x)).to(units.um).value

mean_molecular_mass = 1.2

class SpectrumInterpolator:

    def __init__(self,n,i,t):
        self.n = n
        self.i = i
        self.t = t
        self.gen_depth_freq_lookup()

    def gen_depth_freq_lookup(self):
        """"SpectrumInterpolator.gen_depth_freq_lookup()

        Finds emissivity spectra files in table directory for the particular density (n), intensity (i), and gas temperature (t)
        define during construction
        Then builds a list of all column densities and wavelengths used in these files.

        Then, you can call `SpectrumInterpolator.get_spectrum(depth)` to get the spectrum at that column density (depth)
        Wavelengths are stored in SpectrumInterpolator.wavelengths

        You can also call `SpectrumInterpolator.get_all_spectra()` to get the full matrix of all spectra at all wavelength
        """
        # get list of spectrum files
        data_dir = pr.format_dir(self.n,self.i,self.t)
        spec_files = glob.glob(f'{data_dir}*znu*')
        self.file_base = data_dir+pr.format_file(self.n,self.i,self.t)+".znu"

        freq = []

        self.nfiles = len(spec_files)

        # check filenames are formatted correctly, extract frequency of each file
        for i in range(self.nfiles):
            file = f"{self.file_base}{i}"
            if file not in spec_files:
                raise FileNotFoundError(f"{file} not in .znu* files - filenames may not be contiguous?")
            with open(f"{file}") as f:
                header_line = f.readline()
            freq.append(float(header_line[-18:-4]))
        self.wavelengths = ryd2microns(np.array(freq))[::-1]

        # extract depths and distances from one file (should all be identical!)
        file = f"{self.file_base}0"
        d = np.loadtxt(file,usecols=[0,1],skiprows=1)

        self.radii  = d[:,0]
        self.depths = d[:,1]

    def get_spectrum(self,depth):
        """SpectrumInterpolator.get_spectrum(depth)

        Generates the spectrum at a given column
        This reads one line from each of ~10k files, so is slow for generating multiple spectra
        """
        #TODO: proper interpolation. For now just using lower bound
        idepth = np.searchsorted(self.depths,depth)
        emissivity = []

        # for i in range(self.nfiles) :
        for i in range(self.nfiles-1, -1, -1): # invert order to get increasing wavelength, not increasing frequency
            file = f"{self.file_base}{i}"
            emissivity.append(
                np.loadtxt(file,usecols=2,skiprows=1)[idepth] # load in emmisivity column, select correct depth
                )

        return np.array(emissivity)

    def gather_all_spectra(self,chunksize=32):
        """SpectrumInterpolator.gather_all_spectra()

        Combines all dumped emissivities and opacities (scattering and absorbed) into 3 big matrices, for quick spectrum plotting
        May use a lot of ram
        Returns matrices, also saves them as SpectrumInterpolator.emissivity_matrix,SpectrumInterpolator.opac_abs_matrix,SpectrumInterpolator.opac_scat_matrix
        """
        #TODO: proper interpolation. For now just using lower bound
        spec_temp_file = "workspace/tempspec.tmp"

        print("concatenating files")

        for nchars in range(1,4+1):
            wildcard = "?"*(nchars)
            if nchars==1:
                pipe = ">"
            else:
                pipe = ">>"
            #TODO: replace os.system with newer preferred version
            os.system(f"tail -q -n +2 {self.file_base}{wildcard} | cut -f 3,4,5 {pipe} {spec_temp_file}")

        d = pd.read_csv(spec_temp_file,header=None,sep='\t').values

        # reshape, then flip
        ndepth = d.shape[0]//self.nfiles
        self.emissivity_matrix = np.flip(d[:,0].reshape((self.nfiles,ndepth)).T,axis=1)
        self.opac_abs_matrix = np.flip(d[:,1].reshape((self.nfiles,ndepth)).T,axis=1)
        self.opac_scat_matrix = np.flip(d[:,2].reshape((self.nfiles,ndepth)).T,axis=1)

        return self.emissivity_matrix,self.opac_abs_matrix,self.opac_scat_matrix

    def calculate_depths(self,recalculate=False):
        if recalculate :
            self.gather_all_spectra()
        else:
            try:
                x=self.opac_abs_matrix
            except:
                self.gather_all_spectra()
        dr = np.gradient(self.radii)
        self.opac_abs_tau = np.cumsum(self.opac_abs_matrix*dr[:,np.newaxis],axis=0)
        self.opac_scat_tau = np.cumsum(self.opac_scat_matrix*dr[:,np.newaxis],axis=0)
        print(self.opac_abs_tau.shape)

    def read_incident_transmitted_spectrum(self):
        data_dir = pr.format_dir(self.n,self.i,self.t)
        confilename = data_dir+pr.format_file(self.n,self.i,self.t)+".con" # putain de con

        d = np.loadtxt(confilename,usecols = [0,1,2])

        self.con_wavelengths = np.flip(ryd2microns(d[:,0])) # TODO: check ~equal to self.wavelengths
        self.con_incident = np.flip(d[:,1])
        self.con_trasmitted = np.flip(d[:,2])

    def calculate_extinguished_spectrum(self,recalculate=False):
        if recalculate :
            self.calculate_depths(True)
            self.read_incident_transmitted_spectrum()
        else:
            try:
                x=self.opac_abs_tau
            except:
                self.calculate_depths(True)
            try :
                x = self.con_incident
            except :
                self.read_incident_transmitted_spectrum()
        tau_abs_end = self.opac_abs_tau[-1,:]
        tau_scat_end = self.opac_scat_tau[-1,:]
        self.con_extinguished_abs = self.con_incident*np.exp(-tau_abs_end)
        self.con_extinguished_scat = self.con_incident*np.exp(-tau_scat_end)
        self.con_extinguished_tot = self.con_incident*np.exp(-tau_scat_end-tau_abs_end)
        print(self.con_extinguished_tot.shape)

    def calculate_rad_pressure(self,recalculate=False):
        """SpectrumInterpolator.calculate_rad_pressure()

        Using extinguished continuum
        """
        if recalculate :
            self.weighted_optical_depth(True)
        else:
            try:
                x=self.weighted_opacity
            except:
                self.weighted_optical_depth(True)
        self.arad = self.weighted_opacity * self.con_extinguished_tot_bolometric
        self.arad /= constants.c.to(units.cm/units.s).value
        self.arad /= logn_to_rho(self.n)

    def weighted_optical_depth(self,recalculate=False):
        if recalculate :
            self.calculate_depths(True)
            self.read_incident_transmitted_spectrum()
        else:
            try:
                x=self.opac_abs_tau
            except:
                self.calculate_depths(True)
            try :
                x = self.con_incident
            except :
                self.read_incident_transmitted_spectrum()

        extinguished_continuum_matrix = self.con_incident*np.exp(-self.opac_abs_tau - self.opac_scat_tau)/self.con_wavelengths
        self.con_extinguished_tot_bolometric = np.trapz(
                y= extinguished_continuum_matrix
                ,x=self.con_wavelengths)
        self.weighted_opacity = np.trapz(y=extinguished_continuum_matrix*(self.opac_abs_matrix+self.opac_scat_matrix)
                             ,x=self.con_wavelengths)/self.con_extinguished_tot_bolometric
        dr = np.gradient(self.radii)
        self.weighted_tau = np.cumsum(self.weighted_opacity*dr)



if __name__ == "__main__" :
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt

    profile = True
    if profile:
        prof = cProfile.Profile()
        prof.enable()

    spec = SpectrumInterpolator(1.,10.2,1.)
    print("Initialised!")
    # s=spec.get_spectrum(1.e24)
    # print("One spectrum!")

    # spec.weighted_optical_depth()
    spec.calculate_rad_pressure()
    np.savetxt("data/arad.txt",np.array(
                                        [spec.depths
                                        ,spec.weighted_tau
                                        ,spec.con_extinguished_tot_bolometric
                                        ,spec.arad
                                        ,spec.weighted_opacity/logn_to_rho(spec.n) #convert to mass units
                                         ]
                                        ).T)


    if profile:
        prof.disable()
        s = io.StringIO()
        sortby = SortKey.CUMULATIVE
        ps = pstats.Stats(prof, stream=s).sort_stats(sortby)
        ps.print_stats()
        print(s.getvalue())

    # spec.read_incident_transmitted_spectrum()
    #
    # # print(spec.wavelengths)
    # # print(spec.con_wavelengths)
    #
    # spec.gather_all_spectra()
    # spec.calculate_depths()
    # spec.calculate_extinguished_spectrum()
    #
    # outp = [spec.wavelengths,
    #         spec.con_wavelengths,
    #         spec.con_incident,
    #         spec.con_trasmitted,
    #         spec.con_extinguished_abs,
    #         spec.con_extinguished_scat,
    #         spec.con_extinguished_tot
    #         ]
    # outp = np.array(outp).T
    # np.savetxt("data/specout.txt",outp)
    #
    # print("plotting")
    # fig,sp = plt.subplots(figsize=(9,6))
    # sp.set_xscale('log')
    # sp.set_yscale('log')
    # sp.set_ylim([1.e-10,1e2])
    # # sp.set_ylim([1.e5,1.e11])
    # sp.plot(spec.con_wavelengths,spec.con_incident,label="Incident continuum")
    # sp.plot(spec.con_wavelengths,spec.con_trasmitted,label="Transmitted continuum")
    # sp.plot(spec.con_wavelengths,spec.con_extinguished_abs,label="Extinguished continuum (abs)")
    # sp.plot(spec.con_wavelengths,spec.con_extinguished_scat,label="Extinguished continuum (scat)")
    # sp.plot(spec.con_wavelengths,spec.con_extinguished_tot,label="Extinguished continuum (sum)")
    # sp.legend()
    # fig.savefig(f"pics/con.pdf")


    # alld = spec.gather_all_spectra()
    # #spec.calculate_depths()
    #
    # for irow,d in enumerate(alld):
    #     print("plotting")
    #     fig,sp = plt.subplots(figsize=(9,6))
    #     sp.set_xscale('log')
    #     sp.set_yscale('log')
    #     for i in range(0,d.shape[0],50):
    #         sp.plot(spec.wavelengths,d[i,:],label="{:6.3f}".format(np.log10(spec.depths[i])))
    #     sp.legend()
    #     fig.savefig(f"pics/spec{irow}.pdf")