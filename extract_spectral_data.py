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

import time

import tempfile

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
    # TODO: pare down to bare minimum to run at max efficiency for table

    def __init__(self,run_data):
        self.rd = run_data
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
        spec_files = glob.glob(f'{self.rd.data_dir}*znu*')


        self.nfiles = len(spec_files)

        # check filenames are formatted correctly, extract frequency of each file
        # self.wavelengths is never used?
        # freq = []
        # for i in range(self.nfiles):
        #     file = f"{self.rd.spec_file_base}{i}"
        #     if file not in spec_files:
        #         raise FileNotFoundError(f"{file} not in .znu* files - filenames may not be contiguous?")
        #     with open(f"{file}") as f:
        #         header_line = f.readline()
        #     freq.append(float(header_line[-18:-4]))
        # self.wavelengths = ryd2microns(np.array(freq))[::-1]

        # extract depths and distances from one file (should all be identical!)
        # file = f"{self.rd.spec_file_base}0"
        # d = np.loadtxt(file,usecols=[0,1],skiprows=1)

        # self.radii  = d[:,0]
        # self.depths = d[:,1]

    def get_spectrum(self,depth): # for inspection
        """SpectrumInterpolator.get_spectrum(depth)

        Generates the spectrum at a given column
        This reads one line from each of ~10k files, so is slow for generating multiple spectra
        """
        #TODO: proper interpolation. For now just using lower bound
        idepth = np.searchsorted(self.rd.depths,depth)
        emissivity = []

        # for i in range(self.nfiles) :
        for i in range(self.nfiles-1, -1, -1): # invert order to get increasing wavelength, not increasing frequency
            file = f"{self.rd.spec_file_base}{i}"
            emissivity.append(
                np.loadtxt(file,usecols=2,skiprows=1)[idepth] # load in emmisivity column, select correct depth
                )

        return np.array(emissivity)

    def calculate_extinguished_spectrum(self,recalculate=False):# for inspection
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



    def gather_all_spectra(self):
        """SpectrumInterpolator.gather_all_spectra()

        Combines all dumped emissivities and opacities (scattering and absorbed) into 3 big matrices, for quick spectrum plotting
        May use a lot of ram
        Returns matrices, also saves them as SpectrumInterpolator.emissivity_matrix,SpectrumInterpolator.opac_abs_matrix,SpectrumInterpolator.opac_scat_matrix
        """

        # combine all znu files into a single file, for quicker read_csv processing
        # this is much faster than using read_csv 8000 times
        with tempfile.NamedTemporaryFile() as tf:

            t0 = time.time()
            for nchars in range(1,4+1):
                wildcard = "?"*(nchars)
                if nchars==1:
                    pipe = ">"
                else:
                    pipe = ">>"
                #TODO: replace os.system with newer preferred version
                os.system(f"tail -q -n +2 {self.rd.spec_file_base}{wildcard} | cut -f 3,4,5 {pipe} {tf.name}")
            print(time.time()-t0)

            d = pd.read_csv(tf.name,header=None,sep='\t').values
        print(time.time()-t0)

        # reshape, then flip
        ndepth = d.shape[0]//self.nfiles
        # self.emissivity_matrix = np.flip(d[:,0].reshape((self.nfiles,ndepth)).T,axis=1) # only for check
        self.opac_abs_matrix = np.flip(d[:,1].reshape((self.nfiles,ndepth)).T,axis=1)
        self.opac_scat_matrix = np.flip(d[:,2].reshape((self.nfiles,ndepth)).T,axis=1)

        # return self.emissivity_matrix,self.opac_abs_matrix,self.opac_scat_matrix

    def calculate_depths(self,recalculate=False):
        if recalculate :
            self.gather_all_spectra()
        else:
            try:
                x=self.opac_abs_matrix
            except:
                self.gather_all_spectra()
        dr = np.gradient(self.rd.radii)
        self.opac_abs_tau = np.cumsum(self.opac_abs_matrix*dr[:,np.newaxis],axis=0)
        self.opac_scat_tau = np.cumsum(self.opac_scat_matrix*dr[:,np.newaxis],axis=0)

    def read_incident_transmitted_spectrum(self):
        confilename = self.rd.data_dir+pr.format_file(self.rd)+".con"

        d = np.loadtxt(confilename,usecols = [0,1,2])

        self.con_incident = np.flip(d[:,1])
        self.con_wavelengths = np.flip(ryd2microns(d[:,0]))
        # self.con_trasmitted = np.flip(d[:,2]) # only for checking

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
        self.weighted_abs_opacity = np.trapz(y=extinguished_continuum_matrix*self.opac_abs_matrix
                             ,x=self.con_wavelengths)/self.con_extinguished_tot_bolometric
        self.weighted_scat_opacity = np.trapz(y=extinguished_continuum_matrix*self.opac_scat_matrix
                             ,x=self.con_wavelengths)/self.con_extinguished_tot_bolometric
        self.weighted_opacity = self.weighted_abs_opacity + self.weighted_scat_opacity
        dr = np.gradient(self.rd.radii)
        self.weighted_tau = np.cumsum(self.weighted_opacity*dr)


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
        self.arad /= logn_to_rho(self.rd.n)

    def gen_dataframe(self):

        table_data = pd.DataFrame(np.c_[
                                        self.weighted_tau
                                        ,self.arad
                                        ,self.weighted_abs_opacity
                                        ,self.weighted_scat_opacity
                                       ],
                                       columns=[
                                            "tau"
                                           ,"arad"
                                           ,"kabs"
                                           ,"kscat"
                                           ])
        return table_data


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
        print(pstats
                .Stats(prof, stream=io.StringIO())
                .sort_stats(SortKey.CUMULATIVE)
                .print_stats()
                .getvalue()
              )
