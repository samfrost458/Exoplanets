import numpy as np
from scipy import special

class Telescope:
    '''
    Telescope class

    Attributes
    ----------
    Physical attributes:
        All read-in properties are listed in Telescope.property_dict:
            name          : str
            mzp_V         : float
            mzp_I         : float
            mzp_J         : float
            mzp_H         : float
            mzp_K         : float
            spectral_bins : float
            t_over        : float
            Theta_see     : float
            Theta_DF      : float
            Omega_pix     : float
            Half_Well     : float
            texp_min      : float
    Other (calculated) attributes:
        _C_T_calculable   : bool

    Methods
    -------
    exposure_time_(mag, mzp, DF = False, bgd = 0)
        Calculate an exposure time for a given apparent magnitude, zero point magnitude and background signal
    _C_T_calculable_()
        Determine whether C_T is (in principle) calculable for self
    _tele_coeff_(planet)
        Calculate C_T with a given planet, according to arXiv:1802.05645
    '''
    
    property_dict = {
        'name'           : 0,
        #'Element?'       : 1,
        #'mzp_U'          : 2,
        #'mzp_u'          : 3,
        #'mzp_B'          : 4,
        #'mzp_g'          : 5,
        'mzp_V'          : 6,
        #'mzp_r'          : 7,
        #'mzp_R'          : 8,
        #'mzp_i'          : 9,
        'mzp_I'          : 10,
        #'mzp_z'          : 11,
        'mzp_J'          : 12,
        'mzp_H'          : 13,
        'mzp_K'          : 14,
        #'mzp_HK'         : 15,
        #'mzp_T'          : 16,
        #'mzp_vr'         : 17,
        #'mzp_RISE'       : 18,
        #'mzp_SPRAT_Red'  : 19,
        #'mzp_SPRAT_Blue' : 20,
        #'mzp_600RI+19'   : 21,
        #'mzp_G750M'      : 22,
        #'mzp_GR01'       : 23,
        #'mzp_Band1'      : 24,
        #'mzp_Band2'      : 25,
        #'mzp_Band3'      : 26,
        'spectral_bins'  : 27,
        #'Lambda_Range'   : 28,
        #'Dispersion'     : 29,
        't_over'         : 30,
        'Theta_see'      : 31,
        # 'Theta_DF'       : 32,
        'Omega_pix'      : 33,
        #'Gain'           : 34,
        'Half_Well'      : 35,
        'texp_min'       : 36,
    }
    filter_list = [prop.replace('mzp_', '') for prop in property_dict.keys() if 'mzp' in prop]
    
    def __init__(self):
        '''Telescope initialiser'''
        
        # initialise all values in the dictionary to 0
        for key, val in self.property_dict.items():
            setattr(self, key, 0)
         
        # other properties to initialise
        self._C_T_calculable = 0

    def exposure_time_(self, mag, mzp, DF = False, bgd = 0):
        '''Function to calculate exposure times'''

        # weighting stuff here
        if DF: psf_sigma = self.Theta_DF  / 2.35482
        else:  psf_sigma = self.Theta_see / 2.35482
        pixel_scale = np.sqrt(self.Omega_pix) / psf_sigma # in units of psf_sigma
        
        # gaussian integral is the error function 
        central_pix_weight = special.erf(pixel_scale / (2*np.sqrt(2))) ** 2.0

        # calculate exposure
        t_exp = self.Half_Well * (10 ** (-0.4 * (mzp - mag)) * central_pix_weight + bgd)
        return t_exp

    def _C_T_calculable_(self):
        '''Function to calculate whether C_T is calculable or not'''
        
        missing_parameters, at_least_one_filter = False, False
        if self.spectral_bins == '' or self.t_over == '' or self.Theta_see == '': missing_parameters = True
        for filter_ in self.filter_list:
            if getattr(self, 'mzp_' + filter_) != '' and getattr(self, 'mzp_' + filter_) is not None: at_least_one_filter = True
        if not at_least_one_filter: missing_parameters = True
        
        if missing_parameters: self._C_T_calculable = False
        else: self._C_T_calculable = True
        return self._C_T_calculable
    
    def _tele_coeff_(self, planet):
        '''
        Function to calculate the telescope coefficient (a.k.a. C_T) for a given planet and telescope
        Method described in arXiv:1802.05645

        Parameters
        ----------
        planet : Exoplanet
            An instance of an Exoplanet object

        Returns
        -------
        Tuple:
            C_T : float
                The final telescope coefficient C_T (if incalculable returns 0)
            best_filter : str
                The filter which gives the best C_T, out of those available (if incalculable returns 'nan')
            best_t_exp : float
                The exposure time to fill the telescope to half-well for the calculated best_filter (if incalculable returns 'nan')
        '''
        
        # determine whether necessary values are usable or calculable
        if self._C_T_calculable != True:
            return (0, 'nan', 'nan')
        
        if self.texp_min == '' or self.texp_min is None: self.texp_min = 0
        
        # calculate t_exp and C_T for each filter
        t_exps, C_Ts = [], []
        for filter_ in self.filter_list:
            # calculate t_exp
            mzp = getattr(self, 'mzp_%s' %filter_)
            mag = getattr(planet.host, 'mag_%s' %filter_.lower())
            t_exp = ''
            if mag != '' and (mzp != '' and mzp is not None):
                t_exp = self.exposure_time_(mag, mzp)
                if t_exp >= self.texp_min: t_exps.append(t_exp)
                else: t_exps.append('') # exposure time required does not meet telescope minimum, append blank
            else: t_exps.append('') # both no planet magnitude or telescope mzp for this band, append blank

            # calculate C_T
            if t_exp != '' and self.t_over != '' and t_exp >= self.t_over:
                C_Ts.append( self.spectral_bins**-0.5 * 10**(0.2 * (mzp - mag)) * (t_exp / (t_exp + self.t_over))**0.5 )
            else: C_Ts.append(0) # divide by zero error occured or t_exp < t_over, append 0
        
        # find highest nonzero value and set filter and corresponding t_exp
        best_index = np.argmax(C_Ts)
        C_T_max = C_Ts[best_index]
        best_filter = self.filter_list[best_index]
        best_t_exp = t_exps[best_index]
        if C_T_max == 0:
            best_filter = 'nan'
            best_t_exp  = 'nan'
        
        return (C_T_max, best_filter, best_t_exp)

