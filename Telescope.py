# a class for telescopes
class Telescope:
    
    # class attribute dictionary for telescope properties (all possible properties left in for completeness)
    property_dict = {
        'name'                       : 0,
        #'Element?'                   : 1,
        #'mzp_U'                      : 2,
        #'mzp_u'                      : 3,
        #'mzp_B'                      : 4,
        #'mzp_g'                      : 5,
        'mzp_V'                      : 6,
        #'mzp_r'                      : 7,
        #'mzp_R'                      : 8,
        #'mzp_i'                      : 9,
        'mzp_I'                      : 10,
        #'mzp_z'                      : 11,
        'mzp_J'                      : 12,
        'mzp_H'                      : 13,
        'mzp_K'                      : 14,
        #'mzp_HK'                     : 15,
        #'mzp_T'                      : 16,
        #'mzp_vr'                     : 17,
        #'mzp_RISE'                   : 18,
        #'mzp_SPRAT_Red'              : 19,
        #'mzp_SPRAT_Blue'             : 20,
        #'mzp_600RI+19'               : 21,
        #'mzp_G750M'                  : 22,
        #'mzp_GR01'                   : 23,
        #'mzp_Band1'                  : 24,
        #'mzp_Band2'                  : 25,
        #'mzp_Band3'                  : 26,
        'spectral_bins'              : 27,
        #'Lambda_Range'               : 28,
        #'Dispersion'                 : 29,
        't_over'                     : 30,
        #'Theta_see'                  : 31,
        #'Theta_DF'                   : 32,
        #'Omega_pix'                  : 33,
        #'Gain'                       : 34,
        #'Half_Well'                  : 35,
        #'texp_min'                   : 36
        }
    
    def __init__(self):
        # initialise all values in the dictionary to 0
        for key, val in self.property_dict.items():
            setattr(self, key, 0)
         
        # other properties to initialise
        self.t_exp = 0
        self.filter = 0
        self.C_T = 0
     
    # function to calculate telescope coefficient
    def _tele_coeff(self):
        # determine whether necessary values are usable or calculable
        missing_parameters = False
        if self.spectral_bins == '' or self.t_exp == '' or self.t_over == '': missing_parameters = True
        if self.mzp_V == '' and self.mzp_I == '' and self.mzp_J == '' and self.mzp_H == '' and self.mzp_K == '': missing_parameters = True
        if missing_parameters:
            # cannot calculate coefficient for this telescope, return 0
            self.C_T = 0
            return 0
        
        # determine which filter to use (i.e highest mzp)
        filter_mzps = [self.mzp_V, self.mzp_I, self.mzp_J, self.mzp_H, self.mzp_K]
        mzp = max([x for x in filter_mzps if x !='']) # find highest nonzero value
        if mzp == filter_mzps[0]: self.filter = 'v'
        if mzp == filter_mzps[1]: self.filter = 'i'
        if mzp == filter_mzps[2]: self.filter = 'j'
        if mzp == filter_mzps[3]: self.filter = 'h'
        if mzp == filter_mzps[4]: self.filter = 'k'
         
        # calculate coefficient
        C_T = self.spectral_bins**-0.5 * 10**(0.2 * mzp) * (self.t_exp / (self.t_exp + self.t_over))**0.5
        self.C_T = C_T
        return C_T