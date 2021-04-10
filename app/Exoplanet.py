import numpy as np

# unit conversion constants
R_jup    = 7.1492        * 10**7   # in metres
M_jup    = 1.89818717    * 10**27  # in kilograms
R_sun    = 6.957         * 10**8   # in metres
M_sun    = 1.98847542    * 10**30  # in kilograms
L_sun    = 3.828         * 10**26  # in watts
au       = 1.49597870710 * 10**11  # in metres, exact
G        = 6.67408       * 10**-11 # in metres^3 / (kilograms * seconds^2)
sigma_sb = 5.670367      * 10**-8  # in watts / (kelvin^4 * metres^2)
day      = 60.0 * 60.0 * 24.0      # in seconds

class Exoplanet:
    '''
    Exoplanet class

    Attributes
    ----------
    Physical attributes:
        All physical properties are listed in Exoplanet.property_dict:
            name                                                                     : str
            planet_status                                                            : str
            mass, mass_error_min, mass_error_max                                     : float
            mass_sini, mass_sini_error_min, mass_sini_error_max                      : float
            radius, radius_error_min, radius_error_max                               : float
            orbital_period, orbital_period_error_min, orbital_period_error_max       : float
            semi_major_axis, semi_major_axis_error_min, semi_major_axis_error_max    : float
            eccentricity, eccentricity_error_min, eccentricity_error_max             : float
            inclination, inclination_error_min, inclination_error_max                : float
            angular_distance                                                         : float
            discovered                                                               : float
            updated                                                                  : float
            omega, omega_error_min, omega_error_max                                  : float
            tperi, tperi_error_min, tperi_error_max                                  : float
            tconj, tconj_error_min, tconj_error_max                                  : float
            tzero_tr, tzero_tr_error_min, tzero_tr_error_max                         : float
            tzero_tr_sec, tzero_tr_sec_error_min, tzero_tr_sec_error_max             : float
            lambda_angle, lambda_angle_error_min, lambda_angle_error_max             : float
            impact_parameter, impact_parameter_error_min, impact_parameter_error_max : float
            tzero_vr, tzero_vr_error_min, tzero_vr_error_max                         : float
            k, k_error_min, k_error_max                                              : float
            temp_calculated, temp_calculated_error_min, temp_calculated_error_max    : float
            temp_measured                                                            : float
            hot_point_lon                                                            : float
            geometric_albedo, geometric_albedo_error_min, geometric_albedo_error_max : float
            log_g                                                                    : float
            publication                                                              : str
            detection_type                                                           : str
            mass_detection_type                                                      : str
            radius_detection_type                                                    : str
            alternate_names                                                          : str
            molecules                                                                : str
        Host star properties are contained in a Star object:
            host.name                                                                : str
            host.ra                                                                  : float
            host.dec                                                                 : float
            host.mag_v                                                               : float
            host.mag_i                                                               : float
            host.mag_j                                                               : float
            host.mag_h                                                               : float
            host.mag_k                                                               : float
            host.distance, host.distance_error_min, host.distance_error_max          : float
            host.metallicity, host.metallicity_error_min, host.metallicity_error_max : float
            host.mass, host.mass_error_min, host.mass_error_max                      : float
            host.radius, host.radius_error_min, host.radius_error_max                : float
            host.sp_type                                                             : str
            host.age, host.age_error_min, host.age_error_max                         : float
            host.teff, host.teff_error_min, host.teff_error_max                      : float
            host.detected_disc                                                       : str
            host.magnetic_field                                                      : bool
            host.alternate_names                                                     : str
    Other (calculated) attributes:
        decision_metric    : float
        filter             : float
        t_exp              : float
        C_T                : float
        _metric_calculable : bool

    Methods
    -------
    temp_calculated_()
        Calculate or simply return planet temperature
    orbital_period_()
        Calculate or simply return orbital period
    semi_major_axis_()
        Calculate or simply return semi major axis
    inclination_()
        Calculate or simply return inclination angle
    transit_time_()
        Calculate transit time
    _metric_calculable_()
        Determine whether the metric is (in principle) calculable for self
    decision_metric_(telescope)
        Calculate decision metric with a given telescope according to arXiv:1802.05645
    '''
    
    property_dict = {
        'name'                       : '',
        'planet_status'              : '',
        'mass'                       : 'M_jup',
        'mass_error_min'             : 'M_jup',
        'mass_error_max'             : 'M_jup',
        'mass_sini'                  : 'M_jup',
        'mass_sini_error_min'        : 'M_jup',
        'mass_sini_error_max'        : 'M_jup',
        'radius'                     : 'R_jup',
        'radius_error_min'           : 'R_jup',
        'radius_error_max'           : 'R_jup',
        'orbital_period'             : 'days',
        'orbital_period_error_min'   : 'days',
        'orbital_period_error_max'   : 'days',
        'semi_major_axis'            : 'AU',
        'semi_major_axis_error_min'  : 'AU',
        'semi_major_axis_error_max'  : 'AU',
        'eccentricity'               : 'unitless',
        'eccentricity_error_min'     : 'unitless',
        'eccentricity_error_max'     : 'unitless',
        'inclination'                : 'deg',
        'inclination_error_min'      : 'deg',
        'inclination_error_max'      : 'deg',
        'angular_distance'           : 'arcsec',
        'discovered'                 : 'year',
        'updated'                    : 'date',
        'omega'                      : 'deg',
        'omega_error_min'            : 'deg',
        'omega_error_max'            : 'deg',
        'tperi'                      : 'JD',
        'tperi_error_min'            : 'JD',
        'tperi_error_max'            : 'JD',
        'tconj'                      : 'JD',
        'tconj_error_min'            : 'JD',
        'tconj_error_max'            : 'JD',
        'tzero_tr'                   : 'JD',
        'tzero_tr_error_min'         : 'JD',
        'tzero_tr_error_max'         : 'JD',
        'tzero_tr_sec'               : 'JD',
        'tzero_tr_sec_error_min'     : 'JD',
        'tzero_tr_sec_error_max'     : 'JD',
        'lambda_angle'               : 'deg',
        'lambda_angle_error_min'     : 'deg',
        'lambda_angle_error_max'     : 'deg',
        'impact_parameter'           : '%',
        'impact_parameter_error_min' : '%',
        'impact_parameter_error_max' : '%',
        'tzero_vr'                   : 'JD',
        'tzero_vr_error_min'         : 'JD',
        'tzero_vr_error_max'         : 'JD',
        'k'                          : 'm/s',
        'k_error_min'                : 'm/s',
        'k_error_max'                : 'm/s',
        'temp_calculated'            : 'K',
        'temp_calculated_error_min'  : 'K',
        'temp_calculated_error_max'  : 'K',
        'temp_measured'              : 'K',
        'hot_point_lon'              : 'deg',
        'geometric_albedo'           : 'unitless',
        'geometric_albedo_error_min' : 'unitless',
        'geometric_albedo_error_max' : 'unitless',
        'log_g'                      : 'unitless',
        'publication'                : '',
        'detection_type'             : '',
        'mass_detection_type'        : '',
        'radius_detection_type'      : '',
        'alternate_names'            : '',
        'molecules'                  : '',

        # host starproperties now contained in a star object instead

        # 'star_name'                  : '',
        # 'ra'                         : 'deg',
        # 'dec'                        : 'deg',
        # 'mag_v'                      : 'unitless',
        # 'mag_i'                      : 'unitless',
        # 'mag_j'                      : 'unitless',
        # 'mag_h'                      : 'unitless',
        # 'mag_k'                      : 'unitless',
        # 'star_distance'              : 'pc',
        # 'star_distance_error_min'    : 'pc',
        # 'star_distance_error_max'    : 'pc',
        # 'star_metallicity'           : 'unitless',
        # 'star_metallicity_error_min' : 'unitless',
        # 'star_metallicity_error_max' : 'unitless',
        # 'star_mass'                  : 'M_sun',
        # 'star_mass_error_min'        : 'M_sun',
        # 'star_mass_error_max'        : 'M_sun',
        # 'star_radius'                : 'R_sun',
        # 'star_radius_error_min'      : 'R_sun',
        # 'star_radius_error_max'      : 'R_sun',
        # 'star_sp_type'               : '',
        # 'star_age'                   : 'Gyr',
        # 'star_age_error_min'         : 'Gyr',
        # 'star_age_error_max'         : 'Gyr',
        # 'star_teff'                  : 'K',
        # 'star_teff_error_min'        : 'K',
        # 'star_teff_error_max'        : 'K',
        # 'star_detected_disc'         : '',
        # 'star_magnetic_field'        : '',
        # 'star_alternate_names'       : '',
    }
    
    def __init__(self):
        '''Exoplanet initialiser'''
        # initialiser setting all properties to int 0
        # (to distinguish between new objects and objects with blank values, read in as blank strings)

        # initialise all values in the dictionary to 0
        for prop in self.property_dict.keys():
            setattr(self, prop, 0)
         
        # other properties to initialise
        self.host = Star()
        self.decision_metric = 0
        self.filter = 0
        self.t_exp = 0
        self.C_T = 0
        self._metric_calculable = 0
    
    def temp_calculated_(self):
        '''Function to calculate planet temperature'''

        # define symbols
        T_eq = self.temp_calculated
        A    = self.geometric_albedo
        R_s  = self.host.radius_()
        T_s  = self.host.teff_()
        a    = self.semi_major_axis_()
        
        if R_s != '' and T_s != '' and a != '':
            if T_eq != '':
                # temperature already given
                T_eq = self.temp_calculated
            elif A != '':
                # calculate from albedo
                T_eq = (1.0 - A)**0.25 * (R_s*R_sun / (2.0 * a*au))**0.5 * T_s
            else:
                # calculate with A = 0.3
                A = 0.3
                T_eq = (1.0 - A)**0.25 * (R_s*R_sun / (2.0 * a*au))**0.5 * T_s
            
            # set and return
            self.temp_calculated = T_eq
            return T_eq
        
        else:
            # mising parameters, return measured temperature instead
            return self.temp_measured
    
    def orbital_period_(self):
        '''Function to calculate planet orbital period'''

        # define symbols
        P   = self.orbital_period
        a   = self.semi_major_axis
        M_s = self.host.mass_()
        M_p = self.mass
        
        if a != '' and M_s != '':
            if P != '':
                # period already given
                P = self.orbital_period
            elif M_p != '':
                # calculate with planet mass
                P = (4.0*np.pi**2.0 / (G * (M_s*M_sun + M_p*M_jup)) * (a*au)**3.0) ** 0.5 / day
            else:
                # calculate with star mass only
                P = (4.0*np.pi**2.0 / (G * (M_s*M_sun)) * (a*au)**3.0) ** 0.5 / day
            
            # set and return
            self.orbital_period = P
            return P
        
        else:
            # mising parameters
            return P
    
    def semi_major_axis_(self):
        '''Function to calculate planet semi major axis'''

        # define symbols
        a   = self.semi_major_axis
        P   = self.orbital_period
        M_s = self.host.mass_()
        M_p = self.mass
        
        if P != '' and M_s != '':
            if a != '':
                # period already given
                a = self.semi_major_axis
            elif M_p != '':
                # calculate with planet mass
                a = ( G * (M_s*M_sun + M_p*M_jup) / (4.0*np.pi**2.0) * (P*day)**2.0 ) ** (1.0/3.0) / au
            else:
                # calculate with star mass only
                a = ( G * (M_s*M_sun) / (4.0*np.pi**2.0) * (P*day)**2.0 ) ** (1.0/3.0) / au
            
            # set and return
            self.semi_major_axis = a
            return a
        
        else:
            # mising parameters
            return a
    
    def inclination_(self):
        '''Function to get planet orbit inclination'''

        # define symbols
        i = self.inclination
        
        if i != '':
            # inclination already given
            i = self.inclination
        elif self.detection_type == 'Primary Transit' or self.radius_detection_type == 'Primary Transit':
            # assume 90 degree inclination for transits
            i = 90.0
        
        # set and return
        self.inclination = i
        return i

    def transit_time_(self):
        '''Function to calculate planet transit time'''

        # define symbols
        R_s = self.host.radius
        R_p = self.radius
        a   = self.semi_major_axis
        i   = np.radians(self.inclination_())
        e   = self.eccentricity
        P   = self.orbital_period
        # calculation can be done without omega, so it is not checked for in _metric_calculable_()
        try: w = np.radians(self.omega)
        except: w = ''

        # attempt with omega
        if w != '' and e < 1:
            topfrac = (R_s*R_sun + R_p*R_jup)**2.0 - (a*au * np.cos(i) * (1.0-e**2.0)/(1.0+e*np.sin(w)))**2.0
            if topfrac > 0:
                frac = topfrac ** 0.5 / (a*au * np.sin(i))
            else: return 0 # cannot sqrt a negative number
            
            t_14 = P*day/np.pi * np.arcsin(frac) * np.sqrt(1.0 - e**2.0)/(1.0 + e*np.sin(w))
           
        # attempt without omega 
        else:
            topfrac = (R_s*R_sun + R_p*R_jup)**2.0 - (a*au * np.cos(i))**2.0
            if topfrac > 0:
                frac = topfrac ** 0.5 / (a*au)
            else: return 0 # cannot sqrt a negative number
            
            t_14 = P*day/np.pi * np.arcsin(frac)
        
        if np.isnan(t_14): return 0 # something wrong with np.arcsin(frac), like an input not in range +/- pi/2
        else: return t_14
    
    def _metric_calculable_(self):
        '''Function to see if a metric is calculatable or not'''

        # check for missing parameters
        missing_parameters = False
        # planet radius necessary in several steps
        if self.radius == '': missing_parameters = True
        # at least one magnitude
        if self.host.mag_v == '' and self.host.mag_i == '' and self.host.mag_j == '' and self.host.mag_h == '' and self.host.mag_k == '': missing_parameters = True
        # for t_14
        if self.orbital_period_() == '' or self.semi_major_axis_() == '' or self.host.radius_() == '': missing_parameters = True
        # for T_eq
        if self.temp_calculated_() == '': missing_parameters = True
        # inclination
        if self.inclination_() == '': missing_parameters = True
        
        if missing_parameters:
            # cannot calculate metric for this planet
            self._metric_calculable = False
            return False
        else:
            self._metric_calculable = True
            return True

    def decision_metric_(self, telescope):
        '''
        Function to calculate the decision metric for planet and given telescope
        Method described in arXiv:1802.05645

        Parameters
        ----------
        telescope : Telescope
            An instance of a Telescope object

        Returns
        -------
        float
            The final decision metric score (if incalculable returns 0)
        '''
        
        # thoughout, will return 0 if there is no way to calculate further
        
        # determine whether necessary values are usable or calculable
        if self._metric_calculable != True: return 0
        
        # define symbols
        R_s  = self.host.radius
        R_p  = self.radius
        M_p  = self.mass
        a    = self.semi_major_axis
        i    = self.inclination
        e    = self.eccentricity
        P    = self.orbital_period
        T_eq = self.temp_calculated
            
        # calculate t_14
        t_14 = self.transit_time_()
        if t_14 == 0: return 0
        
        # calculate delta
        delta = (R_p*R_jup / (R_s*R_sun)) ** 2.0
        
        # calculate C_T from telescope input (INCLUDES 10^-0.2mag already)
        (C_T, best_filter, t_exp) = telescope._tele_coeff_(self)
        if C_T == 0: return 0
        self.C_T    = C_T
        self.filter = best_filter
        self.t_exp  = t_exp
        
        # determine whether to use regular metric or known mass metric
        if self.mass != '':
            # calculate Dmass
            Dmass = C_T * t_14**0.5 * T_eq * R_p / M_p * delta * (1.0 - delta)**-0.5
            
            self.decision_metric = Dmass
            return Dmass
        
        # use mass-radius relation
        else:
            # gas giants
            if self.radius >= 0.8: 
                D = C_T * t_14**0.5 * T_eq * R_p * delta
            else:
                # rocky
                if self.radius <= 0.18:
                    a = 0.252721
                    b = -0.919588
                # icy
                elif self.radius <= 0.292:
                    a = 0.085456
                    b = -1.193344
                # transitional
                else:
                    a = 0.466035
                    b = 0.552415

                D = C_T * t_14**0.5 * T_eq * np.exp(b/a) * (R_p)**(1.0 - 1/a) * delta
            
            # set and return values
            self.decision_metric = D
            
            return D


class Star:
    '''
    Star class

    Attributes
    ----------
    Physical attributes:
        All physical attributes are listed in Star.property_dict:
            name                                                      : str
            ra                                                        : float
            dec                                                       : float
            mag_v                                                     : float
            mag_i                                                     : float
            mag_j                                                     : float
            mag_h                                                     : float
            mag_k                                                     : float
            distance, distance_error_min, distance_error_max          : float
            metallicity, metallicity_error_min, metallicity_error_max : float
            mass, mass_error_min, mass_error_max                      : float
            radius, radius_error_min, radius_error_max                : float
            sp_type                                                   : str
            age, age_error_min, age_error_max                         : float
            teff, teff_error_min, teff_error_max                      : float
            detected_disc                                             : str
            magnetic_field                                            : bool
            alternate_names                                           : str

    Methods
    -------
    mass_()
        Calculate or simply return star mass
    radius_()
        Calculate or simply return star radius
    teff_()
        Calculate or simply return star effective temperature
    '''

    property_dict = {
        'name'                  : '',
        'ra'                    : 'deg',
        'dec'                   : 'deg',
        'mag_v'                 : 'unitless',
        'mag_i'                 : 'unitless',
        'mag_j'                 : 'unitless',
        'mag_h'                 : 'unitless',
        'mag_k'                 : 'unitless',
        'distance'              : 'pc',
        'distance_error_min'    : 'pc',
        'distance_error_max'    : 'pc',
        'metallicity'           : 'unitless',
        'metallicity_error_min' : 'unitless',
        'metallicity_error_max' : 'unitless',
        'mass'                  : 'M_sun',
        'mass_error_min'        : 'M_sun',
        'mass_error_max'        : 'M_sun',
        'radius'                : 'R_sun',
        'radius_error_min'      : 'R_sun',
        'radius_error_max'      : 'R_sun',
        'sp_type'               : '',
        'age'                   : 'Gyr',
        'age_error_min'         : 'Gyr',
        'age_error_max'         : 'Gyr',
        'teff'                  : 'K',
        'teff_error_min'        : 'K',
        'teff_error_max'        : 'K',
        'detected_disc'         : '',
        'magnetic_field'        : '',
        'alternate_names'       : '',
    }

    def __init__(self):
        '''Star initialiser'''
        # initialiser setting all properties to int 0
        # (to distinguish between new objects and objects with blank values, read in as blank strings)

        # initialise all values in the dictionary to 0
        for prop in self.property_dict.keys():
            setattr(self, prop, 0)
    
    def mass_(self):
        '''Function to calculate star mass'''

        # define symbols
        M_s = self.mass
        R_s = self.radius
        T_s = self.teff
        
        if R_s != '' and T_s != '':
            if M_s != '':
                # mass already given
                M_s = self.mass
            else:
                # find luminosity first from Stefan-Boltzmann law
                L = 4.0*np.pi*sigma_sb * (R_s*R_sun)**2.0 * T_s**4.0
                
                # use power law to find mass (valid for 0.43 M_sun - 2 M_sun)
                M_s = (L / L_sun)**0.25
                
            # set and return
            self.mass = M_s
            return M_s
            
        else:
            # mising parameters
            return M_s
    
    def radius_(self):
        '''Function to calculate star radius'''

        # define symbols
        R_s = self.radius
        M_s = self.mass
        T_s = self.teff
        
        if M_s != '' and T_s != '':
            if R_s != '':
                # mass already given
                R_s = self.radius
            else:
                # find luminosity first from power law
                if                 M_s < 0.43: L =  0.23 * M_s ** 2.3 * L_sun
                if 0.43 <= M_s and M_s <  2.0: L =         M_s ** 4.0 * L_sun
                if 2.0  <= M_s and M_s < 55.0: L =   1.4 * M_s ** 3.5 * L_sun
                if 55.0 <= M_s               : L = 32000 * M_s        * L_sun
                
                # use Stefan Boltzmann law to find radius
                R_s = (L / (4.0*np.pi*sigma_sb * T_s**4.0)) ** 0.5 / R_sun
                
            # set and return
            self.radius = R_s
            return R_s
            
        else:
            # mising parameters
            return R_s
    
    def teff_(self):
        '''Function to calculate star temperature'''

        # define symbols
        T_s = self.teff
        M_s = self.mass
        R_s = self.radius
        
        if M_s != '' and R_s != '':
            if T_s != '':
                # temperature already given
                T_s = self.teff
            else:
                # find luminosity first from power law
                if                 M_s < 0.43: L =  0.23 * M_s ** 2.3 * L_sun
                if 0.43 <= M_s and M_s <  2.0: L =         M_s ** 4.0 * L_sun
                if 2.0  <= M_s and M_s < 55.0: L =   1.4 * M_s ** 3.5 * L_sun
                if 55.0 <= M_s               : L = 32000 * M_s        * L_sun
                
                # use Stefan Boltzmann law to find radius
                T_s = (L / (4.0*np.pi*sigma_sb * (R_s*R_sun)**2.0)) ** 0.25
                
            # set and return
            self.teff = T_s
            return T_s
            
        else:
            # mising parameters
            return T_s
    
