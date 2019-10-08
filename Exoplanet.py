# a class for exoplanets
import math
import Telescope
class Exoplanet:
    
    # class attribute dictionary for planet properties (same as exoplanet.eu) (all possible properties left in for completeness)
    property_dict = {
        'name'                       : 0,
        #'planet_status'              : 1,
        'mass'                       : 2,
        'mass_error_min'             : 3,
        'mass_error_max'             : 4,
        #'mass_sini'                  : 5,
        #'mass_sini_error_min'        : 6,
        #'mass_sini_error_max'        : 7,
        'radius'                     : 8,
        'radius_error_min'           : 9,
        'radius_error_max'           : 10,
        'orbital_period'             : 11,
        'orbital_period_error_min'   : 12,
        'orbital_period_error_max'   : 13,
        'semi_major_axis'            : 14,
        'semi_major_axis_error_min'  : 15,
        'semi_major_axis_error_max'  : 16,
        'eccentricity'               : 17,
        'eccentricity_error_min'     : 18,
        'eccentricity_error_max'     : 19,
        'inclination'                : 20,
        'inclination_error_min'      : 21,
        'inclination_error_max'      : 22,
        #'angular_distance'           : 23,
        #'discovered'                 : 24,
        #'updated'                    : 25,
        'omega'                      : 26,
        'omega_error_min'            : 27,
        'omega_error_max'            : 28,
        #'tperi'                      : 29,
        #'tperi_error_min'            : 30,
        #'tperi_error_max'            : 31,
        #'tconj'                      : 32,
        #'tconj_error_min'            : 33,
        #'tconj_error_max'            : 34,
        #'tzero_tr'                   : 35,
        #'tzero_tr_error_min'         : 36,
        #'tzero_tr_error_max'         : 37,
        #'tzero_tr_sec'               : 38,
        #'tzero_tr_sec_error_min'     : 39,
        #'tzero_tr_sec_error_max'     : 40,
        #'lambda_angle'               : 41,
        #'lambda_angle_error_min'     : 42,
        #'lambda_angle_error_max'     : 43,
        #'impact_parameter'           : 44,
        #'impact_parameter_error_min' : 45,
        #'impact_parameter_error_max' : 46,
        #'tzero_vr'                   : 47,
        #'tzero_vr_error_min'         : 48,
        #'tzero_vr_error_max'         : 49,
        #'k'                          : 50,
        #'k_error_min'                : 51,
        #'k_error_max'                : 52,
        'temp_calculated'            : 53,
        'temp_calculated_error_min'  : 54,
        'temp_calculated_error_max'  : 55,
        'temp_measured'              : 56,
        #'hot_point_lon'              : 57,
        'geometric_albedo'           : 58,
        'geometric_albedo_error_min' : 59,
        'geometric_albedo_error_max' : 60,
        #'log_g'                      : 61,
        #'publication'                : 62,
        #'detection_type'             : 63,
        #'mass_detection_type'        : 64,
        #'radius_detection_type'      : 65,
        #'alternate_names'            : 66,
        #'molecules'                  : 67,
        #'star_name'                  : 68,
        'ra'                         : 69,
        'dec'                        : 70,
        'mag_v'                      : 71,
        'mag_i'                      : 72,
        'mag_j'                      : 73,
        'mag_h'                      : 74,
        'mag_k'                      : 75,
        #'star_distance'              : 76,
        #'star_distance_error_min'    : 77,
        #'star_distance_error_max'    : 78,
        #'star_metallicity'           : 79,
        #'star_metallicity_error_min' : 80,
        #'star_metallicity_error_max' : 81,
        #'star_mass'                  : 82,
        #'star_mass_error_min'        : 83,
        #'star_mass_error_max'        : 84,
        'star_radius'                : 85,
        'star_radius_error_min'      : 86,
        'star_radius_error_max'      : 87,
        #'star_sp_type'               : 88,
        #'star_age'                   : 89,
        #'star_age_error_min'         : 90,
        #'star_age_error_max'         : 91,
        'star_teff'                  : 92,
        'star_teff_error_min'        : 93,
        'star_teff_error_max'        : 94,
        #'star_detected_disc'         : 95,
        #'star_magnetic_field'        : 96,
        #'star_alternate_names'       : 97
    }
     
    # initialiser setting all properties to 0
    def __init__(self):
        # initialise all values in the dictionary to 0
        for key, val in self.property_dict.items():
            setattr(self, key, 0)
         
        # other properties to initialise
        self.decision_metric = 0
        self.decision_metric_error = 0
        self.filter = 0
     
    # to do: telescope database, telescope coefficient, decision metric function
     
    # function to calculate decision metric
    def _decision_metric(self, telescope):
        # determine whether necessary values are usable or calculable
        missing_parameters = False
        if self.mag_v == '' and self.mag_i == '' and self.mag_j == '' and self.mag_h == '' and self.mag_k == '': missing_parameters = True
        if self.orbital_period == '' or self.star_radius == '' or self.semi_major_axis == '' or self.inclination == '': missing_parameters = True
        if self.temp_measured == '' and self.temp_calculated == '' or self.star_radius == '' and self.star_teff == '': missing_parameters = True
        if self.radius == '': missing_parameters = True
        if missing_parameters:
            # cannot calculate metric for this planet, return 0
            self.decision_metric = 0
            return 0
         
        # determine which filter to use (i.e lowest magnitude)
        filter_mags = [self.mag_v, self.mag_i, self.mag_j, self.mag_h, self.mag_k]
        mag = min([x for x in filter_mags if x !='']) # find lowest nonzero value
        if mag == filter_mags[0]: self.filter = 'v'
        if mag == filter_mags[1]: self.filter = 'i'
        if mag == filter_mags[2]: self.filter = 'j'
        if mag == filter_mags[3]: self.filter = 'h'
        if mag == filter_mags[4]: self.filter = 'k'
         
        # determine which temperature to use
        if self.temp_measured != '': temp = self.temp_measured
        elif self.temp_calculated != '': temp = self.temp_calculated
        elif self.geometric_albedo != '': temp = (1.0 - self.geometric_albedo)**0.25 * (self.star_radius / (2.0 * self.semi_major_axis))**0.5 * self.star_teff
        else: temp = (self.star_radius / (2.0 * self.semi_major_axis))**0.5 * self.star_teff
         
        # calculate t_14
        M_Jup = 1.898 * 10**27 # in kg
        R_Jup = 6.9911 * 10**7 # in metres
        R_sun = 6.955 * 10**8 # in metres
        AU = 1.49597870700 * 10**11 # in metres, exact
        day = 60.0 * 60.0 * 24.0 # in seconds
        G = 6.67408 * 10**11 # in m^3 kg^-1 s^-2
        try:
            sqrt = ( (self.star_radius*R_sun + self.radius*R_Jup)**2.0 - (self.semi_major_axis*AU * math.cos(math.radians(self.inclination)))**2.0 )**0.5
            t_14 = self.orbital_period*day/math.pi * math.asin( sqrt / (self.semi_major_axis*AU * math.sin(math.radians(self.inclination))) )
            
            # calculate delta
            delta = (self.radius*R_Jup / (self.star_radius*R_sun)) ** 2.0
            
            # calculate C_T from telescope input
            #C_T = telescope.C_T()
             
            # determine whether to use regular metric or known mass metric
            if self.mass != '':
                D = 10.0**(-0.2 * mag) * t_14**0.5 * temp * self.radius / (self.mass) * delta * (1.0 - delta)**-0.5
                
                
                #print('Dmass = ' + str(D))
                if self.radius > 0.8 : n = 0.0
                else: n = 2.0
                Dreg = 10.0**(-0.2 * mag) * t_14**0.5 * temp * (self.radius / 0.8)**(1.0 - n) * delta * (1.0 - delta)**-0.5
                #print('Dreg  = ' + str(Dreg))
                #print(D/Dreg)
                
            else:
                # determine n (radius is already in units of Jupiter radii)
                if self.radius > 0.8 : n = 0.0
                else: n = 2.0
                             
                D = 10.0**(-0.2 * mag) * t_14**0.5 * temp * (self.radius / 0.8)**(1.0 - n) * delta * (1.0 - delta)**-0.5
             
            # set and return values
            self.decision_metric = D
            return D
         
        except:
            # some problem occured with t_14 (likely sqrt returned complex value), return 0
            self.decision_metric = 0
            return 0