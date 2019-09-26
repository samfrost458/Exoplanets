# 24/9/19 playing around with web access in python

import urllib
import re

# a class for quantities with error values, initialised with all three values
class quantity:
    def __init__(self, value, error_min, error_max):
        self.value = value
        self.error_min = error_min
        self.error_max = error_max

# a class for exoplanets, initialised by planet name
class exoplanet:
    def __init__(self, name):
        self.name = name
        self.planet_status = 0
        self.mass = 0
        self.mass_error_min = 0
        self.mass_error_max = 0
        self.mass_sini = 0
        self.mass_sini_error_min = 0
        self.mass_sini_error_max = 0
        self.radius = 0
        self.radius_error_min = 0
        self.radius_error_max = 0
        self.orbital_period = 0
        self.orbital_period_error_min = 0
        self.orbital_period_error_max = 0
        self.semi_major_axis = 0
        self.semi_major_axis_error_min = 0
        self.semi_major_axis_error_max = 0
        self.eccentricity = 0
        self.eccentricity_error_min = 0
        self.eccentricity_error_max = 0
        self.inclination = 0
        self.inclination_error_min = 0
        self.inclination_error_max = 0
        self.angular_distance = 0
        self.discovered = 0
        self.updated = 0
        self.omega = 0
        self.omega_error_min = 0
        self.omega_error_max = 0
        self.tperi = 0
        self.tperi_error_min = 0
        self.tperi_error_max = 0
        self.tconj = 0
        self.tconj_error_min = 0
        self.tconj_error_max = 0
        self.tzero_tr = 0
        self.tzero_tr_error_min = 0
        self.tzero_tr_error_max = 0
        self.tzero_tr_sec = 0
        self.tzero_tr_sec_error_min = 0
        self.tzero_tr_sec_error_max = 0
        self.lambda_angle = 0
        self.lambda_angle_error_min = 0
        self.lambda_angle_error_max = 0
        self.impact_parameter = 0
        self.impact_parameter_error_min = 0
        self.impact_parameter_error_max = 0
        self.tzero_vr = 0
        self.tzero_vr_error_min = 0
        self.tzero_vr_error_max = 0
        self.k = 0
        self.k_error_min = 0
        self.k_error_max = 0
        self.temp_calculated = 0
        self.temp_calculated_error_min = 0
        self.temp_calculated_error_max = 0
        self.temp_measured = 0
        self.hot_point_lon = 0
        self.geometric_albedo = 0
        self.geometric_albedo_error_min = 0
        self.geometric_albedo_error_max = 0
        self.log_g = 0
        self.publication = 0
        self.detection_type = 0
        self.mass_detection_type = 0
        self.radius_detection_type = 0
        self.alternate_names = 0
        self.molecules = 0
        self.star_name = 0
        self.ra = 0
        self.dec = 0
        self.mag_v = 0
        self.mag_i = 0
        self.mag_j = 0
        self.mag_h = 0
        self.mag_k = 0
        self.star_distance = 0
        self.star_distance_error_min = 0
        self.star_distance_error_max = 0
        self.star_metallicity = 0
        self.star_metallicity_error_min = 0
        self.star_metallicity_error_max = 0
        self.star_mass = 0
        self.star_mass_error_min = 0
        self.star_mass_error_max = 0
        self.star_radius = 0
        self.star_radius_error_min = 0
        self.star_radius_error_max = 0
        self.star_sp_type = 0
        self.star_age = 0
        self.star_age_error_min = 0
        self.star_age_error_max = 0
        self.star_teff = 0
        self.star_teff_error_min = 0
        self.star_teff_error_max = 0
        self.star_detected_disc = 0
        self.star_magnetic_field = 0
        self.star_alternate_names = 0

# access 1st database and read the data from the csv file
csv_raw = urllib.request.urlopen('http://www.exoplanet.eu/catalog/csv')

# convert the web file to string and split into rows (each row still a csv)
csv_str = str(csv_raw.read())
rows = csv_str.split('\\r\\n')

# remove the first row containing the table headers and create an empty array of exoplanet objects
rows.pop(0)
exoplanet_array = []

# split each row by comma and fill the empty exoplanet array with all assigned properties
# specific to the exoplanets.eu database
for i in range(len(rows) - 1):
    #columns = re.split('[,](?!\s)', rows[i])
    # using imported re to split by comma except when inside quotes
    # adapted from https://stackoverflow.com/questions/18893390/splitting-on-comma-outside-quotes
    columns = re.split(",(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)", rows[i])
    
    exoplanet_array.append(exoplanet(columns[0]))
    exoplanet_array[i].planet_status = columns[1]
    exoplanet_array[i].mass = columns[2]
    exoplanet_array[i].mass_error_min = columns[3]
    exoplanet_array[i].mass_error_max = columns[4]
    exoplanet_array[i].mass_sini = columns[5]
    exoplanet_array[i].mass_sini_error_min = columns[6]
    exoplanet_array[i].mass_sini_error_max = columns[7]
    exoplanet_array[i].radius = columns[8]
    exoplanet_array[i].radius_error_min = columns[9]
    exoplanet_array[i].radius_error_max = columns[10]
    exoplanet_array[i].orbital_period = columns[11]
    exoplanet_array[i].orbital_period_error_min = columns[12]
    exoplanet_array[i].orbital_period_error_max = columns[13]
    exoplanet_array[i].semi_major_axis = columns[14]
    exoplanet_array[i].semi_major_axis_error_min = columns[15]
    exoplanet_array[i].semi_major_axis_error_max = columns[16]
    exoplanet_array[i].eccentricity = columns[17]
    exoplanet_array[i].eccentricity_error_min = columns[18]
    exoplanet_array[i].eccentricity_error_max = columns[19]
    exoplanet_array[i].inclination = columns[20]
    exoplanet_array[i].inclination_error_min = columns[21]
    exoplanet_array[i].inclination_error_max = columns[22]
    exoplanet_array[i].angular_distance = columns[23]
    exoplanet_array[i].discovered = columns[24]
    exoplanet_array[i].updated = columns[25]
    exoplanet_array[i].omega = columns[26]
    exoplanet_array[i].omega_error_min = columns[27]
    exoplanet_array[i].omega_error_max = columns[28]
    exoplanet_array[i].tperi = columns[29]
    exoplanet_array[i].tperi_error_min = columns[30]
    exoplanet_array[i].tperi_error_max = columns[31]
    exoplanet_array[i].tconj = columns[32]
    exoplanet_array[i].tconj_error_min = columns[33]
    exoplanet_array[i].tconj_error_max = columns[34]
    exoplanet_array[i].tzero_tr = columns[35]
    exoplanet_array[i].tzero_tr_error_min = columns[36]
    exoplanet_array[i].tzero_tr_error_max = columns[37]
    exoplanet_array[i].tzero_tr_sec = columns[38]
    exoplanet_array[i].tzero_tr_sec_error_min = columns[39]
    exoplanet_array[i].tzero_tr_sec_error_max = columns[40]
    exoplanet_array[i].lambda_angle = columns[41]
    exoplanet_array[i].lambda_angle_error_min = columns[42]
    exoplanet_array[i].lambda_angle_error_max = columns[43]
    exoplanet_array[i].impact_parameter = columns[44]
    exoplanet_array[i].impact_parameter_error_min = columns[45]
    exoplanet_array[i].impact_parameter_error_max = columns[46]
    exoplanet_array[i].tzero_vr = columns[47]
    exoplanet_array[i].tzero_vr_error_min = columns[48]
    exoplanet_array[i].tzero_vr_error_max = columns[49]
    exoplanet_array[i].k = columns[50]
    exoplanet_array[i].k_error_min = columns[51]
    exoplanet_array[i].k_error_max = columns[52]
    exoplanet_array[i].temp_calculated = columns[53]
    exoplanet_array[i].temp_calculated_error_min = columns[54]
    exoplanet_array[i].temp_calculated_error_max = columns[55]
    exoplanet_array[i].temp_measured = columns[56]
    exoplanet_array[i].hot_point_lon = columns[57]
    exoplanet_array[i].geometric_albedo = columns[58]
    exoplanet_array[i].geometric_albedo_error_min = columns[59]
    exoplanet_array[i].geometric_albedo_error_max = columns[60]
    exoplanet_array[i].log_g = columns[61]
    exoplanet_array[i].publication = columns[62]
    exoplanet_array[i].detection_type = columns[63]
    exoplanet_array[i].mass_detection_type = columns[64]
    exoplanet_array[i].radius_detection_type = columns[65]
    exoplanet_array[i].alternate_names = columns[66]
    exoplanet_array[i].molecules = columns[67]
    exoplanet_array[i].star_name = columns[68]
    exoplanet_array[i].ra = columns[69]
    exoplanet_array[i].dec = columns[70]
    exoplanet_array[i].mag_v = columns[71]
    exoplanet_array[i].mag_i = columns[72]
    exoplanet_array[i].mag_j = columns[73]
    exoplanet_array[i].mag_h = columns[74]
    exoplanet_array[i].mag_k = columns[75]
    exoplanet_array[i].star_distance = columns[76]
    exoplanet_array[i].star_distance_error_min = columns[77]
    exoplanet_array[i].star_distance_error_max = columns[78]
    exoplanet_array[i].star_metallicity = columns[79]
    exoplanet_array[i].star_metallicity_error_min = columns[80]
    exoplanet_array[i].star_metallicity_error_max = columns[81]
    exoplanet_array[i].star_mass = columns[82]
    exoplanet_array[i].star_mass_error_min = columns[83]
    exoplanet_array[i].star_mass_error_max = columns[84]
    exoplanet_array[i].star_radius = columns[85]
    exoplanet_array[i].star_radius_error_min = columns[86]
    exoplanet_array[i].star_radius_error_max = columns[87]
    exoplanet_array[i].star_sp_type = columns[88]
    exoplanet_array[i].star_age = columns[89]
    exoplanet_array[i].star_age_error_min = columns[90]
    exoplanet_array[i].star_age_error_max = columns[91]
    exoplanet_array[i].star_teff = columns[92]
    exoplanet_array[i].star_teff_error_min = columns[93]
    exoplanet_array[i].star_teff_error_max = columns[94]
    exoplanet_array[i].star_detected_disc = columns[95]
    exoplanet_array[i].star_magnetic_field = columns[96]
    exoplanet_array[i].star_alternate_names = columns[97]

# for testing
for i in range(len(exoplanet_array)):
    x = exoplanet_array[i].star_alternate_names
    if x != "": print(x)
