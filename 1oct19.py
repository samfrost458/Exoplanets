# 1/10/19 messing with html and dictionaries, and plotting albedo vs radius etc.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import urllib
import re

# a class for exoplanets
class Exoplanet:
   
    # class attribute dictionary for planet properties (same as exoplanet.eu)
    property_dict = {
        'name'                       : 0,
        'planet_status'              : 1,
        'mass'                       : 2,
        'mass_error_min'             : 3,
        'mass_error_max'             : 4,
        'mass_sini'                  : 5,
        'mass_sini_error_min'        : 6,
        'mass_sini_error_max'        : 7,
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
        'angular_distance'           : 23,
        'discovered'                 : 24,
        'updated'                    : 25,
        'omega'                      : 26,
        'omega_error_min'            : 27,
        'omega_error_max'            : 28,
        'tperi'                      : 29,
        'tperi_error_min'            : 30,
        'tperi_error_max'            : 31,
        'tconj'                      : 32,
        'tconj_error_min'            : 33,
        'tconj_error_max'            : 34,
        'tzero_tr'                   : 35,
        'tzero_tr_error_min'         : 36,
        'tzero_tr_error_max'         : 37,
        'tzero_tr_sec'               : 38,
        'tzero_tr_sec_error_min'     : 39,
        'tzero_tr_sec_error_max'     : 40,
        'lambda_angle'               : 41,
        'lambda_angle_error_min'     : 42,
        'lambda_angle_error_max'     : 43,
        'impact_parameter'           : 44,
        'impact_parameter_error_min' : 45,
        'impact_parameter_error_max' : 46,
        'tzero_vr'                   : 47,
        'tzero_vr_error_min'         : 48,
        'tzero_vr_error_max'         : 49,
        'k'                          : 50,
        'k_error_min'                : 51,
        'k_error_max'                : 52,
        'temp_calculated'            : 53,
        'temp_calculated_error_min'  : 54,
        'temp_calculated_error_max'  : 55,
        'temp_measured'              : 56,
        'hot_point_lon'              : 57,
        'geometric_albedo'           : 58,
        'geometric_albedo_error_min' : 59,
        'geometric_albedo_error_max' : 60,
        'log_g'                      : 61,
        'publication'                : 62,
        'detection_type'             : 63,
        'mass_detection_type'        : 64,
        'radius_detection_type'      : 65,
        'alternate_names'            : 66,
        'molecules'                  : 67,
        'star_name'                  : 68,
        'ra'                         : 69,
        'dec'                        : 70,
        'mag_v'                      : 71,
        'mag_i'                      : 72,
        'mag_j'                      : 73,
        'mag_h'                      : 74,
        'mag_k'                      : 75,
        'star_distance'              : 76,
        'star_distance_error_min'    : 77,
        'star_distance_error_max'    : 78,
        'star_metallicity'           : 79,
        'star_metallicity_error_min' : 80,
        'star_metallicity_error_max' : 81,
        'star_mass'                  : 82,
        'star_mass_error_min'        : 83,
        'star_mass_error_max'        : 84,
        'star_radius'                : 85,
        'star_radius_error_min'      : 86,
        'star_radius_error_max'      : 87,
        'star_sp_type'               : 88,
        'star_age'                   : 89,
        'star_age_error_min'         : 90,
        'star_age_error_max'         : 91,
        'star_teff'                  : 92,
        'star_teff_error_min'        : 93,
        'star_teff_error_max'        : 94,
        'star_detected_disc'         : 95,
        'star_magnetic_field'        : 96,
        'star_alternate_names'       : 97
    }
    
    # initialiser setting all properties in dictionary to 0
    def __init__(self):
        # initialise all values in the dictionary to 0
        for key, val in self.property_dict.items():
            setattr(self, key, 0)
        
        # initialise decision metric score as 0
        self.decision_metric = 0
    
    # to do: telescope database, telescope coefficient, decision metric function


'''
# Pandas database alternate approach - not used yet

# access the databases and put into panda dataframe
df_exoplanets_eu = pd.read_csv('http://www.exoplanet.eu/catalog/csv', delimiter = ',')
# see https://exoplanetarchive.ipac.caltech.edu/docs/program_interfaces.html for details on which url to request
df_exoplanetarchive = pd.read_csv('https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=exoplanets', delimiter = ',')
#print(df_exoplanetarchive)
'''
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
    # using imported re to split by comma except when inside quotes
    columns = re.split(",(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)", rows[i]) # adapted from https://stackoverflow.com/questions/18893390/splitting-on-comma-outside-quotes
    
    # append each new exoplanet and fill with values according to the dictionary
    exoplanet_array.append(Exoplanet())
    for key, val in Exoplanet.property_dict.items():
        setattr(exoplanet_array[i], key, columns[val])

'''
# plotting albedo vs radius for james
exoplanet_array.sort(key=lambda x: x.geometric_albedo)
albedo = [0.06, 0.75, 0.30, 0.29, 0.52, 0.47, 0.51, 0.41,          0.03, 0.038, 0.4, 0.5, 0.58]
radius = [0.034, 0.085, 0.089, 0.48, 1.0, 0.84, 0.36, 0.35]
stellarDistance = [0.386, 0.72, 1.0, 1.52, 5.19, 9.55, 19.1, 29.7, 143.0, 47.0, 19.3, 14.7, 173.0]
for i in range(len(exoplanet_array)):
    x = exoplanet_array[i].geometric_albedo
    if x != "":
        print(x, exoplanet_array[i].star_distance)
        #albedo.append(float(x))
        #radius.append(float(exoplanet_array[i].star_distance))
plt.scatter(stellarDistance, albedo)
plt.xlabel("Stellar distance")
plt.ylabel("Albedo")
plt.show()
'''

# testing output to html file
outfile = open('htmltest.html','w')
outfile.write("<html>")

# page header
outfile.write("<head>\nTest print of database\n")

# table style
style = "<style>\n"
style += "table { border-collapse: collapse; }\n" # collapsed borders
style += "table, th, td { border: 1px solid black; }\n" # black borders
style += "th { background-color: #4CAF50; color: white; }\n" # header in white writing against green background
style += "tr:nth-child(even) {background-color: #f2f2f2;}\n" # every other row shaded grey
style += ".fixed_header thead th { position: sticky; top: 0; }\n" # frozen header
style += "</style>\n"
outfile.write(style)

outfile.write("</head>\n<body>")

# create table
table = "<table class=\"fixed_header\">\n"

# table header
table += "<thead>\n<tr>\n"
for key, val in Exoplanet.property_dict.items():
    table += "<th>" + key + "</th>\n"
table += "</tr>\n</thead>\n"

# table body
table += "<tbody>\n"
for planet in exoplanet_array:
    table += "<tr>\n"
    
    for key, val in Exoplanet.property_dict.items():
        # if printing the name, make it a hyperlink
        if key == 'name':
            link = "http://www.exoplanet.eu/catalog/%s/" %getattr(planet, key).replace(' ', '_')
            display = getattr(planet, key)
            table += "<td> <a href=\"%s\">%s " %(link, display) + "</td>\n"
        else:
            table += "<td>" + getattr(planet, key) + "</td>\n"
    
    table += "</tr>\n"

table += "</tbody>\n</table>\n"

outfile.write(table)
outfile.write("</body>\n</html>")

outfile.close()