# 8/10/19 creating telescope coeff function, attempting to find source of discontinuity in decision metric
 
import os # for calling terminal
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import urllib
import re
 
# a class for exoplanets
class Exoplanet:
    
    # class attribute dictionary for planet properties (same as exoplanet.eu)
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
    def _decision_metric(self):
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
         
# function to convert degrees into ra hh:mm:ss string
def _deg_to_ra_hms(angle_deg, delim = ':', round_seconds = False):
    hours = (angle_deg - (angle_deg % 15)) / 15
    minutes = (angle_deg / 15 - hours) * 60
    seconds = (minutes % 1) * 60
     
    # decide to round the seconds (for display) or leave unrounded (internal)
    if round_seconds: return str(int(hours)) + delim + str(int(minutes)) + delim + str(round(seconds))
    else: return str(int(hours)) + delim + str(int(minutes)) + delim + str(seconds)
 
# function to convert degrees into dec deg:mm:ss string
def _deg_to_dec_dms(angle_deg, delim = ':', round_seconds = False):
    if angle_deg >= 0:
        deg = angle_deg - (angle_deg % 1)
        minutes = (angle_deg % 1) * 60
    else:
        deg = angle_deg + (-angle_deg % 1)
        minutes = (-angle_deg % 1) * 60
 
    seconds = (minutes % 1) * 60
     
    # decide to round the seconds (for display) or leave unrounded (internal)
    if round_seconds: return str(int(deg)) + delim + str(int(minutes)) + delim + str(round(seconds))
    else: return str(int(deg)) + delim + str(int(minutes)) + delim + str(seconds)
 
# function to get image urls from ESO database given ra dec in degrees
def _ra_dec_img_url(ra_deg, dec_deg, imageonly = True):
    # convert to hms/dms
    ra = _deg_to_ra_hms(ra_deg, delim = '%3A')
    dec = _deg_to_dec_dms(dec_deg, delim = '%3A')
     
    # find first url and read off html
    url_full = 'https://archive.eso.org/dss/dss/image?'
    url_full += 'ra=' + ra
    url_full += '&dec=' + dec
    url_full += '&equinox=J2000&name=&x=5&y=5&Sky-Survey=DSS1&mime-type=image%2Fgif&statsmode=WEBFORM'
     
    # search for true image url
    if imageonly:
        url_full_HTML = str( urllib.request.urlopen(url_full).read() )
        
        start = url_full_HTML.find('IMG SRC=/dss/dss') + 8
        end = url_full_HTML.find('.gif ALT') + 4
        url_pic = 'https://archive.eso.org'
        for i in range(start, end):
            url_pic += url_full_HTML[i]
         
        return url_pic
     
    # want full url, not just the image (which expires anyway)
    else: return url_full
 
# function to get commands printed for ESO batch system
def _ra_dec_batch_command(exoplanet_array, filename):
    # open command file
    batchfile = open(filename, 'w')
     
    for planet in exoplanet_array:
        # convert to hms/dms
        ra = _deg_to_ra_hms(planet.ra, delim = ' ')
        dec = _deg_to_dec_dms(planet.dec, delim = ' ')
         
        # write the command line
        batchfile.write(planet.name.replace(' ', '_') + ' ' + ra + ' ' + dec + ' 5 5\n') # 5by5 arcmin image
     
    # close file
    batchfile.close()
 
# a class for telescopes
class Telescope:
    # class attribute dictionary for telescope properties
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
 
# access 1st database and read the data from the csv file
planet_csv = urllib.request.urlopen('http://www.exoplanet.eu/catalog/csv')
 
# convert the web file to string and split into rows (each row still a csv)
planet_csv = str(planet_csv.read())
planet_rows = planet_csv.split('\\r\\n')
 
# remove the first row containing the table headers and create an empty array of exoplanet objects
planet_rows.pop(0)
exoplanet_array = []
 
# split each row by comma and fill the empty exoplanet array with all assigned properties
# specific to the exoplanets.eu database
for i in range(len(planet_rows) - 1):
    # using imported re to split by comma except when inside quotes
    columns = re.split(",(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)", planet_rows[i]) # adapted from https://stackoverflow.com/questions/18893390/splitting-on-comma-outside-quotes
     
    # append each new exoplanet and fill with values according to the dictionary
    exoplanet_array.append(Exoplanet())
    for key, val in Exoplanet.property_dict.items():
        try:
            # attempt to float
            setattr(exoplanet_array[i], key, float(columns[val]))
        except ValueError:
            # could not float, save as string instead
            setattr(exoplanet_array[i], key, columns[val])

# read in telescope data
tele_csv = open('/Users/samfrost/Dropbox/Mphys_project/Scope_Masters.csv')
tele_csv = tele_csv.read()
tele_csv = tele_csv.replace('NaN', '') # replace NaN in csv with blank space
tele_rows = tele_csv.split('\n')

# remove the first row containing the table headers and create an empty array of telescope objects
tele_rows.pop(0)
tele_array = []

# split each row by comma and fill the empty telescope array with all assigned properties
for i in range(len(tele_rows) - 1):
    # split by comma
    columns = re.split(',', tele_rows[i]) 
     
    # append each new telescope and fill with values according to the dictionary
    tele_array.append(Telescope())
    for key, val in Telescope.property_dict.items():
        try:
            # attempt to float
            setattr(tele_array[i], key, float(columns[val]))
        except ValueError:
            # could not float, save as string instead
            setattr(tele_array[i], key, columns[val])

# test print C_Ts
for telescope in tele_array:
    C_T = telescope._tele_coeff()
    print(C_T)
 
# try calculating metrics and ranking
for planet in exoplanet_array:
    planet._decision_metric()
# sort
exoplanet_array.sort(key=lambda x: x.decision_metric, reverse=True)
 
# truncate
init_length = len(exoplanet_array)
for i in range(init_length):
    # since we are popping the array, i will go out of range halfway through unless use a reverse index
    if exoplanet_array[init_length - 1 - i].decision_metric == 0: exoplanet_array.pop()

# plot
rank = []
metric = []
radius = []
mass = []
semi_major_axis = []

for i in range(len(exoplanet_array)):
    rank.append(i)
    metric.append(exoplanet_array[i].decision_metric)
    if exoplanet_array[i].radius == '': radius.append(0)
    else: radius.append(exoplanet_array[i].radius)
    #radius.append(exoplanet_array[i].radius)
    if exoplanet_array[i].mass == '': mass.append(0)
    else: mass.append(exoplanet_array[i].mass)
    semi_major_axis.append(exoplanet_array[i].semi_major_axis)

'''
for i in range(len(radius)):
   if radius[i] > 0.5:
       plt.scatter(rank[i], radius[i], color = 'b')
   else: plt.scatter(rank[i], radius[i], color = 'r')
'''

plt.scatter(mass, radius)
plt.show()
 
 
# fetch images
path = '/Users/samfrost/Dropbox/Mphys_project/Thumbnails/'
batchfilename = 'thumbnail_list.txt'
'''
# create command file
_ra_dec_batch_command(exoplanet_array, path + batchfilename)
# call dss1 executable
os.system("""
          cd %s
          ./dss1 -i %s
          """ %(path, batchfilename))
'''
 
# output to html file
outfile = open('htmltest.html','w')
outfile.write("<html>\n<head>\n")
 
# table style
style = "<style>"
style += """
table {
    border-collapse: collapse; /* collapsed borders */
    width: 100%; /* full width */
    }\n"""
style += "table, th, td { border: 1px solid black; }\n /* black borders */"
style += """
th {
    background-color: #4CAF50; /* green background for header */
    color: white; /* white header text */
    padding: 12px; /* header padding */
    font-size: 18px; /* font size */
    }\n"""
style += "tr:nth-child(even) { background-color: #f2f2f2; }\n /* every other row shaded grey */"
style += ".fixed_header thead th { position: sticky; top: 0; }\n /* frozen header */"
 
# search button styling by https://www.w3schools.com/howto/howto_js_filter_table.asp
style += """
#userInput {
  background-image: url('/css/searchicon.png'); /* Add a search icon to input */
  background-position: 10px 12px; /* Position the search icon */
  background-repeat: no-repeat; /* Do not repeat the icon image */
  width: 100%; /* Full-width */
  font-size: 16px; /* Increase font-size */
  padding: 12px 20px 12px 40px; /* Add some padding */
  border: 1px solid #ddd; /* Add a grey border */
  margin-bottom: 12px; /* Add some space below the input */
  }\n"""
style += "</style>\n"
outfile.write(style)
 
outfile.write("</head>\n<body>")
outfile.write("""<input type="text" id="userInput" onkeyup="searchFunction()" placeholder="Search for names..." title="Type in a name">""")
 
# create table
table = """<table class="fixed_header" id="planetTable">\n"""
 
print_properties = ['name', 'mass', 'radius', 'orbital_period', 'semi_major_axis', 'ra', 'dec']
 
# table header
table += "<thead>\n<tr>\n"
table += "<th>image</th>\n"
for key in print_properties:
    table += "<th>" + key.replace('_', ' ') + "</th>\n"
table += "</tr>\n</thead>\n"
 
# table body
table += "<tbody>\n"
for planet in exoplanet_array:
    table += "<tr>\n"
     
    # old url-based image collection - this expires so regular rerunning of python script is necessary unless we download and store the images
    #table += "<td><img src=\"" + _ra_dec_img_url(planet.ra, planet.dec) + "\" alt = \"%s\">\n" %planet.name + "</td>\n"
     
    # new system, still goes by url if necessary but stores the image for indefinite use (does not use batch tool since this only downloads FITS)
    foundImage = False
    path = '/Users/samfrost/Dropbox/Mphys_project/Thumbnails/'
    #url_pic = _ra_dec_img_url(planet.ra, planet.dec)
    url_full = _ra_dec_img_url(planet.ra, planet.dec, imageonly = False)
     
    # search for image file in directory
    for imagefilename in os.listdir(path):
        if imagefilename == planet.name.replace(' ', '_').replace('/', '_') + '.gif':
            # image found in file, put in table
            foundImage = True
            table += """<td><a href="%s"> <img src="%s%s" alt = "%s" style="border:0"></a></td>\n""" %(url_full, path, imagefilename, planet.name)
             
    # if no image found, download a new one
    if not foundImage:
        # download new image
        url_pic = _ra_dec_img_url(planet.ra, planet.dec)
        imagefilename = path + planet.name.replace(' ', '_').replace('/', '_') + '.gif'
        urllib.request.urlretrieve(url_pic, imagefilename)
         
        # file downloaded, put in table
        table += """<td><a href="%s"><img src="%s%s" alt = "%s" style="border:0"></a></td>\n""" %(url_full, path, imagefilename, planet.name)
     
    '''
    # new batch tool for image collection, images permanently downloaded already, downloads FITS (can't view image)
    foundImage = False
    # search for image file in directory
    for imagefilename in os.listdir(path):
        if imagefilename.startswith(planet.name.replace(' ', '_')) and imagefilename.endswith('.fits'):
            table += "<td><img src=\"" + path + imagefilename + "\" alt = \"%s\">\n" %planet.name + "</td>\n"
            foundImage = True
    # if no image found, download a new one
    if not foundImage:
        table += "<td></td>\n"
        os.system("""
          cd %s
          ./dss1
          %s
          """ %(path, batchfilename))
    '''
     
    # print all other properties
    for key in print_properties:
        # if printing the name, make it a hyperlink
        if key == 'name':
            link = "http://www.exoplanet.eu/catalog/%s/" %planet.name.replace(' ', '_')
            display = planet.name
            table += "<td> <a href=\"%s\">%s " %(link, display) + "</td>\n"
         
        # for ra and dec display in hms/dms
        elif key == 'ra':
            table += "<td>" + _deg_to_ra_hms(planet.ra, round_seconds = True) + "</td>\n"
        elif key == 'dec':
            table += "<td>" + _deg_to_dec_dms(planet.dec, round_seconds = True) + "</td>\n"
         
        # otherwise display as normal
        else:
            table += "<td>" + str(getattr(planet, key)) + "</td>\n"
     
    table += "</tr>\n"
 
table += "</tbody>\n</table>\n"
 
outfile.write(table)
 
# search function from https://www.w3schools.com/howto/howto_js_filter_table.asp
outfile.write("""
<script>
function searchFunction() {
  // Declare variables 
  var input, filter, table, tr, td, i, txtValue;
  input = document.getElementById("userInput");
  filter = input.value.toUpperCase();
  table = document.getElementById("planetTable");
  tr = table.getElementsByTagName("tr");
 
  // Loop through all table rows, and hide those who don't match the search query
  for (i = 0; i < tr.length; i++) {
    td = tr[i].getElementsByTagName("td")[1]; // element 1 for name column
    if (td) {
      txtValue = td.textContent || td.innerText;
      if (txtValue.toUpperCase().indexOf(filter) > -1) {
        tr[i].style.display = "";
      } else {
        tr[i].style.display = "none";
      }
    } 
  }
}
</script>""")
 
outfile.write("</body>\n</html>")
 
outfile.close()