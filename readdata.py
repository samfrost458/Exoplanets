# a script to read in the data
import urllib
import re

import Exoplanet
import Telescope
import htmlwrite

# access database and read the data from the csv file
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
    exoplanet_array.append(Exoplanet.Exoplanet())
    for key, val in Exoplanet.Exoplanet.property_dict.items():
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
    tele_array.append(Telescope.Telescope())
    for key, val in Telescope.Telescope.property_dict.items():
        try:
            # attempt to float
            setattr(tele_array[i], key, float(columns[val]))
        except ValueError:
            # could not float, save as string instead
            setattr(tele_array[i], key, columns[val])
 
# try calculating metrics and ranking
for planet in exoplanet_array:
    metrics = []
    for telescope in tele_array:
        metrics.append(planet._decision_metric(telescope))
    
    # find highest nonzero value
    planet.metric = max(metrics)
    
# sort
exoplanet_array.sort(key=lambda x: x.decision_metric, reverse=True)
 
# truncate
init_length = len(exoplanet_array)
for i in range(init_length):
    # since we are popping the array, i will go out of range halfway through unless use a reverse index
    if exoplanet_array[init_length - 1 - i].decision_metric == 0: exoplanet_array.pop()

# print html page
print_properties = ['name', 'mass', 'radius', 'orbital_period', 'semi_major_axis', 'ra', 'dec'] # list of properties to display
path = '/Users/samfrost/Dropbox/Mphys_project/Thumbnails/' # directory for images
htmlwrite.htmlwrite(exoplanet_array, tele_array, 'htmltest', print_properties, path)