# a document with assorted functions for reading and writing web data

# standard python imports
import urllib
import re
import numpy as np
import copy
import os
import pickle
import time
import datetime
import json

# recursive getattr and setattr tools
import functools
def rsetattr(obj, attr, val):
    '''recursive setattr function'''
    pre, _, post = attr.rpartition('.')
    return setattr(rgetattr(obj, pre) if pre else obj, post, val)
def rgetattr(obj, attr, *args):
    '''recursive setattr function'''
    def _getattr(obj, attr):
        return getattr(obj, attr, *args)
    return functools.reduce(_getattr, [obj] + attr.split('.'))

# graph plotting imports
from bokeh.plotting import figure
from bokeh.models import (ColumnDataSource, Legend, ColorBar, LogTicker, Whisker,
    PanTool, BoxZoomTool, BoxSelectTool, LassoSelectTool, WheelZoomTool, TapTool, OpenURL, UndoTool, RedoTool, CrosshairTool, HoverTool, SaveTool, ResetTool)
from bokeh.palettes import Category10, Turbo256
from bokeh.transform import factor_cmap, linear_cmap, log_cmap
from bokeh.layouts import column, row, gridplot
from bokeh.resources import CDN
from bokeh.embed import file_html
from bokeh.io import export_png

# observability imports
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astroplan import Observer, FixedTarget, AltitudeConstraint, AtNightConstraint, months_observable, is_observable
from astropy.time import Time

# custom exoplanet and telescope classes
from Exoplanet import Exoplanet, Star
from Telescope import Telescope

def read_planet_telescope_data_csv(planet_src, tele_src, path_to_static, onlyConfirmed, to_remove=[]):
    '''
    Function to read in the planet and telescope data, determine if calculable, truncate and return arrays for both, and check if images are present in the static folder

    Parameters
    ----------
    planet_src : str
        A url for the csv containing planet data
    tele_src : str
        A full file path to the csv containing telescope data
    path_to_static : str
        A file path to the static folder, for checking and saving images, saving pickle files as a backup and saving pickle files to check observability

    Returns
    -------
    exoplanet_array : list of Exoplanet
    tele_array      : list of Telescope
    plottingsource  : dict
        A dictionary with Exoplanet properties as keys and a list as values, ready to be converted into a Bokeh plotting source
    '''

    try:
        ################
        # exoplanet read in

        start = time.perf_counter()
        print('Updated data ' + str(datetime.datetime.now()))

        # access database and read the data from the csv file
        if not onlyConfirmed:
            try:
                # submit form data to get unconfirmed planets too
                values = {'Confirmed' : 'on',
                          'Candidate' : 'on',
                          'Other'     : 'on',
                          'query_f'   : '',
                          'search'    : ''}
                data = urllib.parse.urlencode(values).encode('utf-8')
                url = planet_src + '?status_1=on&status_2=on&status_4=on&query_f=&search=/'
                req = urllib.request.Request(url, data=data)
                planet_csv = urllib.request.urlopen(req)
            except:
                # didn't work, just get confirmed ones by default
                planet_csv = urllib.request.urlopen(planet_src)
                onlyConfirmed = True
        else: planet_csv = urllib.request.urlopen(planet_src)
         
        # convert the web file to string and split into rows (each row still a csv)
        planet_csv = str(planet_csv.read())
        planet_rows = planet_csv.split('\\r\\n')
         
        # remove the first row containing the table headers and the last blank row and create an empty array of exoplanet objects
        planet_rows.pop(0)
        planet_rows.pop()
        exoplanet_array = []
        host_star_array = []
         
        # split each row by comma and fill the empty exoplanet array with all assigned properties
        # specific to the exoplanets.eu database, Exoplanet.property_dict has the same headers in the same order
        for row in planet_rows:
            # split by comma (except when inside quotes)
            columns = re.split(",(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)", row)
            
            planet_col = columns[:68]
            host_col   = columns[68:]
            prop_units = list(Exoplanet.property_dict.items()) + [('host.' + prop, val) for prop, val in Star.property_dict.items()]
            
            # planet object (and host object contained inside)
            planet = Exoplanet()
             
            # add info
            for val, (prop, unit) in zip(columns, prop_units):
                try:
                    # attempt to float
                    if not np.isnan(float(val)) and not np.isinf(float(val)):
                        rsetattr(planet, prop, float(val))
                        
                        # years will should be datetime timestamps, for plotting sake
                        if unit == 'year': rsetattr(planet, prop, datetime.datetime.strptime(val, '%Y').timestamp()*1000)
                    
                    else:
                        # could float but the value is nan or inf, set a blank instead
                        rsetattr(planet, prop, '')
                
                # could not float         
                except ValueError:
                    # dates should be datetime timestamps, for plotting sake
                    if unit == 'date': rsetattr(planet, prop, datetime.datetime.strptime(val, '%Y-%m-%d').timestamp()*1000)

                    # yes/no strings converted to bool
                    elif val == 'yes': rsetattr(planet, prop, True)
                    elif val == 'no' : rsetattr(planet, prop, False)
                    
                    # save as string instead
                    else: rsetattr(planet, prop, val.replace('"', ''))
            
            # append to lists
            exoplanet_array.append(planet)
            host_star_array.append(planet.host)

        # manually remove erroneous members of the list
        for planet in exoplanet_array:
            if planet.name in to_remove:
                exoplanet_array.remove(planet)
                host_star_array.remove(planet.host)

        ####################
        # telescope read in
        
        # read in telescope data
        tele_csv = open(tele_src, 'r')
        tele_csv = tele_csv.read()
        tele_csv = tele_csv.replace('NaN', '') # replace NaN in csv with blank space
        tele_rows = tele_csv.split('\n')

        # remove the first row containing the table headers and last blank row and create an empty array of telescope objects
        tele_rows.pop(0)
        tele_rows.pop()
        tele_array = []

        # split each row by comma and fill the empty telescope array with all assigned properties
        for row in tele_rows:
            # split by comma
            columns = re.split(',', row)

            # telescope object
            telescope = Telescope()
             
            # append each new telescope and fill with values according to the dictionary
            for prop, index in Telescope.property_dict.items():
                val = columns[index]
                try:
                    # spectral bins and half well are ints
                    if prop == 'spectral_bins' or prop == 'Half_Well':
                        setattr(telescope, prop, int(val))
                    
                    # otherwise attempt to float
                    else: setattr(telescope, prop, float(val))
                
                # could not float, save as string instead
                except ValueError: setattr(telescope, prop, val)

            # append telescope
            tele_array.append(telescope)

        ###################
        # other jobs to do

        # determine if calculable
        for planet in exoplanet_array:
            planet._metric_calculable_()
        for telescope in tele_array:
            telescope._C_T_calculable_()

        # check images of each planet exist in the static folder, if not download it from the internet
        new_names = []
        for planet in exoplanet_array:
            if not os.path.isfile(path_to_static + planet.name.replace(' ', '_').replace('/', '_') + '.gif'):
                url_pic = ra_dec_img_url(planet.host.ra, planet.host.dec)
                filename = path_to_static + planet.name.replace(' ', '_').replace('/', '_') + '.gif'
                urllib.request.urlretrieve(url_pic, filename)
                new_names.append(planet.name)

        # make bokeh plotting source
        source_dict = plotting_source_dict(exoplanet_array)

        ###################
        # save progress
        
        # save in pickle just in case
        with open(path_to_static + 'exoplanet_array.pkl', 'wb') as f:
            pickle.dump(exoplanet_array, f)
        with open(path_to_static + 'tele_array.pkl', 'wb') as f:
            pickle.dump(tele_array, f)
        with open(path_to_static + 'source_dict.pkl', 'wb') as f:
            pickle.dump(source_dict, f)

        end = time.perf_counter()
        print('Data read-in took ' + str(round(end - start, 3)) + ' seconds for %i planets (reduced to %i) and %i telescopes' %(len(planet_rows), len(exoplanet_array), len(tele_array)))
        if new_names:
            print('New planets since last update:')
            for name in new_names: print(name)

    except:
        # maybe a timeout error occured, just reload from stored pickle arrays instead
        with open(path_to_static + 'exoplanet_array.pkl', 'rb') as f:
            exoplanet_array = pickle.load(f)
        with open(path_to_static + 'tele_array.pkl', 'rb') as f:
            tele_array = pickle.load(f)
        with open(path_to_static + 'source_dict.pkl', 'rb') as f:
            source_dict = pickle.load(f)

        print('Could not load new data, loading data from last update instead')
        print('%i planets and %i telescopes' %(len(exoplanet_array), len(tele_array)))

    return (exoplanet_array, tele_array, source_dict)

def find_observability(exoplanet_array, default_lookahead, longitude, latitude, elevation, start_date, end_date, min_altitude, max_altitude, resolution, flash, display):
    '''Function to calculate'''
    '''# find timezone
    timezone_name = tf.timezone_at(lng=longitude, lat=latitude)
    delta_degree = 1
    while timezone_name is None:
        delta_degree += 1
        timezone_name = tf.closest_timezone_at(lng=longitude, lat=latitude, delta_degree=delta_degree)
    if flash: flash.append('Timezone   : ' + timezone_name)
    print(timezone_name)'''
    
    # create Observer instance
    if elevation is None: elevation = 0
    custom_location = Observer(longitude=longitude*u.deg, latitude=latitude*u.deg, elevation=elevation*u.m)
    
    # make an array of FixedTarget instances
    target_table = []
    for planet in exoplanet_array: target_table.append((planet.name, planet.host.ra, planet.host.dec))
    targets = [FixedTarget(coord=SkyCoord(ra=ra*u.deg, dec=dec*u.deg), name=name)
               for name, ra, dec in target_table]

    # time range
    if start_date is not None: start_date = Time(datetime.combine(start_date, datetime.min.time()), out_subfmt='date')
    else: start_date = Time(datetime.combine(datetime.now(), datetime.min.time()), out_subfmt='date')
    if end_date is not None: end_date = Time(datetime.combine(end_date, datetime.min.time()), out_subfmt='date')
    else: end_date = start_date + default_lookahead*u.day
    time_range = Time([start_date, end_date])
    #time_range_int = # no. days as an int
    if display:
        flash.append('Start date : ' + str(start_date.iso))
        flash.append('End date   : ' + str(end_date.iso))

    # altitude and night time constraints
    if min_altitude is None: min_altitude = 0
    if max_altitude is None: max_altitude = 90
    constraints = [AltitudeConstraint(min_altitude*u.deg, max_altitude*u.deg), AtNightConstraint.twilight_astronomical()]
    
    # calculate if observable within lookahead time from now
    if resolution is None: resolution = 0.5

    # new route, faster by not bothering to recheck planets already found to be observable
    # UPDATE, ASTROPY IS ALREADY OPTIMISED FOR THIS, THIS TAKES 10 TIMES LONGER, LEAVING IN FOR DEMONSTRATION PURPOSES ONLY (commented out lines in 'old route' are also only for demonstrating the difference)
    # start = time.perf_counter()
    # ever_observable = [False] * len(exoplanet_array)
    # day = 0
    # while start_date + day*u.day != end_date:
    #     print(str(start_date + day*u.day) + ', ' + str(end_date))
    #     for i in range(len(ever_observable)):
    #         if ever_observable[i] == False:
    #             ever_observable[i] = is_observable(constraints, custom_location, [targets[i]], time_range=Time([start_date+day*u.day, start_date+(day+1)*u.day]), time_grid_resolution=resolution*u.hour)[0]
    #             print(str(day) + ', ' + str(i) + ', ' + targets[i].name + ', ' + str(ever_observable[i]))
    #     day += 1
    # end = time.perf_counter()
    # print(str(end - start) + ' seconds fast')

    # new route, faster by not bothering to recheck planets already found to be observable
    # UPDATE, ASTROPY IS ALREADY OPTIMISED FOR THIS, THIS TAKES 10 TIMES LONGER, LEAVING IN FOR DEMONSTRATION PURPOSES ONLY (commented out lines in 'old route' are also only for demonstrating the difference)
    # start = time.perf_counter()
    # ever_observable = [False] * len(exoplanet_array)
    # day = 0
    # for i in range(len(ever_observable)):
    #     for day in range(30):
    #         print(str(start_date + day*u.day) + ', ' + str(end_date))
        
    #         if ever_observable[i] == False:
    #             ever_observable[i] = is_observable(constraints, custom_location, [targets[i]], time_range=Time([start_date+day*u.day, start_date+(day+1)*u.day]), time_grid_resolution=resolution*u.hour)[0]
    #             print(str(day) + ', ' + str(i) + ', ' + targets[i].name + ', ' + str(ever_observable[i]))
    #         else: break
    # end = time.perf_counter()
    # print(str(end - start) + ' seconds fast')

    # trying to split array into chunks so we can say the task has been completed to the nearest percent
    # tosplit_target_array = copy.deepcopy(targets)
    # start = time.perf_counter()
    # target_chunks = np.array_split(tosplit_target_array, 10) # split targets into 10 parts (not 100, this makes the calculation time too long)
    # for i in range(len(target_chunks)):
    #     is_observable(constraints, custom_location, target_chunks[i].tolist(), time_range=time_range, time_grid_resolution=resolution*u.hour)
    #     #print('%i calculated' %i)
    # end = time.perf_counter()
    # print(str(end - start) + ' seconds fast')

    # old route, slow, commented out lines are for demonstrating that it is better than a manual method by a factor of 10
    start = time.perf_counter()
    ever_observable = is_observable(constraints, custom_location, targets, time_range=time_range, time_grid_resolution=resolution*u.hour)
    end = time.perf_counter()
    print(str(end - start) + ' seconds slow')

    for i in range(len(exoplanet_array)):
        exoplanet_array[i].observable = ever_observable[i]
    if display: flash.append('Timing resolution = {} hrs'.format(resolution))

def sort_and_truncate_by_property(exoplanet_array, property = 'decision_metric'):
    ''' Function to sort and truncate an array of planets by checking if the property value for each is 0'''
    
    # sort
    exoplanet_array.sort(key=lambda x: rgetattr(x, property), reverse=True)

    # truncate
    init_length = len(exoplanet_array)
    for i in range(init_length):
        # since we are popping the array, i will go out of range halfway through unless use a reverse index
        if rgetattr(exoplanet_array[init_length - 1 - i], property) == 0: exoplanet_array.pop()

def random_planet_selector(sorted_exoplanet_array):
    '''Function to select a random planet from an array, subject to a custom probability distribution and a cutoff'''

    # pass by value so poping later on does not change array outside function
    exoplanet_array = copy.deepcopy(sorted_exoplanet_array)

    # find cumulative scores (assuming array is sorted by metric already)
    cumulativeScore = []
    cumulativeScoreTotal = 0
    for planet in exoplanet_array:
        cumulativeScoreTotal += planet.decision_metric
        cumulativeScore.append(cumulativeScoreTotal)

    # truncate by cutoff percentage of total cumulative
    for i in range(len(cumulativeScore)):
        if cumulativeScore[i] > 0.95 * cumulativeScoreTotal:
            init_length = len(exoplanet_array)
            del exoplanet_array[-(init_length - i):]
            break
    
    # make custom probability distribution
    y = []
    sumY = 0
    for planet in exoplanet_array: y.append(np.log(planet.decision_metric))
    for i in y: sumY += i
    p = np.divide(y, sumY)

    # pick one planet
    randomRank = np.random.choice(a=len(exoplanet_array), p=p)

    # return planet
    return exoplanet_array[randomRank]

def planet_array_to_json_array(exoplanet_array, delim):
    dicts = []
    extra_props = ['decision_metric', 'filter', 't_exp']
    for planet in exoplanet_array:
        dict_ = {delim + prop : getattr(planet.host, prop) for prop in Star.property_dict.keys()}
        dict_.update({prop : getattr(planet, prop) for prop in Exoplanet.property_dict.keys()})
        dict_.update({prop : getattr(planet, prop) for prop in extra_props})
        dicts.append(dict_)
    return json.dumps(dicts)

def deg_to_ra_hms(angle_deg, delim = ':', round_seconds = False):
    '''Function to convert degrees into ra hh:mm:ss string '''
    
    hours = (angle_deg - (angle_deg % 15)) / 15
    minutes = (angle_deg / 15 - hours) * 60
    seconds = (minutes % 1) * 60
     
    # decide to round the seconds (for display) or leave unrounded (internal)
    if round_seconds: return str(int(hours)) + delim + str(int(minutes)) + delim + str(round(seconds))
    else: return str(int(hours)) + delim + str(int(minutes)) + delim + str(seconds)

def deg_to_dec_dms(angle_deg, delim = ':', round_seconds = False):
    '''Function to convert degrees into dec deg:mm:ss string'''
    
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

def ra_dec_img_url(ra_deg, dec_deg, imageonly = True):
    '''Function to get image urls from ESO database, given ra dec in degrees'''

    # convert to hms/dms
    ra = deg_to_ra_hms(ra_deg, delim = '%3A')
    dec = deg_to_dec_dms(dec_deg, delim = '%3A')
     
    # find first url and read off html
    url_full  = 'https://archive.eso.org/dss/dss/image?'
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

def round_sig_figs(number, significant):
    # try for float
    try: return np.format_float_positional(number, significant, unique=False, fractional=False, trim='k')
    # if it doesn't work number is probably a string
    except: return number

def categorical_legend(categories, counts, title='', palette=Category10[10], fig=None, html=False, width=215, height=295):
    if fig is None: fig = figure(toolbar_location=None, width=width, height=height, outline_line_color=None, background_fill_alpha=0.7, border_fill_alpha=0.7,)
    legend = Legend(items=[
        (
            category + ' - %i' %count,
            [fig.scatter(x=[np.nan], y=[np.nan], color=color)]
        ) for (category, count, color) in zip(categories, counts, palette)
    ], background_fill_alpha=0, border_line_alpha=0,)
    if title != '': legend.title = title
    fig.add_layout(legend)

    if html: return file_html(fig, CDN)
    else: return legend

def plotting_source_dict(exoplanet_array):
    '''Function to make bokeh plotting source dictionary'''

    source_dict = {}
    p_properties = list(Exoplanet.property_dict.keys())
    s_properties = ['host.' + prop for prop in Star.property_dict.keys()]
    
    for prop in p_properties + s_properties:
        arr = []
        for planet in exoplanet_array:    
            # errorbars have absolute values
            if '_error_min' in prop:
                eval_min = rgetattr(planet, prop)
                baseval = rgetattr(planet, prop.replace('_error_min', ''))
                try: val = baseval - eval_min
                except: val = baseval
            elif '_error_max' in prop:
                eval_max = rgetattr(planet, prop)
                baseval = rgetattr(planet, prop.replace('_error_max', ''))
                try: val = baseval + eval_max
                except: val = baseval
            # otherwise all other properties have the same values
            else: val = rgetattr(planet, prop)
            
            if val == '': val = np.nan # otherwise bokeh tries to plot blank values with 0
            arr.append(val)

        source_dict.update({prop: arr})
    return source_dict

def linked_scatter_graph(master_source_dict, xproperties, yproperties, cproperties, sproperties, xlogs, ylogs, clogs, slogs, errorbar_xs, errorbar_ys, x_hists, y_hists, tooltip_list):
    '''Funtion to generate html for an interactive bokeh scatter/errorbar graph'''

    # copy source_dict to avoid overmapping
    source_dict = copy.deepcopy(master_source_dict)

    # size mapping (do before making the data source)
    for i, (sproperty, slog) in enumerate(zip(sproperties, slogs)):
        if sproperty == 'None':
            source_dict.update( {'smap%i' %i: [4]*len(source_dict['name'])} ) # size 4 is default
        else:
            arr = source_dict[sproperty]
            if slog: arr = np.log(arr)
            source_dict.update( {'smap%i' %i: np.nan_to_num( np.multiply(np.subtract(arr, np.nanmin(arr)), 15/np.nanmax(arr)) ).tolist()} )

    # make data source
    source = ColumnDataSource(data=source_dict)
    
    plots = []    
    for i in range(len(xproperties)): # don't use i for loops inside this 
        xproperty, yproperty, cproperty, sproperty, xlog, ylog, clog, slog, errorbar_x, errorbar_y, x_hist, y_hist = xproperties[i], yproperties[i], cproperties[i], sproperties[i], xlogs[i], ylogs[i], clogs[i], slogs[i], errorbar_xs[i], errorbar_ys[i], x_hists[i], y_hists[i]

        if xproperty != 'None' and yproperty != 'None':
            # units
            units = [''] * 4
            properties = [xproperty, yproperty, cproperty, sproperty]
            for j, prop in enumerate(properties): # must use enumerate, not zip
                if prop != 'None':
                    if 'host.' in prop: units[j] = Star.property_dict[prop.replace('host.', '')]
                    else: units[j] = Exoplanet.property_dict[prop]
            (xunit, yunit, cunit, sunit) = (unit for unit in units)
            formatters = [''] * 4 # for datetime formatting in tooltips
            for j, unit in enumerate(units):
                if unit == 'date': formatters[j] = '{%F}'
                if unit == 'year': formatters[j] = '{%Y}'

            # axis scale
            if xlog: x_axis_type = 'log'
            else: x_axis_type = 'linear'
            if ylog: y_axis_type = 'log'
            else: y_axis_type = 'linear'
            if xunit == 'date' or xunit == 'year': x_axis_type='datetime'
            if yunit == 'date' or yunit == 'year': y_axis_type='datetime'

            # create tooltip table
            tooltips = [
                ('(%s, %s)' %(xproperty.replace('.', '_'), yproperty.replace('.', '_')), '(@{%s}%s, @{%s}%s)' %(xproperty, formatters[0], yproperty, formatters[1])),
                ('(%s, %s)' %(cproperty.replace('.', '_'), sproperty.replace('.', '_')), '(@{%s}%s, @{%s}%s)' %(cproperty, formatters[2], sproperty, formatters[3])),]
            for prop in tooltip_list: tooltips.append( (prop, '@%s' %prop) )

            # create plot figure
            plot = figure(
                x_axis_type=x_axis_type,
                y_axis_type=y_axis_type,
                sizing_mode='stretch_both',
                tools=[PanTool(), BoxZoomTool(), BoxSelectTool(), LassoSelectTool(), WheelZoomTool(), TapTool(), UndoTool(), RedoTool(), CrosshairTool(), SaveTool(), ResetTool(),
                    HoverTool(tooltips=tooltips, formatters={'@{discovered}': 'datetime', '@{updated}': 'datetime'})],
                toolbar_location='above',
                # title=str(datetime.date.today()),
                background_fill_alpha=0.7,
                border_fill_alpha=0.7,
            )

            # invert scales for magnitudes
            if 'mag_' in xproperty: plot.x_range.flipped = True
            if 'mag_' in yproperty: plot.y_range.flipped = True

            # link points to exoplanet.eu
            url = 'http://www.exoplanet.eu/catalog/@name'
            taptool = plot.select(type=TapTool)
            taptool.callback = OpenURL(url=url.replace(' ', '_'))

            # label axes
            plot.xaxis.axis_label = xproperty.replace('_', ' ').replace('.', ' ').capitalize() + ' (%s)'%xunit
            plot.yaxis.axis_label = yproperty.replace('_', ' ').replace('.', ' ').capitalize() + ' (%s)'%yunit
            
            # font sizes
            fontsize = '12pt'
            plot.xaxis.axis_label_text_font_size = fontsize
            plot.yaxis.axis_label_text_font_size = fontsize

            # color mapping
            if cproperty == 'None': cmap = Category10[10][0]
            elif cproperty == 'detection_type':
                # detection type is categorical so needs special treatment
                detection_methods = ['Primary Transit', 'Radial Velocity', 'Imaging', 'Microlensing', 'Timing', 'Astrometry', 'TTV', 'Default', 'Primary Transit, TTV', 'Controversial']
                cmap = factor_cmap('detection_type', palette=Category10[10], factors=detection_methods)
                
                # counts
                counts = {method: 0 for method in detection_methods}
                for j, method in enumerate(source_dict['detection_type']):
                    if np.isfinite(source_dict[xproperty][j]) and np.isfinite(source_dict[yproperty][j]): counts[method] += 1
                
                # need to create legend manually so we can place it outside the plot area
                # legend placement depends on whether histograms are present
                if y_hist:
                    blank = figure(width_policy='min', height_policy='min', width=205, height=290, outline_line_color=None, background_fill_alpha=0.7, border_fill_alpha=0.7,)
                    legend = categorical_legend(detection_methods, counts.values(), '', Category10[10], fig=blank)
                else:
                    legend = categorical_legend(detection_methods, counts.values(), '', Category10[10], fig=plot)
                    plot.add_layout(legend, 'right')
            else:
                if clog:
                    cmap = log_cmap(cproperty, palette=Turbo256, low=np.nanmin(source_dict[cproperty]), high=np.nanmax(source_dict[cproperty]))
                    color_bar = ColorBar(color_mapper=cmap['transform'], ticker=LogTicker(), width=8, location=(0,0), background_fill_alpha=0)
                else:
                    cmap = linear_cmap(cproperty, palette=Turbo256, low=np.nanmin(source_dict[cproperty]), high=np.nanmax(source_dict[cproperty]))
                    color_bar = ColorBar(color_mapper=cmap['transform'], width=8, location=(0,0), background_fill_alpha=0)
                plot.add_layout(color_bar, 'right')
                # bokeh does not have a native colorbar label item so we will use a plot title as a substitute
                plot.title.text = cproperty.replace('_', ' ').replace('.', ' ').capitalize() + ' (%s)'%cunit
                plot.title_location = 'right'
                plot.title.align = 'center'
                plot.title.text_font_size = fontsize
                plot.title.text_font_style = 'italic'

            # plot graph
            plot.scatter(
                xproperty,
                yproperty,
                color=cmap,
                size='smap%i' %i,
                source=source,
            )

            # errorbars
            if xproperty + '_error_min' in Exoplanet.property_dict or xproperty.replace('host.', '') + '_error_min' in Star.property_dict:
                wx = Whisker(
                    source=source,
                    dimension='width',
                    visible=errorbar_x,
                    base=yproperty, upper=xproperty + '_error_max', lower=xproperty + '_error_min', line_color=cmap, upper_head=None, lower_head=None,
                )
                plot.add_layout(wx)
            if yproperty + '_error_min' in Exoplanet.property_dict or yproperty.replace('host.', '') + '_error_min' in Star.property_dict:
                wy = Whisker(
                    source=source,
                    dimension='height',
                    visible=errorbar_y,
                    base=xproperty, upper=yproperty + '_error_max', lower=yproperty + '_error_min', line_color=cmap, upper_head=None, lower_head=None,
                ) 
                plot.add_layout(wy)

            # histograms
            full_colour   = '#756bb1'
            sample_colour = '#74c476'
            if x_hist:
                # list setup
                full_list = source_dict[xproperty]
                short_list = [xval for xval, yval in zip(source_dict[xproperty], source_dict[yproperty]) if not np.isnan(yval)]
                min_  = np.nanmin(full_list)
                max_  = np.nanmax(full_list)

                # plot setup
                hxplot = figure(x_range=plot.x_range, x_axis_type=x_axis_type, height=300, sizing_mode='stretch_width', background_fill_alpha=0.7, border_fill_alpha=0.7,)    
                hxplot.xaxis.axis_label = xproperty.replace('_', ' ').replace('.', ' ').capitalize() + ' (%s)'%xunit
                hxplot.yaxis.axis_label = 'Frequency'
                hxplot.xaxis.axis_label_text_font_size = fontsize
                hxplot.yaxis.axis_label_text_font_size = fontsize

                # bins
                if xlog: bins = np.logspace(np.log10(min_), np.log10(max_), 100)
                else: bins = np.linspace(min_, max_, 100)

                # full range histogram
                hist, edges = np.histogram(full_list, bins=bins)
                hxplot.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:], color=full_colour, fill_alpha=0.5, legend_label='All')

                # only plotted histogram
                hist, edges = np.histogram(short_list, bins=bins)
                hxplot.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:], color=sample_colour, fill_alpha=0.5, legend_label='Shown')
                
                # legend setup
                hxplot.legend.click_policy = 'hide'
                hxplot.legend.background_fill_alpha = 0
                
            if y_hist:
                # list setup
                full_list = source_dict[yproperty]
                short_list = [yval for xval, yval in zip(source_dict[xproperty], source_dict[yproperty]) if not np.isnan(xval)]
                min_  = np.nanmin(full_list)
                max_  = np.nanmax(full_list)

                # plot setup
                hyplot = figure(y_range=plot.y_range, y_axis_type=y_axis_type, width=300, sizing_mode='stretch_height', background_fill_alpha=0.7, border_fill_alpha=0.7,)    
                hyplot.yaxis.axis_label = yproperty.replace('_', ' ').replace('.', ' ').capitalize() + ' (%s)'%yunit
                hyplot.xaxis.axis_label = 'Frequency'
                hyplot.xaxis.axis_label_text_font_size = fontsize
                hyplot.yaxis.axis_label_text_font_size = fontsize
                
                # bins
                if xlog: bins = np.logspace(np.log10(min_), np.log10(max_), 100)
                else: bins = np.linspace(min_, max_, 100)

                # full range histogram
                hist, edges = np.histogram(full_list, bins=bins)
                hyplot.quad(right=hist, left=0, top=edges[:-1], bottom=edges[1:], color=full_colour, fill_alpha=0.5, legend_label='All')

                # only plotted histogram
                hist, edges = np.histogram(short_list, bins=bins)
                hyplot.quad(right=hist, left=0, top=edges[:-1], bottom=edges[1:], color=sample_colour, fill_alpha=0.5, legend_label='Shown')
                
                # legend setup
                hyplot.legend.click_policy = 'hide'
                hyplot.legend.background_fill_alpha = 0

            # layouts
            if cproperty != 'detection_type': blank = None # don't need this plot if no legend is present
            if x_hist and y_hist: layout = gridplot([[hxplot, blank], [plot, hyplot]], sizing_mode='stretch_both')
            elif x_hist: layout = gridplot([[hxplot], [plot]], sizing_mode='stretch_both')
            elif y_hist: layout = gridplot([[plot, hyplot, blank]], sizing_mode='stretch_both')
            else: layout = gridplot([[plot]], sizing_mode='stretch_both')
            plots.append(layout)

    # layout
    layout = column(row(plots), sizing_mode='stretch_both')
    
    return file_html(layout, CDN)

def metric_rank_bar_graph(exoplanet_array, tooltip_dict=dict()):
    '''Funtion to generate html for an interactive bokeh bar graph, specifically for metric score vs rank'''
    
    # append to tooltip arrays
    for planet in exoplanet_array:
        for prop, vals in tooltip_dict.items():
            vals.append(rgetattr(planet, prop))
    
    # axis scale
    x_axis_type = 'linear'
    y_axis_type = 'log'

    # create tooltip table
    tooltips = [('(rank, metric)', '(@rank, @decision_metric)')]
    for prop in tooltip_dict.keys():
        tooltips.append( (prop, '@%s' %prop) )

    # create plot figure
    plot = figure(
        x_axis_type=x_axis_type,
        y_axis_type=y_axis_type,
        plot_height=400,
        sizing_mode='stretch_width',
        tools=[PanTool(), BoxZoomTool(), WheelZoomTool(), TapTool(), UndoTool(), RedoTool(), CrosshairTool(), HoverTool(mode='vline'), SaveTool(), ResetTool()],
        tooltips=tooltips,
        background_fill_alpha=0.7,
        border_fill_alpha=0.7,
    )

    # link points to exoplanet.eu
    url = 'http://www.exoplanet.eu/catalog/@name'
    taptool = plot.select(type=TapTool)
    taptool.callback = OpenURL(url=url.replace(' ', '_'))

    # label axes
    plot.xaxis.axis_label = 'Rank'
    plot.yaxis.axis_label = 'Decision Metric'
    fontsize = '12pt'
    plot.xaxis.axis_label_text_font_size = fontsize
    plot.yaxis.axis_label_text_font_size = fontsize

    # make data source
    source_dict = dict(
        rank=[i for i in range(1, len(exoplanet_array)+1)],
        decision_metric=[planet.decision_metric for planet in exoplanet_array],
    )
    source_dict.update(tooltip_dict) # add tooltip info to data source
    source = ColumnDataSource(data=source_dict)

    # make colormap
    detection_methods = ['Primary Transit', 'Radial Velocity', 'Imaging', 'Microlensing', 'Timing', 'Astrometry', 'TTV', 'Default', 'Primary Transit, TTV', 'Controversial',]
    cmap = factor_cmap('detection_type', palette=Category10[10], factors=detection_methods)
    
    # plot graph
    plot.vbar(
        x='rank',
        top='decision_metric',
        color=cmap,
        bottom=1, # so that log scale can work,
        width=0.9,
        legend_field='detection_type',
        source=source
    )

    # legend
    plot.legend.title = 'Detection method'
    plot.legend.location = 'top_right'
    plot.legend.background_fill_alpha = 0
    
    return file_html(gridplot([[plot]], sizing_mode='stretch_both'), CDN)

