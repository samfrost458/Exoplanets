'''document containing routes'''

################################################################
# imports

# flask imports
from flask import render_template, url_for, flash, redirect, request, session, jsonify
from app import app
from app.forms import DefaultTelescopeForm, CustomTelescopeForm, GraphForm, GraphFromSelectionSubmit
from apscheduler.schedulers.background import BackgroundScheduler
from celery import Celery
# Celery configuration (leave in even though this is already in settings.py)
app.config['CELERY_BROKER_URL'] = 'redis://localhost:6379/0'
app.config['CELERY_RESULT_BACKEND'] = 'redis://localhost:6379/0'
# Initialize Celery
celery = Celery(app.name, broker=app.config['CELERY_BROKER_URL'])
celery.conf.update(app.config)

# standard python imports
import copy
import json
import os
import numpy as np

# astropy imports (for observability calculations)
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astroplan import Observer, FixedTarget, AltitudeConstraint, AtNightConstraint, months_observable, is_observable
from astropy.time import Time, TimeDelta
from datetime import datetime

# custom imports
from Exoplanet import Exoplanet, Star
from Telescope import Telescope
import functions


################################################################
# initial setup

# setup parameters to be fed into _update_data()
path_to_static = os.path.dirname(os.path.abspath(__file__)) + '/static/'
onlyConfirmed  = True # set to false to include unconfirmed planets 
planet_src     = 'http://www.exoplanet.eu/catalog/csv/'
tele_src       = path_to_static + 'Scope_Masters.csv'
to_remove      = ['1I/\\\'Oumuamua', '2I/Borisov'] # these items are comets, not planets

def _update_data():
    global exoplanet_array, tele_array, source_dict
    (exoplanet_array, tele_array, source_dict) = functions.read_planet_telescope_data_csv(planet_src, tele_src, path_to_static, onlyConfirmed, to_remove)
_update_data() # call _update_data on startup

#################################
# REMOVE WHEN DONE
import pickle
with open(path_to_static + 'exoplanet_array.pkl', 'rb') as f:
    exoplanet_array = pickle.load(f)
with open(path_to_static + 'tele_array.pkl', 'rb') as f:
    tele_array = pickle.load(f)
with open(path_to_static + 'source_dict.pkl', 'rb') as f:
    source_dict = pickle.load(f)
##################################

def _update_IERS_A():
    from astroplan import download_IERS_A
    download_IERS_A()

# schedule data to update itself in the background periodically
scheduler = BackgroundScheduler(daemon=True)  # grace time of an hour to make sure jobs don't get missed if they can't fire on the dot
scheduler.add_job(_update_data, 'cron', day='*', misfire_grace_time = 3600)     # update catalogue every day
scheduler.add_job(_update_IERS_A, 'cron', week='*', misfire_grace_time = 3600)  # update IERS Bulletin A every week
scheduler.start()

# list of properties to display
print_properties = [
    ['name',            'Planet name'],
    ['planet_status',   'Planet status'],
    ['mass',            'Mass in M_jup'],
    ['radius',          'Radius in R_jup'],
    ['orbital_period',  'Period in days'],
    ['semi_major_axis', 'Semi Major Axis in AU'],
    ['temp_calculated', 'Temperature in K'],
    ['ra',              'Right Ascension in h:m:s (J2000)'],
    ['dec',             'Declination in d:m:s (J2000)'],
    ['decision_metric', 'Calculated decision metric score (arbitrary units)'],
    ['filter',          'Filter which gives the best metric score'],
    ['t_exp',           'Exposure time in seconds']
]
if onlyConfirmed: print_properties.pop(1)


################################################################
# flask routes

@app.route('/')
def base(): return redirect(url_for('home'))


@app.route('/home', methods=['GET', 'POST'])
def home():
    '''Homepage. Renders full exoplanet library in alphabetical order and pulldown menu to select telescopes'''

    form = DefaultTelescopeForm(tele_array)
    if form.validate_on_submit():
        # find selected telescope from the tele_array
        for telescope in tele_array:
            if telescope.name == form.teleName.data:
                chosen_telescope = telescope

        # calculate metrics
        chosen_exoplanet_array = copy.deepcopy(exoplanet_array)
        for planet in chosen_exoplanet_array:
            planet.decision_metric_(chosen_telescope)

        # sort and truncate
        functions.sort_and_truncate_by_property(chosen_exoplanet_array, 'decision_metric')

        # choose a random planet
        if len(chosen_exoplanet_array) > 0:
            randomPlanet = functions.random_planet_selector(chosen_exoplanet_array)
        else: randomPlanet = False

        # graph of decision metric vs rank
        tooltip_dict = {
            'name'            : [],
            'mass'            : [],
            'radius'          : [],
            'orbital_period'  : [],
            'semi_major_axis' : [],
            'temp_calculated' : [],
            'detection_type'  : [],
            'decision_metric' : [],
            'filter'          : [],
            't_exp'           : [],
        }
        graph = functions.metric_rank_bar_graph(chosen_exoplanet_array, tooltip_dict)

        # return render of default table
        return render_template('telechoose.html', graph=graph, home=False, custom=False, form=form, print_properties=print_properties, exoplanet_array=chosen_exoplanet_array, tele_array=tele_array, randomPlanet=randomPlanet,
            _ra_dec_img_url=functions.ra_dec_img_url, _deg_to_ra_hms=functions.deg_to_ra_hms, _deg_to_dec_dms=functions.deg_to_dec_dms, _round_sig_figs=functions.round_sig_figs)

    # render the dropdown list for default telescope names above the rendered table
    base_properties = copy.deepcopy(print_properties)
    del base_properties[-3:]

    return render_template('telechoose.html', home=True, custom=False, form=form, print_properties=base_properties, exoplanet_array=exoplanet_array, tele_array=tele_array,
        _ra_dec_img_url=functions.ra_dec_img_url, _deg_to_ra_hms=functions.deg_to_ra_hms, _deg_to_dec_dms=functions.deg_to_dec_dms, _round_sig_figs=functions.round_sig_figs)


@app.route('/teleinput', methods=['GET', 'POST'])
def teleinput():
    '''Page to input custom telescope and observability parameters'''

    form = CustomTelescopeForm()
    location_properties = ['longitude', 'latitude', 'elevation', 'min_altitude', 'max_altitude', 'start_date', 'end_date', 'resolution']#, 'fast_return']
    if form.validate_on_submit():
        # create new telescope instance
        custom_telescope = Telescope()
        flash = [] # array to keep messages to flash the user to confirm their choices
        for prop in Telescope.property_dict.keys():
            if prop != 'name':
                val = getattr(form, prop).data
                setattr(custom_telescope, prop, val)
                # flash a message to the user confirming their choices
                if val is not None: flash.append('{} = {}'.format(prop, val))

        # error checking already done by form.py
        custom_telescope._C_T_calculable = True 

        # make a copy of the exoplanet array so as to not alter the master copy
        custom_exoplanet_array = copy.deepcopy(exoplanet_array)
        
        # calculate new metrics
        for planet in custom_exoplanet_array:
            planet.decision_metric_(custom_telescope)
        # sort and truncate custom_exoplanet_array
        functions.sort_and_truncate_by_property(custom_exoplanet_array, 'decision_metric')
        
        # copy print properties, since location may append it
        observable_properties = copy.deepcopy(print_properties)

        # if telescope location is given, see which planets are visible
        for property in location_properties:
            value = getattr(form, property).data
            if value is not None and property != 'start_date' and property != 'end_date':
                flash.append('{} = {}'.format(property, value))

        longitude    = form.longitude.data
        latitude     = form.latitude.data
        elevation    = form.elevation.data
        start_date   = form.start_date.data
        end_date     = form.end_date.data
        min_altitude = form.min_altitude.data
        max_altitude = form.max_altitude.data
        resolution   = form.resolution.data
        if longitude is not None and latitude is not None:

            # create json for exoplanet array, celery demands
            json_exoplanet_array = functions.planet_array_to_json_array(custom_exoplanet_array, 'host.')

            # start task
            task = observability.apply_async(args=(json_exoplanet_array, longitude, latitude, elevation, start_date, end_date, min_altitude, max_altitude, resolution))
            
            # render template, ajax calls to celery task will fill in the page
            return render_template('observability.html', home=False, custom=True, print_properties=print_properties, flash=flash, exoplanet_array=custom_exoplanet_array, task_id=task.id,
                _ra_dec_img_url=functions.ra_dec_img_url, _deg_to_ra_hms=functions.deg_to_ra_hms, _deg_to_dec_dms=functions.deg_to_dec_dms, _round_sig_figs=functions.round_sig_figs)

        # choose a random planet
        if len(custom_exoplanet_array) > 0:
            randomPlanet = functions.random_planet_selector(custom_exoplanet_array)
        else: randomPlanet = False

        # graph of decision metric vs rank
        tooltip_dict = {
            'name'            : [],
            'mass'            : [],
            'radius'          : [],
            'orbital_period'  : [],
            'semi_major_axis' : [],
            'temp_calculated' : [],
            'detection_type'  : [],
            'decision_metric' : [],
            'filter'          : [],
            't_exp'           : [],
        }
        graph = functions.metric_rank_bar_graph(custom_exoplanet_array, tooltip_dict)
        
        # render new table
        return render_template('table.html', graph=graph, home=False, custom=True, flash=flash, print_properties=print_properties, exoplanet_array=custom_exoplanet_array, randomPlanet=randomPlanet,
            _ra_dec_img_url=functions.ra_dec_img_url, _deg_to_ra_hms=functions.deg_to_ra_hms, _deg_to_dec_dms=functions.deg_to_dec_dms, _round_sig_figs=functions.round_sig_figs)
    
    # render the input form for telescope data
    return render_template('teleinput.html', property_list=Telescope.property_dict.keys(), location_properties=location_properties, form=form)


@app.route('/graph', methods=['GET', 'POST'])
def graph():
    '''Page to plot graphs on'''

    form = GraphForm()
    if form.validate_on_submit():
        # form properties
        xproperties  = [form.xproperty1.data, form.xproperty2.data]
        yproperties  = [form.yproperty1.data, form.yproperty2.data]
        cproperties  = [form.cproperty1.data, form.cproperty2.data]
        sproperties  = [form.sproperty1.data, form.sproperty2.data]
        xlogs        = [form.xlogscale1.data, form.xlogscale2.data]
        ylogs        = [form.ylogscale1.data, form.ylogscale2.data]
        clogs        = [form.clogscale1.data, form.clogscale2.data]
        slogs        = [form.slogscale1.data, form.slogscale2.data]
        xerrs        = [form.showxerr1.data,  form.showxerr2.data]
        yerrs        = [form.showyerr1.data,  form.showyerr2.data]
        x_hists      = [form.xhist1.data,     form.xhist2.data]
        y_hists      = [form.yhist1.data,     form.yhist2.data]
        tooltip_list = [
            # (x, y)
            # (c, s)
            'name',
            'mass',
            'radius',
            'orbital_period',
            'semi_major_axis',
            'temp_calculated',
            'detection_type',
        ]
        
        # plot graph
        graph = functions.linked_scatter_graph(source_dict, xproperties, yproperties, cproperties, sproperties, xlogs, ylogs, clogs, slogs, xerrs, yerrs, x_hists, y_hists, tooltip_list)
        return render_template('graph.html', form=form, graph=graph)

    return render_template('graph.html', form=form)


@app.route('/telelist')
def telelist():
    '''Page to show telescope library'''

    # get the list of properties from the Telescope property dictionary
    tele_properties = []
    for prop in Telescope.property_dict.keys():
        tele_properties.append(prop)
    
    # render the telescope table
    return render_template('telelist.html', tele_properties=tele_properties, tele_array=tele_array)


@app.route('/about')
def about():
    '''About page'''

    # make a legend with detection methods and counts for display purposesp
    detection_methods = ['Primary Transit', 'Radial Velocity', 'Imaging', 'Microlensing', 'Timing', 'Astrometry', 'TTV', 'Default', 'Primary Transit, TTV', 'Controversial']
    counts = {method: 0 for method in detection_methods}
    for method in source_dict['detection_type']: counts[method] += 1
    legend = functions.categorical_legend(detection_methods, counts.values(), 'Table colour coding', html=True, width=205, height=290)

    # render the about page
    return render_template('about.html', legend=legend)


################################################################
# celery tasks

@celery.task(bind=True)
def observability(self, json_exoplanet_array, longitude, latitude, elevation, start_date, end_date, min_altitude, max_altitude, resolution):
    '''background task to (slowly) find planet observabilities'''
    
    # create Observer instance
    if elevation is None: elevation = 0
    custom_location = Observer(longitude=longitude*u.deg, latitude=latitude*u.deg, elevation=elevation*u.m)

    # time range
    if start_date is not None: start_date = Time(datetime.combine(datetime.fromisoformat(start_date), datetime.min.time()), out_subfmt='date')
    else: start_date = Time(datetime.combine(datetime.now(), datetime.min.time()), out_subfmt='date')
    if end_date is not None: end_date = Time(datetime.combine(datetime.fromisoformat(end_date), datetime.min.time()), out_subfmt='date')
    else: end_date = start_date + 30*u.day
    time_range = Time([start_date, end_date])

    # altitude and night time constraints
    if min_altitude is None: min_altitude = 0
    if max_altitude is None: max_altitude = 90
    constraints = [AltitudeConstraint(min_altitude*u.deg, max_altitude*u.deg), AtNightConstraint.twilight_astronomical()]
    
    # calculate if observable within lookahead time from now
    if resolution is None: resolution = 0.5
    
    # since we had to use a json string as argument, recreate the custom_exoplanet_array from the json
    custom_exoplanet_array = []
    json_unload = json.loads(json_exoplanet_array) # unpack json to make a dictionary
    for dict_ in json_unload:
        planet = Exoplanet()
        # copy properties back over from dictionary
        for prop, val in dict_.items(): functions.rsetattr(planet, prop, val)
        custom_exoplanet_array.append(planet)

    # trying to split array into chunks so we can say the task has been completed to the nearest percent
    num_chunks = round(len(custom_exoplanet_array)/800 * (resolution/0.5)**-1 * ((end_date-start_date).jd/30)**2 * 10)
    if num_chunks == 0: num_chunks = 1
    exoplanet_array_chunks = np.array_split(custom_exoplanet_array, num_chunks) # split exoplanet array into num_chunk parts

    planets_done = 0 # planets calculated for so far
    planets_lost = 0 # planets removed so far

    recombined_chunks = []
    for chunk in exoplanet_array_chunks:
        # make an array of FixedTarget instances
        target_table = []
        for planet in chunk:
            target_table.append((planet.name, planet.host.ra, planet.host.dec))
        targets = [FixedTarget(coord=SkyCoord(ra=ra*u.deg, dec=dec*u.deg), name=name)
                   for name, ra, dec in target_table]

        # calculate
        ever_observable = is_observable(constraints, custom_location, targets, time_range=time_range, time_grid_resolution=resolution*u.hour)

        # note unobservables
        to_pop = []
        for i in range(len(chunk)):
            chunk[i].observable = bool(ever_observable[i])
            if chunk[i].observable == False:
                to_pop.append(i)
        # remove unobservables
        chunk = np.delete(chunk, to_pop)

        # record numbers
        planets_done += len(targets)
        planets_lost += len(targets) - len(chunk)

        # write json
        recombined_chunks += chunk.tolist()
        observable_exoplanet_array = [planet for planet in recombined_chunks if planet.observable]
        json_exoplanet_array = functions.planet_array_to_json_array(observable_exoplanet_array, 'host_')

        # graph of decision metric vs rank
        tooltip_dict = {
            'name'            : [],
            'mass'            : [],
            'radius'          : [],
            'orbital_period'  : [],
            'semi_major_axis' : [],
            'temp_calculated' : [],
            'detection_type'  : [],
            'decision_metric' : [],
            'filter'          : [],
            't_exp'           : [],
        }
        graph = functions.metric_rank_bar_graph(observable_exoplanet_array, tooltip_dict)

        # update status
        self.update_state(state='PROGRESS',
                          meta={'current'         : planets_done,
                                'removed'         : planets_lost,
                                'total'           : len(custom_exoplanet_array),
                                'exoplanet_array' : json_exoplanet_array,
                                'graph'           : graph,
                                'finished'        : 'n'})
    
    return {'current': planets_done, 'removed': planets_lost, 'total': len(custom_exoplanet_array), 'exoplanet_array': json_exoplanet_array, 'graph': graph, 'finished': 'y'}

@app.route('/status/<task_id>')
def taskstatus(task_id):
    task = observability.AsyncResult(task_id)
    if task.state == 'PENDING':
        response = {
            'state'          : task.state,
            'current'        : 0,
            'removed'        : 0,
            'total'          : 1,
            'exoplanet_array': 'Pending...',
            'graph'          : 0,
            'finished'       : 'n',
        }
    elif task.state != 'FAILURE':
        response = {
            'state'          : task.state,
            'current'        : task.info.get('current', 0),
            'removed'        : task.info.get('removed', 0),
            'total'          : task.info.get('total', 1),
            'exoplanet_array': task.info.get('exoplanet_array'),
            'graph'          : task.info.get('graph'),
            'finished'       : task.info.get('finished'),
        }
    else:
        # something went wrong in the background job
        response = {
            'state'          : task.state,
            'current'        : 1,
            'removed'        : 0,
            'total'          : 1,
            'exoplanet_array': str(task.info),  # this is the exception raised
            'graph'          : 0,
            'finished'       : 'n',
        }
    return jsonify(response)


################################################################

if __name__ == '__main__':
    app.jinja_env.auto_reload = True
    app.config['TEMPLATES_AUTO_RELOAD'] = True
    app.run(debug=True)

