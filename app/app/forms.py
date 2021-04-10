'''document containing flask forms'''

from flask_wtf import FlaskForm
from wtforms import SelectField, FloatField, IntegerField, StringField, DateField, BooleanField, SubmitField
from wtforms.validators import DataRequired, InputRequired, Optional, NumberRange
from datetime import date
import copy
from Exoplanet import Exoplanet, Star
from Telescope import Telescope

class DefaultTelescopeForm(FlaskForm):
    '''Handles pull down menu for selecting preexisting telescopes from the database'''

    choices = []
    # def __init__ necessary to pass tele_array variable from routes.py
    def __init__(self, tele_array, choices=choices):
        choices.clear() # necessary so updating the form won't keep doubling in length with no new data
        for telescope in tele_array:
            if telescope._C_T_calculable == True:
                choices.append((telescope.name, telescope.name.replace('_', ' ')))
        FlaskForm.__init__(self)
    
    teleName = SelectField( 'Select telescope', choices = choices, validators = [DataRequired()])
    # graph    = BooleanField('Graph plotter with returned planets')
    submit   = SubmitField( 'Submit')


class CustomTelescopeForm(FlaskForm):
    '''Handles custom telescope input page and optional observability inputs'''

    # magnitude zero points
    # for filter_ in Telescope.filter_list:
    #     if 'mag_' + filter_ in Star.property_dict:
    #         setattr(self, 'mag_' + filter_,
    #                 FloatField(   validators = [Optional()], render_kw = {"placeholder": "mzp %s" %filter_, "title": "Zero-point magnitude in the %s band" %filter_}))
    mzp_V 	      = FloatField(   validators = [Optional()], render_kw = {"placeholder": "mzp V", "title": "Zero-point magnitude in the V band"})
    mzp_I 	      = FloatField(   validators = [Optional()], render_kw = {"placeholder": "mzp I", "title": "Zero-point magnitude in the I band"})
    mzp_J 	      = FloatField(   validators = [Optional()], render_kw = {"placeholder": "mzp J", "title": "Zero-point magnitude in the J band"})
    mzp_H     	  = FloatField(   validators = [Optional()], render_kw = {"placeholder": "mzp H", "title": "Zero-point magnitude in the H band"})
    mzp_K 	      = FloatField(   validators = [Optional()], render_kw = {"placeholder": "mzp K", "title": "Zero-point magnitude in the K band"})
    
    # other telescope parameters
    spectral_bins = IntegerField( validators = [DataRequired(),  NumberRange(min = 1)], render_kw = {"placeholder": "No. spectral bins",                           "title": "No. spectral bins: 1 for photometric, >1 for spectroscopic (must be integer)"})
    t_over	      = FloatField(   validators = [InputRequired(), NumberRange(min = 0)], render_kw = {"placeholder": "Overhead time (s)",                           "title": "Time in seconds between exposures (can be 0 for our calculations)"})
    Half_Well 	  = IntegerField( validators = [DataRequired(),  NumberRange(min = 0)], render_kw = {"placeholder": "Half-well electron count",                    "title": "Number of photoelectrons to saturate each pixel to 50% (must be integer)"})
    texp_min 	  = FloatField(   validators = [Optional(),      NumberRange(min = 0)], render_kw = {"placeholder": "Minimum exposure time (s) (optional)",        "title": "Minimum exposure time in seconds (needed for some telescopes)"})
    Theta_see     = FloatField(   validators = [DataRequired(),  NumberRange(min = 0)], render_kw = {"placeholder": "Seeing size (arcsec)",                        "title": "(Average) seeing size at the telescope site in arcseconds, as a FWHM of a Gaussian\n(may of course be different on the night)"})
    Theta_DF      = FloatField(   validators = [Optional(),      NumberRange(min = 0)], render_kw = {"placeholder": "Defocused aperture size (arcsec) (optional)", "title": "Aperture size used for defocused operations in arcseconds"})
    Omega_pix     = FloatField(   validators = [DataRequired(),  NumberRange(min = 0)], render_kw = {"placeholder": "Pixel size (square arcsec)",                  "title": "Pixel area in square arcseconds"})

    # telescope location
    longitude     = FloatField(   validators = [Optional()], render_kw = {"placeholder": "Longitude (deg)", "title": "Longitude at telescope location in degrees (optional)"})
    latitude      = FloatField(   validators = [Optional()], render_kw = {"placeholder": "Latitude (deg)",  "title": "Latitude at telescope location in degrees (optional)"})
    elevation     = FloatField(   validators = [Optional()], render_kw = {"placeholder": "Elevation (m)",   "title": "Telescope elevation in metres (default 0)"})
    
    # further observability parameters
    min_altitude  = FloatField(   validators = [Optional()], render_kw = {"placeholder": "Minimum altitude (deg)",   "title": "Minimum altitude telescope can observe (deg) (default 0)"})
    max_altitude  = FloatField(   validators = [Optional()], render_kw = {"placeholder": "Maximum altitude (deg)",   "title": "Maximum altitude telescope can observe (deg) (default 90)"})
    start_date    = DateField(    validators = [Optional()], render_kw = {"placeholder": "Start date (yyyy-mm-dd)",  "title": "Starting date to check for observability (default today)"})
    end_date      = DateField(    validators = [Optional()], render_kw = {"placeholder": "End date (yyyy-mm-dd)",    "title": "Ending date to check for observability (default 30 days from start date)"})
    resolution    = FloatField(   validators = [Optional()], render_kw = {"placeholder": "Timing resolution (hr)",   "title": "Timing resolution in hours for checking observability (will significantly affect calculation time) (default 0.5)"})

    # radio button and sumbit button
    submit        = SubmitField()
    
    def validate(self):
        '''custom verification on top of standard checks already applied'''

        if not super(CustomTelescopeForm, self).validate():
            return False
        
        fine = True

        # check at least one mzp given
        if not (self.mzp_V.data or self.mzp_I.data or self.mzp_J.data or self.mzp_H.data or self.mzp_K.data):
            msg = 'At least one mzp must be set'
            self.mzp_V.errors.append(msg)
            self.mzp_I.errors.append(msg)
            self.mzp_J.errors.append(msg)
            self.mzp_H.errors.append(msg)
            self.mzp_K.errors.append(msg)
            fine = False
        
        # check spectral_bins and Half_Well are integers
        if not (isinstance(self.spectral_bins.data, int) and isinstance(self.Half_Well.data, int)):
            msg = 'Must be a positive integer'
            self.spectral_bins.errors.append(msg)
            self.Half_Well.errors.append(msg)
            fine = False

        # check if any location parameters given, then both latitude and longitude are given, and check limits
        if self.longitude.data or self.latitude.data or self.elevation.data or self.min_altitude.data or self.max_altitude.data or self.start_date.data or self.end_date.data or self.resolution.data:
            msg = 'Both latitude and longitude needed to check observability (otherwise leave all location fields blank)'
            if self.longitude.data is None:
                self.longitude.errors.append(msg)
                fine = False
            elif not (-180 <= self.longitude.data <= 180):
                self.longitude.errors.append('Longitude must be between -180 and 180 degrees')
                fine = False
            
            if self.latitude.data is None:
                self.latitude.errors.append(msg)
                fine = False
            elif not (-90 <= self.latitude.data <= 90):
                self.latitude.errors.append('Latitude must be between -90 and 90 degrees')
                fine = False

        # check end_date > start_date
        if self.end_date.data is not None:
            msg = 'End date must be after start date'
            if self.start_date.data is not None:
                if self.end_date.data <= self.start_date.data:
                    self.end_date.errors.append(msg)
                    fine = False
            else:
                if self.end_date.data <= date.today():
                    self.end_date.errors.append(msg)
                    fine = False
            
        # check max > min altitude
        if self.min_altitude.data is not None and self.max_altitude.data is not None:
            if self.min_altitude.data >= self.max_altitude.data:
                msg = 'Maximum altitude must be greater than minimum'
                self.max_altitude.errors.append(msg)
                fine = False
        
        # if all fine, return True, else return false
        return fine


class GraphForm(FlaskForm):
    '''Handles pull down menu and radioboxes for graph plotting'''
    
    # dropdown choices
    choices = [('None', 'None')]
    for prop, unit in Exoplanet.property_dict.items():
        if unit != '' and '_error' not in prop:
            choices.append((prop, prop.replace('_', ' ').capitalize()))
    for prop, unit in Star.property_dict.items():
        if unit != '' and '_error' not in prop:
            choices.append(('host.' + prop, 'Host star ' + prop.replace('_', ' ')))
    colour_choices = copy.deepcopy(choices)
    colour_choices.insert(1, ('detection_type', 'Detection type'))
    
    xproperty1 = SelectField('Select x property', default = 'mass',   choices = choices, validators = [DataRequired()])
    yproperty1 = SelectField('Select y property', default = 'radius', choices = choices, validators = [DataRequired()])
    cproperty1 = SelectField('Select colour property', default = 'detection_type',  choices = colour_choices, validators = [DataRequired()])
    sproperty1 = SelectField('Select size property',   default = 'None', choices = choices, validators = [DataRequired()])
    xlogscale1 = BooleanField('Log scale x axis', default = 'checked')
    ylogscale1 = BooleanField('Log scale y axis', default = 'checked')
    clogscale1 = BooleanField('Log scale colour', default = 'checked')
    slogscale1 = BooleanField('Log scale size',   default = 'checked')
    showxerr1  = BooleanField('Errorbars x axis')
    showyerr1  = BooleanField('Errorbars y axis')
    xhist1     = BooleanField('Histogram x')
    yhist1     = BooleanField('Histogram y')

    xproperty2 = SelectField('Select x property', default = 'None', choices = choices, validators = [DataRequired()])
    yproperty2 = SelectField('Select y property', default = 'None', choices = choices, validators = [DataRequired()])
    cproperty2 = SelectField('Select colour property', default = 'detection_type',  choices = colour_choices, validators = [DataRequired()])
    sproperty2 = SelectField('Select size property',   default = 'None', choices = choices, validators = [DataRequired()])
    xlogscale2 = BooleanField('Log scale x axis', default = 'checked')
    ylogscale2 = BooleanField('Log scale y axis', default = 'checked')
    clogscale2 = BooleanField('Log scale colour', default = 'checked')
    slogscale2 = BooleanField('Log scale size',   default = 'checked')
    showxerr2  = BooleanField('Errorbars x axis')
    showyerr2  = BooleanField('Errorbars y axis')
    xhist2     = BooleanField('Histogram x')
    yhist2     = BooleanField('Histogram y')

    submit     = SubmitField('Submit')


class GraphFromSelectionSubmit(FlaskForm):
    '''To redirect user to graph page but with their custom exoplanet array'''
    submit = SubmitField('Graph plotter for just these planets')

