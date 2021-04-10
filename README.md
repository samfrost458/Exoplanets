# Exoplanet project
## A Master's project for The University of Manchester in 2019/2020, automating the selection of already-discovered exoplanets for follow-up observations

### Setup
#### Easy mode
In a terminal window:
```
cd /path/to/app  # OUTER folder, containing run.command
bash run.command # must write here; double-clicking the file does not seem to work
```
This (should) do everything from creating a venv and downloading all the packages to starting each subshell process. They will keep going in the background, so if they need to be stopped, double-click on stop.command

If it doesn't work then...
#### Manual mode
First check all the packages in the venv are version specific:
```
cd path/to/app/app # INNER folder called app, containing routes.py
python3 -m venv venv   		
source venv/bin/activate
pip install -r requirements.txt
```
Following this, open three separate terminal tabs:
###### Terminal 1
```
cd path/to/app/app
source venv/bin/activate
./run-redis.sh
```
###### Terminal 2
```
cd path/to/app/app
source venv/bin/activate
cd ..
celery worker -A app.routes.celery --loglevel=info
```
###### Terminal 3
```
cd path/to/app/app
source venv/bin/activate
cd ..
FLASK_ENV=development
flask run
```

#### Upon success
Once the message
```
 * Running on http://127.0.0.1:5000/ (Press CTRL+C to quit)
```
appears (after a few seconds while the up-to-date database downloads), open a browser and head to `127.0.0.1:5000/home`

(If it takes a while to run it is likely downloading many images from the internet (check in app/app/static for new files); it will not take this long on subsequent startups)

It should be safe to close the terminal window and have everything still run in the background, so double-click on stop.command if you actually need the processes to stop

### About
This browser-based app is the result of an MPhys project undertaken at The University of Manchester in 2019/2020 by [Sam Frost](mailto:sam.frost458@gmail.com) and [James Berezowski](mailto:jamesberezowski93@gmail.com), supervised by [Dr. Eamonn Kerins](http://www.jb.man.ac.uk/~ekerins/current.html).
Hundreds to thousands of exoplanets are discovered every year, and many deserve further observations for study. However, previous observations of this sort were chosen by hand. It is now clear with the increasing numbers of planets being found that a way of automating this process must be developed. The team at [SPEARNET](http://www.spearnet-team.org) have developed a [decision metric](https://arxiv.org/abs/1802.05645) based on the signal-to-noise ratio for a given planet observation. During the project we validated assumptions made by this metric and suggested improvements, including use of a new mass-radius relation, and have produced ranked tables with metrics calculated from the [exoplanet.eu](http://www.exoplanet.eu) archive. Images are provided by [ESO Online Digitized Sky Survey](https://archive.eso.org/dss/dss).

The app returns a list of planets, ranked by decision metric for a given telescope (choose from a list or input details for your own). The metric attempts to show which planets should return a high signal-to-noise ratio when observed, effectively showing which planets would be 'easiest' to observe. If a user wants to plan an observation during a specific window, options exist to eliminate the list to only those planets visible during the window, given the location of the telescope.

As a convenience, we include graph plotting tools, similar to those found on [exoplanet.eu](http://www.exoplanet.eu/diagrams) or [Filtergraph](https://filtergraph.com/nea), with a few extra features they do not offer.

The project owes a debt of gratitude to [Miguel Grinberg](https://github.com/miguelgrinberg), whose tutorials on [flask](https://blog.miguelgrinberg.com/post/the-flask-mega-tutorial-part-i-hello-world-legacy) and [celery](https://blog.miguelgrinberg.com/post/using-celery-with-flask) were invaluable. The [run-redis.sh file](app/app/run-redis.sh) is taken directly from his [flask-celery-example](https://github.com/miguelgrinberg/flask-celery-example/blob/master/run-redis.sh).
