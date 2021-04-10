#!/usr/bin/env bash

# stop previous processes just in case
bash "$(dirname "$0")"/stop.command

# change directory
cd app

# if venv not installed then do so
if [ ! -d "./venv" ]
then
	(
	python3 -m venv ./venv   		
	source ./venv/bin/activate
	pip install -r requirements.txt
	)
else echo "venv found"
fi

# activate venv
source ./venv/bin/activate

# first subshell - redis
(
	./run-redis.sh
) &

# second subshell - celery
(
	cd ..
	celery worker -A app.routes.celery --loglevel=info
) &

# third subshell - flask
(
	cd ..
	flask run
) &