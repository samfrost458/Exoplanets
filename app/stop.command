#!/usr/bin/env bash

kill -9 $(ps aux | grep redis | grep -v grep | awk '{print $2}' | tr '\n' ' ') > /dev/null 2>&1
kill -9 $(ps aux | grep celery | grep -v grep | awk '{print $2}' | tr '\n' ' ') > /dev/null 2>&1
kill -9 $(ps aux | grep flask | grep -v grep | awk '{print $2}' | tr '\n' ' ') > /dev/null 2>&1