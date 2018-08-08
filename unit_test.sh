#!/bin/sh

cd src && python cbuild.py ; cd ..
py.test -vr src/*.py
