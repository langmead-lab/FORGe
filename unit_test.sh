#!/bin/sh

cd src && python squeakr_build.py ; cd ..
py.test -vr src/*.py
