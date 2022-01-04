#!/bin/sh

mkdir -p array_layouts
python create_array_layouts.py

mkdir -p sky_models
python create_sky_models.py

python run_accuracy_test.py
