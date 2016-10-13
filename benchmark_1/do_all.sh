#!/bin/bash

# Do temperature calculation
scripts/setup_effgrain_temperature.py
source initial.sh

# Compute images
scripts/setup_images.py
source images.sh

# Compute SEDs
scripts/setup_seds.py
source seds.sh
