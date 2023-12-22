#!/bin/bash

# Install required Python packages using pip
#pip install skyfield --target "$(pwd)/relative/path/to/installation"
pip install skyfield --target '%cd%\relative\path\to\installation'