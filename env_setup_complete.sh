#!/bin/bash
# filepath: setup.sh

# Set GFAL Python binary path
export GFAL_PYTHONBIN=/usr/bin/python3

# Add CVMFS CMS common tools to PATH
export PATH="/cvmfs/cms.cern.ch/common/:${PATH}"

echo "CMS environment variables configured:"
echo "GFAL_PYTHONBIN=${GFAL_PYTHONBIN}"
echo "PATH now includes /cvmfs/cms.cern.ch/common/"