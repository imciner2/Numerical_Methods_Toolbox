#!/bin/bash
#
# This script will parse a .SIF file into the MATLAB CUTEst framework
# The files required for the mcutest to work will be packaged into a
# .tar.gz archive with the problem name
#
#
# Created by: Ian McInerney
# Created on: November 23, 2017
# Version: 1.0
# Last Modified: November 23, 2017
#
# Revision History:
#   1.0 - Initial Release

# Make sure that the cutest program is on the path
CUTEST_EXEC=$(command -v cutest2matlab)
if [ -z $CUTEST_EXEC ]; then
  echo cutest2matlab not found
  return 1
fi

# Get the problem name (which is the base of the filename)
PROBLEM_PATH=$1
PROBLEM_BASE=$2
ACTIVE_PATH=$3

# Get the problem's SIF filename
PROBLEM_SIF="$PROBLEM_BASE.SIF"

# Navigate to where the problem is located
CURR_DIR=$(pwd)

# Check to see if a tar archive of the problem exists in the directory
PROBLEM_TAR="$PROBLEM_BASE.tar.gz"
if [ ! -e $PROBLEM_PATH/$PROBLEM_TAR ]; then
  # The tar archive isn't here, create it

  # Navigate to the active problem directory
  cd $ACTIVE_PATH

  # Copy the SIF file
  cp $CURR_DIR/$PROBLEM_PATH/$PROBLEM_SIF ./

  # Create the appropriate files
  echo $CUTEST_EXEC
  $CUTEST_EXEC $PROBLEM_SIF
  
  # Create a tar.gz archive of the relevant files
  rm -f $PROBLEM_SIF
  tar -c * -f $PROBLEM_TAR
  cp $PROBLEM_TAR $CURR_DIR/$PROBLEM_PATH

fi
