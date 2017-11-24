#!/bin/bash
#
# This script will fetch the appropriate problem from the Maros and Meszaros
# QP benchmark suite, then parse the .SIF into the mcutest format and place it
# into the active problem folder
#
#
# Created by: Ian McInerney
# Created on: November 23, 2017
# Version: 1.0
# Last Modified: November 23, 2017
#
# Revision History:
#   1.0 - Initial Release

# Get the problem name from the 1st argument
PROBLEM_NAME=$1

# Change the filename to be all uppercase with the .SIF extension

CURR_DIR=$(pwd)
PROBLEM_PATH="../problems/marosmeszaros"
ACTIVE_PATH="../problems/activeCUTEst/"

PROBLEM_BASE=$(echo $PROBLEM_NAME | awk '{print toupper($0)}')
PROBLEM_SIF="$PROBLEM_BASE.SIF"

# Navigate to the problem directory
mkdir -p $PROBLEM_PATH
cd $PROBLEM_PATH

# Check to see if the .SIF file already exists in the directory
if [ ! -e $PROBLEM_SIF ]; then
  echo Problem not downloaded, downloading now
  # The .SIF file isn't here, get it from the web
  URL="ftp://ftp.numerical.rl.ac.uk/pub/cuter/marosmeszaros/$PROBLEM_SIF"
  wget $URL
else
  echo Problem already downloaded
fi

cd $CURR_DIR

# Parse the SIF file into the mcutest format
./parseSIF.sh $PROBLEM_PATH $PROBLEM_BASE $ACTIVE_PATH

# Paths and files
PROBLEM_TAR="$PROBLEM_BASE.tar.gz"

# Move the tar into the active problem directory and remove existing problems
echo Extracting problem into active problem directory
mkdir -p $ACTIVE_PATH
cd $ACTIVE_PATH
rm -f *

# Extract the new problem into the directory
cp $CURR_DIR/$PROBLEM_PATH/$PROBLEM_TAR ./
tar -xf $PROBLEM_TAR
rm -f $PROBLEM_TAR
