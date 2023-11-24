#!/bin/bash


EXE=$MDTROOT/app/application/appWCTESingleEvent

WCSIM_FILE=/disk01/usr5/kmtsui/hyperk_repo/WCSim/build/install/wcsim.root
CONFIG_FILE=$MDTROOT/parameter/MDTParamenter_WCTE.txt
OUT_FILE=out_appWCTESingleEvent2.root
SEED=65457869

$EXE -i $WCSIM_FILE\
	 -p $CONFIG_FILE\
	 -o $OUT_FILE\
	 -s $SEED\
	 -n -1 # to run all events
