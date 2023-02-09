#!/bin/bash

DIR=$1
YEAR=$2 # all, 2018, 2017, 2016APV, 2016nonAPV
DATA=$3 # 0, 1
BKG=$4 # 0, 1
SIG=$5 # 0, 1
SAM=$6
while ! [ -z "$7" ]; do
    FLAGS="$FLAGS $7"; shift;
done
COMMAND="./main.exe $DIR $YEAR $DATA $BKG $SIG $SAM$FLAGS"

echo ""
echo "Arguments:"
echo "----------"
echo "Directory = "$DIR
echo "Year = "$YEAR
echo "Run data = "$DATA
echo "Run bkg = "$BKG
echo "Run signal = "$SIG
echo "Sample = "$SAM
echo "Variation flags = "${FLAGS#", "}
echo ""

eval $COMMAND
