#!/bin/bash -l
NOW=$(date +"%d-%m-%Y-%s")
mkdir "temp_$NOW"

mv ./*.npy ./temp_$NOW/
python many_to_one_dicts.py temp_$NOW
tar czf "data$NOW.tar.gz" "temp_$NOW.npy"

rm -rf "temp_$NOW"
rm -rf "temp_$NOW.npy"