#!/bin/sh
# Just a simple script that will convert our readme file into a VERY basic webpage.
# The converted file needs to be uploaded to lys at: /var/www/html/lys_ucsd_edu/hexblender

INPUT_FILE="Readme.txt"
OUTPUT_FILE="index.html"

# Add a title to the page
echo "<title>Hexblender README</title>" > $OUTPUT_FILE 

# Add html breaks to the end of each line
sed 's|$|<br>|' $INPUT_FILE >> $OUTPUT_FILE
