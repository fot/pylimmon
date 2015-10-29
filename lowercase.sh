#!/bin/bash
#
for filename in $(ls *.csv); do
  echo "Converting $filename to lower case"
  dd if=$filename of=$filename.lower conv=lcase
done

for i in *; do
  echo "Converting $i file name to lower case"
  mv "$i" "$(echo $i|tr A-Z a-z)";
done

# on a mac:
# rename -s csv.lower csv *