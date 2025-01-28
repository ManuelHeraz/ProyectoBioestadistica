#!/bin/bash

# obtain most recent GO-Term description file
if [ ! -f "go-basic.obo" ]; then
    echo "Will download go-basic.obo now."
    wget http://purl.obolibrary.org/obo/go/go-basic.obo
fi

## remove header lines with
# first line with real content
lines=$(grep -n -m 1 "\[Term\]" go-basic.obo |cut -f1 -d: )
tail -n +${lines} go-basic.obo > tmp

# replace newlines with space
sed -i ':a;N;$!ba;s/\n/ /g' tmp

# reintroduce newlines via replacing [Term]
sed -i 's/\[Term\]/\n/g' tmp

# select fields
awk -F':' '{print $2":"$3"#\t\""$4"#\"\t"$5"#"}' tmp > tmp2

# discard unwanted stuff
sed 's/ namespace#//g' tmp2 | sed 's/ def#//g' | sed 's/ alt_id#//g' | sed 's/ name#//g' > tmp

# strip leading whitspaces
sed -i 's/\t /\t/g' tmp
sed -i 's/" /"/g' tmp
sed -i 's/^ //g' tmp

# drop first line
tail -n +2 tmp > tmp2

# finish here
mv tmp2 GOBioProcessAnnotationSimple.csv
rm tmp

gzip GOBioProcessAnnotationSimple.csv


