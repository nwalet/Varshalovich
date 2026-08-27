#!/usr/local/bin/bash
for i in fig11_*  
do     
    j=${i%.*}; 
    if ! grep -qF -- "$j}" ../Chap11.tex
    then
        echo "not found $j"
        mv $i obsolete/
        rm ../images/$j.pdf
    fi
done
