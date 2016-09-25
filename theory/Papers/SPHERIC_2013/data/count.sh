#!/bin/bash

# Collect the average data in each file
for d in Cosmo Sedov SodShock
do
    rm -f ${d}.tot
    for f in `ls $d/perf*.dat | sort -V`
    do
        awk -v c=${f//[^0-9]/} '!/^#/{for(k=1;k<=NF;k++)sum[k]+=$k;}END{printf("%i ",c);for(k=1;k<=NF;k++)printf("%.3e ",sum[k]);printf("\n");}' $f >> ${d}.tot
    done
    rm -f ${d}_fixed.tot
    for f in `ls ${d}_fixed/perf*.dat | sort -V`
    do
        awk -v c=${f//[^0-9]/} '!/^#/{tot+=1;for(k=1;k<=NF;k++)sum[k]+=$k;}END{printf("%i ",c);for(k=1;k<=NF;k++)printf("%.3e ",sum[k]/tot);printf("\n");}' $f >> ${d}_fixed.tot
    done
done

