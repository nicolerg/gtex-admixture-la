#!/usr/bin/python

# convert a bed file of genomic positions to centimorgans using a genetic map 
# source: https://github.com/joepickrell/1000-genomes-genetic-maps/blob/master/scripts/interpolate_maps.py
# adjusted by Nicole Gay 

import sys, os, gzip

infile = open(sys.argv[1]) #.bed
mapfile = open(sys.argv[2]) #input map file, either the HapMap map or the 1000 Genomes OMNI map


posin = list()
rsin = list()
mappos = list()
mapgpos = list()
chrin = list()
mapchr = list()

line = infile.readline()
while line:
    line = line.strip().split()
    chr = line[0]
    pos = int(line[2])
    posin.append(pos)
    chrin.append(chr)
    rs = line[3]
    rsin.append(rs)
    line = infile.readline()

line = mapfile.readline()
line = mapfile.readline()
while line:
    line = line.strip().split()
    pos = int(line[3]) #for BEAGLE genetic map 
    chr = line[0]
    gpos = float(line[2]) #for BEAGLE genetic map 
    mappos.append(pos)
    mapgpos.append(gpos)
    mapchr.append(chr)
    line = mapfile.readline()

with gzip.open(sys.argv[3], "w") as outfile, gzip.open(sys.argv[4], "w") as loc:

    index1 = 0 #index for posin (bed pos in bp)
    index2 = 0 #index for mappos (map pos in bp)
    while index1 < len(posin):
        pos = posin[index1]
        chr = chrin[index1]
        rs = rsin[index1]
        if pos == mappos[index2]:
            #the 1000 Genomes site was genotyped as part of the map
            outfile.write('\t'.join([str(chr), str(pos), str(mapgpos[index2]), rs]) + '\n')
            loc.write(str(mapgpos[index2]) + '\n')
            # print >> outfile, chr, pos, mapgpos[index2], rs
            # print >> loc, mapgpos[index2]
            index1 = index1+1
        elif pos < mappos[index2]:
            #current position in interpolation before marker
            if index2 ==0:
                #before the first site in the map (genetic position = 0)
                outfile.write('\t'.join([str(chr), str(pos), str(mapgpos[index2]), rs]) + '\n')
                loc.write(str(mapgpos[index2]) + '\n')
                # print >> outfile, chr, pos, mapgpos[index2], rs 
                # print >> loc, mapgpos[index2]
                index1 = index1+1
            else:
                #interpolate
                prevg = mapgpos[index2-1]
                prevpos = mappos[index2]
                frac = (float(pos)-float(mappos[index2-1]))/ (float(mappos[index2]) - float(mappos[index2-1]))
                tmpg = prevg + frac* (mapgpos[index2]-prevg)
                outfile.write('\t'.join([str(chr), str(pos), str(tmpg), rs]) + '\n')
                loc.write(str(tmpg) + '\n')
                # print >> outfile, chr, pos, tmpg, rs 
                # print >> loc, tmpg
                index1 = index1+1
        elif pos > mappos[index2]:
            #current position in interpolation after marker
            if index2 == len(mappos)-1:
                #after the last site in the map (genetic position = maximum in map, note could try to extrapolate based on rate instead)
                outfile.write('\t'.join([str(chr), str(pos), str(mapgpos[index2]), rs]) + '\n')
                loc.write(str(mapgpos[index2]) + '\n')
                # print >> outfile, chr, pos, mapgpos[index2], rs
                # print >> loc, mapgpos[index2]
                index1 = index1+1
            else:
                #increment the marker
                index2 = index2+1
