#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 14:42:55 2020

@author: vijendersingh
"""


import scipy
import numpy
import pysam
import csv
import os
import sys as sys
#import subprocess

# BAMFILE
pathTObam=sys.argv[1]

# CSV FILE INFO

file_path='/path/to/CSV/dir'
# List CSV files seperated by comma
csvlist=['SOX2_motifs_ChIP_mm9']
# Specify the class of each CSV file seperated by comma, either TSS and TFB : TSS="Transcript start file", TFB="Transcription factor binding" sites
TFBorTSS=["TFB"]  

# OUTPUTDIR
output_directory="/path/to/output/dir/"

# SCRIPT PARAMETERS
frag_size=500

bp_up=2500

bp_dwn=2500

FrangmentSize_lowerlimit=0

FrangmentSize_Upperlimit=1000

# CLUSTER PLOT : Options:  TRUE or FALSE
clust_plot=False


############################################################
fileprefix=pathTObam.split("/")[-1][:-4]

bamfile=pysam.Samfile(pathTObam,"rb")

total=bamfile.mapped #+ bamfile.unmapped

# This poarameter "frag_size" represents avargae size of fragments that were sequenced
# Its better to set its value to be in range of 500-200bps for chip and nuc seq
# If you see a lack of counts towards left handside of the plot that means this value has to be increased.


cord_correct=int(frag_size)/2

#Number of bps analysed upstream and downstream of TSS or TF binding site
#these values can be set asymmetrically as well
# here we are looking into 5kb region flanking TSS, However in plot we will use only visualise 4KB flanking region.
# Ideally map 500bps extra on either side e.g if we need to look at at 3000bps upstream and down stream than
# set the value as 3500 (500 extra).

region=1+int(bp_up)+int(bp_dwn)+200

region_mid=1+int(bp_up)+int(bp_dwn)

read_count_start_coordinate=int(bp_up+cord_correct)

read_count_end_coordinate=int(bp_dwn+cord_correct)

# This variable is set to verify if the chromosome values are in "NN" or "chrNN" format
#It is used later while reading chromosome values from csv files
refl=int(len(bamfile.references[0]))

# If it is required to go through only one CSV file like CTCF.csv then commentout the above list
# and make csvlist variable with only one element
#   csvlist=['TSS_TTseqGenes_0vs2h']
csvt=0
for x in csvlist:
    cordCounter=500
    chromID=[]
    midpoint=[]
    start=[]
    end=[]
    strand=[]
    clust_plot_index=[]
    filename=str(csvlist[csvt]+'.csv')
    print (filename)
    initials=csvlist[csvt]
    os.chdir(file_path)
    print ('Bamfile Analysed : ' + pathTObam,'\t CSV file analysed : '+filename )
    #openfile=csv.reader(open(filename,'rb'))
    if TFBorTSS[csvt] == "TFB":
        openfile=csv.reader(open(filename))
        next(openfile,None)
        for row in openfile:
            if refl>2:
                chromID.append(row[0])
            else:
                chromID.append(row[0][3:])
            start.append(int(row[1]))
            end.append(int(row[2]))
            midpoint.append((int(row[1])+int(row[2]))/2)
    elif TFBorTSS[csvt] == "TSS":
        openfile=csv.reader(open(filename))
        next(openfile,None)
        for row in openfile:
            if refl>2 :
                chromID.append(row[0])
            else:
                chromID.append(row[0][3:])
            start.append(int(row[1]))
            end.append(int(row[2]))
            strand.append=[row[3]]
            midpoint=start
    else:
        print("The metatdata about CSV file is missing, \
              It is not specified if the files has TFB ot TSS  coocrdinates")

    signal_global=scipy.zeros((int(region_mid)))
    if clust_plot==True:
        clust_matrix=scipy.zeros((1001))

    icounter=0
    print (len(chromID))
    while icounter < len(chromID):
        signal_local=scipy.zeros((int(region_mid)))
        try:
            for reads in bamfile.fetch(chromID[icounter],midpoint[icounter]-read_count_start_coordinate, \
                                       midpoint[icounter]+read_count_end_coordinate):
                if FrangmentSize_lowerlimit < reads.tlen < FrangmentSize_Upperlimit:
                        coord=reads.pos+int((reads.tlen)/2)
                        postn=int(coord-(midpoint[icounter]-bp_up))
                        if postn>=0:
                                 try:
                                     signal_local[postn]+=1

                                 except IndexError:
                                         pass
        except (IndexError, ValueError):
            pass
        
        if TFBorTSS[csvt] == "TSS"  and strand[csvt]=="-1":
            signal_global+=signal_local[::-1]
        else:
            signal_global+=signal_local
            
        if clust_plot==True and numpy.sum(signal_local[1500:2501])>0:
            clust_matrix=numpy.vstack((clust_matrix,signal_local[1500:2501]))
            clust_plot_index.append(chromID[icounter]+":"+str(start[icounter]))                                       
        if icounter==cordCounter:
            print (icounter)
            cordCounter+=500
        icounter+=1


    Signal=((signal_global/total)*1000000)/(len(chromID))

    # If output directory is different then Add path for the output Directory and uncomment the lines below
    if FrangmentSize_lowerlimit==0 and FrangmentSize_Upperlimit==1000:
        output_npy=output_directory+fileprefix+'_'+initials+'.npy'
    else:
        output_npy=output_directory+fileprefix+'_'+initials+'_lower_'+ \
            str(FrangmentSize_lowerlimit)+'_upper_'+str(FrangmentSize_Upperlimit)+'.npy'
   # output_npy=fileprefix+'_'+initials+'.npy'    
    numpy.save(output_npy,Signal)
    del(Signal,signal_global)
    if clust_plot==True:
        import pandas as pd
        clust_matrix=numpy.delete(clust_matrix,0,0)
        df=pd.DataFrame(clust_matrix,index=clust_plot_index)
        df.to_csv("ClusterDF_"+fileprefix+'_'+initials+".csv",index=True)
    csvt+=1
  #  print (output_npy)
