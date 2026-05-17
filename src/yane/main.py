# Description
##YAne: A high-Yielding primer Assessment and database construction tool for eDNA metabarcoding

########################################################################################################################
# region Import modules
import sys
import argparse
import os
import numpy as np
import pandas as pd
import seaborn as sns #may need to update
import matplotlib.pyplot as plt
from Bio import SeqIO
import warnings
warnings.filterwarnings("ignore", "Support for passing numbers through unit converters")
# endregion
########################################################################################################################
# region Function#1: Download NCBI data

def downloadmock(args):

    # region Assign parsed arguments from command flags to python variables
    gateshead = args.i_mock
    taxid = args.p_taxa
    my = args.p_dl_period
    jobs = args.p_jobs
    threads = args.p_threads
    # endregion

    # region Save "taxid" and "my" variables to input argument log file #1 >> use in processprimer function
    farg1 = open("inp_args1.log.txt","w")
    farg1.write("taxa"+"\t"+"dl_period"+"\n")
    farg1.write(taxid+"\t"+my)
    farg1.close()
    # endregion

    # region Create new directories
    newdir = "mkdir "+"NCBI_data"
    os.system(newdir)
    # endregion

    # region Generate QIIME2 commands: Downloading NCBI data
    print("#===================================== Download Mock Community ===========================================#")

    ## Download NCBI data with "qiime rescript get-ncbi-data"
    print("#### Download sequences and taxonomy from NCBI ####")
    f = open(gateshead,"r") #open gene table
    header = f.readline() #read header
    for line in f:
        #Variables
        gene = line.split("\t")[0]
        query = line.split("\t")[1].strip()
        #Files
        out_seq = "NCBI_data/"+"txid"+taxid+"_"+gene+"_"+my+"_seq.qza"
        out_tax = "NCBI_data/"+"txid"+taxid+"_"+gene+"_"+my+"_tax.qza"
        #Command
        cmd = ("qiime rescript get-ncbi-data \\\n"+"--p-query "+"'"+query+"'"+" \\\n"+"--p-rank-propagation \\\n"+"--p-n-jobs "+str(jobs)+" \\\n"+"--o-sequences "+out_seq+" \\\n"+"--o-taxonomy "+out_tax)
        #Running
        print("\n# "+gene+" >>> Running . . . ")
        print(cmd)
        os.system(cmd)
    f.close()

    ## Dereplicate datasets with "qiime rescript dereplicate"
    print("\n#### Dereplicate NCBI sequences and taxonomy ####")
    f = open(gateshead,"r") 
    header = f.readline()
    for line in f:
        #Variables
        gene = line.split("\t")[0]
        query = line.split("\t")[1].strip()
        #Files
        in_seq = "NCBI_data/"+"txid"+taxid+"_"+gene+"_"+my+"_seq.qza"
        in_tax = "NCBI_data/"+"txid"+taxid+"_"+gene+"_"+my+"_tax.qza"
        out_seq = "NCBI_data/"+"txid"+taxid+"_"+gene+"_"+my+"_drp_seq.qza"
        out_tax = "NCBI_data/"+"txid"+taxid+"_"+gene+"_"+my+"_drp_tax.qza"
        #Command
        cmd = "qiime rescript dereplicate \\\n"+"--i-sequences "+in_seq+" \\\n"+"--i-taxa "+in_tax+" \\\n"+"--p-mode "+"uniq"+" \\\n"+"--p-threads "+str(threads)+ " \\\n"+"--o-dereplicated-sequences "+out_seq+" \\\n"+"--o-dereplicated-taxa "+out_tax
        #Running
        print("\n# "+gene+" >>> Running . . . ")
        print(cmd)
        os.system(cmd)
    f.close()

    # endregion

# endregion
########################################################################################################################
# region Function#2: Process NCBI data

def processprimer(args):

    # region Assign parsed arguments from command flags to python variables
    gateshead = args.i_mock
    lowood = args.i_primer
    thornfield = args.o_dir
    jobs = args.p_jobs
    threads = args.p_threads
    identity = args.p_identity
    hmplm = args.p_homopolymer
    degen = args.p_degenerate
    prilab = args.p_label
    observe_len = args.p_observe_length
    level = args.p_length_stat_level
    # endregion

    # region Retrieve "taxid" and "my" variables from the previous command (downloadmock)
    farg1 = open("inp_args1.log.txt","r")
    farg1.readline()
    for line in farg1:
        taxid = line.split("\t")[0]
        my = line.split("\t")[1].strip()
    farg1.close()
    # endregion

#=======================================================================================================================
    # region Extract hypervariable regions from primer pairs

    print("#=================================== Extract hypervariable regions =======================================#")

    # Create new directories
    # region Create new directories
    newdir1 = "mkdir "+thornfield
    newdir2 = "mkdir "+thornfield+"/extract"
    os.system(newdir1)
    os.system(newdir2)
    # endregion

    # Vary parameter #1: one default identity (identity = 0.8)
    if type(identity) == float:
        # region Generate QIIME2 commands: Extract hypervariable regions

        ## Extract hypervariable regions with "qiime feature-classifier extract-reads"
        print("\n#### Extract hypervariable regions ####")
        f = open(lowood,"r") #open primer table
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Files
            in_seq = "NCBI_data/"+"txid"+taxid+"_"+gene+"_"+my+"_drp_seq.qza"
            out_seq = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_seq.qza"
            #Command
            cmd = "qiime feature-classifier extract-reads \\\n"+"--i-sequences "+in_seq+" \\\n"+"--p-f-primer "+forward+" \\\n"+"--p-r-primer "+reverse+" \\\n"+"--p-identity "+str(identity)+" \\\n"+"--p-n-jobs "+str(jobs)+" \\\n"+"--o-reads "+out_seq
            #Running
            print("\n# "+primer+" >>> Running . . .")
            print(cmd)
            os.system(cmd)
        f.close()

        ## Dereplicate datasets with "qiime rescript dereplicate"
        print("\n#### Dereplicate hypervariable regions ####")
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Files
            in_seq = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_seq.qza"
            in_tax = "NCBI_data/"+"txid"+taxid+"_"+gene+"_"+my+"_drp_tax.qza"
            out_seq = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_seq.qza"
            out_tax = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_tax.qza"
            #Command
            cmd = "qiime rescript dereplicate \\\n"+"--i-sequences "+in_seq+" \\\n"+"--i-taxa "+in_tax+" \\\n"+"--p-mode "+"uniq"+" \\\n"+"--p-threads "+str(threads)+ " \\\n"+"--o-dereplicated-sequences "+out_seq+" \\\n"+"--o-dereplicated-taxa "+out_tax
            #Running
            print("\n# "+primer+" >>> Running . . . ")
            print(cmd)
            os.system(cmd)
        f.close()

        # endregion

    # Vary parameter#2: one customized identity (e.g., identity = 0.85)
    elif type(identity) == list and len(identity) == 1:
        # region Generate QIIME2 commands: Extract hypervariable regions

        ## Extract hypervariable regions with "qiime feature-classifier extract-reads"
        print("\n#### Extract hypervariable regions ####")
        f = open(lowood,"r") #open primer table
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Files
            in_seq = "NCBI_data/"+"txid"+taxid+"_"+gene+"_"+my+"_drp_seq.qza"
            out_seq = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_seq.qza"
            #Command
            cmd = "qiime feature-classifier extract-reads \\\n"+"--i-sequences "+in_seq+" \\\n"+"--p-f-primer "+forward+" \\\n"+"--p-r-primer "+reverse+" \\\n"+"--p-identity "+str(identity[0])+" \\\n"+"--p-n-jobs "+str(jobs)+" \\\n"+"--o-reads "+out_seq
            #Running
            print("\n# "+primer+" >>> Running . . .")
            print(cmd)
            os.system(cmd)
        f.close()

        ## Dereplicate datasets with "qiime rescript dereplicate"
        print("\n#### Dereplicate hypervariable regions ####")
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Files
            in_seq = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_seq.qza"
            in_tax = "NCBI_data/"+"txid"+taxid+"_"+gene+"_"+my+"_drp_tax.qza"
            out_seq = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_seq.qza"
            out_tax = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_tax.qza"
            #Command
            cmd = "qiime rescript dereplicate \\\n"+"--i-sequences "+in_seq+" \\\n"+"--i-taxa "+in_tax+" \\\n"+"--p-mode "+"uniq"+" \\\n"+"--p-threads "+str(threads)+ " \\\n"+"--o-dereplicated-sequences "+out_seq+" \\\n"+"--o-dereplicated-taxa "+out_tax
            #Running
            print("\n# "+primer+" >>> Running . . . ")
            print(cmd)
            os.system(cmd)
        f.close()

        # endregion

    # Vary parameter#3: multiple customized identity specific for each primer dataset
    elif type(identity) == list and len(identity) == len(prilab):
        # region Generate QIIME2 commands: Extract hypervariable regions

        ## Extract hypervariable regions with "qiime feature-classifier extract-reads"
        print("\n#### Extract hypervariable regions ####")
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Files
            in_seq = "NCBI_data/"+"txid"+taxid+"_"+gene+"_"+my+"_drp_seq.qza"
            out_seq = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_seq.qza"
            #Customized variables
            for i in range(len(prilab)):
                if primer == prilab[i]:
                    cmd = "qiime feature-classifier extract-reads \\\n"+"--i-sequences "+in_seq+" \\\n"+"--p-f-primer "+forward+" \\\n"+"--p-r-primer "+reverse+" \\\n"+"--p-identity "+str(identity[i])+" \\\n"+"--p-n-jobs "+str(jobs)+" \\\n"+"--o-reads "+out_seq
                    #Running
                    print("\n# "+primer,">>> Running . . .")
                    print(cmd)
                    os.system(cmd)

        ## Dereplicate datasets with "qiime rescript dereplicate"
        print("\n#### Dereplicate hypervariable regions ####")
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Files
            in_seq = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_seq.qza"
            in_tax = "NCBI_data/"+"txid"+taxid+"_"+gene+"_"+my+"_drp_tax.qza"
            out_seq = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_seq.qza"
            out_tax = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_tax.qza"
            #Command
            cmd = "qiime rescript dereplicate \\\n"+"--i-sequences "+in_seq+" \\\n"+"--i-taxa "+in_tax+" \\\n"+"--p-mode "+"uniq"+" \\\n"+"--p-threads "+str(threads)+ " \\\n"+"--o-dereplicated-sequences "+out_seq+" \\\n"+"--o-dereplicated-taxa "+out_tax
            #Running
            print("\n# "+primer,">>> Running . . . ")
            print(cmd)
            os.system(cmd)
        f.close()

        # endregion

    # In case of an ERROR (Did not provide matching number of primer labels with the customized parameters)
    else:
        print("Error: number of customized parameters and primer names are not matched. Please check if all of the primer names are provided.")

    # endregion
#=======================================================================================================================
    # region Filter homopolymers and degenerate base lengths

    print("#====================================== Filter hmplm & degen =============================================#")

    # Create new directories
    # region Create new directories
    newdir1 = "mkdir "+thornfield+"/quality"
    newdir2 = "mkdir "+thornfield+"/quality/hmplm_degen"
    os.system(newdir1)
    os.system(newdir2)
    # endregion
    
    # Vary parameter #1: default degen & default hmplm (degen = 5, hmplm = 8)
    if type(degen) == int and type(hmplm) == int:
        # region Generate QIIME2 commands: Filter hmplm & degen length

        ## Filter homopolymers and degenerate bases with "qiime rescript cull-seqs"
        print("\n#### Filter homopolymers & degenerate bases ####")
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Files
            in_seq = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_seq.qza"
            out_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_seq.qza"
            #Command
            cmd = "qiime rescript cull-seqs \\\n"+"--i-sequences "+in_seq+" \\\n"+"--p-n-jobs "+str(jobs)+" \\\n"+"--p-num-degenerates "+str(degen)+" \\\n"+"--p-homopolymer-length "+str(hmplm)+ " \\\n"+"--o-clean-sequences "+out_seq
            #Running
            print("\n# "+primer+" >>> Running . . . ")
            print(cmd)
            os.system(cmd)
        f.close()

        ## Dereplicate datasets with "qiime rescript dereplicate"
        print("\n#### Dereplicate filtered datasets ####")
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Files
            in_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_seq.qza"
            in_tax = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_tax.qza"
            out_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_seq.qza"
            out_tax = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_tax.qza"
            #Command
            cmd = "qiime rescript dereplicate \\\n"+"--i-sequences "+in_seq+" \\\n"+"--i-taxa "+in_tax+" \\\n"+"--p-mode "+"uniq"+" \\\n"+"--p-threads "+str(threads)+ " \\\n"+"--o-dereplicated-sequences "+out_seq+" \\\n"+"--o-dereplicated-taxa "+out_tax
            #Running
            print("\n# "+primer+" >>> Running . . . ")
            print(cmd)
            os.system(cmd)
        f.close()

        # endregion

    # Vary parameter #2: one-custom degen & default hmplm (e.g., degen = 1, hmplm = 8)
    elif type(degen) == list and type(hmplm) == int and len(degen) == 1:
        # region Generate QIIME2 commands: Filter hmplm & degen length

        ## Filter homopolymers & degenerate bases with "qiime rescript cull-seqs"
        print("\n#### Filter homopolymers & degenerate bases ####")
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Files
            in_seq = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_seq.qza"
            out_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_seq.qza"
            #Command
            cmd = "qiime rescript cull-seqs \\\n"+"--i-sequences "+in_seq+" \\\n"+"--p-n-jobs "+str(jobs)+" \\\n"+"--p-num-degenerates "+str(degen[0])+" \\\n"+"--p-homopolymer-length "+str(hmplm)+ " \\\n"+"--o-clean-sequences "+out_seq
            #Running
            print("\n# "+primer+" >>> Running . . . ")
            print(cmd)
            os.system(cmd)
        f.close()

        ## Dereplicate datasets with "qiime rescript dereplicate"
        print("\n#### Dereplicate filtered datasets ####")
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Files
            in_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_seq.qza"
            in_tax = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_tax.qza"
            out_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_seq.qza"
            out_tax = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_tax.qza"
            #Command
            cmd = "qiime rescript dereplicate \\\n"+"--i-sequences "+in_seq+" \\\n"+"--i-taxa "+in_tax+" \\\n"+"--p-mode "+"uniq"+" \\\n"+"--p-threads "+str(threads)+ " \\\n"+"--o-dereplicated-sequences "+out_seq+" \\\n"+"--o-dereplicated-taxa "+out_tax
            #Running
            print("\n# "+primer+" >>> Running . . . ")
            print(cmd)
            os.system(cmd)
        f.close()

        # endregion

    # Vary parameter #3: many-custom degen & default hmplm (e.g., degen = 1,2,3,5,8,9..., hmplm = 8)
    elif type(degen) == list and type(hmplm) == int and len(degen) == len(prilab):
        # region Generate QIIME2 commands: Filter hmplm & degen length

        ## Filter homopolymers & degenerate bases with "qiime rescript cull-seqs"
        print("\n#### Filter homopolymers & degenerate bases ####")
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Files
            in_seq = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_seq.qza"
            out_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_seq.qza"
            #Customized variables
            for i in range(len(prilab)):
                if primer == prilab[i]:
                    cmd = "qiime rescript cull-seqs \\\n"+"--i-sequences "+in_seq+" \\\n"+"--p-n-jobs "+str(jobs)+" \\\n"+"--p-num-degenerates "+str(degen[i])+" \\\n"+"--p-homopolymer-length "+str(hmplm)+ " \\\n"+"--o-clean-sequences "+out_seq
                    #Running
                    print("\n# "+primer,">>> Running . . . ")
                    print(cmd)
                    os.system(cmd)
        f.close()

        ## Dereplicate datasets with "qiime rescript dereplicate"
        print("\n#### Dereplicate filtered datasets ####")
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Files
            in_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_seq.qza"
            in_tax = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_tax.qza"
            out_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_seq.qza"
            out_tax = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_tax.qza"
            #Command
            cmd = "qiime rescript dereplicate \\\n"+"--i-sequences "+in_seq+" \\\n"+"--i-taxa "+in_tax+" \\\n"+"--p-mode "+"uniq"+" \\\n"+"--p-threads "+str(threads)+ " \\\n"+"--o-dereplicated-sequences "+out_seq+" \\\n"+"--o-dereplicated-taxa "+out_tax
            #Running
            print("\n# "+primer+" >>> Running . . . ")
            print(cmd)
            os.system(cmd)
        f.close()

        # endregion

    # Vary parameter #4: default degen & one-custom hmplm (e.g., degen = 5, hmplm = 10)
    elif type(degen) == int and type(hmplm) == list and len(hmplm) == 1:
        # region Generate QIIME2 commands: Filter hmplm & degen length

        ## Filter homopolymers & degenerate bases with "qiime rescript cull-seqs"
        print("\n#### Filter homopolymers & degenerate bases ####")
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Files
            in_seq = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_seq.qza"
            out_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_seq.qza"
            #Command
            cmd = "qiime rescript cull-seqs \\\n"+"--i-sequences "+in_seq+" \\\n"+"--p-n-jobs "+str(jobs)+" \\\n"+"--p-num-degenerates "+str(degen)+" \\\n"+"--p-homopolymer-length "+str(hmplm[0])+ " \\\n"+"--o-clean-sequences "+out_seq
            #Running
            print("\n# "+primer+" >>> Running . . . ")
            print(cmd)
            os.system(cmd)
        f.close()

        ## Dereplicate datasets with "qiime rescript dereplicate"
        print("\n#### Dereplicate filtered datasets ####")
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Files
            in_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_seq.qza"
            in_tax = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_tax.qza"
            out_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_seq.qza"
            out_tax = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_tax.qza"
            #Command
            cmd = "qiime rescript dereplicate \\\n"+"--i-sequences "+in_seq+" \\\n"+"--i-taxa "+in_tax+" \\\n"+"--p-mode "+"uniq"+" \\\n"+"--p-threads "+str(threads)+ " \\\n"+"--o-dereplicated-sequences "+out_seq+" \\\n"+"--o-dereplicated-taxa "+out_tax
            #Running
            print("\n# "+primer+" >>> Running . . . ")
            print(cmd)
            os.system(cmd)
        f.close()

        # endregion

    # Vary parameter #5: one-custom degen & one-custom hmplm (e.g., degen = 1, hmplm = 10)
    elif type(degen) == list and type(hmplm) == list and len(degen) == 1 and len(hmplm) == 1:
        # region Generate QIIME2 commands: Filter hmplm & degen length

        ## Filter homopolymers & degenerate bases with "qiime rescript cull-seqs"
        print("\n#### Filter homopolymers & degenerate bases ####")
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Files
            in_seq = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_seq.qza"
            out_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_seq.qza"
            #Command
            cmd = "qiime rescript cull-seqs \\\n"+"--i-sequences "+in_seq+" \\\n"+"--p-n-jobs "+str(jobs)+" \\\n"+"--p-num-degenerates "+str(degen[0])+" \\\n"+"--p-homopolymer-length "+str(hmplm[0])+ " \\\n"+"--o-clean-sequences "+out_seq
            #Running
            print("\n# "+primer+" >>> Running . . . ")
            print(cmd)
            os.system(cmd)
        f.close()

        ## Dereplicate datasets with "qiime rescript dereplicate"
        print("\n#### Dereplicate filtered datasets ####")
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Files
            in_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_seq.qza"
            in_tax = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_tax.qza"
            out_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_seq.qza"
            out_tax = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_tax.qza"
            #Command
            cmd = "qiime rescript dereplicate \\\n"+"--i-sequences "+in_seq+" \\\n"+"--i-taxa "+in_tax+" \\\n"+"--p-mode "+"uniq"+" \\\n"+"--p-threads "+str(threads)+ " \\\n"+"--o-dereplicated-sequences "+out_seq+" \\\n"+"--o-dereplicated-taxa "+out_tax
            #Running
            print("\n# "+primer+" >>> Running . . . ")
            print(cmd)
            os.system(cmd)
        f.close()

        # endregion

    # Vary parameter #6: many-custom degen & one-custom hmplm (e.g., degen = 1,2,3,4,5,..., hmplm = 10)
    elif type(degen) == list and type(hmplm) == list and len(degen) == len(prilab) and len(hmplm) == 1:
        # region Generate QIIME2 commands: Filter hmplm & degen length

        ## Filter homopolymers & degenerate bases with "qiime rescript cull-seqs"
        print("\n#### Filter homopolymers & degenerate bases ####")
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Files
            in_seq = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_seq.qza"
            out_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_seq.qza"
            #Customized variables
            for i in range(len(prilab)):
                if primer == prilab[i]:
                    cmd = "qiime rescript cull-seqs \\\n"+"--i-sequences "+in_seq+" \\\n"+"--p-n-jobs "+str(jobs)+" \\\n"+"--p-num-degenerates "+str(degen[i])+" \\\n"+"--p-homopolymer-length "+str(hmplm[0])+ " \\\n"+"--o-clean-sequences "+out_seq
                    #Running
                    print("\n# "+primer+" >>> Running . . . ")
                    print(cmd)
                    os.system(cmd)
        f.close()

        ## Dereplicate datasets with "qiime rescript dereplicate"
        print("\n#### Dereplicate filtered datasets ####")
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Files
            in_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_seq.qza"
            in_tax = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_tax.qza"
            out_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_seq.qza"
            out_tax = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_tax.qza"
            #Command
            cmd = "qiime rescript dereplicate \\\n"+"--i-sequences "+in_seq+" \\\n"+"--i-taxa "+in_tax+" \\\n"+"--p-mode "+"uniq"+" \\\n"+"--p-threads "+str(threads)+ " \\\n"+"--o-dereplicated-sequences "+out_seq+" \\\n"+"--o-dereplicated-taxa "+out_tax
            #Running
            print("\n# "+primer+" >>> Running . . . ")
            print(cmd)
            os.system(cmd)
        f.close()

        # endregion
    
    # Vary parameter #7: default degen & many-custom hmplm (e.g., degen = 8, hmplm = 1,2,3,4,5,...)
    elif type(degen) == int and type(hmplm) == list and len(hmplm) == len (prilab):
        # region Generate QIIME2 commands: Filter hmplm & degen length

        ## Filter homopolymers & degenerate bases with "qiime rescript cull-seqs"
        print("\n#### Filter homopolymers & degenerate bases ####")
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Files
            in_seq = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_seq.qza"
            out_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_seq.qza"
            #Customized variables
            for i in range(len(prilab)):
                if primer == prilab[i]:
                    cmd = "qiime rescript cull-seqs \\\n"+"--i-sequences "+in_seq+" \\\n"+"--p-n-jobs "+str(jobs)+" \\\n"+"--p-num-degenerates "+str(degen)+" \\\n"+"--p-homopolymer-length "+str(hmplm[i])+ " \\\n"+"--o-clean-sequences "+out_seq
                    #Running
                    print("\n# "+primer,">>> Running . . . ")
                    print(cmd)
                    os.system(cmd)
        f.close()

        ## Dereplicate datasets with "qiime rescript dereplicate"
        print("\n#### Dereplicate filtered datasets ####")
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Files
            in_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_seq.qza"
            in_tax = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_tax.qza"
            out_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_seq.qza"
            out_tax = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_tax.qza"
            #Command
            cmd = "qiime rescript dereplicate \\\n"+"--i-sequences "+in_seq+" \\\n"+"--i-taxa "+in_tax+" \\\n"+"--p-mode "+"uniq"+" \\\n"+"--p-threads "+str(threads)+ " \\\n"+"--o-dereplicated-sequences "+out_seq+" \\\n"+"--o-dereplicated-taxa "+out_tax
            #Running
            print("\n# "+primer+" >>> Running . . . ")
            print(cmd)
            os.system(cmd)
        f.close()

        # endregion

    # Vary parameter #8: one-custom degen & many-custom homo (e.g., degen = 1, hmplm = 1,2,3,4,5,...)
    elif type(degen) == list and type(hmplm) == list and len(degen) == 1 and len(hmplm) == len (prilab):
        # region Generate QIIME2 commands: Filter hmplm & degen length

        ## Filter homopolymers & degenerate bases with "qiime rescript cull-seqs"
        print("\n#### Filter homopolymers & degenerate bases ####")
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Files
            in_seq = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_seq.qza"
            out_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_seq.qza"
            #Customized variables
            for i in range(len(prilab)):
                if primer == prilab[i]:
                    cmd = "qiime rescript cull-seqs \\\n"+"--i-sequences "+in_seq+" \\\n"+"--p-n-jobs "+str(jobs)+" \\\n"+"--p-num-degenerates "+str(degen[0])+" \\\n"+"--p-homopolymer-length "+str(hmplm[i])+ " \\\n"+"--o-clean-sequences "+out_seq
                    #Running
                    print("\n# "+primer,">>> Running . . . ")
                    print(cmd)
                    os.system(cmd)
        f.close()

        ## Dereplicate datasets with "qiime rescript dereplicate"
        print("\n#### Dereplicate filtered datasets ####")
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Files
            in_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_seq.qza"
            in_tax = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_tax.qza"
            out_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_seq.qza"
            out_tax = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_tax.qza"
            #Command
            cmd = "qiime rescript dereplicate \\\n"+"--i-sequences "+in_seq+" \\\n"+"--i-taxa "+in_tax+" \\\n"+"--p-mode "+"uniq"+" \\\n"+"--p-threads "+str(threads)+ " \\\n"+"--o-dereplicated-sequences "+out_seq+" \\\n"+"--o-dereplicated-taxa "+out_tax
            #Running
            print("\n# "+primer+" >>> Running . . . ")
            print(cmd)
            os.system(cmd)
        f.close()

        # endregion

    ## Vary parameter #9: many-custom degen & many-custom homo (e.g., degen = 12,3,4,5,..., hmplm = 1,2,3,4,5,...)
    elif type(degen) == list and type(hmplm) == list and len(degen) == len(prilab) and len(hmplm) == len (prilab):
        # region Generate QIIME2 commands: Filter hmplm & degen length

        ## Filter homopolymers & degenerate bases with "qiime rescript cull-seqs"
        print("\n#### Filter homopolymers & degenerate bases ####")
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Files
            in_seq = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_seq.qza"
            out_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_seq.qza"
            #Customized variables
            for i in range(len(prilab)):
                if primer == prilab[i]:
                    cmd = "qiime rescript cull-seqs \\\n"+"--i-sequences "+in_seq+" \\\n"+"--p-n-jobs "+str(jobs)+" \\\n"+"--p-num-degenerates "+str(degen[i])+" \\\n"+"--p-homopolymer-length "+str(hmplm[i])+ " \\\n"+"--o-clean-sequences "+out_seq
                    #Running
                    print("\n# "+primer,">>> Running . . . ")
                    print(cmd)
                    os.system(cmd)
        f.close()

        ## Dereplicate datasets with "qiime rescript dereplicate"
        print("\n#### Dereplicate filtered datasets ####")
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Files
            in_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_seq.qza"
            in_tax = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_tax.qza"
            out_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_seq.qza"
            out_tax = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_tax.qza"
            #Command
            cmd = "qiime rescript dereplicate \\\n"+"--i-sequences "+in_seq+" \\\n"+"--i-taxa "+in_tax+" \\\n"+"--p-mode "+"uniq"+" \\\n"+"--p-threads "+str(threads)+ " \\\n"+"--o-dereplicated-sequences "+out_seq+" \\\n"+"--o-dereplicated-taxa "+out_tax
            #Running
            print("\n# "+primer+" >>> Running . . . ")
            print(cmd)
            os.system(cmd)
        f.close()

        # endregion

    ## In case of an ERROR (Did not provide matching number of primer labels with the customized parameters)
    else:
        print("Error: number of customized parameters and primer names are not matched. Please check if all of the primer names are provided.")

    # endregion
#=======================================================================================================================
    # region Filter sequence length with default criteria or observe length before setting customized criteria

    # Create new directories
    # region Create new directories
    newdir1 = "mkdir "+thornfield+"/export"
    newdir2 = "mkdir "+thornfield+"/quality/seq_len"
    newdir3 = "mkdir "+thornfield+"/clean_datasets"
    os.system(newdir1)
    os.system(newdir2)
    os.system(newdir3)
    # endregion

    # Filter length with default parameters (observe_length == "no")
    if observe_len == "no":
        print("#=========================== Filter sequence length with default criteria ============================#")

        # region Export dereplicated hmplm and degen filtered datasets with "qiime tools export"
        print("\n#### Export dereplicated homopolymer and degenerate base filtered datasets ####")
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            # Files
            in_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_seq.qza"
            in_tax = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_tax.qza"
            out_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_seq"
            out_tax = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_tax"
            # Commands
            cmd_seq = "qiime tools export \\\n"+"--input-path "+in_seq+" \\\n"+"--output-path "+out_seq
            cmd_tax = "qiime tools export \\\n"+"--input-path "+in_tax+" \\\n"+"--output-path "+out_tax
            # Running
            print("\n# "+primer+" >>> Running . . .")
            print(cmd_seq)
            os.system(cmd_seq)
            print(cmd_tax)
            os.system(cmd_tax)
        f.close()
        # endregion

        # region Create sequence length dictionary
        print("\n#### Get sequence length statistics of all records among primer datasets ####")
        # Create sequence length dictionary
        seqlen_dict = {}
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            # Files
            in_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_seq/dna-sequences.fasta"
            in_tax = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_tax/taxonomy.tsv"
            # Create sequence length dictionary
            len_list = []
            records = list(SeqIO.parse(in_seq, "fasta"))
            for record in records:
                len_list.append(len(record.seq))
            seqlen_dict[primer] = len_list
        f.close()
        print("#>>> Created sequence length dictionary for all datasets")
        # endregion

        # region Calculate sequence length statistics
        out_res = thornfield+"/quality/seq_len/seqlen_stat.txt"
        fw = open(out_res,"w")
        fw.write("Primer"+"\t"+"Total records"+"\t"+"Arithmetic mean"+"\t"+"sd"+"\t"+"Min"+"\t"+"Max"+"\t"+"P1"+"\t"+"P2.5"+"\t"+"P5"+"\t"+"P25"+"\t"+"P50"+"\t"+"P75"+"\t"+"P95"+"\t"+"P97.5"+"\t"+"P99"+"\t"+"Lower fence"+"\t"+"Upper fence"+"\n")
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Statistics
            length_array = np.array(seqlen_dict[primer])
            num_rec = len(seqlen_dict[primer])
            arthmean = np.mean(length_array)
            sd = np.std(length_array)
            min_len = np.min(length_array)
            max_len = np.max(length_array)
            p010 = np.percentile(length_array, 1)
            p025 = np.percentile(length_array, 2.5)
            p050 = np.percentile(length_array, 5)
            p250 = np.percentile(length_array, 25)
            p500 = np.percentile(length_array, 50)
            p750 = np.percentile(length_array, 75)
            p950 = np.percentile(length_array, 95)
            p975 = np.percentile(length_array, 97.5)
            p990 = np.percentile(length_array, 99)
            low_fen = p250-1.5*(p750-p250)
            up_fen = p750+1.5*(p750-p250)
            fw.write(primer+"\t"+str(num_rec)+"\t"+str(arthmean)+"\t"+str(sd)+"\t"+str(min_len)+"\t"+str(max_len)+"\t"+str(p010)+"\t"+str(p025)+"\t"+str(p050)+"\t"+str(p250)+"\t"+str(p500)+"\t"+str(p750)+"\t"+str(p950)+"\t"+str(p975)+"\t"+str(p990)+"\t"+str(low_fen)+"\t"+str(up_fen)+"\n")
        fw.close()    
        f.close()
        # endregion

        # region Parse calculated percentiles (min = p1, max = p99) to be used as default filtering criteria
        primer_dict = {}
        fs = open(out_res,"r")
        header = fs.readline()
        for line in fs:
            cri_list = []
            pm = line.split("\t")[0]
            p1 = line.split("\t")[6]
            p99 = line.split("\t")[14]
            cri_list.append(int(round(float(p1))))
            cri_list.append(int(round(float(p99))))
            primer_dict[pm] = cri_list
        fs.close()
        print("#>>> Retrieved sequence length filtering criteria (min=p1, max=p99) for all datasets")
        print("primer"+"\t"+"min"+"\t"+"max")
        for primer in primer_dict.keys():
            print(primer+"\t"+str(primer_dict[primer][0])+"\t"+str(primer_dict[primer][1]))
        # endregion

        # region Transpose statistics table results and save as new file (primers as columns)
        df = pd.read_csv(out_res, sep='\t')
        df_transposed = df.transpose()
        df_transposed.to_csv(out_res, sep='\t', header=False)    
        print("#>>> Saved sequence length statistics table as:",out_res)
        # endregion

        # region Generate QIIME2 commands: Filter sequence length with default criteria
        ## Filter sequence length with "qiime rescript filter-seqs-length"
        print("\n#### Filter sequence length by default parameter (min=p1, max=p99) ####")
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            # Files
            in_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_seq.qza"
            out_seq = thornfield+"/quality/seq_len/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_len_seq.qza"
            out_seq_DISC = thornfield+"/quality/seq_len/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_len_DISC_seq.qza"
            # Command
            cmd = "qiime rescript filter-seqs-length \\\n"+"--i-sequences "+in_seq+" \\\n"+"--p-global-min "+str(primer_dict[primer][0])+" \\\n"+"--p-global-max "+str(primer_dict[primer][1])+" \\\n"+"--o-filtered-seqs "+out_seq+" \\\n"+"--o-discarded-seqs "+out_seq_DISC
            #Running
            print("\n# "+primer+" >>> Running . . .")
            print(cmd)
            os.system(cmd)
        f.close()

        ## Dereplicate datasets with "qiime rescript dereplicate"
        print("\n#### Dereplicate filtered datasets ####")
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Files
            in_seq = thornfield+"/quality/seq_len/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_len_seq.qza" 
            in_tax = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_tax.qza"
            out_seq = thornfield+"/clean_datasets/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_seq.qza"
            out_tax = thornfield+"/clean_datasets/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_tax.qza"
            #Command
            cmd = "qiime rescript dereplicate \\\n"+"--i-sequences "+in_seq+" \\\n"+"--i-taxa "+in_tax+" \\\n"+"--p-mode "+"uniq"+" \\\n"+"--p-threads "+str(threads)+ " \\\n"+"--o-dereplicated-sequences "+out_seq+" \\\n"+"--o-dereplicated-taxa "+out_tax
            #Running
            print("\n# "+primer+" >>> Running . . . ")
            print(cmd)
            os.system(cmd)
        f.close()

        # endregion

        # region Summarize number of records through filtering steps (default length filtering)
        print("#============================= Count number of records (default length) ==============================#")

        # region Export artifact files from QIIME2 with "qiime tools export"
        ## Mock community datasets
        print("\n#### Export artifact files ####")
        fg = open(gateshead,"r")
        header = fg.readline()
        for line in fg:
            #Variables
            gene = line.split("\t")[0]
            query = line.split("\t")[1].strip()
            # Files
            ## Inputs
            in_NCBI_seq = "NCBI_data/"+"txid"+taxid+"_"+gene+"_"+my+"_seq.qza"
            in_NCBI_tax = "NCBI_data/"+"txid"+taxid+"_"+gene+"_"+my+"_tax.qza"
            in_drpNCBI_seq = "NCBI_data/"+"txid"+taxid+"_"+gene+"_"+my+"_drp_seq.qza"
            in_drpNCBI_tax = "NCBI_data/"+"txid"+taxid+"_"+gene+"_"+my+"_drp_tax.qza"
            ## Outputs
            out_NCBI_seq = thornfield+"/export/"+"txid"+taxid+"_"+gene+"_"+my+"_seq"
            out_NCBI_tax = thornfield+"/export/"+"txid"+taxid+"_"+gene+"_"+my+"_tax"
            out_drpNCBI_seq = thornfield+"/export/"+"txid"+taxid+"_"+gene+"_"+my+"_drp_seq"
            out_drpNCBI_tax = thornfield+"/export/"+"txid"+taxid+"_"+gene+"_"+my+"_drp_tax"
            # Commands 
            cmd_NCBI_seq = "qiime tools export \\\n"+"--input-path "+in_NCBI_seq+" \\\n"+"--output-path "+out_NCBI_seq
            cmd_NCBI_tax = "qiime tools export \\\n"+"--input-path "+in_NCBI_tax+" \\\n"+"--output-path "+out_NCBI_tax
            cmd_drpNCBI_seq = "qiime tools export \\\n"+"--input-path "+in_drpNCBI_seq+" \\\n"+"--output-path "+out_drpNCBI_seq
            cmd_drpNCBI_tax = "qiime tools export \\\n"+"--input-path "+in_drpNCBI_tax+" \\\n"+"--output-path "+out_drpNCBI_tax
            #Running
            print("\n# "+gene+" >>> Running . . .")
            print("\n## NCBI data:")
            print(cmd_NCBI_seq)
            os.system(cmd_NCBI_seq)
            print(cmd_NCBI_tax)
            os.system(cmd_NCBI_tax)
            print("\n## Dereplicated NCBI data:")
            print(cmd_drpNCBI_seq)
            os.system(cmd_drpNCBI_seq)
            print(cmd_drpNCBI_tax)
            os.system(cmd_drpNCBI_tax)
        fg.close()

        ## Primer-extracted datasets
        fp = open(lowood,"r")
        header = fp.readline()
        for line in fp:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            # Files
            ## Inputs
            in_hyper_seq = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_seq.qza"
            in_drphyper_seq = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_seq.qza"
            in_drphyper_tax = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_tax.qza"
            in_cull_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_seq.qza"
            in_drpcull_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_seq.qza"
            in_drpcull_tax = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_tax.qza"
            in_fillen_seq = thornfield+"/quality/seq_len/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_len_seq.qza"
            in_drpfillen_seq = thornfield+"/clean_datasets/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_seq.qza"
            in_drpfillen_tax = thornfield+"/clean_datasets/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_tax.qza"
            ## Outputs
            out_hyper_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_seq"
            out_drphyper_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_seq"
            out_drphyper_tax = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_tax"
            out_cull_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_seq"
            out_drpcull_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_seq"
            out_drpcull_tax = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_tax"
            out_fillen_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_len_seq"
            out_drpfillen_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_seq"
            out_drpfillen_tax = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_tax"
            # Commands
            cmd_hyper_seq = "qiime tools export \\\n"+"--input-path "+in_hyper_seq+" \\\n"+"--output-path "+out_hyper_seq
            cmd_drphyper_seq = "qiime tools export \\\n"+"--input-path "+in_drphyper_seq+" \\\n"+"--output-path "+out_drphyper_seq
            cmd_drphyper_tax = "qiime tools export \\\n"+"--input-path "+in_drphyper_tax+" \\\n"+"--output-path "+out_drphyper_tax
            cmd_cull_seq = "qiime tools export \\\n"+"--input-path "+in_cull_seq+" \\\n"+"--output-path "+out_cull_seq
            cmd_drpcull_seq = "qiime tools export \\\n"+"--input-path "+in_drpcull_seq+" \\\n"+"--output-path "+out_drpcull_seq
            cmd_drpcull_tax = "qiime tools export \\\n"+"--input-path "+in_drpcull_tax+" \\\n"+"--output-path "+out_drpcull_tax
            cmd_fillen_seq = "qiime tools export \\\n"+"--input-path "+in_fillen_seq+" \\\n"+"--output-path "+out_fillen_seq
            cmd_drpfillen_seq = "qiime tools export \\\n"+"--input-path "+in_drpfillen_seq+" \\\n"+"--output-path "+out_drpfillen_seq
            cmd_drpfillen_tax = "qiime tools export \\\n"+"--input-path "+in_drpfillen_tax+" \\\n"+"--output-path "+out_drpfillen_tax
            #Running
            print("\n# "+primer,">>> Running . . .")
            print("\n## Extracted hypervariable regions:")
            print(cmd_hyper_seq)
            os.system(cmd_hyper_seq)
            print("\n## Dereplicated hypervariable regions:")
            print(cmd_drphyper_seq)
            os.system(cmd_drphyper_seq)
            print(cmd_drphyper_tax)
            os.system(cmd_drphyper_tax)
            print("\n## Filtered homopolymers (hmplm) & degenerate bases (degen):")
            print(cmd_cull_seq)
            os.system(cmd_cull_seq)
            print("\n## Dereplicated filtered hmplm & degen:")
            print(cmd_drpcull_seq)
            os.system(cmd_drpcull_seq)
            print(cmd_drpcull_tax)
            os.system(cmd_drpcull_tax)
            print("\n## Filtered sequence lengths:")
            print(cmd_fillen_seq)
            os.system(cmd_fillen_seq)
            print("\n## Dereplicated filtered lengths (cleaned datasets):")
            print(cmd_drpfillen_seq)
            os.system(cmd_drpfillen_seq)
            print(cmd_drpfillen_tax)
            os.system(cmd_drpfillen_tax)
        fp.close()
        # endregion

        # region Create new directories
        newdir1 = "mkdir "+thornfield+"/results"
        newdir2 = "mkdir "+thornfield+"/results/number_of_records"
        os.system(newdir1)
        os.system(newdir2)
        # endregion

        # region Count number of records through processing steps
        ## Count records in mock community datasets with "grep"
        print("\n#### Count number of records: Gene datasets ####")
        fg = open(gateshead,"r")
        header = fg.readline()
        out_res1 = thornfield+"/results/number_of_records/grep_mock.log.txt"
        for line in fg:
            #Variables
            gene = line.split("\t")[0]
            query = line.split("\t")[1].strip()
            # Files
            in_NCBI = thornfield+"/export/"+"txid"+taxid+"_"+gene+"_"+my+"_seq/dna-sequences.fasta"
            in_drpNCBI = thornfield+"/export/"+"txid"+taxid+"_"+gene+"_"+my+"_drp_seq/dna-sequences.fasta"
            # Commands
            write_gene = "echo "+"'>>"+gene+"'"+" >> "+out_res1
            write_ns = "echo '#NCBI_data' >> "+out_res1
            count_ns = "grep '>' "+in_NCBI+" | wc -l"+" >> "+out_res1
            write_dns = "echo '#Dereplicated_NCBI_data' >> "+out_res1
            count_dns = "grep '>' "+in_drpNCBI+" | wc -l"+" >> "+out_res1
            #Running
            print("\n# "+gene+" >>> Running . . .")
            os.system(write_gene)
            print("\n## NCBI data:")
            os.system(write_ns)
            print(count_ns)
            os.system(count_ns)
            print("\n## Dereplicated NCBI data:")
            os.system(write_dns)
            print(count_dns)
            os.system(count_dns)
            print("#---------------------------------------------------")
        fg.close()

        ## Count records in primer-extracted datasets with "grep"
        print("\n#### Count number of records: Primer datasets ####")
        fp = open(lowood,"r")
        header = fp.readline()
        out_res2 = thornfield+"/results/number_of_records/grep_primer.log.txt"
        for line in fp:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            # Files   
            in_hyper = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_seq/dna-sequences.fasta"
            in_drphyper = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_seq/dna-sequences.fasta"
            in_cull = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_seq/dna-sequences.fasta"
            in_drpcull = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_seq/dna-sequences.fasta"
            in_fillen = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_len_seq/dna-sequences.fasta"
            in_drpfillen = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_seq/dna-sequences.fasta"
            # Commands
            write_primer = "echo "+"'>>"+primer+"'"+" >> "+out_res2
            write_hs = "echo '#Extracted_hyperregions:' >> "+out_res2
            count_hs = "grep '>' "+in_hyper+" | wc -l"+" >> "+out_res2
            write_dhs = "echo '#Dereplicated_hyperregions:' >> "+out_res2
            count_dhs = "grep '>' "+in_drphyper+" | wc -l"+" >> "+out_res2
            write_cs = "echo '#Filtered_hmplm_degen:' >> "+out_res2
            count_cs = "grep '>' "+in_cull+" | wc -l"+" >> "+out_res2
            write_dcs = "echo '#Dereplicated_filtered_hmplm_degen:' >> "+out_res2
            count_dcs = "grep '>' "+in_drpcull+" | wc -l"+" >> "+out_res2
            write_fs = "echo '#Filtered_lengths:' >> "+out_res2
            count_fs = "grep '>' "+in_fillen+" | wc -l"+" >> "+out_res2
            write_dfs = "echo '#Dereplicated_filtered_lengths:' >> "+out_res2
            count_dfs = "grep '>' "+in_drpfillen+" | wc -l"+" >> "+out_res2
            #Running
            print("\n# "+primer+" >>> Running . . .")
            os.system(write_primer)
            print("\n## Extracted hypervariable regions:")
            os.system(write_hs)
            print(count_hs)
            os.system(count_hs)
            print("\n## Dereplicated hypervariable regions:")
            os.system(write_dhs)
            print(count_dhs)
            os.system(count_dhs)
            print("\n## Filter homopolymers and degenerate bases:")
            os.system(write_cs)
            print(count_cs)
            os.system(count_cs)
            print("\n## Dereplicated filtered homopolymers and degenerate bases:")
            os.system(write_dcs)
            print(count_dcs)
            os.system(count_dcs)
            print("\n## Filtered sequence lengths:")
            os.system(write_fs)
            print(count_fs)
            os.system(count_fs)
            print("\n## Dereplicated filtered sequence lengths (cleaned datasets):")
            os.system(write_dfs)
            print(count_dfs)
            os.system(count_dfs)
            print("#---------------------------------------------------")
        fp.close()
        # endregion

        # region Create result files: tables and visualization
        print("\n#### Write number of record report file & create a bar plot ####")
        
        ## Create number of records dictionaries (e.g., gene_dict -> {gene1:{step1:num_rec ,step2:num_rec}, gene2:...})
        ### Mock community datasets
        out_res1 = thornfield+"/results/number_of_records/grep_mock.log.txt"
        fr1 = open(out_res1,"r")
        gene_dict = {}
        step_dict = {}
        for line in fr1:
            if line.strip()[0] == ">":
                gene = line.strip().split(">>")[1]
                gene_dict[gene] = ''
            elif line.strip()[0] == "#":
                step = line.strip().split("#")[1].replace(":","")
                step_dict[step] = ''
            else:
                num_rec = line.strip()
                step_dict[step] = num_rec
                gene_dict[gene] = step_dict.copy()
        fr1.close()
        ### Primer-extracted datasets
        out_res2 = thornfield+"/results/number_of_records/grep_primer.log.txt"
        fr2 = open(out_res2,"r")
        primer_dict = {}
        step_dict = {}
        for line in fr2:
            if line.strip()[0] == ">":
                primer = line.strip().split(">>")[1]
                primer_dict[primer] = ''
            elif line.strip()[0] == "#":
                step = line.strip().split("#")[1].replace(":","")
                step_dict[step] = ''
            else:
                num_rec = line.strip()
                step_dict[step] = num_rec
                primer_dict[primer] = step_dict.copy()
        fr2.close()

        ## Write number of records report table
        f = open(lowood,"r")
        header = f.readline()
        final_res = thornfield+"/results/number_of_records/num_records_REPORT.txt"
        fw = open(final_res,"w")
        fw.write("primer"+"\t"+"Downloaded NCBI data"+"\t"+"Dereplicated NCBI data"+"\t"+"Extracted hyperregion"+"\t"+"Dereplicated hyperregion"+"\t"+"Filtered hmplm and degen "+"\t"+"Dereplicated filtered hmplm and degen"+"\t"+"Filtered length"+"\t"+"Dereplicated filtered length (cleaned datasets)"+"\n")
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Write file
            fw.write(primer+"\t"+gene_dict[gene]["NCBI_data"]+"\t"+gene_dict[gene]["Dereplicated_NCBI_data"]+"\t"+primer_dict[primer]["Extracted_hyperregions"]+"\t"+primer_dict[primer]["Dereplicated_hyperregions"]+"\t"+primer_dict[primer]["Filtered_hmplm_degen"]+"\t"+primer_dict[primer]["Dereplicated_filtered_hmplm_degen"]+"\t"+primer_dict[primer]["Filtered_lengths"]+"\t"+primer_dict[primer]["Dereplicated_filtered_lengths"]+"\n")
        f.close()
        fw.close()
        print("#>>> Saved report as:", final_res)

        ## Create number of records bar plot
        df = pd.read_csv(final_res, sep='\t')
        out_bar = thornfield+"/results/number_of_records/num_records_BARPLOT.png"
        major_col = df.columns[0]
        df_long = df.melt(id_vars=major_col, var_name='Minor_Category', value_name='Value')
        plt.figure(figsize=(15, 10))
        sns.set_theme(style="darkgrid")
        sns.barplot(data=df_long, x=major_col, y='Value', hue='Minor_Category', edgecolor="black", linewidth=1, palette="Set3")
        plt.title('Number of records through processing steps', fontsize=20)
        plt.xlabel('Primers', fontsize = 18)
        plt.ylabel('Number of records', fontsize = 18)
        plt.xticks(fontsize=16, rotation=90)
        plt.yticks(fontsize=16)
        plt.legend(title='Preprocessing steps', title_fontsize=16, fontsize=14, facecolor='#eaeaf2', edgecolor='#eaeaf2', loc='center left', bbox_to_anchor=(1, 0.5))
        plt.savefig(out_bar, dpi=300, bbox_inches='tight')
        print("#>>> Saved bar plot as:", out_bar)
        # endregion
        # endregion

    # Observe sequence length for customized criteria (observe_length == "yes")
    elif observe_len == "yes":
        print("#======================== Observe length before setting customized criteria ==========================#")

        # region Save variables to input arguments log file #2 >> use in customlength function
        farg2 = open("inp_args2.log.txt","w")
        farg2.write("query"+"\t"+"primer"+"\t"+"out_dir"+"\t"+"taxa"+"\t"+"dl_period"+"\t"+"jobs"+"\t"+"threads"+"\n")
        farg2.write(gateshead+"\t"+lowood+"\t"+thornfield+"\t"+taxid+"\t"+my+"\t"+jobs+"\t"+threads)
        farg2.close()
        # endregion

        # region Export dereplicated hmplm and degen filtered datasets with "qiime tools export"
        print("\n#### Export dereplicated homopolymer and degenerate base filtered dataset ####")
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            # Files
            in_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_seq.qza"
            in_tax = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_tax.qza"
            out_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_seq"
            out_tax = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_tax"
            # Commands
            cmd_seq = "qiime tools export \\\n"+"--input-path "+in_seq+" \\\n"+"--output-path "+out_seq
            cmd_tax = "qiime tools export \\\n"+"--input-path "+in_tax+" \\\n"+"--output-path "+out_tax
            #Running
            print("\n# "+primer,">>> Running . . .")
            print(cmd_seq)
            os.system(cmd_seq)
            print(cmd_tax)
            os.system(cmd_tax)
        f.close()
        # endregion

        # region Get sequence length report of each record in the datasets
        print("\n#### Create sequence length report of each record in the datasets ####")

        ## Create primer-length dictionary and taxonomy-length dictionary
        ### plen_dict[primer][accession] = length --> {primer1:{accession1:len1, acc2:len2,...}, primer2...}
        ### ptax_dict[primer][accession]["Kingdom"] = Kingdom of taxa1 --> {primer1:{accession1:{Kingdom:taxa1, Phylum:taxa1,..., Species:taxa1}}, primer2...}
        plen_dict = {}
        ptax_dict = {}
        f = open(lowood, "r")
        header = f.readline()
        for line in f:
            # Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()

            # Files
            in_seq = thornfield + "/export/" + "txid" + taxid + "_" + primer + "_" + my + "_drp_cll_drp_seq/dna-sequences.fasta"
            in_tax = thornfield + "/export/" + "txid" + taxid + "_" + primer + "_" + my + "_drp_cll_drp_tax/taxonomy.tsv"
            out_len = thornfield + "/quality/seq_len/" + primer + "_prefilter_seqlen_REPORT.txt"

            # Running message
            print("\n# " + primer + " >>> Running . . .")

            # Create sequence length dictionary: seqlen_dict[accession] = length --> {accession1:len1, acc2:len2,...}
            seqlen_dict = {}
            records = list(SeqIO.parse(in_seq, "fasta"))
            for record in records:
                seqlen_dict[record.id] = len(record.seq)
            ## Add seqlen_dict to plen_dict
            plen_dict[primer] = seqlen_dict

            # Create taxonomy dictionary: taxa_dict[accession]["taxonomic level (e.g.,Class)"] = taxa --> {accession1:{Kingdom:taxa1, Phylum:taxa1,..., Species:taxa1}}
            taxa_dict = {}
            f_tax = open(in_tax, "r")
            f_tax.readline()
            for line in f_tax:
                ## subtaxa_dict["taxonomic level"] = taxa --> {Kingdom:taxa1, Phylum:taxa1,..., Species:taxa1}
                subtaxa_dict = {}
                accession = line.split("\t")[0]
                Kingdom = line.split("\t")[1].strip().split(";")[0].split("__")[1]
                subtaxa_dict["Kingdom"] = Kingdom
                Phylum = line.split("\t")[1].strip().split(";")[1].split("__")[1]
                subtaxa_dict["Phylum"] = Phylum
                Class = line.split("\t")[1].strip().split(";")[2].split("__")[1]
                subtaxa_dict["Class"] = Class
                Order = line.split("\t")[1].strip().split(";")[3].split("__")[1]
                subtaxa_dict["Order"] = Order
                Family = line.split("\t")[1].strip().split(";")[4].split("__")[1]
                subtaxa_dict["Family"] = Family
                Genus = line.split("\t")[1].strip().split(";")[5].split("__")[1]
                subtaxa_dict["Genus"] = Genus
                Species = line.split("\t")[1].strip().split(";")[5].split("__")[1] + " " + \
                          line.split("\t")[1].strip().split(";")[6].split("__")[1]
                subtaxa_dict["Species"] = Species
                taxa_dict[accession] = subtaxa_dict
            ## Add taxa_dict to ptax_dict
            ptax_dict[primer] = taxa_dict
            f_tax.close()

            # Write sequence length results
            num_line = 0
            fw_len = open(out_len, "w")
            fw_len.write(
                "Accession" + "\t" + "Kingdom" + "\t" + "Phylum" + "\t" + "Class" + "\t" + "Order" + "\t" + "Family" + "\t" + "Genus" + "\t" + "Species" + "\t" + "Seqlen" + "\n")
            for acc in ptax_dict[primer]:
                accession = acc
                Kingdom = ptax_dict[primer][acc]["Kingdom"]
                Phylum = ptax_dict[primer][acc]["Phylum"]
                Class = ptax_dict[primer][acc]["Class"]
                Order = ptax_dict[primer][acc]["Order"]
                Family = ptax_dict[primer][acc]["Family"]
                Genus = ptax_dict[primer][acc]["Genus"]
                Species = ptax_dict[primer][acc]["Species"]
                seqlen = plen_dict[primer][acc]
                fw_len.write(
                    accession + "\t" + Kingdom + "\t" + Phylum + "\t" + Class + "\t" + Order + "\t" + Family + "\t" + Genus + "\t" + Species + "\t" + str(
                        seqlen) + "\n")
                num_line += 1
            fw_len.close()
            print(">>> Saved as:", out_len)
            print("### Number of written records:", num_line)

        f.close()
        # endregion

        # region Get sequence length statistics and distribution visualization for ALL records
        # region Create sequence length dictionary
        seqlen_dict = {}
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            # Files
            in_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_seq/dna-sequences.fasta"
            in_tax = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_tax/taxonomy.tsv"
            # Create sequence length dictionary
            len_list = []
            records = list(SeqIO.parse(in_seq, "fasta"))
            for record in records:
                len_list.append(len(record.seq))
            seqlen_dict[primer] = len_list
        f.close()
        # endregion

        # region Calculate sequence length statistics
        out_res = thornfield+"/quality/seq_len/allprimer_prefilter_seqlen_STAT.txt"
        fw = open(out_res,"w")
        fw.write("Primer"+"\t"+"Total records"+"\t"+"Arithmetic mean"+"\t"+"sd"+"\t"+"Min"+"\t"+"Max"+"\t"+"P1"+"\t"+"P2.5"+"\t"+"P5"+"\t"+"P25"+"\t"+"P50"+"\t"+"P75"+"\t"+"P95"+"\t"+"P97.5"+"\t"+"P99"+"\t"+"Lower fence"+"\t"+"Upper fence"+"\n")
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            # Files
            length_array = np.array(seqlen_dict[primer])
            num_rec = len(seqlen_dict[primer])
            arthmean = np.mean(length_array)
            sd = np.std(length_array)
            min_len = np.min(length_array)
            max_len = np.max(length_array)
            p010 = np.percentile(length_array, 1)
            p025 = np.percentile(length_array, 2.5)
            p050 = np.percentile(length_array, 5)
            p250 = np.percentile(length_array, 25)
            p500 = np.percentile(length_array, 50)
            p750 = np.percentile(length_array, 75)
            p950 = np.percentile(length_array, 95)
            p975 = np.percentile(length_array, 97.5)
            p990 = np.percentile(length_array, 99)
            low_fen = p250-1.5*(p750-p250)
            up_fen = p750+1.5*(p750-p250)
            fw.write(primer+"\t"+str(num_rec)+"\t"+str(arthmean)+"\t"+str(sd)+"\t"+str(min_len)+"\t"+str(max_len)+"\t"+str(p010)+"\t"+str(p025)+"\t"+str(p050)+"\t"+str(p250)+"\t"+str(p500)+"\t"+str(p750)+"\t"+str(p950)+"\t"+str(p975)+"\t"+str(p990)+"\t"+str(low_fen)+"\t"+str(up_fen)+"\n")
        fw.close()    
        f.close()
        # endregion

        # region Write statistics table
        print("\n#### Report sequence length statistics and distribution before filtering: all records in the datasets ####")
        df = pd.read_csv(out_res, sep='\t')
        df_transposed = df.transpose()
        df_transposed.to_csv(out_res, sep='\t', header=False)
        print("#>>> Saved statistics table as:",out_res)
        # endregion

        # region Create visualization results: Box plot & Density plot
        ## Data
        data = seqlen_dict
        df = pd.DataFrame([{'Group': k, 'Value': v} for k, values in data.items() for v in values])

        ## Boxplot
        out_box = thornfield+"/quality/seq_len/prefilter_seqlen_BOXPLOT.png"
        plt.figure(figsize=(15, 10))
        sns.set_theme(style="darkgrid")
        sns.boxplot(data=df, x='Group', y='Value', palette="Set3", dodge=False)
        plt.title('Box plot of sequence length before filtering across all primers', fontsize=15)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.xlabel('Primer', fontsize=14)
        plt.ylabel('Sequence length', fontsize=14)
        plt.savefig(out_box, dpi=300, bbox_inches='tight')
        print("#>>> Saved boxplot as:",out_box)

        ## Density plot
        out_den = thornfield+"/quality/seq_len/prefilter_seqlen_DENSITY.png"
        plt.figure(figsize=(15, 10))
        sns.set_theme(style="darkgrid")
        ax = sns.kdeplot(data=df, x='Value', hue='Group', fill=True, common_norm=False, linewidth=1.5, alpha=0.7, palette="Set3", warn_singular=False)
        sns.move_legend(ax, "best", title='Primer', title_fontsize=14,fontsize=12,frameon=True, facecolor='white', framealpha=1, edgecolor='white')
        plt.title('Density plot of sequence length before filtering across all primers',fontsize=15)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.xlabel('Sequence length', fontsize=14)
        plt.ylabel('Density', fontsize=14)
        plt.savefig(out_den, dpi=300, bbox_inches='tight')
        print("#>>> Saved density plot as:",out_den)
        # endregion
        # endregion

        # region Get sequence length statistics and distribution visualization grouped by TAXA
        print("\n#### Report sequence length statistics and distribution before filtering: grouped by taxonomic levels ####")

        flen_dict = {}
        for primer in plen_dict.keys():
            # Running message
            print("\n# "+primer+" >>> Running . . .")

            # Create dictionaries
            ## flen_dict[primer]["uniq taxa1"] = [len1, len2, len3, ...] --> {primer:{taxa1:[len1, len2, len3, ...]}}
            ## e.g., flen_dict[primer]["Actinopteri"] = [len1, len2, len3, ...]
            i=0 #number of taxa
            taxlen_dict = {} #{taxa1:[len1, len2, len3, ...]}
            len_list = [] #[len1, len2, len3, ...]
            k_set = set()
            p_set = set()
            c_set = set()
            o_set = set()
            f_set = set()
            g_set = set()
            s_set = set()
            for access in plen_dict[primer].keys():
                k_set.add(ptax_dict[primer][access]["Kingdom"])
                p_set.add(ptax_dict[primer][access]["Phylum"])
                c_set.add(ptax_dict[primer][access]["Class"])
                o_set.add(ptax_dict[primer][access]["Order"])
                f_set.add(ptax_dict[primer][access]["Family"])
                g_set.add(ptax_dict[primer][access]["Genus"])
                s_set.add(ptax_dict[primer][access]["Species"])
            k_list = list(k_set)
            p_list = list(p_set)
            c_list = list(c_set)
            o_list = list(o_set)
            f_list = list(f_set)
            g_list = list(g_set)
            s_list = list(s_set)

            if level == "Kingdom":
                for uniqtax in k_list:
                    for access in plen_dict[primer].keys():
                        if ptax_dict[primer][access]["Kingdom"] == uniqtax:
                            len_list.append(plen_dict[primer][access])
                            taxlen_dict[uniqtax] = len_list
                    len_list = []
                flen_dict[primer] = taxlen_dict

            elif level == "Phylum":
                for uniqtax in p_list:
                    for access in plen_dict[primer].keys():
                        if ptax_dict[primer][access]["Phylum"] == uniqtax:
                            len_list.append(plen_dict[primer][access])
                            taxlen_dict[uniqtax] = len_list
                    len_list = []
                flen_dict[primer] = taxlen_dict

            elif level == "Class":
                for uniqtax in c_list:
                    for access in plen_dict[primer].keys():
                        if ptax_dict[primer][access]["Class"] == uniqtax:
                            len_list.append(plen_dict[primer][access])
                            taxlen_dict[uniqtax] = len_list
                    len_list = []
                flen_dict[primer] = taxlen_dict

            elif level == "Order":
                for uniqtax in o_list:
                    for access in plen_dict[primer].keys():
                        if ptax_dict[primer][access]["Order"] == uniqtax:
                            len_list.append(plen_dict[primer][access])
                            taxlen_dict[uniqtax] = len_list
                    len_list = []
                flen_dict[primer] = taxlen_dict 

            elif level == "Family":
                for uniqtax in f_list:
                    for access in plen_dict[primer].keys():
                        if ptax_dict[primer][access]["Family"] == uniqtax:
                            len_list.append(plen_dict[primer][access])
                            taxlen_dict[uniqtax] = len_list
                    len_list = []
                flen_dict[primer] = taxlen_dict 

            elif level == "Genus":
                for uniqtax in g_list:
                    for access in plen_dict[primer].keys():
                        if ptax_dict[primer][access]["Genus"] == uniqtax:
                            len_list.append(plen_dict[primer][access])
                            taxlen_dict[uniqtax] = len_list
                    len_list = []
                flen_dict[primer] = taxlen_dict 

            elif level == "Species":
                for uniqtax in s_list:
                    for access in plen_dict[primer].keys():
                        if ptax_dict[primer][access]["Species"] == uniqtax:
                            len_list.append(plen_dict[primer][access])
                            taxlen_dict[uniqtax] = len_list
                    len_list = []
                flen_dict[primer] = taxlen_dict 

            # Write statistics tables
            out_stat = thornfield+"/quality/seq_len/"+primer+"_prefilter_seqlen_STAT-"+level+".txt"
            fw_stat = open(out_stat,"w")
            fw_stat.write("Taxa"+"\t"+"Total records"+"\t"+"Seqlen arithmetic mean"+"\t"+"sd"+"\t"+"Min"+"\t"+"Max"+"\t"+"P1"+"\t"+"P2.5"+"\t"+"P5"+"\t"+"P25"+"\t"+"P50"+"\t"+"P75"+"\t"+"P95"+"\t"+"P97.5"+"\t"+"P99"+"\t"+"Lower fence"+"\t"+"Upper fence"+"\n")
            for taxa in flen_dict[primer].keys():
                length_list = flen_dict[primer][taxa]
                length_array = np.array(flen_dict[primer][taxa])
                taxonomy = taxa 
                num_rec = len(length_list)
                arthmean = np.mean(length_array)
                sd = np.std(length_array)
                min_len = np.min(length_array)
                max_len = np.max(length_array)
                p010 = np.percentile(length_array, 1)
                p025 = np.percentile(length_array, 2.5)
                p050 = np.percentile(length_array, 5)
                p250 = np.percentile(length_array, 25)
                p500 = np.percentile(length_array, 50)
                p750 = np.percentile(length_array, 75)
                p950 = np.percentile(length_array, 95)
                p975 = np.percentile(length_array, 97.5)
                p990 = np.percentile(length_array, 99)
                low_fen = p250-1.5*(p750-p250)
                up_fen = p750+1.5*(p750-p250)
                fw_stat.write(taxonomy+"\t"+str(num_rec)+"\t"+str(arthmean)+"\t"+str(sd)+"\t"+str(min_len)+"\t"+str(max_len)+"\t"+str(p010)+"\t"+str(p025)+"\t"+str(p050)+"\t"+str(p250)+"\t"+str(p500)+"\t"+str(p750)+"\t"+str(p950)+"\t"+str(p975)+"\t"+str(p990)+"\t"+str(low_fen)+"\t"+str(up_fen)+"\n")
                i+=1
            fw_stat.close()
            print("#>>> Saved statisitcs table as:",out_stat)
            print("### Number of taxa summarized:",str(i),level)
                
            # Create box plots
            data = flen_dict[primer]
            df = pd.DataFrame([{'Group': k, 'Value': v} for k, values in data.items() for v in values])
            out_box = thornfield+"/quality/seq_len/"+primer+"_prefilter_seqlen_BOXPLOT-"+level+".png"
            plt.figure(figsize=(15, 10))
            sns.set_theme(style="darkgrid")
            sns.boxplot(data=df, x='Group', y='Value', palette="Set3", dodge=False)
            plt.title('Box plot of sequence length before filtering grouped by '+level+" - "+primer, fontsize=15)
            plt.xticks(rotation=90, fontsize=12)
            plt.yticks(fontsize=12)
            plt.xlabel(level, fontsize=14)
            plt.ylabel('Sequence length', fontsize=14)
            plt.savefig(out_box, dpi=300, bbox_inches='tight')
            print("#>>> Saved box plot as:",out_box)
            
            ## Create density plots
            data = flen_dict[primer]
            df = pd.DataFrame([{'Group': k, 'Value': v} for k, values in data.items() for v in values])
            out_den = thornfield+"/quality/seq_len/"+primer+"_prefilter_seqlen_DENSITY-"+level+".png"
            plt.figure(figsize=(15, 10))
            sns.set_theme(style="darkgrid")
            ax = sns.kdeplot(data=df, x='Value', hue='Group', fill=True, common_norm=False, linewidth=1.5, alpha=0.7, palette="Set3", warn_singular=False)
            sns.move_legend(ax, "best", title='Taxa', title_fontsize=14,fontsize=12,frameon=True, facecolor='white', framealpha=1, edgecolor='white')
            plt.title('Density Plot of sequence length before filtering grouped by '+level+" - "+primer, fontsize=15)
            plt.xticks(fontsize=12)
            plt.yticks(fontsize=12)
            plt.xlabel('Sequence length', fontsize=14)
            plt.ylabel('Density', fontsize=14)
            plt.savefig(out_den, dpi=300, bbox_inches='tight')
            print("#>>> Saved density plot as:",out_den)
        # endregion

    # endregion
#=======================================================================================================================

# endregion
########################################################################################################################
# region Function #3: Customize length filtering criteria
                        
def customlength(args):
    # region Assign parsed arguments from command flags to python variables
    prilab = args.p_label
    minlen = args.p_min_length
    maxlen = args.p_max_length
    taxlen = args.p_fillen_by_taxon
    # endregion

    # region Retrieve variables from the previous command (processprimer)
    farg2 = open("inp_args2.log.txt","r")
    farg2.readline()
    for line in farg2:
        gateshead = line.split("\t")[0]
        lowood = line.split("\t")[1]
        thornfield = line.split("\t")[2]
        taxid = line.split("\t")[3]
        my = line.split("\t")[4]
        jobs = line.split("\t")[5]
        threads = line.split("\t")[6].strip()
    farg2.close()
    # endregion

    # Vary parameter #1: Filter length with customized parameters for ALL TAXA
    if minlen!=0 and maxlen!=0 and taxlen=="" and len(minlen)==len(maxlen) and len(maxlen)==len(prilab):
        # region Generate QIIME2 commands: Filter sequence length with customized criteria for ALL TAXA
        print("#============================== Filter length with customized criteria ===============================#")

        ## Filter sequence length with "qiime rescript filter-seqs-length"
        print("\n#### Filter sequence length for all taxa with customized parameter ####")
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            # Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            # Files
            in_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_seq.qza"
            out_seq = thornfield+"/quality/seq_len/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_len_seq.qza"
            out_seq_DISC = thornfield+"/quality/seq_len/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_len_DISC_seq.qza"            
            # Get customized variables from input command arguments
            for i in range(len(prilab)):
                if primer == prilab[i]:
                    cmd = "qiime rescript filter-seqs-length \\\n"+"--i-sequences "+in_seq+" \\\n"+"--p-global-min "+str(minlen[i])+" \\\n"+"--p-global-max "+str(maxlen[i])+" \\\n"+"--o-filtered-seqs "+out_seq+" \\\n"+"--o-discarded-seqs "+out_seq_DISC
                    #Running
                    print("\n# "+primer,">>> Running . . . ")
                    print(cmd)
                    os.system(cmd)
        f.close()

        ## Dereplicate datasets with "qiime rescript dereplicate"
        print("\n#### Dereplicate length filtered datasets (cleaned datasets) ####")
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Files
            in_seq = thornfield+"/quality/seq_len/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_len_seq.qza" 
            in_tax = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_tax.qza"
            out_seq = thornfield+"/clean_datasets/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_seq.qza"
            out_tax = thornfield+"/clean_datasets/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_tax.qza"
            #Command
            cmd = "qiime rescript dereplicate \\\n"+"--i-sequences "+in_seq+" \\\n"+"--i-taxa "+in_tax+" \\\n"+"--p-mode "+"uniq"+" \\\n"+"--p-threads "+str(threads)+ " \\\n"+"--o-dereplicated-sequences "+out_seq+" \\\n"+"--o-dereplicated-taxa "+out_tax
            #Running
            print("\n# "+primer,">>> Running . . . ")
            print(cmd)
            os.system(cmd)
        f.close()

        # endregion

        # region Summarize number of records through filtering steps (customized length filtering - ALL TAXA)
        print("#=========================== Count number of records (customized length) =============================#")

        # region Export artifact files from QIIME2 with "qiime tools export"
        ## Mock community datasets
        print("\n#### Export artifact files ####")
        fg = open(gateshead,"r")
        header = fg.readline()
        for line in fg:
            #Variables
            gene = line.split("\t")[0]
            query = line.split("\t")[1].strip()
            # Files
            ## Inputs
            in_NCBI_seq = "NCBI_data/"+"txid"+taxid+"_"+gene+"_"+my+"_seq.qza"
            in_NCBI_tax = "NCBI_data/"+"txid"+taxid+"_"+gene+"_"+my+"_tax.qza"
            in_drpNCBI_seq = "NCBI_data/"+"txid"+taxid+"_"+gene+"_"+my+"_drp_seq.qza"
            in_drpNCBI_tax = "NCBI_data/"+"txid"+taxid+"_"+gene+"_"+my+"_drp_tax.qza"
            ## Outputs
            out_NCBI_seq = thornfield+"/export/"+"txid"+taxid+"_"+gene+"_"+my+"_seq"
            out_NCBI_tax = thornfield+"/export/"+"txid"+taxid+"_"+gene+"_"+my+"_tax"
            out_drpNCBI_seq = thornfield+"/export/"+"txid"+taxid+"_"+gene+"_"+my+"_drp_seq"
            out_drpNCBI_tax = thornfield+"/export/"+"txid"+taxid+"_"+gene+"_"+my+"_drp_tax"
            # Commands 
            cmd_NCBI_seq = "qiime tools export \\\n"+"--input-path "+in_NCBI_seq+" \\\n"+"--output-path "+out_NCBI_seq
            cmd_NCBI_tax = "qiime tools export \\\n"+"--input-path "+in_NCBI_tax+" \\\n"+"--output-path "+out_NCBI_tax
            cmd_drpNCBI_seq = "qiime tools export \\\n"+"--input-path "+in_drpNCBI_seq+" \\\n"+"--output-path "+out_drpNCBI_seq
            cmd_drpNCBI_tax = "qiime tools export \\\n"+"--input-path "+in_drpNCBI_tax+" \\\n"+"--output-path "+out_drpNCBI_tax
            #Running
            print("\n# "+gene,">>> Running . . .")
            print("\n## NCBI data:")
            print(cmd_NCBI_seq)
            os.system(cmd_NCBI_seq)
            print(cmd_NCBI_tax)
            os.system(cmd_NCBI_tax)
            print("\n## Dereplicated NCBI data:")
            print(cmd_drpNCBI_seq)
            os.system(cmd_drpNCBI_seq)
            print(cmd_drpNCBI_tax)
            os.system(cmd_drpNCBI_tax)
            print("#-------------------------------------------------------------------------------------------------#")
        fg.close()

        ## Primer-extracted datasets
        fp = open(lowood,"r")
        header = fp.readline()
        for line in fp:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            # Files
            ## Inputs
            in_hyper_seq = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_seq.qza"
            in_drphyper_seq = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_seq.qza"
            in_drphyper_tax = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_tax.qza"
            in_cull_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_seq.qza"
            in_drpcull_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_seq.qza"
            in_drpcull_tax = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_tax.qza"
            in_fillen_seq = thornfield+"/quality/seq_len/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_len_seq.qza"
            in_drpfillen_seq = thornfield+"/clean_datasets/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_seq.qza"
            in_drpfillen_tax = thornfield+"/clean_datasets/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_tax.qza"
            ## Outputs
            out_hyper_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_seq"
            out_drphyper_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_seq"
            out_drphyper_tax = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_tax"
            out_cull_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_seq"
            out_drpcull_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_seq"
            out_drpcull_tax = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_tax"
            out_fillen_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_len_seq"
            out_drpfillen_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_seq"
            out_drpfillen_tax = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_tax"
            # Commands
            cmd_hyper_seq = "qiime tools export \\\n"+"--input-path "+in_hyper_seq+" \\\n"+"--output-path "+out_hyper_seq
            cmd_drphyper_seq = "qiime tools export \\\n"+"--input-path "+in_drphyper_seq+" \\\n"+"--output-path "+out_drphyper_seq
            cmd_drphyper_tax = "qiime tools export \\\n"+"--input-path "+in_drphyper_tax+" \\\n"+"--output-path "+out_drphyper_tax
            cmd_cull_seq = "qiime tools export \\\n"+"--input-path "+in_cull_seq+" \\\n"+"--output-path "+out_cull_seq
            cmd_drpcull_seq = "qiime tools export \\\n"+"--input-path "+in_drpcull_seq+" \\\n"+"--output-path "+out_drpcull_seq
            cmd_drpcull_tax = "qiime tools export \\\n"+"--input-path "+in_drpcull_tax+" \\\n"+"--output-path "+out_drpcull_tax
            cmd_fillen_seq = "qiime tools export \\\n"+"--input-path "+in_fillen_seq+" \\\n"+"--output-path "+out_fillen_seq
            cmd_drpfillen_seq = "qiime tools export \\\n"+"--input-path "+in_drpfillen_seq+" \\\n"+"--output-path "+out_drpfillen_seq
            cmd_drpfillen_tax = "qiime tools export \\\n"+"--input-path "+in_drpfillen_tax+" \\\n"+"--output-path "+out_drpfillen_tax
            #Running
            print("\n# "+primer,">>> Running . . .")
            print("\n## Extracted hypervariable regions:")
            print(cmd_hyper_seq)
            os.system(cmd_hyper_seq)
            print("\n## Dereplicated hypervariable regions:")
            print(cmd_drphyper_seq)
            os.system(cmd_drphyper_seq)
            print(cmd_drphyper_tax)
            os.system(cmd_drphyper_tax)
            print("\n## Filtered homopolymers & degenerate bases:")
            print(cmd_cull_seq)
            os.system(cmd_cull_seq)
            print("\n## Dereplicated cull-seqs:")
            print(cmd_drpcull_seq)
            os.system(cmd_drpcull_seq)
            print(cmd_drpcull_tax)
            os.system(cmd_drpcull_tax)
            print("\n## Filtered sequence lengths:")
            print(cmd_fillen_seq)
            os.system(cmd_fillen_seq)
            print("\n## Dereplicated filtered lengths (cleaned datasets):")
            print(cmd_drpfillen_seq)
            os.system(cmd_drpfillen_seq)
            print(cmd_drpfillen_tax)
            os.system(cmd_drpfillen_tax)
            print("#-------------------------------------------------------------------------------------------------#")
        fp.close()
        # endregion

        # region Create new directories
        newdir1 = "mkdir "+thornfield+"/results"
        newdir2 = "mkdir "+thornfield+"/results/number_of_records"
        os.system(newdir1)
        os.system(newdir2)
        # endregion

        # region Count number of records through processing steps
        ## Count records in mock community datasets with "grep"
        print("\n#### Count number of records: Mock community datasets ####")
        fg = open(gateshead,"r")
        header = fg.readline()
        for line in fg:
            #Variables
            gene = line.split("\t")[0]
            query = line.split("\t")[1].strip()
            # Files
            in_NCBI_seq = thornfield+"/export/"+"txid"+taxid+"_"+gene+"_"+my+"_seq/dna-sequences.fasta"
            in_drpNCBI_seq = thornfield+"/export/"+"txid"+taxid+"_"+gene+"_"+my+"_drp_seq/dna-sequences.fasta"
            out_res1 = thornfield+"/results/number_of_records/grep_mock.log.txt"
            # Commands
            write_gene = "echo "+"'>>"+gene+"'"+" >> "+out_res1
            write_ns = "echo '#NCBI_data' >> "+out_res1
            count_ns = "grep '>' "+in_NCBI_seq+" | wc -l"+" >> "+out_res1
            write_dns = "echo '#Dereplicated_NCBI_data' >> "+out_res1
            count_dns = "grep '>' "+in_drpNCBI_seq+" | wc -l"+" >> "+out_res1
            #Running
            print("\n# "+gene,">>> Running . . .")
            os.system(write_gene)
            print("\n## NCBI data:")
            os.system(write_ns)
            print(count_ns)
            os.system(count_ns)
            print("\n## Dereplicated NCBI data:")
            os.system(write_dns)
            print(count_dns)
            os.system(count_dns)
            print("#-------------------------------------------------------------------------------------------------#")
        fg.close()

        ## Count records in primer-extracted datasets with "grep"
        print("\n#### Count number of records: Primer datasets ####")
        fp = open(lowood,"r")
        header = fp.readline()
        for line in fp:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            # Files   
            in_hyper_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_seq/dna-sequences.fasta"
            in_drphyper_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_seq/dna-sequences.fasta"
            in_cull_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_seq/dna-sequences.fasta"
            in_drpcull_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_seq/dna-sequences.fasta"
            in_fillen_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_len_seq/dna-sequences.fasta"
            in_drpfillen_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_seq/dna-sequences.fasta"
            out_res2 = thornfield+"/results/number_of_records/grep_primer.log.txt"
            # Commands
            write_primer = "echo "+"'>>"+primer+"'"+" >> "+out_res2
            write_hs = "echo '#Extracted_hyperregions:' >> "+out_res2
            count_hs = "grep '>' "+in_hyper_seq+" | wc -l"+" >> "+out_res2
            write_dhs = "echo '#Dereplicated_hyperregions:' >> "+out_res2
            count_dhs = "grep '>' "+in_drphyper_seq+" | wc -l"+" >> "+out_res2
            write_cs = "echo '#Filtered_hmplm_degen:' >> "+out_res2
            count_cs = "grep '>' "+in_cull_seq+" | wc -l"+" >> "+out_res2
            write_dcs = "echo '#Dereplicated_filtered_hmplm_degen:' >> "+out_res2
            count_dcs = "grep '>' "+in_drpcull_seq+" | wc -l"+" >> "+out_res2
            write_fs = "echo '#Filtered_lengths:' >> "+out_res2
            count_fs = "grep '>' "+in_fillen_seq+" | wc -l"+" >> "+out_res2
            write_dfs = "echo '#Dereplicated_filtered_lengths:' >> "+out_res2
            count_dfs = "grep '>' "+in_drpfillen_seq+" | wc -l"+" >> "+out_res2
            #Running
            print("\n# "+primer,">>> Running . . .")
            os.system(write_primer)
            print("\n## Extracted hypervariable regions:")
            os.system(write_hs)
            print(count_hs)
            os.system(count_hs)
            print("\n## Dereplicated hypervariable regions:")
            os.system(write_dhs)
            print(count_dhs)
            os.system(count_dhs)
            print("\n## Filter homopolymers and degenerate bases:")
            os.system(write_cs)
            print(count_cs)
            os.system(count_cs)
            print("\n## Dereplicated filtered homopolymers and degenerate bases:")
            os.system(write_dcs)
            print(count_dcs)
            os.system(count_dcs)
            print("\n## Filtered sequence lengths:")
            os.system(write_fs)
            print(count_fs)
            os.system(count_fs)
            print("\n## Dereplicated filtered sequence lengths (cleaned datasets):")
            os.system(write_dfs)
            print(count_dfs)
            os.system(count_dfs)
            print("#-------------------------------------------------------------------------------------------------#")
        fp.close()
        # endregion

        # region Create result files: tables and visualization
        ## Create number of records dictionaries: (gene_dict -> {gene1:{step1:num_rec ,step2:num_rec}, gene2:...})
        ### Mock community datasets
        out_res1 = thornfield+"/results/number_of_records/grep_mock.log.txt"
        fr1 = open(out_res1,"r")
        gene_dict = {}
        step_dict = {}
        for line in fr1:
            if line.strip()[0] == ">":
                gene = line.strip().split(">>")[1]
                gene_dict[gene] = ''
            elif line.strip()[0] == "#":
                step = line.strip().split("#")[1].replace(":","")
                step_dict[step] = ''
            else:
                num_rec = line.strip()
                step_dict[step] = num_rec
                gene_dict[gene] = step_dict.copy()
        fr1.close()
        ### Primer-extracted datasets
        out_res2 = thornfield+"/results/number_of_records/grep_primer.log.txt"
        fr2 = open(out_res2,"r")
        primer_dict = {}
        step_dict = {}
        for line in fr2:
            if line.strip()[0] == ">":
                primer = line.strip().split(">>")[1]
                primer_dict[primer] = ''
            elif line.strip()[0] == "#":
                step = line.strip().split("#")[1].replace(":","")
                step_dict[step] = ''
            else:
                num_rec = line.strip()
                step_dict[step] = num_rec
                primer_dict[primer] = step_dict.copy()
        fr2.close()

        ## Write number of records report table
        f = open(lowood,"r")
        header = f.readline()
        final_res = thornfield+"/results/number_of_records/num_records_REPORT.txt"
        fw = open(final_res,"w")
        fw.write("primer"+"\t"+"Downloaded NCBI data"+"\t"+"Dereplicated NCBI data"+"\t"+"Extracted hyperregion"+"\t"+"Dereplicated hyperregion"+"\t"+"Filtered hmplm and degen "+"\t"+"Dereplicated filtered hmplm and degen"+"\t"+"Filtered length"+"\t"+"Dereplicated filtered length (cleaned datasets)"+"\n")
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Write file
            fw.write(primer+"\t"+gene_dict[gene]["NCBI_data"]+"\t"+gene_dict[gene]["Dereplicated_NCBI_data"]+"\t"+primer_dict[primer]["Extracted_hyperregions"]+"\t"+primer_dict[primer]["Dereplicated_hyperregions"]+"\t"+primer_dict[primer]["Filtered_hmplm_degen"]+"\t"+primer_dict[primer]["Dereplicated_filtered_hmplm_degen"]+"\t"+primer_dict[primer]["Filtered_lengths"]+"\t"+primer_dict[primer]["Dereplicated_filtered_lengths"]+"\n")
        f.close()
        fw.close()

        ## Create number of records bar plot
        df = pd.read_csv(final_res, sep='\t')
        out_bar = thornfield+"/results/number_of_records/num_records_BARPLOT.png"
        major_col = df.columns[0]
        df_long = df.melt(id_vars=major_col, var_name='Minor_Category', value_name='Value')
        plt.figure(figsize=(15, 10))
        sns.set_theme(style="darkgrid")
        sns.barplot(data=df_long, x=major_col, y='Value', hue='Minor_Category', edgecolor="black", linewidth=1, palette="Set3")
        plt.title('Number of records through processing steps', fontsize=20)
        plt.xlabel('Primers', fontsize = 18)
        plt.ylabel('Number of records', fontsize = 18)
        plt.xticks(fontsize=16, rotation=90)
        plt.yticks(fontsize=16)
        plt.legend(title='Preprocessing steps', title_fontsize=16, fontsize=14, facecolor='#eaeaf2', edgecolor='#eaeaf2', loc='center left', bbox_to_anchor=(1, 0.5))
        plt.savefig(out_bar, dpi=300, bbox_inches='tight')
        print("#>>> Saved bar plot as:", out_bar)
        # endregion
        # endregion

    # Vary parameter #2: Filter length with customized parameters SPECIFIC to TAXA
    elif taxlen!="" and minlen==0 and maxlen==0 and len(taxlen)==len(prilab):
        # region Generate QIIME2 commands: Filter sequence length with customized criteria specific to taxa
        print("#==================== Filter length with customized criteria specific to taxa ========================#")

        # region Filter sequence length by taxon with "qiime rescript filter-seqs-length-by-taxon"
        print("\n#### Filter sequence length by taxon with customized parameter ####")
        out_last = {} #Keep final output file paths
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Files
            in_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_seq.qza"
            in_tax = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_tax.qza"
            #Get customized filtering criteria files from input command arguments
            for i in range(len(prilab)):
                if primer == prilab[i]:
                    in_filter = taxlen[i] 
            ff = open(in_filter,"r")
            ff.readline()
            #Generate and execute the command (cascade method)
            print("\n# "+primer,">>> Running . . .")
            j=0
            for line in ff:
                if j == 0: #first round
                    #Variables
                    taxa_name = line.strip().split("\t")[0]
                    taxa_pref = line.strip().split("\t")[0][0:3]
                    input_seq = in_seq
                    input_tax = in_tax
                    min_lens = line.strip().split("\t")[1]
                    max_lens = line.strip().split("\t")[2]
                    #Files
                    output_seq = thornfield+"/quality/seq_len/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_len-"+taxa_pref+"_seq.qza"
                    output_seq_disc = thornfield+"/quality/seq_len/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_len-"+taxa_pref+"_seq_DISC.qza"
                    #Command
                    filter_seq_cmd = "qiime rescript filter-seqs-length-by-taxon \\\n --i-sequences "+input_seq+" \\\n --i-taxonomy "+input_tax+" \\\n --p-labels "+taxa_name+" \\\n --p-min-lens "+min_lens+" \\\n --p-max-lens "+max_lens+" \\\n --o-filtered-seqs "+output_seq+" \\\n --o-discarded-seqs "+output_seq_disc
                    #Running
                    print("\n## "+taxa_name)
                    print(filter_seq_cmd)
                    os.system(filter_seq_cmd)
                    j += 1

                elif j > 0: #cascade filtering using files from previous step as inputs
                    #Variables
                    taxa_name = line.strip().split("\t")[0]
                    taxa_pref = line.strip().split("\t")[0][0:3]
                    input_seq = output_seq
                    input_tax = input_tax
                    min_lens = line.strip().split("\t")[1]
                    max_lens = line.strip().split("\t")[2]
                    #Files
                    output_seq = str(output_seq).replace("_seq.qza","-"+taxa_pref+"_seq.qza")
                    output_seq_disc = str(output_seq).replace("_seq.qza","_seq_DISC.qza")
                    #Command
                    filter_seq_cmd = "qiime rescript filter-seqs-length-by-taxon \\\n --i-sequences "+input_seq+" \\\n --i-taxonomy "+input_tax+" \\\n --p-labels "+taxa_name+" \\\n --p-min-lens "+min_lens+" \\\n --p-max-lens "+max_lens+" \\\n --o-filtered-seqs "+output_seq+" \\\n --o-discarded-seqs "+output_seq_disc
                    #Running
                    print("\n## "+taxa_name)
                    print(filter_seq_cmd)
                    os.system(filter_seq_cmd)
                    out_last[primer] = output_seq #Keeping final filtered files in a dictionary; key = primer, values = file path
            ff.close()
            print("#---------------------------------------------------------------------------------------------------------------#")
        f.close()
        # endregion

        # region Keep final output file paths as a tab-delimited file
        ## Write file
        print("\n# Writing final output file paths")
        fl = open(thornfield+"/quality/seq_len/out_last.txt","w")
        fl.write("primer"+"\t"+"file_path"+"\n")
        for key in out_last.keys():
            fl.write(key+"\t"+out_last[key]+"\n")
        print("## >>> Saved final output file paths as: "+thornfield+"/quality/seq_len/out_last.txt")
        fl.close()
        ## Make dictionary
        out_last = {}
        fl = open(thornfield+"/quality/seq_len/out_last.txt","r")
        fl.readline()
        for line in fl:
            primer = line.split("\t")[0]
            file = line.split("\t")[1].strip()
            out_last[primer] = file
        fl.close()
        print("#-----------------------------------------------------------------------------------------------------#")
        # endregion

        # region Dereplicate datasets with "qiime rescript dereplicate"
        print("\n#### Dereplicate length filtered datasets (cleaned datasets) ####")
        f = open(lowood,"r")
        header = f.readline()
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Files
            in_seq = out_last[primer]
            in_tax = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_tax.qza"
            out_seq = thornfield+"/clean_datasets/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_seq.qza"
            out_tax = thornfield+"/clean_datasets/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_tax.qza"
            #Command
            cmd = "qiime rescript dereplicate \\\n"+"--i-sequences "+in_seq+" \\\n"+"--i-taxa "+in_tax+" \\\n"+"--p-mode "+"uniq"+" \\\n"+"--p-threads "+str(threads)+ " \\\n"+"--o-dereplicated-sequences "+out_seq+" \\\n"+"--o-dereplicated-taxa "+out_tax
            #Running
            print("\n# "+primer,">>> Running . . . ")
            print(cmd)
            os.system(cmd)
        f.close()
        # endregion

        # endregion

        # region Summarize number of records through filtering steps (customized length filtering - TAXA SPECIFIC)
        print("#====================== Count number of records (customized length by taxa) ==========================#")

        # region Export artifact files from QIIME2 with "qiime tools export"
        ## Mock community datasets
        print("\n#### Export artifact files ####")
        fg = open(gateshead,"r")
        header = fg.readline()
        for line in fg:
            #Variables
            gene = line.split("\t")[0]
            query = line.split("\t")[1].strip()
            # Files
            ## Inputs
            in_NCBI_seq = "NCBI_data/"+"txid"+taxid+"_"+gene+"_"+my+"_seq.qza"
            in_NCBI_tax = "NCBI_data/"+"txid"+taxid+"_"+gene+"_"+my+"_tax.qza"
            in_drpNCBI_seq = "NCBI_data/"+"txid"+taxid+"_"+gene+"_"+my+"_drp_seq.qza"
            in_drpNCBI_tax = "NCBI_data/"+"txid"+taxid+"_"+gene+"_"+my+"_drp_tax.qza"
            ## Outputs
            out_NCBI_seq = thornfield+"/export/"+"txid"+taxid+"_"+gene+"_"+my+"_seq"
            out_NCBI_tax = thornfield+"/export/"+"txid"+taxid+"_"+gene+"_"+my+"_tax"
            out_drpNCBI_seq = thornfield+"/export/"+"txid"+taxid+"_"+gene+"_"+my+"_drp_seq"
            out_drpNCBI_tax = thornfield+"/export/"+"txid"+taxid+"_"+gene+"_"+my+"_drp_tax"
            # Commands 
            cmd_NCBI_seq = "qiime tools export \\\n"+"--input-path "+in_NCBI_seq+" \\\n"+"--output-path "+out_NCBI_seq
            cmd_NCBI_tax = "qiime tools export \\\n"+"--input-path "+in_NCBI_tax+" \\\n"+"--output-path "+out_NCBI_tax
            cmd_drpNCBI_seq = "qiime tools export \\\n"+"--input-path "+in_drpNCBI_seq+" \\\n"+"--output-path "+out_drpNCBI_seq
            cmd_drpNCBI_tax = "qiime tools export \\\n"+"--input-path "+in_drpNCBI_tax+" \\\n"+"--output-path "+out_drpNCBI_tax
            #Running
            print("\n# "+gene,">>> Running . . .")
            print("\n## NCBI data:")
            print(cmd_NCBI_seq)
            os.system(cmd_NCBI_seq)
            print(cmd_NCBI_tax)
            os.system(cmd_NCBI_tax)
            print("\n## Dereplicated NCBI data:")
            print(cmd_drpNCBI_seq)
            os.system(cmd_drpNCBI_seq)
            print(cmd_drpNCBI_tax)
            os.system(cmd_drpNCBI_tax)
            print("#----------------------------------------------------------------------------------------------------------------#")
        fg.close()

        ## Primer-extracted datasets
        ### Make dictionary for final output file paths during filter sequence length step
        out_last = {}
        fl = open(thornfield+"/quality/seq_len/out_last.txt","r")
        fl.readline()
        for line in fl:
            primer = line.split("\t")[0]
            file = line.split("\t")[1].strip()
            out_last[primer] = file
        fl.close()
        ### Get parameters & generate exporting commands
        fp = open(lowood,"r")
        header = fp.readline()
        for line in fp:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            # Files
            ## Inputs
            in_hyper_seq = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_seq.qza"
            in_drphyper_seq = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_seq.qza"
            in_drphyper_tax = thornfield+"/extract/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_tax.qza"
            in_cull_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_seq.qza"
            in_drpcull_seq = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_seq.qza"
            in_drpcull_tax = thornfield+"/quality/hmplm_degen/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_tax.qza"
            in_fillen_seq = out_last[primer]
            in_drpfillen_seq = thornfield+"/clean_datasets/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_seq.qza"
            in_drpfillen_tax = thornfield+"/clean_datasets/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_tax.qza"
            ## Outputs
            out_hyper_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_seq"
            out_drphyper_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_seq"
            out_drphyper_tax = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_tax"
            out_cull_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_seq"
            out_drpcull_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_seq"
            out_drpcull_tax = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_tax"
            out_fillen_seq = out_last[primer].replace("quality/seq_len","export").replace(".qza","")
            out_drpfillen_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_seq"
            out_drpfillen_tax = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_tax"
            # Commands
            cmd_hyper_seq = "qiime tools export \\\n"+"--input-path "+in_hyper_seq+" \\\n"+"--output-path "+out_hyper_seq
            cmd_drphyper_seq = "qiime tools export \\\n"+"--input-path "+in_drphyper_seq+" \\\n"+"--output-path "+out_drphyper_seq
            cmd_drphyper_tax = "qiime tools export \\\n"+"--input-path "+in_drphyper_tax+" \\\n"+"--output-path "+out_drphyper_tax
            cmd_cull_seq = "qiime tools export \\\n"+"--input-path "+in_cull_seq+" \\\n"+"--output-path "+out_cull_seq
            cmd_drpcull_seq = "qiime tools export \\\n"+"--input-path "+in_drpcull_seq+" \\\n"+"--output-path "+out_drpcull_seq
            cmd_drpcull_tax = "qiime tools export \\\n"+"--input-path "+in_drpcull_tax+" \\\n"+"--output-path "+out_drpcull_tax
            cmd_fillen_seq = "qiime tools export \\\n"+"--input-path "+in_fillen_seq+" \\\n"+"--output-path "+out_fillen_seq
            cmd_drpfillen_seq = "qiime tools export \\\n"+"--input-path "+in_drpfillen_seq+" \\\n"+"--output-path "+out_drpfillen_seq
            cmd_drpfillen_tax = "qiime tools export \\\n"+"--input-path "+in_drpfillen_tax+" \\\n"+"--output-path "+out_drpfillen_tax
            #Running
            print("\n# "+primer,">>> Running . . .")
            print("\n## Extracted hypervariable regions:")
            print(cmd_hyper_seq)
            os.system(cmd_hyper_seq)
            print("\n## Dereplicated hypervariable regions:")
            print(cmd_drphyper_seq)
            os.system(cmd_drphyper_seq)
            print(cmd_drphyper_tax)
            os.system(cmd_drphyper_tax)
            print("\n## Filtered homopolymers & degenerate bases:")
            print(cmd_cull_seq)
            os.system(cmd_cull_seq)
            print("\n## Dereplicated cull-seqs:")
            print(cmd_drpcull_seq)
            os.system(cmd_drpcull_seq)
            print(cmd_drpcull_tax)
            os.system(cmd_drpcull_tax)
            print("\n## Filtered sequence lengths:")
            print(cmd_fillen_seq)
            os.system(cmd_fillen_seq)
            print("\n## Dereplicated filtered lengths (cleaned datasets):")
            print(cmd_drpfillen_seq)
            os.system(cmd_drpfillen_seq)
            print(cmd_drpfillen_tax)
            os.system(cmd_drpfillen_tax)
            print("#-------------------------------------------------------------------------------------------------#")
        fp.close()
        # endregion

        # region Create new directories
        newdir1 = "mkdir "+thornfield+"/results"
        newdir2 = "mkdir "+thornfield+"/results/number_of_records"
        os.system(newdir1)
        os.system(newdir2)
        # endregion

        # region Count number of records through processing steps
        ## Count records in mock community datasets with "grep"
        print("\n#### Count number of records: Mock community datasets ####")
        fg = open(gateshead,"r")
        header = fg.readline()
        for line in fg:
            #Variables
            gene = line.split("\t")[0]
            query = line.split("\t")[1].strip()
            # Files
            in_NCBI_seq = thornfield+"/export/"+"txid"+taxid+"_"+gene+"_"+my+"_seq/dna-sequences.fasta"
            in_drpNCBI_seq = thornfield+"/export/"+"txid"+taxid+"_"+gene+"_"+my+"_drp_seq/dna-sequences.fasta"
            out_res1 = thornfield+"/results/number_of_records/grep_mock.log.txt"
            # Commands
            write_gene = "echo "+"'>>"+gene+"'"+" >> "+out_res1
            write_ns = "echo '#NCBI_data' >> "+out_res1
            count_ns = "grep '>' "+in_NCBI_seq+" | wc -l"+" >> "+out_res1
            write_dns = "echo '#Dereplicated_NCBI_data' >> "+out_res1
            count_dns = "grep '>' "+in_drpNCBI_seq+" | wc -l"+" >> "+out_res1
            #Running
            print("\n# "+gene,">>> Running . . .")
            os.system(write_gene)
            print("\n## NCBI data:")
            os.system(write_ns)
            print(count_ns)
            os.system(count_ns)
            print("\n## Dereplicated NCBI data:")
            os.system(write_dns)
            print(count_dns)
            os.system(count_dns)
            print("#-------------------------------------------------------------------------------------------------#")
        fg.close()

        ## Count records in primer-extracted datasets with "grep"
        print("\n#### Count number of records: Primer datasets ####")
        ### Make dictionary for output file paths during filter sequence length step
        out_last = {}
        fl = open(thornfield+"/quality/seq_len/out_last.txt","r")
        fl.readline()
        for line in fl:
            primer = line.split("\t")[0]
            file = line.split("\t")[1].strip()
            out_last[primer] = file
        fl.close()
        ### Write results
        fp = open(lowood,"r")
        header = fp.readline()
        for line in fp:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Files
            in_hyper_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_seq/dna-sequences.fasta"
            in_drphyper_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_seq/dna-sequences.fasta"
            in_cull_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_seq/dna-sequences.fasta"
            in_drpcull_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_drp_cll_drp_seq/dna-sequences.fasta"
            in_fillen_seq = out_last[primer].replace("quality/seq_len","export").replace(".qza","/dna-sequences.fasta")
            in_drpfillen_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_seq/dna-sequences.fasta"
            out_res2 = thornfield+"/results/number_of_records/grep_primer.log.txt"
            #Commands
            write_primer = "echo "+"'>>"+primer+"'"+" >> "+out_res2
            write_hs = "echo '#Extracted_hyperregions:' >> "+out_res2
            count_hs = "grep '>' "+in_hyper_seq+" | wc -l"+" >> "+out_res2
            write_dhs = "echo '#Dereplicated_hyperregions:' >> "+out_res2
            count_dhs = "grep '>' "+in_drphyper_seq+" | wc -l"+" >> "+out_res2
            write_cs = "echo '#Filtered_hmplm_degen:' >> "+out_res2
            count_cs = "grep '>' "+in_cull_seq+" | wc -l"+" >> "+out_res2
            write_dcs = "echo '#Dereplicated_filtered_hmplm_degen:' >> "+out_res2
            count_dcs = "grep '>' "+in_drpcull_seq+" | wc -l"+" >> "+out_res2
            write_fs = "echo '#Filtered_lengths:' >> "+out_res2
            count_fs = "grep '>' "+in_fillen_seq+" | wc -l"+" >> "+out_res2
            write_dfs = "echo '#Dereplicated_filtered_lengths:' >> "+out_res2
            count_dfs = "grep '>' "+in_drpfillen_seq+" | wc -l"+" >> "+out_res2
            #Running
            print("\n# "+primer,">>> Running . . .")
            os.system(write_primer)
            print("\n## Extracted hypervariable regions:")
            os.system(write_hs)
            print(count_hs)
            os.system(count_hs)
            print("\n## Dereplicated hypervariable regions:")
            os.system(write_dhs)
            print(count_dhs)
            os.system(count_dhs)
            print("\n## Filter homopolymers and degenerate bases:")
            os.system(write_cs)
            print(count_cs)
            os.system(count_cs)
            print("\n## Dereplicated filtered homopolymers and degenerate bases:")
            os.system(write_dcs)
            print(count_dcs)
            os.system(count_dcs)
            print("\n## Filtered sequence lengths:")
            os.system(write_fs)
            print(count_fs)
            os.system(count_fs)
            print("\n## Dereplicated filtered sequence lengths (cleaned datasets):")
            os.system(write_dfs)
            print(count_dfs)
            os.system(count_dfs)
            print("#-------------------------------------------------------------------------------------------------#")
        fp.close()
        # endregion

        # region Create result files: tables and visualization
        ## Create number of records dictionaries: gene_dict -> {gene1:{step1:num_rec ,step2:num_rec}, gene2:...})
        ### Mock community datasets
        out_res1 = thornfield+"/results/number_of_records/grep_mock.log.txt"
        fr1 = open(out_res1,"r")
        gene_dict = {}
        step_dict = {}
        for line in fr1:
            if line.strip()[0] == ">":
                gene = line.strip().split(">>")[1]
                gene_dict[gene] = ''
            elif line.strip()[0] == "#":
                step = line.strip().split("#")[1].replace(":","")
                step_dict[step] = ''
            else:
                num_rec = line.strip()
                step_dict[step] = num_rec
                gene_dict[gene] = step_dict.copy()
        fr1.close()
        ### Primer-extracted datasets
        out_res2 = thornfield+"/results/number_of_records/grep_primer.log.txt"
        fr2 = open(out_res2,"r")
        primer_dict = {}
        step_dict = {}
        for line in fr2:
            if line.strip()[0] == ">":
                primer = line.strip().split(">>")[1]
                primer_dict[primer] = ''
            elif line.strip()[0] == "#":
                step = line.strip().split("#")[1].replace(":","")
                step_dict[step] = ''
            else:
                num_rec = line.strip()
                step_dict[step] = num_rec
                primer_dict[primer] = step_dict.copy()
        fr2.close()

        ## Write number of records report table
        f = open(lowood,"r")
        header = f.readline()
        final_res = thornfield+"/results/number_of_records/num_records_REPORT.txt"
        fw = open(final_res,"w")
        fw.write("primer"+"\t"+"Downloaded NCBI data"+"\t"+"Dereplicated NCBI data"+"\t"+"Extracted hyperregion"+"\t"+"Dereplicated hyperregion"+"\t"+"Filtered hmplm and degen "+"\t"+"Dereplicated filtered hmplm and degen"+"\t"+"Filtered length"+"\t"+"Dereplicated filtered length (cleaned datasets)"+"\n")
        for line in f:
            #Variables
            primer = line.split("\t")[0]
            gene = line.split("\t")[1]
            forward = line.split("\t")[2]
            reverse = line.split("\t")[3].strip()
            #Write file
            fw.write(primer+"\t"+gene_dict[gene]["NCBI_data"]+"\t"+gene_dict[gene]["Dereplicated_NCBI_data"]+"\t"+primer_dict[primer]["Extracted_hyperregions"]+"\t"+primer_dict[primer]["Dereplicated_hyperregions"]+"\t"+primer_dict[primer]["Filtered_hmplm_degen"]+"\t"+primer_dict[primer]["Dereplicated_filtered_hmplm_degen"]+"\t"+primer_dict[primer]["Filtered_lengths"]+"\t"+primer_dict[primer]["Dereplicated_filtered_lengths"]+"\n")
        f.close()
        fw.close()

        ## Create number of records bar plot
        df = pd.read_csv(final_res, sep='\t')
        out_bar = thornfield+"/results/number_of_records/num_records_BARPLOT.png"
        major_col = df.columns[0]
        df_long = df.melt(id_vars=major_col, var_name='Minor_Category', value_name='Value')
        plt.figure(figsize=(15, 10))
        sns.set_theme(style="darkgrid")
        sns.barplot(data=df_long, x=major_col, y='Value', hue='Minor_Category', edgecolor="black", linewidth=1, palette="Set3")
        plt.title('Number of records through processing steps', fontsize=20)
        plt.xlabel('Primers', fontsize = 18)
        plt.ylabel('Number of records', fontsize = 18)
        plt.xticks(fontsize=16, rotation=90)
        plt.yticks(fontsize=16)
        plt.legend(title='Preprocessing steps', title_fontsize=16, fontsize=14, facecolor='#eaeaf2', edgecolor='#eaeaf2', loc='center left', bbox_to_anchor=(1, 0.5))
        plt.savefig(out_bar, dpi=300, bbox_inches='tight')
        print("#>>> Saved bar plot as:", out_bar)
        # endregion
        # endregion
                                        
    # In case of ERRORS
    ## Occur if number of input files/parameters are not matched with number of primer labels, or both global min max and filter length by taxa files are provided
    else:
        if len(minlen)!=len(maxlen) and len(maxlen)!=len(prilab) or len(taxlen)==len(prilab):
            print("ERROR: number of customized parameters are not matched with the primer labels provided.")
        elif minlen!=0 and maxlen!=0 and taxlen!="":
            print("ERROR: you can choose only one filtering method, either 1.) provide minimum and maximum length for filtering length of all taxa or 2.) provide length filering criteria file for specific taxa.")
        elif minlen==0 and maxlen==0 and taxlen=="":
            print("ERROR: please provide customized sequence length filtering criteria.")

# endregion
########################################################################################################################
# region Function #4: evaluate primer datasets

def evaluateprimer(args):

    # region Assign parsed arguments from command flags to python variables
    gateshead = args.i_mock
    lowood = args.i_primer
    thornfield = args.i_dir
    target = args.i_target_taxa
    themoor = args.o_dir
    level = args.p_tax_level
    p = args.p_vary_confidence
    # endregion

    # region Retrieve "taxid" and "my" variables from the previous command (downloadmock)
    farg1 = open("inp_args1.log.txt","r")
    farg1.readline()
    for line in farg1:
        taxid = line.split("\t")[0]
        my = line.split("\t")[1].strip()
    farg1.close()
    # endregion

#=======================================================================================================================
    # region Specificity of primers

    print("#========================================== Specificity ==================================================#")

    # region Create new directories
    newdir1 = "mkdir "+themoor
    newdir_res = "mkdir "+themoor+"/results"
    newdir2 = "mkdir "+themoor+"/specificity"
    newdir3 = "mkdir "+themoor+"/specificity/temp_files"
    newdir4 = "mkdir "+themoor+"/specificity/mock_community"
    newdir5 = "mkdir "+themoor+"/specificity/primer_extracted"
    newdir6 = "mkdir "+themoor+"/specificity/results"
    os.system(newdir1)
    os.system(newdir_res)
    os.system(newdir2)
    os.system(newdir3)
    os.system(newdir4)
    os.system(newdir5)
    os.system(newdir6)
    # endregion

    # region Get richness grouped by taxonomic levels
    ## Mock community datasets
    print("\n#### Specificity: Get richness grouped by taxonomic levels - mock community datasets ####")
    f = open(gateshead,"r")
    header = f.readline()
    for line in f:
        # region Variables
        gene = line.split("\t")[0]
        query = line.split("\t")[1].strip()
        # endregion

        # region Files
        in_tax = thornfield+"/export/"+"txid"+taxid+"_"+gene+"_"+my+"_drp_tax/taxonomy.tsv"
        out_tax = themoor+"/specificity/temp_files/"+gene+"_pre-PCR_tax.txt"
        out_tax_uniqK = themoor+"/specificity/temp_files/"+gene+"_pre-PCR_tax_RM-dupK.txt"
        out_tax_uniqP = themoor+"/specificity/temp_files/"+gene+"_pre-PCR_tax_RM-dupP.txt"
        out_tax_uniqC = themoor+"/specificity/temp_files/"+gene+"_pre-PCR_tax_RM-dupC.txt"
        out_tax_uniqO = themoor+"/specificity/temp_files/"+gene+"_pre-PCR_tax_RM-dupO.txt"
        out_tax_uniqF = themoor+"/specificity/temp_files/"+gene+"_pre-PCR_tax_RM-dupF.txt"
        out_tax_uniqG = themoor+"/specificity/temp_files/"+gene+"_pre-PCR_tax_RM-dupG.txt"
        out_tax_uniqS = themoor+"/specificity/temp_files/"+gene+"_pre-PCR_tax_RM-dupS.txt"
        out_tax_resK = themoor+"/specificity/mock_community/"+gene+"_pre-PCR_richness_Kingdom.txt"
        out_tax_resP = themoor+"/specificity/mock_community/"+gene+"_pre-PCR_richness_Phylum.txt"
        out_tax_resC = themoor+"/specificity/mock_community/"+gene+"_pre-PCR_richness_Class.txt"
        out_tax_resO = themoor+"/specificity/mock_community/"+gene+"_pre-PCR_richness_Order.txt"
        out_tax_resF = themoor+"/specificity/mock_community/"+gene+"_pre-PCR_richness_Family.txt"
        out_tax_resG = themoor+"/specificity/mock_community/"+gene+"_pre-PCR_richness_Genus.txt"
        out_tax_resS = themoor+"/specificity/mock_community/"+gene+"_pre-PCR_richness_Species.txt"
        # endregion

        # region Running message
        print("\n# "+gene,">>> Running . . .")
        # endregion

        # region Create taxonomy dictionaries
        ## Write 2 levels of dictionaries 1) inferior = {taxonomic level:value,...} and 2) superior = {accession number:{taxonomic level:velue,...},...}
        print("\n## [1.] Create taxonomy dictionary:")
        inf_taxa_dict = {}
        sup_taxa_dict = {}
        f_tax = open(in_tax, "r")
        f_tax.readline() 
        for line in f_tax:
            accss = line.split("\t")[0]
            Kingdom = line.split("\t")[1].strip().split(";")[0].split("__")[1]
            inf_taxa_dict["kingdom"] = Kingdom #inf_taxa_dict = {kingdom:Metazoa}
            #print(Kingdom)
            Phylum = line.split("\t")[1].strip().split(";")[1].split("__")[1]
            inf_taxa_dict["phylum"] = Phylum #inf_taxa_dict = {kingdom:Metazoa, phylum: ...}
            #print(Phylum)
            Class = line.split("\t")[1].strip().split(";")[2].split("__")[1]
            inf_taxa_dict["class"] = Class #inf_taxa_dict = {kingdom:Metazoa, phylum: ..., class:...}
            #print(Class)
            Order = line.split("\t")[1].strip().split(";")[3].split("__")[1]
            inf_taxa_dict["order"] = Order
            #print(Order)
            Family = line.split("\t")[1].strip().split(";")[4].split("__")[1]
            inf_taxa_dict["family"] = Family
            #print(Family)
            Genus = line.split("\t")[1].strip().split(";")[5].split("__")[1]
            inf_taxa_dict["genus"] = Genus
            #print(Genus)
            Species = line.split("\t")[1].strip().split(";")[6].split("__")[1]
            inf_taxa_dict["species"] = Species #inf_taxa_dict = {kingdom:Metazoa, phylum: ..., class:... ,.., .species:...}
            #print(Species)
            sup_taxa_dict[accss] = inf_taxa_dict #sup_taxa_dict = {accession:{kingdom:Metazoa, phylum: ..., class:... ,.., .species:...}}
            inf_taxa_dict = {}
        f_tax.close()
        print(">>> Created dictionary: sup_taxa_dict ---> {Keys = accession[taxonomic level], values = taxonomic information}")
        print("### Number of taxonomic records written in a dictionary:", len(sup_taxa_dict))
        # endregion

        # region Write all taxonomy table >> saved as specificity/temp_files/*_pre-PCR_tax.txt
        print("\n## [2.] Write a new taxonomy file:")
        fw = open(out_tax, "w")
        fw.write("accession"+"\t"+"Kingdom"+"\t"+"Phylum"+"\t"+"Class"+"\t"+"Order"+"\t"+"Family"+"\t"+"Genus"+"\t"+"Species"+"\n")
        write_line = 0
        for accession in sup_taxa_dict.keys():
            fw.write(accession+"\t"+sup_taxa_dict[accession]["kingdom"]+"\t"+sup_taxa_dict[accession]["phylum"]+"\t"+sup_taxa_dict[accession]["class"]+"\t"+sup_taxa_dict[accession]["order"]+"\t"+sup_taxa_dict[accession]["family"]+"\t"+sup_taxa_dict[accession]["genus"]+"\t"+sup_taxa_dict[accession]["genus"]+" "+sup_taxa_dict[accession]["species"]+"\n")
            write_line += 1
        print(">>> Saved as:",out_tax)
        print("### Number of records run:",write_line)
        fw.close()
        ## Check written files
        fc = open(out_tax, "r")
        num_line=0
        for line in fc:
            num_line +=1
        print("### Number of written records:",num_line-1)
        fc.close()
        # endregion

        # region Write unique taxonomy table >> saved as specificity/temp_files/*_pre-PCR_RM-dup*.txt
        print("\n## [3.] Prepare unique taxa files:")

        ## Data
        df = pd.read_table(out_tax)

        ## Prepare unique taxa file: kingdom
        dfK = df.drop_duplicates(subset=["Kingdom"]) 
        select_col = ["Kingdom"]
        dfK_select = dfK[select_col]
        dfK_select.to_csv(out_tax_uniqK, sep='\t', index=False)
        print(">>> Kingdom - saved as:",out_tax_uniqK)
        
        ## Prepare unique taxa file: phylum
        dfP = df.drop_duplicates(subset=["Phylum"])
        select_col = ["Kingdom","Phylum"]
        dfP_select = dfP[select_col]
        dfP_select.to_csv(out_tax_uniqP, sep='\t', index=False)
        print(">>> Phylum - saved as:",out_tax_uniqP)
        
        ## Prepare unique taxa file: class
        dfC = df.drop_duplicates(subset=["Class"])
        select_col = ["Kingdom","Phylum","Class"]
        dfC_select = dfC[select_col]
        dfC_select.to_csv(out_tax_uniqC, sep='\t', index=False)
        print(">>> Class - saved as:",out_tax_uniqC)
        
        ## Prepare unique taxa file: order
        dfO = df.drop_duplicates(subset=["Order"])
        select_col = ["Kingdom","Phylum","Class","Order"]
        dfO_select = dfO[select_col]
        dfO_select.to_csv(out_tax_uniqO, sep='\t', index=False)
        print(">>> Order - saved as:",out_tax_uniqO)
        
        ## Prepare unique taxa file: family
        dfF = df.drop_duplicates(subset=["Family"])
        select_col = ["Kingdom","Phylum","Class","Order","Family"]
        dfF_select = dfF[select_col]
        dfF_select.to_csv(out_tax_uniqF, sep='\t', index=False)
        print(">>> Family - saved as:",out_tax_uniqF)
        
        ## Prepare unique taxa file: genus
        dfG = df.drop_duplicates(subset=["Genus"])
        select_col = ["Kingdom","Phylum","Class","Order","Family","Genus"]
        dfG_select = dfG[select_col]
        dfG_select.to_csv(out_tax_uniqG, sep='\t', index=False)
        print(">>> Genus - saved as:",out_tax_uniqG)
        
        ## Prepare unique taxa file: species
        dfS = df.drop_duplicates(subset=["Species"]) 
        select_col = ["Kingdom","Phylum","Class","Order","Family","Genus","Species"]
        dfS_select = dfS[select_col]
        dfS_select.to_csv(out_tax_uniqS, sep='\t', index=False)
        print(">>> Species - saved as:",out_tax_uniqS)
        # endregion

        # region Write result files: Mock community richness grouped by taxonomic levels
        ## >> saved as specificity/mock_community/*_pre-PCR_richness_*.txt
        print("\n## [4.] Write result files:")

        ## Write result file: kingdom
        fK = open(out_tax_resK,"w")
        fK.write("Kingdom"+"\t"+"num_records"+"\n")
        fK_count = open(out_tax_uniqK,"r")
        fK_count.readline()
        kingdom_inf_dict = {}
        kingdom_sup_dict = {}
        kingdom_count_dict = {}
        for line in fK_count:
            #Assign variables
            Kingdom = line.strip().split("\t")[0]
            #Create kingdom_count_dict 
            dfK_count = pd.read_table(out_tax_uniqS)
            count = (dfK_count["Kingdom"] == Kingdom).sum()
            kingdom_count_dict[Kingdom] = count
            #Write output file
            fK.write(str(Kingdom)+"\t"+str(kingdom_count_dict[Kingdom])+"\n")
        print(">>> Kingdom - saved as:",out_tax_resK)
        fK_count.close()
        fK.close()
        
        ## Write result file: phylum
        fP = open(out_tax_resP,"w")
        fP.write("Kingdom"+"\t"+"Phylum"+"\t"+"num_records"+"\n")
        fP_count = open(out_tax_uniqP,"r")
        fP_count.readline()
        phylum_inf_dict = {}
        phylum_sup_dict = {}
        phylum_count_dict = {}
        for line in fP_count:
            #Assign variables
            Kingdom = line.strip().split("\t")[0]
            Phylum = line.strip().split("\t")[1]
            #Create phylum_sup_dict 
            phylum_inf_dict["Kingdom"] = Kingdom
            phylum_sup_dict[Phylum] = phylum_inf_dict
            #Create phylum_count_dict 
            dfP_count = pd.read_table(out_tax_uniqS)
            count = (dfP_count["Phylum"] == Phylum).sum()
            phylum_count_dict[Phylum] = count
            #Write output file
            fP.write(str(phylum_sup_dict[Phylum]["Kingdom"])+"\t"+str(Phylum)+"\t"+str(phylum_count_dict[Phylum])+"\n")
        print(">>> Phylum - saved as:",out_tax_resP)
        fP_count.close()
        fP.close()
        
        ## Write result file: class
        fC = open(out_tax_resC,"w")
        fC.write("Kingdom"+"\t"+"Phylum"+"\t"+"Class"+"\t"+"num_records"+"\n")
        fC_count = open(out_tax_uniqC,"r")
        fC_count.readline()
        class_inf_dict = {}
        class_sup_dict = {}
        class_count_dict = {}
        for line in fC_count:
            #Assign variables
            Kingdom = line.strip().split("\t")[0]
            Phylum = line.strip().split("\t")[1]
            Class = line.strip().split("\t")[2]
            #Create class_sup_dict --> {Class:{Kingdom:...,Phylum:...}}
            class_inf_dict["Kingdom"] = Kingdom
            class_inf_dict["Phylum"] = Phylum
            class_sup_dict[Class] = class_inf_dict
            #Create class_count_dict --> {Class: Number of unique species}
            dfC_count = pd.read_table(out_tax_uniqS)
            count = (dfC_count["Class"] == Class).sum()
            class_count_dict[Class] = count
            #Write output file
            fC.write(str(class_sup_dict[Class]["Kingdom"])+"\t"+str(class_sup_dict[Class]["Phylum"])+"\t"+str(Class)+"\t"+str(class_count_dict[Class])+"\n")
        print(">>> Class - saved as:",out_tax_resC)
        fC_count.close()
        fC.close()

        ## Write result file: order
        fO = open(out_tax_resO,"w")
        fO.write("Kingdom"+"\t"+"Phylum"+"\t"+"Class"+"\t"+"Order"+"\t"+"num_records"+"\n")
        fO_count = open(out_tax_uniqO,"r")
        fO_count.readline()
        order_inf_dict = {}
        order_sup_dict = {}
        order_count_dict = {}
        for line in fO_count:
            #Assign variables
            Kingdom = line.strip().split("\t")[0]
            Phylum = line.strip().split("\t")[1]
            Class = line.strip().split("\t")[2]
            Order = line.strip().split("\t")[3]
            #Create order_sup_dict
            order_inf_dict["Kingdom"] = Kingdom
            order_inf_dict["Phylum"] = Phylum
            order_inf_dict["Class"] = Class
            order_sup_dict[Order] = order_inf_dict
            #Create order_count_dict
            dfO_count = pd.read_table(out_tax_uniqS)
            count = (dfO_count["Order"] == Order).sum()
            order_count_dict[Order] = count
            #Write output file
            fO.write(str(order_sup_dict[Order]["Kingdom"])+"\t"+str(order_sup_dict[Order]["Phylum"])+"\t"+str(order_sup_dict[Order]["Class"])+"\t"+str(Order)+"\t"+str(order_count_dict[Order])+"\n")
        print(">>> Order - saved as:",out_tax_resO)
        fO_count.close()
        fO.close()

        ## Write result file: family
        fF = open(out_tax_resF,"w")
        fF.write("Kingdom"+"\t"+"Phylum"+"\t"+"Class"+"\t"+"Order"+"\t"+"Family"+"\t"+"num_records"+"\n")
        fF_count = open(out_tax_uniqF,"r")
        fF_count.readline()
        family_inf_dict = {}
        family_sup_dict = {}
        family_count_dict = {}
        for line in fF_count:
            #Assign variables
            Kingdom = line.strip().split("\t")[0]
            Phylum = line.strip().split("\t")[1]
            Class = line.strip().split("\t")[2]
            Order = line.strip().split("\t")[3]
            Family = line.strip().split("\t")[4]
            #Create family_sup_dict
            family_inf_dict["Kingdom"] = Kingdom
            family_inf_dict["Phylum"] = Phylum
            family_inf_dict["Class"] = Class
            family_inf_dict["Order"] = Order
            family_sup_dict[Family] = family_inf_dict
            #Create family_count_dict
            dfF_count = pd.read_table(out_tax_uniqS)
            count = (dfF_count["Family"] == Family).sum()
            family_count_dict[Family] = count
            #Write output file
            fF.write(str(family_sup_dict[Family]["Kingdom"])+"\t"+str(family_sup_dict[Family]["Phylum"])+"\t"+str(family_sup_dict[Family]["Class"])+"\t"+str(family_sup_dict[Family]["Order"])+"\t"+str(Family)+"\t"+str(family_count_dict[Family])+"\n")
        print(">>> Family - saved as:",out_tax_resF)
        fF_count.close()
        fF.close()

        ## Write result file: genus
        fG = open(out_tax_resG,"w")
        fG.write("Kingdom"+"\t"+"Phylum"+"\t"+"Class"+"\t"+"Order"+"\t"+"Family"+"\t"+"Genus"+"\t"+"num_records"+"\n")
        fG_count = open(out_tax_uniqG,"r")
        fG_count.readline()
        genus_inf_dict = {}
        genus_sup_dict = {}
        genus_count_dict = {}
        for line in fG_count:
            #Assign variables
            Kingdom = line.strip().split("\t")[0]
            Phylum = line.strip().split("\t")[1]
            Class = line.strip().split("\t")[2]
            Order = line.strip().split("\t")[3]
            Family = line.strip().split("\t")[4]
            Genus = line.strip().split("\t")[5]
            #Create genus_sup_dict
            genus_inf_dict["Kingdom"] = Kingdom
            genus_inf_dict["Phylum"] = Phylum
            genus_inf_dict["Class"] = Class
            genus_inf_dict["Order"] = Order
            genus_inf_dict["Family"] = Family
            genus_sup_dict[Genus] = genus_inf_dict
            #Create genus_count_dict
            dfG_count = pd.read_table(out_tax_uniqS)
            count = (dfG_count["Genus"] == Genus).sum()
            genus_count_dict[Genus] = count
            #Write output file
            fG.write(str(genus_sup_dict[Genus]["Kingdom"])+"\t"+str(genus_sup_dict[Genus]["Phylum"])+"\t"+str(genus_sup_dict[Genus]["Class"])+"\t"+str(genus_sup_dict[Genus]["Order"])+"\t"+str(genus_sup_dict[Genus]["Family"])+"\t"+str(Genus)+"\t"+str(genus_count_dict[Genus])+"\n")
        print(">>> Genus - saved as:",out_tax_resG)
        fG_count.close()
        fG.close()

        ## Write result file: species
        fS = open(out_tax_resS,"w")
        fS.write("Kingdom"+"\t"+"Phylum"+"\t"+"Class"+"\t"+"Order"+"\t"+"Family"+"\t"+"Genus"+"\t"+"Species"+"\t"+"num_records"+"\n")
        fS_count = open(out_tax_uniqS,"r")
        fS_count.readline()
        species_inf_dict = {}
        species_sup_dict = {}
        species_count_dict = {}
        for line in fS_count:
            #Assign variables
            Kingdom = line.strip().split("\t")[0]
            Phylum = line.strip().split("\t")[1]
            Class = line.strip().split("\t")[2]
            Order = line.strip().split("\t")[3]
            Family = line.strip().split("\t")[4]
            Genus = line.strip().split("\t")[5]
            Species = line.strip().split("\t")[6]
            #Create species_sup_dict
            species_inf_dict["Kingdom"] = Kingdom
            species_inf_dict["Phylum"] = Phylum
            species_inf_dict["Class"] = Class
            species_inf_dict["Order"] = Order
            species_inf_dict["Family"] = Family
            species_inf_dict["Genus"] = Genus
            species_sup_dict[Species] = species_inf_dict
            #Create species_count_dict
            dfS_count = pd.read_table(out_tax)
            count = (dfS_count["Species"] == Species).sum()
            species_count_dict[Species] = count
            #Write output file
            fS.write(str(species_sup_dict[Species]["Kingdom"])+"\t"+str(species_sup_dict[Species]["Phylum"])+"\t"+str(species_sup_dict[Species]["Class"])+"\t"+str(species_sup_dict[Species]["Order"])+"\t"+str(species_sup_dict[Species]["Family"])+"\t"+str(species_sup_dict[Species]["Genus"])+"\t"+str(Species)+"\t"+str(species_count_dict[Species])+"\n")
        print(">>> Species* - saved as:",out_tax_resS)
        print("*Number of all records including intraspecific variations.")
        fS_count.close()
        fS.close()
        # endregion
    f.close()

    ## Primer-extracted datasets
    print("\n#### Specificity: Get richness grouped by taxonomic levels - primer extracted datasets ####")
    f = open(lowood,"r")
    header = f.readline()
    for line in f:
        # region Variables
        primer = line.split("\t")[0]
        gene = line.split("\t")[1]
        forward = line.split("\t")[2]
        reverse = line.split("\t")[3].strip()
        # endregion

        # region Files
        in_tax = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_tax/taxonomy.tsv"
        out_tax = themoor+"/specificity/temp_files/"+primer+"_post-PCR_tax.txt"
        out_tax_uniqK = themoor+"/specificity/temp_files/"+primer+"_post-PCR_tax_RM-dupK.txt"
        out_tax_uniqP = themoor+"/specificity/temp_files/"+primer+"_post-PCR_tax_RM-dupP.txt"
        out_tax_uniqC = themoor+"/specificity/temp_files/"+primer+"_post-PCR_tax_RM-dupC.txt"
        out_tax_uniqO = themoor+"/specificity/temp_files/"+primer+"_post-PCR_tax_RM-dupO.txt"
        out_tax_uniqF = themoor+"/specificity/temp_files/"+primer+"_post-PCR_tax_RM-dupF.txt"
        out_tax_uniqG = themoor+"/specificity/temp_files/"+primer+"_post-PCR_tax_RM-dupG.txt"
        out_tax_uniqS = themoor+"/specificity/temp_files/"+primer+"_post-PCR_tax_RM-dupS.txt"
        out_tax_resK = themoor+"/specificity/primer_extracted/"+primer+"_post-PCR_richness_Kingdom.txt"
        out_tax_resP = themoor+"/specificity/primer_extracted/"+primer+"_post-PCR_richness_Phylum.txt"
        out_tax_resC = themoor+"/specificity/primer_extracted/"+primer+"_post-PCR_richness_Class.txt"
        out_tax_resO = themoor+"/specificity/primer_extracted/"+primer+"_post-PCR_richness_Order.txt"
        out_tax_resF = themoor+"/specificity/primer_extracted/"+primer+"_post-PCR_richness_Family.txt"
        out_tax_resG = themoor+"/specificity/primer_extracted/"+primer+"_post-PCR_richness_Genus.txt"
        out_tax_resS = themoor+"/specificity/primer_extracted/"+primer+"_post-PCR_richness_Species.txt"
        # endregion

        # region Running message
        print("\n# "+primer,">>> Running . . .")
        # endregion

        # region Create taxonomy dictionaries
        ## Write 2 levels of dictionaries 1) inferior = {taxonomic level:value,...} and 2) superior = {accession number:{taxonomic level:velue,...},...}
        print("\n## [1.] Create taxonomy dictionary:")
        inf_taxa_dict = {}
        sup_taxa_dict = {}
        f_tax = open(in_tax, "r")
        f_tax.readline() #read header 
        for line in f_tax:
            accss = line.split("\t")[0]
            Kingdom = line.split("\t")[1].strip().split(";")[0].split("__")[1]
            inf_taxa_dict["kingdom"] = Kingdom #inf_taxa_dict = {kingdom:Metazoa}
            #print(Kingdom)
            Phylum = line.split("\t")[1].strip().split(";")[1].split("__")[1]
            inf_taxa_dict["phylum"] = Phylum #inf_taxa_dict = {kingdom:Metazoa, phylum: ...}
            #print(Phylum)
            Class = line.split("\t")[1].strip().split(";")[2].split("__")[1]
            inf_taxa_dict["class"] = Class #inf_taxa_dict = {kingdom:Metazoa, phylum: ..., class:...}
            #print(Class)
            Order = line.split("\t")[1].strip().split(";")[3].split("__")[1]
            inf_taxa_dict["order"] = Order
            #print(Order)
            Family = line.split("\t")[1].strip().split(";")[4].split("__")[1]
            inf_taxa_dict["family"] = Family
            #print(Family)
            Genus = line.split("\t")[1].strip().split(";")[5].split("__")[1]
            inf_taxa_dict["genus"] = Genus
            #print(Genus)
            Species = line.split("\t")[1].strip().split(";")[6].split("__")[1]
            inf_taxa_dict["species"] = Species #inf_taxa_dict = {kingdom:Metazoa, phylum: ..., class:... ,.., .species:...}
            #print(Species)
            sup_taxa_dict[accss] = inf_taxa_dict #sup_taxa_dict = {accession:{kingdom:Metazoa, phylum: ..., class:... ,.., .species:...}}
            inf_taxa_dict = {}
        f_tax.close()
        print(">>> Created dictionary: sup_taxa_dict ---> {Keys = accession[taxonomic level], values = taxonomic information}")
        print("### Number of taxonomic records written in a dictionary:", len(sup_taxa_dict))
        # endregion

        # region Write all taxonomy table >> saved as specificity/temp_files/*_post-PCR_tax.txt
        print("\n## [2.] Write a new taxonomy file:")
        fw = open(out_tax, "w")
        fw.write("accession"+"\t"+"Kingdom"+"\t"+"Phylum"+"\t"+"Class"+"\t"+"Order"+"\t"+"Family"+"\t"+"Genus"+"\t"+"Species"+"\n")
        write_line = 0
        for accession in sup_taxa_dict.keys():
            fw.write(accession+"\t"+sup_taxa_dict[accession]["kingdom"]+"\t"+sup_taxa_dict[accession]["phylum"]+"\t"+sup_taxa_dict[accession]["class"]+"\t"+sup_taxa_dict[accession]["order"]+"\t"+sup_taxa_dict[accession]["family"]+"\t"+sup_taxa_dict[accession]["genus"]+"\t"+sup_taxa_dict[accession]["genus"]+" "+sup_taxa_dict[accession]["species"]+"\n")
            write_line += 1
        print(">>> Saved as:",out_tax)
        print("### Number of records run:",write_line)
        fw.close()
        ## Check written files
        fc = open(out_tax, "r")
        num_line=0
        for line in fc:
            num_line +=1
        print("### Number of written records:",num_line-1)
        fc.close()
        # endregion

        # region Write unique taxonomy table >> saved as specificity/temp_files/*_post-PCR_tax_RM-dup*.txt
        print("\n## [3.] Prepare unique taxa files:")

        ## Data
        df = pd.read_table(out_tax)

        ## Prepare unique file: kingdom
        dfK = df.drop_duplicates(subset=["Kingdom"]) 
        select_col = ["Kingdom"]
        dfK_select = dfK[select_col]
        dfK_select.to_csv(out_tax_uniqK, sep='\t', index=False)
        print(">>> Kingdom - saved as:",out_tax_uniqK)
        
        ## Prepare unique file: phylum
        dfP = df.drop_duplicates(subset=["Phylum"])
        select_col = ["Kingdom","Phylum"]
        dfP_select = dfP[select_col]
        dfP_select.to_csv(out_tax_uniqP, sep='\t', index=False)
        print(">>> Phylum - saved as:",out_tax_uniqP)

        ## Prepare unique file: class
        dfC = df.drop_duplicates(subset=["Class"])
        select_col = ["Kingdom","Phylum","Class"]
        dfC_select = dfC[select_col]
        dfC_select.to_csv(out_tax_uniqC, sep='\t', index=False)
        print(">>> Class - saved as:",out_tax_uniqC)

        ## Prepare unique file: order
        dfO = df.drop_duplicates(subset=["Order"])
        select_col = ["Kingdom","Phylum","Class","Order"]
        dfO_select = dfO[select_col]
        dfO_select.to_csv(out_tax_uniqO, sep='\t', index=False)
        print(">>> Order - saved as:",out_tax_uniqO)

        ## Prepare unique file: family
        dfF = df.drop_duplicates(subset=["Family"])
        select_col = ["Kingdom","Phylum","Class","Order","Family"]
        dfF_select = dfF[select_col]
        dfF_select.to_csv(out_tax_uniqF, sep='\t', index=False)
        print(">>> Family - saved as:",out_tax_uniqF)   

        ## Prepare unique file: genus
        dfG = df.drop_duplicates(subset=["Genus"])
        select_col = ["Kingdom","Phylum","Class","Order","Family","Genus"]
        dfG_select = dfG[select_col]
        dfG_select.to_csv(out_tax_uniqG, sep='\t', index=False)
        print(">>> Genus - saved as:",out_tax_uniqG)
        
        ## Prepare unique file: species
        dfS = df.drop_duplicates(subset=["Species"]) 
        select_col = ["Kingdom","Phylum","Class","Order","Family","Genus","Species"]
        dfS_select = dfS[select_col]
        dfS_select.to_csv(out_tax_uniqS, sep='\t', index=False)
        print(">>> Species - saved as:",out_tax_uniqS)
        # endregion

        # region Write result files: Primer-extracted richness grouped by taxonomic levels
        ## saved as specificity/primer_extracted/*_post-PCR_richness_*.txt
        print("\n## [4.] Write result files:")

        ## Write result file: kingdom
        fK = open(out_tax_resK,"w")
        fK.write("Kingdom"+"\t"+"num_records"+"\n")
        fK_count = open(out_tax_uniqK,"r")
        fK_count.readline()
        kingdom_inf_dict = {}
        kingdom_sup_dict = {}
        kingdom_count_dict = {}
        for line in fK_count:
            #Assign variables
            Kingdom = line.strip().split("\t")[0]
            #Create kingdom_count_dict 
            dfK_count = pd.read_table(out_tax_uniqS)
            count = (dfK_count["Kingdom"] == Kingdom).sum()
            kingdom_count_dict[Kingdom] = count
            #Write output file
            fK.write(str(Kingdom)+"\t"+str(kingdom_count_dict[Kingdom])+"\n")
        print(">>> Kingdom - saved as:",out_tax_resK)
        fK_count.close()
        fK.close()
        
        ## Write result file: phylum
        fP = open(out_tax_resP,"w")
        fP.write("Kingdom"+"\t"+"Phylum"+"\t"+"num_records"+"\n")
        fP_count = open(out_tax_uniqP,"r")
        fP_count.readline()
        phylum_inf_dict = {}
        phylum_sup_dict = {}
        phylum_count_dict = {}
        for line in fP_count:
            #Assign variables
            Kingdom = line.strip().split("\t")[0]
            Phylum = line.strip().split("\t")[1]
            #Create phylum_sup_dict 
            phylum_inf_dict["Kingdom"] = Kingdom
            phylum_sup_dict[Phylum] = phylum_inf_dict
            #Create phylum_count_dict 
            dfP_count = pd.read_table(out_tax_uniqS)
            count = (dfP_count["Phylum"] == Phylum).sum()
            phylum_count_dict[Phylum] = count
            #Write output file
            fP.write(str(phylum_sup_dict[Phylum]["Kingdom"])+"\t"+str(Phylum)+"\t"+str(phylum_count_dict[Phylum])+"\n")
        print(">>> Phylum - saved as:",out_tax_resP)
        fP_count.close()
        fP.close()
        
        ## Write result file: class
        fC = open(out_tax_resC,"w")
        fC.write("Kingdom"+"\t"+"Phylum"+"\t"+"Class"+"\t"+"num_records"+"\n")
        fC_count = open(out_tax_uniqC,"r")
        fC_count.readline()
        class_inf_dict = {}
        class_sup_dict = {}
        class_count_dict = {}
        for line in fC_count:
            #Assign variables
            Kingdom = line.strip().split("\t")[0]
            Phylum = line.strip().split("\t")[1]
            Class = line.strip().split("\t")[2]
            #Create class_sup_dict --> {Class:{Kingdom:...,Phylum:...}}
            class_inf_dict["Kingdom"] = Kingdom
            class_inf_dict["Phylum"] = Phylum
            class_sup_dict[Class] = class_inf_dict
            #Create class_count_dict --> {Class: Number of unique species}
            dfC_count = pd.read_table(out_tax_uniqS)
            count = (dfC_count["Class"] == Class).sum()
            class_count_dict[Class] = count
            #Write output file
            fC.write(str(class_sup_dict[Class]["Kingdom"])+"\t"+str(class_sup_dict[Class]["Phylum"])+"\t"+str(Class)+"\t"+str(class_count_dict[Class])+"\n")
        print(">>> Class - saved as:",out_tax_resC)
        fC_count.close()
        fC.close()

        ## Write result file: order
        fO = open(out_tax_resO,"w")
        fO.write("Kingdom"+"\t"+"Phylum"+"\t"+"Class"+"\t"+"Order"+"\t"+"num_records"+"\n")
        fO_count = open(out_tax_uniqO,"r")
        fO_count.readline()
        order_inf_dict = {}
        order_sup_dict = {}
        order_count_dict = {}
        for line in fO_count:
            #Assign variables
            Kingdom = line.strip().split("\t")[0]
            Phylum = line.strip().split("\t")[1]
            Class = line.strip().split("\t")[2]
            Order = line.strip().split("\t")[3]
            #Create order_sup_dict
            order_inf_dict["Kingdom"] = Kingdom
            order_inf_dict["Phylum"] = Phylum
            order_inf_dict["Class"] = Class
            order_sup_dict[Order] = order_inf_dict
            #Create order_count_dict
            dfO_count = pd.read_table(out_tax_uniqS)
            count = (dfO_count["Order"] == Order).sum()
            order_count_dict[Order] = count
            #Write output file
            fO.write(str(order_sup_dict[Order]["Kingdom"])+"\t"+str(order_sup_dict[Order]["Phylum"])+"\t"+str(order_sup_dict[Order]["Class"])+"\t"+str(Order)+"\t"+str(order_count_dict[Order])+"\n")
        print(">>> Order - saved as:",out_tax_resO)
        fO_count.close()
        fO.close()

        ## Write result file: family
        fF = open(out_tax_resF,"w")
        fF.write("Kingdom"+"\t"+"Phylum"+"\t"+"Class"+"\t"+"Order"+"\t"+"Family"+"\t"+"num_records"+"\n")
        fF_count = open(out_tax_uniqF,"r")
        fF_count.readline()
        family_inf_dict = {}
        family_sup_dict = {}
        family_count_dict = {}
        for line in fF_count:
            #Assign variables
            Kingdom = line.strip().split("\t")[0]
            Phylum = line.strip().split("\t")[1]
            Class = line.strip().split("\t")[2]
            Order = line.strip().split("\t")[3]
            Family = line.strip().split("\t")[4]
            #Create family_sup_dict
            family_inf_dict["Kingdom"] = Kingdom
            family_inf_dict["Phylum"] = Phylum
            family_inf_dict["Class"] = Class
            family_inf_dict["Order"] = Order
            family_sup_dict[Family] = family_inf_dict
            #Create family_count_dict
            dfF_count = pd.read_table(out_tax_uniqS)
            count = (dfF_count["Family"] == Family).sum()
            family_count_dict[Family] = count
            #Write output file
            fF.write(str(family_sup_dict[Family]["Kingdom"])+"\t"+str(family_sup_dict[Family]["Phylum"])+"\t"+str(family_sup_dict[Family]["Class"])+"\t"+str(family_sup_dict[Family]["Order"])+"\t"+str(Family)+"\t"+str(family_count_dict[Family])+"\n")
        print(">>> Family - saved as:",out_tax_resF)
        fF_count.close()
        fF.close()

        ## Write result file: genus
        fG = open(out_tax_resG,"w")
        fG.write("Kingdom"+"\t"+"Phylum"+"\t"+"Class"+"\t"+"Order"+"\t"+"Family"+"\t"+"Genus"+"\t"+"num_records"+"\n")
        fG_count = open(out_tax_uniqG,"r")
        fG_count.readline()
        genus_inf_dict = {}
        genus_sup_dict = {}
        genus_count_dict = {}
        for line in fG_count:
            #Assign variables
            Kingdom = line.strip().split("\t")[0]
            Phylum = line.strip().split("\t")[1]
            Class = line.strip().split("\t")[2]
            Order = line.strip().split("\t")[3]
            Family = line.strip().split("\t")[4]
            Genus = line.strip().split("\t")[5]
            #Create genus_sup_dict
            genus_inf_dict["Kingdom"] = Kingdom
            genus_inf_dict["Phylum"] = Phylum
            genus_inf_dict["Class"] = Class
            genus_inf_dict["Order"] = Order
            genus_inf_dict["Family"] = Family
            genus_sup_dict[Genus] = genus_inf_dict
            #Create genus_count_dict
            dfG_count = pd.read_table(out_tax_uniqS)
            count = (dfG_count["Genus"] == Genus).sum()
            genus_count_dict[Genus] = count
            #Write output file
            fG.write(str(genus_sup_dict[Genus]["Kingdom"])+"\t"+str(genus_sup_dict[Genus]["Phylum"])+"\t"+str(genus_sup_dict[Genus]["Class"])+"\t"+str(genus_sup_dict[Genus]["Order"])+"\t"+str(genus_sup_dict[Genus]["Family"])+"\t"+str(Genus)+"\t"+str(genus_count_dict[Genus])+"\n")
        print(">>> Genus - saved as:",out_tax_resG)
        fG_count.close()
        fG.close()

        ## Write result file: species
        fS = open(out_tax_resS,"w")
        fS.write("Kingdom"+"\t"+"Phylum"+"\t"+"Class"+"\t"+"Order"+"\t"+"Family"+"\t"+"Genus"+"\t"+"Species"+"\t"+"num_records"+"\n")
        fS_count = open(out_tax_uniqS,"r")
        fS_count.readline()
        species_inf_dict = {}
        species_sup_dict = {}
        species_count_dict = {}
        for line in fS_count:
            #Assign variables
            Kingdom = line.strip().split("\t")[0]
            Phylum = line.strip().split("\t")[1]
            Class = line.strip().split("\t")[2]
            Order = line.strip().split("\t")[3]
            Family = line.strip().split("\t")[4]
            Genus = line.strip().split("\t")[5]
            Species = line.strip().split("\t")[6]
            #Create species_sup_dict
            species_inf_dict["Kingdom"] = Kingdom
            species_inf_dict["Phylum"] = Phylum
            species_inf_dict["Class"] = Class
            species_inf_dict["Order"] = Order
            species_inf_dict["Family"] = Family
            species_inf_dict["Genus"] = Genus
            species_sup_dict[Species] = species_inf_dict
            #Create species_count_dict
            dfS_count = pd.read_table(out_tax)
            count = (dfS_count["Species"] == Species).sum()
            species_count_dict[Species] = count
            #Write output file
            fS.write(str(species_sup_dict[Species]["Kingdom"])+"\t"+str(species_sup_dict[Species]["Phylum"])+"\t"+str(species_sup_dict[Species]["Class"])+"\t"+str(species_sup_dict[Species]["Order"])+"\t"+str(species_sup_dict[Species]["Family"])+"\t"+str(species_sup_dict[Species]["Genus"])+"\t"+str(Species)+"\t"+str(species_count_dict[Species])+"\n")
        print(">>> Species* - saved as:",out_tax_resS)
        print("*Number of all records including intraspecific variations.")
        fS_count.close()
        fS.close()
        # endregion
    f.close()
    # endregion

    # region Summarize richness grouped by taxonomic levels among all datasets
    print("\n#### Specificity: Summarize richness grouped by taxonomic levels among all datasets ####")

    ## Mock community datasets
    # region Create dictionaries containing richness data from the previous step
    all_taxa_dict = {}  #{species1 : {Kingdom:taxa,...,Species:taxa},...,speciesn:{...}}, keep taxonomic classification of every species from every primer
    gene_dict = {} #{gene1: {species:{Kingdom:taxa,..., num_rec:1}, gene2:{...}}
    sup_taxa_dict = {} #{level:{Kingdom:taxa,...,num_rec:1}} --> use as value for primer_dict
    inf_taxa_dict = {} #{Kingdom:taxa,...,num_rec:1} --> use as value for species_dict
    f = open(gateshead,"r")
    header = f.readline()
    for line in f:
        sup_taxa_dict = {}
        #Variables
        gene = line.split("\t")[0]
        query = line.split("\t")[1].strip()
        # Files
        in_tax = themoor+"/specificity/mock_community/"+gene+"_pre-PCR_richness_"+level+".txt"
        # Running
        f_tax = open(in_tax,"r")
        f_tax.readline()
        for line in f_tax:
            inf_taxa_dict = {}
            if level == "Kingdom":
                Kingdom = line.split("\t")[0]
                num_records = line.split("\t")[1].strip()
                inf_taxa_dict["Kingdom"] = Kingdom
                inf_taxa_dict["num_records"] = num_records
                all_taxa_dict[Kingdom] = inf_taxa_dict #add to dict --> {all species(uniq) : {Kingdom:taxa,...,Species:taxa}}
                sup_taxa_dict[Kingdom] = inf_taxa_dict #add to dict --> {species:{Kingdom:taxa,...,num_rec:1}}

            elif level == "Phylum":
                Kingdom = line.split("\t")[0]
                Phylum = line.split("\t")[1]
                num_records = line.split("\t")[2].strip()
                inf_taxa_dict["Kingdom"] = Kingdom
                inf_taxa_dict["Phylum"] = Phylum
                inf_taxa_dict["num_records"] = num_records
                all_taxa_dict[Phylum] = inf_taxa_dict 
                sup_taxa_dict[Phylum] = inf_taxa_dict         

            elif level == "Class":
                Kingdom = line.split("\t")[0]
                Phylum = line.split("\t")[1]
                Class = line.split("\t")[2]
                num_records = line.split("\t")[3].strip()
                inf_taxa_dict["Kingdom"] = Kingdom
                inf_taxa_dict["Phylum"] = Phylum
                inf_taxa_dict["Class"] = Class
                inf_taxa_dict["num_records"] = num_records
                all_taxa_dict[Class] = inf_taxa_dict 
                sup_taxa_dict[Class] = inf_taxa_dict 

            elif level == "Order":
                Kingdom = line.split("\t")[0]
                Phylum = line.split("\t")[1]
                Class = line.split("\t")[2]
                Order = line.split("\t")[3]
                num_records = line.split("\t")[4].strip()
                inf_taxa_dict["Kingdom"] = Kingdom
                inf_taxa_dict["Phylum"] = Phylum
                inf_taxa_dict["Class"] = Class
                inf_taxa_dict["Order"] = Order
                inf_taxa_dict["num_records"] = num_records
                all_taxa_dict[Order] = inf_taxa_dict 
                sup_taxa_dict[Order] = inf_taxa_dict 

            elif level == "Family":
                Kingdom = line.split("\t")[0]
                Phylum = line.split("\t")[1]
                Class = line.split("\t")[2]
                Order = line.split("\t")[3]
                Family = line.split("\t")[4]
                num_records = line.split("\t")[5].strip()
                inf_taxa_dict["Kingdom"] = Kingdom
                inf_taxa_dict["Phylum"] = Phylum
                inf_taxa_dict["Class"] = Class
                inf_taxa_dict["Order"] = Order
                inf_taxa_dict["Family"] = Family
                inf_taxa_dict["num_records"] = num_records
                all_taxa_dict[Family] = inf_taxa_dict 
                sup_taxa_dict[Family] = inf_taxa_dict 

            elif level == "Genus":
                Kingdom = line.split("\t")[0]
                Phylum = line.split("\t")[1]
                Class = line.split("\t")[2]
                Order = line.split("\t")[3]
                Family = line.split("\t")[4]
                Genus = line.split("\t")[5]
                num_records = line.split("\t")[6].strip()
                inf_taxa_dict["Kingdom"] = Kingdom
                inf_taxa_dict["Phylum"] = Phylum
                inf_taxa_dict["Class"] = Class
                inf_taxa_dict["Order"] = Order
                inf_taxa_dict["Family"] = Family
                inf_taxa_dict["Genus"] = Genus
                inf_taxa_dict["num_records"] = num_records
                all_taxa_dict[Genus] = inf_taxa_dict 
                sup_taxa_dict[Genus] = inf_taxa_dict 

            elif level == "Species":
                Kingdom = line.split("\t")[0]
                Phylum = line.split("\t")[1]
                Class = line.split("\t")[2]
                Order = line.split("\t")[3]
                Family = line.split("\t")[4]
                Genus = line.split("\t")[5]
                Species = line.split("\t")[6]
                num_records = line.split("\t")[7].strip()
                inf_taxa_dict["Kingdom"] = Kingdom
                inf_taxa_dict["Phylum"] = Phylum
                inf_taxa_dict["Class"] = Class
                inf_taxa_dict["Order"] = Order
                inf_taxa_dict["Family"] = Family
                inf_taxa_dict["Genus"] = Genus
                inf_taxa_dict["Species"] = Species
                inf_taxa_dict["num_records"] = num_records
                all_taxa_dict[Species] = inf_taxa_dict 
                sup_taxa_dict[Species] = inf_taxa_dict             
        gene_dict[gene] = sup_taxa_dict #add to dict --> {gene1: {species:{Kingdom:taxa,..., num_rec:1}, primer2:{...}}
        f_tax.close()
    f.close()
    # endregion

    # region Create dictionaries to reformat result files & Write results
    ## Create final_gene_dict: key = all species, value = {gene1: num_rec, gene2:num_rec,...}
    final_gene_dict = {} #{species1:{gene1:num_rec1, gene2:num_rec2,...} ,species2:{},...}
    numrec_gene_dict = {} #{gene1:num_rec1, gene2:num_rec2, ..., genen:num_rec3}
    for taxa in all_taxa_dict.keys():
        for gene in gene_dict.keys():
            if taxa in gene_dict[gene]:
                numrec_gene_dict[gene] = gene_dict[gene][taxa]["num_records"] 
            else:
                numrec_gene_dict[gene] = "0"
        final_gene_dict[taxa] = numrec_gene_dict #{species1:{primer1:num_rec1, primer2:num_rec2,...} ,species2:{},...}
        numrec_gene_dict = {}

    ## Create header string 
    header_list1 = []
    header_list2 = []
    for gene in gene_dict.keys():
        header_list1.append('"'+gene+'"'+'+"\\t"+')
    for i in range(len(header_list1)):
        if i < len(header_list1)-1:
            header_list2.append(header_list1[i])
        elif i == len(header_list1)-1:
            header_list2.append(header_list1[i].replace('"\\t"+','"\\n"'))
    header_string = "".join(header_list2)
    
    ## Create values string
    val_list1 = []
    val_list2 = []
    for gene in gene_dict.keys():
        val_list1.append('final_gene_dict[taxa]["'+gene+'"]+"\\t"+')
    for i in range(len(val_list1)):
        if i < len(val_list1)-1:
            val_list2.append(val_list1[i])
        elif i == len(val_list1)-1:
            val_list2.append(val_list1[i].replace('"\\t"+','"\\n"'))
    val_string = "".join(val_list2)

    ## Write results
    out_tax = themoor+"/specificity/results/mock-community_richness_"+level+".txt"
    fw = open(out_tax,"w")
    if level == "Kingdom":
        fw.write("Kingdom"+"\t"+eval(header_string))
        for taxa in final_gene_dict.keys():
            fw.write(all_taxa_dict[taxa]["Kingdom"]+"\t"+eval(val_string))
            
    elif level == "Phylum":
        fw.write("Kingdom"+"\t"+"Phylum"+"\t"+eval(header_string))
        for taxa in final_gene_dict.keys():
            fw.write(all_taxa_dict[taxa]["Kingdom"]+"\t"+all_taxa_dict[taxa]["Phylum"]+"\t"+eval(val_string))

    elif level == "Class":
        fw.write("Kingdom"+"\t"+"Phylum"+"\t"+"Class"+"\t"+eval(header_string))
        for taxa in final_gene_dict.keys():
            fw.write(all_taxa_dict[taxa]["Kingdom"]+"\t"+all_taxa_dict[taxa]["Phylum"]+"\t"+all_taxa_dict[taxa]["Class"]+"\t"+eval(val_string))

    elif level == "Order":
        fw.write("Kingdom"+"\t"+"Phylum"+"\t"+"Class"+"\t"+"Order"+"\t"+eval(header_string))
        for taxa in final_gene_dict.keys():
            fw.write(all_taxa_dict[taxa]["Kingdom"]+"\t"+all_taxa_dict[taxa]["Phylum"]+"\t"+all_taxa_dict[taxa]["Class"]+"\t"+all_taxa_dict[taxa]["Order"]+"\t"+eval(val_string))       

    elif level == "Family":
        fw.write("Kingdom"+"\t"+"Phylum"+"\t"+"Class"+"\t"+"Order"+"\t"+"Family"+"\t"+eval(header_string))
        for taxa in final_gene_dict.keys():
            fw.write(all_taxa_dict[taxa]["Kingdom"]+"\t"+all_taxa_dict[taxa]["Phylum"]+"\t"+all_taxa_dict[taxa]["Class"]+"\t"+all_taxa_dict[taxa]["Order"]+"\t"+all_taxa_dict[taxa]["Family"]+"\t"+eval(val_string))               

    elif level == "Genus":
        fw.write("Kingdom"+"\t"+"Phylum"+"\t"+"Class"+"\t"+"Order"+"\t"+"Family"+"\t"+"Genus"+"\t"+eval(header_string))
        for taxa in final_gene_dict.keys():
            fw.write(all_taxa_dict[taxa]["Kingdom"]+"\t"+all_taxa_dict[taxa]["Phylum"]+"\t"+all_taxa_dict[taxa]["Class"]+"\t"+all_taxa_dict[taxa]["Order"]+"\t"+all_taxa_dict[taxa]["Family"]+"\t"+all_taxa_dict[taxa]["Genus"]+"\t"+eval(val_string))                       

    elif level == "Species":
        fw.write("Kingdom"+"\t"+"Phylum"+"\t"+"Class"+"\t"+"Order"+"\t"+"Family"+"\t"+"Genus"+"\t"+"Species"+"\t"+eval(header_string))
        for taxa in final_gene_dict.keys():
            fw.write(all_taxa_dict[taxa]["Kingdom"]+"\t"+all_taxa_dict[taxa]["Phylum"]+"\t"+all_taxa_dict[taxa]["Class"]+"\t"+all_taxa_dict[taxa]["Order"]+"\t"+all_taxa_dict[taxa]["Family"]+"\t"+all_taxa_dict[taxa]["Genus"]+"\t"+all_taxa_dict[taxa]["Species"]+"\t"+eval(val_string))                               

    fw.close()
    print(">>> Saved results as:",out_tax)
    # endregion

    ## Primer-extracted datasets
    # region Create dictionaries containing richness data from the previous step
    all_taxa_dict = {}  #{species1 : {Kingdom:taxa,...,Species:taxa},...,speciesn:{...}}, keep taxonomic classification of every species from every primer
    primer_dict = {} #{primer1: {species:{Kingdom:taxa,..., num_rec:1}, primer2:{...}}
    sup_taxa_dict = {} #{level:{Kingdom:taxa,...,num_rec:1}} --> use as value for primer_dict
    inf_taxa_dict = {} #{Kingdom:taxa,...,num_rec:1} --> use as value for species_dict
    f = open(lowood,"r")
    header = f.readline()
    for line in f:
        sup_taxa_dict = {}
        # Variables
        primer = line.split("\t")[0]
        gene = line.split("\t")[1]
        forward = line.split("\t")[2]
        reverse = line.split("\t")[3].strip()
        # Files
        in_tax = themoor+"/specificity/primer_extracted/"+primer+"_post-PCR_richness_"+level+".txt"
        # Running
        f_tax = open(in_tax,"r")
        f_tax.readline()
        for line in f_tax:
            inf_taxa_dict = {}
            
            if level == "Kingdom":
                Kingdom = line.split("\t")[0]
                num_records = line.split("\t")[1].strip()
                inf_taxa_dict["Kingdom"] = Kingdom
                inf_taxa_dict["num_records"] = num_records
                all_taxa_dict[Kingdom] = inf_taxa_dict
                sup_taxa_dict[Kingdom] = inf_taxa_dict  

            elif level == "Phylum":
                Kingdom = line.split("\t")[0]
                Phylum = line.split("\t")[1]
                num_records = line.split("\t")[2].strip()
                inf_taxa_dict["Kingdom"] = Kingdom
                inf_taxa_dict["Phylum"] = Phylum
                inf_taxa_dict["num_records"] = num_records
                all_taxa_dict[Phylum] = inf_taxa_dict
                sup_taxa_dict[Phylum] = inf_taxa_dict     

            elif level == "Class":
                Kingdom = line.split("\t")[0]
                Phylum = line.split("\t")[1]
                Class = line.split("\t")[2]
                num_records = line.split("\t")[3].strip()
                inf_taxa_dict["Kingdom"] = Kingdom
                inf_taxa_dict["Phylum"] = Phylum
                inf_taxa_dict["Class"] = Class
                inf_taxa_dict["num_records"] = num_records
                all_taxa_dict[Class] = inf_taxa_dict 
                sup_taxa_dict[Class] = inf_taxa_dict 

            elif level == "Order":
                Kingdom = line.split("\t")[0]
                Phylum = line.split("\t")[1]
                Class = line.split("\t")[2]
                Order = line.split("\t")[3]
                num_records = line.split("\t")[4].strip()
                inf_taxa_dict["Kingdom"] = Kingdom
                inf_taxa_dict["Phylum"] = Phylum
                inf_taxa_dict["Class"] = Class
                inf_taxa_dict["Order"] = Order
                inf_taxa_dict["num_records"] = num_records
                all_taxa_dict[Order] = inf_taxa_dict 
                sup_taxa_dict[Order] = inf_taxa_dict 

            elif level == "Family":
                Kingdom = line.split("\t")[0]
                Phylum = line.split("\t")[1]
                Class = line.split("\t")[2]
                Order = line.split("\t")[3]
                Family = line.split("\t")[4]
                num_records = line.split("\t")[5].strip()
                inf_taxa_dict["Kingdom"] = Kingdom
                inf_taxa_dict["Phylum"] = Phylum
                inf_taxa_dict["Class"] = Class
                inf_taxa_dict["Order"] = Order
                inf_taxa_dict["Family"] = Family
                inf_taxa_dict["num_records"] = num_records
                all_taxa_dict[Family] = inf_taxa_dict 
                sup_taxa_dict[Family] = inf_taxa_dict 

            elif level == "Genus":
                Kingdom = line.split("\t")[0]
                Phylum = line.split("\t")[1]
                Class = line.split("\t")[2]
                Order = line.split("\t")[3]
                Family = line.split("\t")[4]
                Genus = line.split("\t")[5]
                num_records = line.split("\t")[6].strip()
                inf_taxa_dict["Kingdom"] = Kingdom
                inf_taxa_dict["Phylum"] = Phylum
                inf_taxa_dict["Class"] = Class
                inf_taxa_dict["Order"] = Order
                inf_taxa_dict["Family"] = Family
                inf_taxa_dict["Genus"] = Genus
                inf_taxa_dict["num_records"] = num_records
                all_taxa_dict[Genus] = inf_taxa_dict 
                sup_taxa_dict[Genus] = inf_taxa_dict 

            elif level == "Species":
                Kingdom = line.split("\t")[0]
                Phylum = line.split("\t")[1]
                Class = line.split("\t")[2]
                Order = line.split("\t")[3]
                Family = line.split("\t")[4]
                Genus = line.split("\t")[5]
                Species = line.split("\t")[6]
                num_records = line.split("\t")[7].strip()
                inf_taxa_dict["Kingdom"] = Kingdom
                inf_taxa_dict["Phylum"] = Phylum
                inf_taxa_dict["Class"] = Class
                inf_taxa_dict["Order"] = Order
                inf_taxa_dict["Family"] = Family
                inf_taxa_dict["Genus"] = Genus
                inf_taxa_dict["Species"] = Species
                inf_taxa_dict["num_records"] = num_records
                all_taxa_dict[Species] = inf_taxa_dict 
                sup_taxa_dict[Species] = inf_taxa_dict      
                
        primer_dict[primer] = sup_taxa_dict #add to dict --> {primer1: {species:{Kingdom:taxa,..., num_rec:1}, primer2:{...}}
        f_tax.close()
    f.close()
    # endregion

    # region Create dictionaries to reformat result files & Write results
    ## Create final_primer_dict
    final_primer_dict = {} #{species1:{primer1:num_rec1, primer2:num_rec2,...} ,species2:{},...}
    numrec_primer_dict = {} #{primer1:num_rec1, primer2:num_rec2, ..., primern:num_rec3}
    for taxa in all_taxa_dict.keys():
        for primer in primer_dict.keys():
            if taxa in primer_dict[primer]:
                numrec_primer_dict[primer] = primer_dict[primer][taxa]["num_records"] 
            else:
                numrec_primer_dict[primer] = "0"
        final_primer_dict[taxa] = numrec_primer_dict #{species1:{primer1:num_rec1, primer2:num_rec2,...} ,species2:{},...}
        numrec_primer_dict = {}

    ## Create header string 
    header_list1 = []
    header_list2 = []
    for primer in primer_dict.keys():
        header_list1.append('"'+primer+'"'+'+"\\t"+')
    for i in range(len(header_list1)):
        if i < len(header_list1)-1:
            header_list2.append(header_list1[i])
        elif i == len(header_list1)-1:
            header_list2.append(header_list1[i].replace('"\\t"+','"\\n"'))
    header_string = "".join(header_list2)

    ## Create values string
    val_list1 = []
    val_list2 = []
    for primer in primer_dict.keys():
        val_list1.append('final_primer_dict[taxa]["'+primer+'"]+"\\t"+')
    for i in range(len(val_list1)):
        if i < len(val_list1)-1:
            val_list2.append(val_list1[i])
        elif i == len(val_list1)-1:
            val_list2.append(val_list1[i].replace('"\\t"+','"\\n"'))
    val_string = "".join(val_list2)

    ## Write results
    out_tax = themoor+"/specificity/results/primer-extracted_richness_"+level+".txt"
    fw = open(out_tax,"w")
    
    if level == "Kingdom":
        fw.write("Kingdom"+"\t"+eval(header_string))
        for taxa in final_primer_dict.keys():
            fw.write(all_taxa_dict[taxa]["Kingdom"]+"\t"+eval(val_string))
            
    elif level == "Phylum":
        fw.write("Kingdom"+"\t"+"Phylum"+"\t"+eval(header_string))
        for taxa in final_primer_dict.keys():
            fw.write(all_taxa_dict[taxa]["Kingdom"]+"\t"+all_taxa_dict[taxa]["Phylum"]+"\t"+eval(val_string))

    elif level == "Class":
        fw.write("Kingdom"+"\t"+"Phylum"+"\t"+"Class"+"\t"+eval(header_string))
        for taxa in final_primer_dict.keys():
            fw.write(all_taxa_dict[taxa]["Kingdom"]+"\t"+all_taxa_dict[taxa]["Phylum"]+"\t"+all_taxa_dict[taxa]["Class"]+"\t"+eval(val_string))

    elif level == "Order":
        fw.write("Kingdom"+"\t"+"Phylum"+"\t"+"Class"+"\t"+"Order"+"\t"+eval(header_string))
        for taxa in final_primer_dict.keys():
            fw.write(all_taxa_dict[taxa]["Kingdom"]+"\t"+all_taxa_dict[taxa]["Phylum"]+"\t"+all_taxa_dict[taxa]["Class"]+"\t"+all_taxa_dict[taxa]["Order"]+"\t"+eval(val_string))       

    elif level == "Family":
        fw.write("Kingdom"+"\t"+"Phylum"+"\t"+"Class"+"\t"+"Order"+"\t"+"Family"+"\t"+eval(header_string))
        for taxa in final_primer_dict.keys():
            fw.write(all_taxa_dict[taxa]["Kingdom"]+"\t"+all_taxa_dict[taxa]["Phylum"]+"\t"+all_taxa_dict[taxa]["Class"]+"\t"+all_taxa_dict[taxa]["Order"]+"\t"+all_taxa_dict[taxa]["Family"]+"\t"+eval(val_string))               

    elif level == "Genus":
        fw.write("Kingdom"+"\t"+"Phylum"+"\t"+"Class"+"\t"+"Order"+"\t"+"Family"+"\t"+"Genus"+"\t"+eval(header_string))
        for taxa in final_primer_dict.keys():
            fw.write(all_taxa_dict[taxa]["Kingdom"]+"\t"+all_taxa_dict[taxa]["Phylum"]+"\t"+all_taxa_dict[taxa]["Class"]+"\t"+all_taxa_dict[taxa]["Order"]+"\t"+all_taxa_dict[taxa]["Family"]+"\t"+all_taxa_dict[taxa]["Genus"]+"\t"+eval(val_string))                       

    elif level == "Species":
        fw.write("Kingdom"+"\t"+"Phylum"+"\t"+"Class"+"\t"+"Order"+"\t"+"Family"+"\t"+"Genus"+"\t"+"Species"+"\t"+eval(header_string))
        for taxa in final_primer_dict.keys():
            fw.write(all_taxa_dict[taxa]["Kingdom"]+"\t"+all_taxa_dict[taxa]["Phylum"]+"\t"+all_taxa_dict[taxa]["Class"]+"\t"+all_taxa_dict[taxa]["Order"]+"\t"+all_taxa_dict[taxa]["Family"]+"\t"+all_taxa_dict[taxa]["Genus"]+"\t"+all_taxa_dict[taxa]["Species"]+"\t"+eval(val_string))                               
    fw.close()
    print(">>> Saved results as:",out_tax)
    # endregion
    # endregion

    # region Visualize specificity results: Richness bar plot before and after primer extraction
    print("\n#### Specificity: Get specificity report and bar plot ####")

    # region Prepare dictionaries
    ## Mock community datasets
    gene_dict = {}
    f = open(gateshead,"r") 
    header = f.readline()
    for line in f:
        richness_dict = {}
        #Variables
        gene = line.split("\t")[0]
        query = line.split("\t")[1].strip()
        # Files
        in_tax = themoor+"/specificity/mock_community/"+gene+"_pre-PCR_richness_"+level+".txt"
        f_tax = open(in_tax,"r")
        f_tax.readline()
        for line in f_tax:
            if level == "Kingdom":
                taxa = line.split("\t")[0]
                richness = line.split("\t")[1].strip()
                richness_dict[taxa] = richness

            elif level == "Phylum":
                taxa = line.split("\t")[1]
                richness = line.split("\t")[2].strip()
                richness_dict[taxa] = richness

            elif level == "Class":
                taxa = line.split("\t")[2]
                richness = line.split("\t")[3].strip()
                richness_dict[taxa] = richness

            elif level == "Order":
                taxa = line.split("\t")[3]
                richness = line.split("\t")[4].strip()
                richness_dict[taxa] = richness

            elif level == "Family":
                taxa = line.split("\t")[4]
                richness = line.split("\t")[5].strip()
                richness_dict[taxa] = richness

            elif level == "Genus":
                taxa = line.split("\t")[5]
                richness = line.split("\t")[6].strip()
                richness_dict[taxa] = richness

            elif level == "Species":
                taxa = line.split("\t")[6]
                richness = line.split("\t")[7].strip()
                richness_dict[taxa] = richness

        gene_dict[gene] = richness_dict
        f_tax.close()
    f.close()

    ## Primer-extracted datasets
    primer_dict = {}
    f = open(lowood,"r") 
    header = f.readline()
    for line in f:
        richness_dict = {}
        #Variables
        primer = line.split("\t")[0]
        gene = line.split("\t")[1]
        forward = line.split("\t")[2]
        reverse = line.split("\t")[3].strip()
        # Files
        in_tax = themoor+"/specificity/primer_extracted/"+primer+"_post-PCR_richness_"+level+".txt"
        f_tax = open(in_tax,"r")
        f_tax.readline()
        for line in f_tax:
            if level == "Kingdom":
                taxa = line.split("\t")[0]
                richness = line.split("\t")[1].strip()
                richness_dict[taxa] = richness

            elif level == "Phylum":
                taxa = line.split("\t")[1]
                richness = line.split("\t")[2].strip()
                richness_dict[taxa] = richness

            elif level == "Class":
                taxa = line.split("\t")[2]
                richness = line.split("\t")[3].strip()
                richness_dict[taxa] = richness

            elif level == "Order":
                taxa = line.split("\t")[3]
                richness = line.split("\t")[4].strip()
                richness_dict[taxa] = richness

            elif level == "Family":
                taxa = line.split("\t")[4]
                richness = line.split("\t")[5].strip()
                richness_dict[taxa] = richness 

            elif level == "Genus":
                taxa = line.split("\t")[5]
                richness = line.split("\t")[6].strip()
                richness_dict[taxa] = richness 

            elif level == "Species":
                taxa = line.split("\t")[6]
                richness = line.split("\t")[7].strip()
                richness_dict[taxa] = richness 

        primer_dict[gene+"_"+primer] = richness_dict
        f_tax.close()
    f.close()
    # endregion

    # region Write report files
    out_res = themoor+"/specificity/results/specificity_"+level+"_REPORT.txt"
    fw = open(out_res,"w")
    fw.write("gene"+"\t"+"primer"+"\t"+"taxa"+"\t"+"richness"+"\n")
    for gene in gene_dict.keys():
        for taxa in gene_dict[gene]:
            fw.write(gene+"\t"+"pre-extract"+"\t"+taxa+"\t"+gene_dict[gene][taxa]+"\n")
    for primer in primer_dict.keys():
        for taxa in primer_dict[primer]:
            fw.write(primer.split("_")[0]+"\t"+primer.split("_")[1]+"\t"+taxa+"\t"+primer_dict[primer][taxa]+"\n")
    fw.close()
    print(">>> Saved results as:",out_res)

    ##copy result to main result dir
    write_res_cmd1 = "cp "+out_res+" "+themoor+"/results/"+out_res.split("/")[3]
    os.system(write_res_cmd1)
    # endregion

    # region Create bar plots
    ## Output file
    out_bar = themoor+"/specificity/results/specificity_"+level+"_BARPLOT.png"
    ## Apply the ggplot style
    plt.style.use('ggplot')
    ## Load data
    df = pd.read_csv(out_res, sep='\t')
    ## LOCK THE ORDER of the X-axis bars to match the table exactly
    df['gene'] = pd.Categorical(df['gene'], categories=df['gene'].unique(), ordered=True)
    df['primer'] = pd.Categorical(df['primer'], categories=df['primer'].unique(), ordered=True)
    ## Reshape data
    plot_df = df.pivot(index=["gene", "primer"], columns="taxa", values="richness")
    ## Calculate the total richness for each taxon and sort them descending
    taxon_order = plot_df.sum().sort_values(ascending=False).index
    plot_df = plot_df[taxon_order] 
    ## Insert Gaps between Genes for grouping
    new_index = []
    last_gene = None
    for gene, primer in plot_df.index:
        if last_gene is not None and gene != last_gene:
            new_index.append((f"gap_{gene}", "")) 
        new_index.append((gene, primer))
        last_gene = gene
    plot_df = plot_df.reindex(new_index)
    ## Create the Plot
    ### width=0.6 makes the bars thinner
    ax = plot_df.plot(kind='bar', stacked=True, colormap="Set3", 
                      width=0.6, edgecolor='white', figsize=(16, 10))
    ## Format X-axis Ticks (Primers)
    primer_labels = [idx[1] for idx in plot_df.index]
    ax.set_xticklabels(primer_labels, rotation=90, fontsize=12, color='black')
    ## Add Gene Labels (Lowered)
    genes = df['gene'].unique()
    current_pos = 0
    for gene in genes:
        n_primers = len(df[df['gene'] == gene]['primer'].unique())
        center_pos = current_pos + (n_primers - 1) / 2
        # Position lowered to 0.30 to avoid overlap with primer names
        y_position = -ax.get_ylim()[1] * 0.20 
        ax.text(center_pos, y_position, str(gene), 
                ha='center', va='top', fontsize=15, color='black')
        current_pos += n_primers + 1
    ## Axis Label and Tick Customization
    ax.set_ylabel('Species richness', fontsize=18, color='black', labelpad=15)
    ax.set_xlabel('', fontsize=18, color='black')
    ax.tick_params(axis='y', labelsize=14, labelcolor='black')
    ## Grid and Background
    ax.grid(False)
    ax.grid(True, axis='y', color='white', linestyle='-', linewidth=1)
    ax.set_facecolor('#E5E5E5') 
    ## Legend and Title
    plt.legend(title='Taxa', title_fontsize='16', fontsize='14', 
               frameon=True, facecolor='white', edgecolor='grey', 
               framealpha=1, loc='center left', bbox_to_anchor=(1, 0.6))
    plt.title('Species richness before and after primer extraction', 
              fontsize=18, color='black', pad=25)
    ## Remove Borders
    for spine in ax.spines.values():
        spine.set_visible(False)
    ## Final adjustment and export plot
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.3) 
    plt.savefig(out_bar, dpi=300, bbox_inches='tight')
    print(">>> Saved results as:",out_bar)

    ##copy result to main result dir
    write_res_cmd2 = "cp "+out_bar+" "+themoor+"/results/"+out_bar.split("/")[3]
    os.system(write_res_cmd2)
    # endregion
    # endregion

    # region Get proportion of target to all species retrieved by the primer
    print("\n#### Specificity: Calculate proportion of target taxa to all taxa retrieved by primers ####")

    # region Create new directories
    newdir1 = "mkdir "+themoor+"/availability"
    newdir2 = "mkdir "+themoor+"/availability/primer_extracted"
    newdir3 = "mkdir "+themoor+"/availability/results"
    os.system(newdir1)
    os.system(newdir2)
    os.system(newdir3)
    # endregion

    # region Retrieve a list of target taxa
    prefix = target.split(".")[0]
    fx = open(target,"r")
    taxa_lev = str(fx.readline().strip())
    ext_taxa_list = []
    for line in fx:
        ext_taxa_list.append(line.strip())
    fx.close()
    print("\n# Extract taxonomic level:",taxa_lev)
    print("# Extract taxa:")
    for item in ext_taxa_list:
        print(item)
    # endregion

    # region Create dictionaries: Count target species richness & Extract variations for the availability assessment
    f = open(lowood,"r")
    header = f.readline()
    all_line_dict = {}
    target_line_dict = {}
    for line in f:
        # Variables
        primer = line.split("\t")[0]
        gene = line.split("\t")[1]
        forward = line.split("\t")[2]
        reverse = line.split("\t")[3].strip()
        # Files
        in_tax = themoor+"/specificity/primer_extracted/"+primer+"_post-PCR_richness_Species.txt"
        out_tax = themoor+"/availability/primer_extracted/"+primer+"_"+prefix+"_variations.txt"
        # Running
        f_tax = open(in_tax,"r")
        header = f_tax.readline()
        fw = open(out_tax,"w")
        fw.write(header)
        tax_dict = {}
        all_line = 0
        target_line = 0
        for line in f_tax:
            all_line += 1 #count all species richness
            Kingdom = line.strip().split("\t")[0]
            Phylum = line.strip().split("\t")[1]
            Class = line.strip().split("\t")[2]
            Order = line.strip().split("\t")[3]
            Family = line.strip().split("\t")[4]
            Genus = line.strip().split("\t")[5]
            Species = line.strip().split("\t")[6]
            tax_dict["Kingdom"] = Kingdom
            tax_dict["Phylum"] = Phylum
            tax_dict["Class"] = Class
            tax_dict["Order"] = Order
            tax_dict["Family"] = Family
            tax_dict["Genus"] = Genus
            tax_dict["Species"] = Species
            for item in ext_taxa_list:
                if tax_dict[taxa_lev] == item:
                    fw.write(line)
                    target_line += 1 #count target species richness
            tax_dict = {}
        all_line_dict[primer] = all_line
        target_line_dict[primer] = target_line
        f_tax.close()
        fw.close()
    f.close()
    # endregion

    # region Write target species richness proportion table
    out_prop = themoor+"/specificity/results/proportion_"+prefix+"_REPORT.txt"
    fw2 = open(out_prop,"w")
    print("\nprimer"+"\t"+"all_richness"+"\t"+"target_richness"+"\t"+"perc_proportion")
    fw2.write("primer"+"\t"+"all_richness"+"\t"+"target_richness"+"\t"+"perc_proportion"+"\n")
    f = open(lowood,"r")
    header = f.readline()
    for line in f:
        # Variables
        primer = line.split("\t")[0]
        gene = line.split("\t")[1]
        forward = line.split("\t")[2]
        reverse = line.split("\t")[3].strip()
        # Write results
        print(primer+"\t"+str(all_line_dict[primer])+"\t"+str(target_line_dict[primer])+"\t"+str((target_line_dict[primer]/all_line_dict[primer])*100))
        fw2.write(primer+"\t"+str(all_line_dict[primer])+"\t"+str(target_line_dict[primer])+"\t"+str((target_line_dict[primer]/all_line_dict[primer])*100)+"\n")
    fw2.close()
    f.close()
    ## Transpose table
    df = pd.read_csv(out_prop, sep='\t')
    df_transposed = df.transpose()
    df_transposed.to_csv(out_prop, sep='\t', header=False)
    print("\n>>> Saved result as:",out_prop)

    ##copy result to main result dir
    write_res_cmd3 = "cp "+out_prop+" "+themoor+"/results/specificity_"+out_prop.split("/")[3]
    os.system(write_res_cmd3)
    # endregion
    # endregion

    # endregion
#=======================================================================================================================
    # region Reference availability and intraspecific variations

    print("\n#====================== Reference availability & intraspecific variations  ==============================")

    # region Assign empty dictionaries
    all_species_dict = {}  #{species1 : {Kingdom:taxa,...,Species:taxa},...,speciesn:{...}}, keep taxonomic classification of every species from every primer
    primer_dict = {} #{primer1: {species:{Kingdom:taxa,..., num_rec:1}, primer2:{...}}
    species_dict = {} #{species:{Kingdom:taxa,...,num_rec:1}} --> use as value for primer_dict
    taxa_dict = {} #{Kingdom:taxa,...,num_rec:1} --> use as value for species_dict
    # endregion

    # region Create dictionaries containing number of variations of target species
    ## all_species_dict --> {all species(uniq) : {Kingdom:taxa,...,Species:taxa}}
    ## primer_dict --> {primer1: {species:{Kingdom:taxa,..., num_rec:1}, primer2:{...}}
    print("\n#### Summarize target species' availability & variations from all primer-extracted datasets ####")
    f = open(lowood,"r")
    header = f.readline()
    for line in f:
        species_dict = {}
        # Variables
        primer = line.split("\t")[0]
        gene = line.split("\t")[1]
        forward = line.split("\t")[2]
        reverse = line.split("\t")[3].strip()
        # Files
        in_tax = themoor+"/availability/primer_extracted/"+primer+"_"+prefix+"_variations.txt"
        # Running, creating dictionaries
        f_tax = open(in_tax,"r")
        f_tax.readline()
        for line in f_tax:
            taxa_dict = {}
            Species = line.split("\t")[6]
            Kingdom = line.split("\t")[0]
            Phylum = line.split("\t")[1]
            Class = line.split("\t")[2]
            Order = line.split("\t")[3]
            Family = line.split("\t")[4]
            Genus = line.split("\t")[5]
            num_records = line.split("\t")[7].strip()
            taxa_dict["Kingdom"] = Kingdom
            taxa_dict["Phylum"] = Phylum
            taxa_dict["Class"] = Class
            taxa_dict["Order"] = Order
            taxa_dict["Family"] = Family
            taxa_dict["Genus"] = Genus
            taxa_dict["num_records"] = num_records
            all_species_dict[Species] = taxa_dict #add to dict --> {all species(uniq) : {Kingdom:taxa,...,Species:taxa}}
            species_dict[Species] = taxa_dict #add to dict --> {species:{Kingdom:taxa,...,num_rec:1}}
        primer_dict[primer] = species_dict #add to dict --> {primer1: {species:{Kingdom:taxa,..., num_rec:1}, primer2:{...}}
        f_tax.close()
    f.close()
    # endregion

    # region Create dictionaries to reformat result files & Write results
    ## Create final_primer_dict: key = all species, value = {primer2: num_rec, primer2:num_rec,...}
    final_primer_dict = {} #{species1:{primer1:num_rec1, primer2:num_rec2,...} ,species2:{},...}
    numrec_primer_dict = {} #{primer1:num_rec1, primer2:num_rec2, ..., primern:num_rec3}
    for spp in all_species_dict.keys():
        for primer in primer_dict.keys():
            if spp in primer_dict[primer]:
                numrec_primer_dict[primer] = primer_dict[primer][spp]["num_records"] 
            else:
                numrec_primer_dict[primer] = "0"
        final_primer_dict[spp] = numrec_primer_dict #{species1:{primer1:num_rec1, primer2:num_rec2,...} ,species2:{},...}
        numrec_primer_dict = {}

    ## Create header string 
    header_list1 = []
    header_list2 = []
    for primer in primer_dict.keys():
        header_list1.append('"'+primer+'"'+'+"\\t"+')
    for i in range(len(header_list1)):
        if i < len(header_list1)-1:
            header_list2.append(header_list1[i])
        elif i == len(header_list1)-1:
            header_list2.append(header_list1[i].replace('"\\t"+','"\\n"'))
    header_string = "".join(header_list2)

    ## Create values string
    val_list1 = []
    val_list2 = []
    for primer in primer_dict.keys():
        val_list1.append('final_primer_dict[species]["'+primer+'"]+"\\t"+')
    for i in range(len(val_list1)):
        if i < len(val_list1)-1:
            val_list2.append(val_list1[i])
        elif i == len(val_list1)-1:
            val_list2.append(val_list1[i].replace('"\\t"+','"\\n"'))
    val_string = "".join(val_list2)

    ## Write results
    out_avail = themoor+"/availability/results/availability-variations_"+prefix+"_REPORT.txt"
    fw = open(out_avail,"w")
    fw.write("Kingdom"+"\t"+"Phylum"+"\t"+"Class"+"\t"+"Order"+"\t"+"Family"+"\t"+"Genus"+"\t"+"Species"+"\t"+eval(header_string))
    for species in final_primer_dict.keys():
        fw.write(all_species_dict[species]["Kingdom"]+"\t"+all_species_dict[species]["Phylum"]+"\t"+all_species_dict[species]["Class"]+"\t"+all_species_dict[species]["Order"]+"\t"+all_species_dict[species]["Family"]+"\t"+all_species_dict[species]["Genus"]+"\t"+species+"\t"+eval(val_string))
    fw.close()
    print(">>> Saved results as:",out_avail)

    ##copy result to main result dir
    write_res_cmd4 = "cp "+out_avail+" "+themoor+"/results/"+out_avail.split("/")[3]
    os.system(write_res_cmd4)

    # endregion

    # endregion
#=======================================================================================================================
    # region Sequence length

    print("\n#======================================= Sequence length ================================================")

    # region Create new directories
    newdir1 = "mkdir "+themoor+"/length"
    newdir2 = "mkdir "+themoor+"/length/results"
    newdir3 = "mkdir "+themoor+"/length/results/report"
    newdir4 = "mkdir "+themoor+"/length/results/statistics"
    newdir5 = "mkdir "+themoor+"/length/results/visualization"
    newdir6 = "mkdir "+themoor+"/length/results/statistics/"+level
    newdir7 = "mkdir "+themoor+"/length/results/visualization/"+level
    os.system(newdir1)
    os.system(newdir2)
    os.system(newdir3)
    os.system(newdir4)
    os.system(newdir5)
    os.system(newdir6)
    os.system(newdir7)
    # endregion

    # region Get sequence length report of each record in the datasets
    print("\n#### Sequence length report: all records in the datasets ####")

    ## Create primer-length dictionary and taxonomy-length dictionary
    ### plen_dict[primer][accession] = length --> {primer1:{accession1:len1, acc2:len2,...}, primer2...}
    ### ptax_dict[primer][accession]["Kingdom"] = Kingdom of taxa1 --> {primer1:{accession1:{Kingdom:taxa1, Phylum:taxa1,..., Species:taxa1}}, primer2...}
    plen_dict = {}
    ptax_dict = {}
    f = open(lowood, "r")
    header = f.readline()
    for line in f:
        # Variables
        primer = line.split("\t")[0]
        gene = line.split("\t")[1]
        forward = line.split("\t")[2]
        reverse = line.split("\t")[3].strip()

        # Files
        in_seq = thornfield + "/export/" + "txid" + taxid + "_" + primer + "_" + my + "_clean-db_seq/dna-sequences.fasta"
        in_tax = thornfield + "/export/" + "txid" + taxid + "_" + primer + "_" + my + "_clean-db_tax/taxonomy.tsv"
        out_len = themoor + "/length/results/report/" + primer + "_seqlen_REPORT.txt"

        # Running message
        print("\n# " + primer + " >>> Running . . .")

        # Create sequence length dictionary: seqlen_dict[accession] = length --> {accession1:len1, acc2:len2,...}
        seqlen_dict = {}
        records = list(SeqIO.parse(in_seq, "fasta"))
        for record in records:
            seqlen_dict[record.id] = len(record.seq)
        ## Add seqlen_dict to plen_dict
        plen_dict[primer] = seqlen_dict

        # Create taxonomy dictionary: taxa_dict[accession]["taxonomic level (e.g.,Class)"] = taxa --> {accession1:{Kingdom:taxa1, Phylum:taxa1,..., Species:taxa1}}
        taxa_dict = {}
        f_tax = open(in_tax, "r")
        f_tax.readline()  # read header
        for line in f_tax:
            # subtaxa_dict["taxonomic level"] = taxa --> {Kingdom:taxa1, Phylum:taxa1,..., Species:taxa1}
            subtaxa_dict = {}
            accession = line.split("\t")[0]
            Kingdom = line.split("\t")[1].strip().split(";")[0].split("__")[1]
            subtaxa_dict["Kingdom"] = Kingdom
            Phylum = line.split("\t")[1].strip().split(";")[1].split("__")[1]
            subtaxa_dict["Phylum"] = Phylum
            Class = line.split("\t")[1].strip().split(";")[2].split("__")[1]
            subtaxa_dict["Class"] = Class
            Order = line.split("\t")[1].strip().split(";")[3].split("__")[1]
            subtaxa_dict["Order"] = Order
            Family = line.split("\t")[1].strip().split(";")[4].split("__")[1]
            subtaxa_dict["Family"] = Family
            Genus = line.split("\t")[1].strip().split(";")[5].split("__")[1]
            subtaxa_dict["Genus"] = Genus
            Species = line.split("\t")[1].strip().split(";")[5].split("__")[1] + " " + \
                      line.split("\t")[1].strip().split(";")[6].split("__")[1]
            subtaxa_dict["Species"] = Species
            taxa_dict[accession] = subtaxa_dict
        ## Add taxa_dict to ptax_dict
        ptax_dict[primer] = taxa_dict
        f_tax.close()

        # Write sequence length results
        num_line = 0
        fw_len = open(out_len, "w")
        fw_len.write(
            "Accession" + "\t" + "Kingdom" + "\t" + "Phylum" + "\t" + "Class" + "\t" + "Order" + "\t" + "Family" + "\t" + "Genus" + "\t" + "Species" + "\t" + "Seqlen" + "\n")
        for acc in ptax_dict[primer]:
            accession = acc
            Kingdom = ptax_dict[primer][acc]["Kingdom"]
            Phylum = ptax_dict[primer][acc]["Phylum"]
            Class = ptax_dict[primer][acc]["Class"]
            Order = ptax_dict[primer][acc]["Order"]
            Family = ptax_dict[primer][acc]["Family"]
            Genus = ptax_dict[primer][acc]["Genus"]
            Species = ptax_dict[primer][acc]["Species"]
            seqlen = plen_dict[primer][acc]
            fw_len.write(
                accession + "\t" + Kingdom + "\t" + Phylum + "\t" + Class + "\t" + Order + "\t" + Family + "\t" + Genus + "\t" + Species + "\t" + str(
                    seqlen) + "\n")
            num_line += 1
        fw_len.close()
        print(">>> Saved as:", out_len)
        print("### Number of written records:", num_line)
    f.close()
    # endregion

    # region Get sequence length statistics and distribution visualization for ALL records
    print("\n#### Sequence length statistics, box plot, and density plot: all records in the datasets ####")

    # region Create sequence length dictionary
    seqlen_dict = {}
    f = open(lowood,"r")
    header = f.readline()
    for line in f:
        #Variables
        primer = line.split("\t")[0]
        gene = line.split("\t")[1]
        forward = line.split("\t")[2]
        reverse = line.split("\t")[3].strip()
        # Files
        in_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_seq/dna-sequences.fasta"
        in_tax = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_tax/taxonomy.tsv"
        # Create sequence length dictionary
        len_list = []
        records = list(SeqIO.parse(in_seq, "fasta"))
        for record in records:
            len_list.append(len(record.seq))
        seqlen_dict[primer] = len_list
    f.close()
    # endregion

    # region Calculate sequence length statistics
    out_res = themoor+"/length/results/statistics/seqlen_STAT.txt"
    fw = open(out_res,"w")
    fw.write("Primer"+"\t"+"Total records"+"\t"+"Arithmetic mean"+"\t"+"sd"+"\t"+"Min"+"\t"+"Max"+"\t"+"P1"+"\t"+"P2.5"+"\t"+"P5"+"\t"+"P25"+"\t"+"P50"+"\t"+"P75"+"\t"+"P95"+"\t"+"P97.5"+"\t"+"P99"+"\t"+"Lower fence"+"\t"+"Upper fence"+"\n")
    f = open(lowood,"r")
    header = f.readline()
    for line in f:
        #Variables
        primer = line.split("\t")[0]
        gene = line.split("\t")[1]
        forward = line.split("\t")[2]
        reverse = line.split("\t")[3].strip()
        # Files
        length_array = np.array(seqlen_dict[primer])
        num_rec = len(seqlen_dict[primer])
        arthmean = np.mean(length_array)
        sd = np.std(length_array)
        min_len = np.min(length_array)
        max_len = np.max(length_array)
        p010 = np.percentile(length_array, 1)
        p025 = np.percentile(length_array, 2.5)
        p050 = np.percentile(length_array, 5)
        p250 = np.percentile(length_array, 25)
        p500 = np.percentile(length_array, 50)
        p750 = np.percentile(length_array, 75)
        p950 = np.percentile(length_array, 95)
        p975 = np.percentile(length_array, 97.5)
        p990 = np.percentile(length_array, 99)
        low_fen = p250-1.5*(p750-p250)
        up_fen = p750+1.5*(p750-p250)
        fw.write(primer+"\t"+str(num_rec)+"\t"+str(arthmean)+"\t"+str(sd)+"\t"+str(min_len)+"\t"+str(max_len)+"\t"+str(p010)+"\t"+str(p025)+"\t"+str(p050)+"\t"+str(p250)+"\t"+str(p500)+"\t"+str(p750)+"\t"+str(p950)+"\t"+str(p975)+"\t"+str(p990)+"\t"+str(low_fen)+"\t"+str(up_fen)+"\n")
    fw.close()    
    f.close()
    # endregion

    # region Write statistics table
    df = pd.read_csv(out_res, sep='\t')
    df_transposed = df.transpose()
    df_transposed.to_csv(out_res, sep='\t', header=False)
    print("#>>> Saved statistics table as:",out_res)

    ##copy result to main result dir
    write_res_cmd5 = "cp "+out_res+" "+themoor+"/results/"+out_res.split("/")[4].replace("seqlen","length")
    os.system(write_res_cmd5)
    # endregion

    # region Create visualization results: Box plot & Density plot
    ## Data
    data = seqlen_dict
    df = pd.DataFrame([{'Group': k, 'Value': v} for k, values in data.items() for v in values])
    
    ## Boxplot
    out_box = themoor+"/length/results/visualization/seqlen_BOXPLOT.png"
    plt.figure(figsize=(15, 10))
    sns.set_theme(style="darkgrid")
    sns.boxplot(data=df, x='Group', y='Value', palette="Set3", dodge=False)
    plt.title('Box plot of sequence length across all primers',fontsize=15)
    plt.xticks(rotation=90, fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel('Primer', fontsize=14)
    plt.ylabel('Sequence length', fontsize=14)
    plt.savefig(out_box, dpi=300, bbox_inches='tight')
    print("#>>> Saved box plot as:",out_box)
    
    ## Density plot
    out_den = themoor+"/length/results/visualization/seqlen_DENSITY.png"
    plt.figure(figsize=(15, 10))
    sns.set_theme(style="darkgrid")
    ax = sns.kdeplot(data=df, x='Value', hue='Group', fill=True, common_norm=False, linewidth=1.5, alpha=0.7, palette="Set3", warn_singular=False)
    sns.move_legend(ax, "best", title='Primer', title_fontsize=14,fontsize=12,frameon=True, facecolor='white', framealpha=1, edgecolor='white')
    plt.title('Density plot of sequence length across all primers',fontsize=15)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel('Sequence length', fontsize=14)
    plt.ylabel('Density', fontsize=14)
    plt.savefig(out_den, dpi=300, bbox_inches='tight')
    print("#>>> Saved density plot as:",out_den)

    ##copy result to main result dir
    write_res_cmd6 = "cp "+out_box+" "+themoor+"/results/"+out_box.split("/")[4].replace("seqlen","length")
    os.system(write_res_cmd6)
    write_res_cmd7 = "cp "+out_den+" "+themoor+"/results/"+out_den.split("/")[4].replace("seqlen","length")
    os.system(write_res_cmd7)
    # endregion
    # endregion

    # region Get sequence length statistics and distribution visualization grouped by TAXA
    print("\n#### Sequence length statistics, box plot, and desity plot: grouped by taxaonomic level ####")

    flen_dict = {}
    for primer in plen_dict.keys():
        # Running message
        print("\n# "+primer+" >>> Running . . .")

        # Create dictionaries
        ## flen_dict[primer]["uniq taxa1"] = [len1, len2, len3, ...] --> {primer:{taxa1:[len1, len2, len3, ...]}}
        ## e.g., flen_dict[primer]["Actinopteri"] = [len1, len2, len3, ...]
        i=0
        taxlen_dict = {} #{taxa1:[len1, len2, len3, ...]}
        len_list = [] #[len1, len2, len3, ...]
        k_set = set()
        p_set = set()
        c_set = set()
        o_set = set()
        f_set = set()
        g_set = set()
        s_set = set()
        for access in plen_dict[primer].keys():
            k_set.add(ptax_dict[primer][access]["Kingdom"])
            p_set.add(ptax_dict[primer][access]["Phylum"])
            c_set.add(ptax_dict[primer][access]["Class"])
            o_set.add(ptax_dict[primer][access]["Order"])
            f_set.add(ptax_dict[primer][access]["Family"])
            g_set.add(ptax_dict[primer][access]["Genus"])
            s_set.add(ptax_dict[primer][access]["Species"])
        k_list = list(k_set)
        p_list = list(p_set)
        c_list = list(c_set)
        o_list = list(o_set)
        f_list = list(f_set)
        g_list = list(g_set)
        s_list = list(s_set)

        if level == "Kingdom":
            for uniqtax in k_list:
                for access in plen_dict[primer].keys():
                    if ptax_dict[primer][access]["Kingdom"] == uniqtax:
                        len_list.append(plen_dict[primer][access])
                        taxlen_dict[uniqtax] = len_list
                len_list = []
            flen_dict[primer] = taxlen_dict
        elif level == "Phylum":
            for uniqtax in p_list:
                for access in plen_dict[primer].keys():
                    if ptax_dict[primer][access]["Phylum"] == uniqtax:
                        len_list.append(plen_dict[primer][access])
                        taxlen_dict[uniqtax] = len_list
                len_list = []
            flen_dict[primer] = taxlen_dict
        elif level == "Class":
            for uniqtax in c_list:
                for access in plen_dict[primer].keys():
                    if ptax_dict[primer][access]["Class"] == uniqtax:
                        len_list.append(plen_dict[primer][access])
                        taxlen_dict[uniqtax] = len_list
                len_list = []
            flen_dict[primer] = taxlen_dict
        elif level == "Order":
            for uniqtax in o_list:
                for access in plen_dict[primer].keys():
                    if ptax_dict[primer][access]["Order"] == uniqtax:
                        len_list.append(plen_dict[primer][access])
                        taxlen_dict[uniqtax] = len_list
                len_list = []
            flen_dict[primer] = taxlen_dict 
        elif level == "Family":
            for uniqtax in f_list:
                for access in plen_dict[primer].keys():
                    if ptax_dict[primer][access]["Family"] == uniqtax:
                        len_list.append(plen_dict[primer][access])
                        taxlen_dict[uniqtax] = len_list
                len_list = []
            flen_dict[primer] = taxlen_dict 
        elif level == "Genus":
            for uniqtax in g_list:
                for access in plen_dict[primer].keys():
                    if ptax_dict[primer][access]["Genus"] == uniqtax:
                        len_list.append(plen_dict[primer][access])
                        taxlen_dict[uniqtax] = len_list
                len_list = []
            flen_dict[primer] = taxlen_dict 
        elif level == "Species":
            for uniqtax in s_list:
                for access in plen_dict[primer].keys():
                    if ptax_dict[primer][access]["Species"] == uniqtax:
                        len_list.append(plen_dict[primer][access])
                        taxlen_dict[uniqtax] = len_list
                len_list = []
            flen_dict[primer] = taxlen_dict 
            
        # Write statistics tables
        out_stat = themoor+"/length/results/statistics/"+level+"/"+primer+"_seqlen_"+level+"_STAT.txt"
        fw_stat = open(out_stat,"w")
        fw_stat.write("Taxa"+"\t"+"Total records"+"\t"+"Seqlen arithmetic mean"+"\t"+"sd"+"\t"+"Min"+"\t"+"Max"+"\t"+"P1"+"\t"+"P2.5"+"\t"+"P5"+"\t"+"P25"+"\t"+"P50"+"\t"+"P75"+"\t"+"P95"+"\t"+"P97.5"+"\t"+"P99"+"\t"+"Lower fence"+"\t"+"Upper fence"+"\n")
        for taxa in flen_dict[primer].keys():
            length_list = flen_dict[primer][taxa]
            length_array = np.array(flen_dict[primer][taxa])
            taxonomy = taxa 
            num_rec = len(length_list)
            arthmean = np.mean(length_array)
            sd = np.std(length_array)
            min_len = np.min(length_array)
            max_len = np.max(length_array)
            p010 = np.percentile(length_array, 1)
            p025 = np.percentile(length_array, 2.5)
            p050 = np.percentile(length_array, 5)
            p250 = np.percentile(length_array, 25)
            p500 = np.percentile(length_array, 50)
            p750 = np.percentile(length_array, 75)
            p950 = np.percentile(length_array, 95)
            p975 = np.percentile(length_array, 97.5)
            p990 = np.percentile(length_array, 99)
            low_fen = p250-1.5*(p750-p250)
            up_fen = p750+1.5*(p750-p250)
            fw_stat.write(taxonomy+"\t"+str(num_rec)+"\t"+str(arthmean)+"\t"+str(sd)+"\t"+str(min_len)+"\t"+str(max_len)+"\t"+str(p010)+"\t"+str(p025)+"\t"+str(p050)+"\t"+str(p250)+"\t"+str(p500)+"\t"+str(p750)+"\t"+str(p950)+"\t"+str(p975)+"\t"+str(p990)+"\t"+str(low_fen)+"\t"+str(up_fen)+"\n")
            i+=1
        fw_stat.close()
        ## Transpose table
        df = pd.read_csv(out_stat, sep='\t')
        df_transposed = df.transpose()
        df_transposed.to_csv(out_stat, sep='\t', header=False)
        print("#>>> Saved statistics table as:",out_stat)
                
        # Create box plots
        out_box = themoor+"/length/results/visualization/"+level+"/"+primer+"_seqlen_"+level+"_BOXPLOT.png"
        data = flen_dict[primer]
        df = pd.DataFrame([{'Group': k, 'Value': v} for k, values in data.items() for v in values])
        plt.figure(figsize=(15, 10))
        sns.set_theme(style="darkgrid")
        sns.boxplot(data=df, x='Group', y='Value', palette="Set3", dodge=False)
        plt.title('Box plot of sequence length grouped by '+level+" - "+primer, fontsize=15)
        plt.xticks(rotation=90, fontsize=12)
        plt.yticks(fontsize=12)
        plt.xlabel(level, fontsize=14)
        plt.ylabel('Sequence length', fontsize=14)
        plt.savefig(out_box, dpi=300, bbox_inches='tight')
        print("#>>> Saved box plot as:",out_box)

        # Create density plot
        out_den = themoor+"/length/results/visualization/"+level+"/"+primer+"_seqlen_"+level+"_DENSITY.png"
        data = flen_dict[primer]
        df = pd.DataFrame([{'Group': k, 'Value': v} for k, values in data.items() for v in values])
        plt.figure(figsize=(15, 10))
        sns.set_theme(style="darkgrid")
        ax = sns.kdeplot(data=df, x='Value', hue='Group', fill=True, common_norm=False, linewidth=1.5, alpha=0.7, palette="Set3", warn_singular=False)
        sns.move_legend(ax, "best", title='Taxa', title_fontsize=14,fontsize=12,frameon=True, facecolor='white', framealpha=1, edgecolor='white')
        plt.title('Density plot of sequence length grouped by '+level+" - "+primer, fontsize=15)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.xlabel('Sequence length', fontsize=14)
        plt.ylabel('Density', fontsize=14)
        plt.savefig(out_den, dpi=300, bbox_inches='tight')
        print("#>>> Saved density plot as:",out_den)
    # endregion

    # endregion
#=======================================================================================================================
    # region Taxonomic resolution

    print("\n#======================================== Taxonomic resolution ==========================================")

    # region Create new directories
    prefix = target.split(".")[0]
    newdir1 = "mkdir "+themoor+"/resolution"
    newdir2 = "mkdir "+themoor+"/resolution/classifier"
    newdir3 = "mkdir "+themoor+"/resolution/obstax_"+prefix
    newdir4 = "mkdir "+themoor+"/resolution/evaluation_"+prefix
    newdir5 = "mkdir "+themoor+"/resolution/results_"+prefix
    newdir6 = "mkdir "+themoor+"/resolution/results_"+prefix+"/best_accuracy"
    newdir7 = "mkdir "+themoor+"/resolution/results_"+prefix+"/accuracy_comparison"
    os.system(newdir1)
    os.system(newdir2)
    os.system(newdir3)
    os.system(newdir4)
    os.system(newdir5)
    os.system(newdir6)
    os.system(newdir7)
    # endregion

    # region Train classifiers from primer-extracted datasets with "qiime feature-classifier fit-classifier-naive-bayes"
    print("\n#### Train QIIME2 classifiers from primer datasets ####")
    f = open(lowood,"r") #Modify input file here
    header = f.readline()
    for line in f:
        #Variables
        primer = line.split("\t")[0]
        gene = line.split("\t")[1]
        forward = line.split("\t")[2]
        reverse = line.split("\t")[3].strip()
        # Files
        in_seq = thornfield+"/clean_datasets/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_seq.qza"
        in_tax = thornfield+"/clean_datasets/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_tax.qza"
        out_cls = themoor+"/resolution/classifier/"+"txid"+taxid+"_"+primer+"_"+my+"_classifier.qza"
        # Commands
        cmd_cls = "qiime feature-classifier fit-classifier-naive-bayes \\\n"+"--i-reference-reads "+in_seq+" \\\n"+"--i-reference-taxonomy "+in_tax+" \\\n"+"--o-classifier "+out_cls
        #Running
        print("\n# "+primer+" >>> Running . . .")
        print(cmd_cls)
        os.system(cmd_cls)
    f.close()

    ##copy result to main result dir
    write_res_cmd8 = "cp -r "+themoor+"/resolution/classifier"+" " +themoor
    os.system(write_res_cmd8)
    # endregion

    # region Extract target taxa sequences and taxonomy from preprocessed primer-extracted dataset
    print("\n#### Extract target taxa sequence and taxonomy ####")
    ## Get list of target taxa to be extracted
    fx = open(target,"r")
    taxa_lev = str(fx.readline().strip())
    ext_taxa_list = []
    for line in fx:
        ext_taxa_list.append(line.strip())
    fx.close()
    print("\nExtract taxonomic level:",taxa_lev)
    print("Extract taxa:")
    for item in ext_taxa_list:
        print(item)
    print("---------------------------------------------------------------")

    ## Extract taxa
    f = open(lowood,"r")
    header = f.readline()
    for line in f:
        # Variables
        primer = line.split("\t")[0]
        gene = line.split("\t")[1]
        forward = line.split("\t")[2]
        reverse = line.split("\t")[3].strip()

        #Files
        in_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_seq/dna-sequences.fasta"
        in_tax = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_tax/taxonomy.tsv"
        out_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_seq/dna-sequences_"+prefix+".fasta"
        out_tax = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_tax/taxonomy_"+prefix+".tsv"

        #Running message
        print("\n# "+primer+" >>> Running . . . ")

        # 1) Create sequence dictionary
        ## Write sequence and put in a dictionary {accession number:sequence,...,...}
        print("\n# Create sequence dictionary")
        seq_dict = {}
        f_seq = open(in_seq, "r")
        for line in f_seq:
            if line[0] == ">":
                accss_id = line[1:].strip()
                seq_dict[accss_id] = ''
            else:
                sequence = line[0:].strip()
                seq_dict[accss_id] = sequence
        f_seq.close()
        print(">>> Created dictionary: seq_dict ---> {Keys = accession: values = sequences}")
        print("## Number of sequences written in a dictionary:", len(seq_dict))

        # 2) Create taxonomy dictionaries
        ## Write 2 levels of dictionaries 1) inferior = {taxonomic level:value,...} and 2) superior = {accession number:{taxonomic level:velue,...},...}
        print("\n# Create taxonomy dictionary")
        inf_taxa_dict = {}
        sup_taxa_dict = {}
        f_tax = open(in_tax, "r")
        f_tax.readline() #read header
        for line in f_tax:
            accss = line.split("\t")[0]
            Kingdom = line.split("\t")[1].strip().split(";")[0].split("__")[1]
            inf_taxa_dict["Kingdom"] = Kingdom #inf_taxa_dict = {kingdom:Metazoa}
            #print(Kingdom)
            Phylum = line.split("\t")[1].strip().split(";")[1].split("__")[1]
            inf_taxa_dict["Phylum"] = Phylum #inf_taxa_dict = {kingdom:Metazoa, phylum: ...}
            #print(Phylum)
            Class = line.split("\t")[1].strip().split(";")[2].split("__")[1]
            inf_taxa_dict["Class"] = Class #inf_taxa_dict = {kingdom:Metazoa, phylum: ..., class:...}
            #print(Class)
            Order = line.split("\t")[1].strip().split(";")[3].split("__")[1]
            inf_taxa_dict["Order"] = Order
            #print(Order)
            Family = line.split("\t")[1].strip().split(";")[4].split("__")[1]
            inf_taxa_dict["Family"] = Family
            #print(Family)
            Genus = line.split("\t")[1].strip().split(";")[5].split("__")[1]
            inf_taxa_dict["Genus"] = Genus
            #print(Genus)
            Species = line.split("\t")[1].strip().split(";")[6].split("__")[1]
            inf_taxa_dict["Species"] = Species #inf_taxa_dict = {kingdom:Metazoa, phylum: ..., class:... ,.., .species:...}
            #print(Species)
            sup_taxa_dict[accss] = inf_taxa_dict
            inf_taxa_dict = {}
        f_tax.close()
        print(">>> Created dictionary: sup_taxa_dict ---> {Keys = accession[taxonomic level], values = taxonomic information")
        print("## Number of taxonomic records written in a dictionary:", len(sup_taxa_dict))

        # 3) Write results as .fasta files
        print("\n# Write extracted taxa sequence files")
        fw = open(out_seq,"w")
        u=0
        for accession in sup_taxa_dict.keys():
            for item in ext_taxa_list:
                if sup_taxa_dict[accession][taxa_lev] == item:
                    fw.write(">"+accession+"\n"+seq_dict[accession]+"\n")
                    u+=1
        fw.close()
        print(">>> Saved file as:",out_seq)
        print("## Number of lines written:", u)

        # 4) Write results as .tsv files
        print("\n# Write extracted taxa taxonomy files")
        f_tax = open(in_tax, "r")
        header = f_tax.readline()
        fw = open(out_tax,"w")
        fw.write(header)
        u=0
        for line in f_tax:
            accession = line.split("\t")[0]
            for item in ext_taxa_list:
                if sup_taxa_dict[accession][taxa_lev] == item:
                    fw.write(line)
                    u+=1
        fw.close()
        f_tax.close()
        print(">>> Saved file as:",out_tax)
        print("## Number of lines written:", u)
        print("------------------------------------------------------------")
    f.close()
    # endregion

    # region Import target taxa sequence and taxonomy files to QIIME2 with "qiime tools import"
    print("\n#### Import target taxa sequence and taxonomy files to QIIME2 ####")
    f = open(lowood,"r")
    header = f.readline()
    for line in f:
        #Variables
        primer = line.split("\t")[0]
        gene = line.split("\t")[1]
        forward = line.split("\t")[2]
        reverse = line.split("\t")[3].strip()
        #Files
        in_seq = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_seq/dna-sequences_"+prefix+".fasta"
        in_tax = thornfield+"/export/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_tax/taxonomy_"+prefix+".tsv"
        out_seq = thornfield+"/clean_datasets/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_seq_"+prefix+".qza"
        out_tax = thornfield+"/clean_datasets/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_tax_"+prefix+".qza"
        # Commands
        cmd_impseq = "qiime tools import \\\n"+"--input-path "+in_seq+" \\\n"+"--output-path "+out_seq+" \\\n"+"--type "+"'FeatureData[Sequence]'"
        cmd_imptax = "qiime tools import \\\n"+"--input-path "+in_tax+" \\\n"+"--output-path "+out_tax+" \\\n"+"--type "+"'FeatureData[Taxonomy]'"
        #Running
        print("\n# "+primer+" >>> Running . . . ")
        print(cmd_impseq)
        os.system(cmd_impseq)
        print(cmd_imptax)
        os.system(cmd_imptax)
        print("-----------------------------------------------------------")
    f.close()
    # endregion

    # region Classify target taxa - varying confidence threshold of limiting taxonomic depth
    print("\n#### Classify target taxa - varying confidence threshold of limiting taxonomic depth ####")

    f = open(lowood,"r")
    header = f.readline()
    for line in f:
        #Variables
        primer = line.split("\t")[0]
        gene = line.split("\t")[1]
        forward = line.split("\t")[2]
        reverse = line.split("\t")[3].strip()
        # Files
        in_seq = thornfield+"/clean_datasets/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_seq_"+prefix+".qza"
        in_tax = thornfield+"/clean_datasets/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_tax_"+prefix+".qza"
        in_cls = themoor+"/resolution/classifier/"+"txid"+taxid+"_"+primer+"_"+my+"_classifier.qza"
        #Running
        print("\n# "+primer+" >>> Running . . . ")
        for item in p:
            print("\n## p-confidence = "+str(item))
            out_obs = themoor+"/resolution/obstax_"+prefix+"/"+primer+"_obstax-p"+str(item)+".qza"
            # Commands
            cmd_obstax = "qiime feature-classifier classify-sklearn \\\n"+"--i-classifier "+in_cls+" \\\n"+"--i-reads "+in_seq+" \\\n"+"--p-confidence "+str(item)+" \\\n"+"--o-classification "+out_obs
            print(cmd_obstax)
            os.system(cmd_obstax)
        print("#-----------------------------------------------------------------------------------------------------#")
    f.close()
    # endregion

    # region Evaluate classification results
    print("\n#### Evaluate classification results ####")

    ## [1.] Create input files and parameters dictionary
    ### Create lists & dict
    primer_dict = {}
    f = open(lowood,"r")
    header = f.readline()
    for line in f:
        para_dict = {}
        exptax_list = []
        obstax_list = []
        plabel_list = []
        #Variables
        primer = line.split("\t")[0]
        gene = line.split("\t")[1]
        forward = line.split("\t")[2]
        reverse = line.split("\t")[3].strip()
        for item in p:
            plabel = "p"+str(item)
            exptax = thornfield+"/clean_datasets/"+"txid"+taxid+"_"+primer+"_"+my+"_clean-db_tax_"+prefix+".qza"
            obstax = themoor+"/resolution/obstax_"+prefix+"/"+primer+"_obstax-"+plabel+".qza"
            plabel_list.append(plabel)
            exptax_list.append(exptax)
            obstax_list.append(obstax)
        para_dict["label"] = plabel_list
        para_dict["exp"] = exptax_list
        para_dict["obs"] = obstax_list
        primer_dict[primer] = para_dict
    f.close()

    ## [2.] Generate QIIME2 commands to evaluate classifiers
    f = open(lowood,"r")
    header = f.readline()
    for line in f:
        #Variables
        primer = line.split("\t")[0]
        gene = line.split("\t")[1]
        forward = line.split("\t")[2]
        reverse = line.split("\t")[3].strip()
        # File
        in_exptax = " ".join(primer_dict[primer]["exp"])
        in_obstax = " ".join(primer_dict[primer]["obs"])
        in_plabel = " ".join(primer_dict[primer]["label"])
        out_eval = themoor+"/resolution/evaluation_"+prefix+"/"+"txid"+taxid+"_"+primer+"_"+my+"_evaldb.qzv"
        #print(in_exptax)
        #print(in_obstax)
        #print(in_plabel)
        # Commands
        cmd_eval = "qiime rescript evaluate-classifications \\\n"+"--i-expected-taxonomies "+in_exptax+" \\\n"+"--i-observed-taxonomies "+in_obstax+" \\\n"+"--p-labels "+in_plabel+" \\\n"+"--o-evaluation "+out_eval
        # Running
        print("\n# "+primer+" >>> Running . . .")
        print(cmd_eval)
        os.system(cmd_eval)
        print("#-----------------------------------------------------------------------------------------------------#")
    f.close()
    # endregion

    # region Export evaluation results and observed taxonomy with "qiime tools export"
    print("\n#### Export evaluation results and observed taxonomy ####")

    f = open(lowood,"r")
    header = f.readline()
    for line in f:
        #Variables
        primer = line.split("\t")[0]
        gene = line.split("\t")[1]
        forward = line.split("\t")[2]
        reverse = line.split("\t")[3].strip()
        # File - eval
        in_eval = themoor+"/resolution/evaluation_"+prefix+"/"+"txid"+taxid+"_"+primer+"_"+my+"_evaldb.qzv"
        out_eval = themoor+"/resolution/evaluation_"+prefix+"/"+"txid"+taxid+"_"+primer+"_"+my+"_evaldb"
        # Commands - eval
        cmd_expeval = "qiime tools export \\\n"+"--input-path "+in_eval+" \\\n"+"--output-path "+out_eval
        # Running
        print("\n# "+primer+" >>> Running . . . ")
        print("\n## Primer evaluation results")
        print(cmd_expeval)
        os.system(cmd_expeval)
        print("\n## Observed taxonomy results")
        for item in p:
            # Files - obstax
            print("### p-confidence = "+str(item))
            in_obs = themoor+"/resolution/obstax_"+prefix+"/"+primer+"_obstax-p"+str(item)+".qza"
            out_obsvis = themoor+"/resolution/obstax_"+prefix+"/"+primer+"_obstax-p"+str(item)+".qzv"
            out_obsdir = themoor+"/resolution/obstax_"+prefix+"/"+primer+"_obstax-p"+str(item)
            # Commands - obstax
            cmd_obsvis = "qiime metadata tabulate \\\n"+"--m-input-file "+in_obs+" \\\n"+"--o-visualization "+out_obsvis
            print(cmd_obsvis)
            os.system(cmd_obsvis)
            cmd_obsdir = "qiime tools export \\\n"+"--input-path "+out_obsvis+" \\\n"+"--output-path "+out_obsdir
            print(cmd_obsdir+"\n")
            os.system(cmd_obsdir)
        print("#--------------------------------------------------------------------------------------------------------#")
    f.close()
    # endregion

    # region Manipulate results: Target taxa classification accuracy comparison among primers
    print("\n#### Manipulate results: compare accuracy ####")

    # region Make dictionaries
    f = open(lowood,"r")
    header = f.readline()
    primer_dict1 = {}
    for line in f:
        #Dictionary
        conf_dict2 = {}
        level_dict3 = {}
        accuracy_dict4 = {}
        value_list = []
        #Variables
        primer = line.split("\t")[0]
        gene = line.split("\t")[1]
        forward = line.split("\t")[2]
        reverse = line.split("\t")[3].strip()
        #Files
        in_eval = themoor+"/resolution/evaluation_"+prefix+"/"+"txid"+taxid+"_"+primer+"_"+my+"_evaldb/data.tsv"
        #Running
        fe = open(in_eval,"r")
        fe.readline()
        fe.readline()
        for line in fe:
            value_list = []
            level = line.split("\t")[1]
            precision = line.split("\t")[2]
            recall = line.split("\t")[3]
            fmeasure = line.split("\t")[4]
            confident = line.split("\t")[5].strip()
            conf_dict2[confident]={}
        fe.close()
        for c in conf_dict2.keys():
            level_dict3 = {}
            fe = open(in_eval,"r")
            fe.readline()
            fe.readline()
            for line in fe:
                value_list = []
                level = line.split("\t")[1]
                precision = line.split("\t")[2]
                recall = line.split("\t")[3]
                fmeasure = line.split("\t")[4]
                confident = line.split("\t")[5].strip()
                if c == confident:
                    value_list.append(precision)
                    value_list.append(recall)
                    value_list.append(fmeasure)
                    level_dict3[level] = value_list
                    conf_dict2[confident] = level_dict3
            fe.close()
        primer_dict1[primer] = conf_dict2
    f.close()
    # endregion

    # region Write results - Reports
    print("\n## Report classification accuracy compared among primer datasets: precision, recall, and f-measure")
    for item in p:
        item = str(item)
        print("\n# Confidence:",item)
        # make new directories
        outdir = "mkdir "+themoor+"/resolution/results_"+prefix+"/accuracy_comparison/confidence_"+item
        os.system(outdir)
        # files
        out_pre = themoor+"/resolution/results_"+prefix+"/accuracy_comparison/confidence_"+item+"/precision_REPORT.txt"
        out_rec = themoor+"/resolution/results_"+prefix+"/accuracy_comparison/confidence_"+item+"/recall_REPORT.txt"
        out_fme = themoor+"/resolution/results_"+prefix+"/accuracy_comparison/confidence_"+item+"/fmeasure_REPORT.txt"
        # variables
        plabel = "p"+item
        # running
        ## Precision
        ### Write results
        fp = open(out_pre,"w")
        fp.write("Level"+"\t"+"Kingdom"+"\t"+"Phylum"+"\t"+"Class"+"\t"+"Order"+"\t"+"Family"+"\t"+"Genus"+"\t"+"Species"+"\n")
        for primer in primer_dict1.keys():
            precision_K = primer_dict1[primer][plabel]["1"][0]
            precision_P = primer_dict1[primer][plabel]["2"][0]
            precision_C = primer_dict1[primer][plabel]["3"][0]
            precision_O = primer_dict1[primer][plabel]["4"][0]
            precision_F = primer_dict1[primer][plabel]["5"][0]
            precision_G = primer_dict1[primer][plabel]["6"][0]
            precision_S = primer_dict1[primer][plabel]["7"][0]
            fp.write(primer+"\t"+precision_K+"\t"+precision_P+"\t"+precision_C+"\t"+precision_O+"\t"+precision_F+"\t"+precision_G+"\t"+precision_S+"\n")
        fp.close()
        ### Transpose table
        df = pd.read_csv(out_pre, sep='\t')
        df_transposed = df.transpose()
        df_transposed.to_csv(out_pre, sep='\t', header=False)
        print("#>>> Saved precision report as:",out_pre)

        # Recall
        ### Write results
        fr = open(out_rec,"w")
        fr.write("Level"+"\t"+"Kingdom"+"\t"+"Phylum"+"\t"+"Class"+"\t"+"Order"+"\t"+"Family"+"\t"+"Genus"+"\t"+"Species"+"\n")
        for primer in primer_dict1.keys():
            recall_K = primer_dict1[primer][plabel]["1"][1]
            recall_P = primer_dict1[primer][plabel]["2"][1]
            recall_C = primer_dict1[primer][plabel]["3"][1]
            recall_O = primer_dict1[primer][plabel]["4"][1]
            recall_F = primer_dict1[primer][plabel]["5"][1]
            recall_G = primer_dict1[primer][plabel]["6"][1]
            recall_S = primer_dict1[primer][plabel]["7"][1]
            fr.write(primer+"\t"+recall_K+"\t"+recall_P+"\t"+recall_C+"\t"+recall_O+"\t"+recall_F+"\t"+recall_G+"\t"+recall_S+"\n")
        fr.close()
        ### Transpose table
        df = pd.read_csv(out_rec, sep='\t')
        df_transposed = df.transpose()
        df_transposed.to_csv(out_rec, sep='\t', header=False)
        print("#>>> Saved recall report as:",out_rec)

        # F-measure
        ### Write results
        ff = open(out_fme,"w")
        ff.write("Level"+"\t"+"Kingdom"+"\t"+"Phylum"+"\t"+"Class"+"\t"+"Order"+"\t"+"Family"+"\t"+"Genus"+"\t"+"Species"+"\n")
        for primer in primer_dict1.keys():
            fmeasure_K = primer_dict1[primer][plabel]["1"][2]
            fmeasure_P = primer_dict1[primer][plabel]["2"][2]
            fmeasure_C = primer_dict1[primer][plabel]["3"][2]
            fmeasure_O = primer_dict1[primer][plabel]["4"][2]
            fmeasure_F = primer_dict1[primer][plabel]["5"][2]
            fmeasure_G = primer_dict1[primer][plabel]["6"][2]
            fmeasure_S = primer_dict1[primer][plabel]["7"][2]
            ff.write(primer+"\t"+fmeasure_K+"\t"+fmeasure_P+"\t"+fmeasure_C+"\t"+fmeasure_O+"\t"+fmeasure_F+"\t"+fmeasure_G+"\t"+fmeasure_S+"\n")
        ff.close()
        ### Transpose table
        df = pd.read_csv(out_fme, sep='\t')
        df_transposed = df.transpose()
        df_transposed.to_csv(out_fme, sep='\t', header=False)
        print("#>>> Saved f-measure report as:",out_fme)
    # endregion

    # region Write results - Line plots
    print("\n## Visualize classification accuracy compared among primer datasets: precision, recall, and f-measure")

    for item in p:
        item = str(item)
        print("\n# Confidence:",item)
        in_pre = themoor+"/resolution/results_"+prefix+"/accuracy_comparison/confidence_"+item+"/precision_REPORT.txt"
        in_rec = themoor+"/resolution/results_"+prefix+"/accuracy_comparison/confidence_"+item+"/recall_REPORT.txt"
        in_fme = themoor+"/resolution/results_"+prefix+"/accuracy_comparison/confidence_"+item+"/fmeasure_REPORT.txt"
        out_figpre = themoor+"/resolution/results_"+prefix+"/accuracy_comparison/confidence_"+item+"/precision_LINEPLOT.png"
        out_figrec = themoor+"/resolution/results_"+prefix+"/accuracy_comparison/confidence_"+item+"/recall_LINEPLOT.png"
        out_figfme = themoor+"/resolution/results_"+prefix+"/accuracy_comparison/confidence_"+item+"/fmeasure_LINEPLOT.png"

        # Precision
        # 1. Load and Reshape
        df = pd.read_csv(in_pre, sep='\t')
        df_long = df.melt(id_vars='Level', var_name='Primer', value_name='Precision')
        original_order = df_long['Level'].unique()
        df_long['Level'] = pd.Categorical(df_long['Level'], categories=original_order, ordered=True)
        # 2. Set the style
        sns.set_theme(style="darkgrid")
        # 3. Create the plot
        plt.figure(figsize=(15, 10))
        sns.lineplot(data=df_long, x='Level', y='Precision', hue='Primer', palette='Set3', marker='o', markersize=9, linewidth=4)
        # 4. Add Title
        plt.title('Precision compared across all primers - '+item+" confidence", fontsize=18)
        plt.xlabel('Taxonomic level', fontsize=16, labelpad=10)
        plt.ylabel('Precision', fontsize=16, labelpad=10)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.legend(title='Primer', title_fontsize=16, fontsize=14, loc='center left', bbox_to_anchor=(1, 0.5), facecolor='#eaeaf2', edgecolor='#eaeaf2')
        # 5. Save plot
        plt.tight_layout()
        plt.savefig(out_figpre, dpi=300, bbox_inches='tight')
        plt.close()
        print("#>>> Saved precision line plot as:",out_figpre)

        # Recall
        # 1. Load and Reshape
        df = pd.read_csv(in_rec, sep='\t')
        df_long = df.melt(id_vars='Level', var_name='Primer', value_name='Recall')
        original_order = df_long['Level'].unique()
        df_long['Level'] = pd.Categorical(df_long['Level'], categories=original_order, ordered=True)
        # 2. Set the style
        sns.set_theme(style="darkgrid")
        # 3. Create the plot
        plt.figure(figsize=(15, 10))
        sns.lineplot(data=df_long, x='Level', y='Recall', hue='Primer', palette='Set3', marker='o', markersize=9, linewidth=4)
        # 4. Add Title
        plt.title('Recall compared across all primers - '+item+" confidence", fontsize=18)
        plt.xlabel('Taxonomic level', fontsize=16, labelpad=10)
        plt.ylabel('Recall', fontsize=16, labelpad=10)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.legend(title='Primer', title_fontsize=16, fontsize=14, loc='center left', bbox_to_anchor=(1, 0.5), facecolor='#eaeaf2', edgecolor='#eaeaf2')
        # 5. Save plot
        plt.tight_layout()
        plt.savefig(out_figrec, dpi=300, bbox_inches='tight')
        plt.close()
        print("#>>> Saved recall line plot as:",out_figrec)

        # F-measure
        # 1. Load and Reshape
        df = pd.read_csv(in_fme, sep='\t')
        df_long = df.melt(id_vars='Level', var_name='Primer', value_name='F-measure')
        original_order = df_long['Level'].unique()
        df_long['Level'] = pd.Categorical(df_long['Level'], categories=original_order, ordered=True)
        # 2. Set the style
        sns.set_theme(style="darkgrid")
        # 3. Create the plot
        plt.figure(figsize=(15, 10))
        sns.lineplot(data=df_long, x='Level', y='F-measure', hue='Primer', palette='Set3', marker='o', markersize=9, linewidth=4)
        # 4. Add Title
        plt.title('F-measure compared across all primers - '+item+" confidence", fontsize=18)
        plt.xlabel('Taxonomic level', fontsize=16, labelpad=10)
        plt.ylabel('F-measure', fontsize=16, labelpad=10)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.legend(title='Primer', title_fontsize=16, fontsize=14, loc='center left', bbox_to_anchor=(1, 0.5), facecolor='#eaeaf2', edgecolor='#eaeaf2')
        # 5. Save plot
        plt.tight_layout()
        plt.savefig(out_figfme, dpi=300, bbox_inches='tight')
        plt.close()
        print("#>>> Saved F-measure line plot as:",out_figfme)
    # endregion
    # endregion

    # region Manipulate results: Find confidence threshold giving best classification accuracy
    print("\n#### Manipulate results: best accuracy ####")

    # region Make dictionaries
    f = open(lowood,"r")
    header = f.readline()
    primer_dict1 = {}
    for line in f:
        #Dictionary
        conf_dict2 = {}
        level_dict3 = {}
        accuracy_dict4 = {}
        value_list = []
        #Variables
        primer = line.split("\t")[0]
        gene = line.split("\t")[1]
        forward = line.split("\t")[2]
        reverse = line.split("\t")[3].strip()
        #Files
        in_eval = themoor+"/resolution/evaluation_"+prefix+"/"+"txid"+taxid+"_"+primer+"_"+my+"_evaldb/data.tsv"
        #Running
        fe = open(in_eval,"r")
        fe.readline()
        fe.readline()
        for line in fe:
            value_list = []
            level = line.split("\t")[1]
            precision = line.split("\t")[2]
            recall = line.split("\t")[3]
            fmeasure = line.split("\t")[4]
            confident = line.split("\t")[5].strip()
            conf_dict2[confident]={}
        fe.close()

        for c in conf_dict2.keys():
            level_dict3 = {}
            fe = open(in_eval,"r")
            fe.readline()
            fe.readline()
            for line in fe:
                value_list = []
                level = line.split("\t")[1]
                precision = line.split("\t")[2]
                recall = line.split("\t")[3]
                fmeasure = line.split("\t")[4]
                confident = line.split("\t")[5].strip()
                if c == confident:
                    value_list.append(precision)
                    value_list.append(recall)
                    value_list.append(fmeasure)
                    level_dict3[level] = value_list
                    conf_dict2[confident] = level_dict3
            fe.close()
        primer_dict1[primer] = conf_dict2
    f.close()
    # endregion

    # region Best Precision
    precision_dict = {}
    recall_dict = {}
    fmeasure_dict = {}
    confidence_dict = {}
    for primer in primer_dict1.keys():
        precision_list = []
        recall_list = []
        fmeasure_list = []
        confidence_list = []
        for item in p:
            # variables
            plabel = "p"+str(item)
            # running
            precision = float(primer_dict1[primer][plabel]["7"][0])
            recall = float(primer_dict1[primer][plabel]["7"][1])
            fmeasure = float(primer_dict1[primer][plabel]["7"][2])
            confidence = float(item)
            # add to list
            precision_list.append(precision)
            recall_list.append(recall)
            fmeasure_list.append(fmeasure)
            confidence_list.append(confidence)
        precision_dict[primer] = precision_list
        recall_dict[primer] = recall_list
        fmeasure_dict[primer] = fmeasure_list
        confidence_dict[primer] = confidence_list

    new_precision_dict = {}
    new_recall_dict = {}
    new_fmeasure_dict = {}
    new_confidence_dict = {}
    for primer in confidence_dict.keys():
        new_precision_list = []
        new_recall_list = []
        new_fmeasure_list = []
        new_confidence_list = []
        max_precision = max(precision_dict[primer])
        for i in range(len(confidence_dict[primer])):
            if max_precision == precision_dict[primer][i]:
                new_precision_list.append(precision_dict[primer][i])
                new_recall_list.append(recall_dict[primer][i])
                new_fmeasure_list.append(fmeasure_dict[primer][i])
                new_confidence_list.append(confidence_dict[primer][i])
        new_precision_dict[primer] = new_precision_list
        new_recall_dict[primer] = new_recall_list
        new_fmeasure_dict[primer] = new_fmeasure_list
        new_confidence_dict[primer] = new_confidence_list

    pre_out_res = themoor+"/resolution/results_"+prefix+"/best_accuracy/best_precision.txt"
    fw = open(pre_out_res, "w")
    fw.write("accuracy"+"\t"+"confidence"+"\t"+"precision"+"\t"+"recall"+"\t"+"f-measure"+"\n")
    for primer in new_confidence_dict.keys():
        min_confidence = min(new_confidence_dict[primer])
        for i in range(len(new_confidence_dict[primer])):
            if min_confidence == new_confidence_dict[primer][i]:
                fw.write(primer+"\t"+str(new_confidence_dict[primer][i])+"\t"+str(new_precision_dict[primer][i])+"\t"+str(new_recall_dict[primer][i])+"\t"+str(new_fmeasure_dict[primer][i])+"\n")
    fw.close()
    df = pd.read_csv(pre_out_res, sep='\t')
    df_transposed = df.transpose()
    df_transposed.to_csv(pre_out_res, sep='\t', header=False)
    print("#>>> Saved best precision report as:",pre_out_res)
    # endregion

    # region Best Recall
    precision_dict = {}
    recall_dict = {}
    fmeasure_dict = {}
    confidence_dict = {}
    for primer in primer_dict1.keys():
        precision_list = []
        recall_list = []
        fmeasure_list = []
        confidence_list = []
        for item in p:
            # variables
            plabel = "p"+str(item)
            # running
            precision = float(primer_dict1[primer][plabel]["7"][0])
            recall = float(primer_dict1[primer][plabel]["7"][1])
            fmeasure = float(primer_dict1[primer][plabel]["7"][2])
            confidence = float(item)
            # add to list
            precision_list.append(precision)
            recall_list.append(recall)
            fmeasure_list.append(fmeasure)
            confidence_list.append(confidence)
        precision_dict[primer] = precision_list
        recall_dict[primer] = recall_list
        fmeasure_dict[primer] = fmeasure_list
        confidence_dict[primer] = confidence_list

    new_precision_dict = {}
    new_recall_dict = {}
    new_fmeasure_dict = {}
    new_confidence_dict = {}
    for primer in confidence_dict.keys():
        new_precision_list = []
        new_recall_list = []
        new_fmeasure_list = []
        new_confidence_list = []
        max_recall = max(recall_dict[primer])
        for i in range(len(confidence_dict[primer])):
            if max_recall == recall_dict[primer][i]:
                new_precision_list.append(precision_dict[primer][i])
                new_recall_list.append(recall_dict[primer][i])
                new_fmeasure_list.append(fmeasure_dict[primer][i])
                new_confidence_list.append(confidence_dict[primer][i])
        new_precision_dict[primer] = new_precision_list
        new_recall_dict[primer] = new_recall_list
        new_fmeasure_dict[primer] = new_fmeasure_list
        new_confidence_dict[primer] = new_confidence_list

    rec_out_res = themoor+"/resolution/results_"+prefix+"/best_accuracy/best_recall.txt"
    fw = open(rec_out_res, "w")
    fw.write("accuracy"+"\t"+"confidence"+"\t"+"precision"+"\t"+"recall"+"\t"+"f-measure"+"\n")
    for primer in new_confidence_dict.keys():
        min_confidence = min(new_confidence_dict[primer])
        for i in range(len(new_confidence_dict[primer])):
            if min_confidence == new_confidence_dict[primer][i]:
                fw.write(primer+"\t"+str(new_confidence_dict[primer][i])+"\t"+str(new_precision_dict[primer][i])+"\t"+str(new_recall_dict[primer][i])+"\t"+str(new_fmeasure_dict[primer][i])+"\n")
    fw.close()
    df = pd.read_csv(rec_out_res, sep='\t')
    df_transposed = df.transpose()
    df_transposed.to_csv(rec_out_res, sep='\t', header=False)
    print("#>>> Saved best recall report as:",rec_out_res)
    # endregion

    # region Best F-measure
    precision_dict = {}
    recall_dict = {}
    fmeasure_dict = {}
    confidence_dict = {}
    for primer in primer_dict1.keys():
        precision_list = []
        recall_list = []
        fmeasure_list = []
        confidence_list = []
        for item in p:
            # variables
            plabel = "p"+str(item)
            # running
            precision = float(primer_dict1[primer][plabel]["7"][0])
            recall = float(primer_dict1[primer][plabel]["7"][1])
            fmeasure = float(primer_dict1[primer][plabel]["7"][2])
            confidence = float(item)
            # add to list
            precision_list.append(precision)
            recall_list.append(recall)
            fmeasure_list.append(fmeasure)
            confidence_list.append(confidence)
        precision_dict[primer] = precision_list
        recall_dict[primer] = recall_list
        fmeasure_dict[primer] = fmeasure_list
        confidence_dict[primer] = confidence_list

    new_precision_dict = {}
    new_recall_dict = {}
    new_fmeasure_dict = {}
    new_confidence_dict = {}
    for primer in confidence_dict.keys():
        new_precision_list = []
        new_recall_list = []
        new_fmeasure_list = []
        new_confidence_list = []
        max_fmeasure = max(fmeasure_dict[primer])
        for i in range(len(confidence_dict[primer])):
            if max_fmeasure == fmeasure_dict[primer][i]:
                new_precision_list.append(precision_dict[primer][i])
                new_recall_list.append(recall_dict[primer][i])
                new_fmeasure_list.append(fmeasure_dict[primer][i])
                new_confidence_list.append(confidence_dict[primer][i])
        new_precision_dict[primer] = new_precision_list
        new_recall_dict[primer] = new_recall_list
        new_fmeasure_dict[primer] = new_fmeasure_list
        new_confidence_dict[primer] = new_confidence_list

    fme_out_res = themoor+"/resolution/results_"+prefix+"/best_accuracy/best_fmeasure.txt"
    fw = open(fme_out_res, "w")
    fw.write("accuracy"+"\t"+"confidence"+"\t"+"precision"+"\t"+"recall"+"\t"+"f-measure"+"\n")
    for primer in new_confidence_dict.keys():
        min_confidence = min(new_confidence_dict[primer])
        for i in range(len(new_confidence_dict[primer])):
            if min_confidence == new_confidence_dict[primer][i]:
                fw.write(primer+"\t"+str(new_confidence_dict[primer][i])+"\t"+str(new_precision_dict[primer][i])+"\t"+str(new_recall_dict[primer][i])+"\t"+str(new_fmeasure_dict[primer][i])+"\n")
    fw.close()
    df = pd.read_csv(fme_out_res, sep='\t')
    df_transposed = df.transpose()
    df_transposed.to_csv(fme_out_res, sep='\t', header=False)
    print("#>>> Saved best F-measure report as:",fme_out_res)
    # endregion
    # endregion

    ##copy result to main result dir
    write_res_cmd9 = "cp " + out_pre.replace("confidence_1","confidence_0") + " " + themoor + "/results/resolution_" + prefix + "_" + out_pre.split("/")[5]
    os.system(write_res_cmd9)
    write_res_cmd10 = "cp " + out_rec.replace("confidence_1","confidence_0") + " " + themoor + "/results/resolution_" + prefix + "_" + out_rec.split("/")[5]
    os.system(write_res_cmd10)
    write_res_cmd11 = "cp " + out_fme.replace("confidence_1","confidence_0") + " " + themoor + "/results/resolution_" + prefix + "_" + out_fme.split("/")[5]
    os.system(write_res_cmd11)
    write_res_cmd12 = "cp " + out_figpre.replace("confidence_1","confidence_0") + " " + themoor + "/results/resolution_" + prefix + "_" + out_figpre.split("/")[5]
    os.system(write_res_cmd12)
    write_res_cmd13 = "cp " + out_figrec.replace("confidence_1","confidence_0") + " " + themoor + "/results/resolution_" + prefix + "_" + out_figrec.split("/")[5]
    os.system(write_res_cmd13)
    write_res_cmd14 = "cp " + out_figfme.replace("confidence_1","confidence_0") + " " + themoor + "/results/resolution_" + prefix + "_" + out_figfme.split("/")[5]
    os.system(write_res_cmd14)
    write_res_cmd15 = "cp " + pre_out_res + " " + themoor + "/results/confidence_with_" + pre_out_res.split("/")[4]
    os.system(write_res_cmd15)
    write_res_cmd16 = "cp " + rec_out_res + " " + themoor + "/results/confidence_with_" + rec_out_res.split("/")[4]
    os.system(write_res_cmd16)
    write_res_cmd17 = "cp " + fme_out_res + " " + themoor + "/results/confidence_with_" + fme_out_res.split("/")[4]
    os.system(write_res_cmd17)

    #endregion
#=======================================================================================================================

# endregion
########################################################################################################################

def main():
    # Main parser ======================================================================================================
    parser = argparse.ArgumentParser(description="YAne: A high-Yielding primer Assessment and database construction tool for eDNA metabarcoding")
    subparsers = parser.add_subparsers(dest="command", help="Available functions")

    # Def#1 downloadmock ===============================================================================================
    parser_dl = subparsers.add_parser("downloadmock", help="Download and dereplicate mock community sequences and taxonomic information from NCBI")
    ## Required arguments
    parser_dl.add_argument("--i-mock", type=str, required=True, help="Mock community table file path [required]")
    parser_dl.add_argument("--p-taxa", type=str, required=True, help="Mock community Taxa ID [String; required]")
    parser_dl.add_argument("--p-dl-period", type=str, required=True, help="Data downloading period [String; required]")
    ## Optional arguments
    parser_dl.add_argument("--p-jobs", type=int, default=5, help="Number of concurrent download connection [Integer, default = 5]")
    parser_dl.add_argument("--p-threads", type=int, default=4, help="Number of computation threads [Integer, default = 4]")
    parser_dl.set_defaults(func=downloadmock)

    # Def#2 processprimer ==============================================================================================
    parser_pc = subparsers.add_parser("processprimer", help="Conduct in-silico PCR against the mock community and clean the datasets through homopolymer, degenerate base, and sequence lengths filtering steps")
    ## Required arguments
    parser_pc.add_argument("--i-mock", type=str, required=True, help="Mock community table file path [required]")
    parser_pc.add_argument("--i-primer", type=str, required=True, help="Primer table file path [required]")
    parser_pc.add_argument("--o-dir", type=str, required=True, help="Output directory name [required]")
    ## Optional arguments
    parser_pc.add_argument("--p-jobs", type=int, default=5, help="Number of concurrent download connection [integer, default = 5]")
    parser_pc.add_argument("--p-threads", type=int, default=4, help="Number of computation threads [integer, default = 4]")
    parser_pc.add_argument("--p-identity", type=float, default=0.8, nargs = "+", help="Combined primer matching identity (can be set globally for all datasets or customized for each primer dataset). [Float (0-1); default = 0.8]")
    parser_pc.add_argument("--p-homopolymer", type=int, default=8, nargs = "+", help="Homopolymer length filtering criteria (can be set globally for all datasets or customized for each primer dataset). Sequences with homopolymer length equal to or more than the criteria will be discarded. [Integer, default = 8]")
    parser_pc.add_argument("--p-degenerate", type=int, default=5, nargs = "+", help="Degenerate base length filtering criteria (can be set globally for all datasets or customized for each primer dataset). Sequences with degenerate base length equal to or more than the criteria will be discarded. [Integer, default = 5]")
    parser_pc.add_argument("--p-label", type=str, default="", nargs = "+", help="Labels of the primer datasets corresponding to the primer-specific customized parameters (identity/homopolymer/degenerate base). The labels must be similar to that in the primer table and arrange in the same order as the customized parameters.")
    parser_pc.add_argument("--p-observe-length", type=str, default="no", help="Observe amplicon length before customization. Enable this flag before setting customized amplicon length filtering criteria with the function customlength. The package will provide amplicon length report, statistics, and distribution to assist filtering criteria decisions. [yes/no; default = no (use 1st and 99th percentiles length filtering criteria)]")
    parser_pc.add_argument("--p-length-stat-level", type=str, default="Class", help="If --p-observe-length is 'yes', amplicon length statistics will be reported and can be grouped by taxonomic level. The desired taxonomic level grouping can be determined with this flag. [Taxonomic levels: Kingdom, Phylum, Class, Order, Family, Genus, Species; Default = Class]")
    parser_pc.set_defaults(func=processprimer)
    
    # Def#3 customlength ==============================================================================================
    parser_cl = subparsers.add_parser("customlength", help="Filter sequence lengths based on observed length statistics and distribution obtained from the function processprimer (--p-observe-length = yes). Records with lengths less or more than the minimum or maximum criteria are removed from the datasets.")
    parser_cl.add_argument("--p-label", type=str, required=True, nargs = "+", help="Labels of the primer datasets corresponding to the primer-specific customized parameters [required]")
    parser_cl.add_argument("--p-min-length", default = 0, type = int, nargs = "+", help="Minimum length [Integer; required for length filtering of all taxa, but not required for taxa-specific filtering]")
    parser_cl.add_argument("--p-max-length", default = 0, type = int, nargs = "+", help="Maximum length [Integer; required for length filtering of all taxa, but not required for taxa-specific filtering]")
    parser_cl.add_argument("--p-fillen-by-taxon", default = "",type = str, nargs = "+", help="Length filtering criteria file path [required for length filtering specific to taxa, but not required for all-taxa filtering] ")
    parser_cl.set_defaults(func=customlength)

    # Def#4 evaluateprimer ============================================================================================
    parser_ep = subparsers.add_parser("evaluateprimer", help="Analyze the primer-extracted datasets and generates comparative tables and graphical representations for primer evaluation in four aspects; primer specificity, reference database availability and intraspecific genetic variation, amplicon length, and taxonomic resolution.")
    ## Required arguments
    parser_ep.add_argument("--i-mock", type=str, required=True, help="Mock community table file path [required]")
    parser_ep.add_argument("--i-primer", type=str, required=True, help="Primer table file path [required]")
    parser_ep.add_argument("--i-dir", type=str, required=True, help="	Directory name of the output from the data preparation step [required]")
    parser_ep.add_argument("--i-target-taxa", type=str, required=True, help="Target taxa table file path [required]")
    parser_ep.add_argument("--o-dir", type=str, required=True, help="Directory name of output for this step (primer evaluation) [required]")
    ## Optional arguments
    parser_ep.add_argument("--p-tax-level", type=str, default="Class", help="Taxonomic level for grouping in primer evaluation results [Taxonomic levels: Kingdom, Phylum, Class, Order, Family, Genus, and Species; default = Class].")
    parser_ep.add_argument("--p-vary-confidence", type=float, nargs = "+", default=[0,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.96,0.97,0.98,0.99,1], help="Confidence values to be varied for target taxa classification [default = 0,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.96,0.97,0.98,0.99,1]")
    parser_ep.set_defaults(func=evaluateprimer)

    # Parse and Execute ================================================================================================
    args = parser.parse_args()
    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()