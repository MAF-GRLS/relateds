#!/bin/bash
## Pedigree analysis


mamba install -c conda-forge ghostscript
mamba install -c conda-forge r-igraph
mamba install -c conda-forge r-kinship2
# Download KING executable in your path.
wget https://www.kingrelatedness.com/Linux-king.tar.gz
tar -xzvf Linux-king.tar.gz

# Define your input file (in PLINK binary format: .bed, .bim, .fam)
INPUT_FILE="AxiomGT1v2.comp_merge.deDup.sexConfirm.public_ids" ## end of '3.Prep Array Sets'

# prepare more information from dog profiles
AxiomGT1v2.comp_merge.deDup.sexConfirm.public_ids.fam
awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$3]=$4;next}{if(a[$1])print a[$1],$0}' tmp_map/map_id.tab <(tail -n+2 phenotypes/dog_profile.csv | tr ',' '\t') > dog_profile_grls.tab
awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1]=1;next}{if(!a[$3])print}' <(tail -n+2 phenotypes/dog_profile.csv | tr ',' '\t') tmp_map/map_id.tab > dog_profile_grls.missing.tab ## all missing are oldies


# Prepare the Plink input files: change chr format & make each sample in a seperate family
plink --bfile ${INPUT_FILE} --chr-set 38 no-xy --allow-extra-chr '0' \
      --geno 0.05 --maf 0.01 \
      --make-bed --out ${INPUT_FILE}.king
cp ${INPUT_FILE}.king.fam ${INPUT_FILE}.king.fam_temp
awk '{$1=$2;print}' ${INPUT_FILE}.king.fam_temp > ${INPUT_FILE}.king.fam

# Step 1: Run KING to estimate kinship coefficients.
# -b specifies the input PLINK binary file.
#./king -b ${INPUT_FILE}.king.bed --kinship --prefix ${INPUT_FILE}.kingOutput &> king.kinship.log

# Step 2: Predict the pedigree based on kinship estimates.
# This will produce a .ped file with predicted relationships.
./king -b ${INPUT_FILE}.king.bed --related --degree 2 --sexchr 39 --prefix ${INPUT_FILE}.kingRelatedOutput --rplot &> king.related.log
./king -b ${INPUT_FILE}.king.bed --cluster --degree 2 --sexchr 39 --prefix ${INPUT_FILE}.kingClusterOutput --rplot &> king.cluster.log
#./king -b ${INPUT_FILE}.king.bed --build --degree 2 --sexchr 39 --prefix ${INPUT_FILE}.kingBuildOutput --rplot &> king.build.log
#./king -b ${INPUT_FILE}.king.bed --build --degree 2 --sexchr 39 --prefix ${INPUT_FILE}.kingBuildOutput &> king.build.log ## --degree 2 causes core dump
./king -b ${INPUT_FILE}.king.bed --build --sexchr 39 --prefix ${INPUT_FILE}.kingBuildOutput &> king.build.log

## merge outputs from related and cluster
tail -n+2 AxiomGT1v2.comp_merge.deDup.sexConfirm.public_ids.kingRelatedOutput.kin0 | cut -f2,4,14 | \
while read id1 id2 rel;do
  printf "Related $id1 $id2 $rel " | tr ' ' '\t';
  echo "cluster" $(grep -w "$id1" AxiomGT1v2.comp_merge.deDup.sexConfirm.public_ids.kingClusterOutputcluster.kin | grep -w "$id2" | cut -f1-3,15) | tr ' ' '\t';done > Related_cluster.tab
cat Related_cluster.tab | awk '{if($9 && $4!=$9)print}' ## perfect match
cat Related_cluster.tab | awk 'BEGIN{FS=OFS="\t"}{if(!$6) $6="-";print $2,$3,$4,$6}' > Related_cluster_short.tab

':
## Python
import pandas as pd
import networkx as nx
import sys

def read_file(file_path):
    """
    Reads the tab-separated file and returns a list of edges
    using only the first two columns.
    """
    data = pd.read_csv(file_path, sep='\t', header=None, usecols=[0, 1])
    edges = list(data.itertuples(index=False, name=None))
    return edges

def find_connected_components(edges):
    """
    Finds and returns the connected components given a list of edges.
    """
    G = nx.Graph()
    G.add_edges_from(edges)
    connected_components = list(nx.connected_components(G))
    return connected_components

def main(file_path):
    edges = read_file(file_path)
    components = find_connected_components(edges)
    
    for i, component in enumerate(components, start=1):
        print(f"Component {i}: {component}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <path_to_file>")
    else:
        file_path = sys.argv[1]
        main(file_path)
'
## Save this script as 'connected_components.py'
python connected_components.py Related_cluster_short.tab > Related_cluster_all_conn.txt
cat Related_cluster_short.tab | awk 'BEGIN{FS=OFS="\t"}{if($3!="2nd")print $0}' > Related_cluster_closeRel_short.tab
python connected_components.py Related_cluster_closeRel_short.tab > Related_cluster_closeRel_conn.txt

rclone -v --copy-links copy Related_cluster_short.tab remote_UCDavis_GoogleDr:MAF/forSequencing/related/
rclone -v --copy-links copy Related_cluster_all_conn.txt remote_UCDavis_GoogleDr:MAF/forSequencing/related/
rclone -v --copy-links copy Related_cluster_closeRel_conn.txt remote_UCDavis_GoogleDr:MAF/forSequencing/related/


## The build output
# 1. errors
cat AxiomGT1v2.comp_merge.deDup.sexConfirm.public_ids.kingBuildOutputupdateparents.txt | awk '{if($2==$3)print}' ## 6
cat AxiomGT1v2.comp_merge.deDup.sexConfirm.public_ids.kingBuildOutputupdateparents.txt | awk '{if($2==$4)print}' ## 11

# 2. confirm FS && PO relationship && fix the issue of the match between a parent and sample ID
cat AxiomGT1v2.comp_merge.deDup.sexConfirm.public_ids.kingBuildOutputupdateids.txt | sed 's/KING/FAM/' > AxiomGT1v2.comp_merge.deDup.sexConfirm.public_ids.kingBuildOutputupdateids_FAM.txt
cat AxiomGT1v2.comp_merge.deDup.sexConfirm.public_ids.kingBuildOutputupdateparents.txt | sed 's/KING/FAM/' | awk 'BEGIN{FS=OFS="\t"}{if($3==$2)$3=0;if($4==$2)$4=0;print $0}' > AxiomGT1v2.comp_merge.deDup.sexConfirm.public_ids.kingBuildOutputupdateparents_FAM.txt


awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$2]=$0;next}{if(a[$2])print a[$2];else print "-","-","-","-";if(a[$4])print a[$4];else print "-","-","-","-";}' AxiomGT1v2.comp_merge.deDup.sexConfirm.public_ids.kingBuildOutputupdateparents_FAM.txt <(tail -n+2 AxiomGT1v2.comp_merge.deDup.sexConfirm.public_ids.kingRelatedOutput.kin0) | paste - - > Related_families.tab

paste Related_cluster_short.tab Related_families.tab > Related_cluster_families.tab


# 2a) PO relationship
cat AxiomGT1v2.comp_merge.deDup.sexConfirm.public_ids.kingBuildOutputupdateparents_FAM.txt | awk 'BEGIN{FS=OFS="\t"}{if($2!=$3 && $3 !~ /^[0-9]+$/)print $2,$3;}' > build.po1
cat AxiomGT1v2.comp_merge.deDup.sexConfirm.public_ids.kingBuildOutputupdateparents_FAM.txt | awk 'BEGIN{FS=OFS="\t"}{if($2!=$4 && $4 !~ /^[0-9]+$/)print $2,$4;}' > build.po2
cat build.po1 build.po2 | sort | uniq > build.po
grep -w PO Related_cluster_short.tab | while read id1 id2 moreData;do grep -w $id1 build.po | grep -w $id2;done > build.po.confirmed_temp
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1"."$2]=1;next}{if(a[$2"."$1])print}' build.po.confirmed_temp build.po.confirmed_temp ## 2 pairs of bidirectional PO

cat Related_cluster_families.tab | awk 'BEGIN{FS=OFS="\t"}{if($3=="PO" && ($2==$7 || $2==$8) && ($1!=$11 && $1!=$12))print}' > Related_cluster_families.po.confirmed.tab
cat Related_cluster_families.tab | awk 'BEGIN{FS=OFS="\t"}{if($3=="PO" && ($2!=$7 && $2!=$8) && ($1==$11 || $1==$12))print $2,$1,$3,$4,$9,$10,$11,$12,$5,$6,$7,$8}' >> Related_cluster_families.po.confirmed.tab
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1"."$2]=1;next}{if(a[$1"."$2])print;}' build.po.confirmed_temp Related_cluster_families.po.confirmed.tab > build.po.confirmed ## 79

awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1]=$4;next}{if(!a[$1])a[$1]="-";if(!a[$2])a[$2]="-";print $0,a[$1],a[$2]}' dog_profile_grls.tab build.po.confirmed > build.po.confirmed_withProfile
awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$4]=$1;next}{if($13=="-")$13=a[$1];if($14=="-")$14=a[$2];print $0}' tmp_map/map_id.tab build.po.confirmed_withProfile > build.po.confirmed_withProfile2

cat build.po.confirmed_withProfile2 | awk '{if($13=="Oldies")print $1;if($14=="Oldies")print $2;}' | wc -l ## 46
cat build.po.confirmed_withProfile2 | awk '{if($13=="Oldies")print $1;if($14=="Oldies")print $2;}' | sort | uniq | wc -l ## 19


# Function to calculate difference in months between two dates
calculate_month_diff() {
    local start_date=$1
    local end_date=$2

    # Extract year and month from the dates
    local start_year=$(echo $start_date | cut -d'-' -f1)
    local start_month=$(echo $start_date | cut -d'-' -f2)
    local end_year=$(echo $end_date | cut -d'-' -f1)
    local end_month=$(echo $end_date | cut -d'-' -f2)

    # Force Bash to interpret numbers as decimal
    start_year=$((10#$start_year))
    start_month=$((10#$start_month))
    end_year=$((10#$end_year))
    end_month=$((10#$end_month))

    # Calculate difference in months
    local month_diff=$(( (end_year - start_year) * 12 + end_month - start_month ))
    echo $month_diff
}

# Read from the file line by line
while read end_date start_date; do
    # Calculate and print the month difference
    month_diff=$(calculate_month_diff $start_date $end_date)
    echo $month_diff
done < <(cat build.po.confirmed_withProfile2 | sed 's/Oldies/2005-01/' | cut -f13,14) > time_diff
paste build.po.confirmed_withProfile2 time_diff | awk '{if($15>0)print $1,$2,$13,$14}' > build.po.confirmed_final
paste build.po.confirmed_withProfile2 time_diff | awk '{if($15<0)print}' > build.po.confirmed_negative

rclone -v --copy-links copy build.po.confirmed_final remote_UCDavis_GoogleDr:MAF/Pedigree/


# 2b) FS relationship (check if they have the same parents ==> FS
cat Related_cluster_families.tab | awk 'BEGIN{FS=OFS="\t"}{if($3=="FS" && $7!="-" && $8!="-" && (($1==$11 || $1==$12) || ($2==$7 || $2==$8)))print}' > build.fs_po.toBeAssessed ## 117 ## Detected as FS but was also reported as PO
rclone -v --copy-links copy build.fs_po.toBeAssessed remote_UCDavis_GoogleDr:MAF/Pedigree/


cat Related_cluster_families.tab | awk 'BEGIN{FS=OFS="\t"}{if($3=="FS" && $7!="-" && $8!="-" && $7==$11 && $8==$12)print}' > build.fs.confirmed ## 4452
cat build.fs.confirmed | cut -f5-8 > build.fs.confirmed_ind  && cat build.fs.confirmed | cut -f9-12 >> build.fs.confirmed_ind
cat build.fs.confirmed_ind | sort | uniq | wc -l ## 1757

cat build.fs.confirmed | sort -k7,8 > build.fs.confirmed.sorted
cut -f7,8 build.fs.confirmed.sorted | uniq | wc -l  ## 425
mkdir families_ped
cut -f5,7,8 build.fs.confirmed.sorted | uniq | while read fam id1 id2;do
  cat build.fs.confirmed.sorted | awk -v id1="$id1" -v id2="$id2" 'BEGIN{FS=OFS="\t"}{if($7==id1 && $8==id2)print}' | cut -f1-5,7,8 > families_ped/$fam.$id1.$id2.ped
  #cat build.fs.confirmed | grep -w "$id1" | grep -w "$id2" | cut -f1-5,7,8 > families_ped/$fam.$id1.$id2.ped
done

for f in families_ped/*.ped;do
  echo $f
  cat $f | awk 'BEGIN{FS=OFS="\t"}{a[$1]+=1;a[$2]+=1;}END{for (i in a)print i,a[i]}'
done > families.consistancy

awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1]=$0;next}{if(a[$1])print $2,a[$1];else print $1}' dog_profile_grls.tab families.consistancy > families.consistancy_withProfile
awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1]=$0;next}/^families_ped/{print;next}{if(a[$1])print $2,a[$1];else print $2,$1}' dog_profile_grls.tab families.consistancy > families.consistancy_withProfile


# 4. dogs identified as FS by --related/--cluster commands are predicted to have totally different parents by the '--build' command
cat Related_cluster_families.tab | awk 'BEGIN{FS=OFS="\t"}{if($3=="FS" && $7!="-" && $11!="-" && $7!=$11 && $8!=$12)print}'



tail -n+2 AxiomGT1v2.comp_merge.deDup.sexConfirm.public_ids.kingRelatedOutput.kin0 | awk 'BEGIN{FS=OFS="\t"}{if($14=="PO"){print $1,$2;print $3,$4;}}' | sort | uniq > Related_PO.lst
#tail -n+2 AxiomGT1v2.comp_merge.deDup.sexConfirm.public_ids.kingRelatedOutput.kin0 | awk 'BEGIN{FS=OFS="\t"}{if($14!="2nd"){print $1,$2;print $3,$4;}}' | sort | uniq > Related_1st.lst

# Prepare the Plink input files: change chr format & make each sample in a seperate family
plink --bfile ${INPUT_FILE}.king --chr-set 38 no-xy --allow-extra-chr '0' \
      --keep Related_PO.lst --maf 0.01 \
      --recode --out ${INPUT_FILE}.PO.toSequoia
      
awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$3]=$4;next}{if(a[$1])print a[$1],$3}' tmp_map/map_id.tab <(tail -n+2 phenotypes/dog_profile.csv | tr ',' '\t') | cut -d"-" -f1 > sequoia.birthYear
cat ${INPUT_FILE}.king.fam | awk 'BEGIN{OFS="\t"}{if($5==1)print $2,2;else if($5==2)print $2,1;}' > sequoia.gender
echo ID Sex BirthYear | tr ' ' '\t' > sequoia.LHD
awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1]=$0;next}{if(a[$1])print a[$1],$2}' sequoia.gender sequoia.birthYear >> sequoia.LHD


> R
install.packages("sequoia")
library(sequoia)
input="AxiomGT1v2.comp_merge.deDup.sexConfirm.public_ids.PO.toSequoia.ped"
GenoM_name="AxiomGT1v2.comp_merge.deDup.sexConfirm.public_ids.po.sequoia"
GenoConvert(InFile = input, InFormat = "ped", OutFile = GenoM_name, OutFormat = "seq")
GenoM <- read.table(GenoM_name, sep = " ", header = FALSE, row.names = 1)
GenoM <- as.matrix(GenoM)
lhd=read.table("sequoia.LHD", header=TRUE)
SeqOUT <- sequoia(GenoM = GenoM,LifeHistData = lhd,Err = 0.01)



# One PO confimred in a cluster
grep -w KING47 Related_cluster_short.tab | awk 'BEGIN{OFS="\t"}{print $1,$1;print $2,$2;}' > Related_KING47.lst
plink --bfile ${INPUT_FILE}.king --chr-set 38 no-xy --allow-extra-chr '0' \
      --keep Related_KING47.lst --maf 0.05 --snps-only \
      --recode --out ${INPUT_FILE}.PO_KING47.toSequoia


> R
install.packages("sequoia")
library(sequoia)
input="AxiomGT1v2.comp_merge.deDup.sexConfirm.public_ids.PO_KING47.toSequoia.ped"
GenoM_name="AxiomGT1v2.comp_merge.deDup.sexConfirm.public_ids.PO_KING47.sequoia"
GenoConvert(InFile = input, InFormat = "ped", OutFile = GenoM_name, OutFormat = "seq")
GenoM <- read.table(GenoM_name, sep = " ", header = FALSE, row.names = 1)
GenoM <- as.matrix(GenoM)
SeqOUT <- sequoia(GenoM = GenoM,Err = 0.01)

#
cat Related_cluster_families.tab | awk 'BEGIN{FS=OFS="\t"}{if($3=="FS" && $7!="-" && $7 !~ /^[0-9]+$/ && $8!="-" && $8 !~ /^[0-9]+$/ && $7==$11 && $8==$12)print}' | grep FAM24 | awk 'BEGIN{OFS="\t"}{print $1,$1;print $2,$2;print $7,$7;print $8,$8;}' | sort | uniq > Related_FAM24.lst
plink --bfile ${INPUT_FILE}.king --chr-set 38 no-xy --allow-extra-chr '0' \
      --keep Related_FAM24.lst --geno 0 --maf 0.05 --snps-only \
      --recode --out ${INPUT_FILE}.FS_FAM24.toSequoia

> R
install.packages("sequoia")
library(sequoia)
input="AxiomGT1v2.comp_merge.deDup.sexConfirm.public_ids.FS_FAM24.toSequoia.ped"
GenoM_name="AxiomGT1v2.comp_merge.deDup.sexConfirm.public_ids.FS_FAM24.sequoia"
GenoConvert(InFile = input, InFormat = "ped", OutFile = GenoM_name, OutFormat = "seq")
GenoM <- read.table(GenoM_name, sep = " ", header = FALSE, row.names = 1)
GenoM <- as.matrix(GenoM)
SeqOUT <- sequoia(GenoM = GenoM,Err = 0.005)


####
cat AxiomGT1v2.comp_merge.deDup.sexConfirm.public_ids.PO.toSequoia.ped | cut -d" " -f1-100 > temp.ped  ## ie. 47 markers
head -n194 AxiomGT1v2.comp_merge.deDup.sexConfirm.public_ids.PO.toSequoia.map > temp.map
> R
library(sequoia)
input="temp.ped"
GenoConvert(InFile = input, InFormat = "ped", OutFile = "temp.sequoia", OutFormat = "seq")


