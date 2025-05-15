# ### modified pipe
#output_folder=$PWD/viruses/investigation
#viruses_folder=$PWD/viruses	
#assembly=$PWD/assembly/final_assembly.fasta


conda activate virsorter2_env
virsorter run all -i "$assembly" -w virsorter -j 28  --include-groups "dsDNAphage,ssDNA"

conda activate dvf_env
python /mnt/g/stone_meta/software/DeepVirFinder/dvf.py -i "$assembly"  -o dvf -l 5000 -c 6


conda activate dmc_env
DeepMicroClass predict -i "$assembly"   -o dmc 
/mnt/d/Myfile/Research/Script/MetagenomicsScripts/process_deepmicroclass_file.sh   dmc/final_assembly.fasta_pred_one-hot_hybrid.tsv > dmc/class.tsv



#cat vibrant/VIBRANT_final_assembly/VIBRANT_phages_final_assembly/final_assembly.phages_combined.fna | grep ">" | sed "s/_fragment_1//g;s/>//g"   > vibrant_filtered_data.txt
cat virsorter/final-viral-combined.fa  | grep ">" | sed "s/_fragment_1//g;s/>//g" |  cut -f1 -d "|" > virsorter2_filtered_data.txt
cat dvf/final_assembly.fasta_gt5000bp_dvfpred.txt | awk -F'\t' '{ if ( $4 <= 0.05) print }'| cut -f1  > dvf_filtered_data.txt 
cat dmc/class.tsv | grep "Virus" |cut -f1 |sed '1d'   > dmc_filtered_data.txt



mkdir -p dereplication
cat *.txt | sort -u >  dereplication/viral_unique_contigs
##extract viral contig from ../../assembly/final_assembly.fasta  to dereplication/uvigs_1.fa


# ###Dereprelication
cd  dereplication
cat uvigs_*.fa  | seqtk seq  -L 5000 | sed  's/___contigs/__ctgs/' > uvigs.fa 

##remove low quality
export CHECKVDB=/mnt/y/database/assembly/checkv/checkv-db-v1.0
conda activate "checkv_env
checkv end_to_end uvigs.fa vcheck_quality -t 28


grep -w "Not-determined" vcheck_quality/quality_summary.tsv | cut  -f 1 > bad_quality_contigs.txt
#grep -w "Low-quality" vcheck_quality/quality_summary.tsv | cut  -f 1 >> bad_quality_contigs.txt
awk '$10 == "" || $10 < 40 {print $1}'  vcheck_quality/quality_summary.tsv  >> bad_quality_contigs.txt

seqkit grep -v -f bad_quality_contigs.txt uvigs.fa > uvigs-filter.fa
seqkit stat uvigs-filter.fa

mkdir ../../v_otu
#cp uvigs-filter.fa ../../v_otu/v_otu.fa

mmseqs easy-cluster uvigs-filter.fa --min-seq-id 0.95 -c 0.85  uvigs-cluster.fa tmp --threads 16 --cluster-mode 2
seqkit stat uvigs-cluster.fa_rep_seq.fasta
sed 's/ //g'  -i uvigs-cluster.fa_rep_seq.fasta
cp uvigs-cluster.fa_rep_seq.fasta ../../v_otu/v_otu.fa

cd ../..


