
##Basic info
seqkit fx2tab  --name  --gc v_otu/v_otu.fa > v_otu/gc.txt ##GC content
omic_type=""    #  omic_type=_trans
mkdir bam_file$omic_type
fastq_fold=/mnt/h/SJYT_metagenomics/clean_fastq  # fastq_fold=/mnt/y/metatrascriptome/clean_fastq
tmp_location=/mnt/i/tmp;mkdir  $tmp_location/bam_file$omic_type
for var in $( ls $fastq_fold/* -d| xargs -n 1 basename ) 
do  (
/mnt/g/stone_meta/software/bbmap/bbmap.sh ref=v_otu/v_otu.fa  \
in1=$fastq_fold/$var/${var}.R_1.fastq.gz \
in2=$fastq_fold/$var/${var}.R_2.fastq.gz  \
out=stdout.bam -Xmx70g minid=0.9 threads=28   | samtools view -F 0x4 -b - | samtools sort -@ 6  > $tmp_location/bam_file${omic_type}/${var}.bam )
done
mv $tmp_location/bam_file${omic_type}  ./
conda activate test
coverm contig --bam-files  bam_file${omic_type}/*.bam --output-file v_otu/coverage${omic_type}.tpm.stat -m tpm --trim-min 0.10 --trim-max 0.90 --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -t 4  



##Identification of AMGs
conda activate /home/jintao/miniconda3/envs/mudoger_env/dependencies/conda/envs/virsorter2_env 
#virsorter setup -d /mnt/y/database/virsorters -j 4
virsorter run --prep-for-dramv -w dramv_output/virsorter -i v_otu/v_otu.fa   -j 20   --min-score 0   --viral-gene-enrich-off   --provirus-off  --min-score 0 #   --viral-gene-required
#DRAM-setup.py prepare_databases --output_dir DRAM_data --threads 12
#rm -r dramv_output/annotation
conda activate DRAM
DRAM-v.py annotate -i dramv_output/virsorter/for-dramv/final-viral-combined-for-dramv.fa -v dramv_output/virsorter/for-dramv/viral-affi-contigs-for-dramv.tab -o dramv_output/annotation  --threads 24 
DRAM-v.py distill -i dramv_output/annotation/annotations.tsv -o  dramv_output/annotation/distilled



##Annotation of vOTU sequences
conda activate phabox
python /mnt/g/stone_meta/software/PhaBOX/main.py  --contigs  v_otu/v_otu.fa  --rootpth phabox/ --threads  12 --len 3000  --dbdir /mnt/g/stone_meta/software/PhaBOX/database/ --parampth  /mnt/g/stone_meta/software/PhaBOX/parameters/

conda activate phatyp
cp v_otu/v_otu.fa /mnt/g/stone_meta/software/PhaTYP/
cd /mnt/g/stone_meta/software/PhaTYP/
python preprocessing.py --contigs v_otu.fa  --prodigal ./prodigal-gv   --threads 16    --len  3000
python PhaTYP.py --out PhaTYP.final_prediction.csv
cd -
cp /mnt/g/stone_meta/software/PhaTYP/PhaTYP.final_prediction.csv phabox/out/PhaTYP.final_prediction.csv

conda activate phagcn2
cp v_otu/v_otu.fa /mnt/g/stone_meta/software/PhaGCN2.0/
cd /mnt/g/stone_meta/software/PhaGCN2.0/
python /mnt/g/stone_meta/software/PhaGCN2.0/run_Speed_up.py  --contigs v_otu.fa --len 3000
cd -
cp /mnt/g/stone_meta/software/PhaGCN2.0/final_prediction.csv phabox/out/PhaGCN2.final_prediction.csv

conda activate phamer
cp v_otu/v_otu.fa /mnt/g/stone_meta/software/PhaMer/
cd /mnt/g/stone_meta/software/PhaMer/
python PhaMer_preprocessing_gv.py --contigs v_otu.fa  --len 3000   --threads  16
python PhaMer.py --out phamer.final_prediction.csv 
cd -
cp /mnt/g/stone_meta/software/PhaMer/phamer.final_prediction.csv phabox/out/phamer.final_prediction.csv



conda activatevcontact2_env
folder_taxonomy=vcontact2
mkdir vcontact2
cd vcontact2
/mnt/g/stone_meta/software/PhaTYP/prodigal-gv  -i ../v_otu/v_otu.fa -o viral_genomes.genes -a viral_genomes.faa -p meta -q -a v_otu_orf.faa -d  v_otu_orf.fa -o v_otu_orf.gff -f gff 
vcontact2_gene2genome -p viral_genomes.faa -o viral_genomes_g2g.csv -s 'Prodigal-FAA'
vcontact2 -t 20 --raw-proteins viral_genomes.faa --rel-mode 'Diamond' --proteins-fp viral_genomes_g2g.csv --db 'ProkaryoticViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin /mnt/h/analysis/database/vcontact2_env/bin/cluster_one-1.0.jar --output-dir vcontact2





##DS
mkdir defense-finder
conda activate defensefinder
find dereplicated_genomes_fasta -name '*.fasta' | xargs -n 1 basename | parallel -j12 "
mkdir -p defense-finder/{.}
prodigal -i dereplicated_genomes_fasta/{.}.fasta -a dereplicated_genomes_fasta/{.}.faa   -p meta 
defense-finder run dereplicated_genomes_fasta/{.}.faa --workers 1 -o defense-finder/{.}   "  ##can use prokka


##ADS
cat  genes.faa >  /mnt/h/analysis/database/dbAPIS/your_sequence.faa
cd /mnt/h/analysis/database/dbAPIS/
hmmscan --domtblout hmmscan.out --noali  --cpu 16  /mnt/h/analysis/database/dbAPIS/dbAPIS.hmm your_sequence.faa  2>/dev/null
#diamond makedb --in anti_defense.pep -d APIS_db 2>/dev/null
diamond blastp --db /mnt/h/analysis/database/dbAPIS/APIS_db -q your_sequence.faa  -p 12 -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen -o diamond.out --max-target-seqs 10000
bash parse_annotation_result.sh hmmscan.out diamond.out
cd -
cp /mnt/h/analysis/database/dbAPIS/diamond.out.parsed.tsv /mnt/h/analysis/database/dbAPIS/hmmscan.out.parsed.tsv ./v_otu
