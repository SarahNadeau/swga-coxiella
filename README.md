# SWGA Coxiella

This project is to try to reproduce the results in Cocking et al. 2020 "[Selective whole genome amplification and sequencing of Coxiella burnetii directly from environmental samples](https://doi.org/10.1016/j.ygeno.2019.10.022)" using the [swga](https://github.com/eclarke/swga) program.

### Download swga workflow
```
git clone git@github.com:sarahnadeau/wf-swga.git
```

### Get data
The mitochondrial sequences to exclude are included in this repo because a) they're small and b) I couldn't find a link to them.
```
git clone git@github.com:sarahnadeau/swga-coxiella.git
DATADIR=$PWD/swga-coxiella/data

# Background genomes
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/317/765/GCA_000317765.1_CHIR_1.0/GCA_000317765.1_CHIR_1.0_genomic.fna.gz -O $DATADIR/background_goat.fasta.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/003/055/GCA_000003055.5_Bos_taurus_UMD_3.1.1/GCA_000003055.5_Bos_taurus_UMD_3.1.1_genomic.fna.gz -O $DATADIR/background_cow.fasta.gz
gunzip -c $DATADIR/background*.fasta.gz >> $DATADIR/background_concat.fasta
gzip $DATADIR/background_concat.fasta
rm $DATADIR/background_goat.fasta.gz $DATADIR/background_cow.fasta.gz

# Target genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/017/105/GCF_000017105.1_ASM1710v1/GCF_000017105.1_ASM1710v1_genomic.fna.gz -O $DATADIR/target_coxiella.fasta.gz

# Exclude genomes
# Use file already in the downloaded git rep, which was created with:
# gunzip -c $DATADIR/exclude*.fasta.gz >> $DATADIR/exclude_concat.fasta
# gzip $DATADIR/exclude_concat.fasta
# rm $DATADIR/exclude_goat.fasta.gz $DATADIR/exclude_cow.fasta.gz
```

### Run workflow (with down-sampling)
```
CHUNK_SIZE=80000
TARGET_N_CHUNKS=3
BACKGR_N_CHUNKS=3

TARGET_LEN=$(( $CHUNK_SIZE * $TARGET_N_CHUNKS ))
BACKGR_LEN=$(( $CHUNK_SIZE * $BACKGR_N_CHUNKS ))

# min/max k-mer sizes comes from min/max of published primers
# max melting temp comes from published methods
# min_fg/max_bg bind rates comes from published methods
# find_sets_max_size based on published set of 40 primers plus a buffer of 10 to see what happens
# max_sets_search and chunk #s specified to make the analysis run quicker 
nextflow run \
    -profile docker wf-swga/main.nf \
    --outpath swga-coxiella/swga_results \
    --target $DATADIR/target_coxiella.fasta.gz \
    --background $DATADIR/background_concat.fasta.gz \
    --exclude $DATADIR/exclude_concat.fasta.gz \
    --target_length $TARGET_LEN \
    --backgr_length $BACKGR_LEN \
    --min_kmer_size 9 \
    --max_kmer_size 10 \
    --max_tm 30 \
    --min_bg_bind_dist 80000 \
    --max_fg_bind_dist 7000 \
    --find_sets_max_size 50 \
    --max_sets_search 1000 \
    --target-chunk-size $CHUNK_SIZE \
	--backgr-chunk-size $CHUNK_SIZE \
	--target-n-chunks $TARGET_N_CHUNKS \
	--backgr-n-chunks $BACKGR_N_CHUNKS
```
The analysis hung after filtering, exporting primers and before finding, exporting sets.
Re-run manually to see what went wrong.
```
# Need to get target genome used in order to calculate primer distances while finding sets
cp work/b6/b3bf0ece44f846187c1030b73fae09/target.fasta.noblanks.downsampled swga-coxiella/swga_results/swga

# Run swga interactively in docker container
docker run \
    --rm \
    -it \
    --mount type=bind,src=$PWD/swga-coxiella/swga_results/swga,dst=/data \
    snads/swga@sha256:776a2988b0ba727efe0b5c1420242c0309cd8e82bff67e9acf98215bf9f1f418
    
# Put target genome where swga expects it based on nextflow run
WORKDIR=/Users/nadeau/Documents/CDC_ORISE/Projects/work/b6/b3bf0ece44f846187c1030b73fae09
mkdir -p $WORKDIR
mv target.fasta.noblanks.downsampled $WORKDIR

swga find_sets \
    --workers 1 \
    --max_sets 1000 \
    --min_bg_bind_dist 80000 \
    --max_fg_bind_dist 7000 \
    --max_size 50
```
Okay, problem is no sets pass the filters for at least 40min.
Trying with more lenient filters for binding distance didn't help.
Looking at primer output, all exported bound only 2 times in ~240000bp target, 0 times in background --> infinite ratio.
Need to re-run with more background sequence and/or more stringent filter criteria for target binding rate.
Try min target binding rate equivalent to max target bind distance (1/7000 = 0.000142)
```
nextflow run \
    -profile docker wf-swga/main.nf \
    --outpath swga-coxiella/swga_results_min_fg_bind_rate_0.0001 \
    --target $DATADIR/target_coxiella.fasta.gz \
    --background $DATADIR/background_concat.fasta.gz \
    --exclude $DATADIR/exclude_concat.fasta.gz \
    --target_length $TARGET_LEN \
    --backgr_length $BACKGR_LEN \
    --min_kmer_size 9 \
    --max_kmer_size 10 \
    --max_tm 30 \
    --min_bg_bind_dist 80000 \
    --max_fg_bind_dist 7000 \
    --target-chunk-size $CHUNK_SIZE \
	--backgr-chunk-size $CHUNK_SIZE \
	--target-n-chunks $TARGET_N_CHUNKS \
	--backgr-n-chunks $BACKGR_N_CHUNKS \
	--min_fg_bind_rate 0.0001 \
	--run_find_sets "false"
	
cat swga-coxiella/swga_results_min_fg_bind_rate_0.0001/.log/stderr.nextflow.txt
# 52/140256 primers bind the foreground genome >= 24 times
# 52/52 primers bind the background genome <= 2 times
# 52 primers satisfy all filters so far.
# Finding melting temps for 52 primers...
# 14/52 primers have a melting temp between 15.0 and 30.0 C
# Finding binding locations for 14 primers...
# Finding Gini coefficients for 14 primers...
# 12/14 primers have a Gini coefficient <= 0.6
# Marked 12 primers as active.
```
Are all these primers in the final set they used?
Software exported many more than 12 primers. How to get only active ones?
"When you're happy with the primers you've selected (remember, you can check the statistics with swga summary)"

### Run workflow (without down-sampling)
This was run on the Aspen cluster, not locally.
```
qlogin

# Downloaded data as above

# Get genome lengths
TARGET_LEN=$(gunzip -c $DATADIR/target_coxiella.fasta.gz | wc -c | awk '{print $1}')
BACKGR_LEN=$(gunzip -c $DATADIR/background_concat.fasta.gz | wc -c | awk '{print $1}')

nextflow run \
    -profile singularity wf-swga/main.nf \
    --outpath swga-coxiella/swga_results_22-03-16_no_sets \
    --target $DATADIR/target_coxiella.fasta.gz \
    --background $DATADIR/background_concat.fasta.gz \
    --exclude $DATADIR/exclude_concat.fasta.gz \
    --target_length $TARGET_LEN \
    --backgr_length $BACKGR_LEN \
    --min_kmer_size 9 \
    --max_kmer_size 10 \
    --max_tm 30 \
    --min_bg_bind_dist 80000 \
    --max_fg_bind_dist 7000 \
    --target-chunk-size $CHUNK_SIZE \
	--backgr-chunk-size $CHUNK_SIZE \
	--target-n-chunks $TARGET_N_CHUNKS \
	--backgr-n-chunks $BACKGR_N_CHUNKS \
	--min_fg_bind_rate 0.0001 \
	--run_find_sets "false"
```