# Calculating polygenic scores from Nebula Genomics WGS data

## What is this?

This repository contains the commands and important data files that I used for the analysis presented in
[this blog post](https://mfasold.net/blog/calculating-polygenic-risk-scores-for-wgs/).

The following files are provided:

 - the actual reference genome used by Nebula Genomics (download and `gunzip` [this file](https://mfasold.net/download/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz) and the `.tbi` file the subfolder `reference-genome`)
 - The "alleles VCF" file containing the locations and alleles for all DNA variants used in the
   PGSKB CLI tool (as of 2024-05) in the subfolder `PRSKB`
   
The commands to calculate PRS are provided directly in this README.

*Disclaimer*: This is a proof-of-concept analysis that disregards several cases (InDels) and
contains hacky one-liners without any error checking.

## Prerequisites

### Required input files
You will need to download the following raw data files from Nebula Genomics. 

1. The DNA variants file. For me the file was named like
   `NG123456M.mm2.sortdup.bqsr.hc.vcf.gz`. Let's rename this file to `me.vcf.gz`.
2. The DNA variants index file. For me the file was named like
   `NG123456M.mm2.sortdup.bqsr.hc.vcf.gz.tbi`. Let's rename this file to `me.vcf.gz.tbi`.
3. The alignment file. For me this file was names like `NG123456M.mm2.sortdup.bqsr.cram`. Let's
   rename this file to `me.cram`

### Required additional files

It will be very useful if we can annotate our DNA variants with DBsnp identifiers, for wich we need 
a file with all human SNPs. However, the dbSNP database is growing very fast and the files get very big.
I therefore use the following slightly old and preprocessed file, with a size of 3.2GB.

```
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz.tbi
```

### System setup
All the commands work in Linux only (for example, GNU Awk is used) and may need 32 GB memory or more.

I suggest to use the conda package manager for installing required tools in the specific
versions. Please follow the instructions [here](https://docs.anaconda.com/free/miniconda/#quick-command-line-install). 

We need access to bioconda:
```sh
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Once you have the conda environment activated, install the variant calling tool GATK4. A specific version is needed to
work with CRAM files of Nebula Genomics (see [this issue](https://gatk.broadinstitute.org/hc/en-us/community/posts/11440622639387-Unable-to-trim-uncertain-bases-without-flow-order-information)).

```sh
conda install gatk4=4.2.6.1
```

The htslib package is needed for tools like `bgzip` and `tabix`.
```sh
conda install htslib
```

Furthermore, we install some python tools with `pip`. This is best done in a dedicated VirtualEnv

```
pip install virtualenv
python -m virtualenv venv
source ./venv/bin/activate
```

At this point, your command prompt should begin with `(venv) (base)` that shows you that both environments are active.

## Calculating PRS with the PRSKB CLI

### Preparation

Download the PRSKB CLI tool from [here](https://prs.byu.edu/cli_download.html) and install the required dependencies

```
pip install "setuptools<58.0"
pip install filelock requests pyvcf myvariant Bio
```

Then comment out [this line](https://github.com/kauwelab/PolyRiskScore/blob/b41e3538ab3886f0ad173ee9ad8141241560c46a/static/downloadables/connect_to_server.py#L897) from the file `connect_to_server.py`.

We also need the [jq tool](https://jqlang.github.io/jq/download/) to work with the JSON files. On Ubuntu/Debian it can be installed with `sudo apt-get install jq`.


### Create a VCF file with all variants in PRSKB

*Note:* If you just want to run the analysis on your variants, you can skip this step and just use the files `PRSKB_all_variants.dbsnp.vcf.gz` and `PRSKB_all_variants.dbsnp.vcf.gz.tbi` provided in the subfolder `PRSKB` of this repository. 

The easiest way to obtain the variants used by PRSKB CLI is to run it, as it will download the required data and keep it in the folder `.workingFiles`. We can simply use the VCF from Nebula Genomics for this "dummy run"

```
bash runPrsCLI.sh -f me.vcf.gz -o PRS.from_NG_vcf.tsv -r hg38 -p EUR -c 0.05
```
We will look at the parameters later.

The file `.workingFiles/allAssociations_hg38.txt` now contains all associations of SNPs with
study-trait-combinations within one very big JSON file. The `jq` tool allows us to extract the
required information from this JSON file and store it in the tab-separated text file
`allAllelesFromAssociations.tsv`. 

Next, we must turn these into valid VCF records suitable for variant calling. We must obtain the
correct REF alleles, which we do by joining above tsv-file with the big DBsnp table downloaded
earlier on rsIDs using `awk`. Effect alleles are stored in the INFO column of the VCF record. In the last step, we iterate over the VCF records, checking if those effect alleles are among the alleles for that rsID. If not, we add them.

```
cat .workingFiles/allAssociations_hg38.txt | jq -r '.associations | to_entries | .[] | .key as $rsid | .value | objects | .pos as $pos | .traits | .[] | .[] | .[] | [$pos, $rsid, (keys[])] | @tsv' | uniq | sort | uniq >allAllelesFromAssociations.tsv
awk -F'\t' -v OFS='\t' 'FNR==NR {if ($2 in d) d[$2] = d[$2]","$3; else d[$2]=$3; next;} ($3 in d) { print $1,$2,$3,$4,$5,1,".","A="d[$3]," "," "}' allAllelesFromAssociations.tsv <(zcat dbsnp_146.hg38.vcf.gz) >dbSnpWithAllelesFromAssoc.pre.vcf
# There are 1728 cases where the effect allele is not one of the alleles in dbSNP - we need to add them
cat dbSnpWithAllelesFromAssoc.pre.vcf | awk -F '\t' -v OFS='\t' '{i=$8; gsub("A=", "", i); split(i,a,","); for (e in a) if (a[e] != $4 && index($5, a[e]) == 0) $5=$5","a[e]; print $0;}' >dbSnpWithAllelesFromAssoc.pre.fixed.vcf
```

Lastly, we turn this collection of VCF records into a valid VCF file. A simple option is to just use the VCF header from Nebula Genomics' VCF file, and concatenate it with the records. We also compress and index the resulting VCF file.

```
bcftools view -h me.vcf.gz >PRSKB_all_variants.dbsnp.vcf
cat dbSnpWithAllelesFromAssoc.pre.vcf >>PRSKB_all_variants.dbsnp.vcf

bgzip PRSKB_all_variants.dbsnp.vcf
tabix PRSKB_all_variants.dbsnp.vcf.gz
```

### Calling DNA variants for effect alleles in PRSKB

The GATK 4 Haplotype caller can be used to call your variants based on the locations and alleles
from PRSKB. We use a special mode of GATK, activated with the `--alleles`, `-L` and `--output-mode`
options. Important additional input files for this tool are the aligned reads (me.cram, downloaded
from NG), the reference genome used by NG (in this repository in the `ReferenceFromNebula`
folder). For convenience, the rsIDs shall be added to output VCF using the `--dbsnp` option.

```
gatk --java-options '-Xms2g -Xmx4g'  HaplotypeCaller -R ReferenceFromNebula/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna -L PRSKB_all_variants.dbsnp.vcf.gz --alleles PRSKB_all_variants.dbsnp.vcf.gz -I me.cram --output-mode EMIT_VARIANTS_ONLY --dbsnp PRSKB_all_variants.dbsnp.vcf.gz -O my_wgs.PRSKB_all_variants.vcf.gz
```

### Run the PRSKB CLI tool

It is time to calculate the scores which can be done with just one command. We have to provide the
DNA variants calculated in the previous step. The verbose option `-v` is required to get the
percentiles. We need to specify the correct reference genome with `-r hg38`, the super-population
with `-p EUR` (adapt, if required) and the p-value treshold with `-c 0.05`. The final option `-h
0.1` defines the Imputation Threshold, an upper limit for how many SNPs are allowed to be imputed:
"We divide the number of imputed SNPs by the total number of SNPs in the calculation and if that
number exceeds the threshold we do not report that study." Setting it to 10% is a safety measure to
ensure that we considered enough SNPs per study. 

```
bash runPrsCLI.sh -f my_wgs.PRSKB_all_variants.vcf.gz -v -o PRS.PRSKB_all_variants.tsv -r hg38 -p EUR -c 0.05 -h 0.1
```

The output file `PRS.PRSKB_all_variants.tsv` is a big table that has the percentile in column 15.

## Calculating PGS with pgs_calc

### A note on alternative approaches
We could compute the PGS in a similar way as we did with PRSKB CLI. That is

1. Collect all locations and alleles of all relevant PGS variants in one big VCF file
2. Use GATK4 to call DNA variants on this locations
3. Transform the results variant VCF in a form suitable for pgs_calc (Hint: use `bcftools norm  -m -any` to get bi-allelic VCF records)

I tried this, but it turned out that there are several small complications like a more
difficult joining with the DBsnp database and more difficult transformations of the variant VCF. So
I choose an alternative approach that focuses on the locations of PGS variants only.

### Preparation

Install Docker has discribed [here](https://docs.docker.com/desktop/install/linux-install/).

Install Nextflow:

```
conda install nextflow
```

Download the reference database required to compute percentiles
```
wget https://ftp.ebi.ac.uk/pub/databases/spot/pgs/resources/pgsc_HGDP+1kGP_v1.tar.zst
```

### Creation of a combined locations file

We will do this analysis for the six PGS
PGS000662, PGS000192, PGS000043, PGS000210, PGS000882, PGS000904. I tried to do the analysis for *all* scores in
one go, but it did not work out. Please replace these PGS with whatever list of scores you are interested in.

```
mkdir scorefiles
cd scorefiles
pip install pgscatalog-utils

# The following unfortunatly does not seem to work - the resulting file has not all locations
# wget https://ftp.ebi.ac.uk/pub/databases/spot/pgs/pgs_scores_list.txt
# for pgsid in $(cat pgs_scores_list.txt); do echo $pgsid; pgscatalog-download -i $pgsid -o . -b GRCh38; done;

for pgsid in PGS000662 PGS000192 PGS000043 PGS000210 PGS000882 PGS000904; do echo $pgsid; pgscatalog-download -i $pgsid -o . -b GRCh38; done;
pgscatalog-combine -s PGS*.txt.gz -t GRCh38 -o all_scores_combined.txt.gz
cd ..
```

We can now collect a list of variant locations, where we will exclude some special cases (chromosomes not included in Nebula Genomics' reference genome, undefined loci).

```
zcat scorefiles/all_scores_combined.txt.gz | grep -v -e "chr" -e "##" |  awk -F'\t' '{hash = $1":"$2; dict[hash];} END { for (h in dict) print h;}' | grep -v -e "::" -e "MT" -e "alt" | sort -t':' -k1,1 -k2,2n | awk '{print "chr"$0;}' >scorefiles_locations.list
```

### Calling DNA variants on variant locations

The GATK HaplotypeCaller will be used to determine DNA variants on all locations from our list. In this case, two steps are required.

```
gatk --java-options '-Xms2g -Xmx4g'  HaplotypeCaller   -R ReferenceFromNebula/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna   -L scorefiles_locations.list -I me.cram -O my_wgs.scorefile_locations.vcf.gz -ERC BP_RESOLUTION 
gatk --java-options '-Xms2g -Xmx4g'  GenotypeGVCFs -R ReferenceFromNebula/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna -V my_wgs.scorefile_locations.vcf.gz -O my_wgs.scorefile_locations.gt.vcf.gz --dbsnp dbsnp_146.hg38.vcf.gz --include-non-variant-sites true
```

### Run pgsc_calc and prepare required input files 

We now have to transfom this VCF into a suitable form for `pgsc_calc` (which is geared towards input
data from arrays/imputation servers). In particular, this means adding the other effect alleles from
the combined scorefile for homyzgous reference variants because otherwise this very large fraction
of variants would not be matched. At this point we will cut some corners, like ignoring entries
where other_allele is not set (empty ALT is not accepted!). 

The general idea is to first read in the combined scorefile and store all locations together with
the effect allele and other allele in memory. Then go through our variant file, and if we find a
homref variant (ALT is not set) then add entries for all possible effect alleles and their other
alleles. A new record is created for each other allele, as `pgsc_calc` only supports bi-allelic
records. We also deal with the special case where GATK sets the QUAL to "Infinity", which is also
not supported.

```
awk -F'\t' -v OFS='\t' 'NR == FNR {h="chr"$1":"$2; if ($4) d[h][$3]=$4; next;} $6 == "Infinity" { $6 = 10000} $5 == "." {h=$1":"$2; if (h in d) for (ea in d[h]) { if ($4 != ea) { $5 = ea } else { $5 = d[h][ea] }; if (!seen[h":"$4":"$5]++) print; }; next;} {print;}' <(zcat scorefiles/all_scores_combined.txt.gz) <(zcat my_wgs.scorefile_locations.gt.vcf.gz) |  bgzip -c >my_wgs.scorefile_locations.homrefAllelesAdded.vcf.gz
tabix my_wgs.scorefile_locations.homrefAllelesAdded.vcf.gz
```

Now we prepare the folder structure and files required for `pgsc_calc` (as detailed in their [guide](https://pgsc-calc.readthedocs.io/en/latest/getting-started.html)

```
mkdir my-wgs
ln -s my_wgs.scorefile_locations.homrefAllelesAdded.vcf.gz my-wgs/me.vcf.gz
ln -s my_wgs.scorefile_locations.homrefAllelesAdded.vcf.gz.tbi my-wgs/me.vcf.gz.tbi

echo -e "sampleset,path_prefix,chrom,format\nmy-wgs,$PWD/my-wgs/me,,vcf">samplesheet.csv

nextflow run pgscatalog/pgsc_calc -profile docker --input samplesheet.csv --target_build GRCh38 --pgs_id PGS000662,PGS000192,PGS000043,PGS000210,PGS000882,PGS000904 --run_ancestry pgsc_HGDP+1kGP_v1.tar.zst --keep_ambiguous
```

After the analysis, the score will be in the folder `results/my-wgs/score/`. If you run the pipeline again with changed input, remember to delete the `work/` directory. 
