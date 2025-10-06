# deepTumour

The DeepTumour algorithm predicts the tissue of origin of a tumour based on the pattern of passenger mutations identified by Whole Genome Sequencing (WGS).

## Overview

## Dependencies

* [deep-tumour 3.0.5](https://github.com/LincolnSteinLab/DeepTumour)


## Usage

### Cromwell
```
java -jar cromwell.jar run deepTumour.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`inputVcf`|File|The input vcf file
`inputVcfIndex`|File|index of input vcf
`outputFileNamePrefix`|String|Prefix for output files
`reference`|String|The genome reference build. For example: hg19, hg38
`filter_method`|String|method to filter input vcf directly or filter using maf, value can only be either 'maf' or 'vcf'


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`inputMaf`|File?|None|the input maf file


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`filterMaf.t_depth`|Int|1|tumour depth filter threshold
`filterMaf.t_vaf`|Float|0.01|Tumor Variant Allele Frequency threshold
`filterMaf.gnomad_af`|Float|0.001|gnomAD allele frequency threshold
`filterMaf.valid_exonic`|Array[String]|["5'Flank", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Silent", "Splice_Region", "Splice_Site", "Targeted_Region", "Translation_Start_Site"]|list of exonic variants to keep
`filterMaf.exclude_mutations`|Array[String]|["str_contraction", "t_lod_fstar"]|list of exonic variants to exclude
`filterMaf.modules`|String|"pandas/2.1.3 bcftools/1.9"|Required environment modules
`filterMaf.jobMemory`|Int|24|Memory allocated indexing job
`filterMaf.timeout`|Int|2|Hours before task timeout
`filterVcf.t_depth`|Int|1|tumour depth filter threshold
`filterVcf.t_vaf`|Float|0.01|Tumor Variant Allele Frequency threshold
`filterVcf.modules`|String|"bcftools/1.9"|Required environment modules
`filterVcf.jobMemory`|Int|8|Memory allocated indexing job
`filterVcf.timeout`|Int|1|Hours before task timeout
`runDeepTumour.jobMemory`|Int|16|Memory allocated indexing job
`runDeepTumour.timeout`|Int|4|Hours before task timeout


### Outputs

Output | Type | Description | Labels
---|---|---|---
`deepTumourOutputJson`|File|the output json assigns a match probability from 0.0 to 1.0 for each of the 29 tumour types on which it was trained and chooses the tumour type with the highest probability score. The algorithm also calculates a type of confidence score based on the probability scores' distributione. A low entropy (< 2.0) is considered a confident score. HIgher values are unreliable (but might be correct).|vidarr_label: deepTumourOutputJson


## Commands
 This section lists command(s) run by deepTumour workflow
 
 * Running deepTumour
 
 ```
         # --- Step 1: filter MAF with criteria ---
         python3<<CODE
         import pandas as pd
         import os
 
         df = pd.read_csv("~{maf_file}", sep="\t", comment="#", low_memory=False)
 
         def safe_float(x):
             try: return float(x)
             except: return None
 
         df["gnomAD_AF"] = df["gnomAD_AF"].apply(safe_float)
         valid_exonic_list = "~{sep=',' valid_exonic}".split(",")
         exclude_mutations_list = "~{sep=',' exclude_mutations}".split(",")
 
         filtered = df[
             (df["gnomAD_AF"].notna()) &
             (df["gnomAD_AF"] < ~{gnomad_af}) &
             (df["Variant_Classification"].isin(valid_exonic_list)) &
             (~df["Variant_Classification"].isin(exclude_mutations_list)) &
             (df["t_depth"] > ~{t_depth}) &
             (df["t_alt_count"] / df["t_depth"] > ~{t_vaf})
         ]
 
         bed_df = filtered[["Chromosome","Start_Position","End_Position"]].copy()
         bed_df["Start_Position"] = bed_df["Start_Position"] - 1
 
         
         bed_path = os.path.join(os.getcwd(), "filter.bed")
         bed_df.to_csv(bed_path, sep="\t", index=False, header=False)
 
         CODE
 
         # --- Step 2: apply BED filter to original VCF ---
         # Extract tumour/normal names
         grep "^##tumor_sample" ~{vcf_file} | cut -d '=' -f2 > samples.txt
         grep "^##normal_sample" ~{vcf_file} | cut -d '=' -f2 >> samples.txt
 
         # Subset VCF to variants in BED + reorder samples
         bcftools view -f PASS -S samples.txt -R "$PWD/filter.bed" ~{vcf_file} -Oz -o filtered.vcf.gz
 
 ```
 ```
         set -euo pipefail
 
         # extract tumour/normal names from header
         zcat ~{vcf_file}|grep "^##tumor_sample"|cut -d '=' -f2 > samples.txt
         zcat ~{vcf_file}|grep "^##normal_sample"|cut -d '=' -f2 >> samples.txt
 
         # Filter chain:
         #   1. Keep PASS only
         #   2. Reorder samples to Normal, Tumour
         #   3. Filter on tumour depth
         #   4. Filter on tumour VAF
         bcftools view -f PASS -S samples.txt ~{vcf_file} -Ou \
           | bcftools filter -i "(FORMAT/DP[1]) >= ~{t_depth}" -Ou \
           | bcftools filter -i "(FORMAT/AD[1:1])/(FORMAT/DP[1]) >= ~{t_vaf}" -Oz -o filtered.vcf.gz
 ```
 ```
         set -euo pipefail
 
         mkdir out
         source $DEEP_TUMOUR_ROOT/.venv/bin/activate
         python $DEEP_TUMOUR_ROOT/src/DeepTumour.py --vcfFile ~{vcf} --reference $HG19_ROOT/hg19_random.fa ~{liftover} --outDir out --keep_input
         mv out/predictions_DeepTumour.json ~{outputFileNamePrefix}.predictions_DeepTumour.json
 
 ```
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
