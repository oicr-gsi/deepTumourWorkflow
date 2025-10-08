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
        python3 <<CODE
        import csv
        import gzip

        output_bed = "filter.bed"
        open_func = gzip.open if "~{maf_file}".endswith(".gz") else open
        with open_func("~{maf_file}", "rt") as maf, open(output_bed, "w") as out:
            reader = csv.DictReader((l for l in maf if not l.startswith("#")), delimiter="\t")

            for row in reader:
                try:
                    t_depth = int(row.get("t_depth", 0))
                    t_alt = int(row.get("t_alt_count", 0))
                    af = float(row.get("gnomAD_AF", "0") or "0")
                except ValueError:
                    continue

                if (
                    af < ~{gnomad_af}
                    and row["Variant_Classification"] in "~{sep=',' valid_exonic}"
                    and row["Variant_Classification"] not in "~{sep=',' exclude_mutations}"
                    and t_depth > ~{t_depth}
                    and (t_alt / t_depth) > ~{t_vaf}
                ):
                    chrom = row["Chromosome"]
                    start = int(row["Start_Position"]) - 1
                    end = int(row["End_Position"])
                    out.write(f"{chrom}\t{start}\t{end}\n")
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
