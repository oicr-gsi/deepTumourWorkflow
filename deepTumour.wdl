version 1.0


workflow deepTumour {
    input {
        File inputVcf
        File inputVcfIndex
        String outputFileNamePrefix     
        String reference
    }

    parameter_meta {
        inputVcf: "The input vcf file"
        inputVcfIndex: "index of input vcf"
        outputFileNamePrefix: "Prefix for output files"
        reference: "The genome reference build. For example: hg19, hg38"
    }

    call filterVcf {
        input:
        vcf_file = inputVcf,
        vcf_index = inputVcfIndex
    }
    
    call runDeepTumour {
        input:
        vcf = filterVcf.filtered_vcf,
        outputFileNamePrefix = outputFileNamePrefix,
        reference_genome = reference,
        modules = "deep-tumour/3.0.5.1 hg19/p13 bcftools/1.9"
    }

    meta {
        author: "Gavin Peng"
        email: "gpeng@oicr.on.ca"
        description: "The DeepTumour algorithm predicts the tissue of origin of a tumour based on the pattern of passenger mutations identified by Whole Genome Sequencing (WGS)."
        dependencies: [
            {
                name: "deep-tumour/3.0.5",
                url: "https://github.com/LincolnSteinLab/DeepTumour"
            }
        ]
      output_meta: {
        deepTumourOutputJson: {
            description: "the output json assigns a match probability from 0.0 to 1.0 for each of the 29 tumour types on which it was trained and chooses the tumour type with the highest probability score. The algorithm also calculates a type of confidence score based on the probability scores' distributione. A low entropy (< 2.0) is considered a confident score. HIgher values are unreliable (but might be correct).",
            vidarr_label: "deepTumourOutputJson"
        },
        filteredVcFile: {
            description: "the filtered vcf file as input of deepTumour, provision out for inspection",
            vidarr_label: "filteredVcFile"
        }
      }
    }
    output {
        File deepTumourOutputJson = runDeepTumour.outputJson
    }
}

task filterVcf {
    input {
        File vcf_file
        File vcf_index
        File repeat_bed = "/.mounts/labs/gsiprojects/gsi/gsiusers/gpeng/workflow/deepTumour/test/bed_files/hg38.repeat_regions.merged.bgzip.bed.gz"
        File repeat_bed_idx = "/.mounts/labs/gsiprojects/gsi/gsiusers/gpeng/workflow/deepTumour/test/bed_files/hg38.repeat_regions.merged.bgzip.bed.gz.tbi"
        Float t_vaf = 0.15
        Float n_vaf = 0.03
        Int mmq_threshold = 40
        Int cluster_window = 10
        Int cluster_count = 3
        Int indel_proximity = 5
        String modules = "bcftools/1.9 python/3.10.6"
        Int jobMemory = 24
        Int timeout = 4
    }
    parameter_meta {
        vcf_file:        "Input Mutect2 VCF (hg38, gzipped)"
        vcf_index:       "Index of input VCF"
        repeat_bed:      "Merged repeat regions BED file (bgzipped, hg38)"
        repeat_bed_idx:  "Tabix index for repeat_bed"
        t_vaf:           "Minimum tumor VAF (default 0.15)"
        n_vaf:           "Maximum normal VAF — variants above this are excluded (default 0.03)"
        mmq_threshold:   "Minimum alt allele median mapping quality (default 40)"
        cluster_window:  "Window size in bp for clustered SNV exclusion (default 10)"
        cluster_count:   "Min SNVs in window to trigger cluster exclusion (default 3)"
        indel_proximity: "Exclude SNVs within this many bp of an indel (default 5)"
        jobMemory:       "Memory allocated to job"
        modules:         "Required environment modules"
        timeout:         "Hours before task timeout"
    }

    command <<<
       set -euo pipefail

        echo "=== Step 1: bcftools filters ===" >&2

        # Extract sample names safely — grep -m1 causes SIGPIPE with set -o pipefail
        set +o pipefail
        TUMOR=$(zcat ~{vcf_file} | grep -m1 "^##tumor_sample"  | cut -d'=' -f2 | tr -d '\r')
        NORMAL=$(zcat ~{vcf_file} | grep -m1 "^##normal_sample" | cut -d'=' -f2 | tr -d '\r')
        set -o pipefail

        echo "Tumor:  $TUMOR" >&2
        echo "Normal: $NORMAL" >&2

        if [ -z "$TUMOR" ] || [ -z "$NORMAL" ]; then
            echo "ERROR: Could not extract tumor/normal sample names from VCF header" >&2
            exit 1
        fi

        bcftools view -f PASS ~{vcf_file} \
        | bcftools view -s "${TUMOR},${NORMAL}" \
        | bcftools filter \
            -i "FORMAT/AF[0:0] >= ~{t_vaf}
                && FORMAT/AF[1:0] < ~{n_vaf}
                && INFO/MMQ[1] >= 40" \
        | bcftools view -T ^~{repeat_bed} \
        -Oz -o prefiltered.vcf.gz

        bcftools index -t prefiltered.vcf.gz
        
        echo "After bcftools filters:" >&2
        bcftools stats prefiltered.vcf.gz | grep "^SN" >&2

        echo "=== Step 2: clustered SNV + indel proximity filter ===" >&2

        # Extract SNP positions and indel positions separately
        bcftools view -v snps prefiltered.vcf.gz \
          | bcftools query -f '%CHROM\t%POS\n' > snp_positions.txt
        bcftools view -v indels prefiltered.vcf.gz \
          | bcftools query -f '%CHROM\t%POS\n' > indel_positions.txt

        echo "SNPs before clustering/indel filter: $(wc -l < snp_positions.txt)" >&2
        echo "Indels (for proximity filter): $(wc -l < indel_positions.txt)" >&2

        python3 <<CODE
import sys
from collections import defaultdict

clust_win  = ~{cluster_window}
clust_cnt  = ~{cluster_count}
indel_prox = ~{indel_proximity}

# Load SNP positions
snp_positions = defaultdict(list)
with open("snp_positions.txt") as f:
    for line in f:
        chrom, pos = line.strip().split("\t")
        snp_positions[chrom].append(int(pos))

# Load indel positions
indel_positions = defaultdict(list)
with open("indel_positions.txt") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) == 2:
            indel_positions[parts[0]].append(int(parts[1]))

# Clustered SNV filter
def is_clustered(chrom, pos):
    return sum(1 for p in snp_positions[chrom] if abs(p - pos) <= clust_win) >= clust_cnt

# Indel proximity filter
def near_indel(chrom, pos):
    return any(abs(pos - ipos) <= indel_prox for ipos in indel_positions.get(chrom, []))

n_cluster = 0
n_indel   = 0
kept = []
for chrom, positions in snp_positions.items():
    for pos in positions:
        if is_clustered(chrom, pos):
            n_cluster += 1
        elif near_indel(chrom, pos):
            n_indel += 1
        else:
            kept.append((chrom, pos))

print(f"Fail clustered SNV:       {n_cluster}", file=sys.stderr)
print(f"Fail indel proximity:     {n_indel}", file=sys.stderr)
print(f"Pass all filters:         {len(kept)}", file=sys.stderr)

# Write passing positions as BED (0-based)
with open("filter.bed", "w") as out:
    for chrom, pos in kept:
        out.write(f"{chrom}\t{pos-1}\t{pos}\n")
CODE

        echo "=== Step 3: extract final VCF ===" >&2

        # Extract sample names for correct column ordering
        set +o pipefail
        zcat ~{vcf_file} | grep -m1 "^##tumor_sample"  | cut -d '=' -f2 >  samples.txt
        zcat ~{vcf_file} | grep -m1 "^##normal_sample" | cut -d '=' -f2 >> samples.txt
        set -o pipefail
        echo "Samples:" >&2
        cat samples.txt >&2

        bcftools view -S samples.txt -R "$PWD/filter.bed" prefiltered.vcf.gz \
          -Oz -o filtered.vcf.gz
        bcftools index -t filtered.vcf.gz

        echo "=== Final output ===" >&2
        bcftools stats filtered.vcf.gz | grep "^SN" >&2

    >>>

    runtime {
        memory:  "~{jobMemory} GB"
        modules: "~{modules}"
        timeout: "~{timeout}"
    }

    output {
        File filtered_vcf     = "filtered.vcf.gz"
        File filtered_vcf_idx = "filtered.vcf.gz.tbi"
        File filter_bed       = "filter.bed"
    }
}

task runDeepTumour {
    input {
        File vcf
        String outputFileNamePrefix
        String reference_genome
        Int jobMemory = 16
        String modules
        Int timeout = 4
    }
    parameter_meta {
        vcf:  "Input vcf file"
        outputFileNamePrefix: "Prefix for output file"
        reference_genome: "the reference genome fasta"
        jobMemory: "Memory allocated indexing job"
        modules:   "Required environment modules"
        timeout:   "Hours before task timeout"    
    }
    String liftover = if reference_genome == "hg38" then "--hg38" else ""

    command <<<
        set -euo pipefail

        mkdir out
        source $DEEP_TUMOUR_ROOT/.venv/bin/activate
        python $DEEP_TUMOUR_ROOT/src/DeepTumour.py --vcfFile ~{vcf} --reference $HG19_ROOT/hg19_random.fa ~{liftover} --outDir out --keep_input
        mv out/predictions_DeepTumour.json ~{outputFileNamePrefix}.predictions_DeepTumour.json

    >>>

    runtime {
        memory: "~{jobMemory} GB"
        modules: "~{modules}"
        timeout: "~{timeout}"
    }

    output {
        File outputJson = "~{outputFileNamePrefix}.predictions_DeepTumour.json"
    }

    meta {
        output_meta: {
            outputJson: "the output json of run DeepTumour"
        }
    }       
}

