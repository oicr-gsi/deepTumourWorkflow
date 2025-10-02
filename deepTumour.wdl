version 1.0


workflow deepTumour {
    input {
        File inputVcf
        File inputVcfIndex
        File? inputMaf
        String outputFileNamePrefix     
        String reference
        String filter_method
        Array[String]? valid_exonic 
        Array[String]? exclude_mutations 
        Int? t_depth = 1
        Float? t_vaf = 0.01
        Float? gnomad_af = 0.001
    }

    parameter_meta {
        inputVcf: "The input vcf file"
        inputVcfIndex: "index of input vcf"
        outputFileNamePrefix: "Prefix for output files"
        reference: "The genome reference build. For example: hg19, hg38"
        filterVcf: "whether to filter input vcf for rows that has PASS value in FILTER column"
        valid_exonic: "list of exonic variants to keep"
        exclude_mutations: "list of exonic variants to exclude"
        t_depth: "tumour depth filter threshold"
        t_vaf: "Tumor Variant Allele Frequency threshold"
        gnomad_af: "gnomAD allele frequency threshold"
    }
    if (filter_method == "maf" && defined(inputMaf)) {
        call filterMaf {
            input:
            maf_file = select_first([inputMaf]),
            vcf_file = inputVcf,
            vcf_index = inputVcfIndex,
            t_depth = t_depth,
            t_vaf = t_vaf,
            gnomad_af = gnomad_af,
            valid_exonic = valid_exonic,
            exclude_mutations = exclude_mutations
        }
    }
    if (filter_method == "vcf") {
        call filterVcf {
            input:
            vcf_file = inputVcf
        }
    }
    File filteredVcf = select_first([filterMaf.vcf_filtered, filterVcf.vcf_filtered])
    
    call runDeepTumour {
        input:
        vcf = filteredVcf,
        outputFileNamePrefix = outputFileNamePrefix,
        reference_genome = reference,
        modules = "deep-tumour/3.0.4 hg19/p13 bcftools/1.9"
    }

    meta {
        author: "Gavin Peng"
        email: "gpeng@oicr.on.ca"
        description: "The DeepTumour algorithm predicts the tissue of origin of a tumour based on the pattern of passenger mutations identified by Whole Genome Sequencing (WGS)."
        dependencies: [
            {
                name: "deep-tumour/3.0.4",
                url: "https://github.com/LincolnSteinLab/DeepTumour"
            }
        ]
      output_meta: {
        deepTumourOutputJson: {
            description: "the output json assigns a match probability from 0.0 to 1.0 for each of the 29 tumour types on which it was trained and chooses the tumour type with the highest probability score. The algorithm also calculates a type of confidence score based on the probability scores' distributione. A low entropy (< 2.0) is considered a confident score. HIgher values are unreliable (but might be correct).",
            vidarr_label: "deepTumourOutputJson"
        }
      }
    }
    output {
        File deepTumourOutputJson = runDeepTumour.outputJson
    }
}

task filterMaf {
    input {
        File maf_file
        File vcf_file
        File vcf_index
        Int t_depth = 1
        Float t_vaf = 0.1
        Float gnomad_af = 0.001
        Array[String]? valid_exonic = ["5'Flank","Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins",
            "Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Silent",
            "Splice_Region","Splice_Site","Targeted_Region","Translation_Start_Site"]
        Array[String]? exclude_mutations = ["str_contraction", "t_lod_fstar"]
        String modules = "pandas/2.1.3 bcftools/1.9"
        Int jobMemory = 24
        Int timeout = 2
    }
    parameter_meta {
        maf_file:  "Input maf file"
        vcf_file:  "Input vcf file"
        vcf_index: "index of input vcf"
        t_depth: "tumour depth filter threshold"
        t_vaf: "Tumor Variant Allele Frequency threshold"
        gnomad_af: "gnomAD allele frequency threshold"
        jobMemory: "Memory allocated indexing job"
        modules:   "Required environment modules"
        timeout:   "Hours before task timeout"
    }

    command <<<
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

    >>>
    runtime {
        memory: "~{jobMemory} GB"
        modules: "~{modules}"
        timeout: "~{timeout}"
    }
    output {
        File vcf_filtered = "filtered.vcf.gz"
    }
}

task filterVcf {
    input {
        File vcf_file
        Int t_depth = 1
        Float t_vaf = 0.01
        String modules = "bcftools/1.9"
        Int jobMemory = 8
        Int timeout = 1
    }
    parameter_meta {
        vcf_file:  "Input vcf file"
        t_depth: "tumour depth filter threshold"
        t_vaf: "Tumor Variant Allele Frequency threshold"
        jobMemory: "Memory allocated indexing job"
        modules:   "Required environment modules"
        timeout:   "Hours before task timeout"
    }

    command <<<
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
    >>>

    runtime {
        memory: "~{jobMemory} GB"
        modules: "~{modules}"
        timeout: "~{timeout}"
    }

    output {
        File vcf_filtered = "filtered.vcf.gz"
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

