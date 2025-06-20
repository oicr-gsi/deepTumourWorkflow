version 1.0


workflow deepTumour {
    input {
        File vcf_file
        String outputFileNamePrefix     
        String reference
    }

    parameter_meta {
        vcf_file: "The input vcf file"
        outputFileNamePrefix: "Prefix for output files"
        reference: "The genome reference build. For example: hg19, hg38"
    }
    
    call runDeepTumour {
        input:
        vcf = vcf_file,
        outputFileNamePrefix = outputFileNamePrefix,
        reference_genome = reference,
        modules = "deep-tumour/3.0.1 hg19/p13"
    }

    meta {
        author: "Gavin Peng"
        email: "gpeng@oicr.on.ca"
        description: "The DeepTumour algorithm predicts the tissue of origin of a tumour based on the pattern of passenger mutations identified by Whole Genome Sequencing (WGS)."
        dependencies: [
            {
                name: "deep-tumour/3.0.1",
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

task runDeepTumour {
    input {
        File vcf
        String outputFileNamePrefix
        String reference_genome
        Int jobMemory = 24
        String modules
        Int timeout = 24
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
        python $DEEP_TUMOUR_ROOT/src/DeepTumour.py --vcfFile ~{vcf} --reference $HG19_ROOT/hg19_random.fa ~{liftover} --outDir out
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

