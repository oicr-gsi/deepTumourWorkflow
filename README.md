# deepTumour

The DeepTumour algorithm predicts the tissue of origin of a tumour based on the pattern of passenger mutations identified by Whole Genome Sequencing (WGS).

## Overview

## Dependencies

* [deep-tumour 3.0.1](https://github.com/oicr-gsi/DeepTumour)


## Usage

### Cromwell
```
java -jar cromwell.jar run deepTumour.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`vcf_file`|File|The input vcf file
`outputFileNamePrefix`|String|Prefix for output files
`reference`|String|The genome reference build. For example: hg19, hg38


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`runDeepTumour.jobMemory`|Int|24|Memory allocated indexing job
`runDeepTumour.timeout`|Int|24|Hours before task timeout


### Outputs

Output | Type | Description | Labels
---|---|---|---
`deepTumourOutputJson`|File|the output json assigns a match probability from 0.0 to 1.0 for each of the 29 tumour types on which it was trained and chooses the tumour type with the highest probability score. The algorithm also calculates a type of confidence score based on the probability scores' distributione. A low entropy (< 2.0) is considered a confident score. HIgher values are unreliable (but might be correct).|vidarr_label: deepTumourOutputJson


## Commands
This section lists command(s) run by deepTumour workflow

* Running deepTumour

```
        set -euo pipefail
        
        mkdir out
        python $DEEP_TUMOUR_ROOT/src/DeepTumour.py --vcfFile ~{vcf} --reference $HG19_ROOT/hg19_random.fa ~{liftover} --outDir out
        mv out/predictions_DeepTumour.json ~{outputFileNamePrefix}.predictions_DeepTumour.json

```

 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
