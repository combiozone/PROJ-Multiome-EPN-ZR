#!/bin/sh

${CELLRANGER} count --id=${NAME} \
                    --transcriptome=${ref_gex} \
                    --fastqs=${FQ_DIR} \
                    --create-bam=true \
                    --localcores=16 \
                    --localmem=64


