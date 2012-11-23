source $(dirname $0)/../common.sh

BAM=$(dirname $0)/801.2L.10k.bam
REGIONS=$(dirname $0)/fly.2L.10k.regions
FASTA=$(dirname $0)/dmel.2L.fa

#known good files were generated with the following command:
#samtools index ${BAM}
#repeatseq -repeatseq -calls 801.2L.10k.bam dmel.2L.fa fly.2L.10k.regions

#delete index
rm dmel.2L.fa.fai

#delete output
rm ${BAM}.repeatseq ${BAM}.calls ${BAM}.vcf

#delete output logs
rm test.out

tar xf dmel.2L.fa.tgz
[ ! -f dmel.2L.fa ] && err "Failed to decompress FASTA reference"

${OGE} repeatseq ${BAM} -L ${REGIONS} -R ${FASTA} --repeatseq --calls -v 2>test.out

#check that output was generated
[ ! -f test.out ] && err "Failed to find test.out"
[ ! -f ${BAM}.repeatseq ] && err "Failed to find ${BAM}.repeatseq"
[ ! -f ${BAM}.calls ] && err "Failed to find ${BAM}.calls"
[ ! -f ${BAM}.vcf ] && err "Failed to find ${BAM}.vcf"

#check output
diff ${BAM}.repeatseq ${BAM}.good.repeatseq|| err "Repeatseq file differs from known good file"
diff ${BAM}.calls ${BAM}.good.calls|| err "Calls file differs from known good file"
diff ${BAM}.vcf ${BAM}.good.vcf|| err "VCF file differs from known good file"

true
