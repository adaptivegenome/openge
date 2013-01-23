#!/bin/bash
source $(dirname $0)/../common.sh
rm -f *.sam *.bam

function unbam {
	cp $1 /tmp/unbam.gz
	gunzip /tmp/unbam.gz
	cp /tmp/unbam ${1}raw
	rm /tmp/unbam
}

function cmpbam {
	unbam $1
	unbam $2
	cmp ${1}raw ${2}raw
}

# make a sam reference
samtools view -bS $DATA/simple.sam > simple_.bam
samtools view -b simple_.bam > simple.bam
samtools view -h simple.bam > simple.sam
unbam simple.bam

# compare sam and bam conversions with bamtools to verify sam and bam io
# SAM->SAM
$OGE view --nopg simple.sam -o convert_sam_sam.sam
cmp simple.sam convert_sam_sam.sam || err "SAM SAM conversion failed."
# SAM->BAM
$OGE view --nopg simple.sam -o convert_sam_bam.bam
cmpbam simple.bam convert_sam_bam.bam || err "SAM BAM conversion failed."
# BAM->BAM
$OGE view --nopg simple.bam -o convert_bam_bam.bam
cmpbam simple.bam convert_bam_bam.bam || err "BAM BAM conversion failed."
# BAM->SAM
$OGE view --nopg simple.bam -o convert_bam_sam.sam
cmp simple.sam convert_bam_sam.sam || err "BAM SAM conversion failed."

# test console redirection with bam
# in
$OGE view --nopg < simple.bam -o redirect_in.bam
cmpbam simple.bam redirect_in.bam || err "BAM input redirection failed."
# out
$OGE view --nopg simple.bam --format bam > redirect_out.bam
cmpbam simple.bam redirect_out.bam || err "BAM output redirection failed."
# both
$OGE view --nopg --format bam < simple.bam > redirect_both.bam
cmpbam simple.bam redirect_both.bam || err "BAM I/O redirection failed."

# test console redirection with sam
# in
$OGE view --nopg < simple.sam -o redirect_in.sam
cmp simple.sam redirect_in.sam || err "SAM input redirection failed."
# out
$OGE view --nopg simple.sam > redirect_out.sam
cmp simple.sam redirect_out.sam || err "SAM output redirection failed."
# both
$OGE view --nopg < simple.sam > redirect_both.sam
cmp simple.sam redirect_both.sam || err "SAM I/O redirection failed."

# test piping inputs and outputs
# in
cat simple.bam | $OGE view --nopg -o pipe_in.bam
cmpbam simple.bam pipe_in.bam || err "BAM input piping failed."
# out
$OGE view --nopg simple.bam --format bam | tee pipe_out.bam > /dev/null
cmpbam simple.bam pipe_out.bam || err "BAM output piping failed."
# both
cat simple.bam | $OGE view --nopg --format bam | tee pipe_both.bam > /dev/null
cmpbam simple.bam pipe_both.bam || err "BAM I/O piping failed."

true
