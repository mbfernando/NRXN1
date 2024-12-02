# Adapted from https://gist.github.com/aryarm/f001b808d2827f09589e9cdb05d710b3
samtools index ${in_fn}
samtools view --no-PG -bh ${in_fn} GRCh38_chr{1..22} GRCh38_chr{X,Y,M} | \
samtools reheader --no-PG -c 'perl -pe "s/^(@SQ.*)(\tSN:GRCh38_chrM)/\$1\${2}T/"' - | \
samtools reheader --no-PG -c 'perl -pe "s/^(@SQ.*)(\tSN:)GRCh38_chr/\$1\$2/"' - | \
samtools view --no-PG -h | \
grep -Ev '^@SQ' | \
samtools view --no-PG -f2 -F4 -F256 -bhT hg38.fa.gz - 2>/dev/null > out/${in_file}_hg38.bam
samtools index out/${in_file}_hg38.bam
