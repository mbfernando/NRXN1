STAR --genomeDir nrxn1_ref/GRCh38_and_mm10/star \
        --readFilesIn ${in_dir}/${in_file}_1.fq.gz ${in_dir}/${in_file}_2.fq.gz \
        --readFilesCommand zcat \
        --outFileNamePrefix out/${in_file}/${in_file}_ \
        --runThreadN 16 \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMstrandField intronMotif \
        --twopassMode Basic
