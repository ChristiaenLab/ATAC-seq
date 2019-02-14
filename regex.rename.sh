for i in *; do mv $i $(echo $i | sed 's/\.[a-z0-9]\+\.fastq\.gz/.fastq.gz/');done
