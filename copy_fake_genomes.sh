rm samples_fake_genomes.txt

for folder in  /home/ashwini/ash/synthetic-fastq/fake_genome*; do
    echo $folder
    name=${folder##*fake_genome}
    cp $folder/fake_genome-R1.fastq fake_genome_${name}_R1.fastq
    cp $folder/fake_genome-R2.fastq fake_genome_${name}_R2.fastq
    echo fake_genome_${name} fake_genome_${name}_R1.fastq fake_genome_${name}_R2.fastq >> samples_fake_genomes.txt
done
