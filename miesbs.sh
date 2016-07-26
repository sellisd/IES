# miesbs.sh: mobile IES BLAST and SiLiX
# BLAST all against all IESs for each species and cluster results with SiLiX
set -x
blastn -evalue 0.0000001 -query /home/dsellis/data/IES/analysis/mies/fasta/pprIES.fa -db /home/dsellis/data/IES/analysis/mies/bdbs/ppries -outfmt 6 -out /home/dsellis/data/IES/analysis/mies/blastout/ppries.blastout -dust no -task blastn-short -num_threads 7 -max_target_seqs 10000
blastn -evalue 0.0000001 -query /home/dsellis/data/IES/analysis/mies/fasta/pbiIES.fa -db /home/dsellis/data/IES/analysis/mies/bdbs/pbiies -outfmt 6 -out /home/dsellis/data/IES/analysis/mies/blastout/pbiies.blastout -dust no -task blastn-short -num_threads 7 -max_target_seqs 10000
blastn -evalue 0.0000001 -query /home/dsellis/data/IES/analysis/mies/fasta/pteIES.fa -db /home/dsellis/data/IES/analysis/mies/bdbs/pteies -outfmt 6 -out /home/dsellis/data/IES/analysis/mies/blastout/pteies.blastout -dust no -task blastn-short -num_threads 7 -max_target_seqs 10000
blastn -evalue 0.0000001 -query /home/dsellis/data/IES/analysis/mies/fasta/ppeIES.fa -db /home/dsellis/data/IES/analysis/mies/bdbs/ppeies -outfmt 6 -out /home/dsellis/data/IES/analysis/mies/blastout/ppeies.blastout -dust no -task blastn-short -num_threads 7 -max_target_seqs 10000
blastn -evalue 0.0000001 -query /home/dsellis/data/IES/analysis/mies/fasta/pseIES.fa -db /home/dsellis/data/IES/analysis/mies/bdbs/pseies -outfmt 6 -out /home/dsellis/data/IES/analysis/mies/blastout/pseies.blastout -dust no -task blastn-short -num_threads 7 -max_target_seqs 10000
blastn -evalue 0.0000001 -query /home/dsellis/data/IES/analysis/mies/fasta/pocIES.fa -db /home/dsellis/data/IES/analysis/mies/bdbs/pocies -outfmt 6 -out /home/dsellis/data/IES/analysis/mies/blastout/pocies.blastout -dust no -task blastn-short -num_threads 7 -max_target_seqs 10000
blastn -evalue 0.0000001 -query /home/dsellis/data/IES/analysis/mies/fasta/ptrIES.fa -db /home/dsellis/data/IES/analysis/mies/bdbs/ptries -outfmt 6 -out /home/dsellis/data/IES/analysis/mies/blastout/ptries.blastout -dust no -task blastn-short -num_threads 7 -max_target_seqs 10000
blastn -evalue 0.0000001 -query /home/dsellis/data/IES/analysis/mies/fasta/psoIES.fa -db /home/dsellis/data/IES/analysis/mies/bdbs/psoies -outfmt 6 -out /home/dsellis/data/IES/analysis/mies/blastout/psoies.blastout -dust no -task blastn-short -num_threads 7 -max_target_seqs 10000
blastn -evalue 0.0000001 -query /home/dsellis/data/IES/analysis/mies/fasta/pcaIES.fa -db /home/dsellis/data/IES/analysis/mies/bdbs/pcaies -outfmt 6 -out /home/dsellis/data/IES/analysis/mies/blastout/pcaies.blastout -dust no -task blastn-short -num_threads 7 -max_target_seqs 10000

silix /home/dsellis/data/IES/analysis/mies/fasta/pprIES.fa /home/dsellis/data/IES/analysis/mies/blastout/ppries.blastout > /home/dsellis/data/IES/analysis/mies/silixout/ppries.silixout
silix /home/dsellis/data/IES/analysis/mies/fasta/pbiIES.fa /home/dsellis/data/IES/analysis/mies/blastout/pbiies.blastout > /home/dsellis/data/IES/analysis/mies/silixout/pbiies.silixout
silix /home/dsellis/data/IES/analysis/mies/fasta/pteIES.fa /home/dsellis/data/IES/analysis/mies/blastout/pteies.blastout > /home/dsellis/data/IES/analysis/mies/silixout/pteies.silixout
silix /home/dsellis/data/IES/analysis/mies/fasta/ppeIES.fa /home/dsellis/data/IES/analysis/mies/blastout/ppeies.blastout > /home/dsellis/data/IES/analysis/mies/silixout/ppeies.silixout
silix /home/dsellis/data/IES/analysis/mies/fasta/pseIES.fa /home/dsellis/data/IES/analysis/mies/blastout/pseies.blastout > /home/dsellis/data/IES/analysis/mies/silixout/pseies.silixout
silix /home/dsellis/data/IES/analysis/mies/fasta/pocIES.fa /home/dsellis/data/IES/analysis/mies/blastout/pocies.blastout > /home/dsellis/data/IES/analysis/mies/silixout/pocies.silixout
silix /home/dsellis/data/IES/analysis/mies/fasta/ptrIES.fa /home/dsellis/data/IES/analysis/mies/blastout/ptries.blastout > /home/dsellis/data/IES/analysis/mies/silixout/ptries.silixout
silix /home/dsellis/data/IES/analysis/mies/fasta/psoIES.fa /home/dsellis/data/IES/analysis/mies/blastout/psoies.blastout > /home/dsellis/data/IES/analysis/mies/silixout/psoies.silixout
silix /home/dsellis/data/IES/analysis/mies/fasta/pcaIES.fa /home/dsellis/data/IES/analysis/mies/blastout/pcaies.blastout > /home/dsellis/data/IES/analysis/mies/silixout/pcaies.silixout
