conda activate transdecoder

dir="MMETSP0754.transdecoder"
fasta="mmetsp/transcriptome/MMETSP0754.trinity_out_2.2.0.Trinity.fasta.renamed.fasta"
fasta_name="MMETSP0754.trinity_out_2.2.0.Trinity.fasta.renamed.fasta"

mkdir -p $dir
cp $fasta $dir
cd $dir

TransDecoder.LongOrfs -t $fasta_name
TransDecoder.Predict -t  $fasta_name 

