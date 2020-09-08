cd $1

for BEDFILE in *.bed; do
#for BEDFILE in gencode.as1.bed kras.alu.bed kras.mer.bed sz.sx.sx1.y.bed; do
  docker run \
    -u $(id -u ${USER}):$(id -g ${USER})\
    -v '/public/groups/kimlab/aale.kras/data/bulk.rna.seq/exotic/output/exotic.alu.edit.in':'/data/input_files'\
    -v '/public/groups/kimlab/aale.kras/data/bulk.rna.seq/exotic/output/exotic.alu.edit.out':'/data/output_dir'\
    rna_editing:1.0 RNAEditingIndex -d '/data/input_files/'\
    -l '/data/output_dir/log.files' -o '/data/output_dir/cmpileups'\
    -os '/data/output_dir/summary_dir'\
    -f '.bam'\
    --genome hg38 --verbose -rb '/data/input_files/'$BEDFILE
done

    #-v '/public/groups/kimlab/aale.kras/data/single.cell.rna.seq/output/rna.editing.out/hg38.out':'/data/output_dir' \

