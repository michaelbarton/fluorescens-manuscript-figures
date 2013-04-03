gene_density.eps: bin/plot_gene_density gene_count.csv
	$^ $@

gene_count.csv: bin/count genomes/.decompressed
	echo "species,strain,genome_size,gene_count" > $@
	ls $(dir $(lastword $^))**/**/genome.gb | parallel $< >> $@

genomes/.decompressed: genomes.tar.xz
	xz --force --keep --decompress $<
	tar -xf $(basename $<)
	rm $(basename $<)
	touch $@

genomes.tar.xz:
	s3cmd get s3://com-cavescience-r124-manuscript/$@

clean:
	rm -rf genomes/ *.csv *.eps
