all: gene_density.eps dotplot.eps growth.eps

growth.eps: bin/plot_growth growth_data.csv
	$^ $@

dotplot.eps: bin/dotplot synteny.tab
	$^ $@

synteny.tab: bin/align genomes/.fasta
	ls $(dir $(lastword $^))**/**/genome.fna \
		| grep "pf01\|sbw25\|pf5" \
		| parallel $< genomes/fluorescens/r124/genome.fna {}
	sleep 1
	cat *.coords | sort > $@
	rm *.coords *.delta

genomes/.fasta: bin/fasta genomes/.decompressed
	ls $(dir $(lastword $^))**/**/genome.gb | parallel $<
	touch $@

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
