genomes/.decompressed: genomes.tar.xz
	xz --force --keep --decompress $<
	tar -xf $(basename $<)
	rm $(basename $<)
	touch $@

genomes.tar.xz:
	s3cmd get s3://com-cavescience-r124-manuscript/$@
