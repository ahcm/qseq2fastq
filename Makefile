INSTALL_PREFIX?=.
all:
	cd src && make

install:
	mkdir -p $(INSTALL_PREFIX)/bin
	install src/qseq2fastq $(INSTALL_PREFIX)/bin/qseq2fastq
