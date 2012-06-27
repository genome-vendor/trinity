###################################################################

all: 
	cd Inchworm && (test -e configure || autoreconf) \
                && ./configure --prefix=`pwd` && $(MAKE) install
	cd Chrysalis && $(MAKE) UNSUPPORTED=yes
	cd trinity-plugins/rsem && $(MAKE) 
#	cd trinity-plugins/bwa-0.5.7-patched_multi_map && $(MAKE)
	cd trinity-plugins/jellyfish && ./configure --prefix=`pwd` && $(MAKE)
	cd trinity-plugins/fastool && $(MAKE)
	cd trinity-plugins/slclust && $(MAKE)
	cd trinity-plugins/kmer && $(MAKE)

clean:
	cd Inchworm && make clean
	cd Chrysalis && $(MAKE) clean UNSUPPORTED=yes
	cd trinity-plugins/rsem && $(MAKE) clean
#	cd trinity-plugins/bwa-0.5.7-patched_multi_map && $(MAKE) clean
	cd trinity-plugins/jellyfish && $(MAKE) clean 
	cd trinity-plugins/fastool && $(MAKE) clean
	cd trinity-plugins/slclust && $(MAKE) clean
	cd trinity-plugins/kmer && $(MAKE) clean
	find trinity-plugins/kmer \( -name Make.compilers \
                                     -o -name "*.[Cc].d" \) -exec rm -f {} +

###################################################################


