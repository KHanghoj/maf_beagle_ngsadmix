maf_beagle: maf_beagle.cpp
	g++ maf_beagle.cpp -o maf_beagle -O3

test: maf_beagle
	zcat test_beagle.gz | ./maf_beagle 0.05 2> test_log.txt | gzip -c > filtered_beagle.gz 

clean:
	rm -f maf_beagle
