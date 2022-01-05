.PHONY: all clean

files = output.txt spep.csv faint.csv

all : $(files)

output.txt: spep_6sta3n_fix_anno3.txt spep_classifier.R spep_functions.R
	Rscript spep_classifier.R > $@
	perl -pi -e 's/\n/\r\n/g' $@
	dot -Tpng roc_graph.dot > roc_g.png

spep.csv: spep_6sta3n_fix_anno3.txt
	perl -pe 's/,/_/g' $< | perl -pe 's/\t/,/g' > $@

faint.csv: spep.csv
	grep -i faint $< > $@

clean :
	rm -f $(files)
	rm -f Rplots.pdf roc_g.png roc_graph.dot
	rm -f feat_hist_fig3_scale05.tiff feat_hist_fig3_8080.tiff
