

RepeatMasker --lib ../../../REF_DATA/Repeat_libraries/BosTau9_repeat_library.fasta insertion_sequences.sQTL.fa
 cat insertion_sequences.sQTL.fa.out | grep -v "Simple_repeat" | grep -v "Low_complexity" | /cluster/work/pausch/alex/software/RepeatMasker/util/buildSummary.pl - > insertion_sequences.sQTL.fa.tbl
