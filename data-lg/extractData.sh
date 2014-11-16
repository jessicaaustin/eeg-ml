#!/bin/bash

cut  -f2,5,7,8,9 eeg_lex_view_2012_2013.xls > eeg_lex_view_2012_2013_cols.xls
cut  -f2,5,7,8,9 eeg_lex_view_2011_2012.xls > eeg_lex_view_2011_2012_cols.xls
cut  -f2,5,7,8,9 eeg_lex_view_2013_2014.xls > eeg_lex_view_2013_2014_cols.xls

sed -i 's/\\N/0.0000/g' eeg_lex_view_2011_2012_cols.xls
sed -i 's/\\N/0.0000/g' eeg_lex_view_2011_2012_cols.xls
sed -i 's/\\N/0.0000/g' eeg_lex_view_2013_2014_cols.xls 

