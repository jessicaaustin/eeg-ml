#!/bin/sh


#EEG_FILE=listen_eeg_2011_2012.xls
#ASR_FILE=eeg_lex_view_2011_2012.xls
#YEAR=2011
#EEG_FILE=listen_eeg_2012_2013.xls
#ASR_FILE=eeg_lex_view_2012_2013.xls
#YEAR=2012
EEG_FILE=listen_eeg_2013_2014.xls
ASR_FILE=eeg_lex_view_2013_2014.xls
YEAR=2013


SUBJECTS=$(cut -f2 ${ASR_FILE} | sort | uniq)

for s in $SUBJECTS
do
    grep $s $ASR_FILE | cut -f 1,2,3,4 > tmp/asr.xls

    MACHINE_NAMES=$(cut -f1 tmp/asr.xls | sort | uniq)
    for m in $MACHINE_NAMES
    do
        grep $m $EEG_FILE | cut -f 1,3,4,5 > tmp/eeg_mn.xls
        grep $m tmp/asr.xls > tmp/asr_mn.xls

        cut -f3 tmp/asr_mn.xls | cut -c1-19 > tmp/times_mn.xls
        found=false
        while read t
        do
            if grep -q "$t" tmp/eeg_mn.xls; then
                grep "$t" tmp/eeg_mn.xls >> users/eeg_found_${s}.xls
                grep "$t" tmp/asr_mn.xls >> users/asr_found_${s}.xls
                found=true
            fi
        done < tmp/times_mn.xls
        if [ "$found" = true ]; then
            echo "match found for $s"
            subjectid=$(cut -f2 users/asr_found_${s}.xls | sort | uniq)
            userid=$(cut -f2 users/eeg_found_${s}.xls | sort | uniq)
            echo "${subjectid}\t${userid}\n"
            echo "${subjectid}\t${userid}" >> subjectToUserId_${YEAR}.xls
            break
        fi
        
    done
done 
