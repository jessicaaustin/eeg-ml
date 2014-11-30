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
    userid=$(grep $s subjectToUserId.xls | cut -f2)
    echo "${s}\t${userid}\n"

    asrfile=userdata/${s}_asr.xls
    if [ ! -f $asrfile ];
    then
        head -n1 $ASR_FILE |  cut -f 1,2,3,4,5,7,8,9 > $asrfile
    fi
    grep $s $ASR_FILE | cut -f 1,2,3,4,5,7,8,9 >> $asrfile
    sed -i 's/\\N/0.0000/g' $asrfile

    eegfile=userdata/${s}_eeg.xls
    if [ ! -f $eegfile ];
    then
        head -n1 $EEG_FILE |  cut -f 2,3,4,8,9,10,11-18,21 > $eegfile
    fi
    grep -P "\t${userid}\t" $EEG_FILE | cut -f 2,3,4,8,9,10,11-18,21 >> $eegfile
done 
