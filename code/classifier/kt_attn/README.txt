These scripts depend on aligned data in ../../../data-lg

Once you have aligned data, run these in order:
1. findDataWithGoodEEG
2. prepareAllAlignedData

At this point, there are SUBJECTID.mat files in the subjects/ folder. Each 
of these files contains everything needed for running the models.

To prepare data for the HMM, run
1. prepareDataForHMM

This will put data for each subject into word sequences, and save them
to subjects/SUBJECTID_sequences.mat

Then, you can do a few things. To estimate KT params (without any attention 
data), run:
1. estimateParamsWithHMM_KTonly

To investigate how time elapsed and attention are correlated, run:
1. correlateTimeElapsedAndAttention