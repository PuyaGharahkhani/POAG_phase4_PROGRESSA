#!/bin/bash
#PBS -l ncpus=1
#PBS -l walltime=12:00:00
#PBS -l mem=8750mb
#PBS -N clean_GWAS2016
#PBS -j oe
#PBS -o /mnt/lustre/home/matthewL/condor/clean_GWAS2016.log

#Author: Matthew H. Law
#Email: matthew.law@qimrberghofer.edu.au
#Date Created: 16:31 13/Jun/2016
#Date Modified: 12/09/2016

########################################################################
#0.1 Updates
########################################################################
#12/09/2016 script using > 1 cpu, found ~34 instances where didn't add --threads 1 to plink.
#08/09/2016 while I don't have the OCAC control genotype data testing with the imputed data shows they will be suitable (see /mnt/lustre/home/matthewL/Melanoma/MIA/SCRIPTS/script_to_clean_MIA_KK_data_simplified_fork.sh and I have cleaned the MIA data to the OCAC guidelines, so work that in now and add in the control data when ready. That way I can use the current 2016 cleaned case sets for the breslow work. Run through manually seems fine (had to fix up a few indels that had the wrong overlapping SNP ID via the oncoarray remap files) but run formally overnight and continue Monday. >>>worked perfectly for MIA

#25/08/2016 pending incorporating the new extraction of the MIA data I am renaming and saving the imputation results. Is laborious as you have to unpack huge zips, rename them, check for correct line breaks (sometimes the final lines are missing line breaks) and then rezip.

#25/08/2016 some how a section of the lifover functions has been copied multiple time through the script. I think I have cleaned out each instance but script may abort when it hits a copied chunk I missed... likewise sections of the initial cleaning of EPIGENE/QTWINS copied into a couple of places, no idea how this happened. Still if the script suddently fails when it worked likely due to this.

#14/08/2016 Scripts works fine, and AMFS data imputed for all SNPs for Autosomes (chrX not implimented yet). No major errors from a quick scan of the logs, though the QC files do flag a number of SNPs with very different freqs to the ref panel. Not sure if this is just chance or not. Only problem with starting work before I go on leave is it is TOO quick - barely 36 hours to run AMFS, and they only keep the data for a week... 
   #minimac 3 gives a info file like impute2, so should be trivial to flag bad SNPs, though chr22 looks good. Will need a remap for the chr:POS as well through

#09/08/2016 looks to work fine for all current sets. Now just weaving in EPIGENE, and then can start the rest of the cleaning/testing prior to imputation. Really should weave QTWNS into the control comparisons etc, then clean, then recombine but that it not done yet
#08/08/2016 cleaned up and improved the method for finding SNPs duplicated on the same array set (e.g. express + exome) and all the sets eem to merge properly now. Running from start through to generating the merged case control sets 
#03/08/2016 most recent run through fine but still identifying a handful of SNPs when merging the case/control sets that need to be updated or excluded. As such the noted SNP counts etc might be out by 1-2 after a run through but not worth keeping updating the notes until I can get the case and control sets to merge without errors - still hanging on IBD + OAG, see if latest fixes fix, and then merge WAMHS with IBD + OAG then done....
#01-02/08/2016 - kept finding a few more SNPs with incorrect and hard to resolve positions (e.g. two AG SNPs next to each other with the wrong positions, correcting).
#28/07/2016 +++++ re-merging case and controls had IDs ~250 SNPs with wrong positions in the various starting sets, so added an early step to collate these and fix across all sets. Tetesting first on OAG1 and OAG2 as OAG2 has ~228 of them itself. The rest are in omni/amfs, and only showed up when omni merged with SHD, but are wrong in AMFS as well.....
    #++++++first past test with OAG1 failed as I had duplicated some of the SNPs in the remap list, so redo. May as well run all the samples rather than just OAG
#27/07/2016 Initial cleaning and aligning to 1kG runs through for all sets
#26/07/2016 - script looks to run through for all sets but redo IBD (output to a diff log file) as I don't think I cleared temp files from IBD testing before running....
#25/07/2016 after a few false starts and tweaks all the datasets seems to be intergrated into the cleaning, check the output though
#21/07/2016 checking this over, found that for melanoma samples I wasn't actually writing out the SNPSEX, just the whole file (which would have resulted in the PED sex being used). Fixed, need to rerun for all (minimal effect i imagine as only 6 people with changed sex but nice to be accurate)
#21/07/2016 looking good, work in GERMAN and WAMHS now that the core new functions are working.
#20/07/2016 further tweaks - flipscan revealed merge-euqal-pos directly merges regardless of the quality of the data; it merged a SNP with missing data at every position with an overlapping 100% genotyped SNP, and worse used the alleles of the empty SNP - such that the correct SNP ended up with C->T and T->C, breaking LD... wrote and tested a manual version
#13/07/2016 a few more tweaks to how positional overlaps are handled - one to better align SNPs where there is an indel at the same position as a SNP - and two to properly remove SNPS that are ambigious, overlap an ambigious SNP in 1kg by position, but have a MAF too different to safely align (before was silently ignoring them). Can now succesfully merge the three case/controls sets back together, and merge any overlapping SNPs successfully. Still need to update the final SNP N for the latest script run through, and perform QC on the merged case control sets - run a quick GWAS first to find any badlty wrong SNPS....
#12/07/2016 control v control GWAS revealed ridiculous p values - realised needed to drop any ambigious SNP that can't be mapped to 1KG (as can't align strand), and also run some basic HWE checks etc before running control-control GWAS. Post this the results were what you expect for a dummy GWAS - one or 2 SNPs ~ 5e-7, a few more at 5e-6. Nothing of concern, data looks good
#11/07/2016 Fixed the re,xlab="CHROMOSOME",ylab="PERCENTAGE SNPs with freq chisq diff > 300"maining het haploid/non missing non male Y chromsome calls. Did a sex check as a part of this so the sample and SNP N trivially changes (loses 4 people with ambigious gender by chrX calls, and ~8k SNPs that are het in males). Otherwise alignment seems to work
#06/07/2016 found I was missing some classes of triplicate SNPs so expanded that check to cover, I think, all possible permutations of 3 alleles. Also added a miinor if test step update the alleles where I/D indels could be aligned to actual alleles in 1KG. Might make a slight difference to N output, check. Due to not capturing all permutation of triplicates AMFS_controls failed one of the early checks (SNP MAF 0.49 so A1 A2 order swapped in controls, and the case trip was checked for (A C vs C T) but no the controls (C A bs C T). So AMFS controls need to be thoroughly checked while the others just update SNP N output if need be
#13/06/2016 - 30/06/2016 lots of error checking, bug fixing, but aligning to 1kg seems to be working finally and SDH, 610k, 610k_controls run through without error. Next run as a loop with all the samples..
#11/11/2016 - PG started to modify the script for POAG phase 4 and 5 

########################################################################
#0.2 Background
########################################################################

#Main cleaning steps - updated from earlier steps with the new guidelines from MI
#1. Remove ambigious A/T and G/C SNPs , check is hg19 <<<aligned those I could to 1000 genomes
#2. Where required merge with control set (will likely require some flipping and manual checking) <<done
#3. Clean on missingness (--mind 0.03 --geno 0.03 --maf 0.001 ) then HWE (10-4 in controls, 10-10 in cases). <<done
#4 Also filter on extreme het in samples (see the MIA scripts etc; MI suggests >0.05 or < -0.05) <<done
#4. Check sex matches <<done
#5. Test for relationship. This has been decreased in Mark's new protocol -  (PI_HAT>0.15) <<done 
#6 Non european outliers  <<done
#1.6. Do a control-control GWAS between controls to ID any abnormalities in the data <<<aligned those I could to 1000 genomes, done control-control GWAS as well

########################################################################
#0.3 Variables
########################################################################
#working
dirI=/working/lab_stuartma/puyaG/glaucoma/POAG2016genotypes/ANZRAG_phase4-Non-advancedGlaucoma_phase5-Progressa

#write out
dirK=/working/lab_stuartma/puyaG/glaucoma/POAG2016genotypes/ANZRAG_phase4-Non-advancedGlaucoma_phase5-Progressa/cleaned

#liftover
#dirL=/working/lab_stuartma/scratch/matthewL/melanoma_meta_analysis/MANIFEST_CHAIN_QC_FILES

#maf threshold for handling ambigious SNPs
MAF=0.35

#difference between 1kg and data to all aligning/retaining ambigious SNPs
DIFF=0.1

#Raw data Non-advancedGlaucoma
dirJ=/working/lab_stuartma/puyaG/glaucoma/POAG2016genotypes/20150615_Hewitt_Non-advancedGlaucoma_DataRelease/PLINK_201016_0109
raw_data=20150615_Hewitt_Non-advancedGlaucoma

#Progressa
dirJB=/working/lab_stuartma/puyaG/glaucoma/POAG2016genotypes/20150615_Hewitt_Progressa_DataRelease/PLINK_191016_0336
raw_data2=20150615_Hewitt_Progressa

#ENDO cases and controls to be used as controls here:
#first extract ENDO cases and controls from the Release8_Observed_CoreExome data
grep ENDO_CoreExome /reference/genepi/GWAS_release/Release8/Release8_Observed_CoreExome/GWAS_preMZ_projectbysample.txt > ${dirI}/ENDO_CoreExome
awk 'NR==FNR{a[$1];next}($2 in a)' ${dirI}/ENDO_CoreExome /reference/genepi/GWAS_release/Release8/Release8_Observed_CoreExome/PLINK_format/GWAS_b37PlusStrand.fam | awk '{print $1, $2}' > ${dirI}/ENDO_CoreExome_extract
module load plink/1.90b3.40
plink --bfile /reference/genepi/GWAS_release/Release8/Release8_Observed_CoreExome/PLINK_format/GWAS_b37PlusStrand --keep ${dirI}/ENDO_CoreExome_extract --make-bed --out ${dirI}/ENDO_Release8_Observed_CoreExome
raw_data3=${dirI}/ENDO_Release8_Observed_CoreExome

#POAG phase 3:
dirJG=/working/lab_stuartma/puyaG/glaucoma/POAG_cases_recent_011015_0928
raw_data4=ANZRAG_POAG_cases_phase3_QIMR_twins_controls_merged_cleaned5_relatedness_removed2_ancestry_outliers_removed

#OAG phase 1 and 2
dirOAG=/working/lab_stuartma/puyaG/glaucoma
raw_data5=WAMHS_GLUACOMA_IBD_BEACON_final_clean_IBD_PCA_outliers_removed

########################################################################
#0.4 Functions
########################################################################
#some of these functions are clunky to have in the script, and I call them more than once - so define them here

#awk_check_based_on_ID_1KG; takes SNP allele and freq data paired with 1kG data that has the same rsID and works out if they are the same SNP, and aligned, and if their positions (CHR, BP, are different). To restore the original function take the middle - XXXX in BEGIN {} {XXXXX} END {} and drop the \ before each col (e.g. \$3 becomes $3
  #This is ugly but it checks if (1) are SNPs non amibigious and (2) and on the same strand,  and if so prints them to concordant ( CHR and BP the same), or same strand but a difference in chr, bp or both. If not (2) - on the other strand, but still non-ambigious, it flags them as flipped but with matching/non matching chr and bp. For non ambifious SNPs it doesn't care which allele is ref/non ref because these are non ambigious. If not (1) - if ambigious - it additionally checks if the PLINK MAF is low enough to sure of selecting correct allele based on MAF (currently 0.35 - too close to 0.5 and you can't judge which strand). If it is, it then tries to work which allele is being tested relative to 1kG, and it calls it as such if the difference in allele frequency (DIFF) is no more than 10%. If it is within 10% it then assigns them  to flipped or not, again with and without matching CHR/BP. 15/06/2016 checked some of the results, found the allele freq test for ambigious SNPs etc was using the wrong field - should use 12 (alt allele freq from 1kg); was using 13, maf. also found my handling function was wrong for some. Handling ambigious alleles with plink (A1_minor, A2, MAF) data versus 1kg (a0,a1,a1_freq) is a headache but it helps if you think of it this way: 
    #ambigious, same stand (aligned), same allele = A B MAF B A a1_freq=MAF
    #ambigious, same stand (aligned), other allele = A B MAF A B a1_freq=(1 - MAF)
    #ambigious, other stand (flipped), same allele = A B MAF A B a1_freq=MAF
    #ambigious, same stand (flipped), other allele = A B MAF B A a1_freq=(1 - MAF)

  #printing to multiple files awk '{if($2>10) {print > "outfile1"} else {print > "outfile2"}}' infile. Print SNPs out to aligned, not aligned, so I can check total to make sure the function works. Handling indels if fiddly. Not many overlap when you merge based on bp (4 in melanoma.bim) but worth fixing and giving them the right allele. Can't work out how to also handle variables in the multiple outfile format (inside "") so use generic file names and remove at the end of each loop
  #27/06/2016 tweaked the temp_ambigious_MAF_to_high_too_assess to specifically check for ambigious alleles in both the PLINK and 1KG data; previously was just checking in the PLINK data and it was double assigning SNPs that were triallelic and abimigious in PLINK (e.g. A/T in plink and A/G in 1KG). 
  #28/06/2016 was missing some combinations of ambigious high maf alleles to simply exclude (they were being handled correctly in the aligning funciton, but not the search on maf and and drop

(cat <<__EOF__
#!/usr/bin/awk -f
# written by script_to_clean_QC_ausmel_GWAS_data_2016.sh.
# takes SNP allele and freq data paired with 1kG data that has the same rsID and works out if they are the same SNP, and aligned, and if their positions (CHR, BP, are different)
BEGIN {}
{
    if( !(\$3=="A" && \$4=="T") && !(\$3=="T" && \$4=="A") && !(\$3=="C" && \$4=="G") && !(\$3=="G" && \$4=="C") && \$1==\$8 && \$7==\$9 && ( (\$3==\$10 && \$4==\$11) || (\$3==\$11 && \$4==\$10 ) ) )
      {print \$2 > "temp_aligned_concordant"};
    if( !(\$3=="A" && \$4=="T") && !(\$3=="T" && \$4=="A") && !(\$3=="C" && \$4=="G") && !(\$3=="G" && \$4=="C") && \$1!=\$8 && \$7==\$9 && ( (\$3==\$10 && \$4==\$11) || (\$3==\$11 && \$4==\$10 ) ) )
      {print \$2 > "temp_aligned_chr_disc_bp_conc"};
    if( !(\$3=="A" && \$4=="T") && !(\$3=="T" && \$4=="A") && !(\$3=="C" && \$4=="G") && !(\$3=="G" && \$4=="C") && \$1==\$8 && \$7!=\$9 && ( (\$3==\$10 && \$4==\$11) || (\$3==\$11 && \$4==\$10 ) ) )
      {print \$2 > "temp_aligned_chr_conc_bp_disc"};
    if( !(\$3=="A" && \$4=="T") && !(\$3=="T" && \$4=="A") && !(\$3=="C" && \$4=="G") && !(\$3=="G" && \$4=="C") && \$1!=\$8 && \$7!=\$9 && ( (\$3==\$10 && \$4==\$11) || (\$3==\$11 && \$4==\$10 ) ) )
      {print \$2 > "temp_aligned_chr_disc_bp_disc"};
    if( \$1==\$8 && \$7==\$9 && ( (\$3=="G" && \$4=="A" && \$10=="C" && \$11=="T" ) || (\$3=="G" && \$4=="A" && \$10=="T" && \$11=="C" ) || (\$3=="G" && \$4=="T" && \$10=="C" && \$11=="A" ) || (\$3=="G" && \$4=="T" && \$10=="A" && \$11=="C") || (\$3=="C" && \$4=="A" && \$10=="G" && \$11=="T") || (\$3=="C" && \$4=="A" && \$10=="T" && \$11=="G" ) || (\$3=="C" && \$4=="T" && \$10=="G" && \$11=="A") || (\$3=="C" && \$4=="T" && \$10=="A" && \$11=="G") || (\$3=="A" && \$4=="G" && \$10=="C" && \$11=="T" ) || (\$3=="A" && \$4=="G" && \$10=="T" && \$11=="C" ) || (\$3=="T" && \$4=="G" && \$10=="C" && \$11=="A" ) || (\$3=="T" && \$4=="G" && \$10=="A" && \$11=="C" ) || (\$3=="A" && \$4=="C" && \$10=="G" && \$11=="T" ) || (\$3=="A" && \$4=="C" && \$10=="T" && \$11=="G" ) || (\$3=="T" && \$4=="C" && \$10=="G" && \$11=="A" ) || (\$3=="T" && \$4=="C" && \$10=="A" && \$11=="G") ) )
      {print \$2 > "temp_flipped_concordant"};
    if( \$1!=\$8 && \$7==\$9 && ( (\$3=="G" && \$4=="A" && \$10=="C" && \$11=="T" ) || (\$3=="G" && \$4=="A" && \$10=="T" && \$11=="C" ) || (\$3=="G" && \$4=="T" && \$10=="C" && \$11=="A" ) || (\$3=="G" && \$4=="T" && \$10=="A" && \$11=="C") || (\$3=="C" && \$4=="A" && \$10=="G" && \$11=="T") || (\$3=="C" && \$4=="A" && \$10=="T" && \$11=="G" ) || (\$3=="C" && \$4=="T" && \$10=="G" && \$11=="A") || (\$3=="C" && \$4=="T" && \$10=="A" && \$11=="G") || (\$3=="A" && \$4=="G" && \$10=="C" && \$11=="T" ) || (\$3=="A" && \$4=="G" && \$10=="T" && \$11=="C" ) || (\$3=="T" && \$4=="G" && \$10=="C" && \$11=="A" ) || (\$3=="T" && \$4=="G" && \$10=="A" && \$11=="C" ) || (\$3=="A" && \$4=="C" && \$10=="G" && \$11=="T" ) || (\$3=="A" && \$4=="C" && \$10=="T" && \$11=="G" ) || (\$3=="T" && \$4=="C" && \$10=="G" && \$11=="A" ) || (\$3=="T" && \$4=="C" && \$10=="A" && \$11=="G") ) )
      {print \$2 > "temp_flipped_chr_disc_bp_conc"};
    if( \$1==\$8 && \$7!=\$9 && ( (\$3=="G" && \$4=="A" && \$10=="C" && \$11=="T" ) || (\$3=="G" && \$4=="A" && \$10=="T" && \$11=="C" ) || (\$3=="G" && \$4=="T" && \$10=="C" && \$11=="A" ) || (\$3=="G" && \$4=="T" && \$10=="A" && \$11=="C") || (\$3=="C" && \$4=="A" && \$10=="G" && \$11=="T") || (\$3=="C" && \$4=="A" && \$10=="T" && \$11=="G" ) || (\$3=="C" && \$4=="T" && \$10=="G" && \$11=="A") || (\$3=="C" && \$4=="T" && \$10=="A" && \$11=="G") || (\$3=="A" && \$4=="G" && \$10=="C" && \$11=="T" ) || (\$3=="A" && \$4=="G" && \$10=="T" && \$11=="C" ) || (\$3=="T" && \$4=="G" && \$10=="C" && \$11=="A" ) || (\$3=="T" && \$4=="G" && \$10=="A" && \$11=="C" ) || (\$3=="A" && \$4=="C" && \$10=="G" && \$11=="T" ) || (\$3=="A" && \$4=="C" && \$10=="T" && \$11=="G" ) || (\$3=="T" && \$4=="C" && \$10=="G" && \$11=="A" ) || (\$3=="T" && \$4=="C" && \$10=="A" && \$11=="G") ) )
     {print \$2 > "temp_flipped_chr_con_bp_disc"};
    if( \$1!=\$8 && \$7!=\$9 && ( (\$3=="G" && \$4=="A" && \$10=="C" && \$11=="T" ) || (\$3=="G" && \$4=="A" && \$10=="T" && \$11=="C" ) || (\$3=="G" && \$4=="T" && \$10=="C" && \$11=="A" ) || (\$3=="G" && \$4=="T" && \$10=="A" && \$11=="C") || (\$3=="C" && \$4=="A" && \$10=="G" && \$11=="T") || (\$3=="C" && \$4=="A" && \$10=="T" && \$11=="G" ) || (\$3=="C" && \$4=="T" && \$10=="G" && \$11=="A") || (\$3=="C" && \$4=="T" && \$10=="A" && \$11=="G") || (\$3=="A" && \$4=="G" && \$10=="C" && \$11=="T" ) || (\$3=="A" && \$4=="G" && \$10=="T" && \$11=="C" ) || (\$3=="T" && \$4=="G" && \$10=="C" && \$11=="A" ) || (\$3=="T" && \$4=="G" && \$10=="A" && \$11=="C" ) || (\$3=="A" && \$4=="C" && \$10=="G" && \$11=="T" ) || (\$3=="A" && \$4=="C" && \$10=="T" && \$11=="G" ) || (\$3=="T" && \$4=="C" && \$10=="G" && \$11=="A" ) || (\$3=="T" && \$4=="C" && \$10=="A" && \$11=="G") ) )
     {print \$2 > "temp_flipped_chr_disc_bp_disc"};
  if( ( (\$3=="A" && \$4=="T" && \$10=="T" && \$11=="A") || (\$3=="T" && \$4=="A" && \$10=="A" && \$11=="T") || (\$3=="A" && \$4=="T" && \$10=="A" && \$11=="T") || (\$3=="T" && \$4=="A" && \$10=="T" && \$11=="A") || (\$3=="G" && \$4=="C" && \$10=="C" && \$11=="G") || (\$3=="C" && \$4=="G" && \$10=="G" && \$11=="C")   || (\$3=="G" && \$4=="C" && \$10=="G" && \$11=="C") || (\$3=="C" && \$4=="G" && \$10=="C" && \$11=="G") )  && \$5>= MAF)
      {print \$2 > "temp_ambigious_MAF_to_high_too_assess"};
  if( \$1==\$8 && \$7==\$9 && ( (\$3=="A" && \$4=="T" && \$10=="T" && \$11=="A") || (\$3=="T" && \$4=="A" && \$10=="A" && \$11=="T") || (\$3=="G" && \$4=="C" && \$10=="C" && \$11=="G") || (\$3=="C" && \$4=="G" && \$10=="G" && \$11=="C") ) && \$5 < MAF && ((\$5 - \$12) < DIFF ) && ((\$5 - \$12) > -(DIFF) )  )
      {print \$2 > "temp_ambigious_aligned_concordant_same_allele"};
  if( \$1!=\$8 && \$7==\$9 && ( (\$3=="A" && \$4=="T" && \$10=="T" && \$11=="A") || (\$3=="T" && \$4=="A" && \$10=="A" && \$11=="T") || (\$3=="G" && \$4=="C" && \$10=="C" && \$11=="G") || (\$3=="C" && \$4=="G" && \$10=="G" && \$11=="C") ) && \$5 < MAF && ((\$5 - \$12) < DIFF ) && ((\$5 - \$12) > -(DIFF) )  )
      {print \$2 > "temp_ambigious_aligned_chr_disc_bp_conc_same_allele"};
  if( \$1==\$8 && \$7!=\$9 && ( (\$3=="A" && \$4=="T" && \$10=="T" && \$11=="A") || (\$3=="T" && \$4=="A" && \$10=="A" && \$11=="T") || (\$3=="G" && \$4=="C" && \$10=="C" && \$11=="G") || (\$3=="C" && \$4=="G" && \$10=="G" && \$11=="C") ) && \$5 < MAF && ((\$5 - \$12) < DIFF ) && ((\$5 - \$12) > -(DIFF) )  )
      {print \$2 > "temp_ambigious_aligned_chr_conc_bp_disc_same_allele"};
  if( \$1!=\$8 && \$7!=\$9 && ( (\$3=="A" && \$4=="T" && \$10=="T" && \$11=="A") || (\$3=="T" && \$4=="A" && \$10=="A" && \$11=="T") || (\$3=="G" && \$4=="C" && \$10=="C" && \$11=="G") || (\$3=="C" && \$4=="G" && \$10=="G" && \$11=="C") ) && \$5 < MAF && ((\$5 - \$12) < DIFF ) && ((\$5 - \$12) > -(DIFF) )  )
      {print \$2 > "temp_ambigious_aligned_chr_disc_bp_disc_same_allele"};
  if( \$1==\$8 && \$7==\$9 && ( (\$3=="A" && \$4=="T" && \$10=="A" && \$11=="T") || (\$3=="T" && \$4=="A" && \$10=="T" && \$11=="A") || (\$3=="G" && \$4=="C" && \$10=="G" && \$11=="C") || (\$3=="C" && \$4=="G" && \$10=="C" && \$11=="G") ) && \$5 < MAF && ((\$5 - (1 - \$12)) < DIFF ) && ((\$5 - (1-\$12)) > -(DIFF) )    )
      {print \$2 > "temp_ambigious_aligned_concordant_other_allele"};
  if( \$1!=\$8 && \$7==\$9 && ( (\$3=="A" && \$4=="T" && \$10=="A" && \$11=="T") || (\$3=="T" && \$4=="A" && \$10=="T" && \$11=="A") || (\$3=="G" && \$4=="C" && \$10=="G" && \$11=="C") || (\$3=="C" && \$4=="G" && \$10=="C" && \$11=="G") ) && \$5 < MAF && ((\$5 - (1 - \$12)) < DIFF ) && ((\$5 - (1-\$12)) > -(DIFF) )   ) 
      {print \$2 > "temp_ambigious_aligned_chr_disc_bp_conc_other_allele"};
  if( \$1==\$8 && \$7!=\$9 && ( (\$3=="A" && \$4=="T" && \$10=="A" && \$11=="T") || (\$3=="T" && \$4=="A" && \$10=="T" && \$11=="A") || (\$3=="G" && \$4=="C" && \$10=="G" && \$11=="C") || (\$3=="C" && \$4=="G" && \$10=="C" && \$11=="G") ) && \$5 < MAF && ((\$5 - (1 - \$12)) < DIFF ) && ((\$5 - (1-\$12)) > -(DIFF) )    )
      {print \$2 > "temp_ambigious_aligned_chr_conc_bp_disc_other_allele"};
  if( \$1!=\$8 && \$7!=\$9 && ( (\$3=="A" && \$4=="T" && \$10=="A" && \$11=="T") || (\$3=="T" && \$4=="A" && \$10=="T" && \$11=="A") || (\$3=="G" && \$4=="C" && \$10=="G" && \$11=="C") || (\$3=="C" && \$4=="G" && \$10=="C" && \$11=="G") ) && \$5 < MAF && ((\$5 - (1 - \$12)) < DIFF ) && ((\$5 - (1-\$12)) > -(DIFF) )    )
      {print \$2 > "temp_ambigious_aligned_chr_disc_bp_disc_other_allele"};
  if( \$1==\$8 && \$7==\$9 && ( (\$3=="A" && \$4=="T" && \$10=="A" && \$11=="T") || (\$3=="T" && \$4=="A" && \$10=="T" && \$11=="A") || (\$3=="G" && \$4=="C" && \$10=="G" && \$11=="C") || (\$3=="C" && \$4=="G" && \$10=="C" && \$11=="G") )  && \$5 < MAF && ((\$5 - \$12) < DIFF ) && ((\$5 - \$12) > -(DIFF) )   )
      {print \$2 > "temp_ambigious_flipped_concordant_same_allele"};
  if( \$1!=\$8 && \$7==\$9 && ( (\$3=="A" && \$4=="T" && \$10=="A" && \$11=="T") || (\$3=="T" && \$4=="A" && \$10=="T" && \$11=="A") || (\$3=="G" && \$4=="C" && \$10=="G" && \$11=="C") || (\$3=="C" && \$4=="G" && \$10=="C" && \$11=="G") )  && \$5 < MAF && ((\$5 - \$12) < DIFF ) && ((\$5 - \$12) > -(DIFF) )   )
      {print \$2 > "temp_ambigious_flipped_chr_disc_bp_conc_same_allele"};
  if( \$1==\$8 && \$7!=\$9 && ( (\$3=="A" && \$4=="T" && \$10=="A" && \$11=="T") || (\$3=="T" && \$4=="A" && \$10=="T" && \$11=="A") || (\$3=="G" && \$4=="C" && \$10=="G" && \$11=="C") || (\$3=="C" && \$4=="G" && \$10=="C" && \$11=="G") )  && \$5 < MAF && ((\$5 - \$12) < DIFF ) && ((\$5 - \$12) > -(DIFF) )   )
   {print \$2 > "temp_ambigious_flipped_chr_conc_bp_disc_same_allele"};
  if( \$1!=\$8 && \$7!=\$9 && ( (\$3=="A" && \$4=="T" && \$10=="A" && \$11=="T") || (\$3=="T" && \$4=="A" && \$10=="T" && \$11=="A") || (\$3=="G" && \$4=="C" && \$10=="G" && \$11=="C") || (\$3=="C" && \$4=="G" && \$10=="C" && \$11=="G") )  && \$5 < MAF && ((\$5 - \$12) < DIFF ) && ((\$5 - \$12) > -(DIFF) )   )
      {print \$2 > "temp_ambigious_flipped_chr_disc_bp_disc_same_allele"};
  if( \$1==\$8 && \$7==\$9 && ( (\$3=="A" && \$4=="T" && \$10=="T" && \$11=="A") || (\$3=="T" && \$4=="A" && \$10=="A" && \$11=="T") || (\$3=="G" && \$4=="C" && \$10=="C" && \$11=="G") || (\$3=="C" && \$4=="G" && \$10=="G" && \$11=="C") ) && \$5 < MAF && ((\$5 - (1 - \$12)) < DIFF ) && ((\$5 - (1-\$12)) > -(DIFF) )  )
      {print \$2 > "temp_ambigious_flipped_concordant_other_allele"};
  if( \$1!=\$8 && \$7==\$9 && ( (\$3=="A" && \$4=="T" && \$10=="T" && \$11=="A") || (\$3=="T" && \$4=="A" && \$10=="A" && \$11=="T") || (\$3=="G" && \$4=="C" && \$10=="C" && \$11=="G") || (\$3=="C" && \$4=="G" && \$10=="G" && \$11=="C") ) && \$5 < MAF && ((\$5 - (1 - \$12)) < DIFF ) && ((\$5 - (1-\$12)) > -(DIFF) )  )
     {print \$2 > "temp_ambigious_flipped_chr_disc_bp_conc_other_allele"};
  if( \$1==\$8 && \$7!=\$9 && ( (\$3=="A" && \$4=="T" && \$10=="T" && \$11=="A") || (\$3=="T" && \$4=="A" && \$10=="A" && \$11=="T") || (\$3=="G" && \$4=="C" && \$10=="C" && \$11=="G") || (\$3=="C" && \$4=="G" && \$10=="G" && \$11=="C") ) && \$5 < MAF && ((\$5 - (1 - \$12)) < DIFF ) && ((\$5 - (1-\$12)) > -(DIFF) )  )
   {print \$2 > "temp_ambigious_flipped_chr_conc_bp_disc_other_allele"};
  if( \$1!=\$8 && \$7!=\$9 && ( (\$3=="A" && \$4=="T" && \$10=="T" && \$11=="A") || (\$3=="T" && \$4=="A" && \$10=="A" && \$11=="T") || (\$3=="G" && \$4=="C" && \$10=="C" && \$11=="G") || (\$3=="C" && \$4=="G" && \$10=="G" && \$11=="C") ) && \$5 < MAF && ((\$5 - (1 - \$12)) < DIFF ) && ((\$5 - (1-\$12)) > -(DIFF) )  )
     {print \$2 > "temp_ambigious_flipped_chr_disc_bp_disc_other_allele"};
  if ( (\$3=="A" && \$4=="C" && \$10=="C" && \$11=="G") || (\$3=="A" && \$4=="C" && \$10=="G" && \$11=="C") || (\$3=="A" && \$4=="C" && \$10=="C" && \$11=="T") || (\$3=="A" && \$4=="C" && \$10=="T" && \$11=="C") || (\$3=="A" && \$4=="C" && \$10=="G" && \$11=="A") || (\$3=="A" && \$4=="C" && \$10=="A" && \$11=="G") || (\$3=="A" && \$4=="C" && \$10=="T" && \$11=="A") || (\$3=="A" && \$4=="C" && \$10=="A" && \$11=="T") || (\$3=="A" && \$4=="G" && \$10=="A" && \$11=="C") || (\$3=="A" && \$4=="G" && \$10=="C" && \$11=="A") || (\$3=="A" && \$4=="G" && \$10=="C" && \$11=="G") || (\$3=="A" && \$4=="G" && \$10=="G" && \$11=="C") || (\$3=="A" && \$4=="G" && \$10=="G" && \$11=="T") || (\$3=="A" && \$4=="G" && \$10=="T" && \$11=="G") || (\$3=="A" && \$4=="G" && \$10=="A" && \$11=="T") || (\$3=="A" && \$4=="G" && \$10=="T" && \$11=="A") || (\$3=="A" && \$4=="T" && \$10=="T" && \$11=="C") || (\$3=="A" && \$4=="T" && \$10=="C" && \$11=="T") || (\$3=="A" && \$4=="T" && \$10=="A" && \$11=="G") || (\$3=="A" && \$4=="T" && \$10=="G" && \$11=="A") || (\$3=="A" && \$4=="T" && \$10=="C" && \$11=="G") || (\$3=="A" && \$4=="T" && \$10=="G" && \$11=="C") || (\$3=="A" && \$4=="T" && \$10=="A" && \$11=="C") || (\$3=="A" && \$4=="T" && \$10=="C" && \$11=="A") || (\$3=="C" && \$4=="A" && \$10=="A" && \$11=="G") || (\$3=="C" && \$4=="A" && \$10=="G" && \$11=="A") || (\$3=="C" && \$4=="A" && \$10=="A" && \$11=="T") || (\$3=="C" && \$4=="A" && \$10=="T" && \$11=="A") || (\$3=="C" && \$4=="A" && \$10=="T" && \$11=="C") || (\$3=="C" && \$4=="A" && \$10=="C" && \$11=="T") || (\$3=="C" && \$4=="A" && \$10=="G" && \$11=="C") || (\$3=="C" && \$4=="A" && \$10=="C" && \$11=="G") || (\$3=="C" && \$4=="G" && \$10=="G" && \$11=="A") || (\$3=="C" && \$4=="G" && \$10=="A" && \$11=="G") || (\$3=="C" && \$4=="G" && \$10=="C" && \$11=="A") || (\$3=="C" && \$4=="G" && \$10=="A" && \$11=="C") || (\$3=="C" && \$4=="G" && \$10=="C" && \$11=="T") || (\$3=="C" && \$4=="G" && \$10=="T" && \$11=="C") || (\$3=="C" && \$4=="G" && \$10=="G" && \$11=="T") || (\$3=="C" && \$4=="G" && \$10=="T" && \$11=="G") || (\$3=="C" && \$4=="T" && \$10=="A" && \$11=="C") || (\$3=="C" && \$4=="T" && \$10=="C" && \$11=="A") || (\$3=="C" && \$4=="T" && \$10=="T" && \$11=="A") || (\$3=="C" && \$4=="T" && \$10=="A" && \$11=="T") || (\$3=="C" && \$4=="T" && \$10=="T" && \$11=="G") || (\$3=="C" && \$4=="T" && \$10=="G" && \$11=="T") || (\$3=="C" && \$4=="T" && \$10=="C" && \$11=="G") || (\$3=="C" && \$4=="T" && \$10=="G" && \$11=="C") || (\$3=="G" && \$4=="A" && \$10=="A" && \$11=="C") || (\$3=="G" && \$4=="A" && \$10=="C" && \$11=="A") || (\$3=="G" && \$4=="A" && \$10=="A" && \$11=="T") || (\$3=="G" && \$4=="A" && \$10=="T" && \$11=="A") || (\$3=="G" && \$4=="A" && \$10=="G" && \$11=="T") || (\$3=="G" && \$4=="A" && \$10=="T" && \$11=="G") || (\$3=="G" && \$4=="A" && \$10=="G" && \$11=="G") || (\$3=="G" && \$4=="A" && \$10=="C" && \$11=="G") || (\$3=="G" && \$4=="T" && \$10=="A" && \$11=="G") || (\$3=="G" && \$4=="T" && \$10=="G" && \$11=="A") || (\$3=="G" && \$4=="T" && \$10=="T" && \$11=="A") || (\$3=="G" && \$4=="T" && \$10=="A" && \$11=="T") || (\$3=="G" && \$4=="T" && \$10=="T" && \$11=="C") || (\$3=="G" && \$4=="T" && \$10=="C" && \$11=="T") || (\$3=="G" && \$4=="T" && \$10=="T" && \$11=="A") || (\$3=="G" && \$4=="T" && \$10=="A" && \$11=="T") || (\$3=="G" && \$4=="C" && \$10=="G" && \$11=="T") || (\$3=="G" && \$4=="C" && \$10=="T" && \$11=="G") || (\$3=="G" && \$4=="C" && \$10=="G" && \$11=="A") || (\$3=="G" && \$4=="C" && \$10=="A" && \$11=="G") || (\$3=="G" && \$4=="C" && \$10=="C" && \$11=="T") || (\$3=="G" && \$4=="C" && \$10=="T" && \$11=="C") || (\$3=="G" && \$4=="C" && \$10=="A" && \$11=="T") || (\$3=="G" && \$4=="C" && \$10=="T" && \$11=="A") || (\$3=="T" && \$4=="A" && \$10=="A" && \$11=="G") || (\$3=="T" && \$4=="A" && \$10=="G" && \$11=="A") || (\$3=="T" && \$4=="A" && \$10=="C" && \$11=="G") || (\$3=="T" && \$4=="A" && \$10=="G" && \$11=="C") || (\$3=="T" && \$4=="A" && \$10=="T" && \$11=="G") || (\$3=="T" && \$4=="A" && \$10=="G" && \$11=="T") || (\$3=="T" && \$4=="A" && \$10=="T" && \$11=="C") || (\$3=="T" && \$4=="A" && \$10=="C" && \$11=="T") || (\$3=="T" && \$4=="C" && \$10=="A" && \$11=="C") || (\$3=="T" && \$4=="C" && \$10=="C" && \$11=="A") || (\$3=="T" && \$4=="C" && \$10=="C" && \$11=="G") || (\$3=="T" && \$4=="C" && \$10=="G" && \$11=="C") || (\$3=="T" && \$4=="C" && \$10=="G" && \$11=="T") || (\$3=="T" && \$4=="C" && \$10=="T" && \$11=="G") || (\$3=="T" && \$4=="C" && \$10=="A" && \$11=="T") || (\$3=="T" && \$4=="C" && \$10=="T" && \$11=="A") || (\$3=="T" && \$4=="G" && \$10=="C" && \$11=="T") || (\$3=="T" && \$4=="G" && \$10=="T" && \$11=="C") || (\$3=="T" && \$4=="G" && \$10=="G" && \$11=="A") || (\$3=="T" && \$4=="G" && \$10=="A" && \$11=="G") || (\$3=="T" && \$4=="G" && \$10=="G" && \$11=="C") || (\$3=="T" && \$4=="G" && \$10=="C" && \$11=="G") || (\$3=="T" && \$4=="G" && \$10=="T" && \$11=="A") || (\$3=="T" && \$4=="G" && \$10=="A" && \$11=="T") )      
  {print \$2 > "temp_non_biallelic_snps"};
  if ( (\$3== "I" && \$4=="D" && (length(\$10)<length(\$11))   )  || (\$3== "D" && \$4=="I" && (length(\$10)>length(\$11)) ) && ((\$5 - \$12) < DIFF ) && ((\$5 - \$12) > -(DIFF)  )  )
      {print \$2 > "temp_indels_aligned_same_allele"};
  if ( (\$3== "I" && \$4=="D" && (length(\$10)>length(\$11))   )  || (\$3== "D" && \$4=="I" && (length(\$10)<length(\$11)) ) && ((\$5 - (1 - \$12)) < DIFF ) && ((\$5 - (1- \$12)) > -(DIFF)  )  )
      {print \$2 > "temp_indels_aligned_other_allele"};
}
END {}
__EOF__
) > ${dirI}/awk_check_based_on_ID_1KG; chmod 744 ${dirI}/awk_check_based_on_ID_1KG

#the following function compares SNPs matched by position (CHR, BP) and check their IDs, and if different, works out if they are likely the same SNP. Also works out flips etc
  #Next thing to to fix the reverse - same position, and alleles, but different ID. This is not necessary for IMPUTE2 (as it does that for you, since it only cares about position) but might been needed for HRC, and will also highlight duplicate SNPs in the data (same SNP under different names). Should probably test for those first within each set, and will need to do again after merging. Oh wait I remember this is a pain to do as it fails each time it hit an unfixable pair, even if there might be many. Hmmm, in SDH many of these seem to have 0,0 for positions, or for alleles

(cat <<__EOF__
#!/usr/bin/awk -f
# written by script_to_clean_QC_ausmel_GWAS_data_2016.sh.
# takes SNP allele and freq data paired with 1kG data that has the same position and works out if they are the same SNP, and aligned
BEGIN {}
{
    if( \$2==\$14 && !(\$3=="A" && \$4=="T") && !(\$3=="T" && \$4=="A") && !(\$3=="C" && \$4=="G") && !(\$3=="G" && \$4=="C") && ( (\$3==\$10 && \$4==\$11) || (\$3==\$11 && \$4==\$10 ) ) )
      {print \$2,\$14 > "temp_aligned_ID_concordant"};
    if( \$2!=\$14 && !(\$3=="A" && \$4=="T") && !(\$3=="T" && \$4=="A") && !(\$3=="C" && \$4=="G") && !(\$3=="G" && \$4=="C") && ( (\$3==\$10 && \$4==\$11) || (\$3==\$11 && \$4==\$10 ) ) )
      {print \$2,\$14 > "temp_aligned_ID_disc"};
    if( \$2==\$14 && ((\$3=="G" && \$4=="A" && \$10=="C" && \$11=="T" ) || (\$3=="G" && \$4=="A" && \$10=="T" && \$11=="C" ) || (\$3=="G" && \$4=="T" && \$10=="C" && \$11=="A" ) || (\$3=="G" && \$4=="T" && \$10=="A" && \$11=="C") || (\$3=="C" && \$4=="A" && \$10=="G" && \$11=="T") || (\$3=="C" && \$4=="A" && \$10=="T" && \$11=="G" ) || (\$3=="C" && \$4=="T" && \$10=="G" && \$11=="A") || (\$3=="C" && \$4=="T" && \$10=="A" && \$11=="G") || (\$3=="A" && \$4=="G" && \$10=="C" && \$11=="T" ) || (\$3=="A" && \$4=="G" && \$10=="T" && \$11=="C" ) || (\$3=="T" && \$4=="G" && \$10=="C" && \$11=="A" ) || (\$3=="T" && \$4=="G" && \$10=="A" && \$11=="C" ) || (\$3=="A" && \$4=="C" && \$10=="G" && \$11=="T" ) || (\$3=="A" && \$4=="C" && \$10=="T" && \$11=="G" ) || (\$3=="T" && \$4=="C" && \$10=="G" && \$11=="A" ) || (\$3=="T" && \$4=="C" && \$10=="A" && \$11=="G")) ) 
      {print \$2,\$14 > "temp_flipped_ID_concordant"};
    if( \$2!=\$14 && ( (\$3=="G" && \$4=="A" && \$10=="C" && \$11=="T" ) || (\$3=="G" && \$4=="A" && \$10=="T" && \$11=="C" ) || (\$3=="G" && \$4=="T" && \$10=="C" && \$11=="A" ) || (\$3=="G" && \$4=="T" && \$10=="A" && \$11=="C") || (\$3=="C" && \$4=="A" && \$10=="G" && \$11=="T") || (\$3=="C" && \$4=="A" && \$10=="T" && \$11=="G" ) || (\$3=="C" && \$4=="T" && \$10=="G" && \$11=="A") || (\$3=="C" && \$4=="T" && \$10=="A" && \$11=="G") || (\$3=="A" && \$4=="G" && \$10=="C" && \$11=="T" ) || (\$3=="A" && \$4=="G" && \$10=="T" && \$11=="C" ) || (\$3=="T" && \$4=="G" && \$10=="C" && \$11=="A" ) || (\$3=="T" && \$4=="G" && \$10=="A" && \$11=="C" ) || (\$3=="A" && \$4=="C" && \$10=="G" && \$11=="T" ) || (\$3=="A" && \$4=="C" && \$10=="T" && \$11=="G" ) || (\$3=="T" && \$4=="C" && \$10=="G" && \$11=="A" ) || (\$3=="T" && \$4=="C" && \$10=="A" && \$11=="G") )) 
      {print \$2,\$14 > "temp_flipped_ID_disc"};
  if( ( (\$3=="A" && \$4=="T" && \$10=="T" && \$11=="A") || (\$3=="T" && \$4=="A" && \$10=="A" && \$11=="T") || (\$3=="A" && \$4=="T" && \$10=="A" && \$11=="T") || (\$3=="T" && \$4=="A" && \$10=="T" && \$11=="A") || (\$3=="G" && \$4=="C" && \$10=="C" && \$11=="G") || (\$3=="C" && \$4=="G" && \$10=="G" && \$11=="C")   || (\$3=="G" && \$4=="C" && \$10=="G" && \$11=="C") || (\$3=="C" && \$4=="G" && \$10=="C" && \$11=="G") )  && \$5>= MAF)
      {print \$2,\$14 > "temp_ambigious_ID_MAF_to_high_too_assess"};
  if( \$2==\$14 && ( (\$3=="A" && \$4=="T" && \$10=="T" && \$11=="A") || (\$3=="T" && \$4=="A" && \$10=="A" && \$11=="T") || (\$3=="G" && \$4=="C" && \$10=="C" && \$11=="G") || (\$3=="C" && \$4=="G" && \$10=="G" && \$11=="C") ) && \$5 < MAF && ((\$5 - \$12) < DIFF ) && ((\$5 - \$12) > -(DIFF) )  )
      {print \$2,\$14 > "temp_ambigious_aligned_ID_concordant_same_allele"};
  if( \$2!=\$14 && ( (\$3=="A" && \$4=="T" && \$10=="T" && \$11=="A") || (\$3=="T" && \$4=="A" && \$10=="A" && \$11=="T") || (\$3=="G" && \$4=="C" && \$10=="C" && \$11=="G") || (\$3=="C" && \$4=="G" && \$10=="G" && \$11=="C") ) && \$5 < MAF && ((\$5 - \$12) < DIFF ) && ((\$5 - \$12) > -(DIFF) )  )
      {print \$2,\$14 > "temp_ambigious_aligned_ID_disc_same_allele"};
  if( \$2==\$14 && ( (\$3=="A" && \$4=="T" && \$10=="A" && \$11=="T") || (\$3=="T" && \$4=="A" && \$10=="T" && \$11=="A") || (\$3=="G" && \$4=="C" && \$10=="G" && \$11=="C") || (\$3=="C" && \$4=="G" && \$10=="C" && \$11=="G") ) && \$5 < MAF && ((\$5 - (1 - \$12)) < DIFF ) && ((\$5 - (1-\$12)) > -(DIFF) )    )
      {print \$2,\$14 > "temp_ambigious_aligned_ID_concordant_other_allele"};
  if( \$2!=\$14 && ( (\$3=="A" && \$4=="T" && \$10=="A" && \$11=="T") || (\$3=="T" && \$4=="A" && \$10=="T" && \$11=="A") || (\$3=="G" && \$4=="C" && \$10=="G" && \$11=="C") || (\$3=="C" && \$4=="G" && \$10=="C" && \$11=="G") ) && \$5 < MAF && ((\$5 - (1 - \$12)) < DIFF ) && ((\$5 - (1-\$12)) > -(DIFF) )    )
      {print \$2,\$14 > "temp_ambigious_aligned_ID_disc_other_allele"};
  if( \$2==\$14 && ( (\$3=="A" && \$4=="T" && \$10=="A" && \$11=="T") || (\$3=="T" && \$4=="A" && \$10=="T" && \$11=="A") || (\$3=="G" && \$4=="C" && \$10=="G" && \$11=="C") || (\$3=="C" && \$4=="G" && \$10=="C" && \$11=="G") )  && \$5 < MAF && ((\$5 - \$12) < DIFF ) && ((\$5 - \$12) > -(DIFF) )   )
      {print \$2,\$14 > "temp_ambigious_flipped_ID_concordant_same_allele"};
  if( \$2!=\$14 && ( (\$3=="A" && \$4=="T" && \$10=="A" && \$11=="T") || (\$3=="T" && \$4=="A" && \$10=="T" && \$11=="A") || (\$3=="G" && \$4=="C" && \$10=="G" && \$11=="C") || (\$3=="C" && \$4=="G" && \$10=="C" && \$11=="G") )  && \$5 < MAF && ((\$5 - \$12) < DIFF ) && ((\$5 - \$12) > -(DIFF) )   )
      {print \$2,\$14 > "temp_ambigious_flipped_ID_disc_same_allele"};
  if( \$2==\$14 && ( (\$3=="A" && \$4=="T" && \$10=="T" && \$11=="A") || (\$3=="T" && \$4=="A" && \$10=="A" && \$11=="T") || (\$3=="G" && \$4=="C" && \$10=="C" && \$11=="G") || (\$3=="C" && \$4=="G" && \$10=="G" && \$11=="C") ) && \$5 < MAF && ((\$5 - (1 - \$12)) < DIFF ) && ((\$5 - (1-\$12)) > -(DIFF) )  )
      {print \$2,\$14 > "temp_ambigious_flipped_ID_concordant_other_allele"};
  if( \$2!=\$14 && ( (\$3=="A" && \$4=="T" && \$10=="T" && \$11=="A") || (\$3=="T" && \$4=="A" && \$10=="A" && \$11=="T") || (\$3=="G" && \$4=="C" && \$10=="C" && \$11=="G") || (\$3=="C" && \$4=="G" && \$10=="G" && \$11=="C") ) && \$5 < MAF && ((\$5 - (1 - \$12)) < DIFF ) && ((\$5 - (1-\$12)) > -(DIFF) )  )
      {print \$2,\$14 > "temp_ambigious_flipped_ID_disc_other_allele"};
  if ( (\$3=="A" && \$4=="C" && \$10=="C" && \$11=="G") || (\$3=="A" && \$4=="C" && \$10=="G" && \$11=="C") || (\$3=="A" && \$4=="C" && \$10=="C" && \$11=="T") || (\$3=="A" && \$4=="C" && \$10=="T" && \$11=="C") || (\$3=="A" && \$4=="C" && \$10=="G" && \$11=="A") || (\$3=="A" && \$4=="C" && \$10=="A" && \$11=="G") || (\$3=="A" && \$4=="C" && \$10=="T" && \$11=="A") || (\$3=="A" && \$4=="C" && \$10=="A" && \$11=="T") || (\$3=="A" && \$4=="G" && \$10=="A" && \$11=="C") || (\$3=="A" && \$4=="G" && \$10=="C" && \$11=="A") || (\$3=="A" && \$4=="G" && \$10=="C" && \$11=="G") || (\$3=="A" && \$4=="G" && \$10=="G" && \$11=="C") || (\$3=="A" && \$4=="G" && \$10=="G" && \$11=="T") || (\$3=="A" && \$4=="G" && \$10=="T" && \$11=="G") || (\$3=="A" && \$4=="G" && \$10=="A" && \$11=="T") || (\$3=="A" && \$4=="G" && \$10=="T" && \$11=="A") || (\$3=="A" && \$4=="T" && \$10=="T" && \$11=="C") || (\$3=="A" && \$4=="T" && \$10=="C" && \$11=="T") || (\$3=="A" && \$4=="T" && \$10=="A" && \$11=="G") || (\$3=="A" && \$4=="T" && \$10=="G" && \$11=="A") || (\$3=="A" && \$4=="T" && \$10=="C" && \$11=="G") || (\$3=="A" && \$4=="T" && \$10=="G" && \$11=="C") || (\$3=="A" && \$4=="T" && \$10=="A" && \$11=="C") || (\$3=="A" && \$4=="T" && \$10=="C" && \$11=="A") || (\$3=="C" && \$4=="A" && \$10=="A" && \$11=="G") || (\$3=="C" && \$4=="A" && \$10=="G" && \$11=="A") || (\$3=="C" && \$4=="A" && \$10=="A" && \$11=="T") || (\$3=="C" && \$4=="A" && \$10=="T" && \$11=="A") || (\$3=="C" && \$4=="A" && \$10=="T" && \$11=="C") || (\$3=="C" && \$4=="A" && \$10=="C" && \$11=="T") || (\$3=="C" && \$4=="A" && \$10=="G" && \$11=="C") || (\$3=="C" && \$4=="A" && \$10=="C" && \$11=="G") || (\$3=="C" && \$4=="G" && \$10=="G" && \$11=="A") || (\$3=="C" && \$4=="G" && \$10=="A" && \$11=="G") || (\$3=="C" && \$4=="G" && \$10=="C" && \$11=="A") || (\$3=="C" && \$4=="G" && \$10=="A" && \$11=="C") || (\$3=="C" && \$4=="G" && \$10=="C" && \$11=="T") || (\$3=="C" && \$4=="G" && \$10=="T" && \$11=="C") || (\$3=="C" && \$4=="G" && \$10=="G" && \$11=="T") || (\$3=="C" && \$4=="G" && \$10=="T" && \$11=="G") || (\$3=="C" && \$4=="T" && \$10=="A" && \$11=="C") || (\$3=="C" && \$4=="T" && \$10=="C" && \$11=="A") || (\$3=="C" && \$4=="T" && \$10=="T" && \$11=="A") || (\$3=="C" && \$4=="T" && \$10=="A" && \$11=="T") || (\$3=="C" && \$4=="T" && \$10=="T" && \$11=="G") || (\$3=="C" && \$4=="T" && \$10=="G" && \$11=="T") || (\$3=="C" && \$4=="T" && \$10=="C" && \$11=="G") || (\$3=="C" && \$4=="T" && \$10=="G" && \$11=="C") || (\$3=="G" && \$4=="A" && \$10=="A" && \$11=="C") || (\$3=="G" && \$4=="A" && \$10=="C" && \$11=="A") || (\$3=="G" && \$4=="A" && \$10=="A" && \$11=="T") || (\$3=="G" && \$4=="A" && \$10=="T" && \$11=="A") || (\$3=="G" && \$4=="A" && \$10=="G" && \$11=="T") || (\$3=="G" && \$4=="A" && \$10=="T" && \$11=="G") || (\$3=="G" && \$4=="A" && \$10=="G" && \$11=="G") || (\$3=="G" && \$4=="A" && \$10=="C" && \$11=="G") || (\$3=="G" && \$4=="T" && \$10=="A" && \$11=="G") || (\$3=="G" && \$4=="T" && \$10=="G" && \$11=="A") || (\$3=="G" && \$4=="T" && \$10=="T" && \$11=="A") || (\$3=="G" && \$4=="T" && \$10=="A" && \$11=="T") || (\$3=="G" && \$4=="T" && \$10=="T" && \$11=="C") || (\$3=="G" && \$4=="T" && \$10=="C" && \$11=="T") || (\$3=="G" && \$4=="T" && \$10=="T" && \$11=="A") || (\$3=="G" && \$4=="T" && \$10=="A" && \$11=="T") || (\$3=="G" && \$4=="C" && \$10=="G" && \$11=="T") || (\$3=="G" && \$4=="C" && \$10=="T" && \$11=="G") || (\$3=="G" && \$4=="C" && \$10=="G" && \$11=="A") || (\$3=="G" && \$4=="C" && \$10=="A" && \$11=="G") || (\$3=="G" && \$4=="C" && \$10=="C" && \$11=="T") || (\$3=="G" && \$4=="C" && \$10=="T" && \$11=="C") || (\$3=="G" && \$4=="C" && \$10=="A" && \$11=="T") || (\$3=="G" && \$4=="C" && \$10=="T" && \$11=="A") || (\$3=="T" && \$4=="A" && \$10=="A" && \$11=="G") || (\$3=="T" && \$4=="A" && \$10=="G" && \$11=="A") || (\$3=="T" && \$4=="A" && \$10=="C" && \$11=="G") || (\$3=="T" && \$4=="A" && \$10=="G" && \$11=="C") || (\$3=="T" && \$4=="A" && \$10=="T" && \$11=="G") || (\$3=="T" && \$4=="A" && \$10=="G" && \$11=="T") || (\$3=="T" && \$4=="A" && \$10=="T" && \$11=="C") || (\$3=="T" && \$4=="A" && \$10=="C" && \$11=="T") || (\$3=="T" && \$4=="C" && \$10=="A" && \$11=="C") || (\$3=="T" && \$4=="C" && \$10=="C" && \$11=="A") || (\$3=="T" && \$4=="C" && \$10=="C" && \$11=="G") || (\$3=="T" && \$4=="C" && \$10=="G" && \$11=="C") || (\$3=="T" && \$4=="C" && \$10=="G" && \$11=="T") || (\$3=="T" && \$4=="C" && \$10=="T" && \$11=="G") || (\$3=="T" && \$4=="C" && \$10=="A" && \$11=="T") || (\$3=="T" && \$4=="C" && \$10=="T" && \$11=="A") || (\$3=="T" && \$4=="G" && \$10=="C" && \$11=="T") || (\$3=="T" && \$4=="G" && \$10=="T" && \$11=="C") || (\$3=="T" && \$4=="G" && \$10=="G" && \$11=="A") || (\$3=="T" && \$4=="G" && \$10=="A" && \$11=="G") || (\$3=="T" && \$4=="G" && \$10=="G" && \$11=="C") || (\$3=="T" && \$4=="G" && \$10=="C" && \$11=="G") || (\$3=="T" && \$4=="G" && \$10=="T" && \$11=="A") || (\$3=="T" && \$4=="G" && \$10=="A" && \$11=="T") )      
      {print \$2,\$14 > "temp_non_biallelic_snps_ID"};
  if ( (\$3== "I" && \$4=="D" && (length(\$10)<length(\$11))   )  || (\$3== "D" && \$4=="I" && (length(\$10)>length(\$11)) ) && ((\$5 - \$12) < DIFF ) && ((\$5 - \$12) > -(DIFF)  )  )
      {print \$2,\$14 > "temp_indels_ID_aligned_same_allele"};
  if ( (\$3== "I" && \$4=="D" && (length(\$10)>length(\$11))   )  || (\$3== "D" && \$4=="I" && (length(\$10)<length(\$11)) ) && ((\$5 - (1 - \$12)) < DIFF ) && ((\$5 - (1- \$12)) > -(DIFF)  )  )
      {print \$2,\$14 > "temp_indels_ID_aligned_other_allele"};
  }
END {}
__EOF__
) > ${dirI}/awk_check_based_on_pos_1KG; chmod 744 ${dirI}/awk_check_based_on_pos_1KG


########################################################################
#0.5 Manual SNP remapping
########################################################################
   #25/07/2016 found a 2 pairs of SNPs where the liftover positions are wrong. Looks to be where to of the same kind of SNPs (e.g. both A/G) with sim freqs are next to each other; the chain has assigned one to the wrong positon, and dropped the other. This will fix that if in the dataset. 27/07/2016 I wonder if it more that these positions ar wrong in the initial bim, and instead liftover just adds bp to convert from hg18-hg19. For OAG2 there are a bunch of SNPs that are wrong by 1b[ in the initial bim, and wrong by 1 post liftover. Might be better to have this step after the liftover to harmonise all steps. Following list accrued from various attempts to merge files etc. On 13/07/2016 found ~16 multiple position SNPs that aren't in 1KG data I have. Looked them up in dbSNP manually, and where possible gave them the hg19 position. A few were flagged as suspect or had some other oddity - rs11930373 has two positions in dbSNP - as in 66358713:66358714 which I have never seen before, so were removed. I was correcting these just in omni (as they had only shown up when merged with SDH but they are wrong in AMFS as well, so fix it now for all
  #Further run throughs identified additional SNPs on the arrays, but not in 1KG ref files, with different positions between SDH and QMEGA_omni; a number had been zeroed out. Where they could be matched by position/freq to real SNPs in dbSNP added and initial step at the start of the script to remap them otherwise exclude them. Some requiring flipping still (as not in 1kG ref files I have can't align them to a set strand.

##PG can use the folowing to remap as well (even if the map being correct in my dataset, it won't cause a problem since the map will remain unged in such case).

cd $dirI

echo "rs11174810 63540295" > fix_remap    #OAG2 map position wrong by 1bp
echo "rs34475990 63540296" >> fix_remap   #OAG2 map position wrong by 1bp
echo "rs17233804 89804723" >> fix_remap #OAG2 map position wrong by 1bp
echo "rs17227417 89804722" >> fix_remap #OAG2 map position wrong by 1bp
echo "rs1166587 78256947" >> fix_remap #wrong position in QMEGA_omni/AMFS files
echo "rs11936434 103715816" >> fix_remap #wrong position in QMEGA_omni/AMFS files
echo "rs36058833 33689482" >> fix_remap #wrong position in QMEGA_omni/AMFS files
echo "rs17879125 32551892" >> fix_remap #wrong position in QMEGA_omni/AMFS files
echo "rs7382836 32734314" >> fix_remap #wrong position in QMEGA_omni/AMFS files
echo "rs10867964 85861561" >> fix_remap #wrong position in QMEGA_omni/AMFS files
echo "rs11594913 39103082" >> fix_remap #wrong position in QMEGA_omni/AMFS files
echo "rs2155163 126524251" >> fix_remap #wrong position in QMEGA_omni/AMFS files
echo "rs2596105 17772069" >> fix_remap #wrong position in QMEGA_omni/AMFS files
echo "rs2266856 151141532" >> fix_remap #wrong position in QMEGA_omni/AMFS files
echo "rs1052571 15850613" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs17887074 22964176" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs2501428 24256932" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs2501430 24264244" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs7534633 24280063" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs7530595 24281106" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs9424346 24286885" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs6696719 24286951" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs4648925 24292830" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs9424410 24292854" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs4276860 24294249" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs6667686 24300849" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs6674323 24302765" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs4532777 24306054" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs4474201 24335104" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs10493019 24336524" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs7553231 24347103" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs6700106 24357576" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs11583297 26961834" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs11580952 27113923" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs910048 27215784" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs35647061 27217699" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs17362452 28157308" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs648759 29092785" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs6666421 30884114" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs4949495 32644104" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs34149698 32717536" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs12402882 33978546" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs1404350 34764974" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs2001943 37205500" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs577674 72309009" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs4147851 94478784" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs17181981 161015268" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs10917701 163191314" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs3830180 167764925" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs7880493 171061033" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs28363522 171061034" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs12758408 171170971" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs598562 172164499" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs431188 241117850" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs17008384 29867599" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs11466299 70675309" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs9288700 158735265" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs2893113 169815139" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs10206282 178883104" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs13013663 219667194" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs17868342 234673516" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs4135260 12422488" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs3773656 30722653" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs1839778 32165335" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs11712410 57892031" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs12721615 119537223" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs814130 127410959" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs1584576 188911398" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs6856283 3602109" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs11936089 18011864" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs11098602 122052641" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs406542 171510501" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs11942861 179496461" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs326121 7876288" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs4700177 66545138" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs387528 6172680" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs5743651 6649867" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs6938256 20887954" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs6929563 27472223" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs2516703 30229922" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs2524212 30345871" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs1049281 31236567" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs9265797 31308630" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs2596487 31325056" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs4645843 31544562" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs34053064 31781748" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs4947331 31821041" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs28435656 31880637" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs12182351 32201707" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs2187823 32439508" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs7748472 32448763" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs7749057 32448904" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs3830135 32548464" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs9270984 32573991" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs9272143 32600803" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs213207 33241169" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs3025051 43748439" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs9476225 57535593" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs17594297 91470169" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs1084651 161089817" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs4252195 161159619" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs1182136 2966445" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs4722357 24442911" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs28746498 87286985" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs3917561 94932674" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs17882992 100486016" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs7804904 123592946" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs1582264 125717914" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs6138 139572123" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs10247107 150683083" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs3918232 150706640" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs2069444 150754699" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs4149238 27361352" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs4745571 79318367" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs323729 102116594" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs10980418 113273768" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs4979603 118738038" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs4842091 138947167" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs1194492 53966300" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs16928461 54013275" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs2227579 75671289" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs1107688 92495336" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs3824739 92681475" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs8187709 101611277" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs7894236 122236760" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs3793971 1472605" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs12293627 3704770" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs11031961 4389110" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs11023292 14670923" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs12792463 14888625" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs4757587 17760904" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs12786539 34486626" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs11231174 62426102" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs17881532 64956204" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs3218740 94212019" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs11212587 108186610" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs948413 119214771" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs3824995 134158745" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs10878735 68447643" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs776417 70022135" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs17110563 72366306" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs4309196 72740262" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs12821081 72977698" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs1275646 76150495" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs2694820 77815771" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs1994667 77880974" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs7487519 81315248" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs10466971 88496410" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs10745480 88787395" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs10859474 93687153" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs1616091 98890053" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs2374654 107371180" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs1426455 107715409" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs4759383 123215850" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs1296265 133511592" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs2000258 36533992" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs2485287 38519281" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs9603453 39430879" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs4142536 39941628" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs1179957 42077774" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs8015278 64875004" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs17841074 72667671" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs10134356 91365647" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs7142236 93011877" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs7152495 102763731" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs11630521 75645538" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs2305364 85452275" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs4932475 89689583" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs4984951 917120" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs4148375 16212302" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs3743760 77227831" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs17233833 89804330" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs12941873 8068763" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs12051754 8088900" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs17772089 21198118" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs6587138 21214253" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs8075185 36060216" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs4792939 42425313" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs11568589 48761363" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs6502001 67096079" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs6501576 70897142" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs12939354 74106612" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs16970386 75625282" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs4790048 77478072" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs11650413 78228181" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs2446199 5765436" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs5030387 10380282" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs11575934 18186618" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs17884680 36630803" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs12973764 44023153" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs17207376 55248072" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs1654644 55373362" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs2875839 3149729" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs713164 45895799" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs5584 48124462" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs6068812 52774635" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs2833467 32964287" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs17860248 34602261" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs7283856 34959776" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs2190742 17254399" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs16982871 19965563" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs13054858 25424579" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs12330067 25425439" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs5760813 25445405" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs742169 25448604" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs2413140 33051129" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs470033 33071432" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs926755 36134690" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs2157250 36631691" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs738304 37902177" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs6009465 49315035" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs2301416 55028863" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs17875899 70330483" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs17423 153629663" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs36101366 154064493" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs2534636 2657176" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs2058276 2668456" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs2267801 2828196" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs1865680 6868118" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs9306845 6941218" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs7893107 6995523" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs9786714 7173143" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs9785717 7771131" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs9785740 8334875" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs17316729 8353707" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs17842387 8424089" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs9306848 8474189" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs7892893 8590752" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs9786489 9851457" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs2740981 14433100" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs9786043 15472863" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs17307070 15590342" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs34134567 15604899" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs16980360 15809326" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs17316592 18560005" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs35547782 18700150" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs877756 18719565" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs9786197 21117888" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs9786232 21166358" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs2032640 21892572" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs2032652 21917313" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs16980426 22214221" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs16980548 22918577" >> fix_remap #Position wrong by 1bp in OAG2, correct in OAG1, not in 1kG but in dbSNP
echo "rs9341257 38298744" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs17225060 47710314" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs13086387 119844557" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs17808817 187763947" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs10937669 5771285" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs11569705 70723299" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs2910160 156795719" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs3004047 26710391" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs2596454 31436312" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs2262697 152516647" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs34160508 99300485" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs11603505 4400545" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs35676648 69486514" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs28399546 92537466" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs11638512 25476983" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs8051412 55670725" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs12602109 21216602" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs11871216 21236058" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs447802 49116359" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs17207383 55248107" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs3865507 55322376" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs12461010 55332601" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs5744203 36979250" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs28368104 55918846" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs28359542 56137758" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs28624811 42527533" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs34738655 2892150" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs9786139 6753519" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs16981290 7568568" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs895530 7963031" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs3853054 8679843" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs4141564 8685083" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs16981293 8796078" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs16980473 14159846" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs9786712 14577177" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs9786191 14804077" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs2032601 14869076" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs2032600 14888783" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs20320 14898163" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs20321 14902414" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs8179022 15018459" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs8179021 15018582" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs9341290 15020578" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs9786025 15021522" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs9341301 15023364" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs16980601 15415115" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs2032667 15436316" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs17221964 15517851" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs917759 15754313" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs7067279 15935524" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs9786893 16846439" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs16980459 16906683" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs16980499 17256018" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs17307398 17285993" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs9786420 17398598" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs17222419 17508337" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs9786076 17844018" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs17222573 17891241" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs9785905 18097251" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs9785659 18248698" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs9786283 18907236" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs17315758 18964263" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs17315772 19038302" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs13304625 19077409" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs17315821 19166861" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs9786119 19198212" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs9786357 19500107" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs3848982 21717208" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs3865828 22094491" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs34283263 22178569" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs2196155 22665262" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs9341308 22738775" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs1558843 22750583" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs13447357 22751863" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs34126399 23021978" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs17842518 23443971" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs10106770 235832763" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs12496398 115156497" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs2569201 109657233" >> fix_remap # Position wrong in OAG1/OAG2, correct in IBD which matches dbSNP
echo "rs2032629 21865821" >> fix_remap #position wrong in OAG by 1bp, correct in IBD/WAMHS
echo "rs12189802 27714048" >> fix_remap #position wrong in OAG by 1bp, correct in WAMHS
echo "rs200958 27840826" >> fix_remap #position wrong in OAG by 1bp, correct in WAMHS

#plus some ones to remove as in the AMFS/omni data but have wrong position and flagged as unrelaible in dbSNP. Also one in omni/SDH - chr0 in SDH, quad allelic on dbSNP
echo rs11930373 > fix_exclude
echo rs35276574 >> fix_exclude
echo rs1523844 >> fix_exclude
echo rs361359 >> fix_exclude
echo rs2431622 >> fix_exclude
echo rs35288367 >> fix_exclude
echo rs2413140 >> fix_exclude #freq matches between SDH and qmega_omni, but chr0 in SDH. Quad allellic on dbSNP, drop
echo rs12662611 >> fix_exclude #according to dbSNP same SNP as rs114224970; two assays on arrays through with wildly different MAFs
echo rs114224970 >> fix_exclude #according to dbSNP same SNP as rs12662611; two assays on OAG arrays through with wildly different MAFs
echo rs12043679 >> fix_exclude #diff chr positions in OAG, dbSNP, when checked dbSNP was quad allelic
echo rs10493784 >> fix_exclude #same position as exm2275458 but diff allele; site is triplicate in dbSNP so easier to exclude
echo exm2275458 >> fix_exclude #same position as rs10493784 but diff allele; site is triplicate in dbSNP so easier to exclude
echo rs10487391 >> fix_exclude #same position as exm2137387 but diff allele; site is triplicate in dbSNP so easier to exclude
echo exm2137387 >> fix_exclude #same position as rs10487391 but diff allele; site is triplicate in dbSNP so easier to exclude
echo rs7910455 >> fix_exclude #same position as exm2276793 but diff allele; site is triplicate in dbSNP so easier to exclude
echo exm2276793 >> fix_exclude #same position as rs7910455 but diff allele; site is triplicate in dbSNP so easier to exclude
echo exm111624 >> fix_exclude #exm111624,exm2277003,rs140418584 are the same very rare quad SNP in IBD (and probably WAMHS), easier to exclude
echo exm2277003 >> fix_exclude #exm111624,exm2277003,rs140418584 are the same very rare quad SNP in IBD (and probably WAMHS), easier to exclude
echo rs140418584 >> fix_exclude #exm111624,exm2277003,rs140418584 are the same very rare quad SNP in IBD (and probably WAMHS), easier to exclude
echo rs9505819 >> fix_exclude #on most arrays but positions differ in OAG vs others; but no dbSNP validation info so may be suspect, dropping
echo rs5942313 >> fix_exclude #chr23 but giving het in males in WAMHS dataset
echo exm2276927 >> fix_exclude #on exm array, tri allelic
echo exm1669027 >> fix_exclude #on exm array, tri allelic

#few chr errors as well - 0 in SDH and other sets but freq matches case/control/dbSNP so happy with 
echo "rs2190742 22" > fix_remap_chr #freq matches, is chr22 on dbSNP with ~0.3 freq
echo "rs13054858 22" >> fix_remap_chr #freq matches, is chr22 on dbSNP with ~0.2 freq
echo "rs12330067 22" >> fix_remap_chr #freq matches, is chr22 on dbSNP with ~0.35 freq
echo "rs5760813 22" >> fix_remap_chr #freq matches, is chr22 on dbSNP with ~0.35 freq
echo "rs742169 22" >> fix_remap_chr #freq matches, is chr22 on dbSNP with ~0.18 freq
echo "rs470033 22" >> fix_remap_chr #freq matches, is chr22 on dbSNP with ~0.1 freq
echo "rs1296265 12" >> fix_remap_chr #freq matches, is chr22 on dbSNP with ~0.024 freq
echo "rs10106770 2" >> fix_remap_chr #freq matches between OAG and IBD, IBD chr matches dbSNP (changed chr between hapmap and 1KG)
echo "rs12496398 4" >> fix_remap_chr #freq matches between OAG and IBD, IBD chr matches dbSNP (changed chr between hapmap and 1KG)
echo "rs2569201 8" >> fix_remap_chr #freq matches between OAG and IBD, IBD chr matches dbSNP (changed chr between hapmap and 1KG)
echo "exm1625029 25" >> fix_remap_chr #listed as X in IBD/OAG when is X_PAR
echo "exm1625216 25" >> fix_remap_chr #listed as X in IBD/OAG when is X_PAR
echo "exm1667356 25" >> fix_remap_chr #listed as X in IBD/OAG when is X_PAR
echo "rs6615048 23" >> fix_remap_chr #listed as X_PAR in WAMHS when is X
echo "rs6615211 23" >> fix_remap_chr #listed as X_PAR in WAMHS when is X
echo "rs6618904 23" >> fix_remap_chr #listed as X_PAR in WAMHS when is X
echo "rs6618931 23" >> fix_remap_chr #listed as X_PAR in WAMHS when is X
echo "rs6618943 23" >> fix_remap_chr #listed as X_PAR in WAMHS when is X
echo "exm2268606 23" >> fix_remap_chr #listed as X_PAR in WAMHS when is X
echo "rs1198735 23" >> fix_remap_chr #listed as X_PAR in WAMHS when is X

#loop sections 1 and 2 to be able to turn on/off
#for x in #1
#do

  ########################################################################
  #1.0 Fix up the BEACON fam file to make it loopable
  ########################################################################
##I use # to deactivate lines with are not relavant to the current study.
  
#14/06/2016 BEACON data is a bit borked - FID, PID and MID are all NA, which makes PLINK treat them as non founders. So reset them all to IID,0 and 0 respectively. Also fix status to controls
  #echo -e "\nExtracting BEACON controls and fixing their PID, FID and MID (all NA).\n"

  #07/07/2016 these functions are all straight forward (will only change if the source files moves) so add --silent
  #awk '{print $1,$2}' ${dirJB}/${raw_data2_p} > ${dirI}/temp_${raw_data2_p}
  #plink_1.90 --threads 1 --bfile ${dirJB}/${raw_data2} --keep ${dirI}/temp_${raw_data2_p} --make-bed --out ${dirI}/temp_2016_cleaning_SDH --silent
    #Total genotyping rate in remaining samples is 0.879269. 1134514 variants and 570 people pass filters and QC.

  #check the sex
  #plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_SDH --check-sex --out ${dirI}/temp_2016_cleaning_SDH --nonfounders --silent
    #no problems found, disable

  #11/07/2016 getting het haploid warning, and PLINK recommends you fix these asap. I previously tried splitX to fix this but it doesn't (PLINK seems to have a differnt idea about what are the boundaries for the PAR region. 1673 SNPs but many are dups (het in more than one person)
  #echo ""
  #awk '{print $3}' ${dirI}/temp_2016_cleaning_SDH.hh | sort | uniq >  ${dirI}/temp_2016_cleaning_SDH.hh_snps; wc -l ${dirI}/temp_2016_cleaning_SDH.hh_snps
    #563

  #quick check. Some are chrY (17 or so), rest are chromsome X. Looked around a bit and even though GRC etc have slightly different PAR boundaries http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/ this doesn't explain  these. Checked a handful and they have the correct hg19 positions. So probabyly best to just remove them. 
  #awk 'NR==FNR {a[$1];next} $2 in a' ${dirI}/temp_2016_cleaning_SDH.hh_snps ${dirI}/temp_2016_cleaning_SDH.bim | awk '($1==23 && ($4 >= 60001 && $4 <= 2699520)) ||($1==23 && ($4 >= 154931044 && $4 <= 155260560)) || ($1==24 && ($4 >= 10001 && $4 <= 2649520)) || ($1==24 && ($4 >= 59034050 && $4 <= 59363566))' | wc -l
    #0

  #echo  ""
  #plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_SDH --freq --out ${dirI}/temp_2016_cleaning_SDH --nonfounders --silent 

  #echo ""
  #awk '$1!=0 && $5!=0 && $5<0.01 && $5>0.001' ${dirI}/temp_2016_cleaning_SDH.frq | wc -l
    #36538 after filtering .hh SNPs

  #need to fix the IDs/founder status. In fact the whole ID file is a problem  - FID is NA, family status is NA NA, fix to FID=IID and change them to controls
  #awk '{print $2,$2,"0","0",$5,"1"}' ${dirI}/temp_2016_cleaning_SDH.fam > ${dirI}/temp_2016_cleaning_SDH2.fam
  #echo -e "\nRemoving het haploid SNPs that do not fall in any accepted PAR boundary positions"
  #plink_1.90 --threads 1 --bed ${dirI}/temp_2016_cleaning_SDH.bed --fam ${dirI}/temp_2016_cleaning_SDH2.fam --bim ${dirI}/temp_2016_cleaning_SDH.bim --make-bed --out ${dirI}/temp_2016_cleaning_SDH_2 --exclude ${dirI}/temp_2016_cleaning_SDH.hh_snps --silent
    #--exclude: 1133951 variants remaining. 1133951 variants and 570 people pass filters and QC. Among remaining phenotypes, 0 are cases and 570 are controls.

  #intermin clean
  #rm -f ${dirI}/temp_2016_cleaning_SDH_omni_amb_pre ${dirI}/temp_2016_cleaning_SDH.* ${dirI}/temp_${raw_data2_p} ${dirI}/temp_2016_cleaning_SDH.hh_snps 

  ########################################################################
  #2.0 Initial cleaninf for phase 4 , non-advanced POAG
  ########################################################################
  #15/06/2016 I was going to try align/remap the samples before splitting them but I just realised I need to align them as seperate case and controls groups for each analysis set or I wont pick up borked SNPs, especially since I want to try keep in ambigious SNPs (e.g. I may not captue an A/T SNP flipped in QMEGA_610k relative to endo/qtwin 610k. The melanoma dataset is three sets merged without cleaning (provided by SM back when I started). It is the AMFS (omni1M), and QMEGA cases on omni1M, QMEGA cases and 610, and endo + qtwin controls on 610/670k. Need to seperate them out, and give them similar names, so I can then work with loops

  #11/07/2016 I now want to drop .hh SNPs at the start so run a test to generate a file


  echo -e "\nInitial cleaninf for phase 4 , non-advanced POAG by removing het haploid SNPs etc...\n"
  
  # to male bfile format:
  plink_1.90 --threads 1 --ped ${dirJ}/${raw_data}.ped --map ${dirJ}/${raw_data}.map --make-bed --out ${dirJ}/${raw_data}
  plink_1.90 --threads 1 --bfile  ${dirJ}/${raw_data} --freq --out ${dirI}/temp_${raw_data} --silent

  echo ""
  awk '{print $3}' ${dirI}/temp_${raw_data}.hh | sort | uniq > ${dirI}/temp_${raw_data}.hh_snps; wc -l ${dirI}/temp_${raw_data}.hh_snps

  #do any fall in the boundaries http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/
  echo ""
  awk 'NR==FNR {a[$1];next} $2 in a' ${dirI}/temp_${raw_data}.hh_snps ${dirJ}/${raw_data}.bim | awk '($1==23 && ($4 >= 60001 && $4 <= 2699520)) ||($1==23 && ($4 >= 154931044 && $4 <= 155260560)) || ($1==24 && ($4 >= 10001 && $4 <= 2649520)) || ($1==24 && ($4 >= 59034050 && $4 <= 59363566))' | wc -l 
    #0

  #also getting "Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands treat these as missing". This error is less serious; can remove it with --set-hh-missing if you think all the genders are correct (if the gender is incorrect and you do this, you lose data. So first test gender.
  echo ""
  plink_1.90 --threads 1 --bfile  ${dirJ}/${raw_data} --exclude ${dirI}/temp_${raw_data}.hh_snps --check-sex --out ${dirI}/temp_${raw_data} --silent
  echo -e "\nCheck sex results for ${raw_data} : showing non OK results\n"
  grep -v "OK" ${dirI}/temp_${raw_data}.sexcheck
  echo -e "\nCorrect gender based on sex check, and remove people with failed sex check, and removing het haploid SNPs and setting non zero nonmale Y calls to missing\n"
  awk '$4==0 {print $1,$2}' ${dirI}/temp_${raw_data}.sexcheck > ${dirI}/temp_${raw_data}.sexcheck_remove; wc -l ${dirI}/temp_${raw_data}.sexcheck_remove
  echo ""
  #21/07/2016 checking this over, found I wasn't actually writing out the SNPSEX, just the whole file (which would have resulted in the PED sex being used)
  awk '$5=="PROBLEM" && $4!="0" {print $1,$2,$4}' ${dirI}/temp_${raw_data}.sexcheck > ${dirI}/temp_${raw_data}.sexcheck_update; wc -l ${dirI}/temp_${raw_data}.sexcheck_update
    #6
  echo ""
  #so remove them them, and then work from that copy
  plink_1.90 --threads 1 --bfile  ${dirJ}/${raw_data} --exclude ${dirI}/temp_${raw_data}.hh_snps --make-bed --out ${dirI}/temp_${raw_data} --remove ${dirI}/temp_${raw_data}.sexcheck_remove --update-sex ${dirI}/temp_${raw_data}.sexcheck_update --silent

  #still getting "Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands treat these as missing." Running this command oesnt' give a lot of info but the missingness rate, SNP N doesnt' really change. I tried to work out which SNP(s) were causing this but filter-females/filter males then freq found no non male chr 24 SNPs with a valid genotype. #looking at them it doesn't seeem to be systematic; eahc perosn has a different single genotype, and always homozygous, so might just be borderline signals in the cluster plot s? As in rare false positives rather than a bad SNP. So just use hh. Not worth folowing further, just set to missing
    ##echo ""
    plink_1.90 --threads 1 --bfile ${dirI}/temp_${raw_data} --chr 24 --recode --out ${dirI}/temp_${raw_data}_2
    ###use set hh and try find the diff
    plink_1.90 --threads 1 --file  ${dirI}/temp_${raw_data}_2 --set-hh-missing --recode --out ${dirI}/temp_${raw_data}_3
    wc -l ${dirI}/temp_${raw_data}_2.ped
    awk 'NR==FNR {a[$0];next} !($0 in a)' ${dirI}/temp_${raw_data}_3.ped ${dirI}/temp_${raw_data}_2.ped | wc -l
    ###324 changed lines
    rm ${dirI}/temp_${raw_data}_2* ${dirI}/temp_${raw_data}_3*
  #set hh to missing, then overwrite the earlier version
  plink_1.90 --threads 1 --bfile ${dirI}/temp_${raw_data} --set-hh-missing --make-bed --out ${dirI}/temp_${raw_data}_2 --silent
  plink_1.90 --threads 1 --bfile ${dirI}/temp_${raw_data}_2 --make-bed --out ${dirI}/temp_${raw_data} --silent

  echo -e "\nHet haploid SNPs, non missing non male y calls removed\n"
  #interim clean
  rm -f ${dirI}/temp_${raw_data}_2.* ${dirI}/temp_${raw_data}.hh* ${dirI}/temp_${raw_data}.sexcheck*
  #Sample lists
    #978
  #echo ""
  #plink_1.90 --threads 1 --bfile  ${dirI}/temp_${raw_data} --keep ${dirI}/temp_2016_cleaning_AMFS --make-bed --out ${dirI}/temp_2016_cleaning_AMFS --silent
    #--keep: 978 people remaining.Total genotyping rate in remaining samples is 0.781615. (inc lots of 610k/670k SNPs) 1034058 variants and 978 people pass filters and QC.Among remaining phenotypes, 548 are cases and 430 are controls.

  #Now split these into cases and controls for aligning with 1kg
  #plink_1.90 --threads 1 --bfile  ${dirI}/temp_2016_cleaning_AMFS --filter-cases --make-bed --out ${dirI}/temp_2016_cleaning_AMFS_cases --silent
    #Post fixing/removing X/Y calls 1034058 (was 1043006) variants and 548 people pass filters and QC.
  #plink_1.90 --threads 1 --bfile  ${dirI}/temp_2016_cleaning_AMFS --filter-controls --make-bed --out ${dirI}/temp_2016_cleaning_AMFS_controls --silent
    #1034058 variants and 430 people pass filters and QC.

  #get rid of these from the starting file
  #plink_1.90 --threads 1 --bfile  ${dirI}/temp_${raw_data} --remove ${dirI}/temp_2016_cleaning_AMFS --make-bed --out ${dirI}/temp_2016_cleaning_1 --silent
    #6551 phenotype values loaded from .fam.--remove: 5573 people remaining. 1034058 variants and 5573 people pass filters and QC.Among remaining phenotypes, 1618 are cases and 3954 are controls.

  #now all the controls left are on the 610k/670k list; can use this to back determine samples
  #plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_1 --filter-controls --geno 0.1 --make-bed --out ${dirI}/temp_2016_cleaning_2 --silent
    #529966 variants removed due to missing genotype data (--geno). 504092 variants and 3954 people pass filters and QC Among remaining phenotypes, 0 are cases and 3954 are controls 

  #awk '{print $2}' ${dirI}/temp_2016_cleaning_2.bim > ${dirI}/temp_2016_cleaning_2_snps; wc -l ${dirI}/temp_2016_cleaning_2_snps
    #504092
  #now filter back on this
  #plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_1 --extract ${dirI}/temp_2016_cleaning_2_snps --make-bed --out ${dirI}/temp_2016_cleaning_3 --silent
    #1034058 variants loaded from .bim file. 5573 people (1614 males, 3959 females) loaded from .fam. --extract: 504092 variants remaining Total genotyping rate is 0.947502. 504092 variants and 5573 people pass filters and QC. Among remaining phenotypes, 1618 are cases and 3954 are controls.

  #now use this to ID the Q-MEGA cases on the omni chip (as will have high mind
  #plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_3 --mind 0.1 --make-bed --out ${dirI}/temp_2016_cleaning_4 --silent
    #Among remaining phenotypes, 925 are cases and 3954 are controls. <<nice, exaact >> .IDs written to /working/lab_stuartma/scratch/matthewL/melanoma_meta_analysis/working_directory/temp_2016_cleaning_4.irem
    
  #if I did this right the irem will only be cases
  #echo -e "\nChecking phenotype of people in .irem file, should all be cases (2):"
  #awk 'NR==FNR {a[$1,$2];next} ($1,$2) in a {print $6}' ${dirI}/temp_2016_cleaning_4.irem ${dirI}/temp_2016_cleaning_1.fam | sort | uniq -c
    #    693 2
    #      1 -9
  #lets me careful and make sure we haven't lost a few SNPs in this adhoc way, so work back from the sample lists
  #plink_1.90 --threads 1 --bfile ${dirI}/temp_${raw_data} --keep ${dirI}/temp_2016_cleaning_4.irem --make-bed --out ${dirI}/temp_2016_cleaning_QMEGA_omni_cases --silent
    #Among remaining phenotypes, 693 are cases and 0 are controls.  (1 phenotype is missing.)
  #now combine the AMFS and irem list
  #cat ${dirI}/temp_2016_cleaning_4.irem ${dirI}/temp_2016_cleaning_AMFS > ${dirI}/temp_2016_cleaning_non_610k

  #plink_1.90 --threads 1 --bfile ${dirI}/temp_${raw_data} --remove ${dirI}/temp_2016_cleaning_non_610k --make-bed --out ${dirI}/temp_2016_cleaning_QMEGA_610k --silent
    #1034058 variants and 4879 people pass filters and QC. Among remaining phenotypes, 925 are cases and 3954 are controls.
  #15/06/2016 split these into cases and controls
  #plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_QMEGA_610k --filter-cases --make-bed --out ${dirI}/temp_2016_cleaning_QMEGA_610k_cases --silent
    #1034058 variants and 925 people pass filters and QC.
  #plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_QMEGA_610k --filter-controls --make-bed --out ${dirI}/temp_2016_cleaning_QMEGA_610k_controls --silent
    #1034058 variants and 3954 people pass filters and QC.

  #interim clean
  #rm -f ${dirI}/temp_2016_cleaning_2_snps ${dirI}/temp_${raw_data}.*

#done #loop to turn on/off

########################################################################
#2.1 HEIDELBERG
########################################################################
#21/07/2016 as doesn't require additional sets (cleaning, aligning the QMEGA_omni set requires SDH) set this up seperately
for x in #1
do
  #need to fix hh SNPs, sexcheck etc. Fix sex first before fixing hh as many are missing any sex info which will bork hh checks
  echo -e "\nSeperating Heidelberg samples into individual case/controls sets, fixing missing sex and removing het haploid SNPs etc...\n"

  plink_1.90 --threads 1 --bfile  ${dirJG}/${raw_data3} --check-sex --out ${dirI}/temp_${raw_data3} --silent
    #--check-sex: 6958 Xchr and 0 Ychr variant(s) scanned, 373 problems detected.

  echo -e "\nCheck sex results for ${raw_data3} : showing `grep -c "PROBLEM" ${dirI}/temp_${raw_data3}.sexcheck` PROBLEMS. Removing failed, correcting missing\n"
    #373
  awk '$4==0 {print $1,$2}' ${dirI}/temp_${raw_data3}.sexcheck > ${dirI}/temp_${raw_data3}.sexcheck_remove; wc -l ${dirI}/temp_${raw_data3}.sexcheck_remove
  echo ""
    #3
  awk '$5=="PROBLEM" && $4!="0" {print $1,$2,$4}' ${dirI}/temp_${raw_data3}.sexcheck > ${dirI}/temp_${raw_data3}.sexcheck_update; wc -l ${dirI}/temp_${raw_data3}.sexcheck_update
    #370
  echo ""
  plink_1.90 --threads 1 --bfile  ${dirJG}/${raw_data3} --make-bed --out ${dirI}/temp_${raw_data3} --remove ${dirI}/temp_${raw_data3}.sexcheck_remove --update-sex ${dirI}/temp_${raw_data3}.sexcheck_update --silent
    #2441 phenotype values loaded from .fam. --update-sex: 370 people updated --remove: 2438 people remaining. 222482 variants and 2438 people pass filters and QC.
  #not sure when hh SNPs are generated relative to updating sex so regenerate them again to make sure (just run freq)
  plink_1.90 --threads 1 --bfile ${dirI}/temp_${raw_data3} --freq --out ${dirI}/temp_${raw_data3} --silent
    #Warning: 73970 het. haploid genotypes present
  #check if any fall in hg18 NON PAR regions 
  echo ""
  awk '{print $3}' ${dirI}/temp_${raw_data3}.hh | sort | uniq > ${dirI}/temp_${raw_data3}.hh_snps; wc -l ${dirI}/temp_${raw_data3}.hh_snps
    #698
  #do any fall in the hg18 PAR boundaries
  echo ""
  awk 'NR==FNR {a[$1];next} $2 in a' ${dirI}/temp_${raw_data3}.hh_snps ${dirJG}/${raw_data3}.bim | awk '($1==23 && ($4 >= 1 && $4 <= 2709520)) ||($1==23 && ($4 >= 154584238 && $4 <= 154913754)) || ($1==24 && ($4 >= 1 && $4 <= 2709520)) || ($1==24 && ($4 >= 57443438 && $4 <= 57772954))' | wc -l
    #0; so remove
  #as above was still getting the 'non-missing Y errors'; in melanoma.bim found these to be single calls scattered across SNPs, not worth following up, so get rid of those with set-hh-missing
  plink_1.90 --threads 1 --bfile ${dirI}/temp_${raw_data3} --exclude ${dirI}/temp_${raw_data3}.hh_snps --make-bed --out  ${dirI}/temp_${raw_data3}_2 --set-hh-missing --silent
    #--exclude: 221784 variants remaining.
  plink_1.90 --threads 1 --bfile ${dirI}/temp_${raw_data3}_2 --filter-controls  --make-bed --out ${dirI}/temp_2016_cleaning_HEIDELBERG_controls --silent
    #221784 variants and 1221 people pass filters and QC.
  plink_1.90 --threads 1 --bfile ${dirI}/temp_${raw_data3}_2 --filter-cases --make-bed --out ${dirI}/temp_2016_cleaning_HEIDELBERG_cases --silent
    #221784 variants and 1217 people pass filters and QC.

  echo -e "\nHeidelberg extracted to case and control, with hh and sex errors fixed."
  rm -f ${dirI}/temp_${raw_data3}.* ${dirI}/temp_${raw_data3}_2.*

done #loop to turn on/off

########################################################################
#2.2 WAMHS + IBD + Glaucoma
########################################################################
#21/07/2016 set up WAMHS seperately. Not imputed with BEACON this time as Puya is doing that himself.
for x in #1
do
  echo -e "\nExtracting initial IBD data for WAMHS\n"
  #this is a simplified version of the initial steps in /working/lab_stuartma/scratch/matthewL/melanoma_meta_analysis/CIDR_WAMHS/IBD_omni_express_exome/script_to_clean_IBD_GWAS_data.sh
  #fix Gabriel's replacement of all spaces in IDs with underscores, even if just spaces in ID fields
  sed 's/___/_/g' ${dirJI}/${raw_data4}.fam > ${dirI}/temp_${raw_data4}_fam; wc -l ${dirI}/temp_${raw_data4}_fam
    #960
  #there are a bunch of control samples run on each array (e.g. CON_5738_10 is on array 10). I had a way of making these unique. Not sure why, but may as well keep this
  for i in {1..10}
  do
    sed -i "s/_${i} / /g"  ${dirI}/temp_${raw_data4}_fam
  done

  echo ""
  head -5 ${dirJI}/${raw_data4}.fam
  echo ""
  head -5 ${dirI}/temp_${raw_data4}_fam

  for i in CON_216 CON_1490 CON_1652
  do
    var_A=`grep -c ${i} ${dirJI}/${IBD_phenotypes} `
      for k in `seq 1 ${var_A}`
      do
        #replace the duplicate ID with the unique ID FID first
        sed -i "0,/${i} /s/${i} /${i}_${k} /" ${dirI}/temp_${raw_data4}_fam
        #then IID
        sed -i "0,/ ${i} /s/ ${i} / ${i}_${k} /" ${dirI}/temp_${raw_data4}_fam
      done #number of occurances
  done #loop to fix control names
  echo ""
  grep 'CON_216\|CON_1490\|CON_1652' ${dirJI}/${raw_data4}.fam
  echo ""
  grep 'CON_216\|CON_1490\|CON_1652' ${dirI}/temp_${raw_data4}_fam
  echo ""
  awk '{print $1,$2,$3,$4,$5,1}' ${dirI}/temp_${raw_data4}_fam > ${dirI}/temp_${raw_data4}_as_controls.fam; wc -l ${dirI}/temp_${raw_data4}_as_controls.fam  
    #960
  #now I had a fancy process in /working/lab_stuartma/scratch/matthewL/melanoma_meta_analysis/CIDR_WAMHS/IBD_omni_express_exome/script_to_clean_IBD_GWAS_data.sh that used the many, many controls to find discordant SNPs. I don't see a reason to repeat that, just use the list there and drop the duplicate controls used to do that
  echo -e "\nRemoving duplicate control samples from the IBD data, as well as SNPs previously found to be discordant in the duplicate control samples\n"  
  plink_1.90 --threads 1 --bed ${dirJI}/${raw_data4}.bed --bim ${dirJI}/${raw_data4}.bim --fam ${dirI}/temp_${raw_data4}_as_controls.fam --remove ${dirJI}/IBD_duplicate_IDs.txt --exclude ${dirJI}/IBD_snps_discordant_between_replicates.txt --make-bed --out ${dirI}/temp_${raw_data4} --silent
    #951117 variants loaded from .bim file. --exclude: 951101 variants remaining. --remove: 930 people remaining. 951101 variants and 930 people pass filters and QC.
  echo -e "\nChecking sex for IBD data"
  #now check sex
  plink_1.90 --threads 1 --bfile  ${dirI}/temp_${raw_data4} --check-sex --out ${dirI}/temp_${raw_data4} --silent
    #--check-sex: 17013 Xchr and 0 Ychr variant(s) scanned, 1 problem detected.

  echo -e "\nCheck sex results for ${raw_data4} : showing `grep -c "PROBLEM" ${dirI}/temp_${raw_data4}.sexcheck` PROBLEMS. Removing failed, correcting missing\n"
    #1
  awk '$4==0 {print $1,$2}' ${dirI}/temp_${raw_data4}.sexcheck > ${dirI}/temp_${raw_data4}.sexcheck_remove; wc -l ${dirI}/temp_${raw_data4}.sexcheck_remove
  echo ""
    #1
  awk '$5=="PROBLEM" && $4!="0" {print $1,$2,$4}' ${dirI}/temp_${raw_data4}.sexcheck > ${dirI}/temp_${raw_data4}.sexcheck_update; wc -l ${dirI}/temp_${raw_data4}.sexcheck_update
    #0
  echo ""
  plink_1.90 --threads 1 --bfile  ${dirI}/temp_${raw_data4} --make-bed --out ${dirI}/temp_${raw_data4}_2 --remove ${dirI}/temp_${raw_data4}.sexcheck_remove --update-sex ${dirI}/temp_${raw_data4}.sexcheck_update --silent
   #930 phenotype values loaded from .fam --update-sex: 0 people updated. --remove: 929 people remaining.
   #not sure when hh SNPs are generated relative to updating sex so regenerate them again to make sure (just run freq)
   plink_1.90 --threads 1 --bfile ${dirI}/temp_${raw_data4}_2 --freq --out ${dirI}/temp_${raw_data4} --silent
     #Warning: 60547 het. haploid genotypes present
  #check if any fall in hg19 NON PAR regions 
  awk '{print $3}' ${dirI}/temp_${raw_data4}.hh | sort | uniq > ${dirI}/temp_${raw_data4}.hh_snps; wc -l ${dirI}/temp_${raw_data4}.hh_snps
    #296
  #do any fall in the hg19 PAR boundaries
  echo ""
  awk 'NR==FNR {a[$1];next} $2 in a' ${dirI}/temp_${raw_data4}.hh_snps ${dirJI}/${raw_data4}.bim | awk '($1==23 && ($4 >= 60001 && $4 <= 2699520)) ||($1==23 && ($4 >= 154931044 && $4 <= 155260560)) || ($1==24 && ($4 >= 10001 && $4 <= 2649520)) || ($1==24 && ($4 >= 59034050 && $4 <= 59363566))' | wc -l
  echo ""
    #36, this is new. Update them to chr25, exclude the rest
  awk 'NR==FNR {a[$1];next} $2 in a' ${dirI}/temp_${raw_data4}.hh_snps ${dirJI}/${raw_data4}.bim | awk '($1==23 && ($4 >= 60001 && $4 <= 2699520)) ||($1==23 && ($4 >= 154931044 && $4 <= 155260560)) || ($1==24 && ($4 >= 10001 && $4 <= 2649520)) || ($1==24 && ($4 >= 59034050 && $4 <= 59363566)) {print $2,25}' > ${dirI}/temp_${raw_data4}_PAR_snps;wc -l ${dirI}/temp_${raw_data4}_PAR_snps
    #36
  echo ""
  awk 'NR==FNR {a[$1];next} !($1 in a)' ${dirI}/temp_${raw_data4}_PAR_snps ${dirI}/temp_${raw_data4}.hh_snps > ${dirI}/temp_${raw_data4}.hh_snps_exclude; wc -l ${dirI}/temp_${raw_data4}.hh_snps_exclude
    #260
  #as above was still getting the 'non-missing Y errors'; in melanoma.bim found these to be single calls scattered across SNPs, not worth following up. Normally would run this now with set-hh-missing but not sure of order of that with update-chr so run again after
  plink_1.90 --threads 1 --bfile ${dirI}/temp_${raw_data4}_2 --update-chr ${dirI}/temp_${raw_data4}_PAR_snps --exclude ${dirI}/temp_${raw_data4}.hh_snps_exclude --make-bed --out  ${dirI}/temp_${raw_data4} --silent
    #951101 variants loaded from .bim file. --exclude: 950841 variants remaining --update-chr: 36 values updated.
  plink_1.90 --threads 1 --bfile ${dirI}/temp_${raw_data4} --set-hh-missing --make-bed --out ${dirI}/temp_${raw_data4}_2 --silent
    
  #now there are a bunch of SNPs with B alleles - which means gabriel couldn't remap them. Are they all monomorphic?
  echo ""
  plink_1.90 --threads 1 --bfile ${dirI}/temp_${raw_data4}_2 --freq --out ${dirI}/temp_${raw_data4} --silent
  awk '($3=="B" || $4=="B") && $5>0' ${dirI}/temp_${raw_data4}.frq | wc -l
    #3983 - ffs

  #have what I be the update_allele file from http://www.well.ox.ac.uk/~wrayner/strand/ABtoTOPstrand.html
  #zgrep kgp22785494 ${dirI}/HumanOmniExpress-12v1-1_B.update_alleles.txt.gHumanOmniExpressExome-8-v1-2-B.update_alleles.txt.gz 
  #SNPs aren't in there, but messing around with the T/B SNPs that have chr and bp position it looks like (a) most are in 1kg so can remap them. I just realised there may be the reverse - SNPs where the B has been remapped, but not the A - can't readily find these other than them being treated as triplicate SNPs
  awk 'NR==FNR {a[$2]=$3" "$4" "$5" "$6;next} $2 in a {print $0,a[$2]}' <(awk '($3=="B" || $4=="B") && $5>0' ${dirI}/temp_${raw_data4}.frq) <( sed 's/\t/ /g' ${dirI}/temp_${raw_data4}_2.bim | sed 's/ \+/ /g' | sed 's/^ //g' ) > ${dirI}/temp_${raw_data4}_B_allele_remap_1; wc -l ${dirI}/temp_${raw_data4}_B_allele_remap_1 
    #3983
  #now map to 1kg (this will drop those without bp info
  awk 'NR==FNR {a[$1,$4]=$0;next} ($1,$3) in a {print a[$1,$3],$0}' ${dirI}/temp_${raw_data4}_B_allele_remap_1 <(gunzip -dc ${dirI}/ALL_1000G_phase1integrated_v3_ALL_chr_legends_with_chromosome.txt.gz) > ${dirI}/temp_${raw_data4}_B_allele_remap_2; wc -l ${dirI}/temp_${raw_data4}_B_allele_remap_2
   #3797 < expected a few for for overlapping indels, hence merging 1kg onto the data not the other way around (otherwise might miss real pairing)

   cd $dirI
   #same strand, same allele - e.g. T B 0.03983 A T 0.0422. Need to filter out indels in the 1kg data. high maf ambigious, and then try work out from allele matching and freq matching what the B allele is, and may as well flip it to match 1kg while in here (as need to check for that anyway to make sense of the matches)
  awk '{
    if( length($14)==1 && length($15)==1 && $23 < '''${MAF}''' && ( ( $14=="T" && $15=="A" ) || ( $14=="A" && $15=="T" ) || ( $14=="G" && $15=="C" ) || ( $14=="C" && $15=="G" ) ) && ( ( $5==$15) || ( $6==$14) ) && ( ($9 - $19) > -('''${DIFF}''') ) && ( ($9 - $19) < '''${DIFF}''' ) )
      { print $2,$5,$6,$15,$14 > "temp_IBD_B_alleles_ambigious_aligned_same_allele"};
    if( length($14)==1 && length($15)==1 && $23 < '''${MAF}''' && ( ( $14=="T" && $15=="A" ) || ( $14=="A" && $15=="T" ) || ( $14=="G" && $15=="C" ) || ( $14=="C" && $15=="G" ) ) && ( ( $5==$14 ) || ( $6==$15) ) && ( ($9 - ( 1- $19 )) > -('''${DIFF}''') ) && ( ($9 - ( 1- $19 ) ) < '''${DIFF}''' ) )
      { print $2,$5,$6,$14,$15 > "temp_IBD_B_alleles_ambigious_aligned_other_allele"};
    if( length($14)==1 && length($15)==1 && $23 < '''${MAF}''' && ( ( $14=="T" && $15=="A" ) || ( $14=="A" && $15=="T" ) || ( $14=="G" && $15=="C" ) || ( $14=="C" && $15=="G" ) ) && ( ( $5==$14 ) || ( $6==$15) ) && ( ($9 - $19) > -('''${DIFF}''') ) && ( ($9 - $19) < '''${DIFF}''' ) )
      { print $2,$5,$6,$15,$14 > "temp_IBD_B_alleles_ambigious_flipped_same_allele"};
    if( length($14)==1 && length($15)==1 && $23 < '''${MAF}''' && ( ( $14=="T" && $15=="A" ) || ( $14=="A" && $15=="T" ) || ( $14=="G" && $15=="C" ) || ( $14=="C" && $15=="G" ) ) && ( ( $5==$15) || ( $6==$14) ) && ( ($9 - ( 1- $19 )) > -('''${DIFF}''') ) && ( ($9 - ( 1- $19 ) ) < '''${DIFF}''' ) )
      { print $2,$5,$6,$14,$15 > "temp_IBD_B_alleles_ambigious_flipped_other_allele"};
    if( length($14)==1 && length($15)==1 && !( $14=="T" && $15=="A" ) && !( $14=="A" && $15=="T" ) && !( $14=="G" && $15=="C" ) && !( $14=="C" && $15=="G" ) && ( $5==$15 || $6==$14) && ( ($9 - $19) > -('''${DIFF}''') ) && ( ($9 - $19) < '''${DIFF}''' ) )
      { print $2,$5,$6,$15,$14 > "temp_IBD_B_alleles_aligned_same_allele"};
    if( length($14)==1 && length($15)==1 && !( $14=="T" && $15=="A" ) && !( $14=="A" && $15=="T" ) && !( $14=="G" && $15=="C" ) && !( $14=="C" && $15=="G" ) && ( $5==$14 || $6==$15) && ( ($9 - ( 1- $19 )) > -('''${DIFF}''') ) && ( ($9 - ( 1- $19 ) ) < '''${DIFF}''' ) )
      { print $2,$5,$6,$14,$15 > "temp_IBD_B_alleles_aligned_other_allele"};
    if( length($14)==1 && length($15)==1 && !( $14=="T" && $15=="A" ) && !( $14=="A" && $15=="T" ) && !( $14=="G" && $15=="C" ) && !( $14=="C" && $15=="G" ) && ( ($9 - $19) > -('''${DIFF}''') ) && ( ($9 - $19) < '''${DIFF}''' ) && ( ( $5=="T" && $15=="A" && $14=="G" ) || ($6=="C" && $15=="A" && $14=="G" ) ) ) 
      { print $2,$5,$6,"T","C" > "temp_IBD_B_alleles_flipped_same_allele_AG"};
    if( length($14)==1 && length($15)==1 && !( $14=="T" && $15=="A" ) && !( $14=="A" && $15=="T" ) && !( $14=="G" && $15=="C" ) && !( $14=="C" && $15=="G" ) && ( ($9 - $19) > -('''${DIFF}''') ) && ( ($9 - $19) < '''${DIFF}''' ) && ( ( $5=="T" && $15=="A" && $14=="C" ) || ($6=="G" && $15=="A" && $14=="C" ) ) )
      { print $2,$5,$6,"T","G" > "temp_IBD_B_alleles_flipped_same_allele_AC"};
    if( length($14)==1 && length($15)==1 && !( $14=="T" && $15=="A" ) && !( $14=="A" && $15=="T" ) && !( $14=="G" && $15=="C" ) && !( $14=="C" && $15=="G" ) && ( ($9 - $19) > -('''${DIFF}''') ) && ( ($9 - $19) < '''${DIFF}''' ) && ( ( $5=="A" && $15=="T" && $14=="G" ) || ( $6=="C" && $15=="T" && $14=="G" ) ) )
      { print $2,$5,$6,"A","C" > "temp_IBD_B_alleles_flipped_same_allele_TG"};
    if( length($14)==1 && length($15)==1 && !( $14=="T" && $15=="A" ) && !( $14=="A" && $15=="T" ) && !( $14=="G" && $15=="C" ) && !( $14=="C" && $15=="G" ) && ( ($9 - $19) > -('''${DIFF}''') ) && ( ($9 - $19) < '''${DIFF}''' ) && ( ( $5=="A" && $15=="T" && $14=="C" ) || ( $6=="G" && $15=="T" && $14=="C" ) ) )
      { print $2,$5,$6,"A","G" > "temp_IBD_B_alleles_flipped_same_allele_TC"};
    if( length($14)==1 && length($15)==1 && !( $14=="T" && $15=="A" ) && !( $14=="A" && $15=="T" ) && !( $14=="G" && $15=="C" ) && !( $14=="C" && $15=="G" ) && ( ($9 - $19) > -('''${DIFF}''') ) && ( ($9 - $19) < '''${DIFF}''' ) && ( ( $5=="C" && $15=="G" && $14=="A" ) || ( $6=="T" && $15=="G" && $14=="A" ) ) )
      { print $2,$5,$6,"C","T" > "temp_IBD_B_alleles_flipped_same_allele_GA"};
    if( length($14)==1 && length($15)==1 && !( $14=="T" && $15=="A" ) && !( $14=="A" && $15=="T" ) && !( $14=="G" && $15=="C" ) && !( $14=="C" && $15=="G" ) && ( ($9 - $19) > -('''${DIFF}''') ) && ( ($9 - $19) < '''${DIFF}''' ) && ( ( $5=="C" && $15=="G" && $14=="T" ) || ( $6=="A" && $15=="G" && $14=="T" ) ) )
      { print $2,$5,$6,"C","A" > "temp_IBD_B_alleles_flipped_same_allele_GT"};
    if( length($14)==1 && length($15)==1 && !( $14=="T" && $15=="A" ) && !( $14=="A" && $15=="T" ) && !( $14=="G" && $15=="C" ) && !( $14=="C" && $15=="G" ) && ( ($9 - $19) > -('''${DIFF}''') ) && ( ($9 - $19) < '''${DIFF}''' ) && ( ( $5=="G" && $15=="C" && $14=="A" ) || ( $6=="T" && $15=="C" && $14=="A" ) ) )
      { print $2,$5,$6,"G","T" > "temp_IBD_B_alleles_flipped_same_allele_CA"};
   if( length($14)==1 && length($15)==1 && !( $14=="T" && $15=="A" ) && !( $14=="A" && $15=="T" ) && !( $14=="G" && $15=="C" ) && !( $14=="C" && $15=="G" ) && ( ($9 - $19) > -('''${DIFF}''') ) && ( ($9 - $19) < '''${DIFF}''' ) && ( ( $5=="G" && $15=="C" && $14=="T" ) || ( $6=="A" && $15=="C" && $14=="T" ) ) )
      { print $2,$5,$6,"G","A" > "temp_IBD_B_alleles_flipped_same_allele_CT"};
    if( length($14)==1 && length($15)==1 && !( $14=="T" && $15=="A" ) && !( $14=="A" && $15=="T" ) && !( $14=="G" && $15=="C" ) && !( $14=="C" && $15=="G" ) && ( ($9 - ( 1 - $19 ) ) > -('''${DIFF}''') ) && ( ($9 - ( 1 - $19 ) ) < '''${DIFF}''' ) && ( ( $5=="T" && $14=="A" && $15=="G" ) || ( $6=="C" && $14=="A" && $15=="G" ) ) )
      { print $2,$5,$6,"T","C" > "temp_IBD_B_alleles_flipped_other_allele_GA"};
    if( length($14)==1 && length($15)==1 && !( $14=="T" && $15=="A" ) && !( $14=="A" && $15=="T" ) && !( $14=="G" && $15=="C" ) && !( $14=="C" && $15=="G" ) && ( ($9 - ( 1 - $19 ) ) > -('''${DIFF}''') ) && ( ($9 - ( 1 - $19 ) ) < '''${DIFF}''' ) && ( ( $5=="T" && $14=="A" && $15=="C" ) || ( $6=="G" && $14=="A" && $15=="C" ) ) )
      { print $2,$5,$6,"T","G" > "temp_IBD_B_alleles_flipped_other_allele_CA"};
    if( length($14)==1 && length($15)==1 && !( $14=="T" && $15=="A" ) && !( $14=="A" && $15=="T" ) && !( $14=="G" && $15=="C" ) && !( $14=="C" && $15=="G" ) && ( ($9 - ( 1 - $19 ) ) > -('''${DIFF}''') ) && ( ($9 - ( 1 - $19 ) ) < '''${DIFF}''' ) && ( ( $5=="A" && $14=="T" && $15=="G" ) || ( $6=="C" && $14=="T" && $15=="G" ) ) )
      { print $2,$5,$6,"A","C" > "temp_IBD_B_alleles_flipped_other_allele_GT"};
    if( length($14)==1 && length($15)==1 && !( $14=="T" && $15=="A" ) && !( $14=="A" && $15=="T" ) && !( $14=="G" && $15=="C" ) && !( $14=="C" && $15=="G" ) && ( ($9 - ( 1 - $19 ) ) > -('''${DIFF}''') ) && ( ($9 - ( 1 - $19 ) ) < '''${DIFF}''' ) && ( ( $5=="A" && $14=="T" && $15=="C" ) || ( $6=="G" && $14=="T" && $15=="C" ) ) )
      { print $2,$5,$6,"A","G" > "temp_IBD_B_alleles_flipped_other_allele_CT"};
    if( length($14)==1 && length($15)==1 && !( $14=="T" && $15=="A" ) && !( $14=="A" && $15=="T" ) && !( $14=="G" && $15=="C" ) && !( $14=="C" && $15=="G" ) && ( ($9 - ( 1 - $19 ) ) > -('''${DIFF}''') ) && ( ($9 - ( 1 - $19 ) ) < '''${DIFF}''' ) && ( ( $5=="C" && $14=="G" && $15=="A" ) || ( $6=="T" && $14=="G" && $15=="A" ) ) )
      { print $2,$5,$6,"C","T" > "temp_IBD_B_alleles_flipped_other_allele_AG"};
    if( length($14)==1 && length($15)==1 && !( $14=="T" && $15=="A" ) && !( $14=="A" && $15=="T" ) && !( $14=="G" && $15=="C" ) && !( $14=="C" && $15=="G" ) && ( ($9 - ( 1 - $19 ) ) > -('''${DIFF}''') ) && ( ($9 - ( 1 - $19 ) ) < '''${DIFF}''' ) && ( ( $5=="C" && $14=="G" && $15=="T" ) || ( $6=="A" && $14=="G" && $15=="T" ) ) )
      { print $2,$5,$6,"C","A" > "temp_IBD_B_alleles_flipped_other_allele_TG"};
    if( length($14)==1 && length($15)==1 && !( $14=="T" && $15=="A" ) && !( $14=="A" && $15=="T" ) && !( $14=="G" && $15=="C" ) && !( $14=="C" && $15=="G" ) && ( ($9 - ( 1 - $19 ) ) > -('''${DIFF}''') ) && ( ($9 - ( 1 - $19 ) ) < '''${DIFF}''' ) && ( ( $5=="G" && $14=="C" && $15=="A" ) || ( $6=="T" && $14=="C" && $15=="A" ) ) )
      { print $2,$5,$6,"G","T" > "temp_IBD_B_alleles_flipped_other_allele_AC"};
    if( length($14)==1 && length($15)==1 && !( $14=="T" && $15=="A" ) && !( $14=="A" && $15=="T" ) && !( $14=="G" && $15=="C" ) && !( $14=="C" && $15=="G" ) && ( ($9 - ( 1 - $19 ) ) > -('''${DIFF}''') ) && ( ($9 - ( 1 - $19 ) ) < '''${DIFF}''' ) && ( ( $5=="G" && $14=="C" && $15=="T" ) || ( $6=="A" && $14=="C" && $15=="T" ) ) )
      { print $2,$5,$6,"G","A" > "temp_IBD_B_alleles_flipped_other_allele_TC"};
}' ${dirI}/temp_${raw_data4}_B_allele_remap_2
  echo ""
  cat temp_IBD_B_alleles_ambigious_aligned_same_allele temp_IBD_B_alleles_ambigious_aligned_other_allele temp_IBD_B_alleles_ambigious_flipped_same_allele temp_IBD_B_alleles_ambigious_flipped_other_allele temp_IBD_B_alleles_aligned_same_allele temp_IBD_B_alleles_aligned_other_allele temp_IBD_B_alleles_flipped_same_allele_AG temp_IBD_B_alleles_flipped_same_allele_AC temp_IBD_B_alleles_flipped_same_allele_TG temp_IBD_B_alleles_flipped_same_allele_TC temp_IBD_B_alleles_flipped_same_allele_GA temp_IBD_B_alleles_flipped_same_allele_GT temp_IBD_B_alleles_flipped_same_allele_CA temp_IBD_B_alleles_flipped_same_allele_CT temp_IBD_B_alleles_flipped_other_allele_GA temp_IBD_B_alleles_flipped_other_allele_CA temp_IBD_B_alleles_flipped_other_allele_GT temp_IBD_B_alleles_flipped_other_allele_CT temp_IBD_B_alleles_flipped_other_allele_AG temp_IBD_B_alleles_flipped_other_allele_TG temp_IBD_B_alleles_flipped_other_allele_AC temp_IBD_B_alleles_flipped_other_allele_TC 2> /dev/null | wc -l 
   #3419 
  echo ""
  #check for an error - same SNP in more than one file
  cat temp_IBD_B_alleles_ambigious_aligned_same_allele temp_IBD_B_alleles_ambigious_aligned_other_allele temp_IBD_B_alleles_ambigious_flipped_same_allele temp_IBD_B_alleles_ambigious_flipped_other_allele temp_IBD_B_alleles_aligned_same_allele temp_IBD_B_alleles_aligned_other_allele temp_IBD_B_alleles_flipped_same_allele_AG temp_IBD_B_alleles_flipped_same_allele_AC temp_IBD_B_alleles_flipped_same_allele_TG temp_IBD_B_alleles_flipped_same_allele_TC temp_IBD_B_alleles_flipped_same_allele_GA temp_IBD_B_alleles_flipped_same_allele_GT temp_IBD_B_alleles_flipped_same_allele_CA temp_IBD_B_alleles_flipped_same_allele_CT temp_IBD_B_alleles_flipped_other_allele_GA temp_IBD_B_alleles_flipped_other_allele_CA temp_IBD_B_alleles_flipped_other_allele_GT temp_IBD_B_alleles_flipped_other_allele_CT temp_IBD_B_alleles_flipped_other_allele_AG temp_IBD_B_alleles_flipped_other_allele_TG temp_IBD_B_alleles_flipped_other_allele_AC temp_IBD_B_alleles_flipped_other_allele_TC 2> /dev/null | cut -d" " -f1 | sort | uniq -d | wc -l 
  #0
  echo ""
  #those not in there
  awk 'NR==FNR {a[$1];next} !($2 in a )' <( cat temp_IBD_B_alleles_ambigious_aligned_same_allele temp_IBD_B_alleles_ambigious_aligned_other_allele temp_IBD_B_alleles_ambigious_flipped_same_allele temp_IBD_B_alleles_ambigious_flipped_other_allele temp_IBD_B_alleles_aligned_same_allele temp_IBD_B_alleles_aligned_other_allele temp_IBD_B_alleles_flipped_same_allele_AG temp_IBD_B_alleles_flipped_same_allele_AC temp_IBD_B_alleles_flipped_same_allele_TG temp_IBD_B_alleles_flipped_same_allele_TC temp_IBD_B_alleles_flipped_same_allele_GA temp_IBD_B_alleles_flipped_same_allele_GT temp_IBD_B_alleles_flipped_same_allele_CA temp_IBD_B_alleles_flipped_same_allele_CT temp_IBD_B_alleles_flipped_other_allele_GA temp_IBD_B_alleles_flipped_other_allele_CA temp_IBD_B_alleles_flipped_other_allele_GT temp_IBD_B_alleles_flipped_other_allele_CT temp_IBD_B_alleles_flipped_other_allele_AG temp_IBD_B_alleles_flipped_other_allele_TG temp_IBD_B_alleles_flipped_other_allele_AC temp_IBD_B_alleles_flipped_other_allele_TC 2> /dev/null ) ${dirI}/temp_${raw_data4}_B_allele_remap_2 | wc -l
   #376
  echo ""
  #the ones that don't match look like their make sense e.g. A B and G C. It might be that the A B was simply never remapped (as in it not "A" and an unmapped B) but diminishing returns etc, and I am then assuming I can match the two with nothing in commmon; what I am doing is potentially risk as it is (assuming that I can match a SNP with one allele). The following command filters out the ones not matches that shouldn't be matched, e.g. A B and G C. The remaining 13 that aren't matched still fail valid rules - are either pairing with an indel, or have a DIFF > 0.1. So function works. The decent N of A B and G C etc makes me bet that the A and B have both been unmatched but I can't asusme that...
   awk 'NR==FNR {a[$1];next} !($2 in a )' <( cat temp_IBD_B_alleles_ambigious_aligned_same_allele temp_IBD_B_alleles_ambigious_aligned_other_allele temp_IBD_B_alleles_ambigious_flipped_same_allele temp_IBD_B_alleles_ambigious_flipped_other_allele temp_IBD_B_alleles_aligned_same_allele temp_IBD_B_alleles_aligned_other_allele temp_IBD_B_alleles_flipped_same_allele_AG temp_IBD_B_alleles_flipped_same_allele_AC temp_IBD_B_alleles_flipped_same_allele_TG temp_IBD_B_alleles_flipped_same_allele_TC temp_IBD_B_alleles_flipped_same_allele_GA temp_IBD_B_alleles_flipped_same_allele_GT temp_IBD_B_alleles_flipped_same_allele_CA temp_IBD_B_alleles_flipped_same_allele_CT temp_IBD_B_alleles_flipped_other_allele_GA temp_IBD_B_alleles_flipped_other_allele_CA temp_IBD_B_alleles_flipped_other_allele_GT temp_IBD_B_alleles_flipped_other_allele_CT temp_IBD_B_alleles_flipped_other_allele_AG temp_IBD_B_alleles_flipped_other_allele_TG temp_IBD_B_alleles_flipped_other_allele_AC temp_IBD_B_alleles_flipped_other_allele_TC 2> /dev/null ) ${dirI}/temp_${raw_data4}_B_allele_remap_2 | awk '!($5=="B" && $6=="T" && $14=="G" && $15=="C") && !( $5=="B" && $6=="T" && $14=="C" && $15=="G"  ) && !($5=="B" && $6=="A" && $14=="G" && $15=="C") && !( $5=="B" && $6=="A" && $14=="C" && $15=="G"  ) && !($5=="T" && $6=="B" && $14=="G" && $15=="C") && !( $5=="T" && $6=="B" && $14=="C" && $15=="G"  ) && !($5=="A" && $6=="B" && $14=="G" && $15=="C") && !( $5=="A" && $6=="B" && $14=="C" && $15=="G"  ) && $23 < '''${MAF}'''  ' | head
  
  #so allele update list
  cat temp_IBD_B_alleles_ambigious_aligned_same_allele temp_IBD_B_alleles_ambigious_aligned_other_allele temp_IBD_B_alleles_ambigious_flipped_same_allele temp_IBD_B_alleles_ambigious_flipped_other_allele temp_IBD_B_alleles_aligned_same_allele temp_IBD_B_alleles_aligned_other_allele temp_IBD_B_alleles_flipped_same_allele_AG temp_IBD_B_alleles_flipped_same_allele_AC temp_IBD_B_alleles_flipped_same_allele_TG temp_IBD_B_alleles_flipped_same_allele_TC temp_IBD_B_alleles_flipped_same_allele_GA temp_IBD_B_alleles_flipped_same_allele_GT temp_IBD_B_alleles_flipped_same_allele_CA temp_IBD_B_alleles_flipped_same_allele_CT temp_IBD_B_alleles_flipped_other_allele_GA temp_IBD_B_alleles_flipped_other_allele_CA temp_IBD_B_alleles_flipped_other_allele_GT temp_IBD_B_alleles_flipped_other_allele_CT temp_IBD_B_alleles_flipped_other_allele_AG temp_IBD_B_alleles_flipped_other_allele_TG temp_IBD_B_alleles_flipped_other_allele_AC temp_IBD_B_alleles_flipped_other_allele_TC 2> /dev/null > ${dirI}/temp_${raw_data4}_B_allele_remap_file

  #and the flip list (do after)
  cat temp_IBD_B_alleles_*flipped* | cut -d" " -f1 > ${dirI}/temp_${raw_data4}_B_allele_remap_flip
  plink_1.90 --threads 1 --bfile ${dirI}/temp_${raw_data4}_2 --update-alleles ${dirI}/temp_${raw_data4}_B_allele_remap_file --make-bed --out ${dirI}/temp_${raw_data4} --silent
     #--update-alleles: 3419 variants updated. 
  #now get rid of the B alleles you can't align
  echo ""
  awk '$5=="B" || $6=="B" {print $2}' ${dirI}/temp_${raw_data4}.bim > ${dirI}/temp_${raw_data4}_exclude; wc -l ${dirI}/temp_${raw_data4}_exclude
    #2499

  #04/08/2016 there are a number of overlapping SNPs on the IBD set that are flipped relative to each other. Can't make this change for all sets as they are not necessarily fliiped in the WAMHS exm data
    echo "exm-rs6456951" > ${dirI}/temp_${raw_data4}_B_allele_remap_flip #flipped in IBD relative to rs6456951 which is on the same array
    echo "exm-rs9261434" >> ${dirI}/temp_${raw_data4}_B_allele_remap_flip #flipped in IBD relative to rs9261434 which is on the same array
    echo "exm-rs3094116" >> ${dirI}/temp_${raw_data4}_B_allele_remap_flip #flipped in IBD relative to rs3094116 which is on the same array
    echo "exm-rs12138950" >> ${dirI}/temp_${raw_data4}_B_allele_remap_flip #flipped in IBD relative to rs12138950 which is on the same array
    echo "exm-rs349475" >> ${dirI}/temp_${raw_data4}_B_allele_remap_flip #flipped in IBD relative to rs349475 which is on the same array
    echo "exm-rs1150739" >> ${dirI}/temp_${raw_data4}_B_allele_remap_flip #flipped in IBD relative to rs115748975 (alias of rs1150739) which is on the same array
    echo "exm-rs3868082" >> ${dirI}/temp_${raw_data4}_B_allele_remap_flip #flipped in IBD relative to rs114900872 (alias of rs3868082)  which is on the same array
    echo "exm-rs7819412" >> ${dirI}/temp_${raw_data4}_B_allele_remap_flip #flipped in IBD relative to rs7819412 which is on the same array

  plink_1.90 --threads 1 --bfile ${dirI}/temp_${raw_data4} --flip ${dirI}/temp_${raw_data4}_B_allele_remap_flip --make-bed --out ${dirI}/temp_${raw_data4}_2 --exclude ${dirI}/temp_${raw_data4}_exclude --silent

  #write out for next step
  plink_1.90 --threads 1 --bfile ${dirI}/temp_${raw_data4}_2  --make-bed --out ${dirI}/temp_2016_cleaning_IBD_controls --silent
    #948342 variants and 929 people pass filters and QC.
  echo -e "\nFinished extracting IBD data and fixing up A/B alleles"

  rm -f temp_IBD_B_alleles_ambigious_aligned_same_allele temp_IBD_B_alleles_ambigious_aligned_other_allele temp_IBD_B_alleles_ambigious_flipped_same_allele temp_IBD_B_alleles_ambigious_flipped_other_allele temp_IBD_B_alleles_aligned_same_allele temp_IBD_B_alleles_aligned_other_allele temp_IBD_B_alleles_flipped_same_allele_AG temp_IBD_B_alleles_flipped_same_allele_AC temp_IBD_B_alleles_flipped_same_allele_TG temp_IBD_B_alleles_flipped_same_allele_TC temp_IBD_B_alleles_flipped_same_allele_GA temp_IBD_B_alleles_flipped_same_allele_GT temp_IBD_B_alleles_flipped_same_allele_CA temp_IBD_B_alleles_flipped_same_allele_CT temp_IBD_B_alleles_flipped_other_allele_GA temp_IBD_B_alleles_flipped_other_allele_CA temp_IBD_B_alleles_flipped_other_allele_GT temp_IBD_B_alleles_flipped_other_allele_CT temp_IBD_B_alleles_flipped_other_allele_AG temp_IBD_B_alleles_flipped_other_allele_TG temp_IBD_B_alleles_flipped_other_allele_AC temp_IBD_B_alleles_flipped_other_allele_TC ${dirI}/temp_${raw_data4}*

  echo -e "\n++++++Extracting initial Gluacoma data for WAMHS+++++\n"
  #as with IBD I used some internal duplicates (in this case a MZ twinset) to flag discordant SNPs /working/lab_stuartma/scratch/matthewL/melanoma_meta_analysis/CIDR_WAMHS/GLAUCOMA/script_to_clean_glaucoma_data.sh. Script shows nothing else unusual to deal with
  wc -l /working/lab_stuartma/scratch/matthewL/melanoma_meta_analysis/CIDR_WAMHS/GLAUCOMA/Glaucoma_SNPs_discordant_in_MZ_twin_samples.txt
    #114
  #are hg18
  for i in ${OAG1} ${OAG2}
  do
    echo -e "\nSetting ${i} samples as controls\n"
    #set them as controls (to make sure all are used for checks. Note that 35 or so of OAG1 and all of OAG2 are missing phenotypes; my notes from past cleaning don't expand on why this is...
    awk '{print $1,$2,1}' ${dirOAG}/${i}.fam > ${dirI}/temp_${i}.pheno
    plink_1.90 --threads 1 --bfile ${dirOAG}/${i} --check-sex --pheno ${dirI}/temp_${i}.pheno --out ${dirI}/temp_${i} --silent
      #OAG1: 964325 variants loaded. 651 people (0 males, 0 females, 651 ambiguous) 651 phenotypes present after --pheno. NO sex data, so 651 problems
      #OAG2: 733202 variants loaded. 647 people (0 males, 0 females, 647 ambiguous). 647 phenotype values present after --pheno.
    echo -e "\nCheck sex results for ${i} : showing `grep -c "PROBLEM" ${dirI}/temp_${i}.sexcheck` PROBLEMS. Removing failed, correcting missing\n"
      #OAG1 651, OAG2 647
    awk '$4==0 {print $1,$2}' ${dirI}/temp_${i}.sexcheck > ${dirI}/temp_${i}.sexcheck_remove; wc -l ${dirI}/temp_${i}.sexcheck_remove
      #OAG1 0, OAG2 2
    echo ""
    awk '$5=="PROBLEM" && $4!="0" {print $1,$2,$4}' ${dirI}/temp_${i}.sexcheck > ${dirI}/temp_${i}.sexcheck_update; wc -l ${dirI}/temp_${i}.sexcheck_update
      #OAG1 651, OAG2 645
    echo ""
    plink_1.90 --threads 1 --bfile  ${dirOAG}/${i} --pheno ${dirI}/temp_${i}.pheno --make-bed --out ${dirI}/temp_${i} --remove ${dirI}/temp_${i}.sexcheck_remove --update-sex ${dirI}/temp_${i}.sexcheck_update --silent
     #OAG1: 964325 variants 651 people loaded, --update-sex: 651 people updated --remove: 651 people remaining. 964325 variants and 651 people pass filters and QC.
     #OAG2: 733202 variants 647 people loaded, --update-sex: 645 people updated. --remove: 645 people remaining. 733202 variants and 645 people pass filters and QC.
    if `echo ${i} | grep -q "OAGphase1"`
    then
      echo -e "\nOAG1 phase detected, changing "RS"IDs to rsIDs\n"
      head -3 ${dirI}/temp_${i}.bim
      echo ""
      awk '{print $2,$2}' ${dirI}/temp_${i}.bim | sed 's/ RS/ rs/g' > ${dirI}/temp_${i}_update; wc -l ${dirI}/temp_${i}_update
      echo ""
      plink_1.90 --threads 1 --bfile ${dirI}/temp_${i} --update-name ${dirI}/temp_${i}_update --make-bed --out ${dirI}/temp_${i} --silent
        #--update-name: 964325 values updated.
      echo ""
      head -3 ${dirI}/temp_${i}.bim
     fi
   echo ""
   #not sure when hh SNPs are generated relative to updating sex so regenerate them again to make sure (just run freq)
   plink_1.90 --threads 1 --bfile ${dirI}/temp_${i} --freq --out ${dirI}/temp_${i} --silent
       #OAG1: Warning: 53902 het. haploid genotypes present; OAG2 20965 het. haploid;
    #check if any fall in hg18 NON PAR regions 
    echo ""
    awk '{print $3}' ${dirI}/temp_${i}.hh | sort | uniq > ${dirI}/temp_${i}.hh_snps; wc -l ${dirI}/temp_${i}.hh_snps
      #OAG1: 2873; OAG2: 727
    #do any fall in the hg18 PAR boundaries (checked, both are hg18)
    echo ""
    awk 'NR==FNR {a[$1];next} $2 in a' ${dirI}/temp_${i}.hh_snps ${dirI}/temp_${i}.bim | awk '($1==23 && ($4 >= 1 && $4 <= 2709520)) ||($1==23 && ($4 >= 154584238 && $4 <= 154913754)) || ($1==24 && ($4 >= 1 && $4 <= 2709520)) || ($1==24 && ($4 >= 57443438 && $4 <= 57772954))' | wc -l
      #OAG1 1; 24 rs9286410 0 1308324 G	A. Can ignore as Y chr, just remove it. Actually OAG2 has a bunch, may as well do it properly. OAG2 36 
    echo ""
    awk 'NR==FNR {a[$1];next} $2 in a' ${dirI}/temp_${i}.hh_snps ${dirI}/temp_${i}.bim | awk '($1==23 && ($4 >= 1 && $4 <= 2709520)) ||($1==23 && ($4 >= 154584238 && $4 <= 154913754)) || ($1==24 && ($4 >= 1 && $4 <= 2709520)) || ($1==24 && ($4 >= 57443438 && $4 <= 57772954)) {print $2,25}' > ${dirI}/temp_${i}_PAR_update; wc -l ${dirI}/temp_${i}_PAR_update
      #OAG1 1, OAG2 36
    awk 'NR==FNR {a[$1];next} !($1 in a)' ${dirI}/temp_${i}_PAR_update ${dirI}/temp_${i}.hh_snps  > ${dirI}/temp_${i}.hh_snps_exclude; wc -l ${dirI}/temp_${i}.hh_snps_exclude
      #OAG2 691, OAG1 2873
    #as above was still getting the 'non-missing Y errors'; in melanoma.bim found these to be single calls scattered across SNPs, not worth following up, so get rid of those with set-hh-missing; as now updateding chr this has to be done in the second step
    plink_1.90 --threads 1 --bfile ${dirI}/temp_${i} --exclude ${dirI}/temp_${i}.hh_snps_exclude --make-bed --out  ${dirI}/temp_${i}_2 --update-chr ${dirI}/temp_${i}_PAR_update --silent
      #OAG1 --exclude: 961452 variants remaining. --update-chr: 1 value updated
      #OAG2 --exclude: 7332511 variants remaining --update-chr: 36 values updated 
    plink_1.90 --threads 1 --bfile ${dirI}/temp_${i}_2 --set-hh-missing --make-bed --out ${dirI}/temp_2016_cleaning_${i}_controls --silent
      #OAG1 961452 variants and 651 people pass filters and QC. Among remaining phenotypes, 0 are cases and 651 are controls.
      #OAG2 732511 variants and 645 people pass filters and QC. Among remaining phenotypes, 0 are cases and 645 are controls.
    rm -f ${dirI}/temp_${i}*
    echo -e "\nFinished extracting ${i}"
  done #basic handling of the OAG data

  #WAMHS - /working/lab_stuartma/scratch/matthewL/melanoma_meta_analysis/CIDR_WAMHS/working_directory/script_to_clean_WAMHS_GWAS_data.sh for how the initial CIDR cleaned files wer e converted from A/B format to ACGT. No reason to repeat that here
  echo -e "\nExtracting the WAMHS GWAS data...\n"

  plink_1.90 --threads 1 --bfile  ${dirJW}/${raw_data5} --check-sex --out ${dirI}/temp_${raw_data5} --silent
    #--check-sex: 16573 Xchr and 0 Ychr variant(s) scanned, 1 problem detected.
  echo -e "\nCheck sex results for ${raw_data5} : showing `grep -c "PROBLEM" ${dirI}/temp_${raw_data5}.sexcheck` PROBLEMS. Removing failed, correcting missing\n"
    #1
  awk '$4==0 {print $1,$2}' ${dirI}/temp_${raw_data5}.sexcheck > ${dirI}/temp_${raw_data5}.sexcheck_remove; wc -l ${dirI}/temp_${raw_data5}.sexcheck_remove
    #1
  echo ""
  awk '$5=="PROBLEM" && $4!="0" {print $1,$2,$4}' ${dirI}/temp_${raw_data5}.sexcheck > ${dirI}/temp_${raw_data5}.sexcheck_update; wc -l ${dirI}/temp_${raw_data5}.sexcheck_update
    #0
  echo ""
  plink_1.90 --threads 1 --bfile  ${dirJW}/${raw_data5} --make-bed --out ${dirI}/temp_${raw_data5} --remove ${dirI}/temp_${raw_data5}.sexcheck_remove --update-sex ${dirI}/temp_${raw_data5}.sexcheck_update --silent
   #1274 phenotype values loaded from .fam. --update-sex: 0 people updated --remove: 1273 people remaining. 951117 variants and 1273 people pass filters and QC.
   #not sure when hh SNPs are generated relative to updating sex so regenerate them again to make sure (just run freq)
   plink_1.90 --threads 1 --bfile ${dirI}/temp_${raw_data5} --freq --out ${dirI}/temp_${raw_data5} --silent
     #Warning: 380 het. haploid genotypes present
  #check if any fall in hg19 NON PAR regions
  echo ""
  awk '{print $3}' ${dirI}/temp_${raw_data5}.hh | sort | uniq > ${dirI}/temp_${raw_data5}.hh_snps; wc -l ${dirI}/temp_${raw_data5}.hh_snps
    #325
  #do any fall in the hg19 PAR boundaries
  echo ""
  awk 'NR==FNR {a[$1];next} $2 in a' ${dirI}/temp_${raw_data5}.hh_snps ${dirJW}/${raw_data5}.bim | awk '($1==23 && ($4 >= 60001 && $4 <= 2699520)) ||($1==23 && ($4 >= 154931044 && $4 <= 155260560)) || ($1==24 && ($4 >= 10001 && $4 <= 2649520)) || ($1==24 && ($4 >= 59034050 && $4 <= 59363566))' | wc -l
  echo ""
    #0; so remove
  #as above was still getting the 'non-missing Y errors'; in melanoma.bim found these to be single calls scattered across SNPs, not worth following up, so get rid of those with set-hh-missing
  plink_1.90 --threads 1 --bfile ${dirI}/temp_${raw_data5} --exclude ${dirI}/temp_${raw_data5}.hh_snps --make-bed --out  ${dirI}/temp_2016_cleaning_WAMHS_cases --set-hh-missing --silent
    #951117 variants loaded from .bim file --exclude: 950792 variants remaining 950792 variants and 1273 people pass filters and QC

  echo -e "\nWAMHS extracted to case and control, with hh and sex errors fixed."
  rm -f ${dirI}/temp_${raw_data5}*

done #loop to turn on off

########################################################################
#2.3 MIA
########################################################################
 #See /mnt/lustre/home/matthewL/Melanoma/MIA/SCRIPTS/script_to_clean_MIA_KK_data_simplified_fork.sh where I implimented a fix for general oncoarray problems based on OA consortium, and how i verified the UK control sets will be suitable. As of 08/09/2016 don't have the actual contorl genotypes in hand (extracted genotypes from existing control imputed data to test PCA matching but imputation done with --pgs-miss so missing genotypes filled in.

for x in 1
do

  echo -e "\nExtracting MIA case data"
  plink_1.90 --threads 1 --bfile ${dirMIA}/${raw_data6} --check-sex --out ${dirI}/temp_${raw_data6}_1 --silent
    #no problems detected
  echo -e "\nCheck sex results for ${raw_data6} : showing `grep -c "PROBLEM" ${dirI}/temp_${raw_data6}_1.sexcheck` PROBLEMS. Removing failed, correcting missing\n"
    #0
  #not implimenting anything further here re fixing sex check
  awk '{print $3}' ${dirI}/temp_${raw_data6}_1.hh | sort | uniq > ${dirI}/temp_${raw_data6}.hh_snps; wc -l ${dirI}/temp_${raw_data6}.hh_snps
    #1681
  #do any fall in the hg19 PAR boundaries
  echo ""
  awk 'NR==FNR {a[$1];next} $2 in a' ${dirI}/temp_${raw_data6}.hh_snps ${dirMIA}/${raw_data6}.bim | awk '($1==23 && ($4 >= 60001 && $4 <= 2699520)) ||($1==23 && ($4 >= 154931044 && $4 <= 155260560)) || ($1==24 && ($4 >= 10001 && $4 <= 2649520)) || ($1==24 && ($4 >= 59034050 && $4 <= 59363566))' | wc -l
    #3, may as well remove
  echo ""

  #First pass alignment against 1KG found a number of SNPs that were remapped incorrectly by the Oncoarray consurtium files. E.g. it remaps chr2_238556837_A_INDEL to rs189474818, even though they are different SNPs. Where possible I have identified the correct ID - in each case there is a SNP and indel at the same position.
  echo "rs189474818 chr2:238556837:I" > ${dirI}/temp_oncoarray_correction #freq matches existing indel at that position
  echo "rs145123507 chr3:170072132:I" >> ${dirI}/temp_oncoarray_correction #freq matches existing indel at that position
  echo "rs184189083 chr8:76255548:I" >> ${dirI}/temp_oncoarray_correction #freq matches existing indel at that position
  echo "rs182565982 chr8:76654003:I" >> ${dirI}/temp_oncoarray_correction #freq matches existing indel at that position
  echo "rs3786878 chr19:38761733:I" >> ${dirI}/temp_oncoarray_correction #freq matches existing indel at that position

  #the following probably are just a misname but the freq is too different. E.g. 'rs183238081' is MAF 0.0007 in MIA, but expect 15% in EUR for matching chr7:23505149:I. 
  echo "rs183238081" >> ${dirI}/temp_${raw_data6}.hh_snps #matching indel has diff freq
  echo "rs12883217" >> ${dirI}/temp_${raw_data6}.hh_snps #Alleles don't make sense (A/T when should be T/G
  echo "rs186261221" >> ${dirI}/temp_${raw_data6}.hh_snps #no overlapping indel in 1kg phase 1 v3, drop
  echo "rs60682269" >> ${dirI}/temp_${raw_data6}.hh_snps #overlapping indel has diff freq (this is 1% in MIA, expect 47%)

  plink_1.90 --threads 1 --bfile ${dirMIA}/${raw_data6} --exclude ${dirI}/temp_${raw_data6}.hh_snps --set-hh-missing --make-bed --out ${dirI}/temp_2016_cleaning_MIA_KK_cases --silent --update-name ${dirI}/temp_oncoarray_correction
    #493234 variants loaded --exclude: 491549 variants remaining. --update-name: 5 values updated.  Total genotyping rate is 0.995059 491549 variants and 1968 people pass filters and QC.

  rm -f ${dirI}/temp_${raw_data6}.hh_snps ${dirI}/temp_${raw_data6}_1.*  ${dirI}/temp_oncoarray_correction
  echo -e "\nMIA data extracted"

done 

########################################################################
#2.4 EPIGENE + QTWIN
########################################################################
#this is based on /mnt/lustre/home/matthewL/Melanoma/WHITEMAN/EPIGENE/SCRIPTS/script_to_identify_epigene_control_release8.sh but skips what little cleaning/aligning I did there as doing it again here the same way for all sets
  #see /mnt/lustre/home/matthewL/Melanoma/WHITEMAN/EPIGENE/SCRIPTS/script_to_identify_epigene_control_release8.sh for how I identified cases and controls for the EPIGENE data from the most recent imputation batch

for x in #1 
do
  echo -e "\nExtracting EPIGENE and QTWIN controls from release 8 data"
  #work out a N of IDs
  wc -l ${dirJEO}/GWAS_b37PlusStrand_chr10.fam
    #12089
  echo ""
  wc -l ${dirJE}/EPIGENE_hg19_cleaned_no_maf_filter_IBD_PCA_outliers_removed.fam
    #782
  echo ""
  grep -c DWM ${dirJEO}/GWAS_b37PlusStrand_chr10.fam
    #784
  echo ""
  awk 'NR==FNR {a[$1];next} $1 in a' <(cut -f4 ${EPIGENE_IDS} ) <(grep DWM ${dirJEO}/GWAS_b37PlusStrand_chr10.fam | cut -d" " -f2 | sed 's/^DW//g') | wc -l
    #784; good start
  echo ""
  grep DWM ${dirJEO}/GWAS_b37PlusStrand_chr10.fam | cut -d" " -f1,2 > ${dirI}/temp_EPIGENE_controls_IDs; wc -l ${dirI}/temp_EPIGENE_controls_IDs
    #784
  echo ""
  #at the moment want parents of twins (as unrelated and only collected for P193 etc where we have ethics)
  cut -d" " -f2 ${dirJEO}/GWAS_b37PlusStrand_chr10.fam | grep "^8" | grep "03$\|04$" | wc -l
    #2805
  echo ""
  #see the epigene script for how I worked out which samples were twins, and which were suitable for inclusion. Here can just use ${EPIGENE_QTWIN_IDs}
  cat ${EPIGENE_QTWIN_IDs} >> ${dirI}/temp_EPIGENE_controls_IDs; wc -l ${dirI}/temp_EPIGENE_controls_IDs
    #1767

  #first extract; can make a merge list at the same time
  rm -f ${dirI}/temp_merge.list

  echo ""
  for i in {1..22} X
  do
    echo "Extracting chr${i} for EPIGENE/QTWIN"
    plink_1.90 --threads 1 --bfile ${dirJEO}/GWAS_b37PlusStrand_chr${i} --keep ${dirI}/temp_EPIGENE_controls_IDs --make-bed --out ${dirI}/temp_GWAS_b37PlusStrand_chr${i}_EPIGENE_controls_IDs --silent
    echo ${dirI}/temp_GWAS_b37PlusStrand_chr${i}_EPIGENE_controls_IDs >> ${dirI}/temp_merge.list
  done #extract IDs by chromosome
   #example for chr22 - Look like all SNPs, not just the ones I cleaned down to
     #12089 people (4004 males, 8085 females) loaded from .fam. --keep: 1767 people remaining.#Using 1 thread (no multithreaded calculations invoked).  #Before main variant filters, 2103 founders and 30 nonfounders present.#Calculating allele frequencies... done. #Total genotyping rate in remaining samples is 0.867375. #10152 variants and 1767 people pass filters and QC.

  #then merge into a single GWAS file, not by chr
  echo ""
  plink_1.90 --threads 1 --merge-list ${dirI}/temp_merge.list --make-bed --out ${dirI}/temp_${raw_data7}
  echo ""
  awk '{print $1}' ${dirI}/temp_${raw_data7}.bim | sort | uniq -c
    #all chr in there
  echo ""
  #rm the temp files
  rm -f ${dirI}/temp_GWAS_b37PlusStrand_chr*_EPIGENE_controls_IDs.* ${dirI}/temp_merge.list ${dirI}/temp_EPIGENE_controls_IDs

  ########################################################################
  #2.4.1 Fix the IDs, and set the phenotype
  ########################################################################
  #as per http://stackoverflow.com/questions/27252827/awk-replace-in-multiple-columns-specific-character-with-one-command Print an ID update file - without the DW field
  awk '{print $1,$2,$1,$2}' ${dirI}/temp_${raw_data7}.fam | awk ' BEGIN {a[3];a[4]} {for(x in a)gsub(/DW/,"",$x)}7' > ${dirI}/temp_${raw_data7}_ID_update; wc -l ${dirI}/temp_${raw_data7}_ID_update
    #1767
  echo ""
  tail ${dirI}/temp_${raw_data7}_ID_update
  echo ""
  awk '{
  if($2 ~ /DWM/)
    print $3,$4,"2";
  if($2 ~ /^8/)
    print $3,$4,"1";
  }' ${dirI}/temp_${raw_data7}_ID_update > ${dirI}/temp_${raw_data7}_pheno;wc -l ${dirI}/temp_${raw_data7}_pheno
    #1767
  echo ""
  tail ${dirI}/temp_${raw_data7}_pheno
    #looks good
  #fix IDs
  plink_1.90 --threads 1 --bfile ${dirI}/temp_${raw_data7} --update-ids ${dirI}/temp_${raw_data7}_ID_update --make-bed --out ${dirI}/temp_${raw_data7}_v2 --allow-no-sex --silent
  plink_1.90 --threads 1 --bfile ${dirI}/temp_${raw_data7}_v2 --pheno ${dirI}/temp_${raw_data7}_pheno --make-bed --out ${dirI}/temp_${raw_data7}_v3 --allow-no-sex --silent
  #Among remaining phenotypes, 784 are cases and 983 are controls.

  ########################################################################
  #2.4.2 Check for sex problems, hh SNPs, write out as case/controls
  ########################################################################
  plink_1.90 --threads 1 --bfile  ${dirI}/temp_${raw_data7}_v3 --check-sex --out ${dirI}/temp_${raw_data7}_v3 --silent
    #--check-sex: 10039 Xchr and 0 Ychr variant(s) scanned, 2 problem detected.
  echo -e "\nCheck sex results for ${raw_data7} : showing `grep -c "PROBLEM" ${dirI}/temp_${raw_data7}_v3.sexcheck` PROBLEMS. Removing failed, correcting missing\n"
    #2
  awk '$4==0 {print $1,$2}' ${dirI}/temp_${raw_data7}_v3.sexcheck > ${dirI}/temp_${raw_data7}_v3.sexcheck_remove; wc -l ${dirI}/temp_${raw_data7}_v3.sexcheck_remove
    #2
  echo ""
  awk '$5=="PROBLEM" && $4!="0" {print $1,$2,$4}' ${dirI}/temp_${raw_data7}_v3.sexcheck > ${dirI}/temp_${raw_data7}_v3.sexcheck_update; wc -l ${dirI}/temp_${raw_data7}_v3.sexcheck_update
    #0
  plink_1.90 --threads 1 --bfile  ${dirI}/temp_${raw_data7}_v3 --make-bed --out ${dirI}/temp_${raw_data7}_v2 --remove ${dirI}/temp_${raw_data7}_v3.sexcheck_remove --update-sex ${dirI}/temp_${raw_data7}_v3.sexcheck_update --silent
   #1767 phenotype values loaded from .fam. --update-sex: 0 people updated --remove: 1765 people remaining. 577849 variants and 1765 people pass filters and QC.
 
  #Note getting any hh or het Y chromosome errors so write out as it is - 
  plink_1.90 --threads 1 --bfile ${dirI}/temp_${raw_data7}_v2 --filter-cases --make-bed --out ${dirI}/temp_2016_cleaning_EPIGENE_cases --silent
    #577849 variants and 784 people pass filters and QC.
  plink_1.90 --threads 1 --bfile ${dirI}/temp_${raw_data7}_v2 --filter-controls --make-bed --out ${dirI}/temp_2016_cleaning_QTWIN_controls --silent
    #577849 variants and 981 people pass filters and QC.
  rm -f ${dirI}/temp_${raw_data7}* ${dirI}/temp_controls_IDs ${dirI}/EPIGENE_QTWIN_IDs.txt
  echo -e "\nEPIGENE + QTWIN extracted to case and control, with sex errors fixed (no hh errors)."
done #loop to turn on/off EPIGENE/QTWIN data


########################################################################
#3.0 Lift over data, and then fix any mismatches with 1KG
########################################################################
#15/06/2016 I am not quite sure about how to do this. If I seperate them into analysis subsets (AMFS cases, AMFS controls, QMEGA_610k_cases, QMEGA_omni_cases, 610k controls, SDH controls before lift over and aligned to 1kG I have a better chance of catching cryptic ambigious SNPs, since I want to keep those in, (as in if I try align using sam melanoma.bim I wont pick up if it is flipped in one set but not the other.... )
  #not quite sure how the missing SNP (e.g. SNPs not on 610k array thus NA in QMEGA_610k) are going to fall out

#THIS MAY SEEM WASTEFUL to do for the seperate sets taken from the same melanoma.bim file but I can't assume that they were perfectly merged (e.g. ambigious SNPs were properly combined

#_______SHOULD I CHECK FOR non ambigious maf vs 1kg maf - is it worth flagging/dropping SNPs wildly different from 1kg? I don't know....

echo -e "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "Aligning sets with 1KG, including ambigious SNPs where MAF  < ${MAF}"
echo -e "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

for i in ${dirI}/temp_2016_cleaning_MIA_KK_cases  #${dirI}/temp_2016_cleaning_QTWIN_controls ${dirI}/temp_2016_cleaning_EPIGENE_cases ${dirI}/temp_2016_cleaning_IBD_controls ${dirI}/temp_2016_cleaning_WAMHS_cases ${dirI}/temp_2016_cleaning_OAGphase2_controls ${dirI}/temp_2016_cleaning_HEIDELBERG_cases ${dirI}/temp_2016_cleaning_HEIDELBERG_controls ${dirI}/temp_2016_cleaning_SDH_2 ${dirI}/temp_2016_cleaning_QMEGA_610k_cases ${dirI}/temp_2016_cleaning_QMEGA_610k_controls ${dirI}/temp_2016_cleaning_QMEGA_omni_cases ${dirI}/temp_2016_cleaning_AMFS_cases ${dirI}/temp_2016_cleaning_AMFS_controls ${dirI}/temp_2016_cleaning_OAGphase1_controls
do

  out_name=`basename ${i}`

  #SDH is already hg19, as are IBD and WAMHS datasets
  if `echo ${out_name} | grep -q "SDH\|IBD\|WAMHS\|QTWIN\|EPIGENE\|MIA"`
  then
    echo -e "\nHg19 data detected, skipping lift over"
    #write out a version to make it work with the rest of the script
    plink_1.90 --threads 1 --bfile ${i} --silent --make-bed --out ${dirI}/temp_${out_name}_cleaning_4

  else
    echo -e "\n++++++++++++++++++++Hg18 dataset ${out_name} found, using Liftover++++++++++++++\n"

    awk ' {print "chr"$1"\t"$4"\t"$4+1"\t"$2}' ${i}.bim | sed 's/chr23/chrX/g' | sed 's/chr24/chrY/g' > ${dirI}/temp_${out_name}_BED
    head ${dirI}/temp_${out_name}_BED
    echo ""

    # as per wiki liftOver input.bed hg18ToHg19.over.chain.gz output.bed unlifted.bed
    liftOver ${dirI}/temp_${out_name}_BED ${dirL}/hg18ToHg19.over.chain.gz ${dirI}/temp_${out_name}_BED_output ${dirI}/temp_${out_name}_BED_unlifted
    wc -l ${i}.bim
      #1034058 melanoma, 610k
    echo ""
    wc -l ${dirI}/temp_${out_name}_BED
      #1034058 melanoma, 610k
    echo ""
    wc -l ${dirI}/temp_${out_name}_BED_output
      #1034058 melanoma, 610k
    echo ""
    wc -l ${dirI}/temp_${out_name}_BED_unlifted
      #2154; 2 per line so 1077 SNPs dropped, sounds about right, 610k
    echo ""
    grep -c Deleted ${dirI}/temp_${out_name}_BED_unlifted
      #1077 melanoma, 610k
    echo ""
    #make a removal list; keep a copy
    grep -v Deleted ${dirI}/temp_${out_name}_BED_unlifted | cut -f4 > ${dirK}/${out_name}_hg18_snps_not_lifted_exlude.txt; wc -l ${dirK}/${out_name}_hg18_snps_not_lifted_exlude.txt

    awk '{print $4,$2}' ${dirI}/temp_${out_name}_BED_output > ${dirI}/temp_${out_name}_update_map_bp.txt
    awk '{print $4,$1}' ${dirI}/temp_${out_name}_BED_output | sed 's/chrX/chr23/g' | sed 's/chrY/chr24/g'| sed 's/chr//g' > ${dirI}/temp_${out_name}_update_map_chr.txt

    #first drop SNPs and update bp postions. 07/07/2016 these are all straighforward; add --silent
    echo ""
    plink_1.90 --threads 1 --bfile ${i} --exclude ${dirK}/${out_name}_hg18_snps_not_lifted_exlude.txt  --make-bed --out ${dirI}/temp_${out_name}_cleaning_1 --silent
      #--exclude: 1041929 variants remaining, melanoma, 610k etc
    plink_1.90 --threads 1 --bfile ${dirI}/temp_${out_name}_cleaning_1 --update-map  ${dirI}/temp_${out_name}_update_map_bp.txt --make-bed --out ${dirI}/temp_${out_name}_cleaning_2 --silent
      #--update-map: 1041929 values updated.
    #remember you can only update one part of the map file at a time and need to add a flag to specify chr
    plink_1.90 --threads 1 --bfile  ${dirI}/temp_${out_name}_cleaning_2 --update-chr  ${dirI}/temp_${out_name}_update_map_chr.txt --make-bed --out  ${dirI}/temp_${out_name}_cleaning_3 --silent
      #--update-chr: 1041929 values updated.: melanoma, 610k, 610k_controls, omni_cases
    #now split-x to convert X to X_nonpar or will bork later allele freq test. This may be imperfect (may not capture SNPs not lifted over prop and thus not in the right region...).... didn't work when I used hg19 - said no SNPs in this region - suggesting it isn't lifting them over correctly? No it look like plink and 1000 genomes disagree - it is listing SNPs out as far as <5 million as X_nonpar. Use no-fail to force it, may have to fix downstream

    plink_1.90 --threads 1  --bfile  ${dirI}/temp_${out_name}_cleaning_3 --split-x hg19 no-fail --make-bed  --out  ${dirI}/temp_${out_name}_cleaning_4 

  fi #finished remapping with Liftover
  echo -e "\nNow data is hg19 further checking mapping and alignment against 1000 genomes phase 1 v3"

   #19/07/2016 as they keep causing me problems going to remove NA SNPs - SNPs with 0 chromosomes. You can't use maf though - they have no MAF so pass any maf setting
  echo -e "\nDropping SNPs completely missing in the dataset. Also a number of SNPs have incorrect map positions across the different files; correcting them all now"
   
  #need the freq data if I am going to try align ambigious SNPs, then merge with the bp data from the bim file to work out if the positions are wrong (liftover gets some SNPs wrong, or doesn't lift over every SNP. 
  plink_1.90 --threads 1 --bfile ${dirI}/temp_${out_name}_cleaning_4 --freq --out ${dirI}/temp_${out_name}_cleaning_4 --silent 
 
  #combine in the fix exclude SNPs as well (see section 0.5 for how these are defined)
  awk '$6==0 {print $2}' ${dirI}/temp_${out_name}_cleaning_4.frq > ${dirI}/temp_${out_name}_cleaning_4.frq_missing
  cat fix_exclude >> ${dirI}/temp_${out_name}_cleaning_4.frq_missing
    
  echo ""
  plink_1.90 --threads 1 --bfile ${dirI}/temp_${out_name}_cleaning_4 --exclude ${dirI}/temp_${out_name}_cleaning_4.frq_missing --make-bed --out ${dirI}/temp_${out_name}_cleaning_4_temp --update-map ${dirI}/fix_remap --update-chr  ${dirI}/fix_remap_chr
    #SDH: 1133951 variants loaded from .bim file. --update-map: 318 --exclude: 997451 variants remaining. --update-chr: 15 values updated
    #610k_cases: 1032981 variants loaded  --update-map: 223 --exclude: 549288 variants remaining. 549288 variants and 925 people pass filters and QC. --update-chr: 6
    #610k_controls: 1032981 variants loaded from .bim file. --update-map: 223 --exclude: 533839 variants remaining. --update-chr: 6 values updated
    #QMEGA_omni_cases: 1032981 variants loaded from .bim file --update-map: 223 --exclude: 807695 variants remaining. --update-chr: 10 values updated
    #AMFS_cases: 1032981 variants loaded from .bim file --update-map: 223 --exclude: 807695 variants remaining --update-chr: 10 values updated
    #AMFS_controls: 1032981 variants loaded from .bim file --update-map: 223 --exclude: 807695 variants remaining. --update-chr: 10 values updated
    #WAMHS 950792 variants loaded from .bim file --update-map: 302 --exclude: 937937 variants remaining
    #IBD 948342 variants loaded from .bim file --update-map: 302 --exclude: 935488 variants remaining. 
    #Heidelberg_cases: 221764 variants loaded from .bim file --update-map: 76 --exclude: 219349 variants remaining
    #HEI_con: 221764 variants loaded from .bim file. --update-map: 76 --exclude: 221711 variants remaining. --update-chr: 3 values updated
    #OAG1: 961293 variants loaded from .bim file. --update-map 301 --exclude: 932876 variants remaining. --update-chr: 15 
    #OAG2 731194 variants loaded from .bim file. --update-map: 321 --exclude: 731194 variants remaining.
    #QTWIN 577849 variants loaded from .bim file. --update-map: 26 values updated, 299 variant IDs not present. --exclude: 577253 variants remaining. --update-chr: 8 
    #EPIGENE 577849 variants loaded --update-map: 26 values updated --exclude: 500184 variants remaining --update-chr: 8 values updated
    #MIA 491549 variants loaded --update-map: 28 values updated --exclude: 491548 variants remaining --update-chr: 3 values updated
  plink_1.90 --threads 1 --bfile ${dirI}/temp_${out_name}_cleaning_4_temp --make-bed --out ${dirI}/temp_${out_name}_cleaning_4 --silent
  plink_1.90 --threads 1 --bfile ${dirI}/temp_${out_name}_cleaning_4 --freq --out ${dirI}/temp_${out_name}_cleaning_4 --silent
  echo ""
  awk 'NR==FNR {a[$2]=$1" "$3" "$4" "$5" "$9" "$13;next} $2 in a {print $0,a[$2]}' <(gunzip -dc ${dirI}/ALL_1000G_phase1integrated_v3_ALL_chr_legends_with_chromosome.txt.gz) <( awk 'NR==FNR {a[$2]=$4;next} $2 in a {print $0,a[$2]}' <(sed 's/ \+/ /g' ${dirI}/temp_${out_name}_cleaning_4.bim | sed 's/^ //g')  <(sed 's/ \+/ /g' ${dirI}/temp_${out_name}_cleaning_4.frq | sed 's/^ //g' | sed 's/ $//g') ) > ${dirI}/temp_${out_name}_cleaning_4_1KG_merge_rsIDs; wc -l ${dirI}/temp_${out_name}_cleaning_4_1KG_merge_rsIDs
  #SDH 883401 610k_cas 543678, 610k_con 528926, Omni_cases, AMFS_cas and AMFS_con 765318 WAMHS 679031 IBD 679175 Hei_cas 217000 HEI_con 219244 OAG1 895037 OAG2 714247 IBD 679173 253699 EPIGENE 245796 MIA 435523

  echo "" 
  head -4 ${dirI}/temp_${out_name}_cleaning_4_1KG_merge_rsIDs

  cd $dirI
  echo -e "\n+++++++++++Using awk_check_based_on_ID_1KG function to check SNPs paired with 1KG by ID have correct position, and strand. Aligning amibigious SNPs where MAF < $MAF and MAF difference with 1KG is < $DIFF"
  #see 0.4 for the function used here but this will assign the SNPs to various categories - e.g. aligned_concordant - same ID, same position, same strand and so on
  awk -v MAF=${MAF} -v DIFF=${DIFF} -f /working/lab_stuartma/scratch/matthewL/melanoma_meta_analysis/SCRIPTS/awk_check_based_on_ID_1KG ${dirI}/temp_${out_name}_cleaning_4_1KG_merge_rsIDs
  #any duplicates? Suggests I messed up in the function as it is assigning the same SNP to more than one category. Not every iteration generates every file and that is okay, so supress the cat warning about missing files by 2> /dev/null. Tested it by adding a duplicate SNP to one of the files, still works fine. 
  var_dup_1=`cat temp_flipped_concordant temp_flipped_chr_disc_bp_conc temp_aligned_concordant temp_aligned_chr_disc_bp_disc temp_aligned_chr_disc_bp_conc temp_aligned_chr_conc_bp_disc temp_flipped_chr_con_bp_disc temp_ambigious_aligned_concordant_same_allele temp_ambigious_aligned_chr_disc_bp_conc_same_allele temp_ambigious_aligned_chr_conc_bp_disc_same_allele temp_ambigious_aligned_chr_disc_bp_disc_same_allele temp_ambigious_aligned_concordant_other_allele temp_ambigious_aligned_chr_disc_bp_conc_other_allele temp_ambigious_aligned_chr_conc_bp_disc_other_allele temp_ambigious_aligned_chr_disc_bp_disc_other_allele temp_ambigious_flipped_concordant_same_allele temp_ambigious_flipped_chr_disc_bp_conc_same_allele temp_ambigious_flipped_chr_conc_bp_disc_same_allele temp_ambigious_flipped_chr_disc_bp_disc_same_allele temp_ambigious_flipped_concordant_other_allele temp_ambigious_flipped_chr_disc_bp_conc_other_allele temp_ambigious_flipped_chr_conc_bp_disc_other_allele temp_ambigious_flipped_chr_disc_bp_disc_other_allele temp_ambigious_MAF_to_high_too_assess temp_non_biallelic_snps temp_indels_aligned_same_allele temp_indels_aligned_other_allele temp_flipped_chr_disc_bp_disc 2> /dev/null  | sort | uniq -d | wc -l`

  if [[ $var_dup_1 == 0 ]]
  then
    echo -e "\nNo duplicate SNPs found after running awk_check_based_on_ID_1KG, continuing\n"
  else
    echo -e "\nWarning. Running awk_check_based_on_ID_1KG has assigned the same SNP to more than one file. Exiting"
    exit
  fi  #Got various amounts as I found and fixed bugs, now melanoma, SDH give 0 dups

  #how many output in total by the function. supress the cat warning about missing files by 2> /dev/null.
  echo -e "For ${out_name} total of `cat temp_flipped_concordant temp_flipped_chr_disc_bp_conc temp_aligned_concordant temp_aligned_chr_disc_bp_disc temp_aligned_chr_disc_bp_conc temp_aligned_chr_conc_bp_disc temp_flipped_chr_con_bp_disc temp_ambigious_aligned_concordant_same_allele temp_ambigious_aligned_chr_disc_bp_conc_same_allele temp_ambigious_aligned_chr_conc_bp_disc_same_allele temp_ambigious_aligned_chr_disc_bp_disc_same_allele temp_ambigious_aligned_concordant_other_allele temp_ambigious_aligned_chr_disc_bp_conc_other_allele temp_ambigious_aligned_chr_conc_bp_disc_other_allele temp_ambigious_aligned_chr_disc_bp_disc_other_allele temp_ambigious_flipped_concordant_same_allele temp_ambigious_flipped_chr_disc_bp_conc_same_allele temp_ambigious_flipped_chr_conc_bp_disc_same_allele temp_ambigious_flipped_chr_disc_bp_disc_same_allele temp_ambigious_flipped_concordant_other_allele temp_ambigious_flipped_chr_disc_bp_conc_other_allele temp_ambigious_flipped_chr_conc_bp_disc_other_allele temp_ambigious_flipped_chr_disc_bp_disc_other_allele temp_ambigious_MAF_to_high_too_assess temp_non_biallelic_snps temp_indels_aligned_same_allele temp_indels_aligned_other_allele temp_flipped_chr_disc_bp_disc  2> /dev/null | wc -l` SNPs output by awk_check_based_on_ID_1KG\n"
    #SDH - 876722 I wonder if many simply have NA freq in the SDH samples (which is from a few different arrays I think). So melanoma.bim: missing ~58. 9 are autosome ambigious alleles with allele freq diff > 10% between here and 1kg, which is probably just chance. Most are just over 10% but diminishing returns to chase them all up. The rest are non par vs chrX which should be captured by $1!=$8?. Oh wait, wont the allele freq be wrong if they are X_nonPAR? Can't readily fix these though as plink disagrees as to what the X_nonPAR regions are and there isnt much point me trying to resolve as the HRC doesn't handle X yet. Actually no, most are nonPAR and ambigious with >10% freq, so should be dropped anyway. SDH is missing ~ 37781. For 610k_cases 543669, 610k_controls 528913. QMEGA_omni_cases - 765269 AMFS_cases 765268 AMF_controls 765270 WAMHS 679022. IBD: 679166. Hei_cases" 216619. HEI_con: 218858. OAG1: 844351. OAG2 699012 IBD 679164 QTWINS 253572 EPIGENE 245678 MIA 435327

  #which ones are missing? supress the cat warning about missing files by 2> /dev/null.
  awk 'NR==FNR {a[$1];next} !($2 in a) ' <(cat temp_flipped_concordant temp_flipped_chr_disc_bp_conc temp_aligned_concordant temp_aligned_chr_disc_bp_disc temp_aligned_chr_disc_bp_conc temp_aligned_chr_conc_bp_disc temp_flipped_chr_con_bp_disc temp_ambigious_aligned_concordant_same_allele temp_ambigious_aligned_chr_disc_bp_conc_same_allele temp_ambigious_aligned_chr_conc_bp_disc_same_allele temp_ambigious_aligned_chr_disc_bp_disc_same_allele temp_ambigious_aligned_concordant_other_allele temp_ambigious_aligned_chr_disc_bp_conc_other_allele temp_ambigious_aligned_chr_conc_bp_disc_other_allele temp_ambigious_aligned_chr_disc_bp_disc_other_allele temp_ambigious_flipped_concordant_same_allele temp_ambigious_flipped_chr_disc_bp_conc_same_allele temp_ambigious_flipped_chr_conc_bp_disc_same_allele temp_ambigious_flipped_chr_disc_bp_disc_same_allele temp_ambigious_flipped_concordant_other_allele temp_ambigious_flipped_chr_disc_bp_conc_other_allele temp_ambigious_flipped_chr_conc_bp_disc_other_allele temp_ambigious_flipped_chr_disc_bp_disc_other_allele temp_ambigious_MAF_to_high_too_assess temp_non_biallelic_snps temp_indels_aligned_same_allele temp_indels_aligned_other_allele temp_flipped_chr_disc_bp_disc 2> /dev/null  ) ${dirI}/temp_${out_name}_cleaning_4_1KG_merge_rsIDs > ${dirI}/temp_${out_name}_unmappable_snps; wc -l ${dirI}/temp_${out_name}_unmappable_snps 
    #SDH 6679; 610k_cases 9. 610k_control 13.  Omni_cases 49. AMFS_cases 50, AMFS_controls 48, all extreme allele diff, most on chrX. WAMHS, IBD 9. Hei_cases" 381, perhaps shouldn't be surprised is high given the old chip. HEI_con: 386. OAG1: 50686; all look to monoroprh ( e.g. 0 and A). OAG2 15235 all monomorh 0/A etc QTWIN 127 EPIGENE 118. MIA 196
  echo ""
  head ${dirI}/temp_${out_name}_unmappable_snps

  var_unmapple_non_mono_non_amb_1=`awk '!($3==0 || $4==0) && !($3=="A" && $4=="T" && $10=="T" && $11=="A") && !($3=="A" && $4=="T" && $10=="A" && $11=="T") && !($3=="T" && $4=="A" && $10=="T" && $11=="A") && !($3=="T" && $4=="A" && $10=="A" && $11=="T") && !($3=="C" && $4=="G" && $10=="G" && $11=="C") && !($3=="C" && $4=="G" && $10=="C" && $11=="G") && !($3=="G" && $4=="C" && $10=="G" && $11=="C") && !($3=="G" && $4=="C" && $10=="C" && $11=="G")' ${dirI}/temp_${out_name}_unmappable_snps | wc -l`
   
  if [[ $var_unmapple_non_mono_non_amb_1 == 0 ]]
  then
    echo -e "\nAll unmappable SNPs for ${out_name} are either monomorphic or ambigious, continuing with their removal and including additional SNPs to drop (high MAF ambigious, tri-alleleic etc..\n"
  else
    echo -e "\nWarning. For ${out_name} there are $var_unmapple_non_mono_non_amb_1 unmappable SNPs resulting from awk_check_based_on_ID_1KG function that should be mappable (are polymorphic and not ambigious or tri-allelic etc. Check ${dirI}/temp_${out_name}_unmappable_snps, exiting\n" 
    exit
  fi #check for unmappable SNPs that should be mappable

  cut -d" " -f2 ${dirI}/temp_${out_name}_unmappable_snps > ${dirI}/temp_${out_name}_unmappable_snps_final
  cat temp_ambigious_MAF_to_high_too_assess >> ${dirI}/temp_${out_name}_unmappable_snps_final
  cut -d" " -f2 temp_non_biallelic_snps >> ${dirI}/temp_${out_name}_unmappable_snps_final; wc -l ${dirI}/temp_${out_name}_unmappable_snps_final
     #SDH 10609 610k_cas 723 610k_con 670 Omni_cas 3819 AMFS_cas 3827 AMFS_con 3850 WAMHS 583 HEI_cas 404 HEI_con 409 OAG1 54567 OAG2 15818 IBD 579 QTWIN 305 EPIGENE 271 MIA 3402
  echo ""
  #so list of SNPs that need to be flipped; 27/06/2016 found a couple of potential flip files not included
  cat temp_flipped_concordant temp_flipped_chr_disc_bp_conc temp_flipped_chr_con_bp_disc temp_flipped_chr_disc_bp_disc temp_ambigious_flipped_concordant_same_allele temp_ambigious_flipped_chr_disc_bp_conc_same_allele temp_ambigious_flipped_chr_conc_bp_disc_same_allele temp_ambigious_flipped_chr_disc_bp_disc_same_allele temp_ambigious_flipped_concordant_other_allele temp_ambigious_flipped_chr_disc_bp_conc_other_allele  temp_ambigious_flipped_chr_conc_bp_disc_other_allele temp_ambigious_flipped_chr_disc_bp_disc_other_allele 2> /dev/null > temp_${out_name}_flip; wc -l temp_${out_name}_flip
   #SDH 436877 610k_cases 3864. 610k_controls 354. Now this is interesting. I was thinking that doing this (cases and controls sep) might detect ambigious SNPs that were incorrectly merged (merged as if strand compliant but actually ambigous). Could I potentially do the reverse though here? Incorrectly flip ambigious SNPs? Not sure..... QMEGA_omni_cases: 73056. AMFS_cases: 74625. AMFS_control: 73046. WAMHS 339500. HEI_cases:108460. HEI_controls: 109553. OAG1 420771 OAG2 349423 IBD 39 EPIGENE 3 MIA 215757

 #list of SNPs with chromosome wrong; don't run wc -l yet as most of these will just be differences in X chromosome naming convention - 1KG calls X_nonPAR that, where everyone else calls it just X; 27/06/2016 found I was missing some potential combinations
  cat temp_aligned_chr_disc_bp_conc temp_aligned_chr_disc_bp_disc temp_flipped_chr_disc_bp_conc temp_flipped_chr_disc_bp_disc temp_ambigious_aligned_chr_disc_bp_conc_same_allele temp_ambigious_aligned_chr_disc_bp_disc_same_allele temp_ambigious_aligned_chr_disc_bp_conc_other_allele temp_ambigious_aligned_chr_disc_bp_disc_other_allele temp_ambigious_flipped_chr_disc_bp_conc_same_allele temp_ambigious_flipped_chr_disc_bp_disc_same_allele temp_ambigious_flipped_chr_disc_bp_conc_other_allele temp_ambigious_flipped_chr_disc_bp_disc_other_allele 2> /dev/null > temp_${out_name}_chr_disc 

  echo -e "\nFor ${out_name} there are `awk 'NR==FNR {a[$1];next} $2 in a' temp_${out_name}_chr_disc ${dirI}/temp_${out_name}_cleaning_4_1KG_merge_rsIDs | awk '!($1==23 && $8=="X_nonPAR") && !($1==25 && $8=="X_nonPAR")' | wc -l` SNPs with chromosome remapping that is not due to naming conventions for chromosome X\n"
   #SDH 3 out of 21457 - all have 0 CHR and BP but based on MAF and alleles look correct in 1kG so should be remapped. rest chrX. 610k - 0. rs34908637 is in the omni set of data, and its freq suggests it is correct. rs11556134 has the wrong chr and position but the alleles and freqs look correct, but will also need flipping, but that is handled in the flip list above. QMEGA_omni_cases and AMFS_cases: 1. with WAMHS including a few chr25/X_nonPAR add that to the exclusion - was getting 190 for WAMHS and now 1. OAG1 0. OAG1 1
  awk 'NR==FNR {a[$1];next} $2 in a' temp_${out_name}_chr_disc ${dirI}/temp_${out_name}_cleaning_4_1KG_merge_rsIDs | awk '!($1==23 && $8=="X_nonPAR") && !($1==25 && $8=="X_nonPAR") {print $2,$8}' > temp_${out_name}_chr_disc_final
    #melanoma, 610k - two as above SDH just the one (0 rs2192400 A C 0.357 1140 0 2 71106880 A C 0.6174 0.3826) which looks right. WAMHS 1 IBD 1 MIA 0
  #the majority of the SNPs should be aligned and concordant, and further set correct but flipped
  wc -l temp_aligned_concordant
  echo ""
    #SDH 418893. 610k_cas 529341 610k_con 518427 Omni_cas, AMFS_cas 670488 AMFS_con 670489. WAMHS 330319 HEI_cas 104917 HEI_con 106053 OAG1 404623 OAG2 289313 IBD 661111 QTWIN 244966 EPIGENE 237790 MIA 202808
  wc -l temp_flipped_concordant 
    #SDH 419892. 610k_cas 3445, 610k_con 26, Omni_cas, AMFS_cas and AMFS_con 69963. WAMHS 330810 HEI_cas 105195 HEI_con 106278, OAG1 405620 OAG2 289738 IBD 1 QTWIN 0 (SG aligned) EPIGENE 0 MIA 202379
  echo ""
  #same with A/T
  wc -l temp_ambigious_aligned_concordant_same_allele
  echo ""
    #SDH 5303. 610k 1070. 610k_controls 931. QMEGA_omni_cases: 6482 (denser chip). AMFS_cases: 6481 (nice close to omni). AMFS_control: 6473. WAMHS 650. HEI_cases: 0; checked and def no ambigious SNPs in the original files supplies to me... HEI_con: 0. OAG1 4814 OAG2 8 IBD 1203 QTWIN 331 EPIGENE 304 MIA 9246
  wc -l temp_ambigious_aligned_concordant_other_allele
  echo ""
    #SDH 1300. 610k 297. 610k_con 259 Omni_cases 1769 AMFS_cas 1744 AMFS_con 1754. WAMHS 142 HEI_cas 0. HEI_con: 0. OAG1 1270 OAG2 4 IBD 277 QTWIN 106 EPIGENE 103 MIA 1085 
  wc -l temp_ambigious_flipped_concordant_same_allele
  echo ""
    #SDH 5132 610k_cas 201 610k_con 173 Omni_cases 1454 AMFS_cas 1457 AMFS_Con 1443 WAMHS 587 HEI_cas 0 HEI_con 0 OAG1 4725 OAG2 7 IBD 25 QTWIN 0 EPIGENE 0 MIA 9082
  wc -l temp_ambigious_flipped_concordant_other_allele
    #SDH 1331. 610k 85 610k_con 71. QMEGA_omni_cases: 587. AMFS_cases 602, AMFS_control: 588, WAMHS 146. HEI_cas 0. HEI_con 0 OAG1 1287 OAG2 2 IBD 12 QTWIN 0 EPIGENE 0 MIA 1079

  #+++++++++++++NEED TO THINK HOW TO WORK OUT IF DIFFERENTIAL FLIPPING IN CONTROLS VS CASE IS A GOOD THING OR NOT!
  #look at some of these; 15/06/2016 they all look good after some tweaking as above.; same for SDH; same for 610k
  echo -e "\nHead of temp_ambigious_aligned_concordant_same_allele:" 
  awk 'NR==FNR {a[$1];next} $2 in a' temp_ambigious_aligned_concordant_same_allele ${dirI}/temp_${out_name}_cleaning_4_1KG_merge_rsIDs | head
  echo -e "\nHead of temp_ambigious_aligned_concordant_other_allele:" 
  awk 'NR==FNR {a[$1];next} $2 in a' temp_ambigious_aligned_concordant_other_allele ${dirI}/temp_${out_name}_cleaning_4_1KG_merge_rsIDs | head
  echo -e "\nHead of temp_ambigious_flipped_concordant_same_allele:" 
  awk 'NR==FNR {a[$1];next} $2 in a' temp_ambigious_flipped_concordant_same_allele ${dirI}/temp_${out_name}_cleaning_4_1KG_merge_rsIDs | head
  echo -e "\nHead of temp_ambigious_flipped_concordant_other_allele:" 
  awk 'NR==FNR {a[$1];next} $2 in a' temp_ambigious_flipped_concordant_other_allele ${dirI}/temp_${out_name}_cleaning_4_1KG_merge_rsIDs | head

  #so what about the set with just bp differences? 27/06 - missing temp_ambigious_flipped_chr_conc_bp_disc_same_allele. check no more ; was missing some other potential matches
  cat temp_aligned_chr_conc_bp_disc temp_aligned_chr_disc_bp_disc temp_flipped_chr_con_bp_disc temp_flipped_chr_disc_bp_disc temp_ambigious_aligned_chr_conc_bp_disc_same_allele temp_ambigious_aligned_chr_disc_bp_disc_same_allele temp_ambigious_aligned_chr_conc_bp_disc_other_allele temp_ambigious_aligned_chr_disc_bp_disc_other_allele temp_ambigious_flipped_chr_conc_bp_disc_same_allele temp_ambigious_flipped_chr_disc_bp_disc_same_allele temp_ambigious_flipped_chr_conc_bp_disc_other_allele temp_ambigious_flipped_chr_disc_bp_disc_other_allele 2> /dev/null > temp_${out_name}_bp_disc

  echo -e "\nFor ${out_name} there are `sed -n '$=' temp_${out_name}_bp_disc` SNPs with different BP pos to 1KG"
    #SDH: 3 which all look to need remapping. 610k_cases - 153, 610k_controls 135 - weird N is diff, but was different genotyping arrays and runs, some SNPs might imperfectly overlap or be dropped out in one. QMEGA_omni_cases, AMFS_cases and AMFS_controls 399. WAMHS 1 - CIDR remapping I guess is pretty good. HEI_Cases and controls: 38. OAG1 162. OAG has 103958 seems massive! IBD 1 QTWIN 627 EPIGENE 612 MIA 2

  #what about indels? Only care about the ones that can be mapped handle - are non NA in freq
  var_indels=`awk 'NR==FNR {a[$1];next} $2 in a && $5!="NA"' <(cat temp_indels_aligned_same_allele temp_indels_aligned_other_allele  2> /dev/null) ${dirI}/temp_${out_name}_cleaning_4_1KG_merge_rsIDs | wc -l`

  echo -e "\nFor ${out_name} there are ${var_indels} indel SNPs that have non zero MAF, check to make sure these are being handled correctly\n"
  #SDH - none. 610k there are 0. QMEGA_omni_cases: 4, and can asses them. They look to be correct below. AMFS_cases has same 4 unsurprisingly. WAMHS, OAG1 and OAG2 IBD 0 MIA 44

  #if the var indels are in the dataset try fixing their allleles to match the 1kg. Need SNP oldA1 oldA2 newA1 newA2. Given 1kg and plink swap alelle order would need to print reverse order from 1kg if same allele, and 1kg allele order if over allele I think. Only do this for alleles with a non NA freq
  #08/09/2016 this function has worked fine so far for existing sets but the oncoarray file 'remap' has borked a handful of SNPs/indels by not checking if there is a SNP overlapping the indel. So need to add a step here that checks if the freqs make sense - if they do make the remap file, if not delete the SNP. 

  if [[ $var_indels -gt 0 ]]
  then
    echo -e "Indels that can be mapped to 1kg found; replacing I/D alleles with 1kg alleles\n"
    awk 'NR==FNR {a[$1];next} $2 in a && $5!="NA"' <(cat temp_indels_aligned_same_allele temp_indels_aligned_other_allele  2> /dev/null) ${dirI}/temp_${out_name}_cleaning_4_1KG_merge_rsIDs | head
      #QMEGA_omni_cases. These look correct, and again I = larger. So could correct these
      #22 rs11570768 D I 0.05988 1386 38507985 22 38507980 AGGCGG A 0.0752 0.0752
      #22 rs11570368 I D 0.01299 1386 38965330 22 38965330 T TC 0.0119 0.0119
      #22 rs28915376 D I 0.04467 1388 43091219 22 43091216 CAA C 0.0396 0.0396
      #22 rs3216827 I D 0.0879 1388 43091567 22 43091565 T TG 0.1003 0.1003
    if [[ -a temp_indels_aligned_same_allele ]]
    then
      #08/09/2016 old version that just assumed the indels were correct, which up until oncoarray data was the case
      ##awk 'NR==FNR {a[$1];next} $2 in a && $5!="NA" {print $2,$3,$4,$11,$10}' temp_indels_aligned_same_allele ${dirI}/temp_${out_name}_cleaning_4_1KG_merge_rsIDs >> ${dirI}/temp_${out_name}_cleaning_4_update_alleles
     cd $dirI
     rm -f temp_cleaning_4_update_indel_alleles_1 temp_cleaning_4_indels_diff_exclude_1
     awk 'NR==FNR {a[$1];next} {
       if($2 in a && $5!="NA" && (( $5 - $12 ) > (-'''${DIFF}''')) && (( $5 - $12 ) < '''${DIFF}''' ))
         {print  $2,$3,$4,$11,$10 > "temp_cleaning_4_update_indel_alleles_1"}
       if($2 in a && $5!="NA" && !((( $5 - $12 ) > (-'''${DIFF}''')) && (( $5 - $12 ) < '''${DIFF}''' )))
         {print  $2 > "temp_cleaning_4_indels_diff_exclude_1"}
       }' temp_indels_aligned_same_allele ${dirI}/temp_${out_name}_cleaning_4_1KG_merge_rsIDs
     echo ""
     wc -l temp_indels_aligned_same_allele
       #MIA 37
     echo ""
     wc -l temp_cleaning_4_update_indel_alleles_1
       #MIA 35
     echo ""
     wc -l temp_cleaning_4_indels_diff_exclude_1
       #MIA 2
       #22 rs5845861 I D 0.0005192 3852 49371352 22 49371352 C CG 0.3654 0.3654 <<actually looks to be the right indel, but freq too diff
       #22 rs35635786 I D 0.0005097 3924 49548451 22 49548451 C CT 0.3153 0.3153 <<actually looks to be the right indel, but freq too diff
     echo ""
     
    fi
   
    if [[ -a temp_indels_aligned_other_allele ]]
    then
      #08/09/2016 old command
      ##awk 'NR==FNR {a[$1];next} $2 in a && $5!="NA" {print $2,$3,$4,$10,$11}' temp_indels_aligned_other_allele ${dirI}/temp_${out_name}_cleaning_4_1KG_merge_rsIDs >> ${dirI}/temp_${out_name}_cleaning_4_update_alleles
    rm -f temp_cleaning_4_update_indel_alleles_2 temp_cleaning_4_indels_diff_exclude_2
     awk 'NR==FNR {a[$1];next} {
       if($2 in a && $5!="NA" && (( $5 - (1 - $12) ) > (-'''${DIFF}''')) && (( $5 - ( 1 - $12) ) < '''${DIFF}''' ))
         {print  $2,$3,$4,$11,$10 > "temp_cleaning_4_update_indel_alleles_2"}
       if($2 in a && $5!="NA" && !((( $5 - (1 - $12) ) > (-'''${DIFF}''')) && (( $5 - (1 - $12) ) < '''${DIFF}''' )))
         {print  $2 > "temp_cleaning_4_indels_diff_exclude_2"}
       }' temp_indels_aligned_other_allele ${dirI}/temp_${out_name}_cleaning_4_1KG_merge_rsIDs
     wc -l temp_indels_aligned_other_allele
       #MIA 7
     echo ""
     wc -l temp_cleaning_4_update_indel_alleles_2
       #MIA 6
     echo ""
     wc -l temp_cleaning_4_indels_diff_exclude_2
       #MIA 1
       #22 rs57835472 I D 0.03872 3926 26484530 22 26484529 AC A 0.0158 0.0158 <<something badly wrong here. alleles inverted (if is the right SNP major = minor) best to drop
 
    fi
    echo ""
    #08/09/2016 so combine the sets
    cat temp_cleaning_4_update_indel_alleles_1 temp_cleaning_4_update_indel_alleles_2 2> /dev/null > ${dirI}/temp_${out_name}_cleaning_4_update_alleles; wc -l ${dirI}/temp_${out_name}_cleaning_4_update_alleles
      #MIA 41
    echo ""
    cat ${dirI}/temp_${out_name}_cleaning_4_update_alleles | head 
      #QMEGA_omni cases - looks good, in every case I have seen I is larger than D
      #rs11570768 D I A AGGCGG
      #rs11570368 I D TC T
      #rs28915376 D I C CAA
      #rs3216827 I D TG T
    #08/09/2016 combine the new deletion file
    cat temp_cleaning_4_indels_diff_exclude_1 temp_cleaning_4_indels_diff_exclude_2 2> /dev/null > ${dirI}/temp_${out_name}_cleaning_4_exclude_indels
    echo -e "\nThere are `sed -n '$=' ${dirI}/temp_${out_name}_cleaning_4_exclude_indels` indels that differ from 1KG enough to delete\n"
    #update alleles, write to temp, then overwrite original 
    plink_1.90 --threads 1 --bfile ${dirI}/temp_${out_name}_cleaning_4 --update-alleles ${dirI}/temp_${out_name}_cleaning_4_update_alleles --make-bed --out ${dirI}/temp_${out_name}_cleaning_4_temp --exclude ${dirI}/temp_${out_name}_cleaning_4_exclude_indels
      #QMEGA_omni_case, AMFS_cases, AMFS_controls: --update-alleles: 4 variants updated. 807695 variants and 694 people pass filters and QC. MIA --update-alleles: 41
    echo ""
    plink_1.90 --threads 1 --bfile ${dirI}/temp_${out_name}_cleaning_4_temp --make-bed --out ${dirI}/temp_${out_name}_cleaning_4
    rm -f temp_cleaning_4_update_indel_alleles_1 temp_cleaning_4_update_indel_alleles_2 temp_cleaning_4_indels_diff_exclude_1 temp_cleaning_4_indels_diff_exclude_2 ${dirI}/temp_${out_name}_cleaning_4_exclude_indels ${dirI}/temp_${out_name}_cleaning_4_update_alleles
    echo -e "\nFinished updating the alleles for indels\n"
  fi

  #grab the updated positions
  awk 'NR==FNR {a[$1];next} $2 in a {print $2,$9}' temp_${out_name}_bp_disc ${dirI}/temp_${out_name}_cleaning_4_1KG_merge_rsIDs > temp_${out_name}_bp_disc_final

  #so now need to fix chr, fix position, flip. I think you can do it all at once
  plink_1.90 --threads 1 --bfile ${dirI}/temp_${out_name}_cleaning_4 --update-map temp_${out_name}_bp_disc_final --update-chr temp_${out_name}_chr_disc_final --flip temp_${out_name}_flip --make-bed --out ${dirI}/temp_${out_name}_cleaning_5 --exclude  ${dirI}/temp_${out_name}_unmappable_snps_final
    #SDH --update-map: 3. --flip: 436877 --exclude: 986842 variants remaining  #--update-chr: 3 . 986842 variants and 570 people pass filters and QC. --update-chr: 3 
    #610k_cases --update-map: 153 values updated --flip: 3864 --exclude: 548565 variants remaining. 548565 variants and 925 people pass filters and QC. --update-chr: 1
    #610k_controls --update-map: 135 values updated. --flip: 354  --exclude: 533169 variants remaining. --update-chr: 0. 533169 variants and 3954 people pass filters and QC.
    #QMEGA_omni_cases: --update-map: 399 updated --flip: 73056  --exclude: 803876 variants remaining --update-chr: 1. 803876 variants and 694 people pass filters and QC
    #AMFS_cases: --update-map: 399 updated. --flip: 73075 --exclude: 803868 variants remaining. --update-chr: 1. 803868 variants and 548 people pass filters and QC.
    #AMFS_controls: --update-map: 399 values updated. --flip: 73046. --exclude: 803845 variants remaining. --update-chr: 1. 803845 variants and 430 people pass filters and QC.
    #WAMHS --update-map: 1 value updated --flip: 339500 SNPs flipped --exclude: 937354 variants remaining. --update-chr 1. 937354 variants and 1273 people pass filters and QC
    #HEI_cases: --update-map: 38 values updated. --flip: 108460 SNPs flipped. --exclude: 218945 variants remaining.  --update-chr 0
    #HEI_controls: --update-map: 38 values updated. --flip: 109553 SNPs flipped. --exclude: 221302 variants remaining  --update-chr 0
    #OAG1: --update-map: 162 values updated. --flip: 420771 SNPs flipped. --exclude: 878309 variants remaining. --update-chr: 0 
    #OAG2: --update-map: 103958 values updated. --flip: 349423 SNPs flipped. --exclude: 715376 variants remaining. --update-chr: 1 
    #IBD --update-map: 1 value updated. --flip: 39 SNPs flipped. --exclude: 934909 variants remaining. --update-chr: 1
    #QTWIN --update-map: 627 values updated. --flip: 3 SNPs flipped. --exclude: 576948 variants remaining.
    #EPIGENE --update-map: 612 values updated  --flip: 3 SNPs flipped --exclude: 499913 variants remaining. --update-chr: 3 values updated.
    #MIA  --update-map: 2  --flip: 215757 SNPs flipped --exclude: 488143 variants remaining --update-chr: 0 values updated.
  echo ""
  awk '($1==0 && $4==0)  || $5==0 || $6==0 {print $2}' ${dirI}/temp_${out_name}_cleaning_5.bim > ${dirI}/temp_${out_name}_cleaning_overlapping_diff_alleles  ; wc -l ${dirI}/temp_${out_name}_cleaning_overlapping_diff_alleles
    #SDH 46508. 610k, 610k_controls - 0. QMEGA_omni_cases: 0, AMFS_cases and control: 0. WAMHS 799. HEI_cases: 241, HEI_cont: 284. OAG1 22026. OAG2 7839 IBD 1301 QTWIN 117773 (unfiltered exm data) EPIGENE 93143 MIA 0

  echo -e "\nAttempting merge duplicated SNPs (so same SNP on the array under two different names. Some requiring flipping and there is no easy way to do this manually other than extracting errors from plink log and restarting...\n"

  #given the damn exm overlapping content - 10000s of overlapping SNPs - have to have a strucutred way of testing for flips (as some overlapping SNPs are flipped in say WAMHS but not in IBD...)
  #just thought of a kludgy way to handle triplicates, quads etc. Pull them out, sort by missingness, and keep the first two with the lowest missingness and drop them into the standard function. Will sometimes lose SNPs (imagine overlapping indels and SNPs that are in turn duplictaes - keep a SNP and an indel, then lose the indel in the following step - but writing a comparison funciton for three or more pairings will take hours for 50 sites in WAMHS..... Then filter these out, which will just leave duplicates, and those I can handle relatively straightforwardly
  plink_1.90 --threads 1 --bfile ${dirI}/temp_${out_name}_cleaning_5 --freq --out ${dirI}/temp_${out_name}_cleaning_5 --silent
    #so pair triplicate sites with freq, then sort by position (chr then bp), then by chromosome count. If you then get awk to drop every first line will be keeping the 2 most informative SNPs. (e.g. drop 1, 4, 7 etc) 
  awk 'NR==FNR {a[$2]=$0;next} $2 in a {print a[$2],$0}' <( awk 'NR==FNR {a[$1,$2];next} ($1,$4) in a {print $0}' <(awk '{print $1,$4}' ${dirI}/temp_${out_name}_cleaning_5.bim | sort | uniq -c | awk '$1==3 {print $2,$3}') <( sed 's/\t/ /g' ${dirI}/temp_${out_name}_cleaning_5.bim | sed 's/ \+/ /g') ) <( sed 's/ \+/ /g' ${dirI}/temp_${out_name}_cleaning_5.frq | sed 's/^ //g') | sort -k1,1 -k4,4 -k12,12 | head -6
  echo ""
   #IBD 
      #1 exm2250216 0 169513583 T G 1 exm2250216 T G 0.0819 1856
      #1 exm121943 0 169513583 T G 1 exm121943 T G 0.08181 1858
      #1 rs6037 0 169513583 T G 1 rs6037 T G 0.08181 1858
      #1 rs17650204 0 171755071 G A 1 rs17650204 G A 0.09537 1856
      #1 exm123761 0 171755071 G A 1 exm123761 G A 0.09526 1858
      #1 exm2250237 0 171755071 G A 1 exm2250237 G A 0.09526 1858
  #print the first (least informative SNP) for exclusion
  awk 'NR==FNR {a[$2]=$0;next} $2 in a {print a[$2],$0}' <( awk 'NR==FNR {a[$1,$2];next} ($1,$4) in a {print $0}' <(awk '{print $1,$4}' ${dirI}/temp_${out_name}_cleaning_5.bim | sort | uniq -c | awk '$1==3 {print $2,$3}') <( sed 's/\t/ /g' ${dirI}/temp_${out_name}_cleaning_5.bim | sed 's/ \+/ /g') ) <( sed 's/ \+/ /g' ${dirI}/temp_${out_name}_cleaning_5.frq | sed 's/^ //g') | sort -k1,1 -k4,4 -k12,12 | sed -n 1~3p | cut -d" " -f2 >> ${dirI}/temp_${out_name}_cleaning_overlapping_diff_alleles  ; wc -l ${dirI}/temp_${out_name}_cleaning_overlapping_diff_alleles
    #IBD 1350 - 1301 + 49 trip SNPs. WAMHS 849. OAG2 7838 OAG1 22026 HEI_Cas 241 HEI_con 284 SDH 46508 610k_case, 610k_con, Omni_case, AMFS_case, AMFS_con 0 (no more) QTWIN 117797 ~ 24 more or so EPIGENE 93160, 17 more MIA 0
  echo ""
  #first get rid of 0 SNPs (and now the least informative triplicate SNP
 plink_1.90 --threads 1 --bfile ${dirI}/temp_${out_name}_cleaning_5 --exclude ${dirI}/temp_${out_name}_cleaning_overlapping_diff_alleles --make-bed --out ${dirI}/temp_${out_name}_cleaning_6
    #IBD --exclude: 933559 variants remaining WAMHS --exclude: 936505 variants remaining. OAG2 --exclude: 707537 variants remaining OAG1 --exclude: 856283 variants remaining HEI_cas --exclude: 218704 variants remaining HEI_con --exclude: 221018 variants remaining SDH --exclude: 940334 variants remaining 610k_case --exclude: 548565 variants remaining. 610k_con --exclude: 533169 variants remaining. Omni_cases --exclude: 803876 variants remaining AMFS_case --exclude: 803868 variants remaining --exclude: 803845 variants remainingQTWIN --exclude: 459153 variants remaining. EPIGENE --exclude: 406754 variants remaining. EPIGENE - at least 1 duplicate ID i nthe exclude list??? MIA  - no change, 0
  #generate freq
  plink_1.90 --threads 1 --bfile ${dirI}/temp_${out_name}_cleaning_6 --freq --out ${dirI}/temp_${out_name}_cleaning_6 --silent
  #now look for anything more than duplicated sites; if so have messed up
  var_non_dup=`awk '{print $1,$4}' ${dirI}/temp_${out_name}_cleaning_6.bim | sort | uniq -c | awk '$1>2 {print $2,$3}' | wc -l`
  if [[ $var_non_dup -gt 0 ]]
  then
    echo -e "\nWarning there remain sites in ${out_name} where they are more than 2 SNPs overlapping the same site; these cannot be handled, exiting\n"
    exit
  fi 
  #how many dup sites now?
  awk '{print $1,$4}' ${dirI}/temp_${out_name}_cleaning_6.bim | sort | uniq -c | awk '$1==2 {print $2,$3}' | wc -l
    #IBD 18141, WAMHS 18232, OAG2 2, OAG1 0 HEI_cas 0 HEI_con 0 SDH 24 610k_case 2 610k_con 2 OMNI_case 24 AMFS_case 24 AMFS_con 24 QTWIN 7029 EPIGENE 5994 
  echo ""
  #then write the straightforward dups
  awk 'NR==FNR {a[$1,$2];next} ($1,$4) in a {print $0}' <(awk '{print $1,$4}' ${dirI}/temp_${out_name}_cleaning_6.bim | sort | uniq -c | awk '$1==2 {print $2,$3}') <( sed 's/\t/ /g' ${dirI}/temp_${out_name}_cleaning_6.bim | sed 's/ \+/ /g') > ${dirI}/temp_${out_name}_cleaning_6.bim_dups; wc -l ${dirI}/temp_${out_name}_cleaning_6.bim_dups
    #IBD 36282, WAMHS 36464 (doubled). OAG2 4 OAG1 0 HEI_cas 0 HEI_con 0 SDH 48 610k_case 4 610k_con 4 OMNI_case 48 AMFS_case 48 AMFS_con 48 QTWIN 14058 EPIGENE 11998 MIA 0
  echo ""
  #pull out the double matches - reverse the match so we print the first one, (as in the first of the pair that matches) then exclude that and redo it.
  awk 'NR==FNR {a[$1,$4]=$0;next} ($1,$2) in a {print $1,$2,a[$1,$2]}' ${dirI}/temp_${out_name}_cleaning_6.bim_dups  <(awk '{print $1,$4}' ${dirI}/temp_${out_name}_cleaning_6.bim | sort | uniq -c | awk '$1==2 {print $2,$3}') > ${dirI}/temp_${out_name}_cleaning_6_overlap_1; wc -l ${dirI}/temp_${out_name}_cleaning_6_overlap_1
    #IBD 18141, WAMHS 18232 OAG2 2 OAG1 0 HEI_cas 0 HEI_con 0 SDH 24 610k_case 2 610k_con 2 OMNI_case 24 AMFS_case 24 AMFS_con 24 QTWIN 7029 EPIGENE 5994 MIA 0
  echo ""
  #now remove those from the dup set
  awk 'NR==FNR {a[$3,$4,$5,$6,$7,$8];next} !(($1,$2,$3,$4,$5,$6) in a)' ${dirI}/temp_${out_name}_cleaning_6_overlap_1 ${dirI}/temp_${out_name}_cleaning_6.bim_dups > ${dirI}/temp_${out_name}_cleaning_6.bim_dups_remaining; wc -l ${dirI}/temp_${out_name}_cleaning_6.bim_dups_remaining
    #IBD 18141 WAMHS, 18232 OAG2 2 OAG1 0 HEI_cas 0 HEI_con 0 SDH 24 610k_case 2 610k_con 2 OMNI_case 24 AMFS_case 24 AMFS_con 24 QTWIN 7029 EPIGENE 5994 MIA 0
  echo ""
  #and do it again - should now get the duplicates paired up
  awk 'NR==FNR {a[$1,$4]=$0;next} ($1,$2) in a {print $0,a[$1,$2]}' ${dirI}/temp_${out_name}_cleaning_6.bim_dups_remaining ${dirI}/temp_${out_name}_cleaning_6_overlap_1 > ${dirI}/temp_${out_name}_cleaning_6_overlap_2; wc -l ${dirI}/temp_${out_name}_cleaning_6_overlap_2
    #IBD 18141, WAMHS 18232 OAG2 2 OAG1 0 HEI_cas 0 HEI_con 0 SDH 24 610k_case 2 610k_con 2 OMNI_case 24 AMFS_case 24  AMFS_con 24 QTWIN 7029 EPIGENE 5994 MIA 0
  echo ""
  head -4 ${dirI}/temp_${out_name}_cleaning_6_overlap_2
  echo ""
  #now bring in the freq data but merging it in twice
  awk 'NR==FNR {a[$2]=$3" "$4" "$5" "$6} $10 in a {print $0,a[$10]}' <( sed 's/ \+/ /g' ${dirI}/temp_${out_name}_cleaning_6.frq)  <( awk 'NR==FNR {a[$2]=$3" "$4" "$5" "$6} $4 in a {print $0,a[$4]}' <( sed 's/ \+/ /g' ${dirI}/temp_${out_name}_cleaning_6.frq) ${dirI}/temp_${out_name}_cleaning_6_overlap_2 ) > ${dirI}/temp_${out_name}_cleaning_6_overlap_3; wc -l ${dirI}/temp_${out_name}_cleaning_6_overlap_3
   #IBD 18141, WAMHS 18232 OAG2 2 OAG1 0 HEI_cas 0 HEI_con 0 SDH 24 610k_case 2 610k_con 2 OMNI_case 24 AMFS_case 24 AMFS_con 24 QTWIN 7029 MIA 0
  echo ""
  head -4 ${dirI}/temp_${out_name}_cleaning_6_overlap_3
    #WAMHS - first and third looks to the same SNP but fipped, the second the fourth the same SNP and strand. IBD looks sim but obv diff freq etc
    #10 100017453 10 rs1983864 0 100017453 G T 10 exm847519 0 100017453 C A G T 0.3425 2546 C A 0.3425 2546
    #10 100017553 10 rs11189525 0 100017553 A C 10 exm847525 0 100017553 A C A C 0 2546 A C 0 2546
    #10 100144782 10 rs2296441 0 100144782 T C 10 exm847614 0 100144782 A G T C 0.2887 2546 A G 0.2887 2546
    #10 100148176 10 rs2147896 0 100148176 G A 10 exm847632 0 100148176 G A G A 0.34 2544 G A 0.34 2544

  #need to identify the flips, and flip them, so my existing function can then manage comparing them... actually that is silly, don't need to flip them, just work out if they are the same SNP (strand regardless) and drop the extranous one. First drop NA or 0 freq SNPs, then check for same SNP same strand, then same SNP other strang. Could do each in a single command but this way is cleaner to understand/write. Did a quick test where I wrote out both SNPs at each check ($4,$10) and looked for duplicates (as in the same pair was tripping more than one check) and it only happened where both SNPs were monomorphic. Write out the one with the the smallest N of chromosomes. 
  #For WAMHS these setting left 34 SNPs that couldn't be handled. Some were just a bit more than 0.01 diff but 0.02 diff still 22 so not worth changing. The remaining were SNPs near 0.5 MAF and thus flipped, or properly different SNPs. IBD gets 69; allowing 2% diff leaves 38. So set to 2%, but not worth messing about with any more so just drop these ones that don't match. 
  awk '{
    if($17==0 || $17=="NA")
      {print $4 > "temp_overlapping_NA_exclude_1"};
    if($21==0 || $21=="NA")
      {print $10 > "temp_overlapping_NA_exclude_2"};
    if($17!="NA" && $17!=0 && $21!="NA" && $21!=0 && (($17 - $21) > (-0.02) ) && (($17 - $21) < 0.02) && $18<=$22 && $13==$15 && $14==$16)
      {print $4 > "temp_overlapping_same_SNP_1"};
    if($17!="NA" && $17!=0 && $21!="NA" && $21!=0 && (($17 - $21) > (-0.02) ) && (($17 - $21) < 0.02) && $18>$22 && $13==$15 && $14==$16)
      {print $10 > "temp_overlapping_same_SNP_2"};
    if($17!="NA" && $17!=0 && $21!="NA" && $21!=0 && (($17 - $21) > (-0.02) ) && (($17 - $21) < 0.02) && $18<=$22 && (  ( $13=="A" && $14=="T" && $15=="T" && $16=="A" ) ||  ( $13=="T" && $14=="A" && $15=="A" && $16=="T" ) ||  ( $13=="C" && $14=="G" && $15=="G" && $16=="C" ) ||  ( $13=="G" && $14=="C" && $15=="C" && $16=="G" ) ||  ( $13=="A" && $14=="C" && $15=="T" && $16=="G" ) ||  ( $13=="A" && $14=="G" && $15=="T" && $16=="C" ) ||  ( $13=="C" && $14=="A" && $15=="G" && $16=="T" ) ||  ( $13=="C" && $14=="T" && $15=="G" && $16=="A" ) ||  ( $13=="G" && $14=="A" && $15=="C" && $16=="T" ) ||  ( $13=="G" && $14=="T" && $15=="C" && $16=="A" ) ||  ( $13=="T" && $14=="G" && $15=="A" && $16=="C" ) ||  ( $13=="T" && $14=="C" && $15=="A" && $16=="G" ) ) )
      {print $4 > "temp_overlapping_same_SNP_flipped_1"};
    if($17!="NA" && $17!=0 && $21!="NA" && $21!=0 && (($17 - $21) > (-0.02) ) && (($17 - $21) < 0.02) && $18>$22 && (  ( $13=="A" && $14=="T" && $15=="T" && $16=="A" ) ||  ( $13=="T" && $14=="A" && $15=="A" && $16=="T" ) ||  ( $13=="C" && $14=="G" && $15=="G" && $16=="C" ) ||  ( $13=="G" && $14=="C" && $15=="C" && $16=="G" ) ||  ( $13=="A" && $14=="C" && $15=="T" && $16=="G" ) ||  ( $13=="A" && $14=="G" && $15=="T" && $16=="C" ) ||  ( $13=="C" && $14=="A" && $15=="G" && $16=="T" ) ||  ( $13=="C" && $14=="T" && $15=="G" && $16=="A" ) ||  ( $13=="G" && $14=="A" && $15=="C" && $16=="T" ) ||  ( $13=="G" && $14=="T" && $15=="C" && $16=="A" ) ||  ( $13=="T" && $14=="G" && $15=="A" && $16=="C" ) ||  ( $13=="T" && $14=="C" && $15=="A" && $16=="G" ) ) )
      {print $10 > "temp_overlapping_same_SNP_flipped_2"};
   if($17!="NA" && $17!=0 && $21!="NA" && $21!=0 && $13==$16 && $14==$15 && ($17>0.48 && $21>0.48) && $18<=$22 )
      {print $4 > "temp_overlapping_same_SNP_MAF_50_swapped_A1_1" };
   if($17!="NA" && $17!=0 && $21!="NA" && $21!=0 && $13==$16 && $14==$15 && ($17>0.48 && $21>0.48) && $18>$22 )
      {print $10 > "temp_overlapping_same_SNP_MAF_50_swapped_A1_2" };
  }' ${dirI}/temp_${out_name}_cleaning_6_overlap_3
  echo ""
  wc -l temp_overlapping_NA_exclude_1
  echo ""
   #IBD 1198, WAMHS 1360 OAG2 0 OAG1 0 HEI_cas 0 HEI_con 0 SDH 0 610k_case 0 610k_con 0 OMNI_case 0 AMFS_case 0 AMFS_con 0 QTWIN 831 EPIGENE 615 MIA 0
  wc -l temp_overlapping_NA_exclude_2
  echo ""
   #IBD 1194, WAMHS 1362 OAG2 0 OAG1 0 HEI_cas 0 HEI_con 0 SDH 0 610k_case 0 610k_con 0 OMNI_case 0 AMFS_case 0 AMFS_con 0 QTWIN 654 EPIGENE 621 MIA 0
  wc -l temp_overlapping_same_SNP_1
  echo ""
   #IBD 13796, WAMHS 7212 OAG2 0 OAG1 0 HEI_cas 0 HEI_con 0 SDH 8 610k_case 2 610k_con 2 OMNI_case 18 AMFS_case 20 AMFS_con 18 QTWIN 5425 EPIGENE 5190 MIA 0
  wc -l temp_overlapping_same_SNP_2
  echo ""
   #IBD 3091, WAMHS 2195 OAG2 0 OAG1 0 HEI_cas 0 HEI_con 0 SDH 2 610k_case 0 610k_con 0 OMNI_case 5 AMFS_case 3 AMFS_con 5 QTWIN 491 EPIGENE 152 MIA 0
  wc -l temp_overlapping_same_SNP_flipped_1
  echo ""
   #IBD 0, WAMHS 5687 OAG2 0 OAG1 0 HEI_cas 0 HEI_con 0 SDH 9 610k_case 0 610k_con 0 OMNI_case 0 AMFS_case 1 AMFS_con 1 QTWIN 0 EPIGENE 0 MIA 0
  wc -l temp_overlapping_same_SNP_flipped_2
  echo ""
   #IBD 0, WAMHS 1722 OAG2 2 OAG1 0 HEI_cas 0 HEI_con 0 SDH 5 610k_case 0 610k_con 0 OMNI_case 1 AMFS_case 0 AMFS_con 0 QTWIN 0 EPIGENE 0 MIA 0
  wc -l temp_overlapping_same_SNP_MAF_50_swapped_A1_1
  echo ""
   #IBD 4, WAMHS 2 OAG2 0 OAG1 0 HEI_cas 0 HEI_con 0 SDH 0 610k_case 0 610k_con 0 OMNI_case 0 AMFS_case 0 AMFS_con 0 QTWIN 1 EPIGENE 0 MIA 0
  wc -l temp_overlapping_same_SNP_MAF_50_swapped_A1_2
  echo ""
   #IBD 3, WAMHS 1 OAG2 0 OAG1 0 HEI_cas 0 HEI_con 0 SDH 0 610k_case 0 610k_con 0 OMNI_case 0 AMFS_case 0 AMFS_con 0 QTWIN 2 EPIGENE 0 MIA 0
  #what is left?
  awk 'NR==FNR {a[$1];next} !($10 in a)' <(cat temp_overlapping_NA_exclude_2 temp_overlapping_same_SNP_2 temp_overlapping_same_SNP_flipped_2 temp_overlapping_same_SNP_MAF_50_swapped_A1_2  2> /dev/null )  <(awk 'NR==FNR {a[$1];next} !($4 in a)' <(cat temp_overlapping_NA_exclude_1 temp_overlapping_same_SNP_1 temp_overlapping_same_SNP_flipped_1 temp_overlapping_same_SNP_MAF_50_swapped_A1_1  2> /dev/null ) ${dirI}/temp_${out_name}_cleaning_6_overlap_3) | wc -l
  echo ""
     #IBD 34 out of 18k, WAMHS 23 out of 18k OAG2 0 OAG1 0 HEI_cas 0 HEI_con 0 SDH 0 610k_case 0 610k_con 0 OMNI_case 0 AMFS_case 0 AMFS_con 0 QTWIN 223 EPIGENE 12
  #drop these
  awk 'NR==FNR {a[$1];next} !($10 in a)' <(cat temp_overlapping_NA_exclude_2 temp_overlapping_same_SNP_2 temp_overlapping_same_SNP_flipped_2 temp_overlapping_same_SNP_MAF_50_swapped_A1_2  2> /dev/null ) <(awk 'NR==FNR {a[$1];next} !($4 in a)' <(cat temp_overlapping_NA_exclude_1 temp_overlapping_same_SNP_1 temp_overlapping_same_SNP_flipped_1 temp_overlapping_same_SNP_MAF_50_swapped_A1_1  2> /dev/null ) ${dirI}/temp_${out_name}_cleaning_6_overlap_3) | awk '{print $4}' > temp_unpairable_1
  awk 'NR==FNR {a[$1];next} !($10 in a)' <(cat temp_overlapping_NA_exclude_2 temp_overlapping_same_SNP_2 temp_overlapping_same_SNP_flipped_2 temp_overlapping_same_SNP_MAF_50_swapped_A1_2  2> /dev/null )  <(awk 'NR==FNR {a[$1];next} !($4 in a)' <(cat temp_overlapping_NA_exclude_1 temp_overlapping_same_SNP_1 temp_overlapping_same_SNP_flipped_1 temp_overlapping_same_SNP_MAF_50_swapped_A1_1  2> /dev/null ) ${dirI}/temp_${out_name}_cleaning_6_overlap_3) | awk '{print $10}' > temp_unpairable_2
  cat temp_overlapping_NA_exclude_1 temp_overlapping_NA_exclude_2 temp_overlapping_same_SNP_1 temp_overlapping_same_SNP_2 temp_overlapping_same_SNP_flipped_1 temp_overlapping_same_SNP_flipped_2 temp_overlapping_same_SNP_MAF_50_swapped_A1_1 temp_overlapping_same_SNP_MAF_50_swapped_A1_2 temp_unpairable_1 temp_unpairable_2 2> /dev/null | sort | uniq -d | wc -l
    #should be 0; IBD, WAMHS 0 OAG2 0 OAG1 0 HEI_cas 0 HEI_con 0 SDH 0 610k_case 0 610k_con 0 OMNI_case 0 AMFS_case 0 AMFS_con 0 QTWIN 0 EPIGENE 0
  cat temp_overlapping_NA_exclude_1 temp_overlapping_NA_exclude_2 temp_overlapping_same_SNP_1 temp_overlapping_same_SNP_2 temp_overlapping_same_SNP_flipped_1 temp_overlapping_same_SNP_flipped_2 temp_overlapping_same_SNP_MAF_50_swapped_A1_1 temp_overlapping_same_SNP_MAF_50_swapped_A1_2 temp_unpairable_1 temp_unpairable_2 2> /dev/null > temp_overlapping_exclude
  echo -e "\nRemoving the `sed -n '$=' temp_overlapping_exclude` overlapping and redundant SNPs\n" 
     #IBD 19354 WAMHS 19587 OAG2 2 OAG1 0 HEI_cas 0 HEI_con 0 610k_case 2 omni_case 24 AMFS_case 24 QTWIN 7850  EPIGENE 6602 MIA 0
  plink_1.90 --threads 1 --bfile ${dirI}/temp_${out_name}_cleaning_6 --exclude temp_overlapping_exclude --make-bed --out  ${dirI}/temp_${out_name}_cleaning_7
     #IBD 933559 variants loaded --exclude: 914205 variants remaining.
     #WAMHS 936505 variants loaded --exclude: 916918 variants remaining
     #OAG2 707537 variants loaded --exclude: 707535 variants remaining.
     #OAG1 856283 variants loaded --exclude: 856283 variants remaining.
     #HEI_cas 218704 variants loaded from .bim file --exclude: 218704 variants remaining.
     #HEI_con 221018 variants loaded from .bim file. --exclude: 221018 variants remaining
     #SDH 940334 variants loaded from .bim file --exclude: 940310 variants remaining.
     #610k_case 548565 variants loaded  --exclude: 548563 variants remaining.
     #610k_control 533169 variants loaded --exclude: 533167 variants remaining.
     #omni_cases 803876 variants loaded --exclude: 803852 variants remaining.
     #AMFS_cases 803868 variants loaded from .bim file. --exclude: 803844 variants remaining.
     #AMFS_con 803845 variants loaded --exclude: 803821 variants remaining 
     #QTWIN 459153 variants loaded --exclude: 451303 variants remaining.
     #EPIGENE 406754 variants loaded  --exclude: 400152 variants remaining.
     #MIA 488143 variants loaded --exclude: 488143 variants remaining.

  #clean up
  rm -f ${dirI}/temp_overlapping_NA_exclude_1 ${dirI}/temp_overlapping_NA_exclude_2 ${dirI}/temp_overlapping_same_SNP_1 ${dirI}/temp_overlapping_same_SNP_2 ${dirI}/temp_overlapping_same_SNP_flipped_1 ${dirI}/temp_overlapping_same_SNP_flipped_2 ${dirI}/temp_overlapping_same_SNP_MAF_50_swapped_A1_1 ${dirI}/temp_overlapping_same_SNP_MAF_50_swapped_A1_2 ${dirI}/temp_unpairable_1 ${dirI}/temp_unpairable_2 ${dirI}/temp_overlapping_exclude

    #18/07/2016 found the merge equal pos bug when looking at the LD of rs9402684 and rs7745098 following flip-scan - CT/TC alleles should be in LD, but following this step CC and TT became in LD, as merge-equal-pos overwrite the ALLLELES when merging!
      #plink_1.90 --threads 1 --bfile ${dirI}/temp_${out_name}_cleaning_7 --ld rs9402684 rs7745098
      #  Haplotype     Frequency    Expectation under LE
      #   ---------     ---------    --------------------
      #          CT      0.487027                0.237459
      #          TT      0.000541                0.250109
      #          CC      0                       0.249568
      #          TC      0.512432                0.262864
  ##OLD COMMAND - merge-equal-pos does not work as intended
  ##plink_1.90 --threads 1 --bfile ${dirI}/temp_${out_name}_cleaning_6 --bmerge ${dirI}/temp_${out_name}_cleaning_6 --merge-equal-pos --make-bed --out ${dirI}/temp_${out_name}_cleaning_7
    #SDH 944340 markers to be merged --25 same positions warnings (3 + 22) --940315 variants and 570 people pass filters and QC.
    #610k_cases 1018549 markers loaded; 1553 (1550 + 3) same positions warnings; 1016996 variants and 925 people pass filters and QC.
    #610k_controls - no extra errors. 1018323 markers loaded. 1550 (1547 + 3) same-position warnings 1016773 variants and 3954 people pass filters and Q
    #QMEGA_omni_cases - no extra errors. 1551 (3 + 1548) more same-position warnings.  1025732 variants and 694 people pass filters and QC.
    #AMFS_cases: - no extra errors. 1551 (3 + 1548) more same-position warnings.  1025724 variants and 548 people pass filters and QC
    #AMFS_controls:  no extra errors. 1551 (3 + 1548) more same-position warnings. 1025700 variants and 430 people pass filters and QC

  echo -e "\n+++++++Now pairing SNPs with 1KG based on CHR and POS and checking for ID mismatches using awk_check_based_on_pos_1KG function\n"
  plink_1.90 --threads 1 --bfile ${dirI}/temp_${out_name}_cleaning_7 --freq --out ${dirI}/temp_${out_name}_cleaning_7 --silent
  #13/07/2016 this is the older version of this command for posterity. It works EXCEPT where there is a INDEL overlapping the SNP, and is listed first in the 1000 genomes panel before the matching SNP. As each SNP is in PLINK once, only the first match in 1kg is shown, so only the SNP vs INDEL is paired. Reversing this would print both pairings, allowing the SNP to be recovered.
  ##awk 'NR==FNR {a[$1,$3]=$1" "$3" "$4" "$5" "$9" "$13" "$2;next} ($1,$7) in a {print $0,a[$1,$7]}' <(gunzip -dc ${dirI}/ALL_1000G_phase1integrated_v3_ALL_chr_legends_with_chromosome.txt.gz) <( awk 'NR==FNR {a[$2]=$4;next} $2 in a {print $0,a[$2]}' <(sed 's/ \+/ /g' ${dirI}/temp_${out_name}_cleaning_7.bim | sed 's/^ //g')  <(sed 's/ \+/ /g' ${dirI}/temp_${out_name}_cleaning_7.frq | sed 's/^ //g' | sed 's/ $//g') ) > ${dirI}/temp_${out_name}_cleaning_7_1KG_merge_pos; wc -l ${dirI}/temp_${out_name}_cleaning_7_1KG_merge_pos
  ##25/07/2017 Found the reverse - with exm chip data there are multiple overlapping rare variants (e.g. exm1771279 and exm2277033) that 'block' correct matching. exm2277033 looks to be either super rare or a nonsense SNP, whereas exm1771279 is a real 1kg SNP, but exm2277033 is in the bim file first and thus used in the matching with 1KG. But actually his is okay. Such matches are picked up as temp_non_biallelic_snps_ID, and removed, and detected again later in a final test for overlapping SNPs; I just need to add another renaming step there to finally fix these
  awk 'NR==FNR  {a[$1,$7]=$0;next} ($1,$3) in a {print a[$1,$3],$1,$3,$4,$5,$9,$13,$2}' <( awk 'NR==FNR {a[$2]=$4;next} $2 in a {print $0,a[$2]}'  <(sed 's/ \+/ /g' ${dirI}/temp_${out_name}_cleaning_7.bim | sed 's/^ //g')  <(sed 's/ \+/ /g' ${dirI}/temp_${out_name}_cleaning_7.frq | sed 's/^ //g' | sed 's/ $//g') ) <(gunzip -dc ${dirI}/ALL_1000G_phase1integrated_v3_ALL_chr_legends_with_chromosome.txt.gz) > ${dirI}/temp_${out_name}_cleaning_7_1KG_merge_pos; wc -l ${dirI}/temp_${out_name}_cleaning_7_1KG_merge_pos
    #SDH 910865 610k_cas 539245 610k_con 524168 Omni_cas 790304 AMFS_cas 790296 AMFS_con 790270 WAMHS 793779 HEI_cas 211004. HEI_con: 213254. OAG1 834462 OAG2 684914 IBD 791715 QTWIN 372933 EPIGENE 331688 MIA 466762
  echo ""
  wc -l ${dirI}/temp_${out_name}_cleaning_7.bim
    #SDH 940310 610k 548563 610k_con 533167 Omni_cas 803852 AMFS_cas 803844 AMFS_con 803821. WAMHS 916918. HEI_cas 218704 HEI_con 221018 OAG1 856283 OAG2 707535 IBD 914205 QTWIN 451303 EPIGENE 400152 MIA 488143

  echo -e "\nFor ${out_name} there are `  awk '$2!=$14' ${dirI}/temp_${out_name}_cleaning_7_1KG_merge_pos | wc -l` SNPs with mismatching IDs but the same position. This will include true ID differences and overlapping indels/SNPs\n"
    #SDH 59026 610k_cas 4658 610k_con 4150 Omni_cases 39197 AMFS_cas 39198 AMFS_con 39198 WAMHS 143046 (due to exm, and loss of non exm version if it has less info) HEI_cas 856 HEI_con 887 OAG1 11976 OAG2 2911 IBD 141699 (likewise exm) QTWIN 131501 EPIGENE 97347 MIA 41041
  #see 0.4 for this function
  awk -v MAF=${MAF} -v DIFF=${DIFF} -f /working/lab_stuartma/scratch/matthewL/melanoma_meta_analysis/SCRIPTS/awk_check_based_on_pos_1KG  ${dirI}/temp_${out_name}_cleaning_7_1KG_merge_pos

  #13/07/2016 Was not explcitely filtering out SNPs I couldn't align/map to 1kG as I did for the rsID match. This causes some downstream problems (ambigious SNPs that were not aligned, and thus not on the same strand, but retained in the dataset. So explicitely exclude them.
  awk 'NR==FNR {a[$1];next} !($2 in a) ' <(cat temp_aligned_ID_concordant temp_aligned_ID_disc temp_flipped_ID_concordant temp_flipped_ID_disc temp_ambigious_ID_MAF_to_high_too_assess temp_ambigious_aligned_ID_concordant_same_allele temp_ambigious_aligned_ID_disc_same_allele temp_ambigious_aligned_ID_concordant_other_allele temp_ambigious_aligned_ID_disc_other_allele temp_ambigious_flipped_ID_concordant_same_allele temp_ambigious_flipped_ID_disc_same_allele temp_ambigious_flipped_ID_concordant_other_allele temp_ambigious_flipped_ID_disc_other_allele temp_non_biallelic_snps_ID temp_indels_ID_aligned_same_allele temp_indels_ID_aligned_other_allele 2> /dev/null ) ${dirI}/temp_${out_name}_cleaning_7_1KG_merge_pos > ${dirI}/temp_${out_name}_unmappable_snps; wc -l ${dirI}/temp_${out_name}_unmappable_snps
    #SDH 127, 610k_cases 18, 610k_controls 15, QMEGA_omni_cases and AMFS_cases 109, AMFS controls 100, WAMHS 95. HEI_cases and controls: 10. OAG1 55 OAG2 59 IBD 73 QTWIN 165 EPIGENE 73 MIA 397
  echo ""
  wc -l temp_aligned_ID_concordant
  echo ""
    #SDH 838775 610k_cas 532932 610k_con 518583 Omni_cases 740812 AMFS_cas 740811 AMFS_con 740810 WAMHS 649432 HEI_cas 210148 HEI_con 212367 OAG1 810387 OAG2 680472 IBD 648735 QTWIN 241051 EPIGENE 233989 MIA 405229
  wc -l temp_aligned_ID_disc
  echo ""
    #SDH 25875; 610k_cas 4052 610k_con 3761 Omni_cases 31586 AMFS_cases 31588 AMFS_controls 31589, WAMHS 61865 HEI_ca: 388 HEI_con 400 OAG1 5219 OAG2 1265 IBD 124411 QTWIN 118987 EPIGENE 86753 MIA 15771
  wc -l temp_flipped_ID_concordant
  echo ""
    #SDH, 610k, 610k_controls, QMEGA_omni_cases, AMFS_cases, AMFS_control, WAMHS, HEI_cases and controls: 0 OAG1 0 OAG2 0 IBD 0 QTWIN 0 EPIGENE 0 MIA 0
  wc -l temp_flipped_ID_disc
  echo ""
    #SDH 25650. 610k_cases 151 610k_controls 4, Omni_cases 2455 AMFS_cases and AMFS_con 2454 WAMHS 61865. HEI_cases: 398 HEI_con 412 OAG1 5037 OAG2 1324 IBD 46 QTWIN 0 EPIGENE 0 MIA 15338
  wc -l temp_ambigious_ID_MAF_to_high_too_assess
  echo ""
    #SDH 938. 610k_cases 40 610k controls 30 Omni_cases: 997. AMFS_cases 1021. AMFS_control: 1022, WAMHS 653, HEI_cases and con: 0 OAG1 270 OAG2 0 IBD 577 QTWIN 603  EPIGENE 497 MIA 260
  wc -l temp_ambigious_aligned_ID_concordant_same_allele
  echo ""
    #SDH 10435; 610k - 1272. 610k controls 1104. QMEGA_omni_cases: 7939. AMFS_cases 7941. AMFS_control 7919, WAMHS 1050, HEI_cases and con: 0 OAG1 9541 OAG2 1235 IBD 1029 QTWIN 288 EPIGENE 262 MIA 18328
  wc -l temp_ambigious_aligned_ID_disc_same_allele
  echo ""
    #SDH 2681; 610k_cases 128, 610k controls 99. QMEGA_omni_cases: 2869, AMFS_cases 2865. AMFS_control 2858, WAMHS 7795, HEI_cases and con: 0 OAG1 467 OAG2 0 IBD 15242 QTWIN 10675  EPIGENE 9130 MIA 1035
  wc -l temp_ambigious_aligned_ID_concordant_other_allele
  echo ""
    #SDH 2629; 610k 383, 610k control 331. QMEGA_omni_cases: 2356, AMFS_cases 2346. AMFS_control 2343, WAMHS 251, HEI_cases and con: 0 OAG1 2558 OAG2 296 IBD 252 QTWIN 93  EPIGENE 0 90 MIA 2164
  wc -l temp_ambigious_aligned_ID_disc_other_allele
  echo ""
    #SDH 313; 610k_cases 16, 610k controls 17. QMEGA_omni_cases: 492, AMFS_cases 472. AMFS_control 481, WAMHS 209, HEI_cases and con: 0 OAG1 79 OAG2 2 IBD 464 QTWIN 455  EPIGENE 385 MIA 79
  wc -l temp_ambigious_flipped_ID_concordant_same_allele
  echo ""
    #SDH, 610k, 610k controls, QMEGA_omni_cases none, WAMHS none, HEI_cases and con: 0 OAG1 0 OAG2 0 IBD 0 QTWIN 0  EPIGENE 0 MIA 0
  wc -l temp_ambigious_flipped_ID_disc_same_allele
  echo ""
    #SDH 2732; 610k 58, 610k controls 49. QMEGA_omni_cases: 263, AMFS_cases 268. AMFS_control 271. WAMHS 7784, HEI_cases and con: 0 OAG1 476 OAG2 1 IBD 3 QTWIN 0 EPIGENE 0 MIA 964
  wc -l temp_ambigious_flipped_ID_concordant_other_allele
  echo ""
    #SHD, 610k, 610k_controls, QMEGA_omni_cases, WAMHS none, HEI_cases and con: 0 OAG1 0 OAG2 0 IBD 0 QTWIN 0 EPIGENE 0 MIA 0
  wc -l temp_ambigious_flipped_ID_disc_other_allele
  echo ""
    #SDH 296; 610k 12, 610k controls 11. QMEGA_omni_cases: 88, AMFS_cases 83. AMFS_control 85, WAMHS 266, HEI_cases and con: 0 OAG1 78 OAG2 0 IBD 2 QTWIN 0 EPIGENE 0 MIA 85
  wc -l temp_non_biallelic_snps_ID
  echo ""
    #SDH 113; 610k, 610k controls 6, QMEGA_omni_cases, AMFS_cases and AMFS_con 36, WAMHS 671 (exm), HEI_cases and con: 0 OAG1 17 OAG2 32 IBD 657  QTWIN 511 EPIGENE 429 MIA 76

  if [[ -a temp_indels_ID_aligned_same_allele ]]
  then
    echo -e "Indels with non NA freq found, attempting to replace I/D alleles with 1kG format alleles\n"
    wc -l temp_indels_ID_aligned_same_allele
    echo ""
      #SDH - none - checked, all SNPs are alleles or 0/NA. Nothing diff. 610k - 27 but as you can see below are not actually on the 610k array, all NA. will be on omni samles though. 610k_controls same. QMEGA_omni_cases, AMFS_cases, AMFS_controls: 27, WAMHS and IBD 1
    #13/07/2016 now that I have changed the order of pairing PLINK and 1KG indels be matched twice, so make the match more explicit by checking for both SNP ids, not just the first one - not just a[$1], $2 in a;
     #08/09/2016 added a head step - MIA has ~ 6k and all are printed...
    awk 'NR==FNR {a[$1,$2];next} ($2,$14) in a' temp_indels_ID_aligned_same_allele ${dirI}/temp_${out_name}_cleaning_7_1KG_merge_pos | head 
      #for QMEGA_cases_omni all look correct. As lon as the freq is non NA make a list to update alleles; example 2 below. Remember order for alleles swapped in 1kg
      #1 rs4647009 I D 0.05476 1388 59249359 1 59249359 C CTG 0.0449 0.0449 chr1:59249359:I
      #1 rs3783480 I D 0.00938 1386 68153678 1 68153678 A AG 0.0172 0.0172 chr1:68153678:I
      #WAMHS      #16 exm-rs2066847 I D 0.01885 2546 50763778 16 50763778 G GC 0.0106 0.0106 chr16:50763778:I
      #IBD       #16 exm-rs2066847 I D 0.01453 1858 50763778 16 50763778 G GC 0.0106 0.0106 chr16:50763778:I
      #QTWIN  19 rs28371511 D I 0.2095 1924 15771028 19 15771028 TG T 0.1451 0.1451 chr19:15771028:D EPI has same SNP, diff freq obv but v close
    echo ""
    #old command; works for everything by oncoarray which has some borked indel/SNP mix ups    
    ##awk 'NR==FNR {a[$1,$2];next} ($2,$14) in a && $5!="NA" {print $2,$3,$4,$11,$10}' temp_indels_ID_aligned_same_allele ${dirI}/temp_${out_name}_cleaning_7_1KG_merge_pos > ${dirI}/temp_${out_name}_cleaning_7_update_indels; wc -l ${dirI}/temp_${out_name}_cleaning_7_update_indels

     rm -f temp_cleaning_7_update_indel_alleles_1 temp_cleaning_7_indels_diff_exclude_1
     #remember that this includes cases where the same indel is paired to an indel and an overlapping SNP, so pre filter the ${dirI}/temp_${out_name}_cleaning_7_1KG_merge_pos for only indel to indel matches (e.g. cols 10 or 11 more than one byte)
     awk 'NR==FNR {a[$1,$2];next} {
       if(($2,$14) in a && $5!="NA" && (( $5 - $12 ) > (-'''${DIFF}''')) && (( $5 - $12 ) < '''${DIFF}''' ))
         {print  $2,$3,$4,$11,$10 > "temp_cleaning_7_update_indel_alleles_1"}
       if(($2,$14) in a && $5!="NA" && !((( $5 - $12 ) > (-'''${DIFF}''')) && (( $5 - $12 ) < '''${DIFF}''' )))
         {print  $2 > "temp_cleaning_7_indels_diff_exclude_1"}
       }' temp_indels_ID_aligned_same_allele <( awk 'length($10)>1 || length($11)>1' ${dirI}/temp_${out_name}_cleaning_7_1KG_merge_pos)
     wc -l temp_indels_ID_aligned_same_allele; echo ""
       #MIA 5997
     wc -l temp_cleaning_7_update_indel_alleles_1; echo ""
       #MIA 5598
     wc -l temp_cleaning_7_indels_diff_exclude_1; echo ""
       #MIA 399
     #any double matches?
     cat <(cut -d" " -f1 temp_cleaning_7_update_indel_alleles_1) temp_cleaning_7_indels_diff_exclude_1 | sort | uniq -d
      #MIA none
    echo ""
  fi #check for indels and if so examine in the dataset

  if [[ -a temp_indels_ID_aligned_other_allele ]]
  then
    echo -e "Indels with non NA freq testing other allele found, attempting to replace I/D alleles with 1kG format alleles\n"
    wc -l temp_indels_ID_aligned_other_allele
      #QMEGA_omni_cases, AMFS cases, AMFS_controls: 5 QTWIN 1
    echo ""
    #13/07/2016 now that I have changed the order of pairing PLINK and 1KG indels be matched twice, so make the match more explicit by checking for both SNP ids, not just the first one
    ##awk 'NR==FNR {a[$1,$2];next} ($2,$14) in a && $5!="NA" {print $2,$3,$4,$10,$11}' temp_indels_ID_aligned_other_allele ${dirI}/temp_${out_name}_cleaning_7_1KG_merge_pos >> ${dirI}/temp_${out_name}_cleaning_7_update_indels; wc -l ${dirI}/temp_${out_name}_cleaning_7_update_indels
    awk 'NR==FNR {a[$1,$2];next} {
      if(($2,$14) in a && $5!="NA" && (( $5 - (1-$12) ) > (-'''${DIFF}''')) && (( $5 - (1-$12) ) < '''${DIFF}''' ))
        {print  $2,$3,$4,$11,$10 > "temp_cleaning_7_update_indel_alleles_2"}
      if(($2,$14) in a && $5!="NA" && !((( $5 - (1-$12) ) > (-'''${DIFF}''')) && (( $5 - (1-$12) ) < '''${DIFF}''' )))
        {print  $2 > "temp_cleaning_7_indels_diff_exclude_2"}
      }' temp_indels_ID_aligned_other_allele <( awk 'length($10)>1 || length($11)>1' ${dirI}/temp_${out_name}_cleaning_7_1KG_merge_pos)
    wc -l temp_indels_ID_aligned_other_allele; echo ""
      #MIA 766
    wc -l temp_cleaning_7_update_indel_alleles_2; echo ""
      #MIA 710
    wc -l temp_cleaning_7_indels_diff_exclude_2; echo ""
      #MIA 56
    #any double matches?
    cat <(cut -d" " -f1 temp_cleaning_7_update_indel_alleles_2) temp_cleaning_7_indels_diff_exclude_2 | sort | uniq -d
      #MIA 0
  fi
  cat temp_cleaning_7_update_indel_alleles_1 temp_cleaning_7_update_indel_alleles_2 2> /dev/null  > ${dirI}/temp_${out_name}_cleaning_7_update_indels
  cat temp_cleaning_7_indels_diff_exclude_1 temp_cleaning_7_indels_diff_exclude_2 2> /dev/null > ${dirI}/temp_${out_name}_cleaning_7_exclude_indels

  #check for the update file and if it is non zero (if above steps found no indels will be size 0)
  if [[ -s ${dirI}/temp_${out_name}_cleaning_7_update_indels || -s ${dirI}/temp_${out_name}_cleaning_7_exclude_indels  ]]
  then
    echo -e "Found a file to update alleles for indels that could be aligned to 1kG. Using...\n"
    plink_1.90 --threads 1 --bfile ${dirI}/temp_${out_name}_cleaning_7 --update-alleles ${dirI}/temp_${out_name}_cleaning_7_update_indels --make-bed --out ${dirI}/temp_${out_name}_cleaning_7_temp --exclude ${dirI}/temp_${out_name}_cleaning_7_exclude_indels
      #QMEGA_omni: --update-alleles: 32 variants updated. 803856 variants and 694 people pass filters and QC.
      #AMFS_cases: --update-alleles: 32 variants updated. 803844 variants and 548 people pass filters and QC
      #AMFS_controls: --update-alleles: 32 variants updated. 803821 variants and 430 people pass filters and QC.
      #WAMHS --update-alleles: 1 variant updated. 926098 variants and 1273 people pass filters and QC.
      #IBD --update-alleles: 1 variant updated. 916918 variants and 929 people pass filters and QC 
      #EPIGENE --update-alleles: 1, QTWIN --update-alleles: 1
      #MIA --update-alleles: 6308 variants updated. --exclude: 487688 variants remaining. 487688 variants and 1968 people pass filters and QC.
    echo ""
    plink_1.90 --threads 1 --bfile ${dirI}/temp_${out_name}_cleaning_7_temp  --make-bed --out ${dirI}/temp_${out_name}_cleaning_7
    echo -e "\nFinished fixing indel alleles from I/D to match 1KG\n"
  fi

  rm -f ${dirI}/temp_cleaning_7_update_indel_alleles_* ${dirI}/temp_cleaning_7_indels_diff_exclude_* ${dirI}/temp_${out_name}_cleaning_7_update_indels ${dirI}/temp_${out_name}_cleaning_7_exclude_indels
  
  #check line counts
  echo -e "For ${out_name} the function awk_check_based_on_pos_1KG reports `cat temp_aligned_ID_concordant temp_aligned_ID_disc temp_flipped_ID_concordant temp_flipped_ID_disc temp_ambigious_ID_MAF_to_high_too_assess temp_ambigious_aligned_ID_concordant_same_allele temp_ambigious_aligned_ID_disc_same_allele temp_ambigious_aligned_ID_concordant_other_allele temp_ambigious_aligned_ID_disc_other_allele temp_ambigious_flipped_ID_concordant_same_allele temp_ambigious_flipped_ID_disc_same_allele temp_ambigious_flipped_ID_concordant_other_allele temp_ambigious_flipped_ID_disc_other_allele temp_non_biallelic_snps_ID temp_indels_ID_aligned_same_allele temp_indels_ID_aligned_other_allele 2> /dev/null  | wc -l` SNPs"
    #SDH 910437 - some double matches from overlapping indels/SNPs. 610k_cas 539050 610k_con 523993 Omni_cases 789925 AMFS_cas 789917 AMFS_con 789900 WAMHS 793459 HEI_cases: 210934 HEI_con 213179 OAG1 834129 OAG2 684627 IBD 791419 QTWIN 372664 EPIGENE 331536 MIA 466092

  #how many rsID pairs are dups (that is more of an issue)
  var_dup_2=`cat temp_aligned_ID_concordant temp_aligned_ID_disc temp_flipped_ID_concordant temp_flipped_ID_disc temp_ambigious_ID_MAF_to_high_too_assess temp_ambigious_aligned_ID_concordant_same_allele temp_ambigious_aligned_ID_disc_same_allele temp_ambigious_aligned_ID_concordant_other_allele temp_ambigious_aligned_ID_disc_other_allele temp_ambigious_flipped_ID_concordant_same_allele temp_ambigious_flipped_ID_disc_same_allele temp_ambigious_flipped_ID_concordant_other_allele temp_ambigious_flipped_ID_disc_other_allele temp_non_biallelic_snps_ID temp_indels_ID_aligned_same_allele temp_indels_ID_aligned_other_allele 2> /dev/null  | sort | uniq -d | wc -l`

  if [[ $var_dup_2 == 0 ]]
  then
    echo -e "\nNo duplicate pairing of SNPs by position found after running awk_check_based_on_pos_1KG, continuing"
  else
    echo -e "\nWarning. Running awk_check_based_on_pos_1KG has assigned the same SNP pairing by position to more than one file. Exiting"
    exit
  fi  #Got various amounts as I found and fixed bugs, now melanoma, SDH, 610k, 610k_controls, QMEGA_omni_cases give 0 dups

  echo -e "\nFor ${out_name} the function awk_check_based_on_pos_1KG reports `cat temp_aligned_ID_concordant temp_aligned_ID_disc temp_flipped_ID_concordant temp_flipped_ID_disc temp_ambigious_ID_MAF_to_high_too_assess temp_ambigious_aligned_ID_concordant_same_allele temp_ambigious_aligned_ID_disc_same_allele temp_ambigious_aligned_ID_concordant_other_allele temp_ambigious_aligned_ID_disc_other_allele temp_ambigious_flipped_ID_concordant_same_allele temp_ambigious_flipped_ID_disc_same_allele temp_ambigious_flipped_ID_concordant_other_allele temp_ambigious_flipped_ID_disc_other_allele temp_non_biallelic_snps_ID temp_indels_ID_aligned_same_allele temp_indels_ID_aligned_other_allele 2> /dev/null | cut -d" " -f1 |  uniq -d | wc -l` SNPs from input bim file that are paired to more than one 1kG SNP"
    #SDH, 610k, 610k_control, QMEGA_omni_cases, WAMHS, HEI_cases all 0. OAG1 0 OAG2 0 IBD 0 QTWIN EPI 0 MIA 0

   echo -e "\nFor ${out_name} the function awk_check_based_on_pos_1KG reports `cat temp_aligned_ID_concordant temp_aligned_ID_disc temp_flipped_ID_concordant temp_flipped_ID_disc temp_ambigious_ID_MAF_to_high_too_assess temp_ambigious_aligned_ID_concordant_same_allele temp_ambigious_aligned_ID_disc_same_allele temp_ambigious_aligned_ID_concordant_other_allele temp_ambigious_aligned_ID_disc_other_allele temp_ambigious_flipped_ID_concordant_same_allele temp_ambigious_flipped_ID_disc_same_allele temp_ambigious_flipped_ID_concordant_other_allele temp_ambigious_flipped_ID_disc_other_allele temp_non_biallelic_snps_ID temp_indels_ID_aligned_same_allele temp_indels_ID_aligned_other_allele 2> /dev/null | cut -d" " -f2 |  uniq -d | wc -l` SNPs from 1KG that are paired by position to more than one SNP from the .bim file"
    #0 SDH, 610k, 610k_controls. Then how is the N higher than what I started? If there is an indel and a SNP at the same site in one set but not the other you will get a duplicate in the other set. tr " " "\n" would merge the cols into one and I could test for dups that way but that won't help since all the concordant files will be dups in col 1 and 2. Make a combined remap file and see if that will end up with any duplicates. WAMHS, HEI_cases and cont 0 IBD 0 QTWIN EPI MIA 0

  echo -e "\nFor ${out_name} the function awk_check_based_on_pos_1KG reports `cat temp_aligned_ID_disc temp_flipped_ID_disc temp_ambigious_aligned_ID_disc_same_allele temp_ambigious_aligned_ID_disc_other_allele temp_ambigious_flipped_ID_disc_same_allele temp_ambigious_flipped_ID_disc_other_allele 2> /dev/null  | wc -l` SNPs with overlapping positions but different IDs\n"
    #SDH 57547 610k_cas 4417 610k_con 3939 Omni_cases: 37753 AMFS_cas 37730 AMFS_con 37738 WAMHS 141401 HEI_cas 786 HEI_con: 812 OAG1 11356 OAG2 2592 IBD 140168 QTWIN 130117 EPIGENE 96268 MIA 33272

  #lets flip the SNPs before fixing names
  cat temp_flipped_ID_concordant temp_flipped_ID_disc temp_ambigious_flipped_ID_concordant_same_allele temp_ambigious_flipped_ID_disc_same_allele temp_ambigious_flipped_ID_concordant_other_allele temp_ambigious_flipped_ID_disc_other_allele  2> /dev/null | cut -d" " -f1 > temp_ID_flip_list; wc -l temp_ID_flip_list
    #SDH 28678 610k_cas 221 610k_con 62 Omni_case 2806 and AMFS_cases 2805 AMFS_con 2810 WAMHS 69915 HEI_cas 398 HEI_con 412 OAG1 5591 OAG2 1325 IBD 51 QTWIN 0 (SG aligned)  EPIGENE 0 MIA 16387

  #collate SNPs to exclude based on this match - ambigious, too high to match, triallelic, 
  echo ""
  cut -d" " -f1 temp_ambigious_ID_MAF_to_high_too_assess > ${dirI}/temp_${out_name}_ID_excludes_SNPs
  cut -d" " -f1 temp_non_biallelic_snps_ID >> ${dirI}/temp_${out_name}_ID_excludes_SNPs
  #13/07/2016 now add the unmappable SNPs
  cut -d" " -f2 ${dirI}/temp_${out_name}_unmappable_snps >> ${dirI}/temp_${out_name}_ID_excludes_SNPs; wc -l ${dirI}/temp_${out_name}_ID_excludes_SNPs
    #SDH 1178, 610k_cases 64, 610k_controls 51. QMEGA_omni_cases: 1142. AMFS_cases 1166. AMFS_Controls: 1158, WAMHS 1419 HEI_cases and con: 10. OAG1 342 OAG2 91 IBD 1307 QTWIN 1279  EPIGENE 999 MIA 733

  echo -e "\nFlipping SNPs to match 1KG, and excluding any ambigious SNPs with MAF too high to align and triallelic SNPs\n"
  plink_1.90 --threads 1 --bfile ${dirI}/temp_${out_name}_cleaning_7 --flip temp_ID_flip_list --make-bed --out ${dirI}/temp_${out_name}_cleaning_8 --exclude ${dirI}/temp_${out_name}_ID_excludes_SNPs
    #SDH --flip: 28678 SNPs flipped. --exclude: 939132 variants remaining.
    #610k_cases - --flip: 221 SNPs flipped --exclude: 548499 variants remaining. 548499 variants and 925 people pass filters and QC.
    #610k_controls --flip: 62 SNPs flipped. --exclude: 533116 variants remaining. 533116 variants and 3954 people pass filters and QC.
    #QMEGA_omni_cases: --flip: 2806 SNPs flipped. --exclude: 802710 variants remaining. 802710 variants and 694 people pass filters and QC.
    #AMFS_cases: --flip: 2805 SNPs flipped. --exclude: 802678 variants remaining. 802678 variants and 548 people pass filters and QC.
    #AMFS_controls: --flip: 2810 SNPs flipped. --exclude: 802663 variants remaining. 802663 variants and 430 people pass filters and QC.
    #WAMHS --flip: 69915 SNPs flipped --exclude: 915499 variants remaining. 915499 variants and 1273 people pass filters and QC.
    #HEI_cases: --flip: 398 SNPs flipped --exclude: 218694 variants remaining 218694. variants and 1217 people pass filters and QC
    #HEI_con: --flip: 412 SNPs flipped  --exclude: 221008 variants remaining (none excluded) 221008 variants and 1221 people pass filters and QC
    #OAG1 --flip: 5591 SNPs flipped. --exclude: 855941 variants remaining. 855941 variants and 651 people pass filters and QC.
    #OAG2 --flip: 1325 SNPs flipped --exclude: 707444 variants remaining 707444 variants and 645 people pass filters and QC.
    #IBD --flip: 51 SNPs flipped. --exclude: 912898 variants remaining. 912898 variants and 929 people pass filters and QC.
    #QTWIN --flip: 0 SNPs flipped. --exclude: 450024 variants remaining. 450024 variants and 981 people pass filters and QC.
    #EPIGENE --flip: 0 SNPs flipped. --exclude: 399153 variants remaining. 399153 variants and 784 people pass filters and QC.
    #MIA --flip: 16387 --exclude: 486958 variants remaining. 486958 variants and 1968 people pass filters and QC.
  echo ""
  #now make a remapping list that makes sense. if only one is an rsID use that; if not use 1kg ID. 25/07/2016 found a bug - one match wasn't looking for ^rs so some SNPs (exm-rs) were falling into 2 cat
  awk '{
  if($2~/^rs/ && !($1~/^rs/))
    print $1,$2;
  if($2~/^rs/ && $1~/^rs/)
    print $1,$2;
  if(!($2~/^rs/) && !($1~/^rs/))
    print $1,$2;
  }' <( cat temp_aligned_ID_disc temp_flipped_ID_disc temp_ambigious_aligned_ID_disc_same_allele temp_ambigious_aligned_ID_disc_other_allele temp_ambigious_flipped_ID_disc_same_allele temp_ambigious_flipped_ID_disc_other_allele 2> /dev/null  ) > temp_ID_remap; wc -l temp_ID_remap
    #SDH 57099 610k 4403 610k_con 3928 Omni_cas 37703 AMFS_cas 37680 AMFS_con 37688 WAMHS 141273 Lots of exm IDs. HEI_cases: 768. HEI_con 794. OAG1 11090 OAG2 2515 IBD 140042 QTWIN 130109  EPIGENE 96260 MIA 31813
  #check those that weren't remapped
  echo -e "\nExample of the SNPs not being remapped - should be of the form rsXXX chrX:XXXXXX:"
  awk 'NR==FNR {a[$1,$2];next} !(($1,$2) in a)' temp_ID_remap <( cat temp_aligned_ID_disc temp_flipped_ID_disc temp_ambigious_aligned_ID_disc_same_allele temp_ambigious_aligned_ID_disc_other_allele temp_ambigious_flipped_ID_disc_same_allele temp_ambigious_flipped_ID_disc_other_allele 2> /dev/null ) | head
    #SDH; all look right  - rsIDs in illumina, chr1:XXXXX in 1kg
  echo ""
  plink_1.90 --threads 1 --bfile ${dirI}/temp_${out_name}_cleaning_8 --update-name temp_ID_remap --make-bed --out ${dirI}/temp_${out_name}_cleaning_9
    #SDH --update-name: 57099 values updated. 939132 variants and 570 people pass filters and QC.
    #610k_cases --update-name: 4403 values updated. 548499 variants and 925 people pass filters and QC.
    #610k_controls --update-name: 3928 values updated 533116 variants and 3954 people pass filters and QC.
    #QMEGA_omni_cases: --update-name: 37703 values updated. 802710 variants and 694 people pass filters and QC
    #AMFS_cases:  --update-name: 37680 values updated. 802678 variants and 548 people pass filters and QC.
    #AMFS_controls: --update-name: 37688 values updated. 802663 variants and 430 people pass filters and QC.
    #WAMHS --update-name: 141273 values updated. 915499 variants and 1273 people pass filters and QC
    #HEI_cases: --update-name: 768 values updated. 218694 variants and 1217 people pass filters and QC
    #HEI_con: --update-name: 794 values updated. 221008 variants and 1221 people pass filters and QC.
    #OAG1 --update-name: 11090 values updated. 855941 variants and 651 people pass filters and QC.
    #OAG2 --update-name: 2515 values updated. 707444 variants and 645 people pass filters and QC.
    #IBD --update-name: 140042 values updated. 912898 variants and 929 people pass filters and QC.
    #QTWIN --update-name: 130109 values updated. 450024 variants and 981 people pass filters and QC.
    #EPIGENE --update-name: 96260 values updated. 399153 variants and 784 people pass filters and QC.
    #MIA --update-name: 31813 values updated 486958 variants and 1968 people pass filters and QC.

  #check for duplicates - not a terrible problem , just need to merge them on each other but I am not doing that as a default so add a check 
  var_dup_3=`awk '{print $2}' ${dirI}/temp_${out_name}_cleaning_9.bim | sort | uniq -d | wc -l`

  if [[ $var_dup_3 == 0 ]]
  then
    echo -e "\nNo duplicate IDs following updating IDs; continuing"
  else
    echo -e "\nWarning. Running awk_check_based_on_pos_1KG has assigned the same SNP ID to two SNPs. May be resolvable by merging with itself\n"
    exit
  fi

  #interim clean before rerunning the checks
  rm -f ${dirI}/temp_${out_name}_update_map_bp.txt ${dirI}/temp_${out_name}_cleaning_1.* ${dirI}/temp_${out_name}_cleaning_2.* ${dirI}/temp_${out_name}_cleaning_3.* ${dirI}/temp_${out_name}_cleaning_4* ${dirI}/temp_${out_name}_cleaning_5.* ${dirI}/temp_${out_name}_cleaning_6* ${dirI}/temp_${out_name}_cleaning_7* ${dirI}/temp_${out_name}_cleaning_8.* ${dirI}/temp_${out_name}_BED* ${dirI}/temp_${out_name}_cleaning_overlapping* ${dirI}/temp_${out_name}_chr* ${dirI}/temp_${out_name}_bp* ${dirI}/temp_${out_name}_flip* ${dirI}/temp_${out_name}_unmappable* ${dirI}/temp_${out_name}_update_map_chr.txt ${dirI}/temp_${out_name}_ID_excludes_SNPs temp_flipped_concordant temp_flipped_chr_* temp_aligned_concordant temp_aligned_chr_* temp_ambigious_aligned_* temp_ambigious_flipped_* temp_non_biallelic_snps* temp_indels_aligned_* temp_aligned_ID_* temp_flipped_ID_* temp_ambigious* MAF_to_high_too_assess temp_ID_flip_list temp_ID_remap temp_indels_ID_aligned_*_allele
  
  plink_1.90 --threads 1 --bfile ${dirI}/temp_${out_name}_cleaning_9 --freq --out ${dirI}/temp_${out_name}_cleaning_9 --silent
  echo ""
  awk 'NR==FNR {a[$2]=$1" "$3" "$4" "$5" "$9" "$13;next} $2 in a {print $0,a[$2]}' <(gunzip -dc ${dirI}/ALL_1000G_phase1integrated_v3_ALL_chr_legends_with_chromosome.txt.gz) <( awk 'NR==FNR {a[$2]=$4;next} $2 in a {print $0,a[$2]}' <(sed 's/ \+/ /g' ${dirI}/temp_${out_name}_cleaning_9.bim | sed 's/^ //g')  <(sed 's/ \+/ /g' ${dirI}/temp_${out_name}_cleaning_9.frq | sed 's/^ //g' | sed 's/ $//g') ) > ${dirI}/temp_${out_name}_cleaning_9_1KG_merge_rsIDs; wc -l ${dirI}/temp_${out_name}_cleaning_9_1KG_merge_rsIDs
    #SDH 929876 610k_cas 547357 610k_con 532183 Omni_cases 799184 AMFS_cas 799152 AMFS_con 799137 WAMHS 807561 HEI_cas 217364 HEI_con 219629 OAG1 851560 OAG2 700942 IBD 805763 QTWIN 378781 EPIGENE 337213 MIA 463931

  echo -e "\n+++++++++++Re running awk_check_based_on_ID_1KG function to check SNPs paired with 1KG by ID have correct position following updating and aligning. Aligning amibigious SNPs where MAF < $MAF and MAF difference with 1KG is < $DIFF \n"

  #see 0.4 for the function used here but this will assign the SNPs to various categories - e.g. aligned_concordant - same ID, same position, same strand and so on
  awk -v MAF=${MAF} -v DIFF=${DIFF} -f /working/lab_stuartma/scratch/matthewL/melanoma_meta_analysis/SCRIPTS/awk_check_based_on_ID_1KG ${dirI}/temp_${out_name}_cleaning_9_1KG_merge_rsIDs

  #how many are aligned and concordant?
  cat temp_aligned_concordant temp_ambigious_aligned_concordant_same_allele temp_ambigious_aligned_concordant_other_allele temp_indels_aligned_same_allele temp_indels_aligned_other_allele  2> /dev/null | wc -l
    #SDH - 908938 (missing chrX vs X_nonPAR) 610k_cas 538990 610k_con 523946 Omni_cases: 788810 AMFS_cas 788778 AMFS_con 788760 WAMHS 792006 HEI_cas 210916. HEI_con 213161 OAG1 833576 OAG2 684518 IBD 790058 QTWIN 371541 EPIGENE 330601 MIA 457534

  cat temp_flipped_chr_disc_bp_conc temp_aligned_chr_disc_bp_disc temp_aligned_chr_disc_bp_conc temp_ambigious_aligned_chr_disc_bp_conc_same_allele temp_ambigious_aligned_chr_disc_bp_disc_same_allele temp_ambigious_aligned_chr_disc_bp_conc_other_allele temp_ambigious_aligned_chr_disc_bp_disc_other_allele temp_ambigious_flipped_chr_disc_bp_conc_same_allele temp_ambigious_flipped_chr_disc_bp_disc_same_allele temp_ambigious_flipped_chr_disc_bp_conc_other_allele temp_ambigious_flipped_chr_disc_bp_disc_other_allele 2> /dev/null > temp_${out_name}_chr_disc

  #assuming the earlier steps are correct all the SNPs should be concordant once you discount chr_disc solely due to 23 in PLINK and X_nonPAR in 1KG. So tweak earlier comparisons to add the concordant files, and those chr_disc that are just 23 vs X_nonPAR; and flag if the in total doesn't match this out total (that is, there is a SNP that is discordant at position or CHR with 1000 genome. 25/07/2016 - also including chr25/X_nonPAR
  var_in=`sed -n '$=' ${dirI}/temp_${out_name}_cleaning_9_1KG_merge_rsIDs`
  var_out=`cat temp_aligned_concordant temp_ambigious_aligned_concordant_same_allele temp_ambigious_aligned_concordant_other_allele temp_indels_aligned_same_allele temp_indels_aligned_other_allele <(awk 'NR==FNR {a[$1];next} $2 in a' temp_${out_name}_chr_disc ${dirI}/temp_${out_name}_cleaning_9_1KG_merge_rsIDs | awk '($1==23 && $8=="X_nonPAR") || ($1==25 && $8=="X_nonPAR")' )  2> /dev/null  | wc -l`
 
  if [[ $var_in == $var_out ]]
  then
    echo -e "\nRemapping and renaming means all SNPs with a common ID with 1000 genomes are concordant for CHR and BP once you account for CHRX naming convention differences"
  else
    echo -e "\nWarning. After remapping and renaming there are SNPs that overlap in ID between PLINK and 1KG but not CHR and/or BP. May be due to indels having their alleles fixed, and now able to matched by the function. Updating positions to see if this resolves."
    #27/03/2016 610k - at first were getting some but found I was merging all the lists of remap/flip SNPs etc, fixed. 07/07/2016 now getting three for QMEGA_cases, and it is due to three indels that now have the right alleles - this has picked up the positions are wrong. Looks like a disagreement wether you could the indel or not in the bp rang (e.g. 3bp del is 3pb dif between illumina and 1kg. 
      #22 rs11570768 A AGGCGG 0.05988 1386 38507985 22 38507980 AGGCGG A 0.0752 0.0752
      #22 rs28915376 C CAA 0.04467 1388 43091219 22 43091216 CAA C 0.0396 0.0396
      #22 rs3216827 TG T 0.0879 1388 43091567 22 43091565 T TG 0.1003 0.1003
    
      awk 'NR==FNR {a[$1];next} $2 in a {print $2,$9}' temp_aligned_chr_conc_bp_disc ${dirI}/temp_${out_name}_cleaning_9_1KG_merge_rsIDs > ${dirI}/temp_${out_name}_cleaning_9_update_bp
      echo ""
      plink_1.90 --threads 1 --bfile ${dirI}/temp_${out_name}_cleaning_9 --update-map ${dirI}/temp_${out_name}_cleaning_9_update_bp --make-bed --out ${dirI}/temp_${out_name}_cleaning_9_temp
      #QMEGA_omni_cases; AMFS_cases, AMFS_controls, MIA: --update-map: 3 values updated. MIA : 1 update
      echo ""
      #07/07/2016 restore and rerun same check      
      plink_1.90 --threads 1 --bfile ${dirI}/temp_${out_name}_cleaning_9_temp --make-bed --out ${dirI}/temp_${out_name}_cleaning_9
      plink_1.90 --threads 1 --bfile ${dirI}/temp_${out_name}_cleaning_9 --freq --out ${dirI}/temp_${out_name}_cleaning_9 --silent
      echo "" 
      awk 'NR==FNR {a[$2]=$1" "$3" "$4" "$5" "$9" "$13;next} $2 in a {print $0,a[$2]}' <(gunzip -dc ${dirI}/ALL_1000G_phase1integrated_v3_ALL_chr_legends_with_chromosome.txt.gz) <( awk 'NR==FNR {a[$2]=$4;next} $2 in a {print $0,a[$2]}' <(sed 's/ \+/ /g' ${dirI}/temp_${out_name}_cleaning_9.bim | sed 's/^ //g')  <(sed 's/ \+/ /g' ${dirI}/temp_${out_name}_cleaning_9.frq | sed 's/^ //g' | sed 's/ $//g') ) > ${dirI}/temp_${out_name}_cleaning_9_1KG_merge_rsIDs; wc -l ${dirI}/temp_${out_name}_cleaning_9_1KG_merge_rsIDs

      echo -e "\n+++As have had to update positions a second time checking if those now allows all common IDs to line up to 1KG positions"

      rm -f temp_flipped_concordant temp_flipped_chr_* temp_aligned_concordant temp_aligned_chr_* temp_ambigious_aligned_* temp_ambigious_flipped_* temp_non_biallelic_snps* temp_indels_aligned_* temp_aligned_ID_* temp_flipped_ID_* temp_ambigious*MAF_to_high_too_assess temp_ID_flip_list temp_ID_remap temp_indels_ID_aligned_*_allele

      awk -v MAF=${MAF} -v DIFF=${DIFF} -f /working/lab_stuartma/scratch/matthewL/melanoma_meta_analysis/SCRIPTS/awk_check_based_on_ID_1KG ${dirI}/temp_${out_name}_cleaning_9_1KG_merge_rsIDs

      cat temp_flipped_chr_disc_bp_conc temp_aligned_chr_disc_bp_disc temp_aligned_chr_disc_bp_conc temp_ambigious_aligned_chr_disc_bp_conc_same_allele temp_ambigious_aligned_chr_disc_bp_disc_same_allele temp_ambigious_aligned_chr_disc_bp_conc_other_allele temp_ambigious_aligned_chr_disc_bp_disc_other_allele temp_ambigious_flipped_chr_disc_bp_conc_same_allele temp_ambigious_flipped_chr_disc_bp_disc_same_allele temp_ambigious_flipped_chr_disc_bp_conc_other_allele temp_ambigious_flipped_chr_disc_bp_disc_other_allele 2> /dev/null > temp_${out_name}_chr_disc

      var_out2=`cat temp_aligned_concordant temp_ambigious_aligned_concordant_same_allele temp_ambigious_aligned_concordant_other_allele temp_indels_aligned_same_allele temp_indels_aligned_other_allele <(awk 'NR==FNR {a[$1];next} $2 in a' temp_${out_name}_chr_disc ${dirI}/temp_${out_name}_cleaning_9_1KG_merge_rsIDs | awk '$1==23 && $8=="X_nonPAR"' )  2> /dev/null  | wc -l`
     
    #redo check
    if [[ $var_in == $var_out2 ]]
    then
      echo -e "\nSecond round of remapping and renaming means all SNPs with a common ID with 1000 genomes are concordant for CHR and BP once you account for CHRX naming convention differences"
    else
      echo -e "\nWarning. After two rounds of remapping and renaming there are SNPs that overlap in ID between PLINK and 1KG but not CHR and/or BP. Exiting to allow correct as may require another round of aligning post ID renaming."
      exit
    fi #Now fixes the problem for QMEGA_cases and AMFS cases
  fi #checking for any new base pair mismatches introduced by renaming or removing indels.

  #test for unmappable ambigious SNPs - shouldn't be any
  if [[ -a temp_ambigious_MAF_to_high_too_assess ]]
  then
    echo -e "\nWarning. After remapping and renaming there are still ambigious SNPs with too high MAF to be mapped. Exiting."
    exit
  fi  #SDH 0, which is good

  var_dup_4=`cat temp_flipped_concordant temp_flipped_chr_disc_bp_conc temp_aligned_concordant temp_aligned_chr_disc_bp_disc temp_aligned_chr_disc_bp_conc temp_aligned_chr_conc_bp_disc temp_flipped_chr_con_bp_disc temp_ambigious_aligned_concordant_same_allele temp_ambigious_aligned_chr_disc_bp_conc_same_allele temp_ambigious_aligned_chr_conc_bp_disc_same_allele temp_ambigious_aligned_chr_disc_bp_disc_same_allele temp_ambigious_aligned_concordant_other_allele temp_ambigious_aligned_chr_disc_bp_conc_other_allele temp_ambigious_aligned_chr_conc_bp_disc_other_allele temp_ambigious_aligned_chr_disc_bp_disc_other_allele temp_ambigious_flipped_concordant_same_allele temp_ambigious_flipped_chr_disc_bp_conc_same_allele temp_ambigious_flipped_chr_conc_bp_disc_same_allele temp_ambigious_flipped_chr_disc_bp_disc_same_allele temp_ambigious_flipped_concordant_other_allele temp_ambigious_flipped_chr_disc_bp_conc_other_allele temp_ambigious_flipped_chr_conc_bp_disc_other_allele temp_ambigious_flipped_chr_disc_bp_disc_other_allele temp_ambigious_MAF_to_high_too_assess temp_non_biallelic_snps temp_indels_aligned_same_allele temp_indels_aligned_other_allele temp_flipped_chr_disc_bp_disc 2> /dev/null  | sort | uniq -d | wc -l`

  if [[ $var_dup_4 == 0 ]]
  then
    echo -e "\nNo duplicate SNPs found after running awk_check_based_on_ID_1KG on the remapped and aligned set, continuing\n"
  else
    echo -e "\nWarning. Running awk_check_based_on_ID_1KG has assigned the same SNP to more than one file. Exiting"
    exit
  fi  #Got various amounts as I found and fixed bugs, now melanoma, SDH give 0 dups

  #how many output in total by the function. supress the cat warning about missing files by 2> /dev/null.
  echo -e "For remapped and aligned ${out_name} total of `cat temp_flipped_concordant temp_flipped_chr_disc_bp_conc temp_aligned_concordant temp_aligned_chr_disc_bp_disc temp_aligned_chr_disc_bp_conc temp_aligned_chr_conc_bp_disc temp_flipped_chr_con_bp_disc temp_ambigious_aligned_concordant_same_allele temp_ambigious_aligned_chr_disc_bp_conc_same_allele temp_ambigious_aligned_chr_conc_bp_disc_same_allele temp_ambigious_aligned_chr_disc_bp_disc_same_allele temp_ambigious_aligned_concordant_other_allele temp_ambigious_aligned_chr_disc_bp_conc_other_allele temp_ambigious_aligned_chr_conc_bp_disc_other_allele temp_ambigious_aligned_chr_disc_bp_disc_other_allele temp_ambigious_flipped_concordant_same_allele temp_ambigious_flipped_chr_disc_bp_conc_same_allele temp_ambigious_flipped_chr_conc_bp_disc_same_allele temp_ambigious_flipped_chr_disc_bp_disc_same_allele temp_ambigious_flipped_concordant_other_allele temp_ambigious_flipped_chr_disc_bp_conc_other_allele temp_ambigious_flipped_chr_conc_bp_disc_other_allele temp_ambigious_flipped_chr_disc_bp_disc_other_allele temp_ambigious_MAF_to_high_too_assess temp_non_biallelic_snps temp_indels_aligned_same_allele temp_indels_aligned_other_allele temp_flipped_chr_disc_bp_disc  2> /dev/null | wc -l` SNPs output by awk_check_based_on_ID_1KG"
   #SDH 929876 Omni_cas 799184 AMFS_cas 799152 AMFS_con 799144 610k_cas 547357 610k_con 532183 WAMHS 808757. HEI_cas 217364. HEI_con 219629 OAG1 851560 OAG2 700942 IBD 806748 QTWIN 378781 EPIGENE 337213 MIA 463931

  #which ones are missing? supress the cat warning about missing files by 2> /dev/null.
  var_unmappable_SNPs=`awk 'NR==FNR {a[$1];next} !($2 in a) ' <(cat temp_flipped_concordant temp_flipped_chr_disc_bp_conc temp_aligned_concordant temp_aligned_chr_disc_bp_disc temp_aligned_chr_disc_bp_conc temp_aligned_chr_conc_bp_disc temp_flipped_chr_con_bp_disc temp_ambigious_aligned_concordant_same_allele temp_ambigious_aligned_chr_disc_bp_conc_same_allele temp_ambigious_aligned_chr_conc_bp_disc_same_allele temp_ambigious_aligned_chr_disc_bp_disc_same_allele temp_ambigious_aligned_concordant_other_allele temp_ambigious_aligned_chr_disc_bp_conc_other_allele temp_ambigious_aligned_chr_conc_bp_disc_other_allele temp_ambigious_aligned_chr_disc_bp_disc_other_allele temp_ambigious_flipped_concordant_same_allele temp_ambigious_flipped_chr_disc_bp_conc_same_allele temp_ambigious_flipped_chr_conc_bp_disc_same_allele temp_ambigious_flipped_chr_disc_bp_disc_same_allele temp_ambigious_flipped_concordant_other_allele temp_ambigious_flipped_chr_disc_bp_conc_other_allele temp_ambigious_flipped_chr_conc_bp_disc_other_allele temp_ambigious_flipped_chr_disc_bp_disc_other_allele temp_ambigious_MAF_to_high_too_assess temp_non_biallelic_snps temp_indels_aligned_same_allele temp_indels_aligned_other_allele temp_flipped_chr_disc_bp_disc 2> /dev/null  ) ${dirI}/temp_${out_name}_cleaning_9_1KG_merge_rsIDs | wc -l`

  if [[ $var_unmappable_SNPs != 0 ]]
  then
    echo -e "\nWarning. Running awk_check_based_on_ID_1KG a second time has found additional SNPs that cannot be mapped. Exiting"
    exit
  fi  #Got various amounts as I found and fixed bugs, now melanoma, SDH give 0 unammapable as do 60k, 610k controls

  echo -e "\n+++++++++++Re running awk_check_based_on_pos_1KG function to check SNPs paired with 1KG by position have correct IDs following updating and aligning. Aligning amibigious SNPs where MAF < $MAF and MAF difference with 1KG is < $DIFF"

  echo -e "\n+++++++Now pairing SNPs with 1KG based on CHR and POS and checking for ID mismatches using awk_check_based_on_pos_1KG function\n"
  plink_1.90 --threads 1 --bfile ${dirI}/temp_${out_name}_cleaning_9 --freq --out ${dirI}/temp_${out_name}_cleaning_9 --silent

  #13/07/2016 as above; this command works except then there is overlapping indels; replaced with the newer version where match PLINK to 1kg rather than the reverse
  ##awk 'NR==FNR {a[$1,$3]=$1" "$3" "$4" "$5" "$9" "$13" "$2;next} ($1,$7) in a {print $0,a[$1,$7]}' <(gunzip -dc ${dirI}/ALL_1000G_phase1integrated_v3_ALL_chr_legends_with_chromosome.txt.gz) <( awk 'NR==FNR {a[$2]=$4;next} $2 in a {print $0,a[$2]}' <(sed 's/ \+/ /g' ${dirI}/temp_${out_name}_cleaning_9.bim | sed 's/^ //g')  <(sed 's/ \+/ /g' ${dirI}/temp_${out_name}_cleaning_9.frq | sed 's/^ //g' | sed 's/ $//g') ) > ${dirI}/temp_${out_name}_cleaning_9_1KG_merge_pos; wc -l ${dirI}/temp_${out_name}_cleaning_9_1KG_merge_pos
  awk 'NR==FNR  {a[$1,$7]=$0;next} ($1,$3) in a {print a[$1,$3],$1,$3,$4,$5,$9,$13,$2}' <( awk 'NR==FNR {a[$2]=$4;next} $2 in a {print $0,a[$2]}'  <(sed 's/ \+/ /g' ${dirI}/temp_${out_name}_cleaning_9.bim | sed 's/^ //g')  <(sed 's/ \+/ /g' ${dirI}/temp_${out_name}_cleaning_9.frq | sed 's/^ //g' | sed 's/ $//g') ) <(gunzip -dc ${dirI}/ALL_1000G_phase1integrated_v3_ALL_chr_legends_with_chromosome.txt.gz) > ${dirI}/temp_${out_name}_cleaning_9_1KG_merge_pos; wc -l ${dirI}/temp_${out_name}_cleaning_9_1KG_merge_pos
    #SDH: 909685 so ~200 more SNPs; 610k_cas 539181 610k_con 524117 Omni_cases 789163 AMFS_cas 789131, AMFS_con 789113 WAMHS 792360 HEI_cas 210994. HEI_con 213244 OAG1 834120 OAG2 684823 IBD 790408 QTWIN 371654 EPIGENE 330689 MIA 465572
  echo ""
  wc -l ${dirI}/temp_${out_name}_cleaning_9.bim
    #SDH 939132 610k_cas 548499 610k_con 533116 Omni_cases 802710 AMFS_cas 802678 AMFS_con 802667 WAMHS 915499 HEI_cas: 218694 HEI_con 221008 OAG1 855941 OAG2 707444 IBD 914239 QTWIN 450024 EPIGENE 399153 MIA 486958
  echo -e "\nFor ${out_name} there are `awk '$2!=$14' ${dirI}/temp_${out_name}_cleaning_9_1KG_merge_pos | wc -l` SNPs with mismatching IDs but the same position. This will include true ID differences and overlapping indels/SNPs, and SNPs where the 1KG ID was out of date (e.g. PLINK rsXXX, 1KG chrX:XXXX)\n"
    #SDH: 747 - all indels vs SNPs, or 1KG IDs are CHR:POS. 610k_cas 191. 610k_con 171. Omni_cas, AMFS_cas, AMFS_con 350. HEI_cases 78. HEI_con 83 OAG1 544 OAG2 305 IBD 434 WAMHS 354 QTWIN 113  EPIGENE 88 MIA 8037
  #see 0.4 for this function
  awk -v MAF=${MAF} -v DIFF=${DIFF} -f /working/lab_stuartma/scratch/matthewL/melanoma_meta_analysis/SCRIPTS/awk_check_based_on_pos_1KG  ${dirI}/temp_${out_name}_cleaning_9_1KG_merge_pos

  cat temp_aligned_ID_concordant temp_ambigious_aligned_ID_concordant_same_allele temp_ambigious_aligned_ID_concordant_other_allele 2> /dev/null | wc -l
    #SDH 908938 610k_cas 538990 610k_con 523946 Omni_cases 802713 AMFS_cas 788781 AMFS_con 788770 WAMHS 792006 HEI_cas 210916. HEI_con 213161 OAG1 833576 OAG2 684518 IBD 790058 QTWIN 371541  EPIGENE 330601 MIA 457535
  echo ""
  cat temp_aligned_ID_disc temp_flipped_ID_disc temp_ambigious_aligned_ID_disc_same_allele temp_ambigious_aligned_ID_disc_other_allele temp_ambigious_flipped_ID_disc_same_allele temp_ambigious_flipped_ID_disc_other_allele 2> /dev/null | wc -l
    #SDH 448.  610k_cases 14, 610k_controls 11, Omni_cases, AMFS_cases, and AMFS_con 82. WAMHS 129 HEI_cases and con 18. OAG1 266 OAG2 77 IBD 127 QTWIN 9  EPIGENE 9
   #so out of the ones I can line up do any pass my naming test above? MIA 7767
  var_need_rename=`awk '{
  if($2~/^rs/ && !($1~/^rs/))
    print $1,$2;
  if($2~/^rs/ && $1~/^rs/)
    print $1,$2;
  if(!($2~/^rs/) && !($1~/^rs/))
    print $1,$2;
  }' <( cat temp_aligned_ID_disc temp_flipped_ID_disc temp_ambigious_aligned_ID_disc_same_allele temp_ambigious_aligned_ID_disc_other_allele temp_ambigious_flipped_ID_disc_same_allele temp_ambigious_flipped_ID_disc_other_allele 2> /dev/null  ) | wc -l`

  if [[ $var_need_rename != 0 ]]
  then
    echo -e "\nWarning. Running awk_check_based_on_pos_1KG a second time has found additional SNPs that need renaming. This may be due to overlapping exm chip SNPs since removed, rerunning the renaming step\n"
     #25/07/2016 WAMHS has 18, which is the only set so far, and looking at them they are overlapping exm SNPs, where the first one is likely bogus, and removed, allowing the correct one to now be mapped here. 
    awk '{
      if($2~/^rs/ && !($1~/^rs/))
        print $1,$2;
      if($2~/^rs/ && $1~/^rs/)
        print $1,$2;
      if(!($2~/^rs/) && !($1~/^rs/))
        print $1,$2;
      }' <( cat temp_aligned_ID_disc temp_flipped_ID_disc temp_ambigious_aligned_ID_disc_same_allele temp_ambigious_aligned_ID_disc_other_allele temp_ambigious_flipped_ID_disc_same_allele temp_ambigious_flipped_ID_disc_other_allele 2> /dev/null  ) > temp_ID_remap; wc -l temp_ID_remap
        #WAMHS 18 IBD 1 (both exm chip) MIA 595
    echo ""
    plink_1.90 --threads 1 --bfile ${dirI}/temp_${out_name}_cleaning_9 --update-name temp_ID_remap --make-bed --out ${dirI}/temp_${out_name}_cleaning_9_temp
      #WAMHS and IBD --1 value updated MIA --update-name: 595 
    echo ""
    plink_1.90 --threads 1 --bfile ${dirI}/temp_${out_name}_cleaning_9_temp --make-bed --out  ${dirI}/temp_${out_name}_cleaning_9
    echo ""

    #Again check for duplicates - not a terrible problem , just need to merge them on each other but I am not doing that as a default so add a check 
    var_dup_3=`awk '{print $2}' ${dirI}/temp_${out_name}_cleaning_9.bim | sort | uniq -d | wc -l`
    if [[ $var_dup_3 == 0 ]]
    then
      echo -e "No duplicate IDs following second round of updating IDs; continuing\n"
    else
      echo -e "\nWarning. Running awk_check_based_on_pos_1KG has assigned the same SNP ID to two SNPs after a second round of renaming. May be resolvable by merging with itself\n"
      exit
    fi
  fi  #Check for SNPs that can still be renamed.

  #last test - test for any ambigious SNPs with MAF > 0.35; suggests I messed up somewhere (or there are ones in there that don't align to 1000 genomes? Still drop them in case they align to the HRC or other datasets - 54 with SDH,. 

  awk 'NR==FNR {a[$2]=$5;next} $2 in a {print $0,a[$2]}' <(awk ' $5!="NA" && $5> '''${MAF}''' && ( ($3=="A" && $4=="T") || ($3=="T" && $4=="A") || ($3=="G" && $4=="C") || ($3=="C" && $4=="G") )' ${dirI}/temp_${out_name}_cleaning_9.frq)  ${dirI}/temp_${out_name}_cleaning_9.bim > ${dirI}/temp_${out_name}_cleaning_9_ambigious_high_MAF; wc -l ${dirI}/temp_${out_name}_cleaning_9_ambigious_high_MAF
    #SDH 66; 610k_cases and 610k controls 10. Omni_case: 66. AMFS_cas 65. AMFS_con 74. WAMHS 22. HEI_cases and con: 0 OAG1 38 OAG2 51 IBD 18 QTWIN 6 EPIGENE 5 MIA 30
  echo ""
  #do any match with 1KG by ID (should be zero, unless I messed up)
  var_rem_amb_ID=`awk 'NR==FNR {a[$2];next} $2 in a' <(gunzip -dc ${dirI}/ALL_1000G_phase1integrated_v3_ALL_chr_legends_with_chromosome.txt.gz) ${dirI}/temp_${out_name}_cleaning_9_ambigious_high_MAF | wc -l`
  if [[ $var_rem_amb_ID != 0 ]]
  then
    echo -e "\nWarning. There remains ambigious high MAF SNPS that map to 1KG by ID, an earlier step must have failed. Exiting"
    exit
  fi  #

  #do any match with 1KG by pos (should be zero, unless I messed up). 30/06/2016 though the earlier functions werent working as SHD still have 4 overlaps but this check here is too simple. There are 4 high maf ambigious SNPs that share a position with a different INDEL, and shouldn't be counted so added second step to check alleles, which means needs to print out allleles from 1kG
  var_rem_amb_POS=`awk 'NR==FNR {a[$1,$3]=$4" "$5;next} ($1,$4) in a {print $0,a[$1,$4]}'  <(gunzip -dc ${dirI}/ALL_1000G_phase1integrated_v3_ALL_chr_legends_with_chromosome.txt.gz) <(sed 's/\t\+/ /g' ${dirI}/temp_${out_name}_cleaning_9_ambigious_high_MAF) | awk '($5==$8 && $6==$9) || ($5==$9 && $6==$8) ' | wc -l`

  if [[ $var_rem_amb_POS != 0 ]]
  then
    echo -e "\nWarning. There remains ambigious high MAF SNPS that map to 1KG by pos, an earlier step must have failed. Exiting"
    exit
  fi  

  #12/07/2016 Need to filter out any SNP that isn't in 1KG and is ambigious - can't easily align it between sets. Need to filter 1kg for amb as well as overlapping indels are causing a few ambigious SNPs to slip through (as in they were incorrectly passing the positional match)
  awk 'NR==FNR {a[$2];next} !($2 in a) {print $2}'  <(gunzip -dc ${dirI}/ALL_1000G_phase1integrated_v3_ALL_chr_legends_with_chromosome.txt.gz | awk '($4=="A" && $5=="T") || ($4=="T" && $5=="A") || ($4=="G" && $5=="C") || ($4=="C" && $5=="G")' )  <(awk '($3=="A" && $4=="T") || ($3=="T" && $4=="A") || ($3=="G" && $4=="C") || ($3=="C" && $4=="G") ' ${dirI}/temp_${out_name}_cleaning_9.frq)  >  ${dirI}/temp_${out_name}_cleaning_9_non_1KG_amb_rsID; wc -l ${dirI}/temp_${out_name}_cleaning_9_non_1KG_amb_rsID
    #SDH 798, 610k_cases 30, 610k_controls 29, QMEGA_omni_cases, AMFS_cases, AMFS_controls 236. WAMHS 12466. HEI_cases and con 0 OAG1 337 OAG2 347 IBD 12274 QTWIN 7205 EPIGENE 6095 MIA 1585
  echo ""
  awk 'NR==FNR {a[$1,$3];next} !(($1,$4) in a) {print $2}' <(gunzip -dc ${dirI}/ALL_1000G_phase1integrated_v3_ALL_chr_legends_with_chromosome.txt.gz | awk '($4=="A" && $5=="T") || ($4=="T" && $5=="A") || ($4=="G" && $5=="C") || ($4=="C" && $5=="G")' )  <(awk '($5=="A" && $6=="T") || ($5=="T" && $6=="A") || ($5=="G" && $6=="C") || ($5=="C" && $6=="G") ' ${dirI}/temp_${out_name}_cleaning_9.bim) >  ${dirI}/temp_${out_name}_cleaning_9_non_1KG_amb_pos; wc -l ${dirI}/temp_${out_name}_cleaning_9_non_1KG_amb_pos
    #SDH 1019, 610k_cases 55, 610_con 48, Omni_cases and AMFS_cas 353, AMFS_con 356, WAMHS 12491 - all exm. HEI_Cases and con 0 OAG1 514 OAG2 377 IBD 12298 QTWIN 7238 EPIGENE 6100 MIA 1452
  echo ""
  #so need to drop any SNP that is in BOTH files, not just one. 
  cat  ${dirI}/temp_${out_name}_cleaning_9_non_1KG_amb_rsID ${dirI}/temp_${out_name}_cleaning_9_non_1KG_amb_pos | sort | uniq -d > ${dirI}/temp_${out_name}_cleaning_9_non_1KG_amb_exclude_pre; wc -l ${dirI}/temp_${out_name}_cleaning_9_non_1KG_amb_exclude_pre
    #SDH 762, 610k_cases 30, 610k_controls 29, QMEGA_omni_cases, AMFS_cases, and AMFS_controls 236, WAMHS 12465, HEI_Cases and con 0 OAG1 320 OAG2 347 IBD 12273 QTWIN 7204 EPIGENE 6094 MIA 1408

  echo -e "\nTo be safe removing any remaining ambigious SNPs that do not map to 1KG\n"
  awk '{print $2}' ${dirI}/temp_${out_name}_cleaning_9_ambigious_high_MAF >> ${dirI}/temp_${out_name}_cleaning_9_non_1KG_amb_exclude_pre
sort ${dirI}/temp_${out_name}_cleaning_9_non_1KG_amb_exclude_pre | uniq > ${dirI}/temp_${out_name}_cleaning_9_non_1KG_amb_exclude; wc -l ${dirI}/temp_${out_name}_cleaning_9_non_1KG_amb_exclude
  #SDH 762, 610k_cases 30, 610k_controls 29, QMEGA_omni_cases, AMFS_cases, and AMFS_control 236, WAMHS 12465, HEI_Cases and con 0. OAG1 320 OAG2 347 IBD 12273 QTWIN 7204 EPIGENE 6094 MIA 1408
  echo ""
  plink_1.90 --threads 1 --bfile ${dirI}/temp_${out_name}_cleaning_9 --exclude ${dirI}/temp_${out_name}_cleaning_9_non_1KG_amb_exclude --make-bed --out ${dirI}/${out_name}_aligned_1KG
   #SDH: 939132 variants loaded --exclude: 938370 variant remaining. 938370 variants and 570 people pass filters and QC.
   #610k_cases: 548499 variants loaded from .bim file. --exclude: 548469 variants remaining. 548469 variants and 925 people pass filters and QC.
   #610k_controls: 533116 variants loaded from .bim file. --exclude: 533087 variants remaining. 533087 variants and 3954 people pass filters and QC.
   #QMEGA_omni_cases: 802710 variants loaded from .bim file --exclude: 802474 variants remaining 802474 variants and 694 people pass filters and QC
   #AMFS_cases: 802678 variants loaded from .bim file --exclude: 802442 variants remaining. 802442 variants and 548 people pass filters and QC.
   #AMFS_controls: 802663 variants loaded from .bim file --exclude: 802427 variants remaining 802427 variants and 430 people pass filters and QC
   #WAMHS 915499 variants loaded from .bim file. --exclude: 903034 variants remaining. 903034 variants and 1273 people pass filters and QC.
   #HEI_cases: none removed (no ambigious in the set) 218694 variants and 1217 people pass filters and QC
   #HEI_cont: none removed. 221008 variants and 1221 people pass filters and QC
   #OAG1: 855941 variants loaded from .bim file. --exclude: 855621 variants remaining. 855621 variants and 651 people pass filters and QC.
   #OAG2 707444 variants loaded from .bim file. --exclude: 707097 variants remaining. 707097 variants and 645 people pass filters and QC.
   #IBD 912989 variants loaded from .bim file. --exclude: 900625 variants remaining. 900625 variants and 929 people pass filters and QC.
   #QTWIN  450024 variants loaded --exclude: 442820 variants remaining. 442820 variants and 981 people pass filters and QC.
   #EPIGENE 399153 variants loaded --exclude: 393059 variants remaining. 393059 variants and 784 people pass filters and QC
   #MIA 486958 variants loaded --exclude: 485550 variants remaining. 485550 variants and 1968 people pass filters and QC.
  echo -e "\nDataset ${out_name} appears to be aligned to 1000 genomes succesfully\n"

  rm -f ${dirI}/temp_${out_name}_update_map_bp.txt ${dirI}/temp_${out_name}_cleaning_1.* ${dirI}/temp_${out_name}_cleaning_2.* ${dirI}/temp_${out_name}_cleaning_3.* ${dirI}/temp_${out_name}_cleaning_4* ${dirI}/temp_${out_name}_cleaning_5.* ${dirI}/temp_${out_name}_cleaning_6* ${dirI}/temp_${out_name}_cleaning_7* ${dirI}/temp_${out_name}_cleaning_8.* ${dirI}/temp_${out_name}_cleaning_9* ${dirI}/temp_${out_name}_BED* ${dirI}/temp_${out_name}_cleaning_overlapping* ${dirI}/temp_${out_name}_chr* ${dirI}/temp_${out_name}_bp* ${dirI}/temp_${out_name}_flip* ${dirI}/temp_${out_name}_unmappable* temp_flipped_concordant temp_flipped_chr_* temp_aligned_concordant temp_aligned_chr_* temp_ambigious_aligned_* temp_ambigious_flipped_* temp_non_biallelic_snps* temp_indels_aligned_* temp_aligned_ID_* temp_flipped_ID_* temp_ambigious*MAF_to_high_too_assess temp_ID_flip_list temp_ID_remap temp_indels_ID_aligned_*_allele

done #lift over data and align to 1kg. 

echo -e "\nFinito\n"
exit

+++++++++++++++++++done as far as here for MIA data, not worth going futher until have controls for MIA (needed to get MIA ready for the survival work++++++


#clean up interim files
rm -f temp_2016_cleaning_SDH2.* temp_2016_cleaning_SDH_2.* temp_2016_cleaning_AMFS temp_2016_cleaning_AMFS.* temp_2016_cleaning_AMFS_cases.* temp_2016_cleaning_AMFS_controls.* temp_2016_cleaning_1.* temp_2016_cleaning_2.* temp_2016_cleaning_3.* temp_2016_cleaning_4.* temp_2016_cleaning_QMEGA_omni_cases.* temp_2016_cleaning_non_610k temp_2016_cleaning_QMEGA_610k.* temp_2016_cleaning_QMEGA_610k_cases.* temp_2016_cleaning_QMEGA_610k_controls.* temp_2016_cleaning_WAMHS_cases.* temp_2016_cleaning_IBD_controls.* temp_2016_cleaning_OAGphase*_controls.* temp_2016_cleaning_HEIDELBERG_controls.* temp_2016_cleaning_HEIDELBERG_cases.* ${dirI}/temp_2016_cleaning_EPIGENE_cases.* ${dirI}/temp_2016_cleaning_QTWIN_controls.*

##################################################################################
#4.0 Combine the control sets to test for errors in flipping etc
##################################################################################

#loop to turn on/off
for x in #1
do

  echo -e "\nConverting control sets into case sets for test of flipping and other alignment...\n"
  for i in temp_2016_cleaning_QTWIN_controls_aligned_1KG temp_2016_cleaning_QMEGA_610k_controls_aligned_1KG temp_2016_cleaning_SDH_2_aligned_1KG temp_2016_cleaning_AMFS_controls_aligned_1KG temp_2016_cleaning_HEIDELBERG_controls_aligned_1KG temp_2016_cleaning_OAGphase1_controls_aligned_1KG temp_2016_cleaning_OAGphase2_controls_aligned_1KG temp_2016_cleaning_IBD_controls_aligned_1KG
  do
    awk '{print $1,$2,2}' ${dirI}/${i}.fam > ${dirI}/${i}_dummy_case
    plink_1.90 --threads 1 --bfile ${dirI}/${i} --pheno ${dirI}/${i}_dummy_case --make-bed --out ${dirI}/${i}_dummy_case --silent
  done

  #there are some SNPs that are by designed not aligned (e.g. rs10305752 is missing in 610k dataset, and is an indel, so does not have its alleles fixed where as it is in omni sets) so have a file of SNPs to exclude I made by manually flip/merging. only exclude these for control:control merges
  wc -l ${dirI}/SNPs_to_drop_when_testing_controls_vs_controls.txt

  #still about ~20 multiple positions warnings, and ~1k SNPs that need flipping. After some concerning results realised I was carryign through bad SNPs (extreme HWE results) because controls are now cases - clean properly to find actual problems
  k=temp_2016_cleaning_QMEGA_610k_controls_aligned_1KG_dummy_case
  for i in temp_2016_cleaning_QTWIN_controls_aligned_1KG #temp_2016_cleaning_SDH_2_aligned_1KG temp_2016_cleaning_AMFS_controls_aligned_1KG temp_2016_cleaning_HEIDELBERG_controls_aligned_1KG temp_2016_cleaning_OAGphase1_controls_aligned_1KG temp_2016_cleaning_OAGphase2_controls_aligned_1KG temp_2016_cleaning_IBD_controls_aligned_1KG
  do
    echo -e "\nComparing ${k} and ${i}\n"
    plink_1.90 --threads 1 --bfile ${dirI}/${k} --bmerge ${dirI}/${i} --make-bed --out ${dirI}/${i}_${k} --exclude ${dirI}/SNPs_to_drop_when_testing_controls_vs_controls.txt --silent
    plink_1.90 --threads 1 --bfile ${dirI}/${k} --flip ${dirI}/${i}_${k}-merge.missnp --make-bed --out ${dirI}/${k}_flip --exclude ${dirI}/SNPs_to_drop_when_testing_controls_vs_controls.txt --silent
    plink_1.90 --threads 1 --bfile ${dirI}/${k}_flip --bmerge ${dirI}/${i} --make-bed --out ${dirI}/${i}_${k} --exclude ${dirI}/SNPs_to_drop_when_testing_controls_vs_controls.txt --silent
    plink_1.90 --threads 1 --bfile ${dirI}/${i}_${k} --geno 0.03 --hwe 5e-10 midp include-nonctrl --make-bed --out ${dirI}/${i}_${k}_2 --exclude ${dirI}/SNPs_to_drop_when_testing_controls_vs_controls.txt --silent
    plink_1.90 --threads 1 --bfile ${dirI}/${i}_${k}_2 --mind 0.03 --geno 0.03 --hwe 5e-4 midp --assoc --out ${dirI}/${i}_${k} --silent
    echo -e "${k} vs ${i}"
    awk 'NR==1 || $9<5e-6' ${dirI}/${i}_${k}.assoc
    echo ""
    #clean up
    rm -f ${dirI}/${k}_flip.* ${dirI}/${i}_${k}-merge.missnp ${dirI}/${i}_${k}*
  done #compare QMEGA_610k_controls with SDH and AMFs controls

  k=temp_2016_cleaning_SDH_2_aligned_1KG_dummy_case
  for i in temp_2016_cleaning_QTWIN_controls_aligned_1KG #temp_2016_cleaning_AMFS_controls_aligned_1KG temp_2016_cleaning_HEIDELBERG_controls_aligned_1KG temp_2016_cleaning_OAGphase1_controls_aligned_1KG temp_2016_cleaning_OAGphase2_controls_aligned_1KG temp_2016_cleaning_IBD_controls_aligned_1KG
  do
    plink_1.90 --threads 1 --bfile ${dirI}/${k} --bmerge ${dirI}/${i} --make-bed --out ${dirI}/${i}_${k} --exclude ${dirI}/SNPs_to_drop_when_testing_controls_vs_controls.txt --silent
    echo ""
    plink_1.90 --threads 1 --bfile ${dirI}/${k} --flip ${dirI}/${i}_${k}-merge.missnp --make-bed --out ${dirI}/${k}_flip --exclude ${dirI}/SNPs_to_drop_when_testing_controls_vs_controls.txt --silent
    echo ""
    plink_1.90 --threads 1 --bfile ${dirI}/${k}_flip --bmerge ${dirI}/${i} --make-bed --out ${dirI}/${i}_${k} --exclude ${dirI}/SNPs_to_drop_when_testing_controls_vs_controls.txt --silent
    echo ""
    plink_1.90 --threads 1 --bfile ${dirI}/${i}_${k} --geno 0.03 --hwe 5e-10 midp include-nonctrl --make-bed --out ${dirI}/${i}_${k}_2 --exclude ${dirI}/SNPs_to_drop_when_testing_controls_vs_controls.txt --silent
    plink_1.90 --threads 1 --bfile ${dirI}/${i}_${k}_2 --mind 0.03 --geno 0.03 --hwe 5e-4 midp --assoc --out ${dirI}/${i}_${k} --silent
    echo -e "${k} vs ${i}"
    awk 'NR==1 || $9<5e-6' ${dirI}/${i}_${k}.assoc
    echo ""
    #clean up
    rm -f ${dirI}/${k}_flip.* ${dirI}/${i}_${k}-merge.missnp ${dirI}/${i}_${k}*
  done

  k=temp_2016_cleaning_AMFS_controls_aligned_1KG_dummy_case
  for i in temp_2016_cleaning_QTWIN_controls_aligned_1KG #temp_2016_cleaning_HEIDELBERG_controls_aligned_1KG temp_2016_cleaning_OAGphase1_controls_aligned_1KG temp_2016_cleaning_OAGphase2_controls_aligned_1KG temp_2016_cleaning_IBD_controls_aligned_1KG
  do
    plink_1.90 --threads 1 --bfile ${dirI}/${k} --bmerge ${dirI}/${i} --make-bed --out ${dirI}/${i}_${k} --exclude ${dirI}/SNPs_to_drop_when_testing_controls_vs_controls.txt --silent
    echo ""
    plink_1.90 --threads 1 --bfile ${dirI}/${k} --flip ${dirI}/${i}_${k}-merge.missnp --make-bed --out ${dirI}/${k}_flip --exclude ${dirI}/SNPs_to_drop_when_testing_controls_vs_controls.txt --silent
    echo ""
    plink_1.90 --threads 1 --bfile ${dirI}/${k}_flip --bmerge ${dirI}/${i} --make-bed --out ${dirI}/${i}_${k} --exclude ${dirI}/SNPs_to_drop_when_testing_controls_vs_controls.txt --silent
    echo ""
    plink_1.90 --threads 1 --bfile ${dirI}/${i}_${k} --geno 0.03 --hwe 5e-10 midp include-nonctrl --make-bed --out ${dirI}/${i}_${k}_2 --exclude ${dirI}/SNPs_to_drop_when_testing_controls_vs_controls.txt --silent
    plink_1.90 --threads 1 --bfile ${dirI}/${i}_${k}_2 --mind 0.03 --geno 0.03 --hwe 5e-4 midp --assoc --out ${dirI}/${i}_${k} --silent
    echo -e "${k} vs ${i}"
    awk 'NR==1 || $9<5e-6' ${dirI}/${i}_${k}.assoc
    echo ""
    #clean up
    rm -f ${dirI}/${k}_flip.* ${dirI}/${i}_${k}-merge.missnp ${dirI}/${i}_${k}*
  done

  k=temp_2016_cleaning_HEIDELBERG_controls_aligned_1KG_dummy_case
  for i in temp_2016_cleaning_QTWIN_controls_aligned_1KG #temp_2016_cleaning_OAGphase1_controls_aligned_1KG temp_2016_cleaning_OAGphase2_controls_aligned_1KG temp_2016_cleaning_IBD_controls_aligned_1KG
  do
    plink_1.90 --threads 1 --bfile ${dirI}/${k} --bmerge ${dirI}/${i} --make-bed --out ${dirI}/${i}_${k} --exclude ${dirI}/SNPs_to_drop_when_testing_controls_vs_controls.txt --silent
    echo ""
    plink_1.90 --threads 1 --bfile ${dirI}/${k} --flip ${dirI}/${i}_${k}-merge.missnp --make-bed --out ${dirI}/${k}_flip --exclude ${dirI}/SNPs_to_drop_when_testing_controls_vs_controls.txt --silent
    echo ""
    plink_1.90 --threads 1 --bfile ${dirI}/${k}_flip --bmerge ${dirI}/${i} --make-bed --out ${dirI}/${i}_${k} --exclude ${dirI}/SNPs_to_drop_when_testing_controls_vs_controls.txt --silent
    echo ""
    plink_1.90 --threads 1 --bfile ${dirI}/${i}_${k} --geno 0.03 --hwe 5e-10 midp include-nonctrl --make-bed --out ${dirI}/${i}_${k}_2 --exclude ${dirI}/SNPs_to_drop_when_testing_controls_vs_controls.txt --silent
    plink_1.90 --threads 1 --bfile ${dirI}/${i}_${k}_2 --mind 0.03 --geno 0.03 --hwe 5e-4 midp --assoc --out ${dirI}/${i}_${k} --silent
    echo -e "${k} vs ${i}"
    awk 'NR==1 || $9<5e-6' ${dirI}/${i}_${k}.assoc
    echo ""
    #clean up
    rm -f ${dirI}/${k}_flip.* ${dirI}/${i}_${k}-merge.missnp ${dirI}/${i}_${k}*
  done

  k=temp_2016_cleaning_OAGphase1_controls_aligned_1KG_dummy_case
  for i in temp_2016_cleaning_QTWIN_controls_aligned_1KG #temp_2016_cleaning_OAGphase2_controls_aligned_1KG temp_2016_cleaning_IBD_controls_aligned_1KG
  do
    plink_1.90 --threads 1 --bfile ${dirI}/${k} --bmerge ${dirI}/${i} --make-bed --out ${dirI}/${i}_${k} --exclude ${dirI}/SNPs_to_drop_when_testing_controls_vs_controls.txt --silent
    echo ""
    plink_1.90 --threads 1 --bfile ${dirI}/${k} --flip ${dirI}/${i}_${k}-merge.missnp --make-bed --out ${dirI}/${k}_flip --exclude ${dirI}/SNPs_to_drop_when_testing_controls_vs_controls.txt --silent
    echo ""
    plink_1.90 --threads 1 --bfile ${dirI}/${k}_flip --bmerge ${dirI}/${i} --make-bed --out ${dirI}/${i}_${k} --exclude ${dirI}/SNPs_to_drop_when_testing_controls_vs_controls.txt --silent
    echo ""
    plink_1.90 --threads 1 --bfile ${dirI}/${i}_${k} --geno 0.03 --hwe 5e-10 midp include-nonctrl --make-bed --out ${dirI}/${i}_${k}_2 --exclude ${dirI}/SNPs_to_drop_when_testing_controls_vs_controls.txt --silent
    plink_1.90 --threads 1 --bfile ${dirI}/${i}_${k}_2 --mind 0.03 --geno 0.03 --hwe 5e-4 midp --assoc --out ${dirI}/${i}_${k} --silent
    echo -e "${k} vs ${i}"
    awk 'NR==1 || $9<5e-6' ${dirI}/${i}_${k}.assoc
    echo ""
    #clean up
    rm -f ${dirI}/${k}_flip.* ${dirI}/${i}_${k}-merge.missnp ${dirI}/${i}_${k}*
  done

  k=temp_2016_cleaning_OAGphase2_controls_aligned_1KG_dummy_case
  for i in temp_2016_cleaning_QTWIN_controls_aligned_1KG #temp_2016_cleaning_IBD_controls_aligned_1KG
  do
    plink_1.90 --threads 1 --bfile ${dirI}/${k} --bmerge ${dirI}/${i} --make-bed --out ${dirI}/${i}_${k} --exclude ${dirI}/SNPs_to_drop_when_testing_controls_vs_controls.txt --silent
    echo ""
    plink_1.90 --threads 1 --bfile ${dirI}/${k} --flip ${dirI}/${i}_${k}-merge.missnp --make-bed --out ${dirI}/${k}_flip --exclude ${dirI}/SNPs_to_drop_when_testing_controls_vs_controls.txt --silent
    echo ""
    plink_1.90 --threads 1 --bfile ${dirI}/${k}_flip --bmerge ${dirI}/${i} --make-bed --out ${dirI}/${i}_${k} --exclude ${dirI}/SNPs_to_drop_when_testing_controls_vs_controls.txt --silent
    echo ""
    plink_1.90 --threads 1 --bfile ${dirI}/${i}_${k} --geno 0.03 --hwe 5e-10 midp include-nonctrl --make-bed --out ${dirI}/${i}_${k}_2 --exclude ${dirI}/SNPs_to_drop_when_testing_controls_vs_controls.txt --silent
    plink_1.90 --threads 1 --bfile ${dirI}/${i}_${k}_2 --mind 0.03 --geno 0.03 --hwe 5e-4 midp --assoc --out ${dirI}/${i}_${k} --silent
    echo -e "${k} vs ${i}"
    awk 'NR==1 || $9<5e-6' ${dirI}/${i}_${k}.assoc
    echo ""
    #clean up
    rm -f ${dirI}/${k}_flip.* ${dirI}/${i}_${k}-merge.missnp ${dirI}/${i}_${k}*
  done

  #clean up
  rm -f ${dirI}/temp_2016_cleaning_SDH_2_aligned_1KG_dummy_case* ${dirI}/temp_2016_cleaning_QMEGA_610k_controls_aligned_1KG_dummy_case* ${dirI}/temp_2016_cleaning_AMFS_controls_aligned_1KG_dummy_case* ${dirI}/temp_2016_cleaning_OAGphase*_controls_aligned_1KG_dummy_case* ${dirI}/temp_2016_cleaning_IBD_controls_aligned_1KG_dummy_case* ${dirI}/temp_2016_cleaning_HEIDELBERG_controls_aligned_1KG_dummy_case* ${dirI}/temp_2016_cleaning_IBD_controls_aligned_1KG_temp_2016_cleaning_QMEGA_610k_controls_aligned_1KG_dummy_case.* ${dirI}/temp_2016_cleaning_QTWIN_controls_aligned_1KG_dummy_case.* ${dirI}/temp_2016_cleaning_QTWIN_controls_aligned_1KG_dummy_case
done #loop to turn control comparison on/off

  #QTWIN controls vs QMEGA_610k
  #CHR         SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR 
  # 1   rs1482658  199204879    G    0.407   0.3504    T        21.04    4.492e-06        1.273 
  # 1   rs6429457  244422816    G   0.0291  0.05618    A        34.57    4.106e-09       0.5035  << touch worrying. Or an ENDO hit? Doesn't seem to be. Not near wint 4
  # 5   rs7711172  127257406    A 0.009742  0.02345    G        23.97     9.76e-07       0.4098 
  # 7   rs7790316   80933327    T   0.1037  0.06932    C        21.26    4.005e-06        1.554 
  # 8  rs10957495   70445363    T  0.01366  0.03007    C        25.45     4.55e-07       0.4467 
  #10  rs16925541   70405539    G  0.04169  0.06946    A        26.87    2.171e-07       0.5828 
  #12   rs7312252   50744171    T   0.3447   0.4001    C        21.06    4.447e-06       0.7886 
  #12  rs12303082   50754563    T   0.3447   0.4001    G        21.07    4.424e-06       0.7886 
  #17   rs3826543    8132763    T   0.1764   0.2217    C        21.32    3.886e-06       0.7521 
  #18   rs1370449   28762788    T  0.06815   0.1005    C        23.82     1.06e-06       0.6545 
  #23   rs3012642   71628912    G 0.007922  0.02398    A        30.57    3.219e-08        0.325  << touch worrying. Or an ENDO hit?
  #23   rs2369215   71817252    G 0.008203  0.02398    T        28.91    7.601e-08       0.3366  << touch worrying. Or and ENDO hit? 

    #plink_1.90 --threads 1 --bfile temp_2016_cleaning_QTWIN_controls_aligned_1KG --hardy
    #grep "rs6429457\|rs3012642\|rs2369215" plink.hwe
       #nothing big P 0.03. Note the two chrX ones seem to be in LD. 
    #A bit concerning

  #QTWIN controls vs SDH
  # CHR             SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR 
  #   5       rs7732691  157957642    A   0.4395   0.3562    G        21.05    4.471e-06        1.417 
  #  17      rs11078738    8167600    T    0.217   0.2946    C        22.18    2.485e-06       0.6638

  #QTWIN AMFS
  #CHR              SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR 
  # 7        rs7790316   80933327    T   0.1326  0.06932    C        29.62    5.264e-08        2.052 

  #QTWIN HEI - many, but not as many as others, probably only due to limited overlap

  #QTWIN OAG1
  # CHR          SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR 
  # 2   rs34664134  241812461    A  0.02742 0.0005102    G        50.84    1.003e-12        55.23 <<seems to be borked in OAG1, keeps coming up 
  # 4   rs13142179  185781717    T   0.2289   0.1641    C        21.33    3.863e-06        1.512 
  # 6    rs1365711   93127527    T  0.05776 0.005607    C        81.93     1.41e-19        10.87 <keeps coming up, seem to be borked in OAG1; dropped in the last imp
  # 9     rs564398   22029547    C   0.3282   0.4123    T         23.5    1.249e-06       0.6963  <OG locus
  # 9    rs7865618   22031005    G   0.3405    0.422    A        21.82        3e-06       0.7072  <OAAG locus
  #10   rs12411414   51832220    G   0.1359   0.2074    A        26.61    2.483e-07        0.601 
  #19     rs862709   57868336    A  0.01462        0    G        28.73    8.337e-08           NA 
  #23     rs980136  115457731    T   0.2896   0.3798    C        21.73    3.142e-06       0.6658 

  #QTWIN OAG2
 #CHR          SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR 
 #  1   rs11800181  102120197    T   0.1166  0.04434    G        58.69    1.848e-14        2.844  <<icomes up in other matches, borked in OAG. Dropped last imp
 #  7    rs4236302   21211412    A   0.2025   0.2752    G        22.11    2.575e-06       0.6686 
 #  9    rs1063192   22003367    G   0.3349   0.4256    A        26.76    2.305e-07       0.6796 
 #  9     rs564398   22029547    C   0.3128   0.4123    T        32.85    9.952e-09       0.6487  <<OAG  locus
 #  9    rs1412832   22077543    C   0.2214   0.3048    T        27.16    1.875e-07       0.6487 
 # 11   rs10769708    6719118    G   0.3495   0.2737    A        21.03    4.526e-06        1.426 
 # 11    rs1395558    6727468    A   0.3879   0.3078    G        22.15     2.52e-06        1.425 
 # 17   rs12150284   10031090    T   0.3022   0.3853    C        23.51    1.245e-06       0.6908 
 # 17    rs9916496   62750365    A   0.1995 0.0005097    G        419.9    2.602e-93        488.7 <<not show up before, dropped in last imp. 
 # 19     rs862709   57868336    A   0.0234        0    G        46.15    1.094e-11           NA 
      #zgrep rs9916496 ALL_1000G_phase1integrated_v3_ALL_chr_legends_with_chromosome.txt.gz 
      #17 rs9916496 62750365 G A 0.1626 0.0166 0.0000 0.0000 0.1626 0.0166 0.0000 0.0000 <<looks to be incorrect in the "A" which is OAG2

  #SDH vs 610k controls - with cleaning. 296779 variants and 4523 people pass filters and QC. Bonf correction is 1.69e-7, but valid to use normal GWS, so nothing more than chance I think here
  #CHR          SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR 
  #   5    rs6882579  123147064    C  0.04873  0.08246    A        22.58    2.019e-06       0.5701 
  #   6    rs9403856  147759534    T   0.1127   0.1605    C        21.75    3.104e-06       0.6644 
  #   9     rs885290   16203309    A   0.1379   0.1907    G        22.45    2.162e-06       0.6787 
  #  16     rs276626   13156649    G   0.1611   0.2226    A        26.71    2.358e-07       0.6704 

  #610k vs AMFS
  #CHR          SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR
  # 3   rs17035012   36440897    C  0.07806   0.1279    T        25.33    4.827e-07       0.5773
  # 8    rs3110043  110093871    G     0.33   0.4077    A        20.88     4.89e-06       0.7154
  #11     rs978769  110888484    T   0.3079   0.3895    C        23.94    9.946e-07       0.6971
  #12   rs11059421  128436536    G  0.03759  0.07377    A         25.6    4.208e-07       0.4903

  #610k_controls as cases versus HEIDELBERG controls
  #so many GWS hits ~ 520 or So SNP with P < 5e-6, many GWS. Not a good ethnic match for the Aus set but still....
  #As not before This set is so, so weird. A chr24 is weird - 0.05MAF in 610k, 20% in HEI (but again maybe you would expect big freq diff on Y chr?) and a few fairly worrying ones 
  #stopped looking from chr5 to 1, not really useful as I am only justifying my bias (that I think it is stratficiation)
  # 6    rs908024     463753    T   0.2693   0.3357    C        40.29    2.191e-10       0.7294 <near IRF4, but only LD 0.01 with IRF4 SNP. Still suggests strat
  # 6 rs114648420   32209027    C   0.2122   0.3184    T        116.2    4.315e-27       0.5767 <<<<HLA region so may just be stratification
  # 6 rs114205116   32705276    T   0.3794     0.45    C        38.78    4.737e-10       0.7473 <<<<HLA
  # 6 rs114895229   32730940    A   0.2733   0.3643    C        74.28    6.768e-18       0.6563 <<<<<HLA
  # 6 rs115084029   32763514    A   0.3867   0.2902    G        75.02    4.668e-18        1.542 <<HLA
  # 6 rs116812508   32941910    T    0.258   0.3318    G        50.85    9.966e-13       0.7003 <<HLA
  # 9   rs7021746   94176261    A   0.2819    0.219    G         37.6    8.675e-10          1.4
  #9    rs657152  136139265    A   0.3367   0.4066    C        39.81    2.795e-10        0.741
  #13  rs12427790   78055855    T  0.01152  0.04516    C        110.4    7.887e-26       0.2465
  #16  rs12923427   17575065    T    0.232   0.1684    C        44.28    2.846e-11        1.491
  #19    rs504963   49208865    G     0.45   0.5594    A        89.48    3.093e-21       0.6444

  #610k controls vs OAGphase1 < I think these are just  the disease group
  # CHR          SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR
  # 1    rs2797261  215903837    G  0.01025  0.03198    T        39.74    2.902e-10       0.3135
  # 5      rs27183   58588629    A   0.4344   0.5248    C        32.67    1.089e-08       0.6955
  # 7   rs10228341   18920499    G   0.0143  0.04064    A        37.87    7.571e-10       0.3424
  # 7   rs10271398  140868640    C  0.05716   0.1065    T        39.06    4.104e-10       0.5088
  # 8    rs7845936  132352151    A  0.04749  0.09208    G         37.9    7.464e-10       0.4916
  #  8    rs7829952   15809821    T  0.02593   0.1343    C        301.7    1.414e-67       0.1716 <<<<<failed in 610k? Checked the 610k melanoam data and it is 0.03 percent or so in both sets. ~4% n 1kg, failed in OAG? Dropped in the WAMGS/OAG GWAS set. 
  # 9    rs1063192   22003367    G   0.4398   0.3546    A        33.09    8.799e-09        1.429  <<<OAG GWAS hit
  # 9    rs7049105   22028801    G   0.4722   0.5588    A        33.39     7.56e-09       0.7062
  # 9   rs10120688   22056499    A   0.4843   0.5785    G        39.58    3.145e-10       0.6844
  # 9    rs4977756   22068652    G   0.4008   0.3069    A        41.41    1.232e-10         1.51 <<literally a top OAG SNP.
  # 12     rs345684  128061788    T   0.1054   0.1654    G        33.17    8.456e-09       0.5945
  # 14    rs7142344   79250144    T   0.1001   0.1705    C        46.97    7.213e-12       0.5409
  # 17    rs4246444   80038952    T   0.2484   0.3422    G        44.04    3.219e-11       0.6353
  # 19    rs8106771   16089925    T   0.1117   0.1736    C        35.26    2.879e-09       0.5983

  #610k_controls OAG2
  #   1    rs2797261  215903837    G  0.01025  0.03198    T        39.74    2.902e-10       0.3135
  #   4   rs12639641  102040091    C  0.03453  0.06797    T        32.71    1.069e-08       0.4904
  # 6    rs1321204  101468591    G  0.04328  0.08346    A        38.03    6.967e-10       0.4968
  # 7    rs7805156   44730284    T  0.01455  0.03738    C        33.03    9.079e-09       0.3803
  #  7    rs7810984  141761438    T  0.04152  0.08022    C        36.98    1.192e-09       0.4967
  #  8    rs4737534   59963737    G  0.09284   0.1433    A        31.28    2.231e-08       0.6118
  #  9    rs1063192   22003367    G   0.4398   0.3349    A        49.53    1.956e-12        1.559
  #   9    rs4977756   22068652    G   0.4008    0.305    A        42.63    6.625e-11        1.524 <<again OAG hit
  #    9    rs1412832   22077543    C   0.3063   0.2214    T        38.13    6.614e-10        1.552
  #  10    rs3813865  135339244    C  0.02012   0.0466    G        32.97    9.378e-09       0.4201
  #13    rs2992385   29843289    G  0.01202  0.03193    A        30.03    4.249e-08       0.3688
  #13   rs16962524  104137962    G  0.01935  0.04446    A        31.14    2.399e-08       0.4241
  #14    rs1992302   75013174    C  0.02125  0.05218    T        42.41    7.392e-11       0.3944
  #14    rs2359141   75015894    A  0.02112   0.0514    G        40.98    1.536e-10       0.3982
  #  22    rs8137687   37943299    A  0.05148  0.09034    G        31.08    2.471e-08       0.5465

  #610k vs IBD
  #   4  rs10013026  104532648    T  0.05162   0.0863    C        32.95    9.468e-09       0.5763
  #   7   rs2717808  148213724    A   0.3584   0.0625    C        565.4   5.621e-125         8.38 <<<<<<<~35% is the 1kg level, seems to have failed in IBD.
  #17  rs11653068    3408576    T  0.01392   0.0329    G        31.42     2.08e-08       0.414
  #all three are like that in the intial set, not due to renaminimg on my side etc (as in MAF etc same in IBD)

  #grep "rs10013026\|rs2717808\|rs11653068" /working/lab_stuartma/scratch/matthewL/melanoma_meta_analysis/cleaned/WAMHS_GLUACOMA_IBD_BEACON_final_clean_IBD_PCA_outliers_removed.bim
  #4	rs10013026	0	104532648	T	C
  #17	rs11653068	0	3408576	T	G
  #so 7 was dropped for some reason - maybe not in all arrays. While the chr4 SNP is extreme in CIDR_WAMHS, it is even more so in the opposiute direction in HEIDELBERG, and nearly as extreme in FRENCH. I think I am obsessing over nothing.


  #SHD vs AMFS
  # CHR              SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR
  # 2        rs4144323   83841076    A   0.3518    0.457    G        22.67    1.922e-06       0.6448

  #SDH versus HEI
  #also lost of signal, many extreme chr24 hits so again I am thinking strat plus not using the right test there. Cluster of HLA SNPs as well near GWS 
  #of note: rs17774013; seeems to have failed in HEI and shows up in every matching. rs17774013 was dropped in the final melanoma set
  #   1   rs1801274  161479745    G   0.5181   0.4393    A        46.27    1.028e-11        1.372
  #   2  rs17774013   47067214    T   0.2728  0.03654    C        436.4     6.65e-97        9.893 <<<<<< gives similar results in the HEI set sent to me (as in 27% in cases, and 4% in controls, so not something I was given. Not in the cleaned german set so why not? Oh I noted this exact results, and dropped it from the dataset prior to imputation. Need to work out a formal way to flag these SNPs for removal.
  #    2  rs12465220   29559519    C   0.4658   0.5328    T        33.49    7.169e-09       0.7646
  # 2  rs13430738   98839367    G   0.1818   0.2426    A        43.81     3.62e-11       0.6934
  # 2   rs9287442  136522710    A   0.1912   0.2578    G        50.43    1.237e-12       0.6805
  # 2  rs10188066  136539513    G   0.1897   0.2566    A        50.96    9.409e-13       0.6786
  # 2  rs12472293  136648077    C     0.13   0.1828    T        42.48    7.137e-11       0.6682
  # 2   rs4954391  136883823    C   0.2606   0.1996    T        37.42      9.5e-10        1.414
  # 2    rs953387  136907170    C   0.2768   0.3529    A        51.87    5.933e-13        0.702
  
  # again a bunch of HLA hits... and many more highly sig hits. Almost wonder if better off dropping HEIL

  #

  #SDH vs OAG1 - once you ignore the maf 0 SNPs nothin too concerning other than the CDKN2A hit for OAG
  # CHR         SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR
  # 2  rs34664134  241812461    A        0  0.02742    G        31.43    2.062e-08            0  <<<how are these 0 MAF not being removed earlier? 
  # 4   rs4516798  190074821    A        0  0.03285    T        37.64    8.499e-10            0
  # 5  rs10044029   36854542    T        0  0.03296    C        38.23    6.288e-10            0
  # 6 rs147564756   32890939    A        0  0.02701    C        31.18     2.35e-08            0
  # 7  rs10248706  119728147    T 0.002641  0.02989    C        26.55    2.571e-07      0.08595
  # 8   rs4132448   32453838    T        0  0.03074    C        35.62      2.4e-09            0
  # 9   rs1063192   22003367    G   0.4491   0.3546    A        22.63    1.967e-06        1.484 <<OAG CDKN2A hit
  # 9    rs564398   22029547    C   0.4272   0.3282    T        25.39    4.677e-07        1.527 
  # 9   rs7865618   22031005    G   0.4377   0.3405    A        24.19    8.714e-07        1.508
  # 9    rs634537   22032152    G   0.4272   0.3277    T        25.68    4.032e-07         1.53
  # 9   rs4977756   22068652    G   0.4096   0.3069    A        27.99    1.216e-07        1.567
  #11   rs2970332   14360435    G   0.1904    0.274    A        23.27    1.406e-06       0.6228
  #13   rs2298086   47308109    A        0  0.02568    G        29.57    5.387e-08            0
  #14   rs8015216   59295696    G   0.3737   0.4662    A        21.32    3.892e-06       0.6831
  #17    rs936397   73840774    A        0  0.04383    G        51.04    9.037e-13            0
  #19   rs2997212   21512188    A        0  0.02283    G         26.3    2.919e-07            0
  #20   rs6101092   59408388    A        0  0.02922    G        33.83     6.03e-09            0

  #SDH vs OAG2
   #1  rs11800181  102120197    T  0.05088   0.1166    G        32.75    1.047e-08       0.4062
   #6 rs115606363   31074128    T  0.05614   0.3503    C        310.3    1.919e-69       0.1103 <<HLA?
   #9   rs1063192   22003367    G   0.4491   0.3349    A        33.08     8.83e-09        1.619 OAG CDKN2 hit
   #9    rs564398   22029547    C   0.4272   0.3128    T        34.01    5.488e-09        1.639
   #9   rs4977756   22068652    G   0.4096    0.305    A        28.89    7.652e-08        1.581
   #9   rs9632885   22072638    G   0.5272   0.4306    A        22.58    2.012e-06        1.475
   #10   rs8192780  135354125    T   0.1085   0.1968    G         35.2    2.981e-09       0.4969
  #otherwise nothing too concerning. rs11800181 and rs8192780 dropped from WAMHS/OAG. rs11800181 1kg matches SDH, rs8192780 1kg matches OAG. rs8192780 dropped in WAMHS set (where its freq might have been an issue. I wonder why. 


   #SDH vs IBD
   #CHR             SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR
   #2      rs10496917  142967224    T   0.1602   0.1031    G        20.97    4.676e-06         1.66
   #2      rs67254537  179395555    A  0.02018        0    C        37.69    8.272e-10           NA
   #4       rs1489001  117245430    C  0.08187   0.1387    T        21.96    2.777e-06       0.5538
   #6       rs9379861   26370616    C  0.01579        0    G         29.2    6.536e-08           NA
   #6       rs2534657   31472459    T   0.1623   0.1035    C        21.66    3.254e-06        1.678
   #6       rs4707557   90362783    C  0.02368        0    A        42.65    6.558e-11           NA
   #6       rs2504920  160897872    C   0.5202   0.4277    A        24.27    8.379e-07         1.45
   #10        rs293326   53184548    C    0.186   0.1246    T        21.09    4.384e-06        1.605
   #10       rs8192780  135354125    T   0.1085   0.1969    G        39.91    2.654e-10       0.4967
   #11      rs11537993   66109033    C   0.3921   0.3096    T        21.38    3.761e-06        1.438
   #12        rs748125     632922    G  0.03421  0.07767    A        23.23    1.434e-06       0.4206
   #18       rs1521788   32026262    G   0.2947   0.3877    A        26.72    2.352e-07         0.66
   #21       rs2406176   19719426    T   0.1408   0.2438    G        45.86    1.268e-11       0.5085
   #22       rs8135345   35771689    A   0.1592    0.322    G        95.88    1.219e-22       0.3987

   #Ignoring the 0 maf SNPs, rs8192780 - again looks to be wrong in IBD, leaves rs2406176 (~20% in 1kg) and rs8135345 (32% in 1kg). So the latter might be wrong in SDG, and the formaer might just be chance - neither really match the 1kg freq

  #AMFs vs HEI controls - again many <1e-6
   #CHR         SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR
   #HLA signal, though that as extreme as in other sets
   #Only a Y chrom is GWS - maybe just a power issue (as in if AMFS was larger would be as bad as the other set. Still nothing new leaps out.

  #AMFS vs OAG
  #nothing GWAS but a decent signal for SNP all around 16.5Mb - a real OAG hit?  Or just a real chance difference in pops and the SNPs are in LD. 
  # CHR         SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR
  # 1  rs12119764   16552772    G   0.3942   0.2934    A        23.67    1.143e-06        1.567
  # 1   rs2863845   16561179    T   0.4174   0.3065    C        28.02    1.199e-07        1.622
  # 1  rs17369676   16575030    C   0.4418   0.3351    T        24.64    6.923e-07         1.57
  # 1   rs7536390   16576836    G   0.4442   0.3408    A        23.45    1.282e-06        1.546
  # 1   rs4661729   16579070    T   0.4442   0.3387    C        24.44    7.667e-07         1.56
  # 1  rs10907222   16594055    A   0.4452   0.3402    C        24.14    8.943e-07        1.556
  # 1  rs10489598   16619350    A    0.443   0.3372    G        24.65    6.882e-07        1.564
  # 1  rs11804252   16642849    T   0.4442     0.34    C        23.81    1.063e-06        1.551
  # 1  rs10907225   16720112    C    0.436   0.3387    T        20.88    4.883e-06         1.51
  # 8   rs2439322   32520472    G    0.436   0.3379    A        21.22    4.087e-06        1.515
  # 9   rs7040587   16256964    A   0.2988   0.2089    G        22.66    1.931e-06        1.614
  #14   rs7147840   21697858    C    0.157   0.2488    T        26.14    3.176e-07       0.5621
  #18   rs1038670   45898195    A   0.1744   0.1046    G        21.95    2.793e-06        1.808


  #AMFS OAG2
  #   6 rs115606363   31074128    T  0.06163   0.3503    C        237.3    1.514e-53       0.1218
  #  rest are >GWS but do include CDKN2A SNPs

  #AMFS IBD - ignoring the 0 maf SNP, rest aboue whay youtm ight expect by chance. 
  #CHR              SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR
  # 2        rs6712833   12673361    A   0.4266   0.3333    G        22.06    2.646e-06        1.488
  # 2        rs2627765   55883317    G   0.3858   0.4829    T        22.17    2.493e-06       0.6726
  # 2       rs67254537  179395555    A  0.01628        0    C        30.34     3.63e-08           NA
  # 3        rs6769457  180365956    A   0.2698   0.1893    G        22.54    2.056e-06        1.582
  # 6        rs4707557   90362783    C  0.01628        0    A         29.2     6.54e-08           NA
  #11        rs1585962   38664557    T  0.05711   0.1127    C        21.08    4.407e-06       0.4767
  #11        rs1383741   38670343    T  0.05581   0.1112    C        21.25    4.037e-06       0.4723
  #18         rs568388   33829077    T   0.5291   0.4185    C        28.99    7.284e-08        1.561

  #HEI as cases versus OAG
  #again multiple SNPs, including rs17774013 showing a ridiculous signal. Strong HLA sig (GWS, 
  #   6   rs1365711   93127527    T 0.002051  0.05776    C        124.8    5.577e-29      0.03353
  #  22  rs11913963   38028942    T        0  0.04269    C        105.6    9.081e-25            0 <<0 though

  #HEI versus the rest as you might expect

  #OAG1 vs IBD - actually more interesrting
  #CHR         SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR
  # 1   rs4658359   89847372    G  0.03513   0.2657    T        272.6    3.105e-61       0.1006
  # 2  rs34664134  241812461    A  0.02742 0.003236    G        34.05    5.365e-09        8.683
  #   6  rs13195402   26463575    T 0.0007862  0.08054    G        104.8    1.335e-24     0.008982
  # 6   rs4707557   90362783    C  0.02074        0    A        37.32    1.002e-09           NA
  # 6   rs1365711   93127527    T  0.05776   0.0054    C        78.97    6.296e-19        11.29
  # 7   rs4236206   65070538    A   0.2686   0.3699    G        34.41    4.456e-09       0.6257
  # 9   rs2043664  107694245    A   0.3562   0.2643    G        30.58    3.201e-08         1.54
  #21   rs2406176   19719426    T   0.1369   0.2438    G         54.6    1.474e-13       0.4921
  #22   rs9612352   23763969    T    0.306   0.4186    G        41.29     1.31e-10       0.6126
  #22   rs8135345   35771689    A   0.1985    0.322    G        59.03    1.554e-14       0.5213

#temp_2016_cleaning_OAGphase2_controls_aligned_1KG_dummy_case vs temp_2016_cleaning_IBD_controls_aligned_1KG
# CHR         SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR
#   1   rs3795501   36180972    T  0.06555  0.02535    C        29.74    4.945e-08        2.697
#   1  rs11800181  102120197    T   0.1166  0.04746    G         50.7    1.075e-12        2.648
#   5  rs10053872  125532171    T  0.02273 0.001618    G        33.28     7.97e-09        14.35
#   6   rs6930590   50364398    G  0.01869        0    T        34.92    3.433e-09           NA
#   6   rs9498740  102341430    A  0.01194   0.1348    G        145.1    2.051e-33      0.07755
#   6   rs4870310  155033087    C   0.2758   0.4277    T        73.17     1.19e-17       0.5095
#   7   rs4723400   35277426    A  0.03271        0    G        58.44    2.094e-14           NA
#   7   rs6943064  148238976    C   0.0182  0.00054    T        30.64    3.101e-08        34.31
#   8   rs2703260    5775319    A   0.1332   0.3781    G        217.9     2.59e-49       0.2528
#  10   rs7070439   98230489    A   0.1537   0.2508    G        42.44    7.287e-11       0.5426
#  10 rs145397190  121196334    T   0.1147        0    C          223    1.961e-50           NA
#  10   rs7081582  130677124    C   0.2023   0.2988    T        36.68    1.392e-09       0.5953
#  21   rs2406176   19719426    T   0.1285   0.2438    G        63.84    1.352e-15       0.4574



  #All these exrtreme SNPs dropped from the melanoma set. Deliberate or cleaned? 
  #rs4658359 in the seperate input sets.. though I bet the HWE is weird if you combine a SNP with 0.035 and 0.27 MAF - so that might have done the job 


  #AMFS vs 610k - with cleaning 298019 variants and 4383 people pass filters and QC.
  # CHR               SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR 
  #   3        rs17035012   36440897    C  0.07806   0.1279    T        25.33    4.827e-07       0.5773 
  #   8         rs3110043  110093871    G     0.33   0.4077    A        20.88     4.89e-06       0.7154 
  #  11          rs978769  110888484    T   0.3079   0.3895    C        23.94    9.946e-07       0.6971 
  #  12        rs11059421  128436536    G  0.03759  0.07377    A         25.6    4.208e-07       0.4903 

  #AMFS vs SDH  787121 variants and 1000 people pass filters and QC.  (was getting some 1e-22 SNPs that were hwe failures in the 'case' SDH group but not being filtered without --hwe 5e-10 midp include-nonctrl

  #CHR              SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR 
  #  2        rs4144323   83841076    A   0.3518    0.457    G        22.67    1.922e-06       0.6448 

  #12/07/2016 after a couple of false starts and IDing additonal SNPs that needed removing, Nothing here leaps out as worrying

  #OAG1 vs OAG2
 #  CHR         SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR
 #  1  rs11800181  102120197    T  0.04378   0.1166    G        45.86    1.268e-11       0.3469
 #  2   rs4148189   44047530    T  0.08679   0.1503    C           25    5.723e-07       0.5372
 #  5  rs10062847  178179729    C   0.2476   0.3315    T        21.66     3.26e-06       0.6636
 #  6 rs115606363   31074128    T   0.0722   0.3503    C        298.1    8.479e-67       0.1443
 #  6 rs114212145   32671247    T  0.09769   0.1651    C        25.76    3.869e-07       0.5475
 #  7   rs4236206   65070538    A   0.2686   0.3592    G        23.89    1.021e-06       0.6553
 #  9   rs7047933    4659111    T   0.5231    0.429    C        22.89    1.712e-06         1.46
 # 11   rs7111524   25841229    G   0.1828   0.2586    T        21.24    4.054e-06        0.641
 # 14  rs10438117   20398326    T 0.0008117  0.02191    C        24.45    7.632e-07      0.03627
 # 14    rs893298   52690059    T   0.2779   0.3644    G        21.74    3.125e-06       0.6711
 # 18   rs4800685   23501282    T   0.3891   0.3014    C        21.93    2.821e-06        1.476
 # 21   rs9978024   28713216    T   0.3217   0.4165    C        24.61    7.015e-07       0.6644



  #27/07/2016 across the more diverse range of SNPs there is more to be concerned about; HEI, IBD and OAG sets seem to have a number of concerning SNPs. May seem to have been removed during cleaning of the mel + control sets, or in the case of the few extreme HEI SNPs IDd during QC and removed prior to imputation. Still, a bit worrying. 
   #need to come back and see if any of these SNPs are still here post other QC

##################################################################################
#5.0 Minimal cleaning
##################################################################################
#given the unusual SNPs observed I wonder if I should perform so basic cleaning within the seperate case control sets as some extreme examples might be otherwise missed (e.g. a SNP with MAF 0 in cases and 6% in cases might pass a 0.001 maf filter in the combined set
#So minimal cleaning should still retain SNPs with freq higher than 0.001, filter on the most extreme hwe (5e-10, so treat them all as cases) and normal missingness (0.03) because again a SNP might be missing in 5% controls and 0% cases and pass in the combined set
for x in #1
do
  for i in #temp_2016_cleaning_AMFS_cases_aligned_1KG temp_2016_cleaning_AMFS_controls_aligned_1KG temp_2016_cleaning_HEIDELBERG_cases_aligned_1KG temp_2016_cleaning_HEIDELBERG_controls_aligned_1KG temp_2016_cleaning_IBD_controls_aligned_1KG temp_2016_cleaning_OAGphase1_controls_aligned_1KG temp_2016_cleaning_OAGphase2_controls_aligned_1KG temp_2016_cleaning_QMEGA_610k_cases_aligned_1KG temp_2016_cleaning_QMEGA_omni_cases_aligned_1KG temp_2016_cleaning_SDH_2_aligned_1KG temp_2016_cleaning_WAMHS_cases_aligned_1KG
  do
   echo -e  "\nCleaning ${i} to minimal settings - maf 0.001, geno 0.03, mind 0.03, --hwe 5e-10 midp include-nonctrl\n"
   plink_1.90 --threads 1 --bfile ${dirI}/${i} --geno 0.03 --mind 0.03 --maf 0.001 --hwe 5e-10 midp include-nonctrl --make-bed --out ${dirI}/${i}_min_clean
  done #minimal cleaning

    #AMFS_cases: 802442 variants loaded. 548 people (205 males, 343 females) loaded. 0 people removed due to --mind 35 variants removed due to --hwe. 0 variants removed due to minor allele threshold(s) 2877 variants removed due to missing genotype data (--geno). 799530 variants and 548 people pass filters and QC.
    #AMFS_controls: 802427 variants loaded 30 people (179 males, 251 females) loaded 0 people removed due to missing genotype data  2331 variants removed due to missing genotype data (--geno) --hwe: 27 variants removed  1 variant removed due to minor allele threshold(s) 800068 variants and 430 people pass filters and QC
    #HEI_cases: 218694 variants loaded 1217 people (706 males, 511 females) loaded 355 variants removed due to missing genotype data (--geno) --hwe: 86 variants removed due to Hardy-Weinberg exact test 3760 variants removed due to minor allele threshold(s) 214493 variants and 1217 people pass filters and QC 
    #HEI_con: 221008 variants loaded  221 people (610 males, 611 females) loaded 1 person removed due to missing genotype data (--mind) 532 variants removed due to missing genotype data (--geno) --hwe: 111 variants removed due to Hardy-Weinberg 4501 variants removed due to minor allele threshold(s) 215864 variants and 1220 people pass filters and QC
    #IBD 900625 variants loaded 929 people (419 males, 510 females); 2 people removed due to --mind; 5079 variants removed due to --geno; --hwe: 504 variants removed; 197458 variants removed due to --maf; 697584 variants and 927 people pass filters and QC
    #OAG1 855621 variants loaded, 651 people (310 males, 341 females); 1 person removed due to --mind; 19853 variants removed due to --geno; --hwe: 399 variants removed; 24435 variants removed due to maf; 810934 variants and 650 people pass filters and QC.
    #OAG2 707097 variants loaded, 645 people (332 males, 313 females), 3  people removed due --mind; 2933 variants removed due to --geno; --hwe: 328 variants removed; 13161 variants removed due to --maf; 690675 variants and 642 people pass filters and QC.
    #610k_cases 548469 variants loaded, 925 people (436 males, 489 females) loaded, 0 people removed due to --mind, 650 variants removed due to --geno, --hwe: 196 variants removed; 0 variants removed due to minor allele threshold 547623 variants and 925 people pass filters and QC.
    #Omni_cases: 802474 variants loaded 694 people (346 males, 348 females) loaded; 0 people removed due to --mind; 2426 variants removed due to --geno; --hwe: 40 variants removed; 0 variants removed due to --maf, 800008 variants and 694 people pass filters and QC
    #SDH: 938370 variants loaded; 570 people (449 males, 121 females) loaded; 0 people removed due to --mind; 182 variants removed due to --geno; --hwe: 84 variants removed; 103454 variants removed due to  --maf; 834650 variants and 570 people pass filters and QC.
    #WAMHS 903034 variants loaded 1273 people (740 males, 533 females) loaded; 0 people removed due to --mind; 2490 variants removed due to --geno; -hwe: 873 variants removed; 216876 variants removed due to --maf; 682795 variants and 1273 people pass filters and QC.
    #EPIGENE 393059 variants loaded 784 people (520 males, 264 females) loaded; 0 people removed due to missing genotype data (--mind) 97 variants removed due to missing genotype data (--geno). --hwe: 0 . 97028 variants removed due to minor allele threshold(s) 295934 variants and 784 people pass filters and QC.

  #temp_2016_cleaning_QMEGA_610k_controls_aligned_1KG needs to be cleaned seperately; 0.03 missingness dropping 2160 people, I bet due to differences between 610k and 670k arrays. Likewise need to handle QTWIN diff - pulled from a larger set of genotypes so a number of 0 freq SNPs
  for k in temp_2016_cleaning_QTWIN_controls_aligned_1KG temp_2016_cleaning_EPIGENE_cases_aligned_1KG  #temp_2016_cleaning_QMEGA_610k_controls_aligned_1KG
  do
    plink_1.90 --threads 1 --bfile ${dirI}/${k} --geno 0.3 --make-bed --out ${dirI}/${k}_2
    #QMEGA_610k 533087 variants loaded, 3954 people (836 males, 3118 females) loaded. 29637 variants removed due to --geno. 503450 variants and 3954 people pass filters and QC.
    #QTWIN 442820 variants loaded from .bim file. Total genotyping rate is 0.934486. 31052 variants removed due to -geno. 411768 variants and 981 people pass filters and QC.
    echo ""
    plink_1.90 --threads 1 --bfile ${dirI}/${k}_2 --geno 0.03 --mind 0.03 --maf 0.001 --hwe 5e-10 midp include-nonctrl --make-bed --out ${dirI}/${k}_min_clean
    #QMEGA_610k 503450 variants loaded; 3954 people (836 males, 3118 females) loaded; 1 person removed due to --mind; 115 variants removed due to --geno; --hwe: 0 variants removed 0 variants removed due to maf; 503335 variants and 3953 people pass filters and QC.
    #QTWIN 411768 variants loaded from .bim file 0 people removed--mind. Total genotyping rate is 0.999682 157 variants removed due to  (--geno) --hwe: 4 variants removed due to Hardy-Weinberg exact test. 98106 variants removed due to minor allele threshold(s) 313501 variants and 981 people pass filters and QC.
    rm -f ${dirI}/${k}_2
  done #handling 610/QTWIN sep

done #loop on/off

##################################################################################
#5.0 Recombine the control and case sets
##################################################################################

for x in #1
do
  echo -e "\nCombining AMFS cases and controls\n"
  plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_AMFS_cases_aligned_1KG_min_clean --bmerge ${dirI}/temp_2016_cleaning_AMFS_controls_aligned_1KG_min_clean --make-bed --out ${dirI}/temp_2016_cleaning_AMFS_aligned_1KG 
    #Total genotyping rate is 0.996949. 801414 variants and 978 people pass filters and QC
  echo -e "\nCombining QMEGA_610k cases and controls\n"
  plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_QMEGA_610k_cases_aligned_1KG_min_clean --bmerge ${dirI}/temp_2016_cleaning_QMEGA_610k_controls_aligned_1KG_min_clean --make-bed --out ${dirI}/temp_2016_cleaning_QMEGA_610k_aligned_1KG
    #Total genotyping rate is 0.92944. 550293 variants and 4878 people pass filters and QC
  echo -e "\nCombining QMEGA_omni cases and SDH controls\n"
  #A couple of run throughs identified various SNPs on the arrays, but not in 1KG ref files, with different positions between SDH and QMEGA_omni; a number had been zeroed out. Where they could be matched by position/freq to real SNPs in dbSNP added and initial step at the start of the script to remap them otherwise exclude them. Some requiring flipping still (as not in 1kG ref files I have can't align them to a set strand.
  plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_QMEGA_omni_cases_aligned_1KG_min_clean --bmerge ${dirI}/temp_2016_cleaning_SDH_2_aligned_1KG_min_clean  --make-bed --out ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG
    #963 with 3+ alleles
  echo ""
  #check for ambigious (only would happen if something borked above - should have been removed if I couldn't align to 1kG )
  awk 'NR==FNR {a[$1];next} $2 in a' ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG-merge.missnp ${dirI}/temp_2016_cleaning_SDH_2_aligned_1KG_min_clean.bim | awk '($5=="A" && $6=="T") || ($5=="T" && $6=="A") || ($5=="C" && $6=="G") || ($5=="G" && $6=="C")' | wc -l
   #0
  echo ""
  plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_SDH_2_aligned_1KG_min_clean --flip ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG-merge.missnp --make-bed --out ${dirI}/temp_2016_cleaning_SDH_2_aligned_1KG_min_clean_flip
    #Flip 963 SNPs flipped
  echo ""
  plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_QMEGA_omni_cases_aligned_1KG_min_clean --bmerge ${dirI}/temp_2016_cleaning_SDH_2_aligned_1KG_min_clean_flip  --make-bed --out ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG
   #359 same position warning - agin will be SNPs where not in 1KG, or the 1KG is not and rsID
     #Total genotyping rate is 0.960236. 848806 variants and 1264 people pass filters and QC.
  echo ""
  grep Warning ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG.log |  grep "have the same position" | cut -d" " -f3,5 | sed "s/'//g" > ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG_same_pos; wc -l ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG_same_pos
    #359
  echo -e "\nFollowing is head of the files of SNPs that have overlapping positions but different names"
  awk 'NR==FNR {a[$2]=$0;next} $2 in a {print a[$2],$0}' <(awk 'NR==FNR {a[$2]=$0;next} $1 in a {print $0,a[$1]}' <(sed 's/\t/ /g' ${dirI}/temp_2016_cleaning_SDH_2_aligned_1KG_min_clean_flip.bim) ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG_same_pos) <( sed 's/\t/ /g' ${dirI}/temp_2016_cleaning_QMEGA_omni_cases_aligned_1KG_min_clean.bim) > ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG_same_pos_1; head ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG_same_pos_1
    #before 13/07/2016, one of these was ambigious, and looked to be on the wrong strand for one of them. But they should have only been retained if they were alignable to 1kg. : kgp1064009 SNP1-65860417 1 kgp1064009 0 66087829 T A 1 SNP1-65860417 0 66087829 A T
  echo -e "\nNone should be ambigious; checking..."
  awk 'NR==FNR {a[$2]=$0;next} $2 in a {print a[$2],$0}' <(awk 'NR==FNR {a[$2]=$0;next} $1 in a {print $0,a[$1]}' <(sed 's/\t/ /g' ${dirI}/temp_2016_cleaning_SDH_2_aligned_1KG_min_clean_flip.bim) ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG_same_pos) <( sed 's/\t/ /g' ${dirI}/temp_2016_cleaning_QMEGA_omni_cases_aligned_1KG_min_clean.bim) | awk '($7=="A" && $8=="T") || ($7=="T" && $8=="A") || ($7=="G" && $8=="C") || ($7=="C" && $8=="G")' | wc -l
    #0 now I have fixed the alignment process
  echo ""
  #above I found that you can wrongly combine SNPs by assuming the alleles make sense - so should really check here with freq
  plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_QMEGA_omni_cases_aligned_1KG_min_clean --freq --out ${dirI}/temp_2016_cleaning_QMEGA_omni_cases_aligned_1KG_min_clean --silent
  plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_SDH_2_aligned_1KG_min_clean --freq --out ${dirI}/temp_2016_cleaning_SDH_2_aligned_1KG_min_clean --silent

  awk 'NR==FNR {a[$2]=$0;next} $1 in a {print $0,a[$1]}' <(sed 's/ \+/ /g' ${dirI}/temp_2016_cleaning_SDH_2_aligned_1KG_min_clean.frq | sed 's/^ //g') ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG_same_pos_1 > ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG_same_pos_2; wc -l ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG_same_pos_2
    #359
  echo ""
  awk 'NR==FNR {a[$2]=$0;next} $2 in a {print $0,a[$2]}' <(sed 's/ \+/ /g' ${dirI}/temp_2016_cleaning_QMEGA_omni_cases_aligned_1KG_min_clean.frq | sed 's/^ //g') ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG_same_pos_2 > ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG_same_pos_3; wc -l ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG_same_pos_3
    #359 
  echo -e "\nDisplaying SNPs with a more than 5% diff between cases and controls but the same positions"
  awk '($19-$25) < (-0.05) || ($19-$25) > 0.05' ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG_same_pos_3
    #kgp3126359 SNP4-9846076 4 kgp3126359 0 10236978 C A 4 SNP4-9846076 0 10236978 G T 4 kgp3126359 C A 0.3849 1138 4 SNP4-9846076 G T 0.4394 1386
    #kgp9471170 SNP4-9860645 4 kgp9471170 0 10251547 A C 4 SNP4-9860645 0 10251547 T G 4 kgp9471170 A C 0.414 1140 4 SNP4-9860645 T G 0.469 1388
    #kgp12176916 SNP4-10021586 4 kgp12176916 0 10412488 G A 4 SNP4-10021586 0 10412488 C T 4 kgp12176916 G A 0.4009 1140 4 SNP4-10021586 C T 0.4517 1388
    #kgp10197490 SNP20-32543053 20 kgp10197490 0 33079392 C A 20 SNP20-32543053 0 33079392 G T 20 kgp10197490 C A 0.1868 1140 20 SNP20-32543053 G T 0.2478 1380
    #kgp19180985 SNP20-33199732 20 kgp19180985 0 33736071 G A 20 SNP20-33199732 0 33736071 C T 20 kgp19180985 G A 0.2456 1140 20 SNP20-33199732 C T 0.2968 1388

    #are any in 1kG? E.g. is it just that 1kg has an equaly ambigious name? 
    ##zcat ALL_1000G_phase1integrated_v3_ALL_chr_legends_with_chromosome.txt.gz | awk '($1==4 && $3==10236978) || ($1==4 && $3==10251547) || ($1==4 && $3==10412488) || ($1==20 && $3==33079392) || ($1==20 && $3==33736071)'
      #Nope
  echo ""
  #as no ambigious can now use this to work out which ones to flip - do a quick test of $7!=$13; as long as not very close to 0.5 MAF should be testing the same allele. Manually checked a bunch and that looks to be the case. Nope there is one - SNP6-135367351/kgp16958124 with MAF _ 0.48 so have to check both
  awk 'NR==FNR {a[$2]=$0;next} $2 in a {print a[$2],$0}' <(awk 'NR==FNR {a[$2]=$0;next} $1 in a {print $0,a[$1]}' <(sed 's/\t/ /g' ${dirI}/temp_2016_cleaning_SDH_2_aligned_1KG_min_clean_flip.bim) ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG_same_pos) <( sed 's/\t/ /g' ${dirI}/temp_2016_cleaning_QMEGA_omni_cases_aligned_1KG_min_clean.bim) | awk '!(($7==$13 && $8==$14) || ($7==$14 && $8==$13)) {print $1}' >> ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG-merge.missnp; wc -l ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG-merge.missnp
    #1135
  echo -e "\nNow flipping those SNPs with overlapping positions, are non-ambigious to the same strand\n"
  #redo the flip
  plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_SDH_2_aligned_1KG_min_clean --flip ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG-merge.missnp --make-bed --out  ${dirI}/temp_2016_cleaning_SDH_2_aligned_1KG_min_clean_flip
    #--1139 flipped
  #now update the names (may as well match the melanoma set)
  echo ""
  #Before updating names make sure alleles are the same, or alleles swapped A1 and A2 but MAF > 0.5 
  awk '!(($7==$13 && $8==$14) || ($7==$14 && $8==$13 && $19>=0.45) || ($7=="G" && $8=="T" && $13=="C" && $14=="A") || ($7=="A" && $8=="C" && $13=="T" && $14=="G") || ($7=="A" && $8=="G" && $13=="T" && $14=="C") || ($7=="T" && $8=="C" && $13=="A" && $14=="G") || ($7=="G" && $8=="A" && $13=="C" && $14=="T") || ($7=="C" && $8=="A" && $13=="G" && $14=="T") || ($7=="C" && $8=="T" && $13=="G" && $14=="A"))' ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG_same_pos_3
    #kgp1718593 SNP6-135365744 6 kgp1718593 0 135324051 G A 6 SNP6-135365744 0 135324051 T C 6 kgp1718593 G A 0.4842 1136 6 SNP6-135365744 T C 0.4942 1386 << strand flip and MAF close to 0.5 so swapped A1 and A2
    #kgp2808278 SNP20-32070151 20 kgp2808278 0 32606490 A C 20 SNP20-32070151 0 32606490 G T 20 kgp2808278 A C 0.4702 1140 20 SNP20-32070151 G T 0.4971 1386 << strand flip and MAF  close to 0.5 so swapped A1 and A2   
    #okay, so happy to rename
  awk '{print $1,$2}' ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG_same_pos_3 >  ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG_same_pos_3_update
  echo ""
  plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_SDH_2_aligned_1KG_min_clean_flip --update-name ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG_same_pos_3_update --make-bed --out ${dirI}/temp_2016_cleaning_SDH_2_aligned_1KG_min_clean_flip_2
     #--update-name: 359 values updated.
  echo "" 
  plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_QMEGA_omni_cases_aligned_1KG_min_clean --bmerge ${dirI}/temp_2016_cleaning_SDH_2_aligned_1KG_min_clean_flip_2 --make-bed --out ${dirI}/temp_2016_cleaning_QMEGA_omni_SDH_aligned_1KG
    #Total genotyping rate is 0.960647.  848436 variants and 1264 people pass filters and QC.
 
  #interimn clean
  rm -f ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG_same_pos_3 ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG_same_pos_2 ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG_same_pos_1 ${dirI}/temp_2016_cleaning_QMEGA_omni_cases_aligned_1KG_min_clean.frq ${dirI}/temp_2016_cleaning_SDH_2_aligned_1KG_min_clean.frq ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG_same_pos_3_update ${dirI}/temp_2016_cleaning_SDH_2_aligned_1KG_min_clean_flip*.* ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG-merge.missnp 

  #combining the HEI sets
  echo -e "\nCombining the HEI sets\n"

  plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_HEIDELBERG_cases_aligned_1KG_min_clean --bmerge ${dirI}/temp_2016_cleaning_HEIDELBERG_controls_aligned_1KG_min_clean --make-bed --out ${dirI}/temp_2016_cleaning_HEIDELBERG_aligned_1KG
    #no flip warnings, not miss positions warning. Suspiciously easy... Total genotyping rate is 0.989428. 217287 variants and 2437 people pass filters and QC.

  echo -e "\nCombining the control sets for WAMHS\n"

   #originally getting 12 flips but they were just errors in SNP remapping, fixed now
  plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_OAGphase1_controls_aligned_1KG_min_clean --bmerge ${dirI}/temp_2016_cleaning_OAGphase2_controls_aligned_1KG_min_clean --make-bed --out ${dirI}/temp_2016_cleaning_OAG_controls_aligned_1KG
    #Total genotyping rate is 0.858391. 873861 variants and 1292 people pass filters and QC.
  echo ""
    #combine in the IBD set now. Was getting ~76 mutliple position warnings; matching with dbSNP indicates IBD positions correct, intergrated into initial remap step. 4 had different chromosomes and diff bp, checking dbSNP confirms in hapmap were on diff chromosomes, now on the correct chromosome in IBD; freq suggests they are correct.
      #8 rs10106770 76.01 65568328 G A 0.3006 1284 2 rs10106770 0 235832763 C T 0.3069 1854 <<no freq in dbSNP to confirm, but match each other
      #3 rs12496398 216.149 194616354 G A 0.125 1256 4 rs12496398 0 115156497 C T 0.1305 1854 <<no freq in dbSNP to confirm, but match each other
      #5 rs2569201 152.875 149260811 G A 0.06975 1276 8 rs2569201 0 109657233 C T 0.08091 1854 << match each other and dbSNP freq
      #1 rs12043679 145.981 152231425 C A 0.1677 1282 13 rs12043679 0 66066597 C A 0.1753 1854 <<quad allelic in 1kg
  plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_OAG_controls_aligned_1KG --bmerge ${dirI}/temp_2016_cleaning_IBD_controls_aligned_1KG_min_clean --make-bed --out ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG
    # 725 variants with 3+ alleles present
  echo ""
    #IBD has been correct in terms of position so flip OAG
  plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_OAG_controls_aligned_1KG --flip ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG-merge.missnp --make-bed --out ${dirI}/temp_2016_cleaning_OAG_controls_aligned_1KG_flip
      # 725 variants flipped
    #Was getting more same position warnings - added SNPs that are impossible to handle (tri or quad allelic SNPs where two exm arrays assay for each combintation) to fix_exclude
  echo ""
  plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_OAG_controls_aligned_1KG_flip --bmerge ${dirI}/temp_2016_cleaning_IBD_controls_aligned_1KG_min_clean --make-bed --out ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG
    #447 same position warning
  #for some reason I don't understand IBD is the second field in the same position warning, except for one pairing - exm1574697 VG21S23270 - so reverse that if present
  echo -e "\nFor an unknown reason PLINK is reversing the log reporting (file 2 is 2nd) in the overlapping SNP pair exm1574697 VG21S23270 so fixing manually\n"
  grep Warning ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG.log |  grep "have the same position" | cut -d" " -f3,5 | sed "s/'//g" | sed 's/exm1574697 VG21S23270/VG21S23270 exm1574697/g' > ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG_same_pos; wc -l ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG_same_pos
    #447
  echo -e "\nFollowing is head of the files of SNPs that have overlapping positions but different names:"
  awk 'NR==FNR {a[$2]=$0;next} $2 in a {print a[$2],$0}' <(awk 'NR==FNR {a[$2]=$0;next} $1 in a {print $0,a[$1]}' <(sed 's/\t/ /g' ${dirI}/temp_2016_cleaning_OAG_controls_aligned_1KG_flip.bim) ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG_same_pos) <( sed 's/\t/ /g' ${dirI}/temp_2016_cleaning_IBD_controls_aligned_1KG_min_clean.bim) > ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG_same_pos_1; head -4 ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG_same_pos_1
   #rs1052571 chr1:15850613 1 rs1052571 0 15850613 G A 1 chr1:15850613 0 15850613 G A
   #rs17887074 chr1:22964176 1 rs17887074 0 22964176 A G 1 chr1:22964176 0 22964176 A G
   #rs41306580 chr1:23847534 1 rs41306580 0 23847534 A C 1 chr1:23847534 0 23847534 A C
   #rs41263977 chr1:32193185 1 rs41263977 0 32193185 A G 1 chr1:32193185 0 32193185 A G

  echo -e "\nNone should be ambigious; checking..."
   #was getting three more ambigious SNPs with name mismatch, but three could be fixed by flipping the exm SNP in the IBD set, which meant it could be paired with its omni express non exm counterpart (e.g. exm-rs9261434 and rs9261434) and thus removed
  awk 'NR==FNR {a[$2]=$0;next} $2 in a {print a[$2],$0}' <(awk 'NR==FNR {a[$2]=$0;next} $1 in a {print $0,a[$1]}' <(sed 's/\t/ /g' ${dirI}/temp_2016_cleaning_OAG_controls_aligned_1KG_flip.bim) ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG_same_pos) <( sed 's/\t/ /g' ${dirI}/temp_2016_cleaning_IBD_controls_aligned_1KG_min_clean.bim) | awk '($7=="A" && $8=="T") || ($7=="T" && $8=="A") || ($7=="G" && $8=="C") || ($7=="C" && $8=="G")' | wc -l
    #4..... I bet the name in 1000 genomes is chr......
    #rs41273533 chr1:158655129 1 rs41273533 0 158655129 G C 1 chr1:158655129 0 158655129 G C
    #rs35400920 chr4:36163089 4 rs35400920 0 36163089 G C 4 chr4:36163089 0 36163089 G C
    #rs8187788 chr7:87092143 7 rs8187788 0 87092143 C G 7 chr7:87092143 0 87092143 C G
    #rs41279214 chr15:34032188 15 rs41279214 0 34032188 G C 15 chr15:34032188 0 34032188 G C

  #so the four that didn't have their name fixed were aligned to 1KG in terms of allele, just not name fixeed as 1KG name not an rsID.
  ##zcat ALL_1000G_phase1integrated_v3_ALL_chr_legends_with_chromosome.txt.gz | awk '($1==1 && $3==158655129) || ($1==4 && $3==36163089) || ($1==7 && $3==87092143) || ($1==15 && $3==34032188)'
    #1 chr1:158655129 158655129 C G 0.0000 0.0028 0.0000 0.0013 0.0000 0.0028 0.0000 0.0013 <<so neither had their name updated as no better/diff than 1KG
    #4 chr4:36163089 36163089 C G 0.0000 0.0000 0.0000 0.0026 0.0000 0.0000 0.0000 0.0026
    #7 chr7:87092143 87092143 G C 0.0000 0.0000 0.0000 0.0013 0.0000 0.0000 0.0000 0.0013
    #15 chr15:34032188 34032188 C G 0.0000 0.0000 0.0000 0.0013 0.0000 0.0000 0.0000 0.0013
 
  #before go further check the freq but I wonder if my selecting name function is not specific enough but I wonder if my selecting name function is not specific enough..
  echo ""
  #above I found that you can wrongly combine SNPs by assuming the alleles make sense - so should really check here with freq
  plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_OAG_controls_aligned_1KG_flip --freq --out ${dirI}/temp_2016_cleaning_OAG_controls_aligned_1KG_flip --silent
  plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_IBD_controls_aligned_1KG_min_clean --freq --out ${dirI}/temp_2016_cleaning_IBD_controls_aligned_1KG_min_clean --silent

  awk 'NR==FNR {a[$2]=$0;next} $1 in a {print $0,a[$1]}' <(sed 's/ \+/ /g' ${dirI}/temp_2016_cleaning_OAG_controls_aligned_1KG_flip.frq | sed 's/^ //g') ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG_same_pos_1 > ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG_same_pos_2; wc -l ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG_same_pos_2
   #447; was n-1 until I manually reversed the pairing of exm1574697 VG21S23270   
  echo ""
  awk 'NR==FNR {a[$2]=$0;next} $2 in a {print $0,a[$2]}' <(sed 's/ \+/ /g' ${dirI}/temp_2016_cleaning_IBD_controls_aligned_1KG_min_clean.frq | sed 's/^ //g') ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG_same_pos_2 > ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG_same_pos_3; wc -l ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG_same_pos_3
    #295
  echo -e "\nDisplaying SNPs with a more than 5% diff between cases and controls but the same positions"
  #need to check the alleles make sense too - not just for flips but for literally different SNP
  awk '($19-$25) < (-0.05) || ($19-$25) > 0.05' ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG_same_pos_3
    #rs1403543 exm-rs1403543 23 rs1403543 0 115302192 G A 23 exm-rs1403543 0 115302192 A G 23 rs1403543 G A 0.4263 990 23 exm-rs1403543 A G 0.4993 1434
      #that is a reasonably big freq diff - still the exm assays are for the same SNP (checked dbSNP) so update. Note that allele order order has swapped - so it is actually about a 9% diff, not 7%
  echo ""
  awk '$7!=$13' ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG_same_pos_3 | wc -l 
    #35
  echo -e "\nThe following SNP is not a simple flip; check it is only dislaying rs1403543, which has had the minor/major swap"
    #look for ones that don't make sense before renaming
  awk '!(($7==$13 && $8==$14) || ($7==$14 && $8==$13 && $19>=0.45) || ($7=="G" && $8=="T" && $13=="C" && $14=="A") || ($7=="A" && $8=="C" && $13=="T" && $14=="G") || ($7=="A" && $8=="G" && $13=="T" && $14=="C") || ($7=="T" && $8=="C" && $13=="A" && $14=="G") || ($7=="G" && $8=="A" && $13=="C" && $14=="T") || ($7=="C" && $8=="A" && $13=="G" && $14=="T") || ($7=="C" && $8=="T" && $13=="G" && $14=="A"))' ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG_same_pos_3
   #rs1403543 exm-rs1403543 23 rs1403543 0 115302192 G A 23 exm-rs1403543 0 115302192 A G 23 rs1403543 G A 0.4263 990 23 exm-rs1403543 A G 0.4993 1434 <<same as above, allele order swapped as minor allele swapped
   #1434 
  echo ""
   #looks like the OAG content is preferable in terms of names 
  awk '{print $2,$1}' ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG_same_pos_3 > ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG_same_pos_3_update
  cut -d" " -f1 ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG_same_pos_3_update | sort | uniq -d
    #0
  cut -d" " -f2 ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG_same_pos_3_update | sort | uniq -d
    #0
  echo ""
  plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_IBD_controls_aligned_1KG_min_clean --update-name ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG_same_pos_3_update --make-bed --out ${dirI}/temp_2016_cleaning_IBD_controls_aligned_1KG_min_clean_2
   #--update-name: 447 values updated.
  echo ""
  #look for any new flips
  plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_OAG_controls_aligned_1KG --bmerge ${dirI}/temp_2016_cleaning_IBD_controls_aligned_1KG_min_clean_2 --make-bed --out ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG
   #766 - so + ~24 after renaming
  echo ""
  plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_OAG_controls_aligned_1KG --flip ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG-merge.missnp --make-bed --out ${dirI}/temp_2016_cleaning_OAG_controls_aligned_1KG_flip
      # 766  variants flipped
  echo ""
  plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_OAG_controls_aligned_1KG_flip --bmerge ${dirI}/temp_2016_cleaning_IBD_controls_aligned_1KG_min_clean_2 --make-bed --out ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG
    #Total genotyping rate is 0.788449. 923249 variants and 2219 people pass filters and QC.
  #clean up
  rm -f ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG_same_pos ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG-merge.missnp 

  echo -e "\nNow merging WAMHS data into combined control data"
  #now merge in WAMHS. Was getting 4 multiple positin and 10 multiple chromsome warnings. Fixed in main script - 1 SNP was weird (no dbSNP info) so dropped, 3 had wrong OAG positions by 1bp, the 10 chr positions were SNPs assigned to X_PAR or X_nonPAR wrongly.
  plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_WAMHS_cases_aligned_1KG_min_clean --bmerge ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG --make-bed --out ${dirI}/temp_2016_cleaning_WAMHS_OAG_IBD_aligned_1KG
    #Error: 2840 variants with 3+ alleles present.
  echo ""
  #try flipping
  plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_WAMHS_cases_aligned_1KG_min_clean --flip ${dirI}/temp_2016_cleaning_WAMHS_OAG_IBD_aligned_1KG-merge.missnp --make-bed --out ${dirI}/temp_2016_cleaning_WAMHS_cases_aligned_1KG_min_clean_flip
    #2840 flipped
  echo ""
  plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_WAMHS_cases_aligned_1KG_min_clean_flip --bmerge ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG --make-bed --out ${dirI}/temp_2016_cleaning_WAMHS_OAG_IBD_aligned_1KG
   #no more flip errors anymore but 432 (down from ~7400 before I updated handling overlapping SNPS) same position warning. Many look to be exm overlaps within WAMHS data not picked up earlier due to the exm data being flipped
   #Warning: Variants 'rs1052571' and 'chr1:15850613' have the same position.
   #Warning: Variants 'rs41306580' and 'chr1:23847534' have the same position.
   #Warning: Variants 'rs41263977' and 'chr1:32193185' have the same position.
   #432 more same-position warnings: see log file.
   #most of these appear to be duplicated in WAMHs, but not fixed as flipped in WAMHS.....
  echo ""
   grep "Warning" ${dirI}/temp_2016_cleaning_WAMHS_OAG_IBD_aligned_1KG.log | grep "have the same position" | wc -l
  echo ""
     #435
   grep "Warning" ${dirI}/temp_2016_cleaning_WAMHS_OAG_IBD_aligned_1KG.log | grep "have the same position" | cut -d" " -f3,5 | sed "s/'//g" | sed 's/\.//g' | sed 's/exm686274 exm2258773/exm2258773 exm686274/g' | sed 's/exm2251671 exm1134438/exm1134438 exm2251671/g' | sed 's/exm1574697 VG21S23270/VG21S23270 exm1574697/g' | sed 's/rs2032597 exm-rs2032597/exm-rs2032597 rs2032597/g' | sed 's/rs34442126 exm-rs34442126/exm-rs34442126 rs34442126/g' |  sed 's/rs35474563 exm-rs35474563/exm-rs35474563 rs35474563/g' | sed 's/rs700449 exm2273221/exm2273221 rs700449/g'  > ${dirI}/temp_2016_cleaning_WAMHS_OAG_IBD_aligned_1KG_same_pos_1; wc -l ${dirI}/temp_2016_cleaning_WAMHS_OAG_IBD_aligned_1KG_same_pos_1
     #435. How many both from WAMHS?
  echo -e "\nChecking for overlapping SNPs that are both in WAMHS; should be 0:"
   awk 'NR==FNR {a[$2];next} $1 in a && $2 in a' ${dirI}/temp_2016_cleaning_WAMHS_cases_aligned_1KG_min_clean.bim ${dirI}/temp_2016_cleaning_WAMHS_OAG_IBD_aligned_1KG_same_pos_1 | wc -l
  echo -e "\nChecking for overlapping SNPs that are both in OAG + IBD; should be 0:"
    #0 (down from 7298 before I fixed handling overlapping SNPs)
   awk 'NR==FNR {a[$2];next} $1 in a && $2 in a' ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG.bim ${dirI}/temp_2016_cleaning_WAMHS_OAG_IBD_aligned_1KG_same_pos_1 | wc -l
  echo ""
   #0 good; 
   #so lets combine them and make sure they are the right SNPs. We again have 7 SNPs were for some reason OAG+IBD is reported second, and WAMHS first, so manually swap
   awk 'NR==FNR {a[$2]=$0;next} $2 in a {print a[$2],$0}' <(awk 'NR==FNR {a[$2]=$0;next} $1 in a {print $0,a[$1]}' <(sed 's/\t/ /g' ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG.bim) ${dirI}/temp_2016_cleaning_WAMHS_OAG_IBD_aligned_1KG_same_pos_1) <( sed 's/\t/ /g' ${dirI}/temp_2016_cleaning_WAMHS_cases_aligned_1KG_min_clean_flip.bim) > ${dirI}/temp_2016_cleaning_WAMHS_OAG_IBD_aligned_1KG_same_pos_2; wc -l ${dirI}/temp_2016_cleaning_WAMHS_OAG_IBD_aligned_1KG_same_pos_2
    #435 (was 428 before fixing the order of 7 of them)
  echo ""
  #doesn't really matter whose names get updated - which one has more rs rather than exm-rs? Seems neater
  grep -c ^rs ${dirI}/temp_2016_cleaning_WAMHS_OAG_IBD_aligned_1KG_same_pos_2
     #428
  echo ""
  grep -c -w " rs" ${dirI}/temp_2016_cleaning_WAMHS_OAG_IBD_aligned_1KG_same_pos_2
     #0
  echo ""
  grep -v ^rs ${dirI}/temp_2016_cleaning_WAMHS_OAG_IBD_aligned_1KG_same_pos_2
  echo ""
    #exm2258773 exm686274 8 exm2258773 0 19809375 C A 8 exm686274 0 19809375 C A
    #exm1134438 exm2251671 14 exm1134438 0 105622160 A G 14 exm2251671 0 105622160 A G
    #VG21S23270 exm1574697 21 VG21S23270 63.0669 44483184 G A 21 exm1574697 0 44483184 G A
    #exm-rs2032597 rs2032597 24 exm-rs2032597 0 14847792 C A 24 rs2032597 0 14847792 C A
    #exm-rs34442126 rs34442126 24 exm-rs34442126 0 14922583 C T 24 rs34442126 0 14922583 G A
    #exm-rs35474563 rs35474563 24 exm-rs35474563 0 22917995 T C 24 rs35474563 0 22917995 A G
    #exm2273221 rs700449 25 exm2273221 0 155009856 T C 25 rs700449 0 155009856 A G
    #those are mostly the ones I had to manually fix - does it print rs preferentially in the first col?
  awk '{print $2,$1}' ${dirI}/temp_2016_cleaning_WAMHS_OAG_IBD_aligned_1KG_same_pos_2 > ${dirI}/temp_2016_cleaning_WAMHS_OAG_IBD_aligned_1KG_same_pos_2_update
  #update WAMHS
  plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_WAMHS_cases_aligned_1KG_min_clean_flip --update-name ${dirI}/temp_2016_cleaning_WAMHS_OAG_IBD_aligned_1KG_same_pos_2_update --make-bed --out ${dirI}/temp_2016_cleaning_WAMHS_cases_aligned_1KG_min_clean_flip_2
     #435 names updated
  echo ""
  plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_WAMHS_cases_aligned_1KG_min_clean_flip_2 --bmerge ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG --make-bed --out ${dirI}/temp_2016_cleaning_WAMHS_OAG_IBD_aligned_1KG
     #Error: 193 variants with 3+ alleles present.
  echo ""
  plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_WAMHS_cases_aligned_1KG_min_clean_flip_2 --flip ${dirI}/temp_2016_cleaning_WAMHS_OAG_IBD_aligned_1KG-merge.missnp --make-bed --out ${dirI}/temp_2016_cleaning_WAMHS_cases_aligned_1KG_min_clean_flip_3
    #--193 flipped
  echo ""
  plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_WAMHS_cases_aligned_1KG_min_clean_flip_3 --bmerge ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG --make-bed --out ${dirI}/temp_2016_cleaning_WAMHS_OAG_IBD_aligned_1KG
    #Total genotyping rate is 0.765493. 929113 variants and 3492 people pass filters and QC.

  #clean up
  rm -f ${dirI}/temp_2016_cleaning_WAMHS_cases_aligned_1KG_min_clean_flip*.* ${dirI}/temp_2016_cleaning_WAMHS_OAG_IBD_aligned_1KG-merge.missnp ${dirI}/temp_2016_cleaning_WAMHS_OAG_IBD_aligned_1KG_same_pos_* ${dirI}/temp_2016_cleaning_OAG_controls_aligned_1KG_flip.* ${dirI}/temp_2016_cleaning_IBD_controls_aligned_1KG_min_clean_* ${dirI}/temp_2016_cleaning_OAG_IBD_controls_aligned_1KG_same_pos_* 

  echo -e "\nCombining EPIGENE and QTWIN\n"
  plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_EPIGENE_cases_aligned_1KG_min_clean --bmerge ${dirI}/temp_2016_cleaning_QTWIN_controls_aligned_1KG_min_clean --make-bed --out ${dirI}/temp_2016_cleaning_EPIGENE_QTWIN_aligned_1KG
    #No flips, though previous merged by SG. 30 same position warnings
  echo ""
    #7 presented in reverse order, correct  
  grep "Warning" ${dirI}/temp_2016_cleaning_EPIGENE_QTWIN_aligned_1KG.log | grep "have the same position" | cut -d" " -f3,5 | sed "s/'//g" | sed 's/\.//g' | sed 's/rs178715 exm2273099/exm2273099 rs178715/g' | sed 's/rs6525433 exm1644651/exm1644651 rs6525433/g' | sed 's/rs6618889 exm2264764/exm2264764 rs6618889/g' | sed 's/rs7890941 exm2273180/exm2273180 rs7890941/g' |  sed 's/rs2843594 exm2262922/exm2262922 rs2843594/g' | sed 's/rs7885911 exm2273292/exm2273292 rs7885911/g' | sed 's/rs5911653 exm2268622/exm2268622 rs5911653/g' | sed 's/rs4289953 exm1656466/exm1656466 rs4289953/g' | sed 's/rs2761608 exm2273213/exm2273213 rs2761608/g' > ${dirI}/temp_2016_cleaning_EPIGENE_QTWIN_aligned_1KG_same_pos; wc -l ${dirI}/temp_2016_cleaning_EPIGENE_QTWIN_aligned_1KG_same_pos
    #30
  echo -e "\nFollowing is head of the files of SNPs that have overlapping positions but different names:"
  awk 'NR==FNR {a[$2]=$0;next} $2 in a {print a[$2],$0}' <(awk 'NR==FNR {a[$2]=$0;next} $1 in a {print $0,a[$1]}' <(sed 's/\t/ /g' ${dirI}/temp_2016_cleaning_EPIGENE_cases_aligned_1KG_min_clean.bim) ${dirI}/temp_2016_cleaning_EPIGENE_QTWIN_aligned_1KG_same_pos) <( sed 's/\t/ /g' ${dirI}/temp_2016_cleaning_QTWIN_controls_aligned_1KG_min_clean.bim) > ${dirI}/temp_2016_cleaning_EPIGENE_QTWIN_aligned_1KG_same_pos_1; head -4 ${dirI}/temp_2016_cleaning_EPIGENE_QTWIN_aligned_1KG_same_pos_1
    #rs1974522 exm1626530 23 rs1974522 0 3238733 G A 23 exm1626530 0 3238733 G A
    #rs7885458 exm-rs7885458 23 rs7885458 0 6220474 C T 23 exm-rs7885458 0 6220474 C T
    #rs12689910 exm1626935 23 rs12689910 0 7023678 A G 23 exm1626935 0 7023678 A G
    #rs6639094 exm2273087 23 rs6639094 0 11889104 T C 23 exm2273087 0 11889104 T C

  echo -e "\nNone should be ambigious; checking..."
  awk 'NR==FNR {a[$2]=$0;next} $2 in a {print a[$2],$0}' <(awk 'NR==FNR {a[$2]=$0;next} $1 in a {print $0,a[$1]}' <(sed 's/\t/ /g' ${dirI}/temp_2016_cleaning_EPIGENE_cases_aligned_1KG_min_clean.bim) ${dirI}/temp_2016_cleaning_EPIGENE_QTWIN_aligned_1KG_same_pos_1) <( sed 's/\t/ /g' ${dirI}/temp_2016_cleaning_QTWIN_controls_aligned_1KG_min_clean.bim  ) | awk '($7=="A" && $8=="T") || ($7=="T" && $8=="A") || ($7=="G" && $8=="C") || ($7=="C" && $8=="G")' | wc -l
    #0
  #check freq to make sure are the same
  plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_EPIGENE_cases_aligned_1KG_min_clean --freq --out ${dirI}/temp_2016_cleaning_EPIGENE_cases_aligned_1KG_min_clean --silent
  plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_QTWIN_controls_aligned_1KG_min_clean --freq --out ${dirI}/temp_2016_cleaning_QTWIN_controls_aligned_1KG_min_clean --silent

  awk 'NR==FNR {a[$2]=$0;next} $1 in a {print $0,a[$1]}' <(sed 's/ \+/ /g' ${dirI}/temp_2016_cleaning_EPIGENE_cases_aligned_1KG_min_clean.frq | sed 's/^ //g') ${dirI}/temp_2016_cleaning_EPIGENE_QTWIN_aligned_1KG_same_pos_1 > ${dirI}/temp_2016_cleaning_EPIGENE_QTWIN_aligned_1KG_same_pos_2; wc -l ${dirI}/temp_2016_cleaning_EPIGENE_QTWIN_aligned_1KG_same_pos_2
  #30

   #any internal duplicates - shouldn't be any
  awk 'NR==FNR {a[$2]} $1 in a && $2 in a' ${dirI}/temp_2016_cleaning_EPIGENE_cases_aligned_1KG_min_clean.bim ${dirI}/temp_2016_cleaning_EPIGENE_QTWIN_aligned_1KG_same_pos_1
     #0
  awk 'NR==FNR {a[$2]} $1 in a && $2 in a' ${dirI}/temp_2016_cleaning_QTWIN_controls_aligned_1KG_min_clean.bim ${dirI}/temp_2016_cleaning_EPIGENE_QTWIN_aligned_1KG_same_pos_1
     #0
  echo ""
  awk 'NR==FNR {a[$2]=$0;next} $2 in a {print $0,a[$2]}' <(sed 's/ \+/ /g' ${dirI}/temp_2016_cleaning_QTWIN_controls_aligned_1KG_min_clean.frq | sed 's/^ //g') ${dirI}/temp_2016_cleaning_EPIGENE_QTWIN_aligned_1KG_same_pos_2 > ${dirI}/temp_2016_cleaning_EPIGENE_QTWIN_aligned_1KG_same_pos_3; wc -l ${dirI}/temp_2016_cleaning_EPIGENE_QTWIN_aligned_1KG_same_pos_3
    #30
  echo -e "\nDisplaying SNPs with a more than 5% diff between cases and controls but the same positions"
  #need to check the alleles make sense too - not just for flips but for literally different SNP
  awk '($19-$25) < (-0.05) || ($19-$25) > 0.05' ${dirI}/temp_2016_cleaning_EPIGENE_QTWIN_aligned_1KG_same_pos_3
    #None for QTWIN/EPIGENE
  echo ""
  awk '$7!=$13' ${dirI}/temp_2016_cleaning_EPIGENE_QTWIN_aligned_1KG_same_pos_3 | wc -l
    #0
    #look for ones that don't make sense before renaming
  awk '!(($7==$13 && $8==$14) || ($7==$14 && $8==$13 && $19>=0.45) || ($7=="G" && $8=="T" && $13=="C" && $14=="A") || ($7=="A" && $8=="C" && $13=="T" && $14=="G") || ($7=="A" && $8=="G" && $13=="T" && $14=="C") || ($7=="T" && $8=="C" && $13=="A" && $14=="G") || ($7=="G" && $8=="A" && $13=="C" && $14=="T") || ($7=="C" && $8=="A" && $13=="G" && $14=="T") || ($7=="C" && $8=="T" && $13=="G" && $14=="A"))' ${dirI}/temp_2016_cleaning_EPIGENE_QTWIN_aligned_1KG_same_pos_3
   #none

  #mix of rs and exm on one side or the other (due to which had the higher info at random I guess, so do a slightly more nuanced update and rename both
  awk '{
   if($2~/^rs/ && !($1~/^rs/))
    print $1,$2;
  if(!($2~/^rs/) && $1~/^rs/)
    print $2,$1;
  if($2~/^rs/ && $1~/^rs/)
    print $1,$2;
  if(!($2~/^rs/) && !($1~/^rs/))
    print $1,$2;
  }' ${dirI}/temp_2016_cleaning_EPIGENE_QTWIN_aligned_1KG_same_pos_3 > ${dirI}/temp_2016_cleaning_EPIGENE_QTWIN_aligned_1KG_update; wc -l ${dirI}/temp_2016_cleaning_EPIGENE_QTWIN_aligned_1KG_update

  for q in temp_2016_cleaning_EPIGENE_cases_aligned_1KG_min_clean temp_2016_cleaning_QTWIN_controls_aligned_1KG_min_clean
  do
    plink_1.90 --threads 1 --bfile ${dirI}/${q} --update-name ${dirI}/temp_2016_cleaning_EPIGENE_QTWIN_aligned_1KG_update --make-bed --out ${dirI}/${q}_2 --silent
      #EPIGENE --update-name: 9 values updated, 21 variant IDs not present. QTWIN --update-name: 21 values updated, 9 variant IDs not present.
  done
  echo ""
  plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_EPIGENE_cases_aligned_1KG_min_clean_2 --bmerge ${dirI}/temp_2016_cleaning_QTWIN_controls_aligned_1KG_min_clean_2 --make-bed --out ${dirI}/temp_2016_cleaning_EPIGENE_QTWIN_aligned_1KG
   #no more errors Total genotyping rate is 0.943689 323877 variants and 1765 people pass filters and QC.
 
  rm -f ${dirI}/temp_2016_cleaning_EPIGENE_cases_aligned_1KG_min_clean_2.* ${dirI}/temp_2016_cleaning_QTWIN_controls_aligned_1KG_min_clean_2.*

  #clean up
  #rm -f ${dirI}/temp_2016_cleaning_QMEGA_610k_cases_aligned_1KG.* ${dirI}/temp_2016_cleaning_SDH_2_aligned_1KG.* ${dirI}/temp_2016_cleaning_QMEGA_610k_controls_aligned_1KG.* ${dirI}/temp_2016_cleaning_AMFS_cases_aligned_1KG.* ${dirI}/temp_2016_cleaning_AMFS_controls_aligned_1KG.* ${dirI}/temp_2016_cleaning_QMEGA_omni_cases_aligned_1KG.* ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG-merge.missnp ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG_same_pos ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG_update_exclude ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG_update_bp ${dirI}/temp_2016_cleaning_SDH_2_aligned_1KG_flip.* ${dirI}/temp_2016_cleaning_SDH_2_aligned_1KG_2.* ${dirI}/temp_2016_cleaning_QMEGA_omni_cases_aligned_1KG_2.* ${dirI}/temp_2016_cleaning_QMEGA_omni_cases_aligned_1KG_3.* ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG_same_pos* ${dirI}/fix_remap ${dirI}/fix_exclude ${dirI}/fix_remap_chr ${dirI}/temp_2016_cleaning_OAGphase*_controls_aligned_1KG_min_clean.* ${dirI}/temp_2016_cleaning_WAMHS_cases_aligned_1KG_min_clean.* ${dirI}/temp_2016_cleaning_IBD_controls_aligned_1KG_min_clean.* ${dirI}/temp_2016_cleaning_SDH_2_aligned_1KG_min_clean.* ${dirI}/temp_2016_cleaning_QMEGA_omni_cases_aligned_1KG_min_clean.* ${dirI}/temp_2016_cleaning_QMEGA_610k_controls_aligned_1KG_min_clean.* ${dirI}/temp_2016_cleaning_QMEGA_610k_controls_aligned_1KG_2.* ${dirI}/temp_2016_cleaning_QMEGA_610k_cases_aligned_1KG_min_clean.* ${dirI}/temp_2016_cleaning_HEIDELBERG_controls_aligned_1KG_min_clean.* ${dirI}/temp_2016_cleaning_AMFS_controls_aligned_1KG_min_clean.* ${dirI}/temp_2016_cleaning_HEIDELBERG_cases_aligned_1KG_min_clean.* ${dirI}/temp_2016_cleaning_AMFS_cases_aligned_1KG_min_clean.*

done #merging sets

##################################################################################
#6.0 Basic tests of the re-merged data
##################################################################################
#flip scan; loop to turn on/off

for x in #1 
do
  for i in temp_2016_cleaning_AMFS_aligned_1KG temp_2016_cleaning_QMEGA_610k_aligned_1KG temp_2016_cleaning_QMEGA_omni_aligned_1KG temp_2016_cleaning_WAMHS_OAG_IBD_aligned_1KG temp_2016_cleaning_HEIDELBERG_aligned_1KG temp_2016_cleaning_EPIGENE_QTWIN_aligned_1KG 
  do
    plink_1.90 --threads 1 --bfile ${dirI}/${i} --flip-scan verbose --out ${dirI}/${i}_flipscan
      #AMFS --flip-scan: 7 variants with at least one negative LD match. Looking at the results for AMFS I am not that concertaned by this. The SNPs in question are rare so different LD patterns for a single person might be driving this (and might be real, or genotyping errors, no way to know). rs7250407/rs10405660 have no POS SNPs (no SNPs with positive LD). A sign of a flipped SNP would be a collection of SNPs each with POS LD (so say POS = 5) and a single SNP in NEG with all of them (see example on plink website). Still try flip scan verbose to get a better feel/
      #QMEGA_610k was 2 variants, now is 0 with my improved cleaning/aligning/merging
      #QMEGA_omni was 4, now 0 with my improved cleaning/aligning/merging
      #WAMHS/OAG/IBD --flip-scan verbose: 34 variants with at least one negative LD match
      #HEI --flip-scan verbose: 2 variants with at least one negative LD match. Was expecting more
      #EPIGENE/QTWIN  - 73 variants with at least one negaritive LD match

    echo ""
    awk '$9!=0' ${dirI}/${i}_flipscan.flipscan

    #still may as well play it safe and remove these SNPs. That said there are three omni SNPs that are okay, and only flagged due to one SNP - rs326150 - so want to keep rs326153\|rs6883729\|rs326162 (see below for working). WAMHS similar - 4 SNPs are flagged due to rs9818314 which is all over the place; keep them. Likewise rs17003186 impacting 4 others. In all other cases it is a pair of SNPs flagging each other, can't tell which is which so drop both
      #QTWIN/EPIGENE a quick eyeball; there are 73 SNp and at most 5 or 6 would be saved where there seems to be a single bad SNP. Weirdlt there is a combined set where 2 SNPs each disagree with two other SNPs - rs1969818|rs7919863|rs546640 - drop them all. Still in EPIGENE retained the few that looked to be caused only by a single bad SNP. rs57928695 and rs1026992
    awk 'NR>1 && $9!=0 {print $2}' ${dirI}/${i}_flipscan.flipscan | grep -v "rs326153\|rs6883729\|rs326162\|rs6763046\|rs4473559\|rs9831516\|rs9836305\|rs12006821\|rs4142629\|rs5915311\|rs6521881\|rs6757850\|rs6708255\|rs10193618\|rs8071181\|rs34146848\|rs9896466\|rs16944458" > ${dirI}/${i}_flipscan.flipscan_exclude; wc -l ${dirI}/${i}_flipscan.flipscan_exclude
     #WAMHS 26 AMFS 7 QMEGA 610k 0 QMEGA_omni 0 HEI 2. QTWIN/EPIGENE 66
    echo ""
    plink_1.90 --threads 1 --bfile ${dirI}/${i} --exclude ${dirI}/${i}_flipscan.flipscan_exclude --make-bed --out ${dirI}/temp_${i}
      #AMFS 801414 variants loaded from .bim file. --exclude: 801407 variants remaining. 802722 variants and 978 people pass filters and QC.
      #QMEGA_610k - none removed.
      #QMEGA_omni - none removed.
      #WAMHS 929113 variants loaded from .bim file. --exclude: 929087 variants remaining.
      #HEI 217287 variants loaded  --exclude: 217285 variants remaining.
      #QTWIN/EPIGENE 323877 variants loaded from .bim file. --exclude: 323811 variants remaining.
    echo ""
    plink_1.90 --threads 1 --bfile ${dirI}/temp_${i} --make-bed --out ${dirI}/${i}
    echo ""
    #double check
    plink_1.90 --threads 1  --bfile ${dirI}/${i} --flip-scan --out ${dirI}/${i}_flipscan
      #AMFS, QMEGA_610k QMEGA_omni WAMHS HEI --flip-scan: 0 variants. EPIGENE 0 now
    echo ""
    rm -f  ${dirI}/${i}_flipscan.* ${dirI}/temp_${i}.*
  done #flipscan removal
nk_1.90 --threads 1 --bfile ${dirI}/temp${i}_7 --extract ${dirI}/temp${i}_7_IBD.prune.in --genome --min 0.1 --out ${dirI}/temp${i}_7_IBD 
done
    ############################################################
    #AMFS
    ############################################################
      #CHR               SNP           BP   A1   A2        F    POS    R_POS    NEG    R_NEG NEGSNPS
      #19         rs7250407     53831824    C    T   0.0158      0       NA      1    0.279 rs10405660
      #19        rs10405660     53843999    G    A   0.0123      0       NA      1    0.279 rs7250407
      #23        rs12556807     91511898    G    T   0.0121      2    0.661      1    0.267 rs12014732
      #23        rs12014732     91520376    A    G  0.00762      3    0.748      1    0.267 rs12556807
      #23        rs12387794    130024275    C    T  0.00826      3    0.659      2    0.298 rs7063174|rs7057262
      #23         rs7063174    130168050    A    G   0.0108      4    0.817      1    0.298 rs12387794
      #23         rs7057262    130178487    C    T   0.0108      3    0.927      1    0.298 rs12387794
    ############################################################
    #610k_670k
    ############################################################
    #before I fixed the merge equal positon bug found a bunch of errors - now just these two - update - these no longer appear with the newly cleaned data 09/08 but useful to keep as it revealed an error in the merge equal position function in PLINK
      #   CHR         SNP           BP   A1   A2        F    POS    R_POS    NEG    R_NEG NEGSNPS
      #3   rs7647349    140848571    T    C   0.0811      0       NA      1    0.287 rs7614592
      #3   rs7614592    140863405    T    C    0.194      0       NA      1    0.287 rs7647349
      #3 rs7614592 140863405 0.195 1.000 1.000 2 0.668 0.875 0.664   <<while not terrible this is far from perfect; moderately well imputed (0.67) but only 0.87 concordant, so probably right to remove it
    ############################################################
    #QMEGA_omni
    ############################################################
    # <<<09/08 these no longer appear with the new cleaning, just keep a historical record of it
    #   CHR               SNP           BP   A1   A2        F    POS    R_POS    NEG    R_NEG NEGSNPS
    # 5          rs326150      7814593    G    A  0.00715      0       NA      3      nan rs326153|rs6883729|rs326162
    # 5          rs326153      7815809    C    T   0.0142      2    0.926      1      nan rs326150
    # 5         rs6883729      7819569    G    A   0.0123      2    0.886      1      nan rs326150
    # 5          rs326162      7819850    C    A   0.0146      2    0.941      1      nan rs326150
    #now that is interesting. Here it is not so clear cut; rs326153, rs6883729, rs326162 look okay and probably shouldn't be removed 

    #zgrep -a "rs326150\|rs326153\|rs6883729\|rs326162" /working/lab_stuartma/scratch/matthewL/melanoma_meta_analysis/cleaned/IMPUTATION_RESULTS/further_imputing/Archive_QMEGA_BEACON_chromosome_5_imputation.tar.gz  
      #--- rs326150 7814593 0.013 0.910 0.997 0 -1 -1 -1 <<< and this is the one that looks weird with the others. Pretty rare mind so would only need to be a couple of bad calls to mess up the LD
      #5 rs326153 7815809 0.014 1.000 1.000 2 0.873 0.997 0.928 <<well imputed all, and highly concordant. 
      #5 rs6883729 7819569 0.012 1.000 1.000 2 0.906 0.999 0.969
      #5 rs326162 7819850 0.015 0.982 0.999 2 0.924 0.998 0.964
 
    #awk '{print $1,$2}' ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG.fam > ${dirI}/temp_keep
    #plink_1.90 --threads 1 --bfile /working/lab_stuartma/scratch/matthewL/melanoma_meta_analysis/melanoma_GWAS_original/melanoma --keep ${dirI}/temp_keep --snps rs326150,rs326153,rs6883729,rs326162  --make-bed --out ${dirI}/plink
    #plink_1.90 --threads 1 --bfile ${dirI}/plink --ld rs326150 rs326153
    
   #   --ld rs326150 rs326153:
   #Multiple haplotype phasing solutions; sample size, HWE, or random mating assumption may be violated.
   #HWE exact test p-values
   #-----------------------
   #rs326150: 1
   #rs326153: 0.169135
   #Solution #1:
   #R-sq = 0.000174659    D' = 1
   #Haplotype     Frequency    Expectation under LE
   #---------     ---------    --------------------
   #       GC      0                       0.000170
   #       AC      0.013043                0.012873
   #       GT      0.013043                0.012873
   #       AT      0.973913                0.974083
   #In phase alleles are GT/AC

    #Solution #2:
   #R-sq = 1              D' = 1
   #Haplotype     Frequency    Expectation under LE
   #---------     ---------    --------------------
   #       GC      0.013043                0.000170
   #       AC      0                       0.012873
   #       GT      0                       0.012873
   #       AT      0.986957                0.974083
   #In phase alleles are GC/AT

   #plink_1.90 --threads 1 --bfile /working/lab_stuartma/scratch/matthewL/melanoma_meta_analysis/BEACON_BACKUP/beagess_australia_matthew --ld rs326150 rs326153 --make-founders
     #Error: rs326150 is monomorphic across all valid observations.
   
   #plink_1.90 --threads 1 --bfile ${dirI}/temp_2016_cleaning_QMEGA_omni_aligned_1KG --filter-controls --ld rs326150 rs326153
     #Error: rs326150 is monomorphic across all valid observations.

    #So I think I can safely drop this SNP, but the other three look okay, just due to rs326150 so kludge them back in

   #20/07/2016 the old error due to merge-multiple-position was found by testing ld for rs9402684 rs7745098 in 610k -shown here is that the LD should be but was getting TT and CC in phase as the alleles got swapped for one of the SNPs - see above but deleted most of the non-applicable notes no3w
    #plink_1.90  --bfile ${dirI}/${i} --filter-controls --ld rs9402684 rs7745098
         #           Rsq = 0.994915       D' = 0.99847
         #Haplotype     Frequency    Expectation under LE
    #   ---------     ---------    --------------------
     #     TC      0.486521                0.237322
     #     CC      0.000382                0.249581
     #     TT      0.000890                0.250089
     #     CT      0.512207                0.263008
    ###########################################################
    #WAMHS/OAG/IBD
    ###########################################################
    #  CHR                 SNP           BP   A1   A2        F    POS    R_POS    NEG    R_NEG NEGSNPS
    # 3           rs6763046     69220207    A    C  0.00531     11    0.489      1    0.323 rs9818314
    # 3           rs4473559     69222863    T    G  0.00861      6    0.625      1    0.327 rs9818314
    # 3           rs9818314     69230026    T    C  0.00264      1    0.418      4    0.314 rs6763046|rs4473559|rs9831516|rs9836305 <<looks to be the problem rs9818314
    # 3           rs9831516     69230061    G    A   0.0093      7    0.566      1    0.291 rs9818314
    # 3           rs9836305     69230801    G    A    0.011      5    0.625      1    0.316 rs9818314
    #13           rs9586420    104855694    A    G  0.00258      4    0.561      1     0.26 rs2146645 <<can't tell which one it is, drop both
    #13          rs16964351    104867556    G    A  0.00243      5    0.516      1    0.318 rs7329271 <<can't tell which one it is, drop both
    #13           rs2146645    104877036    G    A   0.0076     10    0.581      1     0.26 rs9586420 <<can't tell which one it is, drop both
    #13           rs7329271    104888924    C    T  0.00645      8    0.578      1    0.318 rs16964351 <<can't tell which one it is, drop both
    #15           rs4776147     53709744    T    C   0.0113      3    0.557      1    0.266 rs1697452 <<can't tell which one it is, drop both
    #15           rs1697452     53719045    G    A  0.00917      4    0.563      1    0.266 rs4776147 <<can't tell which one it is, drop both
    #15           rs7171438    101401700    T    C  0.00157      1    0.254      1    0.257 rs7182439 <<can't tell which one it is, drop both
    #15           rs7182439    101402323    G    T   0.0123      4      0.7      1    0.257 rs7171438 <<can't tell which one it is, drop both
    #19          rs10402010     22086049    G    T  0.00183      1      0.5      1    0.282 rs10424424 <<can't tell which one it is, drop both
    #19          rs10424424     22110219    T    C   0.0103      2    0.614      1    0.282 rs10402010 <<can't tell which one it is, drop both
    #20           rs7261749     17731502    T    C  0.00301      0       NA      1    0.335 rs6034900 <<can't tell which one it is, drop both
    #20           rs6034900     17746155    T    C  0.00136      0       NA      1    0.335 rs7261749 <<can't tell which one it is, drop both
    #20           rs6023927     53781235    C    T    0.002      1    0.463      1     0.26 rs6092078 <<can't tell which one it is, drop both
    #20           rs6092078     53801129    A    G  0.00818      1    0.959      1     0.26 rs6023927 <<can't tell which one it is, drop both
    #23           rs7049877     17610299    A    G  0.00668      3    0.476      1    0.281 rs16980660
    #23          rs16980660     17648165    T    C  0.00407      0       NA      1    0.281 rs7049877
    #23          rs17148207     48329702    C    T   0.0024      2    0.466      1    0.336 rs7066831
    #23           rs7066831     48341462    C    T  0.00219      0       NA      1    0.336 rs17148207
    #23          rs12006821     50424688    T    C  0.00174      1    0.408      1    0.289 rs17003186
    #23          rs17003186     50449578    T    C  0.00145      1    0.353      4    0.339 rs12006821|rs4142629|rs5915311|rs6521881
    #23           rs4142629     50462442    C    T  0.00268      5    0.679      1    0.289 rs17003186
    #23           rs5915311     50502078    A    C   0.0029      3    0.811      1    0.388 rs17003186
    #23           rs6521881     50507516    C    T   0.0029      3    0.811      1    0.388 rs17003186
    #23           rs6568263    108778290    T    C  0.00239      1    0.288      1    0.285 rs1924829
    #23           rs1924829    108792848    A    C  0.00406      1    0.267      1    0.285 rs6568263
    #23           rs5980580    147778631    T    C  0.00145      0       NA      1    0.307 rs16994698
    #23          rs16994698    147828662    G    A   0.0056      0       NA      1    0.307 rs5980580
    #26          exm2216447        14318    C    T  0.00182      1    0.353      1    0.354 exm2216458
    #26          exm2216458        14978    G    A  0.00182      0       NA      1    0.354 exm2216447  
   #zgrep -a rs9818314 /working/lab_stuartma/scratch/matthewL/melanoma_meta_analysis/cleaned/IMPUTATION_RESULTS/CIDR_WAMHS_IMPUTATION/Archive_CIDR_WAMHS_chromosome_3_imputation.tar.gz
    #not in the imputation?
  #zgrep -w rs9818314 ${dirI}/ALL_1000G_phase1integrated_v3_ALL_chr_legends_with_chromosome.txt.gz
    #3 rs9818314 69230026 T C 0.7358 0.9945 0.9283 1.0000 0.2642 0.0055 0.0717 0.0000
  #grep rs9818314 /working/lab_stuartma/scratch/matthewL/melanoma_meta_analysis/cleaned/WAMHS_GLUACOMA_IBD_BEACON_final_clean_IBD_PCA_outliers_removed.bim
    #not in there 

  #for i in ${dirOAG}/${OAG1} ${dirOAG}/${OAG2} ${dirJI}/${raw_data4} ${dirJW}/${raw_data5}
  #do
    #plink_1.90 --threads 1 --bfile ${i} --ld rs6763046 rs9818314
    #echo ""
    #plink_1.90 --threads 1 --bfile ${i} --ld rs6763046 rs4473559
    #echo
    #plink_1.90 --threads 1 --bfile ${i} --ld rs9831516 rs9836305
  #done

  #OAG phase 1 
     #one or the other not in there - which is probably why this never came up before
    #grep rs9818314 ${dirOAG}/${OAG1}.bim
       #not in there
   #OAG phase 2
     #--ld rs6763046 rs9818314:
     # R-sq = 0.273818       D' = 0.830297
     # Haplotype     Frequency    Expectation under LE
     # ---------     ---------    --------------------
     #        AA      0.003883                0.000054
     #        CA      0.000783                0.004611
     #        AG      0.007781                0.011610
     #        CG      0.987553                0.983725
     #In phase alleles are AA/CG
  #IBD
    # --ld rs6763046 rs9818314:
    # R-sq = 0.498431       D' = 1
    #Haplotype     Frequency    Expectation under LE
    #---------     ---------    --------------------
    #      AT      0.003128                0.000020
    #      CT      0                       0.003109
    #      AC      0.003128                0.006237
    #      CC      0.993743                0.990635
    #In phase alleles are AT/CC
  #WAMHS 
   #--ld rs6763046 rs9818314:
   #R-sq = 3.09635e-06    D' = 1
   #Haplotype     Frequency    Expectation under LE
   #---------     ---------    --------------------
   #       AA      0                       0.000003
   #       CA      0.001964                0.001961
   #       AG      0.001571                0.001568
   #       CG      0.996465                0.996468
   #In phase alleles are AG/CA  #hah ha what. << dif in each set, but see 1kg - monomorphic, so need to drop rs9818314. Rest might be alright

#for the LD estimates differ (one SNP is sig rare than the other) but for rs6763046 rs4473559 In phase alleles are AA/CC (note second SNP flipped in OAG2)
  #thousand genomes  
  #plink_1.90 --threads 1 --bfile /reference/genepi/public_reference_panels/1000G_20101123_v3/derived_plinkformat/chr3.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR --ld rs6763046 rs9818314
    #Error: rs9818314 is monomorphic across all valid observations. Would suggest the rs9818314 is producing nonsense and should be removed.

  #for i in ${dirOAG}/${OAG1} ${dirOAG}/${OAG2} ${dirJI}/${raw_data4} ${dirJW}/${raw_data5} 
  #do
    #plink_1.90 --threads 1 --bfile ${i} --ld rs6763046 rs9818314
    #echo ""
    #plink_1.90 --threads 1 --bfile ${i} --ld rs6763046 rs4473559
    #echo
    #plink_1.90 --threads 1 --bfile ${i} --ld rs9831516 rs9836305
    #plink_1.90 --threads 1 --bfile ${i} --ld rs17003186 rs12006821
    #WAMHS  R-sq = 2.76241e-06    D' = 1    In phase alleles are AG/GA <<< 
    #IBD    R-sq = 0.443399       D' = 0.665882    In phase alleles are TT/CC
    #OAG2    R-sq = 0.39907        D' = 1    In phase alleles are AA/GG
    #OAG1 one or more is missing.
    #plink_1.90 --threads 1 --bfile ${i} --ld rs9831516 rs4473559
  #done

  #plink_1.90 --threads 1 --bfile /reference/genepi/public_reference_panels/1000G_20101123_v3/derived_plinkformat/chr3.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR --ld rs17003186 rs4142629
    #no chr X in 1000G_20101123_v3, need giant with the EUR set
  #awk '{print $1,$2}' /reference/genepi/public_reference_panels/1000G_20101123_v3/derived_plinkformat/chr3.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.fam > ${dirI}/temp_keep
  #plink_1.90 --threads 1 --bfile /reference/genepi/public_reference_panels/1000G_20101123_v3_GIANT/derived_plinkformat/chrX.no.auto.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL --keep ${dirI}/temp_keep --ld rs17003186 rs12006821
   #Error: rs4142629 is monomorphic across all valid observations.rs5915311 is monomorphic across all valid observations. rs6521881 is mono as is rs12006821
  ##################################################################
  #HEIDELBERG set ${dirJG}/${raw_data3}
  ####################################################################
  # CHR         SNP           BP   A1   A2        F    POS    R_POS    NEG    R_NEG NEGSNPS
  #  11   rs2011233     19638535    G    A  0.00144      0       NA      1    0.266 rs714585
  #  11    rs714585     19668438    G    A  0.00349      0       NA      1    0.266 rs2011233
  #Drop both


########################################################################
#6.0 QUICK GWAS for major problems
########################################################################
 #run a quick GWAS to look for any weird problems before cleaning properly
for x in #1
do
  for i in #temp_2016_cleaning_WAMHS_OAG_IBD_aligned_1KG temp_2016_cleaning_AMFS_aligned_1KG temp_2016_cleaning_QMEGA_610k_aligned_1KG temp_2016_cleaning_QMEGA_omni_aligned_1KG temp_2016_cleaning_EPIGENE_QTWIN_aligned_1KG temp_2016_cleaning_HEIDELBERG_aligned_1KG
  do
    echo -e "\nquick GWAS for ${i}\n"
    plink_1.90 --threads 1  --bfile ${dirI}/${i} --assoc --out ${dirI}/${i}
    echo "near GWS hits"
    awk 'NR==1 || $9<5e-7' ${dirI}/${i}.assoc
    rm -f ${dirI}/${i}.assoc
  done #loop through sets
done #loop to turn on/off

  #AMFS 
  # CHR               SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR 
  #16          rs258322   89755903    A   0.1734  0.08256    G        34.36    4.584e-09         2.33 
   #doesn't suggest any major problems

  #QMEGA_610k
  #CHR         SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR 
  # 2  rs17278676  124849389    G   0.3692   0.2985    T        34.54    4.178e-09        1.376 <<had same OR in previous 610k data, just chance
  # 3  rs10048973   20310573    T   0.3598   0.4328    C        31.88    1.643e-08       0.7365 << cleaned out of the other 610k data
  # 4   rs7682929   44297909    A   0.1141  0.07764    C        25.68    4.027e-07        1.529 <<slightly worse p value in old data looks like a couple of people carrying alleles removed
  # 4   rs7663379  171634933    C   0.2616   0.3328    T         33.8    6.115e-09       0.7101 << cleaned out of the other 610k data
  # 7  rs17157992   83349707    T  0.05493  0.02314    C        52.26    4.865e-13        2.454 << cleaned out of the other 610k data
  # 8   rs1440788   53388439    A   0.1405  0.08749    G        35.97        2e-09        1.704 << cleaned out of the other 610k data
  #15   rs2305252   28267338    T   0.3027   0.2356    C        35.09    3.155e-09        1.408  << cleaned out of the other 610k data
  #16    rs352935   89648580    T   0.4313   0.5073    C        34.67     3.91e-09       0.7364   <<MC1R
  #16    rs164741   89692298    A   0.3755   0.3134    G        26.33    2.878e-07        1.317 
  #16   rs7188458   89726484    A   0.5227   0.4503    G        31.61    1.887e-08        1.337 
  #16    rs258322   89755903    A   0.1676     0.11    G        46.72    8.209e-12        1.628 
  #16   rs1800359   89805261    A   0.3341   0.4006    G        27.97    1.235e-07       0.7505 
  #16   rs1800286   89869761    T   0.3341   0.4012    C        28.49    9.436e-08       0.7486 
  #16  rs11861084   89875710    A   0.3378   0.4047    C        28.08    1.162e-07       0.7506 
  #16   rs8049897   90024202    A   0.2286   0.1624    G        45.64    1.425e-11        1.529 
  #16   rs4408545   90044028    T   0.4065   0.4917    C         43.6     4.02e-11       0.7081 
  #16   rs4238833   90050689    G    0.447   0.3818    T        26.71    2.369e-07        1.309 
  #16   rs4785763   90066936    A   0.4114   0.3463    C         27.6    1.494e-07        1.319 

  #nothing leaps out at me as weird here

  #QMEGA_omni
  # CHR               SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR 
  # 9        rs10811595   21715801    G   0.4094   0.5185    T        29.86     4.64e-08       0.6439 <<CDKN2A
  # 9         rs6475564   21728683    T   0.5447    0.443    C         25.9    3.588e-07        1.505
  # 9         rs7871345   21729895    T    0.544    0.443    C        25.54    4.337e-07          1.5
  # 9         rs7854222   21734348    G   0.5441   0.4428    A        25.61    4.187e-07        1.502
  # 9        rs10811600   21734802    A    0.544    0.443    G        25.54    4.337e-07          1.5
  # 9         rs4636294   21747803    A    0.544   0.4421    G        25.98    3.443e-07        1.506
  # 9         rs1965153   21752105    T   0.5447   0.4421    G        26.35    2.844e-07         1.51
  # 9         rs7852710   21760254    T   0.5441   0.4421    C           26    3.415e-07        1.506
  # 9         rs7866787   21760639    A    0.544   0.4421    G        25.98    3.443e-07        1.506
  # 9         rs7856941   21761322    T    0.544    0.442    C        26.01    3.395e-07        1.506
  # 9         rs7850446   21763742    G   0.5455   0.4421    A        26.72    2.346e-07        1.514
  # 9         rs6475576   21765601    C   0.5447   0.4421    T        26.35    2.844e-07         1.51
  # 9         rs3849929   21769412    G   0.5448   0.4421    A        26.37     2.82e-07         1.51
  # 9         rs7045768   21770709    T    0.544   0.4421    C        25.98    3.443e-07        1.506
  # 9         rs6475579   21771756    A    0.544    0.443    G        25.54    4.337e-07          1.5
  # 9         rs7037577   21772037    G    0.544    0.442    A        26.01    3.395e-07        1.506
  # 9        rs10757253   21772267    T    0.544    0.443    C        25.54    4.337e-07          1.5
  # 9         rs1414228   21773541    A    0.542     0.44    C         25.9    3.589e-07        1.506
  # 9         rs7038708   21788081    T   0.5397   0.4386    G        25.58    4.238e-07        1.501
  # 9        rs10811617   21790067    A    0.539   0.4377    G        25.65    4.092e-07        1.502
  #16         rs7195043   90020861    C   0.4488   0.5526    T        26.99    2.045e-07       0.6591 <<MC1R
  #20         rs6059655   32665748    A   0.1522  0.08509    G        26.34    2.858e-07        1.931 <<ASIP, literally an ASIP color SNP

  #EPIGENE TQWIN
 #  CHR                   SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR
 #  1             rs9662633    3649562    A 0.001923  0.03619    G         49.6    1.883e-12      0.05132
 #  1            rs41306613   33235789    T 0.001276  0.02905    C        40.91    1.591e-10      0.04268
 #  2           rs143730975   25523096    G   0.9981 0.003568    T         3490            0    1.457e+05 <<cleaned
 #  2             rs4670855   38504680    A 0.001276   0.1998    G        346.4    2.597e-77     0.005115 
 #  3            rs11711450  178207097    T 0.002551   0.6391    C         1540            0     0.001444
 #  4             rs4493574   64073599    A 0.001913   0.0316    G         42.5    7.086e-11      0.05874
 #  6           rs114206289   29431698    G 0.001276   0.5757    A         1319    7.05e-289    0.0009411
 #  6             rs9443079   75553805    T 0.001276   0.1435    C        238.6     7.79e-54     0.007622
 #  6            rs61731700  109771691    A 0.001279  0.01937    G        25.39    4.682e-07      0.06483
 #  6             rs9403839  147282434    T 0.001913  0.02497    C        31.93    1.595e-08      0.07484
 #  8            rs13265018   22886020    A   0.0102   0.1081    G        138.1    6.901e-32       0.0851
 #  8            rs10107655   57994818    T 0.001276   0.0841    G        132.6    1.084e-30      0.01391
 #  9             exm779484  125486542    G 0.001276  0.02345    A        31.94    1.593e-08       0.0532
 # 10            rs61866610   95326676    T 0.003831  0.04638    C        58.97    1.602e-14      0.07908
 # 10              rs942674  131922572    A 0.001928     0.37    C        715.8    1.12e-157     0.003289
 # 13        chr13:37584732   37584732    C   0.9962 0.002041    A         3488            0    1.273e+05
 # 14             rs2415664   41577085    G 0.001913    0.402    A        800.9    3.37e-176     0.002852
 # 15            exm1142547   25924865    C   0.9962 0.003058    T         3476            0    8.454e+04
 # 16            rs74017766   50703016    T  0.02934 0.007136    C        25.71    3.969e-07        4.205
 # 17            exm1345519   62028843    T  0.01722 0.001019    C        28.07    1.169e-07        17.17
 # 20              rs910873   33171772    A   0.1422  0.08563    G        28.34    1.019e-07         1.77 <<ASIP
 # 21             rs7276462   31429636    T 0.001913  0.05505    C        80.78    2.519e-19      0.03291 
 # 23             rs2799952   93634938    G   0.0126  0.05574    A        31.13    2.414e-08       0.2161
 # 23             rs1174258  144642152    G  0.00194  0.03435    A        31.04    2.523e-08      0.05464

 #9 missing based on IDs (IDs changed during alignment. All are unique positional SNPs (so I haven't merged SNPs etc)
 #grep "33235789\|25523096\|29431698\|109771691\|22886020\|95326676\|37584732\|50703016\|33171772" temp_2016_cleaning_QTWIN_controls.bim
   #1	exm41165	0	33235789	T	C
   #2	exm177625	0	25523096	G	T
   #6	exm-rs3128854	0	29431698	A	G
   #6	exm570865	0	109771691	A	G
   #8	exm689164	0	22886020	A	G
   #10	exm842865	0	95326676	T	C
   #13	exm1063738	0	37584732	C	A
   #16	exm1239770	0	50703016	T	C
   #20	exm-rs910873	0	33171772	A	G

 #grep "33235789\|25523096\|29431698\|109771691\|22886020\|95326676\|37584732\|50703016\|33171772" temp_2016_cleaning_EPIGENE_cases.bim
   #1	exm41165	0	33235789	T	C
   #2	exm177625	0	25523096	T	G
   #6	exm-rs3128854	0	29431698	G	A
   #6	exm570865	0	109771691	A	G
   #8	exm689164	0	22886020	A	G
   #10	exm842865	0	95326676	T	C
   #13	exm1063738	0	37584732	A	C
   #16	exm1239770	0	50703016	T	C
   #20	exm-rs910873	0	33171772	A	G

  #need to work out if any were caused by my merging.
  #awk 'NR==1 || $9<5e-7 {print $2}' ${dirI}/${i}.assoc > ${dirI}/${i}.assoc_snps; wc -l ${dirI}/${i}.assoc_snps
  #grep "33235789\|25523096\|29431698\|109771691\|22886020\|95326676\|37584732\|50703016\|33171772" temp_2016_cleaning_EPIGENE_cases.bim | awk '{print $2}' >> ${dirI}/${i}.assoc_snps; wc -l ${dirI}/${i}.assoc_snps
   #33 (some in there twice under two IDs
  #plink_1.90 --threads 1 --bfile temp_2016_cleaning_QTWIN_controls --extract ${dirI}/${i}.assoc_snps --freq
  # CHR         SNP   A1   A2          MAF  NCHROBS
  # 1   rs9662633    A    G      0.03619     1934 <<same in QTWIN and EPI original data
  # 1        exm41165    T    C      0.02947     1934 <<same in QTWIN and EPI original data
  #   2       exm177625    G    T     0.003619     1934 <<same in QTWIN and EPI original data (as in really is minor-major in the other array.
  # 2   rs4670855    A    G        0.198     1934 <<same in QTWIN and EPI original data
  # 3  rs11711450    C    T       0.3594     1934 <<same in QTWIN original data NA in EPI and now has a freq
  # 4   rs4493574    A    G      0.03154     1934 <<same in QTWIN and EPI original data
  # 6   exm-rs3128854    A    G       0.4221     1926 << really is ~0.42 here and 0.001 in EPI
  # 6   rs9443079    T    C       0.1435     1930 <<same in QTWIN and EPI original data
  # 6       exm570865    A    G      0.01965     1934 <<same in QTWIN and EPI original data
  # 6   rs9403839    T    C      0.02534     1934 <<same in QTWIN and EPI original data
  # 8       exm689164    A    G       0.1075     1934 <<same in QTWIN and EPI original data
  # 8  rs10107655    T    G       0.0848     1934 <<same in QTWIN and EPI original data
  # 9   exm779484    G    A      0.02327     1934 <<same in QTWIN and EPI original data
  # 10       exm842865    T    C      0.04654     1934 <<same in QTWIN and EPI original dat
  #10    rs942674    A    C       0.3713     1934 <<same in QTWIN and EPI original data
  # 13      exm1063738    C    A      0.00207     1932 <<same in QTWIN and EPI original data - major and minor really are swapped
  #14   rs2415664    G    A       0.4028     1924 <<same in QTWIN and EPI original data
  #15  exm1142547    C    T     0.003102     1934 <<same in QTWIN and EPI original data
  #17  exm1345519    T    C     0.001034     1934 <<same in QTWIN and EPI original data
  # 20    exm-rs910873    A    G       0.0848     1934 <<real locus
  #21   rs7276462    T    C      0.05533     1934 <<same in QTWIN and EPI original data
  #23   rs2799952    G    A      0.05585     1522 <<same in QTWIN and EPI original data
  #23   rs1174258    G    A      0.03351     1522 <<same in QTWIN and EPI original data
  #15 only
  #plink_1.90 --threads 1 --bfile temp_2016_cleaning_EPIGENE_cases --extract ${dirI}/${i}.assoc_snps --freq
  #15 only
  # CHR         SNP   A1   A2          MAF  NCHROBS
  # 1   rs9662633    A    G     0.001923     1560
  # 1        exm41165    T    C     0.001276     1568
  # 2       exm177625    T    G     0.001913     1568 
  # 2   rs4670855    A    G     0.001276     1568
  # 3  rs11711450    C    T           NA        0 <<?
  # 4   rs4493574    A    G     0.001913     1568 
  # 6   exm-rs3128854    G    A     0.001276     1568 
  # 6   rs9443079    T    C     0.001276     1568
  # 6       exm570865    A    G     0.001279     1564
  # 6   rs9403839    T    C     0.001913     1568
  # 8       exm689164    A    G       0.0102     1568
  # 8  rs10107655    T    G     0.001276     1568
  # 9   exm779484    G    A     0.001276     1568
  # 10       exm842865    T    C     0.003831     1566
  #10    rs942674    A    C     0.001928     1556
  # 13      exm1063738    A    C     0.003827     1568
  #14   rs2415664    G    A     0.001913     1568
  #15  exm1142547    T    C     0.003841     1562
  # 16      exm1239770    T    C      0.02934     1568
  #17  exm1345519    T    C      0.01722     1568
  # 20    exm-rs910873    A    G       0.1422     1568
  #21   rs7276462    T    C     0.001913     1568
  #23   rs2799952    G    A       0.0126     1032
  #23   rs1174258    G    A      0.00194     1031

  #SO EPIGENE - QTWIN - no systematic differences but a bunch of concerning SNPs



 #HEI
 #CHR         SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR
 #  1   rs3795501   36180972    T  0.02099  0.04899    C        28.07    1.168e-07       0.4161 
 #  2  rs17774013   47067214    T   0.2806  0.03654    C        543.7   2.924e-120        10.29 <<< removed in last imp
 #  3   rs3774262  186571814    A   0.1253  0.07869    G        28.92    7.542e-08        1.677 <<present in last run. Less extreme post cleaning, and in the meta set while OR was 1.28 it was almost as extreme in other sets (1.16 in QMEGA beacon)
 #  5     rs31489    1342714    A   0.4663   0.3946    C        25.56    4.291e-07        1.341 
 #  7    rs725092  125467854    G    0.514   0.4418    T        25.43     4.58e-07        1.336
 #  8    rs870585    5615643    A   0.1499  0.08794    C        44.19    2.979e-11        1.828  <<< removed in last imp
 #  9    rs644893   79293847    C  0.02054 0.004098    T         27.1    1.933e-07        5.096
 #  9   rs2031852  129755509    A   0.2679   0.2056    G        26.11     3.22e-07        1.414
 # 13    rs283957   77117462    T   0.4337   0.3624    C        25.81    3.772e-07        1.348
 # 15   rs7172696  101991311    A   0.2814   0.3539    G         29.1     6.89e-08       0.7151 <<< removed in last imp
  
 #this suggests a lot of the extreme diff in control vs HEI control differences might be actual ethnic differences (suggested by preponderance of HLA hits)

 #Other than monitoring EPI/QTWIN post cleaning nothing super concerning, and if I can't remove some of the extreme SNPs through cleaning imputation might shed light

 #WAMHS vs OAG/IBD - decent amount, filtering down to <8e08 or so for clarity 
 #CHR                 SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR
 #  1          rs10915285   29782051    A 0.001181  0.01422    G        29.27    6.304e-08      0.08197 <<not in prev imp
 #  1          rs16834715  152844585    T 0.001965    0.016    C        29.58    5.378e-08       0.1211 <<not in prev imp
 #  2           rs6730157  135907088    G   0.2675   0.3327    A        32.23    1.369e-08       0.7323 <<
 #  2           rs1900741  136002500    T   0.1313   0.1817    C        30.06    4.186e-08       0.6807
 #  2          rs12465802  136381348    A    0.262   0.3256    G        31.04    2.529e-08       0.7353 <<none of these chrom 2 SNPs eem to be in the final set. SNPs have been associated with BMI at e-6
 #  2            rs313518  136439090    T    0.188   0.2452    C        30.28    3.732e-08       0.7131
 #  2           rs1438307  136499166    G   0.1843   0.2397    T        28.84    7.853e-08       0.7163
 #  2           rs9287442  136522710    A   0.1651   0.2189    G        29.32    6.126e-08       0.7056
 #  2           rs4988235  136608646    G   0.2553   0.3186    A        31.15    2.392e-08       0.7332
 #  2            rs182549  136616754    C    0.249   0.3146    T        33.72    6.372e-09       0.7224
 #  3           rs2700372  123633220    G 0.007492   0.0281    T        32.12    1.448e-08       0.2611 <<not in prev imp, not in meta at al
 #  4          rs10033656   19607089    T   0.2467    0.314    C        28.74    8.265e-08       0.7153 
 #  4          rs13104573   77225313    C   0.1517   0.2112    T        37.21    1.061e-09        0.668 << in prev imp. Extreme in WAMHS> High het SNP in meta... not sure. 
 #  4           rs7658585   77252196    T   0.1386   0.1955    C        36.29    1.704e-09       0.6622
 #  5             rs35407   33946571    A  0.01257  0.05385    G        56.79     4.84e-14       0.2237 SLC45A2 SNP
 #  5            rs250416   33947544    A 0.005106  0.02298    C        31.93    1.601e-08       0.2182
 #  5             rs35391   33955673    T  0.01021  0.03969    C        50.02    1.522e-12       0.2496
 #  5             rs28777   33958959    C  0.01181  0.05308    A         58.3    2.251e-14       0.2132
 #  5             rs28117   33962770    G  0.01453  0.04304    A        41.63    1.105e-10       0.3279
 #  5            rs183671   33964210    T  0.01414  0.05478    G         52.4    4.524e-13       0.2475
 #  6           rs2256965   31555130    A   0.4493   0.3837    G        28.88    7.719e-08        1.311  
 #  6           rs9498740  102341430    A   0.1354  0.08521    G        36.67    1.398e-09        1.682 <<not in previous imp, much less extreme following imp ((OR ~1.1) 
 # 10          rs11816696   25362890    C   0.2913   0.1712    T        64.41     1.01e-15         1.99 <<not in previous imp, much less extreme following imp OR 1.06
 # 12          rs11066009  112140815    A 0.001967  0.01633    G        29.53    5.515e-08       0.1188 
 # 12           rs2339840  112254941    T 0.001969  0.01625    C        29.33    6.108e-08       0.1194
 # 12          rs34708049  132879672    C   0.5622   0.1764    T        819.6   2.962e-180        5.996 < not in previous imp. 
 # 14          rs11157646   49198386    C 0.001574  0.01533    T        29.22     6.45e-08       0.1013
 # 15           rs1129038   28356859    C   0.2054    0.263    T         29.2    6.533e-08       0.7246 
  #bunch of MC1R GWS hits not showing
 # 16           rs1805007   89986117    T   0.1546   0.0836    C        49.53    1.957e-12        2.005 <<MC1R
 # 18           rs2578189   43680700    T    0.132  0.06641    C        37.39    9.668e-10        2.137 < not in previous imp.
 # 22            rs138000   40354981    T   0.1081    0.313    G        370.6    1.388e-82        0.266 < not in previous imp.
    #again nothing to concerning, other than watching them post proper cleaning, but looks like most will be filtered out.
  

########################################################################
#7.0 Basic cleaning 
########################################################################
#3. Clean on missingness (--mind 0.03 --geno 0.03 --maf 0.001 ) then HWE (10-4 in controls, 10-10 in cases). In this case don't do the maf cleaning yet. do that later
#4 Also filter on extreme het in samples (see the MIA scripts etc; MI suggests >0.05 or < -0.05)

for x in #1
do

  for i in #temp_2016_cleaning_WAMHS_OAG_IBD_aligned_1KG temp_2016_cleaning_AMFS_aligned_1KG temp_2016_cleaning_QMEGA_610k_aligned_1KG temp_2016_cleaning_QMEGA_omni_aligned_1KG temp_2016_cleaning_EPIGENE_QTWIN_aligned_1KG temp_2016_cleaning_HEIDELBERG_aligned_1KG
  do

    #need to filter any SNPs missing in one array or the other for merged sets
    plink_1.90 --threads 1 --bfile ${dirI}/${i} --geno 0.2 --make-bed --out ${dirI}/${i}_geno02
      #610k 550293 variants loaded from .bim file. 47382 variants removed due to missing genotype data (--geno). 502911 variants and 4878 people pass filters and QC.
      #EPI/QTWIN - 323811 variants loaded from .bim file. 38319 variants removed due to --geno 285492 variants and 1765 people pass filters and QC.
      #AMFS 801407 variants loaded from .bim file. 3230 variants removed due to missing genotype data (--geno). 798177 variants and 978 people pass filters and QC.
      #QMEGA_Omni  848795 variants loaded 62932 variants removed due to --geno 785863 variants and 1264 people pass filters and QC.
      #HEI 217285 variants loaded 4217 variants removed due to missing genotype data 213068 variants and 2437 people pass filters and QC.
      #WAMHS 929087 variants loaded 293235 variants removed due to --geno 635852 variants and 3492 people pass filters and QC. <<doesn't seem to be right for WMAHS, next stip is removing all the OAG1 samples. 

    if `echo ${i} | grep -q WAMHS `
    then
      echo -e "\nWAMHS data found; need an diff cleaning step or lose all OAG1 sample due to smaller 610 array\n"
        #3500 people, 650 0AG1 < 20% 
      plink_1.90 --threads 1 --bfile ${dirI}/${i} --geno 0.1 --make-bed --out ${dirI}/${i}_geno02
        #929087 variants loaded from .bim file. 337363 variants removed due to --geno 591724 variants and 3492 people pass filters and QC.

    fi
   
    data_in=${dirI}/${i}_geno02
    plink_1.90 --threads 1 --bfile ${data_in} --missing --out ${data_in}_miss --silent
    plink_1.90 --threads 1 --bfile ${data_in} --het --out ${data_in}_het --silent

  #Following is based off the GWAS cleaning paper UY has been using; is a nice method so I see no need to remove it
  cd $dirI

  module load R/3.2.2

R --no-save <<EOF

#read in the per person (individual) missingness data = imiss. h=T means there is a header - the first line is column names. dim is dimensions; show us the size of the newly made object
imiss=read.table("${data_in}_miss.imiss",h=T);dim(imiss)

#look at the file
head(imiss)

#plotting will be easier on a log scale
imiss\$logF_MISS = log10(imiss[,6])

het=read.table("${data_in}_het.het",h=T);dim(het)
het\$meanHet = (het\$N.NM. - het\$O.HOM.)/het\$N.NM.
#I think geneplotter doesn't work in R anymore so exit out
#library("geneplotter")
colors  <- densCols(imiss\$logF_MISS,het\$meanHet)
#non fatal warning: Warning message: In KernSmooth::bkde2D(x, bandwidth = bandwidth, gridsize = nbin,  :  Binning grid too coarse for current (small) bandwidth: consider increasing 'gridsize'

pdf("${data_in}_imiss-vs-het.pdf")
plot(imiss\$logF_MISS,het\$meanHet, col=colors, xlim=c(-3,0),ylim=c(0,0.5),pch=20, xlab="Proportion of missing genotypes", ylab="Heterozygosity rate",axes=F)
axis(2,at=c(0,0.05,0.10,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5),tick=T)
axis(1,at=c(-3,-2,-1,0),labels=c(0.001,0.01,0.1,1))

abline(h=mean(het\$meanHet)-(2*sd(het\$meanHet)),col="RED",lty=2)
abline(h=mean(het\$meanHet)+(2*sd(het\$meanHet)),col="RED",lty=2)
#drawing 2sd not 3sd
abline(h=mean(het\$meanHet)-(3*sd(het\$meanHet)),col="GREEN",lty=2)
abline(h=mean(het\$meanHet)+(3*sd(het\$meanHet)),col="GREEN",lty=2)
abline(v=-1.522879, col="RED", lty=2)
dev.off()

#write out IDs that have high missingness or >3sd het
imiss[imiss\$F_MISS>=0.03,c(1,2)]->imiss_missing;dim(imiss_missing)
  #610k 0, EPI/QTWIN 0, AMFS 0, QMEGA_omni 0, HEI 0. WAMHS 0 (was all 650 OAG1 at 0.2)

write.table(imiss_missing,"${data_in}_miss.imiss.missing",row.names=F,quote=F)
 
#write out IDs with meanHet more than 3 sd from the mean
het_min<-mean(het\$meanHet)-(3*sd(het\$meanHet))
het_max<-mean(het\$meanHet)+(3*sd(het\$meanHet))
het[het\$meanHet<het_min | het\$meanHet>het_max,c(1,2)]->het_extreme;dim(het_extreme)
  #610k 27, most look to be just outside the boundaries. QTWIN/EPI 37. AMFS 7. QMEGA_omni 12. HEI 18. WAMHS 68

write.table(het_extreme,"${data_in}_het.het.extreme",row.names=F,quote=F)

q()
EOF
     echo ""
   cat ${data_in}_het.het.extreme ${data_in}_miss.imiss.missing | sort | uniq | grep -v "FID IID" > ${data_in}_remove; wc -l ${data_in}_remove
      #610k 27. EPI/QTWIN 37. Plot looks okay as in doesn't seem to be that many outliers but all the values are in the range 0.28-0.36 so would have been plotted (I was concerned some extreme values might not be plotted, but not the case. AMFS 7. QMEGA_omni 12 HEI 18. WAMHS 68
    echo ""
    #now clean properly - hwe first at 5e-4, then 5e-10 in cases
    plink_1.90 --threads 1 --bfile ${data_in} --geno 0.03 --mind 0.03 --hwe 5e-4 midp --make-bed --out ${data_in}_geno3_mind3_hweCo --remove ${data_in}_remove
       #610k: 502911 variants loaded from .bim file. 4878 people (1271 males, 3607 females) loaded from .fam. --remove: 4851 people remaining. 0 people removed due to missing genotype data (--mind). 2246 variants removed due to missing genotype data (--geno). --hwe: 410 variants removed. 500255 variants and 4851 people pass filters and QC.
      #EPIGENE/QTWIN 285492 variants loaded from .bim file. 1765 people (939 males, 826 females) loaded from .fam. --remove: 1728 people remaining. 0 people removed due to --mind; 0 variants removed due to --geno --hwe: 206 variants removed. Total genotyping rate in remaining samples is 0.999821. 285286 variants and 1728 people pass filters and QC.
      #AMFS 798177 variants loaded 978 people (384 males, 594 females) loaded --remove: 971 people remaining. 0 people removed due to --mind 0 variants removed due to --geno --hwe: 423 variants removed. 797754 variants and 971 people pass filters and QC.
      #QMEGA_omni 785863 variants loaded 1264 people (795 males, 469 females) loaded --remove: 1252 people remaining. 0 people removed due to --mind 0 variants removed due to --geno --hwe: 533 variants removed 785330 variants and 1252 people pass filters and QC.
      #HEI 213068 variants loaded 2437 people (1316 males, 1121 females) loaded  --remove: 2419 people remaining 0 people removed due to --mind 0 variants removed due to --geno --hwe: 206 variants removed  212862 variants and 2419 people pass filters and QC.
      #WAMHS 591724 variants loaded 3492 people (1801 males, 1691 females) loaded --remove: 3424 people remaining. 0 people removed due to --mind 0 variants removed due to --geno --hwe: 1269 variants removed 590455 variants and 3424 people pass filters and QC
    echo ""
    plink_1.90 --threads 1 --bfile ${data_in}_geno3_mind3_hweCo --hwe 5e-10 include-nonctrl midp --make-bed --out ${data_in}_geno3_mind3_hweCo_hweCa
      #610k --hwe: 0 variant removed due to Hardy-Weinberg exact test. 500255 variants and 4851 people pass filters and QC.
      #EPIGENE/QTWIN --hwe: 7 variants removed 285279 variants and 1728 people pass filters and QC.
      #AMFS --hwe: 1 variant removed 797753 variants and 971 people pass filters and QC
      #QMEGA_omni --hwe: 0 variants removed 785330 variants and 1252 people pass filters and QC.
      #HEI --hwe: 1 variant removed 212861 variants and 2419 people pass filters and QC.
      #WAMHS --hwe: 3 variants removed 590452 variants and 3424 people pass filters and QC.

    #intermin clean
    rm -f ${data_in}_miss.* ${data_in}_het.* ${data_in}_remove ${data_in}.*

    #based on /mnt/lustre/home/matthewL/Melanoma/WHITEMAN/EPIGENE/SCRIPTS/script_to_identify_epigene_control_release8.sh  

    echo -e "\nLooking for SNPs with call rates that differ between genders\n"
    awk '{print $1,$2,$3,$4,$5,$5}' ${data_in}_geno3_mind3_hweCo_hweCa.fam > ${data_in}_geno3_mind3_hweCo_hweCa_genders_case.fam
    plink_1.90 --threads 1 --bed ${data_in}_geno3_mind3_hweCo_hweCa.bed --bim ${data_in}_geno3_mind3_hweCo_hweCa.bim --fam ${data_in}_geno3_mind3_hweCo_hweCa_genders_case.fam --test-missing midp perm --out ${data_in}_geno3_mind3_hweCo_hweCa
    echo ""   
    #not sure how I setlled on this p value - probably N SNPs in MIA
    awk 'NR==1 || $5<1.6e-7' ${data_in}_geno3_mind3_hweCo_hweCa.missing
      #EPI/QTWIN - none. AMFS Bunch on chrX with 0 missing in A and ~5% in controls. 610k - bunch 0 missing in 'controls' and ~3% in 'cases'. QMEGA_omni 3 on chrX with 0 missingness in 'unnafected' amd near 4% in 'cases'. HEI - 3 SNPs again 3% missing in 'controls' nad 0% in 'case'. WAMHS - all but 1 chrX and as the others - complete in ;Affected;, high missing (3%) in UNAF
    echo ""
    #again not sure on this p value
    awk '$3<5e-4 {print $2}' ${data_in}_geno3_mind3_hweCo_hweCa.missing.perm > ${dirI}/temp_SNP_diff_missing_by_gender; wc -l ${dirI}/temp_SNP_diff_missing_by_gender
      #EPIGENE/QTWIN 7. AMFS 373. 610k 1026. QMEGA_omni 73. HEI 32. WAMHS 185

    echo -e "\nLooking for SNPs with call rates that differ between cases and controls\n"
    plink_1.90 --threads 1 --bfile ${data_in}_geno3_mind3_hweCo_hweCa --test-missing midp perm --out ${data_in}_geno3_mind3_hweCo_hweCa
    echo ""
    #again no idea on these p values, presumably just analogous to what we use in HWE
    awk 'NR==1 || $3<5e-4' ${data_in}_geno3_mind3_hweCo_hweCa.missing.perm | wc -l
      #EPI/QTWIN - 389. AMFS 10, 610k 4082. QMEGA_omni 2674 HEI 4532 WAMHS 6889
    echo ""
    awk 'NR==1 || $5<5e-4' ${data_in}_geno3_mind3_hweCo_hweCa.missing | wc -l
       #EPI/QTWIN - 403. AMFS 10. 610k 4180. QMEGA_omni 2758 HEI 4542 WAMHS 7109
    echo ""
    awk '$3<5e-4 {print $2}' ${data_in}_geno3_mind3_hweCo_hweCa.missing.perm > ${dirI}/temp_SNP_diff_missing_by_status; wc -l ${dirI}/temp_SNP_diff_missing_by_status
       #EPI/QTIN 388 (drop header). AMFS 9. 410k 4081 QMEGA_omni 2673 HEI 4531 WAMHS 6888
    echo ""
    cat ${dirI}/temp_SNP_diff_missing_by_gender ${dirI}/temp_SNP_diff_missing_by_status | sort | uniq  > ${dirI}/temp_SNP_diff_exclude; wc -l ${dirI}/temp_SNP_diff_exclude
       #EPI/QTWIN 391 (4 failed both test). AMFS 382 - all uniq. 610k 4540 ~ 500 overlap (note though that gender it a bit confoudned with status here). QMEGA_omni 2716; some overlap. HEI 4556 ~10 overlap WAMHS 7046 - most overlap
    echo ""
    plink_1.90 --threads 1 --bfile ${data_in}_geno3_mind3_hweCo_hweCa --exclude ${dirI}/temp_SNP_diff_exclude --make-bed --out ${data_in}_geno3_mind3_hweCo_hweCa_diff
      #EPI/QTWIN 285279 variants loaded --exclude: 284888 variants Total genotyping rate is 0.999834. 284888 variants and 1728 people pass filters and QC.
      #AMFS 797753 variants loaded --exclude: 797371 variants remaining Total genotyping rate is 0.999041. 797371 variants and 971 people pass filters and QC. 
      #610k 500255 variants loaded --exclude: 495715 variants remaining. Total genotyping rate is 0.99958. 495715 variants and 4851 people pass filters and QC.
      #QMEGA_omni 785330 variants loaded --exclude: 782614 variants remaining. Total genotyping rate is 0.999374. 782614 variants and 1252 people pass filters and QC.
      #HEI 212861 variants loaded --exclude: 208305 variants remaining. Total genotyping rate is 0.999284. 208305 variants and 2419 people pass filters and QC
      #WAMHS 590452 variants loaded --exclude: 583406 variants remaining. Total genotyping rate is 0.999112. 583406 variants and 3424 people pass filters and QC.
    echo ""
    #clean
    rm -f ${dirI}/temp_SNP_diff_exclude ${dirI}/temp_SNP_diff_missing_by_status ${dirI}/temp_SNP_diff_missing_by_gender ${data_in}_geno3_mind3_hweCo_hweCa.* ${data_in}_geno3_mind3_hweCo_hweCa_genders_case.fam ${data_in}_geno3_mind3_hweCo

  done #cleaning
done #turn on/off

########################################################################
#8.0 Combine with each other, and with 1KG then run PCA and IBD
########################################################################
#From what I have seen elsewhere the easiest thing to do is combine the GWAS sets, get the minimal list of clean SNPs, merge in 1KG, clean, and then LD prune etc.

#temp_2016_cleaning_AMFS_aligned_1KG temp_2016_cleaning_QMEGA_610k_aligned_1KG temp_2016_cleaning_QMEGA_omni_aligned_1KG temp_2016_cleaning_EPIGENE_QTWIN_aligned_1KG temp_2016_cleaning_HEIDELBERG_aligned_1KG

#no fancy way to do this off IDs so
  #awk '{print $1,$2,"OAG1_controls"}' temp_2016_cleaning_OAGphase1_controls_aligned_1KG.fam > OAG_controls_IDs.txt
  #awk '{print $1,$2,"OAG2_controls"}' temp_2016_cleaning_OAGphase2_controls_aligned_1KG.fam >> OAG_controls_IDs.txt

#takes a couple of hours, don't really need to do each  time
for x in #1
do
  echo -e "\nRunning IBD/PCA\n"
  echo "FID IID GROUP" > ${dirI}/temp_PCA_group

  set1=temp_2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff
  var1=`echo $set1 | sed 's/temp_2016_cleaning_//g' | sed 's/_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff//g'`
    awk '$6==1 {print $1,$2,"'''${var1}'''_control"}' ${dirI}/${set1}.fam >> ${dirI}/temp_PCA_group
    awk '$6==2 {print $1,$2,"'''${var1}'''_case"}' ${dirI}/${set1}.fam  >> ${dirI}/temp_PCA_group 
    cut -d" " -f3 ${dirI}/temp_PCA_group | sort | uniq -c
  echo ""
  set2=temp_2016_cleaning_QMEGA_610k_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff
  var2=`echo $set2 | sed 's/temp_2016_cleaning_//g' | sed 's/_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff//g'` 
    awk '$6==1 && $1~/E/ {print $1,$2,"ENDO_control"}' ${dirI}/${set2}.fam  >> ${dirI}/temp_PCA_group
    awk '$6==1 && !($1~/E/) {print $1,$2,"610k_QTWIN_control"}' ${dirI}/${set2}.fam  >> ${dirI}/temp_PCA_group 
    awk '$6==2 {print $1,$2,"'''${var2}'''_case"}' ${dirI}/${set2}.fam  >> ${dirI}/temp_PCA_group
    cut -d" " -f3 ${dirI}/temp_PCA_group | sort | uniq -c
      #   1794 610k_QTWIN_control
      # 544 AMFS_case
      #427 AMFS_control
      #2136 ENDO_control
      #1 group
      #921 QMEGA_610k_case
  echo ""
  plink_1.90 --threads 1 --bfile ${dirI}/${set1} --bmerge ${dirI}/${set2} --make-bed --out ${dirI}/temp_merge_1
     #AMFS and QMEGA_610k - rs7961667' and 'SNP12-132128064 have the same position, only pair, not worth fixing. Total genotyping rate is 0.546202 999157 variants and 5822 people pass filters and QC.
  echo ""
  set1=temp_merge_1
  set2=temp_2016_cleaning_QMEGA_omni_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff
   var2=`echo $set2 | sed 's/temp_2016_cleaning_//g' | sed 's/_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff//g'`
    awk '$6==1 {print $1,$2,"SDH_control"}' ${dirI}/${set2}.fam  >> ${dirI}/temp_PCA_group
    awk '$6==2 {print $1,$2,"'''${var2}'''_case"}' ${dirI}/${set2}.fam  >> ${dirI}/temp_PCA_group
    cut -d" " -f3 ${dirI}/temp_PCA_group | sort | uniq -c
  echo ""
  plink_1.90 --threads 1 --bfile ${dirI}/${set1} --bmerge ${dirI}/${set2} --make-bed --out ${dirI}/temp_merge_2
     #again warning rs7961667' and 'SNP12-132128064 . Total genotyping rate is 0.586701. 1001534 variants and 7074 people pass filters and QC.
  echo ""
  set1=temp_merge_2
  set2=temp_2016_cleaning_EPIGENE_QTWIN_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff
   var2=`echo $set2 | sed 's/temp_2016_cleaning_//g' | sed 's/_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff//g'`
    awk '$6==1 {print $1,$2,"'''${var2}'''_control"}' ${dirI}/${set2}.fam  >> ${dirI}/temp_PCA_group
    awk '$6==2 {print $1,$2,"'''${var2}'''_case"}' ${dirI}/${set2}.fam  >> ${dirI}/temp_PCA_group
    cut -d" " -f3 ${dirI}/temp_PCA_group | sort | uniq -c
  echo""
  plink_1.90 --threads 1 --bfile ${dirI}/${set1} --bmerge ${dirI}/${set2} --make-bed --out ${dirI}/temp_merge_1
     #now have 11 multiple position warnings which is annoying, but have to be SNPs not in 1kG. Not enough to fix given this is a quick compariosn. 11 3+ alleles.
  #First pass testing found 2 indels not in phase 1 1kg (so left as I/D) but have the correct alleles in the QTWIN/EPIGENE so can't be resolved. Checked the other nine are just flips not in 1kG so couldn't be aligned
  echo rs28369942 > ${dirI}/temp_merge_1-merge.exclude
  echo rs17884215 >> ${dirI}/temp_merge_1-merge.exclude
  plink_1.90 --threads 1 --bfile ${dirI}/${set2} --flip ${dirI}/temp_merge_1-merge.missnp --make-bed --out ${dirI}/${set2}_flip --exclude ${dirI}/temp_merge_1-merge.exclude
  echo ""
  plink_1.90 --threads 1 --bfile ${dirI}/${set1} --bmerge ${dirI}/${set2}_flip --make-bed --out ${dirI}/temp_merge_1
     #263 same position warning, 9 multiple position warnings. Not sure if may as well remap them properly to have it done? Total genotyping rate is 0.505295 1045237 variants and 8802 people pass filters and QC.
  #interimn clean
  rm -f ${dirI}/temp_merge_1-merge.exclude ${dirI}/temp_merge_1-merge.missnp ${dirI}/${set2}_flip
  echo ""
  set1=temp_merge_1
  set2=temp_2016_cleaning_HEIDELBERG_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff
   var2=`echo $set2 | sed 's/temp_2016_cleaning_//g' | sed 's/_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff//g'`
    awk '$6==1 {print $1,$2,"'''${var2}'''_control"}' ${dirI}/${set2}.fam  >> ${dirI}/temp_PCA_group
    awk '$6==2 {print $1,$2,"'''${var2}'''_case"}' ${dirI}/${set2}.fam  >> ${dirI}/temp_PCA_group
    cut -d" " -f3 ${dirI}/temp_PCA_group | sort | uniq -c
  echo ""
  plink_1.90 --threads 1 --bfile ${dirI}/${set1} --bmerge ${dirI}/${set2} --make-bed --out ${dirI}/temp_merge_2
    #221 3+ alleles
  echo ""
  plink_1.90 --threads 1 --bfile ${dirI}/${set2} --flip ${dirI}/temp_merge_2-merge.missnp --make-bed --out ${dirI}/${set2}_flip
     #--flip: 221 SNPs flipped.
  echo ""
  plink_1.90 --threads 1 --bfile ${dirI}/${set1} --bmerge ${dirI}/${set2}_flip --make-bed --out ${dirI}/temp_merge_2
     #322 same positions warnings, Total genotyping rate is 0.437729 1048978 variants and 11221 people pass filters and QC.
  rm -f ${dirI}/${set2}_flip ${dirI}/temp_merge_2-merge.missnp
  echo ""
  set1=temp_merge_2
  set2=temp_2016_cleaning_WAMHS_OAG_IBD_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff
  var2=`echo $set2 | sed 's/temp_2016_cleaning_//g' | sed 's/_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff//g'`
  cat ${dirI}/OAG_controls_IDs.txt  >> ${dirI}/temp_PCA_group
  awk 'NR==FNR {a[$1,$2];next} !($1,$2) in a {print $1,$2,"IBD_controls"}' ${dirI}/OAG_controls_IDs.txt <( awk '$6==1 {print $1,$2}' ${dirI}/${set2}.fam ) >> ${dirI}/temp_PCA_group
  awk '$6==2 {print $1,$2,"'''${var2}'''_case"}' ${dirI}/${set2}.fam  >> ${dirI}/temp_PCA_group
  cut -d" " -f3 ${dirI}/temp_PCA_group | sort | uniq -c
    # 1794 610k_QTWIN_control
    #  544 AMFS_case
    #  427 AMFS_control
    # 2136 ENDO_control
    #  779 EPIGENE_QTWIN_case
    #  949 EPIGENE_QTWIN_control
    #    1 GROUP
    # 1202 HEIDELBERG_case
    # 1217 HEIDELBERG_control
    #  904 IBD_controls
    #  651 OAG1_controls
    #  645 OAG2_controls
    #  921 QMEGA_610k_case
    #  689 QMEGA_omni_case
    #  562 SDH_control
    # 1268 WAMHS_OAG_IBD_case
  plink_1.90 --threads 1 --bfile ${dirI}/${set1} --bmerge ${dirI}/${set2} --make-bed --out ${dirI}/temp_merge_1
     #108 3+ alleles
  plink_1.90 --threads 1 --bfile ${dirI}/${set2} --flip ${dirI}/temp_merge_1-merge.missnp --make-bed --out ${dirI}/${set2}_flip
     #--flip: 108 SNPs flipped.
  echo ""
  plink_1.90 --threads 1 --bfile ${dirI}/${set1} --bmerge ${dirI}/${set2}_flip --make-bed --out ${dirI}/temp_merge_1
     #408 same position warnings
  rm -f ${dirI}/temp_merge_1-merge.missnp ${dirI}/${set2}_flip 
  echo
  wc -l ${dirI}/temp_merge_1.fam
    #14645
  echo ""
  wc -l ${dirI}/temp_PCA_group
    #14689 <<about right, grabeed the OAG IDs pre cleaning. 
  echo ""
  #get the 1kG group data  
  awk '{print $1,$2,$3,$4,$5,$6}' /mnt/lustre/home/matthewL/QUINN_PCA/DATA/ALL_data_merged_1000G_numeric_codes_unix.txt > ${dirI}/temp_ALL_data_merged_1000G_numeric_codes.txt
  \cp /mnt/lustre/home/matthewL/QUINN_PCA/DATA/final_1000G_codes_w_ceph_unix.txt ${dirI}/temp_final_1000G_codes_w_ceph.txt
  sed -i 's/\t/ /g' ${dirI}/temp_ALL_data_merged_1000G_numeric_codes.txt
  sed -i 's/ \+/ /g' ${dirI}/temp_ALL_data_merged_1000G_numeric_codes.txt

  #make a headed file to allow awk to merge headers
  echo "PHENOTYPE GROUP" > ${dirI}/temp_final_1000G_codes_w_ceph_2.txt
  sed 's/\t/ /g' ${dirI}/temp_final_1000G_codes_w_ceph.txt >> ${dirI}/temp_final_1000G_codes_w_ceph_2.txt
  echo ""
  head ${dirI}/temp_final_1000G_codes_w_ceph_2.txt
  echo ""
  head ${dirI}/temp_ALL_data_merged_1000G_numeric_codes.txt
  echo ""
  #append GROUP names to IDs
  awk 'NR==FNR {a[$1]=$2;next} $6 in a {print $1,$2,a[$6]}' ${dirI}/temp_final_1000G_codes_w_ceph_2.txt ${dirI}/temp_ALL_data_merged_1000G_numeric_codes.txt > ${dirI}/temp_ALL_data_merged_1000G_numeric_codes_groups.txt; wc -l ${dirI}/temp_ALL_data_merged_1000G_numeric_codes_groups.txt
    #1995
  echo ""
  head ${dirI}/temp_ALL_data_merged_1000G_numeric_codes_groups.txt
  echo ""
  echo "FID IID GROUP SOURCE" > ${dirI}/temp_ALL_data_merged_1000G_numeric_codes_groups_GWAS.txt  
  awk 'NR>1 {print $1,$2,$3,"GWAS"}' ${dirI}/temp_PCA_group >> ${dirI}/temp_ALL_data_merged_1000G_numeric_codes_groups_GWAS.txt
  awk 'NR>1 {print $1,$2,$3,"1000G"}' ${dirI}/temp_ALL_data_merged_1000G_numeric_codes_groups.txt >> ${dirI}/temp_ALL_data_merged_1000G_numeric_codes_groups_GWAS.txt; wc -l ${dirI}/temp_ALL_data_merged_1000G_numeric_codes_groups_GWAS.txt

  rm -f ${dirI}/temp_ALL_data_merged_1000G_numeric_codes.txt ${dirI}/temp_ALL_data_merged_1000G_numeric_codes_groups.txt ${dirI}/temp_PCA_group  ${dirI}/temp_final_1000G_codes_w_ceph.txt ${dirI}/temp_final_1000G_codes_w_ceph_2.txt

  G1000_pheno=${dirI}/temp_ALL_data_merged_1000G_numeric_codes_groups_GWAS.txt

  #extract 1kg on SNPLIST
  awk '{print $2}' ${dirI}/temp_merge_1.bim > ${dirI}/temp_merge_1_snplist
  for j in `seq 1 22`
  do
    dirV=/reference/genepi/public_reference_panels/1000G_20101123_v3/derived_plinkformat
    LD_FILE=/working/lab_stuartma/scratch/matthewL/melanoma_meta_analysis/MANIFEST_CHAIN_QC_FILES/eigenstrat_excluded_regions_hg19.txt
    echo -e "\nExtracting chr ${j} of the 1000 genomes data for temp_merge_1 SNPLIST"
    plink_1.90 --threads 1 --bfile ${dirV}/chr${j}.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL --extract ${dirI}/temp_merge_1_snplist --make-bed --out ${dirI}/temp_chr${j}.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL --silent
  done
  #make a merge list
  echo "${dirI}/temp_chr2.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.bed ${dirI}/temp_chr2.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.bim ${dirI}/temp_chr2.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.fam" > ${dirI}/temp_merge_1_HAPMAP_mergelist_chromsomes.txt

  for z in `seq 3 22`
  do
    echo "${dirI}/temp_chr${z}.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.bed ${dirI}/temp_chr${z}.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.bim ${dirI}/temp_chr${z}.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.fam" >> ${dirI}/temp_merge_1_HAPMAP_mergelist_chromsomes.txt
  done #merge loop

  awk '{print $1,$2,1}' ${dirI}/temp_chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.fam > ${dirI}/temp_chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.pheno

  LD_FILE=/working/lab_stuartma/scratch/matthewL/melanoma_meta_analysis/MANIFEST_CHAIN_QC_FILES/eigenstrat_excluded_regions_hg19.txt

  plink_1.90 --threads 1 --bfile ${dirI}/temp_chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL --merge-list ${dirI}/temp_merge_1_HAPMAP_mergelist_chromsomes.txt --exclude range ${LD_FILE} --make-bed --out ${dirI}/temp_merge_1_ALL.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL_pre --pheno ${dirI}/temp_chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.pheno
   #991299 variants and 1092 people pass filters and QC.

  #interimn clean
  rm -f ${dirI}/temp_chr*.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.* ${dirI}/temp_merge_1_HAPMAP_mergelist_chromsomes.txt ${dirI}/temp_chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.pheno

  echo ""
  #as 1kG aligned hopefully wont need to flip any SNPS
  plink_1.90 --threads 1 --bfile ${dirI}/temp_merge_1 --bmerge ${dirI}/temp_merge_1_ALL.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL_pre --make-bed --out ${dirI}/temp_merge_2
    #408 same position warnings. Total genotyping rate is 0.493379. 1060058 variants and 15737 people pass filters and QC. Among remaining phenotypes, 5403 are cases and 10333 are controls.
  echo ""
  #clean
  plink_1.90 --threads 1 --bfile ${dirI}/temp_merge_2 --geno 0.03 --make-bed --out  ${dirI}/temp_merge_1
     #1004866 variants removed due to missing genotype data (--geno) 55192 variants and 15737 people pass filters and QC. (!!!!)
  #+++++might need a different approach, losing too many SNPs
  plink_1.90 --threads 1 --bfile ${dirI}/temp_merge_1 --indep 50 5 2 --out ${dirI}/temp_merge_1_LD --exclude range ${LD_FILE}
    #Pruning complete.  12995 of 55192 variants removed.
  echo ""
  plink_1.90 --threads 1 --bfile ${dirI}/temp_merge_1 --extract ${dirI}/temp_merge_1_LD.prune.in --make-bed --out ${dirI}/temp_merge_pruned
    #42197 variants and 15737 people pass filters and QC.

  ###rough PCA ---takes about 30-40 min to run
  echo -e "\nRunning PCA with --silent otherwise floods log"
  plink_1.90 --threads 1 --bfile ${dirI}/temp_merge_pruned --pca 20 header --out ${dirI}/temp_merge_pruned_PCA --silent
  echo -e "PCA finished...\n"
  awk 'NR==FNR {a[$2]=$3" "$4;next} $2 in a {print $0,a[$2]}' ${G1000_pheno} ${dirI}/temp_merge_pruned_PCA.eigenvec > ${dirI}/temp_merge_pruned_PCA.eigenvec_IDs; wc -l ${dirI}/temp_merge_pruned_PCA.eigenvec_IDs
    #15737
  echo ""
  wc -l ${dirI}/temp_merge_pruned_PCA.eigenvec
    #15738 inc header
  echo ""
  head -4 ${dirI}/temp_merge_pruned_PCA.eigenvec_IDs
  echo -e "\nPlotting PCA results across sets; also backing up eigenvec data"
  rm -f ${dirK}/Melanoma_clean_2016_merged_PCA.eigenvec_IDs
  cp ${dirI}/temp_merge_pruned_PCA.eigenvec_IDs ${dirK}/Melanoma_clean_2016_merged_PCA.eigenvec_IDs
  

  data_set_R=temp_merge_pruned_PCA.eigenvec_IDs
  #PCA_out=MEGA_analysis
  PCA_out=Melanoma_clean_2016_merged
  module load R/3.2.2
R --no-save <<EOF

#16/03/2016 getting x11 errors with png, this might be a fix
options(bitmapType='cairo')
library(ggplot2)

read.table("${dirI}/${data_set_R}",header=T)->PCA_ALL;dim(PCA_ALL)
names(PCA_ALL)
table(PCA_ALL\$GROUP)

#################################
#CEU western european refence panel
#################################
PCA_ALL[PCA_ALL\$GROUP %in% "CEU",]->HAPMAP_CEU;dim(HAPMAP_CEU)
#85 23

table(HAPMAP_CEU\$GROUP)

mP1_CEU<-mean(HAPMAP_CEU\$PC1)
mP2_CEU<-mean(HAPMAP_CEU\$PC2)
sdP1_CEU<-sd(HAPMAP_CEU\$PC1)
sdP2_CEU<-sd(HAPMAP_CEU\$PC2)

#rect(mP1_CEU-6*sdP1_CEU,mP2_CEU-6*sdP2_CEU,mP1_CEU+6*sdP1_CEU,mP2_CEU+6*sdP2_CEU)

PCA_ALL[(PCA_ALL\$PC1>(mP1_CEU+6*sdP1_CEU) | PCA_ALL\$PC1<(mP1_CEU-6*sdP1_CEU) | PCA_ALL\$PC2>(mP2_CEU+6*sdP2_CEU) | PCA_ALL\$PC2<(mP2_CEU-6*sdP2_CEU) ),]->PCA_ALL_CEU_outliers6; dim(PCA_ALL_CEU_outliers6)

table(PCA_ALL_CEU_outliers6\$GROUP)

#################################
#GBR reference panel
#################################
PCA_ALL[PCA_ALL\$GROUP %in% "GBR",]->HAPMAP_GBR;dim(HAPMAP_GBR)

table(HAPMAP_GBR\$GROUP)

mP1_GBR<-mean(HAPMAP_GBR\$PC1)
mP2_GBR<-mean(HAPMAP_GBR\$PC2)
sdP1_GBR<-sd(HAPMAP_GBR\$PC1)
sdP2_GBR<-sd(HAPMAP_GBR\$PC2)

#rect(mP1_GBR-6*sdP1_GBR,mP2_GBR-6*sdP2_GBR,mP1_GBR+6*sdP1_GBR,mP2_GBR+6*sdP2_GBR)

PCA_ALL[(PCA_ALL\$PC1>(mP1_GBR+6*sdP1_GBR) | PCA_ALL\$PC1<(mP1_GBR-6*sdP1_GBR) | PCA_ALL\$PC2>(mP2_GBR+6*sdP2_GBR) | PCA_ALL\$PC2<(mP2_GBR-6*sdP2_GBR) ),]->PCA_ALL_GBR_outliers6; dim(PCA_ALL_GBR_outliers6)

table(PCA_ALL_GBR_outliers6\$GROUP)

#################################
#CEU+GBR+FIN northern/western european reference panel
#################################

PCA_ALL[PCA_ALL\$GROUP %in% c("GBR","CEU","FIN"),]->PCA_NW;dim(PCA_NW)

table(PCA_NW\$GROUP)

mP1_NW<-mean(PCA_NW\$PC1)
sdP1_NW<-sd(PCA_NW\$PC1)
mP2_NW<-mean(PCA_NW\$PC2)
sdP2_NW<-sd(PCA_NW\$PC2)

#24/04/2015 for historical reasons was called PCA_ALL_NW_outlier2; change N to reflect SD now

PCA_ALL[(PCA_ALL\$PC1>(mP1_NW+6*sdP1_NW) | PCA_ALL\$PC1<(mP1_NW-6*sdP1_NW) | PCA_ALL\$PC2>(mP2_NW+6*sdP2_NW) | PCA_ALL\$PC2<(mP2_NW-6*sdP2_NW) ),]->PCA_ALL_NW_outliers6; dim(PCA_ALL_NW_outliers6)

#24/04/2014 was calling this PCA_ALL_NW_outliers2, make it PCA_ALL_NW_outliers6 to reflect SD
table(PCA_ALL_NW_outliers6\$GROUP)

PCA_ALL[(PCA_ALL\$PC1>(mP1_NW+3*sdP1_NW) | PCA_ALL\$PC1<(mP1_NW-3*sdP1_NW) | PCA_ALL\$PC2>(mP2_NW+3*sdP2_NW) | PCA_ALL\$PC2<(mP2_NW-3*sdP2_NW) ),]->PCA_ALL_NW_outliers3; dim(PCA_ALL_NW_outliers3)

table(PCA_ALL_NW_outliers3\$GROUP)

#filters out some more MDACC and WAMHS which is not a problem - doesn't get all the CLM and MXL pops though, removes the PUR though.
#as always it is a balance - too strict etc

PCA_ALL[(PCA_ALL\$PC1>(mP1_NW+2*sdP1_NW) | PCA_ALL\$PC1<(mP1_NW-2*sdP1_NW) | PCA_ALL\$PC2>(mP2_NW+2*sdP2_NW) | PCA_ALL\$PC2<(mP2_NW-2*sdP2_NW) ),]->PCA_ALL_NW_outliers2; dim(PCA_ALL_NW_outliers2)

table(PCA_ALL_NW_outliers2\$GROUP)
#################################
#PLOTS
#################################
#24/04/2015 update rect naming to include SD info, and for NW include different SD for comparison

#need to call rect slightly differently to get geom_rect to work; assign the data to a data_frame
rect_NW6<-data.frame(xmin=mP1_NW-6*sdP1_NW,ymin=mP2_NW-6*sdP2_NW,xmax=mP1_NW+6*sdP1_NW,ymax=mP2_NW+6*sdP2_NW)
rect_NW3<-data.frame(xmin=mP1_NW-3*sdP1_NW,ymin=mP2_NW-3*sdP2_NW,xmax=mP1_NW+3*sdP1_NW,ymax=mP2_NW+3*sdP2_NW)
rect_NW2<-data.frame(xmin=mP1_NW-2*sdP1_NW,ymin=mP2_NW-2*sdP2_NW,xmax=mP1_NW+2*sdP1_NW,ymax=mP2_NW+2*sdP2_NW)
rect_GBR6<-data.frame(xmin=mP1_GBR-6*sdP1_GBR,ymin=mP2_GBR-6*sdP2_GBR,xmax=mP1_GBR+6*sdP1_GBR,ymax=mP2_GBR+6*sdP2_GBR)
rect_CEU6<-data.frame(xmin=mP1_CEU-6*sdP1_CEU,ymin=mP2_CEU-6*sdP2_CEU,xmax=mP1_CEU+6*sdP1_CEU,ymax=mP2_CEU+6*sdP2_CEU)

#for geom_rect to work I have to assign the plot to an object, then call them together

plot_NW6<-qplot(PC1,PC2,data=PCA_ALL,colour=GROUP,main="MELANOMA vs. 1KG, bounded on CEU+GBR+FIN, 6sd")
plot_NW3<-qplot(PC1,PC2,data=PCA_ALL,colour=GROUP,main="MELANOMA vs. 1KG, bounded on CEU+GBR+FIN, 3sd")
plot_NW2<-qplot(PC1,PC2,data=PCA_ALL,colour=GROUP,main="MELANOMA vs. 1KG, bounded on CEU+GBR+FIN, 2sd")
plot_GBR6<-qplot(PC1,PC2,data=PCA_ALL,colour=GROUP,main="MELANOMA vs. 1KG, bounded on GBR, 6sd")
plot_CEU6<-qplot(PC1,PC2,data=PCA_ALL,colour=GROUP,main="MELANOMA vs. 1KG, bounded on CEU, 6sd")


#NW outliers, 6SD
png("${dirK}/${PCA_out}_PCA_plot_NW_reference.png")
plot_NW6 + geom_rect(data=rect_NW6,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

#NW outliers, 3SD
png("${dirK}/${PCA_out}_PCA_plot_NW_reference_3sd.png")
plot_NW3 + geom_rect(data=rect_NW3,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

#NW outliers, 2SD
png("${dirK}/${PCA_out}_PCA_plot_NW_reference_2sd.png")
plot_NW2 + geom_rect(data=rect_NW2,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

#GBR outliers, 6SD
png("${dirK}/${PCA_out}_PCA_plot_GBR_reference.png")
plot_GBR6 + geom_rect(data=rect_GBR6,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

#CEU outliers, 6SD
png("${dirK}/${PCA_out}_PCA_plot_CEU_reference.png")
plot_CEU6 + geom_rect(data=rect_CEU6,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

ptions(bitmapType='cairo')
library(ggplot2)

#write out study outlier just for record keeping, may as well exclude 1000G samples
write.table(PCA_ALL_NW_outliers6[PCA_ALL_NW_outliers6\$SOURCE %in% "GWAS",c(1,2)],"${dirK}/${PCA_out}_GBR_CEU_FIN_outlier_6sd.txt",row.names=F,col.names=T,quote=F)
write.table(PCA_ALL_NW_outliers3[PCA_ALL_NW_outliers3\$SOURCE %in% "GWAS",c(1,2)],"${dirK}/${PCA_out}_GBR_CEU_FIN_outlier_3sd.txt",row.names=F,col.names=T,quote=F)
write.table(PCA_ALL_NW_outliers2[PCA_ALL_NW_outliers2\$SOURCE %in% "GWAS",c(1,2)],"${dirK}/${PCA_out}_GBR_CEU_FIN_outlier_2sd.txt",row.names=F,col.names=T,quote=F)

#################################
#14.7 PLOTS - EUR only
#################################
#15/03/2016 would be useful to focus in on the EUR plots as well - I developed this for the MIA analysis but makes sense to use here as the full world plot is a bit hard to work out

#Since i am using NW the most look at some other PCs
mP3_NW<-mean(PCA_NW\$PC3)
sdP3_NW<-sd(PCA_NW\$PC3)
mP4_NW<-mean(PCA_NW\$PC4)
sdP4_NW<-sd(PCA_NW\$PC4)

#15/03/2016 change this to exclude non EUR ( in the source it included EUR + MIA but the inclusion will change as I add more sets
subset(PCA_ALL,!(PCA_ALL\$GROUP %in% c("ASW","CHB","CHS","CLM","JPT","LWK","MXL","PUR","YRI")))->PCA_EUR; dim(PCA_EUR)
  #[1] 6922   24

table(PCA_EUR\$GROUP)
  #     AMFS       ASW      CEPH       CEU       CHB       CHS       CLM   EPIGENE 
  #      547         0         2        83         0         0         0       717 
  #      FIN       GBR       IBS       JPT       LWK     MDACC       MIA       MXL 
  #       93        88        14         0         0      1806      1191         0 
  #      PUR     QMEGA QMEGA_610       TSI     WAMHS       YRI 
  #        0       481       551        98      1251         0

plot_NW6_EUR<-qplot(PC1,PC2,data=PCA_EUR,colour=GROUP,main="MELANOMA vs. 1KG EUR, bounded on CEU+GBR+FIN, 6sd")
plot_NW3_EUR<-qplot(PC1,PC2,data=PCA_EUR,colour=GROUP,main="MELANOMA vs. 1KG EUR, bounded on CEU+GBR+FIN, 3sd")
plot_NW2_EUR<-qplot(PC1,PC2,data=PCA_EUR,colour=GROUP,main="MELANOMA vs. 1KG EUR, bounded on CEU+GBR+FIN, 2sd")
plot_GBR6_EUR<-qplot(PC1,PC2,data=PCA_EUR,colour=GROUP,main="MELANOMA vs. 1KG EUR, bounded on GBR, 6sd")
plot_CEU6_EUR<-qplot(PC1,PC2,data=PCA_EUR,colour=GROUP,main="MELANOMA vs. 1KG EUR, bounded on CEU, 6sd")

#NW outliers, 6SD
png("${dirK}/${PCA_out}_PCA_plot_NW_reference_EUR.png")
plot_NW6_EUR + geom_rect(data=rect_NW6,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

#NW outliers, 3SD
png("${dirK}/${PCA_out}_PCA_plot_NW_reference_3sd_EUR.png")
plot_NW3_EUR + geom_rect(data=rect_NW3,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

#NW outliers, 2SD
png("${dirK}/${PCA_out}_PCA_plot_NW_reference_2sd_EUR.png")
plot_NW2_EUR + geom_rect(data=rect_NW2,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

#GBR outliers, 6SD
png("${dirK}/${PCA_out}_PCA_plot_GBR_reference_EUR.png")
plot_GBR6_EUR + geom_rect(data=rect_GBR6,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

#CEU outliers, 6SD
png("${dirK}/${PCA_out}_PCA_plot_CEU_reference_EUR.png")
plot_CEU6_EUR + geom_rect(data=rect_CEU6,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()


#####################################################################
#PLOTS per GWAS
#####################################################################
 #Big plots not very good for seeing within study differences so spin off per GWAS plots

 subset(PCA_ALL,!(PCA_ALL\$GROUP %in% c("AMFS_case", "AMFS_control", "EPIGENE_QTWIN_case", "EPIGENE_QTWIN_control", "HEIDELBERG_case", "HEIDELBERG_control", "WAMHS_OAG_IBD_case", "IBD_controls", "OAG1_controls", "OAG2_controls", "QMEGA_omni_case","SDH_control")))->PCA_610k

plot_NW6_610k<-qplot(PC1,PC2,data=PCA_610k,colour=GROUP,main="QMEGA_610k vs. 1KG, bounded on CEU+GBR+FIN, 6sd")
plot_NW3_610k<-qplot(PC1,PC2,data=PCA_610k,colour=GROUP,main="QMEGA_610k vs. 1KG, bounded on CEU+GBR+FIN, 3sd")
plot_NW2_610k<-qplot(PC1,PC2,data=PCA_610k,colour=GROUP,main="QMEGA_610k vs. 1KG, bounded on CEU+GBR+FIN, 2sd")
plot_GBR6_610k<-qplot(PC1,PC2,data=PCA_610k,colour=GROUP,main="QMEGA_610k vs. 1KG, bounded on GBR, 6sd")
plot_CEU6_610k<-qplot(PC1,PC2,data=PCA_610k,colour=GROUP,main="QMEGA_610k vs. 1KG, bounded on CEU, 6sd")

#at the moment only plot 6 and 3 NW
#NW outliers, 6SD
png("${dirK}/${PCA_out}_PCA_plot_NW_reference_610k.png")
plot_NW6_610k + geom_rect(data=rect_NW6,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

#NW outliers, 3SD
png("${dirK}/${PCA_out}_PCA_plot_NW_reference_3sd_610k.png")
plot_NW6_610k + geom_rect(data=rect_NW3,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

 subset(PCA_ALL,!(PCA_ALL\$GROUP %in% c("QMEGA_610k_case","610k_QTWIN_control", "ENDO_control", "EPIGENE_QTWIN_case", "EPIGENE_QTWIN_control", "HEIDELBERG_case", "HEIDELBERG_control", "WAMHS_OAG_IBD_case", "IBD_controls", "OAG1_controls", "OAG2_controls", "QMEGA_omni_case","SDH_control")))->PCA_AMFS

plot_NW6_AMFS<-qplot(PC1,PC2,data=PCA_AMFS,colour=GROUP,main="AMFS vs. 1KG, bounded on CEU+GBR+FIN, 6sd")
plot_NW3_AMFS<-qplot(PC1,PC2,data=PCA_AMFS,colour=GROUP,main="AMFS vs. 1KG, bounded on CEU+GBR+FIN, 3sd")
plot_NW2_AMFS<-qplot(PC1,PC2,data=PCA_AMFS,colour=GROUP,main="AMFS vs. 1KG, bounded on CEU+GBR+FIN, 2sd")
plot_GBR6_AMFS<-qplot(PC1,PC2,data=PCA_AMFS,colour=GROUP,main="AMFS vs. 1KG, bounded on GBR, 6sd")
plot_CEU6_AMFS<-qplot(PC1,PC2,data=PCA_AMFS,colour=GROUP,main="AMFS vs. 1KG, bounded on CEU, 6sd")

#at the moment only plot 6 and 3 NW
#NW outliers, 6SD
png("${dirK}/${PCA_out}_PCA_plot_NW_reference_AMFS.png")
plot_NW6_AMFS + geom_rect(data=rect_NW6,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

#NW outliers, 3SD
png("${dirK}/${PCA_out}_PCA_plot_NW_reference_3sd_AMFS.png")
plot_NW6_AMFS + geom_rect(data=rect_NW3,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

 subset(PCA_ALL,!(PCA_ALL\$GROUP %in% c("QMEGA_610k_case","610k_QTWIN_control", "ENDO_control", "AMFS_case", "AMFS_control", "HEIDELBERG_case", "HEIDELBERG_control", "WAMHS_OAG_IBD_case", "IBD_controls", "OAG1_controls", "OAG2_controls", "QMEGA_omni_case","SDH_control")))->PCA_EPIGENE_QTWIN

plot_NW6_EPIGENE_QTWIN<-qplot(PC1,PC2,data=PCA_EPIGENE_QTWIN,colour=GROUP,main="EPIGENE_QTWIN vs. 1KG, bounded on CEU+GBR+FIN, 6sd")
plot_NW3_EPIGENE_QTWIN<-qplot(PC1,PC2,data=PCA_EPIGENE_QTWIN,colour=GROUP,main="EPIGENE_QTWIN vs. 1KG, bounded on CEU+GBR+FIN, 3sd")
plot_NW2_EPIGENE_QTWIN<-qplot(PC1,PC2,data=PCA_EPIGENE_QTWIN,colour=GROUP,main="EPIGENE_QTWIN vs. 1KG, bounded on CEU+GBR+FIN, 2sd")
plot_GBR6_EPIGENE_QTWIN<-qplot(PC1,PC2,data=PCA_EPIGENE_QTWIN,colour=GROUP,main="EPIGENE_QTWIN vs. 1KG, bounded on GBR, 6sd")
plot_CEU6_EPIGENE_QTWIN<-qplot(PC1,PC2,data=PCA_EPIGENE_QTWIN,colour=GROUP,main="EPIGENE_QTWIN vs. 1KG, bounded on CEU, 6sd")

#at the moment only plot 6 and 3 NW
#NW outliers, 6SD
png("${dirK}/${PCA_out}_PCA_plot_NW_reference_EPIGENE_QTWIN.png")
plot_NW6_EPIGENE_QTWIN + geom_rect(data=rect_NW6,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

#NW outliers, 3SD
png("${dirK}/${PCA_out}_PCA_plot_NW_reference_3sd_EPIGENE_QTWIN.png")
plot_NW6_EPIGENE_QTWIN + geom_rect(data=rect_NW3,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

 subset(PCA_ALL,!(PCA_ALL\$GROUP %in% c("QMEGA_610k_case","610k_QTWIN_control", "ENDO_control", "AMFS_case", "AMFS_control", "EPIGENE_QTWIN_case", "EPIGENE_QTWIN_control", "WAMHS_OAG_IBD_case", "IBD_controls", "OAG1_controls", "OAG2_controls", "QMEGA_omni_case","SDH_control")))->PCA_HEI

plot_NW6_HEI<-qplot(PC1,PC2,data=PCA_HEI,colour=GROUP,main="HEIDELBERG vs. 1KG, bounded on CEU+GBR+FIN, 6sd")
plot_NW3_HEI<-qplot(PC1,PC2,data=PCA_HEI,colour=GROUP,main="HEIDELBERG vs. 1KG, bounded on CEU+GBR+FIN, 3sd")
plot_NW2_HEI<-qplot(PC1,PC2,data=PCA_HEI,colour=GROUP,main="HEIDELBERG vs. 1KG, bounded on CEU+GBR+FIN, 2sd")
plot_GBR6_HEI<-qplot(PC1,PC2,data=PCA_HEI,colour=GROUP,main="HEIDELBERG vs. 1KG, bounded on GBR, 6sd")
plot_CEU6_HEI<-qplot(PC1,PC2,data=PCA_HEI,colour=GROUP,main="HEIDELBERG vs. 1KG, bounded on CEU, 6sd")

#at the moment only plot 6 and 3 NW
#NW outliers, 6SD
png("${dirK}/${PCA_out}_PCA_plot_NW_reference_HEI.png")
plot_NW6_HEI + geom_rect(data=rect_NW6,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

#NW outliers, 3SD
png("${dirK}/${PCA_out}_PCA_plot_NW_reference_3sd_HEI.png")
plot_NW6_HEI + geom_rect(data=rect_NW3,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

 subset(PCA_ALL,!(PCA_ALL\$GROUP %in% c("QMEGA_610k_case","610k_QTWIN_control", "ENDO_control", "AMFS_case", "AMFS_control", "EPIGENE_QTWIN_case", "EPIGENE_QTWIN_control", "HEIDELBERG_case", "HEIDELBERG_control", "QMEGA_omni_case","SDH_control")))->PCA_WAMHS

plot_NW6_WAMHS<-qplot(PC1,PC2,data=PCA_WAMHS,colour=GROUP,main="WAMHS vs. 1KG, bounded on CEU+GBR+FIN, 6sd")
plot_NW3_WAMHS<-qplot(PC1,PC2,data=PCA_WAMHS,colour=GROUP,main="WAMHS vs. 1KG, bounded on CEU+GBR+FIN, 3sd")
plot_NW2_WAMHS<-qplot(PC1,PC2,data=PCA_WAMHS,colour=GROUP,main="WAMHS vs. 1KG, bounded on CEU+GBR+FIN, 2sd")
plot_GBR6_WAMHS<-qplot(PC1,PC2,data=PCA_WAMHS,colour=GROUP,main="WAMHS vs. 1KG, bounded on GBR, 6sd")
plot_CEU6_WAMHS<-qplot(PC1,PC2,data=PCA_WAMHS,colour=GROUP,main="WAMHS vs. 1KG, bounded on CEU, 6sd")

#at the moment only plot 6 and 3 NW
#NW outliers, 6SD
png("${dirK}/${PCA_out}_PCA_plot_NW_reference_WAMHS.png")
plot_NW6_WAMHS + geom_rect(data=rect_NW6,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

#NW outliers, 3SD
png("${dirK}/${PCA_out}_PCA_plot_NW_reference_3sd_WAMHS.png")
plot_NW6_WAMHS + geom_rect(data=rect_NW3,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

 subset(PCA_ALL,!(PCA_ALL\$GROUP %in% c("QMEGA_610k_case","610k_QTWIN_control", "ENDO_control", "AMFS_case", "AMFS_control", "EPIGENE_QTWIN_case", "EPIGENE_QTWIN_control", "HEIDELBERG_case", "HEIDELBERG_control", "WAMHS_OAG_IBD_case", "IBD_controls", "OAG1_controls", "OAG2_controls")))->PCA_QMEGA_omni

plot_NW6_QMEGA_omni<-qplot(PC1,PC2,data=PCA_QMEGA_omni,colour=GROUP,main="QMEGA_omni vs. 1KG, bounded on CEU+GBR+FIN, 6sd")
plot_NW3_QMEGA_omni<-qplot(PC1,PC2,data=PCA_QMEGA_omni,colour=GROUP,main="QMEGA_omni vs. 1KG, bounded on CEU+GBR+FIN, 3sd")
plot_NW2_QMEGA_omni<-qplot(PC1,PC2,data=PCA_QMEGA_omni,colour=GROUP,main="QMEGA_omni vs. 1KG, bounded on CEU+GBR+FIN, 2sd")
plot_GBR6_QMEGA_omni<-qplot(PC1,PC2,data=PCA_QMEGA_omni,colour=GROUP,main="QMEGA_omni vs. 1KG, bounded on GBR, 6sd")
plot_CEU6_QMEGA_omni<-qplot(PC1,PC2,data=PCA_QMEGA_omni,colour=GROUP,main="QMEGA_omni vs. 1KG, bounded on CEU, 6sd")

#at the moment only plot 6 and 3 NW
#NW outliers, 6SD
png("${dirK}/${PCA_out}_PCA_plot_NW_reference_QMEGA_omni.png")
plot_NW6_QMEGA_omni + geom_rect(data=rect_NW6,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

#NW outliers, 3SD
png("${dirK}/${PCA_out}_PCA_plot_NW_reference_3sd_QMEGA_omni.png")
plot_NW6_QMEGA_omni + geom_rect(data=rect_NW3,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()


#####################################################################
#PLOTS per GWAS EUR only
#####################################################################
 #Big plots not very good for seeing within study differences to spin of per GWAS plots - agin for EUR set

 subset(PCA_EUR,!(PCA_EUR\$GROUP %in% c("AMFS_case", "AMFS_control", "EPIGENE_QTWIN_case", "EPIGENE_QTWIN_control", "HEIDELBERG_case", "HEIDELBERG_control", "WAMHS_OAG_IBD_case", "IBD_controls", "OAG1_controls", "OAG2_controls", "QMEGA_omni_case","SDH_control")))->PCA_EUR_610k

plot_NW6_610k_EUR<-qplot(PC1,PC2,data=PCA_EUR_610k,colour=GROUP,main="QMEGA_610k vs. 1KG EUR, bounded on CEU+GBR+FIN, 6sd")
plot_NW3_610k_EUR<-qplot(PC1,PC2,data=PCA_EUR_610k,colour=GROUP,main="QMEGA_610k vs. 1KG EUR, bounded on CEU+GBR+FIN, 3sd")
plot_NW2_610k_EUR<-qplot(PC1,PC2,data=PCA_EUR_610k,colour=GROUP,main="QMEGA_610k vs. 1KG EUR, bounded on CEU+GBR+FIN, 2sd")
plot_GBR6_610k_EUR<-qplot(PC1,PC2,data=PCA_EUR_610k,colour=GROUP,main="QMEGA_610k vs. 1KG EUR, bounded on GBR, 6sd")
plot_CEU6_610k_EUR<-qplot(PC1,PC2,data=PCA_EUR_610k,colour=GROUP,main="QMEGA_610k vs. 1KG EUR, bounded on CEU, 6sd")

#at the moment only plot 6 and 3 NW
#NW outliers, 6SD
png("${dirK}/${PCA_out}_PCA_plot_NW_reference_610k_EUR.png")
plot_NW6_610k_EUR + geom_rect(data=rect_NW6,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

#NW outliers, 3SD
png("${dirK}/${PCA_out}_PCA_plot_NW_reference_3sd_610k_EUR.png")
plot_NW6_610k_EUR + geom_rect(data=rect_NW3,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

 subset(PCA_EUR,!(PCA_EUR\$GROUP %in% c("QMEGA_610k_case","610k_QTWIN_control", "ENDO_control", "EPIGENE_QTWIN_case", "EPIGENE_QTWIN_control", "HEIDELBERG_case", "HEIDELBERG_control", "WAMHS_OAG_IBD_case", "IBD_controls", "OAG1_controls", "OAG2_controls", "QMEGA_omni_case","SDH_control")))->PCA_EUR_AMFS

plot_NW6_AMFS_EUR<-qplot(PC1,PC2,data=PCA_EUR_AMFS,colour=GROUP,main="AMFS vs. 1KG EUR, bounded on CEU+GBR+FIN, 6sd")
plot_NW3_AMFS_EUR<-qplot(PC1,PC2,data=PCA_EUR_AMFS,colour=GROUP,main="AMFS vs. 1KG EUR, bounded on CEU+GBR+FIN, 3sd")
plot_NW2_AMFS_EUR<-qplot(PC1,PC2,data=PCA_EUR_AMFS,colour=GROUP,main="AMFS vs. 1KG EUR, bounded on CEU+GBR+FIN, 2sd")
plot_GBR6_AMFS_EUR<-qplot(PC1,PC2,data=PCA_EUR_AMFS,colour=GROUP,main="AMFS vs. 1KG EUR, bounded on GBR, 6sd")
plot_CEU6_AMFS_EUR<-qplot(PC1,PC2,data=PCA_EUR_AMFS,colour=GROUP,main="AMFS vs. 1KG EUR, bounded on CEU, 6sd")

#at the moment only plot 6 and 3 NW
#NW outliers, 6SD
png("${dirK}/${PCA_out}_PCA_plot_NW_reference_AMFS_EUR.png")
plot_NW6_AMFS_EUR + geom_rect(data=rect_NW6,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

#NW outliers, 3SD
png("${dirK}/${PCA_out}_PCA_plot_NW_reference_3sd_AMFS_EUR.png")
plot_NW6_AMFS_EUR + geom_rect(data=rect_NW3,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

 subset(PCA_EUR,!(PCA_EUR\$GROUP %in% c("QMEGA_610k_case","610k_QTWIN_control", "ENDO_control", "AMFS_case", "AMFS_control", "HEIDELBERG_case", "HEIDELBERG_control", "WAMHS_OAG_IBD_case", "IBD_controls", "OAG1_controls", "OAG2_controls", "QMEGA_omni_case","SDH_control")))->PCA_EPIGENE_QTWIN_EUR

plot_NW6_EPIGENE_QTWIN_EUR<-qplot(PC1,PC2,data=PCA_EPIGENE_QTWIN_EUR,colour=GROUP,main="EPIGENE_QTWIN vs. 1KG EUR, bounded on CEU+GBR+FIN, 6sd")
plot_NW3_EPIGENE_QTWIN_EUR<-qplot(PC1,PC2,data=PCA_EPIGENE_QTWIN_EUR,colour=GROUP,main="EPIGENE_QTWIN vs. 1KG EUR, bounded on CEU+GBR+FIN, 3sd")
plot_NW2_EPIGENE_QTWIN_EUR<-qplot(PC1,PC2,data=PCA_EPIGENE_QTWIN_EUR,colour=GROUP,main="EPIGENE_QTWIN vs. 1KG EUR, bounded on CEU+GBR+FIN, 2sd")
plot_GBR6_EPIGENE_QTWIN_EUR<-qplot(PC1,PC2,data=PCA_EPIGENE_QTWIN_EUR,colour=GROUP,main="EPIGENE_QTWIN vs. 1KG EUR, bounded on GBR, 6sd")
plot_CEU6_EPIGENE_QTWIN_EUR<-qplot(PC1,PC2,data=PCA_EPIGENE_QTWIN_EUR,colour=GROUP,main="EPIGENE_QTWIN vs. 1KG EUR, bounded on CEU, 6sd")
m=c(0,0.5)
#at the moment only plot 6 and 3 NW
#NW outliers, 6SD
png("${dirK}/${PCA_out}_PCA_plot_NW_reference_EPIGENE_QTWIN_EUR.png")
plot_NW6_EPIGENE_QTWIN_EUR + geom_rect(data=rect_NW6,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

#NW outliers, 3SD
png("${dirK}/${PCA_out}_PCA_plot_NW_reference_3sd_EPIGENE_QTWIN_EUR.png")
plot_NW6_EPIGENE_QTWIN_EUR + geom_rect(data=rect_NW3,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

 subset(PCA_EUR,!(PCA_EUR\$GROUP %in% c("QMEGA_610k_case","610k_QTWIN_control", "ENDO_control", "AMFS_case", "AMFS_control", "EPIGENE_QTWIN_case", "EPIGENE_QTWIN_control", "WAMHS_OAG_IBD_case", "IBD_controls", "OAG1_controls", "OAG2_controls", "QMEGA_omni_case","SDH_control")))->PCA_HEI_EUR

plot_NW6_HEI_EUR<-qplot(PC1,PC2,data=PCA_HEI_EUR,colour=GROUP,main="HEIDELBERG vs. 1KG EUR, bounded on CEU+GBR+FIN, 6sd")
plot_NW3_HEI_EUR<-qplot(PC1,PC2,data=PCA_HEI_EUR,colour=GROUP,main="HEIDELBERG vs. 1KG EUR, bounded on CEU+GBR+FIN, 3sd")
plot_NW2_HEI_EUR<-qplot(PC1,PC2,data=PCA_HEI_EUR,colour=GROUP,main="HEIDELBERG vs. 1KG EUR, bounded on CEU+GBR+FIN, 2sd")
plot_GBR6_HEI_EUR<-qplot(PC1,PC2,data=PCA_HEI_EUR,colour=GROUP,main="HEIDELBERG vs. 1KG EUR, bounded on GBR, 6sd")
plot_CEU6_HEI_EUR<-qplot(PC1,PC2,data=PCA_HEI_EUR,colour=GROUP,main="HEIDELBERG vs. 1KG EUR, bounded on CEU, 6sd")

#at the moment only plot 6 and 3 NW
#NW outliers, 6SD
png("${dirK}/${PCA_out}_PCA_plot_NW_reference_HEI_EUR.png")
plot_NW6_HEI_EUR + geom_rect(data=rect_NW6,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

#NW outliers, 3SD
png("${dirK}/${PCA_out}_PCA_plot_NW_reference_3sd_HEI_EUR.png")
plot_NW6_HEI_EUR + geom_rect(data=rect_NW3,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

 subset(PCA_EUR,!(PCA_EUR\$GROUP %in% c("QMEGA_610k_case","610k_QTWIN_control", "ENDO_control", "AMFS_case", "AMFS_control", "EPIGENE_QTWIN_case", "EPIGENE_QTWIN_control", "HEIDELBERG_case", "HEIDELBERG_control", "QMEGA_omni_case","SDH_control")))->PCA_WAMHS_EUR

plot_NW6_WAMHS_EUR<-qplot(PC1,PC2,data=PCA_WAMHS_EUR,colour=GROUP,main="WAMHS vs. 1KG EUR, bounded on CEU+GBR+FIN, 6sd")
plot_NW3_WAMHS_EUR<-qplot(PC1,PC2,data=PCA_WAMHS_EUR,colour=GROUP,main="WAMHS vs. 1KG EUR, bounded on CEU+GBR+FIN, 3sd")
plot_NW2_WAMHS_EUR<-qplot(PC1,PC2,data=PCA_WAMHS_EUR,colour=GROUP,main="WAMHS vs. 1KG EUR, bounded on CEU+GBR+FIN, 2sd")
plot_GBR6_WAMHS_EUR<-qplot(PC1,PC2,data=PCA_WAMHS_EUR,colour=GROUP,main="WAMHS vs. 1KG EUR, bounded on GBR, 6sd")
plot_CEU6_WAMHS_EUR<-qplot(PC1,PC2,data=PCA_WAMHS_EUR,colour=GROUP,main="WAMHS vs. 1KG EUR, bounded on CEU, 6sd")

#at the moment only plot 6 and 3 NW
#NW outliers, 6SD
png("${dirK}/${PCA_out}_PCA_plot_NW_reference_WAMHS_EUR.png")
plot_NW6_WAMHS_EUR + geom_rect(data=rect_NW6,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

#NW outliers, 3SD
png("${dirK}/${PCA_out}_PCA_plot_NW_reference_3sd_WAMHS_EUR.png")
plot_NW6_WAMHS_EUR + geom_rect(data=rect_NW3,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

 subset(PCA_EUR,!(PCA_EUR\$GROUP %in% c("QMEGA_610k_case","610k_QTWIN_control", "ENDO_control", "AMFS_case", "AMFS_control", "EPIGENE_QTWIN_case", "EPIGENE_QTWIN_control", "HEIDELBERG_case", "HEIDELBERG_control", "WAMHS_OAG_IBD_case", "IBD_controls", "OAG1_controls", "OAG2_controls")))->PCA_QMEGA_omni_EUR

plot_NW6_QMEGA_omni_EUR<-qplot(PC1,PC2,data=PCA_QMEGA_omni_EUR,colour=GROUP,main="QMEGA_omni vs. 1KG EUR, bounded on CEU+GBR+FIN, 6sd")
plot_NW3_QMEGA_omni_EUR<-qplot(PC1,PC2,data=PCA_QMEGA_omni_EUR,colour=GROUP,main="QMEGA_omni vs. 1KG EUR, bounded on CEU+GBR+FIN, 3sd")
plot_NW2_QMEGA_omni_EUR<-qplot(PC1,PC2,data=PCA_QMEGA_omni_EUR,colour=GROUP,main="QMEGA_omni vs. 1KG EUR, bounded on CEU+GBR+FIN, 2sd")
plot_GBR6_QMEGA_omni_EUR<-qplot(PC1,PC2,data=PCA_QMEGA_omni_EUR,colour=GROUP,main="QMEGA_omni vs. 1KG EUR, bounded on GBR, 6sd")
plot_CEU6_QMEGA_omni_EUR<-qplot(PC1,PC2,data=PCA_QMEGA_omni_EUR,colour=GROUP,main="QMEGA_omni vs. 1KG EUR, bounded on CEU, 6sd")

#at the moment only plot 6 and 3 NW
#NW outliers, 6SD
png("${dirK}/${PCA_out}_PCA_plot_NW_reference_QMEGA_omni_EUR.png")
plot_NW6_QMEGA_omni_EUR + geom_rect(data=rect_NW6,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

#NW outliers, 3SD
png("${dirK}/${PCA_out}_PCA_plot_NW_reference_3sd_QMEGA_omni_EUR.png")
plot_NW6_QMEGA_omni_EUR + geom_rect(data=rect_NW3,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()


#################################
#14.8 PLOTS - NW, some other PCs
#################################
#since I am focusing on NW 3SD look at some other PCs
rect_NW3_2<-data.frame(xmin=mP1_NW-3*sdP1_NW,ymin=mP3_NW-3*sdP3_NW,xmax=mP1_NW+3*sdP1_NW,ymax=mP3_NW+3*sdP3_NW)
rect_NW3_3<-data.frame(xmin=mP2_NW-3*sdP2_NW,ymin=mP3_NW-3*sdP3_NW,xmax=mP2_NW+3*sdP2_NW,ymax=mP3_NW+3*sdP3_NW)
rect_NW3_4<-data.frame(xmin=mP1_NW-3*sdP1_NW,ymin=mP4_NW-3*sdP4_NW,xmax=mP1_NW+3*sdP1_NW,ymax=mP4_NW+3*sdP4_NW)

plot_NW3_2<-qplot(PC1,PC3,data=PCA_ALL,colour=GROUP,main="MELANOMA vs. 1KG EUR, bounded on CEU+GBR+FIN, 3sd")
plot_NW3_3<-qplot(PC2,PC3,data=PCA_ALL,colour=GROUP,main="MELANOMA vs. 1KG EUR, bounded on CEU+GBR+FIN, 3sd")
plot_NW3_4<-qplot(PC1,PC4,data=PCA_ALL,colour=GROUP,main="MELANOMA vs. 1KG EUR, bounded on CEU+GBR+FIN, 3sd")

#all samples, PC1 and PC3, NW 3sd
 #does not really help with the new samples - seperates out more CLM, MXL samples
png("${dirK}/${PCA_out}_PCA_plot_NW_reference_3sd_EUR_PC1_PC3.png")
plot_NW3_2 + geom_rect(data=rect_NW3_2,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

#all samples, PC2 and PC3, NW 3sd; again doesn't really help much
png("${dirK}/${PCA_out}_PCA_plot_NW_reference_3sd_EUR_PC2_PC3.png")
plot_NW3_3 + geom_rect(data=rect_NW3_3,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()

#all samples, PC1 and PC4, NW 3sd. Huh that is interesting. PC4 looks almost like Nth/Sth EUR division - Fin ++, TSI --, and again the slighly outlier clump for MIA more -ve than TSI. I wonder again if this is GRK admixture or something? TSI is toscani italian is central italian )well more northern, so it might suggest we have more southern italian. Well without greek or sicilian I can't resolve it, but hardly concerning.
png("${dirK}/${PCA_out}_PCA_plot_NW_reference_3sd_EUR_PC1_PC4.png")
plot_NW3_4 + geom_rect(data=rect_NW3_4,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), color="grey20", alpha=0.5, inherit.aes = FALSE)
dev.off()


q()
EOF
 
  #11/08/2016 Summary - PCA looks fine. Detail results below
  #CEU looks about right. Handful of AMFS, 610k, EPI etc. Is still pretty strict (as I know it is). More EPIGENE than the others but it never quite sat over the others. 
  #table(PCA_ALL_CEU_outliers6$GROUP)
    #610k_QTWIN_control             AMFS_case          AMFS_control 
    #                18                     6                     6 
    #               ASW                  CEPH                   CEU 
    #                61                     0                     0 
    #               CHB                   CHS                   CLM 
    #                97                   100                    59 
    #      ENDO_control    EPIGENE_QTWIN_case EPIGENE_QTWIN_control 
    #                24                    10                    51 
    #               FIN                   GBR       HEIDELBERG_case 
    #                57                     0                    25 
    #HEIDELBERG_control          IBD_controls                   IBS 
    #                 0                    37                     0 
    #               JPT                   LWK                   MXL 
    #                89                    97                    65 
    #     OAG1_controls         OAG2_controls                   PUR 
    #                19                    24                    55 
    #   QMEGA_610k_case       QMEGA_omni_case           SDH_control 
    #                 2                     4                    13 
    #               TSI    WAMHS_OAG_IBD_case                   YRI 
    #                 0                    17                    88 

    #GBR even more strict so as expected more losses
    #    610k_QTWIN_control             AMFS_case          AMFS_control 
    #                  37                    16                    12 
    #                 ASW                  CEPH                   CEU 
    #                  61                     0                     0 
    #                CHB                   CHS                   CLM 
    #                  97                   100                    60 
    #        ENDO_control    EPIGENE_QTWIN_case EPIGENE_QTWIN_control 
    #                  55                    14                    72 
    #                 FIN                   GBR       HEIDELBERG_case 
    #                  72                     0                    36   <<loses 36 HEI case which are def EUROPEAN. 
    #  HEIDELBERG_control          IBD_controls                   IBS 
    #                   0                    54                     0 
    #                 JPT                   LWK                   MXL 
    #                  89                    97                    66 
    #       OAG1_controls         OAG2_controls                   PUR 
    #                  27                    30                    55 <<loses PUR
    #     QMEGA_610k_case       QMEGA_omni_case           SDH_control 
    #                   3                     6                    19 
    #                 TSI    WAMHS_OAG_IBD_case                   YRI 
    #                   0                    27                    88 

  #2sd from NW group 
    #   610k_QTWIN_control             AMFS_case          AMFS_control 
    #                   15                     5                     3 
    #                  ASW                  CEPH                   CEU 
    #                   61                     0                     0 
    #                  CHB                   CHS                   CLM 
    #                   97                   100                    59 
    #         ENDO_control    EPIGENE_QTWIN_case EPIGENE_QTWIN_control 
    #                   17                     7                    48 
    #                  FIN                   GBR       HEIDELBERG_case 
    #                    8                     0                    28 
    #   HEIDELBERG_control          IBD_controls                   IBS 
    #                    0                    36                     1 
    #                  JPT                   LWK                   MXL 
    #                   89                    97                    65 
    #        OAG1_controls         OAG2_controls                   PUR 
    #                   19                    23                    55 
    #      QMEGA_610k_case       QMEGA_omni_case           SDH_control 
    #                    3                     3                    13 
    #                  TSI    WAMHS_OAG_IBD_case                   YRI 
    #                    1                    14                    88  <<starting to lose TSI, IBS, 
    
  #3s from NW group
    #   610k_QTWIN_control             AMFS_case          AMFS_control 
    #                    0                     0                     0  <<EUR 610k and AMFS only plot support these. Given past cleaning really discrete with EUR pops
    #                  ASW                  CEPH                   CEU 
    #                   61                     0                     0 
    #                  CHB                   CHS                   CLM 
    #                   97                   100                    58 
    #         ENDO_control    EPIGENE_QTWIN_case EPIGENE_QTWIN_control 
    #                    0                     2                    18  <<EUR EPIGENE only plot matches this. Most overlp, but there are a few proper EUR outl. Matches breslow
    #                  FIN                   GBR       HEIDELBERG_case 
    #                    0                     0                     8 <<HEI look better than I remember, to be honest. Overlap each other. As noted last time there seems to abe smattering of controls matching more with IBS/TSI so might be southern vs northern germany
    #   HEIDELBERG_control          IBD_controls                   IBS 
    #                    0                    18                     0 
    #                  JPT                   LWK                   MXL 
    #                   89                    97                    64 
    #        OAG1_controls         OAG2_controls                   PUR 
    #                   13                    15                    53 
    #      QMEGA_610k_case       QMEGA_omni_case           SDH_control 
    #                    0                     1                     5 <<QMEGA and SDH overlap each other and GBR/CEU very well. the few outliers are really outliers - looks like asian or middle eastern admixture.
    #                  TSI    WAMHS_OAG_IBD_case                   YRI 
    #                    0                     5                    88 <<there were a handful of WAMHS outliers - middle eastern background. Plot not as helpful due to the few OAG people that look to be near 100% asian or african, and thus screw with te plot boundaries. IBD look to be more close to TSI/IBS than WAMHS etc

  #6sd from NW
    #   610k_QTWIN_control             AMFS_case          AMFS_control 
    #                    0                     0                     0  <<previously cleaned on EUR so not surprising 0
    #                  ASW                  CEPH                   CEU 
    #                   61                     0                     0 
    #                  CHB                   CHS                   CLM 
    #                   97                   100                    48 
    #         ENDO_control    EPIGENE_QTWIN_case EPIGENE_QTWIN_control 
    #                    0                     0                     7 
    #                  FIN                   GBR       HEIDELBERG_case 
    #                    0                     0                     1 
    #   HEIDELBERG_control          IBD_controls                   IBS 
    #                    0                    13                     0 
    #                  JPT                   LWK                   MXL 
    #                   89                    97                    63 
    #        OAG1_controls         OAG2_controls                   PUR 
    #                    9                    11                    41 
    #      QMEGA_610k_case       QMEGA_omni_case           SDH_control 
    #                    0                     0                     1 
    #                  TSI    WAMHS_OAG_IBD_case                   YRI 
    #                    0                     1                    88  <<there is one extreme WAMHS outlier form EUR - Australian Aboriginal

  #in the breslow work, and after a bunch of back of forth using 3sd NW; Stuart thinks could easily use 6sd so it makes little difference, other than not having papers with different requirements, so use 3sd
  echo -e "\nUsing NW group (GBR_CEU_FIN_) 3sd to determine non EUR samples\n"
  wc -l ${dirK}/${PCA_out}_GBR_CEU_FIN_outlier_3sd.txt
  #86 /working/lab_stuartma/scratch/matthewL/melanoma_meta_analysis/cleaned/Melanoma_clean_2016_merged_GBR_CEU_FIN_outlier_3sd.tx

  echo -e "\nPerforming IBD\n"
  plink_1.90 --threads 1 --bfile ${dirI}/temp_merge_pruned --genome --min 0.1 --out ${dirI}/temp_merge_pruned --silent
  echo -e "\nIBD checks complete...\n"
  wc -l ${dirI}/temp_merge_pruned.genome
     #604
  awk '$10>=0.15' ${dirI}/temp_merge_pruned.genome | wc -l
     #370
  #look for anything weird - someone related to lots of people - manually filter these to stop excess people being removed
  echo -e "\mDisplaying people in col1,2 who match mutliple other people at IBD > 0.15; manually will filter these\n"
  awk 'NR==FNR {a[$2];next} $2 in a' <(awk '$10>=0.15 {print $2}' ${dirI}/temp_merge_pruned.genome | sort | uniq -c | awk '$1>1') ${dirI}/temp_merge_pruned.genome
    #  07MH01800H      07MH01800H            434         GACT2-2 UN    NA  0.1582  0.5249  0.3169  0.5794   0  0.857989  1.0000 16.0497 
    #  07MH01800H      07MH01800H            498         GACT2-7 UN    NA  0.2917  0.4720  0.2363  0.4723   0  0.828706  1.0000  9.0363 <<
    #  07MH01800H      07MH01800H            599        GACT2-26 UN    NA  0.7029  0.2971  0.0000  0.1486   0  0.740696  1.0000  3.0350
    #         294     KH09_04_A12            302     KH09_04_B08 UN    NA  0.5536  0.4464  0.0000  0.2232   1  0.756978  1.0000  3.9598
    #         294     KH09_04_A12            729     KH09_08_H06 UN    NA  0.5205  0.4721  0.0073  0.2434   1  0.761943  1.0000  4.9183
    #         358        GTas2-57            359        GTas2-68 UN    NA  0.7301  0.2699  0.0000  0.1350  -1  0.736608  1.0000  2.7730
    #         358        GTas2-57            438        GTas2-19 UN    NA  0.7695  0.2276  0.0029  0.1167  -1  0.734006  1.0000  2.6833
    #         358        GTas2-57            466        GTas2-21 UN    NA  0.6913  0.3000  0.0087  0.1587  -1  0.743575  1.0000  3.2247
    #         358        GTas2-57            487         GTas2-5 UN    NA  0.0003  0.9916  0.0080  0.5039  -1  0.818770  1.0000      NA
    #         358        GTas2-57            507        GTas2-53 UN    NA  0.2021  0.4994  0.2985  0.5482  -1  0.849829  1.0000 13.8641
    #         358        GTas2-57            562        GIST2-23 UN    NA  0.5553  0.3622  0.0825  0.2636  -1  0.771889  1.0000  3.9950
    #         358        GTas2-57            563        GIST2-55 UN    NA  0.4274  0.5566  0.0160  0.2943  -1  0.773673  1.0000  5.6023
    #         359        GTas2-68            438        GTas2-19 UN    NA  0.7046  0.2869  0.0085  0.1519  -1  0.742088  1.0000  3.0160
    #         359        GTas2-68            466        GTas2-21 UN    NA  0.7499  0.2370  0.0131  0.1316  -1  0.737991  1.0000  2.9140
    #         359        GTas2-68            507        GTas2-53 UN    NA  0.7024  0.2870  0.0106  0.1541  -1  0.742720  1.0000  2.9495
    #         359        GTas2-68            562        GIST2-23 UN    NA  0.7643  0.2275  0.0082  0.1219  -1  0.735520  1.0000  2.7761
    #         360        GTas3-11            565        GIST3-13 UN    NA  0.2474  0.5556  0.1970  0.4748  -1  0.826357  1.0000 11.3704
    #         360        GTas3-11         AG-323          AG-323 UN    NA  0.6440  0.3518  0.0042  0.1801  -1  0.747901  1.0000  3.3741
    #         361         GTas4-4            362        GTas4-10 UN    NA  0.0013  0.9908  0.0080  0.5033  -1  0.818651  1.0000 1408.5000
    #         361         GTas4-4            549        GTas4-30 UN    NA  0.3195  0.4524  0.2281  0.4543  -1  0.824170  1.0000  7.7585
    #         361         GTas4-4            550        GTas4-34 UN    NA  0.7503  0.2462  0.0035  0.1266  -1  0.736194  1.0000  2.8720
    #         361         GTas4-4            551        GTas4-41 UN    NA  0.5018  0.4872  0.0111  0.2547  -1  0.764674  1.0000  4.9225
    #         361         GTas4-4            566        GIST4-16 UN    NA  0.2596  0.4400  0.3004  0.5204  -1  0.843926  1.0000 11.3043
    #         361         GTas4-4         AG-324          AG-324 UN    NA  0.0003  0.9849  0.0148  0.5072  -1  0.819997  1.0000      NA
    #         362        GTas4-10            549        GTas4-30 UN    NA  0.5989  0.3950  0.0061  0.2036  -1  0.753168  1.0000  3.5806
    #         362        GTas4-10            551        GTas4-41 UN    NA  0.7396  0.2571  0.0033  0.1319  -1  0.737334  1.0000  2.8939
    #         362        GTas4-10            566        GIST4-16 UN    NA  0.4806  0.5108  0.0086  0.2640  -1  0.766534  1.0000  4.8115
    #         362        GTas4-10         AG-324          AG-324 UN    NA  0.2320  0.4738  0.2942  0.5311  -1  0.845798  1.0000 12.8009
    #         363        GTas12-2            439        GTas12-5 UN    NA  0.2170  0.5307  0.2522  0.5176  -1  0.839758  1.0000 11.5429
    #         363        GTas12-2            567        GIST12-1 UN    NA  0.0006  0.9862  0.0132  0.5063  -1  0.819674  1.0000 2813.0000
    #         366       GTas37-24            367       GTas37-29 UN    NA  0.2812  0.4918  0.2269  0.4729  -1  0.828138  1.0000  9.8085
    #         366       GTas37-24            463       GTas417-1 UN    NA  0.7114  0.2765  0.0121  0.1503  -1  0.741999  1.0000  3.1799
    #         366       GTas37-24         AG-332          AG-332 UN    NA  0.3845  0.6106  0.0049  0.3102  -1  0.776316  1.0000  6.1400
    #         371        GTas53-1            372        GTas53-7 UN    NA  0.2803  0.4358  0.2839  0.5018  -1  0.838645  1.0000  8.8350
    #         371        GTas53-1            373       GTas53-11 UN    NA  0.7642  0.2209  0.0150  0.1254  -1  0.736781  1.0000  2.7962
    #         371        GTas53-1            488        GTas53-4 UN    NA  0.4917  0.5018  0.0065  0.2574  -1  0.764925  1.0000  5.2154
    #         376        GTas73-1            377        GTas73-6 UN    NA  0.2561  0.4867  0.2572  0.5006  -1  0.836418  1.0000 10.8069
    #         376        GTas73-1         AG-348          AG-348 UN    NA  0.2087  0.4545  0.3368  0.5640  -1  0.856112  1.0000 13.6150
    #         382       GTas149-4            383       GTas149-7 UN    NA  0.2625  0.4949  0.2426  0.4901  -1  0.833045  1.0000 10.5602
    #         382       GTas149-4         AG-370          AG-370 UN    NA  0.2184  0.4418  0.3398  0.5607  -1  0.855609  1.0000 12.8009
    #         384       GTas156-1            449       GTas156-2 UN    NA  0.1924  0.5438  0.2638  0.5357  -1  0.844551  1.0000 13.3009
    #         384       GTas156-1          E8382       E83820001 UN    NA  0.6730  0.3241  0.0029  0.1650  -1  0.744515  1.0000  3.1939
    #         393       GTas230-6         AG-180          AG-180 UN    NA  0.4311  0.5654  0.0035  0.2862  -1  0.770989  1.0000  5.4978
    #         393       GTas230-6         AG-392          AG-392 UN    NA  0.0006  0.9922  0.0071  0.5032  -1  0.818570  1.0000 1403.0000
    #         397       GTas305-5            467       GTas305-2 UN    NA  0.4864  0.5136  0.0000  0.2568  -1  0.763491  1.0000  4.4954
    #         397       GTas305-5         AG-419          AG-419 UN    NA  0.1842  0.5081  0.3077  0.5617  -1  0.853466  1.0000 15.0207
    #         398      GTas309-10            489       GTas309-6 UN    NA  0.1862  0.5070  0.3067  0.5603  -1  0.853069  1.0000 15.5187
    #         398      GTas309-10         AG-315          AG-315 UN    NA  0.1815  0.5365  0.2820  0.5503  -1  0.849068  1.0000 14.0000
    #         427        Gvic1-37            435         Gvic1-8 UN    NA  0.4912  0.5007  0.0081  0.2585  -1  0.765286  1.0000  4.8235
    #         427        Gvic1-37            506        GVic1-38 UN    NA  0.2496  0.4642  0.2863  0.5184  -1  0.842430  1.0000 11.5911
    #         427        Gvic1-37            602        Gvic1-39 UN    NA  0.2786  0.4442  0.2773  0.4994  -1  0.837627  1.0000  9.8551
    #         427        Gvic1-37            603        Gvic1-40 UN    NA  0.2602  0.4646  0.2752  0.5075  -1  0.839249  1.0000 10.7045
    #         435         Gvic1-8            506        GVic1-38 UN    NA  0.3143  0.6857  0.0000  0.3429  -1  0.781722  1.0000  7.7380
    #         435         Gvic1-8            602        Gvic1-39 UN    NA  0.5785  0.4150  0.0066  0.2140  -1  0.755487  1.0000  3.9570
    #         435         Gvic1-8            603        Gvic1-40 UN    NA  0.0003  0.9893  0.0104  0.5050  -1  0.819202  1.0000 2817.0000
    #         478          AG0318         AG-020          AG-020 UN    NA  0.2299  0.4551  0.3150  0.5425  -1  0.849817  1.0000 11.3730
    #         478          AG0318         AG-063          AG-063 UN    NA  0.0000  1.0000  0.0000  0.5000  -1  0.816581  1.0000      NA
    #     480.013         480.013        480.014         480.014 UN    NA  0.0000  0.0000  1.0000  1.0000  -1  1.000000  1.0000      NA
    #     480.013         480.013          81470         8147001 UN    NA  0.0000  0.0000  1.0000  1.0000  -1  1.000000  1.0000      NA
    #         506        GVic1-38            602        Gvic1-39 UN    NA  0.2325  0.4772  0.2903  0.5289  -1  0.845021  1.0000 11.7418
    #         506        GVic1-38            603        Gvic1-40 UN    NA  0.1597  0.4510  0.3893  0.6148  -1  0.871035  1.0000 16.8427
    #         507        GTas2-53            562        GIST2-23 UN    NA  0.5471  0.4045  0.0484  0.2506  -1  0.766538  1.0000  4.0906
    #         507        GTas2-53            563        GIST2-55 UN    NA  0.0006  0.9791  0.0202  0.5098  -1  0.820963  1.0000      NA
    #         549        GTas4-30            550        GTas4-34 UN    NA  0.6877  0.3066  0.0057  0.1590  -1  0.743421  1.0000  2.9048
    #         549        GTas4-30            551        GTas4-41 UN    NA  0.6008  0.3902  0.0090  0.2041  -1  0.753489  1.0000  3.6827
    #         549        GTas4-30            566        GIST4-16 UN    NA  0.1717  0.5168  0.3115  0.5699  -1  0.855519  1.0000 16.7874
    #         549        GTas4-30         AG-324          AG-324 UN    NA  0.5256  0.4744  0.0000  0.2372  -1  0.759347  1.0000  4.1278
    #     657.001         657.001          83800         8380003 UN    NA  0.0003  0.9972  0.0025  0.5011  -1  0.817752  1.0000 2770.0000
    #     657.001         657.001          83800         8380004 UN    NA  0.0003  0.9840  0.0156  0.5077  -1  0.820158  1.0000      NA
    #        2382         NA19672           m008         NA19660 UN    NA  0.2400  0.4719  0.2882  0.5241  -1  0.843816  1.0000 11.9915
    #        2382         NA19672           m011         NA19685 UN    NA  0.4299  0.5701  0.0000  0.2850  -1  0.766607  1.0000  5.2790
    #        2382         NA19672           m012         NA19664 UN    NA  0.4757  0.5243  0.0000  0.2622  -1  0.762436  1.0000  4.7609
    #    6247.001        6247.001          83890         8389003 UN    NA  0.0000  0.9868  0.0132  0.5066  -1  0.819746  1.0000      NA
    #    6247.001        6247.001          83890         8389004 UN    NA  0.0000  0.9994  0.0006  0.5003  -1  0.817452  1.0000      NA
    #       61589           61589          88833         8883303 UN    NA  0.0000  0.9974  0.0026  0.5013   0  0.817813  1.0000      NA
    #       61589           61589          88833         8883304 UN    NA  0.0000  0.9987  0.0013  0.5007   0  0.817577  1.0000      NA
    #      AG-103          AG-103         AG-104          AG-104 UN    NA  0.5379  0.4525  0.0097  0.2359  -1  0.760481  1.0000  4.1429
    #      AG-103          AG-103         AG-105          AG-105 UN    NA  0.0000  0.9965  0.0035  0.5018  -1  0.817986  1.0000      NA
    #      AG-107          AG-107       AG-107.4        AG-107.4 UN    NA  0.2065  0.4849  0.3086  0.5510  -1  0.851193  1.0000 12.9276
    #      AG-107          AG-107       AG-107.5        AG-107.5 UN    NA  0.0000  0.9823  0.0177  0.5088  -1  0.820565  1.0000      NA
    #      AG-107          AG-107       AG-107.6        AG-107.6 UN    NA  0.0000  0.9991  0.0009  0.5004  -1  0.817493  1.0000      NA
    #    AG-107.4        AG-107.4       AG-107.5        AG-107.5 UN    NA  0.4425  0.5389  0.0186  0.2881  -1  0.772513  1.0000  5.6810
    #    AG-107.4        AG-107.4       AG-107.6        AG-107.6 UN    NA  0.4512  0.5488  0.0000  0.2744  -1  0.766928  1.0000  5.1970
    #        m004         NA19675           m009         NA19678 UN    NA  0.0000  1.0000  0.0000  0.5000  -1  0.812890  1.0000 1419.5000
    #        m004         NA19675           m009         NA19679 UN    NA  0.0000  1.0000  0.0000  0.5000  -1  0.811503  1.0000 883.6667
    #        m008         NA19660           m011         NA19685 UN    NA  0.0000  1.0000  0.0000  0.5000  -1  0.815864  1.0000 1379.5000
    #        m008         NA19660           m012         NA19664 UN    NA  0.5996  0.3948  0.0055  0.2029  -1  0.752992  1.0000  3.9630
    #     NA19313         NA19313        NA19331         NA19331 UN    NA  0.0019  0.9444  0.0537  0.5259  -1  0.826931  1.0000 557.8000
    #     NA19313         NA19313        NA19334         NA19334 UN    NA  0.5134  0.4240  0.0625  0.2746  -1  0.772804  1.0000  4.2205
    #     NA19355         NA19355        NA19434         NA19434 UN    NA  0.7267  0.2298  0.0435  0.1584  -1  0.746084  1.0000  2.8204
    #     NA19355         NA19355        NA19444         NA19444 UN    NA  0.6884  0.2588  0.0527  0.1822  -1  0.751937  1.0000  2.9238
    #     NA19380         NA19380        NA19381         NA19381 UN    NA  0.6231  0.3239  0.0530  0.2149  -1  0.759094  1.0000  3.5642
    #     NA19380         NA19380        NA19382         NA19382 UN    NA  0.4150  0.5303  0.0547  0.3199  -1  0.782105  1.0000  5.6508
    #     NA19434         NA19434        NA19444         NA19444 UN    NA  0.2528  0.4841  0.2631  0.5051  -1  0.837832  1.0000  9.5156
    #     NA19434         NA19434        NA19453         NA19453 UN    NA  0.5781  0.3607  0.0612  0.2415  -1  0.765505  1.0000  3.5889
    #     NA19443         NA19443        NA19469         NA19469 UN    NA  0.5044  0.4405  0.0551  0.2753  -1  0.772425  1.0000  4.7529
    #     NA19443         NA19443        NA19470         NA19470 UN    NA  0.2371  0.4731  0.2898  0.5264  -1  0.844432  1.0000 11.6214
    #       SH028         HG00501          SH032         HG00512 UN    NA  0.2895  0.4373  0.2731  0.4918  -1  0.835676  1.0000  8.4089
    #       SH028         HG00501          SH036         HG00524 UN    NA  0.2300  0.4898  0.2802  0.5251  -1  0.843449  1.0000 10.1654
    #       SH052         HG00577          SH054         HG00584 UN    NA  0.2056  0.5151  0.2793  0.5369  -1  0.845949  1.0000 13.2297
    #       SH052         HG00577          SH058         HG00595 UN    NA  0.7367  0.2085  0.0549  0.1591  -1  0.747067  1.0000  2.7308
    #       SH052         HG00578          SH053         HG00581 UN    NA  0.1969  0.4987  0.3044  0.5538  -1  0.851482  1.0000 12.7880
    #       SH052         HG00578          SH062         HG00607 UN    NA  0.6952  0.2514  0.0534  0.1791  -1  0.751321  1.0000  3.0385
    #       SH052         HG00578          SH071         HG00635 UN    NA  0.2252  0.5120  0.2628  0.5188  -1  0.840806  1.0000 11.5897
    #       SH053         HG00581          SH062         HG00607 UN    NA  0.7463  0.2004  0.0533  0.1535  -1  0.745728  1.0000  2.8008
    #       SH053         HG00581          SH071         HG00635 UN    NA  0.2387  0.5394  0.2219  0.4916  -1  0.831848  1.0000 12.1300
 
  #Manually used these above list over extreme overlapping people to make a list of people related to multiple people that can be more efficiently removed (parents etc) then use this to further filter the data; should limit the number of unexessary removals
  echo "07MH01800H 07MH01800H" > ${dirI}/IBD_common_people_to_remove
  echo "294 KH09_04_A12" >> ${dirI}/IBD_common_people_to_remove
  echo "358 GTas2-57" >> ${dirI}/IBD_common_people_to_remove
  echo "361 GTas4-4" >> ${dirI}/IBD_common_people_to_remove
  echo "362 GTas4-10" >> ${dirI}/IBD_common_people_to_remove
  echo "363 GTas12-2" >> ${dirI}/IBD_common_people_to_remove
  echo "371 GTas53-1" >> ${dirI}/IBD_common_people_to_remove
  echo "376 GTas73-1" >> ${dirI}/IBD_common_people_to_remove
  echo "393 GTas230-6" >> ${dirI}/IBD_common_people_to_remove
  echo "397 GTas305-5" >> ${dirI}/IBD_common_people_to_remove
  echo "398 GTas309-10" >> ${dirI}/IBD_common_people_to_remove
  echo "427 Gvic1-37" >> ${dirI}/IBD_common_people_to_remove
  echo "435 Gvic1-8" >> ${dirI}/IBD_common_people_to_remove
  echo "478 AG0318" >> ${dirI}/IBD_common_people_to_remove
  echo "480.013 480.013" >> ${dirI}/IBD_common_people_to_remove
  echo "480.014 480.014" >> ${dirI}/IBD_common_people_to_remove
  echo "506 GVic1-38" >> ${dirI}/IBD_common_people_to_remove
  echo "507 GTas2-53" >> ${dirI}/IBD_common_people_to_remove
  echo "549 GTas4-30" >> ${dirI}/IBD_common_people_to_remove
  echo "657.001 657.001" >> ${dirI}/IBD_common_people_to_remove
  echo "2382 NA19672" >> ${dirI}/IBD_common_people_to_remove
  echo "6247.001 6247.001" >> ${dirI}/IBD_common_people_to_remove
  echo "61589 61589" >> ${dirI}/IBD_common_people_to_remove
  echo "AG-103 AG-103" >> ${dirI}/IBD_common_people_to_remove
  echo "AG-107 AG-107" >> ${dirI}/IBD_common_people_to_remove
  echo "AG-107.4 AG-107.4" >> ${dirI}/IBD_common_people_to_remove
  echo "m004 NA19675" >> ${dirI}/IBD_common_people_to_remove
  echo "m008 NA19660" >> ${dirI}/IBD_common_people_to_remove
  echo "NA19313 NA19313" >> ${dirI}/IBD_common_people_to_remove
  echo "NA19355 NA19355" >> ${dirI}/IBD_common_people_to_remove
  echo "NA19380 NA19380" >> ${dirI}/IBD_common_people_to_remove
  echo "NA19434 NA19434" >> ${dirI}/IBD_common_people_to_remove
  echo "NA19443 NA19443" >> ${dirI}/IBD_common_people_to_remove
  echo "SH028 HG00501" >> ${dirI}/IBD_common_people_to_remove
  echo "SH052 HG00577" >> ${dirI}/IBD_common_people_to_remove
  echo "SH053 HG00581" >> ${dirI}/IBD_common_people_to_remove
  echo ""
  #then do the same for the second set of IDs, after removing these people
  awk 'NR==FNR {a[$1,$2];next} !($1,$2) in a' ${dirI}/IBD_common_people_to_remove <( sed 's/ \+/ /g' ${dirI}/temp_merge_pruned.genome | sed 's/^ //g' ) | awk '$10>=0.15 {print $4}' | sort | uniq -c | awk '$1>1' | wc -l
    #10
  echo -e "\nThe following people are from col3,4 and match more than one person at IBD > 0.15, will manually filter\n" #need to sort on col 4 to make eyeballing easier
  awk 'NR==FNR {a[$2];next} $4 in a'   <(awk 'NR==FNR {a[$1,$2];next} !($1,$2) in a' ${dirI}/IBD_common_people_to_remove <( sed 's/ \+/ /g' ${dirI}/temp_merge_pruned.genome | sed 's/^ //g' ) | awk '$10>=0.15 {print $4}' | sort | uniq -c | awk '$1>1' ) ${dirI}/temp_merge_pruned.genome | sed 's/ \+/ /g' | sed 's/^ //g' | sort -k4,4
    #27 
    #366 GTas37-24 AG-332 AG-332 UN NA 0.3845 0.6106 0.0049 0.3102 -1 0.776316 1.0000 6.1400
    #367 GTas37-29 AG-332 AG-332 UN NA 0.5091 0.4862 0.0047 0.2478 -1 0.762704 1.0000 4.4695
    #382 GTas149-4 AG-370 AG-370 UN NA 0.2184 0.4418 0.3398 0.5607 -1 0.855609 1.0000 12.8009
    #383 GTas149-7 AG-370 AG-370 UN NA 0.2061 0.5819 0.2120 0.5030 -1 0.833598 1.0000 12.8721
    #636 AG-1461 AG-585 AG-585 UN NA 0.0010 0.9906 0.0085 0.5038 -1 0.818781 1.0000 NA
    #AG-156 AG-156 AG-585 AG-585 UN NA 0.0006 0.9886 0.0108 0.5051 -1 0.819234 1.0000 NA
    #88816 8881603 AMI010013 AMI010013 UN NA 0.0000 1.0000 0.0000 0.5000 0 0.816582 1.0000 2798.0000
    #88816 8881604 AMI010013 AMI010013 UN NA 0.0000 0.9852 0.0148 0.5074 0 0.820050 1.0000 NA
    #384 GTas156-1 E8382 E83820001 UN NA 0.6730 0.3241 0.0029 0.1650 -1 0.744515 1.0000 3.1939
    #449 GTas156-2 E8382 E83820001 UN NA 0.6653 0.3347 0.0000 0.1673 -1 0.743534 1.0000 3.1498
    #361 GTas4-4 566 GIST4-16 UN NA 0.2596 0.4400 0.3004 0.5204 -1 0.843926 1.0000 11.3043
    #362 GTas4-10 566 GIST4-16 UN NA 0.4806 0.5108 0.0086 0.2640 -1 0.766534 1.0000 4.8115
    #549 GTas4-30 566 GIST4-16 UN NA 0.1717 0.5168 0.3115 0.5699 -1 0.855519 1.0000 16.7874
    #550 GTas4-34 566 GIST4-16 UN NA 0.7101 0.2709 0.0190 0.1544 -1 0.743404 1.0000 3.0199
    #551 GTas4-41 566 GIST4-16 UN NA 0.6239 0.3616 0.0145 0.1953 -1 0.751975 1.0000 3.6312
    #358 GTas2-57 507 GTas2-53 UN NA 0.2021 0.4994 0.2985 0.5482 -1 0.849829 1.0000 13.8641
    #359 GTas2-68 507 GTas2-53 UN NA 0.7024 0.2870 0.0106 0.1541 -1 0.742720 1.0000 2.9495
    #438 GTas2-19 507 GTas2-53 UN NA 0.7456 0.2438 0.0106 0.1325 -1 0.738012 1.0000 2.7364
    #466 GTas2-21 507 GTas2-53 UN NA 0.7596 0.2126 0.0277 0.1341 -1 0.739610 1.0000 2.9947
    #487 GTas2-5 507 GTas2-53 UN NA 0.4521 0.5479 0.0000 0.2739 -1 0.767362 1.0000 5.3153
    #SH074 HG00656 SH089 HG00702 UN NA 0.0026 0.9665 0.0310 0.5142 -1 0.822713 1.0000 908.6667
    #SH074 HG00657 SH089 HG00702 UN NA 0.0019 0.9655 0.0325 0.5153 -1 0.823068 1.0000 1352.0000
    #09MH03766H 09MH03766H M30242 M30242 UN NA 0.6010 0.3718 0.0272 0.2131 1 0.756797 1.0000 3.7240
    #11734 11734 M30242 M30242 UN NA 0.6899 0.2970 0.0131 0.1616 1 0.744536 1.0000 3.1393
    #NA19434 NA19434 NA19453 NA19453 UN NA 0.5781 0.3607 0.0612 0.2415 -1 0.765505 1.0000 3.5889
    #NA19444 NA19444 NA19453 NA19453 UN NA 0.4388 0.5065 0.0547 0.3080 -1 0.779510 1.0000 4.9526

  echo "AG-332 AG-332" >> ${dirI}/IBD_common_people_to_remove
  echo "AG-370 AG-370" >> ${dirI}/IBD_common_people_to_remove
  echo "AG-585 AG-585" >> ${dirI}/IBD_common_people_to_remove
  echo "AMI010013 AMI010013" >> ${dirI}/IBD_common_people_to_remove
  echo "E8382 E8382000" >> ${dirI}/IBD_common_people_to_remove
  echo "566 GIST4-16" >> ${dirI}/IBD_common_people_to_remove
  echo "507 GTas2-53" >> ${dirI}/IBD_common_people_to_remove
  echo "SH089 HG00702" >> ${dirI}/IBD_common_people_to_remove
  echo "M30242 M30242" >> ${dirI}/IBD_common_people_to_remove
  echo "NA19434 NA19434" >> ${dirI}/IBD_common_people_to_remove
  wc -l ${dirI}/IBD_common_people_to_remove
    #46
  #need the missingness amount to handle the rest
  echo -e "\nDropping common samples by missingness\n"
  plink_1.90 --threads 1 --bfile ${dirI}/temp_merge_pruned --missing --out ${dirI}/temp_merge_pruned --silent
  echo ""
  #number of people in the file after you drop the people who match multiple people
  awk 'NR==FNR {a[$1,$2];next} !(($1,$2) in a) && !(($3,$4) in a)' ${dirI}/IBD_common_people_to_remove  ${dirI}/temp_merge_pruned.genome | awk '$10>0.15' | wc -l
     #264; down from 370
  awk 'NR==FNR {a[$2]=$6;next} $2 in a {print $0,a[$2]}' ${dirI}/temp_merge_pruned.imiss <( awk 'NR==FNR {a[$1,$2];next} !(($1,$2) in a) && !(($3,$4) in a)' ${dirI}/IBD_common_people_to_remove  ${dirI}/temp_merge_pruned.genome) > ${dirI}/temp_merge_pruned.genome_1
  awk 'NR==FNR {a[$2]=$6;next} $4 in a {print $0,a[$4]}' ${dirI}/temp_merge_pruned.imiss ${dirI}/temp_merge_pruned.genome_1 > ${dirI}/temp_merge_pruned.genome_2
  head ${dirI}/temp_merge_pruned.genome_2
  echo ""
  #make a drop list. 17/03/2016 aargh, drop list formula is wrong - both times was printing $2, $15<$16 should print $4
  awk '{
    if($10>0.15 && $15>=$16)
      print $2;
    if($10>0.15 && $15<$16)
      print $4;
    }' ${dirI}/temp_merge_pruned.genome_2 > ${dirI}/temp_merge_pruned.genome_2_IID; wc -l ${dirI}/temp_merge_pruned.genome_2_IID
       #263

  awk 'NR==FNR {a[$1];next} $2 in a {print $1,$2}' ${dirI}/temp_merge_pruned.genome_2_IID ${dirI}/temp_merge_pruned.fam > ${dirI}/temp_merge_pruned.genome_3_IID
  cat ${dirI}/temp_merge_pruned.genome_3_IID ${dirI}/IBD_common_people_to_remove | sort | uniq >  ${dirK}/2016_cleaning_Melanoma_IBD_drop_ID; wc -l ${dirK}/2016_cleaning_Melanoma_IBD_drop_ID
    #302 samples (note 0.15 more extreme than used before) note this includes a handful of hapmap samples. 

  #quick once of test to make sure this worked
    ##plink_1.90 --threads 1 --bfile ${dirI}/temp_merge_pruned --remove ${dirK}/2016_cleaning_Melanoma_IBD_drop_ID --genome --min 0.1 --out ${dirI}/temp_merge_pruned_v2
    #--remove: 15436 people remaining.
    #awk '$10>0.15' ${dirI}/temp_merge_pruned_v2.genome
      #none

done #loop to turn on/off

########################################################################
#9.0 Drop IBD/PCA outliers
########################################################################

  #see section 11.0 but these were identified as having bad cluster plots in EPIGENE
  echo "rs9662633" > ${dirI}/EPIGENE_failed_cluster
  echo "rs41306613" >> ${dirI}/EPIGENE_failed_cluster
  echo "rs4670855" >> ${dirI}/EPIGENE_failed_cluster
  echo "rs4493574" >> ${dirI}/EPIGENE_failed_cluster
  echo "rs9443079" >> ${dirI}/EPIGENE_failed_cluster
  echo "rs61731700" >> ${dirI}/EPIGENE_failed_cluster
  echo "rs9403839" >> ${dirI}/EPIGENE_failed_cluster
  echo "rs13265018" >> ${dirI}/EPIGENE_failed_cluster
  echo "rs10107655" >> ${dirI}/EPIGENE_failed_cluster
  echo "exm779484" >> ${dirI}/EPIGENE_failed_cluster
  echo "rs61866610" >> ${dirI}/EPIGENE_failed_cluster
  echo "rs74017766" >> ${dirI}/EPIGENE_failed_cluster
  echo "exm1345519" >> ${dirI}/EPIGENE_failed_cluster
  echo "rs7276462" >> ${dirI}/EPIGENE_failed_cluster

  #likewise see 11.0 - extreme ORs and poor past imputation concordance
  echo rs17774013 > ${dirI}/HEIDELBERG_poor_past_imputation
  echo rs870585 >> ${dirI}/HEIDELBERG_poor_past_imputation
  echo rs644893 >> ${dirI}/HEIDELBERG_poor_past_imputation


#will need these for analysis but don't exclude prior to imputation as may as well impute them
awk 'NR==FNR {a[$1,$2];next} ($1,$2) in a'  /mnt/lustre/home/matthewL/Melanoma/Meta_analysis/IDs_selfreportedmel_exclude.txt ${dirI}/temp_2016_cleaning_QMEGA_610k_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff.fam | wc -l
  #128 - sounds about right, thought it was ~120 or so

for x in 1
do
  for i in temp_2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff temp_2016_cleaning_QMEGA_610k_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff temp_2016_cleaning_QMEGA_omni_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff temp_2016_cleaning_WAMHS_OAG_IBD_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff
  do
    echo -e "\nRemoving IBD and PCA outliers for ${i}\n"
    cat ${dirK}/Melanoma_clean_2016_merged_GBR_CEU_FIN_outlier_3sd.txt ${dirK}/2016_cleaning_Melanoma_IBD_drop_ID | sort | uniq > ${dirI}/temp_remove
    plink_1.90 --threads 1 --bfile ${dirI}/${i} --remove ${dirI}/temp_remove --make-bed --out ${dirI}/${i}_PCA_IBD
    rm -f ${dirI}/temp_remove
  done 

  #EPIGENE sep
  i=temp_2016_cleaning_EPIGENE_QTWIN_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff
  cat ${dirK}/Melanoma_clean_2016_merged_GBR_CEU_FIN_outlier_3sd.txt ${dirK}/2016_cleaning_Melanoma_IBD_drop_ID | sort | uniq > ${dirI}/temp_remove
  plink_1.90 --threads 1 --bfile ${dirI}/${i} --remove ${dirI}/temp_remove --make-bed --out ${dirI}/${i}_PCA_IBD --exclude ${dirI}/EPIGENE_failed_cluster

  #HEI
  i=temp_2016_cleaning_HEIDELBERG_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff
  cat ${dirK}/Melanoma_clean_2016_merged_GBR_CEU_FIN_outlier_3sd.txt ${dirK}/2016_cleaning_Melanoma_IBD_drop_ID | sort | uniq > ${dirI}/temp_remove
  plink_1.90 --threads 1 --bfile ${dirI}/${i} --remove ${dirI}/temp_remove --make-bed --out ${dirI}/${i}_PCA_IBD --exclude ${dirI}/HEIDELBERG_poor_past_imputation

done

  #AMFS: 971 phenotype values loaded from .fam. Among remaining phenotypes, 535 are cases and 427 are controls. (Must have been IBD - I have since noted there are some sibs across the QMEGA and AMFS set
  #610k 4851 phenotype values loaded from .fam. --remove: 4826 people remaining. Among remaining phenotypes, 912 are cases and 3914 are controls.. Again IBD, and there are a bunch of overlaps between the control set and other studies (people acted as controls in multiple studies in QLD etc)
  #QMEGA_omni_SDH: 1251 phenotype values loaded from .fam --remove: 1196 people remaining. 656 are cases and 540 are controls  Mix of PCA and IBD across studies
  #EPIGENE_QTWIN: 1728 phenotype values loaded from .fam --remove: 1694 people remaining 773 are cases and 921 are controls. - mix of internal and cross IBD, + PCA 284889 variants loaded --exclude: 284875 variants remaining.
  #HEIDELBERG: 2419 phenotype values loaded --remove: 2404 people remaining 1189 are cases and 1215 are controls. All internal IBD and some PCA - no cross study IBD 208305 variants loaded --exclude: 208302 variants remaining
  #WAMHS_OAG_IBD: 3424 phenotype values loaded --remove: 3217 people remaining 1237 are cases and 1980 are controls


  #previous meta figures (no cross study IBD check and high IBD threshold, more generous PCA, plus now look for abnormal het scores which removed 18 HEI etc). So lost a few more but nothing very concerning.
  #Study	Cases	Controls
  #AMFS	549	430
  #Q-MEGA_610k	926	3,789
  #Q-MEGA_omni	691	553
  #WAMHS	1,253	2,053
  #Essen-Heidelberg	1,212	1,219  

########################################################################
#10.0 HRC or 1000G Imputation preparation and checking 
########################################################################
  #Michigin server suggests you use a tool from http://www.well.ox.ac.uk/~wrayner/tools/ to check your data. it
    #Checks: strand, alleles, position, Ref/Alt assignments and frequency differences. In addition to the reference file v4 and above require the plink .bim and (from the plink --freq command) .frq files. Produces: set of plink commands to update or remove SNPs (see below and changes to V4.2.2) based on the checks  as well as a file (FreqPlot) of cohort allele frequency vs reference panel allele frequency.
      #Seem to remove SNPS with more than 20% MAF difference to ref panel, which makes sense. Use that default setting
  #downloaded the v4.2.5, is in /mnt/lustre/home/matthewL/bin/HRC-1000G-check-bim.pl. Checked, HRC ref panel is hg19

for x in 1
do
  for i in temp_2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD temp_2016_cleaning_QMEGA_610k_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD temp_2016_cleaning_QMEGA_omni_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD temp_2016_cleaning_EPIGENE_QTWIN_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD temp_2016_cleaning_HEIDELBERG_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD temp_2016_cleaning_WAMHS_OAG_IBD_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD
  do
    echo -e "\nAligning ${i} to HRC using HRC-1000G-check-bim.pl\n"
    plink_1.90 --threads 1 --bfile ${dirI}/${i} --freq --out ${dirI}/${i} --silent
    cd $dirI  #note the perl function makes files in the current directory if testing...
    perl /mnt/lustre/home/matthewL/bin/HRC-1000G-check-bim.pl -b ${dirI}/${i}.bim -f ${dirI}/${i}.frq -r /working/lab_stuartma/puyaG/HRC.r1-1.GRCh37.wgs.mac5.sites.tab  -h -t 0.2
    #check the exclusion list
    echo ""
    wc -l ${dirI}/Exclude-${i}-HRC.txt
    echo ""
    #looks good - just need to fix the command in the plink file
    sed 's/plink --bfile/plink_1.90 --threads 1 --bfile/g' ${dirI}/Run-plink.sh >  ${dirI}/Run-plink_v2.sh
    cd $dirI; chmod 744 Run-plink_v2.sh
    ./Run-plink_v2.sh
    #rename and move to a discrete folder. Looks like it makes final chr files; all else can be removed. Actually looks like it does a full conversion (out ${i}-updated.bim etc) then a filal chromosome split
    outname=`echo ${i} | sed 's/temp_//g'`
    echo ""
    plink_1.90 --threads 1 --bfile ${dirI}/${i}-updated --make-bed --out ${dirK}/HRC_upload_cleaned_files_lower_MAF/${outname}-updated --silent
    rm -f ${dirI}/${i}-updated.*
    echo ""
    for w in {1..23}
    do
      plink_1.90 --threads 1 --bfile ${dirI}/${i}-updated-chr${w} --make-bed --out ${dirK}/HRC_upload_cleaned_files_lower_MAF/${outname}-updated-chr${w} --silent
      rm -f ${dirI}/${i}-updated-chr${w}.*
    done
    echo -e "\nAligned and cleaned files moved to ${dirK}/HRC_upload_cleaned_files_lower_MAF. Now making a second send with MAF >=0.01 SNPs only."
    #check check if this set is worth doing a low maf version
    var_rare=`awk 'NR>1 && $5<0.01' ${dirI}/${i}.frq | wc -l`
      
    echo -e "\nFor ${i} there are ${var_rare} SNPs with MAF < 0.01 and > 0.001\n"

    #in the form ${i}-updated-chr23.fam. The default files is MAF 0.001 so just make 0.01 
    for q in {1..23}
    do
      plink_1.90 --threads 1 --bfile ${dirK}/HRC_upload_cleaned_files_lower_MAF/${outname}-updated-chr${q} --maf 0.01 --make-bed --out ${dirK}/HRC_upload_cleaned_files/${outname}_maf01-updated-chr${q} --silent
    done
    echo -e "\n${i} has been aligned to 1KG and final clean files have been made"
    #clean up
    rm -f ${dirI}/Run-plink.sh ${dirI}/Run-plink_v2.sh ${dirI}/Force-Allele1-${i}-HRC.txt ${dirI}/Chromosome-${i}-HRC.txt ${dirI}/Exclude-${i}-HRC.txt ${dirI}/FreqPlot-${i}-HRC.txt ${dirI}/ID-${i}-HRC.txt ${dirI}/LOG-${i}-HRC.txt ${dirI}/Position-${i}-HRC.txt ${dirI}/Strand-Flip-${i}-HRC.txt ${dirI}/${i}.frq
  done #loop through sets
done #loop to turn on/off

  #AMFS 
    #3624 Exclude-${i}-HRC.txt
    #Position Matches  ID matches HRC 784330  ID Doesn't match HRC 11790  Total Position Matches 796120
    #ID Match  Different position to HRC 8 No Match to HRC 834 Skipped (X, XY, Y, MT) 282 Total in bim file 797371 Total processed 797244
    #Indels (ignored in r1) 126
    #SNPs not changed 229166 SNPs to change ref alt 565026 Strand ok 794067 Total Strand ok 794192
    #Strand to change 221 Total checked 796128 Total checked Strand 794288
    #Total removed for allele Frequency diff > 0.2 541
    #Palindromic SNPs with Freq > 0.4 11
    #Non Matching alleles 1829 ID and allele mismatching 1007; where HRC is . 987
    #Duplicates removed 1
    #834 282 126 541 11 1829 1 = 3624. So "No Match to HRC", "Skipped (X, XY, Y, MT):, "Indels" "allele Frequency diff > 0.2", Palindromic SNPs with Freq > 0.4, "Non Matching alleles 1829", "Duplicates removed" all seem to be excluded
    
    #1799 rare variants
  
  #610k
    #Position Matches  ID matches HRC 492863  ID Doesn't match HRC 2740  Total Position Matches 495603
    #ID Match  Different position to HRC 2 No Match to HRC 87 Skipped (X, XY, Y, MT) 23 Total in bim file 495715 Total processed 495715 Indels (ignored in r1) 0
    #SNPs not changed 150772 SNPs to change ref alt 343745 Strand ok 494517 Total Strand ok 494517
    #Strand to change 1 Total checked 495605 Total checked Strand 494518
    #Total removed for allele Frequency diff > 0.2 137
    #Palindromic SNPs with Freq > 0.4 2
    #Non Matching alleles 1085 ID and allele mismatching 604; where HRC is . 603 Duplicates removed 0
      #1085 + 2 + 137 + 87 + 23  = 1334
      #exclude list is 1334, nice
    #36 rare SNPs!

  #QMEGA_omni_SDH
    #Position Matches  ID matches HRC 771380  ID Doesn't match HRC 10732  Total Position Matches 782112
    #ID Match  Different position to HRC 1 No Match to HRC 258 Skipped (X, XY, Y, MT) 243 Total in bim file 782614 Total processed 782614 Indels (ignored in r1) 0
    #SNPs not changed 225705 SNPs to change ref alt 554547 Strand ok 780165 Total Strand ok 780252
    #Strand to change 146 Total checked 782113 Total checked Strand 780311 
    #Total removed for allele Frequency diff > 0.2 367 Palindromic SNPs with Freq > 0.4 11
    #Non Matching alleles 1791 ID and allele mismatching 991; where HRC is . 976 Duplicates removed 0
    #1791 + 11 + 367 + 243 + 258  = 2670
      #Exclude list is 2670
    #2291 rare SNPs

  #EPIGENE/QTWIN (updated for removeing the ~14 bad cluster SNPs but miised the values after total processed = 284676 so will be out by ~14
    #Position Matches ID matches HRC 274409  ID Doesn't match HRC 9963  Total Position Matches 284372
    #ID Match  Different position to HRC 5 No Match to HRC 298 Skipped (X, XY, Y, MT) 1 Total in bim file 284875 Total processed 284676 Indels (ignored in r1) 199
    #SNPs not changed 82360 SNPs to change ref alt 201411 Strand ok 283771 Total Strand ok 283771
    #Strand to change 0 Total checked 284390 Total checked Strand 283771 
    #Total removed for allele Frequency diff > 0.2 52 Palindromic SNPs with Freq > 0.4 2
    #Non Matching alleles 617 ID and allele mismatching 359; where HRC is . 348 Duplicates removed 0
    #617 + 2 + 52 + 199 + 1 + 299 = 1170
      #1170 exclude list
    #22673 rare SNPs



  #HEIDELBERG
    #getting a weird pythin error -- see if there is something that explains it at line 208305 of the bim/freq file? also look at lines 555 and 570 in the .pl
    #Position Matches  ID matches HRC 207273  ID Doesn't match HRC 891  Total Position Matches 208164
    #ID Match  Different position to HRC 0 No Match to HRC 49 Skipped (X, XY, Y, MT) 92  Total in bim file 208305 Total processed 208305 Indels (ignored in r1) 0
    #SNPs not changed 63721 SNPs to change ref alt 143930 Strand ok 207544 Total Strand ok 207651 Strand to change 170 Total checked 208164 Total checked Strand 207714
    #Total removed for allele Frequency diff > 0.2 62  
    ####Use of uninitialized value $palin in concatenation (.) or string at /mnt/lustre/home/matthewL/bin/HRC-1000G-check-bim.pl line 555, <IN> line 208305.
    #Palindromic SNPs with Freq > 0.4
   #Non Matching alleles 450 ID and allele mismatching 254; where HRC is . 254 Duplicates removed 0
   #######Use of uninitialized value $palin in concatenation (.) or string at /mnt/lustre/home/matthewL/bin/HRC-1000G-check-bim.pl line 570, <IN> line 208305.
    #450 + 62 + 92 49 = 653.
      #exclude list 653 so exclude list okay? Script ran so suggests python didn't fail. Does "Use of uninitialized value $palin in concatenation" just mean it tried to call an non existant variable, as there are 0 of that class? Quick google suggests that is the case - this error happens when a variable is undefined. Looking at the script $palin isn't given a defualt value where as indel and duplicate say are given a defaul value of 0, and dont throw an error. IGNORE.
    #2754 rare SNPs

  #WAMHS OAG IBS
    #Position Matches  ID matches HRC 580954  ID Doesn't match HRC 2230  Total Position Matches 583184
    #ID Match  Different position to HRC 1 No Match to HRC 128 Skipped (X, XY, Y, MT) 93 Total in bim file 583406 Total processed 583406 Indels (ignored in r1) 0
    #SNPs not changed 170867 SNPs to change ref alt 411028 Strand ok 581894 Total Strand ok 581895
    #Strand to change 3 Total checked 583185 Total checked Strand 581897 
    #Total removed for allele Frequency diff > 0.2 166 Palindromic SNPs with Freq > 0.4 8
    #Non Matching alleles 1280 ID and allele mismatching 708; where HRC is . 706 Duplicates removed 0
       #exclude list 1675
    #13893 rare SNPs

########################################################################
#11.0 Genotype GWAS check. 
########################################################################
  #Not entirely sure abou the step in the HRC alignment that filters SNPs with an 0.2 MAF difference to EUR. This could potentially drop real SNPs, but will also get rid of nonsense...

for x in 1
do
  for i in temp_2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD temp_2016_cleaning_QMEGA_610k_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD temp_2016_cleaning_QMEGA_omni_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD temp_2016_cleaning_EPIGENE_QTWIN_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD temp_2016_cleaning_HEIDELBERG_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD temp_2016_cleaning_WAMHS_OAG_IBD_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD
  do
    outname=`echo ${i} | sed 's/temp_//g'`
    echo -e "\nRunning GWAS with data before alignment to HRC via HRC-1000G-check-bim\n"
    plink_1.90 --threads 1 --bfile ${dirI}/${i} --assoc --out ${dirI}/${i} --silent
    plink_1.90 --threads 1 --bfile ${dirI}/${i} --logistic --covar ${dirK}/Melanoma_clean_2016_merged_PCA.eigenvec_IDs --covar-name PC1 PC2 PC3 PC4 PC5 PC6 --hide-covar --out ${dirI}/${i} --silent 
    echo -e "\nUncorrected GWAS results for ${i} that are < 5e-7\n"
    awk 'NR==1 || $9<5e-7' ${dirI}/${i}.assoc
    echo -e "\n${i} results corrected for first 6 PCs. Showing 1e-6 SNPs as p values may have decreased\n"
    awk 'NR==1 ||  $9<1e-6' ${dirI}/${i}.assoc.logistic
    
    cd $dirI
    awk 'NR==1 || $9!="NA" {print $9}' ${dirI}/${i}.assoc > ${dirI}/${i}.assoc_t
    awk 'NR==1 || $9!="NA" {print $9}' ${dirI}/${i}.assoc.logistic > ${dirI}/${i}.assoc.logistic_t

    echo -e "\nRunning GWAS for the HRC aligned - some extra SNPs dropped (extreme MAF differences)"
    plink_1.90 --threads 1 --bfile ${dirK}/HRC_upload_cleaned_files_lower_MAF/${outname}-updated --assoc --out ${dirI}/${i}_v2 --silent
    plink_1.90 --threads 1 --bfile ${dirK}/HRC_upload_cleaned_files_lower_MAF/${outname}-updated --logistic --covar ${dirK}/Melanoma_clean_2016_merged_PCA.eigenvec_IDs --covar-name PC1 PC2 PC3 PC4 PC5 PC6 --hide-covar  --out ${dirI}/${i}_v2 --silent
 
    echo -e "\nUncorrected GWAS results for HRC aligned version of ${i} that are <5e-7\n"
    awk 'NR==1 || $9<5e-7' ${dirI}/${i}_v2.assoc
    echo -e "\nHRC aligned version of ${i} results corrected for first 6 PCs. Showing 1e-6 SNPs as p values may have decreased\n"
    awk 'NR==1 ||  $9<1e-6' ${dirI}/${i}_v2.assoc.logistic

    awk 'NR==1 || $9!="NA" {print $9}' ${dirI}/${i}_v2.assoc > ${dirI}/${i}_v2.assoc_t
    awk 'NR==1 || $9!="NA" {print $9}' ${dirI}/${i}_v2.assoc.logistic > ${dirI}/${i}_v2.assoc.logistic_t


module load R/3.2.2
R --no-save <<EOF
options(bitmapType='cairo')
library(ggplot2)
read.table("${dirI}/${i}.assoc_t",header=T)->dd;dim(dd)
read.table("${dirI}/${i}.assoc.logistic_t",header=T)->ee;dim(ee)
read.table("${dirI}/${i}_v2.assoc_t",header=T)->ff;dim(ff)
read.table("${dirI}/${i}_v2.assoc.logistic_t",header=T)->gg;dim(gg)

round(median(qchisq(dd\$P,1,lower.tail=FALSE))/0.456,4)

png("${dirI}/${i}_QQ_no_PCs.png"); qqcustom(dd\$P); dev.off()

round(median(qchisq(ee\$P,1,lower.tail=FALSE))/0.456,4)

png("${dirI}/${i}_QQ_6_PCs.png"); qqcustom(ee\$P) ; dev.off()

round(median(qchisq(ff\$P,1,lower.tail=FALSE))/0.456,4)

png("${dirI}/${i}_HRC_aligned_QQ_no_PCs.png") ; qqcustom(ff\$P) ; dev.off()

round(median(qchisq(gg\$P,1,lower.tail=FALSE))/0.456,4)

png("${dirI}/${i}_HRC_aligned_QQ_6_PCs.png") ; qqcustom(gg\$P) ; dev.off()

q()
EOF

  rm -f ${dirI}/${i}.assoc* ${dirI}/${i}_v2.assoc* 

  done #
done #loop to turn on/off

#AMFS
  #AMFS assoc - GIF 1.0164
  # CHR               SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR 
  #  16          rs258322   89755903    A   0.1729  0.08197    G        34.16    5.087e-09        2.341 
  #  18          rs568388   33829077    T    0.414   0.5293    C        25.36     4.76e-07       0.6284 
  #AMFS + PCs - GIF 1.0159
  # CHR               SNP         BP   A1       TEST    NMISS         OR         STAT            P 
  #  16          rs258322   89755903    A        ADD      962      2.251        5.405     6.48e-08
  #  18          rs568388   33829077    T        ADD      962     0.6258        -4.94    7.797e-
  #AMFS HRC aligned uncorrected. So not impacting the most different SNPs. GIF 1.0159 - so GIF decreases a touch. Not sure how to intepret, but any high maf diff will be likelt to have extreme values... so don't know. That said GIF dif at 4th decimal place, and would be the same with rounding.
  # CHR               SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR 
  #  16          rs258322   89755903    A   0.1729  0.08197    G        34.16    5.087e-09        2.341 
  #  18          rs568388   33829077    T    0.414   0.5293    C        25.36     4.76e-07       0.6284 
  #AMFS HRC aligned +6 PCs.  - GIF [1] 1.0155 < again 4th decimal place
  # CHR               SNP         BP   A1       TEST    NMISS         OR         STAT            P 
  #  16          rs258322   89755903    A        ADD      962      2.251        5.405     6.48e-08
  #  18          rs568388   33829077    T        ADD      962     0.6258        -4.94    7.797e-07

#610k 
  #ASSOC my cleaning alone GIF 1.0321
   # CHR          SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR 
   #4    rs7682929   44297909    A   0.1135  0.07716    C        25.37    4.742e-07        1.531 
   #9   rs10757257   21806564    A   0.3607   0.4251    G        25.29    4.939e-07       0.7631  <<CDKN2A
  #16     rs352935   89648580    T   0.4314   0.5079    C        34.63     3.99e-09        0.735 <MC1R
  #16    rs7188458   89726484    A   0.5219   0.4497    G        31.08    2.482e-08        1.336 
  #16     rs258322   89755903    A   0.1645   0.1095    G        42.41    7.411e-11        1.601 
  #16    rs1800359   89805261    A    0.335    0.401    G        27.14    1.894e-07       0.7524 
  #16    rs1800286   89869761    T    0.335   0.4019    C        27.86    1.304e-07       0.7496 
  #16   rs11861084   89875710    A   0.3388   0.4051    C        27.24    1.794e-07       0.7526 
  #16    rs8049897   90024202    A   0.2275   0.1622    G        43.75    3.725e-11        1.521 
  #16    rs4408545   90044028    T   0.4073   0.4918    C        42.33    7.722e-11       0.7102 
  #my cleaing + first 6PCs GIF 1.015
  # CHR          SNP         BP   A1       TEST    NMISS         OR         STAT            P 
  # 4    rs7682929   44297909    A        ADD     4826      1.541        4.956    7.195e-07
  # 9   rs10123637   21721939    T        ADD     4824     0.7651       -4.978    6.434e-07 
  # 9   rs10757257   21806564    A        ADD     4825     0.7586       -5.128    2.921e-07
  # 9   rs10811629   21840298    G        ADD     4824      0.767       -4.979    6.383e-07
  #16     rs352935   89648580    T        ADD     4821     0.7414       -5.626    1.845e-08
  #16    rs7188458   89726484    A        ADD     4825      1.309          5.1    3.388e-07
  #16     rs258322   89755903    A        ADD     4826      1.576        6.135    8.535e-10
  #16    rs1800286   89869761    T        ADD     4826     0.7601       -4.928    8.301e-07
  #16    rs8049897   90024202    A        ADD     4826      1.494        6.144    8.041e-10
  #16    rs4408545   90044028    T        ADD     4826     0.7168       -6.237    4.473e-10
  #HRC cleaning - getting the same results for assoc and assoc.logistic. HRC uncorrected GIF 1.0321 (identical). PC6 1.0145, diff at 4 dec. 
  
#Omni and SDH
  #my cleaning alone, assoc GIF 1.0117
   #CHR              SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR 
  #   9       rs10811595   21715801    G   0.4096   0.5186    T        28.21     1.09e-07       0.6443 
  #   9        rs6475564   21728683    T    0.545   0.4407    C        25.74    3.903e-07         1.52 
  #   9        rs7871345   21729895    T   0.5442   0.4407    C        25.37     4.74e-07        1.515 
  #   9        rs7854222   21734348    G   0.5443   0.4406    A        25.41    4.634e-07        1.516 
  #   9       rs10811600   21734802    A   0.5442   0.4407    G        25.37     4.74e-07        1.515 
  #   9        rs4636294   21747803    A   0.5442   0.4398    G        25.82    3.742e-07        1.521 
  #   9        rs1965153   21752105    T    0.545   0.4398    G         26.2    3.076e-07        1.525 
  #   9        rs7852710   21760254    T   0.5443   0.4398    C        25.84    3.712e-07        1.521 
  #   9        rs7866787   21760639    A   0.5442   0.4398    G        25.82    3.742e-07        1.521 
  #   9        rs7856941   21761322    T   0.5442   0.4397    C        25.85    3.686e-07        1.521 
  #   9        rs7850446   21763742    G   0.5457   0.4398    A        26.58    2.525e-07         1.53 
  #   9        rs6475576   21765601    C    0.545   0.4398    T         26.2    3.076e-07        1.525 
  #   9        rs3849929   21769412    G    0.545   0.4398    A        26.22     3.05e-07        1.526 
  #   9        rs7045768   21770709    T   0.5442   0.4398    C        25.82    3.742e-07        1.521 
  #   9        rs6475579   21771756    A   0.5442   0.4407    G        25.37     4.74e-07        1.515 
  #   9        rs7037577   21772037    G   0.5442   0.4397    A        25.85    3.686e-07        1.521 
  #   9       rs10757253   21772267    T   0.5442   0.4407    C        25.37     4.74e-07        1.515 
  #   9        rs1414228   21773541    A   0.5421   0.4376    C        25.74    3.898e-07        1.521 
  #  13        rs9538228   59416119    A  0.05496    0.113    G         26.7    2.375e-07       0.4567 
  #my cleaning, assoc.logistic GIF 1.0141
  # CHR              SNP         BP   A1       TEST    NMISS         OR         STAT            P 
  #   9       rs10811595   21715801    G        ADD     1192     0.6398       -5.261     1.43e-07
  #   9        rs6475564   21728683    T        ADD     1196      1.533        5.018    5.209e-07
  #   9        rs7871345   21729895    T        ADD     1196      1.529        4.987    6.125e-07
  #   9        rs7854222   21734348    G        ADD     1194      1.529        4.987    6.122e-07
  #   9       rs10811600   21734802    A        ADD     1196      1.529        4.987    6.125e-07
  #   9        rs4636294   21747803    A        ADD     1196      1.535        5.036    4.749e-07
  #   9        rs1965153   21752105    T        ADD     1196      1.539        5.069    4.008e-07
  #   9        rs7852710   21760254    T        ADD     1195       1.53        5.008    5.488e-07
  #   9        rs7866787   21760639    A        ADD     1196       1.53        5.009    5.461e-07
  #   9        rs7856941   21761322    T        ADD     1195      1.531         5.01     5.43e-07
  #   9        rs7850446   21763742    G        ADD     1196      1.537        5.068     4.01e-07
  #   9        rs6475576   21765601    C        ADD     1196      1.532        5.034    4.795e-07
  #   9        rs3849929   21769412    G        ADD     1195      1.535        5.049    4.441e-07
  #   9        rs7045768   21770709    T        ADD     1196      1.529        5.009    5.468e-07
  #   9        rs6475579   21771756    A        ADD     1196      1.526        4.973    6.594e-07
  #   9        rs7037577   21772037    G        ADD     1195       1.53            5    5.723e-07
  #   9       rs10757253   21772267    T        ADD     1196      1.526        4.973    6.594e-07
  #   9        rs2152272   21773416    A        ADD     1196      1.516          4.9    9.563e-07
  #   9        rs1414228   21773541    A        ADD     1190      1.536        5.024    5.049e-07
  #   9       rs10738599   21774467    T        ADD     1196       1.52        4.936     7.98e-07
  #   9        rs7021012   21775957    G        ADD     1196      1.518        4.902    9.498e-07
  #   9        rs7038708   21788081    T        ADD     1195      1.522         4.91    9.111e-07
  #   9        rs7022856   21788523    A        ADD     1196      1.521        4.899    9.631e-07
  #   9       rs10811617   21790067    A        ADD     1196      1.524        4.924    8.468e-07
  #chr 13 fallen off, rest still there - all CDKN2A hits though note that I think is actually the gluacoma site.... oh no, 610 down thre as well
  #HRC top SNP results for assoc identical, as are logistic.assoc. GIF 1.0122 for assoc, 1.0141 for HRC + 6PCs which is the same as above.

 #EPIGENE/QTWIN
    #assoc my cleaning GIF 0.9842 (same post cleaning bad SNPs)
    # CHR                   SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR 
    #  1             rs9662633    3649562    A 0.001951  0.03366    G        44.68     2.32e-11      0.05611  <<<GONE POST CLEANING
    #  1            rs41306613   33235789    T 0.001294  0.03094    C        43.19    4.975e-11      0.04056 <<<GONE POST CLEANING
    #  2             rs4670855   38504680    A 0.001294   0.1971    G        335.1    7.304e-75     0.005278 <<<GONE POST CLEANING 
    #  4             rs4493574   64073599    A 0.001294   0.0266    G        36.31    1.683e-09       0.0474  <<<GONE POST CLEANING
    #  6             rs9443079   75553805    T 0.001294    0.136    C        221.1    5.123e-50     0.008228  <<<GONE POST CLEANING
    #  6            rs61731700  109771691    A 0.001297  0.02009    G        26.01    3.391e-07      0.06336  <<<GONE POST CLEANING
    #  6             rs9403839  147282434    T  0.00194  0.02226    C        27.07    1.963e-07      0.08541  <<<GONE POST CLEANING
    #  8            rs13265018   22886020    A  0.01035   0.1064    G        132.2    1.342e-30      0.08782  <<<GONE POST CLEANING
    #  8            rs10107655   57994818    T 0.001294  0.08523    G        132.4    1.218e-30       0.0139  <<<GONE POST CLEANING
    #  9             exm779484  125486542    G 0.001294  0.02443    A        32.89    9.755e-09      0.05173  <<<GONE POST CLEANING
    # 10            rs61866610   95326676    T 0.003886   0.0418    C         50.5    1.192e-12      0.08942  <<<GONE POST CLEANING
    # 16            rs74017766   50703016    T  0.02975 0.006515    C        26.98    2.057e-07        4.677  <<<GONE POST CLEANING
    # 17            exm1345519   62028843    T  0.01746 0.001086    C        26.57    2.542e-07        16.35  <<<GONE POST CLEANING
    # 20              rs910873   33171772    A   0.1429   0.0874    G        25.92    3.554e-07        1.741  same
    # 21             rs7276462   31429636    T  0.00194  0.05429    C        78.14    9.597e-19      0.03387  <<<GONE POST CLEANING
    #assoc my cleaning with 6PC GIF 0.9682 <same post cleaning)
    #CHR                   SNP         BP   A1       TEST    NMISS         OR         STAT            P 
    #2             rs4670855   38504680    A        ADD     1694   0.005064       -7.427    1.107e-13<<maybe maf filtered out of Upekha's?
    #6             rs9443079   75553805    T        ADD     1692    0.00797       -6.784    1.172e-11 <<maybe maf filtered out of Upekha's?
    #8            rs13265018   22886020    A        ADD     1694    0.09017       -9.065    1.246e-19 <<extreme in UPEKHAS
    #8            rs10107655   57994818    T        ADD     1694    0.01411       -5.968      2.4e-09  <<maybe maf filtered out of Upekha's?
    #10            rs61866610   95326676    T        ADD     1693    0.09007       -5.623    1.878e-08 <<maybe maf filtered out of Upekha's?
    #21             rs7276462   31429636    T        ADD     1694    0.03128       -5.882    4.041e-09 <<maybe maf filtered out of Upekha's?
    #HRC cleaning, assoc GIF 0.9847 (same post cleaning)
      #following SNP  missing    # 17            exm1345519   62028843    T  0.01746 0.001086    C        26.57    2.542e-07        16.35 
    #HERC cleaning, logistic + 6PC ius same as above. GIF 0.9687
    #REALLY NEED TO CHECK THESE CLUSTER PLOTS - or see what happens with imutation? 

    #zgrep -w "afr\|rs4670855\|rs9443079\|rs13265018\|rs10107655\|rs61866610\|rs7276462" /working/lab_stuartma/scratch/matthewL/melanoma_meta_analysis/working_directory/ALL_1000G_phase1integrated_v3_ALL_chr_legends_with_chromosome.txt.gz | cut -d" " -f1-4,9,13
      #chr id position a0 eur.aaf eur.maf
      #2 rs4670855 38504680 G 0.1702 0.1702 <<aff freq nearl 0.001, con = 1KG <<cluster blot is borked. All three cluster are close together and in major allele plot. Drop
      #6 rs9443079 75553805 C 0.1504 0.1504 <<aff freq near 0.001, con = 1KG  << maybe two clusters close together  called as hom +some noise. Drop
      #8 rs13265018 22886020 A 0.9050 0.0950 << aff freq near 0.001, con 1KG<< not on EPIGENE as that ID, try exm689164. Yep two close hom + het called as hom, looks like hom minor called as het. hom/het to close to seperate, drop
      #8 rs10107655 57994818 G 0.0739 0.0739 << aff freq near 0.001, con 1KG << one cluster, plus some nonsense of which 2 falls in her/ Drop 
      #10 rs61866610 95326676 C 0.0475 0.0475 <<aff freq near 0.003, con 1KG <<not on EPIGENE as that ID, try exm842865... wow, trash. No signal at all. justa smear. Drop.
      #21 rs7276462 31429636 C 0.0383 0.0383 << aff freq near 0.002, con 1KG  <<maybe a hom and het cluster but too cloe to call, and falls under hom major. Drop

   
  #SINCE ALL are borked in EPIGENE -----
    #logged onto genome01, hell, still have EPIGENE_QC_v3 on there from last time. Actually no sure if that will load? sure the data isn't local? Or maybe it is. 
       #Filtered the table by these SNP IDs. Two not in there, must be under exome ID prior to conversion, need to get (got them)
          #all six have clearly borked plots, usually only 1-2 clusters close to each other falling under the hom boundary. Can exclide - need to save picture. 
             #SHould expand this to all the weird ones without PCs that also have extreme values while in here.


       #assoc my cleaning GIF 0.9842
    # CHR                   SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR 
    #  1             rs9662633    3649562    A 0.001951  0.03366    G        44.68     2.32e-11      0.05611 << cluster overlap under hom major, drop
    #  1            rs41306613   33235789    T 0.001294  0.03094    C        43.19    4.975e-11      0.04056 <<as exm41165 in raw data. Nonsense, one big cluster
    #  2             rs4670855   38504680    A 0.001294   0.1971    G        335.1    7.304e-75     0.005278   <<see above. cluster blot is borked. All three cluster are close together and in major allele plot. Drop
    #  4             rs4493574   64073599    A 0.001294   0.0266    G        36.31    1.683e-09       0.0474 <<may have a het cluster, but too close to hom to call, drop
    #  6             rs9443079   75553805    T 0.001294    0.136    C        221.1    5.123e-50     0.008228 << maybe two clusters close together  called as hom +some noise. Drop
    #  6            rs61731700  109771691    A 0.001297  0.02009    G        26.01    3.391e-07      0.06336 <<epi raw data as exm570865. One big smeared cluser, drop.
    #  6             rs9403839  147282434    T  0.00194  0.02226    C        27.07    1.963e-07      0.08541   <<might actually be three clusters, but not being called correctly (het falling between hom major and her, and looks like hom-minor falling in het). Not really worth recalling at this stage, just drop

    #  8            rs13265018   22886020    A  0.01035   0.1064    G        132.2    1.342e-30      0.08782   << in raw epi data as exm689164. two close hom + het called as hom, looks like hom minor called as het. hom/het to close to seperate, drop
    #  8            rs10107655   57994818    T 0.001294  0.08523    G        132.4    1.218e-30       0.0139  < one cluster, plus some nonsense of which 2 falls in her/ Drop 
    #  9             exm779484  125486542    G 0.001294  0.02443    A        32.89    9.755e-09      0.05173  <<one smeareed cluster, drop
    # 10            rs61866610   95326676    T 0.003886   0.0418    C         50.5    1.192e-12      0.08942 << one cluster, plus some nonsense of which 2 falls in her/ Drop
    # 16            rs74017766   50703016    T  0.02975 0.006515    C        26.98    2.057e-07        4.677 << exm1239770 in raw data. Smear. Drop
    # 17            exm1345519   62028843    T  0.01746 0.001086    C        26.57    2.542e-07        16.35 << smear, drop
    # 20              rs910873   33171772    A   0.1429   0.0874    G        25.92    3.554e-07        1.741 <<<ASIP, keep in as a refernce. exm-rs910873 in raw data. 3 nice clusters falling under the expected positions
    # 21             rs7276462   31429636    T  0.00194  0.05429    C        78.14    9.597e-19      0.03387 <<maybe a hom and het cluster but too cloe to call, and falls under hom major. Drop
 
      #So all these SNPs BUT rs910873 can be dropped. All images saved as EPIGENE...XXXX_cluster_plot.png in Document/Windows_transfer
  

    #HEI - no covariates, my cleaning.. INF 1.41 which is what I noted before. Control array is mad
   #  CHR         SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR 
   #2  rs17774013   47067214    T   0.2775  0.03669    C        529.4    3.88e-117        10.09  <<no cluster plots to check but I think we can drop that one.... GONE POST CLEANING
   #3   rs3774262  186571814    A   0.1245  0.07901    G        27.24    1.799e-07        1.657
   #8    rs870585    5615643    A   0.1491  0.08748    C         43.4    4.458e-11        1.827  <<can't presume is false, left in last time. GONE POST CLEANING
   #9    rs644893   79293847    C  0.02061 0.004115    T        26.96    2.073e-07        5.091  <<GONE POST CLEANING
   #9    rs671261  129659918    C   0.2717   0.2091    T        25.82    3.742e-07        1.411 
   #9   rs2031852  129755509    A   0.2696    0.206    G        26.74    2.329e-07        1.422 
  #13    rs283957   77117462    T   0.4348   0.3622    C        26.36    2.835e-07        1.354 

    #with PCs - inflation is 1.084 which is still pretty bad.  - no SNPs with p values < 1e-6 with PC correction. 
    #CHR         SNP         BP   A1       TEST    NMISS         OR         STAT            P 

   #check HEI past imputation
   #zgrep -a rs870585 ../../GERMAN_MELANOMA/IMPUTATION_RESULTS/Archive_GERMAN_MELANOMA_chromosome_8_imputation.tar.gz 
      #8 rs870585 5615643 0.123 0.968 0.993 2 0.872 0.891 0.405 <<imputed well but low concordance, drop it. (below 0.9 for a well imputed SNP is werid
   #zgrep -a rs3774262 ../../GERMAN_MELANOMA/IMPUTATION_RESULTS/Archive_GERMAN_MELANOMA_chromosome_3_imputation.tar.gz
      #3 rs3774262 186571814 0.102 1.000 1.000 2 0.849 0.980 0.916 < well imputed and concordant, leave
   #zgrep -a "rs644893\|rs671261\|rs2031852" ../../GERMAN_MELANOMA/IMPUTATION_RESULTS/Archive_GERMAN_MELANOMA_chromosome_9_imputation.tar.gz
      #9 rs671261 129659918 0.761 1.000 1.000 2 0.941 0.981 0.959 <<fine, keep, Two indep SNps that look to be in LD, can't both be wrong
      #9 rs2031852 129755509 0.763 0.997 0.999 2 0.982 0.995 0.990 <<fine, keep,  Two indep SNps that look to be in LD, can't both be wrong 
      #9 rs644893 79293847 0.988 1.000 1.000 2 0.675 0.991 0.644 <<concord is low but so is masked imputed (0.675). Drop.

   #zgrep -a rs283957 ../../GERMAN_MELANOMA/IMPUTATION_RESULTS/Archive_GERMAN_MELANOMA_chromosome_13_imputation.tar.gz
      #13 rs283957 77117462 0.399 0.997 0.999 2 0.859 0.920 0.848 <<decent imputation, concordance is not that bad, keep
   
  #WAMHS/ETC
  #no covariates, my cleaning - GIF 1.041
  # CHR         SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR 
  # 4   rs2244291   53468479    C   0.1289  0.08838    T        26.88    2.161e-07        1.527 <<near rs2898681, optic nerve SNP 
  # 6   rs2256965   31555130    A   0.4519   0.3866    G        26.84    2.214e-07        1.308  <<HLA? Near a chron's loci rs1799964
  #16    rs164741   89692298    A   0.3816   0.3126    G        32.29    1.327e-08        1.357  <<MC1R
  #16    rs258322   89755903    A   0.1645   0.1096    G        40.41    2.056e-10          1.6 
  #16   rs7195043   90020861    T   0.5356   0.4654    C        29.99    4.335e-08        1.325 
  #16   rs8051733   90024206    G   0.4147   0.3343    A        42.42     7.35e-11        1.411 
  #16   rs4785751   90029417    A   0.4292   0.4977    G        28.68    8.556e-08       0.7588 
  #16   rs4785752   90035141    A   0.5121   0.4386    G        33.05    8.965e-09        1.344 
  #16   rs4238833   90050689    G   0.4587   0.3804    T        38.49    5.491e-10         1.38 
  #16   rs4785759   90050880    A   0.5287   0.4495    C        38.27    6.159e-10        1.374 
  #16   rs9939542   90053048    C   0.3618   0.2831    A        43.81    3.613e-11        1.435 
  #16   rs4785763   90066936    A   0.4182   0.3442    C        35.66    2.352e-09         1.37 
  #16  rs11076650   90067941    G    0.504   0.4293    A        34.26    4.823e-09        1.351 
  #16  rs10852628   90079927    C   0.3541   0.2852    T        33.62    6.687e-09        1.374 

  #PC6, my cleaning - GIF 1.017
  # CHR         SNP         BP   A1       TEST    NMISS         OR         STAT            P 
  #   4   rs2244291   53468479    C        ADD     3217      1.499        4.923    8.533e-07
  #   6   rs2256965   31555130    A        ADD     3216      1.317        5.248    1.539e-07
  #  16    rs164741   89692298    A        ADD     3217      1.333        5.277    1.315e-07
  #  16    rs258322   89755903    A        ADD     3217       1.54        5.738    9.557e-09
  #  16   rs7195043   90020861    T        ADD     3217      1.296        5.018     5.23e-07
  #  16   rs8051733   90024206    G        ADD     3217      1.382        5.995    2.032e-09
  #  16   rs4785752   90035141    A        ADD     3215      1.344        5.557    2.738e-08
  #  16   rs4238833   90050689    G        ADD     3213      1.366        5.854    4.806e-09
  #  16   rs4785759   90050880    A        ADD     3217      1.369         5.89    3.858e-09
  #  16   rs9939542   90053048    C        ADD     3217      1.431        6.374    1.845e-10
  #  16   rs4785763   90066936    A        ADD     3215      1.358         5.67    1.432e-08
  #  16  rs11076650   90067941    G        ADD     3216      1.339        5.467    4.566e-08
  #  16  rs10852628   90079927    C        ADD     3216      1.361        5.489    4.032e-08

  #HRC cleaning - no diff in SNPs/values. No PC6 1.041, with PC6 1.017

  zgrep -a rs2244291 ../cleaned/IMPUTATION_RESULTS/CIDR_WAMHS_IMPUTATION/Archive_CIDR_WAMHS_chromosome_4_imputation.tar.gz 
     #4 rs2244291 53468479 0.105 1.000 1.000 2 0.935 0.991 0.960 <<well mask imputed, highly concordant

  ##plink_1.90 --threads 1 --bfile /reference/genepi/public_reference_panels/1000G_20101123_v3/derived_plinkformat/chr4.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR --ld rs2244291  rs2898681
      #--ld rs2244291 rs2898681:
      # R-sq = 0.0198602      D' = 0.804617
      # Still, pretty big coincidence. 


  zgrep -a rs2256965 ../cleaned/IMPUTATION_RESULTS/CIDR_WAMHS_IMPUTATION/Archive_CIDR_WAMHS_chromosome_6_imputation.tar.gz 
      #6 rs2256965 31555130 0.588 1.000 1.000 2 0.953 0.993 0.982 <<well imputed, hihgly concordant.
plink_1.90 --threads 1 --bfile /reference/genepi/public_reference_panels/1000G_20101123_v3/derived_plinkformat/chr6.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR --ld rs1799964 rs2256965
     #--ld rs1799964 rs2256965:    R-sq = 0.199124       D' = 0.978956  <<suggests residual from the IBD set. 

   #no further action needed on SNPs for WAMHS. Problem with using cases as controls though....

########################################################################
#11.0 HRC VCF files 
########################################################################
#never used recode to vcf but seems to work okay. As per https://imputationserver.sph.umich.edu/start.html#!pages/help also need to create a sorted *.vcf.gz file using VCFtools and tabix (including bgzip):
  #they recommend running checkVCF - downloaded from https://github.com/zhanxw/checkVCF and put 
      #Currently turned vcfCHECK off. Imputation server doesn't seem to require strand alignment, and it is disagreeing with the recommended tool by Will Raynor which also strand aligned explicitely to the HRC. Here Ref genome is hs37d5.fa, dont' know if that makes a difference - not sure if explicitely the same.
            #AMFS chr1 - tells me 18790 inconsistant sites (as in different strand). but aligned by will's function... tempted to drop this step. No other errors. 
            #AMFS 22 - 3072 inconsistant. Again, should already be aligned
   

#from HEI conversion - "Warning: Underscore(s) present in sample IDs." may need to fix?

for x in 1
do
  for i in 2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD 2016_cleaning_QMEGA_610k_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD 2016_cleaning_QMEGA_omni_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD 2016_cleaning_EPIGENE_QTWIN_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD 2016_cleaning_HEIDELBERG_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD 2016_cleaning_WAMHS_OAG_IBD_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD
  do
    module load vcftools/0.1.14
    #make vcf of all SNPs
    echo -e  "\nConverting ${i} files into sorted .vcf.gz"
      for RR in {1..23}
      do
        plink_1.90 --threads 1 --bfile ${dirK}/HRC_upload_cleaned_files_lower_MAF/${i}-updated-chr${RR} --recode vcf --out ${dirK}/HRC_upload_cleaned_files_lower_MAF/${i}-updated-chr${RR} --silent
        #gzip normally doesn't like working across directories, I think, so be safe
        cd ${dirK}/HRC_upload_cleaned_files_lower_MAF
        vcf-sort ${i}-updated-chr${RR}.vcf | bgzip -c > ${i}-updated-chr${RR}.vcf.gz
          ##python /mnt/lustre/home/matthewL/bin/checkVCF/checkVCF.py -r /mnt/lustre/home/matthewL/bin/checkVCF/hs37d5.fa -o test ${i}-updated-chr${RR}.vcf.gz
            #AMFS chr1 - tells me 18790 inconsistant sites (as in different strand). but aligned by will's function... tempted to drop this step. No other errors. 
            #AMFS 22 - 3072 inconsistant. 

        rm -f ${i}-updated-chr${RR}.vcf
        echo ""
        plink_1.90 --threads 1 --bfile ${dirK}/HRC_upload_cleaned_files/${i}_maf01-updated-chr${RR} --recode vcf --out ${dirK}/HRC_upload_cleaned_files/${i}_maf01-updated-chr${RR} --silent
        cd ${dirK}/HRC_upload_cleaned_files
        vcf-sort ${i}_maf01-updated-chr${RR}.vcf | bgzip -c > ${i}_maf01-updated-chr${RR}.vcf.gz
          ##python /mnt/lustre/home/matthewL/bin/checkVCF/checkVCF.py -r /mnt/lustre/home/matthewL/bin/checkVCF/hs37d5.fa -o test ${i}_maf01-updated-chr${RR}.vcf.gz

        rm -f ${i}_maf01-updated-chr${RR}.vcf   
        echo -e "\n${i} chr ${RR} done\n"
      done
  done #loop through sets
done #loop to turn on/off

#quick test have the new version of files

grep -w "rs9662633\|rs41306613\|rs4670855\|rs4493574\|rs9443079\|rs61731700\|rs9403839\|rs13265018\|rs10107655\|exm779484\|rs61866610\|rs74017766\|exm1345519\|rs7276462" ${dirK}/HRC_upload_cleaned_files_lower_MAF/2016_cleaning_EPIGENE_QTWIN_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr*.bim
  #not in there, good

grep -w "rs17774013\|rs870585\|rs644893" ${dirK}/HRC_upload_cleaned_files_lower_MAF/2016_cleaning_HEIDELBERG_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr*.bim
  #not in there, good

 #16/08/2016 did a test run of AMFS in a few batches (kept failing until I realised chrX not working for HRC yet. However i chose HRC r1.1 2016, (only other HRC depreciated). Phasing via EAGLE, EUR pop. QC and imputation whereas we need phasing by SHAPEIT. took only ~48hrs to run though so deleted it and started again


########################################################################
#12.0 Submitting and saving HRC files
########################################################################
#STUART POINTED OUT  - save the output to a folder with the runN. HRC r1.1 2016, SHAPEIT
#15/08/2016 - no MAF filter
  #AMFS - job-20160815-014613 - KircxIgXIO <<<<moved to /working 24/08/2016
  #EPIGENE - job-20160815-015124 - zXDzOpfILl <<moved to /working 25/08/2016
  #HEi - job-20160815-015621 - ULVnSkVXkc <<moved to working 30/08/2016
  #610k - job-20160815-020126 - kVWqFZReRc <<moved 31/08 but borked it, need to rerun 3-9
  #Omni: job-20160815-021257 - tWJNholRxV <<moved 01/09
  #WAMHS: ..... job refused - "Error Bad Request {"success":false,"message":"Only 5 jobs per user can be executed simultaneously."} "
     #well that scuppers that. Could log in remotely and upload more. 

#23/08/2016 started downloading past jobs into job ID filders. Running
  #WAMHS no MAF filter -  job-20160823-000628 - npte job 51, under heavy load << 73GrTP3ISNNhp #moved 05/09/2016


#23/08/2016 with 0.01 MAF Filter
    #AMFS- 1799 rare variants - run job-20160823-010530 - already down to 44 jobs lUjWE521jcnP9 <<moved 06/09/2016
    #610k  - 36 rare SNPs, not worth an entire imp run <<<SKIP 
    #EPIGENE/QTWIN - 22673 rare SNPs - job-20160823-011125 << douCejRNBJJlJ <<moved 07/09/2016 
    #HEI - 2754 rare SNPs - job-20160823-013344   <<  bt6GLbxANQC7O # moved 12/09/2016. 
    #QMEGA_omni - 2291 rare SNPs - run  job-20160823-013754 <<iI4KAhYtJcA1V <<moved 07/09/2016
    #WAMHS - 13893 rare SNPs <<<not enough free job spaces (5 max) - Run 30/08/2016. Think I might have uploaded it twice (first one seem to have failed but looks lke it was still uploading even though the upload window had closed)/ job-20160830-023412 <<M2y68BJoe8T3D #moved 10/09

 #01/09/2016 re-imputing the unfiltered 610k data for chrom 3-9; screwed up moving to working and deleted them to be safe. job-20160831-195838 <<1X4A4WqkCWCI2 #moved 
   #

  #Files are literally just named chr_1.zip etc. Download each to a folder named after the job name, then use this to correct the naming and move to /working
  jobname=job-20160831-195838 #job-20160823-013344 #job-20160830-023412 #job-20160823-013754 #job-20160823-011125 #job-20160823-010530 #job-20160823-000628 #job-20160815-021257 #job-20160815-020126 #job-20160815-015621 #job-20160815-015124  #job-20160815-014613 
  cd ~/Downloads/${jobname}
  new_name=2016_cleaning_QMEGA_610k_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr #2016_cleaning_HEIDELBERG_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD_maf01-updated-chr  #2016_cleaning_WAMHS_OAG_IBD_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD_maf01-updated-chr #2016_cleaning_QMEGA_omni_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD_maf01-updated-chr #2016_cleaning_EPIGENE_QTWIN_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD_maf01-updated-chr  #2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD_maf01-updated-chr #2016_cleaning_WAMHS_OAG_IBD_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr #2016_cleaning_QMEGA_omni_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr #2016_cleaning_QMEGA_610k_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr #2016_cleaning_HEIDELBERG_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr #2016_cleaning_EPIGENE_QTWIN_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr #2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr
  imp_pw=1X4A4WqkCWCI2 #bt6GLbxANQC7O #M2y68BJoe8T3D #iI4KAhYtJcA1V #douCejRNBJJlJ #lUjWE521jcnP9 #73GrTP3ISNNhp #tWJNholRxV #kVWqFZReRc #ULVnSkVXkc #zXDzOpfILl #KircxIgXIO 
  dirQQ=/working/lab_stuartma/matthewL/MELANOMA/IMPUTATION_2016
  #mv statistics_chr14_21.txt ${dirQQ}/statistics_${new_name}14_21.txt
  #mv statistics_chr1_5.txt ${dirQQ}/statistics_${new_name}1_5.txt
  #mv statistics_chr22.txt ${dirQQ}/statistics_${new_name}22.txt
  #mv statistics_chr6_13.txt ${dirQQ}/statistics_${new_name}6_13.txt
  #mv QC-Report_AMFS_chr14_21.html ${dirQQ}/QC-Report_${new_name}14_21.html
  #mv QC-Report_AMFS_chr1_5.html ${dirQQ}/QC-Report_${new_name}1_5.html
  #mv QC-Report_AMFS_chr22.html ${dirQQ}/QC-Report_${new_name}22.html
  #mv QC-Report_AMFS_chr6_13.html ${dirQQ}/QC-Report_${new_name}6_13.html
  
  mv qcreport.html ${dirQQ}/QC-Report_${new_name}.html
  mv statistics.txt ${dirQQ}/Statistics_${new_name}.txt
  #+++++++testing new version for chr22  new AMFS data in ~/Downloads/job-20160815-014613, seems to run okay for chr22, try for rest before the results are deleted from the website
  for i in {3..9} #{1..22}
  do
    mv chr_${i}.log ${new_name}${i}.log
    unzip -P ${imp_pw} chr_${i}.zip
    mv chr${i}.dose.vcf.gz.tbi ${new_name}${i}.dose.vcf.gz.tbi
    mv chr${i}.info.gz ${new_name}${i}.info.gz
    #files may be missing line break on final line based on  files may be missing line break on final line. This test http://stackoverflow.com/questions/10082204/add-a-newline-only-if-it-doesnt-exist should check for a new line and fix it if not there. AMFS chr22 seemed to be okay. Tested it on a set of dummy files and seemed to work okay
    if [ "$(zcat "chr${i}.dose.vcf.gz" | tail -c1 ;echo x)" != $'\nx' ]
      then
        echo -e "\nMissing end of line break found, fixing\n"
        zcat chr${i}.dose.vcf.gz |  sed  -e '$a\' > ${new_name}${i}.dose.vcf; gzip ${new_name}${i}.dose.vcf
        mv chr${i}.dose.vcf.gz ${new_name}${i}.dose.vcf.gz
      else
        echo -e "\nvcf has correct line breaks on each line\n"
        mv chr${i}.dose.vcf.gz ${new_name}${i}.dose.vcf.gz
      fi
    #files may be missing line break on final line - this command should fix that (didn't run for AMFS) BUT NOT TESTED YET!
      #zcat chr${i}.dose.vcf.gz |  sed  -e '$a\' > ${new_name}${i}.dose2.vcf; gzip ${new_name}${i}.dose2.vcf
      #AMFS chr22 seemed to be okay? Not sure if this command is working - testing on new AMFS chr22 data. Flagged .gz as lacking the line break even when it seems to have it
    zip Archive_${new_name}${i}.zip  ${new_name}${i}.*
    rm -f ${new_name}${i}.dose.vcf.gz ${new_name}${i}.dose.vcf.gz.tbi ${new_name}${i}.info.gz ${new_name}${i}.log chr_${i}.zip
    mv Archive_${new_name}${i}.zip ${dirQQ}/ 

  done
 

#####06/09/2016 Further runs. The QC files (*.html) suggest there is a systematic difference in allele freq between my data and the HRC (plot of my MAF vs HRC has 2 lines - most SNPS MAF = MAF but for a subset the HRC MAFS are +20%. Chatted to SG and he said this was a known problem with HRC v1.0, (~10% of the genotype calls were borked) and should have been fixed with the 1.1 release I used. However he checked a log he got last week and indeed he is seeing the same problem with his GWAS set. MQ used 1000 genomes phase 3 and is not seeing this problem, so as a parallel check as we investigate I am going to run AMFS against 1000 genomes
   #So AMFS, all SNPs, SHAPEIT, 1000 genomes phase 3 v5, EUR set.  Didn't note job, I think it is job-20160906-021717 2sjBIdMGLfF9J 
   #19/10/2016 never moved this to scratch, doing so now.

  new_name=2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD_1000G_imp-updated-chr
  imp_pw=2sjBIdMGLfF9J
  dirQQ=/working/lab_stuartma/matthewL/MELANOMA/IMPUTATION_2016

  mv qcreport.html ${dirQQ}/QC-Report_${new_name}.html
  mv statistics.txt ${dirQQ}/Statistics_${new_name}.txt
  for i in {1..21} #22
  do
    mv chr_${i}.log ${new_name}${i}.log
    unzip -P ${imp_pw} chr_${i}.zip
    mv chr${i}.dose.vcf.gz.tbi ${new_name}${i}.dose.vcf.gz.tbi
    mv chr${i}.info.gz ${new_name}${i}.info.gz
    #files may be missing line break on final line based on  files may be missing line break on final line. This test http://stackoverflow.com/questions/10082204/add-a-newline-only-if-it-doesnt-exist should check for a new line and fix it if not there. AMFS chr22 seemed to be okay. Tested it on a set of dummy files and seemed to work okay
    if [ "$(zcat "chr${i}.dose.vcf.gz" | tail -c1 ;echo x)" != $'\nx' ]
      then
        echo -e "\nMissing end of line break found, fixing\n"
        zcat chr${i}.dose.vcf.gz |  sed  -e '$a\' > ${new_name}${i}.dose.vcf; gzip ${new_name}${i}.dose.vcf
        mv chr${i}.dose.vcf.gz ${new_name}${i}.dose.vcf.gz
      else
        echo -e "\nvcf has correct line breaks on each line\n"
        mv chr${i}.dose.vcf.gz ${new_name}${i}.dose.vcf.gz
      fi
    #files may be missing line break on final line - this command should fix that (didn't run for AMFS) BUT NOT TESTED YET!
      #zcat chr${i}.dose.vcf.gz |  sed  -e '$a\' > ${new_name}${i}.dose2.vcf; gzip ${new_name}${i}.dose2.vcf
      #AMFS chr22 seemed to be okay? Not sure if this command is working - testing on new AMFS chr22 data. Flagged .gz as lacking the line break even when it seems to have it
    zip Archive_${new_name}${i}.zip  ${new_name}${i}.*
    rm -f ${new_name}${i}.dose.vcf.gz ${new_name}${i}.dose.vcf.gz.tbi ${new_name}${i}.info.gz ${new_name}${i}.log chr_${i}.zip
    mv Archive_${new_name}${i}.zip ${dirQQ}/

  done


 



#+++++still to do - check the assoc for anything super weid for HEI and WAMHS, filter out the last SNPS from EPIGENE and regenerate, then upload.
  #

########################################################################
#12.0 Initial look at AMFS files
########################################################################
#.log not that useful, just logs N of SNPs imputed
dirQQ=/working/lab_stuartma/matthewL/MELANOMA/IMPUTATION_2016
cd $dirQQ
head 2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr22.log
#[attempt_201606111211_11579_m_000000_0] INFO  Imputation for chunk chunk_22_0040000001_0060000000 successful.
#[attempt_201606111211_11579_m_000000_0] INFO    chunk_22_0040000001_0060000000 Snps in info chunk: 188754
#[attempt_201606111211_11579_m_000001_0] INFO  Imputation for chunk chunk_22_0020000001_0040000000 successful.
#[attempt_201606111211_11579_m_000001_0] INFO    chunk_22_0020000001_0040000000 Snps in info chunk: 288490
#[attempt_201606111211_11579_m_000002_0] INFO  Imputation for chunk chunk_22_0000000001_0020000000 successful.
#[attempt_201606111211_11579_m_000002_0] INFO    chunk_22_0000000001_0020000000 Snps in info chunk: 47300

#the statistics files just seems to be alleles where the ref SNP is different in my data than HRC (e.g. A/G bs G/A)

wc -l statistics_2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr6_13.txt
  #92692 statistics_2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr6_13.txt
grep -c "Allele switch" statistics_2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr6_13.txt
  #92693. Huh? Oh, files aren't formatted correctly, last line break borked.

wc -l statistics_2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr*.txt
#   48395 statistics_2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr14_21.txt
#   79791 statistics_2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr1_5.txt
#    3041 statistics_2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr22.txt
#   92692 statistics_2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr6_13.txt
  #223919

grep -c "Allele switch" statistics_2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr*.txt
#statistics_2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr14_21.txt:48396
#statistics_2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr1_5.txt:79792
#statistics_2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr22.txt:3042
#statistics_2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr6_13.txt:92693
  #looks good

#from the zip you get three files (unpacked chr22 to test
  #2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr22.info.gz
  #2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr22.dose.vcf.gz
  #2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr22.dose.vcf.gz.tbi

 wc -l 2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr22.dose.vcf.gz.tbi
  #120

 #tbi is a vcf specialised format - see https://genome.ucsc.edu/goldenpath/help/vcf.html

zcat 2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr22.info.gz | wc -l
   #524544
  #end is broken - what is wrong with their software? See above for a fix I can pipeline but this will do it

 gunzip 2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr22.info.gz
 sed -i -e '$a\' 2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr22.info
 gzip 2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr22.info

 #vcf seems to be okay
 #gunzip 2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr22.dose.vcf.gz #<might not be smart, a gig or so.
  # sed -i -e '$a\ 2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr22.dose.vcf
 #gzip 2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr22.dose.vcf

 zcat 2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr22.info.gz | grep -c "Genotyped"
   #11549
  #http://genome.sph.umich.edu/wiki/Minimac3_Info_File SNps are in chr:post dammit
    #LooRsq This statistic can only be provided for genotyped sites. This is similar to the estimated Rsq above, but the imputed dosages value used to compare are calculated by hiding all known genotypes for the given SNP (see LooDosage).
  #EmpR, EmpRsq While the LooRsq statistic completely ignores experimental genotypes, EmpR is calculated by calculating the correlation between the true genotyped values and the imputed dosages that were calculated by hiding all known genotyped for the given SNP (see LooDosage). A negative correlation between imputed and experimental genotypes can indicate allele flips. This statistic also can only be provided for genotyped sites. EmpRsq is the square of this correlation.

  #so I want to look for high Loo, say >0.8 and low EmpRsq and also look for -ve EmpR

  zcat 2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr22.info.gz | awk '$9>0.8 && $11<0.9' | wc -l 
    #414

  zcat 2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr22.info.gz | awk '$9>0.8 && $11<0.85' | wc -l
    #162

  zcat 2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr22.info.gz | awk '$9>0.8 && $11<0.8' | wc -l
    #59

  zcat 2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr22.info.gz | awk 'NR==1 || $9>0.8 && $11<0.75' 
#SNP	REF(0)	ALT(1)	ALT_Frq	MAF	AvgCall	Rsq	Genotyped	LooRsq	EmpR	EmpRsq	Dose0	Dose1
#22:17712611	G	A	0.31525	0.31525	0.99320	0.97483	Genotyped	0.812	0.865	0.74797	0.87413	0.08384
#22:19717320	G	A	0.35012	0.35012	0.99311	0.97684	Genotyped	0.817	0.859	0.73705	0.86826	0.08711
#22:25419524	C	T	0.43037	0.43037	0.99626	0.98852	Genotyped	0.830	0.845	0.71432	0.91062	0.13418
#22:25423062	G	A	0.17647	0.17647	0.99234	0.96025	Genotyped	0.827	0.808	0.65222	0.63597	0.00940 <<
#22:25425439	A	C	0.39928	0.39928	0.98995	0.96830	Genotyped	0.840	0.840	0.70627	0.86140	0.09158
#22:25428780	C	T	0.14578	0.14578	0.99540	0.97163	Genotyped	0.850	0.865	0.74881	0.68359	0.00140
#22:25959183	A	C	0.04829	0.04829	0.99727	0.95638	Genotyped	0.818	0.814	0.66188	0.56636	0.00060 <<
#22:27929853	T	C	0.93660	0.06340	0.99999	0.99977	Genotyped	1.000	0.851	0.72483	0.98252	0.00000
#22:31533225	G	A	0.56588	0.43412	0.99363	0.98122	Genotyped	0.815	0.815	0.66375	0.88868	0.15370 <<
#22:32339782	G	A	0.07264	0.07264	0.99863	0.98444	Genotyped	0.803	0.858	0.73587	0.82056	0.02125
#22:33898427	A	G	0.23596	0.23596	0.99545	0.98272	Genotyped	0.810	0.844	0.71245	0.82765	0.06140
#22:35532059	T	C	0.41592	0.41592	0.99522	0.98445	Genotyped	0.842	0.859	0.73765	0.87753	0.08928
#22:37259273	C	T	0.47055	0.47055	0.99833	0.99439	Genotyped	0.813	0.863	0.74475	0.87474	0.09743
#22:38864460	G	A	0.10982	0.10982	0.99859	0.98918	Genotyped	0.802	0.862	0.74379	0.81427	0.02869
#22:39525436	C	T	0.42475	0.42475	0.99000	0.96733	Genotyped	0.876	0.866	0.74999	0.79948	0.01298
#22:44496917	A	G	0.15443	0.15443	0.99863	0.99146	Genotyped	0.806	0.866	0.74926	0.85889	0.04664
#22:44499558	G	A	0.07799	0.07799	0.99618	0.95986	Genotyped	0.804	0.781	0.60982	0.70386	0.02133 <<
#22:44499643	T	C	0.08721	0.08721	0.99869	0.98745	Genotyped	0.802	0.860	0.74012	0.82820	0.02565
#22:45258495	A	G	0.55006	0.44994	0.99143	0.97391	Genotyped	0.927	0.866	0.74950	0.84425	0.00149
#22:49489917	A	G	0.09799	0.09799	0.99585	0.96577	Genotyped	0.842	0.855	0.73124	0.64365	0.00015 
#22:50793661	T	C	0.30208	0.30208	0.98657	0.95466	Genotyped	0.915	0.816	0.66533	0.72079	0.00461 <<
    #21, 22 inc header

#hard to be too concerned - that said so easy to filter and rerun... maybe just for new regions? A few are properly low Info > 0.8, con <0.66 or so

  zcat 2016_cleaning_AMFS_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr22.info.gz | awk 'NR==1 ||$8=="Genotyped" && $10 < 0' | wc -l
  #1 - header only.

  #Still would be an easy QC - get the whole imp, flag low concord SNP, remove, rerun. 


#05/10/2016 looking at the weird results from HRC. Log report in the zip isn't informative; it literally the list of chunks run and nothing more. The QC-Report....html however lists SNPs with extreme p values, and shows the plot of input vs ref MAF. the statistic file seems to be most be allele switches (different A1/A0 order) rather than a problem
  #looking at the QC-Report. The unusual band is a set of SNPs where the minimum MAF is minimum is ~0.20, yet MAF in my data look to be normally distributed.  
  #The table shows for most chr score or so SNPs with big freq differences, e.g. ## Mismatched frequencies for '1:180146385' f[C,T] = [0.8491,0.1509] vs [0.9814,0.01861], chisq 1449
  #For QMEGA_omni there is a huge run on chrom 15 where the the maf in the HRC set is always ~20%, so start looking at chr15. Note though the problem may be more severe than shown by the extreme chisq values - there might be more bad SNPS where the freq difference isn't as sever (is 25 when it should be MAF 20)


cd $dirQQ
mkdir temp
cd temp
cp ../QC-Report_2016_cleaning_QMEGA_omni_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr.html .
cp ../Archive_2016_cleaning_QMEGA_omni_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr16.zip .
cp ../Statistics_2016_cleaning_QMEGA_omni_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr.txt .

  ##grep "Mismatched" QC-Report_2016_cleaning_QMEGA_omni_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr.html | sed 's/## //g' | sed 's/&#39;//g' | grep -v frequencies
    #so all are Mismatched frequencies

#look at a few chr other than 15 to see if there is a systematic difference (e.g. always omni or always HRC diff to 1KG. If random probably just sampling between the sets)
k=EPIGENE_QTWIN # HEIDELBERG # QMEGA_omni
dirV=/working/lab_stuartma/scratch/matthewL/ALL_1000G_phase1integrated_v3_impute

grep "Mismatched" ${dirQQ}/QC-Report_2016_cleaning_${k}_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr.html | sed 's/## //g' | sed 's/&#39;//g' | sed 's/Mismatched frequencies for //g' | sed 's/:/ /g' | sed 's/f\[//g' | sed 's/\] = \[/ /g' | sed 's/\] vs \[/ /g' | sed 's/,/ /g' | sed 's/\]  chisq / /g' | wc -l
  #QMEGA_omni 632 EPIGENE_QTWIN=260.

for i in 15 # 1 4 7 16
do
  wc -l ${dirK}/HRC_upload_cleaned_files_lower_MAF/2016_cleaning_${k}_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr${i}.bim ; echo ""
  #QMEGA_OMNI: 1=63929; 4=46354; 7=41678, 16=24364, 15 = 23189
  #EPIGENE_QTWIN: 1=22777; 4=17447, 15 = 
  grep "Mismatched" ${dirQQ}/QC-Report_2016_cleaning_${k}_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr.html | sed 's/## //g' | sed 's/&#39;//g' | sed 's/Mismatched frequencies for //g' | sed 's/:/ /g' | sed 's/f\[//g' | sed 's/\] = \[/ /g' | sed 's/\] vs \[/ /g' | sed 's/,/ /g' | sed 's/\]  chisq / /g' | awk '$1=='''${i}''' {print $2}' | wc -l ; echo ""
  #QMEGA_OMNI: 1=18; 4=6; 7=10; 16=13; 15=333. So ~50% on chr15.
  #EPIGENE_QTWIN: 1=5; 4=1 (checked, is just 1), 7=15195, 15=160 so ~60% on chr15
  #

  #To allow a more formal comparison, merge 1KG data with this set, and add a step to filter out incorrect matches (indels/SNP overlaps). May as well insert a header with printf while in there. Checked bim freqs and input data (my data) has its freq reported first in the mismatch report
  awk 'NR==FNR {a[$2]=$0;next} $2 in a {print a[$2],$0}' <(grep "Mismatched" ${dirQQ}/QC-Report_2016_cleaning_${k}_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr.html | sed 's/## //g' | sed 's/&#39;//g' | sed 's/Mismatched frequencies for //g' | sed 's/:/ /g' | sed 's/f\[//g' | sed 's/\] = \[/ /g' | sed 's/\] vs \[/ /g' | sed 's/,/ /g' | sed 's/\]  chisq / /g' | awk '$1=='''${i}''' ') <(awk 'NR==FNR {a[$1];next} $2 in a {print $1,$2,$3,$4,(1-$8),$8}' <(grep "Mismatched" ${dirQQ}/QC-Report_2016_cleaning_${k}_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr.html | sed 's/## //g' | sed 's/&#39;//g' | sed 's/Mismatched frequencies for //g' | sed 's/:/ /g' | sed 's/f\[//g' | sed 's/\] = \[/ /g' | sed 's/\] vs \[/ /g' | sed 's/,/ /g' | sed 's/\]  chisq / /g' | awk '$1=='''${i}''' {print $2}') <(gunzip -dc ${dirV}/ALL_1000G_phase1integrated_v3_chr${i}_impute.legend.gz ) ) | awk '($3==$12 && $4==$13) || ($4==$12 && $3==$13)' | awk 'BEGIN{ printf("CHR BP A1 A2 F1_in F2_in F1_HRC F2_HRC chisq 1kg_ID BP a0 a1 Fa0 Fa1\n" )} {print$0}' ; echo ""

  #chr15 has 330 or so mismtachs, a lot more than expect. So have a slightly more complex step after to explore. Are they all testing the same allele?  yep
  
    


done

#while 15 seems to be the worst, looks for other chromosomes with proportionally more. count number of SNPs per CHR in bim file, then count number of times a chr SNP is flagged as mismatch using uniq -c and sort to arrange by chr, and merge on said chromosome, then insert a header


for k in EPIGENE_QTWIN HEIDELBERG QMEGA_omni WAMHS_OAG_IBD AMFS
do
  echo -e "\n${k}"
  awk 'NR==FNR {a[$2]=$1;next} $2 in a {print $2,$1,a[$2]}' <(cut -f1 ${dirK}/HRC_upload_cleaned_files_lower_MAF/2016_cleaning_${k}_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated.bim | awk '$1<23 && $1>0' | sort | uniq -c | sed 's/ \+/ /g' | sed 's/^ //g' | sort -n -k2,2 ) <( grep "Mismatched" ${dirQQ}/QC-Report_2016_cleaning_${k}_aligned_1KG_geno02_geno3_mind3_hweCo_hweCa_diff_PCA_IBD-updated-chr.html | sed 's/## //g' | sed 's/&#39;//g' | sed 's/Mismatched frequencies for //g' | sed 's/:/ /g' | sed 's/f\[//g' | sed 's/\] = \[/ /g' | sed 's/\] vs \[/ /g' | sed 's/,/ /g' | sed 's/\]  chisq / /g' | cut -d" " -f1 | sort | uniq -c | sed 's/ \+/ /g' | sed 's/^ //g' | sort -n -k2,2) |     awk ' BEGIN{ printf("CHR N_miss N_SNP PER\n")} {print $0,(100*($2/$3))}' > ${dirI}/temp_${k}_missN; wc -l  ${dirI}/temp_${k}_missN
   #some 20 not 23 (22 + header) as some chr have no errors. Will make plotting a bit fiddlier

done

echo "CHR N_miss N_SNP PER SET" > ${dirI}/temp_HRC_miss

for k in EPIGENE_QTWIN HEIDELBERG QMEGA_omni WAMHS_OAG_IBD AMFS
do
  for i in {1..22}
  do
    if `cut -d" " -f1 ${dirI}/temp_${k}_missN | grep -q -w ${i}`
    then
      awk '$1=='''${i}''' {print $1,$2,$3,$4,"'''${k}'''" }' ${dirI}/temp_${k}_missN  >> ${dirI}/temp_HRC_miss
    else
      echo "${i} 0 0 0 ${k}" >> ${dirI}/temp_HRC_miss
    fi
  done #fill in missing chr (no mismatches) with 0
 rm -f ${dirI}/temp_${k}_missN
done




cd $dirI
module load R/3.2.2
R --no-save <<EOF
options(bitmapType='cairo')
library(ggplot2)
dd<-read.table("temp_HRC_miss",header=T)

#work out an informative range
summary(dd\$PER)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.01880 0.03162 0.19360 0.05440 4.58300

png("temp_HRC_miss.png")
qplot(CHR,PER,data=dd,colour=SET,ylim=c(0,5),xlab="CHROMOSOME",ylab="PERCENTAGE SNPs with freq chisq diff > 300")
dev.off()

EOF

#okay the plot suggests that while chr15 is definately weird - 5% of WAMHS SNPs are diff - 14 is also enriched, and to a lesser extent are 11,12,13. Chr 6 might have a slight excess but I bet that is real MHC differencs. chr22 is enriched just for HEIDELBERG...






  #quick summary. On balance from checking chr 1, 4, 7, 16 HRC tends to be the one differnt to both 1KG and omni, though not enough for me to think there is a systematic problem. Certainly omni differs from the other two many times in the list, and there are a number cases when none of them agree. So could be safe and drop these SNPs and re impute but nothing looks to be systematically wrong, and that wont fix the cases where the freq MIGHT be wrong in HRC. 


#EPIGENE/TWIN - only one on chrom1 in all three sets, and EPI matches 1KG, HRC diff. only one on chr4, in all three, HRC and 1KG the same (0.255), EPI 0.49 so EPIQTWIN extreme

#quick look at Omni chr1 SNPs that overlap my data, 1kg and HRC Following is mismtachtch report from HRC merged with 1KG data where avaliable
  #CHR BP A1 A2 F1_in F2_in F1_HRC F2_HRC chisq 1kg_ID BP a0 a1 Fa0 Fa1
  #1 36822024 T C 0.6874 0.3126 0.8725 0.1275 632.7 rs10796883 36822024 T C 0.752 0.2480  <<1KG C = 0.248 so neither match but mine is closer
  #1 107090200 G T 0.8303 0.1697 0.9458 0.05422 510 rs17434983 107090200 G T 0.8377 0.1623 << 1KG T is 0.1623; mine is 0.1697, HE is 0.05. HRC most different
  #1 119282172 G A 0.87 0.13 0.9999 0.0001385 4191 rs17022367 119282172 G A 0.8668 0.1332 << 1KG A 0.1332; mine is 0.13, HRC is 0.0001385. HRC most different 
  #1 152609199 C A 0.9227 0.07731 0.9946 0.005387 1125 rs57146904 152609199 C A 0.9129 0.0871 << 1KG A 0.0871, mine 0.077, HRC 0.005 - HRC most different
  #1 162047745 C T 0.3002 0.6998 0.1077 0.8923 776.2 rs4233385 162047745 C T 0.0818 0.9182 << 1KG T 0.91, mine 0.7, HRC 0.89 - mine most different
  #1 168460878 C T 0.9097 0.0903 0.9808 0.01918 480 rs5003114 168460878 C T 0.8892 0.1108 << KG T 0.11; mine 0.09, HRC 0.0192 -  HRC most different
  #1 240223235 T C 0.8261 0.1739 0.9611 0.03888 885.1 rs12091593 240223235 T C 0.8166 0.1834 << 1KG C 0.18, mine 0.17, HRC 0.04.  -  HRC most different
  #1 248738991 G A 0.8744 0.1256 0.9909 0.009126 1786 rs2802074 248738991 G A 0.9459 0.0541 << 1KG A  0.0541, mine is 0.13, HRC 0.009 - neither match
  #1 249148769 A T 0.9352 0.0648 0.9995 0.0005079 1881 rs12136211 249148769 A T 0.9472 0.0528 <<  1KG T 0.053, mine 0.065, HRC 0.0005 -  HRC most different



#chr4
  #CHR BP A1 A2 F1_in F2_in F1_HRC F2_HRC chisq 1kg_ID BP a0 a1 Fa0 Fa1
  #4 1291246 C T 0.8457 0.1543 0.9809 0.01907 1479 rs6823741 1291246 C T 0.8615 0.1385 << 1KG T 0.1385; omni 0.1543; HRC 0.02 < HRC most differen
  #4 10081603 T C 0.06814 0.9319 0.004602 0.9954 1013 rs4697916 10081603 T C 0.0106 0.9894 << 1KG C 0.9894, omni 0.931, HRC 0.99 <omni most different
  #4 99262877 A T 0.9341 0.06591 0.9962 0.003755 1080 rs10003517 99262877 A T 0.9763 0.0237 <<1KG T 0.0237, omni 0.066, HTC 0.003 - both different, omni closer
  #4 123183343 A G 0.8131 0.1869 0.9207 0.07931 325.6 rs79830901 123183343 A G 0.7968 0.2032 << 1KG G 0.20322, omni 0.187, HRC 0.079, HRC most different 
  #4 153995742 G A 0.9322 0.06778 0.9965 0.00354 1166 rs74909745 153995742 G A 0.9842 0.0158 << 1KG A 0.0158, omni 0.068, HRC 0.0035, both differnt, HRC closer.

#chr7
  #CHR BP A1 A2 F1_in F2_in F1_HRC F2_HRC chisq 1kg_ID BP a0 a1 Fa0 Fa1
  #7 6384653 G A 0.9531 0.0469 0.9934 0.006649 386.8 rs10235108 6384653 G A 0.9987 0.0013 <<1000 genomes ~ 0.13%, my data 4.6%, HRC 0.6%; my data most diff
  #7 12963865 T C 0.1883 0.8117 0.3724 0.6276 327.6 rs246831 12963865 T C 0.4208 0.5792 <<1000 genomes 42%, HRC 37%, my data 19%, my data most diff
  #7 27224963 T C 0.9262 0.07376 0.9964 0.003601 1313 rs17437495 27224963 T C 0.9538 0.0462 <<my data 7.4%, 1000 genomes 4.6%, HRC 0.3% ; both diff 
  #7 27889391 C T 0.7638 0.2362 0.9136 0.08642 572.3 rs10951185 27889391 C T 0.7836 0.2164 << My data 23.6%, 1KG genome 21.6, HRC 8.6 ; HRC most diff
  #7 35054045 G T 0.9174 0.08263 0.7296 0.2704 407.6 rs2893491 35054045 G T 0.719 0.2810 << N data T 0.08, 1kG 0.28, HRC 0.27 < my data most diff

#chr16
  #CHR BP A1 A2 F1_in F2_in F1_HRC F2_HRC chisq 1kg_ID BP a0 a1 Fa0 Fa1
  #16 1002648 T G 0.8135 0.1865 0.9244 0.07557 359.2 rs9933730 1002648 T G 0.7744 0.2256 <<1KG G 0.2256, omni 0.1865, HRC 0.075, HRC different
  #16 3613207 G A 0.9169 0.08312 0.9789 0.0211 344.3 rs9940099 3613207 G A 0.9921 0.0079 <<1KG A 0.0079, omni 0.083, HRC 0.02, HRC different
  #16 9325775 G T 0.8163 0.1837 0.9321 0.06787 426.1 rs13337680 9325775 G T 0.814 0.1860 <<1KG T 0.186, omni 0.183, HRC 0.06, HRC different
  #16 10893887 T C 0.9381 0.06192 0.9987 0.001293 1510 rs9930801 10893887 T C 0.9446 0.0554 <<1KG C 0.0554, omni 0.062, HRC 0.001, HRC different
  #16 16687563 G A 0.9306 0.0694 0.9997 0.0003232 2110 rs8051688 16687563 G A 1 0.0000 << 1KG A 0, omni 0.07, HRC 0.0003, Omni diff
  #16 17473650 A C 0.8649 0.1351 0.9997 0.0002924 4285 rs79614102 17473650 A C 0.8456 0.1544 << 1KG C 0.154, omni 0.135, HREC 0.00029, HRC different 
  #16 51096971 G A 0.769 0.231 0.9347 0.06533 862.5 rs8054651 51096971 G A 0.8971 0.1029 << 1KG A 0.1029, omni 0.23, HRC 0.065, Both differnt, omni most extreme
  #16 68017356 A G 0.8482 0.1518 0.9719 0.02813 965.2 rs255054 68017356 A G 0.8483 0.1517 <<1KG G 0.1517, omni 0.1518, HRC 0.028, HRC different 
  #16 77972790 C T 0.5255 0.4745 0.7039 0.2961 331.8 rs3886360 77972790 C T 0.558 0.442 0 <<1KG T 0.4420, omni 0.475, HRC 0.30 - HRC most differen
  #16 81190601 T C 0.02926 0.9707 0.002078 0.9979 419.2 rs1453325 81190601 T C 0 1.0000 <<1KG C 100%; omni 0.97,HRC 0.9979, mine different
 







 
########################################################################
#X.0 Clean up
########################################################################
${dirK}/Melanoma_clean_2016_merged_PCA.eigenvec_IDs




echo -e "\nfinito\n"
exit

########################################################################
#A.0 Appendix
########################################################################






