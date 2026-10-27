The enclosed Matlab-scripts are intented to guide the determination of retention times and selection of mass transitions of a combinatorial library. Mass data in uncompressed mzML-files from samples containing mixtures of products are used as input. Our mass-data has been recorded on a Agilent QTOF-system and the scripts have not been tested on data recorded on other systems.

If, as in our case, mass transition data for all combinations of a library of 30 azides (named Z1-Z30) and 25 alkynes (named A1-A25) should be established, 30 + 25 samples should to be prepared. In the first 30 samples each azide is combined with a mixture of all 25 alkynes. These samples are named Z1A, Z2A, ..., Z30A. In the next 25 samples each alkyne is combined with a mixture of all 30 azides. These samples are named A1Z, A2Z, ..., A25Z.

To execute the scrips, we have used Matlab R2024a.

How to use the script:

1. Convert your mass data to uncompressed mzML-files using e.g. Proteowizards MSconvert. The data and filenames should follow the pattern described above.

2. Specifiy the exact mass of all builiding blocks in the file "Builidingblock masses.xlsx".

3. Open the script "findRetetionTimes.m"

4. Optionally, set the variable "massTolerance" (delta m/z value between theoretical and recored MS1-mass)

5. Specify the number of different collision energies that have been recorded by setting the variable "noCollisionEnergies".

6. Specify the mass-adduct that should be searched for by setting the variable "massAdduct". (eg. 1.0073 for positive mode ionization and -1.0073 for negative mode ionization)

7. Give the file path to where the recorded mzML-files is located by setting the variable "mzMLpath"

8. Run the script.


This give two output files: "RTAnalysis.pdf" and "RT.txt". Manually inspect the file RTAnalysis.pdf that containes extracted mass traces from the recorded data - one compound per page. Since each product should be formed twice, e.g. product Z5A3 should be present in both the file A3Z.mzML and Z5A.mzML. Each page in the pdf-file contains two overlayed mass traces - one from each file. Suggested retention times are printed in the file RT.txt.

9. Manually curate the retention times and update the file "RT.txt" before proceeding with the next step.


Finding mass transitions:

10. Optionally specify unwanted fragment losses in the file "Excluded fragment losses.xlsx" (e.g. 18 for water)

11. Run the script "findMassTransitions.m"

This gives as output the file "MSMSAnalysis.pdf", which gives the fragmentation data for each product at the recorded collision energies, one panel for each of the up to two times each product has been recorded and found.

The script also outputs a file "missingMSMSData.txt" where all products without any found MSMS-data are listed.

The extracted MSMS-data is found in the files "MSMSData.txt" and "MSMSDataCompr.txt"

In the file "Transitions.txt" four suggested mass transitions for each recorded product is given.

12. All transitions should be manually curated.

