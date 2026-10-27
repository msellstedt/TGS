The script "integrateQQQmzML.m" integrates and plots QQQ-MSMS data. It allows for plotting of primary and secondary transitions for each compound. The MSMS-data should be stored in *uncompressed* mzML-format, recorded as centroid data.

The script searches for mass-transitions that are specified in the file "CompoundTransitions.xlsx".

The file "Plotgroups.xlsx", specifies the file names of the mzML-file that will be plotted, as well as plotting styles, and data-labels. Also information about compound groups are stored in the "Plotgroups.xlsx"-file. - For improved analytical purposes, datasets with a                                                                                                                                                                                                                           large number of mass transitions can be divided into several groups allowing detection of a subset of compounds in each datafile.

How to use the script:

1. Specify MSMS-transitions in the file "CompoundTransitions.xlsx"
2. Specify plotting information and compound grouping in the file "Plotgroups.xlsx"
3. Run the script "integrateQQQmzML.m"

A semicolon-separated output txt-file that contains the integrals for all found compounds is created. A pdf-file containing the recorded chromatograms and integration windows for up to three peaks is also generated. This graphical data should be be used to manually verify that correct integrations have been made.

