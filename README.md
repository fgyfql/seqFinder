# seqFinder
****************
If you have 
    a better than ~3.5 Å resolution density map,
    a well built C-alpha model
    and protein sequence cantidates for the C-alpha model, 
then you could try this program for finding the target sequence of the C-alpha model.

If this program is useful in your work, pleas cite our paper: XXXXX.

If you have questions when using this program, please contact Qianglin Fang via fang163@purdue.edu (Email).




****************
1. Prepare input files:

1) a file containing information on protein sequence candidates. (an example: inputfile1.txt)

2) a file containing information on the observed side-chain size distribution of the C-alpha model (an example: inputfile2.txt)

i) open the C-alpha model and the density map in the program COOT;

ii) go through each residue of the C-alpha model from its N terminus to its C terminus:

a) When encountering a residue that has prominent side-chain density (Such residues normally locate in buried regions of the structure.):
    try to mutate the residue to all the 20 different kinds of amino acids in COOT to esitimate the side-chain size;
    write one line to inputfile2.txt to record the residue number and the residue name that have similar size of the estimated size of this residue;

b) When encoutering a residue that has reliable C-alpha atom position but does not have good side-chain density:
    write one line to inputfile2.txt to record its residue number and residue name (use "X");

c) When encountering a residue that has continuous main-chain density but does not show reliable C-alpha atom position:
    do nothing and move to the next residue;

d) When encountering a very poor region whose density does not allow C-alpha atom placement or is missing:
    use information before or after this region to search the target sequence.
    
    
    

****************
2. run the program

An example command: seqFinder.py --annotation_file inputfile2.txt --fasta_file inputfile1.txt --output_file out.txt

After running the program, you will get a file named "out.txt" which contains the alignment results and the fitting scores between the target protein sequence and the side-chain size distribution pattern derived from the density map. The value of the fitting scores is between 0.0 and 1.0. The higher the fitting score, the better the sequence fits to the density map.




****************
3. result validation

Local chemical environment, disulfide bond formation, secondary structure prediction, amino acid modification (such as glycosylation) and mass spectrometry experiments can be used to validate the resultant model.
