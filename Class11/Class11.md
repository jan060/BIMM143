Structural Bioinformatics
================
Julie Nguyen
November 5, 2019

Introduction to the RCSB Protein Data Bank (PDB)
================================================

Download CSV file from PDB website and anaylze the PDB statistics.

``` r
data <- read.csv("Data Export Summary.csv")
```

Determine the percentage of structures solved by X-Ray and Electron Microscopy.

``` r
percent_structure <- round((data$Total/sum(data$Total)) * 100, 2)
percent_structure
```

    ## [1] 89.07  8.14  2.50  0.19  0.10

Also can you determine what proportion of structures are protein?

``` r
proportion_protein <- round(sum(data$Proteins)/sum(data$Total) * 100, 2)
proportion_protein
```

    ## [1] 92.71

Visualizing the HIV-1 protease structure
========================================

In this section we will visualize the X-ray crystal structure of HIV-1 protease with a bound drug molecule indinavir (PDB ID: 1HSG). We will use the VMD molecular viewer to visually inspect the protein, the binding site, and the drug molecule.

Water molecules normally have 3 atoms. Why do we see just one atom per water molecule in this structure?

-   In VMD, atoms can be colored by molecule instead of individual atoms.

There is a conserved water molecule in the binding site. Can you identify this water molecule? What residue number does this water molecule have (see note below)?

-   Residue number 308.

As you have hopefully observed HIV protease is a homodimer (i.e. it is composed of two identical chains). With the aid of the graphic display and the sequence viewer extension can you identify secondary structure elements that are likely to only form in the dimer rather than the monomer?

-   Many turns are associated at the center of the molecule, where the two subunits come togeher.

Introduction to Bio3D in R
==========================

Bio3D is an R package for structural bioinformatics. Features include the ability to read, write, and analyze biomolecular structure, sequence, and dynamic trajectory data.

Protein Data Bank files (or PDB files) are the most common format for the distribution and storage of high-resolution biomolecular coordinate data.

Read and inspect the on-line file with PDB ID 1HSG. How many amino acid residues are there in this pdb object and what are the two non-protein residues?

``` r
library(bio3d)
protease <- read.pdb("1HSG")
```

    ##   Note: Accessing on-line PDB file

-   There are 198 residues and the two non-protein residues are water and MK1, the ligand.

What type of R object is `protease$atom`?

-   This outputs a dataframe.

Atom selection
==============

The Bio3D `atom.select()` function is arguably one of the most challenging for newcomers to master. It is however central to PDB structure manipulation and analysis. At its most basic, this function operates on PDB structure objects (as created by `read.pdb()`) and returns the numeric indices of a selected atom subset. These indices can then be used to access the `$atom` and `$xyz` attributes of PDB structure related objects.

Use the Bio3D write.pdb() function to write out a protein only PDB file for viewing in VMD. Also write out a second separate PDB file for the ligand with residue name MK1.

### Protein only

``` r
protein_only <- atom.select(protease, "protein", value = TRUE)
write.pdb(protein_only, file = "1hsg_protein.pdb")
trim.pdb(protease, inds = atom.select(protease, "protein"))
```

    ## 
    ##  Call:  trim.pdb(pdb = protease, inds = atom.select(protease, "protein"))
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1514,  XYZs#: 4542  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 0  (residues: 0)
    ##      Non-protein/nucleic resid values: [ none ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, helix, sheet, seqres, xyz,
    ##         calpha, call

### Ligand

``` r
ligand <- atom.select(protease, "ligand", value = TRUE)
write.pdb(ligand, file = "1hsg_ligand.pdb")
trim.pdb(protease, inds = atom.select(protease, "ligand"))
```

    ## 
    ##  Call:  trim.pdb(pdb = protease, inds = atom.select(protease, "ligand"))
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 45,  XYZs#: 135  Chains#: 1  (values: B)
    ## 
    ##      Protein Atoms#: 0  (residues/Calpha atoms#: 0)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 45  (residues: 1)
    ##      Non-protein/nucleic resid values: [ MK1 (1) ]
    ## 
    ## + attr: atom, helix, sheet, seqres, xyz,
    ##         calpha, call

3D structure viewing in R
=========================

Working with multiple PDB files
===============================
