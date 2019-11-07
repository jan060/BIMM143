Bioinformatics in Drug Discovery and Design
================

In-silico docking of drugs to HIV-1 protease
============================================

We will continue with the structure of the HIV-1 protease and its inhibitor, indinavir (PDB ID: 1HSG), and computationally dock drug molecules into the binding site of HIV-1 protease to see how well bioinformatics-based docking can reproduce the crystallographically observed binding pose.

``` r
library(bio3d)
protease <- get.pdb("1hsg")
```

    ## Warning in get.pdb("1hsg"): ./1hsg.pdb exists. Skipping download

``` r
hiv <- read.pdb(protease)
hiv
```

    ## 
    ##  Call:  read.pdb(file = protease)
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 172  (residues: 128)
    ##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, xyz, seqres, helix, sheet,
    ##         calpha, remark, call

What is the name of the two non protein `resid` values in this structure? What does `resid` correspond to and how would you get a listing of all `resid` values in this structure?

-   The two non-protein resid values are HOH and MK1. `resid` corresponds to "residue id" and you can access all `resid` with `pdb$atom$resid`.

We want to trim the PDB file so that we produce a protein-only PDB file and a drug (ligand) only PDB file.

``` r
protein <- atom.select(hiv, "protein", value = TRUE)
write.pdb(protein, file="1hsg_protein.pdb")
```

``` r
ligand <- atom.select(hiv, "ligand", value = TRUE)
write.pdb(ligand, file="1hsg_ligand.pdb")
```
