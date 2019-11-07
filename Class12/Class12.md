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

We want to trim the PDB file so that we produce a protein-only PDB file and a drug (ligand) only PDB file. The `trim.pdb` function or `atom.select` function allows you to extract protein-only and ligand-only objects. By using `write.pdb`, these files will automatically be downloaded onto your R studio workspace.

``` r
protein <- atom.select(hiv, "protein", value = TRUE)
write.pdb(protein, file="1hsg_protein.pdb")
```

``` r
ligand <- atom.select(hiv, "ligand", value = TRUE)
write.pdb(ligand, file="1hsg_ligand.pdb")
```

Docking algorithms require each atom to have a charge and an atom type that describes its properties. However, typical PDB structures don’t contain this information. We therefore have to ‘prep’ the protein and ligand files to include these values along with their atomic coordinates. All this will be done in a tool called AutoDock Tools (ADT for short).

Can you locate the binding site visually? Note that crystal structures normally lack hydrogen atoms, why?

-   We are able to see a gap in the protein only structure, which suggests that the ligand binds within this gap. Crystal structures lack hydrogen atoms because visualization would become too difficult if each hydrogen atom was displayed.

Look at the charges. Does it make sense (e.g. based on your knowledge of the physiochemical properties of amino acids)?

-   The charges for each amino acid make sense. Charge values differ depending on the atom in the amino acid.

Docking ligands into HIV-1 protease
===================================

For this section, we will use the program called Autodock Vina \[4\]. Autodock Vina is a fast docking program that requires minimal user intervention and is often employed for highthroughput virtual screening. We will run it from the command line.

In order to visualize the docks and compare to the crystal conformation of the ligand we will process the `all.pdbqt` to a PDB format file that can be loaded into VMD.

``` r
res <- read.pdb("all.pdbqt", multi=TRUE)
write.pdb(res, "results.pdb")
```

Now load both the original 1hsg.pdb file and your new results.pdb file into VMD.

Qualitatively, how good are the docks? Is the crystal binding mode reproduced? Is it the best conformation according to AutoDock Vina?

-   The docks are well-fitted inside the binding site of the HIV-1 protease. The crystal binding mode is reproduced with few rotational errors. This is the best conformation according to AutoDock Vina.

To assess the results quantitatively we will calculate the RMSD (root mean square distance) between each of the docking results and the known crystal structure using the bio3d package.

``` r
ori <- read.pdb("ligand.pdbqt")
rmsd(ori, res)
```

    ##  [1]  0.590 11.163 10.531  4.364 11.040  3.682  5.741  3.864  5.442 10.920
    ## [11]  4.318  6.249 11.084  8.929

Quantitatively how good are the docks? Is the crystal binding mode reproduced within 1Å RMSD for all atoms?

-   The first dock has an RMSD of 0.590.

How would you determine the RMSD for heavy atoms only (i.e. non hydrogen atoms)?

\*Using `atom.select()` with the selection string "noh" will choose non hydrogen atoms.
