# chem_scripts

This contains various scripts to make pseudo-potential generation and use easier. It works with Turbomole.

## Placing pseudo-potentials

Use the coord.py script to place pseudo-potentials on a standard geometry. Keywords are:

 - `guess` This attempts to place sp2 and sp3 potentials throughout the **entire** molecule. In order to use it you will need:
   - A Turbomole geometry file, named `coord`, with **no internal coordinates**
   - No control file (`guess` will make you a control file and assign the correct ecps and basis sets by itself)
   
The guess function also takes the arguments `incl` and `excl`. Using `incl` followed by carbon indices means the guess function will guess the potentialisation for those atoms *only*. Using `excl` followed by carbon indices means the guess function will guess the potentialisation for all atoms *except* those specified. Example usage: `coord.py guess incl 1-6,23` or `coord.py guess excl 30-49,1,6`.
   
   The reference basis file currently uses a combination of potentials optimised for sp2 1e, sp2 2e and sp3 1e Carbon atoms. The dummy atoms are a mixture of Helium and Lithium (this makes it easier to distinguish them from any remaining all-electron Hydrogens. The Carbon basis used is def-SV(P) with the s-functions removed. Example usage: `coord.py guess`.

 - `sp2 int1,int2,int3...` This places the 6 dummy atoms around the specified atoms in the sp2 pattern. The integers represent the atom index (starting from 1). Example: `coord.py sp2 2-4,5`.
 
 - `sp3 in1,int2,int3...` This places 3 dummy atoms around the specified atoms in the sp3 (methyl-group) pattern. The integers represent the atom index (starting from 1). Example: `coord.py sp3 2-4,5`.
 
 - `del int1,int2,int3...` This removes the specified atoms and is chained with the sp2 or sp3 options. Use it to remove the all-electron Hydrogen atoms. For example, `coord.py sp3 1 del 2-4` will add dummy atoms in the methyl pattern around atom 1, and removes the atoms 2, 3 and 4. This is useful because the atom indices often change during pseudo-potentialisation.  
 
 ## Optimising potentials
 
 ```
#######################################
The Multi-Orbital Optimiser  (___)
                             <(o o)>______
         (or MOO)              ../ ` # # \`;   
                                 \ ,___, /
                                  ||   ||  
                                  ^^   ^^  
#######################################
 ```
 
 This is done with the MOO script. Run `multivariate_optimisation.py` to have MOO guide you through the setup. This is a work in progress, a better menu system is on the way...
 
 The `opt.moo` settings file governs the optimisation parameters.
