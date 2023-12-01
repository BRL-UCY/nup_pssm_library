## Profiles

Here we provide the PSSM and HMM profiles that were used in this study to perform identification of nucleoporin sequences in proteomes.

1) [PSSM profiles](./nup-pssm/)\
These profiles were created using PSI-BLAST and our library of collected nucleoporin sequences.
2) [PFAM HMM profiles](./pfam-hmm/)\
These profiles were collected from PFAM, where available, for domains identifying some of the nucleopori families. Nomeclature follows the PFAM ID, refer to [mapping](./pfam-hmm/profile_mapping.csv) file to match the profile to the nucleoporin family.
3) [Nucleoporin HMM profiles](./nup-hmm/)\
These profiles were created using our library of collected nucleoporin sequences and the **hmmbuild**[^1] software.



[^1]: Eddy,S.R. (2011) Accelerated Profile HMM Searches. PLoS Comput. Biol., 7, e1002195.
