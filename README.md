# TCR Analysis with TCRdist3
TCR sequence analysis of scRNA sequencing from T cells isolated from PBMCs and tissue.

Script to comput a pairwise amino acid sequence similarity of TCR β and TCR α sequences from tissue and PBMC samples. This was achieved by computing weighted multi-CDR distances between TCRs using tcrdist3, a Python3 package for TCR repertoire analysis, following the procedure described in Mayer-Blackwell et al., 2021 (Mayer-Blackwell et al., 2021).
The aim of this code is to determine if an specific TCR β and TCR α sequence has similar TCRs in diferent datasets from tissue and PBMCs. The results were plotted in a network to visually if there were group of TCRs connected. The TCR distances were also plotted with ggplot as histograms to compare two different group of TCR and see if they presented a different distribution.
