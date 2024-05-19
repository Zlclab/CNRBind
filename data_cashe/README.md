# RNA Binding Site Prediction

## Datasets

TR60 and TE18 were constructed in RNAsite [1]. By applying pairwise structure similarity clustering with a TM-scoreRNA cutoff of 0.3 [2], 634 redundant structures were removed. To further avoid redundancy between the training and test sets, sequence-based similarity clustering was used [3]. This process resulted in 78 remaining RNAs, which were divided into 57 clusters with 0.3 sequence similarity. The training set (denoted as TR60) comprises 42 diverse clusters with 60 RNAs, while the test set consists of 15 clusters with 18 RNAs (denoted as TE18).

The RB9 dataset is adapted from the RB19 dataset in the RBind method [4]. RBind initially retained 19 RNAs by filtering out RNA structures with multistrand or pseudoknot interactions. However, due to different data collection processes, there were 10 duplicate RNAs between TR60 and RB19. Consequently, RLBind [5] removed these 10 duplicates to re-evaluate model performance, resulting in the RB9 dataset. In this study, the dataset TR60 is utilized as the training dataset, the others are utilized to test the performance of our model.

## References

[1] H. Su, Z. Peng, J. Yang, "Recognition of small molecule-RNA binding sites using RNA sequence and structure," *Bioinformatics*, 2021; 37(1): 36-42. [https://doi.org/10.1093/bioinformatics/btaa1092](https://doi.org/10.1093/bioinformatics/btaa1092)

[2] S. Gong, C. Zhang, Y. Zhang, "RNA-align: quick and accurate alignment of RNA 3D structures based on size-independent TM-scoreRNA," *Bioinformatics*, 2019; 35(21): 4459-4461. [https://doi.org/10.1093/bioinformatics/btz282](https://doi.org/10.1093/bioinformatics/btz282)

[3] W. Li, A. Godzik, "Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences," *Bioinformatics*, 2006; 22(13): 1658-1659. [https://doi.org/10.1093/bioinformatics/btl158](https://doi.org/10.1093/bioinformatics/btl158)

[4] K. Wang, Y. Jian, H. Wang, C. Zeng, Y. Zhao, "RBind: computational network method to predict RNA binding sites," *Bioinformatics*, 2018; 34(18): 3131-3136. [https://doi.org/10.1093/bioinformatics/bty345](https://doi.org/10.1093/bioinformatics/bty345)

[5] K. Wang, R. Zhou, Y. Wu, M. Li, "RLBind: a deep learning method to predict RNA-ligand binding sites," *Brief Bioinform*, 2023; 24(1): bbac486. [https://doi.org/10.1093/bib/bbac486](https://doi.org/10.1093/bib/bbac486)
