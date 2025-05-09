# Code for Inferring Cancer Mutation Causal Networks

## File Structure

### requirements.txt

This file contains the list of packages required to run the code.

### main.py

This file contiains the main code for the project. The command format to run it is as follows:

```
python main.py StudyName --max-genes --num-edges-drawn
```

- *StudyName* corresponds to the type of cancer genes you want to investigate. The valid types of strings for this variable are *BRCA, COADREAD, BLCA*

- *--max-genes* corresponds to the integer number of genes to be included in the analysis 

- *--num-edges-drawn* corresponds to the integer number of edges to be drawn for the final plot for clarity. The strongest *--num-edges-drawn* edges are the ones that are drawn.

On running the code the results would be stored in a new folder as follows:
```
results/causal_network.png, results/top_edges.csv, results/adjacency_matrix.csv
```

They correspond to the resulting image plotting the discovered causal network, listing the top edge waits and the adjacency matrix corresponding to the edge weights between different nodes.

The default set of genes chosen in the code corresponds to Breast Cancer and the Blind 40 set as discussed in the report. These include the genes:

```
genes = [
    "EGFL7", "GATA6", "ASXL2", "RALGDS", "ZNF24", "H3C13", "PDGFB", "HOOK3",
    "SFRP2", "PTPRC", "GEN1", "ZRSR2", "STAG2", "IGH", "MRE11", "H3C2",
    "HOXA13", "STAT4", "CREB3L1", "SOX17", "GPS2", "MAFB", "PALB2", "MN1",
    "MLF1", "DIS3", "WDR90", "FOLR1", "SESN3", "DNM2", "CLTCL1", "DIS3L2",
    "AFF1", "PTPRO", "ELF3", "TCEA1", "ATF1", "CRLF2", "CTDNEP1", "MLH1"
]
```

This code can be easily modified for other studies via minimal changes.

Authors: Noa Kalfus, Jeffrey Wu, and Shreyas Havaldar
