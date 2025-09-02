# Operational Modes

STRmie-HD provides two main CLI modes: **Complete_Pipeline** and **Index_Calculation**.

## Complete Pipeline
This mode executes the **entire workflow** starting directly from raw sequencing data (`.fastq.gz` or `.fasta.gz`).  
It automatically performs:
- Histogram-based **CAG allele peak calling**.  
- **CCG repeat assignment** for each allele.  
- Calculation of **Instability and Expansion indices (II, EI)**.  
- Detection of **Loss of Interruption (LOI)** and **Duplication of Interruption (DOI)** events.  
- Generation of both an **Excel summary report** and an **interactive HTML report** for visual inspection.  

Optional graphical outputs (CAG/CCG histograms) and raw read-level counts are also produced for each sample if requested.



## Index Calculation
This mode is intended for situations where **automatic allele calling requires manual adjustment**.  
When inspecting the interactive HTML report, the user may decide that the automatically identified allele peaks are not accurate.  
Through the HTML interface, it is possible to:
- Visually explore **CAG and CCG histograms** for each sample.  
- **Manually adjust allele peak values** if necessary.  
- Export a curated Excel table with the corrected allele definitions.  

The `Index_Calculation` mode takes this adjusted table as input and **recomputes all instability and expansion indices (II, EI)** accordingly.  
This ensures that downstream results are based on manually validated allele assignments.

