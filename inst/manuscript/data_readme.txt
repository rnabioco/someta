Numeric data underlying Figure 2 and analyses of Puntambekar et al. 2021

fig2a.tsv:
id : GEO accession
year : year of publication
usable : single cell dataset has metadata
usable2 : format the metadata stored in

fig2b.tsv:
id : GEO accession
journal	: published journal
journal_group : subgrouping of journal

fig2c.tsv:
id : GEO accession
usable : single cell dataset has metadata
journal	: published journal
cite : number of citations
year : year of publication
ifs : journal impact factor

fig2d.tsv:
id : GEO accession
usable : single cell dataset has metadata
usable2 : format the metadata stored in
software_author : whether an author is involved with scRNA-seq software development

geo_vs_Svensson.tsv:
id : GEO accession
useable : single cell dataset has metadata. Entries in Svensson but not found in automated search are labeled NA

manual_check.tsv:
id : GEO accession
is_sc : manual check if entry is truly single cell data
called_usable : orginal automated call if single cell dataset has metadata
correct_call : manual check if automated call was correct
has_type : manual check if dataset metadata actually contains cell type column/info
