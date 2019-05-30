---
author: Akshara Balachandra, A11933675
title: Find a Gene Project (BGGN 213)
date: May 15, 2019
---

# Question 1: Name of the protein

**Name**: Sonic Hedgehog protein isoform 1 preprotein (SHH)

**Accession**: NP_000184.1

**Species**: Homo sapiens

# Question 2: BLAST search

**Method**: TBLASTN search against Xenopus ESTs

**Database**: Expressed sequence tags (est)

**Organism**: Xenopus (Taxid: 8353)

## Alignment details

**Chosen match:** Accession BX704643.1, a 874bp clone from Xenopus *tropicalis*.

```
>BX704643.1 BX704643 XGC-tadpole Xenopus tropicalis cDNA
clone TTpA011g13 5', mRNA sequence

cDNA clone TTpA011g13, mRNA sequence.
Length: 874

Score           Expect   Method                         Identities
387 bits(994)   3e-132   Compositional matrix adjust.   197/285(69%)

Positives      Gaps         Frame
220/285(77%)   23/285(8%)   +2


Query  68   RNSERFKELTPNYNPDIIFKDEENTGADRLMTQRCKDKLNALAISVMNQWPGVKLRVTEG  127
            RNS+RFKELTPNYNPDIIFKDEENTGADRLMTQRCKDKLNALAISVMNQWPGVKLRVTEG
Sbjct  41   RNSDRFKELTPNYNPDIIFKDEENTGADRLMTQRCKDKLNALAISVMNQWPGVKLRVTEG  220

Query  128  WDEDGHHSEESLHYEGRAVDITTSDRDRSKYGMLARLAVEAGFDWVYYESKAHIHCSVKA  187
            WDEDGHHSEESLHYEGRAVDITTSDRDRSKYGMLARLAVEAGFDWVY+ESKAHIHCSVKA
Sbjct  221  WDEDGHHSEESLHYEGRAVDITTSDRDRSKYGMLARLAVEAGFDWVYFESKAHIHCSVKA  400

Query  188  ENSVAAKSGGCFPGSATVHLEQGGTKLVKDLSPGDRVLAADDQGRLLYSDFLTFLDRDDG  247
            ENSVAAKSGGCFPGSA V +E GGTK V++L PGDRVL++D QG L+YSDFL F+D++
Sbjct  401  ENSVAAKSGGCFPGSARVMVEPGGTKAVRELRPGDRVLSSDPQGNLIYSDFLLFIDKEHD  580

Query  248  AKKVFYVIETREPRERLLLTAAHLLFVAPHNDSATXXXXXXXXXXXXXXXXXXXRALFAS  307
            KK++YVI+T + R R  +TAAHLLFVA  N + +                   +++FAS
Sbjct  581  VKKLYYVIQTSQTRIR--MTAAHLLFVAQSNGTGS------------------FKSVFAS  700

Query  308  RVRPGQRVYVVAERDGDRRLLPAAVHSVTLSEEAAGAYAPLTAQG  352
            VRPG  +Y    R  D  L  A V  V L EE  G + P    G
Sbjct  701  NVRPGDVIYSADRR--DMTLREAMVEKVDL-EEDIGGFCPRNCPG  826
```

![BLAST parameters](img/tblastn_params.png){ width=5in }

![BLAST search results](img/q2_results.png){ width=5in }

![Alignment details](img/q2_alignment.png){ width=4.5in }

\pagebreak

# Question 3: Novel Protein Information

**Chosen Sequence:**
```
>Xenopus tropicalis protein (from BLAST result)
RNSDRFKELTPNYNPDIIFKDEENTGADRLMTQRCKDKLNALAISVMNQWPGVKLRVTEG
WDEDGHHSEESLHYEGRAVDITTSDRDRSKYGMLARLAVEAGFDWVYFESKAHIHCSVKA
ENSVAAKSGGCFPGSARVMVEPGGTKAVRELRPGDRVLSSDPQGNLIYSDFLLFIDKEHD
VKKLYYVIQTSQTRIRMTAAHLLFVAQSNGTGSFKSVFASNVRPGDVIYSADRRDMTLRE
AMVEKVDLEEDIGGFCPRNCPG
```

**Name:** *Xenopus* sonic hedgehog

**Species:** Xenopus (Silurana) tropicalis (Common name: western clawed frog, Taxid: 8364)

Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;
Amphibia; Batrachia; Anura; Pipoidea; Pipidae; Xenopodinae; Xenopus; Silurana

# Question 4: Proof that gene/protein are novel

A BLASTP search against the NR database yielded a top hit to a protein from
*Xenopus tropicalis* with 100% coverage and 97.71% identity. The alignment
details are given in the screenshots below.

![BLASTP parameters](img/q4_blastp.png){width=6in}

![BLASTP results page](img/q4_results.png){width=6in}

![Alignment for top result](img/q4_alignment.png){width=6in}
