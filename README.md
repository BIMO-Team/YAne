**YAne**: A high-**Y**ielding primer **A**ssessment a**n**d database construction tool for **e**DNA metabarcoding

# Introduction
YAne is a python package for eDNA metabarcoding primer evaluation and database construction.
This package is equiped with functions to perform data retrieval from NCBI (mock community), simulate *in-silico* PCR based on primer pairs, filter sequence quality, and evaluate primer performance on specificity, length, reference database availability and taxonomic resolution. Simultaneously, primer-specific databases are constructed based on primer-extracted datasets in a format of Naive-bayes classifiers to be used for taxonomic assignment with QIIME2.

# Package installation
**1. Download package**: https://anonymous.4open.science/r/YAne-1406

**2. Install package dependencies and create environment**

We used a Conda environment configuration file (`environment.yml`) to manage all software dependencies. In the package's main directory, execute the following command to download all dependencies and create a new conda environment called `yane-2026.1`

```bash
conda env create -f environment.yml
```

**3. Activate conda environment**

```bash
conda activate yane-2026.1
```

# Overview of the workflow & functions
## Workflow

![Overview of YAne workflow](yane_workflow_final.png)

# Prepare input tables

User-defined input tables are required to provide information on **mock community** and **primer pairs** for *in-silico* PCR, and **target taxa** for primer evaluation.

As demonstration, we used vertebrates as our mock community, and eight primer pairs to conduct *in-silico* PCR against the mock community. Whales, dolphins and dugong are our target taxa we aimed to detect with eDNA metabarcoding.

| Data | Table file path (demo) |
| :--- | :--- |
| **Mock community Entrez** | `mock_community.txt` |
| **Primer sequences** | `primer.txt` |
| **Target taxa labels** | `cetacean.txt` |

## 1. Mock community information
* To cunstruct a **Mock community table**, a header is assigned as **mock label** and **NCBI Entrez** as the first and second colums of the first row, and contents are provided from the second row onwards (See **Table 1**)

* [NCBI Entrez](https://www.ncbi.nlm.nih.gov/books/NBK3837/) is used to indicate which mock community and genetic compartments will be downloaded from the NCBI database. 

* Large download (more than 100 requests) need to be done on weekends or between 9 pm and 5 am Eastern Time weekdays according to [NCBI Policy](https://www.ncbi.nlm.nih.gov/home/about/policies/). 

* Many mock communities can be downloaded at the same time by adding more rows to the table and specify their labels. *Primer table* will be used to indicate the mock community labels that primers will be extracted from.

* The table must be saved **tab-delimited** file format (See: `mock_community.txt`), which are used as an input for the functions `downloadmock`, `processprimer` and `evaluateprimer`. 

**Table 1** A mock community table used in our study as a demonstration. We used vertebrate (*txid7742*) mitochondrial DNA of all lengths and regions (*mitochondrion[filter]*) as our mock community and filtered out the records that are from environmental samples and are unclassified ones (*NOT uncultured[Title] NOT unclassified[Title] NOT unidentified[Title] NOT unverified[Title]*).

| mock label | NCBI Entrez |
| :--- | :--- |
| mtDNA | txid7742[Organism] AND mitochondrion[filter] NOT "environmental samples"[Title] NOT "environment"[Title] NOT uncultured[Title] NOT unclassified[Title] NOT unidentified[Title] NOT unverified[Title] |

## 2. Primer information
* A Primer table is consisted of:
    * **primer labels**
    * **mock community labels** that the primer will be extracted from (must be similar to that of a *Mock community table*), 
    * **forward primer sequence** (5' to 3')
    * **reverse primer sequence** (5' to 3') 
   
    The table must have header and content, and the columns must be arranged as shown in **Table 2**.
* If primers conduct in-silico PCR against different mock community datasets, labels of the mock can be indicated in the 2nd column.
* This table must be in **tab-delimited** file format (e.g. `primer.txt`) and will be used as an input for the functions `processprimer` and `evaluateprimer`.

*In this study, we used 8 candidate primers to do in-silico PCR against a single vertebrate mock community called mtDNA.*

**Table 2** A Primer table used to provide information on primer pairs to do in-silico PCR against a single vertebrate mock community called mtDNA
    
| primer | mock | forward | reverse |
| :--- | :--- | :--- | :--- |
| MiFishU | mtDNA | GTCGGTAAAACTCGTGCCAGC | CATAGTGGGGTATCTAATCCCAGTTTG |
| MarVer1 | mtDNA | CGTGCCAGCCACCGCG | GGGTATCTAATCCYAGTTTG |
| V12SU | mtDNA | GTGCCAGCNRCCGCGGTYANAC | ATAGTRGGGTATCTAATCCYAGT |
| Vert01 | mtDNA | TTAGATACCCCACTATGC | TAGAACAGGCTCCTCTAG |
| Mamm01 | mtDNA | CCGCCCGTCACCCTCCT | GTAYRCTTACCWTGTTACGAC |
| MarVer3 | mtDNA | AGACGAGAAGACCCTRTG | GGATTGCGCTGTTATCCC |
| V16SU | mtDNA | ACGAGAAGACCCYRYGRARCTT | TCTHRRANAGGATTGCGCTGTTA |
| dlp1.5dlp4 | mtDNA | TCACCCAAAGCTGRARTTCTA | GCGGGWTRYTGRTTTCACG |


## 3. Target taxa labels
* A target taxa table is used to specify targets that you aim to detect with candidate primers based on eDNA metabarcoding. 
* This table consisted of a **header** indicating **taxonomic level** (Kingdom, Phylum, Class, Order, Family, Genus, or Species) that the taxa names belong to as the first row. The table **contents** indicate target taxa names based on [NCBI Taxonomy Database](https://www.ncbi.nlm.nih.gov/taxonomy) , which can be searched through [NCBI Taxonomy Browser](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/tree/).

**Table 3** A target taxa table of thr demo 

| Family | 
| :--- |
| Balaenopteridae |
| Delphinidae |
| Phocoenidae |
| Monodontidae |
| Ziphiidae |
| Balaenidae |
| Platanistidae |
| Trichechidae |
| Physeteridae |
| Lipotidae |
| Eschrichtiidae |
| Iniidae |
| Pontoporiidae |
| Dugongidae |
| Neobalaenidae |

# Download mock community
The function `downloadmock` downloads and dereplicates mock community sequences and taxonomic information from NCBI. 

## Command-Line Options
* For description and usage of the function and its flags, execute `yane downloadmock --help`.
* Outputs from this step are downloaded sequence and taxonomic information from NCBI database in QIIME2's artifact file format (`sequence.qza` and `taxonomy.qza`) located in the automatically generated directory called `NCBI_data` 

| Category | Flag | Description | Remarks |
| :--- | :--- | :--- | :--- |
| **Inputs** | `--i-mock` | Mock community table file path | Required | 
|| `--p-taxa` | Mock community Taxa ID | Required |
|| `--p-dl-period` | Data downloading period | Required |
| **Parameters** | `--p-jobs` | Number of concurrent download connection | Integer; Default = 5 |
|| `--p-threads` | Number of computation threads | Integer; Default = 4 |

## Usage Examples
* Mock community information table file path: `mock_community.txt`
* Mock community Taxa ID: `7742` (Vertebrates)
* NCBI data downloading period: `1.26` (January 2026)

### Example #1: Default jobs and threads

```bash
yane downloadmock \
--i-mock mock_community_info.txt \
--m-taxa 7742 \
--m-dl-period 1.26
```

### Example #2: Customized jobs = 3 and threads = 6

```bash
yane downloadmock \
--i-mock mock_community_info.txt \
--m-taxa 7742 \
--m-dl-period 1.26 \
--p-jobs 3 \
--p-threads 6
```

# In-silico PCR and process extracted primer datasets
The function `processprimer` conducts in-silico PCR against the mock community and cleans the datasets. 

## Command-Line Options
* For description and usage of the function and its flags, execute `yane processprimer --help`.

| Category | Flag | Description | Remarks |
| :--- | :--- | :--- | :--- |
| **Inputs** | `--i-mock` | Mock community table file path | Required | 
|| `--i-primer` | Primer table file path | Required |
| **Output** | `--o-dir` | Output directory name | Required |
| **Parameters** | `--p-jobs` | Number of concurrent download connection | Integer; Default = 5 |
|| `--p-threads` | Number of computation threads | Integer; Default = 4 |
|| `--p-identity` | Primer matching idenitity* | float (0-1); Default = 0.8  |
|| `--p-homopolymer` | Homopolymer length filtering criteria* | Integer; Default = 8 |
|| `--p-degenerate` | Degenerate base length filtering criteria* | Integer;  Default = 5 |
|| `--p-label` | Primer names of customized parameters (identity/homopolymer/degenerate base) for each primer dataset. Must be similar to the labels in primer table and in the same order of the customized parameters | Required if customized parameters are set |
|| `--p-observe-length` | Observe length before customization | yes/no; Default = no |
|| `--p-length-stat-level` | If `--p-observe-length` is yes. Then the length statistics will report length grouped by taxonomic level as additional information, this can be customized to any taxonomic level. | Taxonomic levels: Kingdom, Phylum, Class, Order, Family, Genus, Species; Default = Class  |

**can be set globally for all datasets or customized to each dataset*

## Usage Examples
* Mock community information table file path: `mock_community.txt`
* Primer information table: `primer.txt`
* Output directory: `data_prep`

### Example #1: Default parameters
* Default parameters; jobs = 5, threads = 4, identity = 0.8, homopolymer = 8, degenerate base = 5, sequence length = 1st and 99th percentiles

```bash
yane processprimer \
--i-mock mock_community.txt \
--i-primer primer.txt \
--o-dir data_prep
```

### Example #2: Customized identity, homopolymer, and/or degenerate base for all primer-extracted datasets
* Customized identity = 0.7 for all primer datasets

```bash
yane processprimer \
--i-mock mock_community_info.txt \
--i-primer primer_info.txt \
--o-dir data_prep_iden07 \
--p-identity 0.7
```

* Customized homopolymer length = 10 bp for all primer datasets

```bash
yane processprimer \
--i-mock mock_community_info.txt \
--i-primer primer_info.txt \
--o-dir data_prep_hmplm10 \
--p-homopolymer 10
```

* Customized degenerate base length = 1 bp for all primer datasets

```bash
yane processprimer \
--i-mock mock_community_info.txt \
--i-primer primer_info.txt \
--o-dir data_prep_degen1 \
--p-degenerate 1
```

### Example #3: Customized identity, homopolymer, and/or degenerate base specific to primer-extracted datasets
*Homopolymer length customization is shown as an example, for identity and degenerate base, they can be applied in the same manner*

| Dataset | Homopolymer length |
| :--- | :---: |
| MiFishU | 7 |
| MarVer1 | 7 |
| V12SU | 7 |
| Vert01 | 7 |
| Mamm01 | 7 |
| MarVer3 | 8 |
| V16SU | 8 |
| dlp.15dlp4 | 8 |

* Customized homopolymer length specific to primer-extracted datasets

```bash
yane processprimer \
--i-mock mock_community_info.txt \
--i-primer primer_info.txt \
--o-dir data_prep_varyhmplm \
--p-homopolymer 7 7 7 7 7 8 8 8 \
--p-label MiFishU MarVer1 V12SU Vert01 Mamm01 MarVer3 V16SU dlp1.5dlp4
```

### Example #4: Observe sequence length before filtering to prepare for customized length filtering criteria

```bash
yane processprimer \
--i-mock mock_community_info.txt \
--i-primer primer_info.txt \
--o-dir data_prep_customlen \
--p-observe-length yes \
--p-length-stat-level Order
```

# Customize sequence length filtering criteria
## Command-Line Options

| Category | Flag | Description | Remarks |
| :--- | :--- | :--- | :--- |
| **Input** | `--p-label` | Primer label for custom length filtering | Required |
| **Parameters** | `--p-min-length` | Minimum length for all taxa filtering | Required* |
| | `--p-max-length` | Maximum length for all taxa filtering | Required* |
| | `--p-fillen-by-taxon` | Length filtering criteria file path for taxa-specific filtering | Requires* |

**User chooses only one way of length filtering, 1.) all taxa or 2.) taxa-specific. If you would like to filter length of all taxa, `--p-min-length` and `--p-max-length` are required. If specific-taxa length filtering criteria are used, `--p-fillen-by-taxon` is required.*

## Usage Examples

### Example #1: Customized length for all taxa

Example of desired length filtering criteria

| Dataset | Min | Max |
| :--- | :--- | :--- | 
| MiFishU | 163 | 183 |
| MarVer1 | 154 | 176 |
| V12SU | 150 | 172 |
| Vert01 | 82 | 102 |
| Mamm01 | 50 | 63 |
| MarVer3 | 150 | 236 |
| V16SU | 148 | 234 |
| dlp1.5dlp4 | 144 | 679 |

**Example of command** 
```bash
yane customlength \
--p-label MiFishU MarVer1 V12SU Vert01 Mamm01 MarVer3 V16SU dlp1.5dlp4 \
--p-min-length 163 154 150 82 50 150 148 144 \
--p-max-length 183 176 172 102 63 236 234 679
```

### Example #2: Customized length specific to taxa

Example of filter length by taxon criteria table of MiFishU dataset: `fillen_MiFishU.txt`

| Taxa | Min | Max |
| :--- | :--- | :--- | 
| Actinopteri | 167 | 178 |
| Amphibia | 158 | 179 |
| Aves | 175 | 189 
| Chondrichthyes | 180 | 185 |
| Cladistia | 171 | 172 |
| Hyperoartia | 176 | 182 |
| Lepidosauria | 159 | 178 |
| Mammalia | 164 | 174 |
| Myxini | 172 | 172 |
| Sarcopterygii | 166 | 184 |

**Example of command** 
```bash
yane customlength \
--p-label MiFishU MarVer1 V12SU Vert01 Mamm01 MarVer3 V16SU dlp1.5dlp4 \
--p-fillen-by-taxon fillen_MiFishU.txt fillen_MarVer1.txt fillen_V12SU.txt fillen_Vert01.txt fillen_Mamm01.txt fillen_MarVer3.txt fillen_V16SU.txt fillen_dlp1.5dlp4.txt
```

# Primer evaluation
The function `evaluateprimer` analyzes the primer-extracted datasets and generates comparative tables and graphical representations for primer evaluation in four aspects.
* Primer specificity
* Reference database availability and intraspecific genetic variation
* Amplicon length
* Taxonomic resolution

## Command-Line options

| Category | Flag | Description | Remarks |
| :--- | :--- | :--- | :--- |
| **Input** | `--i-mock` | Mock community table file path | Required |
| | `--i-primer` | Primer table file path | Required |
| | `--i-dir` | Preprocessed primer dataset directory path | Required |
| | `--m-target-taxa` | Target taxa table file path |  Required |
| **Output** | `--o-dir` | Evaluated primer datasets directory name |  Required |
| **Parameters** | `--p-tax-level` | Taxonomic level for grouping in primer evaluation results | Taxonomic levels: Kingdom, Phylum, Class, Order, Family, Genus, Species; Default = Class |
| | `--p-vary-confidence` | Confidence values to be varied for target taxa classification | 0-1; Default = 0, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 1 |

## Usage Examples

### Example #1: Default parameters

```bash
yane evaluateprimer \
--i-mock mock_community.txt \
--i-primer primer.txt \
--i-dir data_prep \
--i-target-taxa cetacean.txt \
--o-dir primer_eval
```

### Example #2: Customize taxonomic level

```bash
yane evaluateprimer \
--i-mock mock_info.txt \
--i-primer primer_info.txt \
--i-dir processprimer_ds \
--i-target-taxa cetacean.txt \
--o-dir primer_eval \
--p-tax-level Phylum
```

### Example #3: Customize confident threshold

```bash
yane evaluateprimer \
--i-mock mock_info.txt \
--i-primer primer_info.txt \
--i-dir processprimer_ds \
--i-target-taxa cetacean.txt \
--o-dir primer_eval \
--p-vary-confidence 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1
```
