# Radiation Protection Traditional Chinese Medicine Component Screening System

## 1. System Overview

This system is a radiation protection traditional Chinese medicine (TCM) component screening platform developed based on the Flask framework. It aims to help researchers quickly screen TCM components and formulas with radiation protection effects. The system integrates multiple functional modules, including compound screening, network analysis, KEGG GO analysis, etc., providing strong tool support for radiation protection research.

### 1.1 Technical Architecture

- **Backend**: Flask framework
- **Frontend**: HTML, CSS, JavaScript, Bootstrap
- **API Integration**: BATMAN2 API (target prediction), STRING API (protein-protein interaction network construction)

### 1.2 System Directory Structure

```
website/
├── Compound screening results/    # Screening results storage directory
├── templates/                      # Frontend templates
│   ├── web_database_files/         # Static resource files
│   ├── Compound_screening.html     # Compound screening page
│   ├── index.html                  # Home page
│   ├── Herb_miRNA_data.html        # TCM miRNA data page
│   ├── Network_Analysis.html       # Network analysis page
│   ├── KEGG_GO.html                # KEGG GO analysis page
│   └── AI_Analysis.html            # AI analysis page
├── KEGG_GO_Analysis.py             # KEGG GO analysis module
├── Compound_screening_normal.py    # Compound screening in normal mode
├── Network_Analysis.py             # Network analysis module
├── AI_Analysis.py                  # AI analysis module
├── app.py                          # Main application
└── Test-file Zisu.csv              # Test case file
```

## 2. Functional Modules

### 2.1 Radiation Protection Compound Screening (RadioProtect Compound Screening)

**Function Description**:
- Upload TCM component CSV files to screen compounds with radiation protection effects
- Support normal mode screening
- Generate matched and unmatched result files
- Provide result download functionality

**Technical Implementation**:
- Use `Compound_screening_normal.py` to handle screening
- Perform data processing and matching through pandas
- Call BATMAN2 API for target prediction

### 2.2 Network Analysis

**Function Description**:
- Upload CSV results from 2.1 Radiation Protection Compound Screening to analyze protein-protein interaction networks
- Perform protein-protein interaction network analysis, screening core targets and their corresponding radiation protection compounds
- AI mining potential protein-protein interaction relationships and core targets, screening corresponding radiation protection compounds, generating professional analysis reports
- Provide API configuration options, supporting multiple AI models

**Technical Implementation**:
- Use `Network_Analysis.py` to implement analysis logic
- Construct protein-protein interaction networks through STRING API
- Screen core targets and their corresponding radiation protection compounds

### 2.3 KEGG GO Analysis

**Function Description**:
- Upload CSV results from 2.1 Radiation Protection Compound Screening for KEGG GO analysis

**Technical Implementation**:
- Use `KEGG_GO_Analysis.py` to implement

### 2.4 TCM miRNA Data (Herb miRNA data)

**Function Description**:
- Display TCM miRNA related data

**Technical Implementation**:
- Display miRNA data through frontend pages

### 2.5 AI Analysis of Screening Results

**Function Description**:
- Intelligently analyze compound screening results and generate professional analysis reports
- Support analysis of different cell/tissue types
- Provide API configuration options, supporting multiple AI models
- Allow custom system prompts to adjust analysis direction

**Technical Implementation**:
- Use `AI_Analysis.py` to implement backend analysis logic
- Provide user interaction interface through frontend page `AI_Analysis.html`
- Integrate OpenAI SDK to call AI models for analysis
- Read screening result files from `Compound screening results` directory
- Group by medicinal materials, count main genes, pathways, compound types, etc.
- Generate structured analysis reports, including compound type enrichment, main gene pathways, drug likeness and oral bioavailability, protection mechanisms, etc.

## 3. Operation Guide

### 3.1 System Startup

1. Ensure that the required dependency packages are installed:
   ```
   pip install flask pandas requests openai numpy matplotlib seaborn gseapy networkx
   ```

   **Note**: If using the AI analysis function, you need to install the openai library additionally.

2. Enter the website directory and run the application:
   ```
   python app.py
   ```

### 3.2 Radiation Protection Compound Screening Operation Steps

1. Click "RadioProtect Compound Screening" in the function area
2. Click the "Select File" button to upload a CSV file containing TCM components, or enter a medicinal material name in the search box and click the Search button
3. Click the "Run screening" button to submit
4. Wait for the system to complete processing
5. View the screening results, including matched and unmatched compounds
6. Click the "Download" button to download the result files (or view in the "Compound screening results" folder)

### 3.3 KEGG GO Analysis Operation Steps

1. Upload CSV results from 2.1 Radiation Protection Compound Screening
2. Select cell or tissue type
3. Select KEGG GO analysis type and species
4. View analysis results

### 3.4 Network Analysis Operation Steps

1. Upload CSV results from 2.1 Radiation Protection Compound Screening
2. Select cell or tissue type
3. Select network analysis type
4. Select AI model (optional)
5. Click the "Run network analysis" button to submit

### 3.5 TCM miRNA Data Operation Steps

1. Click "Plant microRNA Targets in Human" in the function area
2. Enter the external link page for miRNA target prediction

### 3.6 AI Analysis of Screening Results Operation Steps

1. Upload CSV results from 2.1 Radiation Protection Compound Screening
2. Select cell or tissue type
3. Configure API settings:
   - Enter API Key
   - Set API Base URL (default is https://api.deepseek.com)
   - Select AI model
4. Modify the system prompt as needed to adjust the analysis direction
5. Click the "Run AI Analysis" button to start analysis
6. Wait for the system to complete processing and view the AI-generated analysis report

## 4. Data Format Requirements

### 4.1 Radiation Protection Compound Screening Data Format

CSV files should contain the following columns (any one):
- CID columns: `cid`, `CID`, `Cid`, `cID`, `pubchem_id`, `PubChem_id`, `PubChem_ID`
- Gene name columns: `gene_name`, `Gene_name`

**Example file**: `Test-file Zisu.csv`

## 5. Result Interpretation

### 5.1 Radiation Protection Compound Screening Results

- **Matching results**: Contain compounds that match radiation protection targets
- **Unmapped results**: Unmatched compounds

Result files are saved in the `Compound screening results` directory, named in the format: `Matching normal results of [filename] ingredients.csv`

### 5.2 KEGG GO Analysis Results

- **KEGG analysis results**: Contain KEGG analysis results of radiation protection compounds, including pathway names, pathway IDs, enrichment p-values, etc.
- **GO analysis results**: Contain GO analysis results of radiation protection compounds, including gene names, gene IDs, GO IDs, GO categories, GO category p-values, etc.

### 5.3 Network Analysis Results

- **Network analysis results**: Contain network analysis results of radiation protection compounds, including protein-protein interaction networks, core genes, core gene-TCM compound relationships, etc.
- **AI model results**: Contain AI model analysis results, AI mining potential protein-protein interaction relationships, scoring and screening core genes, core gene-TCM compound relationships, etc.

## 6. System Configuration and Maintenance

### 6.1 Configuration Files

- **app.py**: Main configuration file, containing routes and core functions
- **BATMAN2_API_URL**: BATMAN2 API address (default is 'http://batman2api.cloudna.cn/queryTarget')

### 6.2 Data Files

- **DEGs_normal.csv**: Differentially expressed genes in normal mode
- **herb_ingredient.csv**: TCM component data
- **results of target matching normal.csv**: Target matching results in normal mode

## 7. Common Issues and Solutions

### 7.1 Upload File Failure

- **Issue**: The system prompts an error after uploading a file
- **Solution**: Check if the file format is correct and ensure the file contains the required columns (such as CID)

### 7.2 Empty Screening Results

- **Issue**: No matched compounds after screening
- **Solution**: Check if the CID or gene names in the input file are correct and ensure the data format meets the requirements

### 7.3 Slow System Response

- **Issue**: The system responds slowly when processing large files
- **Solution**: Reduce file size or process data in batches

### 7.4 API Call Failure

- **Issue**: BATMAN2 API call failure
- **Solution**: Check network connection and ensure the API address is correct

## 8. Example Usage Flow

### 8.1 Using Test File for Radiation Protection Compound Screening

1. After starting the system, enter the "RadioProtect Compound Screening" page
2. Select "normal" mode
3. Upload "Test-file Zisu.csv"
4. Click the "Submit" button
5. View the screening results and download "Matching normal results of Test-file Zisu ingredients.csv"

## 9. Technical Support

- **Contact**: Professor Jiahua Ma Email: jiahuama@swust.edu.cn
- **Developers**:
            Jiahua Ma (PhD, Professor, School of Life Science and Agronomy, Southwest University of Science and Technology);
            Shilong Qi (Master's student, School of Life Science and Agronomy, Southwest University of Science and Technology);
            Yang Zhou (Undergraduate student, School of Life Science and Agronomy, Southwest University of Science and Technology).
- **Institution**: College of Life Sciences and Agri-forestry, Southwest University of Science and Technology
- **Project**: Radiation Protection TCM Component Screening System

## 10. Update Log

- **2026/4/13**: Added AI analysis of screening results function, supporting intelligent analysis of compound screening results and generating professional reports
- **2026/3/30**: Modified Gene name-based mapping logic to adapt to BATMAN API data
- **2026/3/16**: Added radiation protection formula screening function
- **2026/3/14**: Added radiation protection TCM compound screening function
