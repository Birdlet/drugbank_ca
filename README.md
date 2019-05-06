# drugbank_ca
Drugbank.ca xml file converted to rational tables for SQL database

This is a fork of Michael Yourshaw's work
> drugbankxml2db 1.0β1 ©2014 Michael Yourshaw all rights reserved


Writes relational database text files for drugs, drug_target, drug_target_action, and targets tables.

  Input: the path of a DrugBank xml file.
  Output: text files that can be used as inputs to SQL tables.


## Requirements
  * lxml
  * codecs
  
  
## Method
  1. Read and parse xml file.
  2. Extract data and save records as key, value pairs.
  3. Write output files.
  
  
## Usage:
```bash
drug_bank_xml2db.py [-h] --input INPUT [--output OUTPUT]

writes relational database text files for drugs, drug_target,
drug_target_action, and targetsdrug_enzyme, drug_enzyme_action and enzyme

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT, -i INPUT
                        drug bank xml file downloaded from http://www.drugbank
                        .ca/system/downloads/current/drugbank.xml.zip
  --output OUTPUT, -o OUTPUT
                        output directory for write relation tables, default is
                        input directory
```
