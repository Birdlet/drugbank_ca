#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys
import os
import argparse
import codecs
from collections import Counter
from lxml import etree
# from xml import etree #untested

# -i /scratch1/tmp/myourshaw/drugbank/drugbank.xml
# -i /Users/myourshaw/Downloads/drugbank.xml


def run(input, output):
    """Writes relational database text files for drugs, drug_target, drug_target_action, and targets tables.
    Input: the path of a DrugBank xml file.
    Output: text files that can be used as inputs to SQL tables.
    Method:
        1. Read and parse xml file.
        2. Extract data and save records as key, value pairs.
        3. Write output files."""
    # output file names
    if output == "default":
        output = os.path.dirname(input)

    drugs_out = os.path.join(output, 'drugs.txt')
    drug_target_out = os.path.join(output, 'drug_target.txt')
    drug_target_action_out = os.path.join(output, 'drug_target_action.txt')
    targets_out = os.path.join(output, 'targets.txt')
    drug_enzyme_out = os.path.join(output, 'drug_enzyme.txt')
    drug_enzyme_action_out = os.path.join(output, 'drug_enzyme_action.txt')
    enzymes_out = os.path.join(output, 'enzymes.txt')

    # counter for number of records in each file
    record_counts = Counter()

    # open input file and parse xml
    # get drugbank namspace
    print('Reading and parsing xml file.')
    tree = etree.ElementTree(file=input)
    # a few namespace tricks to make the code more readable
    ns = tree.getroot().nsmap
    ns['db'] = ns[None]
    del ns[None]

    drugs = tree.xpath('db:drug', namespaces=ns)

    # define tables
    # key and value name tuples for each table
    drugs_key_names = ('drug_id',)
    drugs_value_names = ['drug_name', 'groups', 'drug_description', 'pathways', 'general_references',
                         'type', 'formula', 'inchi', 'smiles', 'chembl']

    drug_target_key_names = ('drug_id', 'target_id', )
    drug_target_value_names = ['drug_target_references', ]
    drug_target_action_key_names = ('drug_id', 'target_id', 'action', )
    drug_target_action_value_names = []
    targets_key_names = ('target_id', )
    targets_value_names = ['target_name', 'UniProtKB']

    drug_enzyme_key_names = ('drug_id', 'target_id', )
    drug_enzyme_value_names = ['drug_enzyme_references', ]
    drug_enzyme_action_key_names = ('drug_id', 'target_id', 'action', )
    drug_enzyme_action_value_names = []
    enzymes_key_names = ('enzyme_id', )
    enzymes_value_names = ['enzyme_name', 'UniProtKB']

    # dicts for each table
    # in the form {<key tuple>: <value dict>}
    drugs_records = {}

    drug_target_records = {}
    drug_target_action_records = {}
    targets_records = {}

    drug_enzyme_records = {}
    drug_enzyme_action_records = {}
    enzymes_records = {}
    # drug_transporter_records = {}
    # drug_transporter_action_records = {}
    # transporters_records = {}
    # drug_carrier_records = {}
    # drug_carrier_action_records = {}
    # carriers_records = {}

    # process drug records
    print('Processing records.')
    for drug in drugs:
        # initialize dict to save unique records to print at end of run
        # in the form {(<key tuple>): {<value dict>}]}
        drugs_record = {}
        drug_id = [i for i in drug.xpath('db:drugbank-id', namespaces=ns) if i.attrib.get('primary') == 'true'][0].text
        drugs_record['drug_name'] = drug.xpath('db:name', namespaces=ns)[0].text
        # Basic Property for Drugs
        # description
        drug_description = drug.xpath('db:description', namespaces=ns)[0].text
        # deal with drugs that have no description or stray newlines, linefeeds, and tabs
        drugs_record['drug_description'] = drug_description.strip().replace('\n', '').replace('\r', '')\
            .replace('\t', ' ') if drug_description else ''
        # restrictions/groups
        drugs_record['groups'] = ','.join([g.text for g in drug.xpath('db:groups/db:group', namespaces=ns)])
        # pathways
        pathway = drug.xpath('db:pathways/db:pathway/db:name', namespaces=ns)
        drugs_record['pathways'] = ','.join([p.text for p in pathway]) if pathway else ''
        # general_references
        general_references = drug.xpath('db:general-references', namespaces=ns)[0].text
        drugs_record['general_references'] = general_references.strip().replace('\n', '').replace('\r', '')\
            .replace('\t', ' ') if general_references else ''
        # drug type
        drugs_record['type'] = drug.attrib["type"]
        # Calculated properties:
        # Formula, InChi, \InChiKey, SMILES, \MW
        drugs_record["formula"] = ''
        drugs_record["inchi"] = ''
        drugs_record["smiles"] = ''
        # drugs_record["molecular_weight"] = ''
        # drugs_record["inchikey"] = ''
        if drug.xpath("db:calculated-properties", namespaces=ns):
            flag = ""
            for property in drug.xpath("db:calculated-properties/*/*/text()", namespaces=ns):
                if flag == "formula":
                    drugs_record["formula"] = property
                    flag = ""
                elif property == "Molecular Formula":
                    flag = "formula"
                elif flag == "smiles":
                    drugs_record["smiles"] = property
                    flag = ""
                elif property == "SMILES":
                    flag = "smiles"
                elif flag == "inchi":
                    drugs_record["inchi"] = property
                    flag = ""
                elif property == "InChI":
                    flag = "inchi"
        # ChEmBL ID
        external_links = drug.xpath('db:external-identifiers/*/*/text()', namespaces=ns)
        drugs_record["chembl"] = ''
        flag = ""
        if external_links:
            for link in external_links:
                if link == "ChEMBL":
                    flag = "chembl"
                elif flag == "chembl":
                    drugs_record["chembl"] = link
                    flag = ""
        # ADMET properties
        # metabolism
        drugs_record['metabolism'] = drug.xpath("db:metabolism/text()", namespaces=ns)[0] \
            if drug.xpath("db:metabolism/text()", namespaces=ns) else ''
        # clearance
        drugs_record['clearance'] = drug.xpath('db:clearance/text()', namespaces=ns)[0] \
            if drug.xpath('db:clearance/text()', namespaces=ns) else ''
        # half-life
        drugs_record['half_life'] = '. '.join(drug.xpath('db:half-life/text()', namespaces=ns)) \
            if drug.xpath('db:half-life/text()', namespaces=ns) else ''
        # absorption
        drugs_record['absorption'] = '. '.join(drug.xpath('db:absorption/text()', namespaces=ns)) \
            if drug.xpath('db:half-life/text()', namespaces=ns) else ''
        # add this record to the table
        drugs_records[(drug_id,)] = drugs_record
        record_counts['drugs'] += 1

        # process the targets of each drug
        for target in drug.xpath('db:targets/db:target', namespaces=ns):
            # target reference
            targets_record = {}
            drug_target_record = {}
            target_id = target.xpath('db:id', namespaces=ns)[0].text
            targets_record['target_name'] = target.xpath('db:name', namespaces=ns)[0].text
            targets_record['UniProtKB'] = ""
            uniprotid_flag = 0
            for extdb in target.xpath("db:polypeptide/db:external-identifiers/db:external-identifier", namespaces=ns):
                for tdbid in extdb.xpath("db:*/text()", namespaces=ns):
                    if uniprotid_flag:
                        targets_record['UniProtKB'] = tdbid
                        uniprotid_flag = 0
                        break
                    if tdbid == 'UniProtKB':
                        uniprotid_flag = 1
            targets_records[(target_id,)] = targets_record
            record_counts['targets'] += 1

            # drug_target_references
            drug_target_references = target.xpath('db:references', namespaces=ns)[0].text
            drug_target_record['drug_target_references'] = drug_target_references.strip().replace('\n', '')\
                .replace('\r', '').replace('\t', ' ') if drug_target_references else ''
            drug_target_records[(drug_id, target_id,)] = drug_target_record
            record_counts['drug_target'] += 1

            # process the actions of the drug on the target
            for action in drug.xpath('db:targets/db:target/db:actions/db:action', namespaces=ns):
                drug_target_action_record = {}
                action = action.text
                drug_target_action_records[(drug_id, target_id, action,)] = drug_target_action_record
                record_counts['drug_target_action'] += 1

        # process the enzymes of each drug
        for enzyme in drug.xpath('db:enzymes/db:enzyme', namespaces=ns):
            # enzyme reference
            enzymes_record = {}
            drug_enzyme_record = {}
            enzyme_id = enzyme.xpath('db:id', namespaces=ns)[0].text
            enzymes_record['enzyme_name'] = enzyme.xpath('db:name', namespaces=ns)[0].text
            enzymes_record['UniProtKB'] = ""
            uniprotid_flag = 0
            for extdb in enzyme.xpath("db:polypeptide/db:external-identifiers/db:external-identifier", namespaces=ns):
                for tdbid in extdb.xpath("db:*/text()", namespaces=ns):
                    if uniprotid_flag:
                        enzymes_record['UniProtKB'] = tdbid
                        uniprotid_flag = 0
                        break
                    if tdbid == 'UniProtKB':
                        uniprotid_flag = 1
            enzymes_records[(enzyme_id,)] = enzymes_record
            record_counts['enzymes'] += 1

            # drug_enzyme_references
            drug_enzyme_references = enzyme.xpath('db:references', namespaces=ns)[0].text
            drug_enzyme_record['drug_enzyme_references'] = drug_enzyme_references.strip().replace('\n', '') \
                .replace('\r', '').replace('\t', ' ') if drug_enzyme_references else ''
            drug_enzyme_records[(drug_id, enzyme_id,)] = drug_enzyme_record
            record_counts['drug_enzyme'] += 1

            # process the actions of the drug on the enzyme
            for action in drug.xpath('db:enzymes/db:enzyme/db:actions/db:action', namespaces=ns):
                drug_enzyme_action_record = {}
                action = action.text
                drug_enzyme_action_records[(drug_id, enzyme_id, action,)] = drug_enzyme_action_record
                record_counts['drug_enzyme_action'] += 1

    # open output files files, using codecs to handle data such as 'Dihomo-γ-linolenic acid'
    print('Writing output files.')
    with codecs.open(drugs_out, 'w', encoding='utf-8') as drugs_f, \
            codecs.open(drug_target_out, 'w', encoding='utf-8') as drug_target_f, \
            codecs.open(drug_target_action_out, 'w', encoding='utf-8') as drug_target_action_f, \
            codecs.open(targets_out, 'w', encoding='utf-8') as targets_f, \
            codecs.open(drug_enzyme_out, 'w', encoding='utf-8') as drug_enzyme_f, \
            codecs.open(drug_enzyme_action_out, 'w', encoding='utf-8') as drug_enzyme_action_f, \
            codecs.open(enzymes_out, 'w', encoding='utf-8') as enzymes_f:

        # write output files from the dicts by iteration over each table, output file pair
        for table, file, key_names, value_names in (
                (drugs_records, drugs_f, drugs_key_names, drugs_value_names),
                (targets_records, targets_f, targets_key_names, targets_value_names),
                (drug_target_records, drug_target_f, drug_target_key_names, drug_target_value_names),
                (drug_target_action_records, drug_target_action_f, drug_target_action_key_names, drug_target_action_value_names),
                (enzymes_records, enzymes_f, enzymes_key_names, enzymes_value_names),
                (drug_enzyme_records, drug_enzyme_f, drug_enzyme_key_names, drug_enzyme_value_names),
                (drug_enzyme_action_records, drug_enzyme_action_f, drug_enzyme_action_key_names, drug_enzyme_action_value_names)
        ):
            # write output file column header
            k = '\t'.join(key_names)
            tab = '\t' if value_names else ''
            v = '\t'.join(value_names)
            file.write(u'{keys}{tab}{values}\n'.format(keys=k,  tab=tab,  values=v))

            # write a key, value record
            for key, value in table.items():
                k = '\t'.join(key)
                tab = '\t' if value_names else ''
                v = '\t'.join([value[n] for n in value_names])
                # the u is necessary to deal with non-ascii characters in the data, such as in 'Dihomo-γ-linolenic acid'
                file.write(u'{keys}{tab}{values}\n'.format(keys=k, tab=tab, values=v))

    print('Got:\n{}'.format('\n'.join(['\t{}: {} records'.format(c, record_counts[c]) for c in record_counts.keys()])))
    print('Done.')
    sys.exit()


def main():

    # command line arguments
    parser = argparse.ArgumentParser(
        description='writes relational database text files for drugs, drug_target, drug_target_action, and targets' +
                    'drug_enzyme, drug_enzyme_action and enzyme',
        epilog='drugbankxml2db 1.0β1 ©2014 Michael Yourshaw all rights reserved',
    )
    parser.add_argument('--input', '-i', required=True,
                        help='drug bank xml file downloaded from http://www.drugbank.ca/system' +
                        '/downloads/current/drugbank.xml.zip',)
    parser.add_argument('--output', '-o', required=False, default="default",
                        help='output directory for write relation tables, default is input directory',)
    args = parser.parse_args()

    run(input=args.input, output=args.output)


if __name__ == "__main__":
    main()
    print("Exit.")
