import os
import re
import pandas as pd
import json, jsonify
os.system('/media/ub/New_vol/project_miRNA/post_lockdown/target/miRanda-1.0b/bin/miranda temp1.fasta DM_tar/DM_3L.fasta -en -50 -sc 150 > res1.txt')
with open('res1.txt') as file:
    res = file.read()
    z = re.findall('^\>\t(\w+)\t(\d+.\d+)\t(\-?\d+.\d+)\t\d+.\d+\t(\d+\s\d+)\t(\d+\s\d+)\t(\d+)', res,
                   flags=re.IGNORECASE | re.MULTILINE | re.DOTALL)
    y = re.findall('(^.+Query:.+\n.+\n.+3\')', res, flags=re.IGNORECASE | re.MULTILINE)
    df = pd.DataFrame(z)
    df1 = pd.DataFrame(y)
    df = pd.concat([df1, df], axis=1)
    df.columns = ['algn','name', "score", "energy", "Query_range", "Ref_range", "aln_len"]
    df["name"] = df["name"].astype("str")
    df["algn"] = df["algn"].astype("str")
    df['urlF'] = "http://flybase.org/reports/" + df.name

    pr = pd.read_csv("mirbase_ann.csv", sep=',')

    df = pd.merge(df,
                         pr,
                         on='name',
                         how='left')
    df['urlm'] = "http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=" + df.mirna
    df = df.sort_values('score', ascending=False)
    result = df.to_json(orient="records")
    parsed = json.loads(result)

    #############################################################
    ##########################################################

    fn = gffutils.example_filename('/media/ub/New_vol/project_miRNA/post_lockdown/mipipe/flask_mi/web/bmori_tar.gff3')
    print(open(fn).read())

    db = gffutils.create_db(fn, dbfn='test.db', force=True, keep_order=True, merge_strategy = 'merge', sort_attribute_values = True)
    gene = db.features_of_type("miRNA")
    ge = list(gene)

    http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=MI0012997

result = df.to_json(orient="records")
parsed = json.loads(result)

json.dumps(parsed, indent=1)