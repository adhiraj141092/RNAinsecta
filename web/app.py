from flask import request, render_template, url_for, redirect, flash
import os
import pandas as pd
import sys
import json, jsonify

sys.path.append('/usr/local/lib/python3.6/site-packages')
import re
import RNA
from flask import Flask, session
import pickle

app = Flask(__name__)


app.secret_key = 'xyzzyspoon!'


@app.route('/')
def home():
    return render_template('home.html')



@app.route('/predict', methods=['POST', 'GET'])
def predict():
    rna = request.form["rna"]
    session["rna"] = rna
    fc = RNA.fold_compound(rna)
    (ss, mfe) = fc.mfe()
    session["struc"] = ss
    clf1 = request.form.get('clf1')

    if clf1 == '2':
        classifier = pickle.load(open('../ML_models/SVM_Final.pkl', 'rb'))
    if clf1 == '1':
        classifier = pickle.load(open('../ML_models/RF_reg1.pkl', 'rb'))


    if mfe == 0:
        return render_template('predict.html', display='This is not a pre-miRNA')
    else:
        f = open("temp1.fasta", "w+")
        f.write(rna)
        f.close()
        os.system('RNAfold < temp1.fasta > rnafold.txt')
        os.system('perl trip.pl < rnafold.txt > triplet.csv')  # make it subprocess pipe
        os.system('perl stats.pl < rnafold.txt > stats.csv')
        os.system("gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r300   -sOutputFile=rna.png rna.ps")
        os.system("cp rna.ps static/rna.ps")
        os.system("cp rna.png static/rna.png")
        s1 = pd.read_csv("triplet.csv", sep='\s+')
        s2 = pd.read_csv("stats.csv", sep='\t')
        s1 = s1.fillna(0)
        s2 = s2.fillna(0)
        s2 = s2.iloc[:, :-1]

        qi = open("rnafold.txt", "r")
        text = qi.readlines()
        qi.close()
        match1 = re.findall("[(.][.)(]+[.)]", str(text))
        match2 = re.findall(r"[AUGC]{15,}", str(text))

        df01 = pd.DataFrame(match1, match2)
        df01 = df01.reset_index()
        df01.columns = ["sequence", "fold"]
        scon = pd.concat([s1, s2, df01], axis=1)

        scon.rename(columns={"Nmfe": "dG"}, inplace=True)
        scon["stem"] = scon["fold"].str.extract('(\(.*\))')
        scon["nstem"] = scon["stem"].apply(len)
        scon["stemte"] = scon["fold"].str.findall('(\(\.+\))')
        scon["nstemte"] = scon["stemte"].apply(len)
        scon["stemte"] = scon["fold"].str.extract('(\(\.+\))')
        scon["lenstemte"] = scon["stemte"].apply(len)

        scon["MFE1"] = scon["dG"] / scon["%G+C"]
        scon["stem0"] = scon["fold"].str.findall('(\({3,})')
        scon["n_stems"] = scon["stem0"].apply(len)
        scon["stem1"] = scon["fold"].str.findall('(\()')
        scon["total_base"] = scon["stem1"].apply(len)
        scon["MFE2"] = scon["dG"] / scon["n_stems"]
        scon["MFE3"] = scon["dG"] / scon["nstemte"]
        scon["MFE4"] = scon["dG"] / scon["total_base"]
        scon["avg_bp"] = scon["total_base"] / scon["n_stems"]

        scon1 = scon[
            ["A...", "A..(", "A.(.", "A.((", "A(..", "A(.(", "A((.", "A(((", "G...", "G..(", "G.(.", "G.((", "G(..",
             "G(.(",
             "G((.", "G(((", "C...", "C..(", "C.(.", "C.((", "C(..", "C(.(", "C((.", "C(((", "U...", "U..(", "U.(.",
             "U.((",
             "U(..", "U(.(", "U((.", "U(((", "Len", "A", "C", "G", "U", "G+C", "A+U", "AA", "AC", "AG", "AU", "CA",
             "CC",
             "CG",
             "CU", "GA", "GC", "GG", "GU", "UA", "UC", "UG", "UU", "%A", "%C", "%G", "%U", "%G+C", "%A+U", "%AA", "%AC",
             "%AG",
             "%AU", "%CA", "%CC", "%CG", "%CU", "%GA", "%GC", "%GG", "%GU", "%UA", "%UC", "%UG", "%UU", "pb", "Npb",
             "mfe",
             "dG", "Q", "NQ", "D", "ND", "nstem", "MFE1", "MFE2", "MFE3", "MFE4", "total_base", "n_stems", "avg_bp"]]

        b = scon[["sequence", "Len", "mfe", "fold", "avg_bp"]]
        scon1 = scon1.abs()


        scon1.to_csv('static/input.csv', header=True, index=False, float_format="%.10f")
        app.config['UPLOAD_FOLDER'] = "/static/"
        full_filename = os.path.join(app.config['UPLOAD_FOLDER'], 'rna.png')
        df = pd.read_csv("static/input.csv")
        X = df.iloc[:, :].values
        y_pred = classifier.predict(X)
        prob_w_br = classifier.predict_proba(X)[:, 1]
        prob = str(prob_w_br)[1:-1]
        proba = (classifier.predict_proba(X)[:, 1] >= 0.50).astype(int)

        p = re.compile("(\(\.+\))")
        m01 = []

        for m in p.finditer(ss):
            m01 += m.span()

        seq3 = []
        seq5 = []

        for i in range(len(m01)):
            if i % 2 == 0:
                seq3 += [m01[i]]
            if i % 2 == 1:
                seq5 += [m01[i]]
        fseq3 = json.dumps(seq3)
        fseq5 = json.dumps(seq5)

        if proba == 0:

            return render_template('predict01.html', display='This is not an insect pre-miRNA', mfe=mfe, str=ss,
                                   plot=full_filename, rna=rna, prob=prob, GC=df['%G+C'].values[0],
                                   Len=df['Len'].values[0])
        else:
            return render_template('predict002.html', display='This looks like an insect pre-miRNA', mfe=mfe, str=ss,
                                   plot=full_filename, rna=rna, prob=prob, GC=df['%G+C'].values[0],
                                   Len=df['Len'].values[0], seq3=fseq3, seq5=fseq5)


@app.route('/predict2', methods=['POST', 'GET'])
def predict2():
    rna = request.form["rna1"]
    clf = request.form.get('clf')

    if clf == '2':
        classifier = pickle.load(open('../ML_models/SVM_Final.pkl', 'rb'))
    if clf == '1':
        classifier = pickle.load(open('../ML_models/RF_reg1.pkl', 'rb'))
    # fc = RNA.fold_compound(rna)
    # (ss, mfe) = fc.mfe()
    # session["struc"] = ss
    # if mfe == 0:
    #     return render_template('predict.html', display='This is not a pre-miRNA')
    # else:
    f = open("temp1.fasta", "w+")
    f.write(rna)
    f.close()
    os.system('RNAfold < temp1.fasta > rnafold.txt')
    os.system('perl trip.pl < rnafold.txt > triplet.csv')  # make it subprocess pipe
    os.system('perl stats.pl < rnafold.txt > stats.csv')
    os.system("gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r300   -sOutputFile=rna.png rna.ps")
    os.system("cp rna.ps static/rna.ps")
    os.system("cp rna.png static/rna.png")
    os.system("rm *.ps")
    s1 = pd.read_csv("triplet.csv", sep='\s+')
    s2 = pd.read_csv("stats.csv", sep='\t')
    s1 = s1.fillna(0)
    s2 = s2.fillna(0)
    s2 = s2.iloc[:, :-1]
    qi = open("rnafold.txt", "r")
    text = qi.readlines()
    qi.close()
    match1 = re.findall("[(.][.)(]+[.)]", str(text))
    match2 = re.findall(r"[AUGC]{15,}", str(text))

    df01 = pd.DataFrame(match1, match2)
    df01 = df01.reset_index()
    df01.columns = ["sequence","fold"]
    scon = pd.concat([s1, s2, df01], axis=1)


    scon.rename(columns={"Nmfe": "dG"}, inplace=True)
    scon["stem"] = scon["fold"].str.extract('(\(.*\))')
    scon["nstem"] = scon["stem"].apply(len)
    scon["stemte"] = scon["fold"].str.findall('(\(\.+\))')
    scon["nstemte"] = scon["stemte"].apply(len)
    scon["stemte"] = scon["fold"].str.extract('(\(\.+\))')
    scon["lenstemte"] = scon["stemte"].apply(len)

    scon["MFE1"] = scon["dG"] / scon["%G+C"]
    scon["stem0"] = scon["fold"].str.findall('(\({3,})')
    scon["n_stems"] = scon["stem0"].apply(len)
    scon["stem1"] = scon["fold"].str.findall('(\()')
    scon["total_base"] = scon["stem1"].apply(len)
    scon["MFE2"] = scon["dG"] / scon["n_stems"]
    scon["MFE3"] = scon["dG"] / scon["nstemte"]
    scon["MFE4"] = scon["dG"] / scon["total_base"]
    scon["avg_bp"] = scon["total_base"] / scon["n_stems"]

    scon1 = scon[
        ["A...", "A..(", "A.(.", "A.((", "A(..", "A(.(", "A((.", "A(((", "G...", "G..(", "G.(.", "G.((", "G(..", "G(.(",
         "G((.", "G(((", "C...", "C..(", "C.(.", "C.((", "C(..", "C(.(", "C((.", "C(((", "U...", "U..(", "U.(.", "U.((",
         "U(..", "U(.(", "U((.", "U(((", "Len", "A", "C", "G", "U", "G+C", "A+U", "AA", "AC", "AG", "AU", "CA", "CC",
         "CG",
         "CU", "GA", "GC", "GG", "GU", "UA", "UC", "UG", "UU", "%A", "%C", "%G", "%U", "%G+C", "%A+U", "%AA", "%AC",
         "%AG",
         "%AU", "%CA", "%CC", "%CG", "%CU", "%GA", "%GC", "%GG", "%GU", "%UA", "%UC", "%UG", "%UU", "pb", "Npb", "mfe",
         "dG", "Q", "NQ", "D", "ND", "nstem", "MFE1", "MFE2", "MFE3", "MFE4", "total_base", "n_stems", "avg_bp"]]

    b = scon[["sequence","Len", "mfe", "fold", "avg_bp"]]

    scon1 = scon1.abs()

    scon1.to_csv('static/input.csv', header=True, index=False, float_format="%.10f")
    app.config['UPLOAD_FOLDER'] = "/static/"
    full_filename = os.path.join(app.config['UPLOAD_FOLDER'], 'rna.png')
    X = scon1.iloc[:, :].values

    #classifier = pickle.load(open('SVM_Final.pkl', 'rb'))
    y_pred = classifier.predict(X)
    prob_w_br = classifier.predict_proba(X)[:, 1]
    prob = str(prob_w_br)[1:-1]
    proba = (classifier.predict_proba(X)[:, 1] >= 0.50).astype(int)

    a = pd.DataFrame([proba, prob_w_br])
    a = a.T
    a.columns=["pred","proba"]

    df = pd.concat([b, a], axis=1)
    df.loc[df['pred'] > 0, 'Prediction'] = 'Pre-microRNA'
    df.loc[df['pred'] < 1, 'Prediction'] = 'Not a Pre-miRNA'

    result = df.to_json(orient="records")
    parsed = json.loads(result)
    return render_template('predict002_batch.html', results=json.dumps(parsed))


@app.route('/val', methods=['GET', 'POST'])
def val():
    s3 = request.form.get('getseq3')
    s5 = request.form.get('getseq5')
    ori = request.form.get('ori')
    rna = session['rna']
    se = ''
    lenr = len(rna)

    if s5 != '0':
        var = int(s5)
    else:
        var = int(s3)

    if ori == '3':
        for element in range(0, var):
            se += rna[element]

    if ori == '5':
        for element in range(var, lenr):
            se += rna[element]

    return f"3' is {s3}<br>5' is {s5}<br> The value is {var}<br>original RNA {rna}<br>Directionality: {ori}<br>Cropped RNA is {se}"
    # return (str(s5))


# com = "/media/ub/New_vol/project_miRNA/post_lockdown/target/miRanda-1.0b/bin/miranda tarRNA.fasta DM_tar/DM_" + chrom + ".fasta -en -50 -sc 150 > res1.txt"


@app.route('/mirtar', methods=['POST', 'GET'])
def mirtar():
    s3 = request.form.get('getseq3')
    s5 = request.form.get('getseq5')
    ori = request.form.get('ori')
    chrom = request.form.get('chrom')
    rna = session['rna']
    se = ''
    lenr = len(rna)

    if s5 != '0':
        var = int(s5)
    else:
        var = int(s3)

    if ori == '3':
        for element in range(0, var):
            se += rna[element]

    if ori == '5':
        for element in range(var, lenr):
            se += rna[element]

    f = open("tarRNA.fasta", "w+")
    f.write(se)
    f.close()

    com = "miRanda-1.0b/miRanda-1.0b/bin/miranda tarRNA.fasta DM_tar/DM_" + chrom + ".fasta > res1.txt"

    os.system(com)

    with open('res1.txt') as file:
        res = file.read()
        z = re.findall('^\>\t(\w+)\t(\d+.\d+)\t(\-?\d+.\d+)\t\d+.\d+\t(\d+\s\d+)\t(\d+\s\d+)\t(\d+)', res,
                       flags=re.IGNORECASE | re.MULTILINE | re.DOTALL)
        y = re.findall('(^.+Query:.+\n.+\n.+3\')', res, flags=re.IGNORECASE | re.MULTILINE)

        if len(z) == 0:
            return render_template('tarno.html')
        else:
            df = pd.DataFrame(z)
            df1 = pd.DataFrame(y)
            df = pd.concat([df1, df], axis=1)
            df = df[:20]
            df.columns = ['algn', 'name', "score", "energy", "Query_range", "Ref_range", "aln_len"]
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
            return render_template('mirtar.html', targets=json.dumps(parsed), mi=se, ori=ori, chrom=chrom)


@app.after_request
def add_header(response):
    # response.cache_control.no_store = True
    response.headers['Cache-Control'] = 'no-store, no-cache, must-revalidate, post-check=0, pre-check=0, max-age=0'
    response.headers['Pragma'] = 'no-cache'
    response.headers['Expires'] = '-1'
    return response


if __name__ == "__main__":
    app.run(debug=True)
