# RNAinsecta: A tool for prediction of pre-microRNA in insects using machine learning algorithms.

Pre-MicroRNAs are the hairpin loops which produces microRNAs that negatively regulate gene expression in several organisms. In insects, microRNAs participate in several biological processes including metamorphosis, reproduction, immune response, etc. 
In this work, we trained machine learning classifiers such as Random Forest, Support Vector Machine, Logistic Regression and k-Nearest Neighbours to predict pre-microRNA hairpin loops in insects while using Synthetic Minority Over-sampling Technique and Near-Miss to handle the class imbalance. The trained model on Support Vector Machine achieved accuracy of 92.19% while the Random Forest attained an accuracy of 82.4% on our validation dataset. These models are hosted online as web application called RNAinsecta. Further, searching target for the predicted pre-microRNA in insect model organism <i>Drosophila melanogaster</i> has been provided in RNAinsecta using miRanda at the backend where experimentally validated genes regulated by microRNA are collected from miRTarBase as target sites. 
RNAinsecta is currently hosted at https://rnainsecta.in
<br>
This repository consist of the source code for hosting the webserver as well as testing the Machine Learning models to replicate the results.
<br><br>

<b>Pre- Requisites:</b>
<li><a href="https://www.tbi.univie.ac.at/RNA/">ViennaRNA Package</a></li>
<li><a href="https://www.python.org/downloads/">Python3</a></li>
<li>Virtual Environment:</li>
<pre>python3 -m pip install --user virtualenv</pre>
<br>

<b>Installation:</b>
<li>Create a virtual environment in python3 engine:</li>
<pre>python3 -m venv my_project_env</pre>
<li>Activate Virtual Environment:</li>
<pre>source my_project_env/bin/activate</pre>
<li>Install the required packages using requirements.txt:</li>
<pre>pip install requirements.txt</pre>
<br>

<b>Testing:</b><br>
For testing, the sequences along with their secondary structure are provided in the Dataset directory. The test dataset consist of true insect pre-microRNA (pos.fold) and pseudo insect pre-microRNA (neg.fold) which are hairpin loops found in insects that closely resembles true pre-microRNA.
Each prediction for true pre-microRNA can give either True Positive (TP) or False Positive (FP) and likewise, pseudo pre-microRNA gives either True Negative (TN) or False Negative (FN). Using these parameters the accuracy, sensitivity, specificity, MCC and F1 scores are calculated. 
<br>
To test the model enter:
<pre>python3 testing.py</pre>

Results are stored in the newly created results directory.

<b>Deploying the Web-Server</b><br>
To run the web-server locally, execute the following commands:
<pre>cd web
python3 app.py</pre>
<br><br>
Acknowledgements:
<br>
RNAinsecta Uses<br>
<li>Perl scripts developed by Kwang Loong Stanley Ng and Santosh K. Mishra for feature calculation.</li><br>
<ol>Kwang Loong Stanley Ng, Santosh K. Mishra, De novo SVM classification of precursor microRNAs from genomic pseudo hairpins using global and intrinsic folding measures, Bioinformatics, Volume 23, Issue 11, June 2007, Pages 1321–1330, https://doi.org/10.1093/bioinformatics/btm026</ol><br><BR>
<li>miRanda: A microRNA target searching tool.</li><br>
<ol>John,B., Enright,A.J., Aravin,A., Tuschl,T., Sander,C. and Marks,D.S. (2005) Correction: Human MicroRNA Targets. PLoS Biol., 3, e264.</ol>

<b>Pre-print</b>
RNAinsecta: A tool for prediction of pre-microRNA in insects using machine learning algorithms.
Adhiraj Nath, Utpal Bora*
bioRxiv 2022.03.31.486617; doi: <a href="https://doi.org/10.1101/2022.03.31.486617">https://doi.org/10.1101/2022.03.31.486617
</a>
<br><br>
The web-server is created and maintained under GNU/GPL v3 license by:<br>
Adhiraj Nath<br>
IIT Guwahati<br>
PhD candidate<br>
email: adhiraj@iitg.ac.in<br>
Mobile: +91 87230 13467<br>
