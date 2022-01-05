SPEP classifier
========

Natural language processing of serum protein electrophoresis reports,
using support vector machine and rules-based approaches.

PDF: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7477876/pdf/CCI.19.00167.pdf

PMC: https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/32813561/

PubMed: https://pubmed.ncbi.nlm.nih.gov/32813561/

Journal: https://www.doi.org/10.1200/CCI.19.00167

Ryu JH, Zimolzak AJ. Natural Language Processing of Serum Protein
Electrophoresis Reports in the Veterans Affairs Health Care System.
JCO Clin Cancer Inform. 2020 Aug;4:749-756. doi: 10.1200/CCI.19.00167.

Code is public domain, because created as work for hire by an employee
of US Federal Government.


Details about iterations of rules-based
========

First attempt
--------
    [1] 0.82
       yhat
    y     0   1
      0 126  41
      1  13 120

Second
--------
       yhat
    y     0   1
      0 166   1
      1  14 119
    [1] 0.95

Third
--------
       yhat
    y     0   1
      0 166   1
      1   9 124
    [1] 0.9666667

Fourth
--------
       yhat
    y     0   1
      0 166   1
      1   5 128
    [1] 0.98

After fixing 3 annotation errors
--------
       yhat
    y     0   1
      0 169   1
      1   2 128
    
    [1] 0.99
