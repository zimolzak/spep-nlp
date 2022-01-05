SPEP classifier
========

Natural language processing of serum protein electrophoresis reports,
using support vector machine and rules-based approaches.


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
