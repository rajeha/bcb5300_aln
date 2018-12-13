## BCB5300 Implementation: Global Alignment

Perl script to globally align two sequences with a custom scoring matrix.



### Example:
```
$ global.pl ADQEGHILKPSTWYV ANDCQEGHILKMFPSTWYVBZX  scores.tsv
A-D-QEGHILK--PSTWYV---
ANDCQEGHILKMFPSTWYVBZX  
```

It takes a sequence (arg1), another sequence (arg2), and a scoring matrix space-delimited file (arg3) as input.
The output is the global alignment of the sequences.

The method is based on dynamic programming. It starts by calculating the cost of aligning every pair i,j in the matrix from neighboring cells and the score of aligning i,j (accounting for a gap). To return and build the alignment result, a trace matrix is kept track of. Then, the global alignment is built recursively from the tracing matrix.


