from __future__ import division
import iotbx.pdb
from cctbx import  miller
from cctbx import maptbx
from libtbx.test_utils import approx_equal
from cctbx.maptbx import resolution_from_map_and_model

pdb_str_box = """
CRYST1   46.014   46.587   51.149  90.00  90.00  90.00 P 1
SCALE1      0.021733  0.000000  0.000000        0.00000
SCALE2      0.000000  0.021465  0.000000        0.00000
SCALE3      0.000000  0.000000  0.019551        0.00000
ATOM      1  N   MET A   1       8.989  25.311  35.844  1.00 55.86           N
ATOM      2  CA  MET A   1      10.212  24.965  36.614  1.00 51.72           C
ATOM      3  C   MET A   1      10.661  26.174  37.411  1.00 46.27           C
ATOM      4  O   MET A   1      10.536  27.310  36.941  1.00 45.16           O
ATOM      5  CB  MET A   1      11.351  24.546  35.668  1.00 53.97           C
ATOM      6  CG  MET A   1      11.560  23.056  35.505  1.00 60.99           C
ATOM      7  SD  MET A   1      13.341  22.661  35.464  1.00 68.34           S
ATOM      8  CE  MET A   1      13.263  20.867  35.618  1.00 67.72           C
ATOM      9  N   THR A   2      11.186  25.930  38.612  1.00 43.01           N
ATOM     10  CA  THR A   2      11.594  27.021  39.490  1.00 40.71           C
ATOM     11  C   THR A   2      12.921  27.574  38.997  1.00 40.79           C
ATOM     12  O   THR A   2      13.738  26.869  38.391  1.00 39.70           O
ATOM     13  CB  THR A   2      11.702  26.615  41.004  1.00 41.37           C
ATOM     14  OG1 THR A   2      12.738  25.639  41.197  1.00 42.89           O
ATOM     15  CG2 THR A   2      10.358  26.062  41.530  1.00 42.18           C
ATOM     16  N   VAL A   3      13.110  28.862  39.223  1.00 39.13           N
ATOM     17  CA  VAL A   3      14.389  29.465  38.976  1.00 40.22           C
ATOM     18  C   VAL A   3      15.412  28.627  39.747  1.00 39.03           C
ATOM     19  O   VAL A   3      16.456  28.254  39.190  1.00 41.08           O
ATOM     20  CB  VAL A   3      14.376  30.996  39.344  1.00 41.84           C
ATOM     21  CG1 VAL A   3      15.777  31.541  39.575  1.00 45.27           C
ATOM     22  CG2 VAL A   3      13.658  31.801  38.239  1.00 42.05           C
ATOM     23  N   GLU A   4      15.093  28.285  40.994  1.00 38.38           N
ATOM     24  CA  GLU A   4      16.000  27.484  41.821  1.00 40.46           C
ATOM     25  C   GLU A   4      16.451  26.224  41.057  1.00 40.45           C
ATOM     26  O   GLU A   4      17.648  26.028  40.825  1.00 40.05           O
ATOM     27  CB  GLU A   4      15.369  27.101  43.170  1.00 41.73           C
ATOM     28  CG  GLU A   4      16.361  26.345  44.075  1.00 48.81           C
ATOM     29  CD  GLU A   4      15.918  26.181  45.520  1.00 55.24           C
ATOM     30  OE1 GLU A   4      16.040  25.045  46.035  1.00 61.84           O
ATOM     31  OE2 GLU A   4      15.486  27.176  46.149  1.00 57.03           O
ATOM     32  N   GLN A   5      15.483  25.407  40.645  1.00 39.24           N
ATOM     33  CA  GLN A   5      15.787  24.163  39.936  1.00 41.50           C
ATOM     34  C   GLN A   5      16.425  24.396  38.567  1.00 41.35           C
ATOM     35  O   GLN A   5      17.343  23.659  38.201  1.00 40.60           O
ATOM     36  CB  GLN A   5      14.543  23.275  39.802  1.00 42.26           C
ATOM     37  CG  GLN A   5      14.061  22.661  41.141  1.00 47.54           C
ATOM     38  CD  GLN A   5      14.893  21.447  41.602  1.00 51.54           C
ATOM     39  OE1 GLN A   5      15.936  21.117  41.025  1.00 49.52           O
ATOM     40  NE2 GLN A   5      14.420  20.787  42.648  1.00 58.93           N
ATOM     41  N   MET A   6      15.973  25.422  37.830  1.00 41.22           N
ATOM     42  CA  MET A   6      16.519  25.678  36.488  1.00 44.06           C
ATOM     43  C   MET A   6      17.996  26.061  36.530  1.00 45.57           C
ATOM     44  O   MET A   6      18.796  25.628  35.679  1.00 48.61           O
ATOM     45  CB  MET A   6      15.735  26.780  35.777  1.00 46.25           C
ATOM     46  CG  MET A   6      14.421  26.312  35.209  1.00 51.71           C
ATOM     47  SD  MET A   6      13.513  27.625  34.359  1.00 64.51           S
ATOM     48  CE  MET A   6      13.144  28.725  35.720  1.00 58.26           C
HETATM   49  N   SME A   7      18.343  26.883  37.513  1.00 44.69           N
HETATM   50  CA  SME A   7      19.702  27.367  37.690  1.00 47.86           C
HETATM   51  C   SME A   7      20.552  26.208  38.098  1.00 47.33           C
HETATM   52  O   SME A   7      21.688  26.091  37.678  1.00 48.73           O
HETATM   53  CB  SME A   7      19.753  28.375  38.828  1.00 48.79           C
HETATM   54  CG  SME A   7      18.930  29.616  38.490  1.00 51.80           C
HETATM   55  CE  SME A   7      21.476  30.890  38.486  1.00 66.41           C
HETATM   56  OE  SME A   7      19.892  30.890  36.374  1.00 65.87           O
HETATM   57  S   SME A   7      19.876  30.883  37.865  1.00 64.46           S
ATOM     58  N   LYS A   8      19.983  25.358  38.942  1.00 45.03           N
ATOM     59  CA  LYS A   8      20.643  24.160  39.468  1.00 46.88           C
ATOM     60  C   LYS A   8      21.037  23.191  38.363  1.00 46.66           C
ATOM     61  O   LYS A   8      22.153  22.689  38.341  1.00 49.15           O
ATOM     62  CB  LYS A   8      19.688  23.464  40.441  1.00 45.97           C
ATOM     63  CG  LYS A   8      20.287  22.433  41.349  1.00 51.13           C
ATOM     64  CD  LYS A   8      19.193  21.950  42.323  1.00 55.04           C
ATOM     65  CE  LYS A   8      19.483  20.592  42.888  1.00 59.32           C
ATOM     66  NZ  LYS A   8      18.353  20.108  43.739  1.00 64.12           N
ATOM     67  N   SER A   9      20.102  22.936  37.452  1.00 45.05           N
ATOM     68  CA  SER A   9      20.318  22.019  36.330  1.00 46.72           C
ATOM     69  C   SER A   9      21.331  22.575  35.333  1.00 49.01           C
ATOM     70  O   SER A   9      22.195  21.850  34.838  1.00 50.72           O
ATOM     71  CB  SER A   9      18.989  21.720  35.616  1.00 46.76           C
ATOM     72  OG  SER A   9      19.207  21.170  34.321  1.00 50.67           O
ATOM     73  N   GLY A  10      21.221  23.864  35.035  1.00 36.96           N
ATOM     74  CA  GLY A  10      22.147  24.515  34.097  1.00 34.97           C
ATOM     75  C   GLY A  10      23.553  24.538  34.661  1.00 34.57           C
ATOM     76  O   GLY A  10      24.529  24.357  33.934  1.00 31.75           O
ATOM     77  N   GLU A  11      23.648  24.730  35.971  1.00 39.42           N
ATOM     78  CA  GLU A  11      24.915  24.561  36.658  1.00 42.55           C
ATOM     79  C   GLU A  11      25.414  23.126  36.611  1.00 42.05           C
ATOM     80  O   GLU A  11      26.581  22.892  36.335  1.00 41.32           O
ATOM     81  CB  GLU A  11      24.832  25.019  38.120  1.00 47.91           C
ATOM     82  CG  GLU A  11      25.325  26.463  38.332  1.00 54.32           C
ATOM     83  CD  GLU A  11      24.245  27.413  38.798  1.00 60.17           C
ATOM     84  OE1 GLU A  11      23.801  27.271  39.965  1.00 62.75           O
ATOM     85  OE2 GLU A  11      23.861  28.305  38.005  1.00 62.45           O
ATOM     86  N   MET A  12      24.529  22.172  36.914  1.00 43.46           N
ATOM     87  CA  MET A  12      24.862  20.738  36.866  1.00 43.08           C
ATOM     88  C   MET A  12      25.434  20.386  35.506  1.00 39.76           C
ATOM     89  O   MET A  12      26.497  19.781  35.399  1.00 40.17           O
ATOM     90  CB  MET A  12      23.618  19.867  37.117  1.00 44.16           C
ATOM     91  CG  MET A  12      23.846  18.330  36.987  1.00 44.29           C
ATOM     92  SD  MET A  12      23.573  17.627  35.315  1.00 40.29           S
ATOM     93  CE  MET A  12      21.773  17.620  35.210  1.00 39.82           C
ATOM     94  N   ILE A  13      24.709  20.780  34.470  1.00 37.62           N
ATOM     95  CA  ILE A  13      25.084  20.460  33.098  1.00 35.09           C
ATOM     96  C   ILE A  13      26.439  21.072  32.753  1.00 33.74           C
ATOM     97  O   ILE A  13      27.267  20.436  32.121  1.00 33.71           O
ATOM     98  CB  ILE A  13      24.015  20.912  32.103  1.00 33.70           C
ATOM     99  CG1 ILE A  13      22.795  20.012  32.238  1.00 36.47           C
ATOM    100  CG2 ILE A  13      24.550  20.839  30.663  1.00 31.03           C
ATOM    101  CD1 ILE A  13      21.527  20.571  31.689  1.00 36.99           C
ATOM    102  N   ARG A  14      26.666  22.309  33.178  1.00 35.83           N
ATOM    103  CA  ARG A  14      27.939  22.966  32.934  1.00 34.63           C
ATOM    104  C   ARG A  14      29.025  22.222  33.671  1.00 35.93           C
ATOM    105  O   ARG A  14      30.079  22.011  33.103  1.00 34.08           O
ATOM    106  CB  ARG A  14      27.882  24.424  33.383  1.00 37.11           C
ATOM    107  CG  ARG A  14      29.163  25.264  33.267  1.00 39.17           C
ATOM    108  CD  ARG A  14      28.913  26.636  33.966  1.00 44.41           C
ATOM    109  NE  ARG A  14      30.015  27.595  33.842  1.00 49.92           N
ATOM    110  CZ  ARG A  14      30.150  28.517  32.873  1.00 50.55           C
ATOM    111  NH1 ARG A  14      29.255  28.650  31.887  1.00 46.69           N
ATOM    112  NH2 ARG A  14      31.208  29.327  32.892  1.00 51.05           N
ATOM    113  N   SER A  15      28.776  21.846  34.936  1.00 39.02           N
ATOM    114  CA  SER A  15      29.757  21.091  35.714  1.00 40.95           C
ATOM    115  C   SER A  15      30.134  19.807  34.999  1.00 39.30           C
ATOM    116  O   SER A  15      31.300  19.503  34.854  1.00 42.27           O
ATOM    117  CB  SER A  15      29.250  20.735  37.120  1.00 45.17           C
ATOM    118  OG  SER A  15      28.850  21.877  37.861  1.00 50.58           O
ATOM    119  N   VAL A  16      29.144  19.036  34.568  1.00 37.98           N
ATOM    120  CA  VAL A  16      29.422  17.746  33.929  1.00 37.22           C
ATOM    121  C   VAL A  16      30.258  17.940  32.707  1.00 35.88           C
ATOM    122  O   VAL A  16      31.251  17.251  32.505  1.00 37.43           O
ATOM    123  CB  VAL A  16      28.132  16.957  33.616  1.00 36.91           C
ATOM    124  CG1 VAL A  16      28.426  15.792  32.690  1.00 40.10           C
ATOM    125  CG2 VAL A  16      27.517  16.432  34.946  1.00 38.88           C
ATOM    126  N   CYS A  17      29.888  18.938  31.916  1.00 33.28           N
ATOM    127  CA  CYS A  17      30.457  19.112  30.606  1.00 32.19           C
ATOM    128  C   CYS A  17      31.839  19.730  30.633  1.00 32.65           C
ATOM    129  O   CYS A  17      32.705  19.319  29.845  1.00 32.70           O
ATOM    130  CB  CYS A  17      29.484  19.882  29.729  1.00 29.99           C
ATOM    131  SG  CYS A  17      28.090  18.841  29.331  1.00 29.57           S
ATOM    132  N   LEU A  18      32.061  20.689  31.541  1.00 32.85           N
ATOM    133  CA  LEU A  18      33.419  21.192  31.770  1.00 33.50           C
ATOM    134  C   LEU A  18      34.334  20.057  32.259  1.00 37.40           C
ATOM    135  O   LEU A  18      35.526  20.029  31.940  1.00 38.44           O
ATOM    136  CB  LEU A  18      33.446  22.308  32.810  1.00 34.41           C
ATOM    137  CG  LEU A  18      32.901  23.703  32.531  1.00 35.31           C
ATOM    138  CD1 LEU A  18      32.910  24.557  33.814  1.00 36.98           C
ATOM    139  CD2 LEU A  18      33.689  24.372  31.360  1.00 34.44           C
ATOM    140  N   GLY A  19      33.774  19.153  33.063  1.00 40.28           N
ATOM    141  CA  GLY A  19      34.506  18.011  33.590  1.00 44.13           C
ATOM    142  C   GLY A  19      34.947  17.026  32.529  1.00 44.20           C
ATOM    143  O   GLY A  19      36.065  16.507  32.601  1.00 46.42           O
ATOM    144  N   LYS A  20      34.078  16.758  31.558  1.00 41.93           N
ATOM    145  CA  LYS A  20      34.367  15.761  30.513  1.00 42.83           C
ATOM    146  C   LYS A  20      35.376  16.242  29.482  1.00 42.66           C
ATOM    147  O   LYS A  20      36.147  15.443  28.952  1.00 42.05           O
ATOM    148  CB  LYS A  20      33.089  15.330  29.792  1.00 41.52           C
ATOM    149  CG  LYS A  20      32.245  14.346  30.623  1.00 46.85           C
ATOM    150  CD  LYS A  20      32.516  12.892  30.192  1.00 53.93           C
ATOM    151  CE  LYS A  20      32.342  11.903  31.321  1.00 56.45           C
ATOM    152  NZ  LYS A  20      33.639  11.254  31.736  1.00 64.84           N
ATOM    153  N   THR A  21      35.374  17.544  29.217  1.00 40.77           N
ATOM    154  CA  THR A  21      36.176  18.117  28.142  1.00 40.50           C
ATOM    155  C   THR A  21      37.378  18.992  28.602  1.00 42.78           C
ATOM    156  O   THR A  21      38.249  19.315  27.793  1.00 44.01           O
ATOM    157  CB  THR A  21      35.255  18.954  27.232  1.00 37.48           C
ATOM    158  OG1 THR A  21      34.872  20.144  27.911  1.00 35.35           O
ATOM    159  CG2 THR A  21      33.989  18.173  26.837  1.00 35.67           C
ATOM    160  N   LYS A  22      37.421  19.391  29.876  1.00 44.63           N
ATOM    161  CA  LYS A  22      38.514  20.235  30.393  1.00 45.73           C
ATOM    162  C   LYS A  22      38.768  21.527  29.585  1.00 44.09           C
ATOM    163  O   LYS A  22      39.916  21.969  29.460  1.00 48.05           O
ATOM    164  CB  LYS A  22      39.815  19.414  30.482  1.00 47.71           C
TER
END
"""

pdb_str_5 = """
CRYST1   46.014   46.587   51.149  90.00  90.00  90.00 P 1
SCALE1      0.021733  0.000000  0.000000        0.00000
SCALE2      0.000000  0.021465  0.000000        0.00000
SCALE3      0.000000  0.000000  0.019551        0.00000
ATOM      1  N   MET A   1      12.998  26.982  38.701  1.00 55.86           N
ATOM      2  CA  MET A   1      12.064  26.130  39.371  1.00 51.72           C
ATOM      3  C   MET A   1      12.541  25.824  40.765  1.00 46.27           C
ATOM      4  O   MET A   1      13.020  26.753  41.457  1.00 45.16           O
ATOM      5  CB  MET A   1      11.841  24.824  38.573  1.00 53.97           C
ATOM      6  CG  MET A   1      10.387  24.301  38.632  1.00 60.99           C
ATOM      7  SD  MET A   1      10.319  22.528  38.121  1.00 68.34           S
ATOM      8  CE  MET A   1      11.465  22.592  36.754  1.00 67.72           C
ATOM      9  N   THR A   2      12.430  24.571  41.223  1.00 43.01           N
ATOM     10  CA  THR A   2      12.760  24.214  42.608  1.00 40.71           C
ATOM     11  C   THR A   2      14.263  24.195  42.844  1.00 40.79           C
ATOM     12  O   THR A   2      15.054  23.955  41.944  1.00 39.70           O
ATOM     13  CB  THR A   2      12.183  22.840  42.995  1.00 41.37           C
ATOM     14  OG1 THR A   2      12.764  21.861  42.123  1.00 42.89           O
ATOM     15  CG2 THR A   2      10.699  22.856  42.865  1.00 42.18           C
ATOM     16  N   VAL A   3      14.641  24.412  44.125  1.00 39.13           N
ATOM     17  CA  VAL A   3      16.065  24.373  44.476  1.00 40.22           C
ATOM     18  C   VAL A   3      16.650  23.058  44.005  1.00 39.03           C
ATOM     19  O   VAL A   3      17.804  22.963  43.597  1.00 41.08           O
ATOM     20  CB  VAL A   3      16.274  24.578  46.008  1.00 41.84           C
ATOM     21  CG1 VAL A   3      17.646  24.086  46.466  1.00 45.27           C
ATOM     22  CG2 VAL A   3      16.040  26.047  46.395  1.00 42.05           C
ATOM     23  N   GLU A   4      15.842  22.007  44.046  1.00 38.38           N
ATOM     24  CA  GLU A   4      16.317  20.733  43.565  1.00 40.46           C
ATOM     25  C   GLU A   4      16.749  20.835  42.111  1.00 40.45           C
ATOM     26  O   GLU A   4      17.915  20.605  41.763  1.00 40.05           O
ATOM     27  CB  GLU A   4      15.245  19.679  43.726  1.00 41.73           C
ATOM     28  CG  GLU A   4      15.722  18.284  43.201  1.00 48.81           C
ATOM     29  CD  GLU A   4      15.580  17.239  44.253  1.00 55.24           C
ATOM     30  OE1 GLU A   4      16.517  17.091  45.034  1.00 61.84           O
ATOM     31  OE2 GLU A   4      14.519  16.594  44.350  1.00 57.03           O
ATOM     32  N   GLN A   5      15.802  21.158  41.231  1.00 39.24           N
ATOM     33  CA  GLN A   5      16.045  21.191  39.789  1.00 41.50           C
ATOM     34  C   GLN A   5      17.171  22.167  39.405  1.00 41.35           C
ATOM     35  O   GLN A   5      17.958  21.904  38.491  1.00 40.60           O
ATOM     36  CB  GLN A   5      14.722  21.505  39.053  1.00 42.26           C
ATOM     37  CG  GLN A   5      13.640  20.538  39.435  1.00 47.54           C
ATOM     38  CD  GLN A   5      14.002  19.109  39.010  1.00 51.54           C
ATOM     39  OE1 GLN A   5      14.738  18.953  38.045  1.00 49.52           O
ATOM     40  NE2 GLN A   5      13.449  18.086  39.660  1.00 58.93           N
ATOM     41  N   MET A   6      17.241  23.302  40.078  1.00 41.22           N
ATOM     42  CA  MET A   6      18.256  24.283  39.690  1.00 44.06           C
ATOM     43  C   MET A   6      19.628  23.887  40.196  1.00 45.57           C
ATOM     44  O   MET A   6      20.616  23.937  39.409  1.00 48.61           O
ATOM     45  CB  MET A   6      17.943  25.686  40.233  1.00 46.25           C
ATOM     46  CG  MET A   6      16.602  26.222  39.848  1.00 51.71           C
ATOM     47  SD  MET A   6      15.994  27.653  40.808  1.00 64.51           S
ATOM     48  CE  MET A   6      16.951  29.028  40.150  1.00 58.26           C
HETATM   49  N   SME A   7      19.727  23.474  41.452  1.00 44.69           N
HETATM   50  CA  SME A   7      21.022  23.098  42.112  1.00 47.86           C
HETATM   51  C   SME A   7      21.613  21.859  41.495  1.00 47.33           C
HETATM   52  O   SME A   7      22.746  21.693  41.056  1.00 48.73           O
HETATM   53  CB  SME A   7      20.891  22.855  43.631  1.00 48.79           C
HETATM   54  CG  SME A   7      21.759  23.780  44.433  1.00 51.80           C
HETATM   55  CE  SME A   7      23.144  23.218  46.623  1.00 66.41           C
HETATM   56  OE  SME A   7      23.168  21.603  44.662  1.00 65.87           O
HETATM   57  S   SME A   7      23.313  23.075  44.898  1.00 64.46           S
ATOM     58  N   LYS A   8      20.731  20.865  41.492  1.00 45.03           N
ATOM     59  CA  LYS A   8      21.070  19.546  40.953  1.00 46.88           C
ATOM     60  C   LYS A   8      21.261  19.557  39.429  1.00 46.66           C
ATOM     61  O   LYS A   8      22.108  18.821  38.908  1.00 49.15           O
ATOM     62  CB  LYS A   8      19.969  18.606  41.300  1.00 45.97           C
ATOM     63  CG  LYS A   8      20.352  17.127  41.054  1.00 51.13           C
ATOM     64  CD  LYS A   8      19.846  16.355  42.252  1.00 55.04           C
ATOM     65  CE  LYS A   8      19.359  15.012  41.841  1.00 59.32           C
ATOM     66  NZ  LYS A   8      18.937  14.155  42.951  1.00 64.12           N
ATOM     67  N   SER A   9      20.455  20.367  38.737  1.00 45.05           N
ATOM     68  CA  SER A   9      20.578  20.448  37.300  1.00 46.72           C
ATOM     69  C   SER A   9      21.860  21.145  36.916  1.00 49.01           C
ATOM     70  O   SER A   9      22.612  20.643  36.083  1.00 50.72           O
ATOM     71  CB  SER A   9      19.369  21.146  36.679  1.00 46.76           C
ATOM     72  OG  SER A   9      19.636  21.471  35.301  1.00 50.67           O
ATOM     73  N   GLY A  10      22.149  22.269  37.591  1.00 36.96           N
ATOM     74  CA  GLY A  10      23.376  22.991  37.285  1.00 34.97           C
ATOM     75  C   GLY A  10      24.607  22.198  37.667  1.00 34.57           C
ATOM     76  O   GLY A  10      25.600  22.163  36.948  1.00 31.75           O
ATOM     77  N   GLU A  11      24.560  21.578  38.828  1.00 39.42           N
ATOM     78  CA  GLU A  11      25.646  20.707  39.311  1.00 42.55           C
ATOM     79  C   GLU A  11      25.867  19.512  38.405  1.00 42.05           C
ATOM     80  O   GLU A  11      27.003  19.115  38.124  1.00 41.32           O
ATOM     81  CB  GLU A  11      25.293  20.265  40.733  1.00 47.91           C
ATOM     82  CG  GLU A  11      26.384  19.458  41.398  1.00 54.32           C
ATOM     83  CD  GLU A  11      27.716  20.215  41.418  1.00 60.17           C
ATOM     84  OE1 GLU A  11      27.680  21.418  41.701  1.00 62.75           O
ATOM     85  OE2 GLU A  11      28.758  19.607  41.203  1.00 62.45           O
ATOM     86  N   MET A  12      24.801  18.925  37.938  1.00 43.46           N
ATOM     87  CA  MET A  12      24.927  17.820  36.992  1.00 43.08           C
ATOM     88  C   MET A  12      25.546  18.283  35.671  1.00 39.76           C
ATOM     89  O   MET A  12      26.514  17.688  35.206  1.00 40.17           O
ATOM     90  CB  MET A  12      23.541  17.186  36.751  1.00 44.16           C
ATOM     91  CG  MET A  12      23.467  15.758  36.225  1.00 44.29           C
ATOM     92  SD  MET A  12      21.833  15.045  35.960  1.00 40.29           S
ATOM     93  CE  MET A  12      22.311  13.286  35.717  1.00 39.82           C
ATOM     94  N   ILE A  13      24.984  19.335  35.038  1.00 37.62           N
ATOM     95  CA  ILE A  13      25.466  19.762  33.721  1.00 35.09           C
ATOM     96  C   ILE A  13      26.898  20.240  33.800  1.00 33.74           C
ATOM     97  O   ILE A  13      27.756  19.815  33.008  1.00 33.71           O
ATOM     98  CB  ILE A  13      24.535  20.872  33.187  1.00 33.70           C
ATOM     99  CG1 ILE A  13      23.089  20.349  33.073  1.00 36.47           C
ATOM    100  CG2 ILE A  13      25.065  21.380  31.842  1.00 31.03           C
ATOM    101  CD1 ILE A  13      22.120  21.142  32.100  1.00 36.99           C
ATOM    102  N   ARG A  14      27.184  21.101  34.798  1.00 35.83           N
ATOM    103  CA  ARG A  14      28.542  21.611  35.053  1.00 34.63           C
ATOM    104  C   ARG A  14      29.486  20.441  35.298  1.00 35.93           C
ATOM    105  O   ARG A  14      30.409  20.244  34.514  1.00 34.08           O
ATOM    106  CB  ARG A  14      28.500  22.648  36.199  1.00 37.11           C
ATOM    107  CG  ARG A  14      29.810  22.934  36.869  1.00 39.17           C
ATOM    108  CD  ARG A  14      29.598  23.433  38.323  1.00 44.41           C
ATOM    109  NE  ARG A  14      30.848  23.570  39.045  1.00 49.92           N
ATOM    110  CZ  ARG A  14      31.706  24.587  38.876  1.00 50.55           C
ATOM    111  NH1 ARG A  14      31.419  25.497  37.982  1.00 46.69           N
ATOM    112  NH2 ARG A  14      32.804  24.741  39.586  1.00 51.05           N
ATOM    113  N   SER A  15      29.184  19.566  36.265  1.00 39.02           N
ATOM    114  CA  SER A  15      30.148  18.527  36.626  1.00 40.95           C
ATOM    115  C   SER A  15      30.376  17.566  35.480  1.00 39.30           C
ATOM    116  O   SER A  15      31.528  17.249  35.111  1.00 42.27           O
ATOM    117  CB  SER A  15      29.602  17.741  37.825  1.00 45.17           C
ATOM    118  OG  SER A  15      29.743  18.374  39.100  1.00 50.58           O
ATOM    119  N   VAL A  16      29.306  17.091  34.890  1.00 37.98           N
ATOM    120  CA  VAL A  16      29.452  16.039  33.860  1.00 37.22           C
ATOM    121  C   VAL A  16      30.109  16.586  32.595  1.00 35.88           C
ATOM    122  O   VAL A  16      31.077  16.065  32.082  1.00 37.43           O
ATOM    123  CB  VAL A  16      28.113  15.374  33.553  1.00 36.91           C
ATOM    124  CG1 VAL A  16      28.291  14.431  32.371  1.00 40.10           C
ATOM    125  CG2 VAL A  16      27.613  14.599  34.763  1.00 38.88           C
ATOM    126  N   CYS A  17      29.572  17.693  32.090  1.00 33.28           N
ATOM    127  CA  CYS A  17      30.128  18.251  30.858  1.00 32.19           C
ATOM    128  C   CYS A  17      31.535  18.829  31.061  1.00 32.65           C
ATOM    129  O   CYS A  17      32.359  18.795  30.135  1.00 32.70           O
ATOM    130  CB  CYS A  17      29.203  19.340  30.353  1.00 29.99           C
ATOM    131  SG  CYS A  17      27.835  18.668  29.321  1.00 29.57           S
ATOM    132  N   LEU A  18      31.796  19.377  32.240  1.00 32.85           N
ATOM    133  CA  LEU A  18      33.128  19.914  32.524  1.00 33.50           C
ATOM    134  C   LEU A  18      34.169  18.799  32.677  1.00 37.40           C
ATOM    135  O   LEU A  18      35.349  18.974  32.398  1.00 38.44           O
ATOM    136  CB  LEU A  18      33.083  20.741  33.779  1.00 34.41           C
ATOM    137  CG  LEU A  18      34.371  21.301  34.321  1.00 35.31           C
ATOM    138  CD1 LEU A  18      35.057  22.222  33.327  1.00 36.98           C
ATOM    139  CD2 LEU A  18      34.173  21.979  35.744  1.00 34.44           C
ATOM    140  N   GLY A  19      33.779  17.627  33.161  1.00 40.28           N
ATOM    141  CA  GLY A  19      34.635  16.474  33.045  1.00 44.13           C
ATOM    142  C   GLY A  19      34.730  15.975  31.615  1.00 44.20           C
ATOM    143  O   GLY A  19      35.681  15.250  31.269  1.00 46.42           O
ATOM    144  N   LYS A  20      33.735  16.258  30.773  1.00 41.93           N
ATOM    145  CA  LYS A  20      33.819  15.832  29.374  1.00 42.83           C
ATOM    146  C   LYS A  20      34.745  16.718  28.549  1.00 42.66           C
ATOM    147  O   LYS A  20      35.427  16.243  27.638  1.00 42.05           O
ATOM    148  CB  LYS A  20      32.431  15.789  28.727  1.00 41.52           C
ATOM    149  CG  LYS A  20      31.613  14.516  28.974  1.00 46.85           C
ATOM    150  CD  LYS A  20      32.440  13.297  28.518  1.00 53.93           C
ATOM    151  CE  LYS A  20      32.197  12.935  27.029  1.00 56.45           C
ATOM    152  NZ  LYS A  20      33.346  13.503  26.261  1.00 64.84           N
ATOM    153  N   THR A  21      34.805  18.007  28.797  1.00 40.77           N
ATOM    154  CA  THR A  21      35.593  18.896  27.954  1.00 40.50           C
ATOM    155  C   THR A  21      36.842  19.399  28.654  1.00 42.78           C
ATOM    156  O   THR A  21      37.823  19.657  27.983  1.00 44.01           O
ATOM    157  CB  THR A  21      34.758  20.107  27.476  1.00 37.48           C
ATOM    158  OG1 THR A  21      33.429  19.990  28.020  1.00 35.35           O
ATOM    159  CG2 THR A  21      34.672  20.088  25.972  1.00 35.67           C
ATOM    160  N   LYS A  22      36.815  19.486  29.974  1.00 44.63           N
ATOM    161  CA  LYS A  22      37.951  19.848  30.817  1.00 45.73           C
ATOM    162  C   LYS A  22      38.543  21.187  30.359  1.00 44.09           C
ATOM    163  O   LYS A  22      39.660  21.265  29.818  1.00 48.05           O
ATOM    164  CB  LYS A  22      39.014  18.761  30.772  1.00 47.71           C
TER
END
"""

pdb_str_5bb = """
CRYST1   46.014   46.587   51.149  90.00  90.00  90.00 P 1
SCALE1      0.021733  0.000000  0.000000        0.00000
SCALE2      0.000000  0.021465  0.000000        0.00000
SCALE3      0.000000  0.000000  0.019551        0.00000
ATOM      1  N   MET A   1      12.998  26.982  38.701  1.00 55.86           N
ATOM      2  CA  MET A   1      12.064  26.130  39.371  1.00 51.72           C
ATOM      3  C   MET A   1      12.541  25.824  40.765  1.00 46.27           C
ATOM      4  O   MET A   1      13.020  26.753  41.457  1.00 45.16           O
ATOM      5  N   THR A   2      12.430  24.571  41.223  1.00 43.01           N
ATOM      6  CA  THR A   2      12.760  24.214  42.608  1.00 40.71           C
ATOM      7  C   THR A   2      14.263  24.195  42.844  1.00 40.79           C
ATOM      8  O   THR A   2      15.054  23.955  41.944  1.00 39.70           O
ATOM      9  N   VAL A   3      14.641  24.412  44.125  1.00 39.13           N
ATOM     10  CA  VAL A   3      16.065  24.373  44.476  1.00 40.22           C
ATOM     11  C   VAL A   3      16.650  23.058  44.005  1.00 39.03           C
ATOM     12  O   VAL A   3      17.804  22.963  43.597  1.00 41.08           O
ATOM     13  N   GLU A   4      15.842  22.007  44.046  1.00 38.38           N
ATOM     14  CA  GLU A   4      16.317  20.733  43.565  1.00 40.46           C
ATOM     15  C   GLU A   4      16.749  20.835  42.111  1.00 40.45           C
ATOM     16  O   GLU A   4      17.915  20.605  41.763  1.00 40.05           O
ATOM     17  N   GLN A   5      15.802  21.158  41.231  1.00 39.24           N
ATOM     18  CA  GLN A   5      16.045  21.191  39.789  1.00 41.50           C
ATOM     19  C   GLN A   5      17.171  22.167  39.405  1.00 41.35           C
ATOM     20  O   GLN A   5      17.958  21.904  38.491  1.00 40.60           O
ATOM     21  N   MET A   6      17.241  23.302  40.078  1.00 41.22           N
ATOM     22  CA  MET A   6      18.256  24.283  39.690  1.00 44.06           C
ATOM     23  C   MET A   6      19.628  23.887  40.196  1.00 45.57           C
ATOM     24  O   MET A   6      20.616  23.937  39.409  1.00 48.61           O
HETATM   25  N   SME A   7      19.727  23.474  41.452  1.00 44.69           N
HETATM   26  CA  SME A   7      21.022  23.098  42.112  1.00 47.86           C
HETATM   27  C   SME A   7      21.613  21.859  41.495  1.00 47.33           C
HETATM   28  O   SME A   7      22.746  21.693  41.056  1.00 48.73           O
ATOM     29  N   LYS A   8      20.731  20.865  41.492  1.00 45.03           N
ATOM     30  CA  LYS A   8      21.070  19.546  40.953  1.00 46.88           C
ATOM     31  C   LYS A   8      21.261  19.557  39.429  1.00 46.66           C
ATOM     32  O   LYS A   8      22.108  18.821  38.908  1.00 49.15           O
ATOM     33  N   SER A   9      20.455  20.367  38.737  1.00 45.05           N
ATOM     34  CA  SER A   9      20.578  20.448  37.300  1.00 46.72           C
ATOM     35  C   SER A   9      21.860  21.145  36.916  1.00 49.01           C
ATOM     36  O   SER A   9      22.612  20.643  36.083  1.00 50.72           O
ATOM     37  N   GLY A  10      22.149  22.269  37.591  1.00 36.96           N
ATOM     38  CA  GLY A  10      23.376  22.991  37.285  1.00 34.97           C
ATOM     39  C   GLY A  10      24.607  22.198  37.667  1.00 34.57           C
ATOM     40  O   GLY A  10      25.600  22.163  36.948  1.00 31.75           O
ATOM     41  N   GLU A  11      24.560  21.578  38.828  1.00 39.42           N
ATOM     42  CA  GLU A  11      25.646  20.707  39.311  1.00 42.55           C
ATOM     43  C   GLU A  11      25.867  19.512  38.405  1.00 42.05           C
ATOM     44  O   GLU A  11      27.003  19.115  38.124  1.00 41.32           O
ATOM     45  N   MET A  12      24.801  18.925  37.938  1.00 43.46           N
ATOM     46  CA  MET A  12      24.927  17.820  36.992  1.00 43.08           C
ATOM     47  C   MET A  12      25.546  18.283  35.671  1.00 39.76           C
ATOM     48  O   MET A  12      26.514  17.688  35.206  1.00 40.17           O
ATOM     49  N   ILE A  13      24.984  19.335  35.038  1.00 37.62           N
ATOM     50  CA  ILE A  13      25.466  19.762  33.721  1.00 35.09           C
ATOM     51  C   ILE A  13      26.898  20.240  33.800  1.00 33.74           C
ATOM     52  O   ILE A  13      27.756  19.815  33.008  1.00 33.71           O
ATOM     53  N   ARG A  14      27.184  21.101  34.798  1.00 35.83           N
ATOM     54  CA  ARG A  14      28.542  21.611  35.053  1.00 34.63           C
ATOM     55  C   ARG A  14      29.486  20.441  35.298  1.00 35.93           C
ATOM     56  O   ARG A  14      30.409  20.244  34.514  1.00 34.08           O
ATOM     57  N   SER A  15      29.184  19.566  36.265  1.00 39.02           N
ATOM     58  CA  SER A  15      30.148  18.527  36.626  1.00 40.95           C
ATOM     59  C   SER A  15      30.376  17.566  35.480  1.00 39.30           C
ATOM     60  O   SER A  15      31.528  17.249  35.111  1.00 42.27           O
ATOM     61  N   VAL A  16      29.306  17.091  34.890  1.00 37.98           N
ATOM     62  CA  VAL A  16      29.452  16.039  33.860  1.00 37.22           C
ATOM     63  C   VAL A  16      30.109  16.586  32.595  1.00 35.88           C
ATOM     64  O   VAL A  16      31.077  16.065  32.082  1.00 37.43           O
ATOM     65  N   CYS A  17      29.572  17.693  32.090  1.00 33.28           N
ATOM     66  CA  CYS A  17      30.128  18.251  30.858  1.00 32.19           C
ATOM     67  C   CYS A  17      31.535  18.829  31.061  1.00 32.65           C
ATOM     68  O   CYS A  17      32.359  18.795  30.135  1.00 32.70           O
ATOM     69  N   LEU A  18      31.796  19.377  32.240  1.00 32.85           N
ATOM     70  CA  LEU A  18      33.128  19.914  32.524  1.00 33.50           C
ATOM     71  C   LEU A  18      34.169  18.799  32.677  1.00 37.40           C
ATOM     72  O   LEU A  18      35.349  18.974  32.398  1.00 38.44           O
ATOM     73  N   GLY A  19      33.779  17.627  33.161  1.00 40.28           N
ATOM     74  CA  GLY A  19      34.635  16.474  33.045  1.00 44.13           C
ATOM     75  C   GLY A  19      34.730  15.975  31.615  1.00 44.20           C
ATOM     76  O   GLY A  19      35.681  15.250  31.269  1.00 46.42           O
ATOM     77  N   LYS A  20      33.735  16.258  30.773  1.00 41.93           N
ATOM     78  CA  LYS A  20      33.819  15.832  29.374  1.00 42.83           C
ATOM     79  C   LYS A  20      34.745  16.718  28.549  1.00 42.66           C
ATOM     80  O   LYS A  20      35.427  16.243  27.638  1.00 42.05           O
ATOM     81  N   THR A  21      34.805  18.007  28.797  1.00 40.77           N
ATOM     82  CA  THR A  21      35.593  18.896  27.954  1.00 40.50           C
ATOM     83  C   THR A  21      36.842  19.399  28.654  1.00 42.78           C
ATOM     84  O   THR A  21      37.823  19.657  27.983  1.00 44.01           O
ATOM     85  N   LYS A  22      36.815  19.486  29.974  1.00 44.63           N
ATOM     86  CA  LYS A  22      37.951  19.848  30.817  1.00 45.73           C
ATOM     87  C   LYS A  22      38.543  21.187  30.359  1.00 44.09           C
ATOM     88  O   LYS A  22      39.660  21.265  29.818  1.00 48.05           O
TER
END
"""

def run(pdb_str, d_min, b):
  #
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str_box)
  cs = pdb_inp.crystal_symmetry()
  ph = pdb_inp.construct_hierarchy()
  xrs = ph.extract_xray_structure(crystal_symmetry=cs)
  xrs = xrs.set_b_iso(value=b)
  fc = xrs.structure_factors(d_min=d_min).f_calc()
  cg = maptbx.crystal_gridding(
    unit_cell         = cs.unit_cell(),
    space_group_info  = cs.space_group_info(),
    d_min             = d_min)
  fft_map = miller.fft_map(
    crystal_gridding     = cg,
    fourier_coefficients = fc)
  map_obs = fft_map.real_map_unpadded()
  #
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  xrs = ph.extract_xray_structure(crystal_symmetry=cs)
  o = maptbx.resolution_from_map_and_model.run(map_data=map_obs,
    xray_structure=xrs)
  return o.d_min, o.b_iso, o.d_fsc_model

if (__name__ == "__main__"):
  for pdb_str in [pdb_str_5bb, pdb_str_5]:
    for d_min in [2., 3., 4., 5.]:
      for b in [20,80,100,200]:
        if((d_min==2. or d_min==3.) and b>50.): continue
        result, b_result, d_fsc_model = run(pdb_str=pdb_str, d_min=d_min, b=b)
        print b, d_min, "<>", b_result, result, d_fsc_model
        #assert approx_equal(d_min, result)
        #XXX assert cc>0.9
        #XXX if(not randomize):
        #XXX   assert approx_equal(b, b_result,5.+1)
