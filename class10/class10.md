# Class 10: Structural Bioinformatics Pt. 1
Steven Studdard (PID: A19123470)

## Background

The main repository of of high-resolution structural data on
biomolecules is called the **Protein Data Bank** (PDB).

## PDB Statistics

What is in the PDB in terms of molecule type and structure determination
method?

``` r
library(readr)
pdb<- read_csv("Data Export Summary.csv")
```

    Rows: 6 Columns: 9
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: ","
    chr (1): Molecular Type
    dbl (4): Integrative, Multiple methods, Neutron, Other
    num (4): X-ray, EM, NMR, Total

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
head(pdb)
```

    # A tibble: 6 × 9
      `Molecular Type`    `X-ray`    EM   NMR Integrative `Multiple methods` Neutron
      <chr>                 <dbl> <dbl> <dbl>       <dbl>              <dbl>   <dbl>
    1 Protein (only)       180758 23111 12813         348                229      84
    2 Protein/Oligosacch…   10488  3741    34           8                 11       1
    3 Protein/NA             9205  6751   287          26                  8       0
    4 Nucleic acid (only)    3154   250  1578           3                 15       3
    5 Other                   178    27    35           4                  0       0
    6 Oligosaccharide (o…      11     0     6           0                  1       0
    # ℹ 2 more variables: Other <dbl>, Total <dbl>

> Q1: What percentage of structures in the PDB are solved by X-Ray and
> Electron Microscopy.

93.87% of structures in the PDB are solved by X-Ray and Electron
Microscopy.

``` r
x <- sum(pdb$`X-ray`) + sum(pdb$`EM`) 
y <- sum(pdb$Total)
x/y * 100
```

    [1] 93.86623

> Q2: What proportion of structures in the PDB are protein?

The total number of protein sequences in UniProt is 202556314. There are
217375 in the dataset. Leading to a ~0.1%.

``` r
pdb$Total[1]/202556314 * 100
```

    [1] 0.1073158

> **Key-Point**: We have a very, very small structural coverage of known
> proteins (~0.1%). Most structues we know are about (~%80) come from
> one method (X-Ray Crystalography).

> Q3: Type HIV in the PDB website search box on the home page and
> determine how many HIV-1 protease structures are in the current PDB?

After typing HIV into the search box in PDB, I saw that there were 4998
Structures that match the query.

## Visualizing PDB Data with Mol-Star

Main stand alone web version with all features is at
https://molstar.org/viewer/.

![The HIV-Protease Enzyme is a homodimer with two chains.](1HSG.png)

![Surface display showing the binding cleft site for the inhibitor
(drug) molecule.](1HSG_surf.png)

![Spacefill display of catalytic ASP25 amino acids and key binding site
water molecule](1HSG_label.png)

# GEtting started with the Bio3D package

Bio3D is on R package from CRAN for structural bioinformatics

> Q7

> Q8

> Q9

``` r
library(bio3d)

pdb <- read.pdb("1hsg")
```

      Note: Accessing on-line PDB file

``` r
pdb
```


     Call:  read.pdb(file = "1hsg")

       Total Models#: 1
         Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)

         Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
         Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)

         Non-protein/nucleic Atoms#: 172  (residues: 128)
         Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]

       Protein sequence:
          PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
          QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
          ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
          VNIIGRNLLTQIGCTLNF

    + attr: atom, xyz, seqres, helix, sheet,
            calpha, remark, call

``` r
pdb$atom
```

           type eleno elety  alt resid chain resno insert      x       y       z o
    1      ATOM     1     N <NA>   PRO     A     1   <NA> 29.361  39.686   5.862 1
    2      ATOM     2    CA <NA>   PRO     A     1   <NA> 30.307  38.663   5.319 1
    3      ATOM     3     C <NA>   PRO     A     1   <NA> 29.760  38.071   4.022 1
    4      ATOM     4     O <NA>   PRO     A     1   <NA> 28.600  38.302   3.676 1
    5      ATOM     5    CB <NA>   PRO     A     1   <NA> 30.508  37.541   6.342 1
    6      ATOM     6    CG <NA>   PRO     A     1   <NA> 29.296  37.591   7.162 1
    7      ATOM     7    CD <NA>   PRO     A     1   <NA> 28.778  39.015   7.019 1
    8      ATOM     8     N <NA>   GLN     A     2   <NA> 30.607  37.334   3.305 1
    9      ATOM     9    CA <NA>   GLN     A     2   <NA> 30.158  36.492   2.199 1
    10     ATOM    10     C <NA>   GLN     A     2   <NA> 30.298  35.041   2.643 1
    11     ATOM    11     O <NA>   GLN     A     2   <NA> 31.401  34.494   2.763 1
    12     ATOM    12    CB <NA>   GLN     A     2   <NA> 30.970  36.738   0.926 1
    13     ATOM    13    CG <NA>   GLN     A     2   <NA> 30.625  35.783  -0.201 1
    14     ATOM    14    CD <NA>   GLN     A     2   <NA> 31.184  36.217  -1.549 1
    15     ATOM    15   OE1 <NA>   GLN     A     2   <NA> 32.006  35.518  -2.156 1
    16     ATOM    16   NE2 <NA>   GLN     A     2   <NA> 30.684  37.339  -2.061 1
    17     ATOM    17     N <NA>   ILE     A     3   <NA> 29.160  34.436   2.919 1
    18     ATOM    18    CA <NA>   ILE     A     3   <NA> 29.123  33.098   3.397 1
    19     ATOM    19     C <NA>   ILE     A     3   <NA> 28.968  32.155   2.198 1
    20     ATOM    20     O <NA>   ILE     A     3   <NA> 28.088  32.330   1.368 1
    21     ATOM    21    CB <NA>   ILE     A     3   <NA> 27.977  32.995   4.409 1
    22     ATOM    22   CG1 <NA>   ILE     A     3   <NA> 28.341  33.820   5.652 1
    23     ATOM    23   CG2 <NA>   ILE     A     3   <NA> 27.692  31.548   4.745 1
    24     ATOM    24   CD1 <NA>   ILE     A     3   <NA> 27.264  33.884   6.696 1
    25     ATOM    25     N <NA>   THR     A     4   <NA> 29.891  31.210   2.066 1
    26     ATOM    26    CA <NA>   THR     A     4   <NA> 29.774  30.143   1.062 1
    27     ATOM    27     C <NA>   THR     A     4   <NA> 28.986  28.975   1.658 1
    28     ATOM    28     O <NA>   THR     A     4   <NA> 28.690  28.948   2.875 1
    29     ATOM    29    CB <NA>   THR     A     4   <NA> 31.165  29.618   0.634 1
    30     ATOM    30   OG1 <NA>   THR     A     4   <NA> 31.866  29.209   1.815 1
    31     ATOM    31   CG2 <NA>   THR     A     4   <NA> 31.980  30.688  -0.085 1
    32     ATOM    32     N <NA>   LEU     A     5   <NA> 28.641  28.019   0.803 1
    33     ATOM    33    CA <NA>   LEU     A     5   <NA> 27.644  27.003   1.144 1
    34     ATOM    34     C <NA>   LEU     A     5   <NA> 28.204  25.559   1.071 1
    35     ATOM    35     O <NA>   LEU     A     5   <NA> 27.446  24.583   0.969 1
    36     ATOM    36    CB <NA>   LEU     A     5   <NA> 26.411  27.139   0.226 1
    37     ATOM    37    CG <NA>   LEU     A     5   <NA> 25.676  28.479   0.352 1
    38     ATOM    38   CD1 <NA>   LEU     A     5   <NA> 24.624  28.624  -0.753 1
    39     ATOM    39   CD2 <NA>   LEU     A     5   <NA> 25.088  28.590   1.745 1
    40     ATOM    40     N <NA>   TRP     A     6   <NA> 29.528  25.436   1.146 1
    41     ATOM    41    CA <NA>   TRP     A     6   <NA> 30.177  24.150   1.279 1
    42     ATOM    42     C <NA>   TRP     A     6   <NA> 29.837  23.488   2.611 1
    43     ATOM    43     O <NA>   TRP     A     6   <NA> 29.706  22.271   2.673 1
    44     ATOM    44    CB <NA>   TRP     A     6   <NA> 31.685  24.301   1.109 1
    45     ATOM    45    CG <NA>   TRP     A     6   <NA> 32.152  24.955  -0.189 1
    46     ATOM    46   CD1 <NA>   TRP     A     6   <NA> 32.681  26.216  -0.345 1
    47     ATOM    47   CD2 <NA>   TRP     A     6   <NA> 32.274  24.314  -1.478 1
    48     ATOM    48   NE1 <NA>   TRP     A     6   <NA> 33.102  26.385  -1.655 1
    49     ATOM    49   CE2 <NA>   TRP     A     6   <NA> 32.864  25.258  -2.369 1
    50     ATOM    50   CE3 <NA>   TRP     A     6   <NA> 31.949  23.035  -1.986 1
    51     ATOM    51   CZ2 <NA>   TRP     A     6   <NA> 33.093  24.968  -3.717 1
    52     ATOM    52   CZ3 <NA>   TRP     A     6   <NA> 32.195  22.755  -3.294 1
    53     ATOM    53   CH2 <NA>   TRP     A     6   <NA> 32.754  23.722  -4.169 1
    54     ATOM    54     N <NA>   GLN     A     7   <NA> 29.667  24.280   3.667 1
    55     ATOM    55    CA <NA>   GLN     A     7   <NA> 29.141  23.799   4.960 1
    56     ATOM    56     C <NA>   GLN     A     7   <NA> 27.747  24.395   5.208 1
    57     ATOM    57     O <NA>   GLN     A     7   <NA> 27.349  25.330   4.547 1
    58     ATOM    58    CB <NA>   GLN     A     7   <NA> 30.072  24.227   6.100 1
    59     ATOM    59    CG <NA>   GLN     A     7   <NA> 31.512  23.694   5.995 1
    60     ATOM    60    CD <NA>   GLN     A     7   <NA> 32.521  24.750   5.469 1
    61     ATOM    61   OE1 <NA>   GLN     A     7   <NA> 32.666  25.860   6.038 1
    62     ATOM    62   NE2 <NA>   GLN     A     7   <NA> 33.268  24.374   4.419 1
    63     ATOM    63     N <NA>   ARG     A     8   <NA> 26.992  23.877   6.169 1
    64     ATOM    64    CA <NA>   ARG     A     8   <NA> 25.757  24.566   6.593 1
    65     ATOM    65     C <NA>   ARG     A     8   <NA> 26.029  26.025   6.996 1
    66     ATOM    66     O <NA>   ARG     A     8   <NA> 26.947  26.291   7.775 1
    67     ATOM    67    CB <NA>   ARG     A     8   <NA> 25.087  23.849   7.776 1
    68     ATOM    68    CG <NA>   ARG     A     8   <NA> 24.646  22.409   7.505 1
    69     ATOM    69    CD <NA>   ARG     A     8   <NA> 23.728  21.896   8.637 1
    70     ATOM    70    NE <NA>   ARG     A     8   <NA> 22.952  20.730   8.230 1
    71     ATOM    71    CZ <NA>   ARG     A     8   <NA> 22.367  19.871   9.064 1
    72     ATOM    72   NH1 <NA>   ARG     A     8   <NA> 22.376  20.074  10.370 1
    73     ATOM    73   NH2 <NA>   ARG     A     8   <NA> 21.776  18.789   8.589 1
    74     ATOM    74     N <NA>   PRO     A     9   <NA> 25.123  26.955   6.645 1
    75     ATOM    75    CA <NA>   PRO     A     9   <NA> 25.491  28.352   6.938 1
    76     ATOM    76     C <NA>   PRO     A     9   <NA> 25.127  28.763   8.364 1
    77     ATOM    77     O <NA>   PRO     A     9   <NA> 24.136  29.472   8.578 1
    78     ATOM    78    CB <NA>   PRO     A     9   <NA> 24.719  29.176   5.916 1
    79     ATOM    79    CG <NA>   PRO     A     9   <NA> 23.625  28.254   5.407 1
    80     ATOM    80    CD <NA>   PRO     A     9   <NA> 24.096  26.855   5.591 1
    81     ATOM    81     N <NA>   LEU     A    10   <NA> 25.905  28.285   9.330 1
    82     ATOM    82    CA <NA>   LEU     A    10   <NA> 25.653  28.510  10.750 1
    83     ATOM    83     C <NA>   LEU     A    10   <NA> 26.383  29.770  11.208 1
    84     ATOM    84     O <NA>   LEU     A    10   <NA> 27.567  29.927  10.938 1
    85     ATOM    85    CB <NA>   LEU     A    10   <NA> 26.120  27.284  11.573 1
    86     ATOM    86    CG <NA>   LEU     A    10   <NA> 25.161  26.082  11.544 1
    87     ATOM    87   CD1 <NA>   LEU     A    10   <NA> 25.895  24.743  11.662 1
    88     ATOM    88   CD2 <NA>   LEU     A    10   <NA> 24.206  26.196  12.696 1
    89     ATOM    89     N <NA>   VAL     A    11   <NA> 25.667  30.672  11.872 1
    90     ATOM    90    CA <NA>   VAL     A    11   <NA> 26.267  31.854  12.497 1
    91     ATOM    91     C <NA>   VAL     A    11   <NA> 25.818  31.957  13.955 1
    92     ATOM    92     O <NA>   VAL     A    11   <NA> 24.929  31.184  14.402 1
    93     ATOM    93    CB <NA>   VAL     A    11   <NA> 25.824  33.131  11.791 1
    94     ATOM    94   CG1 <NA>   VAL     A    11   <NA> 26.270  33.089  10.323 1
    95     ATOM    95   CG2 <NA>   VAL     A    11   <NA> 24.333  33.275  11.879 1
    96     ATOM    96     N <NA>   THR     A    12   <NA> 26.397  32.913  14.700 1
    97     ATOM    97    CA <NA>   THR     A    12   <NA> 26.001  33.143  16.102 1
    98     ATOM    98     C <NA>   THR     A    12   <NA> 24.915  34.200  16.204 1
    99     ATOM    99     O <NA>   THR     A    12   <NA> 25.010  35.279  15.610 1
    100    ATOM   100    CB <NA>   THR     A    12   <NA> 27.201  33.565  16.998 1
    101    ATOM   101   OG1 <NA>   THR     A    12   <NA> 28.330  32.709  16.751 1
    102    ATOM   102   CG2 <NA>   THR     A    12   <NA> 26.827  33.430  18.450 1
    103    ATOM   103     N <NA>   ILE     A    13   <NA> 23.848  33.868  16.909 1
    104    ATOM   104    CA <NA>   ILE     A    13   <NA> 22.842  34.875  17.206 1
    105    ATOM   105     C <NA>   ILE     A    13   <NA> 22.770  35.114  18.707 1
    106    ATOM   106     O <NA>   ILE     A    13   <NA> 23.328  34.363  19.500 1
    107    ATOM   107    CB <NA>   ILE     A    13   <NA> 21.413  34.460  16.661 1
    108    ATOM   108   CG1 <NA>   ILE     A    13   <NA> 20.878  33.229  17.431 1
    109    ATOM   109   CG2 <NA>   ILE     A    13   <NA> 21.510  34.194  15.162 1
    110    ATOM   110   CD1 <NA>   ILE     A    13   <NA> 19.353  33.201  17.603 1
    111    ATOM   111     N <NA>   LYS     A    14   <NA> 22.106  36.199  19.087 1
    112    ATOM   112    CA <NA>   LYS     A    14   <NA> 21.894  36.545  20.492 1
    113    ATOM   113     C <NA>   LYS     A    14   <NA> 20.442  36.943  20.615 1
    114    ATOM   114     O <NA>   LYS     A    14   <NA> 19.960  37.808  19.873 1
    115    ATOM   115    CB <NA>   LYS     A    14   <NA> 22.777  37.724  20.896 1
    116    ATOM   116    CG <NA>   LYS     A    14   <NA> 22.727  38.056  22.383 1
    117    ATOM   117    CD <NA>   LYS     A    14   <NA> 23.270  39.450  22.678 1
    118    ATOM   118    CE <NA>   LYS     A    14   <NA> 24.814  39.490  22.755 1
    119    ATOM   119    NZ <NA>   LYS     A    14   <NA> 25.394  40.891  22.572 1
    120    ATOM   120     N <NA>   ILE     A    15   <NA> 19.739  36.267  21.512 1
    121    ATOM   121    CA <NA>   ILE     A    15   <NA> 18.345  36.563  21.813 1
    122    ATOM   122     C <NA>   ILE     A    15   <NA> 18.224  36.327  23.316 1
    123    ATOM   123     O <NA>   ILE     A    15   <NA> 18.886  35.449  23.864 1
    124    ATOM   124    CB <NA>   ILE     A    15   <NA> 17.380  35.592  21.022 1
    125    ATOM   125   CG1 <NA>   ILE     A    15   <NA> 15.935  35.812  21.435 1
    126    ATOM   126   CG2 <NA>   ILE     A    15   <NA> 17.745  34.137  21.266 1
    127    ATOM   127   CD1 <NA>   ILE     A    15   <NA> 14.929  35.116  20.526 1
    128    ATOM   128     N <NA>   GLY     A    16   <NA> 17.446  37.139  24.012 1
    129    ATOM   129    CA <NA>   GLY     A    16   <NA> 17.356  36.968  25.459 1
    130    ATOM   130     C <NA>   GLY     A    16   <NA> 18.711  36.871  26.160 1
    131    ATOM   131     O <NA>   GLY     A    16   <NA> 18.866  36.162  27.153 1
    132    ATOM   132     N <NA>   GLY     A    17   <NA> 19.671  37.659  25.697 1
    133    ATOM   133    CA <NA>   GLY     A    17   <NA> 20.970  37.660  26.340 1
    134    ATOM   134     C <NA>   GLY     A    17   <NA> 21.680  36.316  26.278 1
    135    ATOM   135     O <NA>   GLY     A    17   <NA> 22.785  36.163  26.794 1
    136    ATOM   136     N <NA>   GLN     A    18   <NA> 21.093  35.361  25.572 1
    137    ATOM   137    CA <NA>   GLN     A    18   <NA> 21.780  34.106  25.263 1
    138    ATOM   138     C <NA>   GLN     A    18   <NA> 22.500  34.159  23.907 1
    139    ATOM   139     O <NA>   GLN     A    18   <NA> 21.937  34.624  22.915 1
    140    ATOM   140    CB <NA>   GLN     A    18   <NA> 20.776  32.957  25.228 1
    141    ATOM   141    CG <NA>   GLN     A    18   <NA> 19.599  33.116  26.176 1
    142    ATOM   142    CD <NA>   GLN     A    18   <NA> 19.556  31.997  27.179 1
    143    ATOM   143   OE1 <NA>   GLN     A    18   <NA> 20.393  31.944  28.082 1
    144    ATOM   144   NE2 <NA>   GLN     A    18   <NA> 18.647  31.035  26.975 1
    145    ATOM   145     N <NA>   LEU     A    19   <NA> 23.733  33.672  23.848 1
    146    ATOM   146    CA <NA>   LEU     A    19   <NA> 24.334  33.365  22.552 1
    147    ATOM   147     C <NA>   LEU     A    19   <NA> 23.896  31.963  22.106 1
    148    ATOM   148     O <NA>   LEU     A    19   <NA> 23.975  31.020  22.863 1
    149    ATOM   149    CB <NA>   LEU     A    19   <NA> 25.869  33.432  22.625 1
    150    ATOM   150    CG <NA>   LEU     A    19   <NA> 26.561  34.761  22.968 1
    151    ATOM   151   CD1 <NA>   LEU     A    19   <NA> 28.007  34.629  22.620 1
    152    ATOM   152   CD2 <NA>   LEU     A    19   <NA> 25.983  35.913  22.194 1
    153    ATOM   153     N <NA>   LYS     A    20   <NA> 23.416  31.855  20.876 1
    154    ATOM   154    CA <NA>   LYS     A    20   <NA> 23.006  30.584  20.266 1
    155    ATOM   155     C <NA>   LYS     A    20   <NA> 23.626  30.463  18.874 1
    156    ATOM   156     O <NA>   LYS     A    20   <NA> 24.024  31.460  18.283 1
    157    ATOM   157    CB <NA>   LYS     A    20   <NA> 21.494  30.523  20.107 1
    158    ATOM   158    CG <NA>   LYS     A    20   <NA> 20.778  29.875  21.264 1
    159    ATOM   159    CD <NA>   LYS     A    20   <NA> 19.868  30.857  21.939 1
    160    ATOM   160    CE <NA>   LYS     A    20   <NA> 19.112  30.168  23.043 1
    161    ATOM   161    NZ <NA>   LYS     A    20   <NA> 18.467  28.892  22.571 1
    162    ATOM   162     N <NA>   GLU     A    21   <NA> 23.725  29.250  18.342 1
    163    ATOM   163    CA <NA>   GLU     A    21   <NA> 24.053  29.117  16.931 1
    164    ATOM   164     C <NA>   GLU     A    21   <NA> 22.822  28.761  16.150 1
    165    ATOM   165     O <NA>   GLU     A    21   <NA> 21.879  28.136  16.672 1
    166    ATOM   166    CB <NA>   GLU     A    21   <NA> 25.197  28.130  16.679 1
    167    ATOM   167    CG <NA>   GLU     A    21   <NA> 25.035  26.716  17.168 1
    168    ATOM   168    CD <NA>   GLU     A    21   <NA> 25.878  25.743  16.334 1
    169    ATOM   169   OE1 <NA>   GLU     A    21   <NA> 27.022  26.130  15.972 1
    170    ATOM   170   OE2 <NA>   GLU     A    21   <NA> 25.379  24.639  15.983 1
    171    ATOM   171     N <NA>   ALA     A    22   <NA> 22.778  29.268  14.927 1
    172    ATOM   172    CA <NA>   ALA     A    22   <NA> 21.553  29.189  14.165 1
    173    ATOM   173     C <NA>   ALA     A    22   <NA> 21.870  29.183  12.682 1
    174    ATOM   174     O <NA>   ALA     A    22   <NA> 22.975  29.578  12.252 1
    175    ATOM   175    CB <NA>   ALA     A    22   <NA> 20.625  30.359  14.524 1
    176    ATOM   176     N <NA>   LEU     A    23   <NA> 20.893  28.726  11.903 1
    177    ATOM   177    CA <NA>   LEU     A    23   <NA> 21.047  28.473  10.476 1
    178    ATOM   178     C <NA>   LEU     A    23   <NA> 20.381  29.596   9.664 1
    179    ATOM   179     O <NA>   LEU     A    23   <NA> 19.231  29.943   9.912 1
    180    ATOM   180    CB <NA>   LEU     A    23   <NA> 20.382  27.135  10.174 1
    181    ATOM   181    CG <NA>   LEU     A    23   <NA> 20.532  26.573   8.786 1
    182    ATOM   182   CD1 <NA>   LEU     A    23   <NA> 21.939  26.039   8.621 1
    183    ATOM   183   CD2 <NA>   LEU     A    23   <NA> 19.490  25.490   8.627 1
    184    ATOM   184     N <NA>   LEU     A    24   <NA> 21.122  30.163   8.715 1
    185    ATOM   185    CA <NA>   LEU     A    24   <NA> 20.617  31.144   7.775 1
    186    ATOM   186     C <NA>   LEU     A    24   <NA> 19.940  30.412   6.617 1
    187    ATOM   187     O <NA>   LEU     A    24   <NA> 20.567  29.833   5.740 1
    188    ATOM   188    CB <NA>   LEU     A    24   <NA> 21.767  32.023   7.262 1
    189    ATOM   189    CG <NA>   LEU     A    24   <NA> 22.647  32.673   8.359 1
    190    ATOM   190   CD1 <NA>   LEU     A    24   <NA> 23.698  33.581   7.738 1
    191    ATOM   191   CD2 <NA>   LEU     A    24   <NA> 21.797  33.496   9.368 1
    192    ATOM   192     N <NA>   ASP     A    25   <NA> 18.626  30.444   6.627 1
    193    ATOM   193    CA <NA>   ASP     A    25   <NA> 17.853  29.516   5.837 1
    194    ATOM   194     C <NA>   ASP     A    25   <NA> 16.945  30.292   4.886 1
    195    ATOM   195     O <NA>   ASP     A    25   <NA> 15.843  30.678   5.237 1
    196    ATOM   196    CB <NA>   ASP     A    25   <NA> 17.047  28.642   6.811 1
    197    ATOM   197    CG <NA>   ASP     A    25   <NA> 16.316  27.513   6.146 1
    198    ATOM   198   OD1 <NA>   ASP     A    25   <NA> 16.236  27.458   4.905 1
    199    ATOM   199   OD2 <NA>   ASP     A    25   <NA> 15.762  26.696   6.882 1
    200    ATOM   200     N <NA>   THR     A    26   <NA> 17.364  30.439   3.645 1
    201    ATOM   201    CA <NA>   THR     A    26   <NA> 16.548  31.148   2.684 1
    202    ATOM   202     C <NA>   THR     A    26   <NA> 15.302  30.382   2.289 1
    203    ATOM   203     O <NA>   THR     A    26   <NA> 14.412  30.939   1.615 1
    204    ATOM   204    CB <NA>   THR     A    26   <NA> 17.328  31.419   1.447 1
    205    ATOM   205   OG1 <NA>   THR     A    26   <NA> 17.693  30.177   0.863 1
    206    ATOM   206   CG2 <NA>   THR     A    26   <NA> 18.601  32.197   1.773 1
    207    ATOM   207     N <NA>   GLY     A    27   <NA> 15.213  29.111   2.702 1
    208    ATOM   208    CA <NA>   GLY     A    27   <NA> 14.043  28.331   2.349 1
    209    ATOM   209     C <NA>   GLY     A    27   <NA> 12.958  28.456   3.394 1
    210    ATOM   210     O <NA>   GLY     A    27   <NA> 11.832  28.015   3.171 1
    211    ATOM   211     N <NA>   ALA     A    28   <NA> 13.301  28.967   4.569 1
    212    ATOM   212    CA <NA>   ALA     A    28   <NA> 12.324  29.192   5.642 1
    213    ATOM   213     C <NA>   ALA     A    28   <NA> 11.629  30.567   5.554 1
    214    ATOM   214     O <NA>   ALA     A    28   <NA> 12.303  31.618   5.504 1
    215    ATOM   215    CB <NA>   ALA     A    28   <NA> 13.031  29.084   6.978 1
    216    ATOM   216     N <NA>   ASP     A    29   <NA> 10.296  30.560   5.587 1
    217    ATOM   217    CA <NA>   ASP     A    29   <NA>  9.512  31.798   5.651 1
    218    ATOM   218     C <NA>   ASP     A    29   <NA>  9.632  32.481   6.994 1
    219    ATOM   219     O <NA>   ASP     A    29   <NA>  9.671  33.706   7.056 1
    220    ATOM   220    CB <NA>   ASP     A    29   <NA>  8.029  31.534   5.402 1
    221    ATOM   221    CG <NA>   ASP     A    29   <NA>  7.752  31.004   4.015 1
    222    ATOM   222   OD1 <NA>   ASP     A    29   <NA>  8.591  31.192   3.104 1
    223    ATOM   223   OD2 <NA>   ASP     A    29   <NA>  6.661  30.410   3.833 1
    224    ATOM   224     N <NA>   ASP     A    30   <NA>  9.698  31.685   8.062 1
    225    ATOM   225    CA <NA>   ASP     A    30   <NA>  9.718  32.168   9.444 1
    226    ATOM   226     C <NA>   ASP     A    30   <NA> 10.988  31.819  10.163 1
    227    ATOM   227     O <NA>   ASP     A    30   <NA> 11.818  31.072   9.679 1
    228    ATOM   228    CB <NA>   ASP     A    30   <NA>  8.549  31.585  10.214 1
    229    ATOM   229    CG <NA>   ASP     A    30   <NA>  7.254  31.916   9.579 1
    230    ATOM   230   OD1 <NA>   ASP     A    30   <NA>  6.951  33.118   9.473 1
    231    ATOM   231   OD2 <NA>   ASP     A    30   <NA>  6.561  31.008   9.099 1
    232    ATOM   232     N <NA>   THR     A    31   <NA> 11.161  32.408  11.326 1
    233    ATOM   233    CA <NA>   THR     A    31   <NA> 12.248  32.053  12.215 1
    234    ATOM   234     C <NA>   THR     A    31   <NA> 11.707  31.128  13.318 1
    235    ATOM   235     O <NA>   THR     A    31   <NA> 10.660  31.408  13.910 1
    236    ATOM   236    CB <NA>   THR     A    31   <NA> 12.896  33.338  12.795 1
    237    ATOM   237   OG1 <NA>   THR     A    31   <NA> 13.451  34.082  11.707 1
    238    ATOM   238   CG2 <NA>   THR     A    31   <NA> 14.027  32.992  13.816 1
    239    ATOM   239     N <NA>   VAL     A    32   <NA> 12.390  30.005  13.537 1
    240    ATOM   240    CA <NA>   VAL     A    32   <NA> 11.893  28.983  14.419 1
    241    ATOM   241     C <NA>   VAL     A    32   <NA> 13.036  28.655  15.292 1
    242    ATOM   242     O <NA>   VAL     A    32   <NA> 14.067  28.221  14.821 1
    243    ATOM   243    CB <NA>   VAL     A    32   <NA> 11.528  27.683  13.690 1
    244    ATOM   244   CG1 <NA>   VAL     A    32   <NA> 10.656  26.825  14.592 1
    245    ATOM   245   CG2 <NA>   VAL     A    32   <NA> 10.805  27.963  12.423 1
    246    ATOM   246     N <NA>   LEU     A    33   <NA> 12.899  28.904  16.576 1
    247    ATOM   247    CA <NA>   LEU     A    33   <NA> 13.996  28.594  17.500 1
    248    ATOM   248     C <NA>   LEU     A    33   <NA> 13.571  27.454  18.396 1
    249    ATOM   249     O <NA>   LEU     A    33   <NA> 12.363  27.234  18.612 1
    250    ATOM   250    CB <NA>   LEU     A    33   <NA> 14.337  29.799  18.375 1
    251    ATOM   251    CG <NA>   LEU     A    33   <NA> 14.849  31.061  17.691 1
    252    ATOM   252   CD1 <NA>   LEU     A    33   <NA> 15.091  32.156  18.733 1
    253    ATOM   253   CD2 <NA>   LEU     A    33   <NA> 16.139  30.718  16.927 1
    254    ATOM   254     N <NA>   GLU     A    34   <NA> 14.568  26.722  18.889 1
    255    ATOM   255    CA <NA>   GLU     A    34   <NA> 14.390  25.650  19.874 1
    256    ATOM   256     C <NA>   GLU     A    34   <NA> 13.739  26.125  21.160 1
    257    ATOM   257     O <NA>   GLU     A    34   <NA> 13.982  27.222  21.618 1
    258    ATOM   258    CB <NA>   GLU     A    34   <NA> 15.727  25.007  20.190 1
    259    ATOM   259    CG <NA>   GLU     A    34   <NA> 16.297  24.250  18.988 1
    260    ATOM   260    CD <NA>   GLU     A    34   <NA> 17.726  23.801  19.191 1
    261    ATOM   261   OE1 <NA>   GLU     A    34   <NA> 18.134  23.657  20.375 1
    262    ATOM   262   OE2 <NA>   GLU     A    34   <NA> 18.443  23.614  18.182 1
    263    ATOM   263     N <NA>   GLU     A    35   <NA> 12.865  25.288  21.703 1
    264    ATOM   264    CA <NA>   GLU     A    35   <NA> 12.183  25.482  22.981 1
    265    ATOM   265     C <NA>   GLU     A    35   <NA> 12.971  26.329  23.986 1
    266    ATOM   266     O <NA>   GLU     A    35   <NA> 13.981  25.861  24.497 1
    267    ATOM   267    CB <NA>   GLU     A    35   <NA> 11.941  24.114  23.581 1
    268    ATOM   268    CG <NA>   GLU     A    35   <NA> 10.800  24.049  24.516 1
    269    ATOM   269    CD <NA>   GLU     A    35   <NA>  9.489  24.067  23.809 1
    270    ATOM   270   OE1 <NA>   GLU     A    35   <NA>  9.134  23.066  23.120 1
    271    ATOM   271   OE2 <NA>   GLU     A    35   <NA>  8.758  25.047  24.035 1
    272    ATOM   272     N <NA>   MET     A    36   <NA> 12.495  27.556  24.269 1
    273    ATOM   273    CA <NA>   MET     A    36   <NA> 13.101  28.479  25.261 1
    274    ATOM   274     C <NA>   MET     A    36   <NA> 12.013  29.437  25.755 1
    275    ATOM   275     O <NA>   MET     A    36   <NA> 10.969  29.535  25.133 1
    276    ATOM   276    CB <NA>   MET     A    36   <NA> 14.216  29.306  24.621 1
    277    ATOM   277    CG <NA>   MET     A    36   <NA> 13.741  30.341  23.590 1
    278    ATOM   278    SD <NA>   MET     A    36   <NA> 15.123  31.051  22.667 1
    279    ATOM   279    CE <NA>   MET     A    36   <NA> 15.783  32.224  23.828 1
    280    ATOM   280     N <NA>   SER     A    37   <NA> 12.249  30.179  26.852 1
    281    ATOM   281    CA <NA>   SER     A    37   <NA> 11.163  31.036  27.460 1
    282    ATOM   282     C <NA>   SER     A    37   <NA> 11.347  32.442  26.960 1
    283    ATOM   283     O <NA>   SER     A    37   <NA> 12.475  32.929  27.007 1
    284    ATOM   284    CB <NA>   SER     A    37   <NA> 11.300  31.071  28.976 1
    285    ATOM   285    OG <NA>   SER     A    37   <NA>  9.948  31.148  29.559 1
    286    ATOM   286     N <NA>   LEU     A    38   <NA> 10.330  33.112  26.459 1
    287    ATOM   287    CA <NA>   LEU     A    38   <NA> 10.566  34.476  26.033 1
    288    ATOM   288     C <NA>   LEU     A    38   <NA>  9.594  35.375  26.756 1
    289    ATOM   289     O <NA>   LEU     A    38   <NA>  8.616  34.906  27.343 1
    290    ATOM   290    CB <NA>   LEU     A    38   <NA> 10.409  34.626  24.500 1
    291    ATOM   291    CG <NA>   LEU     A    38   <NA> 11.559  34.187  23.577 1
    292    ATOM   292   CD1 <NA>   LEU     A    38   <NA> 11.171  34.399  22.132 1
    293    ATOM   293   CD2 <NA>   LEU     A    38   <NA> 12.807  34.964  23.875 1
    294    ATOM   294     N <NA>   PRO     A    39   <NA>  9.929  36.666  26.880 1
    295    ATOM   295    CA <NA>   PRO     A    39   <NA>  8.980  37.700  27.301 1
    296    ATOM   296     C <NA>   PRO     A    39   <NA>  7.760  37.785  26.410 1
    297    ATOM   297     O <NA>   PRO     A    39   <NA>  7.866  37.883  25.194 1
    298    ATOM   298    CB <NA>   PRO     A    39   <NA>  9.778  38.989  27.220 1
    299    ATOM   299    CG <NA>   PRO     A    39   <NA> 11.021  38.637  26.370 1
    300    ATOM   300    CD <NA>   PRO     A    39   <NA> 11.291  37.226  26.725 1
    301    ATOM   301     N <NA>   GLY     A    40   <NA>  6.601  37.811  27.029 1
    302    ATOM   302    CA <NA>   GLY     A    40   <NA>  5.419  38.218  26.302 1
    303    ATOM   303     C <NA>   GLY     A    40   <NA>  4.430  37.094  26.331 1
    304    ATOM   304     O <NA>   GLY     A    40   <NA>  4.591  36.107  27.055 1
    305    ATOM   305     N <NA>   ARG     A    41   <NA>  3.289  37.341  25.729 1
    306    ATOM   306    CA <NA>   ARG     A    41   <NA>  2.382  36.252  25.509 1
    307    ATOM   307     C <NA>   ARG     A    41   <NA>  2.606  35.713  24.096 1
    308    ATOM   308     O <NA>   ARG     A    41   <NA>  3.225  36.383  23.273 1
    309    ATOM   309    CB <NA>   ARG     A    41   <NA>  0.956  36.719  25.748 1
    310    ATOM   310    CG <NA>   ARG     A    41   <NA>  0.288  36.021  26.959 1
    311    ATOM   311    CD <NA>   ARG     A    41   <NA>  0.118  36.953  28.169 1
    312    ATOM   312    NE <NA>   ARG     A    41   <NA>  1.356  37.143  28.933 1
    313    ATOM   313    CZ <NA>   ARG     A    41   <NA>  1.830  38.332  29.318 1
    314    ATOM   314   NH1 <NA>   ARG     A    41   <NA>  1.161  39.450  29.040 1
    315    ATOM   315   NH2 <NA>   ARG     A    41   <NA>  2.973  38.402  29.994 1
    316    ATOM   316     N <NA>   TRP     A    42   <NA>  2.145  34.493  23.834 1
    317    ATOM   317    CA <NA>   TRP     A    42   <NA>  2.295  33.862  22.533 1
    318    ATOM   318     C <NA>   TRP     A    42   <NA>  0.934  33.419  21.959 1
    319    ATOM   319     O <NA>   TRP     A    42   <NA>  0.031  33.049  22.694 1
    320    ATOM   320    CB <NA>   TRP     A    42   <NA>  3.207  32.645  22.642 1
    321    ATOM   321    CG <NA>   TRP     A    42   <NA>  2.946  31.783  23.787 1
    322    ATOM   322   CD1 <NA>   TRP     A    42   <NA>  3.473  31.896  25.041 1
    323    ATOM   323   CD2 <NA>   TRP     A    42   <NA>  1.989  30.710  23.857 1
    324    ATOM   324   NE1 <NA>   TRP     A    42   <NA>  2.882  30.973  25.884 1
    325    ATOM   325   CE2 <NA>   TRP     A    42   <NA>  1.966  30.246  25.193 1
    326    ATOM   326   CE3 <NA>   TRP     A    42   <NA>  1.129  30.108  22.909 1
    327    ATOM   327   CZ2 <NA>   TRP     A    42   <NA>  1.117  29.219  25.618 1
    328    ATOM   328   CZ3 <NA>   TRP     A    42   <NA>  0.313  29.091  23.334 1
    329    ATOM   329   CH2 <NA>   TRP     A    42   <NA>  0.305  28.651  24.686 1
    330    ATOM   330     N <NA>   LYS     A    43   <NA>  0.781  33.465  20.639 1
    331    ATOM   331    CA <NA>   LYS     A    43   <NA> -0.305  32.776  19.928 1
    332    ATOM   332     C <NA>   LYS     A    43   <NA>  0.220  31.412  19.477 1
    333    ATOM   333     O <NA>   LYS     A    43   <NA>  1.400  31.270  19.145 1
    334    ATOM   334    CB <NA>   LYS     A    43   <NA> -0.739  33.603  18.700 1
    335    ATOM   335    CG <NA>   LYS     A    43   <NA> -1.311  34.967  19.027 1
    336    ATOM   336    CD <NA>   LYS     A    43   <NA> -1.066  35.945  17.886 1
    337    ATOM   337    CE <NA>   LYS     A    43   <NA> -1.726  37.319  18.156 1
    338    ATOM   338    NZ <NA>   LYS     A    43   <NA> -0.979  38.292  19.067 1
    339    ATOM   339     N <NA>   PRO     A    44   <NA> -0.601  30.363  19.541 1
    340    ATOM   340    CA <NA>   PRO     A    44   <NA> -0.088  29.114  18.966 1
    341    ATOM   341     C <NA>   PRO     A    44   <NA> -0.275  29.085  17.454 1
    342    ATOM   342     O <NA>   PRO     A    44   <NA> -1.185  29.721  16.929 1
    343    ATOM   343    CB <NA>   PRO     A    44   <NA> -0.893  28.021  19.667 1
    344    ATOM   344    CG <NA>   PRO     A    44   <NA> -2.170  28.683  20.012 1
    345    ATOM   345    CD <NA>   PRO     A    44   <NA> -1.825  30.129  20.325 1
    346    ATOM   346     N <NA>   LYS     A    45   <NA>  0.586  28.336  16.762 1
    347    ATOM   347    CA <NA>   LYS     A    45   <NA>  0.634  28.302  15.290 1
    348    ATOM   348     C <NA>   LYS     A    45   <NA>  1.025  26.869  14.873 1
    349    ATOM   349     O <NA>   LYS     A    45   <NA>  1.711  26.163  15.608 1
    350    ATOM   350    CB <NA>   LYS     A    45   <NA>  1.693  29.299  14.799 1
    351    ATOM   351    CG <NA>   LYS     A    45   <NA>  1.495  29.822  13.396 1
    352    ATOM   352    CD <NA>   LYS     A    45   <NA>  2.628  30.791  12.986 1
    353    ATOM   353    CE <NA>   LYS     A    45   <NA>  2.662  31.103  11.471 1
    354    ATOM   354    NZ <NA>   LYS     A    45   <NA>  1.837  32.301  11.181 1
    355    ATOM   355     N <NA>   MET     A    46   <NA>  0.570  26.438  13.707 1
    356    ATOM   356    CA <NA>   MET     A    46   <NA>  1.091  25.232  13.081 1
    357    ATOM   357     C <NA>   MET     A    46   <NA>  1.899  25.638  11.874 1
    358    ATOM   358     O <NA>   MET     A    46   <NA>  1.385  26.292  10.965 1
    359    ATOM   359    CB <NA>   MET     A    46   <NA> -0.047  24.319  12.624 1
    360    ATOM   360    CG <NA>   MET     A    46   <NA> -0.970  23.867  13.745 1
    361    ATOM   361    SD <NA>   MET     A    46   <NA> -0.348  22.459  14.681 1
    362    ATOM   362    CE <NA>   MET     A    46   <NA> -0.585  21.081  13.426 1
    363    ATOM   363     N <NA>   ILE     A    47   <NA>  3.167  25.258  11.849 1
    364    ATOM   364    CA <NA>   ILE     A    47   <NA>  3.978  25.464  10.655 1
    365    ATOM   365     C <NA>   ILE     A    47   <NA>  4.420  24.140  10.114 1
    366    ATOM   366     O <NA>   ILE     A    47   <NA>  4.667  23.248  10.887 1
    367    ATOM   367    CB <NA>   ILE     A    47   <NA>  5.234  26.259  10.953 1
    368    ATOM   368   CG1 <NA>   ILE     A    47   <NA>  5.959  25.628  12.127 1
    369    ATOM   369   CG2 <NA>   ILE     A    47   <NA>  4.898  27.703  11.148 1
    370    ATOM   370   CD1 <NA>   ILE     A    47   <NA>  7.369  26.170  12.291 1
    371    ATOM   371     N <NA>   GLY     A    48   <NA>  4.567  24.042   8.795 1
    372    ATOM   372    CA <NA>   GLY     A    48   <NA>  4.886  22.777   8.145 1
    373    ATOM   373     C <NA>   GLY     A    48   <NA>  6.265  22.720   7.521 1
    374    ATOM   374     O <NA>   GLY     A    48   <NA>  6.723  23.689   6.900 1
    375    ATOM   375     N <NA>   GLY     A    49   <NA>  7.015  21.688   7.884 1
    376    ATOM   376    CA <NA>   GLY     A    49   <NA>  8.313  21.471   7.276 1
    377    ATOM   377     C <NA>   GLY     A    49   <NA>  8.221  20.379   6.232 1
    378    ATOM   378     O <NA>   GLY     A    49   <NA>  7.177  20.197   5.583 1
    379    ATOM   379     N <NA>   ILE     A    50   <NA>  9.309  19.619   6.155 1
    380    ATOM   380    CA <NA>   ILE     A    50   <NA>  9.537  18.544   5.194 1
    381    ATOM   381     C <NA>   ILE     A    50   <NA>  8.802  17.261   5.658 1
    382    ATOM   382     O <NA>   ILE     A    50   <NA>  8.143  16.589   4.863 1
    383    ATOM   383    CB <NA>   ILE     A    50   <NA> 11.095  18.362   5.046 1
    384    ATOM   384   CG1 <NA>   ILE     A    50   <NA> 11.553  18.874   3.682 1
    385    ATOM   385   CG2 <NA>   ILE     A    50   <NA> 11.521  16.945   5.317 1
    386    ATOM   386   CD1 <NA>   ILE     A    50   <NA> 10.910  18.225   2.567 1
    387    ATOM   387     N <NA>   GLY     A    51   <NA>  8.865  16.952   6.945 1
    388    ATOM   388    CA <NA>   GLY     A    51   <NA>  8.174  15.771   7.405 1
    389    ATOM   389     C <NA>   GLY     A    51   <NA>  6.812  16.062   7.983 1
    390    ATOM   390     O <NA>   GLY     A    51   <NA>  6.408  15.342   8.870 1
    391    ATOM   391     N <NA>   GLY     A    52   <NA>  6.141  17.132   7.563 1
    392    ATOM   392    CA <NA>   GLY     A    52   <NA>  4.855  17.480   8.157 1
    393    ATOM   393     C <NA>   GLY     A    52   <NA>  4.884  18.624   9.170 1
    394    ATOM   394     O <NA>   GLY     A    52   <NA>  5.873  19.342   9.280 1
    395    ATOM   395     N <NA>   PHE     A    53   <NA>  3.806  18.788   9.925 1
    396    ATOM   396    CA <NA>   PHE     A    53   <NA>  3.593  19.996  10.731 1
    397    ATOM   397     C <NA>   PHE     A    53   <NA>  4.015  19.881  12.194 1
    398    ATOM   398     O <NA>   PHE     A    53   <NA>  3.930  18.810  12.781 1
    399    ATOM   399    CB <NA>   PHE     A    53   <NA>  2.121  20.351  10.670 1
    400    ATOM   400    CG <NA>   PHE     A    53   <NA>  1.760  21.152   9.484 1
    401    ATOM   401   CD1 <NA>   PHE     A    53   <NA>  1.725  20.567   8.216 1
    402    ATOM   402   CD2 <NA>   PHE     A    53   <NA>  1.556  22.518   9.607 1
    403    ATOM   403   CE1 <NA>   PHE     A    53   <NA>  1.500  21.332   7.075 1
    404    ATOM   404   CE2 <NA>   PHE     A    53   <NA>  1.327  23.302   8.496 1
    405    ATOM   405    CZ <NA>   PHE     A    53   <NA>  1.290  22.718   7.212 1
    406    ATOM   406     N <NA>   ILE     A    54   <NA>  4.483  20.969  12.792 1
    407    ATOM   407    CA <NA>   ILE     A    54   <NA>  4.689  21.005  14.248 1
    408    ATOM   408     C <NA>   ILE     A    54   <NA>  3.921  22.179  14.858 1
    409    ATOM   409     O <NA>   ILE     A    54   <NA>  3.575  23.139  14.182 1
    410    ATOM   410    CB <NA>   ILE     A    54   <NA>  6.199  21.155  14.625 1
    411    ATOM   411   CG1 <NA>   ILE     A    54   <NA>  6.796  22.408  13.939 1
    412    ATOM   412   CG2 <NA>   ILE     A    54   <NA>  6.967  19.888  14.203 1
    413    ATOM   413   CD1 <NA>   ILE     A    54   <NA>  8.110  22.821  14.465 1
    414    ATOM   414     N <NA>   LYS     A    55   <NA>  3.632  22.095  16.145 1
    415    ATOM   415    CA <NA>   LYS     A    55   <NA>  2.968  23.196  16.823 1
    416    ATOM   416     C <NA>   LYS     A    55   <NA>  4.038  24.093  17.449 1
    417    ATOM   417     O <NA>   LYS     A    55   <NA>  4.949  23.610  18.157 1
    418    ATOM   418    CB <NA>   LYS     A    55   <NA>  2.021  22.661  17.895 1
    419    ATOM   419    CG <NA>   LYS     A    55   <NA>  0.974  23.665  18.300 1
    420    ATOM   420    CD <NA>   LYS     A    55   <NA>  0.006  23.101  19.304 1
    421    ATOM   421    CE <NA>   LYS     A    55   <NA> -0.580  24.217  20.149 1
    422    ATOM   422    NZ <NA>   LYS     A    55   <NA>  0.439  24.751  21.104 1
    423    ATOM   423     N <NA>   VAL     A    56   <NA>  3.953  25.391  17.185 1
    424    ATOM   424    CA <NA>   VAL     A    56   <NA>  4.927  26.319  17.754 1
    425    ATOM   425     C <NA>   VAL     A    56   <NA>  4.225  27.378  18.556 1
    426    ATOM   426     O <NA>   VAL     A    56   <NA>  3.023  27.557  18.455 1
    427    ATOM   427    CB <NA>   VAL     A    56   <NA>  5.769  27.009  16.663 1
    428    ATOM   428   CG1 <NA>   VAL     A    56   <NA>  6.791  26.027  16.120 1
    429    ATOM   429   CG2 <NA>   VAL     A    56   <NA>  4.881  27.551  15.567 1
    430    ATOM   430     N <NA>   ARG     A    57   <NA>  4.978  28.079  19.377 1
    431    ATOM   431    CA <NA>   ARG     A    57   <NA>  4.459  29.271  20.026 1
    432    ATOM   432     C <NA>   ARG     A    57   <NA>  5.038  30.502  19.335 1
    433    ATOM   433     O <NA>   ARG     A    57   <NA>  6.242  30.596  19.138 1
    434    ATOM   434    CB <NA>   ARG     A    57   <NA>  4.824  29.245  21.505 1
    435    ATOM   435    CG <NA>   ARG     A    57   <NA>  4.168  28.102  22.260 1
    436    ATOM   436    CD <NA>   ARG     A    57   <NA>  4.656  28.068  23.693 1
    437    ATOM   437    NE <NA>   ARG     A    57   <NA>  6.032  27.573  23.790 1
    438    ATOM   438    CZ <NA>   ARG     A    57   <NA>  7.027  28.273  24.317 1
    439    ATOM   439   NH1 <NA>   ARG     A    57   <NA>  6.825  29.532  24.678 1
    440    ATOM   440   NH2 <NA>   ARG     A    57   <NA>  8.223  27.723  24.467 1
    441    ATOM   441     N <NA>   GLN     A    58   <NA>  4.171  31.431  18.958 1
    442    ATOM   442    CA <NA>   GLN     A    58   <NA>  4.570  32.596  18.172 1
    443    ATOM   443     C <NA>   GLN     A    58   <NA>  4.681  33.818  19.085 1
    444    ATOM   444     O <NA>   GLN     A    58   <NA>  3.694  34.242  19.683 1
    445    ATOM   445    CB <NA>   GLN     A    58   <NA>  3.539  32.859  17.094 1
    446    ATOM   446    CG <NA>   GLN     A    58   <NA>  3.736  34.104  16.321 1
    447    ATOM   447    CD <NA>   GLN     A    58   <NA>  2.500  34.473  15.541 1
    448    ATOM   448   OE1 <NA>   GLN     A    58   <NA>  1.530  33.703  15.489 1
    449    ATOM   449   NE2 <NA>   GLN     A    58   <NA>  2.508  35.651  14.940 1
    450    ATOM   450     N <NA>   TYR     A    59   <NA>  5.883  34.378  19.196 1
    451    ATOM   451    CA <NA>   TYR     A    59   <NA>  6.097  35.658  19.896 1
    452    ATOM   452     C <NA>   TYR     A    59   <NA>  6.304  36.752  18.835 1
    453    ATOM   453     O <NA>   TYR     A    59   <NA>  6.923  36.513  17.800 1
    454    ATOM   454    CB <NA>   TYR     A    59   <NA>  7.354  35.588  20.765 1
    455    ATOM   455    CG <NA>   TYR     A    59   <NA>  7.213  34.624  21.955 1
    456    ATOM   456   CD1 <NA>   TYR     A    59   <NA>  7.479  33.237  21.805 1
    457    ATOM   457   CD2 <NA>   TYR     A    59   <NA>  6.795  35.087  23.223 1
    458    ATOM   458   CE1 <NA>   TYR     A    59   <NA>  7.345  32.365  22.871 1
    459    ATOM   459   CE2 <NA>   TYR     A    59   <NA>  6.638  34.224  24.268 1
    460    ATOM   460    CZ <NA>   TYR     A    59   <NA>  6.926  32.869  24.102 1
    461    ATOM   461    OH <NA>   TYR     A    59   <NA>  6.967  32.076  25.200 1
    462    ATOM   462     N <NA>   ASP     A    60   <NA>  5.767  37.940  19.049 1
    463    ATOM   463    CA <NA>   ASP     A    60   <NA>  6.022  39.045  18.125 1
    464    ATOM   464     C <NA>   ASP     A    60   <NA>  7.025  40.015  18.725 1
    465    ATOM   465     O <NA>   ASP     A    60   <NA>  7.340  39.951  19.900 1
    466    ATOM   466    CB <NA>   ASP     A    60   <NA>  4.719  39.777  17.832 1
    467    ATOM   467    CG <NA>   ASP     A    60   <NA>  3.699  38.899  17.148 1
    468    ATOM   468   OD1 <NA>   ASP     A    60   <NA>  3.989  38.368  16.050 1
    469    ATOM   469   OD2 <NA>   ASP     A    60   <NA>  2.570  38.799  17.672 1
    470    ATOM   470     N <NA>   GLN     A    61   <NA>  7.529  40.913  17.896 1
    471    ATOM   471    CA <NA>   GLN     A    61   <NA>  8.337  42.060  18.331 1
    472    ATOM   472     C <NA>   GLN     A    61   <NA>  9.535  41.630  19.179 1
    473    ATOM   473     O <NA>   GLN     A    61   <NA>  9.777  42.191  20.264 1
    474    ATOM   474    CB <NA>   GLN     A    61   <NA>  7.471  43.051  19.131 1
    475    ATOM   475    CG <NA>   GLN     A    61   <NA>  7.718  44.555  18.814 1
    476    ATOM   476    CD <NA>   GLN     A    61   <NA>  7.182  45.552  19.907 1
    477    ATOM   477   OE1 <NA>   GLN     A    61   <NA>  7.936  46.461  20.398 1
    478    ATOM   478   NE2 <NA>   GLN     A    61   <NA>  5.892  45.377  20.306 1
    479    ATOM   479     N <NA>   ILE     A    62   <NA> 10.283  40.645  18.676 1
    480    ATOM   480    CA <NA>   ILE     A    62   <NA> 11.484  40.115  19.328 1
    481    ATOM   481     C <NA>   ILE     A    62   <NA> 12.745  40.584  18.614 1
    482    ATOM   482     O <NA>   ILE     A    62   <NA> 12.830  40.485  17.396 1
    483    ATOM   483    CB <NA>   ILE     A    62   <NA> 11.465  38.545  19.309 1
    484    ATOM   484   CG1 <NA>   ILE     A    62   <NA> 10.152  38.035  19.930 1
    485    ATOM   485   CG2 <NA>   ILE     A    62   <NA> 12.688  37.973  20.027 1
    486    ATOM   486   CD1 <NA>   ILE     A    62   <NA>  9.966  38.380  21.363 1
    487    ATOM   487     N <NA>   LEU     A    63   <NA> 13.722  41.086  19.369 1
    488    ATOM   488    CA <NA>   LEU     A    63   <NA> 15.038  41.476  18.822 1
    489    ATOM   489     C <NA>   LEU     A    63   <NA> 16.033  40.304  18.862 1
    490    ATOM   490     O <NA>   LEU     A    63   <NA> 16.195  39.661  19.897 1
    491    ATOM   491    CB <NA>   LEU     A    63   <NA> 15.631  42.666  19.619 1
    492    ATOM   492    CG <NA>   LEU     A    63   <NA> 16.776  43.426  18.914 1
    493    ATOM   493   CD1 <NA>   LEU     A    63   <NA> 16.560  44.922  18.993 1
    494    ATOM   494   CD2 <NA>   LEU     A    63   <NA> 18.103  43.062  19.558 1
    495    ATOM   495     N <NA>   ILE     A    64   <NA> 16.686  40.036  17.738 1
    496    ATOM   496    CA <NA>   ILE     A    64   <NA> 17.760  39.039  17.653 1
    497    ATOM   497     C <NA>   ILE     A    64   <NA> 18.991  39.753  17.116 1
    498    ATOM   498     O <NA>   ILE     A    64   <NA> 18.862  40.632  16.288 1
    499    ATOM   499    CB <NA>   ILE     A    64   <NA> 17.390  37.875  16.634 1
    500    ATOM   500   CG1 <NA>   ILE     A    64   <NA> 16.127  37.137  17.110 1
    501    ATOM   501   CG2 <NA>   ILE     A    64   <NA> 18.551  36.857  16.518 1
    502    ATOM   502   CD1 <NA>   ILE     A    64   <NA> 16.194  35.643  16.885 1
    503    ATOM   503     N <NA>   GLU     A    65   <NA> 20.181  39.387  17.573 1
    504    ATOM   504    CA <NA>   GLU     A    65   <NA> 21.406  39.983  17.036 1
    505    ATOM   505     C <NA>   GLU     A    65   <NA> 22.192  39.019  16.135 1
    506    ATOM   506     O <NA>   GLU     A    65   <NA> 22.866  38.128  16.631 1
    507    ATOM   507    CB <NA>   GLU     A    65   <NA> 22.307  40.439  18.180 1
    508    ATOM   508    CG <NA>   GLU     A    65   <NA> 21.987  41.843  18.676 1
    509    ATOM   509    CD <NA>   GLU     A    65   <NA> 23.228  42.558  19.187 1
    510    ATOM   510   OE1 <NA>   GLU     A    65   <NA> 24.094  42.925  18.348 1
    511    ATOM   511   OE2 <NA>   GLU     A    65   <NA> 23.348  42.730  20.429 1
    512    ATOM   512     N <NA>   ILE     A    66   <NA> 22.148  39.215  14.823 1
    513    ATOM   513    CA <NA>   ILE     A    66   <NA> 22.818  38.327  13.863 1
    514    ATOM   514     C <NA>   ILE     A    66   <NA> 24.123  38.945  13.303 1
    515    ATOM   515     O <NA>   ILE     A    66   <NA> 24.114  39.903  12.521 1
    516    ATOM   516    CB <NA>   ILE     A    66   <NA> 21.831  37.964  12.694 1
    517    ATOM   517   CG1 <NA>   ILE     A    66   <NA> 20.481  37.516  13.278 1
    518    ATOM   518   CG2 <NA>   ILE     A    66   <NA> 22.444  36.880  11.799 1
    519    ATOM   519   CD1 <NA>   ILE     A    66   <NA> 19.263  38.072  12.560 1
    520    ATOM   520     N <NA>   CYS     A    67   <NA> 25.251  38.396  13.731 1
    521    ATOM   521    CA <NA>   CYS     A    67   <NA> 26.590  38.871  13.351 1
    522    ATOM   522     C <NA>   CYS     A    67   <NA> 26.798  40.370  13.560 1
    523    ATOM   523     O <NA>   CYS     A    67   <NA> 27.461  41.045  12.763 1
    524    ATOM   524    CB <NA>   CYS     A    67   <NA> 26.907  38.497  11.905 1
    525    ATOM   525    SG <NA>   CYS     A    67   <NA> 27.274  36.734  11.705 1
    526    ATOM   526     N <NA>   GLY     A    68   <NA> 26.239  40.887  14.650 1
    527    ATOM   527    CA <NA>   GLY     A    68   <NA> 26.319  42.307  14.938 1
    528    ATOM   528     C <NA>   GLY     A    68   <NA> 24.980  42.946  14.665 1
    529    ATOM   529     O <NA>   GLY     A    68   <NA> 24.371  43.518  15.561 1
    530    ATOM   530     N <NA>   HIS     A    69   <NA> 24.449  42.689  13.479 1
    531    ATOM   531    CA <NA>   HIS     A    69   <NA> 23.251  43.368  13.013 1
    532    ATOM   532     C <NA>   HIS     A    69   <NA> 22.009  43.025  13.836 1
    533    ATOM   533     O <NA>   HIS     A    69   <NA> 21.626  41.865  13.958 1
    534    ATOM   534    CB <NA>   HIS     A    69   <NA> 22.975  43.009  11.568 1
    535    ATOM   535    CG <NA>   HIS     A    69   <NA> 24.157  43.172  10.639 1
    536    ATOM   536   ND1 <NA>   HIS     A    69   <NA> 25.294  42.398  10.762 1
    537    ATOM   537   CD2 <NA>   HIS     A    69   <NA> 24.283  43.886   9.518 1
    538    ATOM   538   CE1 <NA>   HIS     A    69   <NA> 26.075  42.630   9.723 1
    539    ATOM   539   NE2 <NA>   HIS     A    69   <NA> 25.505  43.532   8.946 1
    540    ATOM   540     N <NA>   LYS     A    70   <NA> 21.364  44.050  14.376 1
    541    ATOM   541    CA <NA>   LYS     A    70   <NA> 20.030  43.883  14.945 1
    542    ATOM   542     C <NA>   LYS     A    70   <NA> 18.929  43.698  13.882 1
    543    ATOM   543     O <NA>   LYS     A    70   <NA> 19.000  44.197  12.761 1
    544    ATOM   544    CB <NA>   LYS     A    70   <NA> 19.667  45.075  15.840 1
    545    ATOM   545    CG <NA>   LYS     A    70   <NA> 20.422  45.141  17.154 1
    546    ATOM   546    CD <NA>   LYS     A    70   <NA> 19.885  46.281  18.036 1
    547    ATOM   547    CE <NA>   LYS     A    70   <NA> 20.905  46.699  19.095 1
    548    ATOM   548    NZ <NA>   LYS     A    70   <NA> 20.376  47.817  19.946 1
    549    ATOM   549     N <NA>   ALA     A    71   <NA> 17.879  43.010  14.305 1
    550    ATOM   550    CA <NA>   ALA     A    71   <NA> 16.703  42.737  13.499 1
    551    ATOM   551     C <NA>   ALA     A    71   <NA> 15.600  42.546  14.556 1
    552    ATOM   552     O <NA>   ALA     A    71   <NA> 15.880  42.175  15.705 1
    553    ATOM   553    CB <NA>   ALA     A    71   <NA> 16.937  41.464  12.683 1
    554    ATOM   554     N <NA>   ILE     A    72   <NA> 14.366  42.881  14.219 1
    555    ATOM   555    CA <NA>   ILE     A    72   <NA> 13.233  42.663  15.111 1
    556    ATOM   556     C <NA>   ILE     A    72   <NA> 12.196  41.969  14.237 1
    557    ATOM   557     O <NA>   ILE     A    72   <NA> 12.083  42.266  13.025 1
    558    ATOM   558    CB <NA>   ILE     A    72   <NA> 12.635  44.019  15.609 1
    559    ATOM   559   CG1 <NA>   ILE     A    72   <NA> 13.664  44.777  16.446 1
    560    ATOM   560   CG2 <NA>   ILE     A    72   <NA> 11.353  43.782  16.409 1
    561    ATOM   561   CD1 <NA>   ILE     A    72   <NA> 13.437  46.286  16.451 1
    562    ATOM   562     N <NA>   GLY     A    73   <NA> 11.457  41.034  14.821 1
    563    ATOM   563    CA <NA>   GLY     A    73   <NA> 10.359  40.420  14.090 1
    564    ATOM   564     C <NA>   GLY     A    73   <NA>  9.750  39.238  14.835 1
    565    ATOM   565     O <NA>   GLY     A    73   <NA> 10.019  39.007  16.009 1
    566    ATOM   566     N <NA>   THR     A    74   <NA>  8.853  38.524  14.181 1
    567    ATOM   567    CA <NA>   THR     A    74   <NA>  8.224  37.361  14.786 1
    568    ATOM   568     C <NA>   THR     A    74   <NA>  9.169  36.194  14.895 1
    569    ATOM   569     O <NA>   THR     A    74   <NA>  9.964  35.908  14.007 1
    570    ATOM   570    CB <NA>   THR     A    74   <NA>  6.961  36.929  14.010 1
    571    ATOM   571   OG1 <NA>   THR     A    74   <NA>  6.030  38.013  14.001 1
    572    ATOM   572   CG2 <NA>   THR     A    74   <NA>  6.287  35.760  14.688 1
    573    ATOM   573     N <NA>   VAL     A    75   <NA>  9.085  35.534  16.025 1
    574    ATOM   574    CA <NA>   VAL     A    75   <NA>  9.944  34.420  16.291 1
    575    ATOM   575     C <NA>   VAL     A    75   <NA>  9.018  33.325  16.795 1
    576    ATOM   576     O <NA>   VAL     A    75   <NA>  8.104  33.576  17.594 1
    577    ATOM   577    CB <NA>   VAL     A    75   <NA> 11.002  34.791  17.336 1
    578    ATOM   578   CG1 <NA>   VAL     A    75   <NA> 11.638  33.561  17.884 1
    579    ATOM   579   CG2 <NA>   VAL     A    75   <NA> 12.070  35.661  16.648 1
    580    ATOM   580     N <NA>   LEU     A    76   <NA>  9.200  32.127  16.258 1
    581    ATOM   581    CA <NA>   LEU     A    76   <NA>  8.353  30.995  16.628 1
    582    ATOM   582     C <NA>   LEU     A    76   <NA>  9.229  30.084  17.483 1
    583    ATOM   583     O <NA>   LEU     A    76   <NA> 10.409  29.923  17.194 1
    584    ATOM   584    CB <NA>   LEU     A    76   <NA>  7.937  30.247  15.365 1
    585    ATOM   585    CG <NA>   LEU     A    76   <NA>  7.222  30.919  14.207 1
    586    ATOM   586   CD1 <NA>   LEU     A    76   <NA>  6.909  29.914  13.139 1
    587    ATOM   587   CD2 <NA>   LEU     A    76   <NA>  5.961  31.526  14.713 1
    588    ATOM   588     N <NA>   VAL     A    77   <NA>  8.689  29.507  18.538 1
    589    ATOM   589    CA <NA>   VAL     A    77   <NA>  9.519  28.691  19.405 1
    590    ATOM   590     C <NA>   VAL     A    77   <NA>  8.852  27.309  19.499 1
    591    ATOM   591     O <NA>   VAL     A    77   <NA>  7.621  27.241  19.608 1
    592    ATOM   592    CB <NA>   VAL     A    77   <NA>  9.648  29.353  20.800 1
    593    ATOM   593   CG1 <NA>   VAL     A    77   <NA> 10.306  28.419  21.773 1
    594    ATOM   594   CG2 <NA>   VAL     A    77   <NA> 10.426  30.686  20.648 1
    595    ATOM   595     N <NA>   GLY     A    78   <NA>  9.637  26.227  19.375 1
    596    ATOM   596    CA <NA>   GLY     A    78   <NA>  9.027  24.918  19.185 1
    597    ATOM   597     C <NA>   GLY     A    78   <NA> 10.005  23.782  18.984 1
    598    ATOM   598     O <NA>   GLY     A    78   <NA> 11.215  24.012  19.065 1
    599    ATOM   599     N <NA>   PRO     A    79   <NA>  9.538  22.542  18.758 1
    600    ATOM   600    CA <NA>   PRO     A    79   <NA> 10.442  21.388  18.731 1
    601    ATOM   601     C <NA>   PRO     A    79   <NA> 11.136  21.242  17.389 1
    602    ATOM   602     O <NA>   PRO     A    79   <NA> 10.748  20.397  16.559 1
    603    ATOM   603    CB <NA>   PRO     A    79   <NA>  9.515  20.215  18.995 1
    604    ATOM   604    CG <NA>   PRO     A    79   <NA>  8.239  20.616  18.339 1
    605    ATOM   605    CD <NA>   PRO     A    79   <NA>  8.170  22.144  18.383 1
    606    ATOM   606     N <NA>   THR     A    80   <NA> 12.138  22.078  17.162 1
    607    ATOM   607    CA <NA>   THR     A    80   <NA> 12.910  22.043  15.936 1
    608    ATOM   608     C <NA>   THR     A    80   <NA> 14.280  21.503  16.325 1
    609    ATOM   609     O <NA>   THR     A    80   <NA> 14.692  21.666  17.473 1
    610    ATOM   610    CB <NA>   THR     A    80   <NA> 13.031  23.467  15.319 1
    611    ATOM   611   OG1 <NA>   THR     A    80   <NA> 13.971  23.419  14.244 1
    612    ATOM   612   CG2 <NA>   THR     A    80   <NA> 13.542  24.478  16.300 1
    613    ATOM   613     N <NA>   PRO     A    81   <NA> 14.915  20.725  15.449 1
    614    ATOM   614    CA <NA>   PRO     A    81   <NA> 16.233  20.170  15.766 1
    615    ATOM   615     C <NA>   PRO     A    81   <NA> 17.352  21.205  15.795 1
    616    ATOM   616     O <NA>   PRO     A    81   <NA> 18.306  21.085  16.567 1
    617    ATOM   617    CB <NA>   PRO     A    81   <NA> 16.467  19.139  14.669 1
    618    ATOM   618    CG <NA>   PRO     A    81   <NA> 15.477  19.489  13.607 1
    619    ATOM   619    CD <NA>   PRO     A    81   <NA> 14.331  20.114  14.238 1
    620    ATOM   620     N <NA>   VAL     A    82   <NA> 17.220  22.219  14.948 1
    621    ATOM   621    CA <NA>   VAL     A    82   <NA> 18.236  23.243  14.762 1
    622    ATOM   622     C <NA>   VAL     A    82   <NA> 17.549  24.626  14.801 1
    623    ATOM   623     O <NA>   VAL     A    82   <NA> 16.328  24.735  14.569 1
    624    ATOM   624    CB <NA>   VAL     A    82   <NA> 18.969  23.017  13.388 1
    625    ATOM   625   CG1 <NA>   VAL     A    82   <NA> 18.009  22.769  12.263 1
    626    ATOM   626   CG2 <NA>   VAL     A    82   <NA> 19.843  24.144  13.080 1
    627    ATOM   627     N <NA>   ASN     A    83   <NA> 18.292  25.671  15.167 1
    628    ATOM   628    CA <NA>   ASN     A    83   <NA> 17.731  27.051  15.193 1
    629    ATOM   629     C <NA>   ASN     A    83   <NA> 17.799  27.628  13.799 1
    630    ATOM   630     O <NA>   ASN     A    83   <NA> 18.794  27.456  13.113 1
    631    ATOM   631    CB <NA>   ASN     A    83   <NA> 18.508  27.950  16.161 1
    632    ATOM   632    CG <NA>   ASN     A    83   <NA> 18.236  27.600  17.636 1
    633    ATOM   633   OD1 <NA>   ASN     A    83   <NA> 17.102  27.381  18.051 1
    634    ATOM   634   ND2 <NA>   ASN     A    83   <NA> 19.287  27.475  18.398 1
    635    ATOM   635     N <NA>   ILE     A    84   <NA> 16.729  28.283  13.370 1
    636    ATOM   636    CA <NA>   ILE     A    84   <NA> 16.588  28.651  11.955 1
    637    ATOM   637     C <NA>   ILE     A    84   <NA> 16.190  30.132  11.811 1
    638    ATOM   638     O <NA>   ILE     A    84   <NA> 15.178  30.600  12.384 1
    639    ATOM   639    CB <NA>   ILE     A    84   <NA> 15.504  27.755  11.242 1
    640    ATOM   640   CG1 <NA>   ILE     A    84   <NA> 16.049  26.365  11.040 1
    641    ATOM   641   CG2 <NA>   ILE     A    84   <NA> 15.081  28.324   9.932 1
    642    ATOM   642   CD1 <NA>   ILE     A    84   <NA> 14.957  25.345  11.321 1
    643    ATOM   643     N <NA>   ILE     A    85   <NA> 16.991  30.863  11.053 1
    644    ATOM   644    CA <NA>   ILE     A    85   <NA> 16.710  32.257  10.821 1
    645    ATOM   645     C <NA>   ILE     A    85   <NA> 16.155  32.252   9.420 1
    646    ATOM   646     O <NA>   ILE     A    85   <NA> 16.857  31.950   8.464 1
    647    ATOM   647    CB <NA>   ILE     A    85   <NA> 17.999  33.119  10.844 1
    648    ATOM   648   CG1 <NA>   ILE     A    85   <NA> 18.724  32.983  12.207 1
    649    ATOM   649   CG2 <NA>   ILE     A    85   <NA> 17.653  34.600  10.517 1
    650    ATOM   650   CD1 <NA>   ILE     A    85   <NA> 17.830  33.204  13.418 1
    651    ATOM   651     N <NA>   GLY     A    86   <NA> 14.888  32.608   9.302 1
    652    ATOM   652    CA <NA>   GLY     A    86   <NA> 14.238  32.655   8.011 1
    653    ATOM   653     C <NA>   GLY     A    86   <NA> 14.181  34.034   7.393 1
    654    ATOM   654     O <NA>   GLY     A    86   <NA> 14.714  34.991   7.894 1
    655    ATOM   655     N <NA>   ARG     A    87   <NA> 13.398  34.131   6.334 1
    656    ATOM   656    CA <NA>   ARG     A    87   <NA> 13.370  35.302   5.456 1
    657    ATOM   657     C <NA>   ARG     A    87   <NA> 12.931  36.548   6.216 1
    658    ATOM   658     O <NA>   ARG     A    87   <NA> 13.429  37.633   5.958 1
    659    ATOM   659    CB <NA>   ARG     A    87   <NA> 12.459  35.018   4.263 1
    660    ATOM   660    CG <NA>   ARG     A    87   <NA> 13.030  34.044   3.318 1
    661    ATOM   661    CD <NA>   ARG     A    87   <NA> 12.310  34.101   1.966 1
    662    ATOM   662    NE <NA>   ARG     A    87   <NA> 10.903  33.713   2.088 1
    663    ATOM   663    CZ <NA>   ARG     A    87   <NA>  9.899  34.575   2.137 1
    664    ATOM   664   NH1 <NA>   ARG     A    87   <NA> 10.126  35.860   1.922 1
    665    ATOM   665   NH2 <NA>   ARG     A    87   <NA>  8.657  34.143   2.300 1
    666    ATOM   666     N <NA>   ASN     A    88   <NA> 12.063  36.352   7.198 1
    667    ATOM   667    CA <NA>   ASN     A    88   <NA> 11.458  37.455   7.919 1
    668    ATOM   668     C <NA>   ASN     A    88   <NA> 12.511  38.266   8.653 1
    669    ATOM   669     O <NA>   ASN     A    88   <NA> 12.419  39.481   8.723 1
    670    ATOM   670    CB <NA>   ASN     A    88   <NA> 10.359  36.959   8.863 1
    671    ATOM   671    CG <NA>   ASN     A    88   <NA> 10.892  36.309  10.125 1
    672    ATOM   672   OD1 <NA>   ASN     A    88   <NA> 11.577  35.273  10.065 1
    673    ATOM   673   ND2 <NA>   ASN     A    88   <NA> 10.446  36.809  11.283 1
    674    ATOM   674     N <NA>   LEU     A    89   <NA> 13.551  37.601   9.136 1
    675    ATOM   675    CA <NA>   LEU     A    89   <NA> 14.633  38.337   9.729 1
    676    ATOM   676     C <NA>   LEU     A    89   <NA> 15.802  38.487   8.783 1
    677    ATOM   677     O <NA>   LEU     A    89   <NA> 16.580  39.394   8.941 1
    678    ATOM   678    CB <NA>   LEU     A    89   <NA> 15.096  37.687  11.029 1
    679    ATOM   679    CG <NA>   LEU     A    89   <NA> 14.146  37.555  12.221 1
    680    ATOM   680   CD1 <NA>   LEU     A    89   <NA> 14.939  37.022  13.434 1
    681    ATOM   681   CD2 <NA>   LEU     A    89   <NA> 13.508  38.913  12.527 1
    682    ATOM   682     N <NA>   LEU     A    90   <NA> 15.910  37.651   7.761 1
    683    ATOM   683    CA <NA>   LEU     A    90   <NA> 17.078  37.773   6.865 1
    684    ATOM   684     C <NA>   LEU     A    90   <NA> 17.048  39.068   6.049 1
    685    ATOM   685     O <NA>   LEU     A    90   <NA> 18.098  39.613   5.652 1
    686    ATOM   686    CB <NA>   LEU     A    90   <NA> 17.200  36.561   5.913 1
    687    ATOM   687    CG <NA>   LEU     A    90   <NA> 17.754  35.226   6.440 1
    688    ATOM   688   CD1 <NA>   LEU     A    90   <NA> 17.798  34.207   5.340 1
    689    ATOM   689   CD2 <NA>   LEU     A    90   <NA> 19.138  35.441   7.022 1
    690    ATOM   690     N <NA>   THR     A    91   <NA> 15.835  39.547   5.808 1
    691    ATOM   691    CA <NA>   THR     A    91   <NA> 15.653  40.738   4.985 1
    692    ATOM   692     C <NA>   THR     A    91   <NA> 16.137  41.999   5.772 1
    693    ATOM   693     O <NA>   THR     A    91   <NA> 16.848  42.870   5.252 1
    694    ATOM   694    CB <NA>   THR     A    91   <NA> 14.157  40.860   4.594 1
    695    ATOM   695   OG1 <NA>   THR     A    91   <NA> 13.342  40.810   5.786 1
    696    ATOM   696   CG2 <NA>   THR     A    91   <NA> 13.740  39.709   3.681 1
    697    ATOM   697     N <NA>   GLN     A    92   <NA> 15.842  42.001   7.064 1
    698    ATOM   698    CA <NA>   GLN     A    92   <NA> 16.170  43.099   7.967 1
    699    ATOM   699     C <NA>   GLN     A    92   <NA> 17.650  43.315   8.126 1
    700    ATOM   700     O <NA>   GLN     A    92   <NA> 18.071  44.440   8.459 1
    701    ATOM   701    CB <NA>   GLN     A    92   <NA> 15.580  42.837   9.347 1
    702    ATOM   702    CG <NA>   GLN     A    92   <NA> 14.081  42.603   9.325 1
    703    ATOM   703    CD <NA>   GLN     A    92   <NA> 13.318  43.811   8.800 1
    704    ATOM   704   OE1 <NA>   GLN     A    92   <NA> 12.836  43.807   7.683 1
    705    ATOM   705   NE2 <NA>   GLN     A    92   <NA> 13.349  44.892   9.543 1
    706    ATOM   706     N <NA>   ILE     A    93   <NA> 18.444  42.266   7.913 1
    707    ATOM   707    CA <NA>   ILE     A    93   <NA> 19.887  42.417   8.035 1
    708    ATOM   708     C <NA>   ILE     A    93   <NA> 20.530  42.555   6.677 1
    709    ATOM   709     O <NA>   ILE     A    93   <NA> 21.744  42.639   6.577 1
    710    ATOM   710    CB <NA>   ILE     A    93   <NA> 20.557  41.231   8.872 1
    711    ATOM   711   CG1 <NA>   ILE     A    93   <NA> 20.472  39.878   8.113 1
    712    ATOM   712   CG2 <NA>   ILE     A    93   <NA> 19.926  41.182  10.272 1
    713    ATOM   713   CD1 <NA>   ILE     A    93   <NA> 21.520  38.796   8.533 1
    714    ATOM   714     N <NA>   GLY     A    94   <NA> 19.713  42.563   5.636 1
    715    ATOM   715    CA <NA>   GLY     A    94   <NA> 20.226  42.870   4.327 1
    716    ATOM   716     C <NA>   GLY     A    94   <NA> 20.843  41.689   3.611 1
    717    ATOM   717     O <NA>   GLY     A    94   <NA> 21.680  41.866   2.723 1
    718    ATOM   718     N <NA>   CYS     A    95   <NA> 20.294  40.507   3.844 1
    719    ATOM   719    CA <NA>   CYS     A    95   <NA> 20.890  39.296   3.297 1
    720    ATOM   720     C <NA>   CYS     A    95   <NA> 20.407  39.005   1.897 1
    721    ATOM   721     O <NA>   CYS     A    95   <NA> 19.217  39.101   1.606 1
    722    ATOM   722    CB <NA>   CYS     A    95   <NA> 20.585  38.134   4.215 1
    723    ATOM   723    SG <NA>   CYS     A    95   <NA> 21.505  36.663   3.848 1
    724    ATOM   724     N <NA>   THR     A    96   <NA> 21.339  38.698   1.005 1
    725    ATOM   725    CA <NA>   THR     A    96   <NA> 20.969  38.276  -0.346 1
    726    ATOM   726     C <NA>   THR     A    96   <NA> 21.721  37.024  -0.758 1
    727    ATOM   727     O <NA>   THR     A    96   <NA> 22.753  36.645  -0.178 1
    728    ATOM   728    CB <NA>   THR     A    96   <NA> 21.298  39.331  -1.459 1
    729    ATOM   729   OG1 <NA>   THR     A    96   <NA> 22.688  39.698  -1.364 1
    730    ATOM   730   CG2 <NA>   THR     A    96   <NA> 20.405  40.560  -1.303 1
    731    ATOM   731     N <NA>   LEU     A    97   <NA> 21.183  36.413  -1.805 1
    732    ATOM   732    CA <NA>   LEU     A    97   <NA> 21.718  35.228  -2.473 1
    733    ATOM   733     C <NA>   LEU     A    97   <NA> 22.314  35.753  -3.779 1
    734    ATOM   734     O <NA>   LEU     A    97   <NA> 21.616  36.396  -4.577 1
    735    ATOM   735    CB <NA>   LEU     A    97   <NA> 20.518  34.349  -2.810 1
    736    ATOM   736    CG <NA>   LEU     A    97   <NA> 20.318  32.897  -2.464 1
    737    ATOM   737   CD1 <NA>   LEU     A    97   <NA> 21.285  32.388  -1.422 1
    738    ATOM   738   CD2 <NA>   LEU     A    97   <NA> 18.896  32.782  -2.010 1
    739    ATOM   739     N <NA>   ASN     A    98   <NA> 23.577  35.469  -4.026 1
    740    ATOM   740    CA <NA>   ASN     A    98   <NA> 24.207  35.904  -5.257 1
    741    ATOM   741     C <NA>   ASN     A    98   <NA> 24.863  34.748  -5.999 1
    742    ATOM   742     O <NA>   ASN     A    98   <NA> 25.635  33.982  -5.396 1
    743    ATOM   743    CB <NA>   ASN     A    98   <NA> 25.229  36.991  -4.938 1
    744    ATOM   744    CG <NA>   ASN     A    98   <NA> 24.609  38.190  -4.212 1
    745    ATOM   745   OD1 <NA>   ASN     A    98   <NA> 24.263  38.112  -3.019 1
    746    ATOM   746   ND2 <NA>   ASN     A    98   <NA> 24.354  39.252  -4.955 1
    747    ATOM   747     N <NA>   PHE     A    99   <NA> 24.539  34.602  -7.287 1
    748    ATOM   748    CA <NA>   PHE     A    99   <NA> 25.376  33.754  -8.172 1
    749    ATOM   749     C <NA>   PHE     A    99   <NA> 25.726  34.383  -9.544 1
    750    ATOM   750     O <NA>   PHE     A    99   <NA> 24.797  34.890 -10.226 1
    751    ATOM   751    CB <NA>   PHE     A    99   <NA> 24.743  32.354  -8.375 1
    752    ATOM   752    CG <NA>   PHE     A    99   <NA> 23.328  32.389  -8.856 1
    753    ATOM   753   CD1 <NA>   PHE     A    99   <NA> 22.303  32.588  -7.961 1
    754    ATOM   754   CD2 <NA>   PHE     A    99   <NA> 23.029  32.118 -10.181 1
    755    ATOM   755   CE1 <NA>   PHE     A    99   <NA> 21.000  32.515  -8.357 1
    756    ATOM   756   CE2 <NA>   PHE     A    99   <NA> 21.730  32.028 -10.613 1
    757    ATOM   757    CZ <NA>   PHE     A    99   <NA> 20.700  32.221  -9.700 1
    758    ATOM   759     N <NA>   PRO     B     1   <NA> 22.659  36.727 -10.823 1
    759    ATOM   760    CA <NA>   PRO     B     1   <NA> 21.708  37.741 -10.269 1
    760    ATOM   761     C <NA>   PRO     B     1   <NA> 21.931  37.939  -8.779 1
    761    ATOM   762     O <NA>   PRO     B     1   <NA> 22.755  37.283  -8.190 1
    762    ATOM   763    CB <NA>   PRO     B     1   <NA> 20.263  37.289 -10.512 1
    763    ATOM   764    CG <NA>   PRO     B     1   <NA> 20.385  35.811 -10.891 1
    764    ATOM   765    CD <NA>   PRO     B     1   <NA> 21.753  35.755 -11.555 1
    765    ATOM   766     N <NA>   GLN     B     2   <NA> 21.203  38.873  -8.191 1
    766    ATOM   767    CA <NA>   GLN     B     2   <NA> 21.156  39.043  -6.744 1
    767    ATOM   768     C <NA>   GLN     B     2   <NA> 19.698  38.882  -6.389 1
    768    ATOM   769     O <NA>   GLN     B     2   <NA> 18.850  39.538  -6.975 1
    769    ATOM   770    CB <NA>   GLN     B     2   <NA> 21.625  40.447  -6.329 1
    770    ATOM   771    CG <NA>   GLN     B     2   <NA> 21.353  40.777  -4.865 1
    771    ATOM   772    CD <NA>   GLN     B     2   <NA> 22.139  41.975  -4.358 1
    772    ATOM   773   OE1 <NA>   GLN     B     2   <NA> 21.577  42.881  -3.752 1
    773    ATOM   774   NE2 <NA>   GLN     B     2   <NA> 23.450  41.976  -4.586 1
    774    ATOM   775     N <NA>   ILE     B     3   <NA> 19.405  38.008  -5.448 1
    775    ATOM   776    CA <NA>   ILE     B     3   <NA> 18.037  37.742  -5.100 1
    776    ATOM   777     C <NA>   ILE     B     3   <NA> 17.832  38.164  -3.660 1
    777    ATOM   778     O <NA>   ILE     B     3   <NA> 18.457  37.631  -2.746 1
    778    ATOM   779    CB <NA>   ILE     B     3   <NA> 17.694  36.224  -5.238 1
    779    ATOM   780   CG1 <NA>   ILE     B     3   <NA> 17.788  35.772  -6.692 1
    780    ATOM   781   CG2 <NA>   ILE     B     3   <NA> 16.284  35.967  -4.750 1
    781    ATOM   782   CD1 <NA>   ILE     B     3   <NA> 18.183  34.327  -6.802 1
    782    ATOM   783     N <NA>   THR     B     4   <NA> 16.960  39.132  -3.461 1
    783    ATOM   784    CA <NA>   THR     B     4   <NA> 16.635  39.571  -2.117 1
    784    ATOM   785     C <NA>   THR     B     4   <NA> 15.555  38.634  -1.627 1
    785    ATOM   786     O <NA>   THR     B     4   <NA> 15.066  37.784  -2.372 1
    786    ATOM   787    CB <NA>   THR     B     4   <NA> 16.147  41.074  -2.110 1
    787    ATOM   788   OG1 <NA>   THR     B     4   <NA> 15.093  41.256  -3.079 1
    788    ATOM   789   CG2 <NA>   THR     B     4   <NA> 17.283  42.019  -2.472 1
    789    ATOM   790     N <NA>   LEU     B     5   <NA> 15.157  38.761  -0.379 1
    790    ATOM   791    CA <NA>   LEU     B     5   <NA> 14.466  37.636   0.259 1
    791    ATOM   792     C <NA>   LEU     B     5   <NA> 13.125  38.065   0.821 1
    792    ATOM   793     O <NA>   LEU     B     5   <NA> 12.585  37.446   1.733 1
    793    ATOM   794    CB <NA>   LEU     B     5   <NA> 15.340  37.047   1.374 1
    794    ATOM   795    CG <NA>   LEU     B     5   <NA> 16.622  36.365   0.892 1
    795    ATOM   796   CD1 <NA>   LEU     B     5   <NA> 17.455  35.953   2.080 1
    796    ATOM   797   CD2 <NA>   LEU     B     5   <NA> 16.248  35.136  -0.006 1
    797    ATOM   798     N <NA>   TRP     B     6   <NA> 12.567  39.127   0.262 1
    798    ATOM   799    CA <NA>   TRP     B     6   <NA> 11.260  39.582   0.682 1
    799    ATOM   800     C <NA>   TRP     B     6   <NA> 10.196  38.601   0.218 1
    800    ATOM   801     O <NA>   TRP     B     6   <NA>  9.192  38.404   0.903 1
    801    ATOM   802    CB <NA>   TRP     B     6   <NA> 11.004  40.992   0.135 1
    802    ATOM   803    CG <NA>   TRP     B     6   <NA> 12.065  42.014   0.478 1
    803    ATOM   804   CD1 <NA>   TRP     B     6   <NA> 13.157  42.366  -0.279 1
    804    ATOM   805   CD2 <NA>   TRP     B     6   <NA> 12.209  42.716   1.739 1
    805    ATOM   806   NE1 <NA>   TRP     B     6   <NA> 13.979  43.196   0.470 1
    806    ATOM   807   CE2 <NA>   TRP     B     6   <NA> 13.433  43.441   1.686 1
    807    ATOM   808   CE3 <NA>   TRP     B     6   <NA> 11.443  42.805   2.913 1
    808    ATOM   809   CZ2 <NA>   TRP     B     6   <NA> 13.914  44.211   2.785 1
    809    ATOM   810   CZ3 <NA>   TRP     B     6   <NA> 11.903  43.585   3.953 1
    810    ATOM   811   CH2 <NA>   TRP     B     6   <NA> 13.148  44.273   3.896 1
    811    ATOM   812     N <NA>   GLN     B     7   <NA> 10.396  38.008  -0.958 1
    812    ATOM   813    CA <NA>   GLN     B     7   <NA>  9.518  36.960  -1.516 1
    813    ATOM   814     C <NA>   GLN     B     7   <NA> 10.321  35.670  -1.433 1
    814    ATOM   815     O <NA>   GLN     B     7   <NA> 11.546  35.694  -1.298 1
    815    ATOM   816    CB <NA>   GLN     B     7   <NA>  9.215  37.238  -3.000 1
    816    ATOM   817    CG <NA>   GLN     B     7   <NA>  8.278  38.433  -3.326 1
    817    ATOM   818    CD <NA>   GLN     B     7   <NA>  7.629  38.332  -4.744 1
    818    ATOM   819   OE1 <NA>   GLN     B     7   <NA>  8.319  38.100  -5.754 1
    819    ATOM   820   NE2 <NA>   GLN     B     7   <NA>  6.307  38.532  -4.814 1
    820    ATOM   821     N <NA>   ARG     B     8   <NA>  9.661  34.528  -1.525 1
    821    ATOM   822    CA <NA>   ARG     B     8   <NA> 10.385  33.251  -1.642 1
    822    ATOM   823     C <NA>   ARG     B     8   <NA> 11.348  33.305  -2.780 1
    823    ATOM   824     O <NA>   ARG     B     8   <NA> 10.964  33.682  -3.906 1
    824    ATOM   825    CB <NA>   ARG     B     8   <NA>  9.435  32.061  -1.862 1
    825    ATOM   826    CG <NA>   ARG     B     8   <NA>  8.623  31.716  -0.644 1
    826    ATOM   827    CD <NA>   ARG     B     8   <NA>  7.828  30.471  -0.844 1
    827    ATOM   828    NE <NA>   ARG     B     8   <NA>  7.143  30.102   0.378 1
    828    ATOM   829    CZ <NA>   ARG     B     8   <NA>  6.523  28.938   0.561 1
    829    ATOM   830   NH1 <NA>   ARG     B     8   <NA>  6.476  28.028  -0.411 1
    830    ATOM   831   NH2 <NA>   ARG     B     8   <NA>  5.925  28.688   1.719 1
    831    ATOM   832     N <NA>   PRO     B     9   <NA> 12.533  32.711  -2.598 1
    832    ATOM   833    CA <NA>   PRO     B     9   <NA> 13.443  32.743  -3.746 1
    833    ATOM   834     C <NA>   PRO     B     9   <NA> 13.174  31.592  -4.703 1
    834    ATOM   835     O <NA>   PRO     B     9   <NA> 13.897  30.587  -4.698 1
    835    ATOM   836    CB <NA>   PRO     B     9   <NA> 14.813  32.658  -3.125 1
    836    ATOM   837    CG <NA>   PRO     B     9   <NA> 14.564  31.823  -1.917 1
    837    ATOM   838    CD <NA>   PRO     B     9   <NA> 13.238  32.289  -1.373 1
    838    ATOM   839     N <NA>   LEU     B    10   <NA> 12.134  31.727  -5.504 1
    839    ATOM   840    CA <NA>   LEU     B    10   <NA> 11.816  30.740  -6.534 1
    840    ATOM   841     C <NA>   LEU     B    10   <NA> 12.459  31.075  -7.877 1
    841    ATOM   842     O <NA>   LEU     B    10   <NA> 12.274  32.150  -8.406 1
    842    ATOM   843    CB <NA>   LEU     B    10   <NA> 10.303  30.637  -6.738 1
    843    ATOM   844    CG <NA>   LEU     B    10   <NA>  9.483  30.307  -5.497 1
    844    ATOM   845   CD1 <NA>   LEU     B    10   <NA>  8.028  30.334  -5.876 1
    845    ATOM   846   CD2 <NA>   LEU     B    10   <NA>  9.845  28.975  -4.951 1
    846    ATOM   847     N <NA>   VAL     B    11   <NA> 13.225  30.141  -8.420 1
    847    ATOM   848    CA <NA>   VAL     B    11   <NA> 13.759  30.227  -9.768 1
    848    ATOM   849     C <NA>   VAL     B    11   <NA> 13.103  29.153 -10.641 1
    849    ATOM   850     O <NA>   VAL     B    11   <NA> 12.381  28.285 -10.135 1
    850    ATOM   851    CB <NA>   VAL     B    11   <NA> 15.253  29.988  -9.735 1
    851    ATOM   852   CG1 <NA>   VAL     B    11   <NA> 15.898  31.108  -8.939 1
    852    ATOM   853   CG2 <NA>   VAL     B    11   <NA> 15.573  28.635  -9.104 1
    853    ATOM   854     N <NA>   THR     B    12   <NA> 13.346  29.214 -11.949 1
    854    ATOM   855    CA <NA>   THR     B    12   <NA> 12.809  28.220 -12.873 1
    855    ATOM   856     C <NA>   THR     B    12   <NA> 13.951  27.253 -13.210 1
    856    ATOM   857     O <NA>   THR     B    12   <NA> 15.089  27.656 -13.442 1
    857    ATOM   858    CB <NA>   THR     B    12   <NA> 12.259  28.907 -14.158 1
    858    ATOM   859   OG1 <NA>   THR     B    12   <NA> 11.693  27.933 -15.028 1
    859    ATOM   860   CG2 <NA>   THR     B    12   <NA> 13.341  29.639 -14.925 1
    860    ATOM   861     N <NA>   ILE     B    13   <NA> 13.684  25.961 -13.124 1
    861    ATOM   862    CA <NA>   ILE     B    13   <NA> 14.708  24.966 -13.446 1
    862    ATOM   863     C <NA>   ILE     B    13   <NA> 14.230  24.165 -14.639 1
    863    ATOM   864     O <NA>   ILE     B    13   <NA> 13.014  24.057 -14.918 1
    864    ATOM   865    CB <NA>   ILE     B    13   <NA> 14.993  23.953 -12.269 1
    865    ATOM   866   CG1 <NA>   ILE     B    13   <NA> 13.699  23.190 -11.870 1
    866    ATOM   867   CG2 <NA>   ILE     B    13   <NA> 15.564  24.705 -11.094 1
    867    ATOM   868   CD1 <NA>   ILE     B    13   <NA> 13.900  22.077 -10.834 1
    868    ATOM   869     N <NA>   LYS     B    14   <NA> 15.186  23.630 -15.378 1
    869    ATOM   870    CA <NA>   LYS     B    14   <NA> 14.828  22.733 -16.464 1
    870    ATOM   871     C <NA>   LYS     B    14   <NA> 15.482  21.387 -16.141 1
    871    ATOM   872     O <NA>   LYS     B    14   <NA> 16.690  21.305 -15.886 1
    872    ATOM   873    CB <NA>   LYS     B    14   <NA> 15.340  23.256 -17.814 1
    873    ATOM   874    CG <NA>   LYS     B    14   <NA> 14.868  22.447 -18.992 1
    874    ATOM   875    CD <NA>   LYS     B    14   <NA> 14.687  23.295 -20.194 1
    875    ATOM   876    CE <NA>   LYS     B    14   <NA> 15.979  23.453 -20.922 1
    876    ATOM   877    NZ <NA>   LYS     B    14   <NA> 15.739  23.451 -22.389 1
    877    ATOM   878     N <NA>   ILE     B    15   <NA> 14.660  20.354 -16.136 1
    878    ATOM   879    CA <NA>   ILE     B    15   <NA> 15.108  18.999 -15.906 1
    879    ATOM   880     C <NA>   ILE     B    15   <NA> 14.188  18.067 -16.702 1
    880    ATOM   881     O <NA>   ILE     B    15   <NA> 12.968  18.199 -16.719 1
    881    ATOM   882    CB <NA>   ILE     B    15   <NA> 15.090  18.663 -14.360 1
    882    ATOM   883   CG1 <NA>   ILE     B    15   <NA> 15.694  17.283 -14.101 1
    883    ATOM   884   CG2 <NA>   ILE     B    15   <NA> 13.682  18.760 -13.780 1
    884    ATOM   885   CD1 <NA>   ILE     B    15   <NA> 16.011  17.036 -12.625 1
    885    ATOM   886     N <NA>   GLY     B    16   <NA> 14.799  17.159 -17.438 1
    886    ATOM   887    CA <NA>   GLY     B    16   <NA> 14.024  16.279 -18.286 1
    887    ATOM   888     C <NA>   GLY     B    16   <NA> 13.311  17.029 -19.391 1
    888    ATOM   889     O <NA>   GLY     B    16   <NA> 12.318  16.549 -19.924 1
    889    ATOM   890     N <NA>   GLY     B    17   <NA> 13.887  18.145 -19.823 1
    890    ATOM   891    CA <NA>   GLY     B    17   <NA> 13.243  18.938 -20.850 1
    891    ATOM   892     C <NA>   GLY     B    17   <NA> 12.009  19.638 -20.345 1
    892    ATOM   893     O <NA>   GLY     B    17   <NA> 11.392  20.369 -21.084 1
    893    ATOM   894     N <NA>   GLN     B    18   <NA> 11.676  19.486 -19.073 1
    894    ATOM   895    CA <NA>   GLN     B    18   <NA> 10.572  20.249 -18.515 1
    895    ATOM   896     C <NA>   GLN     B    18   <NA> 11.056  21.466 -17.699 1
    896    ATOM   897     O <NA>   GLN     B    18   <NA> 12.196  21.494 -17.175 1
    897    ATOM   898    CB <NA>   GLN     B    18   <NA>  9.708  19.388 -17.619 1
    898    ATOM   899    CG <NA>   GLN     B    18   <NA>  8.936  18.297 -18.309 1
    899    ATOM   900    CD <NA>   GLN     B    18   <NA>  9.088  16.973 -17.566 1
    900    ATOM   901   OE1 <NA>   GLN     B    18   <NA>  9.813  16.078 -18.014 1
    901    ATOM   902   NE2 <NA>   GLN     B    18   <NA>  8.488  16.884 -16.371 1
    902    ATOM   903     N <NA>   LEU     B    19   <NA> 10.186  22.475 -17.590 1
    903    ATOM   904    CA <NA>   LEU     B    19   <NA> 10.450  23.647 -16.749 1
    904    ATOM   905     C <NA>   LEU     B    19   <NA>  9.631  23.511 -15.483 1
    905    ATOM   906     O <NA>   LEU     B    19   <NA>  8.432  23.277 -15.551 1
    906    ATOM   907    CB <NA>   LEU     B    19   <NA> 10.010  24.932 -17.474 1
    907    ATOM   908    CG <NA>   LEU     B    19   <NA> 10.775  25.419 -18.723 1
    908    ATOM   909   CD1 <NA>   LEU     B    19   <NA> 10.226  26.727 -19.190 1
    909    ATOM   910   CD2 <NA>   LEU     B    19   <NA> 12.241  25.543 -18.395 1
    910    ATOM   911     N <NA>   LYS     B    20   <NA> 10.273  23.619 -14.329 1
    911    ATOM   912    CA <NA>   LYS     B    20   <NA>  9.576  23.583 -13.044 1
    912    ATOM   913     C <NA>   LYS     B    20   <NA> 10.018  24.784 -12.260 1
    913    ATOM   914     O <NA>   LYS     B    20   <NA> 10.998  25.429 -12.605 1
    914    ATOM   915    CB <NA>   LYS     B    20   <NA>  9.970  22.364 -12.236 1
    915    ATOM   916    CG <NA>   LYS     B    20   <NA> 10.281  21.123 -13.051 1
    916    ATOM   917    CD <NA>   LYS     B    20   <NA>  9.037  20.289 -13.305 1
    917    ATOM   918    CE <NA>   LYS     B    20   <NA>  9.400  18.836 -13.595 1
    918    ATOM   919    NZ <NA>   LYS     B    20   <NA>  9.672  18.102 -12.328 1
    919    ATOM   920     N <NA>   GLU     B    21   <NA>  9.324  25.048 -11.162 1
    920    ATOM   921    CA <NA>   GLU     B    21   <NA>  9.705  26.075 -10.199 1
    921    ATOM   922     C <NA>   GLU     B    21   <NA> 10.311  25.408  -8.990 1
    922    ATOM   923     O <NA>   GLU     B    21   <NA>  9.847  24.351  -8.547 1
    923    ATOM   924    CB <NA>   GLU     B    21   <NA>  8.468  26.848  -9.767 1
    924    ATOM   925    CG <NA>   GLU     B    21   <NA>  8.676  28.343  -9.687 1
    925    ATOM   926    CD <NA>   GLU     B    21   <NA>  7.419  29.111 -10.041 1
    926    ATOM   927   OE1 <NA>   GLU     B    21   <NA>  6.374  28.877  -9.383 1
    927    ATOM   928   OE2 <NA>   GLU     B    21   <NA>  7.461  29.911 -11.004 1
    928    ATOM   929     N <NA>   ALA     B    22   <NA> 11.348  26.008  -8.443 1
    929    ATOM   930    CA <NA>   ALA     B    22   <NA> 12.025  25.431  -7.271 1
    930    ATOM   931     C <NA>   ALA     B    22   <NA> 12.472  26.569  -6.331 1
    931    ATOM   932     O <NA>   ALA     B    22   <NA> 12.709  27.701  -6.765 1
    932    ATOM   933    CB <NA>   ALA     B    22   <NA> 13.248  24.643  -7.708 1
    933    ATOM   934     N <NA>   LEU     B    23   <NA> 12.560  26.257  -5.054 1
    934    ATOM   935    CA <NA>   LEU     B    23   <NA> 13.017  27.151  -4.001 1
    935    ATOM   936     C <NA>   LEU     B    23   <NA> 14.518  27.005  -3.822 1
    936    ATOM   937     O <NA>   LEU     B    23   <NA> 15.008  25.915  -3.588 1
    937    ATOM   938    CB <NA>   LEU     B    23   <NA> 12.330  26.721  -2.741 1
    938    ATOM   939    CG <NA>   LEU     B    23   <NA> 12.592  27.472  -1.469 1
    939    ATOM   940   CD1 <NA>   LEU     B    23   <NA> 12.006  28.877  -1.562 1
    940    ATOM   941   CD2 <NA>   LEU     B    23   <NA> 11.917  26.658  -0.379 1
    941    ATOM   942     N <NA>   LEU     B    24   <NA> 15.266  28.090  -3.963 1
    942    ATOM   943    CA <NA>   LEU     B    24   <NA> 16.672  28.114  -3.511 1
    943    ATOM   944     C <NA>   LEU     B    24   <NA> 16.809  28.098  -1.992 1
    944    ATOM   945     O <NA>   LEU     B    24   <NA> 16.417  29.019  -1.338 1
    945    ATOM   946    CB <NA>   LEU     B    24   <NA> 17.416  29.342  -4.065 1
    946    ATOM   947    CG <NA>   LEU     B    24   <NA> 17.444  29.528  -5.585 1
    947    ATOM   948   CD1 <NA>   LEU     B    24   <NA> 18.151  30.843  -5.884 1
    948    ATOM   949   CD2 <NA>   LEU     B    24   <NA> 18.170  28.385  -6.270 1
    949    ATOM   950     N <NA>   ASP     B    25   <NA> 17.407  27.054  -1.437 1
    950    ATOM   951    CA <NA>   ASP     B    25   <NA> 17.227  26.751  -0.026 1
    951    ATOM   952     C <NA>   ASP     B    25   <NA> 18.555  26.446   0.653 1
    952    ATOM   953     O <NA>   ASP     B    25   <NA> 19.003  25.309   0.692 1
    953    ATOM   954    CB <NA>   ASP     B    25   <NA> 16.258  25.572   0.084 1
    954    ATOM   955    CG <NA>   ASP     B    25   <NA> 15.759  25.336   1.493 1
    955    ATOM   956   OD1 <NA>   ASP     B    25   <NA> 16.399  25.780   2.453 1
    956    ATOM   957   OD2 <NA>   ASP     B    25   <NA> 14.731  24.675   1.645 1
    957    ATOM   958     N <NA>   THR     B    26   <NA> 19.163  27.455   1.257 1
    958    ATOM   959    CA <NA>   THR     B    26   <NA> 20.441  27.290   1.920 1
    959    ATOM   960     C <NA>   THR     B    26   <NA> 20.319  26.393   3.168 1
    960    ATOM   961     O <NA>   THR     B    26   <NA> 21.316  25.867   3.637 1
    961    ATOM   962    CB <NA>   THR     B    26   <NA> 21.063  28.678   2.282 1
    962    ATOM   963   OG1 <NA>   THR     B    26   <NA> 20.188  29.407   3.146 1
    963    ATOM   964   CG2 <NA>   THR     B    26   <NA> 21.279  29.499   1.024 1
    964    ATOM   965     N <NA>   GLY     B    27   <NA> 19.106  26.199   3.688 1
    965    ATOM   966    CA <NA>   GLY     B    27   <NA> 18.957  25.372   4.876 1
    966    ATOM   967     C <NA>   GLY     B    27   <NA> 18.845  23.903   4.520 1
    967    ATOM   968     O <NA>   GLY     B    27   <NA> 18.660  23.054   5.417 1
    968    ATOM   969     N <NA>   ALA     B    28   <NA> 18.819  23.600   3.217 1
    969    ATOM   970    CA <NA>   ALA     B    28   <NA> 18.721  22.211   2.738 1
    970    ATOM   971     C <NA>   ALA     B    28   <NA> 20.089  21.661   2.292 1
    971    ATOM   972     O <NA>   ALA     B    28   <NA> 20.749  22.243   1.403 1
    972    ATOM   973    CB <NA>   ALA     B    28   <NA> 17.682  22.117   1.578 1
    973    ATOM   974     N <NA>   ASP     B    29   <NA> 20.536  20.559   2.918 1
    974    ATOM   975    CA <NA>   ASP     B    29   <NA> 21.779  19.912   2.496 1
    975    ATOM   976     C <NA>   ASP     B    29   <NA> 21.693  19.374   1.107 1
    976    ATOM   977     O <NA>   ASP     B    29   <NA> 22.642  19.502   0.361 1
    977    ATOM   978    CB <NA>   ASP     B    29   <NA> 22.169  18.827   3.447 1
    978    ATOM   979    CG <NA>   ASP     B    29   <NA> 22.272  19.337   4.844 1
    979    ATOM   980   OD1 <NA>   ASP     B    29   <NA> 22.714  20.482   5.029 1
    980    ATOM   981   OD2 <NA>   ASP     B    29   <NA> 21.836  18.647   5.778 1
    981    ATOM   982     N <NA>   ASP     B    30   <NA> 20.524  18.868   0.719 1
    982    ATOM   983    CA <NA>   ASP     B    30   <NA> 20.355  18.227  -0.584 1
    983    ATOM   984     C <NA>   ASP     B    30   <NA> 19.212  18.807  -1.371 1
    984    ATOM   985     O <NA>   ASP     B    30   <NA> 18.383  19.523  -0.846 1
    985    ATOM   986    CB <NA>   ASP     B    30   <NA> 20.084  16.745  -0.413 1
    986    ATOM   987    CG <NA>   ASP     B    30   <NA> 21.023  16.088   0.586 1
    987    ATOM   988   OD1 <NA>   ASP     B    30   <NA> 22.233  15.970   0.265 1
    988    ATOM   989   OD2 <NA>   ASP     B    30   <NA> 20.531  15.682   1.672 1
    989    ATOM   990     N <NA>   THR     B    31   <NA> 19.152  18.428  -2.643 1
    990    ATOM   991    CA <NA>   THR     B    31   <NA> 18.113  18.843  -3.603 1
    991    ATOM   992     C <NA>   THR     B    31   <NA> 17.019  17.743  -3.682 1
    992    ATOM   993     O <NA>   THR     B    31   <NA> 17.342  16.586  -3.904 1
    993    ATOM   994    CB <NA>   THR     B    31   <NA> 18.810  19.051  -4.967 1
    994    ATOM   995   OG1 <NA>   THR     B    31   <NA> 19.740  20.123  -4.831 1
    995    ATOM   996   CG2 <NA>   THR     B    31   <NA> 17.844  19.330  -6.078 1
    996    ATOM   997     N <NA>   VAL     B    32   <NA> 15.750  18.102  -3.440 1
    997    ATOM   998    CA <NA>   VAL     B    32   <NA> 14.628  17.162  -3.514 1
    998    ATOM   999     C <NA>   VAL     B    32   <NA> 13.618  17.726  -4.436 1
    999    ATOM  1000     O <NA>   VAL     B    32   <NA> 13.169  18.861  -4.263 1
    1000   ATOM  1001    CB <NA>   VAL     B    32   <NA> 13.781  17.005  -2.245 1
    1001   ATOM  1002   CG1 <NA>   VAL     B    32   <NA> 13.297  15.592  -2.184 1
    1002   ATOM  1003   CG2 <NA>   VAL     B    32   <NA> 14.518  17.455  -1.007 1
    1003   ATOM  1004     N <NA>   LEU     B    33   <NA> 13.199  16.926  -5.401 1
    1004   ATOM  1005    CA <NA>   LEU     B    33   <NA> 12.141  17.335  -6.327 1
    1005   ATOM  1006     C <NA>   LEU     B    33   <NA> 10.876  16.500  -6.065 1
    1006   ATOM  1007     O <NA>   LEU     B    33   <NA> 10.948  15.389  -5.534 1
    1007   ATOM  1008    CB <NA>   LEU     B    33   <NA> 12.618  17.139  -7.766 1
    1008   ATOM  1009    CG <NA>   LEU     B    33   <NA> 13.889  17.846  -8.247 1
    1009   ATOM  1010   CD1 <NA>   LEU     B    33   <NA> 13.942  17.794  -9.731 1
    1010   ATOM  1011   CD2 <NA>   LEU     B    33   <NA> 13.897  19.266  -7.782 1
    1011   ATOM  1012     N <NA>   GLU     B    34   <NA>  9.719  17.083  -6.361 1
    1012   ATOM  1013    CA <NA>   GLU     B    34   <NA>  8.442  16.392  -6.346 1
    1013   ATOM  1014     C <NA>   GLU     B    34   <NA>  8.514  15.172  -7.224 1
    1014   ATOM  1015     O <NA>   GLU     B    34   <NA>  9.413  15.013  -8.040 1
    1015   ATOM  1016    CB <NA>   GLU     B    34   <NA>  7.316  17.305  -6.819 1
    1016   ATOM  1017    CG <NA>   GLU     B    34   <NA>  6.914  18.362  -5.808 1
    1017   ATOM  1018    CD <NA>   GLU     B    34   <NA>  6.205  19.552  -6.439 1
    1018   ATOM  1019   OE1 <NA>   GLU     B    34   <NA>  6.323  19.742  -7.666 1
    1019   ATOM  1020   OE2 <NA>   GLU     B    34   <NA>  5.613  20.369  -5.715 1
    1020   ATOM  1021     N <NA>   GLU     B    35   <NA>  7.526  14.309  -7.044 1
    1021   ATOM  1022    CA <NA>   GLU     B    35   <NA>  7.425  13.006  -7.682 1
    1022   ATOM  1023     C <NA>   GLU     B    35   <NA>  7.528  13.141  -9.172 1
    1023   ATOM  1024     O <NA>   GLU     B    35   <NA>  6.660  13.711  -9.819 1
    1024   ATOM  1025    CB <NA>   GLU     B    35   <NA>  6.100  12.361  -7.297 1
    1025   ATOM  1026    CG <NA>   GLU     B    35   <NA>  5.907  10.953  -7.838 1
    1026   ATOM  1027    CD <NA>   GLU     B    35   <NA>  7.182  10.093  -7.854 1
    1027   ATOM  1028   OE1 <NA>   GLU     B    35   <NA>  7.743   9.853  -6.766 1
    1028   ATOM  1029   OE2 <NA>   GLU     B    35   <NA>  7.521   9.561  -8.946 1
    1029   ATOM  1030     N <NA>   MET     B    36   <NA>  8.627  12.651  -9.705 1
    1030   ATOM  1031    CA <NA>   MET     B    36   <NA>  8.791  12.578 -11.145 1
    1031   ATOM  1032     C <NA>   MET     B    36   <NA>  9.583  11.322 -11.483 1
    1032   ATOM  1033     O <NA>   MET     B    36   <NA> 10.100  10.637 -10.616 1
    1033   ATOM  1034    CB <NA>   MET     B    36   <NA>  9.546  13.808 -11.654 1
    1034   ATOM  1035    CG <NA>   MET     B    36   <NA> 10.867  14.095 -11.014 1
    1035   ATOM  1036    SD <NA>   MET     B    36   <NA> 11.575  15.547 -11.778 1
    1036   ATOM  1037    CE <NA>   MET     B    36   <NA> 11.710  15.108 -13.551 1
    1037   ATOM  1038     N <NA>   SER     B    37   <NA>  9.657  11.016 -12.763 1
    1038   ATOM  1039    CA <NA>   SER     B    37   <NA> 10.411   9.858 -13.218 1
    1039   ATOM  1040     C <NA>   SER     B    37   <NA> 11.673  10.374 -13.825 1
    1040   ATOM  1041     O <NA>   SER     B    37   <NA> 11.636  11.272 -14.685 1
    1041   ATOM  1042    CB <NA>   SER     B    37   <NA>  9.573   9.049 -14.221 1
    1042   ATOM  1043    OG <NA>   SER     B    37   <NA>  8.330   8.594 -13.579 1
    1043   ATOM  1044     N <NA>   LEU     B    38   <NA> 12.793   9.884 -13.330 1
    1044   ATOM  1045    CA <NA>   LEU     B    38   <NA> 14.091  10.261 -13.857 1
    1045   ATOM  1046     C <NA>   LEU     B    38   <NA> 14.818   9.018 -14.330 1
    1046   ATOM  1047     O <NA>   LEU     B    38   <NA> 14.416   7.898 -13.995 1
    1047   ATOM  1048    CB <NA>   LEU     B    38   <NA> 14.866  10.942 -12.759 1
    1048   ATOM  1049    CG <NA>   LEU     B    38   <NA> 14.480  12.376 -12.556 1
    1049   ATOM  1050   CD1 <NA>   LEU     B    38   <NA> 15.159  12.900 -11.300 1
    1050   ATOM  1051   CD2 <NA>   LEU     B    38   <NA> 14.955  13.131 -13.766 1
    1051   ATOM  1052     N <NA>   PRO     B    39   <NA> 15.767   9.161 -15.261 1
    1052   ATOM  1053    CA <NA>   PRO     B    39   <NA> 16.525   8.031 -15.798 1
    1053   ATOM  1054     C <NA>   PRO     B    39   <NA> 17.366   7.241 -14.777 1
    1054   ATOM  1055     O <NA>   PRO     B    39   <NA> 17.943   7.817 -13.847 1
    1055   ATOM  1056    CB <NA>   PRO     B    39   <NA> 17.407   8.673 -16.857 1
    1056   ATOM  1057    CG <NA>   PRO     B    39   <NA> 17.532  10.100 -16.423 1
    1057   ATOM  1058    CD <NA>   PRO     B    39   <NA> 16.150  10.405 -15.968 1
    1058   ATOM  1059     N <NA>   GLY     B    40   <NA> 17.477   5.932 -14.999 1
    1059   ATOM  1060    CA <NA>   GLY     B    40   <NA> 18.494   5.171 -14.302 1
    1060   ATOM  1061     C <NA>   GLY     B    40   <NA> 18.048   4.556 -12.995 1
    1061   ATOM  1062     O <NA>   GLY     B    40   <NA> 16.865   4.438 -12.729 1
    1062   ATOM  1063     N <NA>   ARG     B    41   <NA> 19.000   3.939 -12.313 1
    1063   ATOM  1064    CA <NA>   ARG     B    41   <NA> 18.722   3.282 -11.042 1
    1064   ATOM  1065     C <NA>   ARG     B    41   <NA> 18.615   4.306  -9.916 1
    1065   ATOM  1066     O <NA>   ARG     B    41   <NA> 19.168   5.390 -10.018 1
    1066   ATOM  1067    CB <NA>   ARG     B    41   <NA> 19.852   2.333 -10.662 1
    1067   ATOM  1068    CG <NA>   ARG     B    41   <NA> 20.359   1.417 -11.726 1
    1068   ATOM  1069    CD <NA>   ARG     B    41   <NA> 20.905   0.188 -11.023 1
    1069   ATOM  1070    NE <NA>   ARG     B    41   <NA> 19.927  -0.864 -11.169 1
    1070   ATOM  1071    CZ <NA>   ARG     B    41   <NA> 19.354  -1.538 -10.177 1
    1071   ATOM  1072   NH1 <NA>   ARG     B    41   <NA> 19.839  -1.513  -8.939 1
    1072   ATOM  1073   NH2 <NA>   ARG     B    41   <NA> 18.333  -2.314 -10.472 1
    1073   ATOM  1074     N <NA>   TRP     B    42   <NA> 17.989   3.918  -8.810 1
    1074   ATOM  1075    CA <NA>   TRP     B    42   <NA> 17.920   4.766  -7.634 1
    1075   ATOM  1076     C <NA>   TRP     B    42   <NA> 18.295   3.969  -6.380 1
    1076   ATOM  1077     O <NA>   TRP     B    42   <NA> 18.163   2.742  -6.361 1
    1077   ATOM  1078    CB <NA>   TRP     B    42   <NA> 16.525   5.338  -7.459 1
    1078   ATOM  1079    CG <NA>   TRP     B    42   <NA> 15.444   4.372  -7.312 1
    1079   ATOM  1080   CD1 <NA>   TRP     B    42   <NA> 14.681   3.834  -8.299 1
    1080   ATOM  1081   CD2 <NA>   TRP     B    42   <NA> 14.840   3.957  -6.079 1
    1081   ATOM  1082   NE1 <NA>   TRP     B    42   <NA> 13.640   3.121  -7.756 1
    1082   ATOM  1083   CE2 <NA>   TRP     B    42   <NA> 13.719   3.166  -6.402 1
    1083   ATOM  1084   CE3 <NA>   TRP     B    42   <NA> 15.154   4.180  -4.715 1
    1084   ATOM  1085   CZ2 <NA>   TRP     B    42   <NA> 12.881   2.589  -5.411 1
    1085   ATOM  1086   CZ3 <NA>   TRP     B    42   <NA> 14.300   3.625  -3.745 1
    1086   ATOM  1087   CH2 <NA>   TRP     B    42   <NA> 13.168   2.842  -4.106 1
    1087   ATOM  1088     N <NA>   LYS     B    43   <NA> 18.801   4.689  -5.365 1
    1088   ATOM  1089    CA <NA>   LYS     B    43   <NA> 19.180   4.182  -4.032 1
    1089   ATOM  1090     C <NA>   LYS     B    43   <NA> 18.127   4.736  -3.089 1
    1090   ATOM  1091     O <NA>   LYS     B    43   <NA> 17.442   5.729  -3.400 1
    1091   ATOM  1092    CB <NA>   LYS     B    43   <NA> 20.561   4.731  -3.581 1
    1092   ATOM  1093    CG <NA>   LYS     B    43   <NA> 21.777   4.400  -4.445 1
    1093   ATOM  1094    CD <NA>   LYS     B    43   <NA> 22.996   5.295  -4.048 1
    1094   ATOM  1095    CE <NA>   LYS     B    43   <NA> 24.193   5.280  -5.063 1
    1095   ATOM  1096    NZ <NA>   LYS     B    43   <NA> 25.251   6.324  -4.725 1
    1096   ATOM  1097     N <NA>   PRO     B    44   <NA> 18.053   4.208  -1.878 1
    1097   ATOM  1098    CA <NA>   PRO     B    44   <NA> 17.102   4.804  -0.946 1
    1098   ATOM  1099     C <NA>   PRO     B    44   <NA> 17.754   5.853  -0.023 1
    1099   ATOM  1100     O <NA>   PRO     B    44   <NA> 18.929   5.769   0.330 1
    1100   ATOM  1101    CB <NA>   PRO     B    44   <NA> 16.596   3.610  -0.171 1
    1101   ATOM  1102    CG <NA>   PRO     B    44   <NA> 17.803   2.676  -0.117 1
    1102   ATOM  1103    CD <NA>   PRO     B    44   <NA> 18.649   2.962  -1.335 1
    1103   ATOM  1104     N <NA>   LYS     B    45   <NA> 16.974   6.847   0.381 1
    1104   ATOM  1105    CA <NA>   LYS     B    45   <NA> 17.443   7.812   1.371 1
    1105   ATOM  1106     C <NA>   LYS     B    45   <NA> 16.334   8.328   2.257 1
    1106   ATOM  1107     O <NA>   LYS     B    45   <NA> 15.192   8.470   1.828 1
    1107   ATOM  1108    CB <NA>   LYS     B    45   <NA> 18.177   8.988   0.701 1
    1108   ATOM  1109    CG <NA>   LYS     B    45   <NA> 19.183   9.659   1.670 1
    1109   ATOM  1110    CD <NA>   LYS     B    45   <NA> 20.095  10.640   1.011 1
    1110   ATOM  1111    CE <NA>   LYS     B    45   <NA> 20.751  11.491   2.069 1
    1111   ATOM  1112    NZ <NA>   LYS     B    45   <NA> 21.413  12.705   1.509 1
    1112   ATOM  1113     N <NA>   MET     B    46   <NA> 16.672   8.585   3.514 1
    1113   ATOM  1114    CA <NA>   MET     B    46   <NA> 15.755   9.281   4.404 1
    1114   ATOM  1115     C <NA>   MET     B    46   <NA> 16.373  10.584   4.732 1
    1115   ATOM  1116     O <NA>   MET     B    46   <NA> 17.542  10.636   5.104 1
    1116   ATOM  1117    CB <NA>   MET     B    46   <NA> 15.562   8.530   5.708 1
    1117   ATOM  1118    CG <NA>   MET     B    46   <NA> 14.763   7.266   5.540 1
    1118   ATOM  1119    SD <NA>   MET     B    46   <NA> 13.357   7.367   6.566 1
    1119   ATOM  1120    CE <NA>   MET     B    46   <NA> 14.146   6.922   8.167 1
    1120   ATOM  1121     N <NA>   ILE     B    47   <NA> 15.582  11.636   4.604 1
    1121   ATOM  1122    CA <NA>   ILE     B    47   <NA> 16.003  12.986   4.955 1
    1122   ATOM  1123     C <NA>   ILE     B    47   <NA> 15.018  13.507   5.961 1
    1123   ATOM  1124     O <NA>   ILE     B    47   <NA> 13.822  13.222   5.884 1
    1124   ATOM  1125    CB <NA>   ILE     B    47   <NA> 16.040  13.915   3.699 1
    1125   ATOM  1126   CG1 <NA>   ILE     B    47   <NA> 14.745  13.761   2.918 1
    1126   ATOM  1127   CG2 <NA>   ILE     B    47   <NA> 17.229  13.534   2.782 1
    1127   ATOM  1128   CD1 <NA>   ILE     B    47   <NA> 14.742  14.561   1.698 1
    1128   ATOM  1129     N <NA>   GLY     B    48   <NA> 15.542  14.239   6.941 1
    1129   ATOM  1130    CA <NA>   GLY     B    48   <NA> 14.714  14.799   8.014 1
    1130   ATOM  1131     C <NA>   GLY     B    48   <NA> 14.762  16.314   8.114 1
    1131   ATOM  1132     O <NA>   GLY     B    48   <NA> 15.803  16.952   7.888 1
    1132   ATOM  1133     N <NA>   GLY     B    49   <NA> 13.583  16.896   8.239 1
    1133   ATOM  1134    CA <NA>   GLY     B    49   <NA> 13.484  18.319   8.459 1
    1134   ATOM  1135     C <NA>   GLY     B    49   <NA> 12.647  18.565   9.672 1
    1135   ATOM  1136     O <NA>   GLY     B    49   <NA> 12.880  17.976  10.721 1
    1136   ATOM  1137     N <NA>   ILE     B    50   <NA> 11.850  19.611   9.584 1
    1137   ATOM  1138    CA <NA>   ILE     B    50   <NA> 10.858  19.974  10.594 1
    1138   ATOM  1139     C <NA>   ILE     B    50   <NA>  9.707  18.988  10.393 1
    1139   ATOM  1140     O <NA>   ILE     B    50   <NA>  9.341  18.670   9.251 1
    1140   ATOM  1141    CB <NA>   ILE     B    50   <NA> 10.365  21.453  10.316 1
    1141   ATOM  1142   CG1 <NA>   ILE     B    50   <NA> 11.556  22.421  10.394 1
    1142   ATOM  1143   CG2 <NA>   ILE     B    50   <NA>  9.299  21.862  11.270 1
    1143   ATOM  1144   CD1 <NA>   ILE     B    50   <NA> 11.936  22.850  11.773 1
    1144   ATOM  1145     N <NA>   GLY     B    51   <NA>  9.164  18.455  11.474 1
    1145   ATOM  1146    CA <NA>   GLY     B    51   <NA>  8.011  17.583  11.313 1
    1146   ATOM  1147     C <NA>   GLY     B    51   <NA>  8.360  16.130  11.054 1
    1147   ATOM  1148     O <NA>   GLY     B    51   <NA>  7.494  15.283  11.167 1
    1148   ATOM  1149     N <NA>   GLY     B    52   <NA>  9.638  15.842  10.818 1
    1149   ATOM  1150    CA <NA>   GLY     B    52   <NA> 10.123  14.474  10.792 1
    1150   ATOM  1151     C <NA>   GLY     B    52   <NA> 10.860  14.080   9.524 1
    1151   ATOM  1152     O <NA>   GLY     B    52   <NA> 11.419  14.930   8.826 1
    1152   ATOM  1153     N <NA>   PHE     B    53   <NA> 10.878  12.788   9.221 1
    1153   ATOM  1154    CA <NA>   PHE     B    53   <NA> 11.638  12.302   8.079 1
    1154   ATOM  1155     C <NA>   PHE     B    53   <NA> 10.739  11.914   6.924 1
    1155   ATOM  1156     O <NA>   PHE     B    53   <NA>  9.601  11.543   7.137 1
    1156   ATOM  1157    CB <NA>   PHE     B    53   <NA> 12.458  11.126   8.531 1
    1157   ATOM  1158    CG <NA>   PHE     B    53   <NA> 13.464  11.471   9.564 1
    1158   ATOM  1159   CD1 <NA>   PHE     B    53   <NA> 13.092  11.697  10.886 1
    1159   ATOM  1160   CD2 <NA>   PHE     B    53   <NA> 14.789  11.625   9.189 1
    1160   ATOM  1161   CE1 <NA>   PHE     B    53   <NA> 14.036  12.076  11.825 1
    1161   ATOM  1162   CE2 <NA>   PHE     B    53   <NA> 15.753  12.001  10.078 1
    1162   ATOM  1163    CZ <NA>   PHE     B    53   <NA> 15.392  12.225  11.421 1
    1163   ATOM  1164     N <NA>   ILE     B    54   <NA> 11.204  12.078   5.695 1
    1164   ATOM  1165    CA <NA>   ILE     B    54   <NA> 10.538  11.431   4.563 1
    1165   ATOM  1166     C <NA>   ILE     B    54   <NA> 11.513  10.453   3.866 1
    1166   ATOM  1167     O <NA>   ILE     B    54   <NA> 12.727  10.529   4.052 1
    1167   ATOM  1168    CB <NA>   ILE     B    54   <NA>  9.923  12.446   3.500 1
    1168   ATOM  1169   CG1 <NA>   ILE     B    54   <NA> 10.968  13.414   2.964 1
    1169   ATOM  1170   CG2 <NA>   ILE     B    54   <NA>  8.754  13.195   4.090 1
    1170   ATOM  1171   CD1 <NA>   ILE     B    54   <NA> 10.571  14.020   1.616 1
    1171   ATOM  1172     N <NA>   LYS     B    55   <NA> 10.983   9.503   3.111 1
    1172   ATOM  1173    CA <NA>   LYS     B    55   <NA> 11.816   8.478   2.482 1
    1173   ATOM  1174     C <NA>   LYS     B    55   <NA> 11.862   8.988   1.074 1
    1174   ATOM  1175     O <NA>   LYS     B    55   <NA> 10.827   9.256   0.525 1
    1175   ATOM  1176    CB <NA>   LYS     B    55   <NA> 11.062   7.136   2.489 1
    1176   ATOM  1177    CG <NA>   LYS     B    55   <NA> 11.699   5.963   3.273 1
    1177   ATOM  1178    CD <NA>   LYS     B    55   <NA> 13.070   5.502   2.689 1
    1178   ATOM  1179    CE <NA>   LYS     B    55   <NA> 12.949   4.923   1.253 1
    1179   ATOM  1180    NZ <NA>   LYS     B    55   <NA> 13.964   5.445   0.291 1
    1180   ATOM  1181     N <NA>   VAL     B    56   <NA> 13.024   9.137   0.474 1
    1181   ATOM  1182    CA <NA>   VAL     B    56   <NA> 13.072   9.652  -0.897 1
    1182   ATOM  1183     C <NA>   VAL     B    56   <NA> 13.885   8.719  -1.786 1
    1183   ATOM  1184     O <NA>   VAL     B    56   <NA> 14.547   7.817  -1.279 1
    1184   ATOM  1185    CB <NA>   VAL     B    56   <NA> 13.757  11.033  -0.959 1
    1185   ATOM  1186   CG1 <NA>   VAL     B    56   <NA> 12.766  12.134  -0.523 1
    1186   ATOM  1187   CG2 <NA>   VAL     B    56   <NA> 15.032  11.016  -0.119 1
    1187   ATOM  1188     N <NA>   ARG     B    57   <NA> 13.889   8.977  -3.102 1
    1188   ATOM  1189    CA <NA>   ARG     B    57   <NA> 14.697   8.220  -4.083 1
    1189   ATOM  1190     C <NA>   ARG     B    57   <NA> 15.894   9.023  -4.591 1
    1190   ATOM  1191     O <NA>   ARG     B    57   <NA> 15.735  10.066  -5.210 1
    1191   ATOM  1192    CB <NA>   ARG     B    57   <NA> 13.862   7.854  -5.298 1
    1192   ATOM  1193    CG <NA>   ARG     B    57   <NA> 12.767   6.845  -5.041 1
    1193   ATOM  1194    CD <NA>   ARG     B    57   <NA> 12.224   6.316  -6.354 1
    1194   ATOM  1195    NE <NA>   ARG     B    57   <NA> 10.944   6.937  -6.617 1
    1195   ATOM  1196    CZ <NA>   ARG     B    57   <NA> 10.717   7.776  -7.614 1
    1196   ATOM  1197   NH1 <NA>   ARG     B    57   <NA> 11.534   7.792  -8.658 1
    1197   ATOM  1198   NH2 <NA>   ARG     B    57   <NA>  9.555   8.403  -7.678 1
    1198   ATOM  1199     N <NA>   GLN     B    58   <NA> 17.095   8.516  -4.388 1
    1199   ATOM  1200    CA <NA>   GLN     B    58   <NA> 18.306   9.218  -4.819 1
    1200   ATOM  1201     C <NA>   GLN     B    58   <NA> 18.742   8.833  -6.226 1
    1201   ATOM  1202     O <NA>   GLN     B    58   <NA> 19.157   7.694  -6.438 1
    1202   ATOM  1203    CB <NA>   GLN     B    58   <NA> 19.465   8.920  -3.861 1
    1203   ATOM  1204    CG <NA>   GLN     B    58   <NA> 20.738   9.622  -4.271 1
    1204   ATOM  1205    CD <NA>   GLN     B    58   <NA> 21.825   9.463  -3.248 1
    1205   ATOM  1206   OE1 <NA>   GLN     B    58   <NA> 21.554   9.424  -2.048 1
    1206   ATOM  1207   NE2 <NA>   GLN     B    58   <NA> 23.045   9.365  -3.692 1
    1207   ATOM  1208     N <NA>   TYR     B    59   <NA> 18.705   9.773  -7.167 1
    1208   ATOM  1209    CA <NA>   TYR     B    59   <NA> 19.361   9.587  -8.464 1
    1209   ATOM  1210     C <NA>   TYR     B    59   <NA> 20.689  10.330  -8.544 1
    1210   ATOM  1211     O <NA>   TYR     B    59   <NA> 20.860  11.368  -7.943 1
    1211   ATOM  1212    CB <NA>   TYR     B    59   <NA> 18.472  10.082  -9.563 1
    1212   ATOM  1213    CG <NA>   TYR     B    59   <NA> 17.116   9.383  -9.609 1
    1213   ATOM  1214   CD1 <NA>   TYR     B    59   <NA> 16.157   9.612  -8.615 1
    1214   ATOM  1215   CD2 <NA>   TYR     B    59   <NA> 16.814   8.484 -10.637 1
    1215   ATOM  1216   CE1 <NA>   TYR     B    59   <NA> 14.959   8.977  -8.640 1
    1216   ATOM  1217   CE2 <NA>   TYR     B    59   <NA> 15.647   7.851 -10.673 1
    1217   ATOM  1218    CZ <NA>   TYR     B    59   <NA> 14.704   8.066  -9.679 1
    1218   ATOM  1219    OH <NA>   TYR     B    59   <NA> 13.561   7.307  -9.711 1
    1219   ATOM  1220     N <NA>   ASP     B    60   <NA> 21.665   9.797  -9.258 1
    1220   ATOM  1221    CA <NA>   ASP     B    60   <NA> 22.959  10.470  -9.336 1
    1221   ATOM  1222     C <NA>   ASP     B    60   <NA> 23.303  10.921 -10.737 1
    1222   ATOM  1223     O <NA>   ASP     B    60   <NA> 22.793  10.396 -11.707 1
    1223   ATOM  1224    CB <NA>   ASP     B    60   <NA> 24.042   9.554  -8.834 1
    1224   ATOM  1225    CG <NA>   ASP     B    60   <NA> 23.843   9.184  -7.407 1
    1225   ATOM  1226   OD1 <NA>   ASP     B    60   <NA> 23.463  10.074  -6.620 1
    1226   ATOM  1227   OD2 <NA>   ASP     B    60   <NA> 24.107   8.012  -7.044 1
    1227   ATOM  1228     N <NA>   GLN     B    61   <NA> 24.189  11.897 -10.837 1
    1228   ATOM  1229    CA <NA>   GLN     B    61   <NA> 24.622  12.478 -12.111 1
    1229   ATOM  1230     C <NA>   GLN     B    61   <NA> 23.474  12.903 -13.022 1
    1230   ATOM  1231     O <NA>   GLN     B    61   <NA> 23.492  12.665 -14.229 1
    1231   ATOM  1232    CB <NA>   GLN     B    61   <NA> 25.596  11.549 -12.869 1
    1232   ATOM  1233    CG <NA>   GLN     B    61   <NA> 26.892  12.283 -13.369 1
    1233   ATOM  1234    CD <NA>   GLN     B    61   <NA> 28.007  11.392 -13.985 1
    1234   ATOM  1235   OE1 <NA>   GLN     B    61   <NA> 28.747  11.845 -14.875 1
    1235   ATOM  1236   NE2 <NA>   GLN     B    61   <NA> 28.197  10.173 -13.442 1
    1236   ATOM  1237     N <NA>   ILE     B    62   <NA> 22.493  13.589 -12.452 1
    1237   ATOM  1238    CA <NA>   ILE     B    62   <NA> 21.380  14.141 -13.224 1
    1238   ATOM  1239     C <NA>   ILE     B    62   <NA> 21.710  15.575 -13.686 1
    1239   ATOM  1240     O <NA>   ILE     B    62   <NA> 22.247  16.373 -12.924 1
    1240   ATOM  1241    CB <NA>   ILE     B    62   <NA> 20.090  14.147 -12.331 1
    1241   ATOM  1242   CG1 <NA>   ILE     B    62   <NA> 19.709  12.702 -11.952 1
    1242   ATOM  1243   CG2 <NA>   ILE     B    62   <NA> 18.950  14.826 -13.045 1
    1243   ATOM  1244   CD1 <NA>   ILE     B    62   <NA> 19.068  11.933 -13.057 1
    1244   ATOM  1245     N <NA>   LEU     B    63   <NA> 21.400  15.900 -14.931 1
    1245   ATOM  1246    CA <NA>   LEU     B    63   <NA> 21.600  17.263 -15.431 1
    1246   ATOM  1247     C <NA>   LEU     B    63   <NA> 20.386  18.085 -15.052 1
    1247   ATOM  1248     O <NA>   LEU     B    63   <NA> 19.260  17.703 -15.355 1
    1248   ATOM  1249    CB <NA>   LEU     B    63   <NA> 21.769  17.265 -16.962 1
    1249   ATOM  1250    CG <NA>   LEU     B    63   <NA> 21.792  18.587 -17.759 1
    1250   ATOM  1251   CD1 <NA>   LEU     B    63   <NA> 22.903  19.529 -17.300 1
    1251   ATOM  1252   CD2 <NA>   LEU     B    63   <NA> 21.997  18.246 -19.205 1
    1252   ATOM  1253     N <NA>   ILE     B    64   <NA> 20.626  19.203 -14.381 1
    1253   ATOM  1254    CA <NA>   ILE     B    64   <NA> 19.548  20.111 -14.029 1
    1254   ATOM  1255     C <NA>   ILE     B    64   <NA> 20.089  21.523 -14.223 1
    1255   ATOM  1256     O <NA>   ILE     B    64   <NA> 21.175  21.858 -13.763 1
    1256   ATOM  1257    CB <NA>   ILE     B    64   <NA> 19.107  19.880 -12.540 1
    1257   ATOM  1258   CG1 <NA>   ILE     B    64   <NA> 18.216  20.990 -12.036 1
    1258   ATOM  1259   CG2 <NA>   ILE     B    64   <NA> 20.311  19.766 -11.655 1
    1259   ATOM  1260   CD1 <NA>   ILE     B    64   <NA> 17.324  20.544 -10.930 1
    1260   ATOM  1261     N <NA>   GLU     B    65   <NA> 19.327  22.330 -14.953 1
    1261   ATOM  1262    CA <NA>   GLU     B    65   <NA> 19.661  23.719 -15.251 1
    1262   ATOM  1263     C <NA>   GLU     B    65   <NA> 18.938  24.686 -14.338 1
    1263   ATOM  1264     O <NA>   GLU     B    65   <NA> 17.700  24.745 -14.345 1
    1264   ATOM  1265    CB <NA>   GLU     B    65   <NA> 19.282  24.017 -16.688 1
    1265   ATOM  1266    CG <NA>   GLU     B    65   <NA> 20.180  25.011 -17.326 1
    1266   ATOM  1267    CD <NA>   GLU     B    65   <NA> 19.960  25.126 -18.805 1
    1267   ATOM  1268   OE1 <NA>   GLU     B    65   <NA> 19.601  24.094 -19.445 1
    1268   ATOM  1269   OE2 <NA>   GLU     B    65   <NA> 20.214  26.241 -19.314 1
    1269   ATOM  1270     N <NA>   ILE     B    66   <NA> 19.709  25.427 -13.547 1
    1270   ATOM  1271    CA <NA>   ILE     B    66   <NA> 19.159  26.431 -12.654 1
    1271   ATOM  1272     C <NA>   ILE     B    66   <NA> 19.519  27.864 -13.127 1
    1272   ATOM  1273     O <NA>   ILE     B    66   <NA> 20.684  28.254 -13.167 1
    1273   ATOM  1274    CB <NA>   ILE     B    66   <NA> 19.663  26.199 -11.203 1
    1274   ATOM  1275   CG1 <NA>   ILE     B    66   <NA> 19.566  24.717 -10.848 1
    1275   ATOM  1276   CG2 <NA>   ILE     B    66   <NA> 18.824  27.018 -10.232 1
    1276   ATOM  1277   CD1 <NA>   ILE     B    66   <NA> 20.510  24.280  -9.745 1
    1277   ATOM  1278     N <NA>   CYS     B    67   <NA> 18.504  28.630 -13.516 1
    1278   ATOM  1279    CA <NA>   CYS     B    67   <NA> 18.684  29.971 -14.104 1
    1279   ATOM  1280     C <NA>   CYS     B    67   <NA> 19.685  29.990 -15.245 1
    1280   ATOM  1281     O <NA>   CYS     B    67   <NA> 20.565  30.852 -15.294 1
    1281   ATOM  1282    CB <NA>   CYS     B    67   <NA> 19.124  30.970 -13.037 1
    1282   ATOM  1283    SG <NA>   CYS     B    67   <NA> 17.736  31.542 -12.037 1
    1283   ATOM  1284     N <NA>   GLY     B    68   <NA> 19.562  29.010 -16.143 1
    1284   ATOM  1285    CA <NA>   GLY     B    68   <NA> 20.485  28.891 -17.256 1
    1285   ATOM  1286     C <NA>   GLY     B    68   <NA> 21.899  28.603 -16.811 1
    1286   ATOM  1287     O <NA>   GLY     B    68   <NA> 22.843  28.970 -17.493 1
    1287   ATOM  1288     N <NA>   HIS     B    69   <NA> 22.059  28.002 -15.633 1
    1288   ATOM  1289    CA <NA>   HIS     B    69   <NA> 23.354  27.479 -15.197 1
    1289   ATOM  1290     C <NA>   HIS     B    69   <NA> 23.178  25.984 -15.087 1
    1290   ATOM  1291     O <NA>   HIS     B    69   <NA> 22.307  25.523 -14.354 1
    1291   ATOM  1292    CB <NA>   HIS     B    69   <NA> 23.711  27.993 -13.810 1
    1292   ATOM  1293    CG <NA>   HIS     B    69   <NA> 23.976  29.473 -13.740 1
    1293   ATOM  1294   ND1 <NA>   HIS     B    69   <NA> 25.177  29.994 -13.361 1
    1294   ATOM  1295   CD2 <NA>   HIS     B    69   <NA> 23.114  30.512 -13.920 1
    1295   ATOM  1296   CE1 <NA>   HIS     B    69   <NA> 25.058  31.325 -13.286 1
    1296   ATOM  1297   NE2 <NA>   HIS     B    69   <NA> 23.849  31.652 -13.613 1
    1297   ATOM  1298     N <NA>   LYS     B    70   <NA> 23.995  25.240 -15.820 1
    1298   ATOM  1299    CA <NA>   LYS     B    70   <NA> 23.935  23.791 -15.800 1
    1299   ATOM  1300     C <NA>   LYS     B    70   <NA> 24.749  23.243 -14.652 1
    1300   ATOM  1301     O <NA>   LYS     B    70   <NA> 25.875  23.676 -14.403 1
    1301   ATOM  1302    CB <NA>   LYS     B    70   <NA> 24.423  23.212 -17.115 1
    1302   ATOM  1303    CG <NA>   LYS     B    70   <NA> 23.463  23.445 -18.279 1
    1303   ATOM  1304    CD <NA>   LYS     B    70   <NA> 24.261  23.524 -19.576 1
    1304   ATOM  1305    CE <NA>   LYS     B    70   <NA> 23.377  23.737 -20.816 1
    1305   ATOM  1306    NZ <NA>   LYS     B    70   <NA> 22.863  22.443 -21.404 1
    1306   ATOM  1307     N <NA>   ALA     B    71   <NA> 24.104  22.353 -13.909 1
    1307   ATOM  1308    CA <NA>   ALA     B    71   <NA> 24.689  21.612 -12.802 1
    1308   ATOM  1309     C <NA>   ALA     B    71   <NA> 24.391  20.156 -13.159 1
    1309   ATOM  1310     O <NA>   ALA     B    71   <NA> 23.339  19.865 -13.735 1
    1310   ATOM  1311    CB <NA>   ALA     B    71   <NA> 23.991  21.992 -11.485 1
    1311   ATOM  1312     N <NA>   ILE     B    72   <NA> 25.330  19.253 -12.902 1
    1312   ATOM  1313    CA <NA>   ILE     B    72   <NA> 25.048  17.816 -13.016 1
    1313   ATOM  1314     C <NA>   ILE     B    72   <NA> 25.312  17.246 -11.637 1
    1314   ATOM  1315     O <NA>   ILE     B    72   <NA> 26.442  17.315 -11.167 1
    1315   ATOM  1316    CB <NA>   ILE     B    72   <NA> 26.029  17.094 -13.983 1
    1316   ATOM  1317   CG1 <NA>   ILE     B    72   <NA> 26.092  17.805 -15.333 1
    1317   ATOM  1318   CG2 <NA>   ILE     B    72   <NA> 25.615  15.649 -14.171 1
    1318   ATOM  1319   CD1 <NA>   ILE     B    72   <NA> 27.241  17.224 -16.197 1
    1319   ATOM  1320     N <NA>   GLY     B    73   <NA> 24.303  16.690 -10.975 1
    1320   ATOM  1321    CA <NA>   GLY     B    73   <NA> 24.545  16.248  -9.616 1
    1321   ATOM  1322     C <NA>   GLY     B    73   <NA> 23.470  15.337  -9.095 1
    1322   ATOM  1323     O <NA>   GLY     B    73   <NA> 22.674  14.852  -9.881 1
    1323   ATOM  1324     N <NA>   THR     B    74   <NA> 23.517  15.041  -7.794 1
    1324   ATOM  1325    CA <NA>   THR     B    74   <NA> 22.568  14.162  -7.116 1
    1325   ATOM  1326     C <NA>   THR     B    74   <NA> 21.280  14.889  -6.829 1
    1326   ATOM  1327     O <NA>   THR     B    74   <NA> 21.299  16.005  -6.275 1
    1327   ATOM  1328    CB <NA>   THR     B    74   <NA> 23.128  13.638  -5.748 1
    1328   ATOM  1329   OG1 <NA>   THR     B    74   <NA> 24.323  12.880  -5.961 1
    1329   ATOM  1330   CG2 <NA>   THR     B    74   <NA> 22.071  12.763  -5.030 1
    1330   ATOM  1331     N <NA>   VAL     B    75   <NA> 20.177  14.227  -7.192 1
    1331   ATOM  1332    CA <NA>   VAL     B    75   <NA> 18.832  14.773  -7.068 1
    1332   ATOM  1333     C <NA>   VAL     B    75   <NA> 17.989  13.688  -6.360 1
    1333   ATOM  1334     O <NA>   VAL     B    75   <NA> 17.993  12.504  -6.743 1
    1334   ATOM  1335    CB <NA>   VAL     B    75   <NA> 18.229  15.113  -8.501 1
    1335   ATOM  1336   CG1 <NA>   VAL     B    75   <NA> 16.760  15.421  -8.434 1
    1336   ATOM  1337   CG2 <NA>   VAL     B    75   <NA> 18.967  16.317  -9.099 1
    1337   ATOM  1338     N <NA>   LEU     B    76   <NA> 17.295  14.086  -5.303 1
    1338   ATOM  1339    CA <NA>   LEU     B    76   <NA> 16.390  13.196  -4.575 1
    1339   ATOM  1340     C <NA>   LEU     B    76   <NA> 14.972  13.423  -5.073 1
    1340   ATOM  1341     O <NA>   LEU     B    76   <NA> 14.652  14.514  -5.518 1
    1341   ATOM  1342    CB <NA>   LEU     B    76   <NA> 16.450  13.476  -3.067 1
    1342   ATOM  1343    CG <NA>   LEU     B    76   <NA> 17.787  13.595  -2.330 1
    1343   ATOM  1344   CD1 <NA>   LEU     B    76   <NA> 17.591  13.640  -0.845 1
    1344   ATOM  1345   CD2 <NA>   LEU     B    76   <NA> 18.678  12.463  -2.705 1
    1345   ATOM  1346     N <NA>   VAL     B    77   <NA> 14.135  12.391  -5.049 1
    1346   ATOM  1347    CA <NA>   VAL     B    77   <NA> 12.749  12.509  -5.485 1
    1347   ATOM  1348     C <NA>   VAL     B    77   <NA> 11.826  11.894  -4.472 1
    1348   ATOM  1349     O <NA>   VAL     B    77   <NA> 12.052  10.766  -3.999 1
    1349   ATOM  1350    CB <NA>   VAL     B    77   <NA> 12.502  11.835  -6.868 1
    1350   ATOM  1351   CG1 <NA>   VAL     B    77   <NA> 11.065  12.148  -7.398 1
    1351   ATOM  1352   CG2 <NA>   VAL     B    77   <NA> 13.593  12.322  -7.843 1
    1352   ATOM  1353     N <NA>   GLY     B    78   <NA> 10.778  12.616  -4.125 1
    1353   ATOM  1354    CA <NA>   GLY     B    78   <NA> 10.004  12.229  -2.965 1
    1354   ATOM  1355     C <NA>   GLY     B    78   <NA>  8.832  13.128  -2.672 1
    1355   ATOM  1356     O <NA>   GLY     B    78   <NA>  8.614  14.117  -3.393 1
    1356   ATOM  1357     N <NA>   PRO     B    79   <NA>  8.032  12.814  -1.646 1
    1357   ATOM  1358    CA <NA>   PRO     B    79   <NA>  6.887  13.664  -1.350 1
    1358   ATOM  1359     C <NA>   PRO     B    79   <NA>  7.292  14.915  -0.550 1
    1359   ATOM  1360     O <NA>   PRO     B    79   <NA>  7.007  15.036   0.638 1
    1360   ATOM  1361    CB <NA>   PRO     B    79   <NA>  5.951  12.742  -0.594 1
    1361   ATOM  1362    CG <NA>   PRO     B    79   <NA>  6.838  11.781   0.040 1
    1362   ATOM  1363    CD <NA>   PRO     B    79   <NA>  8.096  11.665  -0.739 1
    1363   ATOM  1364     N <NA>   THR     B    80   <NA>  7.997  15.816  -1.220 1
    1364   ATOM  1365    CA <NA>   THR     B    80   <NA>  8.324  17.137  -0.702 1
    1365   ATOM  1366     C <NA>   THR     B    80   <NA>  7.227  18.114  -1.090 1
    1366   ATOM  1367     O <NA>   THR     B    80   <NA>  6.528  17.896  -2.080 1
    1367   ATOM  1368    CB <NA>   THR     B    80   <NA>  9.677  17.594  -1.299 1
    1368   ATOM  1369   OG1 <NA>   THR     B    80   <NA>  9.952  18.951  -0.924 1
    1369   ATOM  1370   CG2 <NA>   THR     B    80   <NA>  9.688  17.449  -2.825 1
    1370   ATOM  1371     N <NA>   PRO     B    81   <NA>  6.922  19.074  -0.214 1
    1371   ATOM  1372    CA <NA>   PRO     B    81   <NA>  5.896  20.066  -0.556 1
    1372   ATOM  1373     C <NA>   PRO     B    81   <NA>  6.244  20.969  -1.727 1
    1373   ATOM  1374     O <NA>   PRO     B    81   <NA>  5.343  21.294  -2.509 1
    1374   ATOM  1375    CB <NA>   PRO     B    81   <NA>  5.694  20.874   0.729 1
    1375   ATOM  1376    CG <NA>   PRO     B    81   <NA>  6.274  19.954   1.831 1
    1376   ATOM  1377    CD <NA>   PRO     B    81   <NA>  7.387  19.243   1.190 1
    1377   ATOM  1378     N <NA>   VAL     B    82   <NA>  7.520  21.355  -1.868 1
    1378   ATOM  1379    CA <NA>   VAL     B    82   <NA>  7.990  22.207  -2.983 1
    1379   ATOM  1380     C <NA>   VAL     B    82   <NA>  9.255  21.636  -3.581 1
    1380   ATOM  1381     O <NA>   VAL     B    82   <NA>  9.973  20.911  -2.929 1
    1381   ATOM  1382    CB <NA>   VAL     B    82   <NA>  8.375  23.616  -2.484 1
    1382   ATOM  1383   CG1 <NA>   VAL     B    82   <NA>  7.122  24.513  -2.361 1
    1383   ATOM  1384   CG2 <NA>   VAL     B    82   <NA>  9.101  23.485  -1.163 1
    1384   ATOM  1385     N <NA>   ASN     B    83   <NA>  9.588  21.960  -4.812 1
    1385   ATOM  1386    CA <NA>   ASN     B    83   <NA> 10.914  21.521  -5.319 1
    1386   ATOM  1387     C <NA>   ASN     B    83   <NA> 11.922  22.373  -4.576 1
    1387   ATOM  1388     O <NA>   ASN     B    83   <NA> 11.716  23.574  -4.429 1
    1388   ATOM  1389    CB <NA>   ASN     B    83   <NA> 11.068  21.763  -6.823 1
    1389   ATOM  1390    CG <NA>   ASN     B    83   <NA> 10.096  20.950  -7.647 1
    1390   ATOM  1391   OD1 <NA>   ASN     B    83   <NA> 10.013  19.733  -7.494 1
    1391   ATOM  1392   ND2 <NA>   ASN     B    83   <NA>  9.305  21.627  -8.498 1
    1392   ATOM  1393     N <NA>   ILE     B    84   <NA> 12.983  21.755  -4.066 1
    1393   ATOM  1394    CA <NA>   ILE     B    84   <NA> 13.979  22.449  -3.233 1
    1394   ATOM  1395     C <NA>   ILE     B    84   <NA> 15.345  22.255  -3.870 1
    1395   ATOM  1396     O <NA>   ILE     B    84   <NA> 15.779  21.111  -4.072 1
    1396   ATOM  1397    CB <NA>   ILE     B    84   <NA> 14.041  21.833  -1.837 1
    1397   ATOM  1398   CG1 <NA>   ILE     B    84   <NA> 12.754  22.151  -1.088 1
    1398   ATOM  1399   CG2 <NA>   ILE     B    84   <NA> 15.323  22.264  -1.122 1
    1399   ATOM  1400   CD1 <NA>   ILE     B    84   <NA> 12.438  21.217   0.060 1
    1400   ATOM  1401     N <NA>   ILE     B    85   <NA> 16.044  23.346  -4.167 1
    1401   ATOM  1402    CA <NA>   ILE     B    85   <NA> 17.441  23.243  -4.575 1
    1402   ATOM  1403     C <NA>   ILE     B    85   <NA> 18.305  23.504  -3.345 1
    1403   ATOM  1404     O <NA>   ILE     B    85   <NA> 18.262  24.613  -2.802 1
    1404   ATOM  1405    CB <NA>   ILE     B    85   <NA> 17.805  24.332  -5.644 1
    1405   ATOM  1406   CG1 <NA>   ILE     B    85   <NA> 16.809  24.299  -6.827 1
    1406   ATOM  1407   CG2 <NA>   ILE     B    85   <NA> 19.271  24.191  -6.067 1
    1407   ATOM  1408   CD1 <NA>   ILE     B    85   <NA> 16.672  22.956  -7.581 1
    1408   ATOM  1409     N <NA>   GLY     B    86   <NA> 19.121  22.530  -2.936 1
    1409   ATOM  1410    CA <NA>   GLY     B    86   <NA> 19.857  22.642  -1.690 1
    1410   ATOM  1411     C <NA>   GLY     B    86   <NA> 21.324  22.800  -1.940 1
    1411   ATOM  1412     O <NA>   GLY     B    86   <NA> 21.750  22.958  -3.056 1
    1412   ATOM  1413     N <NA>   ARG     B    87   <NA> 22.117  22.755  -0.887 1
    1413   ATOM  1414    CA <NA>   ARG     B    87   <NA> 23.533  23.126  -0.985 1
    1414   ATOM  1415     C <NA>   ARG     B    87   <NA> 24.413  22.356  -1.973 1
    1415   ATOM  1416     O <NA>   ARG     B    87   <NA> 25.398  22.903  -2.461 1
    1416   ATOM  1417    CB <NA>   ARG     B    87   <NA> 24.171  23.147   0.403 1
    1417   ATOM  1418    CG <NA>   ARG     B    87   <NA> 23.646  24.256   1.283 1
    1418   ATOM  1419    CD <NA>   ARG     B    87   <NA> 24.429  24.306   2.576 1
    1419   ATOM  1420    NE <NA>   ARG     B    87   <NA> 24.362  23.068   3.371 1
    1420   ATOM  1421    CZ <NA>   ARG     B    87   <NA> 25.357  22.185   3.495 1
    1421   ATOM  1422   NH1 <NA>   ARG     B    87   <NA> 26.467  22.303   2.795 1
    1422   ATOM  1423   NH2 <NA>   ARG     B    87   <NA> 25.255  21.190   4.368 1
    1423   ATOM  1424     N <NA>   ASN     B    88   <NA> 24.074  21.101  -2.275 1
    1424   ATOM  1425    CA <NA>   ASN     B    88   <NA> 24.950  20.312  -3.132 1
    1425   ATOM  1426     C <NA>   ASN     B    88   <NA> 24.980  20.893  -4.527 1
    1426   ATOM  1427     O <NA>   ASN     B    88   <NA> 26.015  20.853  -5.202 1
    1427   ATOM  1428    CB <NA>   ASN     B    88   <NA> 24.512  18.849  -3.197 1
    1428   ATOM  1429    CG <NA>   ASN     B    88   <NA> 23.126  18.649  -3.777 1
    1429   ATOM  1430   OD1 <NA>   ASN     B    88   <NA> 22.167  19.280  -3.351 1
    1430   ATOM  1431   ND2 <NA>   ASN     B    88   <NA> 23.020  17.775  -4.767 1
    1431   ATOM  1432     N <NA>   LEU     B    89   <NA> 23.863  21.490  -4.949 1
    1432   ATOM  1433    CA <NA>   LEU     B    89   <NA> 23.811  22.125  -6.273 1
    1433   ATOM  1434     C <NA>   LEU     B    89   <NA> 24.018  23.645  -6.231 1
    1434   ATOM  1435     O <NA>   LEU     B    89   <NA> 24.321  24.242  -7.236 1
    1435   ATOM  1436    CB <NA>   LEU     B    89   <NA> 22.457  21.815  -6.962 1
    1436   ATOM  1437    CG <NA>   LEU     B    89   <NA> 22.219  20.372  -7.436 1
    1437   ATOM  1438   CD1 <NA>   LEU     B    89   <NA> 20.937  20.300  -8.243 1
    1438   ATOM  1439   CD2 <NA>   LEU     B    89   <NA> 23.408  19.901  -8.273 1
    1439   ATOM  1440     N <NA>   LEU     B    90   <NA> 23.819  24.255  -5.075 1
    1440   ATOM  1441    CA <NA>   LEU     B    90   <NA> 24.020  25.701  -4.954 1
    1441   ATOM  1442     C <NA>   LEU     B    90   <NA> 25.511  26.005  -5.072 1
    1442   ATOM  1443     O <NA>   LEU     B    90   <NA> 25.907  26.959  -5.732 1
    1443   ATOM  1444    CB <NA>   LEU     B    90   <NA> 23.430  26.234  -3.624 1
    1444   ATOM  1445    CG <NA>   LEU     B    90   <NA> 21.900  26.309  -3.475 1
    1445   ATOM  1446   CD1 <NA>   LEU     B    90   <NA> 21.487  26.779  -2.081 1
    1446   ATOM  1447   CD2 <NA>   LEU     B    90   <NA> 21.358  27.241  -4.506 1
    1447   ATOM  1448     N <NA>   THR     B    91   <NA> 26.336  25.143  -4.491 1
    1448   ATOM  1449    CA <NA>   THR     B    91   <NA> 27.785  25.304  -4.546 1
    1449   ATOM  1450     C <NA>   THR     B    91   <NA> 28.270  25.184  -5.969 1
    1450   ATOM  1451     O <NA>   THR     B    91   <NA> 29.168  25.903  -6.378 1
    1451   ATOM  1452    CB <NA>   THR     B    91   <NA> 28.501  24.261  -3.669 1
    1452   ATOM  1453   OG1 <NA>   THR     B    91   <NA> 27.898  22.990  -3.894 1
    1453   ATOM  1454   CG2 <NA>   THR     B    91   <NA> 28.366  24.579  -2.208 1
    1454   ATOM  1455     N <NA>   GLN     B    92   <NA> 27.619  24.323  -6.741 1
    1455   ATOM  1456    CA <NA>   GLN     B    92   <NA> 28.009  24.110  -8.150 1
    1456   ATOM  1457     C <NA>   GLN     B    92   <NA> 27.823  25.361  -8.999 1
    1457   ATOM  1458     O <NA>   GLN     B    92   <NA> 28.719  25.701  -9.794 1
    1458   ATOM  1459    CB <NA>   GLN     B    92   <NA> 27.226  22.956  -8.793 1
    1459   ATOM  1460    CG <NA>   GLN     B    92   <NA> 27.720  21.588  -8.406 1
    1460   ATOM  1461    CD <NA>   GLN     B    92   <NA> 27.313  20.496  -9.421 1
    1461   ATOM  1462   OE1 <NA>   GLN     B    92   <NA> 27.138  20.746 -10.620 1
    1462   ATOM  1463   NE2 <NA>   GLN     B    92   <NA> 27.230  19.281  -8.941 1
    1463   ATOM  1464     N <NA>   ILE     B    93   <NA> 26.683  26.043  -8.820 1
    1464   ATOM  1465    CA <NA>   ILE     B    93   <NA> 26.362  27.233  -9.606 1
    1465   ATOM  1466     C <NA>   ILE     B    93   <NA> 26.904  28.524  -8.963 1
    1466   ATOM  1467     O <NA>   ILE     B    93   <NA> 26.574  29.627  -9.385 1
    1467   ATOM  1468    CB <NA>   ILE     B    93   <NA> 24.838  27.349  -9.857 1
    1468   ATOM  1469   CG1 <NA>   ILE     B    93   <NA> 24.103  27.648  -8.559 1
    1469   ATOM  1470   CG2 <NA>   ILE     B    93   <NA> 24.300  26.092 -10.460 1
    1470   ATOM  1471   CD1 <NA>   ILE     B    93   <NA> 22.672  28.017  -8.785 1
    1471   ATOM  1472     N <NA>   GLY     B    94   <NA> 27.741  28.369  -7.949 1
    1472   ATOM  1473    CA <NA>   GLY     B    94   <NA> 28.481  29.498  -7.421 1
    1473   ATOM  1474     C <NA>   GLY     B    94   <NA> 27.749  30.439  -6.488 1
    1474   ATOM  1475     O <NA>   GLY     B    94   <NA> 28.177  31.555  -6.261 1
    1475   ATOM  1476     N <NA>   CYS     B    95   <NA> 26.808  29.899  -5.748 1
    1476   ATOM  1477    CA <NA>   CYS     B    95   <NA> 25.895  30.718  -5.004 1
    1477   ATOM  1478     C <NA>   CYS     B    95   <NA> 26.408  30.993  -3.598 1
    1478   ATOM  1479     O <NA>   CYS     B    95   <NA> 26.769  30.065  -2.870 1
    1479   ATOM  1480    CB <NA>   CYS     B    95   <NA> 24.578  29.989  -4.978 1
    1480   ATOM  1481    SG <NA>   CYS     B    95   <NA> 23.221  30.929  -4.410 1
    1481   ATOM  1482     N <NA>   THR     B    96   <NA> 26.473  32.277  -3.244 1
    1482   ATOM  1483    CA <NA>   THR     B    96   <NA> 26.794  32.734  -1.882 1
    1483   ATOM  1484     C <NA>   THR     B    96   <NA> 25.672  33.544  -1.205 1
    1484   ATOM  1485     O <NA>   THR     B    96   <NA> 24.760  34.079  -1.852 1
    1485   ATOM  1486    CB <NA>   THR     B    96   <NA> 28.051  33.660  -1.858 1
    1486   ATOM  1487   OG1 <NA>   THR     B    96   <NA> 27.888  34.689  -2.857 1
    1487   ATOM  1488   CG2 <NA>   THR     B    96   <NA> 29.316  32.870  -2.141 1
    1488   ATOM  1489     N <NA>   LEU     B    97   <NA> 25.759  33.617   0.119 1
    1489   ATOM  1490    CA <NA>   LEU     B    97   <NA> 24.902  34.468   0.963 1
    1490   ATOM  1491     C <NA>   LEU     B    97   <NA> 25.714  35.689   1.398 1
    1491   ATOM  1492     O <NA>   LEU     B    97   <NA> 26.854  35.558   1.870 1
    1492   ATOM  1493    CB <NA>   LEU     B    97   <NA> 24.489  33.718   2.236 1
    1493   ATOM  1494    CG <NA>   LEU     B    97   <NA> 23.211  32.901   2.344 1
    1494   ATOM  1495   CD1 <NA>   LEU     B    97   <NA> 23.114  32.358   3.719 1
    1495   ATOM  1496   CD2 <NA>   LEU     B    97   <NA> 22.037  33.773   2.076 1
    1496   ATOM  1497     N <NA>   ASN     B    98   <NA> 25.121  36.868   1.264 1
    1497   ATOM  1498    CA <NA>   ASN     B    98   <NA> 25.870  38.101   1.449 1
    1498   ATOM  1499     C <NA>   ASN     B    98   <NA> 25.102  39.038   2.370 1
    1499   ATOM  1500     O <NA>   ASN     B    98   <NA> 23.889  39.124   2.261 1
    1500   ATOM  1501    CB <NA>   ASN     B    98   <NA> 26.140  38.756   0.086 1
    1501   ATOM  1502    CG <NA>   ASN     B    98   <NA> 27.048  37.921  -0.814 1
    1502   ATOM  1503   OD1 <NA>   ASN     B    98   <NA> 28.268  37.895  -0.630 1
    1503   ATOM  1504   ND2 <NA>   ASN     B    98   <NA> 26.455  37.192  -1.754 1
    1504   ATOM  1505     N <NA>   PHE     B    99   <NA> 25.809  39.706   3.283 1
    1505   ATOM  1506    CA <NA>   PHE     B    99   <NA> 25.267  40.855   4.034 1
    1506   ATOM  1507     C <NA>   PHE     B    99   <NA> 26.351  41.742   4.659 1
    1507   ATOM  1508     O <NA>   PHE     B    99   <NA> 27.448  41.208   5.013 1
    1508   ATOM  1509    CB <NA>   PHE     B    99   <NA> 24.284  40.418   5.127 1
    1509   ATOM  1510    CG <NA>   PHE     B    99   <NA> 24.859  39.479   6.130 1
    1510   ATOM  1511   CD1 <NA>   PHE     B    99   <NA> 25.061  38.150   5.808 1
    1511   ATOM  1512   CD2 <NA>   PHE     B    99   <NA> 25.071  39.890   7.436 1
    1512   ATOM  1513   CE1 <NA>   PHE     B    99   <NA> 25.450  37.240   6.756 1
    1513   ATOM  1514   CE2 <NA>   PHE     B    99   <NA> 25.473  38.988   8.409 1
    1514   ATOM  1515    CZ <NA>   PHE     B    99   <NA> 25.658  37.663   8.073 1
    1515 HETATM  1517    N1 <NA>   MK1     B   902   <NA>  9.280  23.763   3.004 1
    1516 HETATM  1518    C1 <NA>   MK1     B   902   <NA>  9.498  23.983   4.459 1
    1517 HETATM  1519    C2 <NA>   MK1     B   902   <NA> 10.591  24.905   4.962 1
    1518 HETATM  1520    C3 <NA>   MK1     B   902   <NA> 10.591  24.864   6.466 1
    1519 HETATM  1521    O1 <NA>   MK1     B   902   <NA> 10.937  23.849   7.057 1
    1520 HETATM  1522    N2 <NA>   MK1     B   902   <NA> 10.193  25.953   7.094 1
    1521 HETATM  1523    C4 <NA>   MK1     B   902   <NA> 10.145  26.250   8.490 1
    1522 HETATM  1524    C5 <NA>   MK1     B   902   <NA>  9.379  27.577   8.641 1
    1523 HETATM  1525    C6 <NA>   MK1     B   902   <NA> 11.398  26.347   9.074 1
    1524 HETATM  1526    C7 <NA>   MK1     B   902   <NA>  9.364  25.283   9.268 1
    1525 HETATM  1527    N3 <NA>   MK1     B   902   <NA> 11.819  24.282   4.355 1
    1526 HETATM  1528    C8 <NA>   MK1     B   902   <NA> 11.753  23.776   2.961 1
    1527 HETATM  1529    C9 <NA>   MK1     B   902   <NA> 10.440  23.182   2.493 1
    1528 HETATM  1530   C10 <NA>   MK1     B   902   <NA> 13.083  24.963   4.552 1
    1529 HETATM  1531   C11 <NA>   MK1     B   902   <NA> 14.203  24.064   5.078 1
    1530 HETATM  1532    O2 <NA>   MK1     B   902   <NA> 15.242  24.884   4.634 1
    1531 HETATM  1533   C12 <NA>   MK1     B   902   <NA> 14.440  23.761   6.569 1
    1532 HETATM  1534   C13 <NA>   MK1     B   902   <NA> 15.573  22.821   7.005 1
    1533 HETATM  1535   C14 <NA>   MK1     B   902   <NA> 15.644  22.664   8.534 1
    1534 HETATM  1536   C15 <NA>   MK1     B   902   <NA> 16.733  21.750   8.961 1
    1535 HETATM  1537   C16 <NA>   MK1     B   902   <NA> 18.058  21.916   8.553 1
    1536 HETATM  1538   C17 <NA>   MK1     B   902   <NA> 19.037  21.016   8.947 1
    1537 HETATM  1539   C18 <NA>   MK1     B   902   <NA> 18.673  19.939   9.758 1
    1538 HETATM  1540   C19 <NA>   MK1     B   902   <NA> 17.347  19.773  10.176 1
    1539 HETATM  1541   C20 <NA>   MK1     B   902   <NA> 16.374  20.687   9.772 1
    1540 HETATM  1542   C21 <NA>   MK1     B   902   <NA> 15.447  21.440   6.373 1
    1541 HETATM  1543    O3 <NA>   MK1     B   902   <NA> 14.367  20.831   6.397 1
    1542 HETATM  1544    N4 <NA>   MK1     B   902   <NA> 16.583  20.913   5.924 1
    1543 HETATM  1545   C22 <NA>   MK1     B   902   <NA> 16.692  19.500   5.604 1
    1544 HETATM  1546   C23 <NA>   MK1     B   902   <NA> 18.067  18.945   5.936 1
    1545 HETATM  1547    O4 <NA>   MK1     B   902   <NA> 19.061  19.938   5.729 1
    1546 HETATM  1548   C24 <NA>   MK1     B   902   <NA> 18.226  17.726   5.057 1
    1547 HETATM  1549   C25 <NA>   MK1     B   902   <NA> 17.476  17.904   3.760 1
    1548 HETATM  1550   C26 <NA>   MK1     B   902   <NA> 17.500  17.363   2.496 1
    1549 HETATM  1551   C27 <NA>   MK1     B   902   <NA> 16.613  17.872   1.541 1
    1550 HETATM  1552   C28 <NA>   MK1     B   902   <NA> 15.722  18.906   1.865 1
    1551 HETATM  1553   C29 <NA>   MK1     B   902   <NA> 15.683  19.479   3.129 1
    1552 HETATM  1554   C30 <NA>   MK1     B   902   <NA> 16.504  19.061   4.128 1
    1553 HETATM  1555   C31 <NA>   MK1     B   902   <NA>  8.033  23.100   2.604 1
    1554 HETATM  1556   C32 <NA>   MK1     B   902   <NA>  6.666  23.739   2.876 1
    1555 HETATM  1557   C33 <NA>   MK1     B   902   <NA>  6.158  24.808   2.124 1
    1556 HETATM  1558    N5 <NA>   MK1     B   902   <NA>  4.911  25.430   2.300 1
    1557 HETATM  1559   C34 <NA>   MK1     B   902   <NA>  4.207  24.839   3.348 1
    1558 HETATM  1560   C35 <NA>   MK1     B   902   <NA>  4.654  23.774   4.136 1
    1559 HETATM  1561   C36 <NA>   MK1     B   902   <NA>  5.905  23.211   3.897 1
    1560 HETATM  1562     O <NA>   HOH     A   305   <NA> 20.857  43.192  21.450 1
    1561 HETATM  1563     O <NA>   HOH     A   307   <NA> 14.076  19.789  19.440 1
    1562 HETATM  1564     O <NA>   HOH     A   309   <NA> 28.075  21.177   7.222 1
    1563 HETATM  1565     O <NA>   HOH     A   314   <NA> 16.759  40.274   1.287 1
    1564 HETATM  1566     O <NA>   HOH     A   315   <NA> 13.997  22.233  21.468 1
    1565 HETATM  1567     O <NA>   HOH     A   324   <NA> 11.282  30.738   1.625 1
    1566 HETATM  1568     O <NA>   HOH     A   325   <NA> 16.774  42.740   2.296 1
    1567 HETATM  1569     O <NA>   HOH     A   327   <NA> 14.623  29.552  28.235 1
    1568 HETATM  1570     O <NA>   HOH     A   328   <NA>  1.651  36.463  19.459 1
    1569 HETATM  1571     O <NA>   HOH     A   329   <NA> 14.435  44.966  11.913 1
    1570 HETATM  1572     O <NA>   HOH     A   330   <NA> 19.877  40.160  21.917 1
    1571 HETATM  1573     O <NA>   HOH     A   331   <NA> 17.126  45.934  10.801 1
    1572 HETATM  1574     O <NA>   HOH     A   332   <NA>  8.840  28.026   4.860 1
    1573 HETATM  1575     O <NA>   HOH     A   335   <NA> 10.341  18.480  14.477 1
    1574 HETATM  1576     O <NA>   HOH     A   341   <NA> 19.233  16.711   9.027 1
    1575 HETATM  1577     O <NA>   HOH     A   342   <NA> 23.799  21.928  12.391 1
    1576 HETATM  1578     O <NA>   HOH     A   344   <NA>  9.953  37.934   4.548 1
    1577 HETATM  1579     O <NA>   HOH     A   345   <NA>  8.478  35.995   5.789 1
    1578 HETATM  1580     O <NA>   HOH     A   357   <NA>  3.960  19.389  17.384 1
    1579 HETATM  1581     O <NA>   HOH     A   373   <NA> 27.561  43.155  19.015 1
    1580 HETATM  1582     O <NA>   HOH     A   384   <NA>  1.245  19.292  18.124 1
    1581 HETATM  1583     O <NA>   HOH     A   386   <NA> 31.402  27.051   3.335 1
    1582 HETATM  1584     O <NA>   HOH     A   389   <NA> 32.446  31.200   4.417 1
    1583 HETATM  1585     O <NA>   HOH     A   391   <NA> 25.480  38.468  17.938 1
    1584 HETATM  1586     O <NA>   HOH     A   394   <NA> 23.940  41.721   0.346 1
    1585 HETATM  1587     O <NA>   HOH     A   401   <NA>  5.912  15.727   3.369 1
    1586 HETATM  1588     O <NA>   HOH     A   406   <NA>  9.272  33.891  12.681 1
    1587 HETATM  1589     O <NA>   HOH     A   408   <NA> 21.185  25.233  16.048 1
    1588 HETATM  1590     O <NA>   HOH     A   416   <NA> 18.474  26.012  21.664 1
    1589 HETATM  1591     O <NA>   HOH     A   420   <NA>  9.469  16.910  17.371 1
    1590 HETATM  1592     O <NA>   HOH     A   422   <NA> 13.074  17.786  16.615 1
    1591 HETATM  1593     O <NA>   HOH     A   439   <NA> 28.821  29.338   7.342 1
    1592 HETATM  1594     O <NA>   HOH     A   457   <NA> 23.284  23.107  15.132 1
    1593 HETATM  1595     O <NA>   HOH     A   468   <NA>  3.114  26.260   6.773 1
    1594 HETATM  1596     O <NA>   HOH     A   501   <NA>  6.382  26.424   5.893 1
    1595 HETATM  1597     O <NA>   HOH     A   503   <NA> 35.293  43.006   5.212 1
    1596 HETATM  1598     O <NA>   HOH     A   510   <NA> 21.891  49.715   7.192 1
    1597 HETATM  1599     O <NA>   HOH     A   524   <NA> 34.085  32.735   2.849 1
    1598 HETATM  1600     O <NA>   HOH     A   529   <NA> 31.491  41.524   6.678 1
    1599 HETATM  1601     O <NA>   HOH     A   553   <NA>  5.943  34.223   6.748 1
    1600 HETATM  1602     O <NA>   HOH     A   561   <NA>  0.934  40.259  19.405 1
    1601 HETATM  1603     O <NA>   HOH     A   567   <NA> 29.539  25.486  13.281 1
    1602 HETATM  1604     O <NA>   HOH     A   572   <NA> 24.552  17.352  10.295 1
    1603 HETATM  1605     O <NA>   HOH     A   575   <NA> 23.112  15.510   8.776 1
    1604 HETATM  1606     O <NA>   HOH     B   301   <NA> 20.445   8.036 -12.631 1
    1605 HETATM  1607     O <NA>   HOH     B   303   <NA> 20.044  14.822   4.638 1
    1606 HETATM  1608     O <NA>   HOH     B   304   <NA> 21.538   6.875 -10.099 1
    1607 HETATM  1609     O <NA>   HOH     B   306   <NA> 22.449  23.958   5.252 1
    1608 HETATM  1610     O <NA>   HOH     B   308   <NA> 11.720  21.289   7.190 1
    1609 HETATM  1611     O <NA>   HOH     B   312   <NA> 14.097   5.111 -11.638 1
    1610 HETATM  1612     O <NA>   HOH     B   313   <NA> 20.998  21.834   6.561 1
    1611 HETATM  1613     O <NA>   HOH     B   316   <NA> 22.659  14.583  -2.196 1
    1612 HETATM  1614     O <NA>   HOH     B   317   <NA> 28.724  15.629 -11.660 1
    1613 HETATM  1615     O <NA>   HOH     B   318   <NA> 16.539  45.207   0.079 1
    1614 HETATM  1616     O <NA>   HOH     B   319   <NA> 23.678  14.931   2.680 1
    1615 HETATM  1617     O <NA>   HOH     B   321   <NA> 20.718  15.976  -3.657 1
    1616 HETATM  1618     O <NA>   HOH     B   323   <NA> 31.249  26.796  -9.595 1
    1617 HETATM  1619     O <NA>   HOH     B   326   <NA> 28.813  28.445  -2.106 1
    1618 HETATM  1620     O <NA>   HOH     B   333   <NA> 12.251  39.551  -2.672 1
    1619 HETATM  1621     O <NA>   HOH     B   334   <NA> 25.465  12.592  -8.670 1
    1620 HETATM  1622     O <NA>   HOH     B   338   <NA> 12.998  36.205  -3.972 1
    1621 HETATM  1623     O <NA>   HOH     B   339   <NA> 17.541  17.060 -17.194 1
    1622 HETATM  1624     O <NA>   HOH     B   340   <NA>  5.321  14.325  -4.866 1
    1623 HETATM  1625     O <NA>   HOH     B   346   <NA>  9.314  17.330  -9.801 1
    1624 HETATM  1626     O <NA>   HOH     B   347   <NA>  7.435  26.652 -14.854 1
    1625 HETATM  1627     O <NA>   HOH     B   348   <NA>  4.405  16.704  -3.635 1
    1626 HETATM  1628     O <NA>   HOH     B   349   <NA> 19.414   7.026   4.428 1
    1627 HETATM  1629     O <NA>   HOH     B   350   <NA>  6.718  34.538  -1.322 1
    1628 HETATM  1630     O <NA>   HOH     B   354   <NA> 15.041  31.743 -13.235 1
    1629 HETATM  1631     O <NA>   HOH     B   355   <NA> 27.404  32.078 -10.860 1
    1630 HETATM  1632     O <NA>   HOH     B   356   <NA> 27.673  18.789  -6.155 1
    1631 HETATM  1633     O <NA>   HOH     B   358   <NA> 21.289  -1.161  -5.102 1
    1632 HETATM  1634     O <NA>   HOH     B   359   <NA>  6.973  36.523   1.489 1
    1633 HETATM  1635     O <NA>   HOH     B   360   <NA> 27.602  21.234  -0.635 1
    1634 HETATM  1636     O <NA>   HOH     B   362   <NA>  3.902   9.376  -0.027 1
    1635 HETATM  1637     O <NA>   HOH     B   364   <NA> 28.498  36.632  -7.529 1
    1636 HETATM  1638     O <NA>   HOH     B   366   <NA> 18.572  40.567 -10.042 1
    1637 HETATM  1639     O <NA>   HOH     B   367   <NA> 25.658  18.970   0.428 1
    1638 HETATM  1640     O <NA>   HOH     B   369   <NA> 20.843   1.263  -7.014 1
    1639 HETATM  1641     O <NA>   HOH     B   370   <NA> 13.975  15.741  12.070 1
    1640 HETATM  1642     O <NA>   HOH     B   374   <NA>  7.661  23.876  -6.324 1
    1641 HETATM  1643     O <NA>   HOH     B   375   <NA> 10.125   5.706  -1.458 1
    1642 HETATM  1644     O <NA>   HOH     B   376   <NA> 18.450  20.497 -18.728 1
    1643 HETATM  1645     O <NA>   HOH     B   377   <NA> 29.267  20.487  -3.497 1
    1644 HETATM  1646     O <NA>   HOH     B   379   <NA>  6.685  26.541  -5.608 1
    1645 HETATM  1647     O <NA>   HOH     B   381   <NA> 25.810  26.789 -19.106 1
    1646 HETATM  1648     O <NA>   HOH     B   383   <NA> 21.144  -4.428 -11.331 1
    1647 HETATM  1649     O <NA>   HOH     B   387   <NA> 16.904  27.594 -15.938 1
    1648 HETATM  1650     O <NA>   HOH     B   388   <NA> 23.926  45.612  -4.998 1
    1649 HETATM  1651     O <NA>   HOH     B   390   <NA> 25.300  17.493   3.076 1
    1650 HETATM  1652     O <NA>   HOH     B   392   <NA>  6.618  28.079  -3.427 1
    1651 HETATM  1653     O <NA>   HOH     B   393   <NA> 19.795  13.651 -16.606 1
    1652 HETATM  1654     O <NA>   HOH     B   395   <NA>  7.202   9.982  -4.103 1
    1653 HETATM  1655     O <NA>   HOH     B   400   <NA>  8.474  34.203  -4.893 1
    1654 HETATM  1656     O <NA>   HOH     B   405   <NA> 16.659  15.866  11.446 1
    1655 HETATM  1657     O <NA>   HOH     B   410   <NA> 26.400  10.057  -3.287 1
    1656 HETATM  1658     O <NA>   HOH     B   414   <NA>  9.503   3.489  -4.419 1
    1657 HETATM  1659     O <NA>   HOH     B   419   <NA> 15.438  12.973 -18.484 1
    1658 HETATM  1660     O <NA>   HOH     B   425   <NA> 11.428  19.956 -24.551 1
    1659 HETATM  1661     O <NA>   HOH     B   430   <NA> 18.725  43.171  -5.575 1
    1660 HETATM  1662     O <NA>   HOH     B   436   <NA> 32.141  29.620  -8.580 1
    1661 HETATM  1663     O <NA>   HOH     B   443   <NA>  8.811  13.667 -20.256 1
    1662 HETATM  1664     O <NA>   HOH     B   444   <NA>  4.071  26.169  -0.230 1
    1663 HETATM  1665     O <NA>   HOH     B   461   <NA> 11.425  44.636  -3.033 1
    1664 HETATM  1666     O <NA>   HOH     B   469   <NA>  6.902  23.686 -10.066 1
    1665 HETATM  1667     O <NA>   HOH     B   471   <NA>  5.749  25.785 -19.792 1
    1666 HETATM  1668     O <NA>   HOH     B   500   <NA> 25.592  16.404  -5.805 1
    1667 HETATM  1669     O <NA>   HOH     B   502   <NA>  4.040  15.516  -7.200 1
    1668 HETATM  1670     O <NA>   HOH     B   505   <NA> 28.640  34.232  -5.637 1
    1669 HETATM  1671     O <NA>   HOH     B   506   <NA>  8.979  11.173  11.112 1
    1670 HETATM  1672     O <NA>   HOH     B   509   <NA> 19.882   3.986 -18.136 1
    1671 HETATM  1673     O <NA>   HOH     B   514   <NA> 27.409  15.355   2.200 1
    1672 HETATM  1674     O <NA>   HOH     B   515   <NA> 17.222  39.766 -23.774 1
    1673 HETATM  1675     O <NA>   HOH     B   517   <NA> 28.742  24.158 -16.641 1
    1674 HETATM  1676     O <NA>   HOH     B   525   <NA> 22.694  -2.192 -12.589 1
    1675 HETATM  1677     O <NA>   HOH     B   526   <NA> 17.901  43.157 -14.082 1
    1676 HETATM  1678     O <NA>   HOH     B   531   <NA> 18.192   8.914  11.344 1
    1677 HETATM  1679     O <NA>   HOH     B   532   <NA> 19.507  45.215   1.709 1
    1678 HETATM  1680     O <NA>   HOH     B   548   <NA>  1.442  14.700  -6.128 1
    1679 HETATM  1681     O <NA>   HOH     B   549   <NA> 19.908   8.718 -19.215 1
    1680 HETATM  1682     O <NA>   HOH     B   556   <NA> 21.499  44.884  -1.280 1
    1681 HETATM  1683     O <NA>   HOH     B   564   <NA> 10.031   8.593 -22.052 1
    1682 HETATM  1684     O <NA>   HOH     B   568   <NA>  2.817  28.133   2.191 1
    1683 HETATM  1685     O <NA>   HOH     B   591   <NA> 15.835  40.105  -5.971 1
    1684 HETATM  1686     O <NA>   HOH     B   595   <NA>  4.515  36.451  -4.499 1
    1685 HETATM  1687     O <NA>   HOH     B   613   <NA> 24.127 -10.994  -0.982 1
    1686 HETATM  1688     O <NA>   HOH     B   617   <NA> 30.112  17.912  -4.791 1
             b segid elesy charge
    1    38.10  <NA>     N   <NA>
    2    40.62  <NA>     C   <NA>
    3    42.64  <NA>     C   <NA>
    4    43.40  <NA>     O   <NA>
    5    37.87  <NA>     C   <NA>
    6    38.40  <NA>     C   <NA>
    7    38.74  <NA>     C   <NA>
    8    41.76  <NA>     N   <NA>
    9    41.30  <NA>     C   <NA>
    10   41.38  <NA>     C   <NA>
    11   43.09  <NA>     O   <NA>
    12   40.81  <NA>     C   <NA>
    13   46.61  <NA>     C   <NA>
    14   50.36  <NA>     C   <NA>
    15   53.89  <NA>     O   <NA>
    16   51.46  <NA>     N   <NA>
    17   37.80  <NA>     N   <NA>
    18   34.13  <NA>     C   <NA>
    19   33.19  <NA>     C   <NA>
    20   32.74  <NA>     O   <NA>
    21   34.34  <NA>     C   <NA>
    22   33.95  <NA>     C   <NA>
    23   33.06  <NA>     C   <NA>
    24   32.50  <NA>     C   <NA>
    25   31.65  <NA>     N   <NA>
    26   30.14  <NA>     C   <NA>
    27   29.74  <NA>     C   <NA>
    28   27.88  <NA>     O   <NA>
    29   29.24  <NA>     C   <NA>
    30   27.60  <NA>     O   <NA>
    31   26.26  <NA>     C   <NA>
    32   29.27  <NA>     N   <NA>
    33   30.12  <NA>     C   <NA>
    34   32.98  <NA>     C   <NA>
    35   32.32  <NA>     O   <NA>
    36   26.21  <NA>     C   <NA>
    37   24.56  <NA>     C   <NA>
    38   20.74  <NA>     C   <NA>
    39   21.87  <NA>     C   <NA>
    40   31.52  <NA>     N   <NA>
    41   30.82  <NA>     C   <NA>
    42   31.90  <NA>     C   <NA>
    43   33.26  <NA>     O   <NA>
    44   28.66  <NA>     C   <NA>
    45   26.42  <NA>     C   <NA>
    46   28.49  <NA>     C   <NA>
    47   28.47  <NA>     C   <NA>
    48   29.88  <NA>     N   <NA>
    49   28.86  <NA>     C   <NA>
    50   25.06  <NA>     C   <NA>
    51   32.14  <NA>     C   <NA>
    52   27.41  <NA>     C   <NA>
    53   30.03  <NA>     C   <NA>
    54   33.01  <NA>     N   <NA>
    55   32.37  <NA>     C   <NA>
    56   29.95  <NA>     C   <NA>
    57   26.55  <NA>     O   <NA>
    58   38.12  <NA>     C   <NA>
    59   48.91  <NA>     C   <NA>
    60   59.75  <NA>     C   <NA>
    61   61.83  <NA>     O   <NA>
    62   59.99  <NA>     N   <NA>
    63   27.34  <NA>     N   <NA>
    64   29.87  <NA>     C   <NA>
    65   31.94  <NA>     C   <NA>
    66   33.83  <NA>     O   <NA>
    67   28.16  <NA>     C   <NA>
    68   27.47  <NA>     C   <NA>
    69   25.45  <NA>     C   <NA>
    70   23.06  <NA>     N   <NA>
    71   28.80  <NA>     C   <NA>
    72   30.95  <NA>     N   <NA>
    73   26.02  <NA>     N   <NA>
    74   30.21  <NA>     N   <NA>
    75   29.48  <NA>     C   <NA>
    76   29.45  <NA>     C   <NA>
    77   28.45  <NA>     O   <NA>
    78   27.88  <NA>     C   <NA>
    79   28.12  <NA>     C   <NA>
    80   31.64  <NA>     C   <NA>
    81   28.83  <NA>     N   <NA>
    82   31.57  <NA>     C   <NA>
    83   30.48  <NA>     C   <NA>
    84   31.00  <NA>     O   <NA>
    85   31.09  <NA>     C   <NA>
    86   35.91  <NA>     C   <NA>
    87   40.15  <NA>     C   <NA>
    88   40.51  <NA>     C   <NA>
    89   30.80  <NA>     N   <NA>
    90   30.14  <NA>     C   <NA>
    91   33.13  <NA>     C   <NA>
    92   34.48  <NA>     O   <NA>
    93   27.12  <NA>     C   <NA>
    94   28.48  <NA>     C   <NA>
    95   26.29  <NA>     C   <NA>
    96   32.34  <NA>     N   <NA>
    97   32.56  <NA>     C   <NA>
    98   33.07  <NA>     C   <NA>
    99   33.62  <NA>     O   <NA>
    100  33.81  <NA>     C   <NA>
    101  40.47  <NA>     O   <NA>
    102  34.22  <NA>     C   <NA>
    103  31.84  <NA>     N   <NA>
    104  32.26  <NA>     C   <NA>
    105  33.69  <NA>     C   <NA>
    106  30.43  <NA>     O   <NA>
    107  32.80  <NA>     C   <NA>
    108  31.81  <NA>     C   <NA>
    109  27.69  <NA>     C   <NA>
    110  32.46  <NA>     C   <NA>
    111  36.78  <NA>     N   <NA>
    112  38.20  <NA>     C   <NA>
    113  37.51  <NA>     C   <NA>
    114  33.78  <NA>     O   <NA>
    115  43.07  <NA>     C   <NA>
    116  50.67  <NA>     C   <NA>
    117  56.97  <NA>     C   <NA>
    118  62.89  <NA>     C   <NA>
    119  69.50  <NA>     N   <NA>
    120  40.26  <NA>     N   <NA>
    121  46.34  <NA>     C   <NA>
    122  49.77  <NA>     C   <NA>
    123  52.38  <NA>     O   <NA>
    124  45.04  <NA>     C   <NA>
    125  46.91  <NA>     C   <NA>
    126  47.78  <NA>     C   <NA>
    127  50.24  <NA>     C   <NA>
    128  53.37  <NA>     N   <NA>
    129  56.32  <NA>     C   <NA>
    130  56.91  <NA>     C   <NA>
    131  55.53  <NA>     O   <NA>
    132  57.90  <NA>     N   <NA>
    133  59.71  <NA>     C   <NA>
    134  60.57  <NA>     C   <NA>
    135  63.20  <NA>     O   <NA>
    136  59.37  <NA>     N   <NA>
    137  58.98  <NA>     C   <NA>
    138  56.87  <NA>     C   <NA>
    139  56.84  <NA>     O   <NA>
    140  63.37  <NA>     C   <NA>
    141  67.60  <NA>     C   <NA>
    142  72.39  <NA>     C   <NA>
    143  76.16  <NA>     O   <NA>
    144  74.55  <NA>     N   <NA>
    145  53.88  <NA>     N   <NA>
    146  49.32  <NA>     C   <NA>
    147  46.38  <NA>     C   <NA>
    148  43.18  <NA>     O   <NA>
    149  47.14  <NA>     C   <NA>
    150  46.02  <NA>     C   <NA>
    151  46.45  <NA>     C   <NA>
    152  46.45  <NA>     C   <NA>
    153  44.76  <NA>     N   <NA>
    154  44.52  <NA>     C   <NA>
    155  42.87  <NA>     C   <NA>
    156  41.39  <NA>     O   <NA>
    157  43.50  <NA>     C   <NA>
    158  46.90  <NA>     C   <NA>
    159  49.38  <NA>     C   <NA>
    160  54.07  <NA>     C   <NA>
    161  60.37  <NA>     N   <NA>
    162  41.75  <NA>     N   <NA>
    163  40.02  <NA>     C   <NA>
    164  35.18  <NA>     C   <NA>
    165  30.96  <NA>     O   <NA>
    166  45.21  <NA>     C   <NA>
    167  55.20  <NA>     C   <NA>
    168  63.88  <NA>     C   <NA>
    169  67.73  <NA>     O   <NA>
    170  60.14  <NA>     O   <NA>
    171  30.04  <NA>     N   <NA>
    172  26.79  <NA>     C   <NA>
    173  26.36  <NA>     C   <NA>
    174  24.88  <NA>     O   <NA>
    175  22.80  <NA>     C   <NA>
    176  25.58  <NA>     N   <NA>
    177  25.08  <NA>     C   <NA>
    178  23.67  <NA>     C   <NA>
    179  25.94  <NA>     O   <NA>
    180  27.31  <NA>     C   <NA>
    181  27.99  <NA>     C   <NA>
    182  29.29  <NA>     C   <NA>
    183  31.05  <NA>     C   <NA>
    184  21.67  <NA>     N   <NA>
    185  18.03  <NA>     C   <NA>
    186  19.77  <NA>     C   <NA>
    187  20.40  <NA>     O   <NA>
    188  15.36  <NA>     C   <NA>
    189  18.86  <NA>     C   <NA>
    190  16.70  <NA>     C   <NA>
    191  18.96  <NA>     C   <NA>
    192  18.83  <NA>     N   <NA>
    193  18.91  <NA>     C   <NA>
    194  16.83  <NA>     C   <NA>
    195  17.64  <NA>     O   <NA>
    196  16.55  <NA>     C   <NA>
    197  21.74  <NA>     C   <NA>
    198  22.87  <NA>     O   <NA>
    199  24.00  <NA>     O   <NA>
    200  15.44  <NA>     N   <NA>
    201  15.16  <NA>     C   <NA>
    202  15.12  <NA>     C   <NA>
    203  12.07  <NA>     O   <NA>
    204  17.69  <NA>     C   <NA>
    205  16.60  <NA>     O   <NA>
    206  13.90  <NA>     C   <NA>
    207  17.61  <NA>     N   <NA>
    208  14.06  <NA>     C   <NA>
    209  13.84  <NA>     C   <NA>
    210  14.08  <NA>     O   <NA>
    211  12.94  <NA>     N   <NA>
    212  13.20  <NA>     C   <NA>
    213  17.23  <NA>     C   <NA>
    214  15.43  <NA>     O   <NA>
    215  12.92  <NA>     C   <NA>
    216  17.65  <NA>     N   <NA>
    217  19.68  <NA>     C   <NA>
    218  21.65  <NA>     C   <NA>
    219  26.39  <NA>     O   <NA>
    220  18.61  <NA>     C   <NA>
    221  23.38  <NA>     C   <NA>
    222  23.84  <NA>     O   <NA>
    223  24.58  <NA>     O   <NA>
    224  23.67  <NA>     N   <NA>
    225  22.55  <NA>     C   <NA>
    226  20.55  <NA>     C   <NA>
    227  22.98  <NA>     O   <NA>
    228  24.31  <NA>     C   <NA>
    229  25.77  <NA>     C   <NA>
    230  37.03  <NA>     O   <NA>
    231  33.64  <NA>     O   <NA>
    232  21.65  <NA>     N   <NA>
    233  19.28  <NA>     C   <NA>
    234  20.65  <NA>     C   <NA>
    235  21.27  <NA>     O   <NA>
    236  21.21  <NA>     C   <NA>
    237  17.85  <NA>     O   <NA>
    238  20.08  <NA>     C   <NA>
    239  19.73  <NA>     N   <NA>
    240  20.18  <NA>     C   <NA>
    241  18.49  <NA>     C   <NA>
    242  21.86  <NA>     O   <NA>
    243  20.27  <NA>     C   <NA>
    244  21.42  <NA>     C   <NA>
    245  20.29  <NA>     C   <NA>
    246  21.84  <NA>     N   <NA>
    247  24.48  <NA>     C   <NA>
    248  25.37  <NA>     C   <NA>
    249  25.11  <NA>     O   <NA>
    250  24.32  <NA>     C   <NA>
    251  28.97  <NA>     C   <NA>
    252  31.03  <NA>     C   <NA>
    253  30.18  <NA>     C   <NA>
    254  29.09  <NA>     N   <NA>
    255  32.79  <NA>     C   <NA>
    256  33.96  <NA>     C   <NA>
    257  34.70  <NA>     O   <NA>
    258  36.48  <NA>     C   <NA>
    259  43.49  <NA>     C   <NA>
    260  49.89  <NA>     C   <NA>
    261  52.23  <NA>     O   <NA>
    262  53.58  <NA>     O   <NA>
    263  32.88  <NA>     N   <NA>
    264  33.85  <NA>     C   <NA>
    265  34.72  <NA>     C   <NA>
    266  34.03  <NA>     O   <NA>
    267  35.16  <NA>     C   <NA>
    268  39.66  <NA>     C   <NA>
    269  45.60  <NA>     C   <NA>
    270  50.36  <NA>     O   <NA>
    271  47.94  <NA>     O   <NA>
    272  33.77  <NA>     N   <NA>
    273  32.08  <NA>     C   <NA>
    274  33.50  <NA>     C   <NA>
    275  33.82  <NA>     O   <NA>
    276  31.77  <NA>     C   <NA>
    277  33.78  <NA>     C   <NA>
    278  38.76  <NA>     S   <NA>
    279  34.69  <NA>     C   <NA>
    280  34.23  <NA>     N   <NA>
    281  33.81  <NA>     C   <NA>
    282  32.96  <NA>     C   <NA>
    283  32.62  <NA>     O   <NA>
    284  34.02  <NA>     C   <NA>
    285  38.20  <NA>     O   <NA>
    286  35.12  <NA>     N   <NA>
    287  38.10  <NA>     C   <NA>
    288  42.61  <NA>     C   <NA>
    289  41.25  <NA>     O   <NA>
    290  35.60  <NA>     C   <NA>
    291  33.68  <NA>     C   <NA>
    292  33.91  <NA>     C   <NA>
    293  34.05  <NA>     C   <NA>
    294  45.17  <NA>     N   <NA>
    295  45.32  <NA>     C   <NA>
    296  44.13  <NA>     C   <NA>
    297  44.57  <NA>     O   <NA>
    298  48.51  <NA>     C   <NA>
    299  47.48  <NA>     C   <NA>
    300  48.21  <NA>     C   <NA>
    301  43.26  <NA>     N   <NA>
    302  44.04  <NA>     C   <NA>
    303  46.02  <NA>     C   <NA>
    304  48.40  <NA>     O   <NA>
    305  47.22  <NA>     N   <NA>
    306  48.14  <NA>     C   <NA>
    307  45.12  <NA>     C   <NA>
    308  43.73  <NA>     O   <NA>
    309  53.98  <NA>     C   <NA>
    310  61.76  <NA>     C   <NA>
    311  64.66  <NA>     C   <NA>
    312  67.95  <NA>     N   <NA>
    313  69.74  <NA>     C   <NA>
    314  69.56  <NA>     N   <NA>
    315  68.51  <NA>     N   <NA>
    316  42.00  <NA>     N   <NA>
    317  40.42  <NA>     C   <NA>
    318  40.51  <NA>     C   <NA>
    319  41.61  <NA>     O   <NA>
    320  42.08  <NA>     C   <NA>
    321  42.80  <NA>     C   <NA>
    322  46.14  <NA>     C   <NA>
    323  44.40  <NA>     C   <NA>
    324  48.15  <NA>     N   <NA>
    325  45.80  <NA>     C   <NA>
    326  44.14  <NA>     C   <NA>
    327  47.06  <NA>     C   <NA>
    328  45.02  <NA>     C   <NA>
    329  47.76  <NA>     C   <NA>
    330  39.07  <NA>     N   <NA>
    331  38.61  <NA>     C   <NA>
    332  38.09  <NA>     C   <NA>
    333  37.97  <NA>     O   <NA>
    334  40.22  <NA>     C   <NA>
    335  40.71  <NA>     C   <NA>
    336  46.25  <NA>     C   <NA>
    337  53.77  <NA>     C   <NA>
    338  55.67  <NA>     N   <NA>
    339  36.78  <NA>     N   <NA>
    340  36.73  <NA>     C   <NA>
    341  35.93  <NA>     C   <NA>
    342  35.84  <NA>     O   <NA>
    343  37.26  <NA>     C   <NA>
    344  38.05  <NA>     C   <NA>
    345  40.15  <NA>     C   <NA>
    346  34.68  <NA>     N   <NA>
    347  31.98  <NA>     C   <NA>
    348  29.00  <NA>     C   <NA>
    349  27.87  <NA>     O   <NA>
    350  31.87  <NA>     C   <NA>
    351  32.85  <NA>     C   <NA>
    352  39.51  <NA>     C   <NA>
    353  38.60  <NA>     C   <NA>
    354  45.48  <NA>     N   <NA>
    355  28.78  <NA>     N   <NA>
    356  28.75  <NA>     C   <NA>
    357  26.81  <NA>     C   <NA>
    358  28.48  <NA>     O   <NA>
    359  30.70  <NA>     C   <NA>
    360  34.30  <NA>     C   <NA>
    361  42.21  <NA>     S   <NA>
    362  43.40  <NA>     C   <NA>
    363  25.63  <NA>     N   <NA>
    364  25.64  <NA>     C   <NA>
    365  23.79  <NA>     C   <NA>
    366  25.39  <NA>     O   <NA>
    367  25.82  <NA>     C   <NA>
    368  27.41  <NA>     C   <NA>
    369  29.47  <NA>     C   <NA>
    370  33.28  <NA>     C   <NA>
    371  23.56  <NA>     N   <NA>
    372  28.25  <NA>     C   <NA>
    373  30.81  <NA>     C   <NA>
    374  30.06  <NA>     O   <NA>
    375  31.96  <NA>     N   <NA>
    376  33.92  <NA>     C   <NA>
    377  37.71  <NA>     C   <NA>
    378  36.71  <NA>     O   <NA>
    379  39.37  <NA>     N   <NA>
    380  40.56  <NA>     C   <NA>
    381  39.02  <NA>     C   <NA>
    382  40.72  <NA>     O   <NA>
    383  40.75  <NA>     C   <NA>
    384  39.18  <NA>     C   <NA>
    385  42.58  <NA>     C   <NA>
    386  39.03  <NA>     C   <NA>
    387  35.64  <NA>     N   <NA>
    388  37.36  <NA>     C   <NA>
    389  38.38  <NA>     C   <NA>
    390  42.43  <NA>     O   <NA>
    391  35.83  <NA>     N   <NA>
    392  34.71  <NA>     C   <NA>
    393  32.77  <NA>     C   <NA>
    394  33.65  <NA>     O   <NA>
    395  31.50  <NA>     N   <NA>
    396  30.76  <NA>     C   <NA>
    397  31.86  <NA>     C   <NA>
    398  36.37  <NA>     O   <NA>
    399  28.77  <NA>     C   <NA>
    400  30.51  <NA>     C   <NA>
    401  27.94  <NA>     C   <NA>
    402  30.74  <NA>     C   <NA>
    403  29.27  <NA>     C   <NA>
    404  32.36  <NA>     C   <NA>
    405  28.01  <NA>     C   <NA>
    406  28.86  <NA>     N   <NA>
    407  28.83  <NA>     C   <NA>
    408  29.53  <NA>     C   <NA>
    409  28.42  <NA>     O   <NA>
    410  28.74  <NA>     C   <NA>
    411  27.76  <NA>     C   <NA>
    412  28.78  <NA>     C   <NA>
    413  23.65  <NA>     C   <NA>
    414  31.86  <NA>     N   <NA>
    415  31.37  <NA>     C   <NA>
    416  31.40  <NA>     C   <NA>
    417  30.39  <NA>     O   <NA>
    418  34.83  <NA>     C   <NA>
    419  41.25  <NA>     C   <NA>
    420  48.65  <NA>     C   <NA>
    421  52.68  <NA>     C   <NA>
    422  56.46  <NA>     N   <NA>
    423  27.17  <NA>     N   <NA>
    424  26.79  <NA>     C   <NA>
    425  25.98  <NA>     C   <NA>
    426  28.45  <NA>     O   <NA>
    427  23.58  <NA>     C   <NA>
    428  24.00  <NA>     C   <NA>
    429  19.88  <NA>     C   <NA>
    430  27.52  <NA>     N   <NA>
    431  29.29  <NA>     C   <NA>
    432  26.53  <NA>     C   <NA>
    433  26.51  <NA>     O   <NA>
    434  33.25  <NA>     C   <NA>
    435  39.60  <NA>     C   <NA>
    436  44.47  <NA>     C   <NA>
    437  50.18  <NA>     N   <NA>
    438  49.84  <NA>     C   <NA>
    439  51.65  <NA>     N   <NA>
    440  52.43  <NA>     N   <NA>
    441  25.48  <NA>     N   <NA>
    442  27.73  <NA>     C   <NA>
    443  30.15  <NA>     C   <NA>
    444  33.75  <NA>     O   <NA>
    445  26.15  <NA>     C   <NA>
    446  29.09  <NA>     C   <NA>
    447  37.13  <NA>     C   <NA>
    448  40.89  <NA>     O   <NA>
    449  39.09  <NA>     N   <NA>
    450  32.83  <NA>     N   <NA>
    451  34.57  <NA>     C   <NA>
    452  35.86  <NA>     C   <NA>
    453  38.78  <NA>     O   <NA>
    454  33.06  <NA>     C   <NA>
    455  36.56  <NA>     C   <NA>
    456  37.53  <NA>     C   <NA>
    457  33.33  <NA>     C   <NA>
    458  34.69  <NA>     C   <NA>
    459  33.09  <NA>     C   <NA>
    460  33.21  <NA>     C   <NA>
    461  40.16  <NA>     O   <NA>
    462  37.21  <NA>     N   <NA>
    463  36.50  <NA>     C   <NA>
    464  35.62  <NA>     C   <NA>
    465  33.30  <NA>     O   <NA>
    466  43.82  <NA>     C   <NA>
    467  49.03  <NA>     C   <NA>
    468  53.15  <NA>     O   <NA>
    469  54.18  <NA>     O   <NA>
    470  35.15  <NA>     N   <NA>
    471  37.27  <NA>     C   <NA>
    472  36.83  <NA>     C   <NA>
    473  39.02  <NA>     O   <NA>
    474  41.54  <NA>     C   <NA>
    475  53.40  <NA>     C   <NA>
    476  61.08  <NA>     C   <NA>
    477  65.46  <NA>     O   <NA>
    478  58.46  <NA>     N   <NA>
    479  33.10  <NA>     N   <NA>
    480  33.95  <NA>     C   <NA>
    481  35.70  <NA>     C   <NA>
    482  38.95  <NA>     O   <NA>
    483  33.50  <NA>     C   <NA>
    484  31.68  <NA>     C   <NA>
    485  33.86  <NA>     C   <NA>
    486  35.62  <NA>     C   <NA>
    487  38.14  <NA>     N   <NA>
    488  41.38  <NA>     C   <NA>
    489  42.34  <NA>     C   <NA>
    490  43.91  <NA>     O   <NA>
    491  41.40  <NA>     C   <NA>
    492  42.58  <NA>     C   <NA>
    493  46.33  <NA>     C   <NA>
    494  45.30  <NA>     C   <NA>
    495  44.21  <NA>     N   <NA>
    496  46.55  <NA>     C   <NA>
    497  46.30  <NA>     C   <NA>
    498  45.60  <NA>     O   <NA>
    499  48.42  <NA>     C   <NA>
    500  51.56  <NA>     C   <NA>
    501  45.78  <NA>     C   <NA>
    502  52.26  <NA>     C   <NA>
    503  48.38  <NA>     N   <NA>
    504  52.65  <NA>     C   <NA>
    505  53.34  <NA>     C   <NA>
    506  57.14  <NA>     O   <NA>
    507  53.39  <NA>     C   <NA>
    508  58.82  <NA>     C   <NA>
    509  61.85  <NA>     C   <NA>
    510  62.51  <NA>     O   <NA>
    511  63.85  <NA>     O   <NA>
    512  54.18  <NA>     N   <NA>
    513  53.69  <NA>     C   <NA>
    514  55.63  <NA>     C   <NA>
    515  55.17  <NA>     O   <NA>
    516  53.26  <NA>     C   <NA>
    517  50.35  <NA>     C   <NA>
    518  51.36  <NA>     C   <NA>
    519  49.54  <NA>     C   <NA>
    520  57.47  <NA>     N   <NA>
    521  59.30  <NA>     C   <NA>
    522  58.98  <NA>     C   <NA>
    523  59.83  <NA>     O   <NA>
    524  60.45  <NA>     C   <NA>
    525  66.26  <NA>     S   <NA>
    526  58.28  <NA>     N   <NA>
    527  56.53  <NA>     C   <NA>
    528  54.84  <NA>     C   <NA>
    529  55.20  <NA>     O   <NA>
    530  52.63  <NA>     N   <NA>
    531  49.83  <NA>     C   <NA>
    532  49.49  <NA>     C   <NA>
    533  50.24  <NA>     O   <NA>
    534  48.37  <NA>     C   <NA>
    535  48.55  <NA>     C   <NA>
    536  46.92  <NA>     N   <NA>
    537  47.07  <NA>     C   <NA>
    538  45.30  <NA>     C   <NA>
    539  49.32  <NA>     N   <NA>
    540  48.04  <NA>     N   <NA>
    541  44.00  <NA>     C   <NA>
    542  40.53  <NA>     C   <NA>
    543  37.44  <NA>     O   <NA>
    544  46.06  <NA>     C   <NA>
    545  48.65  <NA>     C   <NA>
    546  52.10  <NA>     C   <NA>
    547  52.81  <NA>     C   <NA>
    548  57.49  <NA>     N   <NA>
    549  37.00  <NA>     N   <NA>
    550  33.57  <NA>     C   <NA>
    551  31.56  <NA>     C   <NA>
    552  34.29  <NA>     O   <NA>
    553  28.40  <NA>     C   <NA>
    554  28.94  <NA>     N   <NA>
    555  26.93  <NA>     C   <NA>
    556  25.42  <NA>     C   <NA>
    557  23.83  <NA>     O   <NA>
    558  34.11  <NA>     C   <NA>
    559  39.50  <NA>     C   <NA>
    560  34.28  <NA>     C   <NA>
    561  44.26  <NA>     C   <NA>
    562  23.40  <NA>     N   <NA>
    563  22.07  <NA>     C   <NA>
    564  24.48  <NA>     C   <NA>
    565  24.29  <NA>     O   <NA>
    566  23.67  <NA>     N   <NA>
    567  24.46  <NA>     C   <NA>
    568  25.63  <NA>     C   <NA>
    569  25.67  <NA>     O   <NA>
    570  24.94  <NA>     C   <NA>
    571  29.61  <NA>     O   <NA>
    572  19.06  <NA>     C   <NA>
    573  25.14  <NA>     N   <NA>
    574  26.87  <NA>     C   <NA>
    575  27.97  <NA>     C   <NA>
    576  27.87  <NA>     O   <NA>
    577  25.74  <NA>     C   <NA>
    578  27.82  <NA>     C   <NA>
    579  21.53  <NA>     C   <NA>
    580  27.21  <NA>     N   <NA>
    581  25.80  <NA>     C   <NA>
    582  24.05  <NA>     C   <NA>
    583  23.89  <NA>     O   <NA>
    584  24.57  <NA>     C   <NA>
    585  24.94  <NA>     C   <NA>
    586  21.47  <NA>     C   <NA>
    587  22.70  <NA>     C   <NA>
    588  22.25  <NA>     N   <NA>
    589  25.77  <NA>     C   <NA>
    590  25.41  <NA>     C   <NA>
    591  21.46  <NA>     O   <NA>
    592  26.64  <NA>     C   <NA>
    593  24.10  <NA>     C   <NA>
    594  24.74  <NA>     C   <NA>
    595  23.37  <NA>     N   <NA>
    596  25.28  <NA>     C   <NA>
    597  24.69  <NA>     C   <NA>
    598  27.34  <NA>     O   <NA>
    599  24.96  <NA>     N   <NA>
    600  26.27  <NA>     C   <NA>
    601  27.30  <NA>     C   <NA>
    602  30.31  <NA>     O   <NA>
    603  26.15  <NA>     C   <NA>
    604  27.52  <NA>     C   <NA>
    605  27.77  <NA>     C   <NA>
    606  27.82  <NA>     N   <NA>
    607  27.04  <NA>     C   <NA>
    608  28.49  <NA>     C   <NA>
    609  33.89  <NA>     O   <NA>
    610  24.34  <NA>     C   <NA>
    611  28.59  <NA>     O   <NA>
    612  17.97  <NA>     C   <NA>
    613  27.56  <NA>     N   <NA>
    614  28.35  <NA>     C   <NA>
    615  32.18  <NA>     C   <NA>
    616  35.74  <NA>     O   <NA>
    617  28.05  <NA>     C   <NA>
    618  24.27  <NA>     C   <NA>
    619  24.21  <NA>     C   <NA>
    620  31.78  <NA>     N   <NA>
    621  30.36  <NA>     C   <NA>
    622  29.01  <NA>     C   <NA>
    623  28.83  <NA>     O   <NA>
    624  27.70  <NA>     C   <NA>
    625  25.82  <NA>     C   <NA>
    626  29.05  <NA>     C   <NA>
    627  28.44  <NA>     N   <NA>
    628  23.14  <NA>     C   <NA>
    629  20.40  <NA>     C   <NA>
    630  22.72  <NA>     O   <NA>
    631  23.96  <NA>     C   <NA>
    632  19.36  <NA>     C   <NA>
    633  24.75  <NA>     O   <NA>
    634  23.07  <NA>     N   <NA>
    635  19.46  <NA>     N   <NA>
    636  18.70  <NA>     C   <NA>
    637  15.24  <NA>     C   <NA>
    638  16.01  <NA>     O   <NA>
    639  18.01  <NA>     C   <NA>
    640  21.46  <NA>     C   <NA>
    641  18.09  <NA>     C   <NA>
    642  21.77  <NA>     C   <NA>
    643  17.41  <NA>     N   <NA>
    644  17.84  <NA>     C   <NA>
    645  18.00  <NA>     C   <NA>
    646  16.33  <NA>     O   <NA>
    647  19.87  <NA>     C   <NA>
    648  21.72  <NA>     C   <NA>
    649  19.84  <NA>     C   <NA>
    650  23.27  <NA>     C   <NA>
    651  17.67  <NA>     N   <NA>
    652  19.47  <NA>     C   <NA>
    653  20.59  <NA>     C   <NA>
    654  19.31  <NA>     O   <NA>
    655  19.04  <NA>     N   <NA>
    656  21.43  <NA>     C   <NA>
    657  20.84  <NA>     C   <NA>
    658  24.40  <NA>     O   <NA>
    659  15.37  <NA>     C   <NA>
    660  18.83  <NA>     C   <NA>
    661  21.67  <NA>     C   <NA>
    662  21.37  <NA>     N   <NA>
    663  22.28  <NA>     C   <NA>
    664  20.00  <NA>     N   <NA>
    665  21.02  <NA>     N   <NA>
    666  21.72  <NA>     N   <NA>
    667  19.12  <NA>     C   <NA>
    668  21.99  <NA>     C   <NA>
    669  24.63  <NA>     O   <NA>
    670  17.71  <NA>     C   <NA>
    671  23.05  <NA>     C   <NA>
    672  24.51  <NA>     O   <NA>
    673  21.17  <NA>     N   <NA>
    674  21.73  <NA>     N   <NA>
    675  21.46  <NA>     C   <NA>
    676  21.01  <NA>     C   <NA>
    677  22.81  <NA>     O   <NA>
    678  24.37  <NA>     C   <NA>
    679  26.43  <NA>     C   <NA>
    680  27.76  <NA>     C   <NA>
    681  29.10  <NA>     C   <NA>
    682  18.04  <NA>     N   <NA>
    683  21.90  <NA>     C   <NA>
    684  22.24  <NA>     C   <NA>
    685  20.42  <NA>     O   <NA>
    686  21.44  <NA>     C   <NA>
    687  19.50  <NA>     C   <NA>
    688  20.48  <NA>     C   <NA>
    689  15.96  <NA>     C   <NA>
    690  20.88  <NA>     N   <NA>
    691  25.31  <NA>     C   <NA>
    692  25.69  <NA>     C   <NA>
    693  29.29  <NA>     O   <NA>
    694  21.94  <NA>     C   <NA>
    695  20.91  <NA>     O   <NA>
    696  24.44  <NA>     C   <NA>
    697  23.55  <NA>     N   <NA>
    698  23.15  <NA>     C   <NA>
    699  25.63  <NA>     C   <NA>
    700  29.05  <NA>     O   <NA>
    701  17.48  <NA>     C   <NA>
    702  12.58  <NA>     C   <NA>
    703  23.01  <NA>     C   <NA>
    704  24.98  <NA>     O   <NA>
    705  15.18  <NA>     N   <NA>
    706  25.00  <NA>     N   <NA>
    707  25.37  <NA>     C   <NA>
    708  24.18  <NA>     C   <NA>
    709  26.26  <NA>     O   <NA>
    710  24.44  <NA>     C   <NA>
    711  24.46  <NA>     C   <NA>
    712  22.23  <NA>     C   <NA>
    713  24.68  <NA>     C   <NA>
    714  24.04  <NA>     N   <NA>
    715  25.09  <NA>     C   <NA>
    716  27.08  <NA>     C   <NA>
    717  28.80  <NA>     O   <NA>
    718  27.02  <NA>     N   <NA>
    719  26.34  <NA>     C   <NA>
    720  26.35  <NA>     C   <NA>
    721  28.34  <NA>     O   <NA>
    722  27.62  <NA>     C   <NA>
    723  32.32  <NA>     S   <NA>
    724  27.53  <NA>     N   <NA>
    725  27.51  <NA>     C   <NA>
    726  27.27  <NA>     C   <NA>
    727  27.30  <NA>     O   <NA>
    728  28.47  <NA>     C   <NA>
    729  28.16  <NA>     O   <NA>
    730  26.73  <NA>     C   <NA>
    731  28.68  <NA>     N   <NA>
    732  26.62  <NA>     C   <NA>
    733  25.49  <NA>     C   <NA>
    734  23.96  <NA>     O   <NA>
    735  24.94  <NA>     C   <NA>
    736  28.32  <NA>     C   <NA>
    737  27.79  <NA>     C   <NA>
    738  26.48  <NA>     C   <NA>
    739  26.29  <NA>     N   <NA>
    740  29.12  <NA>     C   <NA>
    741  29.85  <NA>     C   <NA>
    742  26.89  <NA>     O   <NA>
    743  33.75  <NA>     C   <NA>
    744  40.13  <NA>     C   <NA>
    745  43.34  <NA>     O   <NA>
    746  42.07  <NA>     N   <NA>
    747  31.14  <NA>     N   <NA>
    748  35.14  <NA>     C   <NA>
    749  34.93  <NA>     C   <NA>
    750  36.66  <NA>     O   <NA>
    751  32.92  <NA>     C   <NA>
    752  31.47  <NA>     C   <NA>
    753  29.83  <NA>     C   <NA>
    754  28.13  <NA>     C   <NA>
    755  27.98  <NA>     C   <NA>
    756  25.49  <NA>     C   <NA>
    757  27.25  <NA>     C   <NA>
    758  48.12  <NA>     N   <NA>
    759  43.36  <NA>     C   <NA>
    760  39.59  <NA>     C   <NA>
    761  37.70  <NA>     O   <NA>
    762  46.58  <NA>     C   <NA>
    763  48.47  <NA>     C   <NA>
    764  50.98  <NA>     C   <NA>
    765  36.85  <NA>     N   <NA>
    766  37.15  <NA>     C   <NA>
    767  36.43  <NA>     C   <NA>
    768  39.41  <NA>     O   <NA>
    769  38.60  <NA>     C   <NA>
    770  39.92  <NA>     C   <NA>
    771  44.52  <NA>     C   <NA>
    772  48.57  <NA>     O   <NA>
    773  45.25  <NA>     N   <NA>
    774  32.18  <NA>     N   <NA>
    775  30.91  <NA>     C   <NA>
    776  28.84  <NA>     C   <NA>
    777  28.43  <NA>     O   <NA>
    778  33.01  <NA>     C   <NA>
    779  35.76  <NA>     C   <NA>
    780  31.17  <NA>     C   <NA>
    781  36.59  <NA>     C   <NA>
    782  24.52  <NA>     N   <NA>
    783  26.46  <NA>     C   <NA>
    784  23.26  <NA>     C   <NA>
    785  23.16  <NA>     O   <NA>
    786  30.32  <NA>     C   <NA>
    787  32.67  <NA>     O   <NA>
    788  30.57  <NA>     C   <NA>
    789  21.95  <NA>     N   <NA>
    790  21.45  <NA>     C   <NA>
    791  22.43  <NA>     C   <NA>
    792  22.26  <NA>     O   <NA>
    793  18.69  <NA>     C   <NA>
    794  18.41  <NA>     C   <NA>
    795  16.32  <NA>     C   <NA>
    796  16.07  <NA>     C   <NA>
    797  21.66  <NA>     N   <NA>
    798  21.58  <NA>     C   <NA>
    799  21.21  <NA>     C   <NA>
    800  22.05  <NA>     O   <NA>
    801  19.59  <NA>     C   <NA>
    802  18.40  <NA>     C   <NA>
    803  17.80  <NA>     C   <NA>
    804  14.37  <NA>     C   <NA>
    805  17.13  <NA>     N   <NA>
    806  17.08  <NA>     C   <NA>
    807  17.33  <NA>     C   <NA>
    808  16.59  <NA>     C   <NA>
    809  15.43  <NA>     C   <NA>
    810  18.86  <NA>     C   <NA>
    811  23.30  <NA>     N   <NA>
    812  24.55  <NA>     C   <NA>
    813  23.64  <NA>     C   <NA>
    814  23.66  <NA>     O   <NA>
    815  32.55  <NA>     C   <NA>
    816  41.57  <NA>     C   <NA>
    817  50.90  <NA>     C   <NA>
    818  55.87  <NA>     O   <NA>
    819  51.09  <NA>     N   <NA>
    820  21.29  <NA>     N   <NA>
    821  20.70  <NA>     C   <NA>
    822  19.16  <NA>     C   <NA>
    823  22.46  <NA>     O   <NA>
    824  21.59  <NA>     C   <NA>
    825  27.03  <NA>     C   <NA>
    826  28.48  <NA>     C   <NA>
    827  38.82  <NA>     N   <NA>
    828  45.44  <NA>     C   <NA>
    829  47.55  <NA>     N   <NA>
    830  49.71  <NA>     N   <NA>
    831  18.31  <NA>     N   <NA>
    832  16.67  <NA>     C   <NA>
    833  18.20  <NA>     C   <NA>
    834  17.79  <NA>     O   <NA>
    835  15.13  <NA>     C   <NA>
    836  13.81  <NA>     C   <NA>
    837  15.31  <NA>     C   <NA>
    838  18.74  <NA>     N   <NA>
    839  24.75  <NA>     C   <NA>
    840  28.33  <NA>     C   <NA>
    841  34.15  <NA>     O   <NA>
    842  22.30  <NA>     C   <NA>
    843  26.19  <NA>     C   <NA>
    844  26.68  <NA>     C   <NA>
    845  25.72  <NA>     C   <NA>
    846  29.04  <NA>     N   <NA>
    847  25.94  <NA>     C   <NA>
    848  28.64  <NA>     C   <NA>
    849  26.28  <NA>     O   <NA>
    850  26.71  <NA>     C   <NA>
    851  26.27  <NA>     C   <NA>
    852  23.17  <NA>     C   <NA>
    853  29.16  <NA>     N   <NA>
    854  29.59  <NA>     C   <NA>
    855  27.37  <NA>     C   <NA>
    856  28.24  <NA>     O   <NA>
    857  33.38  <NA>     C   <NA>
    858  41.16  <NA>     O   <NA>
    859  37.24  <NA>     C   <NA>
    860  27.04  <NA>     N   <NA>
    861  26.55  <NA>     C   <NA>
    862  28.05  <NA>     C   <NA>
    863  25.03  <NA>     O   <NA>
    864  25.84  <NA>     C   <NA>
    865  28.87  <NA>     C   <NA>
    866  20.54  <NA>     C   <NA>
    867  29.32  <NA>     C   <NA>
    868  27.45  <NA>     N   <NA>
    869  26.92  <NA>     C   <NA>
    870  26.01  <NA>     C   <NA>
    871  24.75  <NA>     O   <NA>
    872  25.93  <NA>     C   <NA>
    873  27.04  <NA>     C   <NA>
    874  32.12  <NA>     C   <NA>
    875  34.02  <NA>     C   <NA>
    876  41.05  <NA>     N   <NA>
    877  26.74  <NA>     N   <NA>
    878  30.24  <NA>     C   <NA>
    879  33.16  <NA>     C   <NA>
    880  31.38  <NA>     O   <NA>
    881  31.34  <NA>     C   <NA>
    882  32.33  <NA>     C   <NA>
    883  28.02  <NA>     C   <NA>
    884  33.84  <NA>     C   <NA>
    885  36.75  <NA>     N   <NA>
    886  39.96  <NA>     C   <NA>
    887  40.86  <NA>     C   <NA>
    888  42.44  <NA>     O   <NA>
    889  40.79  <NA>     N   <NA>
    890  37.08  <NA>     C   <NA>
    891  35.82  <NA>     C   <NA>
    892  37.78  <NA>     O   <NA>
    893  32.85  <NA>     N   <NA>
    894  35.03  <NA>     C   <NA>
    895  34.53  <NA>     C   <NA>
    896  33.51  <NA>     O   <NA>
    897  40.01  <NA>     C   <NA>
    898  47.43  <NA>     C   <NA>
    899  53.35  <NA>     C   <NA>
    900  55.74  <NA>     O   <NA>
    901  54.03  <NA>     N   <NA>
    902  32.63  <NA>     N   <NA>
    903  29.37  <NA>     C   <NA>
    904  28.10  <NA>     C   <NA>
    905  29.04  <NA>     O   <NA>
    906  28.40  <NA>     C   <NA>
    907  26.42  <NA>     C   <NA>
    908  22.91  <NA>     C   <NA>
    909  23.99  <NA>     C   <NA>
    910  28.58  <NA>     N   <NA>
    911  25.77  <NA>     C   <NA>
    912  26.35  <NA>     C   <NA>
    913  25.15  <NA>     O   <NA>
    914  29.54  <NA>     C   <NA>
    915  36.66  <NA>     C   <NA>
    916  41.68  <NA>     C   <NA>
    917  45.08  <NA>     C   <NA>
    918  46.33  <NA>     N   <NA>
    919  24.88  <NA>     N   <NA>
    920  25.39  <NA>     C   <NA>
    921  23.31  <NA>     C   <NA>
    922  25.34  <NA>     O   <NA>
    923  32.42  <NA>     C   <NA>
    924  47.17  <NA>     C   <NA>
    925  57.11  <NA>     C   <NA>
    926  62.83  <NA>     O   <NA>
    927  64.02  <NA>     O   <NA>
    928  16.20  <NA>     N   <NA>
    929  14.52  <NA>     C   <NA>
    930  15.55  <NA>     C   <NA>
    931  17.95  <NA>     O   <NA>
    932  11.99  <NA>     C   <NA>
    933  14.37  <NA>     N   <NA>
    934  15.39  <NA>     C   <NA>
    935  15.08  <NA>     C   <NA>
    936  14.05  <NA>     O   <NA>
    937  16.68  <NA>     C   <NA>
    938  20.94  <NA>     C   <NA>
    939  18.77  <NA>     C   <NA>
    940  21.23  <NA>     C   <NA>
    941  13.16  <NA>     N   <NA>
    942  17.92  <NA>     C   <NA>
    943  18.17  <NA>     C   <NA>
    944  19.32  <NA>     O   <NA>
    945  17.42  <NA>     C   <NA>
    946  19.17  <NA>     C   <NA>
    947  20.12  <NA>     C   <NA>
    948  19.21  <NA>     C   <NA>
    949  19.08  <NA>     N   <NA>
    950  18.20  <NA>     C   <NA>
    951  19.09  <NA>     C   <NA>
    952  18.68  <NA>     O   <NA>
    953  18.15  <NA>     C   <NA>
    954  20.90  <NA>     C   <NA>
    955  23.13  <NA>     O   <NA>
    956  22.18  <NA>     O   <NA>
    957  14.92  <NA>     N   <NA>
    958  13.68  <NA>     C   <NA>
    959  15.84  <NA>     C   <NA>
    960  18.69  <NA>     O   <NA>
    961  15.43  <NA>     C   <NA>
    962  15.52  <NA>     O   <NA>
    963  15.01  <NA>     C   <NA>
    964  13.19  <NA>     N   <NA>
    965  13.54  <NA>     C   <NA>
    966  17.44  <NA>     C   <NA>
    967  18.85  <NA>     O   <NA>
    968  17.91  <NA>     N   <NA>
    969  18.89  <NA>     C   <NA>
    970  20.66  <NA>     C   <NA>
    971  22.38  <NA>     O   <NA>
    972  13.77  <NA>     C   <NA>
    973  20.64  <NA>     N   <NA>
    974  14.85  <NA>     C   <NA>
    975  14.01  <NA>     C   <NA>
    976  17.76  <NA>     O   <NA>
    977  15.99  <NA>     C   <NA>
    978  22.55  <NA>     C   <NA>
    979  27.81  <NA>     O   <NA>
    980  29.43  <NA>     O   <NA>
    981  16.70  <NA>     N   <NA>
    982  15.10  <NA>     C   <NA>
    983  16.38  <NA>     C   <NA>
    984  17.11  <NA>     O   <NA>
    985  21.76  <NA>     C   <NA>
    986  25.41  <NA>     C   <NA>
    987  25.58  <NA>     O   <NA>
    988  27.89  <NA>     O   <NA>
    989  18.44  <NA>     N   <NA>
    990  15.45  <NA>     C   <NA>
    991  19.23  <NA>     C   <NA>
    992  18.95  <NA>     O   <NA>
    993  12.28  <NA>     C   <NA>
    994  13.54  <NA>     O   <NA>
    995  10.95  <NA>     C   <NA>
    996  19.53  <NA>     N   <NA>
    997  18.71  <NA>     C   <NA>
    998  19.01  <NA>     C   <NA>
    999  19.90  <NA>     O   <NA>
    1000 20.49  <NA>     C   <NA>
    1001 26.38  <NA>     C   <NA>
    1002 20.87  <NA>     C   <NA>
    1003 22.78  <NA>     N   <NA>
    1004 23.84  <NA>     C   <NA>
    1005 27.19  <NA>     C   <NA>
    1006 26.71  <NA>     O   <NA>
    1007 21.02  <NA>     C   <NA>
    1008 20.79  <NA>     C   <NA>
    1009 21.64  <NA>     C   <NA>
    1010 23.27  <NA>     C   <NA>
    1011 26.46  <NA>     N   <NA>
    1012 26.40  <NA>     C   <NA>
    1013 26.66  <NA>     C   <NA>
    1014 26.75  <NA>     O   <NA>
    1015 24.89  <NA>     C   <NA>
    1016 31.34  <NA>     C   <NA>
    1017 32.62  <NA>     C   <NA>
    1018 39.77  <NA>     O   <NA>
    1019 39.65  <NA>     O   <NA>
    1020 30.87  <NA>     N   <NA>
    1021 32.05  <NA>     C   <NA>
    1022 32.29  <NA>     C   <NA>
    1023 37.80  <NA>     O   <NA>
    1024 35.17  <NA>     C   <NA>
    1025 44.03  <NA>     C   <NA>
    1026 44.48  <NA>     C   <NA>
    1027 38.16  <NA>     O   <NA>
    1028 46.13  <NA>     O   <NA>
    1029 31.85  <NA>     N   <NA>
    1030 33.53  <NA>     C   <NA>
    1031 33.93  <NA>     C   <NA>
    1032 36.78  <NA>     O   <NA>
    1033 35.97  <NA>     C   <NA>
    1034 35.75  <NA>     C   <NA>
    1035 43.16  <NA>     S   <NA>
    1036 42.89  <NA>     C   <NA>
    1037 33.54  <NA>     N   <NA>
    1038 35.84  <NA>     C   <NA>
    1039 34.42  <NA>     C   <NA>
    1040 33.84  <NA>     O   <NA>
    1041 38.59  <NA>     C   <NA>
    1042 41.47  <NA>     O   <NA>
    1043 35.38  <NA>     N   <NA>
    1044 34.72  <NA>     C   <NA>
    1045 36.10  <NA>     C   <NA>
    1046 35.71  <NA>     O   <NA>
    1047 30.77  <NA>     C   <NA>
    1048 26.92  <NA>     C   <NA>
    1049 29.33  <NA>     C   <NA>
    1050 27.40  <NA>     C   <NA>
    1051 39.43  <NA>     N   <NA>
    1052 40.91  <NA>     C   <NA>
    1053 41.38  <NA>     C   <NA>
    1054 42.07  <NA>     O   <NA>
    1055 41.75  <NA>     C   <NA>
    1056 41.42  <NA>     C   <NA>
    1057 42.37  <NA>     C   <NA>
    1058 35.62  <NA>     N   <NA>
    1059 33.54  <NA>     C   <NA>
    1060 29.12  <NA>     C   <NA>
    1061 28.38  <NA>     O   <NA>
    1062 27.66  <NA>     N   <NA>
    1063 28.95  <NA>     C   <NA>
    1064 30.49  <NA>     C   <NA>
    1065 31.54  <NA>     O   <NA>
    1066 28.94  <NA>     C   <NA>
    1067 34.11  <NA>     C   <NA>
    1068 34.17  <NA>     C   <NA>
    1069 38.61  <NA>     N   <NA>
    1070 35.13  <NA>     C   <NA>
    1071 27.84  <NA>     N   <NA>
    1072 27.82  <NA>     N   <NA>
    1073 27.94  <NA>     N   <NA>
    1074 28.18  <NA>     C   <NA>
    1075 30.79  <NA>     C   <NA>
    1076 28.97  <NA>     O   <NA>
    1077 24.28  <NA>     C   <NA>
    1078 18.74  <NA>     C   <NA>
    1079 16.49  <NA>     C   <NA>
    1080 18.74  <NA>     C   <NA>
    1081 16.06  <NA>     N   <NA>
    1082 14.09  <NA>     C   <NA>
    1083 20.09  <NA>     C   <NA>
    1084 18.75  <NA>     C   <NA>
    1085 22.86  <NA>     C   <NA>
    1086 20.90  <NA>     C   <NA>
    1087 31.89  <NA>     N   <NA>
    1088 30.18  <NA>     C   <NA>
    1089 27.53  <NA>     C   <NA>
    1090 25.86  <NA>     O   <NA>
    1091 33.16  <NA>     C   <NA>
    1092 40.56  <NA>     C   <NA>
    1093 49.90  <NA>     C   <NA>
    1094 52.53  <NA>     C   <NA>
    1095 52.90  <NA>     N   <NA>
    1096 25.82  <NA>     N   <NA>
    1097 25.30  <NA>     C   <NA>
    1098 27.71  <NA>     C   <NA>
    1099 25.94  <NA>     O   <NA>
    1100 25.59  <NA>     C   <NA>
    1101 26.07  <NA>     C   <NA>
    1102 27.58  <NA>     C   <NA>
    1103 29.96  <NA>     N   <NA>
    1104 28.39  <NA>     C   <NA>
    1105 27.43  <NA>     C   <NA>
    1106 26.93  <NA>     O   <NA>
    1107 27.50  <NA>     C   <NA>
    1108 27.69  <NA>     C   <NA>
    1109 25.34  <NA>     C   <NA>
    1110 32.56  <NA>     C   <NA>
    1111 32.92  <NA>     N   <NA>
    1112 29.02  <NA>     N   <NA>
    1113 31.16  <NA>     C   <NA>
    1114 29.75  <NA>     C   <NA>
    1115 29.99  <NA>     O   <NA>
    1116 38.22  <NA>     C   <NA>
    1117 45.23  <NA>     C   <NA>
    1118 52.02  <NA>     S   <NA>
    1119 49.41  <NA>     C   <NA>
    1120 30.87  <NA>     N   <NA>
    1121 29.85  <NA>     C   <NA>
    1122 27.34  <NA>     C   <NA>
    1123 26.79  <NA>     O   <NA>
    1124 29.30  <NA>     C   <NA>
    1125 25.99  <NA>     C   <NA>
    1126 28.73  <NA>     C   <NA>
    1127 27.41  <NA>     C   <NA>
    1128 31.53  <NA>     N   <NA>
    1129 31.93  <NA>     C   <NA>
    1130 33.10  <NA>     C   <NA>
    1131 30.88  <NA>     O   <NA>
    1132 32.08  <NA>     N   <NA>
    1133 34.90  <NA>     C   <NA>
    1134 35.58  <NA>     C   <NA>
    1135 36.70  <NA>     O   <NA>
    1136 35.37  <NA>     N   <NA>
    1137 35.93  <NA>     C   <NA>
    1138 35.51  <NA>     C   <NA>
    1139 34.79  <NA>     O   <NA>
    1140 36.57  <NA>     C   <NA>
    1141 38.14  <NA>     C   <NA>
    1142 34.90  <NA>     C   <NA>
    1143 36.99  <NA>     C   <NA>
    1144 32.76  <NA>     N   <NA>
    1145 34.23  <NA>     C   <NA>
    1146 37.16  <NA>     C   <NA>
    1147 40.75  <NA>     O   <NA>
    1148 36.36  <NA>     N   <NA>
    1149 35.14  <NA>     C   <NA>
    1150 35.06  <NA>     C   <NA>
    1151 33.99  <NA>     O   <NA>
    1152 33.56  <NA>     N   <NA>
    1153 33.93  <NA>     C   <NA>
    1154 33.07  <NA>     C   <NA>
    1155 35.98  <NA>     O   <NA>
    1156 35.25  <NA>     C   <NA>
    1157 39.11  <NA>     C   <NA>
    1158 40.72  <NA>     C   <NA>
    1159 40.77  <NA>     C   <NA>
    1160 42.79  <NA>     C   <NA>
    1161 41.95  <NA>     C   <NA>
    1162 44.36  <NA>     C   <NA>
    1163 30.50  <NA>     N   <NA>
    1164 28.82  <NA>     C   <NA>
    1165 28.68  <NA>     C   <NA>
    1166 30.25  <NA>     O   <NA>
    1167 25.49  <NA>     C   <NA>
    1168 24.37  <NA>     C   <NA>
    1169 27.89  <NA>     C   <NA>
    1170 16.71  <NA>     C   <NA>
    1171 29.03  <NA>     N   <NA>
    1172 31.48  <NA>     C   <NA>
    1173 28.45  <NA>     C   <NA>
    1174 28.39  <NA>     O   <NA>
    1175 38.52  <NA>     C   <NA>
    1176 45.52  <NA>     C   <NA>
    1177 51.68  <NA>     C   <NA>
    1178 55.32  <NA>     C   <NA>
    1179 51.46  <NA>     N   <NA>
    1180 25.88  <NA>     N   <NA>
    1181 26.87  <NA>     C   <NA>
    1182 26.06  <NA>     C   <NA>
    1183 29.12  <NA>     O   <NA>
    1184 27.47  <NA>     C   <NA>
    1185 27.66  <NA>     C   <NA>
    1186 29.79  <NA>     C   <NA>
    1187 24.26  <NA>     N   <NA>
    1188 21.61  <NA>     C   <NA>
    1189 20.76  <NA>     C   <NA>
    1190 24.14  <NA>     O   <NA>
    1191 20.35  <NA>     C   <NA>
    1192 26.15  <NA>     C   <NA>
    1193 24.10  <NA>     C   <NA>
    1194 34.58  <NA>     N   <NA>
    1195 31.59  <NA>     C   <NA>
    1196 33.77  <NA>     N   <NA>
    1197 29.56  <NA>     N   <NA>
    1198 19.42  <NA>     N   <NA>
    1199 21.62  <NA>     C   <NA>
    1200 23.18  <NA>     C   <NA>
    1201 26.34  <NA>     O   <NA>
    1202 20.63  <NA>     C   <NA>
    1203 18.81  <NA>     C   <NA>
    1204 20.16  <NA>     C   <NA>
    1205 28.16  <NA>     O   <NA>
    1206 21.91  <NA>     N   <NA>
    1207 21.35  <NA>     N   <NA>
    1208 19.01  <NA>     C   <NA>
    1209 21.82  <NA>     C   <NA>
    1210 25.32  <NA>     O   <NA>
    1211 18.59  <NA>     C   <NA>
    1212 15.73  <NA>     C   <NA>
    1213 15.43  <NA>     C   <NA>
    1214 19.10  <NA>     C   <NA>
    1215 16.84  <NA>     C   <NA>
    1216 19.50  <NA>     C   <NA>
    1217 22.02  <NA>     C   <NA>
    1218 23.01  <NA>     O   <NA>
    1219 26.17  <NA>     N   <NA>
    1220 25.69  <NA>     C   <NA>
    1221 27.15  <NA>     C   <NA>
    1222 26.71  <NA>     O   <NA>
    1223 30.50  <NA>     C   <NA>
    1224 36.31  <NA>     C   <NA>
    1225 37.41  <NA>     O   <NA>
    1226 40.37  <NA>     O   <NA>
    1227 27.15  <NA>     N   <NA>
    1228 29.27  <NA>     C   <NA>
    1229 28.21  <NA>     C   <NA>
    1230 29.57  <NA>     O   <NA>
    1231 35.98  <NA>     C   <NA>
    1232 45.09  <NA>     C   <NA>
    1233 50.40  <NA>     C   <NA>
    1234 49.27  <NA>     O   <NA>
    1235 55.52  <NA>     N   <NA>
    1236 26.57  <NA>     N   <NA>
    1237 24.82  <NA>     C   <NA>
    1238 25.46  <NA>     C   <NA>
    1239 25.12  <NA>     O   <NA>
    1240 24.26  <NA>     C   <NA>
    1241 23.74  <NA>     C   <NA>
    1242 21.48  <NA>     C   <NA>
    1243 23.09  <NA>     C   <NA>
    1244 25.95  <NA>     N   <NA>
    1245 29.84  <NA>     C   <NA>
    1246 30.58  <NA>     C   <NA>
    1247 31.36  <NA>     O   <NA>
    1248 32.16  <NA>     C   <NA>
    1249 33.80  <NA>     C   <NA>
    1250 36.66  <NA>     C   <NA>
    1251 38.83  <NA>     C   <NA>
    1252 31.59  <NA>     N   <NA>
    1253 31.52  <NA>     C   <NA>
    1254 33.29  <NA>     C   <NA>
    1255 33.56  <NA>     O   <NA>
    1256 29.35  <NA>     C   <NA>
    1257 30.67  <NA>     C   <NA>
    1258 34.19  <NA>     C   <NA>
    1259 34.50  <NA>     C   <NA>
    1260 35.37  <NA>     N   <NA>
    1261 34.26  <NA>     C   <NA>
    1262 33.11  <NA>     C   <NA>
    1263 32.35  <NA>     O   <NA>
    1264 37.05  <NA>     C   <NA>
    1265 48.15  <NA>     C   <NA>
    1266 54.00  <NA>     C   <NA>
    1267 54.40  <NA>     O   <NA>
    1268 59.11  <NA>     O   <NA>
    1269 34.12  <NA>     N   <NA>
    1270 36.05  <NA>     C   <NA>
    1271 37.57  <NA>     C   <NA>
    1272 36.67  <NA>     O   <NA>
    1273 33.70  <NA>     C   <NA>
    1274 32.55  <NA>     C   <NA>
    1275 33.60  <NA>     C   <NA>
    1276 35.00  <NA>     C   <NA>
    1277 41.68  <NA>     N   <NA>
    1278 44.87  <NA>     C   <NA>
    1279 44.44  <NA>     C   <NA>
    1280 42.31  <NA>     O   <NA>
    1281 48.45  <NA>     C   <NA>
    1282 54.18  <NA>     S   <NA>
    1283 43.24  <NA>     N   <NA>
    1284 41.28  <NA>     C   <NA>
    1285 41.98  <NA>     C   <NA>
    1286 46.85  <NA>     O   <NA>
    1287 40.00  <NA>     N   <NA>
    1288 36.63  <NA>     C   <NA>
    1289 34.62  <NA>     C   <NA>
    1290 32.72  <NA>     O   <NA>
    1291 39.25  <NA>     C   <NA>
    1292 44.38  <NA>     C   <NA>
    1293 46.88  <NA>     N   <NA>
    1294 45.80  <NA>     C   <NA>
    1295 50.18  <NA>     C   <NA>
    1296 47.74  <NA>     N   <NA>
    1297 33.27  <NA>     N   <NA>
    1298 31.17  <NA>     C   <NA>
    1299 29.41  <NA>     C   <NA>
    1300 29.56  <NA>     O   <NA>
    1301 32.71  <NA>     C   <NA>
    1302 34.39  <NA>     C   <NA>
    1303 44.19  <NA>     C   <NA>
    1304 48.61  <NA>     C   <NA>
    1305 54.35  <NA>     N   <NA>
    1306 24.69  <NA>     N   <NA>
    1307 23.65  <NA>     C   <NA>
    1308 24.51  <NA>     C   <NA>
    1309 24.89  <NA>     O   <NA>
    1310 16.65  <NA>     C   <NA>
    1311 25.50  <NA>     N   <NA>
    1312 24.06  <NA>     C   <NA>
    1313 22.15  <NA>     C   <NA>
    1314 23.87  <NA>     O   <NA>
    1315 24.38  <NA>     C   <NA>
    1316 26.24  <NA>     C   <NA>
    1317 23.69  <NA>     C   <NA>
    1318 23.26  <NA>     C   <NA>
    1319 19.12  <NA>     N   <NA>
    1320 19.31  <NA>     C   <NA>
    1321 19.37  <NA>     C   <NA>
    1322 21.87  <NA>     O   <NA>
    1323 19.61  <NA>     N   <NA>
    1324 19.07  <NA>     C   <NA>
    1325 22.17  <NA>     C   <NA>
    1326 19.62  <NA>     O   <NA>
    1327 19.04  <NA>     C   <NA>
    1328 24.28  <NA>     O   <NA>
    1329 17.43  <NA>     C   <NA>
    1330 18.83  <NA>     N   <NA>
    1331 20.46  <NA>     C   <NA>
    1332 22.85  <NA>     C   <NA>
    1333 20.15  <NA>     O   <NA>
    1334 22.69  <NA>     C   <NA>
    1335 19.10  <NA>     C   <NA>
    1336 18.68  <NA>     C   <NA>
    1337 21.32  <NA>     N   <NA>
    1338 19.82  <NA>     C   <NA>
    1339 21.74  <NA>     C   <NA>
    1340 25.18  <NA>     O   <NA>
    1341 16.45  <NA>     C   <NA>
    1342 15.15  <NA>     C   <NA>
    1343 14.83  <NA>     C   <NA>
    1344 20.59  <NA>     C   <NA>
    1345 23.03  <NA>     N   <NA>
    1346 17.63  <NA>     C   <NA>
    1347 17.51  <NA>     C   <NA>
    1348 16.69  <NA>     O   <NA>
    1349 18.16  <NA>     C   <NA>
    1350 20.02  <NA>     C   <NA>
    1351 16.66  <NA>     C   <NA>
    1352 16.55  <NA>     N   <NA>
    1353 20.06  <NA>     C   <NA>
    1354 21.31  <NA>     C   <NA>
    1355 24.13  <NA>     O   <NA>
    1356 19.47  <NA>     N   <NA>
    1357 21.78  <NA>     C   <NA>
    1358 21.89  <NA>     C   <NA>
    1359 22.16  <NA>     O   <NA>
    1360 17.28  <NA>     C   <NA>
    1361 20.54  <NA>     C   <NA>
    1362 20.64  <NA>     C   <NA>
    1363 24.38  <NA>     N   <NA>
    1364 22.15  <NA>     C   <NA>
    1365 20.84  <NA>     C   <NA>
    1366 20.99  <NA>     O   <NA>
    1367 21.20  <NA>     C   <NA>
    1368 22.05  <NA>     O   <NA>
    1369 13.16  <NA>     C   <NA>
    1370 24.57  <NA>     N   <NA>
    1371 24.68  <NA>     C   <NA>
    1372 26.51  <NA>     C   <NA>
    1373 32.43  <NA>     O   <NA>
    1374 22.89  <NA>     C   <NA>
    1375 24.08  <NA>     C   <NA>
    1376 21.46  <NA>     C   <NA>
    1377 23.91  <NA>     N   <NA>
    1378 24.09  <NA>     C   <NA>
    1379 21.54  <NA>     C   <NA>
    1380 26.47  <NA>     O   <NA>
    1381 26.47  <NA>     C   <NA>
    1382 26.39  <NA>     C   <NA>
    1383 24.12  <NA>     C   <NA>
    1384 16.61  <NA>     N   <NA>
    1385 16.90  <NA>     C   <NA>
    1386 19.55  <NA>     C   <NA>
    1387 18.36  <NA>     O   <NA>
    1388 16.55  <NA>     C   <NA>
    1389 22.00  <NA>     C   <NA>
    1390 23.24  <NA>     O   <NA>
    1391 22.85  <NA>     N   <NA>
    1392 17.31  <NA>     N   <NA>
    1393 13.51  <NA>     C   <NA>
    1394 19.87  <NA>     C   <NA>
    1395 20.54  <NA>     O   <NA>
    1396 10.82  <NA>     C   <NA>
    1397 13.42  <NA>     C   <NA>
    1398 11.48  <NA>     C   <NA>
    1399 20.08  <NA>     C   <NA>
    1400 19.47  <NA>     N   <NA>
    1401 19.69  <NA>     C   <NA>
    1402 22.78  <NA>     C   <NA>
    1403 24.74  <NA>     O   <NA>
    1404 18.51  <NA>     C   <NA>
    1405 19.44  <NA>     C   <NA>
    1406 14.76  <NA>     C   <NA>
    1407 22.90  <NA>     C   <NA>
    1408 19.72  <NA>     N   <NA>
    1409 14.65  <NA>     C   <NA>
    1410 13.18  <NA>     C   <NA>
    1411 13.07  <NA>     O   <NA>
    1412 11.71  <NA>     N   <NA>
    1413 14.23  <NA>     C   <NA>
    1414 19.73  <NA>     C   <NA>
    1415 17.28  <NA>     O   <NA>
    1416 15.02  <NA>     C   <NA>
    1417  8.76  <NA>     C   <NA>
    1418 16.02  <NA>     C   <NA>
    1419 21.61  <NA>     N   <NA>
    1420 21.15  <NA>     C   <NA>
    1421 18.57  <NA>     N   <NA>
    1422 26.79  <NA>     N   <NA>
    1423 18.56  <NA>     N   <NA>
    1424 17.41  <NA>     C   <NA>
    1425 19.59  <NA>     C   <NA>
    1426 20.40  <NA>     O   <NA>
    1427 19.90  <NA>     C   <NA>
    1428 15.28  <NA>     C   <NA>
    1429 22.61  <NA>     O   <NA>
    1430 16.65  <NA>     N   <NA>
    1431 20.98  <NA>     N   <NA>
    1432 21.49  <NA>     C   <NA>
    1433 20.99  <NA>     C   <NA>
    1434 21.25  <NA>     O   <NA>
    1435 20.00  <NA>     C   <NA>
    1436 17.87  <NA>     C   <NA>
    1437 20.88  <NA>     C   <NA>
    1438 17.17  <NA>     C   <NA>
    1439 20.01  <NA>     N   <NA>
    1440 22.17  <NA>     C   <NA>
    1441 19.95  <NA>     C   <NA>
    1442 21.08  <NA>     O   <NA>
    1443 17.71  <NA>     C   <NA>
    1444 16.96  <NA>     C   <NA>
    1445 11.88  <NA>     C   <NA>
    1446 18.70  <NA>     C   <NA>
    1447 23.70  <NA>     N   <NA>
    1448 23.96  <NA>     C   <NA>
    1449 24.81  <NA>     C   <NA>
    1450 26.72  <NA>     O   <NA>
    1451 21.26  <NA>     C   <NA>
    1452 25.49  <NA>     O   <NA>
    1453 17.64  <NA>     C   <NA>
    1454 25.95  <NA>     N   <NA>
    1455 27.31  <NA>     C   <NA>
    1456 27.97  <NA>     C   <NA>
    1457 30.31  <NA>     O   <NA>
    1458 24.09  <NA>     C   <NA>
    1459 25.97  <NA>     C   <NA>
    1460 27.56  <NA>     C   <NA>
    1461 26.23  <NA>     O   <NA>
    1462 27.50  <NA>     N   <NA>
    1463 28.73  <NA>     N   <NA>
    1464 28.11  <NA>     C   <NA>
    1465 30.61  <NA>     C   <NA>
    1466 32.40  <NA>     O   <NA>
    1467 23.33  <NA>     C   <NA>
    1468 25.74  <NA>     C   <NA>
    1469 21.29  <NA>     C   <NA>
    1470 27.65  <NA>     C   <NA>
    1471 30.93  <NA>     N   <NA>
    1472 31.94  <NA>     C   <NA>
    1473 30.20  <NA>     C   <NA>
    1474 31.24  <NA>     O   <NA>
    1475 27.82  <NA>     N   <NA>
    1476 26.76  <NA>     C   <NA>
    1477 27.20  <NA>     C   <NA>
    1478 27.85  <NA>     O   <NA>
    1479 27.09  <NA>     C   <NA>
    1480 36.34  <NA>     S   <NA>
    1481 26.04  <NA>     N   <NA>
    1482 25.60  <NA>     C   <NA>
    1483 24.26  <NA>     C   <NA>
    1484 23.28  <NA>     O   <NA>
    1485 27.82  <NA>     C   <NA>
    1486 31.12  <NA>     O   <NA>
    1487 24.73  <NA>     C   <NA>
    1488 24.54  <NA>     N   <NA>
    1489 27.92  <NA>     C   <NA>
    1490 30.25  <NA>     C   <NA>
    1491 29.52  <NA>     O   <NA>
    1492 27.18  <NA>     C   <NA>
    1493 22.23  <NA>     C   <NA>
    1494 23.52  <NA>     C   <NA>
    1495 20.87  <NA>     C   <NA>
    1496 32.20  <NA>     N   <NA>
    1497 34.86  <NA>     C   <NA>
    1498 34.26  <NA>     C   <NA>
    1499 36.18  <NA>     O   <NA>
    1500 31.24  <NA>     C   <NA>
    1501 33.93  <NA>     C   <NA>
    1502 37.15  <NA>     O   <NA>
    1503 33.06  <NA>     N   <NA>
    1504 36.76  <NA>     N   <NA>
    1505 36.49  <NA>     C   <NA>
    1506 35.50  <NA>     C   <NA>
    1507 37.49  <NA>     O   <NA>
    1508 34.88  <NA>     C   <NA>
    1509 36.75  <NA>     C   <NA>
    1510 37.04  <NA>     C   <NA>
    1511 38.13  <NA>     C   <NA>
    1512 37.02  <NA>     C   <NA>
    1513 37.11  <NA>     C   <NA>
    1514 36.24  <NA>     C   <NA>
    1515 28.25  <NA>     N   <NA>
    1516 30.30  <NA>     C   <NA>
    1517 27.27  <NA>     C   <NA>
    1518 28.85  <NA>     C   <NA>
    1519 29.59  <NA>     O   <NA>
    1520 22.29  <NA>     N   <NA>
    1521 23.47  <NA>     C   <NA>
    1522 27.66  <NA>     C   <NA>
    1523 21.71  <NA>     C   <NA>
    1524 22.75  <NA>     C   <NA>
    1525 28.91  <NA>     N   <NA>
    1526 26.24  <NA>     C   <NA>
    1527 27.47  <NA>     C   <NA>
    1528 20.86  <NA>     C   <NA>
    1529 21.68  <NA>     C   <NA>
    1530 15.87  <NA>     O   <NA>
    1531 21.49  <NA>     C   <NA>
    1532 26.89  <NA>     C   <NA>
    1533 28.67  <NA>     C   <NA>
    1534 26.89  <NA>     C   <NA>
    1535 29.22  <NA>     C   <NA>
    1536 29.22  <NA>     C   <NA>
    1537 30.97  <NA>     C   <NA>
    1538 29.25  <NA>     C   <NA>
    1539 29.96  <NA>     C   <NA>
    1540 29.35  <NA>     C   <NA>
    1541 32.66  <NA>     O   <NA>
    1542 31.19  <NA>     N   <NA>
    1543 29.22  <NA>     C   <NA>
    1544 28.82  <NA>     C   <NA>
    1545 28.32  <NA>     O   <NA>
    1546 32.05  <NA>     C   <NA>
    1547 31.29  <NA>     C   <NA>
    1548 32.00  <NA>     C   <NA>
    1549 28.00  <NA>     C   <NA>
    1550 29.01  <NA>     C   <NA>
    1551 27.70  <NA>     C   <NA>
    1552 31.86  <NA>     C   <NA>
    1553 36.25  <NA>     C   <NA>
    1554 42.75  <NA>     C   <NA>
    1555 47.41  <NA>     C   <NA>
    1556 51.38  <NA>     N   <NA>
    1557 50.60  <NA>     C   <NA>
    1558 49.34  <NA>     C   <NA>
    1559 44.71  <NA>     C   <NA>
    1560 63.07  <NA>     O   <NA>
    1561 63.34  <NA>     O   <NA>
    1562 66.96  <NA>     O   <NA>
    1563 36.09  <NA>     O   <NA>
    1564 64.67  <NA>     O   <NA>
    1565 21.55  <NA>     O   <NA>
    1566 26.65  <NA>     O   <NA>
    1567 60.45  <NA>     O   <NA>
    1568 25.82  <NA>     O   <NA>
    1569 32.52  <NA>     O   <NA>
    1570 41.02  <NA>     O   <NA>
    1571 41.93  <NA>     O   <NA>
    1572 27.94  <NA>     O   <NA>
    1573 51.87  <NA>     O   <NA>
    1574 66.74  <NA>     O   <NA>
    1575 65.58  <NA>     O   <NA>
    1576 67.74  <NA>     O   <NA>
    1577 43.98  <NA>     O   <NA>
    1578 37.23  <NA>     O   <NA>
    1579 69.15  <NA>     O   <NA>
    1580 70.78  <NA>     O   <NA>
    1581 21.93  <NA>     O   <NA>
    1582 46.57  <NA>     O   <NA>
    1583 63.81  <NA>     O   <NA>
    1584 47.08  <NA>     O   <NA>
    1585 63.52  <NA>     O   <NA>
    1586 31.73  <NA>     O   <NA>
    1587 49.24  <NA>     O   <NA>
    1588 65.44  <NA>     O   <NA>
    1589 75.86  <NA>     O   <NA>
    1590 67.42  <NA>     O   <NA>
    1591 57.13  <NA>     O   <NA>
    1592 60.42  <NA>     O   <NA>
    1593 75.52  <NA>     O   <NA>
    1594 38.21  <NA>     O   <NA>
    1595 50.02  <NA>     O   <NA>
    1596 53.78  <NA>     O   <NA>
    1597 61.00  <NA>     O   <NA>
    1598 73.78  <NA>     O   <NA>
    1599 61.39  <NA>     O   <NA>
    1600 49.60  <NA>     O   <NA>
    1601 71.88  <NA>     O   <NA>
    1602 66.74  <NA>     O   <NA>
    1603 70.97  <NA>     O   <NA>
    1604 63.94  <NA>     O   <NA>
    1605 73.81  <NA>     O   <NA>
    1606 42.37  <NA>     O   <NA>
    1607 51.24  <NA>     O   <NA>
    1608 18.18  <NA>     O   <NA>
    1609 53.13  <NA>     O   <NA>
    1610 47.68  <NA>     O   <NA>
    1611 65.44  <NA>     O   <NA>
    1612 38.53  <NA>     O   <NA>
    1613 32.25  <NA>     O   <NA>
    1614 61.86  <NA>     O   <NA>
    1615 22.69  <NA>     O   <NA>
    1616 59.93  <NA>     O   <NA>
    1617 33.99  <NA>     O   <NA>
    1618 79.22  <NA>     O   <NA>
    1619 31.58  <NA>     O   <NA>
    1620 47.41  <NA>     O   <NA>
    1621 46.59  <NA>     O   <NA>
    1622 48.25  <NA>     O   <NA>
    1623 48.73  <NA>     O   <NA>
    1624 54.68  <NA>     O   <NA>
    1625 37.86  <NA>     O   <NA>
    1626 68.44  <NA>     O   <NA>
    1627 42.81  <NA>     O   <NA>
    1628 60.62  <NA>     O   <NA>
    1629 61.36  <NA>     O   <NA>
    1630 35.03  <NA>     O   <NA>
    1631 62.75  <NA>     O   <NA>
    1632 71.64  <NA>     O   <NA>
    1633 57.53  <NA>     O   <NA>
    1634 50.97  <NA>     O   <NA>
    1635 73.30  <NA>     O   <NA>
    1636 62.30  <NA>     O   <NA>
    1637 65.69  <NA>     O   <NA>
    1638 61.76  <NA>     O   <NA>
    1639 67.21  <NA>     O   <NA>
    1640 61.89  <NA>     O   <NA>
    1641 74.72  <NA>     O   <NA>
    1642 48.75  <NA>     O   <NA>
    1643 60.17  <NA>     O   <NA>
    1644 43.92  <NA>     O   <NA>
    1645 70.16  <NA>     O   <NA>
    1646 22.10  <NA>     O   <NA>
    1647 27.84  <NA>     O   <NA>
    1648 65.78  <NA>     O   <NA>
    1649 67.04  <NA>     O   <NA>
    1650 53.99  <NA>     O   <NA>
    1651 54.21  <NA>     O   <NA>
    1652 62.03  <NA>     O   <NA>
    1653 63.64  <NA>     O   <NA>
    1654 42.47  <NA>     O   <NA>
    1655 65.50  <NA>     O   <NA>
    1656 65.50  <NA>     O   <NA>
    1657 73.55  <NA>     O   <NA>
    1658 63.48  <NA>     O   <NA>
    1659 52.97  <NA>     O   <NA>
    1660 72.75  <NA>     O   <NA>
    1661 75.75  <NA>     O   <NA>
    1662 38.25  <NA>     O   <NA>
    1663 68.43  <NA>     O   <NA>
    1664 54.20  <NA>     O   <NA>
    1665 63.96  <NA>     O   <NA>
    1666 23.98  <NA>     O   <NA>
    1667 52.93  <NA>     O   <NA>
    1668 58.06  <NA>     O   <NA>
    1669 64.79  <NA>     O   <NA>
    1670 55.54  <NA>     O   <NA>
    1671 61.69  <NA>     O   <NA>
    1672 69.12  <NA>     O   <NA>
    1673 78.93  <NA>     O   <NA>
    1674 71.37  <NA>     O   <NA>
    1675 78.14  <NA>     O   <NA>
    1676 54.05  <NA>     O   <NA>
    1677 72.78  <NA>     O   <NA>
    1678 58.40  <NA>     O   <NA>
    1679 58.78  <NA>     O   <NA>
    1680 68.40  <NA>     O   <NA>
    1681 64.90  <NA>     O   <NA>
    1682 67.95  <NA>     O   <NA>
    1683 53.68  <NA>     O   <NA>
    1684 49.41  <NA>     O   <NA>
    1685 64.49  <NA>     O   <NA>
    1686 54.09  <NA>     O   <NA>

``` r
attributes(pdb)
```

    $names
    [1] "atom"   "xyz"    "seqres" "helix"  "sheet"  "calpha" "remark" "call"  

    $class
    [1] "pdb" "sse"

There are lots of functions that can work with these `pdb` objects.

``` r
head( pdbseq(pdb))
```

      1   2   3   4   5   6 
    "P" "Q" "I" "T" "L" "W" 

We can have a quick interactive view of any of these `pdb` objects:

``` r
#install.packages("pak")

#pak::pkg_install(c("bioboot/bio3dview",
                   #"NGLVieweR", 
                   #"bioc::msa")
                 #)
```

``` r
library(bio3dview)

#view.pdb(pdb)
```

## Predict the flexiblity of a given structure

Let’s do a Normal Mode Analysis (NMA) to predict the felxibility of a
give `pdb` object:

``` r
adk <- read.pdb("6s36")
```

      Note: Accessing on-line PDB file
       PDB has ALT records, taking A only, rm.alt=TRUE

``` r
m <- nma(adk)
```

     Building Hessian...        Done in 0.04 seconds.
     Diagonalizing Hessian...   Done in 0.35 seconds.

``` r
plot(m)
```

![](class10_files/figure-commonmark/unnamed-chunk-10-1.png)

View the results with an interactive structure view

``` r
#view.nma(m)
```

``` r
mktrj(m, file="nma.pdb")
```

> Q. Create a custom view highlighting the active site ASP (`resno=25`),
> the two chains (in your choice of colors) and the ligand all on a
> custom color background.

``` r
#library(NGLVieweR)
#active.site <- atom.select(pdb, resno=25)

#view.pdb(pdb,
         #cols = c("red", "blue"),
         #highlight = active.site,
         #highlight.style = "spacefill",
         #backgroundColor = "pink") |>
  #setRock()
```

## Comparative analysis of the ADK family

Our first step is find a sequence for this family. We will use the
database ID “1ake_A” here:

``` r
id <- "1ake_A"

aa <- get.seq(id)
```

    Warning in get.seq(id): Removing existing file: seqs.fasta

    Fetching... Please wait. Done.

``` r
aa
```

                 1        .         .         .         .         .         60 
    pdb|1AKE|A   MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVT
                 1        .         .         .         .         .         60 

                61        .         .         .         .         .         120 
    pdb|1AKE|A   DELVIALVKERIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFDVPDELIVDRI
                61        .         .         .         .         .         120 

               121        .         .         .         .         .         180 
    pdb|1AKE|A   VGRRVHAPSGRVYHVKFNPPKVEGKDDVTGEELTTRKDDQEETVRKRLVEYHQMTAPLIG
               121        .         .         .         .         .         180 

               181        .         .         .   214 
    pdb|1AKE|A   YYSKEAEAGNTKYAKVDGTKPVAEVRADLEKILG
               181        .         .         .   214 

    Call:
      read.fasta(file = outfile)

    Class:
      fasta

    Alignment dimensions:
      1 sequence rows; 214 position columns (214 non-gap, 0 gap) 

    + attr: id, ali, call

Search for related sequences in the database

``` r
#blast <- blast.pdb(aa)
```

``` r
#head(blast$hit.tbl)
```

``` r
#hits <- plot(blast) 
```

``` r
#hits$pdb.id
```

``` r
#f <- get.pdb(hits$pdb.id, path="pdbs", split = T, gzip = T)
```

``` r
#pdbs <- pdbaln(f, fit= TRUE, exefile = "msa")
#pdbs
```

``` r
#view.pdbs(pdbs)
```

PCA of all this structural data (x, y, and z atom coordinates):

``` r
#pc <- pca(pdbs)
#plot(pc) 
```

``` r
#plot(pc, 1:2)
```

Interactive view of the PC1 captured structural

``` r
#view.pca(pc)
```

``` r
#mktrj(pc, file="pca.pdb")
```
