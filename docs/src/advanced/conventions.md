# Output and Interface Conventions

Every analysis in NMRAnalysis.jl writes the same set of files, with the same column rules,
whatever kind of experiment it was. This page is the specification those files are written
against. It is aimed at anyone adding or changing an analysis; the per-experiment pages
describe what each one means.

The point is that a results folder should be readable without knowing which routine
produced it, and that a script written against one analysis should work against another.

## The output folder

```
out/
  summary.txt          human-readable record of the whole analysis
  results.csv          fitted and derived parameters, one row per entity
  series.csv           the measurements, one row per entity per coordinate point
  global.csv           parameters fitted once for the whole analysis (when there are any)
  <overview>.pdf       fit.pdf in 1D, summary.pdf in 2D
  regions/             1D: one file pair per region
    signal.csv
    signal.pdf
  peaks/               2D: one file pair per peak
    L23N.csv
    L23N.pdf
  cluster_*.pdf        2D only: one plot per group of overlapping peaks
```

A per-entity CSV and its plot share a basename and differ only in extension, so the data
behind any plot is beside it. The per-entity CSVs are the rows of `series.csv` filtered to
that entity; the duplication is deliberate, since opening one peak's data should not
require filtering a file of several thousand rows.

`global.csv` is written only when the analysis fits something globally: a titration `Kd`,
an exchange `kex`. A plain relaxation fit has nothing global and the file is absent.

Cluster plots have no CSV of their own: a cluster is a group of overlapping peaks rather
than an entity with its own parameters, so it stays at the top level under its existing
`cluster_LABEL.pdf` name.

## Entities, coordinates and scope

An **entity** is whatever the fitted parameters are indexed by: a peak in 2D, a region in
1D, an experiment in a joint fit. Every entity has a `label`.

A **coordinate** is anything that distinguishes one measurement of an entity from another:
a relaxation delay, a gradient strength, a spin-lock field, a concentration, or a
categorical tag such as TROSY versus anti-TROSY. Coordinates may be numeric or
categorical, and there may be several.

A parameter has a **scope**: it belongs to a single series (one entity at one setting of
the grouping coordinates), to an entity as a whole, or to the analysis as a whole. The
first two live in `results.csv` and the third in `global.csv`.

## Column rules

These hold in every CSV.

Column names are **ASCII**, so that `df.tau_c` works: `eta_xy`, `tau_c`, `pulse90`,
`R1rho`. The typeset names belong in the GUI and in `summary.txt`, not in a file header.

A column carrying a physical quantity names its **unit in parentheses** after a space, in
ASCII: `R (s-1)`, `tau_c (ns)`, `D (1e-10 m2/s)`, `pulse90 (us)`. A dimensionless quantity
has no parentheses.

A value column `X` may be accompanied by an **uncertainty** column `X_err` and a **fitted
value** column `X_fit`, each repeating the unit: `R (s-1)`, `R_err (s-1)`. Repeating it
keeps the two symmetrical for anything reading them.

Every other column is a **key**. So the rule for reading one of these files generically is:
strip the parenthesised unit from each header, then any column whose name is `X`, `X_err`
or `X_fit` for some `X` is a value, and everything else identifies the row.

A **blank key** means the row applies to every value of that key. TRACT's `tau_c` describes
a region rather than either component of its TROSY/anti-TROSY pair, so it is written with
`label` filled in and the `which` key blank.

Numbers are written at **full precision**. These are machine files and rounding is
irreversible; `summary.txt` is where numbers are rounded for reading. A value that does not
exist for a row is written `NA`, which is distinct from a blank key.

## `series.csv`

One row per entity per coordinate point, in long form, because an analysis may have
several coordinates and several datasets.

```
source,label,which,time (s),I,I_err,I_fit
12/pdata/1,amide,trosy,0.000,9421.3,12.4,9430.1
12/pdata/1,amide,trosy,0.005,8684.0,12.4,8687.5
13/pdata/1,amide,anti,0.000,9388.2,12.4,9401.7
```

`source` names the dataset each row came from and is always present, even when it is
constant. It is provenance, not meaning: two files may be replicates with identical
coordinates, so whatever physical variable distinguishes datasets (a concentration, a
spin-lock field, a TROSY tag) gets its own coordinate column as well.

`I_fit` is the model evaluated at the measured coordinates, not on a fine grid, so that
residuals are a subtraction. It is `NA` where nothing was fitted. A smooth curve for a
figure is recoverable from the parameters.

## `results.csv`

One row per entity, wide, because this is the table you sort by residue number and plot
against it.

```
label,which,A,A_err,R (s-1),R_err (s-1),eta_xy (s-1),eta_xy_err (s-1),tau_c (ns),tau_c_err (ns)
amide,trosy,9421.3,12.4,16.02,0.41,NA,NA,NA,NA
amide,anti,9388.2,12.4,54.10,0.93,NA,NA,NA,NA
amide,,NA,NA,NA,NA,19.04,0.51,14.86,0.42
```

The key columns are `label` plus the experiment's grouping coordinates. 2D adds `resnum`,
`resname` and `atom`, derived from the label. 1D adds the region bounds `lo (ppm)` and
`hi (ppm)`, which is also how a saved region list is restored.

## `global.csv`

Long, because these parameter sets are small and structurally unlike each other.

```
parameter,value,error,unit
Kd,12.4,0.8,uM
```

## `summary.txt`

The human-readable record, and the one file where numbers are rounded. It carries the
package version and the date, the input filenames and titles, the sample information, the
region or peak definitions with the noise position and integration width, the acquisition
parameters actually used, and the key results formatted with units.

It should also carry the Julia call that would repeat the analysis, together with where
each resolved parameter came from, which makes the annotation lookup auditable:

```
Reproduce:
    relaxation1d("11";
                 tau=[0.01, 0.03, 0.06, 0.10, 0.20, 0.40],
                 model=:exponential,
                 integration=(peakppm=8.21, noiseppm=-1.00, ppmwidth=0.60))

tau from vdlist; model from annotation relaxation.model; region selected interactively.
```

## Entry points

The data is positional; everything else is a keyword. A coordinate list is a keyword even
where it is required, because it may instead be resolved from the experiment, and an
optional positional argument that is sometimes inferred reads badly.

Keyword names are chosen to be informative rather than to match annotation keys, so
`relaxationtimes`, `Trelax` and `Tsat` keep the names a spectroscopist would use.

Where a parameter can be resolved rather than typed, it is looked for in a fixed order:
an explicit argument, then a pulse-sequence annotation, then sample metadata, then a Bruker
acquisition parameter, and finally a question asked before any window opens. `prompt=false`
(the default outside an interactive session) turns that last step into either a documented
default or an error naming the argument to pass.

Not every analysis resolves everything. An experiment combining many datasets, each with
its own lists of offsets, powers and delays, is not something a user can reasonably type,
so `exchange1d` reads those from annotations and does not offer a fallback. Its questions
are about the fit (which model, which molecules, which parameters to fix) rather than about
what the experiment was.

## Implementation status

| Module | `summary.txt` | `results.csv` | `series.csv` | `global.csv` | per-entity files |
|---|---|---|---|---|---|
| Analysis1D | yes | yes | yes | yes | yes |
| GUI2D | no | old format | no | no | plots only |
| R1rho | old format | no | no | no | no |
| Exchange1D | old format | no | no | no | no |

The `Reproduce:` line records the resolved arguments, which each entry point hands to the
writer. It does not yet say *where* each one came from (an annotation, the `vdlist`, a
question), which needs the resolution chain to report its own winner rather than just its
result.
