# Alexs-Cellular-Burden-Model
## Preface
This model is aimed to capture the burden effects due to the shared resources within a cell, and is derived from the work in [this nature article](https://www.nature.com/articles/s41467-020-18392-x) by Velia Siciliano.

The original model derived an effective rate constant for a coupled series of reactions. For example, these reactions:
* A1 + R <=> C1 -> B1 + R
* A2 + R <=> C1 -> B2 + R

Can be approximated with:
* A1 -> B1
* A2 -> B2

By calculating that the new reaction rates as follows:
