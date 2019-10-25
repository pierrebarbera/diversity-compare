Folder structure:
- `datasets`: contains the different data sets used in the evaluations. Should be treated as read-only.
- `metrics`: contains one executable script per tested metric, producing a standardized diversity metric output
- `programs`: contains the programs/dependencies needed to run all scripts in `metrics`. Should aim to be "drop in" directory.
- `workdir`: will contain intermediate data produced during computation of diversity metrics from the raw data sets.
- `viz`: to contain visualizations produced along the way