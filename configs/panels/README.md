# Gene panels

`hvg2000_v1.txt` is the **frozen** 2000-gene highly-variable-gene panel that
defines the model input `X`. It is produced **once** by
`cellfate.data.normalize.fit_gene_panel` on a reference atlas, then committed
here and never edited. Its content hash (`GenePanel.hash()`) is written into
every shard's scalers and into the deployment bundle; a mismatch is a hard error.

Format: one gene symbol per line; lines starting with `#` are ignored.

`example_panel.txt` is a tiny placeholder used by the unit tests only.
