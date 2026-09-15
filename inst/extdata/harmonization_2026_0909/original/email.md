# emails

-2026-09-09

```
For jemeF and hiC_PO, as explained earlier, "id" is the label used in the raw datasets, and "symbol" is the updated/final label that should be used instead; you can get a sense of the scale of the issue by counting how many discordances of "id" vs "symbol" ~ around n=1400-1500; note, these two files are simply KEYS, not data

For hiC_PP, I retained the actual raw data columns since it's just tissue,Promoter(gene1),B(gene2), and added Pro_sym and B_sym; you can see it's on the order 900 genes or so that had discrepancies;

Note, for all of these files I expanded/reshaped to long, when there were multiple genes separated by ";" or "|" or whatever character they used (e.g for JEME, when a given enhancer was predicted to target more than one gene

```
