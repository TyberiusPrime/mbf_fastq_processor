---
weight: 5
---
# Report

Capture data for the final report (see [..](output section)).

You can add multiple reports, at any stage of your transformation chain
to get e.g. before/after filtering reports.

```toml
[[step]]
    action = 'Report'
    label = "report" # Key that the report will be listed under. Must be distinct
    count = true # count reads at this position
    base_statistics = false # include base distribution at each read position, q20, q30, total, gc bases
    length_distribution = false # capture read length distribution
    duplicate_count = false # count duplicates using Cukoo filter
```

Statistics available (for each 'segment'):

- read counts
- total base count
- bases count q20 or better
- bases count q30 or better
- read length distribution
- AGTCN counts at each position
- expected error rate at each position
- duplicate count (if each read occurs twice, duplicate count = read count / 2)

