### FilterDuplicates

```toml
[[steps]]
    action = "FilterDuplicates"
    false_positive_rate = float#
            # the false positive rate of the filter.
            # 0..1
    seed = 59 # required!
    target = "All"|"Read1"|"Read2"|"Index1"|"Index2"
    invert = false # bool, if true, keep only duplicates
```

Remove duplicates from the stream using a [Cuckoo filter](https://en.wikipedia.org/wiki/Cuckoo_filter).

That's a probabilistic data structure, accordingly there's a false positive rate,
and a tunable memory requirement.

Needs a seed for the random number generator, and a target
to know which reads to consider for deduplication (filters the complete molecule, like
all other filters of course).

The lower you set the false positive rate, the higher your memory requirements will be.
0.00001 might be a good place to start. 

If you set the false positive rate to 0.0, a HashSet will be used instead,
which will produce exact results, albeit at the expense of keeping a copy of *all* reads in memory! 

Note that chaining these probably is not what you want (the second filter wouldn't see all fragments!),
therefore we have a special 'All' target here, which will only filter molecules where all segments are duplicated.
