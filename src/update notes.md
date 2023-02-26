The recent commit is mostly tuning

The biggest improvement is combining like terms, which cut run time by 30%.
This reduced time from 36.9s to 26.5s
Also, memory usage reduced from 72.3M allocations to 10.5M, and memory used reduced
from 7.2GiB to 2.5GiB

Then, there are a variety of small improvements, mostly related to memory.
Most notably, created a new data structure that creates arrays for repeated use, which saved some time
that's related to allocation.
This saved about 8%, reducing time from 26.5 s to 24.3s.
Memory Usage reduced from 10.5M to 3.23M, and memory used reduced from 2.5GiB to 334.3MiB
The downside is, code readability becomes horrible.

This commit will be called "Improvements and Benchmarks"
Will do a systematic cleaning in a later commit.
