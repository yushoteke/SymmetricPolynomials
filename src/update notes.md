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

In commit "polynomial use two containers",
I observed that a lot of time was spent on push! new items into polynomial.
Original implementation uses a sorted dict, which is internally represented as a tree, so
every push has O(log(n)) time.
Since sometimes inside the push, we don't add new values, but only modify the values, instead keep
two copies inside the polynomial, a dict for monomial->coefficient pairs, and a tree for sorting
monomials.
This way, when the update doesn't involve adding values, the complexity is O(1).
The runtime reduced from 24.3s to 21s, which is another 12% reduction. Unfortunately, because
an extra container is involved, the memory allocated increated from 3.23M to 17.58M, and memory usage
increase from 334MiB to 2.45GiB

In commit "use no sort containers",
Since most of the time is spent sorting, elliminating that improves most performance.
Notice the order in which we pick the polynomial term is fixed, so I wrote a seperate iteration to
do that. The polynomial containers are now only responsible for holding coefficients and grouping
like terms.
Now, the polynomial structure is holding two dictionaries instead of a dict and a tree.
Also, I noticed there's some redundant information when one of the dicts hold elementary terms.
By elliminating that, no need for the vector to resize anymore.
Also, remove redundant keys from the dictionaries while iterating reduced the likelihood for a resize
to happen.
These improvements collectively reduced runtime from 21s to 18.5s. Memory allocated increated from
17.58M to 38.22M, memory usage increased from 2.45GiB to 5.57GiB

In commit "reduce usage of sets"
An inquiry into memory usage suggests that most of allocated memory is spent on creating new sets
and pushing terms into dictionary. The second most allocated is collecting coefficients into a container.
These improvements collectively reduced fruntime from 18.5s to 16.5s.
Memory allocated decreased from 38.22M to 5.42M, and memory usage decreased from 5.57GiB to 2.623GiB

In commit "further investigations"
A further analysis shows that in decompose12, the creation of new_term is type unstable. Used some
gimmick to make it type stable, but did not result in performance gains.
Memory analysis shows that most of the time is spent to check if key is in dictionary, get value
corresponding to key, and modify value corresponding to key. If I cannot find anymore mathematical
improvements, then the next step is to find a better dictionary.

Some other thoughts:
Since all tuples are fixed length, maybe think in trie direction?
Since we sometimes group all terms which looks like S(a,x), where a is fixed, together, maybe
think of this as a prefix?
Inside push!, sometimes there's double lookup. Maybe tokens could fix this?

In commit "array pool"
By doing memory analysis, most of the memories allocations occur in push!. The suspected reasons
are, when the arrays grow out of size, it needs to be dynamically resized, which results in more
allocation. Since there are many allocation of small and medium sized arrays, the solution to this
is to use an array pool. When the polynomial needs an array, it "borrows" an array from the pool.
When the pool is out of array, it creates a new one. When the polynomial finishes using the array,
it "returns" the array to the pool.
Made two changes. First minor change in utility.jl file avoids allocating extra arrays.
Then, the basic array pool implementation reduced memory allocation from 5.42M to 3.9M and memory
usage from 2.623GiB to 1.91GiB
