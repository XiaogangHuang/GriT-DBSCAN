# GriT-DBSCAN

The C++ code for "[GriT-DBSCAN: A spatial clustering algorithm for very large databases](https://arxiv.org/pdf/2210.07580.pdf)"

Let $\kappa$ be the maximum number of iterations in the merging step.
Then, GriT-DBSCAN runs in $O(d(2\lceil \sqrt{d}\rceil + 1)^{d}\kappa n + d\eta)$ expected time, regardless of the value of $MinPts$.

<!-- In the paper, there is a small bug in the description of "Constructing the Grids" and a small bug in the description of "Non-empty Neighboring Grids Query": -->
<!-- - For any two grids $g_i, g_j \in G_s$ with $i < j$, there exists an integer $z\in [1, d]$ such that $g_{iz} < g_{jz}$ and $g_{iw} = g_{jw}$, for each $w\leq z-1$. (page 4) -->
<!-- - $t^\prime.offset = t.offset+\max\{|KEY(t^\prime) - g_{i2}|-1,0\}^2$. (Eq. (2)) -->


## DATA FORMAT
The input dataset needs to be preprocessed into a text file with the following format:

* In each line: its coordinates and the point's id, where the id is an integer in the range [0,  n-1].

* In total, there are n lines where the numbers at each line are separated by a space. 

For example, a 2-dimensional data set containing 4 points:
 
9.3498     56.7408     17.0527     0

9.3501     56.7406     17.6148     1

9.3505     56.7405     18.0835     2

9.3508     56.7404     18.2794     3
