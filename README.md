# partial_order_alignment

Author: Bansho Masutani

Email: ban-m@g.ecc.u-tokyo.ac.jp

# Short introduction

This repositry implements the Lee(2004)'s partial order alignment. Also, by adding `poa_hmm = {features=["poa_simd"]}` to your `Cargo.toml`, the POA would use SIMD instructions in some part of the alignment.

Compared to the default mode, it runs 2-3 times faster. Note that the SIMD version is entirely optional, and you should opt-in these feature and use the nightly compiler to activate this feature.

Here are some recommendations and remarks;

- The alignment is a semi-global one. In other words, it allows to leave some leading/trailing node that does not align to any query sequence.
- The POA graph is not tuned to generate consensus currently. I recommend to benchmark by some synthetic datasets.
- Unless noted, the alignment is the optimal one. It allocates and fills all the cells of the DP matrix. If it is so slow, please use `banded` functions with an appropriate diameter parameter. If you do not use any noisy dataset such as PacBio's or ONT's long reads, it is the case.
- There is no way to obtain 2nd or 3rd optimal alignment. 
- This crate not only provides you to construct a POA graph but also enables you to calculate the likelihood of a sequence from the given POA graph. See `./src/forward.rs` for more details. Currently, I use this functionality for clustering long-reads into each haplotype.

# License

MIT license
