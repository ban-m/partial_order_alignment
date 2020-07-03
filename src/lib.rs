#[allow(unused_imports)]
#[macro_use]
extern crate log;
#[cfg(feature = "poa_simd")]
extern crate packed_simd;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
mod config;
pub use config::*;
pub mod base;
mod base_table;
use base::Base;
mod add_sequence;
mod formatters;
pub mod forward;
pub mod gen_sample;
mod remove_nodes;
const SMALL: f64 = 0.000_000_001;
const LAMBDA_INS: f64 = 0.09;
const LAMBDA_MATCH: f64 = 0.09;
const THR: f64 = 0.4;
pub const DEFAULT_LK: f64 = -150.;
pub mod generate;
#[cfg(test)]
extern crate rayon;
#[cfg(test)]
mod tests;

/// Edit operations.
/// The inner fields in the `EditOp::Match` and `EditOp::Deltion` enums
/// are to indicate the next position of the POA graph.
/// Note that, when we encounter an `EditOp::Insertion` operation,
/// the position of the graph would never change. Thus,
/// there is no inner fileds in an `EditOp::Insertion`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EditOp {
    /// Match between the query and the graph.
    Match(usize),
    /// Deletion in the graph. In other words, it proceeds the position of the
    /// graph.
    Deletion(usize),
    /// Insertion to the graph. In other words, it proceeds the position of the query.
    Insertion,
    /// Stop operation. It is only used as a sentinel.
    Stop,
}

impl std::fmt::Display for EditOp {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match &self {
            EditOp::Match(_) => write!(f, "M"),
            EditOp::Deletion(_) => write!(f, "D"),
            EditOp::Insertion => write!(f, "I"),
            EditOp::Stop => write!(f, "S"),
        }
    }
}

pub type POA = PartialOrderAlignment;
/// A struct to represent an Partial Order Alignment.
/// To see the examples, please see `./tests.ts`.
#[derive(Clone, Default)]
pub struct PartialOrderAlignment {
    nodes: Vec<Base>,
    weight: f64,
    dfs_flag: Vec<usize>,
    dfs_stack: Vec<usize>,
    mapping: Vec<usize>,
}

type TraceBack = Vec<EditOp>;

impl PartialOrderAlignment {
    /// Return the nodes of the graph.
    /// Note that there's no way to modify each node directly.
    pub fn nodes(&self) -> &[Base] {
        &self.nodes
    }
    /// Return the number of the nodes. Same as `self.nodes().len()`.
    pub fn num_nodes(&self) -> usize {
        self.nodes.len()
    }
    /// Return the number of edges.
    pub fn num_edges(&self) -> usize {
        self.nodes
            .iter()
            .map(|n| n.edges.len() + n.ties.len())
            .sum::<usize>()
    }
    /// Return the weight of the graph.
    /// Usually, it returns the number of sequences added to the graph.
    pub fn weight(&self) -> f64 {
        self.weight
    }
    /// View alignment.
    /// # Example
    /// ```
    /// use poa_hmm::POA;
    /// let query = b"CAGTGTACGTCA";
    /// let poa = POA::default().add_default(query, 1.);
    /// let (_,aln) = poa.align(query, (-1,-1,|x,y|if x == y { 1 } else { -1 }));
    /// let (qry ,rfr) = poa.view(query, &aln);
    /// eprintln!("{}\t{}", qry, rfr);
    /// ```
    pub fn view(&self, seq: &[u8], traceback: &[EditOp]) -> (String, String) {
        let mut q_pos = 0;
        let (mut q, mut g) = (String::new(), String::new());
        for &op in traceback {
            match op {
                EditOp::Deletion(g_pos) => {
                    q.push('-');
                    g.push(self.nodes()[g_pos].base() as char);
                }
                EditOp::Insertion => {
                    g.push('-');
                    q.push(seq[q_pos] as char);
                    q_pos += 1;
                }
                EditOp::Match(g_pos) => {
                    g.push(self.nodes()[g_pos].base() as char);
                    q.push(seq[q_pos] as char);
                    q_pos += 1;
                }
                EditOp::Stop => {}
            }
        }
        (q, g)
    }
    /// Construct a POA graph from a given slice of vectors. Note that
    /// `ws` and `seqs` should have the same length.
    pub fn from_vec<F>(seqs: &[Vec<u8>], ws: &[f64], parameters: (i32, i32, &F)) -> POA
    where
        F: Fn(u8, u8) -> i32,
    {
        let seqs: Vec<_> = seqs.iter().map(|e| e.as_slice()).collect();
        Self::from_slice(&seqs, ws, parameters)
    }
    /// Construct a POA graph from a given slice of slices.
    /// Note that `seqs` and `ws` should have the same length.
    pub fn from_slice<F>(seqs: &[&[u8]], ws: &[f64], parameters: (i32, i32, &F)) -> POA
    where
        F: Fn(u8, u8) -> i32,
    {
        assert_eq!(seqs.len(), ws.len());
        Self::default().update(seqs, ws, parameters)
    }
    /// Construct a POA graph from a given slice of slices.
    /// In this function, the alignment would be computed in a 'banded'-manner.
    /// The diameter is `d` and the computational time is approximately O(Lnd),
    /// where L is the maximum length of seqs and n is the number of sequences.
    /// The weight of each sequences are fixed to 1.
    pub fn from_slice_banded<F>(seqs: &[&[u8]], ps: (i32, i32, &F), d: usize) -> POA
    where
        F: Fn(u8, u8) -> i32,
    {
        // Determine a seed.
        let seed = seqs
            .iter()
            .enumerate()
            .map(|(idx, s)| s.iter().map(|&x| x as u64).sum::<u64>() * idx as u64)
            .map(|x: u64| x % 179132109)
            .sum::<u64>();
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
        if seqs.is_empty() {
            return Self::default();
        }
        let max_len = seqs.iter().map(|s| s.len()).max().unwrap_or(0);
        let node_num_thr = (max_len as f64 * 1.05).floor() as usize;
        rand::seq::index::sample(&mut rng, seqs.len(), seqs.len())
            .into_iter()
            .map(|idx| &seqs[idx])
            .fold(Self::default(), |x, y| {
                if x.nodes.len() > node_num_thr {
                    x.add_banded(y, ps, d).remove_node(0.3)
                } else {
                    x.add_banded(y, ps, d)
                }
            })
            .remove_node(0.3)
            .finalize()
    }
    /// Construct a POA graph from a given slice of slices.
    /// Panics if `seqs.len() != ws.len()`.
    /// Note that it doew not remove any nodes after alignment,
    /// thus potentially it is very slow.
    pub fn from_slice_no_trim<F>(seqs: &[&[u8]], ws: &[f64], param: (i32, i32, &F)) -> POA
    where
        F: Fn(u8, u8) -> i32,
    {
        use rand::seq::SliceRandom;
        let seed = 99_111 * ((ws.iter().sum::<f64>().floor()) as u64);
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
        seqs.choose_multiple(&mut rng, seqs.len())
            .fold(POA::default(), |x, y| x.add(y, 1., param))
            .finalize()
    }
    /// Construct a POA graph from a given slice of slices,
    /// with default parameters and uniformal weights.
    pub fn from_slice_default(seqs: &[&[u8]]) -> POA {
        let ws = vec![1.; seqs.len()];
        fn score(x: u8, y: u8) -> i32 {
            if x == y {
                1
            } else {
                -1
            }
        };
        POA::from_slice(seqs, &ws, (-1, -1, &score))
    }
    /// Construct a POA graph from a given slice of vectors,
    /// with default parameters and uniformal weights.
    pub fn from_vec_default(seqs: &[Vec<u8>]) -> POA {
        let ws = vec![1.; seqs.len()];
        fn score(x: u8, y: u8) -> i32 {
            if x == y {
                1
            } else {
                -1
            }
        };
        let seqs: Vec<_> = seqs.iter().map(|e| e.as_slice()).collect();
        POA::from_slice(&seqs, &ws, (-1, -1, &score))
    }
    /// Update `self` by aligning all `seqs`.
    /// Panics if `seqs.len() != ws.len()`.
    /// Note that it automatically "prunes" the graph.
    /// In other words, if the number of the nodes becames too many,
    /// it removes weak nodes and edges sufficiently.
    /// If you do not want to this behavior, use `update_thr` instead.
    pub fn update<F>(self, seqs: &[&[u8]], ws: &[f64], params: (i32, i32, &F)) -> POA
    where
        F: Fn(u8, u8) -> i32,
    {
        self.update_thr(seqs, ws, params, THR, 1.5)
    }
    /// Update `self` by aligning all `seqs`.
    /// Panics if `seqs.len() != ws.len()` or `thr > 1.0` or `coef < 1.0`
    /// Note that it automatically "prunes" the graph.
    /// In other words, if the number of the nodes becames more than `coef*max`,
    /// where max is the maximum length of seqs,
    /// it removes weak nodes and edges sufficiently.
    /// The `thr` parameters tunes how hard we prune edges and nodes.
    /// If you want to prune edges very hard, please set coef = 0.9.
    pub fn update_thr<F: Fn(u8, u8) -> i32>(
        self,
        seqs: &[&[u8]],
        ws: &[f64],
        params: (i32, i32, &F),
        thr: f64,
        coef: f64,
    ) -> POA {
        // Determine a seed.
        let seed = ws
            .iter()
            .enumerate()
            .map(|(idx, w)| idx as f64 * w)
            .sum::<f64>()
            .floor() as u64
            % 1291021111;
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
        if seqs.is_empty() || ws.iter().all(|&w| w <= 0.001) {
            return self;
        }
        let max_len = seqs
            .iter()
            .zip(ws)
            .filter(|&(_, &w)| w > 0.001)
            .map(|(s, _)| s.len())
            .max()
            .unwrap_or(0);
        let node_num_thr = (max_len as f64 * coef).floor() as usize;
        rand::seq::index::sample(&mut rng, seqs.len(), seqs.len())
            .into_iter()
            .map(|idx| (&seqs[idx], ws[idx]))
            .filter(|&(_, w)| w > 0.001)
            .fold(self, |x, (y, w)| {
                if x.nodes.len() > node_num_thr {
                    x.add(y, w, params).remove_node(thr)
                } else {
                    x.add(y, w, params)
                }
            })
            .remove_node(thr)
            .finalize()
    }
    /// Construct a new POA graph from a given `seq` with `w` weight.
    pub fn new(seq: &[u8], w: f64) -> Self {
        let weight = w;
        let mut nodes: Vec<_> = seq
            .windows(2)
            .enumerate()
            .map(|(idx, pair)| {
                let mut n = Base::new(pair[0]);
                n.add(pair[1], w, idx + 1);
                if idx == 0 {
                    n.is_head = true;
                    n.head_weight += w;
                }
                n
            })
            .collect();
        let mut last = Base::new(*seq.last().unwrap());
        last.is_tail = true;
        nodes.push(last);
        nodes.iter_mut().for_each(|n| n.add_weight(w));
        assert!(nodes.iter().all(|n| n.weight() > 0.));
        let (dfs_flag, dfs_stack, mapping) = (vec![], vec![], vec![]);
        Self {
            weight,
            nodes,
            dfs_stack,
            dfs_flag,
            mapping,
        }
    }
    /// Align `seq` to `self` with `ins del score` parameters and `d` diameters.
    /// Note that the alignment is "semi" global.
    /// In other words, there can be unaligned trainling/leading nodes in the POA graph.
    pub fn align_banded<F>(
        &self,
        seq: &[u8],
        (ins, del, score): (i32, i32, F),
        d: usize,
    ) -> (i32, TraceBack)
    where
        F: Fn(u8, u8) -> i32 + Clone,
    {
        // Search for the starting node. use the first `d` rows to determine the starting point.
        let edges: Vec<Vec<_>> = {
            let mut edges = self.reverse_edges();
            edges.iter_mut().for_each(|eds| {
                if !eds.is_empty() {
                    eds.iter_mut().for_each(|p| *p += 1);
                } else {
                    eds.push(0);
                }
            });
            edges
        };
        let starting_node = {
            let ps = (ins, del, score.clone());
            let seq = &seq[..d / 2];
            let (_, tb) = self.align(seq, ps);
            let node: usize = tb
                .into_iter()
                .take_while(|e| match &e {
                    EditOp::Match(_) => false,
                    _ => true,
                })
                .fold(0, |acc, e| match &e {
                    EditOp::Deletion(_) => acc + 1,
                    _ => acc,
                });
            node
        };
        // Initialize large DP table. It requires O(mn) actually.
        // -----> query position ---->
        // 0 8 8 8 88 8 8 8 88
        // 0
        // 0
        // |
        // |
        // Graph position
        // |
        // v
        let (column, row) = (seq.len() + 1, self.nodes.len() + 1);
        let small = -100_000;
        let mut dp = vec![small; column * row];
        for (j, elm) in dp.iter_mut().enumerate().take(column) {
            *elm = ins * j as i32;
        }
        for i in 0..row {
            dp[i * column] = 0;
        }
        // Filling  DP matrix sparsely.
        let starting_node = starting_node.max(d / 2) - d / 2;
        let filled_cells = self.determine_cells(starting_node, d, seq.len());
        for (g, q) in filled_cells {
            let pos = (g + 1) * column + (q + 1);
            let mut max = dp[pos - 1] + ins;
            for &prev in edges[g].iter() {
                let prev_pos = prev * column + (q + 1);
                let deletion = dp[prev_pos] + del;
                let mat_s = score(seq[q], self.nodes[g].base());
                let mat_s = dp[prev_pos - 1] + mat_s;
                max = max.max(deletion).max(mat_s);
            }
            dp[pos] = max;
        }
        // Traceback.
        let mut q_pos = seq.len();
        let (mut g_pos, opt_score) = (0..row)
            .map(|r| dp[r * column + column - 1])
            .enumerate()
            .max_by(|(a_i, a), (b_i, b)| match a.cmp(&b) {
                std::cmp::Ordering::Equal => a_i.cmp(&b_i),
                x => x,
            })
            .unwrap_or_else(|| panic!("{}", line!()));
        assert_eq!(dp[g_pos * column + q_pos], opt_score);
        let mut operations = vec![];
        'outer: while q_pos > 0 && g_pos > 0 {
            // Current score
            let c_score = dp[g_pos * column + q_pos];
            for &p in edges[g_pos - 1].iter() {
                let pos = p * column + q_pos;
                let del = dp[pos] + del;
                if del == c_score {
                    operations.push(EditOp::Deletion(g_pos - 1));
                    g_pos = p;
                    continue 'outer;
                }
                let mat = dp[pos - 1] + score(seq[q_pos - 1], self.nodes[g_pos - 1].base());
                if mat == c_score {
                    operations.push(EditOp::Match(g_pos - 1));
                    g_pos = p;
                    q_pos -= 1;
                    continue 'outer;
                }
            }
            // Insertion
            let ins = dp[g_pos * column + q_pos - 1] + ins;
            if ins == c_score {
                q_pos -= 1;
                operations.push(EditOp::Insertion);
                continue 'outer;
            }
            panic!("error. none of choices match the current trace table.");
        }
        while q_pos > 0 {
            operations.push(EditOp::Insertion);
            q_pos -= 1;
        }
        operations.reverse();
        (opt_score, operations)
    }
    fn determine_cells(&self, g_start: usize, d: usize, qlen: usize) -> Vec<(usize, usize)> {
        let mut filled_range = vec![(qlen, 0); self.nodes.len()];
        let (mut start, mut end) = (0, d / 2 + 1);
        let mut stack = vec![g_start];
        'search: while !stack.is_empty() {
            let node = *stack.last().unwrap();
            filled_range[node].0 = (filled_range[node].0).min(start);
            filled_range[node].1 = (filled_range[node].1).max(end.min(qlen));
            for &to in self.nodes[node].edges() {
                if start + 1 < filled_range[to].0 || filled_range[to].1 < (end + 1).min(qlen) {
                    stack.push(to);
                    start = if end < d { 0 } else { start + 1 };
                    end += 1;
                    continue 'search;
                }
            }
            stack.pop().unwrap();
            end -= 1;
            start = start.max(1) - 1;
        }
        let mut poss = vec![];
        for (g, &(s, e)) in filled_range.iter().enumerate().skip(g_start) {
            for i in s..e {
                poss.push((g, i));
            }
        }
        poss
    }
    /// Align `seq` to `self` with `param` parameters.
    /// Note that the alignment is "semi" global.
    /// In other words, there can be unaligned trainling/leading nodes in the POA graph.
    #[cfg(feature = "poa_simd")]
    pub fn align<F>(&self, seq: &[u8], (ins, del, score): (i32, i32, F)) -> (i32, TraceBack)
    where
        F: Fn(u8, u8) -> i32,
    {
        let edges = {
            let mut edges = self.reverse_edges();
            edges.iter_mut().for_each(|eds| {
                if !eds.is_empty() {
                    eds.iter_mut().for_each(|p| *p += 1);
                } else {
                    eds.push(0);
                }
            });
            edges
        };
        // -----> query position ---->
        // 0 8 8 8 88 8 8 8 88
        // 0
        // 0
        // |
        // |
        // Graph position
        // |
        // v
        use packed_simd::f32x8 as f32s;
        use packed_simd::i32x8 as i32s;
        const LANE: usize = 8;
        let row = self.nodes.len() + 1;
        let column = if seq.len() % LANE == 0 {
            seq.len() + 1
        } else {
            seq.len() + (LANE - (seq.len() % LANE)) + 1
        };
        let (mut profile, mut dp) = (vec![], vec![]);
        assert!((column - 1) % LANE == 0);
        for &base in b"ACGT" {
            for i in 0..column {
                if i < seq.len() {
                    profile.push(score(base, seq[i]))
                } else {
                    profile.push(0)
                }
            }
        }
        // Initialazation.
        for _ in 0..column * row {
            dp.push(std::i32::MIN);
        }
        // The last elements in each row. Used at the beggining of the trace-back.
        let mut last_elements = vec![std::i32::MIN];
        // Weight to break tie. [i*column + j] is the total weight to (i,j) element.
        let mut route_weight = vec![0.; column * row];
        for j in 0..column {
            dp[j] = ins * j as i32;
            route_weight[j] = j as f32;
        }
        for i in 0..row {
            dp[i * column] = 0;
        }
        let deletions = i32s::splat(del);
        for i in 1..row {
            let bs = base_table::BASE_TABLE[self.nodes[i - 1].base() as usize];
            let nw = self.nodes[i - 1].weight() as f32;
            let node_weight = f32s::splat(nw);
            // Update by SIMD instructions.
            for &p in &edges[i - 1] {
                for lj in (0..column - 1).step_by(LANE) {
                    let start = i * column + lj + 1;
                    let end = start + LANE;
                    let prev_start = p * column + lj + 1;
                    let prev_end = prev_start + LANE;
                    let current = i32s::from_slice_unaligned(&dp[start..end]);
                    let current_weight = f32s::from_slice_unaligned(&route_weight[start..end]);
                    // Update for deletion state.
                    let del_score =
                        i32s::from_slice_unaligned(&dp[prev_start..prev_end]) + deletions;
                    let del_weight =
                        f32s::from_slice_unaligned(&route_weight[prev_start..prev_end]);
                    let mask = current.gt(del_score)
                        | (current.eq(del_score) & current_weight.ge(del_weight));
                    let current = mask.select(current, del_score);
                    let current_weight = mask.select(current_weight, del_weight);
                    // Update for match state.
                    let (seq_s, seq_e) = (column * bs + lj, column * bs + lj + LANE);
                    let mat_score = i32s::from_slice_unaligned(&profile[seq_s..seq_e])
                        + i32s::from_slice_unaligned(&dp[prev_start - 1..prev_end - 1]);
                    let mat_weight = node_weight
                        + f32s::from_slice_unaligned(&route_weight[prev_start - 1..prev_end - 1]);
                    let mask = current.gt(mat_score)
                        | (current.eq(mat_score) & current_weight.gt(mat_weight));
                    mask.select(current, mat_score)
                        .write_to_slice_unaligned(&mut dp[start..end]);
                    mask.select(current_weight, mat_weight)
                        .write_to_slice_unaligned(&mut route_weight[start..end]);
                }
            }
            // Insertions would be updated by usual updates due to dependencies.
            for j in 0..seq.len() {
                let pos = i * column + j + 1;
                let current = dp[pos];
                let current_weight = route_weight[pos];
                let ins = dp[pos - 1] + ins;
                let ins_weight = route_weight[pos - 1] + 1.;
                if ins > current || (ins == current && ins_weight > current_weight) {
                    dp[pos] = ins;
                    route_weight[pos] = ins_weight;
                }
            }
            last_elements.push(dp[i * column + seq.len()]);
        }
        // Traceback.
        let mut q_pos = seq.len();
        let (mut g_pos, &score) = last_elements
            .iter()
            .enumerate()
            .max_by(|(a_i, a), (b_i, b)| match a.cmp(&b) {
                std::cmp::Ordering::Equal => a_i.cmp(&b_i),
                x => x,
            })
            .unwrap_or_else(|| panic!("{}", line!()));
        assert_eq!(dp[g_pos * column + q_pos], score);
        let mut operations = vec![];
        'outer: while q_pos > 0 && g_pos > 0 {
            let w = self.nodes[g_pos - 1].weight() as f32;
            let score = dp[g_pos * column + q_pos];
            let weight = route_weight[g_pos * column + q_pos];
            // Deletion.
            for &p in &edges[g_pos - 1] {
                let pos = p * column + q_pos;
                let (del, del_w) = (dp[pos] + del, route_weight[pos]);
                if del == score && ((del_w - weight) / del_w.max(weight)).abs() < 0.000_1 {
                    operations.push(EditOp::Deletion(g_pos - 1));
                    g_pos = p;
                    continue 'outer;
                }
            }
            // Insertion
            let ins = dp[g_pos * column + q_pos - 1] + ins;
            let ins_w = route_weight[g_pos * column + q_pos - 1] + 1.;
            if ins == score && ((ins_w - weight) / ins_w.max(weight)).abs() < 0.000_1 {
                q_pos -= 1;
                operations.push(EditOp::Insertion);
                continue 'outer;
            }
            // Match/Mismatch
            let bs = base_table::BASE_TABLE[self.nodes[g_pos - 1].base() as usize];
            for &p in &edges[g_pos - 1] {
                let mat = dp[p * column + q_pos - 1] + profile[bs * column + q_pos - 1];
                let mat_w = route_weight[p * column + q_pos - 1] + w;
                if mat == score && ((mat_w - weight) / mat_w.max(weight)).abs() < 0.000_1 {
                    operations.push(EditOp::Match(g_pos - 1));
                    g_pos = p;
                    q_pos -= 1;
                    continue 'outer;
                }
            }
            panic!("error. none of choices match the current trace table.");
        }
        while q_pos > 0 {
            operations.push(EditOp::Insertion);
            q_pos -= 1;
        }
        operations.reverse();
        (score, operations)
    }
    /// Align `seq` to `self` with `param` parameters.
    /// Note that the alignment is "semi" global.
    /// In other words, there can be unaligned trainling/leading nodes in the POA graph.
    #[cfg(not(feature = "poa_simd"))]
    pub fn align<F>(&self, seq: &[u8], (ins, del, score): (i32, i32, F)) -> (i32, TraceBack)
    where
        F: Fn(u8, u8) -> i32,
    {
        let edges = {
            let mut edges = self.reverse_edges();
            edges.iter_mut().for_each(|eds| {
                if !eds.is_empty() {
                    eds.iter_mut().for_each(|p| *p += 1);
                } else {
                    eds.push(0);
                }
            });
            edges
        };

        // -----> query position ---->
        // 0 8 8 8 88 8 8 8 88
        // 0
        // 0
        // |
        // |
        // Graph position
        // |
        // v
        let row = self.nodes.len() + 1;
        let column = seq.len() + 1;
        // Initialazation.
        let mut dp = vec![0; column * row];
        // Weight to break tie. [i*column + j] is the total weight to (i,j) element.
        let mut route_weight = vec![0.; column * row];
        for j in 0..column {
            dp[j] = ins * j as i32;
            route_weight[j] = j as f64;
        }
        for q in 1..column {
            for r in 1..row {
                // Updating position.
                let pos = r * column + q;
                let mut max = dp[pos - 1] + ins;
                let mut max_weight = route_weight[pos - 1] + 1.;
                let mat_score = score(seq[q - 1], self.nodes[r - 1].base());
                for &p in &edges[r - 1] {
                    let prev = p * column + q;
                    let del_score = dp[prev] + del;
                    let del_weight = route_weight[prev];
                    if del_score > max {
                        max = del_score;
                        max_weight = del_weight;
                    }
                    // Update for match state.
                    let mat_score = mat_score + dp[prev - 1];
                    let mat_weight = self.nodes[r - 1].weight() + route_weight[prev - 1];
                    if mat_score > max {
                        max = mat_score;
                        max_weight = mat_weight;
                    }
                }
                dp[pos] = max;
                route_weight[pos] = max_weight;
            }
        }
        // Traceback.
        let mut q_pos = seq.len();
        let (mut g_pos, max_score) = (0..row)
            .map(|g| (g, dp[column * g + column - 1]))
            .max_by(|(a_i, a), (b_i, b)| match a.cmp(&b) {
                std::cmp::Ordering::Equal => a_i.cmp(&b_i),
                x => x,
            })
            .unwrap_or_else(|| panic!("{}", line!()));
        assert_eq!(dp[g_pos * column + q_pos], max_score);
        let mut operations = vec![];
        'outer: while q_pos > 0 && g_pos > 0 {
            let w = self.nodes[g_pos - 1].weight() as f32;
            let c_score = dp[g_pos * column + q_pos];
            let weight = route_weight[g_pos * column + q_pos];
            // Deletion.
            for &p in &edges[g_pos - 1] {
                let pos = p * column + q_pos;
                let (del, del_w) = (dp[pos] + del, route_weight[pos]);
                if del == c_score && ((del_w - weight) / del_w.max(weight)).abs() < 0.000_1 {
                    operations.push(EditOp::Deletion(g_pos - 1));
                    g_pos = p;
                    continue 'outer;
                }
            }
            // Insertion
            let ins = dp[g_pos * column + q_pos - 1] + ins;
            let ins_w = route_weight[g_pos * column + q_pos - 1] + 1.;
            if ins == c_score && ((ins_w - weight) / ins_w.max(weight)).abs() < 0.000_1 {
                q_pos -= 1;
                operations.push(EditOp::Insertion);
                continue 'outer;
            }
            // Match/Mismatch
            for &p in &edges[g_pos - 1] {
                let mat = score(seq[q_pos - 1], self.nodes[g_pos - 1].base())
                    + dp[p * column + q_pos - 1];
                let mat_w = route_weight[p * column + q_pos - 1] + w as f64;
                if mat == c_score && ((mat_w - weight) / mat_w.max(weight)).abs() < 0.000_1 {
                    operations.push(EditOp::Match(g_pos - 1));
                    g_pos = p;
                    q_pos -= 1;
                    continue 'outer;
                }
            }
            panic!("error. none of choices match the current trace table.");
        }
        while q_pos > 0 {
            operations.push(EditOp::Insertion);
            q_pos -= 1;
        }
        operations.reverse();
        (max_score, operations)
    }
    /// Integrate a single sequence to `self` with default parameter and `d` diameter.
    pub fn add_default_banded(self, seq: &[u8], d: usize) -> Self {
        self.add_banded(seq, (-2, -2, &|x, y| if x == y { 1 } else { -1 }), d)
    }
    /// Integrate a single sequence with default parameter.
    pub fn add_default(self, seq: &[u8], w: f64) -> Self {
        self.add(seq, w, (-2, -2, &|x, y| if x == y { 1 } else { -1 }))
    }
    /// Integrate a single sequence with the given parameter and diameter.
    pub fn add_banded<F>(self, seq: &[u8], ps: (i32, i32, &F), d: usize) -> Self
    where
        F: Fn(u8, u8) -> i32,
    {
        if self.weight < SMALL || self.nodes.is_empty() {
            return Self::new(seq, 1.);
        }
        let (_, traceback) = self.align_banded(seq, ps, d);
        self.integrate_alignment(seq, 1., traceback)
    }
    /// Integrate a signle seqnece with `w` weight with the given parameter.
    pub fn add<F>(self, seq: &[u8], w: f64, parameters: (i32, i32, &F)) -> Self
    where
        F: Fn(u8, u8) -> i32,
    {
        if self.weight < SMALL || self.nodes.is_empty() {
            return Self::new(seq, w);
        }
        let (_, traceback) = self.align(seq, parameters);
        self.integrate_alignment(seq, w, traceback)
    }
    fn integrate_alignment(mut self, seq: &[u8], w: f64, traceback: TraceBack) -> Self {
        let mut q_pos = 0;
        let mut previous: Option<usize> = None;
        for operation in traceback {
            match operation {
                EditOp::Match(to) => {
                    let base = seq[q_pos];
                    let position = if self.nodes[to].base() == base {
                        to
                    } else {
                        let mut new_node = Base::new(seq[q_pos]);
                        new_node.ties.push(to);
                        let position = self.nodes.len();
                        self.nodes.push(new_node);
                        self.nodes[to].ties.push(position);
                        position
                    };
                    self.nodes[position].add_weight(w);
                    if q_pos == seq.len() - 1 {
                        self.nodes[position].is_tail = true;
                        self.nodes[position].tail_weight += w;
                    } else if q_pos == 0 {
                        self.nodes[position].is_head = true;
                        self.nodes[position].head_weight += w;
                    }
                    if let Some(p) = previous {
                        self.nodes[p].add(base, w, position);
                    };
                    previous = Some(position);
                    q_pos += 1;
                }
                EditOp::Insertion => {
                    let base = seq[q_pos];
                    let mut new_node = Base::new(base);
                    new_node.add_weight(w);
                    if q_pos == seq.len() - 1 {
                        new_node.is_tail = true;
                        new_node.tail_weight += w;
                    } else if q_pos == 0 {
                        new_node.is_head = true;
                        new_node.head_weight += w;
                    }
                    self.nodes.push(new_node);
                    if let Some(p) = previous {
                        let position = self.nodes.len() - 1;
                        self.nodes[p].add(base, w, position);
                    }
                    previous = Some(self.nodes.len() - 1);
                    q_pos += 1;
                }
                EditOp::Deletion(_) | EditOp::Stop => {}
            }
        }
        assert_eq!(q_pos, seq.len());
        self.weight += w;
        // assert!(self.nodes.iter().all(|node| node.weight() > 0.));
        self.topological_sort()
    }
    /// Return the edges.
    /// Note that there are two kind of edges. One is usual edges, i.e., the edges link consective bases. The other is "tie"-edge, which ties two different bases that have the same child node. In other words, if there is a mis-match like
    /// AACATT
    /// ||XX||
    /// AAGTTT,
    /// the C and G would be tied by a "tie"-edge.
    /// This "tie"-edge is necessary, as I want to align AACT to the resulting graph with NO MISMATCH.
    /// Of course, this choice is arbitrary, and it may be changed in the future.
    pub fn edges(&self) -> Vec<Vec<usize>> {
        let mut edges = vec![vec![]; self.nodes.len()];
        for (from, n) in self.nodes.iter().enumerate() {
            for &to in n.edges.iter() {
                edges[from].push(to);
            }
            for &tied in n.ties.iter() {
                for &to in self.nodes[tied].edges.iter() {
                    if !edges[from].contains(&to) {
                        edges[from].push(to);
                    }
                }
            }
        }
        edges
    }
    /// Return the edges to the reverse directions.
    /// It is useful when you want to back-tracing from some node in the graph.
    pub fn reverse_edges(&self) -> Vec<Vec<usize>> {
        let mut edges = vec![vec![]; self.nodes.len()];
        for (from, n) in self.nodes.iter().enumerate() {
            for &to in n.edges.iter() {
                edges[to].push(from);
            }
            for &tied in n.ties.iter() {
                for &to in self.nodes[tied].edges.iter() {
                    if !edges[to].contains(&from) {
                        edges[to].push(from);
                    }
                }
            }
        }
        edges
    }
    fn finalize(mut self) -> Self {
        let bases: Vec<_> = self.nodes.iter().map(|e| e.base).collect();
        self.nodes.iter_mut().for_each(|e| e.finalize(&bases));
        self
    }
}
