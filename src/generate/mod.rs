use crate::PartialOrderAlignment;
impl PartialOrderAlignment {
    /// Generate consensus from POA graph.
    /// This function is totally immature.
    /// Please read the code, and give me some recommendations.
    pub fn consensus(&self) -> Vec<u8> {
        let (_, mut node) = self.start_node().unwrap_or_else(|| (0, &self.nodes[0]));
        let mut res = vec![node.base()];
        let end_node = self.end_node().map(|e| e.0);
        while node.has_edge() {
            let (&idx, _) = node
                .edges()
                .iter()
                .zip(node.weights.iter())
                .max_by(|a, b| (a.1).partial_cmp(&(b.1)).unwrap())
                .unwrap_or_else(|| panic!("{}", line!()));
            node = &self.nodes[idx];
            res.push(node.base());
            if Some(idx) == end_node {
                break;
            }
        }
        res
    }
    fn end_node(&self) -> Option<(usize, &super::Base)> {
        self.nodes
            .iter()
            .enumerate()
            .filter(|e| e.1.is_tail)
            .max_by(|a, b| (a.1).tail_weight.partial_cmp(&(b.1).tail_weight).unwrap())
    }
    fn start_node(&self) -> Option<(usize, &super::Base)> {
        self.nodes
            .iter()
            .enumerate()
            .filter(|e| e.1.is_head)
            .max_by(|a, b| a.1.head_weight().partial_cmp(&b.1.head_weight()).unwrap())
    }
    /// Generate consensus from POA graph.
    /// This function cares homopolymer-region
    /// by first compress the homopolymer region in the POA graph.
    pub fn consensus_homopolymer(&self) -> Vec<u8> {
        // Grouping nodes, i.e., compress repetitive nodes.
        let group: Vec<_> = {
            let mut group: Vec<_> = (0..self.nodes.len()).collect();
            let mut members: Vec<_> = (0..self.nodes.len()).map(|e| vec![e]).collect();
            for (from, node) in self.nodes.iter().enumerate() {
                for &to in node.edges.iter() {
                    if group[from] != group[to]
                        && node.base() == self.nodes[to].base()
                        && !self.reachable(from, to, &members[group[from]], &members[group[to]])
                    {
                        let pop: Vec<_> = members[group[to]].iter().copied().collect();
                        members[group[to]].clear();
                        for &m in pop.iter() {
                            group[m] = group[from];
                        }
                        members[group[from]].extend(pop);
                    }
                }
            }
            use std::collections::HashMap;
            let deduped_group: HashMap<_, _> = {
                let mut group = group.clone();
                group.sort();
                group.dedup();
                group
                    .into_iter()
                    .enumerate()
                    .map(|(idx, g)| (g, idx))
                    .collect()
            };
            group.iter().map(|g| deduped_group[g]).collect()
        };
        // Never panic.
        let max_group = group.iter().max().unwrap();
        let mut coverages = vec![0.; max_group + 1];
        let mut bases = vec![0; max_group + 1];
        for (node, &g) in self.nodes.iter().zip(group.iter()) {
            coverages[g] += node.weight;
            bases[g] = node.base();
        }
        let median = select_nth_by(&self.nodes, self.nodes.len() / 2, |n| n.weight).unwrap();
        // let median = select_nth_by(&coverages, coverages.len() / 2, |&c| c).unwrap();
        if log_enabled!(log::Level::Trace) {
            debug!("MEDIAN\t{}", median);
            for (idx, (n, &g)) in self.nodes.iter().zip(group.iter()).enumerate() {
                let cov = coverages[g];
                let rep = ((cov + median / 2.) / median).floor();
                let base = n.base() as char;
                let edges: Vec<_> = n.edges.iter().map(|x| format!("{}", x)).collect();
                let edges = edges.join(",");
                debug!("G\t{}\t{}\t{}\t{}\t{}\t{}", idx, g, base, cov, rep, edges);
            }
        }
        let repetitive: Vec<_> = coverages
            .iter()
            .map(|c| (c / median + 0.5).floor() as usize)
            .collect();
        let mut current_group = self.start_node().map(|x| group[x.0]).unwrap_or(0);
        let end_group = self.end_node().map(|x| group[x.0]);
        let mut res = vec![bases[current_group]; repetitive[current_group]];
        let edges = {
            let mut edges: Vec<Vec<(usize, f64)>> = vec![vec![]; self.nodes.len()];
            for (from, node) in self.nodes.iter().enumerate() {
                for (&to, &w) in node.edges.iter().zip(node.weights.iter()) {
                    let from = group[from];
                    let to = group[to];
                    if from == to {
                        continue;
                    }
                    if let Some(res) = edges[from].iter_mut().find(|x| x.0 == to) {
                        res.1 += w;
                    } else {
                        edges[from].push((to, w))
                    }
                }
            }
            for es in edges.iter_mut() {
                let sum = es.iter().map(|x| x.1).sum::<f64>();
                es.iter_mut().for_each(|x| x.1 /= sum);
            }
            edges
        };
        if log_enabled!(log::Level::Trace) {
            for (from, es) in edges.iter().enumerate() {
                for (to, w) in es.iter() {
                    debug!("EDGE\t{}\t{}\t{}", from, to, w);
                }
            }
        }
        while !edges[current_group].is_empty() {
            let &(idx, _) = edges[current_group]
                .iter()
                .max_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap())
                .unwrap_or_else(|| panic!("{}", line!()));
            current_group = idx;
            res.extend(std::iter::repeat(bases[current_group]).take(repetitive[current_group]));
            if Some(current_group) == end_group {
                break;
            }
        }
        res
    }
    // Determine wether there is a path from `source` to `sink` without using
    // the (from, to) edge.
    // Note that, as the self is topolocally sorted,
    // we can early drop if we reached a node with the index more than to.
    fn reachable(&self, from: usize, to: usize, source: &[usize], sink: &[usize]) -> bool {
        // Never panic.
        let max = *sink.iter().max().unwrap();
        let mut stack: Vec<_> = source.iter().filter(|&&e| e < max).copied().collect();
        while !stack.is_empty() {
            let node = stack.pop().unwrap();
            assert!(node < max);
            for &child in self.nodes[node].edges.iter() {
                // Direct edge is not allowed.
                if sink.contains(&child) && (node, child) != (from, to) {
                    return true;
                } else if child < max {
                    stack.push(child);
                }
            }
        }
        false
    }
}

// Return the n-th smallest element. O(log xs.len()) usually.
// If you force it to run in O(log xs.len()) almost always, use `rand` crate
// and randomize the pivot selection like `let pivot = xs.choose(&mut rng).unwrap();`.
use std::cmp::{PartialEq, PartialOrd};
fn select_nth_by<T: Clone, F: Fn(&T) -> K, K>(xs: &[T], n: usize, f: F) -> Option<K>
where
    K: PartialOrd + PartialEq,
{
    if xs.len() <= n {
        return None;
    }
    let pivot = f(&xs[xs.len() / 2]);
    let small = xs.iter().filter(|x| f(x) < pivot).count();
    let same = xs.iter().filter(|x| f(x) == pivot).count();
    // Recursive call.
    if n < small {
        // We can remove elements more than `pivot` from `xs`.
        let xs: Vec<_> = xs.iter().filter(|x| f(&x) < pivot).cloned().collect();
        select_nth_by(&xs, n, f)
    } else if small + same <= n {
        let xs: Vec<_> = xs.iter().filter(|x| f(&x) > pivot).cloned().collect();
        select_nth_by(&xs, n - small - same, f)
    } else {
        assert!(small <= n && n < small + same);
        Some(pivot)
    }
}

#[cfg(test)]
mod tests {
    use crate::*;
    use rand::seq::SliceRandom;
    #[test]
    fn short_consensus() {
        let answer = b"AACGACATGCTAGTCA";
        let seqs = vec![
            b"AACGACATGCTAGTCA".to_vec(),
            b"AAACGACTATGCTAGTCA".to_vec(),
            b"AACGAATGCTAGTCA".to_vec(),
            b"AACGACATGCAGTCA".to_vec(),
            b"AACGACATGCTAGTA".to_vec(),
            b"ACACGACATGCTATGTCA".to_vec(),
        ];
        let res = POA::from_vec_default(&seqs);
        let consensus = res.consensus();
        assert_eq!(
            &consensus,
            answer,
            "{}\t{}",
            String::from_utf8_lossy(&consensus),
            String::from_utf8_lossy(answer)
        );
    }
    #[test]
    fn long_consensus() {
        let answer = b"CAGCTAGTCAGTCGATCGATGGGCAGTGATGCATGATCGGACTTACGATCGACTGACTAGTCA".to_vec();
        let seqs = vec![
            b"CAGCTAGTCAGTCGATCGATGGGCAGTGATGCATGATCGGACTTACGATCGATGACTAGTCA".to_vec(),
            b"CAGCTAGTCAGTCGATCGATGGGCAGTGATGCATGATCGGACTTACGATCGACTGACTAGTCA".to_vec(),
            b"CAGCTAGTCTGTCGATCGATGGGAGTGAGCATGTCGGACTACGATCGACTGACTAGAATCA".to_vec(),
            b"CAGCTAGTCAGTCCATCGATGGGCAGTGATGCATGATCGGACTACGATCGACTGACTAGTCA".to_vec(),
            b"CAGCTAGTCAGTCGATCGATGGCAGTGATGCATGATCGGACTTACGTCGACTGACTAGTCA".to_vec(),
            b"CAGCTAGTCAGTCGATCGATGGGCAGTGATGCATGATCGGACTTACGATCGACTGACTAGGCA".to_vec(),
            b"CAGCTAGTCAGTCGATCATGGGCAGTGATGCATGATCGGACTTACGATCGACTGATAGTCA".to_vec(),
            b"CAGCTAGTCATCGATCGATGGGCAGTGATGCATGATCGGACTTACGTCGACTGACTTAGTCA".to_vec(),
            b"CAGCTAGTCACGATGATGGGCAGTGATGCATGATCGGACTTACGATCGACTACTAGTCA".to_vec(),
            b"CAGCTAGTCAGTCGATCGATGGGCAGTGATGCATGATCGGACTTACGATCGATGACTAGTCA".to_vec(),
        ];
        let res = POA::from_vec_default(&seqs);
        let consensus = res.consensus();
        assert_eq!(
            &consensus,
            &answer,
            "{}\t{}",
            String::from_utf8_lossy(&consensus),
            String::from_utf8_lossy(&answer)
        );
    }
    #[test]
    fn short_consensus_rand() {
        let bases = b"ACTG";
        let coverage = 120;
        let start = 20;
        let len = 10;
        let result = (start..coverage)
            .filter(|&cov| {
                let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(1_234_567);
                let template1: Vec<_> = (0..len)
                    .filter_map(|_| bases.choose(&mut rng))
                    .copied()
                    .collect();
                let seqs: Vec<_> = (0..cov)
                    .map(|_| {
                        gen_sample::introduce_randomness(&template1, &mut rng, &gen_sample::PROFILE)
                    })
                    .collect();
                let consensus = POA::from_vec_default(&seqs).consensus();
                consensus == template1
            })
            .count();
        assert!(result > 80, "{}", result);
    }
    #[test]
    fn long_consensus_rand() {
        let bases = b"ACTG";
        let coverage = 120;
        let start = 40;
        let len = 100;
        let result = (start..coverage)
            .filter(|&cov| {
                let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(1_234_567);
                let template1: Vec<_> = (0..len)
                    .filter_map(|_| bases.choose(&mut rng))
                    .copied()
                    .collect();
                let seqs: Vec<_> = (0..cov)
                    .map(|_| {
                        gen_sample::introduce_randomness(&template1, &mut rng, &gen_sample::PROFILE)
                    })
                    .collect();
                let consensus = POA::from_vec_default(&seqs).consensus();
                consensus == template1
            })
            .count();
        assert!(result > 60, "{}", result);
    }
    #[test]
    fn super_long_consensus_rand() {
        use rayon::prelude::*;
        let bases = b"ACTG";
        let coverage = 100;
        let start = 20;
        let len = 2000;
        let result = (start..coverage)
            .into_par_iter()
            .filter(|&cov| {
                let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(1_234_567);
                let template1: Vec<_> = (0..len)
                    .filter_map(|_| bases.choose(&mut rng))
                    .copied()
                    .collect();
                let seqs: Vec<_> = (0..cov)
                    .map(|_| {
                        gen_sample::introduce_randomness(&template1, &mut rng, &gen_sample::PROFILE)
                    })
                    .collect();
                let consensus = POA::from_vec_default(&seqs).consensus();
                let dist = edit_dist(&consensus, &template1);
                dist <= 2
            })
            .count();
        assert!(result > 40, "{}", result);
    }
    fn edit_dist(x1: &[u8], x2: &[u8]) -> u32 {
        let mut dp = vec![vec![0; x2.len() + 1]; x1.len() + 1];
        for (i, row) in dp.iter_mut().enumerate() {
            row[0] = i as u32;
        }
        for j in 0..=x2.len() {
            dp[0][j] = j as u32;
        }
        for (i, x1_b) in x1.iter().enumerate() {
            for (j, x2_b) in x2.iter().enumerate() {
                let m = if x1_b == x2_b { 0 } else { 1 };
                dp[i + 1][j + 1] = (dp[i][j + 1] + 1).min(dp[i + 1][j] + 1).min(dp[i][j] + m);
            }
        }
        dp[x1.len()][x2.len()]
    }
}
