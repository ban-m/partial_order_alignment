use crate::PartialOrderAlignment;
impl PartialOrderAlignment {
    /// Generate consensus from POA graph.
    /// This function is totally immature.
    /// Please read the code, and give me some recommendations.
    pub fn consensus(&self) -> Vec<u8> {
        let mut node = self
            .nodes
            .iter()
            .filter(|e| e.is_head)
            .max_by(|a, b| a.head_weight().partial_cmp(&b.head_weight()).unwrap())
            .unwrap_or_else(|| &self.nodes[0]);
        let mut res = vec![node.base()];
        let end_nodes = self
            .nodes
            .iter()
            .enumerate()
            .filter(|e| e.1.is_tail)
            .max_by(|a, b| (a.1).tail_weight.partial_cmp(&(b.1).tail_weight).unwrap())
            .map(|x| x.0);
        while node.has_edge() {
            let (&idx, _) = node
                .edges()
                .iter()
                .zip(node.weights.iter())
                .max_by(|a, b| (a.1).partial_cmp(&(b.1)).unwrap())
                .unwrap_or_else(|| panic!("{}", line!()));
            node = &self.nodes[idx];
            res.push(node.base());
            if Some(idx) == end_nodes {
                break;
            }
        }
        res
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
        let start = 20;
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
        assert!(result > 80, "{}", result);
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
