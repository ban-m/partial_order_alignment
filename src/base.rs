use super::base_table::BASE_TABLE;
use super::{LAMBDA_INS, LAMBDA_MATCH};
use std::fmt;
/// Node of the POA graph.
#[derive(Default, Clone)]
pub struct Base {
    /// Base (b'A', b'C', b'G', or b'T').
    pub base: u8,
    /// Edges. Each element is the index of the parent `nodes`.
    pub edges: Vec<usize>,
    /// Mismatch tie. To see the information about ties, see "POA::edges()".
    pub ties: Vec<usize>,
    /// Weights of each edges. After `POA::finalize()`, it is normalized so that
    /// the sum equals to one.
    pub weights: Vec<f64>,
    /// Raw base count of the out-edges.
    pub base_count: [f64; 4],
    /// Weight of the node.
    pub weight: f64,
    /// Weight of the node, especially the situation where
    /// this node is the first character of the sequence.
    pub head_weight: f64,
    /// Weight of the node, especially when
    /// this node is the last chatacter of the sequence.
    pub tail_weight: f64,
    /// Flag. If true, there was a sequence where this node was the first character.
    pub is_tail: bool,
    /// Flag. If true, there was a sequence where this node was the last character.
    pub is_head: bool,
    /// Heviest base.
    pub heaviest: u8,
}

impl Base {
    pub fn new(base: u8) -> Self {
        Self {
            base,
            edges: vec![],
            ties: vec![],
            weights: vec![],
            weight: 0.,
            head_weight: 0.,
            tail_weight: 0.,
            base_count: [0.; 4],
            is_tail: false,
            is_head: false,
            heaviest: 0,
        }
    }
    pub fn weight(&self) -> f64 {
        self.weight
    }
    pub fn head_weight(&self) -> f64 {
        self.head_weight
    }
    /// Remove edges to `node`
    pub fn remove(&mut self, node: usize) {
        if let Some(idx) = self.edges.iter().position(|&to| to == node) {
            self.edges.remove(idx);
            self.weights.remove(idx);
        }
        if let Some(idx) = self.ties.iter().position(|&to| to == node) {
            self.ties.remove(idx);
        }
    }
    pub fn remove_if(&mut self, mapping: &[(usize, bool)]) {
        let retain: Vec<_> = self.edges.iter().map(|&idx| mapping[idx].1).collect();
        {
            let mut idx = 0;
            self.weights.retain(|_| {
                idx += 1;
                retain[idx - 1]
            });
        }
        self.edges.retain(|&idx| mapping[idx].1);
        self.edges.iter_mut().for_each(|idx| *idx = mapping[*idx].0);
        self.ties.retain(|&idx| mapping[idx].1);
        self.ties.iter_mut().for_each(|idx| *idx = mapping[*idx].0);
        assert_eq!(self.edges.len(), self.weights.len());
    }
    pub fn remove_edges(&mut self, e: &[bool]) {
        {
            let mut idx = 0;
            self.weights.retain(|_| {
                idx += 1;
                e[idx - 1]
            });
        }
        {
            let mut idx = 0;
            self.edges.retain(|_| {
                idx += 1;
                e[idx - 1]
            });
        }
    }
    pub fn finalize(&mut self, bases: &[u8]) {
        let tot = self.base_count.iter().sum::<f64>();
        if tot > 0.001 {
            self.base_count.iter_mut().for_each(|e| *e /= tot);
        } else {
            self.base_count.iter_mut().for_each(|e| *e = 0.25);
        }
        let tot = self.weights.iter().sum::<f64>();
        self.weights.iter_mut().for_each(|e| *e /= tot);
        if let Some((_, &argmax)) = self
            .weights
            .iter()
            .zip(self.edges())
            .max_by(|a, b| (a.0).partial_cmp(&b.0).unwrap())
        {
            self.heaviest = bases[argmax];
        }
    }
    pub fn add(&mut self, b: u8, w: f64, idx: usize) {
        let pos = match self.edges.iter().enumerate().find(|&(_, &to)| to == idx) {
            Some((pos, _)) => pos,
            None => {
                self.edges.push(idx);
                self.weights.push(0.);
                self.edges.len() - 1
            }
        };
        self.weights[pos] += w;
        self.base_count[BASE_TABLE[b as usize]] += w;
    }
    pub fn add_weight(&mut self, w: f64) {
        self.weight += w;
    }
    pub fn add_head_weight(&mut self, w: f64) {
        self.head_weight += w;
    }
    pub fn rename_by(&mut self, map: &[usize]) {
        self.edges.iter_mut().for_each(|e| *e = map[*e]);
        self.ties.iter_mut().for_each(|e| *e = map[*e]);
    }
    pub fn base(&self) -> u8 {
        self.base
    }
    pub fn to(&self, to: usize) -> f64 {
        *self
            .edges
            .iter()
            .zip(self.weights.iter())
            .find(|&(&idx, _)| idx == to)
            .unwrap()
            .1
    }
    pub fn prob_with(&self, base: u8, config: &super::Config, from: &Self) -> f64 {
        let p = from.base_count[BASE_TABLE[base as usize]];
        let q = if self.base == base {
            1. - config.mismatch
        } else {
            config.mismatch / 3.
        };
        p * LAMBDA_MATCH + (1. - LAMBDA_MATCH) * q
    }
    pub fn prob(&self, base: u8, config: &super::Config) -> f64 {
        let p = self.base_count[BASE_TABLE[base as usize]];
        let q = if base == self.heaviest {
            1. - config.mismatch
        } else {
            config.mismatch / 3.
        };
        p * LAMBDA_MATCH + (1. - LAMBDA_MATCH) * q
    }
    #[inline]
    pub fn insertion(&self, base: u8) -> f64 {
        let p = self.base_count[BASE_TABLE[base as usize]];
        let q = 0.25;
        p * LAMBDA_INS + (1. - LAMBDA_INS) * q
    }
    #[inline]
    pub fn has_edge(&self) -> bool {
        !self.edges.is_empty()
    }
    pub fn edges(&self) -> &[usize] {
        &self.edges
    }
    pub fn weights_except<'a>(&'a self, e: &'a [bool]) -> impl Iterator<Item = &'a f64> {
        assert!(e.len() == self.edges.len() && e.len() == self.weights.len());
        self.weights
            .iter()
            .zip(e.iter())
            .filter_map(move |(w, b)| if !b { Some(w) } else { None })
    }
    pub fn weights(&self) -> &[f64] {
        &self.weights
    }
}

impl fmt::Display for Base {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "Base\t{}", self.base as char)?;
        let weights: Vec<_> = self.weights.iter().map(|x| format!("{:.3}", x)).collect();
        write!(f, "{}", weights.join("\t"))?;
        for to in self.edges.iter() {
            writeln!(f, "Edge\t{}", to)?;
        }
        writeln!(f, "Is tail\t{}", self.is_tail)?;
        write!(f, "Is head\t{}", self.is_head)
    }
}

impl fmt::Debug for Base {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(
            f,
            "{}\t{}/{}\t{:.2}\t{:.2}",
            self.base as char, self.is_head, self.is_tail, self.weight, self.head_weight,
        )?;
        for (w, to) in self.weights.iter().zip(self.edges.iter()) {
            writeln!(f, "E\t{}\t{:.3}", to, w)?;
        }
        for &b in b"ATGC" {
            let count = self.base_count[BASE_TABLE[b as usize]];
            if count > 0.001 {
                write!(f, "{}:{:.3} ", b as char, count)?;
            }
        }
        Ok(())
    }
}
