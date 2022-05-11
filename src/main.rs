use itertools::iproduct;
use itertools::Itertools;
use na::*;
use nalgebra as na;
use scheme_grid::{print_grid, print_grid_schemes};
use std::collections::HashMap;

mod compression_schemes;
mod scheme_grid;

const EPSILON: f64 = 0.0001;

#[derive(Debug, Clone)]
pub enum Operation {
    E,
    D,
}

#[derive(Debug, Clone)]
pub struct Constraint<const BASE: usize> {
    pub op: Operation,
    pub k: RowSVector<u8, BASE>,
    pub x: RowSVector<u8, BASE>,
    pub y: RowSVector<u8, BASE>,
}

#[derive(Debug)]
pub struct Linicrypt<const OUT: usize, const BASE: usize, const N: usize> {
    m: SMatrix<u8, OUT, BASE>,
    constraints: [Constraint<BASE>; N],
}

type RawConstraint<const BASE: usize> = (Operation, [u8; BASE], [u8; BASE], [u8; BASE]);
impl<const BASE: usize, const N: usize> Linicrypt<1, BASE, N> {
    pub fn new(m: [u8; BASE], cs: [RawConstraint<BASE>; N]) -> Self {
        let constraints = cs.map(|(op, k, x, y)| Constraint {
            op,
            k: RowSVector::from_row_slice(&k),
            x: RowSVector::from_row_slice(&x),
            y: RowSVector::from_row_slice(&y),
        });
        Linicrypt::<1, BASE, N> {
            m: SMatrix::from_row_slice(&m),
            constraints,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Direction {
    F,
    B,
    N,
}

use std::fmt;
impl fmt::Display for Direction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Direction::F => write!(f, "F"),
            Direction::B => write!(f, "B"),
            Direction::N => write!(f, " "),
        }
    }
}

// #[derive(Debug, Clone, Copy)]
// struct SimpleCs<const SAME: usize, const DIFF: usize> {
//     same: [usize; SAME],
//     different: [(usize, CsDirection); DIFF],
// }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CollisionStructure<const N: usize> {
    permutation: [usize; N],
    i_star: usize,
    // first i_star-1 directions are ignored.
    // This is to avoid adding a second const type param
    cs_type: [Direction; N],
}

impl<const N: usize> CollisionStructure<N> {
    fn same(&self) -> &[usize] {
        &self.permutation[..self.i_star]
    }
    fn different(&self) -> impl Iterator<Item = (usize, Direction)> + '_ {
        let diff = &self.permutation[self.i_star..];
        let cs_type = &self.cs_type[self.i_star..];
        diff.iter().copied().zip(cs_type.iter().copied())
    }
    fn i_star(&self) -> (usize, Direction) {
        (self.permutation[self.i_star], self.cs_type[self.i_star])
    }
}
impl CollisionStructure<2> {
    fn id(&self) -> String {
        let perm = format!("{}{}", self.permutation[0], self.permutation[1]);
        let cs_type = format!("{}{}", self.cs_type[0], self.cs_type[1]);
        format!("{perm},{},{cs_type}", self.i_star)
    }
}

fn generate_all_cs_2() -> Vec<CollisionStructure<2>> {
    use Direction::*;
    let mut css = vec![];
    for i_star in 0..2 {
        let perms = (0..2).permutations(2);
        let types = (0..i_star)
            .map(|_| vec![N])
            .chain((i_star..2).map(|_| vec![F, B]))
            .multi_cartesian_product();
        for (permutation, cs_type) in iproduct!(perms, types) {
            let cs = CollisionStructure::<2> {
                permutation: permutation.try_into().unwrap(),
                i_star,
                cs_type: cs_type.try_into().unwrap(),
            };
            css.push(cs);
        }
    }
    css
}

// fn is_in_span(v: RowVector5<u8>, fixed: &[RowVector5<u8>]) -> bool {
//     let matrix = na::OMatrix::<u8, Dynamic, U5>::from_rows(fixed)
//         .cast::<f64>()
//         .transpose();
//     println!("matrix: {matrix}");
//     let v = v.cast().transpose();
//     let linear_combination = matrix.clone().svd(true, true).solve(&v, EPSILON).unwrap();
//     let check = matrix * &linear_combination;
//     println!("check: {check}\nv: {v}\nl: {linear_combination}");
//     if (check - v).norm() < 0.1 {
//         return true;
//     }
//     false
// }
fn is_in_span(v: RowVector5<u8>, fixed: &[RowVector5<u8>]) -> bool {
    let matrix = na::OMatrix::<u8, Dynamic, U5>::from_rows(fixed).cast::<f64>();
    let rows_with_v: Vec<_> = fixed.iter().copied().chain([v].into_iter()).collect();
    let matrix_with_v = na::OMatrix::<u8, Dynamic, U5>::from_rows(&rows_with_v).cast::<f64>();
    if matrix.svd(false, false).rank(EPSILON) == matrix_with_v.svd(false, false).rank(EPSILON) {
        return true;
    }
    false
}
fn full_rank(rows: &[RowVector5<u8>]) -> bool {
    // println!("fixed: {:?}", fixed);
    // println!("free:  {:?}", v);
    let matrix = na::OMatrix::<u8, Dynamic, U5>::from_rows(rows)
        .cast::<f64>()
        .transpose();
    let rank = matrix.svd(false, false).rank(EPSILON);
    rank < min(rows.len(), 5)
}

impl Linicrypt<1, 5, 2> {
    fn has_cs_3_2_1(&self, cs: CollisionStructure<2>) -> bool {
        let same = cs
            .same()
            // .permutation
            .iter()
            // .take(cs.i_star)
            .map(|i| &self.constraints[*i]);
        let mut fixed: Vec<_> = same.into_iter().flat_map(|c| [c.k, c.x, c.y]).collect();
        fixed.push(self.m);
        // Check 2: the i^* query is unconstraint on both sides
        let (i_star, dir_star) = cs.i_star();
        let c_star = &self.constraints[i_star];
        let (free_1, free_2) = match dir_star {
            Direction::F => (c_star.k, c_star.x),
            Direction::B => (c_star.k, c_star.y),
            Direction::N => unreachable!(),
        };
        if is_in_span(free_1, &fixed) && is_in_span(free_2, &fixed) {
            // println!("Cond 2 not fulfilled");
            return false;
        }

        // Check 3: Every query is onconstrained on one side
        for (i, dir) in cs.different() {
            // println!("{i}");
            let c = &self.constraints[i];
            let (should_be_free, fixed_1, fixed_2) = match dir {
                Direction::F => (c.y, c.k, c.x),
                Direction::B => (c.x, c.k, c.y),
                Direction::N => unreachable!(),
            };
            fixed.push(fixed_1);
            fixed.push(fixed_2);
            if is_in_span(should_be_free, &fixed) {
                // println!("Cond 3 not fulfilled at {i}");
                return false;
            }
            fixed.push(should_be_free);
        }

        true
    }
    fn is_degenerate(&self) -> bool {
        let mut vecs: Vec<_> = self
            .constraints
            .iter()
            .flat_map(|c| [c.k, c.x, c.y])
            .collect();
        vecs.push(self.m);
        full_rank(&vecs)
    }
}

fn generate_all_vecs<const BASE: usize, const DIM: usize>(
    last_entries: [u8; DIM],
) -> impl Iterator<Item = RowVector5<u8>> {
    (0..(BASE - DIM))
        .map(|_| 0..=1)
        .multi_cartesian_product()
        .map(move |v| RowVector5::<u8>::from_iterator(v.into_iter().chain(last_entries)))
}

fn generate_3_2_1_programs() -> Vec<Linicrypt<1, 5, 2>> {
    let ms = generate_all_vecs::<5, 0>([]);
    let ks1: Vec<_> = generate_all_vecs::<5, 2>([0, 0]).collect();
    let xs1: Vec<_> = generate_all_vecs::<5, 2>([0, 0]).collect();
    let y1 = RowVector5::<u8>::new(0, 0, 0, 1, 0);
    let cs1: Vec<_> = iproduct!(ks1, xs1)
        .map(|(k, x)| Constraint {
            op: Operation::E,
            k,
            x,
            y: y1,
        })
        .collect();

    let ks2: Vec<_> = generate_all_vecs::<5, 1>([0]).collect();
    let xs2: Vec<_> = generate_all_vecs::<5, 1>([0]).collect();
    let y2 = RowVector5::<u8>::new(0, 0, 0, 0, 1);
    let cs2: Vec<_> = iproduct!(ks2, xs2)
        .map(|(k, x)| Constraint {
            op: Operation::E,
            k,
            x,
            y: y2,
        })
        .collect();

    iproduct!(ms, cs1, cs2)
        .map(|(m, c1, c2)| Linicrypt {
            m,
            constraints: [c1, c2],
        })
        .collect()
}

fn compression_functions() {
    let schemes = compression_schemes::generate_all_schemes();
    let mut counter = HashMap::new();

    for f in &schemes {
        let scheme_type = f.collision_structure_type();
        let c = counter.entry(scheme_type).or_insert(0);
        *c += 1;
    }

    print_grid_schemes(&schemes, compression_schemes::scheme_to_lines, 4);
    println!("Counter {:?}", counter);
}

fn linicrypt_to_lines(p: &Linicrypt<1, 5, 2>) -> Vec<String> {
    let m_line = format!(" M={}{}{}{}{}", p.m[0], p.m[1], p.m[2], p.m[3], p.m[4]);
    let c0 = &p.constraints[0];
    let c1 = &p.constraints[1];
    let k1_line = format!("1k={}{}{}{}{}", c0.k[0], c0.k[1], c0.k[2], c0.k[3], c0.k[4]);
    let x1_line = format!("1x={}{}{}{}{}", c0.x[0], c0.x[1], c0.x[2], c0.x[3], c0.x[4]);
    let y1_line = format!("1y={}{}{}{}{}", c0.y[0], c0.y[1], c0.y[2], c0.y[3], c0.y[4]);
    let k2_line = format!("2k={}{}{}{}{}", c1.k[0], c1.k[1], c1.k[2], c1.k[3], c1.k[4]);
    let x2_line = format!("2x={}{}{}{}{}", c1.x[0], c1.x[1], c1.x[2], c1.x[3], c1.x[4]);
    let y2_line = format!("2y={}{}{}{}{}", c1.y[0], c1.y[1], c1.y[2], c1.y[3], c1.y[4]);
    vec![m_line, k1_line, x1_line, y1_line, k2_line, x2_line, y2_line]
}

pub fn print_linicrypt(p: &Linicrypt<1, 5, 2>) {
    let lines = linicrypt_to_lines(p);
    for line in lines {
        println!("{line}");
    }
}

pub fn linicrypt_to_lines_infos(p: &Linicrypt<1, 5, 2>) -> Vec<String> {
    let all_cs = generate_all_cs_2();
    let mut cs_infos = all_cs
        .into_iter()
        .map(|cs| {
            if p.has_cs_3_2_1(cs) {
                format!("Y{}", cs.id())
            } else {
                format!(" {}", cs.id())
            }
        })
        .collect();

    let mut lines = linicrypt_to_lines(p);
    lines.append(&mut cs_infos);
    lines
}

fn collision_structure_examples() {
    let programs = generate_3_2_1_programs();
    let total_programs = programs.len();

    // idea:
    // generate all (16) possible collision structures for 3_2_1
    // check each of p for all collision structures
    let css = generate_all_cs_2();

    // for p in &programs[0..10000] {
    //     print_linicrypt(p);
    //     for cs in &css {
    //         println!("{}", cs.id());
    //         if p.has_cs_3_2_1(*cs) {
    //             println!("P has cs!")
    //         }
    //     }
    // }

    let interesting_ps: Vec<_> = programs
        .into_iter()
        .filter(|p| {
            !p.is_degenerate()
                && css
                    .iter()
                    .map(|cs| if p.has_cs_3_2_1(*cs) { 1 } else { 0 })
                    .sum::<usize>()
                    <= 2
        })
        .collect();

    print_grid(&interesting_ps, linicrypt_to_lines_infos, 7);
    println!("{} of {total_programs}", interesting_ps.len());
}

fn main() {
    compression_functions();

    println!("\n-------------------------------------------------------------------\n");

    collision_structure_examples();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn check_all_cs_2() {
        use super::Direction::*;
        let manual = vec![
            CollisionStructure::<2> {
                permutation: [0, 1],
                i_star: 0,
                cs_type: [F, F],
            },
            CollisionStructure::<2> {
                permutation: [0, 1],
                i_star: 0,
                cs_type: [F, B],
            },
            CollisionStructure::<2> {
                permutation: [0, 1],
                i_star: 0,
                cs_type: [B, F],
            },
            CollisionStructure::<2> {
                permutation: [0, 1],
                i_star: 0,
                cs_type: [B, B],
            },
            CollisionStructure::<2> {
                permutation: [1, 0],
                i_star: 0,
                cs_type: [F, F],
            },
            CollisionStructure::<2> {
                permutation: [1, 0],
                i_star: 0,
                cs_type: [F, B],
            },
            CollisionStructure::<2> {
                permutation: [1, 0],
                i_star: 0,
                cs_type: [B, F],
            },
            CollisionStructure::<2> {
                permutation: [1, 0],
                i_star: 0,
                cs_type: [B, B],
            },
            CollisionStructure::<2> {
                permutation: [0, 1],
                i_star: 1,
                cs_type: [N, F],
            },
            CollisionStructure::<2> {
                permutation: [0, 1],
                i_star: 1,
                cs_type: [N, B],
            },
            CollisionStructure::<2> {
                permutation: [1, 0],
                i_star: 1,
                cs_type: [N, F],
            },
            CollisionStructure::<2> {
                permutation: [1, 0],
                i_star: 1,
                cs_type: [N, B],
            },
        ];
        let automatic = generate_all_cs_2();
        for (m, a) in manual.iter().zip(automatic.iter()) {
            println!("{:?}", a);
            println!("{:?}", m);
            assert_eq!(a, m);
        }
    }

    #[test]
    fn check_cs_split() {
        use super::Direction::*;

        let cs0 = CollisionStructure::<2> {
            permutation: [0, 1],
            i_star: 0,
            cs_type: [F, B],
        };
        assert_eq!(cs0.same(), &[]);
        let mut different = cs0.different();
        assert_eq!(different.next(), Some((0, F)));
        assert_eq!(different.next(), Some((1, B)));
        assert_eq!(different.next(), None);

        let cs1 = CollisionStructure::<2> {
            permutation: [0, 1],
            i_star: 1,
            cs_type: [N, B],
        };
        assert_eq!(cs1.same(), &[0]);
        let mut different = cs1.different();
        assert_eq!(different.next(), Some((1, B)));
        assert_eq!(different.next(), None);

        let cs1 = CollisionStructure::<2> {
            permutation: [1, 0],
            i_star: 1,
            cs_type: [N, F],
        };
        assert_eq!(cs1.same(), &[1]);
        let mut different = cs1.different();
        assert_eq!(different.next(), Some((0, F)));
        assert_eq!(different.next(), None);
    }

    #[test]
    fn check_linicrypt_3_2_1_cs_1() {
        use super::Direction::*;
        use super::Operation::*;

        let cs = CollisionStructure::<2> {
            permutation: [0, 1],
            i_star: 0,
            cs_type: [B, F],
        };
        let p = Linicrypt::<1, 5, 2>::new(
            [0, 1, 0, 0, 1],
            [
                (E, [1, 0, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0]),
                (E, [0, 0, 1, 1, 0], [0, 1, 0, 0, 0], [0, 0, 0, 0, 1]),
            ],
        );

        assert!(!p.has_cs_3_2_1(cs));
    }

    #[test]
    fn check_linicrypt_3_2_1_cs_2() {
        use super::Direction::*;
        use super::Operation::*;

        let cs = CollisionStructure::<2> {
            permutation: [0, 1],
            i_star: 0,
            cs_type: [B, F],
        };
        let p = Linicrypt::<1, 5, 2>::new(
            [0, 0, 1, 0, 1],
            [
                (E, [0, 1, 0, 0, 0], [1, 0, 0, 0, 0], [0, 0, 0, 1, 0]),
                (E, [1, 1, 1, 1, 0], [0, 1, 0, 1, 0], [0, 0, 0, 0, 1]),
            ],
        );

        assert!(!p.has_cs_3_2_1(cs));
    }

    #[test]
    fn check_linicrypt_3_2_1_cs_3() {
        use super::Direction::*;
        use super::Operation::*;

        let cs = CollisionStructure::<2> {
            permutation: [0, 1],
            i_star: 0,
            cs_type: [F, B],
        };
        let p = Linicrypt::<1, 5, 2>::new(
            [0, 0, 0, 0, 1],
            [
                (E, [0, 0, 1, 0, 0], [0, 1, 0, 0, 0], [0, 0, 0, 1, 0]),
                (E, [1, 1, 1, 0, 0], [0, 1, 0, 1, 0], [0, 0, 0, 0, 1]),
            ],
        );

        assert!(!p.has_cs_3_2_1(cs));
    }
}
