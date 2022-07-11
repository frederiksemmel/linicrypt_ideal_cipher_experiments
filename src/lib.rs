use na::*;
use nalgebra as na;

pub mod next;
pub mod print_grid;

const EPSILON: f64 = 0.0001;

#[derive(Debug, Clone, PartialEq)]
pub enum Operation {
    E,
    D,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Constraint<const BASE: usize> {
    pub op: Operation,
    pub k: RowSVector<u8, BASE>,
    pub x: RowSVector<u8, BASE>,
    pub y: RowSVector<u8, BASE>,
}

#[derive(Debug)]
pub struct AlgebraicRepresentation<const BASE: usize, const N: usize, const OUT: usize> {
    pub m: SMatrix<u8, OUT, BASE>,
    pub constraints: [Constraint<BASE>; N],
}

type RawConstraint<const BASE: usize> = (Operation, [u8; BASE], [u8; BASE], [u8; BASE]);
impl<const BASE: usize, const N: usize> AlgebraicRepresentation<BASE, N, 1> {
    pub fn new(m: [u8; BASE], cs: [RawConstraint<BASE>; N]) -> Self {
        let constraints = cs.map(|(op, k, x, y)| Constraint {
            op,
            k: RowSVector::from_row_slice(&k),
            x: RowSVector::from_row_slice(&x),
            y: RowSVector::from_row_slice(&y),
        });
        AlgebraicRepresentation::<BASE, N, 1> {
            m: SMatrix::from_row_slice(&m),
            constraints,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Direction {
    F,
    B,
}

use std::fmt;
impl fmt::Display for Direction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Direction::F => write!(f, "F"),
            Direction::B => write!(f, "B"),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CollisionStructure<const N: usize, const DIFF: usize> {
    pub permutation: [usize; N],
    pub cs_type: [Direction; DIFF],
}

// TODO This will probably be needed to generate collision structures dynamically
pub trait CollisionStructureTrait {
    fn same(&self) -> &[usize];
    fn different(&self) -> DifferentIter;
    fn types(&self) -> &[Direction];
    fn i_star(&self) -> (usize, Direction);
}

use std::iter::Copied;
use std::iter::Zip;
type DifferentIter<'a> =
    Zip<Copied<std::slice::Iter<'a, usize>>, Copied<std::slice::Iter<'a, Direction>>>;

impl<const N: usize, const DIFF: usize> CollisionStructure<N, DIFF> {
    pub fn same(&self) -> &[usize] {
        &self.permutation[..N - DIFF]
    }
    pub fn different(&self) -> DifferentIter {
        self.permutation[N - DIFF..]
            .iter()
            .copied()
            .zip(self.cs_type.iter().copied())
    }
    pub fn directions(&self) -> &[Direction] {
        &self.cs_type
    }
    pub fn i_star(&self) -> (usize, Direction) {
        (self.permutation[N - DIFF], self.cs_type[0])
    }
}

fn repr_slice<'a>(row: impl IntoIterator<Item = &'a (impl std::fmt::Display + 'a)>) -> String {
    row.into_iter().map(|entry| format!("{}", entry)).collect()
}

impl<const N: usize, const DIFF: usize> CollisionStructure<N, DIFF> {
    pub fn id(&self) -> String {
        let perm = repr_slice(&self.permutation);
        // let perm = format!("{}{}", self.permutation[0], self.permutation[1]);
        let cs_type = repr_slice(&self.cs_type);
        format!("{perm},{},{cs_type}", N - DIFF)
    }
}
// impl CollisionStructure<1, 1> {
//     pub fn id(&self) -> String {
//         let perm = format!("{}", self.permutation[0]);
//         format!("{perm},{},{}", 0, self.cs_type[0])
//     }
// }

fn is_in_span<const BASE: usize>(v: RowSVector<u8, BASE>, fixed: &[RowSVector<u8, BASE>]) -> bool {
    let matrix = na::OMatrix::<u8, Dynamic, Const<BASE>>::from_rows(fixed).cast::<f64>();
    let rows_with_v: Vec<_> = fixed.iter().copied().chain([v].into_iter()).collect();
    let matrix_with_v =
        na::OMatrix::<u8, Dynamic, Const<BASE>>::from_rows(&rows_with_v).cast::<f64>();
    if matrix.svd(false, false).rank(EPSILON) == matrix_with_v.svd(false, false).rank(EPSILON) {
        return true;
    }
    false
}
fn full_rank<const BASE: usize>(rows: &[RowSVector<u8, BASE>]) -> bool {
    // println!("fixed: {:?}", fixed);
    // println!("free:  {:?}", v);
    let matrix = na::OMatrix::<u8, Dynamic, Const<BASE>>::from_rows(rows)
        .cast::<f64>()
        .transpose();
    let rank = matrix.svd(false, false).rank(EPSILON);
    rank < min(rows.len(), BASE)
}

impl<const BASE: usize, const N: usize> AlgebraicRepresentation<BASE, N, 1> {
    pub fn has_cs<const I_STAR: usize>(&self, cs: &CollisionStructure<N, I_STAR>) -> bool {
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
    pub fn is_degenerate(&self) -> bool {
        let mut vecs: Vec<_> = self
            .constraints
            .iter()
            .flat_map(|c| [c.k, c.x, c.y])
            .collect();
        vecs.push(self.m);
        full_rank(&vecs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn check_cs_split() {
        use super::Direction::*;

        let cs0 = CollisionStructure::<2, 2> {
            permutation: [0, 1],
            cs_type: [F, B],
        };
        assert_eq!(cs0.same(), &[]);
        let mut different = cs0.different();
        assert_eq!(different.next(), Some((0, F)));
        assert_eq!(different.next(), Some((1, B)));
        assert_eq!(different.next(), None);

        let cs1 = CollisionStructure::<2, 1> {
            permutation: [0, 1],
            cs_type: [B],
        };
        assert_eq!(cs1.same(), &[0]);
        let mut different = cs1.different();
        assert_eq!(different.next(), Some((1, B)));
        assert_eq!(different.next(), None);

        let cs1 = CollisionStructure::<2, 1> {
            permutation: [1, 0],
            cs_type: [F],
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

        let cs = CollisionStructure::<2, 2> {
            permutation: [0, 1],
            cs_type: [B, F],
        };
        let p = AlgebraicRepresentation::<5, 2, 1>::new(
            [0, 1, 0, 0, 1],
            [
                (E, [1, 0, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0]),
                (E, [0, 0, 1, 1, 0], [0, 1, 0, 0, 0], [0, 0, 0, 0, 1]),
            ],
        );

        assert!(!p.has_cs(&cs));
    }

    #[test]
    fn check_linicrypt_3_2_1_cs_2() {
        use super::Direction::*;
        use super::Operation::*;

        let cs = CollisionStructure::<2, 2> {
            permutation: [0, 1],
            cs_type: [B, F],
        };
        let p = AlgebraicRepresentation::<5, 2, 1>::new(
            [0, 0, 1, 0, 1],
            [
                (E, [0, 1, 0, 0, 0], [1, 0, 0, 0, 0], [0, 0, 0, 1, 0]),
                (E, [1, 1, 1, 1, 0], [0, 1, 0, 1, 0], [0, 0, 0, 0, 1]),
            ],
        );

        assert!(!p.has_cs(&cs));
    }

    #[test]
    fn check_linicrypt_3_2_1_cs_3() {
        use super::Direction::*;
        use super::Operation::*;

        let cs = CollisionStructure::<2, 2> {
            permutation: [0, 1],
            cs_type: [F, B],
        };
        let p = AlgebraicRepresentation::<5, 2, 1>::new(
            [0, 0, 0, 0, 1],
            [
                (E, [0, 0, 1, 0, 0], [0, 1, 0, 0, 0], [0, 0, 0, 1, 0]),
                (E, [1, 1, 1, 0, 0], [0, 1, 0, 1, 0], [0, 0, 0, 0, 1]),
            ],
        );

        assert!(!p.has_cs(&cs));
    }
}
