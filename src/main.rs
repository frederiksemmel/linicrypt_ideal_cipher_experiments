use itertools::iproduct;
use itertools::Itertools;
use na::*;
use nalgebra as na;
use scheme_grid::SchemeGrid;
use std::collections::HashMap;

type Row3 = RowVector3<u8>;
type Row5 = RowVector3<u8>;

mod scheme_grid;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
enum SchemeType {
    Degenerate,
    A,
    B,
    Secure,
}

pub struct SingleQueryScheme {
    m: Row3,
    k: Row3,
    x: Row3,
    y: Row3,
}

#[derive(Debug, Clone)]
enum Operation {
    E,
    D,
}

#[derive(Debug, Clone)]
enum Direction {
    Forward,
    Backward,
}

struct Constraint5 {
    op: Operation,
    k: Row5,
    x: Row5,
    y: Row5,
}

struct Out1In3Query2 {
    m: Row5,
    c1: Constraint5,
    c2: Constraint5,
}

struct CollisionStructure<const BASE: usize> {
    i_star: usize,
    same: Vec<Constraint<BASE>>,
    different: Vec<(Constraint<BASE>, Direction)>,
}

#[derive(Debug, Clone)]
struct Constraint<const BASE: usize> {
    op: Operation,
    k: RowSVector<u8, BASE>,
    x: RowSVector<u8, BASE>,
    y: RowSVector<u8, BASE>,
}

#[derive(Debug)]
struct Linicrypt<const OUT: usize, const BASE: usize> {
    m: SMatrix<u8, OUT, BASE>,
    constraints: Vec<Constraint<BASE>>,
}

impl SingleQueryScheme {
    fn is_y_unconstrained(&self) -> bool {
        let matrix = na::Matrix3::from_rows(&[self.m, self.k, self.x]).cast::<f64>();
        let linear_combination = matrix.lu().solve(&self.y.cast().transpose());
        linear_combination.is_none()
    }
    fn is_x_unconstrained(&self) -> bool {
        let matrix = na::Matrix3::from_rows(&[self.m, self.k, self.y]).cast::<f64>();
        let linear_combination = matrix.lu().solve(&self.x.cast().transpose());
        linear_combination.is_none()
    }
    fn collision_structure_type(&self) -> SchemeType {
        let y_free = self.is_y_unconstrained();
        let x_free = self.is_x_unconstrained();
        match (y_free, x_free) {
            (true, true) => SchemeType::Degenerate,
            (true, false) => SchemeType::A,
            (false, true) => SchemeType::B,
            (false, false) => SchemeType::Secure,
        }
    }
}

fn generate_all_schemes() -> Vec<SingleQueryScheme> {
    let ms = [(0, 0, 1), (0, 1, 1), (1, 0, 1), (1, 1, 1)].map(|(a, b, c)| Row3::new(a, b, c));
    let ks = [(0, 0, 0), (0, 1, 0), (1, 0, 0), (1, 1, 0)].map(|(a, b, c)| Row3::new(a, b, c));
    let xs = [(0, 0, 0), (0, 1, 0), (1, 0, 0), (1, 1, 0)].map(|(a, b, c)| Row3::new(a, b, c));
    let y = Row3::new(0, 0, 1);

    iproduct!(ms, ks, xs)
        .map(|(m, k, x)| SingleQueryScheme { m, k, x, y })
        .collect()
}

fn generate_non_constant_schemes() -> Vec<SingleQueryScheme> {
    let ms = [(0, 1, 1), (1, 0, 1), (1, 1, 1)].map(|(a, b, c)| Row3::new(a, b, c));
    let ks = [(0, 1, 0), (1, 0, 0), (1, 1, 0)].map(|(a, b, c)| Row3::new(a, b, c));
    let xs = [(0, 1, 0), (1, 0, 0), (1, 1, 0)].map(|(a, b, c)| Row3::new(a, b, c));
    let y = Row3::new(0, 0, 1);

    iproduct!(ms, ks, xs)
        .map(|(m, k, x)| SingleQueryScheme { m, k, x, y })
        .collect()
}

fn generate_all_vecs<const BASE: usize, const DIM: usize>(
    last_entries: [u8; DIM],
) -> impl Iterator<Item = RowVector5<u8>> {
    (0..(BASE - DIM))
        .map(|_| 0..=1)
        .multi_cartesian_product()
        .map(move |v| RowVector5::<u8>::from_iterator(v.into_iter().chain(last_entries)))
}

fn generate_3_input_2_query_schemes() -> Vec<Linicrypt<1, 5>> {
    let ms = generate_all_vecs::<5, 1>([1]);
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
    let y2 = RowVector5::<u8>::new(0, 0, 0, 1, 0);
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
            constraints: vec![c1, c2],
        })
        .collect()
}

fn compression_functions() {
    let schemes = generate_all_schemes();
    let _non_constan_schemes = generate_non_constant_schemes();
    let mut counter = HashMap::new();

    for f in &schemes {
        let scheme_type = f.collision_structure_type();
        let c = counter.entry(scheme_type).or_insert(0);
        *c += 1;
    }

    let scheme_grid = SchemeGrid {
        columns: 4,
        schemes,
    };

    println!("{}", scheme_grid);
    println!("Counter {:?}", counter);
}

fn main() {
    compression_functions();

    println!("\n-------------------------------------------------------------------\n");

    for vec in generate_3_input_2_query_schemes() {
        println!("{vec:?}");
    }
}
