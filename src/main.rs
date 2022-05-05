use itertools::iproduct;
use itertools::Itertools;
use na::*;
use nalgebra as na;
use scheme_grid::SchemeGrid;
use std::collections::HashMap;

type Row3 = RowVector3<u8>;

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

#[derive(Debug, Clone, Copy)]
enum CsDirection {
    Forward,
    Backward,
}

#[derive(Debug, Clone, Copy)]
struct SimpleCS<const SAME: usize, const DIFF: usize> {
    same: [usize; SAME],
    different: [(usize, CsDirection); DIFF],
}

trait CollisionStructure {
    fn same(&self) -> &[usize];
    fn different(&self) -> &[(usize, CsDirection)];
}

impl<const SAME: usize, const DIFF: usize> CollisionStructure for SimpleCS<SAME, DIFF> {
    fn same(&self) -> &[usize] {
        &self.same
    }
    fn different(&self) -> &[(usize, CsDirection)] {
        &self.different
    }
}

fn generate_all_cs_2() -> Vec<Box<dyn CollisionStructure>> {
    let cs1 = SimpleCS::<0, 2> {
        same: [],
        different: [(1, CsDirection::Backward), (1, CsDirection::Backward)],
    };
    let cs2 = SimpleCS::<0, 2> {
        same: [],
        different: [(1, CsDirection::Forward), (1, CsDirection::Backward)],
    };
    vec![Box::new(cs1), Box::new(cs2)]
}

fn is_free(v: RowVector5<u8>, fixed: &[RowVector5<u8>]) -> bool {
    // println!("fixed: {:?}", fixed);
    // println!("free:  {:?}", v);
    let matrix = na::OMatrix::<u8, Dynamic, U5>::from_rows(fixed)
        .cast::<f64>()
        .transpose();
    let v = v.cast().transpose();
    let linear_combination = matrix.clone().svd(true, true).solve(&v, 0.0001).unwrap();
    // println!("solution: {linear_combination:?}");
    let check = matrix * linear_combination;
    // println!("check: {check:?}\nv: {v:?}");
    if (check - v).norm() < 0.01 {
        return false;
    }
    true
}
use std::ops::Deref;

impl Linicrypt<1, 5> {
    fn check_for_cs_3_2_1(&self, cs: impl Deref<Target = dyn CollisionStructure>) -> bool {
        let same = cs.same().iter().map(|i| &self.constraints[*i]);
        let mut fixed: Vec<_> = same.into_iter().flat_map(|c| [c.k, c.x, c.y]).collect();
        fixed.push(self.m);
        // Check 2: the i^* query is unconstraint on both sides
        let (i_star, dir_star) = cs.different()[0];
        let c_star = &self.constraints[i_star];
        let (free_1, free_2) = match dir_star {
            CsDirection::Forward => (c_star.k, c_star.x),
            CsDirection::Backward => (c_star.k, c_star.y),
        };
        if !is_free(free_1, &fixed) || !is_free(free_2, &fixed) {
            return false;
        }

        // Check 3: Every query is onconstrained on one side
        for (i, dir) in cs.different() {
            let c = &self.constraints[*i];
            let (should_be_free, fixed_1, fixed_2) = match dir {
                CsDirection::Forward => (c.y, c.k, c.x),
                CsDirection::Backward => (c.x, c.k, c.y),
            };
            fixed.push(fixed_1);
            fixed.push(fixed_2);
            if !is_free(should_be_free, &fixed) {
                return false;
            }
        }

        true
    }
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

fn generate_3_2_1_programs() -> Vec<Linicrypt<1, 5>> {
    let ms = generate_all_vecs::<5, 4>([0, 0, 0, 1]);
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

fn linicrypt_to_lines(p: &Linicrypt<1, 5>) -> Vec<String> {
    let m_line = format!(" M={}{}{}{}{}", p.m[0], p.m[1], p.m[2], p.m[3], p.m[4]);
    let c0 = &p.constraints[0];
    let c1 = &p.constraints[1];
    let k1_line = format!("1k={}{}{}{}{}", c0.k[0], c0.k[1], c0.k[2], c0.k[3], c0.k[4]);
    let x1_line = format!("1x={}{}{}{}{}", c0.x[0], c0.x[1], c0.x[2], c0.x[3], c0.x[4]);
    let y1_line = format!("1y={}{}{}{}{}", c0.y[0], c0.y[1], c0.y[2], c0.y[3], c0.y[4]);
    let k2_line = format!("2k={}{}{}{}{}", c1.k[0], c1.k[1], c1.k[2], c1.k[3], c1.k[4]);
    let x2_line = format!("2x={}{}{}{}{}", c1.x[0], c1.x[1], c1.x[2], c1.x[3], c1.x[4]);
    let y2_line = format!("2y={}{}{}{}{}", c1.y[0], c1.y[1], c1.y[2], c1.y[3], c1.y[4]);

    let cs = SimpleCS::<0, 2> {
        same: [],
        different: [(1, CsDirection::Backward), (1, CsDirection::Backward)],
    };

    let all_cs = generate_all_cs_2();
    let mut cs_infos = all_cs
        .into_iter()
        .map(|cs| {
            if p.check_for_cs_3_2_1(cs) {
                "CS".into()
            } else {
                "No CS".into()
            }
        })
        .collect();

    let lines = vec![m_line, k1_line, x1_line, y1_line, k2_line, x2_line, y2_line];
    lines.append(&mut cs_infos);
    lines
}

fn collision_structure_examples() {
    let programs = generate_3_2_1_programs();

    // idea:
    // generate all (16) possible collision structures for 3_2_1
    // check each of p for all collision structures
    let cs = SimpleCS::<1, 1> {
        same: [0],
        different: [(1, CsDirection::Forward)],
    };

    let programs_with_cs: Vec<_> = programs
        .into_iter()
        .filter(|p| {
            if p.check_for_cs_3_2_1(&cs) {
                println!("CS works for p");
                println!("{:?}", p);
                println!("{:?}", cs);
                return true;
            }
            false
        })
        .collect();
    // compute statistics
    print_grid(&programs_with_cs, linicrypt_to_lines, 5);
}

fn main() {
    compression_functions();

    println!("\n-------------------------------------------------------------------\n");

    collision_structure_examples();
}

use term_grid::{Cell, Direction, Filling, Grid, GridOptions};
fn print_grid<const BASE: usize, const OUT: usize>(
    programs: &[Linicrypt<BASE, OUT>],
    to_lines: impl Fn(&Linicrypt<BASE, OUT>) -> Vec<String> + Copy,
    columns: usize,
) {
    let mut grid = Grid::new(GridOptions {
        filling: Filling::Spaces(3),
        direction: Direction::LeftToRight,
    });

    for row in programs.chunks(columns) {
        let mut row_of_lines: Vec<Vec<String>> = row.iter().map(to_lines).collect();
        let num_lines = row_of_lines[0].len();
        // fill with emtpy blocks
        while row_of_lines.len() < columns {
            row_of_lines.push(vec!["".into(); num_lines])
        }
        // empty line above
        for _i in 0..columns {
            grid.add(Cell::from(""))
        }
        for i in 0..num_lines {
            for lines in &row_of_lines {
                grid.add(Cell::from(lines[i].clone()));
            }
        }
    }

    println!("{}", grid.fit_into_columns(columns))
}
