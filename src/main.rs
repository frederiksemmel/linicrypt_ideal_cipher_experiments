use grid::{Cell, Grid};
use itertools::iproduct;
use na::*;
use nalgebra as na;
use std::collections::HashMap;

type Row = RowVector3<u8>;

struct CompressionFunction {
    m: Row,
    k: Row,
    x: Row,
    y: Row,
}

impl CompressionFunction {
    fn to_cell(&self) -> Cell {
        let m_line = format!("M={}{}{}", self.m[0], self.m[1], self.m[2]);
        let k_line = format!("k={}{}{}", self.k[0], self.k[1], self.k[2]);
        let x_line = format!("x={}{}{}", self.x[0], self.x[1], self.x[2]);
        let y_line = format!("y={}{}{}", self.y[0], self.y[1], self.y[2]);
        let lines = vec![m_line, k_line, x_line, y_line];
        Cell { lines }
    }
}

fn generate_all_schemes() -> Vec<CompressionFunction> {
    let ms = [(0, 0, 1), (0, 1, 1), (1, 0, 1), (1, 1, 1)].map(|(a, b, c)| Row::new(a, b, c));
    let ks = [(0, 0, 1), (0, 1, 1), (1, 0, 1), (1, 1, 1)].map(|(a, b, c)| Row::new(a, b, c));
    let xs = [(0, 0, 1), (0, 1, 1), (1, 0, 1), (1, 1, 1)].map(|(a, b, c)| Row::new(a, b, c));
    let y = Row::new(0, 0, 1);

    iproduct!(ms, ks, xs)
        .map(|(m, k, x)| CompressionFunction { m, k, x, y })
        .collect()
}

fn is_y_unconstrained(f: &CompressionFunction) -> bool {
    let matrix = na::Matrix3::from_rows(&[f.m, f.k, f.x]).cast::<f64>();
    let linear_combination = matrix.lu().solve(&f.y.cast().transpose());
    linear_combination.is_none()
}
fn is_x_unconstrained(f: &CompressionFunction) -> bool {
    let matrix = na::Matrix3::from_rows(&[f.m, f.k, f.y]).cast::<f64>();
    let linear_combination = matrix.lu().solve(&f.x.cast().transpose());
    linear_combination.is_none()
}

#[derive(Debug, PartialEq, Eq, Hash)]
enum SchemeType {
    Degenerate,
    A,
    B,
    Secure,
}

fn compute_scheme_type(f: &CompressionFunction) -> SchemeType {
    let x_free = is_x_unconstrained(f);
    let y_free = is_y_unconstrained(f);
    match (x_free, y_free) {
        (true, true) => SchemeType::Degenerate,
        (true, false) => SchemeType::A,
        (false, true) => SchemeType::B,
        (false, false) => SchemeType::Secure,
    }
}

fn main() {
    let schemes = generate_all_schemes();
    let mut counter = HashMap::new();

    for f in schemes {
        let scheme_type = compute_scheme_type(&f);
        let mut cell = f.to_cell();
        cell.lines.push(format!("{:?}", scheme_type));

        let c = counter.entry(scheme_type).or_insert(0);
        *c += 1;
    }
}

mod grid {
    pub struct Grid {
        pub columns: usize,
        pub cells: Vec<Cell>,
    }

    pub struct Cell {
        pub lines: Vec<String>,
    }
}
