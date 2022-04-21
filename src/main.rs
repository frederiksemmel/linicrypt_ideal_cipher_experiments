use itertools::iproduct;
use na::*;
use nalgebra as na;
use scheme_grid::SchemeGrid;
use std::collections::HashMap;

type Row = RowVector3<u8>;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SchemeType {
    Degenerate,
    A,
    B,
    Secure,
}

pub struct Scheme {
    m: Row,
    k: Row,
    x: Row,
    y: Row,
}

impl Scheme {
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

fn generate_all_schemes() -> Vec<Scheme> {
    let ms = [(0, 0, 1), (0, 1, 1), (1, 0, 1), (1, 1, 1)].map(|(a, b, c)| Row::new(a, b, c));
    let ks = [(0, 0, 0), (0, 1, 0), (1, 0, 0), (1, 1, 0)].map(|(a, b, c)| Row::new(a, b, c));
    let xs = [(0, 0, 0), (0, 1, 0), (1, 0, 0), (1, 1, 0)].map(|(a, b, c)| Row::new(a, b, c));
    let y = Row::new(0, 0, 1);

    iproduct!(ms, ks, xs)
        .map(|(m, k, x)| Scheme { m, k, x, y })
        .collect()
}

fn generate_non_constant_schemes() -> Vec<Scheme> {
    let ms = [(0, 1, 1), (1, 0, 1), (1, 1, 1)].map(|(a, b, c)| Row::new(a, b, c));
    let ks = [(0, 1, 0), (1, 0, 0), (1, 1, 0)].map(|(a, b, c)| Row::new(a, b, c));
    let xs = [(0, 1, 0), (1, 0, 0), (1, 1, 0)].map(|(a, b, c)| Row::new(a, b, c));
    let y = Row::new(0, 0, 1);

    iproduct!(ms, ks, xs)
        .map(|(m, k, x)| Scheme { m, k, x, y })
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
}

mod scheme_grid {
    use super::Scheme;
    use std::fmt;
    use term_grid::{Cell, Direction, Filling, Grid, GridOptions};

    pub struct SchemeGrid {
        pub columns: usize,
        pub schemes: Vec<Scheme>,
    }

    impl fmt::Display for SchemeGrid {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            let mut grid = Grid::new(GridOptions {
                filling: Filling::Spaces(3),
                direction: Direction::LeftToRight,
            });

            for row in self.schemes.chunks(self.columns) {
                let mut row_of_lines: Vec<Vec<String>> = row.iter().map(scheme_to_lines).collect();
                let num_lines = row_of_lines[0].len();
                // fill with emtpy blocks
                while row_of_lines.len() < self.columns {
                    row_of_lines.push(vec!["".into(); num_lines])
                }
                // empty line above
                for _i in 0..self.columns {
                    grid.add(Cell::from(""))
                }
                for i in 0..num_lines {
                    for lines in &row_of_lines {
                        grid.add(Cell::from(lines[i].clone()));
                    }
                }
            }

            write!(f, "{}", grid.fit_into_columns(self.columns))
        }
    }

    fn scheme_to_lines(scheme: &Scheme) -> Vec<String> {
        let m_line = format!("M={}{}{}", scheme.m[0], scheme.m[1], scheme.m[2]);
        let k_line = format!("k={}{}{}", scheme.k[0], scheme.k[1], scheme.k[2]);
        let x_line = format!("x={}{}{}", scheme.x[0], scheme.x[1], scheme.x[2]);
        let y_line = format!("y={}{}{}", scheme.y[0], scheme.y[1], scheme.y[2]);
        let type_line = format!("{:?}", scheme.collision_structure_type());
        vec![m_line, k_line, x_line, y_line, type_line]
    }
}
