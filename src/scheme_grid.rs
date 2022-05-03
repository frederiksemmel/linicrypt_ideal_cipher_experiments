use super::SingleQueryScheme;
use std::fmt;
use term_grid::{Cell, Direction, Filling, Grid, GridOptions};

pub struct SchemeGrid {
    pub columns: usize,
    pub schemes: Vec<SingleQueryScheme>,
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

fn scheme_to_lines(scheme: &SingleQueryScheme) -> Vec<String> {
    let m_line = format!("M={}{}{}", scheme.m[0], scheme.m[1], scheme.m[2]);
    let k_line = format!("k={}{}{}", scheme.k[0], scheme.k[1], scheme.k[2]);
    let x_line = format!("x={}{}{}", scheme.x[0], scheme.x[1], scheme.x[2]);
    let y_line = format!("y={}{}{}", scheme.y[0], scheme.y[1], scheme.y[2]);
    let type_line = format!("{:?}", scheme.collision_structure_type());
    vec![m_line, k_line, x_line, y_line, type_line]
}
