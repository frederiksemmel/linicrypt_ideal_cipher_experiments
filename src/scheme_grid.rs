use crate::compression_schemes::SingleQueryScheme;
use crate::Linicrypt;
use term_grid::{Cell, Direction, Filling, Grid, GridOptions};

pub fn print_grid<const BASE: usize, const OUT: usize>(
    programs: &[Linicrypt<BASE, OUT, 2>],
    to_lines: impl Fn(&Linicrypt<BASE, OUT, 2>) -> Vec<String> + Copy,
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
pub fn print_grid_schemes(
    programs: &[SingleQueryScheme],
    to_lines: impl Fn(&SingleQueryScheme) -> Vec<String> + Copy,
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
