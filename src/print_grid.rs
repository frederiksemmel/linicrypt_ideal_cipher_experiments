use term_grid::{Cell, Direction, Filling, Grid, GridOptions};

pub fn print_grid(cells: Vec<Vec<String>>, columns: usize) {
    let mut grid = Grid::new(GridOptions {
        filling: Filling::Spaces(3),
        direction: Direction::LeftToRight,
    });

    for row in cells.chunks(columns) {
        let mut row: Vec<_> = row.into();
        let num_lines = row[0].len();
        // fill with emtpy blocks
        while row.len() < columns {
            row.push(vec!["".into(); num_lines])
        }
        // empty line above
        for _i in 0..columns {
            grid.add(Cell::from(""))
        }

        for i in 0..num_lines {
            for lines in &row {
                grid.add(Cell::from(lines[i].clone()));
            }
        }
    }

    println!("{}", grid.fit_into_columns(columns))
}
