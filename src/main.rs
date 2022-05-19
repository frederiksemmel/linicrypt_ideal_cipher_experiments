use itertools::iproduct;
use itertools::Itertools;
use linicrypt::print_grid::print_grid;
use na::*;
use nalgebra as na;
use std::collections::HashMap;

use linicrypt::{AlgebraicRepresentation, CollisionStructure, Constraint, Direction, Operation};

fn generate_all_cs_2<const DIFF: usize>() -> impl Iterator<Item = CollisionStructure<2, DIFF>> {
    use Direction::*;
    let perms = (0..2).permutations(2);
    let types = (0..DIFF).map(|_| vec![F, B]).multi_cartesian_product();
    iproduct!(perms, types).map(|(p, t)| CollisionStructure::<2, DIFF> {
        permutation: p.try_into().unwrap(),
        cs_type: t.try_into().unwrap(),
    })
}

fn generate_all_cs_1<const DIFF: usize>() -> impl Iterator<Item = CollisionStructure<1, DIFF>> {
    use Direction::*;
    let perms = (0..1).permutations(1);
    let types = (0..DIFF).map(|_| vec![F, B]).multi_cartesian_product();
    iproduct!(perms, types).map(|(p, t)| CollisionStructure::<1, DIFF> {
        permutation: p.try_into().unwrap(),
        cs_type: t.try_into().unwrap(),
    })
}

fn generate_all_vecs<const BASE: usize, const DIM: usize>(
    last_entries: [u8; DIM],
) -> impl Iterator<Item = RowSVector<u8, BASE>> {
    (0..(BASE - DIM))
        .map(|_| 0..=1)
        .multi_cartesian_product()
        .map(move |v| RowSVector::<u8, BASE>::from_iterator(v.into_iter().chain(last_entries)))
}

fn generate_all_constraints<const BASE: usize, const ZEROS: usize>(
) -> impl Iterator<Item = Constraint<BASE>> {
    let ks1: Vec<_> = generate_all_vecs::<BASE, ZEROS>([0; ZEROS]).collect();
    let xs1: Vec<_> = generate_all_vecs::<BASE, ZEROS>([0; ZEROS]).collect();
    let mut y1 = RowSVector::<u8, BASE>::zeros();
    y1[BASE - ZEROS] = 1;
    iproduct!(ks1, xs1).map(move |(k, x)| Constraint {
        op: Operation::E,
        k,
        x,
        y: y1,
    })
}

fn generate_i_2_1_programs<const BASE: usize>() -> Vec<AlgebraicRepresentation<BASE, 2, 1>> {
    let ms = generate_all_vecs::<BASE, 1>([1]);
    let c1s: Vec<_> = generate_all_constraints::<BASE, 2>().collect();
    let c2s: Vec<_> = generate_all_constraints::<BASE, 1>().collect();

    iproduct!(ms, c1s, c2s)
        .map(|(m, c1, c2)| AlgebraicRepresentation {
            m,
            constraints: [c1, c2],
        })
        .collect()
}

fn generate_2_1_1_programs<const BASE: usize>() -> Vec<AlgebraicRepresentation<BASE, 1, 1>> {
    let ms = generate_all_vecs::<BASE, 1>([1]);
    let css: Vec<_> = generate_all_constraints::<BASE, 1>().collect();

    iproduct!(ms, css)
        .map(|(m, cs)| AlgebraicRepresentation {
            m,
            constraints: [cs],
        })
        .collect()
}

fn repr_vector<const BASE: usize>(row: RowSVector<u8, BASE>) -> String {
    row.iter().map(|entry| format!("{}", entry)).collect()
}

fn linicrypt_to_lines<const BASE: usize, const N: usize>(
    p: &AlgebraicRepresentation<BASE, N, 1>,
) -> Vec<String> {
    let m_line = format!(" M={}", repr_vector(p.m));
    let mut lines = vec![m_line];
    for i in 0..(N) {
        lines.push(format!("{i}k={}", repr_vector(p.constraints[i].k)));
        lines.push(format!("{i}x={}", repr_vector(p.constraints[i].x)));
        lines.push(format!("{i}y={}", repr_vector(p.constraints[i].y)));
    }
    lines
}

pub fn print_linicrypt<const BASE: usize, const N: usize>(p: &AlgebraicRepresentation<BASE, N, 1>) {
    let lines = linicrypt_to_lines(p);
    for line in lines {
        println!("{line}");
    }
}

pub fn linicrypt_to_lines_infos(p: &AlgebraicRepresentation<5, 2, 1>) -> Vec<String> {
    let all_cs = generate_all_cs_2::<2>();
    let mut cs_infos = all_cs
        .map(|cs| {
            if p.has_cs(&cs) {
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

// TODO clean this up, need some struct to hold cs analysis data and print cs_id's
fn check_css<const BASE: usize, const N: usize, const DIFF: usize>(
    p: &AlgebraicRepresentation<BASE, N, 1>,
    css: &[CollisionStructure<N, DIFF>],
    counter: &mut HashMap<String, usize>,
) -> (Vec<usize>, Vec<String>) {
    css.iter()
        .map(|cs| {
            let cs_id = cs.id();
            if p.has_cs(cs) {
                let out = (1, format!("Y{cs_id}"));
                *counter.entry(cs_id).or_insert(0) += 1;
                out
            } else {
                (0, format!(" {cs_id}"))
            }
        })
        .unzip()
}

fn compression_functions() {
    println!();
    println!("Analyzing all 64 compression schemes with 2 input, 1 queries and 1 output.");
    let ps = generate_2_1_1_programs::<{ 2 + 1 }>();
    let css: Vec<_> = generate_all_cs_1::<1>().collect();

    let mut counter: HashMap<String, usize> = HashMap::new();
    let mut cells = vec![];

    for p in &ps {
        let (mut _cs, mut info) = check_css(p, &css, &mut counter);
        let mut cell = linicrypt_to_lines(p);
        cell.append(&mut info);
        // print_linicrypt(p);
        cells.push(cell);
    }

    print_grid(cells, 4);
    println!("{:?}", counter);
}

fn repr_slice<'a>(row: impl IntoIterator<Item = &'a (impl std::fmt::Display + 'a)>) -> String {
    row.into_iter().map(|entry| format!("{}", entry)).collect()
}

fn print_comb_counter(
    css1: &[CollisionStructure<2, 2>],
    css2: &[CollisionStructure<2, 1>],
    counter: HashMap<Vec<usize>, usize>,
) {
    println!("These combinations of types occured.");
    println!("This is the order of types used in the binary representation of the combination.");
    for cs in css1 {
        println!("{}", cs.id());
    }
    for cs in css2 {
        println!("{}", cs.id());
    }
    let mut counter: Vec<(Vec<usize>, usize)> = counter.into_iter().collect();
    counter.sort_by_cached_key(|(comb, _count)| comb.iter().fold(0, |acc, &b| acc * 2 + b as u32));
    for (comb, count) in counter.iter() {
        if *count != 0 {
            let comb_str = repr_slice(comb);
            println!("{}: {count}", comb_str);
        }
    }
}

fn collision_structure_examples() {
    println!();
    println!("Finding interesting examples with 3 input, 2 queries and 1 output.");
    let programs = generate_i_2_1_programs::<{ 3 + 2 }>();

    let css2: Vec<_> = generate_all_cs_2::<2>().collect();
    let css1: Vec<_> = generate_all_cs_2::<1>().collect();

    let mut counter: HashMap<String, usize> = HashMap::new();
    let all_ids = css2
        .iter()
        .map(|cs| cs.id())
        .chain(css1.iter().map(|cs| cs.id()));
    for id in all_ids {
        counter.insert(id, 0);
    }

    let mut combination_counter: HashMap<Vec<usize>, usize> = HashMap::new();
    let all_combinations = (0..12).into_iter().map(|_| 0..=1).multi_cartesian_product();
    for comb in all_combinations {
        combination_counter.insert(comb, 0);
    }

    let non_degenerate: Vec<_> = programs
        .into_iter()
        .filter(|p| !p.is_degenerate())
        .collect();

    let mut cells = vec![];
    for p in &non_degenerate {
        let mut cell = linicrypt_to_lines(p);
        let (mut cs_2, mut info_2) = check_css(p, &css2, &mut counter);
        cell.append(&mut info_2);
        let (mut cs_1, mut info_1) = check_css(p, &css1, &mut counter);
        cell.append(&mut info_1);

        let num_cs_2: usize = cs_2.iter().sum();
        let num_cs_1: usize = cs_1.iter().sum();

        cs_2.append(&mut cs_1);
        let combination_of_cs = cs_2;
        *combination_counter.get_mut(&combination_of_cs).unwrap() += 1;

        if num_cs_2 + num_cs_1 <= 2 {
            cells.push(cell);
        }
    }

    print_grid(cells, 8);
    println!("{:#?}", counter);
    print_comb_counter(&css2, &css1, combination_counter);
}

fn secure_4_2_1() {
    println!();
    println!(
        "Finding a program with 4 inputs, making 2 queries to E without a collision structure"
    );
    let ps = generate_i_2_1_programs::<{ 4 + 2 }>();
    let css2: Vec<_> = generate_all_cs_2::<2>().collect();
    let css1: Vec<_> = generate_all_cs_2::<1>().collect();
    for p in &ps {
        if p.is_degenerate() {
            continue;
        }
        let has_cs_1 = css1.iter().any(|cs| p.has_cs(cs));
        let has_cs_2 = css2.iter().any(|cs| p.has_cs(cs));
        if !has_cs_1 && !has_cs_2 {
            print_linicrypt(p);
            break;
        }
    }
}

fn main() {
    compression_functions();
    collision_structure_examples();
    secure_4_2_1();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn check_all_cs_2_2() {
        use super::Direction::*;
        let manual = vec![
            CollisionStructure::<2, 2> {
                permutation: [0, 1],
                cs_type: [F, F],
            },
            CollisionStructure::<2, 2> {
                permutation: [0, 1],
                cs_type: [F, B],
            },
            CollisionStructure::<2, 2> {
                permutation: [0, 1],
                cs_type: [B, F],
            },
            CollisionStructure::<2, 2> {
                permutation: [0, 1],
                cs_type: [B, B],
            },
            CollisionStructure::<2, 2> {
                permutation: [1, 0],
                cs_type: [F, F],
            },
            CollisionStructure::<2, 2> {
                permutation: [1, 0],
                cs_type: [F, B],
            },
            CollisionStructure::<2, 2> {
                permutation: [1, 0],
                cs_type: [B, F],
            },
            CollisionStructure::<2, 2> {
                permutation: [1, 0],
                cs_type: [B, B],
            },
        ];
        let automatic = generate_all_cs_2::<2>();
        for (m, a) in manual.iter().zip(automatic) {
            println!("{:?}", a);
            println!("{:?}", m);
            assert_eq!(a, *m);
        }
    }

    #[test]
    fn check_all_cs_2_1() {
        use super::Direction::*;
        let manual = vec![
            CollisionStructure::<2, 1> {
                permutation: [0, 1],
                cs_type: [F],
            },
            CollisionStructure::<2, 1> {
                permutation: [0, 1],
                cs_type: [B],
            },
            CollisionStructure::<2, 1> {
                permutation: [1, 0],
                cs_type: [F],
            },
            CollisionStructure::<2, 1> {
                permutation: [1, 0],
                cs_type: [B],
            },
        ];
        let automatic = generate_all_cs_2::<1>();
        for (m, a) in manual.iter().zip(automatic) {
            println!("{:?}", a);
            println!("{:?}", m);
            assert_eq!(a, *m);
        }
    }

    #[test]
    fn check_generate_constraints_5() {
        let mut cs = generate_all_constraints::<5, 1>();
        let manual_c = Constraint {
            op: Operation::E,
            k: RowVector5::new(0, 0, 0, 0, 0),
            x: RowVector5::new(0, 0, 0, 0, 0),
            y: RowVector5::new(0, 0, 0, 0, 1),
        };
        assert_eq!(cs.next(), Some(manual_c));
        let manual_c = Constraint {
            op: Operation::E,
            k: RowVector5::new(0, 0, 0, 0, 0),
            x: RowVector5::new(0, 0, 0, 1, 0),
            y: RowVector5::new(0, 0, 0, 0, 1),
        };
        assert_eq!(cs.next(), Some(manual_c));
        let mut cs = generate_all_constraints::<5, 2>();
        let manual_c = Constraint {
            op: Operation::E,
            k: RowVector5::new(0, 0, 0, 0, 0),
            x: RowVector5::new(0, 0, 0, 0, 0),
            y: RowVector5::new(0, 0, 0, 1, 0),
        };
        assert_eq!(cs.next(), Some(manual_c));
        let manual_c = Constraint {
            op: Operation::E,
            k: RowVector5::new(0, 0, 0, 0, 0),
            x: RowVector5::new(0, 0, 1, 0, 0),
            y: RowVector5::new(0, 0, 0, 1, 0),
        };
        assert_eq!(cs.next(), Some(manual_c));
    }
    #[test]
    fn check_generate_constraints_6() {
        let mut cs = generate_all_constraints::<6, 1>();
        let manual_c = Constraint {
            op: Operation::E,
            k: RowVector6::new(0, 0, 0, 0, 0, 0),
            x: RowVector6::new(0, 0, 0, 0, 0, 0),
            y: RowVector6::new(0, 0, 0, 0, 0, 1),
        };
        assert_eq!(cs.next(), Some(manual_c));
        let manual_c = Constraint {
            op: Operation::E,
            k: RowVector6::new(0, 0, 0, 0, 0, 0),
            x: RowVector6::new(0, 0, 0, 0, 1, 0),
            y: RowVector6::new(0, 0, 0, 0, 0, 1),
        };
        assert_eq!(cs.next(), Some(manual_c));
        let mut cs = generate_all_constraints::<6, 2>();
        let manual_c = Constraint {
            op: Operation::E,
            k: RowVector6::new(0, 0, 0, 0, 0, 0),
            x: RowVector6::new(0, 0, 0, 0, 0, 0),
            y: RowVector6::new(0, 0, 0, 0, 1, 0),
        };
        assert_eq!(cs.next(), Some(manual_c));
        let manual_c = Constraint {
            op: Operation::E,
            k: RowVector6::new(0, 0, 0, 0, 0, 0),
            x: RowVector6::new(0, 0, 0, 1, 0, 0),
            y: RowVector6::new(0, 0, 0, 0, 1, 0),
        };
        assert_eq!(cs.next(), Some(manual_c));
    }
}
