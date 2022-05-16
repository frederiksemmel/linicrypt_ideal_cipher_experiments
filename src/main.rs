use itertools::iproduct;
use itertools::Itertools;
use linicrypt::scheme_grid::{print_grid, print_grid_schemes};
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

fn generate_all_vecs<const BASE: usize, const DIM: usize>(
    last_entries: [u8; DIM],
) -> impl Iterator<Item = RowVector5<u8>> {
    (0..(BASE - DIM))
        .map(|_| 0..=1)
        .multi_cartesian_product()
        .map(move |v| RowVector5::<u8>::from_iterator(v.into_iter().chain(last_entries)))
}

fn generate_3_2_1_programs() -> Vec<AlgebraicRepresentation<1, 5, 2>> {
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
        .map(|(m, c1, c2)| AlgebraicRepresentation {
            m,
            constraints: [c1, c2],
        })
        .collect()
}

fn compression_functions() {
    let schemes = linicrypt::compression_schemes::generate_all_schemes();
    let mut counter = HashMap::new();

    for f in &schemes {
        let scheme_type = f.collision_structure_type();
        let c = counter.entry(scheme_type).or_insert(0);
        *c += 1;
    }

    print_grid_schemes(&schemes, linicrypt::compression_schemes::scheme_to_lines, 4);
    println!("Counter {:?}", counter);
}

fn linicrypt_to_lines(p: &AlgebraicRepresentation<1, 5, 2>) -> Vec<String> {
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

pub fn print_linicrypt(p: &AlgebraicRepresentation<1, 5, 2>) {
    let lines = linicrypt_to_lines(p);
    for line in lines {
        println!("{line}");
    }
}

pub fn linicrypt_to_lines_infos(p: &AlgebraicRepresentation<1, 5, 2>) -> Vec<String> {
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

fn check_css<const DIFF: usize>(
    p: &AlgebraicRepresentation<1, 5, 2>,
    css: &[CollisionStructure<2, DIFF>],
    counter: &mut HashMap<String, usize>,
) {
    for cs in css {
        let cs_id = cs.id();
        if p.has_cs(cs) {
            println!("Y{cs_id}");
            *counter.entry(cs_id).or_insert(0) += 1;
        } else {
            println!(" {cs_id}");
        }
    }
}

fn collision_structure_examples() {
    let programs = generate_3_2_1_programs();

    let css2: Vec<_> = generate_all_cs_2::<2>().collect();
    let css1: Vec<_> = generate_all_cs_2::<1>().collect();

    let mut counter: HashMap<String, usize> = HashMap::new();

    for p in &programs[0..10000] {
        print_linicrypt(p);
        check_css(p, &css2, &mut counter);
        check_css(p, &css1, &mut counter);
    }

    // print_grid(&programs, linicrypt_to_lines_infos, 7);
    println!("{:?}", counter);
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
}
