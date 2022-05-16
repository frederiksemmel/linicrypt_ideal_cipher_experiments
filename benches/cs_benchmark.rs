use criterion::{black_box, criterion_group, criterion_main, Criterion};
use linicrypt::{AlgebraicRepresentation, CollisionStructure, Direction, Operation};

pub fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("single cs check", |b| {
        use Direction::*;
        use Operation::*;
        let p = AlgebraicRepresentation::new(
            [0, 0, 0, 0, 1],
            [
                (E, [1, 0, 1, 0, 0], [1, 0, 0, 0, 0], [0, 0, 0, 1, 0]),
                (E, [0, 1, 0, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]),
            ],
        );
        let cs = CollisionStructure {
            permutation: [1, 0],
            cs_type: [B, B],
        };
        b.iter(|| p.has_cs(black_box(&cs)))
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
