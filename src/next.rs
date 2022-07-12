use na::*;
use nalgebra as na;

pub trait Constraint<const BASE: usize> {}
pub trait DynConstraint {}

type Matrix<const OUT: usize, const BASE: usize> = SMatrix<u8, OUT, BASE>;
type DynMatrix = DMatrix<u8>;

// With compile time checked arrays
#[derive(Debug)]
pub struct Linicrypt<
    const BASE: usize,
    const N: usize,
    const IN: usize,
    const OUT: usize,
    C: Constraint<BASE>,
> {
    pub i: Matrix<IN, BASE>,
    pub m: Matrix<OUT, BASE>,
    pub constraints: [C; N],
}

impl<const BASE: usize, const N: usize, const IN: usize, const OUT: usize, C>
    Linicrypt<BASE, N, IN, OUT, C>
where
    C: Constraint<BASE>,
{
    fn try_invert(self) -> Result<Linicrypt<BASE, N, OUT, IN, C>, ()> {
        Err(())
    }
}

// compile time with macro?

// Nope, not allowed by syntax

// Everything is dynamic
#[derive(Debug)]
pub struct DynLinicrypt<C: DynConstraint> {
    pub i: DynMatrix,
    pub m: DynMatrix,
    pub constraints: Vec<C>,
}

impl<C: DynConstraint> DynLinicrypt<C> {}

// Somethings are dynamic
#[derive(Debug)]
pub struct Dyn2Linicrypt<const BASE: usize, C: Constraint<BASE>> {
    pub i: DynMatrix,
    pub m: DynMatrix,
    pub constraints: Vec<C>,
}

impl<const BASE: usize, C: Constraint<BASE>> Dyn2Linicrypt<BASE, C> {}

// This is an experiment to see if it is possible to clean up the const generics stuff.
// It seems like it is not ready yet...
#[derive(Debug)]
pub struct PLinicrypt<P: Parameters, C: Constraint<{ P::BASE }>>
where
    [(); P::OUT]:,
    [(); P::IN]:,
    [(); P::QUERIES]:,
{
    pub i: Matrix<{ P::IN }, { P::BASE }>,
    pub m: Matrix<{ P::OUT }, { P::BASE }>,
    pub constraints: [C; P::QUERIES],
}

pub trait Parameters {
    const BASE: usize;
    const QUERIES: usize;
    const IN: usize;
    const OUT: usize;
}

impl<P: Parameters, C> PLinicrypt<P, C>
where
    C: Constraint<{ P::BASE }>,
    [(); P::OUT]:,
    [(); P::IN]:,
    [(); P::QUERIES]:,
{
}
