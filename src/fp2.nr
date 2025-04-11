use bignum::{BigNum, BigNumTrait};
use bignum::fields::bls12_381Fq::BLS12_381_Fq_Params;
use std::ops::{Add, Mul, Neg, Sub};

// Parameter reference https://github.com/noir-lang/noir_bigcurve/blob/main/src/curves/bls12_381.nr#L60
pub type Fp = BigNum<4, 381, BLS12_381_Fq_Params>;

// Code reference https://github.com/zkcrypto/bls12_381/blob/main/src/fp2.rs

// c0 + c1*u
// F_p / u^2 + 1
pub struct Fp2 {
    pub c0: Fp,
    pub c1: Fp,
}

impl Fp2 {
    pub fn new(c0: [u128; 4], c1: [u128; 4]) -> Self {
        Fp2 { c0: BigNum { limbs: c0 }, c1: BigNum { limbs: c1 } }
    }

    pub fn zero() -> Fp2 {
        Fp2 { c0: BigNum::zero(), c1: BigNum::zero() }
    }

    pub fn one() -> Fp2 {
        Fp2 { c0: BigNum::one(), c1: BigNum::zero() }
    }

    fn add(self: Self, other: Self) -> Self {
        Self { c0: self.c0 + other.c0, c1: self.c1 + other.c1 }
    }

    fn sub(self: Self, other: Self) -> Self {
        Self { c0: self.c0 - other.c0, c1: self.c1 - other.c1 }
    }

    fn neg(self: Self) -> Self {
        Self { c0: self.c0.neg(), c1: self.c1.neg() }
    }

    pub fn eq(self: Self, other: Self) -> bool {
        (self.c0 == other.c0) & (self.c1 == other.c1)
    }

    pub fn is_zero(self: Self) -> bool {
        (self.c0 == BigNum::zero()) & (self.c1 == BigNum::zero())
    }

    fn mul(self: Self, other: Self) -> Self {
        // Safety: validated using `evaluate_quadratic_expression`
        unsafe {
            // Calculate unconstrained
            // self.c0 * other.c0
            let temp = self.c0.__mul(other.c0);
            // self.c1 * other.c1
            let temp1 = self.c1.__mul(other.c1);
            // c0 = self.c0 * other.c0 - self.c1 * other.c1
            let c0 = temp.__sub(temp1);

            // (self.c0 + self.c1)
            let temp2 = self.c0.__add(self.c1);
            // (other.c0 + other.c1)
            let temp3 = other.c0.__add(other.c1);
            let temp4 = temp2.__mul(temp3);
            // self.c0 * other.c0
            let temp5 = self.c0.__mul(other.c0);
            // self.c1 * other.c1
            let temp6 = self.c1.__mul(other.c1);
            // c1 = (self.c0 + self.c1) * (other.c0 + other.c1) - self.c0 * other.c0 - self.c1 * other.c1
            let c1 = temp4.__sub(temp5).__sub(temp6);

            // Constrain c0
            // c0 = self.c0 * other.c0 - self.c1 * other.c1
            let lhs_terms = [[self.c0], [self.c1]];
            let lhs_flags = [[false], [true]];
            let rhs_terms = [[other.c0], [other.c1]];
            let rhs_flags = [[false], [false]];
            let res = [c0];
            let res_flag = [true];
            BigNum::evaluate_quadratic_expression(
                lhs_terms,
                lhs_flags,
                rhs_terms,
                rhs_flags,
                res,
                res_flag,
            );

            // Constrain c1
            // c1 = (self.c0 + self.c1) * (other.c0 + other.c1) - self.c0 * other.c0 - self.c1 * other.c1
            // (adding zeroes for compilation)
            let lhs_terms1 =
                [[self.c0, self.c1], [self.c0, BigNum::zero()], [self.c1, BigNum::zero()]];
            let lhs_flags1 = [[false, false], [true, false], [true, false]];
            let rhs_terms1 =
                [[other.c0, other.c1], [other.c0, BigNum::zero()], [other.c1, BigNum::zero()]];
            let rhs_flags1 = [[false, false], [false, false], [false, false]];
            let res1 = [c1];
            let res_flag1 = [true];
            BigNum::evaluate_quadratic_expression(
                lhs_terms1,
                lhs_flags1,
                rhs_terms1,
                rhs_flags1,
                res1,
                res_flag1,
            );

            // c0 = self.c0 * other.c0 - self.c1 * other.c1
            // c1 = (self.c0 + self.c1) * (other.c0 + other.c1) - self.c0 * other.c0 - self.c1 * other.c1
            Self { c0, c1 }
        }
    }

    pub fn square(self: Self) -> Self {
        // Safety: validated using `evaluate_quadratic_expression`
        unsafe {
            // a = c0 + c1
            // b = c0 - c1
            // c = c0 + c0
            let a = self.c0.__add(self.c1);
            let b = self.c0.__sub(self.c1);
            let c = self.c0.__add(self.c0);

            // res_c0 = (c0 + c1) * (c0 - c1)
            let res_c0 = a.__mul(b);
            // res_c1 = (c0 + c0) * c1
            let res_c1 = c.__mul(self.c1);

            // Constrain res_c0
            let lhs_terms = [[self.c0, self.c1]];
            let lhs_flags = [[false, false]];
            let rhs_terms = [[self.c0, self.c1]];
            let rhs_flags = [[false, true]];
            let res = [res_c0];
            let res_flag = [true];
            BigNum::evaluate_quadratic_expression(
                lhs_terms,
                lhs_flags,
                rhs_terms,
                rhs_flags,
                res,
                res_flag,
            );

            // Constrain res_c1
            let lhs_terms1 = [[self.c0, self.c0]];
            let lhs_flags1 = [[false, false]];
            let rhs_terms1 = [[self.c1, BigNum::zero()]];
            let rhs_flags1 = [[false, false]];
            let res1 = [res_c1];
            let res_flag1 = [true];
            BigNum::evaluate_quadratic_expression(
                lhs_terms1,
                lhs_flags1,
                rhs_terms1,
                rhs_flags1,
                res1,
                res_flag1,
            );

            Self { c0: res_c0, c1: res_c1 }
        }
    }

    // Returns inverse and errors if there is no inverse
    pub fn invert(self: Self) -> Self {
        // Safety: validated using `evaluate_quadratic_expression`
        unsafe {
            // (self.c0 * self.c0 + self.c1 * self.c1).invert()
            let c0_square = self.c0.__mul(self.c0);
            let c1_square = self.c1.__mul(self.c1);
            let temp = (c0_square.__add(c1_square));

            // Constrain temp = self.c0 * self.c0 + self.c1 * self.c1
            let lhs_terms = [[self.c0], [self.c1]];
            let lhs_flags = [[false], [false]];
            let rhs_terms = [[self.c0], [self.c1]];
            let rhs_flags = [[false], [false]];
            let res = [temp];
            let res_flag = [true];
            BigNum::evaluate_quadratic_expression(
                lhs_terms,
                lhs_flags,
                rhs_terms,
                rhs_flags,
                res,
                res_flag,
            );

            // __invmod asserts that multiplicative inverse exists according to parameters
            // in this case that should pass
            let inv = temp.__invmod();
            // Constrain inverse
            // TODO can this be done more efficiently?
            let check_identity = inv * temp;
            assert(check_identity == BigNum::one());

            // Computing c0, c1 constrained.
            // TODO how can this be improved?
            // c0 = self.c0 * inv
            // c1 = self.c1 * (-inv)
            let res_c0 = self.c0 * inv;
            let res_c1 = self.c1 * inv.neg();

            Self { c0: res_c0, c1: res_c1 }
        }
    }

    // (u+1) is a quadratric residue in Fp2
    // returns self * (u+1)
    pub fn mul_by_nonresidue(self: Self) -> Self {
        // Multiply c0 + c1*u by u + 1:
        // (c0 + c1*u)(u + 1)
        // = c0 * u + c0 + c1 * u^2 + c1 * u
        // and because u^2 = -1, we get
        // = c0 * u + c0 - c1 + c1 * u
        // (c0 - c1) + (c0 + c1)u
        Self { c0: self.c0 - self.c1, c1: self.c0 + self.c1 }
    }

    pub fn conjugate(self: Self) -> Self {
        Fp2 { c0: self.c0, c1: self.c1.neg() }
    }

    // returns self^p
    pub fn frobenius_map(self: Self) -> Self {
        // (c0 + c1*u)^p = c0^p + (c1*u)^p
        // and c0^p=c0, c1^p=c1, and u^p = -u, so the result is:
        // c0 - c1
        self.conjugate()
    }
}

impl Add for Fp2 {
    fn add(self: Self, other: Self) -> Self {
        self.add(other)
    }
}

impl Sub for Fp2 {
    fn sub(self: Self, other: Self) -> Self {
        self.sub(other)
    }
}

impl Mul for Fp2 {
    fn mul(self: Self, other: Self) -> Self {
        self.mul(other)
    }
}

impl Neg for Fp2 {
    fn neg(self: Self) -> Self {
        self.neg()
    }
}

// #[test]
fn test_mul_1() {
    let x = [
        0x55e83ff97a1aeffb3af00adb22c6bb,
        0x8c4f9774b905a14e3a3f171bac586c,
        0xa73197d7942695638c4fa9ac0fc368,
        0x17f1d3,
    ];
    let y = [0, 0, 0, 0];

    let zero = Fp2::zero();
    let test_val = Fp2::new(x, y);

    // test multiplication by 1
    let mult_by_1 = test_val.mul(Fp2::one());

    assert(mult_by_1.c0 == test_val.c0);
    assert(mult_by_1.c1 == test_val.c1);

    // test multiplication by 0
    let mult_by_0 = test_val.mul(zero);
    assert(mult_by_0.c0 == BigNum::zero());
    assert(mult_by_0.c1 == BigNum::zero());
}

#[test]
fn test_addition() {
    let a = Fp2::new(
        [
            0xae4a11560fd2b96a09231943c6dab,
            0xd18b53078ed5d3f5cadea0c65f2734,
            0xff3ed67578ddeae5439da2b02eab31,
            0x8ccec,
        ],
        [
            0x38f87f1c9644538b4bf7ec8dee1be9,
            0x5fb2ff8fddcf6fbb096adcc2b04bfc,
            0x56170fab6f325e416abe497d9e055d,
            0xbd355,
        ],
    );
    let b = Fp2::new(
        [
            0x4579a52cfa8b4f01e3fc00c5962be5,
            0x45b3857ffd05157321b22d5a04c013,
            0xb7e3b626ce947b122db3425e68ae97,
            0xfc9d8,
        ],
        [
            0x8fd0966dcadba0b5f86e4c5ffa1da2,
            0x18864f9cef0cf179640918b4d31a73,
            0xf2af61549b350ae4025242be1008a6,
            0x325f3,
        ],
    );
    let c = Fp2::new(
        [
            0x505e46425b887a98848e3259d29990,
            0x173ed8878bdae968ec90ce2063e747,
            0xb7228c9c477265f77150e50e9759c9,
            0x1896c5,
        ],
        [
            0xc8c9158a611ff441446638ede8398b,
            0x78394f2cccdc61346d73f57783666f,
            0x48c671000a6769256d108c3bae0e03,
            0xef949,
        ],
    );
    let res = a.add(b);
    assert(res.c0 == c.c0);
    assert(res.c1 == c.c1);
}

#[test]
fn test_subtraction() {
    let a = Fp2::new(
        [
            0xae4a11560fd2b96a09231943c6dab,
            0xd18b53078ed5d3f5cadea0c65f2734,
            0xff3ed67578ddeae5439da2b02eab31,
            0x8ccec,
        ],
        [
            0x38f87f1c9644538b4bf7ec8dee1be9,
            0x5fb2ff8fddcf6fbb096adcc2b04bfc,
            0x56170fab6f325e416abe497d9e055d,
            0xbd355,
        ],
    );
    let b = Fp2::new(
        [
            0x4579a52cfa8b4f01e3fc00c5962be5,
            0x45b3857ffd05157321b22d5a04c013,
            0xb7e3b626ce947b122db3425e68ae97,
            0xfc9d8,
        ],
        [
            0x8fd0966dcadba0b5f86e4c5ffa1da2,
            0x18864f9cef0cf179640918b4d31a73,
            0xf2af61549b350ae4025242be1008a6,
            0x325f3,
        ],
    );
    let c = Fp2::new(
        [
            0x716afa99ba71dc4ebb9630cea5ec71,
            0xd75cc10ca49025b37bcd6a1d508b3f,
            0x3194a03544948b7acc2dabfe9d6111,
            0x130426,
        ],
        [
            0xa927e8aecb68b2d55389a02df3fe47,
            0x472caff2eec27e41a561c40ddd3188,
            0x6367ae56d3fd535d686c06bf8dfcb7,
            0x8ad61,
        ],
    );
    let res = a.sub(b);
    assert(res.c0 == c.c0);
    assert(res.c1 == c.c1);
}

#[test]
fn test_multiplication() {
    // zkcrypto test
    // https://github.com/zkcrypto/bls12_381/blob/main/src/fp2.rs#L465
    let a = Fp2::new(
        [
            0xae4a11560fd2b96a09231943c6dab,
            0xd18b53078ed5d3f5cadea0c65f2734,
            0xff3ed67578ddeae5439da2b02eab31,
            0x8ccec,
        ],
        [
            0x38f87f1c9644538b4bf7ec8dee1be9,
            0x5fb2ff8fddcf6fbb096adcc2b04bfc,
            0x56170fab6f325e416abe497d9e055d,
            0xbd355,
        ],
    );
    let b = Fp2::new(
        [
            0x4579a52cfa8b4f01e3fc00c5962be5,
            0x45b3857ffd05157321b22d5a04c013,
            0xb7e3b626ce947b122db3425e68ae97,
            0xfc9d8,
        ],
        [
            0x8fd0966dcadba0b5f86e4c5ffa1da2,
            0x18864f9cef0cf179640918b4d31a73,
            0xf2af61549b350ae4025242be1008a6,
            0x325f3,
        ],
    );
    let c = Fp2::new(
        [
            0x77b92dc10d7dc4df06bbf3f0fe3a47,
            0x68be63dff75f16e1ac619556ce0f96,
            0xeade1ad5acb2d25bcc90c425869481,
            0x80690,
        ],
        [
            0xca7635cc60dd923e8c7db21115127e,
            0x66780a57f3ae3477f44f30c19df036,
            0xb168c449a974a86ef04b08f869126,
            0x6e753,
        ],
    );

    let res = a.mul(b);
    assert(res.c0 == c.c0);
    assert(res.c1 == c.c1);
}

#[test]
fn test_square() {
    // https://github.com/zkcrypto/bls12_381/blob/main/src/fp2.rs#L423
    let a = Fp2::new(
        [
            0xae4a11560fd2b96a09231943c6dab,
            0xd18b53078ed5d3f5cadea0c65f2734,
            0xff3ed67578ddeae5439da2b02eab31,
            0x8ccec,
        ],
        [
            0x38f87f1c9644538b4bf7ec8dee1be9,
            0x5fb2ff8fddcf6fbb096adcc2b04bfc,
            0x56170fab6f325e416abe497d9e055d,
            0xbd355,
        ],
    );
    let c = Fp2::new(
        [
            0x4579a52cfa8b4f01e3fc00c5962be5,
            0x45b3857ffd05157321b22d5a04c013,
            0xb7e3b626ce947b122db3425e68ae97,
            0xfc9d8,
        ],
        [
            0x8fd0966dcadba0b5f86e4c5ffa1da2,
            0x18864f9cef0cf179640918b4d31a73,
            0xf2af61549b350ae4025242be1008a6,
            0x325f3,
        ],
    );
    let res = a.square();
    assert(res.c0 == c.c0);
    assert(res.c1 == c.c1);
}

#[test]
fn test_inversion() {
    let a = Fp2::new(
        [
            0xb77accf6a4e5e337a85d2a62baf070,
            0x5e5a337d7db6344484087a18586c7c,
            0xb3f92a8e29eabd22358787933cd419,
            0xa4aa,
        ],
        [
            0xb33c802bede6916fa9f3d1d08a643b,
            0x4537aef2966b539d7364999eb35707,
            0xde1e413ba79349fbe6cc88298b9f18,
            0x1776be,
        ],
    );
    let c = Fp2::new(
        [
            0x71825b46349b307298025eea77770c,
            0x225cbdc399bb960534036f625a3592,
            0xb32ff3d468eb1aba35968a2bbdb9d0,
            0x66e47,
        ],
        [
            0x3fc9a610ece44f3c738a358aa60150,
            0x602cf426b135e0e0d837a8cf5d9ec7,
            0x4aee9cb5468ad40ee5dccbf0264407,
            0xdc722,
        ],
    );

    let res = a.invert();
    assert(res.c0 == c.c0);
    assert(res.c1 == c.c1);

    let identity_check = res.mul(a);
    assert(identity_check.eq(Fp2::one()));
}
