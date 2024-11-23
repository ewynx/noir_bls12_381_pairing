use crate::fp2::{Fp2, Fp};
use crate::fp6::Fp6;
use crate::fp12::Fp12;
use bignum::BigNum;

// Ref: https://github.com/zkcrypto/bls12_381/blob/main/src/pairings.rs

pub struct G2Projective {
    x: Fp2,
    y: Fp2,
    z: Fp2,
}

fn main(g: G1Affine, h: G2Affine) {
    pairing(g, h);
}

impl G2Projective {
    fn from_affine(p: G2Affine) -> Self {
        let mut z = Fp2::one();
        if (p.infinity) {
            z = Fp2::zero();
        }
        G2Projective { x: p.x, y: p.y, z }
    }
}

pub struct G2Affine {
    x: Fp2,
    y: Fp2,
    infinity: bool,
}

pub struct G1Affine {
    x: Fp,
    y: Fp,
    infinity: bool,
}

impl G2Affine {
    pub fn neg(self: Self) -> Self {
        let mut y = -self.y;
        if self.infinity {
            y = Fp2::one();
        }
        G2Affine { x: self.x, y, infinity: self.infinity }
    }

    pub fn generator() -> Self {
        G2Affine {
            x: Fp2::new(
                [
                    0xac0326a805bbefd48056c8c121bdb8,
                    0x7ad4fa403b02b4510b647ae3d1770b,
                    0xb2f08f0a91260805272dc51051c6e4,
                    0x24aa2,
                ],
                [
                    0x4cf11213945d57e5ac7d055d042b7e,
                    0xd0d09920b61ab5da61bbdc7f504933,
                    0x6052719f607dacd3a088274f65596b,
                    0x13e02b,
                ],
            ),
            y: Fp2::new(
                [
                    0x3ac9cc3baca289e193548608b82801,
                    0x9baa8cbdd3a76d429a695160d12c92,
                    0x27727d6e118cc9cdc6da2e351aadfd,
                    0xce5d5,
                ],
                [
                    0x370d275cec1da1aaa9075ff05f79be,
                    0x287e85a763af267492ab572e99ab3f,
                    0xa02ea734cc32acd2b02bc28b99cb3e,
                    0x606c4,
                ],
            ),
            infinity: false,
        }
    }
}

impl G1Affine {
    pub fn generator() -> Self {
        Self {
            x: BigNum {
                limbs: [
                    0x55e83ff97a1aeffb3af00adb22c6bb,
                    0x8c4f9774b905a14e3a3f171bac586c,
                    0xa73197d7942695638c4fa9ac0fc368,
                    0x17f1d3,
                ],
            },
            y: BigNum {
                limbs: [
                    0x3cc744a2888ae40caa232946c5e7e1,
                    0xe095d5d00af600db18cb2c04b3edd0,
                    0x81e3aaa0f1a09e30ed741d8ae4fcf5,
                    0x8b3f4,
                ],
            },
            infinity: false,
        }
    }

    fn neg(self: Self) -> Self {
        let mut y = self.y.neg();
        if self.infinity {
            y = BigNum::one();
        }
        Self { x: self.x, y, infinity: self.infinity }
    }
}

// https://github.com/zkcrypto/bls12_381/blob/main/src/pairings.rs#L709
pub fn doubling_step(r: G2Projective) -> (Fp2, Fp2, Fp2, G2Projective) {
    let mut tmp0 = r.x.square();
    let mut tmp1 = r.y.square();
    let mut tmp2 = tmp1.square();
    let mut tmp3 = (tmp1 + r.x).square() - tmp0 - tmp2;
    tmp3 = tmp3 + tmp3;
    let tmp4 = tmp0 + tmp0 + tmp0;
    let mut tmp6 = r.x + tmp4;
    let tmp5 = tmp4.square();
    let zsquared = r.z.square();
    let new_r_x = tmp5 - tmp3 - tmp3;
    let new_r_z = (r.z + r.y).square() - tmp1 - zsquared;
    let mut new_r_y = (tmp3 - new_r_x) * tmp4;
    tmp2 = tmp2 + tmp2;
    tmp2 = tmp2 + tmp2;
    tmp2 = tmp2 + tmp2;
    new_r_y = new_r_y - tmp2;
    tmp3 = tmp4 * zsquared;
    tmp3 = tmp3 + tmp3;
    tmp3 = -tmp3;
    tmp6 = (tmp6.square()) - tmp0 - tmp5;
    tmp1 = tmp1 + tmp1;
    tmp1 = tmp1 + tmp1;
    tmp6 = tmp6 - tmp1;
    tmp0 = new_r_z * zsquared;
    tmp0 = tmp0 + tmp0;

    (tmp0, tmp3, tmp6, G2Projective { x: new_r_x, y: new_r_y, z: new_r_z })
}

pub fn addition_step(r: G2Projective, q: G2Affine) -> (Fp2, Fp2, Fp2, G2Projective) {
    let zsquared = r.z.square();
    let ysquared = q.y.square();
    let mut t0 = zsquared * q.x;
    let mut t1 = ((q.y + r.z).square() - ysquared - zsquared) * zsquared;
    let t2 = t0 - r.x;
    let t3 = t2.square();

    let mut t4 = t3 + t3;
    t4 = t4 + t4;
    let t5 = t4 * t2;
    let mut t6 = t1 - r.y - r.y;
    let mut t9 = t6 * q.x;
    let t7 = t4 * r.x;

    let new_r_x = t6.square() - t5 - t7 - t7;

    let new_r_z = (r.z + t2).square() - zsquared - t3;

    let mut t10 = q.y + new_r_z;
    let t8 = (t7 - new_r_x) * t6;
    t0 = r.y * t5;
    t0 = t0 + t0;
    let new_r_y = t8 - t0;

    t10 = t10.square() - ysquared;
    let ztsquared = new_r_z.square();
    t10 = t10 - ztsquared;
    t9 = t9 + t9 - t10;

    let t10_double = new_r_z + new_r_z;
    let t6_neg = -t6;
    let t1_final = t6_neg + t6_neg;
    (t10_double, t1_final, t9, G2Projective { x: new_r_x, y: new_r_y, z: new_r_z })
}

pub fn ell(f: Fp12, coeffs: (Fp2, Fp2, Fp2), p: G1Affine) -> Fp12 {
    let c0c0 = coeffs.0.c0 * p.y;
    let c0c1 = coeffs.0.c1 * p.y;

    let c1c0 = coeffs.1.c0 * p.x;
    let c1c1 = coeffs.1.c1 * p.x;

    f.mul_by_014(coeffs.2, Fp2 { c0: c1c0, c1: c1c1 }, Fp2 { c0: c0c0, c1: c0c1 })
}

pub fn miller_loop(p: G1Affine, q: G2Affine) -> Fp12 {
    // https://github.com/zkcrypto/bls12_381/blob/main/src/pairings.rs#L672
    // Got array from zkcrypto impl
    let booleans = [
        false, true, true, false, true, false, false, true, false, false, false, false, false,
        false, false, false, true, false, false, false, false, false, false, false, false, false,
        false, false, false, false, false, false, false, false, false, false, false, false, false,
        false, false, false, false, false, false, false, false, false, true, false, false, false,
        false, false, false, false, false, false, false, false, false, false, false, false,
    ];

    let mut adder_cur = G2Projective::from_affine(q);
    let mut adder_base = q;
    let mut adder_p = p;
    let mut f = Fp12::one();
    let mut found_one = false;

    for i in 0..64 {
        let i_bool = booleans[i];
        if !found_one {
            found_one = i_bool;
        } else {
            // Doubling step receives f. Does 2 things: doubling_step function & ell function
            let (coeffs_c0, coeffs_c1, coeffs_c2, cur_updated) = doubling_step(adder_cur);
            adder_cur = cur_updated;
            // ell updates f
            f = ell(f, (coeffs_c0, coeffs_c1, coeffs_c2), adder_p);

            if i_bool {
                // addition_step receives f & does 2 things: addition_step & ell function
                let (coeffs_c0, coeffs_c1, coeffs_c2, cur_updated) =
                    addition_step(adder_cur, adder_base);
                adder_cur = cur_updated;
                // ell updates f
                f = ell(f, (coeffs_c0, coeffs_c1, coeffs_c2), adder_p);
            }
            f = f.square();
        }
    }
    // Doubling step receives f. Does 2 things: doubling_step function & ell function
    let (coeffs_c0, coeffs_c1, coeffs_c2, cur_updated) = doubling_step(adder_cur);
    adder_cur = cur_updated;
    // ell updates f
    f = ell(f, (coeffs_c0, coeffs_c1, coeffs_c2), adder_p);

    f = f.conjugate();
    f
}

// https://github.com/zkcrypto/bls12_381/blob/main/src/pairings.rs#L50
pub fn fp4_square(a: Fp2, b: Fp2) -> (Fp2, Fp2) {
    let t0 = a.square();
    let t1 = b.square();
    let mut t2 = t1.mul_by_nonresidue();
    let c0 = t2 + t0;
    t2 = a + b;
    t2 = t2.square();
    t2 = t2 - t0;
    let c1 = t2 - t1;

    (c0, c1)
}

// https://github.com/zkcrypto/bls12_381/blob/main/src/pairings.rs#L66
pub fn cyclotomic_square(f: Fp12) -> Fp12 {
    let mut z0 = f.c0.c0;
    let mut z4 = f.c0.c1;
    let mut z3 = f.c0.c2;
    let mut z2 = f.c1.c0;
    let mut z1 = f.c1.c1;
    let mut z5 = f.c1.c2;

    let (t0, t1) = fp4_square(z0, z1);

    z0 = t0 - z0;
    z0 = z0 + z0 + t0;
    z1 = t1 + z1;
    z1 = z1 + z1 + t1;

    let (mut t0, t1) = fp4_square(z2, z3);
    let (t2, t3) = fp4_square(z4, z5);

    z4 = t0 - z4;
    z4 = z4 + z4 + t0;

    z5 = t1 + z5;
    z5 = z5 + z5 + t1;

    t0 = t3.mul_by_nonresidue();
    z2 = t0 + z2;
    z2 = z2 + z2 + t0;

    z3 = t2 - z3;
    z3 = z3 + z3 + t2;

    Fp12 { c0: Fp6 { c0: z0, c1: z4, c2: z3 }, c1: Fp6 { c0: z2, c1: z1, c2: z5 } }
}

pub fn cycolotomic_exp(f: Fp12) -> Fp12 {
    // https://github.com/zkcrypto/bls12_381/blob/main/src/pairings.rs#L119
    // Got array from zkcrypto impl
    let booleans = [
        true, true, false, true, false, false, true, false, false, false, false, false, false,
        false, false, true, false, false, false, false, false, false, false, false, false, false,
        false, false, false, false, false, false, false, false, false, false, false, false, false,
        false, false, false, false, false, false, false, false, true, false, false, false, false,
        false, false, false, false, false, false, false, false, false, false, false, false,
    ];

    let mut tmp = Fp12::one();
    let mut found_one = false;
    for i in 0..64 {
        if found_one {
            tmp = cyclotomic_square(tmp);
        } else {
            found_one = booleans[i];
        }

        if booleans[i] {
            tmp = tmp * f;
        }
    }

    tmp.conjugate()
}

impl Fp12 {
    pub fn final_exponentiation(self: Self) -> Self {
        let mut f = self;
        let mut t0 = f
            .frobenius_map()
            .frobenius_map()
            .frobenius_map()
            .frobenius_map()
            .frobenius_map()
            .frobenius_map();
        // This will error if inverse doesn't exist
        let mut t1 = f.invert();
        let mut t2 = t0 * t1;
        t1 = t2;
        t2 = t2.frobenius_map().frobenius_map();
        t2 = t2 * t1;
        t1 = cyclotomic_square(t2).conjugate();
        let mut t3 = cycolotomic_exp(t2);
        let mut t4 = cyclotomic_square(t3);
        let mut t5 = t1 * t3;
        t1 = cycolotomic_exp(t5);
        t0 = cycolotomic_exp(t1);
        let mut t6 = cycolotomic_exp(t0);
        t6 = t6 * t4;
        t4 = cycolotomic_exp(t6);
        t5 = t5.conjugate();
        t4 = t4 * t5 * t2;
        t5 = t2.conjugate();
        t1 = t1 * t2;
        t1 = t1.frobenius_map().frobenius_map().frobenius_map();
        t6 = t6 * t5;
        t6 = t6.frobenius_map();
        t3 = t3 * t0;
        t3 = t3.frobenius_map().frobenius_map();
        t3 = t3 * t1;
        t3 = t3 * t6;
        f = t3 * t4;
        f
    }
}

pub fn pairing(p: G1Affine, q: G2Affine) -> Fp12 {
    // If any of the inputs is infinity, the result is 1
    if p.infinity | q.infinity {
        Fp12::one()
    } else {
        let tmp = miller_loop(p, q);
        tmp.final_exponentiation()
    }
}

// #[test]
fn test_doubling_step() {
    let generator = G2Projective {
        x: Fp2::new(
            [
                0xac0326a805bbefd48056c8c121bdb8,
                0x7ad4fa403b02b4510b647ae3d1770b,
                0xb2f08f0a91260805272dc51051c6e4,
                0x24aa2,
            ],
            [
                0x4cf11213945d57e5ac7d055d042b7e,
                0xd0d09920b61ab5da61bbdc7f504933,
                0x6052719f607dacd3a088274f65596b,
                0x13e02b,
            ],
        ),
        y: Fp2::new(
            [
                0x3ac9cc3baca289e193548608b82801,
                0x9baa8cbdd3a76d429a695160d12c92,
                0x27727d6e118cc9cdc6da2e351aadfd,
                0xce5d5,
            ],
            [
                0x370d275cec1da1aaa9075ff05f79be,
                0x287e85a763af267492ab572e99ab3f,
                0xa02ea734cc32acd2b02bc28b99cb3e,
                0x606c4,
            ],
        ),
        z: Fp2::one(),
    };

    let a = Fp2::new(
        [
            0x3f27323d5e8a27cc4e521822e0f559,
            0x23253f723bde4dd997044ed24e8e2a,
            0xb39075d1abe80b8f65256d2793537f,
            0x199642,
        ],
        [
            0xdc349d73b07686aaa41d7fc17de6f8,
            0xa1fa169d8ebc99d24aad5cba66acfc,
            0x80ba9cd330cab34ac0af0a2e672cf8,
            0x181b12,
        ],
    );
    let b = Fp2::new(
        [
            0xb8f2795d21163ce39975990c57da6d,
            0x1f7a219c6ac2d020fff505604caf82,
            0x2c22ea6825bb8178176c3a6ed7df18,
            0x291a,
        ],
        [
            0xe781898edb7994dd769c7c2c6a8447,
            0x690e2cdf9728c6c7a4c90f9ad207a7,
            0x4f2c59babfdd1438e64394468e1eaf,
            0x8a764,
        ],
    );
    let c = Fp2::new(
        [
            0x3cd6e1c4403fdca07f3855ab48b7d,
            0xe23363105a3956fc734cd655c20465,
            0x4989c84ac56fcf81ef6ef0c51408b3,
            0x116916,
        ],
        [
            0x2136e9d5f30d406e4145ecc62c1f99,
            0xed1365d17357536a2d488e2a528552,
            0xac0820ccbd5704388bf7fc5530c4d9,
            0x1b696,
        ],
    );
    let doubled = G2Projective {
        x: Fp2::new(
            [
                0xd403df5b1ac137922b915f1f38ec37,
                0xf87ae23fe110c5aa48bc208fbd1a0e,
                0xdd65245f8908a03288adcd396f489e,
                0x177d52,
            ],
            [
                0x8e26f1e5febc57b750d1c3dc4261d8,
                0x5f169f0d4ddf46e23f15bb4f24304a,
                0x8b9e689553d58194f79863fe02902c,
                0xcf815,
            ],
        ),
        y: Fp2::new(
            [
                0x730c92c30a69d0602a8804c2a7c65f,
                0x9cdb85d9080b2c9b4766434c33d103,
                0x215806e8a55fd2d9ec4af9d2d44859,
                0x65ae9,
            ],
            [
                0x361c973962955ba6405f0e51beb98b,
                0xe41a04569ae33184419e66bbc0271c,
                0x2d8a6d4b3a1cbd02c7bdc0e0aa304d,
                0xe9c34,
            ],
        ),
        z: Fp2::new(
            [
                0x75939877594513c326a90c11705002,
                0x3755197ba74eda8534d2a2c1a25924,
                0x4ee4fadc2319939b8db45c6a355bfb,
                0x19cbaa,
            ],
            [
                0x6e1a4eb9d83b4355520ebfe0bef37c,
                0x50fd0b4ec75e4ce92556ae5d33567e,
                0x405d4e69986559a56057851733967c,
                0xc0d89,
            ],
        ),
    };

    // Result are 3 Fp2's and the adjusted projective point
    // (in Rust this is changed in place, but this is not supported in Noir so it is an additional return value)
    let res = doubling_step(generator);

    // Check changed point
    assert(res.3.x.eq(doubled.x));
    assert(res.3.y.eq(doubled.y));
    assert(res.3.z.eq(doubled.z));

    // Check Fp2's
    assert(res.0.eq(a));
    assert(res.1.eq(b));
    assert(res.2.eq(c));
}

// #[test]
fn test_addition_step() {
    /*
  Values from zkcrypto impl
  #[test]
  fn test_addition_step() {
    let a = Scalar::from_raw([1, 2, 3, 4]).invert().unwrap().square();
    let test_value = G2Affine::generator()*a;
    let mut test_val_proj = G2Projective::from(test_value);
    println!("{}\n", test_val_proj);
    let t = test_value.into();
    println!("{}\n", t);
    let dbl = addition_step(&mut test_val_proj, &t);
    println!("{:?}\n", dbl);
    println!("{}\n", test_val_proj);
  }
  */
    let point1 = G2Projective {
        x: Fp2::new(
            [
                0x335fee64bd06db6573fdcb9ee98319,
                0x593fac9c1b9b2480ae1c7f4e45f9bb,
                0x6da041aadfe26730abaac871996971,
                0x7c7d4,
            ],
            [
                0xf3f398d2eec5613d7e8fb8eea72f72,
                0xdbaafceeb0accc8abf24d95553a45c,
                0x15344fea28e0fd9ab57be98b054672,
                0x2531e,
            ],
        ),
        y: Fp2::new(
            [
                0x69c8bd9c9d71502f1318fe1be21cc0,
                0xc079d0edfc279c0052d38038f28fc5,
                0xadb0d4365ff48c2cd32891ac2615a8,
                0x75e8e,
            ],
            [
                0xc3230bf99a869a5f4a6627fe3851b9,
                0xa833d8fb5f0e5559cd8e03356c9e4d,
                0x889e30f914fe674802ed8803f0f0d7,
                0x590bf,
            ],
        ),
        z: Fp2::new(
            [
                0x91320d77ab2c01f0cf4b79f8644bc8,
                0x375ffeade6cac1c62cf9ab9052de7,
                0x706f3ced1794bc6fd493745e2340a5,
                0xda79,
            ],
            [
                0xc39955cb1638df156cbab123f66dbc,
                0xac805f7513b6a5a8e1e1258d57f404,
                0x8940c37a963a03b135428745f2340c,
                0x18c088,
            ],
        ),
    };

    let affine = G2Affine {
        x: Fp2::new(
            [
                0xb140828914d4659e103008b0f99ce2,
                0xab1484bd130dc52e9a993b56bfa96e,
                0xad2ce8990405d4fafad4edbce84d24,
                0x71593,
            ],
            [
                0x7866ec1320541b8e914140cada9eac,
                0xecc253bae630f295cb97d79685e9ee,
                0x2a956ba1f0c00f7ac188eeda032acb,
                0x6322e,
            ],
        ),
        y: Fp2::new(
            [
                0x688f932215126dacf1612804d3dac3,
                0xd4fb713edca1dbb4cbfa81de69b297,
                0x423f261bb8f9441937551b8c3cc029,
                0xcacc,
            ],
            [
                0x3840064ad1355cd762cc26266d678f,
                0x5ff57854ba0a2e170720dc73b35053,
                0xbff73fc2bd237da821284216f8dc9f,
                0xaab8c,
            ],
        ),
        infinity: false,
    };

    let a = Fp2::new(
        [
            0x1bf08d8f6c261980aa6ec9377eeb5a,
            0x560c420d56d1257cad1c8daa6ffd6a,
            0x15ac001d9750dc97595e05a3ab81a1,
            0x19e61c,
        ],
        [
            0x24330a557692fa7b0c56585c2c0696,
            0x5eb50c857006e9d130313da2ccbbc6,
            0xfaac085ce0a7842b2aa3044ca02fe7,
            0x10048e,
        ],
    );
    let b = Fp2::new(
        [
            0x589dad22b7c810c06691568df9e5c5,
            0xb21119bda48ce10170e15e1d6b435c,
            0xff859681f3e36d0ca9803e70a898f7,
            0x8e196,
        ],
        [
            0xe9c4e582ffd6d9076c4caeefea9722,
            0x2895927e77470516516482f3284d77,
            0x85492b6009b98ac5a4990e3225a290,
            0x80188,
        ],
    );
    let c = Fp2::new(
        [
            0x6be25edc8c984a0727703e81c7369e,
            0x556f550564310a5e0521371de6e410,
            0xeed02b82a7cfb8bbe3547aed326bb7,
            0x8e06c,
        ],
        [
            0xe8e4c2f187c9ba393c2f72314a8bba,
            0x573eb0a205bc6a8eb23577b3f87754,
            0x14d1f9e61982c4efa5fdc5d9ae9afe,
            0x153622,
        ],
    );
    let addition_res_expected = G2Projective {
        x: Fp2::new(
            [
                0xe60e8ea5133cd50e9bbd821aebf8b0,
                0x264be53d3aba569e1414ee1cbdeed6,
                0xf27b9d49940bc78483c46b66300132,
                0x18fe28,
            ],
            [
                0xd85e3b56371e75a3e1ab6f4de4684c,
                0x645f74f93b077a6d50bb261c4fb3f8,
                0x7a0d54da0e220533084ab75e195f31,
                0x78ff4,
            ],
        ),
        y: Fp2::new(
            [
                0xdc24946be6b95b618be1503266315e,
                0x477f1fdb9da707a951559beccf63d3,
                0xdfb35908d8a5d1f5de994c0ab7dd28,
                0x123c72,
            ],
            [
                0x61b1cac14443e5fc724d1f17117692,
                0x3df3f01f65da37a4816ee5269ad530,
                0x8744f8a6346d7242623e875b47e064,
                0x17ab85,
            ],
        ),
        z: Fp2::new(
            [
                0xdf846c7b6130cc05537649bbf75ad,
                0xab062106ab6892be568e46d537feb5,
                0xad6000ecba86e4bacaf02d1d5c0d0,
                0xcf30e,
            ],
            [
                0x1219852abb497d3d862b2c2e16034b,
                0xaf5a8642b80374e898189ed1665de3,
                0x7d56042e7053c215955182265017f3,
                0x80247,
            ],
        ),
    };

    // Result are 3 Fp2's and the adjusted projective point
    // (in Rust this is changed in place, but this is not supported in Noir so it is an additional return value)
    let res: (Fp2, Fp2, Fp2, G2Projective) = addition_step(point1, affine);

    // Check changed point
    assert(res.3.x.eq(addition_res_expected.x));
    assert(res.3.y.eq(addition_res_expected.y));
    assert(res.3.z.eq(addition_res_expected.z));

    // Check Fp2's
    assert(res.0.eq(a));
    assert(res.1.eq(b));
    assert(res.2.eq(c));
}

// #[test]
fn test_final_exponentiation_identity() {
    let t = Fp12::one();
    let res = t.final_exponentiation();
    assert(res.eq(Fp12::one()));
}

// #[test]
fn test_pairing_identity() {
    let p = G1Affine {
        // x, y are irrelevant because infinity is set to true
        x: BigNum::one(),
        y: BigNum::one(),
        infinity: true,
    };

    let expected_identity = pairing(p, G2Affine::generator());
    assert(Fp12::one().eq(expected_identity));
}

// #[test]
fn test_miller_loop_gh() {
    // Testing intermediate step, expected value from zkcrypto impl
    // test values from https://github.com/zkcrypto/bls12_381/blob/main/src/pairings.rs#L835
    let g = G1Affine {
        x: BigNum {
            limbs: [
                0xfabb7be712b08862b6c9f3e9ff6ccc,
                0xdb4476eb5e63dc7314183a309224d,
                0x22fc235f054631617d5f134c6544f0,
                0x845e2,
            ],
        },
        y: BigNum {
            limbs: [
                0xbd5303645f9b22e005f9e2409b286e,
                0xce946b2e0489fbfdff1abb380f4d22,
                0xe2ff06bfa6a8e879aed13e82a17f86,
                0x362e6,
            ],
        },
        infinity: false,
    };

    let h = G2Affine {
        x: Fp2::new(
            [
                0xfb29f4c9b9f7aee559363b9545755c,
                0xbc1390ea1411adbee0eec731b72471,
                0x652b881bca2fe0b80fb53e5c6a8eee,
                0x154358,
            ],
            [
                0x82d763ca4de753d2b3ee304205157,
                0xfe44d89efc569c4bec0ee56c95428f,
                0xac4740bb4f874bdd05a017cb1eb950,
                0x1c656,
            ],
        ),
        y: Fp2::new(
            [
                0x23853393c1100055c5d34c423b3fd2,
                0xf4f22fff37029a585bb1038166e5ce,
                0xe508777c4881dfcf5df11b0425555f,
                0xe3d19,
            ],
            [
                0x1ec46236ee424e72d540ab93e7627a,
                0xf27c49db63114121eaae2ab56dcfad,
                0xd24956daaba822362ec9157039e792,
                0x193caf,
            ],
        ),
        infinity: false,
    };
    let res = miller_loop(g, h);

    let c0 = Fp6 {
        c0: Fp2::new(
            [
                0xeca6fe6fd61761a14d7a4efd04dce8,
                0xae1ce267c3152f809bd0e24cdb3439,
                0x2055adc45db96dbfaad283399f930c,
                0x4b55b,
            ],
            [
                0xfb51edc493eb47d2548b6d3deff2e3,
                0x5f8dc81fe0303d5202e21b27c7cae8,
                0xcd4c031604fe969028f05a22133e4d,
                0x4856c,
            ],
        ),
        c1: Fp2::new(
            [
                0xcea05ef042affb81c007bc1a924e38,
                0x661614d648d4b85aa6ea490e04879e,
                0xf2eb8aef36c8c4aa56791451a956cd,
                0x63e7a,
            ],
            [
                0x658a3de2190a0d68b88af6f2e7f882,
                0x10b4e6591964d7b99fcfd825619448,
                0xab2bf19ab61e92e3e35d137c109f7a,
                0x4c12d,
            ],
        ),
        c2: Fp2::new(
            [
                0x42e14b78b887013477f422fd0b7c53,
                0xfebb9e4384e1d2a5684bf7368e3dcd,
                0xfd5a50a50de7aa9c7f5da287b41cb2,
                0xbe166,
            ],
            [
                0x53a80d17ea07851bfa580d73a5699b,
                0x52445b814699a564bfad2b397a93d3,
                0x475c60246bf1fab412c52f0316c8dd,
                0xa2632,
            ],
        ),
    };

    let c1 = Fp6 {
        c0: Fp2::new(
            [
                0x18568bf1932586e1fc82be25e3e7f,
                0x774bb8fc902e14236184b14831b338,
                0x28a3c9c5d0748adb11a251626eeff4,
                0xcd36f,
            ],
            [
                0xc7d7e9f7df994f4b44b22d9a1a2e7b,
                0xfda5e2c74b5247766e5c23a4adc42,
                0xa10a56360b523bd04edae16f2ddfc9,
                0x19c7db,
            ],
        ),
        c1: Fp2::new(
            [
                0xfbaca33e71d006994e65a453a49d6,
                0xa5f540fa4f2a422f155396fca0f56f,
                0x3f2179791e5d2a8af4375f4a0909ed,
                0x16d48d,
            ],
            [
                0xf27b4ff0c0f8186e52642c37e1399d,
                0xa79376982b5e7582cb7d2408b1581,
                0x7b276530ea982e455da02f5455c66a,
                0x19e623,
            ],
        ),
        c2: Fp2::new(
            [
                0xd57d65cd164deb89f543ba49e5b965,
                0x571046036743b54a03348205caf2a7,
                0xa10301fa6b0ddcbd6bfd20031c35de,
                0x14d6b3,
            ],
            [
                0x6d9fe64f15e04d5f1416e6019dce4f,
                0x5aa69dd12b41d55e2c9b2f62d4d4bf,
                0x12dc81356efd9ef477f8fd46dc9fef,
                0xdb4c1,
            ],
        ),
    };

    let expected_res = Fp12 { c0, c1 };

    assert(res.eq(expected_res));
}

// #[test]
fn test_miller_loop_expected_generator() {
    // Testing intermediate step, expected value from zkcrypto impl
    // test values from https://github.com/zkcrypto/bls12_381/blob/main/src/pairings.rs#L835
    let expected = G1Affine {
        x: BigNum {
            limbs: [
                0xd7db4700715c7fc0a50d36e4d66c19,
                0xab5afe0efc1e6faf7476e6b42bafdc,
                0xf96908d5627e3dd1d3bb6f5476401,
                0x9ee21,
            ],
        },
        y: BigNum {
            limbs: [
                0x6485cf916eac99502622e4752b19f4,
                0x95a16b09abd95b4cbfedb2cbe87f7b,
                0x981c52a89257acb922bf3ffd74a55f,
                0x15878d,
            ],
        },
        infinity: false,
    };

    let res = miller_loop(expected, G2Affine::generator());

    let c0 = Fp6 {
        c0: Fp2::new(
            [
                0x34eeb75b409ba2acf17a3f23b7f67,
                0x3e3856778d39180b6e8509edb2228b,
                0x2d3fea4c1e59809ef2f7df0c003956,
                0xd6a50,
            ],
            [
                0x144d9aa4cc68911bb7c3cb3d953fcf,
                0xc4c6a25d9a086f4b1211e971bbd6d7,
                0x93673ab1a1903a2ece3860a46c2f19,
                0xcc276,
            ],
        ),
        c1: Fp2::new(
            [
                0x42904787d7a0e04698792b05024e06,
                0xedf6ab40f7310f66eba8fd5292d88f,
                0x9642147a44201ad94f7dd01049d74f,
                0xed53,
            ],
            [
                0xb7bed18cc0e883b0f62c7409f075ff,
                0xf6af2d659a89bf3ddd9ad7f6725957,
                0x63f17de251bc34bfccc6934316ab59,
                0x79db6,
            ],
        ),
        c2: Fp2::new(
            [
                0x354121c6f05c68bd3a41bdc5ef2094,
                0x1786120adeb65ec46be891279b345b,
                0x9f4ee0835af6442df8fa981f94d951,
                0x18893d,
            ],
            [
                0x60962d5e5ddc05ce5f5f73d3ebb25b,
                0x79f5056cdb1cae2c4d317286c2dc30,
                0x5f562105010d8b3356904f645be303,
                0x153cd8,
            ],
        ),
    };

    let c1 = Fp6 {
        c0: Fp2::new(
            [
                0x36ece658926e4a9cd0fc4e3c3e47ad,
                0xafbc7e7fecb53133555f636358dbd8,
                0x31d52bdd18050116126e2ddd436a75,
                0x10e412,
            ],
            [
                0xd34536fc5d22bc482ee234cad5a171,
                0xbd277ab51d947b3ca5b106779f5b5c,
                0xd7032a3faf5d557b644837159a2dc3,
                0x15cd26,
            ],
        ),
        c1: Fp2::new(
            [
                0x6936b660ccaab1e62694232ea1477a,
                0xf25be6740181f1e5e295cda62f289e,
                0xc31c6d5cd0dbab08c665c15ffb2d8f,
                0xe2303,
            ],
            [
                0xfadb875e950672048d1654ef16cc0b,
                0xd45a52d793a9bccc8a5ca077be7e75,
                0x8a106d89280e74785b2dad9b16a38e,
                0x142e38,
            ],
        ),
        c2: Fp2::new(
            [
                0x8dbae642d34613bfe1f68da26418af,
                0x2810ad92f6202621761b1ae8624ceb,
                0xabe1701c629b2ae2278b1148caedc7,
                0x6764d,
            ],
            [
                0xee83310f804608f88bae1efbffc05d,
                0x67d140367548c07d8707fd22ad0c12,
                0x14d21e748396b16562ab7a599405ae,
                0x1dac1,
            ],
        ),
    };

    let expected_result = Fp12 { c0, c1 };
    assert(res.eq(expected_result));
}

// #[test]
fn test_pairings_1() {
    // g and h from this test https://github.com/zkcrypto/bls12_381/blob/main/src/pairings.rs#L835
    let g = G1Affine::generator();

    let h = G2Affine::generator();

    let first_pairing = pairing(g.neg(), h);
    let second_pairing = pairing(g, h.neg());
    assert(first_pairing.eq(second_pairing)); // passess!
}

#[test]
fn test_bilinearity() {
    // https://github.com/zkcrypto/bls12_381/blob/main/src/pairings.rs#L835
    let g = G1Affine {
        x: BigNum {
            limbs: [
                0xfabb7be712b08862b6c9f3e9ff6ccc,
                0xdb4476eb5e63dc7314183a309224d,
                0x22fc235f054631617d5f134c6544f0,
                0x845e2,
            ],
        },
        y: BigNum {
            limbs: [
                0xbd5303645f9b22e005f9e2409b286e,
                0xce946b2e0489fbfdff1abb380f4d22,
                0xe2ff06bfa6a8e879aed13e82a17f86,
                0x362e6,
            ],
        },
        infinity: false,
    };

    let h = G2Affine {
        x: Fp2::new(
            [
                0xfb29f4c9b9f7aee559363b9545755c,
                0xbc1390ea1411adbee0eec731b72471,
                0x652b881bca2fe0b80fb53e5c6a8eee,
                0x154358,
            ],
            [
                0x82d763ca4de753d2b3ee304205157,
                0xfe44d89efc569c4bec0ee56c95428f,
                0xac4740bb4f874bdd05a017cb1eb950,
                0x1c656,
            ],
        ),
        y: Fp2::new(
            [
                0x23853393c1100055c5d34c423b3fd2,
                0xf4f22fff37029a585bb1038166e5ce,
                0xe508777c4881dfcf5df11b0425555f,
                0xe3d19,
            ],
            [
                0x1ec46236ee424e72d540ab93e7627a,
                0xf27c49db63114121eaae2ab56dcfad,
                0xd24956daaba822362ec9157039e792,
                0x193caf,
            ],
        ),
        infinity: false,
    };

    let expected = G1Affine {
        x: BigNum {
            limbs: [
                0xd7db4700715c7fc0a50d36e4d66c19,
                0xab5afe0efc1e6faf7476e6b42bafdc,
                0xf96908d5627e3dd1d3bb6f5476401,
                0x9ee21,
            ],
        },
        y: BigNum {
            limbs: [
                0x6485cf916eac99502622e4752b19f4,
                0x95a16b09abd95b4cbfedb2cbe87f7b,
                0x981c52a89257acb922bf3ffd74a55f,
                0x15878d,
            ],
        },
        infinity: false,
    };

    let pairing_g_h = pairing(g, h);
    // Expected value from zkcrypto
    // Gt(0x15cea6a78b89349dcd1c03695a74e72d3657af4ee2cf267337bc96363ef4a1c5d5d7a673cc3a3c1a1350043f99537d62 + 0x001701f60f3d4352c4d388c0f198854908c0e939719709c1b3f82d2a25cc7156a3838bc141e041c259849326fbd0839f*u + (0x171a6e4e0d608e6e6c8ac5ddcf8d6d2a950d06051e6b6c4d3feb6dc8dac2acadd345cadfb890454a2101a112f7471f0e + 0x03c68b88c25be6fd3611cc21c9968e4f87e541beeccc5170b8696a439bb666ad8a6608ab30ebc7dfe56eaf0dd9ab8439*u)*v + (0x1728d688caa8a671f90b55488e272f48c7c5ae32526d3635a5343eb02640358d9ac445c76a5d8f52f653bbaee04ba5ce + 0x05bfe708d381b066b83a883ba8251ce2ea6772cbde51e1322d82b2c8a026a2153f4822e20cb69b8b05003ee74e09cb48*u)*v^2 + (0x17b4a68dedcc04674e5f5739cf0d845ba801263f712ed4ddda59c1d9909148e3f28124ae770682c9b19233bf0bcfa00d + 0x164a54600cb2e24fec168f82542fbf98234dbb9ddf06503dc3c497da88b73db584ba19e685b1b398b51f40160e6c8f09*u + (0x0c214542fcc0cf10e3001fa5e82d398f6ba1ddd1ccdf133bfd75e033eae50aec66bd5e884b8c74d4c1c6ac7c01278ac5 + 0x0d7fcfdf3c1552ada42b7e0521bc2e7389436660c352ecbf2eedf30b77b6b501df302399e6240473af47abe56fc97478*u)*v + (0x130e7ce13260d601c89ee717acfd3d4e1d80f409dd2a5c38b176f0b64d3d0a224c502717270dfecf2b825ac24608215c + 0x14fd52fe9bfd08bbe23fcdf1d3bc5390c62e75a8786a72f8a343123a30a7c5f8d18508a21a2bf902f4db2c068913bc1c*u)*v^2)*w)
    // Converted to Noir
    let c0 = Fp6 {
        c0: Fp2::new(
            [
                0xd7a673cc3a3c1a1350043f99537d62,
                0xaf4ee2cf267337bc96363ef4a1c5d5,
                0xa78b89349dcd1c03695a74e72d3657,
                0x15cea6,
            ],
            [
                0x838bc141e041c259849326fbd0839f,
                0xe939719709c1b3f82d2a25cc7156a3,
                0xf60f3d4352c4d388c0f198854908c0,
                0x1701,
            ],
        ),
        c1: Fp2::new(
            [
                0x45cadfb890454a2101a112f7471f0e,
                0x6051e6b6c4d3feb6dc8dac2acadd3,
                0x4e0d608e6e6c8ac5ddcf8d6d2a950d,
                0x171a6e,
            ],
            [
                0x6608ab30ebc7dfe56eaf0dd9ab8439,
                0x41beeccc5170b8696a439bb666ad8a,
                0x88c25be6fd3611cc21c9968e4f87e5,
                0x3c68b,
            ],
        ),
        c2: Fp2::new(
            [
                0xc445c76a5d8f52f653bbaee04ba5ce,
                0xae32526d3635a5343eb02640358d9a,
                0x88caa8a671f90b55488e272f48c7c5,
                0x1728d6,
            ],
            [
                0x4822e20cb69b8b05003ee74e09cb48,
                0x72cbde51e1322d82b2c8a026a2153f,
                0x8d381b066b83a883ba8251ce2ea67,
                0x5bfe7,
            ],
        ),
    };

    let c1 = Fp6 {
        c0: Fp2::new(
            [
                0x8124ae770682c9b19233bf0bcfa00d,
                0x263f712ed4ddda59c1d9909148e3f2,
                0x8dedcc04674e5f5739cf0d845ba801,
                0x17b4a6,
            ],
            [
                0xba19e685b1b398b51f40160e6c8f09,
                0xbb9ddf06503dc3c497da88b73db584,
                0x600cb2e24fec168f82542fbf98234d,
                0x164a54,
            ],
        ),
        c1: Fp2::new(
            [
                0xbd5e884b8c74d4c1c6ac7c01278ac5,
                0xddd1ccdf133bfd75e033eae50aec66,
                0x42fcc0cf10e3001fa5e82d398f6ba1,
                0xc2145,
            ],
            [
                0x302399e6240473af47abe56fc97478,
                0x6660c352ecbf2eedf30b77b6b501df,
                0xdf3c1552ada42b7e0521bc2e738943,
                0xd7fcf,
            ],
        ),
        c2: Fp2::new(
            [
                0x502717270dfecf2b825ac24608215c,
                0xf409dd2a5c38b176f0b64d3d0a224c,
                0xe13260d601c89ee717acfd3d4e1d80,
                0x130e7c,
            ],
            [
                0x8508a21a2bf902f4db2c068913bc1c,
                0x75a8786a72f8a343123a30a7c5f8d1,
                0xfe9bfd08bbe23fcdf1d3bc5390c62e,
                0x14fd52,
            ],
        ),
    };

    let expected_pairing_g_h = Fp12 { c0, c1 };

    assert(pairing_g_h.eq(expected_pairing_g_h));
    let pairing_expected_id = pairing(expected, G2Affine::generator());
    assert(pairing_g_h.eq(pairing_expected_id)); // passess!
}
