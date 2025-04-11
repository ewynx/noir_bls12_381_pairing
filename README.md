# Pairing over BLS12-381

Implementation of pairing over BLS12-381 in Noir. This uses the new [BigNum](https://github.com/noir-lang/noir-bignum) library. 

## Add dependency

This library uses [BigNum](https://github.com/noir-lang/noir-bignum) library version `v0.6.1`. 

```toml
[dependencies]
noir_bls12_381_pairing = { tag = "v0.2", git = "https://github.com/ewynx/noir_bls12_381_pairing" }
```
## Usage

To do a pairing define the 2 inputs of types `G1Affine` and `G2Affine` and apply the pairing. Example:
```rust
    let g = G1Affine::generator();

    let h = G2Affine::generator();

    let first_pairing = pairing(g.neg(), h);
    let second_pairing = pairing(g, h.neg());
    assert(first_pairing.eq(second_pairing));
```

## Implementation

This codebase uses the zkcrypto repository as a referece: https://github.com/zkcrypto/bls12_381. This has also been used for tests and test values. 
The Noir implementation passes pairing & bilinearity tests obtained from the zkcrypto repo. 

The BLS12-381 Fq field parameters comes from [BigNum](https://github.com/noir-lang/noir-bignum) library and the type `BLS12_381Fq` in `Fp2` follows the [definition](https://github.com/noir-lang/noir_bigcurve/blob/main/src/curves/bls12_381.nr#L60) in the BigCurve library:
```rust
pub type Fp = BigNum<4, BLS12_381_Fq_Params>;
```

We define it here to not have to import the full BigCurve library. 

## Test

Run all tests
```
nargo test
```

Note that a good amount of the tests are commented out because they take a fair amount of time (20-30 min) to run. For example `test_pairings_1` and `test_bilinearity` in `pairings.nr`. 

## Benchmarks

`nargo 0.35.0` and BigNum `0.4.2`:
One pairing:
- "acir_opcodes": 2.441.154
- "circuit_size": 3.210.964

`nargo 1.0.0-beta.3` and BigNum `0.6.1`:
One pairing:
- "acir_opcodes": 1.876.717
- "circuit_size": 5.403.322

For code snippet:
```rust
fn main(p: G1Affine, q: G2Affine) {
    let res = pairing(p, q);
}
```


## Related Noir work
- BLS12_381 Elliptic Curve Pairing and Signature Verification Library by @onurinanc: [repo](https://github.com/onurinanc/noir-bls-signature)
- Noir BigCurve library: [repo](https://github.com/noir-lang/noir_bigcurve)

This repo was forked and more curves are being added in [this repo](https://github.com/iAmMichaelConnor/noir_bls12_381_pairing) by @iAmMichaelConnor. 



