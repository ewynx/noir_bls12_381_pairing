use crate::fp6::Fp6;
use crate::fp2::Fp2;

// c0 + c1*w
// F_{p^12} = F_{p^6} / w^2 - v
pub struct Fp12 {
    pub c0: Fp6,
    pub c1: Fp6,
}

impl Fp12 {
    pub fn new(c0: Fp6, c1: Fp6) -> Self {
        Fp12 { c0, c1 }
    }

    pub fn zero() -> Self {
        Fp12 { c0: Fp6::zero(), c1: Fp6::zero() }
    }

    pub fn one() -> Self {
        Fp12 { c0: Fp6::one(), c1: Fp6::zero() }
    }

    pub fn eq(self: Self, other: Self) -> bool {
        self.c0.eq(other.c0) & self.c1.eq(other.c1)
    }

    pub fn add(self: Self, other: Self) -> Self {
        Self { c0: self.c0.add(other.c0), c1: self.c1.add(other.c1) }
    }

    pub fn sub(self: Self, other: Self) -> Self {
        Self { c0: self.c0.sub(other.c0), c1: self.c1.sub(other.c1) }
    }

    pub fn neg(self: Self) -> Self {
        Self { c0: self.c0.neg(), c1: self.c1.neg() }
    }

    pub fn conjugate(self: Self) -> Self {
        Self { c0: self.c0, c1: self.c1.neg() }
    }

    pub fn square(self: Fp12) -> Fp12 {
        // ab = self.c0 * self.c1
        let ab = self.c0.mul(self.c1);

        // c0_plus_c1 = self.c0 + self.c1
        let c0_plus_c1 = self.c0.add(self.c1);

        // temp = self.c1 * nonresidue
        let temp = self.c1.mul_by_nonresidue();

        // temp_with_c0 = temp + self.c0
        // temp_with_c0 = (self.c1 * nonresidue) + self.c0
        let temp_with_c0 = temp.add(self.c0);

        // step 1: c0 = temp_with_c0 * c0_plus_c1 - ab
        // c0 = (self.c1 * nonresidue + self.c0) * (self.c0 + self.c1) - ab
        let mut c0 = (temp_with_c0.mul(c0_plus_c1)).sub(ab);
        // step 2: c0 = ((self.c1 * nonresidue + self.c0) * (self.c0 + self.c1) - ab) - (ab * nonresidue)
        c0 = c0.sub(ab.mul_by_nonresidue());

        // c1 = 2 * ab
        let c1 = ab.add(ab);

        // Return the final Fp12 structure
        Fp12 { c0, c1 }
    }

    pub fn mul(self: Fp12, other: Fp12) -> Fp12 {
        // aa = self.c0 * other.c0
        let aa = self.c0.mul(other.c0);

        // bb = self.c1 * other.c1
        let bb = self.c1.mul(other.c1);

        // o = other.c0 + other.c1
        let o = other.c0.add(other.c1);

        // c1 = (self.c1 + self.c0) * (other.c0 + other.c1) - aa - bb
        let mut c1 = self.c1.add(self.c0);
        c1 = c1.mul(o);
        c1 = c1.sub(aa);
        c1 = c1.sub(bb);

        // c0 = bb * nonresidue + aa
        let mut c0 = bb.mul_by_nonresidue();
        c0 = c0.add(aa);

        // Return the final Fp12 structure
        Fp12 { c0, c1 }
    }

    pub fn frobenius_map(self: Self) -> Self {
        let c0 = self.c0.frobenius_map();
        let mut c1 = self.c1.frobenius_map();

        // c1 = c1 * (u + 1)^((p - 1) / 6)
        // (u + 1)^((p - 1) / 6)
        let constant = Fp6 {
            c0: Fp2::new(
                [
                    0x7ea53d63e7813d8d0775ed92235fb8,
                    0x3fd3cbd5f4f7b2443d784bab9c4f6,
                    0xbf02bb0667c231beb4202c0d1f0fd6,
                    0x1904d3,
                ],
                [
                    0x2d5ac14d6c7ec22cf78a126ddc4af3,
                    0x4787b6c7b36fec0c8ec971f63c5f28,
                    0x2b36c4e03288e9e902231f9fb854a1,
                    0xfc3e,
                ],
            ),
            c1: Fp2::zero(),
            c2: Fp2::zero(),
        };

        c1 = c1.mul(constant);

        Fp12 { c0, c1 }
    }

    // Returns inverse and errors if there is no inverse
    pub fn invert(self: Self) -> Self {
        let mut t = (self.c0.square()).sub((self.c1.square().mul_by_nonresidue()));
        // This will error if there's no inverse
        t = t.invert();
        Fp12 { c0: self.c0.mul(t), c1: self.c1.mul(t.neg()) }
    }

    pub fn mul_by_014(self: Self, c0: Fp2, c1: Fp2, c4: Fp2) -> Fp12 {
        let aa = self.c0.mul_by_01(c0, c1);
        let bb = self.c1.mul_by_1(c4);
        let o = c1.add(c4);
        let mut c1 = self.c1.add(self.c0);
        c1 = c1.mul_by_01(c0, o);
        c1 = c1.sub(aa).sub(bb);
        let mut c0 = bb;
        c0 = c0.mul_by_nonresidue();
        c0 = c0.add(aa);

        Fp12 { c0, c1 }
    }
}

// #[test]
fn test_arithmetic() {
    let a = Fp12 {
        c0: Fp6 {
            c0: Fp2::new(
                [
                    0xb3af6c44e48eb9b7d1b4b1c41e97b8,
                    0x2b689fdd5da6619b2e16bf2482973d,
                    0x4eb1b4adb3916c37153fa843bc1724,
                    0x1ddb59,
                ],
                [
                    0x687898c89dac3938647a0337390e91,
                    0xfb34d9652fc9a5b7632f105f949e70,
                    0xd4f2477096ad91b69a7a646bec7844,
                    0xed52c,
                ],
            ),
            c1: Fp2::new(
                [
                    0x8d42a6eee24c21ad6d25f3b895a378,
                    0xaf34b84ed0165f4bc59f666809a305,
                    0xa9a43882ae8a5cba9a9f45ace2559c,
                    0x1a2632,
                ],
                [
                    0x836efb730e4bf48a0ec0a552727f8f,
                    0x6de2efd5fb2c80577e407d56e3584d,
                    0x10275ac5fff269bdd4506932e43405,
                    0x72a61,
                ],
            ),
            c2: Fp2::new(
                [
                    0x85d961dd2944c554f41dbdfa55dc3e,
                    0x94dd5f8f62f58e54520ec7ab9f7bdc,
                    0xd76c86b02647531ee625f52f154f17,
                    0x12bebd,
                ],
                [
                    0x2f83f04cbabb9a42140b1cb7871cff,
                    0x2d6ae61a88834a22708ea42961c300,
                    0xcdd27a7dae8e0f9a9612d0a86d41d6,
                    0x915aa,
                ],
            ),
        },
        c1: Fp6 {
            c0: Fp2::new(
                [
                    0x175068bbfe0b1a2bbff8b7a5351255,
                    0x1fa519cb5837aa2e1c4671ef9d5374,
                    0xa3023daea697df7b912350025858bd,
                    0x1676c5,
                ],
                [
                    0x89a7261e177c96b86bb084e8844c12,
                    0x23900e93f79297dc34656842864f9f,
                    0x95b4e58f05d199c2c7760527c09bbe,
                    0x16b0ae,
                ],
            ),
            c1: Fp2::new(
                [
                    0x8350cdccab87c3593f38f0df87efc1,
                    0x8925cda7f91b443d571d7e2bbf8b10,
                    0x6f39338c7c39fb7434af47f709ae2b,
                    0xe95c5,
                ],
                [
                    0xa5ff1ac8f2aa1d45a9bd44f53561ce,
                    0xe52df8692146afee27747a9ec2d33e,
                    0xe3165b35d21f915f3332d58ce673c7,
                    0x12cfd8,
                ],
            ),
            c2: Fp2::new(
                [
                    0x38fb1613304b1b18a981ee7ea62d12,
                    0xce7d54f7127829d3b51d6a0d90c481,
                    0xd46daeed9ca25aa97c964fe9d14b69,
                    0x163a1a,
                ],
                [
                    0x823f42cc9582b3199546dbec99ad49,
                    0xff1210846645f54b39d0cfc0742bd9,
                    0x5b0e86ab5bdd5ada02591ced1e0280,
                    0x19ece0,
                ],
            ),
        },
    };

    let b = Fp12 {
        c0: Fp6 {
            c0: Fp2::new(
                [
                    0x6418e17bb08a26d334754a88af9ed8,
                    0x5c4d4bf9a95873660b982f3b389486,
                    0xe6b9301b716296eceb90902265e9c2,
                    0x344eb,
                ],
                [
                    0x2c4ff8f95c0d905b6460185b3881b2,
                    0x48281080ec95757ca0e3cd5ce42ddf,
                    0xe53b49ce6208454794fc46ac3383f3,
                    0x4c61f,
                ],
            ),
            c1: Fp2::new(
                [
                    0x1a0f38f91126c4385f50f0d87a3bf9,
                    0xdfb2b1a80a69779abdccea85ff97c9,
                    0xe9fc25a78a348ae518826015446871,
                    0x1c49d1,
                ],
                [
                    0x640bec9596ee7f23bf29bda216eb07,
                    0x673020119d960490f956dd027725dd,
                    0xb27c4fe60439d127a1c7a41318dcad,
                    0x1bd151,
                ],
            ),
            c2: Fp2::new(
                [
                    0xba255b575ff9286c9388bca4019c61,
                    0x4c51a6b0b47e7d7598a652f0a73bd8,
                    0xfdc9e2d0587d113f13f4826beb11f9,
                    0x1761ea,
                ],
                [
                    0xea9ac885f1c1478407e854749e0022,
                    0xac5c8ba611a851b9d14acd557e2174,
                    0x37319d961520599c35245529ac93b5,
                    0x11f751,
                ],
            ),
        },
        c1: Fp6 {
            c0: Fp2::new(
                [
                    0xd12f1144e528d5cce5f7eb98ea5d43,
                    0xad4ce3734b537475320ace5bf8e001,
                    0xc0be42a2a9066aa38e777bbd4497a8,
                    0x10951b,
                ],
                [
                    0x270a62aab9b6a4fa161969069340a,
                    0x71b7c596b5fe7ea5a04c14f539c618,
                    0xc0e27d563e716251e497b1851e5443,
                    0x1c3fae,
                ],
            ),
            c1: Fp2::new(
                [
                    0xb1b0ff6adab8c07a0309ac928eab17,
                    0x202bb2a36371c30dee11befaad53ec,
                    0x5c6128c31815fa8f11adff21bc0e93,
                    0x5931,
                ],
                [
                    0x38340422d00d5a2fd80c18732c1166,
                    0xe56145364add39e1433e7666bec7ff,
                    0x7bd111cb2e1ffd57621bfa33f6baec,
                    0x3e8ca,
                ],
            ),
            c2: Fp2::new(
                [
                    0x9951a3096171aa6b2b8f17ab3f4e97,
                    0x4c6a214832ad14421d6205a401a29,
                    0x746a9493f6e699f74f16cb44af029,
                    0x1e06a7,
                ],
                [
                    0x41bca809ec9e7f1fc87ac6e3a3cce9,
                    0xb8a9fd6b10a1fa0a41984d71297189,
                    0x2606b9b5501e92e16fba99ef7b27aa,
                    0x143564,
                ],
            ),
        },
    };

    let c = Fp12 {
        c0: Fp6 {
            c0: Fp2::new(
                [
                    0x720c1a648a6dc4377626d77b03a764,
                    0x6e3db0d99300ef2265f0995ef97dc,
                    0x435d78cbd6c7416bb0370f21ff21cf,
                    0x13a44e,
                ],
                [
                    0x84bfec07a5800c9c7c55d23ec6671b,
                    0x2fdc513cdeffe29adcc3cb947507b5,
                    0x76eb72e11ce785c99eeb4fb5d8040b,
                    0x121e34,
                ],
            ),
            c1: Fp2::new(
                [
                    0xbe7fbf6e2353e21c543e4bc02e0a9c,
                    0xf46431d46c4a56f11e303698d2eebf,
                    0xd496eb7c1c4f6e1d05f6939c63a421,
                    0xd0b56,
                ],
                [
                    0xcd31244bc979628056b1d106ecc8d8,
                    0xfdaf068872d24f31805084a57f6c62,
                    0x62991a2ebf541d9d6ef7c91dc50f97,
                    0x14ab7a,
                ],
            ),
            c2: Fp2::new(
                [
                    0x5f048867ad692fa8210d94fabfa2fd,
                    0x2e62e766035568f1a9721a764287ab,
                    0xbab9bfcac9c257038f6a9f9396dbb2,
                    0x1da83b,
                ],
                [
                    0x53ac5d922508c72ca4362d04e32dfa,
                    0x166c20653130bad0fb9cad59fc12fc,
                    0x655683e232b0eb233b495982cc79f0,
                    0x42738,
                ],
            ),
        },
        c1: Fp6 {
            c0: Fp2::new(
                [
                    0x571bacfab9112ef16ee344f0c16268,
                    0x7a3f5e2e9db974f0d87ea80fa1d428,
                    0xb0104e0b6a7ef0c318eacf111802b0,
                    0xa9736,
                ],
                [
                    0xe41318acb7a8c0e6a3458d84fd1c44,
                    0xaa72d9c749ee1ace125dbefcd38c5b,
                    0xc3217370529d5d92bcb222724d76fd,
                    0x1812f3,
                ],
            ),
            c1: Fp2::new(
                [
                    0x7cd66ddcbf855f85e9692aeabff5d8,
                    0xef4f784c404b42a9f697767fd1aaa7,
                    0x6fdc056bd66cf95fedc12f069ba3e1,
                    0x86933,
                ],
                [
                    0x35205828d34eae50579abdc513cbee,
                    0xa36ac47cdfab5279bd233f26346f95,
                    0x1a3ab7932780db066622c871503529,
                    0xe5913,
                ],
            ),
            c2: Fp2::new(
                [
                    0xb1f4a2c5ab1fde30be6b78791d983c,
                    0xf9cb088c4af65ede4cf0b54f5f5430,
                    0xd49efbdab23637124bda021af5f7b8,
                    0x1b50ba,
                ],
                [
                    0x866177018d092807ccdc7d8e9a1ec7,
                    0x69fcac9ef500eb8a513c902d90ad92,
                    0x1cec57a784c5e8abfb1a4223227d2f,
                    0x74c86,
                ],
            ),
        },
    };

    let a_square = a.square();
    let a_times_a = a.mul(a);
    assert(a_square.eq(a_times_a));

    let b_square = b.square();
    let b_times_b = b.mul(b);
    assert(b_square.eq(b_times_b));

    // (a+b)*c = a*c + b*c
    let lhs = (a.add(b)).mul(c);
    let rhs = (a.mul(c)).add(b.mul(c));
    assert(lhs.eq(rhs));
}

// #[test]
fn test_invert() {
    let a = Fp12 {
        c0: Fp6 {
            c0: Fp2::new(
                [
                    0xcd5d36546bea81490d5309111ddceb,
                    0x60fb0821e8477be42000281597c30d,
                    0xc5bceea125200ead53b9846d34e496,
                    0xe6106,
                ],
                [
                    0x8ebdbdbc10648d74bccd3a998f7d14,
                    0xbc1f6ba6171b9872df351363aff1fe,
                    0xae277a7948eb47594c31af60354995,
                    0x138b40,
                ],
            ),
            c1: Fp2::new(
                [
                    0x93de59c75d21d5a02ed285a867278d,
                    0x68743f57ffe21ee7ce991a169d868c,
                    0xc5eaf10f2f2c5b130bafcb25ceee18,
                    0x14cf6e,
                ],
                [
                    0x27d65314e3e04fe2e42598308f46f6,
                    0xcc9f7ad74ebbce368ea95ae986a3dc,
                    0x891b4de53b5c393493ba1f272679ad,
                    0x143ce9,
                ],
            ),
            c2: Fp2::new(
                [
                    0x5589d21d49d8839b50d88f7990521a,
                    0xfb74891851fe5934d176112797d7af,
                    0x7b93b06afebc4aa2b67b0c035d2567,
                    0x91a36,
                ],
                [
                    0xfe74e25003c99c8d3ed0d02c49d9bc,
                    0xdb28b3c2cf0f4cfbad2e238c6c6c0d,
                    0xafac7905021bc99cfb6dca725c7127,
                    0x7ef0f,
                ],
            ),
        },
        c1: Fp6 {
            c0: Fp2::new(
                [
                    0x364af51cbd5ff7e6819c1d34723deb,
                    0xbf5560c1cad36ff82b7ae6b4e9b669,
                    0xb8b43b7067d5aebbfaf30af941fb6b,
                    0x5d606,
                ],
                [
                    0xd6283df9238fd8892faacf218aab13,
                    0x778fff2b53389c2467832ee0c98fae,
                    0xbe763b594c7de966e4aee7235e03f1,
                    0xded4f,
                ],
            ),
            c1: Fp2::new(
                [
                    0xb1b1dc8139b0ea6c55d4787a32c7fd,
                    0xb0f3f1c9b39a58c54a1d41048049b6,
                    0xa721d13a399308f0e7f083ded43cb8,
                    0x1c010c,
                ],
                [
                    0x973a2a6b834b9ec30d3fa2125d635b,
                    0xcb60de306415b58bb2b74169b3b22d,
                    0x2d5df2d5a06c53ff0f6738cd85426e,
                    0x15dbd0,
                ],
            ),
            c2: Fp2::new(
                [
                    0x3a5579ef31cee3bf72ede9a4ab4751,
                    0xd6e693b6dc6b8074b6832bc70ba85,
                    0x6d2789051356b780ff5ff6d005f5dc,
                    0x1460d,
                ],
                [
                    0xd2c6fb516a060d3f63a44ecc0b8532,
                    0x4d8aea5422fec4e9c93bd8b65553e9,
                    0xd7333204133e13cdb7faf3bcab8632,
                    0x187373,
                ],
            ),
        },
    };

    let a_inv = a.invert();
    let check = a.mul(a_inv);
    assert(check.eq(Fp12::one()));
}

// #[test]
fn test_frobenius_map() {
    let a = Fp12 {
        c0: Fp6 {
            c0: Fp2::new(
                [
                    0x7105ed1df19d3a68cae6ca0a7b1fb9,
                    0x9aa950b7187117843a307a2068c46e,
                    0x838cd18c35a8ff5d6f680a8a229ec5,
                    0xdb672,
                ],
                [
                    0x12dc9536b0cb9bbb6aa1510a73c620,
                    0x42fc8e4eecd96a4da199d5cc569798,
                    0x6a705750eaa4278ced7d8a60ccbdb,
                    0x18cab0,
                ],
            ),
            c1: Fp2::new(
                [
                    0x6eb3cd88be4d899d6b1a927b258785,
                    0xc4e54ca9d34a94ca4ac0fcbb572fa5,
                    0xb8e42b654c384e8d268700b81bbd8c,
                    0xc0211,
                ],
                [
                    0x3ac0dc4d7ea4b1fa45bb455b7ed2ab,
                    0xa2444aca03e12734ee13890a36233b,
                    0xeadea5baebd4b1d4f4422c3b86f1ae,
                    0x42bf4,
                ],
            ),
            c2: Fp2::new(
                [
                    0x5e3ee68784da522e9178953fdc5fc5,
                    0xc15262afda0bd7c9a33a2ed786cb1c,
                    0x428af65ec0fb6b0b786b355ca3a4d7,
                    0xf9909,
                ],
                [
                    0xf338346a22198aed1e5e7a7c98885b,
                    0xaaa01918fd42c50b1e8d42f182d03e,
                    0x6bc534e1ba7eeb3feb9710da6262de,
                    0x12de03,
                ],
            ),
        },
        c1: Fp6 {
            c0: Fp2::new(
                [
                    0x3cead8477cfee8abd857fd3e23e35f,
                    0xef8f21ca80b6310b2c9a43e736abc0,
                    0x9bfc772cbcdac05e480ef3b56ab31e,
                    0x168188,
                ],
                [
                    0x42a29fdc357592f79632d329095662,
                    0x4cd2ca25ce13cdd8431740aa9fd6b4,
                    0x7ebe3a75a2b05294341ffc22d802fd,
                    0x42024,
                ],
            ),
            c1: Fp2::new(
                [
                    0x1a40e8bdd817a157c61810d8c8fe18,
                    0xd432430c44178ed2047356b54c9ceb,
                    0xeb1711f8a4c1a7e0c7356098aecfd,
                    0x48013,
                ],
                [
                    0x2d3a96c4b6957c1783cc12bbdc7628,
                    0xf097be530a9903ea8d8c4fb563d42f,
                    0x31b267ade8c1b22a9215e7665244d8,
                    0x12a1a7,
                ],
            ),
            c2: Fp2::new(
                [
                    0x6b2e550f54623de5d37328d7e87e0,
                    0xd375db0ce7999be6e147ef61cb0e0f,
                    0xed6eaeb86d0031a6f5adfa108d042,
                    0x5396a,
                ],
                [
                    0xf499b0450ac95711376f9a9c54425e,
                    0xcb2d7b1c9acbf9d56704c582b5feca,
                    0x92ac3d19d606a58711c34169b324b2,
                    0xaaa4f,
                ],
            ),
        },
    };
    let res = a
        .frobenius_map()
        .frobenius_map()
        .frobenius_map()
        .frobenius_map()
        .frobenius_map()
        .frobenius_map()
        .frobenius_map()
        .frobenius_map()
        .frobenius_map()
        .frobenius_map()
        .frobenius_map()
        .frobenius_map();

    assert(res.eq(a));
}
