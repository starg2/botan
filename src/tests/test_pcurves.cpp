/*
* (C) 2014,2015,2019 Jack Lloyd
*
* Botan is released under the Simplified BSD License (see license.txt)
*/

#include "tests.h"

#if defined(BOTAN_HAS_PCURVES)
   #include "test_rng.h"
   #include <botan/mem_ops.h>
   #include <botan/internal/pcurves.h>
#endif

namespace Botan_Tests {

#if defined(BOTAN_HAS_PCURVES)

class Pcurve_Basemul_Tests final : public Text_Based_Test {
   public:
      Pcurve_Basemul_Tests() : Text_Based_Test("pubkey/ecc_base_point_mul.vec", "k,P") {}

      Test::Result run_one_test(const std::string& group_id, const VarMap& vars) override {
         Test::Result result("Pcurves base point multiply " + group_id);

         const auto k_bytes = vars.get_req_bin("k");
         const auto P_bytes = vars.get_req_bin("P");

         auto& rng = Test::rng();

         if(auto curve = Botan::PCurve::PrimeOrderCurve::from_name(group_id)) {
            if(auto scalar = curve->deserialize_scalar(k_bytes)) {
               auto pt2 = curve->mul_by_g(scalar.value(), rng).to_affine().serialize();
               result.test_eq("mul_by_g correct", pt2, P_bytes);

               auto g = curve->generator();
               auto pt3 = curve->mul(g, scalar.value(), rng).to_affine().serialize();
               result.test_eq("mul correct", pt3, P_bytes);
            } else {
               result.test_failure("Curve rejected scalar input");
            }
         }

         return result;
      }
};

BOTAN_REGISTER_TEST("pcurves", "pcurves_basemul", Pcurve_Basemul_Tests);

class Pcurve_Ecdsa_Sign_Tests final : public Text_Based_Test {
   public:
      Pcurve_Ecdsa_Sign_Tests() : Text_Based_Test("pubkey/ecdsa_pcurves.vec", "X,E,K,Sig") {}

      bool clear_between_callbacks() const override { return false; }

      static std::vector<uint8_t> ecdsa_sign(const Botan::PCurve::PrimeOrderCurve& curve,
                                             const Botan::PCurve::PrimeOrderCurve::Scalar& x,
                                             std::span<const uint8_t> msg,
                                             Botan::RandomNumberGenerator& nonce_rng,
                                             Botan::RandomNumberGenerator& rng) {
         const auto e = curve.scalar_from_bits_with_trunc(msg);
         const auto k = curve.random_scalar(nonce_rng);
         const auto r = curve.base_point_mul_x_mod_order(k, rng);
         const auto k_inv = k.invert();

         /*
         * Blind the input message and compute x*r+e as (b*x*r + b*e)/b
         */
         auto b = curve.random_scalar(rng);
         auto b_inv = b.invert();

         b = b.square();
         b_inv = b_inv.square();

         const auto be = b * e;
         const auto bx = b * x;

         const auto bxr_e = (bx * r) + be;

         const auto s = (k_inv * bxr_e) * b_inv;

         // With overwhelming probability, a bug rather than actual zero r/s
         if(r.is_zero() || s.is_zero()) {
            throw Botan::Internal_Error("During ECDSA signature generated zero r/s");
         }

         const auto r_bytes = r.serialize();
         const auto s_bytes = s.serialize();

         const size_t sig_len = s_bytes.size() + r_bytes.size();

         std::vector<uint8_t> sig(sig_len);
         Botan::copy_mem(&sig[0], r_bytes.data(), r_bytes.size());
         Botan::copy_mem(&sig[r_bytes.size()], s_bytes.data(), s_bytes.size());
         return sig;
      }

      Test::Result run_one_test(const std::string& group_id, const VarMap& vars) override {
         Test::Result result("Pcurves ECDSA sign " + group_id);

         const auto sk = vars.get_req_bin("X");
         const auto msg = vars.get_req_bin("E");
         const auto nonce = vars.get_req_bin("K");
         const auto expected_sig = vars.get_req_bin("Sig");

         Fixed_Output_RNG nonce_rng(nonce);

         auto curve = Botan::PCurve::PrimeOrderCurve::from_name(group_id);

         if(!curve) {
            result.test_note("Skipping test due to missing pcurve " + group_id);
            return result;
         }

         if(auto x = curve->deserialize_scalar(sk)) {
            auto sig = ecdsa_sign(*curve, x.value(), msg, nonce_rng, rng());
            result.test_eq("correct signature generated", sig, expected_sig);
         } else {
            result.test_failure("Curve rejected scalar input");
         }

         return result;
      }
};

BOTAN_REGISTER_TEST("pcurves", "pcurves_ecdsa_sign", Pcurve_Ecdsa_Sign_Tests);

class Pcurve_Point_Tests final : public Test {
   public:
      std::vector<Test::Result> run() override {
         std::vector<Test::Result> results;

         auto& rng = Test::rng();

         for(auto id : Botan::PCurve::PrimeOrderCurveId::all()) {
            Test::Result result("Pcurves point operations " + id.to_string());

            result.start_timer();

            auto curve = Botan::PCurve::PrimeOrderCurve::from_id(id);

            if(!curve) {
               result.test_failure("ID exists but curve does not");
               continue;
            }

            const auto zero = curve->scalar_zero();
            const auto one = curve->scalar_one();
            const auto g = curve->generator();
            const auto g_bytes = g.serialize();

            const auto inf = curve->mul_by_g(zero, rng);
            result.confirm("g*zero is point at infinity", inf.is_identity());

            std::vector<uint8_t> inf_bytes(g_bytes.size());
            inf_bytes[0] = 0x04;

            result.test_eq("infinity has expected encoding", inf.to_affine().serialize(), inf_bytes);

            const auto inf2 = inf.dbl();
            result.test_eq("infinity * 2 is infinity", inf2.to_affine().serialize(), inf_bytes);

            const auto inf3 = inf2 + inf;
            result.test_eq("infinity plus itself is infinity", inf3.to_affine().serialize(), inf_bytes);

            const auto g_one = curve->mul_by_g(one, rng);
            result.test_eq("g*one == generator", g_one.to_affine().serialize(), g_bytes);

            const auto g_plus_inf = g_one + inf;
            result.test_eq("g + inf == g", g_plus_inf.to_affine().serialize(), g_bytes);

            const auto g_plus_infa = g_one + inf.to_affine();
            result.test_eq("g + inf (affine) == g", g_plus_infa.to_affine().serialize(), g_bytes);

            const auto inf_plus_g = inf + g_one;
            result.test_eq("inf + g == g", inf_plus_g.to_affine().serialize(), g_bytes);

            const auto inf_plus_ga = inf + g_one.to_affine();
            result.test_eq("inf + g (affine) == g", inf_plus_ga.to_affine().serialize(), g_bytes);

            const auto g_neg_one = curve->mul_by_g(one.negate(), rng);

            const auto inf_from_g = g_one + g_neg_one;
            result.test_eq("g - g is infinity", inf_from_g.to_affine().serialize(), inf_bytes);

            const auto g_two = curve->mul_by_g(one + one, rng);
            const auto g_plus_g = g_one + g_one;
            result.test_eq("2*g == g+g", g_two.to_affine().serialize(), g_plus_g.to_affine().serialize());

            result.confirm("Scalar::zero is zero", curve->scalar_zero().is_zero());
            result.confirm("Scalar::one is not zero", !curve->scalar_one().is_zero());

            for(size_t i = 0; i != 16; ++i) {
               const auto pt = curve->mul_by_g(curve->random_scalar(rng), rng).to_affine();
               const auto a = curve->random_scalar(rng);
               const auto b = curve->random_scalar(rng);
               const auto c = a + b;

               const auto Pa = curve->mul(pt, a, rng);
               const auto Pb = curve->mul(pt, b, rng);
               const auto Pc = curve->mul(pt, c, rng);

               const auto Pc_bytes = Pc.to_affine().serialize();

               const auto Pab = Pa + Pb;
               result.test_eq("Pa + Pb == Pc", Pab.to_affine().serialize(), Pc_bytes);

               const auto Pba = Pb + Pa;
               result.test_eq("Pb + Pa == Pc", Pba.to_affine().serialize(), Pc_bytes);

               const auto Pabm = Pa + Pb.to_affine();
               result.test_eq("Pa + Pb == Pc (mixed)", Pabm.to_affine().serialize(), Pc_bytes);
               const auto Pbam = Pb + Pa.to_affine();
               result.test_eq("Pb + Pa == Pc (mixed)", Pbam.to_affine().serialize(), Pc_bytes);
            }

            for(size_t i = 0; i != 16; ++i) {
               const auto pt1 = curve->mul_by_g(curve->random_scalar(rng), rng).to_affine();
               const auto pt2 = curve->mul_by_g(curve->random_scalar(rng), rng).to_affine();

               const auto s1 = curve->random_scalar(rng);
               const auto s2 = curve->random_scalar(rng);

               const auto ref = (curve->mul(pt1, s1, rng) + curve->mul(pt2, s2, rng)).to_affine();
               const auto mul2 = curve->mul2_vartime(pt1, s1, pt2, s2).to_affine();

               result.test_eq("ref == mul2", ref.serialize(), mul2.serialize());
            }

            // Test cases where the two points have a linear relation
            for(size_t i = 0; i != 16; ++i) {
               const auto pt1 = curve->generator();

               auto pt2 = [&]() {
                  const auto lo = curve->scalar_from_u32(static_cast<uint32_t>(i / 2));
                  auto x = curve->mul_by_g(lo, rng);
                  if(i % 2 == 0) {
                     x = x.negate();
                  }
                  return x.to_affine();
               }();

               const auto s1 = curve->random_scalar(rng);
               const auto s2 = curve->random_scalar(rng);

               const auto ref = (curve->mul(pt1, s1, rng) + curve->mul(pt2, s2, rng)).to_affine();
               const auto mul2 = curve->mul2_vartime(pt1, s1, pt2, s2).to_affine();

               result.test_eq("ref == mul2 (linear relation)", ref.serialize(), mul2.serialize());
            }

            result.end_timer();

            results.push_back(result);
         }

         return results;
      }
};

BOTAN_REGISTER_TEST("pcurves", "pcurves_points", Pcurve_Point_Tests);

#endif

}  // namespace Botan_Tests
