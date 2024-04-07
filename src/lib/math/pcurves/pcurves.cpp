/*
* (C) 2024 Jack Lloyd
*
* Botan is released under the Simplified BSD License (see license.txt)
*/

#include <botan/internal/pcurves.h>

#include <botan/internal/pcurves_impl.h>

namespace Botan::PCurve {

// clang-format off

namespace secp256r1 {

class Params final : public EllipticCurveParameters<
   "FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF",
   "FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFC",
   "5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B",
   "FFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551",
   "6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296",
   "4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5",
   -10> {
};

class Curve final : public EllipticCurve<Params> {};

}  // namespace secp256r1

namespace secp384r1 {

class Params final : public EllipticCurveParameters<
   "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFF0000000000000000FFFFFFFF",
   "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFF0000000000000000FFFFFFFC",
   "B3312FA7E23EE7E4988E056BE3F82D19181D9C6EFE8141120314088F5013875AC656398D8A2ED19D2A85C8EDD3EC2AEF",
   "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFC7634D81F4372DDF581A0DB248B0A77AECEC196ACCC52973",
   "AA87CA22BE8B05378EB1C71EF320AD746E1D3B628BA79B9859F741E082542A385502F25DBF55296C3A545E3872760AB7",
   "3617DE4A96262C6F5D9E98BF9292DC29F8F41DBD289A147CE9DA3113B5F0B8C00A60B1CE1D7E819D7A431D7C90EA0E5F",
   -12> {
};

class Curve final : public EllipticCurve<Params> {};

}

namespace secp521r1 {

class Params final : public EllipticCurveParameters<
   "1FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
   "1FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFC",
   "51953EB9618E1C9A1F929A21A0B68540EEA2DA725B99B315F3B8B489918EF109E156193951EC7E937B1652C0BD3BB1BF073573DF883D2C34F1EF451FD46B503F00",
   "1FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFA51868783BF2F966B7FCC0148F709A5D03BB5C9B8899C47AEBB6FB71E91386409",
   "C6858E06B70404E9CD9E3ECB662395B4429C648139053FB521F828AF606B4D3DBAA14B5E77EFE75928FE1DC127A2FFA8DE3348B3C1856A429BF97E7E31C2E5BD66",
   "11839296A789A3BC0045C8A5FB42C7D1BD998F54449579B446817AFBD17273E662C97EE72995EF42640C550B9013FAD0761353C7086A272C24088BE94769FD16650",
   -4> {
};

class Curve final : public EllipticCurve<Params, P521Rep> {};

}

namespace secp256k1 {

class Params final : public EllipticCurveParameters<
   "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F",
   "0",
   "7",
   "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141",
   "79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798",
   "483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8"> {
};

class Curve final : public EllipticCurve<Params> {};

}

namespace brainpool256r1 {

class Params final : public EllipticCurveParameters<
   "A9FB57DBA1EEA9BC3E660A909D838D726E3BF623D52620282013481D1F6E5377",
   "7D5A0975FC2C3057EEF67530417AFFE7FB8055C126DC5C6CE94A4B44F330B5D9",
   "26DC5C6CE94A4B44F330B5D9BBD77CBF958416295CF7E1CE6BCCDC18FF8C07B6",
   "A9FB57DBA1EEA9BC3E660A909D838D718C397AA3B561A6F7901E0E82974856A7",
   "8BD2AEB9CB7E57CB2C4B482FFC81B7AFB9DE27E1E3BD23C23A4453BD9ACE3262",
   "547EF835C3DAC4FD97F8461A14611DC9C27745132DED8E545C1D54C72F046997"> {
};

class Curve final : public EllipticCurve<Params> {};

}

// clang-format on

namespace {

template <PrimeOrderCurveId::Id ID, typename C>
class PrimeOrderCurveImpl final : public PrimeOrderCurve {
   public:
      static_assert(C::OrderBits <= PrimeOrderCurve::MaximumBitLength);
      static_assert(C::PrimeFieldBits <= PrimeOrderCurve::MaximumBitLength);

      std::optional<PrimeOrderCurveId> curve_id() const override { return ID; }

      ProjectivePoint mul_by_g(const Scalar& scalar, RandomNumberGenerator& rng) const override {
         return stash(m_mul_by_g.mul(from_stash(scalar), rng));
      }

      ProjectivePoint mul(const AffinePoint& pt, const Scalar& scalar, RandomNumberGenerator& rng) const override {
         auto tbl = WindowedMulTable<C>(from_stash(pt));
         return stash(tbl.mul(from_stash(scalar), rng));
      }

      Scalar base_point_mul_x_mod_order(const Scalar& scalar, RandomNumberGenerator& rng) const override {
         auto pt = m_mul_by_g.mul(from_stash(scalar), rng);
         const auto x_bytes = pt.to_affine().x().serialize();
         return stash(C::Scalar::from_wide_bytes(std::span{x_bytes}));
      }

      AffinePoint generator() const override { return stash(C::G); }

      AffinePoint to_affine(const ProjectivePoint& pt) const override { return stash(from_stash(pt).to_affine()); }

      std::vector<uint8_t> serialize_point(const AffinePoint& pt, bool compress) const override {
         return from_stash(pt).serialize_to_vec(compress);
      }

      std::vector<uint8_t> serialize_scalar(const Scalar& scalar) const override {
         return from_stash(scalar).serialize_to_vec();
      }

      std::optional<Scalar> deserialize_scalar(std::span<const uint8_t> bytes) const override {
         if(auto scalar = C::Scalar::deserialize(bytes)) {
            return stash(*scalar);
         } else {
            return {};
         }
      }

      std::optional<AffinePoint> deserialize_point(std::span<const uint8_t> bytes) const override {
         if(auto pt = C::AffinePoint::deserialize(bytes)) {
            return stash(*pt);
         } else {
            return {};
         }
      }

      ProjectivePoint hash_to_curve(std::string_view hash,
                                    std::span<const uint8_t> input,
                                    std::span<const uint8_t> domain_sep,
                                    bool random_oracle) const override {
         if constexpr(C::ValidForSswuHash) {
            return stash(hash_to_curve_sswu<C>(hash, random_oracle, input, domain_sep));
         } else {
            throw Not_Implemented("Hash to curve is not implemented for this curve");
         }
      }

      Scalar scalar_add(const Scalar& a, const Scalar& b) const override {
         return stash(from_stash(a) + from_stash(b));
      }

      Scalar scalar_sub(const Scalar& a, const Scalar& b) const override {
         return stash(from_stash(a) - from_stash(b));
      }

      Scalar scalar_mul(const Scalar& a, const Scalar& b) const override {
         return stash(from_stash(a) * from_stash(b));
      }

      Scalar scalar_square(const Scalar& s) const override { return stash(from_stash(s).square()); }

      Scalar scalar_invert(const Scalar& s) const override { return stash(from_stash(s).invert()); }

      Scalar scalar_negate(const Scalar& s) const override { return stash(from_stash(s).negate()); }

      bool scalar_is_zero(const Scalar& s) const override { return from_stash(s).is_zero(); }

      PrimeOrderCurveImpl() : m_mul_by_g(C::G) {}

      static std::shared_ptr<const PrimeOrderCurve> instance() {
         static auto g_curve = std::make_shared<const PrimeOrderCurveImpl<ID, C>>();
         return g_curve;
      }

   private:
      static Scalar stash(const typename C::Scalar& s) {
         return Scalar::_make(instance(), s.template stash_value<StorageWords>());
      }

      static typename C::Scalar from_stash(const Scalar& s) {
         if(s._curve() != instance()) {
            throw Invalid_Argument("Curve mismatch");
         }
         return C::Scalar::from_stash(s._value());
      }

      static AffinePoint stash(const typename C::AffinePoint& pt) {
         auto x_w = pt.x().template stash_value<StorageWords>();
         auto y_w = pt.y().template stash_value<StorageWords>();
         return AffinePoint::_make(instance(), x_w, y_w);
      }

      static typename C::AffinePoint from_stash(const AffinePoint& pt) {
         if(pt._curve() != instance()) {
            throw Invalid_Argument("Curve mismatch");
         }
         auto x = C::FieldElement::from_stash(pt._x());
         auto y = C::FieldElement::from_stash(pt._y());
         return typename C::AffinePoint(x, y);
      }

      static ProjectivePoint stash(const typename C::ProjectivePoint& pt) {
         auto x_w = pt.x().template stash_value<StorageWords>();
         auto y_w = pt.y().template stash_value<StorageWords>();
         auto z_w = pt.z().template stash_value<StorageWords>();
         return ProjectivePoint::_make(instance(), x_w, y_w, z_w);
      }

      static typename C::ProjectivePoint from_stash(const ProjectivePoint& pt) {
         if(pt._curve() != instance()) {
            throw Invalid_Argument("Curve mismatch");
         }
         auto x = C::FieldElement::from_stash(pt._x());
         auto y = C::FieldElement::from_stash(pt._y());
         auto z = C::FieldElement::from_stash(pt._z());
         return typename C::ProjectivePoint(x, y, z);
      }

   private:
      const PrecomputedMulTable<C> m_mul_by_g;
};

}  // namespace

std::shared_ptr<const PrimeOrderCurve> PrimeOrderCurve::from_id(PrimeOrderCurveId id) {
   switch(id.code()) {
      case PrimeOrderCurveId::secp256r1:
         return PrimeOrderCurveImpl<PrimeOrderCurveId::secp256r1, secp256r1::Curve>::instance();
      case PrimeOrderCurveId::secp384r1:
         return PrimeOrderCurveImpl<PrimeOrderCurveId::secp384r1, secp384r1::Curve>::instance();
      case PrimeOrderCurveId::secp521r1:
         return PrimeOrderCurveImpl<PrimeOrderCurveId::secp521r1, secp521r1::Curve>::instance();
      case PrimeOrderCurveId::secp256k1:
         return PrimeOrderCurveImpl<PrimeOrderCurveId::secp256k1, secp256k1::Curve>::instance();
      case PrimeOrderCurveId::brainpool256r1:
         return PrimeOrderCurveImpl<PrimeOrderCurveId::brainpool256r1, brainpool256r1::Curve>::instance();
      default:
         return {};
   }
}

}  // namespace Botan::PCurve
