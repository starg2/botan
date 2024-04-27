/*
* (C) 2024 Jack Lloyd
*
* Botan is released under the Simplified BSD License (see license.txt)
*/

#ifndef BOTAN_PCURVES_WRAP_H_
#define BOTAN_PCURVES_WRAP_H_

#include <botan/internal/pcurves.h>
#include <botan/internal/pcurves_impl.h>

namespace Botan::PCurve {

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

      ProjectivePoint mul2_vartime(const AffinePoint& pt1,
                                   const Scalar& s1,
                                   const AffinePoint& pt2,
                                   const Scalar& s2) const override {
         // Doesn't make sense to use a large window when we throw away the table
         // Even so, W=2 seems slightly better than W=1 here
         auto tbl = WindowedMul2Table<C, 2>(from_stash(pt1), from_stash(pt2));
         return stash(tbl.mul2_vartime(from_stash(s1), from_stash(s2)));
      }

      Scalar base_point_mul_x_mod_order(const Scalar& scalar, RandomNumberGenerator& rng) const override {
         auto pt = m_mul_by_g.mul(from_stash(scalar), rng);
         const auto x_bytes = pt.to_affine().x().serialize();
         return stash(C::Scalar::from_wide_bytes(std::span{x_bytes}));
      }

      AffinePoint generator() const override { return stash(C::G); }

      AffinePoint point_to_affine(const ProjectivePoint& pt) const override {
         return stash(from_stash(pt).to_affine());
      }

      ProjectivePoint point_to_projective(const AffinePoint& pt) const override {
         return stash(C::ProjectivePoint::from_affine(from_stash(pt)));
      }

      ProjectivePoint point_double(const ProjectivePoint& pt) const override { return stash(from_stash(pt).dbl()); }

      ProjectivePoint point_add(const ProjectivePoint& a, const ProjectivePoint& b) const override {
         return stash(from_stash(a) + from_stash(b));
      }

      ProjectivePoint point_add_mixed(const ProjectivePoint& a, const AffinePoint& b) const override {
         return stash(from_stash(a) + from_stash(b));
      }

      ProjectivePoint point_negate(const ProjectivePoint& pt) const override { return stash(from_stash(pt).negate()); }

      bool affine_point_is_identity(const AffinePoint& pt) const override { return from_stash(pt).is_identity(); }

      bool proj_point_is_identity(const ProjectivePoint& pt) const override { return from_stash(pt).is_identity(); }

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

      Scalar scalar_from_bits_with_trunc(std::span<const uint8_t> bytes) const override {
         return stash(C::Scalar::from_bits_with_trunc(bytes));
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

      bool scalar_equal(const Scalar& a, const Scalar& b) const override { return from_stash(a) == from_stash(b); }

      Scalar scalar_zero() const override { return stash(C::Scalar::zero()); }

      Scalar scalar_one() const override { return stash(C::Scalar::one()); }

      Scalar scalar_from_u32(uint32_t x) const override { return stash(C::Scalar::from_word(x)); }

      Scalar random_scalar(RandomNumberGenerator& rng) const override { return stash(C::Scalar::random(rng)); }

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

}  // namespace Botan::PCurve

#endif
