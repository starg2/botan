/*
* (C) 2024 Jack Lloyd
*
* Botan is released under the Simplified BSD License (see license.txt)
*/

#ifndef BOTAN_PCURVES_H_
#define BOTAN_PCURVES_H_

#include <botan/types.h>
#include <array>
#include <optional>
#include <span>
#include <string_view>
#include <vector>

namespace Botan {

class RandomNumberGenerator;

#if defined(BOTAN_HAS_ASN1)
class OID;
#endif

}  // namespace Botan

namespace Botan::PCurve {

/// Identifier for a named prime order curve
class BOTAN_TEST_API PrimeOrderCurveId {
   public:
      enum class Id : uint8_t {
         /// secp256r1 aka P-256
         secp256r1,
         /// secp384r1 aka P-384
         secp384r1,
         /// secp521r1 aka P-521
         secp521r1,
         /// secp256k1
         secp256k1,
         /// brainpool256r1
         brainpool256r1,
      };

      using enum Id;

      Id code() const { return m_id; }

      static std::vector<PrimeOrderCurveId> all();

      std::string to_string() const;

      PrimeOrderCurveId(Id id) : m_id(id) {}

      /// Map a string to a curve identifier
      static std::optional<PrimeOrderCurveId> from_string(std::string_view name);

#if defined(BOTAN_HAS_ASN1)
      /// Map an OID to a curve identifier
      ///
      /// Uses the internal OID table
      static std::optional<PrimeOrderCurveId> from_oid(const OID& oid);
#endif

   private:
      const Id m_id;
};

class BOTAN_TEST_API PrimeOrderCurve {
   public:
      /// Somewhat arbitrary maximum size for a field or scalar
      ///
      /// Sized to fit at least P-521
      static const size_t MaximumBitLength = 521;

      static const size_t MaximumByteLength = (MaximumBitLength + 7) / 8;

      /// Number of words used to store MaximumByteLength
      static const size_t StorageWords = (MaximumByteLength + sizeof(word) - 1) / sizeof(word);

      static std::shared_ptr<const PrimeOrderCurve> from_id(PrimeOrderCurveId id);

      /// Creates a generic non-optimized version
      //static std::shared_ptr<const PrimeOrderCurve> from_params(...);

      typedef std::array<word, StorageWords> StorageUnit;
      typedef std::shared_ptr<const PrimeOrderCurve> CurvePtr;

      class Scalar final {
         public:
            Scalar(const Scalar& other) = default;
            Scalar(Scalar&& other) = default;
            Scalar& operator=(const Scalar& other) = default;
            Scalar& operator=(Scalar&& other) = default;

            friend Scalar operator*(const Scalar& a, const Scalar& b) { return a.m_curve->scalar_mul(a, b); }

            friend Scalar operator+(const Scalar& a, const Scalar& b) { return a.m_curve->scalar_add(a, b); }

            friend Scalar operator-(const Scalar& a, const Scalar& b) { return a.m_curve->scalar_sub(a, b); }

            Scalar negate() const { return m_curve->scalar_negate(*this); }

            Scalar square() const { return m_curve->scalar_square(*this); }

            Scalar invert() const { return m_curve->scalar_invert(*this); }

            bool is_zero() const { return m_curve->scalar_is_zero(*this); }

            const auto& _curve() const { return m_curve; }

            const auto& _value() const { return m_value; }

            static Scalar _make(CurvePtr curve, StorageUnit v) { return Scalar(curve, v); }

         private:
            Scalar(CurvePtr curve, StorageUnit v) : m_curve(curve), m_value(v) {}

            CurvePtr m_curve;
            StorageUnit m_value;
      };

      class AffinePoint final {
         public:
            AffinePoint(const AffinePoint& other) = default;
            AffinePoint(AffinePoint&& other) = default;
            AffinePoint& operator=(const AffinePoint& other) = default;
            AffinePoint& operator=(AffinePoint&& other) = default;

            std::vector<uint8_t> serialize() const { return m_curve->serialize_point(*this, false); }

            std::vector<uint8_t> serialize_compressed() const { return m_curve->serialize_point(*this, true); }

            // AKA the point at infinity
            bool is_identity() const { return m_curve->affine_point_is_identity(*this); }

            const auto& _curve() const { return m_curve; }

            const auto& _x() const { return m_x; }

            const auto& _y() const { return m_y; }

            static AffinePoint _make(CurvePtr curve, StorageUnit x, StorageUnit y) { return AffinePoint(curve, x, y); }

         private:
            AffinePoint(CurvePtr curve, StorageUnit x, StorageUnit y) : m_curve(curve), m_x(x), m_y(y) {}

            CurvePtr m_curve;
            StorageUnit m_x;
            StorageUnit m_y;
      };

      class ProjectivePoint final {
         public:
            ProjectivePoint(const ProjectivePoint& other) = default;
            ProjectivePoint(ProjectivePoint&& other) = default;
            ProjectivePoint& operator=(const ProjectivePoint& other) = default;
            ProjectivePoint& operator=(ProjectivePoint&& other) = default;

            static ProjectivePoint from_projective(const AffinePoint& pt) {
               return pt._curve()->point_to_projective(pt);
            }

            AffinePoint to_affine() const { return m_curve->point_to_affine(*this); }

            bool is_identity() const { return m_curve->proj_point_is_identity(*this); }

            ProjectivePoint dbl() const { return m_curve->point_double(*this); }

            friend ProjectivePoint operator+(const ProjectivePoint& x, const ProjectivePoint& y) {
               return x.m_curve->point_add(x, y);
            }

            friend ProjectivePoint operator+(const ProjectivePoint& x, const AffinePoint& y) {
               return x.m_curve->point_add_mixed(x, y);
            }

            const auto& _curve() const { return m_curve; }

            const auto& _x() const { return m_x; }

            const auto& _y() const { return m_y; }

            const auto& _z() const { return m_z; }

            static ProjectivePoint _make(CurvePtr curve, StorageUnit x, StorageUnit y, StorageUnit z) {
               return ProjectivePoint(curve, x, y, z);
            }

         private:
            ProjectivePoint(CurvePtr curve, StorageUnit x, StorageUnit y, StorageUnit z) :
                  m_curve(curve), m_x(x), m_y(y), m_z(z) {}

            CurvePtr m_curve;
            StorageUnit m_x;
            StorageUnit m_y;
            StorageUnit m_z;
      };

      virtual ~PrimeOrderCurve() = default;

      virtual std::optional<PrimeOrderCurveId> curve_id() const = 0;

      virtual ProjectivePoint mul_by_g(const Scalar& scalar, RandomNumberGenerator& rng) const = 0;

      virtual Scalar base_point_mul_x_mod_order(const Scalar& scalar, RandomNumberGenerator& rng) const = 0;

      virtual ProjectivePoint mul(const AffinePoint& pt, const Scalar& scalar, RandomNumberGenerator& rng) const = 0;

      virtual AffinePoint generator() const = 0;

      virtual AffinePoint point_to_affine(const ProjectivePoint& pt) const = 0;

      virtual ProjectivePoint point_to_projective(const AffinePoint& pt) const = 0;

      virtual bool affine_point_is_identity(const AffinePoint& pt) const = 0;

      virtual bool proj_point_is_identity(const ProjectivePoint& pt) const = 0;

      virtual ProjectivePoint point_double(const ProjectivePoint& pt) const = 0;

      virtual ProjectivePoint point_add(const ProjectivePoint& a, const ProjectivePoint& b) const = 0;

      virtual ProjectivePoint point_add_mixed(const ProjectivePoint& a, const AffinePoint& b) const = 0;

      virtual std::vector<uint8_t> serialize_point(const AffinePoint& pt, bool compress) const = 0;

      virtual std::vector<uint8_t> serialize_scalar(const Scalar& scalar) const = 0;

      virtual std::optional<AffinePoint> deserialize_point(std::span<const uint8_t> bytes) const = 0;

      virtual std::optional<Scalar> deserialize_scalar(std::span<const uint8_t> bytes) const = 0;

      virtual Scalar scalar_zero() const = 0;
      virtual Scalar scalar_one() const = 0;
      virtual Scalar scalar_add(const Scalar& a, const Scalar& b) const = 0;
      virtual Scalar scalar_sub(const Scalar& a, const Scalar& b) const = 0;
      virtual Scalar scalar_mul(const Scalar& a, const Scalar& b) const = 0;
      virtual Scalar scalar_square(const Scalar& s) const = 0;
      virtual Scalar scalar_invert(const Scalar& s) const = 0;
      virtual Scalar scalar_negate(const Scalar& s) const = 0;
      virtual bool scalar_is_zero(const Scalar& s) const = 0;
      virtual bool scalar_equal(const Scalar& a, const Scalar& b) const = 0;

      virtual Scalar random_scalar(RandomNumberGenerator& rng) const = 0;

      virtual ProjectivePoint hash_to_curve(std::string_view hash,
                                            std::span<const uint8_t> input,
                                            std::span<const uint8_t> domain_sep,
                                            bool random_oracle) const = 0;
};

}  // namespace Botan::PCurve

#endif
