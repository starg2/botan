/*
* SM2 Signatures
* (C) 2017,2018 Ribose Inc
* (C) 2018,2024 Jack Lloyd
*
* Botan is released under the Simplified BSD License (see license.txt)
*/

#include <botan/sm2.h>

#include <botan/hash.h>
#include <botan/numthry.h>
#include <botan/internal/keypair.h>
#include <botan/internal/loadstor.h>
#include <botan/internal/parsing.h>
#include <botan/internal/pk_ops_impl.h>

namespace Botan {

std::string SM2_PublicKey::algo_name() const {
   return "SM2";
}

std::unique_ptr<Public_Key> SM2_PrivateKey::public_key() const {
   return std::make_unique<SM2_Signature_PublicKey>(domain(), public_point());
}

bool SM2_PrivateKey::check_key(RandomNumberGenerator& rng, bool strong) const {
   if(!EC_PrivateKey::check_key(rng, strong)) {
      return false;
   }

   // SM2 has an oddity in private key generation when compared to
   // other EC*DSA style signature algorithms described in ISO14888-3:
   // the private key x MUST be in ]0, q-1[ instead of ]0, q[.
   if(m_private_key < 1 || m_private_key >= m_domain_params.get_order() - 1) {
      return false;
   }

   if(!strong) {
      return true;
   }

   return KeyPair::signature_consistency_check(rng, *this, "user@example.com,SM3");
}

SM2_PrivateKey::SM2_PrivateKey(const AlgorithmIdentifier& alg_id, std::span<const uint8_t> key_bits) :
      EC_PrivateKey(alg_id, key_bits) {
   m_da_inv = domain().inverse_mod_order(m_private_key + 1);
}

SM2_PrivateKey::SM2_PrivateKey(RandomNumberGenerator& rng, const EC_Group& domain, const BigInt& x) :
      EC_PrivateKey(rng, domain, x) {
   m_da_inv = domain.inverse_mod_order(m_private_key + 1);
}

std::vector<uint8_t> sm2_compute_za(HashFunction& hash,
                                    std::string_view user_id,
                                    const EC_Group& domain,
                                    const EC_Point& pubkey) {
   if(user_id.size() >= 8192) {
      throw Invalid_Argument("SM2 user id too long to represent");
   }

   const uint16_t uid_len = static_cast<uint16_t>(8 * user_id.size());

   hash.update(get_byte<0>(uid_len));
   hash.update(get_byte<1>(uid_len));
   hash.update(user_id);

   const size_t p_bytes = domain.get_p_bytes();

   hash.update(BigInt::encode_1363(domain.get_a(), p_bytes));
   hash.update(BigInt::encode_1363(domain.get_b(), p_bytes));
   hash.update(BigInt::encode_1363(domain.get_g_x(), p_bytes));
   hash.update(BigInt::encode_1363(domain.get_g_y(), p_bytes));
   hash.update(pubkey.xy_bytes());

   std::vector<uint8_t> za(hash.output_length());
   hash.final(za.data());

   return za;
}

namespace {

/**
* SM2 signature operation
*/
class SM2_Signature_Operation final : public PK_Ops::Signature {
   public:
      SM2_Signature_Operation(const SM2_PrivateKey& sm2, std::string_view ident, std::string_view hash) :
            m_group(sm2.domain()),
            m_x(EC_Scalar::from_bigint(m_group, sm2.private_value())),
            m_da_inv(EC_Scalar::from_bigint(m_group, sm2.get_da_inv())) {
         if(hash == "Raw") {
            // m_hash is null, m_za is empty
         } else {
            m_hash = HashFunction::create_or_throw(hash);
            // ZA=H256(ENTLA || IDA || a || b || xG || yG || xA || yA)
            m_za = sm2_compute_za(*m_hash, ident, m_group, sm2.public_point());
            m_hash->update(m_za);
         }
      }

      size_t signature_length() const override { return 2 * m_group.get_order_bytes(); }

      void update(const uint8_t msg[], size_t msg_len) override {
         if(m_hash) {
            m_hash->update(msg, msg_len);
         } else {
            m_digest.insert(m_digest.end(), msg, msg + msg_len);
         }
      }

      secure_vector<uint8_t> sign(RandomNumberGenerator& rng) override;

      std::string hash_function() const override { return m_hash ? m_hash->name() : "Raw"; }

   private:
      const EC_Group m_group;
      const EC_Scalar m_x;
      const EC_Scalar m_da_inv;

      std::vector<uint8_t> m_za;
      secure_vector<uint8_t> m_digest;
      std::unique_ptr<HashFunction> m_hash;
      std::vector<BigInt> m_ws;
};

secure_vector<uint8_t> SM2_Signature_Operation::sign(RandomNumberGenerator& rng) {
   const auto e = [&]() {
      if(m_hash) {
         auto ie = EC_Scalar::from_bytes_mod_order(m_group, m_hash->final());
         // prepend ZA for next signature if any
         m_hash->update(m_za);
         return ie;
      } else {
         auto ie = EC_Scalar::from_bytes_mod_order(m_group, m_digest);
         m_digest.clear();
         return ie;
      }
   }();

   const auto k = EC_Scalar::random(m_group, rng);

   const auto r = EC_Scalar::gk_x_mod_order(k, rng, m_ws) + e;
   const auto s = (k - r * m_x) * m_da_inv;

   const auto rs = EC_Scalar::serialize_pair(r, s);
   return secure_vector<uint8_t>(rs.begin(), rs.end());
}

/**
* SM2 verification operation
*/
class SM2_Verification_Operation final : public PK_Ops::Verification {
   public:
      SM2_Verification_Operation(const SM2_PublicKey& sm2, std::string_view ident, std::string_view hash) :
            m_group(sm2.domain()), m_gy_mul(m_group, sm2.public_point()) {
         if(hash == "Raw") {
            // m_hash is null, m_za is empty
         } else {
            m_hash = HashFunction::create_or_throw(hash);
            // ZA=H256(ENTLA || IDA || a || b || xG || yG || xA || yA)
            m_za = sm2_compute_za(*m_hash, ident, m_group, sm2.public_point());
            m_hash->update(m_za);
         }
      }

      void update(const uint8_t msg[], size_t msg_len) override {
         if(m_hash) {
            m_hash->update(msg, msg_len);
         } else {
            m_digest.insert(m_digest.end(), msg, msg + msg_len);
         }
      }

      bool is_valid_signature(const uint8_t sig[], size_t sig_len) override;

      std::string hash_function() const override { return m_hash ? m_hash->name() : "Raw"; }

   private:
      const EC_Group m_group;
      const EC_Group::Mul2Table m_gy_mul;
      secure_vector<uint8_t> m_digest;
      std::vector<uint8_t> m_za;
      std::unique_ptr<HashFunction> m_hash;
};

bool SM2_Verification_Operation::is_valid_signature(const uint8_t sig[], size_t sig_len) {
   const auto e = [&]() {
      if(m_hash) {
         auto ie = EC_Scalar::from_bytes_mod_order(m_group, m_hash->final());
         // prepend ZA for next signature if any
         m_hash->update(m_za);
         return ie;
      } else {
         auto ie = EC_Scalar::from_bytes_mod_order(m_group, m_digest);
         m_digest.clear();
         return ie;
      }
   }();

   if(auto rs = EC_Scalar::deserialize_pair(m_group, std::span{sig, sig_len})) {
      const auto& [r, s] = rs.value();

      if(r.is_nonzero() && s.is_nonzero()) {
         const auto t = r + s;
         if(t.is_nonzero()) {
            if(const auto v = m_gy_mul.mul2_x_mod_order(s, t)) {
               return (v.value() + e) == r;
            }
         }
      }
   }
   return false;
}

void parse_sm2_param_string(std::string_view params, std::string& userid, std::string& hash) {
   // GM/T 0009-2012 specifies this as the default userid
   const std::string default_userid = "1234567812345678";

   // defaults:
   userid = default_userid;
   hash = "SM3";

   /*
   * SM2 parameters have the following possible formats:
   * Ident [since 2.2.0]
   * Ident,Hash [since 2.3.0]
   */

   auto comma = params.find(',');
   if(comma == std::string::npos) {
      userid = params;
   } else {
      userid = params.substr(0, comma);
      hash = params.substr(comma + 1, std::string::npos);
   }
}

}  // namespace

std::unique_ptr<Private_Key> SM2_PublicKey::generate_another(RandomNumberGenerator& rng) const {
   return std::make_unique<SM2_PrivateKey>(rng, domain());
}

std::unique_ptr<PK_Ops::Verification> SM2_PublicKey::create_verification_op(std::string_view params,
                                                                            std::string_view provider) const {
   if(provider == "base" || provider.empty()) {
      std::string userid, hash;
      parse_sm2_param_string(params, userid, hash);
      return std::make_unique<SM2_Verification_Operation>(*this, userid, hash);
   }

   throw Provider_Not_Found(algo_name(), provider);
}

std::unique_ptr<PK_Ops::Signature> SM2_PrivateKey::create_signature_op(RandomNumberGenerator& /*rng*/,
                                                                       std::string_view params,
                                                                       std::string_view provider) const {
   if(provider == "base" || provider.empty()) {
      std::string userid, hash;
      parse_sm2_param_string(params, userid, hash);
      return std::make_unique<SM2_Signature_Operation>(*this, userid, hash);
   }

   throw Provider_Not_Found(algo_name(), provider);
}

}  // namespace Botan
