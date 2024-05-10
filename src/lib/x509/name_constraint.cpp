/*
* X.509 Name Constraint
* (C) 2015 Kai Michaelis
*     2024 Jack Lloyd
*
* Botan is released under the Simplified BSD License (see license.txt)
*/

#include <botan/pkix_types.h>

#include <botan/ber_dec.h>
#include <botan/x509cert.h>
#include <botan/internal/fmt.h>
#include <botan/internal/loadstor.h>
#include <botan/internal/parsing.h>
#include <functional>
#include <sstream>

namespace Botan {

class DER_Encoder;

std::string GeneralName::type() const {
   switch(m_type) {
      case NameType::Empty:
         throw Encoding_Error("Could not convert empty NameType to string");
      case NameType::RFC822:
         return "RFC822";
      case NameType::DNS:
         return "DNS";
      case NameType::URI:
         return "URI";
      case NameType::DN:
         return "DN";
      case NameType::IP:
         return "IP";
   }

   BOTAN_ASSERT_UNREACHABLE();
}

std::string GeneralName::name() const {
   const size_t index = m_names.index();

   if(index == 0) {
      return std::get<0>(m_names);
   } else if(index == 1) {
      return std::get<1>(m_names);
   } else if(index == 2) {
      return std::get<2>(m_names);
   } else if(index == 3) {
      return std::get<3>(m_names).to_string();
   } else if(index == 4) {
      auto [net, mask] = std::get<4>(m_names);
      return fmt("{}/{}", ipv4_to_string(net), ipv4_to_string(mask));
   } else {
      BOTAN_ASSERT_UNREACHABLE();
   }
}

void GeneralName::encode_into(DER_Encoder& /*to*/) const {
   throw Not_Implemented("GeneralName encoding");
}

void GeneralName::decode_from(BER_Decoder& ber) {
   BER_Object obj = ber.get_next_object();

   if(obj.is_a(1, ASN1_Class::ContextSpecific)) {
      m_type = NameType::RFC822;
      m_names.emplace<0>(ASN1::to_string(obj));
   } else if(obj.is_a(2, ASN1_Class::ContextSpecific)) {
      m_type = NameType::DNS;
      // Store it in case insensitive form so we don't have to do it
      // again while matching
      m_names.emplace<1>(tolower_string(ASN1::to_string(obj)));
   } else if(obj.is_a(6, ASN1_Class::ContextSpecific)) {
      m_type = NameType::URI;
      m_names.emplace<2>(ASN1::to_string(obj));
   } else if(obj.is_a(4, ASN1_Class::ContextSpecific | ASN1_Class::Constructed)) {
      X509_DN dn;
      BER_Decoder dec(obj);
      dn.decode_from(dec);

      m_type = NameType::DN;
      m_names.emplace<3>(dn);
   } else if(obj.is_a(7, ASN1_Class::ContextSpecific)) {
      if(obj.length() == 8) {
         const uint32_t net = load_be<uint32_t>(obj.bits(), 0);
         const uint32_t mask = load_be<uint32_t>(obj.bits(), 1);

         m_type = NameType::IP;
         m_names.emplace<4>(std::make_pair(net, mask));
      } else if(obj.length() == 32) {
         throw Decoding_Error("Unsupported IPv6 name constraint");
      } else {
         throw Decoding_Error("Invalid IP name constraint size " + std::to_string(obj.length()));
      }
   } else {
      throw Decoding_Error("Found unknown GeneralName type");
   }
}

GeneralName::MatchResult GeneralName::matches(const X509_Certificate& cert) const {
   class MatchScore final {
      public:
         MatchScore() : m_any(false), m_some(false), m_all(true) {}

         void add(bool m) {
            m_any = true;
            m_some |= m;
            m_all &= m;
         }

         MatchResult result() const {
            if(!m_any) {
               return MatchResult::NotFound;
            } else if(m_all) {
               return MatchResult::All;
            } else if(m_some) {
               return MatchResult::Some;
            } else {
               return MatchResult::None;
            }
         }

      private:
         bool m_any;
         bool m_some;
         bool m_all;
   };

   const X509_DN& dn = cert.subject_dn();
   const AlternativeName& alt_name = cert.subject_alt_name();

   MatchScore score;

   if(m_type == NameType::DNS) {
      const auto& constraint = std::get<1>(m_names);

      const auto& alt_names = alt_name.dns();

      for(const std::string& dns : alt_names) {
         score.add(matches_dns(dns, constraint));
      }

      if(alt_names.empty()) {
         // Check CN instead...
         for(const std::string& cn : dn.get_attribute("CN")) {
            score.add(matches_dns(cn, constraint));
         }
      }
   } else if(m_type == NameType::DN) {
      const X509_DN& constraint = std::get<3>(m_names);
      score.add(matches_dn(dn, constraint));

      for(const auto& alt_dn : alt_name.directory_names()) {
         score.add(matches_dn(alt_dn, constraint));
      }
   } else if(m_type == NameType::IP) {
      auto [net, mask] = std::get<4>(m_names);

      for(uint32_t ipv4 : alt_name.ipv4_address()) {
         bool match = (ipv4 & mask) == net;
         score.add(match);
      }
   } else {
      // URI and email name constraint matching not implemented
      return MatchResult::UnknownType;
   }

   return score.result();
}

//static
bool GeneralName::matches_dns(const std::string& name, const std::string& constraint) {
   // constraint is assumed already tolower
   if(name.size() == constraint.size()) {
      return tolower_string(name) == constraint;
   } else if(constraint.size() > name.size()) {
      // The constraint is longer than the issued name: not possibly a match
      return false;
   } else {
      // constraint.size() < name.size()
      // constr is suffix of constraint
      const std::string constr = constraint.front() == '.' ? constraint : "." + constraint;
      BOTAN_ASSERT_NOMSG(name.size() >= constr.size());
      const std::string substr = name.substr(name.size() - constr.size(), constr.size());
      return tolower_string(substr) == constr;
   }
}

//static
bool GeneralName::matches_dn(const X509_DN& name, const X509_DN& constraint) {
   const auto attr = name.get_attributes();
   bool ret = true;
   size_t trys = 0;

   for(const auto& c : constraint.dn_info()) {
      auto i = attr.equal_range(c.first);

      if(i.first != i.second) {
         trys += 1;
         ret = ret && (i.first->second == c.second.value());
      }
   }

   return trys > 0 && ret;
}

std::ostream& operator<<(std::ostream& os, const GeneralName& gn) {
   os << gn.type() << ":" << gn.name();
   return os;
}

void GeneralSubtree::encode_into(DER_Encoder& /*to*/) const {
   throw Not_Implemented("General Subtree encoding");
}

void GeneralSubtree::decode_from(BER_Decoder& ber) {
   ber.start_sequence()
      .decode(m_base)
      .decode_optional(m_minimum, ASN1_Type(0), ASN1_Class::ContextSpecific, size_t(0))
      .end_cons();

   if(m_minimum != 0) {
      throw Decoding_Error("GeneralSubtree minimum must be 0");
   }

   m_maximum = std::numeric_limits<std::size_t>::max();
}

std::ostream& operator<<(std::ostream& os, const GeneralSubtree& gs) {
   os << gs.minimum() << "," << gs.maximum() << "," << gs.base();
   return os;
}
}  // namespace Botan
