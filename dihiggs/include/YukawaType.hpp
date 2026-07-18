#pragma once

#include "THDM.h"

#include <stdexcept>

namespace dihiggs {

inline bool supported_yukawa_type(int type) { return type >= 1 && type <= 4; }

inline void install_yukawa_type(THDM& model, int requested) {
    if (!supported_yukawa_type(requested)) {
        throw std::runtime_error("unsupported_yukawa_type");
    }
    model.set_yukawas_type(requested);
    if (model.get_yukawas_type() != requested) {
        throw std::runtime_error("yukawa_type_installation_mismatch");
    }
}

}  // namespace dihiggs
