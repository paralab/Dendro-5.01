#pragma once

#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

#include "base.h"

namespace dendro_aeh {

// base91 (de)serialization of a flat POD vector, for compact, high-accuracy
// checkpointing. Shared by the AH finder and the BHHistory tracker.
template <typename T>
inline std::string b91_encode(const std::vector<T>& vec) {
    return base<91>::encode(std::string(
        reinterpret_cast<const char*>(vec.data()), vec.size() * sizeof(T)));
}

template <typename T>
inline std::vector<T> b91_decode(const std::string& input) {
    std::string temp         = base<91>::decode(input);
    const size_t num_entries = temp.size() / sizeof(T);
    std::vector<T> output(num_entries);

    std::memcpy(output.data(), temp.data(), sizeof(T) * num_entries);

    return output;
}

template <typename T>
inline void restore_vector(std::vector<T>& original,
                           const std::vector<T>& restored,
                           const size_t expected_size) {
    if (restored.size() != expected_size) {
        throw std::runtime_error(
            "ERROR when restoring a vector! The size is "
            "different from the expected size!");
    }

    // otherwise just copy straight in
    // NOTE: it MUST copy, due to pointer information for prev_horizon
    for (size_t i = 0; i < expected_size; ++i) {
        original[i] = restored[i];
    }
}

}  // namespace dendro_aeh
