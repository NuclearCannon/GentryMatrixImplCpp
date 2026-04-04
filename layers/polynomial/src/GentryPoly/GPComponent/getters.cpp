#include "GentryPoly.hpp"
#include "modops.hpp"

std::vector<int64_t> GPComponent::to_signed() const
{
    std::vector<int64_t> result(data_.size());
    for(int i=0; i<data_.size(); i++)result[i] = mod_centered(data_[i], q_);
    return result;
}