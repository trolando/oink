// SPDX-License-Identifier: Apache-2.0
#ifndef OINK_ENERGY_CREDIT_HPP
#define OINK_ENERGY_CREDIT_HPP

#include <memory>

#include "ep_search.hpp"
#include "gpm.hpp"

namespace pg {

// The alpha table stores the credit needed by the other player. Alpha takes
// the maximum. Top is absorbing and means alpha wins.
class EnergyCredits final : public Measures
{
public:
    EnergyCredits(const Game& game, int player);

    const ep::Integer& value(int slot) const { return values[index(slot)]; }
    ep::Integer progress(const ep::Integer& x, int priority) const;
    bool raise(int slot, const ep::Integer& candidate);

    void copy(int from, int to) override;
    void copy_bot(int to) override;
    void copy_top(int to) override;
    void copy_from(const Measures& other, int from, int to) override;
    void see(int slot, int priority) override;
    void inc(int slot) override;
    int compare(int a, int b) const override;
    bool eq(int a, int b) const override;
    void stream(std::ostream& out, int slot) const override;
    bool is_top(int slot) const override;
    Measures* clone() const override { return new EnergyCredits(*this); }
    void set(const Measures& other) override;
    std::vector<std::byte> store(int slot) const override;
    void load(int slot, const std::vector<std::byte>& bytes) override;
    int code_size() const override { return 1; } // byte digits in a serialized credit

private:
    static std::size_t index(int slot) { return std::size_t(std::int64_t(slot) + 4); }
    struct Domain;
    std::shared_ptr<const Domain> domain;
    ep::Prices values;
};

} // namespace pg

#endif
