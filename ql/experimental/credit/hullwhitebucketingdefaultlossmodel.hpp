/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <https://www.quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#ifndef quantlib_hull_white_bucketing_default_loss_model_hpp
#define quantlib_hull_white_bucketing_default_loss_model_hpp

#include <ql/experimental/credit/lossdistribution.hpp>
#include <ql/experimental/credit/basket.hpp>
#include <ql/experimental/credit/constantlosslatentmodel.hpp>
#include <ql/experimental/credit/defaultlossmodel.hpp>

#include <numeric>

namespace QuantLib {


    template<class copula_policy_type>
    class HullWhiteBucketingDefaultLossModel : public DefaultLossModel {
    public:
        enum UpperBoundStrategy { One, Detachment };

    private:
        const ext::shared_ptr<ConstantLossLatentmodel<copula_policy_type>> copula_;
        Real recoveryScaling_;  // see Amraoui/Hitier, 2008
        Size nBuckets_;  // see Hull/White, 2004
        Real max_;
        Real min_;
        Size nSteps_;
        Real delta_;
        UpperBoundStrategy upperBoundStrategy_;
        bool enforceDistributionProperty_;


        void resetModel() override {
            copula_->resetBasket(basket_.currentLink());
        }

        // a helper implementing the minimal required interface of a Distribution
        struct BucketedDistribution{
            std::vector<Real> probabilities;
            std::vector<Real> means;

            Real cumulativeExcessProbability(const Real a, const Real d) const {
                Real prob = 0.0;
                for (Size k = 0; k < probabilities.size(); ++k) {
                    prob += (std::max(means[k] - a, 0.0) - std::max(means[k] - d, 0.0)) * probabilities[k];
                }
                return prob;
            }
        };

    public:

        typedef copula_policy_type copulaType;  // for base correlation

        HullWhiteBucketingDefaultLossModel(
            const ext::shared_ptr<ConstantLossLatentmodel<copula_policy_type> >& copula,
            const Real recoveryScaling,
            const Size nBuckets,
            const Real max = 5.,
            const Real min = -5.,
            const Size nSteps = 50,
            const UpperBoundStrategy upperBoundStrategy = One,
            const bool enforceDistributionProperty = false
        ):
            copula_(copula),
            recoveryScaling_(recoveryScaling),
            nBuckets_(nBuckets),
            max_(max),
            min_(min),
            nSteps_(nSteps),
            delta_((max - min) / nSteps),
            upperBoundStrategy_(upperBoundStrategy),
            enforceDistributionProperty_(enforceDistributionProperty)
        {
            QL_REQUIRE(copula->numFactors() == 1, "copula->numFactors() = " << copula->numFactors() << " != 1.");
            QL_REQUIRE(recoveryScaling_ > 0.0, "recoveryScaling_ = " << recoveryScaling_ << " <= 0.0");
            QL_REQUIRE(recoveryScaling_ <= 1.0, "recoveryScaling_ = " << recoveryScaling_ << " > 1.0");
        }

        // required DefaultLossModel interface for Baskets and tranche pricing
        Real expectedTrancheLoss(const Date& d) const override {
            // follow the pattern of InhomogeneousPoolLossModel
            return lossDistrib(d).cumulativeExcessProbability(basket_->attachmentAmount(), basket_->detachmentAmount());
        }

    protected:

        BucketedDistribution lossDistrib(const Date& d) const {
            // collect inputs
            std::vector<Real> notionals = basket_->notionals();
            std::vector<Real> asset_probabilities = basket_->probabilities(d);
            std::vector<Real> recoveries = copula_->recoveries();

            QL_REQUIRE(asset_probabilities.size() == notionals.size(), "asset_probabilities mismatch.");
            QL_REQUIRE(recoveries.size() == notionals.size(), "recoveries.size() mismatch.");

            // specify \tilde R_i via recovery markdown
            std::vector<Real> adjusted_recoveries(recoveries.size());
            for (Size i = 0; i < adjusted_recoveries.size(); ++i) {
                adjusted_recoveries[i] = recoveries[i] * recoveryScaling_;
            }

            // calculate \tilde p_i such that expected loss is preserved
            // see Amraoui/Hitier, 2008, p. 6, bottom
            std::vector<Real> adjusted_asset_probabilities(asset_probabilities.size());
            for (Size i = 0; i < adjusted_asset_probabilities.size(); ++i) {
                adjusted_asset_probabilities[i] =
                    (1.0 - recoveries[i]) / (1.0 - adjusted_recoveries[i]) * asset_probabilities[i];
            }

            // z(p_i) via copula inverse cumulative distribution function
            std::vector<Real> default_thesholds(asset_probabilities.size());
            for (Size i = 0; i < default_thesholds.size(); ++i) {
                default_thesholds[i] = copula_->inverseCumulativeY(asset_probabilities[i], i);
            }

            // z(\tilde p_i)
            std::vector<Real> adjusted_default_thesholds(adjusted_asset_probabilities.size());
            for (Size i = 0; i < adjusted_default_thesholds.size(); ++i) {
                adjusted_default_thesholds[i] =
                    copula_->inverseCumulativeY(adjusted_asset_probabilities[i], i);
            }

            // initialize the distribution
            BucketedDistribution dist;
            dist.probabilities = std::vector<Real>(nBuckets_ + 2, 0.0);
            dist.means = std::vector<Real>(nBuckets_ + 2, 0.0);

            Real upperBoundBuckets;
            if (upperBoundStrategy_ == One) {
                upperBoundBuckets = std::accumulate(notionals.begin(), notionals.end(), 0.0);
                // keep bucketing unchanged independent of tranche
            } else if (upperBoundStrategy_ == Detachment) {
                upperBoundBuckets = basket_->detachmentAmount();
                // changes bucketing depending on tranche
                // see comment on 'u' in Hull/White
            } else {
                QL_FAIL("Unknown upperBoundStrategy_");
            }

            // integrate over the market factor
            for (Size j = 0; j < nSteps_; ++j) {
                // place the market factor in the middle of the interval
                // note that the copula interface takes a vector of factors, so we need to wrap the market factor in a vector
                std::vector<Real> market_factor(1, min_ + (j + 0.5) * delta_);
                //
                std::vector<Real> conditional_probabilities(notionals.size());
                std::vector<Real> adjusted_lgds(notionals.size());
                for (Size i = 0; i < notionals.size(); ++i) {
                    // Amraoui/Hitier, 2008, eq. (2)
                    // assume all names are alive, no actual defaults
                    conditional_probabilities[i] = copula_->conditionalDefaultProbabilityInvP(default_thesholds[i], i, market_factor);
                    Real adjusted_conditional_probability = copula_->conditionalDefaultProbabilityInvP(adjusted_default_thesholds[i], i, market_factor);
                    adjusted_lgds[i] = notionals[i] * (1.0 - adjusted_recoveries[i]) * adjusted_conditional_probability / conditional_probabilities[i];
                    // adjusted_lgds is something like 'local LGD'; random but deterministic link to market factor
                }

                // calculate bucketing for adjusted default thresholds and adjusted LGDs

                // see Hull/White, 2004, Appendix B

                std::vector<Real> b(nBuckets_ + 2, 0.0); // lower boundary of buckets, does not change; leave here for clarity
                std::vector<Real> a(nBuckets_ + 2, 0.0); // mean loss in buckets
                std::vector<Real> q(nBuckets_ + 2, 0.0); // probability mass in buckets

                b[0] = 0.0;
                b[1] = 0.0;
                b.back() = upperBoundBuckets;

                Real dx = upperBoundBuckets / nBuckets_;
                for (Size k = 2; k < b.size() - 1; ++k) {
                    b[k] = (k - 1) * dx;
                }
                // b = [0, 0, dx, 2*dx, ..., (nBuckets-1)*dx, detachmentAmount]

                for (Size k = 0; k < a.size() - 1; ++k) {
                    a[k] = 0.5 * (b[k] + b[k + 1]);
                }
                a.back() = b.back();
                // a = [0, 0.5*dx, 1.5*dx, ..., (nBuckets-0.5)*dx, detachmentAmount]

                q[0] = 1.0; // start with no losses, all probability mass in the first bucket
                for (Size i = 0; i < notionals.size(); i++) {
                    Real L = adjusted_lgds[i]; // abbreviate for readability
                    Real P = conditional_probabilities[i];
                    for (Size k = b.size() - 1; k > 0; --k) {
                        // note, we cannot use condition k>= 0 due to Size beeing unsigned
                        // also, do not update the last bucket
                        if (close(q[k - 1], 0.0)) {
                            continue; // no probability to shift
                        }
                        if (close(P, 0.0) || close(L, 0.0)) {
                            // better use some epsilon her to avoid overflow below
                            continue; // no shift in probabilities and losses
                        }
                        QL_REQUIRE(P * L > 0.0, "P * L > 0.0 required, P:" << P << ", L: " << L); // should not throw
                        // find u = u(k) index
                        Size u = b.size() - 1;
                        while (a[k - 1] + L < b[u]) {
                            --u;
                            if (u == 0) {
                                break; // avoid running over index boundaries
                            }
                        }
                        QL_REQUIRE(u >= k - 1, "u >= k - 1 required");  // should not throw
                        if (u == k - 1) { // simple case
                            a[k - 1] += P * L;
                        } else {
                            // shift mean loss
                            Real f = 1.0 / (1.0 + (q[u] / q[k - 1]) / P);
                            a[u] = (1.0 - f) * a[u] + f * (a[k - 1] + L);
                            // shift probabilities
                            q[u] += P * q[k - 1];
                            q[k - 1] -= P * q[k - 1];
                        }
                    }
                }

                // finally build conditional distribution
                Real densitydm = delta_ * copula_->density(market_factor);
                for (Size k = 0; k < nBuckets_ + 2; ++k) {
                    dist.probabilities[k] += q[k] * densitydm;
                    dist.means[k] += a[k] * densitydm;
                }
            }  // market factor

            if (enforceDistributionProperty_) {
                // helps for [0, 1] tranche, but sum of tranches may not match
                //
                // Ensure E[1] == 1
                Real sum_probs = 0.0;
                for (Size k = 0; k < dist.probabilities.size(); ++k)
                    sum_probs += dist.probabilities[k];
                for (Size k = 0; k < dist.probabilities.size(); ++k)
                    dist.probabilities[k] *= 1.0 / sum_probs;
                // Ensure E[ E[ L ] | bucket ] == E[ L ]
                Real exp_loss_dist = 0.0;
                for (Size k = 0; k < dist.probabilities.size(); ++k) {
                    exp_loss_dist += dist.probabilities[k] * dist.means[k];
                }
                Real exp_loss_idx = 0.0;
                for (Size i = 0; i < notionals.size(); ++i) {
                    exp_loss_idx += notionals[i] * (1.0 - recoveries[i]) * asset_probabilities[i];
                } // apply multiplicative scaling
                for (Size k = 0; k < dist.probabilities.size(); ++k) {
                    dist.means[k] *= exp_loss_idx / exp_loss_dist;
                }
            }


            return dist;
        } // lossDistrib

    }; // HullWhiteBucketingDefaultLossModel

    typedef HullWhiteBucketingDefaultLossModel<GaussianCopulaPolicy> GaussianHullWhiteBucketingDefaultLossModel;

}  // namespace QuantLib

#endif
