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

#include <ql/experimental/credit/choekwondefaultlossmodel.hpp>
#include <ql/math/integrals/kronrodintegral.hpp>
#include <ql/math/integrals/gaussianquadratures.hpp>
#include <ql/quotes/simplequote.hpp>

namespace QuantLib {

    // cumulative normal distribution functors
    const CumulativeNormalDistribution ChoeKwonDefaultLossModel::Phi_ = CumulativeNormalDistribution();
    const InverseCumulativeNormal ChoeKwonDefaultLossModel::PhiInv_ = InverseCumulativeNormal();

    // boundaries for inverse normal calculation
    const Real ChoeKwonDefaultLossModel::invNormalMin_ = 1.0e-12;
    const Real ChoeKwonDefaultLossModel::invNormalMax_ = 1.0 - 1.0e-12;

    // GaussKronrodNonAdaptive parameters
    const Real ChoeKwonDefaultLossModel::absoluteAccuracy_ = 1.0e-8;
    const Real ChoeKwonDefaultLossModel::relativeAccuracy_ = 1.0e-8;
    const Size ChoeKwonDefaultLossModel::maxEvaluations_ = 100;

    // GaussHermiteIntegration parameter
    const Size ChoeKwonDefaultLossModel::numberOfIntegrationPoints_ = 10;


    ChoeKwonDefaultLossModel::Intervall ChoeKwonDefaultLossModel::I1(const ModelParameter& p) const {
        Real l = (1.0 - p.R) * p_ * p.F;
        Real u = (1.0 - p.R) * (1.0 - p_ + p_ * p.F);
        return {l, u};
    }

    ChoeKwonDefaultLossModel::Intervall ChoeKwonDefaultLossModel::I2(const ModelParameter& p) const {
        Real l = (1.0 - p.R) * (1.0 - p_ + p_ * p.F);
        Real u = (1.0 - p.R);
        return {l, u};
    }

    Real ChoeKwonDefaultLossModel::h(Real v, const ModelParameter& p) const {
        return Phi_(p.a + b_ * v);
    }

    ChoeKwonDefaultLossModel::Intervall ChoeKwonDefaultLossModel::J1(Real v, const ModelParameter& p) const {
        Real l = (1.0 - h(v, p)) * p_ * p.F;
        Real u = (1.0 - h(v, p)) * (1.0 - p_ + p_ * p.F);
        return {l, u};
    }

    ChoeKwonDefaultLossModel::Intervall ChoeKwonDefaultLossModel::J2(Real v, const ModelParameter& p) const {
        Real l = (1.0 - h(v, p)) * (1.0 - p_ + p_ * p.F);
        Real u = (1.0 - h(v, p));
        return {l, u};
    }

    ChoeKwonDefaultLossModel::Intervall ChoeKwonDefaultLossModel::J3(Real v, const ModelParameter& p) const {
        Real l = 0.0;
        Real u = (1.0 - h(v, p));
        return {l, u};
    }

    ChoeKwonDefaultLossModel::Intervall ChoeKwonDefaultLossModel::J4(Real v, const ModelParameter& p) const {
        Real l = (1.0 - h(v, p));
        Real u = 1.0;
        return {l, u};
    }

    ChoeKwonDefaultLossModel::Intervall ChoeKwonDefaultLossModel::intersect(const Intervall& left, const Intervall& right) const {
        Real l = std::max(left.first, right.first);
        Real u = std::min(left.second, right.second);
        if (l <= u) {
            return {l, u};
        }
        // we pick the mid-point between non-intersecting intervalls
        Real m = 0.5 * (u + l);
        return {m, m}; // or better throw
    }

    Real ChoeKwonDefaultLossModel::length(const Intervall& intervall) const {
        return intervall.second - intervall.first;
    }

    Real ChoeKwonDefaultLossModel::C1(Real x, const ModelParameter& p) const {
        Real t1 = x / (1.0 - p.R);
        Real t2 = PhiInv_(cut(t1));
        Real t3 = sqrt(1.0 - p.rho * p.rho) * t2;
        return 1.0 / p.rho * (p.PhiInv_F - t3);
    }

    Real ChoeKwonDefaultLossModel::C2(Real x, const ModelParameter& p) const {
        Real t1 = x / (1.0 - p.R) - p_ * p.F;
        Real t2 = 1.0 / (1.0 - p_) * t1;
        Real t3 = sqrt(1.0 - p.rho * p.rho) * PhiInv_(cut(t2));
        return 1.0 / p.rho * (p.PhiInv_F - t3);
    }

    Real ChoeKwonDefaultLossModel::C3(Real x, Real v, const ModelParameter& p) const {
        Real t1 = x / (1.0 - h(v, p)) - p_ * p.F;
        Real t2 = 1.0 / (1.0 - p_) * t1;
        Real t3 = sqrt(1.0 - p.rho * p.rho) * PhiInv_(cut(t2));
        return 1.0 / p.rho * (p.PhiInv_F - t3);
    }

    Real ChoeKwonDefaultLossModel::Psi1(Real v, const ModelParameter& p) const {
        return Phi_((mu_ * v - p.PhiInv_F) / sqrt(1.0 - mu_ * mu_));
    }

    Real ChoeKwonDefaultLossModel::Psi2(Real x, Real v, const ModelParameter& p) const {
        return Phi_((mu_ * v - C3(x, v, p)) / sqrt(1.0 - mu_ * mu_));
    }

    Real ChoeKwonDefaultLossModel::xi(Real v, const ModelParameter& p) const {
        Intervall AB = {p.A, p.B};
        Real t2 = (1.0 - q_) * length(intersect(AB, J2(v, p)));
        Real t3 = q_ * Psi1(v, p) * length(intersect(AB, J3(v, p)));
        Real t4 = length(intersect(AB, J4(v, p)));
        //
        Intervall LU = intersect(AB, J1(v, p));
        auto Psi2_v = [v, p, this](Real x) { return Psi2(x, v, p); };
        GaussKronrodNonAdaptive I(absoluteAccuracy_, maxEvaluations_, relativeAccuracy_);
        Real t1 = (1.0 - q_) * I(Psi2_v, LU.first, LU.second);
        //
        return t1 + t2 + t3 + t4;
    }

    Real ChoeKwonDefaultLossModel::xiConstCorr(Real v, const ModelParameter& p) const {
        Intervall AB = {p.A, p.B};
        Intervall LU = intersect(AB, J3(v, p));
        auto Psi2_v = [v, p, this](Real x) { return Psi2(x, v, p); };
        GaussKronrodNonAdaptive I(absoluteAccuracy_, maxEvaluations_, relativeAccuracy_);
        Real integral = I(Psi2_v, LU.first, LU.second);
        Real l = length(intersect(AB, J4(v, p)));
        return integral + l;
    }

    Real ChoeKwonDefaultLossModel::expectedTrancheLossLhp(const ModelParameter& p) const {
        auto Phi_C1 = [p, this](Real x) { return Phi_(C1(x, p)); };
        GaussKronrodNonAdaptive I(absoluteAccuracy_, maxEvaluations_, relativeAccuracy_);
        Real integral = I(Phi_C1, p.A, p.B);
        Real etl = 1.0 / (p.B - p.A) * integral;
        return etl;
    }

    Real ChoeKwonDefaultLossModel::expectedTrancheLossStochasticCorrelation(const ModelParameter& p) const {
        if (p.A >= 1 - p.R) {
            // no loss possible due to deterministic recovery
            return 0.0;
        }
        if (p.A < 1 - p.R && 1 - p.R < p.B) {
            // we need to split the integral from [A, 1-R] and [1-R, B]
            ModelParameter p1 = {p.F, p.R, p.rho, p.a, p.PhiInv_F, p.A, 1 - p.R};
            Real t1 = expectedTrancheLossStochasticCorrelation(p1);  // E[ L(A, 1-R) ]
            Real t2 = (1 - p.R - p.A) / (p.B - p.A);
            Real t3 = (p.B - (1 - p.R)) / (p.B - p.A);
            return 1.0 - (t2 * (1.0 - t1) + t3);
        }
        // A < B <= 1 - R
        Real t1 = (1.0 - q_) / (p.B - p.A);
        //
        Intervall AB = {p.A, p.B};
        Intervall LU = intersect(AB, I1(p));
        auto Phi_minus_C2 = [p, this](Real x) { return Phi_(-C2(x, p)); };
        GaussKronrodNonAdaptive I(absoluteAccuracy_, maxEvaluations_, relativeAccuracy_);
        Real integral = I(Phi_minus_C2, LU.first, LU.second);
        //
        Real t2 = integral + length(intersect(AB, I2(p)));
        Real t3 = q_ * (1.0 - p.F);
        Real etl = 1.0 - t1 * t2 - t3;
        return etl;
    }

    Real ChoeKwonDefaultLossModel::expectedTrancheLossRandomRecovery(const ModelParameter& p) const {
        auto xi_AB = [p, this](Real u) {
            const Real v = M_SQRT2 * u;
            return xiConstCorr(v, p) * M_1_SQRTPI * exp(-u * u);
        };
        GaussHermiteIntegration I(numberOfIntegrationPoints_);
        Real integral = I(xi_AB);
        Real etl = 1.0 - 1.0 / (p.B - p.A) * integral;
        return etl;
    }

    Real ChoeKwonDefaultLossModel::expectedTrancheLossStochasticCorrelationRandomRecovery(const ModelParameter& p) const {
        auto xi_AB = [p, this](Real u) {
            const Real v = M_SQRT2 * u;
            return xi(v, p) * M_1_SQRTPI * exp(-u * u);
        };
        GaussHermiteIntegration I(numberOfIntegrationPoints_);
        Real integral = I(xi_AB);
        Real etl = 1.0 - 1.0 / (p.B - p.A) * integral;
        return etl;
    }


    Real ChoeKwonDefaultLossModel::expectedTrancheLossImpl(
        Real remainingNotional,
        Real averDefaultProbability,
        Real averRecovery,
        Real attachmentRatio,
        Real detachmentRatio) const {
        //
        Real F = averDefaultProbability;
        Real R = averRecovery;
        Real rho = std::sqrt(correlation_->value());
        Real a = PhiInv_(cut(R)) * sqrt(1.0 + b_*b_);
        Real PhiInv_F = PhiInv_(cut(F));
        Real A = attachmentRatio;
        Real B = detachmentRatio;
        // TODO: double-check correlated case, \mu != 0
        QL_REQUIRE(0.0 <= F && F <= 1.0, "0.0 <= F <= 1.0 required.");
        QL_REQUIRE(0 <= R && R < 1.0, "0<= R <1.0 required.");
        QL_REQUIRE(0<rho && rho<=1.0, "0<rho && rho<=1.0 required.");
        QL_REQUIRE(0<=A && A<=1.0, "0<= A<=1.0 required." );
        QL_REQUIRE(0<=B && B<=1.0, "0<= B<=1.0 required." );
        QL_REQUIRE(A < B, "A < B required.");
        QL_REQUIRE(remainingNotional > 0.0, "remainingNotional > 0.0 required.");
        //
        ModelParameter p = {F, R, rho, a, PhiInv_F, A, B};
        Real etl = 0.0;
        if (close(p_, 0.0) && close(q_, 0.0) && close(mu_, 0.0) && close(b_, 0.0)) {
            etl = expectedTrancheLossLhp(p);
        } else if (close(mu_, 0.0) && close(b_, 0.0)) {
            etl = expectedTrancheLossStochasticCorrelation(p);
        } else if (close(p_, 0.0) && close(q_, 0.0)) {
            etl = expectedTrancheLossRandomRecovery(p);
        }else {
            etl = expectedTrancheLossStochasticCorrelationRandomRecovery(p);
        }
        return remainingNotional * (B - A) * etl;

    }



    Probability ChoeKwonDefaultLossModel::averageDefaultProbability(const Date& d) const {
        const std::vector<Probability> probabilities = basket_->remainingProbabilities(d);
        const std::vector<Real> notionals = basket_->remainingNotionals(d);
        Real sum_probs = 0.0;
        Real sum_notionals = 0.0;
        for (Size k = 0; k < probabilities.size(); ++k) {
            if (probabilities[k] < defaultThreshold_) {
                sum_probs += probabilities[k] * notionals[k];
                sum_notionals += notionals[k];
            }
        }
        if (close(sum_notionals,0.0)) {
            return 0.0;
        }
        return sum_probs / sum_notionals;
    }

    Probability ChoeKwonDefaultLossModel::averageRecovery(const Date& d) const {
        const std::vector<Probability> probabilities = basket_->remainingProbabilities(d);
        std::vector<Real> recoveries;
        for (Size k=0; k<recoveries_.size(); ++k) {
            recoveries.push_back(recoveries_[k]->value());
        }
        std::vector<Real> notionals = basket_->remainingNotionals(d);
        Real numerator = 0.0;
        Real denominator = 0.0;
        for (Size k=0; k<notionals.size(); ++k) {
            if (probabilities[k] < defaultThreshold_) {
                numerator += recoveries[k] * notionals[k] * probabilities[k];
                denominator += notionals[k] * probabilities[k];
            }
        }
        if (close(denominator, 0.0)) {
            return 0.0;
        }
        return numerator / denominator;
    }

    void ChoeKwonDefaultLossModel::checkModelInputs() const {
        QL_REQUIRE(0.0 <= correlation_->value() && correlation_->value() <= 1.0, "0.0 <= correlation <= 1.0 required.");
        for (const auto& recovery : recoveries_) {
            QL_REQUIRE(0.0 <= recovery->value() && recovery->value() <= 1.0, "0.0 <= recovery <= 1.0 required.");
        }
        //
        QL_REQUIRE(0.0 <= p_ && p_ < 1.0, "0.0 <= p_ < 1.0 required.");
        QL_REQUIRE(0.0 <= q_ && q_ <= 1.0, "0.0 <= q_ <= 1.0 required.");
        QL_REQUIRE(-1.0 < mu_ && mu_ < 1.0, "-1.0 < mu_ < 1.0 required.");
        QL_REQUIRE(0.0 <= b_, "0.0 <= b_ required.");
        //
        QL_REQUIRE(0.0 < defaultThreshold_ && defaultThreshold_ <= 1.0, "0.0 < defaultThreshold_ <= 1.0 required.");
    }

    ChoeKwonDefaultLossModel::ChoeKwonDefaultLossModel(
        const Handle<Quote>& correlation,
        const std::vector<Handle<RecoveryRateQuote>>& recoveries,
        const Real p,
        const Real q,
        const Real mu,
        const Real b,
        const Real defaultThreshold)
    : LatentModel<GaussianCopulaPolicy>(sqrt(correlation->value()), recoveries.size(), GaussianCopulaPolicy::initTraits()),
      correlation_(correlation), recoveries_(recoveries), p_(p), q_(q), mu_(mu), b_(b),
      defaultThreshold_(defaultThreshold)
    {
        checkModelInputs();
        registerWith(correlation_);
        for (const auto& recovery : recoveries_) {
            registerWith(recovery);
        }
    }


    ChoeKwonDefaultLossModel::ChoeKwonDefaultLossModel(
        Real correlation,
        const std::vector<Real>& recoveries,
        const Real p,
        const Real q,
        const Real mu,
        const Real b,
        const Real defaultThreshold)
    : LatentModel<GaussianCopulaPolicy>(
          sqrt(correlation), recoveries.size(), GaussianCopulaPolicy::initTraits()),
      correlation_(Handle<Quote>(ext::make_shared<SimpleQuote>(correlation))),
      recoveries_(), p_(p), q_(q), mu_(mu), b_(b),
      defaultThreshold_(defaultThreshold)
    {
        for (auto recovery : recoveries) {
            recoveries_.emplace_back(ext::make_shared<RecoveryRateQuote>(recovery));
        }
        checkModelInputs();
        registerWith(correlation_);
    }


    ChoeKwonDefaultLossModel::ChoeKwonDefaultLossModel(
        const Handle<Quote>& correlation,
        const std::vector<Real>& recoveries,
        const Real p,
        const Real q,
        const Real mu,
        const Real b,
        const Real defaultThreshold)
    : LatentModel<GaussianCopulaPolicy>(
          sqrt(correlation->value()), recoveries.size(), GaussianCopulaPolicy::initTraits()),
      correlation_(correlation), recoveries_(), p_(p), q_(q), mu_(mu), b_(b),
      defaultThreshold_(defaultThreshold) 
    {
        for (auto recovery : recoveries) {
            recoveries_.emplace_back(ext::make_shared<RecoveryRateQuote>(recovery));
        }
        checkModelInputs();
        registerWith(correlation_);
    }


    Real ChoeKwonDefaultLossModel::expectedTrancheLoss(const Date& d) const {
        Real basketNotional = basket_->remainingNotional(d);
        Probability averDefaultProbability = averageDefaultProbability(d);
        Real averRecovery = averageRecovery(d);
        Real attachmentAmount = basket_->remainingAttachmentAmount();
        Real detachmentAmount = basket_->remainingDetachmentAmount();
        Real attachmentRatio = attachmentAmount / basketNotional;
        Real detachmentRatio = detachmentAmount / basketNotional;
        //
        const std::vector<Probability> probabilities = basket_->remainingProbabilities(d);
        const std::vector<Real> notionals = basket_->remainingNotionals(d);
        //
        Real lossNotional = 0.0;
        Real realisedLoss = 0.0;
        for (Size k = 0; k < probabilities.size(); ++k) {
            if (probabilities[k] >= defaultThreshold_) {
                lossNotional += notionals[k];
                realisedLoss += (1.0 - recoveries_[k]->value()) * notionals[k];
            }
        }
        Real adjAttachRatio = 0.0;
        Real adjDetachRatio = 0.0;
        if (lossNotional < basketNotional) {
            // typical use case
            adjAttachRatio = (attachmentAmount - realisedLoss) / (basketNotional - realisedLoss);
            adjDetachRatio = (detachmentAmount - realisedLoss) / (basketNotional - realisedLoss);
        }
        else { // all names defaulted, loss determined by recoveries
            return std::min(realisedLoss, detachmentAmount) - std::min(realisedLoss, attachmentAmount);
        }
        if (realisedLoss <= attachmentAmount) {
            // case one
            return expectedTrancheLossImpl(basketNotional-lossNotional, averDefaultProbability, averRecovery, adjAttachRatio, adjDetachRatio);
        }
        if (realisedLoss <= detachmentAmount) {
            // case two
            return realisedLoss - attachmentAmount + expectedTrancheLossImpl(basketNotional-lossNotional, averDefaultProbability, averRecovery, 0.0, adjDetachRatio);
        }
        // full realised loss in tranche
        return detachmentAmount - attachmentAmount;
    }

} // namespace QuantLib