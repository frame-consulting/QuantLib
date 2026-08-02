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

#ifndef quantlib_choe_kwon_default_lossmodel_hpp
#define quantlib_choe_kwon_default_lossmodel_hpp

#include <ql/experimental/credit/defaultlossmodel.hpp>
#include <ql/experimental/math/latentmodel.hpp>

#include <utility>

namespace QuantLib {

    // See Choe/Kwon, The Large Homogeneous Portfolio Approximation..., 2014
    //
    class ChoeKwonDefaultLossModel : public DefaultLossModel, public LatentModel<GaussianCopulaPolicy> {
    public:
        using GaussianCopulaPolicy = copulaType;  // for base correlation model

    private:

        // market inputs
        Handle<Quote> correlation_;
        std::vector<Handle<RecoveryRateQuote>> recoveries_;

        // treat effctively defauted names separately
        Real defaultThreshold_;

        // stochastic correlation model inputs
        Real p_;  // P[ corr == 0 | \eta=0 ]
        Real q_;  // P[ corr == 1 ]

        // random recovery model inputs
        Real mu_; // correlation param. of default and random recovery
        Real b_;  // random recovery input parameter, see eq. (4.1)

        // auxilliary parameters
        struct ModelParameter {
            const Real F;   // (average) default probability
            const Real R;   // (average) recovery rate
            const Real rho; // flat correlation parameter; not market correlation!
            const Real a;   // random recovery function shift
            const Real PhiInv_F;  // \Phi^{-1}(F)
            const Real A;   // attachment point
            const Real B;   // detachment point
        };

        // auxilliary functors
        const static CumulativeNormalDistribution Phi_;
        const static InverseCumulativeNormal PhiInv_;

        // numerical details and parameters

        // boundaries for inverse normal calculation
        const static Real invNormalMin_;
        const static Real invNormalMax_;

        // ensure we can calculate \Phi^{x}
        static Real cut(const Real x) {
            return std::min(invNormalMax_, std::max(invNormalMin_, x));
        }

        // GaussKronrodNonAdaptive parameters
        const static Real absoluteAccuracy_;
        const static Real relativeAccuracy_;
        const static Size maxEvaluations_;

        // GaussHermiteIntegration parameter
        const static Size numberOfIntegrationPoints_;

        // intervall calculations

        using Intervall = std::pair<Real, Real>;  // (min, max)

        // stochastic correlation, constant recovery intervall functions
        Intervall I1(const ModelParameter& p) const;
        Intervall I2(const ModelParameter& p) const;

        // random recovery function
        Real h(Real v, const ModelParameter& p) const; // random recovery function

        // stochastic correlation, random recovery intervall functions
        Intervall J1(Real v, const ModelParameter& o) const;
        Intervall J2(Real v, const ModelParameter& o) const;
        Intervall J3(Real v, const ModelParameter& o) const;
        Intervall J4(Real v, const ModelParameter& o) const;

        // intersection of two intervals
        Intervall intersect(const Intervall& left, const Intervall& right) const;
        // intervall length
        Real length(const Intervall& intervall) const;

        // integrand function

        // constant correlation, constant recovery (standard LHP)
        Real C1(Real x, const ModelParameter& p) const;

        // stochastic correlation, constant recovery
        Real C2(Real x, const ModelParameter& o) const;

        // stochastic correlation, random recovery intervall functions
        Real C3(Real x, Real v, const ModelParameter& p) const;
        Real Psi1(Real v, const ModelParameter& p) const;
        Real Psi2(Real x, Real v, const ModelParameter& p) const;
        Real xi(Real v, const ModelParameter& p) const;
        Real xiConstCorr(Real v, const ModelParameter& p) const;
        Real expectedTrancheLossLhp(const ModelParameter& p) const;
        Real expectedTrancheLossStochasticCorrelation(const ModelParameter& p) const;
        Real expectedTrancheLossRandomRecovery(const ModelParameter& p) const;
        Real expectedTrancheLossStochasticCorrelationRandomRecovery(const ModelParameter& p) const;

        Real expectedTrancheLossImpl(
            Real remainingNotional,
            Real averDefaultProbability,
            Real averRecovery,
            Real attachmentRatio,
            Real detachmentRatio
        ) const;


        // homogeneous pool mapping, see Gaussian LHP model
        Probability averageDefaultProbability(const Date& d) const;
        Real averageRecovery(const Date& d) const;

        // check model inputs at construction
        void checkModelInputs() const;

        // required for DefaultLossModel interface
        void resetModel() override {}

    public:

        ChoeKwonDefaultLossModel(
            const Handle<Quote>& correlation,
            const std::vector<Handle<RecoveryRateQuote>>& recoveries,
            const Real p,
            const Real q,
            const Real mu,
            const Real b,
            const Real defaultThreshold = 1.0-1.0e-8
        );

        ChoeKwonDefaultLossModel(
            Real correlation,
            const std::vector<Real>& recoveries,
            const Real p,
            const Real q,
            const Real mu,
            const Real b,
            const Real defaultThreshold = 1.0-1.0e-8
        );

        ChoeKwonDefaultLossModel(
            const Handle<Quote>& correlation,
            const std::vector<Real>& recoveries,
            const Real p,
            const Real q,
            const Real mu,
            const Real b,
            const Real defaultThreshold = 1.0-1.0e-8
        );

        // Observer interface via LatentModel, see LHP model
        void update() override {
            if (!basket_.empty()) {
                basket_->notifyObservers();
            }
        }

        // actual DefaultLossModel interface
        Real expectedTrancheLoss(const Date& d) const override;


    }; // class ChoeKwonDefaultLossModel



} // namespace QuantLib

#endif  // quantlib_choe_kwon_default_lossmodel_hpp
