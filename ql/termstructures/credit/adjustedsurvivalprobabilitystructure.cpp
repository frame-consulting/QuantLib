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

#include <ql/termstructures/credit/adjustedsurvivalprobabilitystructure.hpp>

namespace QuantLib {

	AdjustedSurvivalProbabilityStructure::AdjustedSurvivalProbabilityStructure(
        const Handle<DefaultProbabilityTermStructure>& defaultProbabilityTermStructureBase,
        const Handle<DefaultProbabilityTermStructure>& adjusterTermStructure)
    : SurvivalProbabilityStructure(defaultProbabilityTermStructureBase->referenceDate(),
                                   defaultProbabilityTermStructureBase->calendar(),
                                   defaultProbabilityTermStructureBase->dayCounter()),
      defaultProbabilityTermStructureBase_(defaultProbabilityTermStructureBase),
      adjusterTermStructure_(ext::dynamic_pointer_cast<adjuster_type>(adjusterTermStructure.currentLink())) {
        QL_REQUIRE(!defaultProbabilityTermStructureBase_.empty(),
                   "defaultProbabilityTermStructureBase_ must not be empty.");
        QL_REQUIRE(adjusterTermStructure_, "adjusterTermStructure must link to a HazardRateCurve.");
        //
        QL_REQUIRE(defaultProbabilityTermStructureBase_->referenceDate() ==
                       adjusterTermStructure_->referenceDate(),
                   "reference dates of defaultProbabilityTermStructureBase and "
                   "adjusterTermStructure must coincide.");
        QL_REQUIRE(adjusterTermStructure_->times().size()>0, "adjusterTermStructure_ must have at least one node");
    }

    Probability AdjustedSurvivalProbabilityStructure::survivalProbabilityImpl(Time t) const {
        QL_REQUIRE(t >= 0.0, "Survival probability calculation requires non-negative time.");
        if (t == 0.0) {
            // do not adjust initial survival probability
            return defaultProbabilityTermStructureBase_->survivalProbability(t);  // should be 1.0, but may differ for defaulted curves
        }
        // find t0 such that t0 = T[k] < t <= T[k+1]
        Time t0 = 0.0;
        for (auto Tk : adjusterTermStructure_->times()) {
            if (Tk < t) {
                t0 = Tk;
            } else {
                break;
            }
        }
        const bool extrapolate = true;
        const Probability Qk = survivalProbabilityImpl(t0); // recursive call
        const Probability Qfwd =
            defaultProbabilityTermStructureBase_->survivalProbability(t, extrapolate) /
            defaultProbabilityTermStructureBase_->survivalProbability(t0, extrapolate);
        const double adjuster = adjusterTermStructure_->hazardRate(t, extrapolate);  // we (mis-)use piece-wise flat hazard rates as adjuster values
        //
        return Qk * pow(Qfwd, adjuster);
    }
}
