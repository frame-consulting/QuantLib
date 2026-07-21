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

#include <ql/experimental/credit/tranchesurvivalprobabilitytermstructure.hpp>

#include <cmath>

namespace QuantLib {

    TrancheSurvivalProbabilityTermStructure::TrancheSurvivalProbabilityTermStructure(
        const Date& referenceDate,
        const Calendar& cal,
        const DayCounter& dc,
        const ext::shared_ptr<Basket>& basket)
    : DefaultProbabilityTermStructure(referenceDate, cal, dc), basket_(basket) {
        registerWith(basket_);
    }

    Probability TrancheSurvivalProbabilityTermStructure::survivalProbabilityImpl(Date d) const {
        Real etl = basket_->expectedTrancheLoss(d);
        // some models struggle to calculate ETL for t=0
        if (d == referenceDate() && isnan(etl)) {
            etl = basket_->expectedTrancheLoss(d + 1);
            // better fix the model...
        }
        Real attachmentAmount = basket_->attachmentAmount();
        Real detachmentAmount = basket_->detachmentAmount();
        //
        return 1.0 - etl / (detachmentAmount - attachmentAmount);
    }

    Probability TrancheSurvivalProbabilityTermStructure::survivalProbabilityImpl(Time t) const {
        Real approxDays = 365.0 * t;
        Size days = lround(approxDays);
        QL_REQUIRE(fabs(days - approxDays) < 3.2e-8,  // 1 sec.
                   "Time t must correspond to a day using act/365 day count.");
        return survivalProbabilityImpl(referenceDate() + days);
    }

    Real TrancheSurvivalProbabilityTermStructure::defaultDensityImpl(Time t) const {
        QL_FAIL("defaultDensityImpl not implemented");
    }

    Real TrancheSurvivalProbabilityTermStructure::hazardRateImpl(Time t) const {
        QL_FAIL("hazardRateImpl not implemented");
    }

} // namespace QuantLib