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


#ifndef quantlib_tranche_survival_probability_termstructure_hpp
#define quantlib_tranche_survival_probability_termstructure_hpp

#include <ql/experimental/credit/basket.hpp>
#include <ql/termstructures/defaulttermstructure.hpp>

namespace QuantLib {

    class TrancheSurvivalProbabilityTermStructure : public DefaultProbabilityTermStructure {
      public:
        TrancheSurvivalProbabilityTermStructure(
            const Date& referenceDate,
            const Calendar& cal,
            const DayCounter& dc,
            const ext::shared_ptr<Basket>& basket);

        Date maxDate() const override { return Date::maxDate(); }
      private:
        Probability survivalProbabilityImpl(Date d) const;

      protected:

        Probability survivalProbabilityImpl(Time) const override;
        Real defaultDensityImpl(Time) const override;
        Real hazardRateImpl(Time) const override;

      private:
        ext::shared_ptr<Basket> basket_;
    };
}

#endif
