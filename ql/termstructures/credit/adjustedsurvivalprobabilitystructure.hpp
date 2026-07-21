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

/*! \file adjustedsurvivalprobabilitystructure.hpp
    \brief A DefaultProbabilityTermStructure with a multiplicative adjustment of hazard rates.
*/

#ifndef quantlib_adjusted_survival_probability_structure_hpp
#define quantlib_adjusted_survival_probability_structure_hpp

#include <ql/math/interpolations/backwardflatinterpolation.hpp>
#include <ql/termstructures/credit/survivalprobabilitystructure.hpp>
#include <ql/termstructures/credit/interpolatedhazardratecurve.hpp>


namespace QuantLib {

    //! Flat hazard-rate curve
    /*! \ingroup defaultprobabilitytermstructures */
    class AdjustedSurvivalProbabilityStructure : public SurvivalProbabilityStructure {
      public:
        typedef InterpolatedHazardRateCurve<BackwardFlat> adjuster_type;

        //! \name Constructors
        //@{
        AdjustedSurvivalProbabilityStructure(
            const Handle<DefaultProbabilityTermStructure>& defaultProbabilityTermStructureBase,
            const Handle<DefaultProbabilityTermStructure>& adjusterTermStructure);

        //@}
        //! \name TermStructure interface
        //@{
        Date maxDate() const override { return defaultProbabilityTermStructureBase_->maxDate(); }
        //@}
      private:
        //! \name SurvivalProbabilityStructure interface
        //@{
        Probability survivalProbabilityImpl(Time) const override;
        //@}

        const Handle<DefaultProbabilityTermStructure> defaultProbabilityTermStructureBase_;
        const ext::shared_ptr<adjuster_type> adjusterTermStructure_;
    };

}

#endif
