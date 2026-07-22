/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2008, 2009 Jose Aparicio
 Copyright (C) 2008 Chris Kenyon
 Copyright (C) 2008 Roland Lichters
 Copyright (C) 2008 StatPro Italia srl
  Copyright (C) 2023 Andrea Pellegatta

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

#include <ql/pricingengines/credit/isdacdsengine.hpp>
#include <ql/pricingengines/credit/midpointcdsengine.hpp>
#include <ql/termstructures/credit/adjustedsurvivalprobabilitystructure.hpp>
#include <ql/termstructures/credit/defaultprobabilityhelpers.hpp>
#include <ql/utilities/null_deleter.hpp>
#include <numeric>
#include <utility>

namespace QuantLib {

    CdsHelper::CdsHelper(const std::variant<Rate, Handle<Quote>>& quote,
                         const Period& tenor,
                         Integer settlementDays,
                         Calendar calendar,
                         Frequency frequency,
                         BusinessDayConvention paymentConvention,
                         DateGeneration::Rule rule,
                         DayCounter dayCounter,
                         Real recoveryRate,
                         const Handle<YieldTermStructure>& discountCurve,
                         bool settlesAccrual,
                         bool paysAtDefaultTime,
                         const Date& startDate,
                         DayCounter lastPeriodDayCounter,
                         const bool rebatesAccrual,
                         const CreditDefaultSwap::PricingModel model)
    : RelativeDateDefaultProbabilityHelper(quote), tenor_(tenor), settlementDays_(settlementDays),
      calendar_(std::move(calendar)), frequency_(frequency), paymentConvention_(paymentConvention),
      rule_(rule), dayCounter_(std::move(dayCounter)), recoveryRate_(recoveryRate),
      discountCurve_(discountCurve), settlesAccrual_(settlesAccrual),
      paysAtDefaultTime_(paysAtDefaultTime), lastPeriodDC_(std::move(lastPeriodDayCounter)),
      rebatesAccrual_(rebatesAccrual), model_(model), startDate_(startDate) {

        CdsHelper::initializeDates();

        registerWith(discountCurve);
    }

    void CdsHelper::setTermStructure(DefaultProbabilityTermStructure* ts) {
        RelativeDateDefaultProbabilityHelper::setTermStructure(ts);

        probability_.linkTo(
            ext::shared_ptr<DefaultProbabilityTermStructure>(ts, null_deleter()),
            false);

        resetEngine();
    }

    void CdsHelper::update() {
        RelativeDateDefaultProbabilityHelper::update();
        resetEngine();
    }

    void CdsHelper::initializeDates() {

        protectionStart_ = evaluationDate_ + settlementDays_;

        Date startDate = startDate_ == Date() ? protectionStart_ : startDate_;
        // Only adjust start date if rule is not CDS or CDS2015. Unsure about OldCDS.
        if (rule_ != DateGeneration::CDS && rule_ != DateGeneration::CDS2015) {
            startDate = calendar_.adjust(startDate, paymentConvention_);
        }

        Date endDate;
        if (rule_ == DateGeneration::CDS2015 || rule_ == DateGeneration::CDS || rule_ == DateGeneration::OldCDS) {
            Date refDate = startDate_ == Date() ? evaluationDate_ : startDate_;
            endDate = cdsMaturity(refDate, tenor_, rule_);
        } else {
            // Keep the old logic here
            Date refDate = startDate_ == Date() ? protectionStart_ : startDate_ + settlementDays_;
            endDate = refDate + tenor_;
        }

        schedule_ =
            MakeSchedule().from(startDate)
                          .to(endDate)
                          .withFrequency(frequency_)
                          .withCalendar(calendar_)
                          .withConvention(paymentConvention_)
                          .withTerminationDateConvention(Unadjusted)
                          .withRule(rule_);
        earliestDate_ = schedule_.dates().front();
        latestDate_   = calendar_.adjust(schedule_.dates().back(),
                                         paymentConvention_);
        if (model_ == CreditDefaultSwap::ISDA)
            ++latestDate_;
    }

    SpreadCdsHelper::SpreadCdsHelper(
                              const std::variant<Rate, Handle<Quote>>& runningSpread,
                              const Period& tenor,
                              Integer settlementDays,
                              const Calendar& calendar,
                              Frequency frequency,
                              BusinessDayConvention paymentConvention,
                              DateGeneration::Rule rule,
                              const DayCounter& dayCounter,
                              Real recoveryRate,
                              const Handle<YieldTermStructure>& discountCurve,
                              bool settlesAccrual,
                              bool paysAtDefaultTime,
                              const Date& startDate,
                              const DayCounter& lastPeriodDayCounter,
                              const bool rebatesAccrual,
                              const CreditDefaultSwap::PricingModel model)
    : CdsHelper(runningSpread, tenor, settlementDays, calendar,
                frequency, paymentConvention, rule, dayCounter,
                recoveryRate, discountCurve, settlesAccrual, paysAtDefaultTime,
                startDate, lastPeriodDayCounter, rebatesAccrual, model) {}

    Real SpreadCdsHelper::impliedQuote() const {
        swap_->recalculate();
        return swap_->fairSpread();
    }

    void SpreadCdsHelper::resetEngine() {
        swap_ = ext::make_shared<CreditDefaultSwap>(
            Protection::Buyer, 100.0, 0.01, schedule_, paymentConvention_,
            dayCounter_, settlesAccrual_, paysAtDefaultTime_, protectionStart_,
            ext::shared_ptr<Claim>(), lastPeriodDC_, rebatesAccrual_, evaluationDate_);

        switch (model_) {
          case CreditDefaultSwap::ISDA:
            swap_->setPricingEngine(ext::make_shared<IsdaCdsEngine>(
                probability_, recoveryRate_, discountCurve_, false,
                IsdaCdsEngine::Taylor, IsdaCdsEngine::HalfDayBias,
                IsdaCdsEngine::Piecewise));
            break;
          case CreditDefaultSwap::Midpoint:
            swap_->setPricingEngine(ext::make_shared<MidPointCdsEngine>(
                probability_, recoveryRate_, discountCurve_));
            break;
          default:
            QL_FAIL("unknown CDS pricing model: " << model_);
        }
    }

    UpfrontCdsHelper::UpfrontCdsHelper(
                              const std::variant<Rate, Handle<Quote>>& upfront,
                              Rate runningSpread,
                              const Period& tenor,
                              Integer settlementDays,
                              const Calendar& calendar,
                              Frequency frequency,
                              BusinessDayConvention paymentConvention,
                              DateGeneration::Rule rule,
                              const DayCounter& dayCounter,
                              Real recoveryRate,
                              const Handle<YieldTermStructure>& discountCurve,
                              Natural upfrontSettlementDays,
                              bool settlesAccrual,
                              bool paysAtDefaultTime,
                              const Date& startDate,
                              const DayCounter& lastPeriodDayCounter,
                              const bool rebatesAccrual,
                              const CreditDefaultSwap::PricingModel model)
    : CdsHelper(upfront, tenor, settlementDays, calendar,
                frequency, paymentConvention, rule, dayCounter,
                recoveryRate, discountCurve, settlesAccrual, paysAtDefaultTime,
                startDate, lastPeriodDayCounter, rebatesAccrual, model),
      upfrontSettlementDays_(upfrontSettlementDays),
      upfrontDate_(upfrontDate()),
      runningSpread_(runningSpread) {}

    Date UpfrontCdsHelper::upfrontDate() {
        return calendar_.advance(evaluationDate_, upfrontSettlementDays_, Days, paymentConvention_);
    }

    void UpfrontCdsHelper::initializeDates() {
        CdsHelper::initializeDates();
        upfrontDate_ = upfrontDate();
    }

    void UpfrontCdsHelper::resetEngine() {
        swap_ = ext::make_shared<CreditDefaultSwap>(
            Protection::Buyer, 100.0, 0.01, runningSpread_, schedule_,
            paymentConvention_, dayCounter_, settlesAccrual_,
            paysAtDefaultTime_, protectionStart_, upfrontDate_,
            ext::shared_ptr<Claim>(), lastPeriodDC_, rebatesAccrual_,
            evaluationDate_);

        switch (model_) {
          case CreditDefaultSwap::ISDA:
            swap_->setPricingEngine(ext::make_shared<IsdaCdsEngine>(
                probability_, recoveryRate_, discountCurve_, false,
                IsdaCdsEngine::Taylor, IsdaCdsEngine::HalfDayBias,
                IsdaCdsEngine::Piecewise));
            break;
          case CreditDefaultSwap::Midpoint:
            swap_->setPricingEngine(ext::make_shared<MidPointCdsEngine>(
                probability_, recoveryRate_, discountCurve_));
            break;
          default:
            QL_FAIL("unknown CDS pricing model: " << model_);
        }
    }

    Real UpfrontCdsHelper::impliedQuote() const {
        SavedSettings backup;
        Settings::instance().includeTodaysCashFlows() = true;
        swap_->recalculate();
        return swap_->fairUpfront();
    }


    SpreadCdsIndexHelper::SpreadCdsIndexHelper(
        const Handle<Quote>& fairIndexSpread,
        const Rate runningSpread,
        const Period& tenor,
        Integer settlementDays,
        const Calendar& calendar,
        Frequency frequency,
        BusinessDayConvention paymentConvention,
        DateGeneration::Rule rule,
        const DayCounter& dayCounter,
        const Handle<YieldTermStructure>& discountCurve,
        const std::vector<Handle<DefaultProbabilityTermStructure>>& baseTermStructures,
        const std::vector<Real>& recoveryRates,
        const std::vector<Real>& weights,
        bool settlesAccrual,
        bool paysAtDefaultTime,
        const Date& startDate,
        const DayCounter& lastPeriodDayCounter,
        bool rebatesAccrual,
        CreditDefaultSwap::PricingModel model)
    : SpreadCdsHelper(fairIndexSpread,
                      tenor,
                      settlementDays,
                      calendar,
                      frequency,
                      paymentConvention,
                      rule,
                      dayCounter,
                      std::accumulate(recoveryRates.begin(), recoveryRates.end(), 0.0) / recoveryRates.size(),  // not used
                      discountCurve,
                      settlesAccrual,
                      paysAtDefaultTime,
                      startDate,
                      lastPeriodDayCounter,
                      rebatesAccrual),
      runningSpread_(runningSpread), baseTermStructures_(baseTermStructures),
      recoveryRates_(recoveryRates), adjustedTermStructures_(baseTermStructures_.size()),
      swaps_(baseTermStructures_.size()), weights_(weights) {
        QL_REQUIRE(baseTermStructures_.size() > 0, "baseTermStructures_.size() > 0 required.");
        QL_REQUIRE(baseTermStructures_.size() == recoveryRates_.size(),
                   "baseTermStructures_.size()== recoveryRates_.size() required.");
        QL_REQUIRE(baseTermStructures_.size() == weights_.size(),
                   "baseTermStructures_.size()== weights_.size() required.");
    }

    SpreadCdsIndexHelper::SpreadCdsIndexHelper(
        Rate fairIndexSpread,
        Rate runningSpread,
        const Period& tenor,
        Integer settlementDays, // ISDA: 1
        const Calendar& calendar,
        Frequency frequency,                     // ISDA: Quarterly
        BusinessDayConvention paymentConvention, // ISDA:Following
        DateGeneration::Rule rule,               // ISDA: CDS
        const DayCounter& dayCounter,            // ISDA: Actual/360
        const Handle<YieldTermStructure>& discountCurve,
        const std::vector<Handle<DefaultProbabilityTermStructure>>& baseTermStructures,
        const std::vector<Real>& recoveryRates,
        const std::vector<Real>& weights,
        bool settlesAccrual,
        bool paysAtDefaultTime,
        const Date& startDate,
        const DayCounter& lastPeriodDayCounter, // ISDA: Actual/360(inc)
        bool rebatesAccrual,                    // ISDA: true
        CreditDefaultSwap::PricingModel model)
    : SpreadCdsHelper(fairIndexSpread,
                      tenor,
                      settlementDays,
                      calendar,
                      frequency,
                      paymentConvention,
                      rule,
                      dayCounter,
                      std::accumulate(recoveryRates.begin(), recoveryRates.end(), 0.0) / recoveryRates.size(),  // not used
                      discountCurve,
                      settlesAccrual,
                      paysAtDefaultTime,
                      startDate,
                      lastPeriodDayCounter,
                      rebatesAccrual),
      runningSpread_(runningSpread), baseTermStructures_(baseTermStructures),
      recoveryRates_(recoveryRates), adjustedTermStructures_(baseTermStructures_.size()),
      swaps_(baseTermStructures_.size()), weights_(weights) {
        QL_REQUIRE(baseTermStructures_.size() > 0, "baseTermStructures_.size() > 0 required.");
        QL_REQUIRE(baseTermStructures_.size() == recoveryRates_.size(),
                   "baseTermStructures_.size()== recoveryRates_.size() required.");
        QL_REQUIRE(baseTermStructures_.size() == weights_.size(),
                   "baseTermStructures_.size()== weights_.size() required.");
    }

    void SpreadCdsIndexHelper::resetEngine() {
        // we use swap_ as a single representative constituent instrument for inspection purposes,
        // but it is not used for pricing
        swap_ = ext::make_shared<CreditDefaultSwap>(
            Protection::Buyer, 100.0, runningSpread_, schedule_, paymentConvention_,
            dayCounter_, settlesAccrual_, paysAtDefaultTime_, protectionStart_,
            ext::shared_ptr<Claim>(), lastPeriodDC_, rebatesAccrual_, evaluationDate_);

        for (Size k = 0; k < adjustedTermStructures_.size(); ++k) {
            // instrument with weight as notional and 100bp spread
            swaps_[k] = ext::shared_ptr<CreditDefaultSwap>(new CreditDefaultSwap(
                Protection::Buyer, 100.0, runningSpread_, schedule_, paymentConvention_,
                dayCounter_, settlesAccrual_, paysAtDefaultTime_, protectionStart_, ext::shared_ptr<Claim>(),
                lastPeriodDC_, rebatesAccrual_, evaluationDate_));

            // credit curve with adjustment
            adjustedTermStructures_[k].linkTo(ext::shared_ptr<DefaultProbabilityTermStructure>(
                new AdjustedSurvivalProbabilityStructure(baseTermStructures_[k], probability_)));

            // engine with adjusted curve
            switch (model_) {
                case CreditDefaultSwap::ISDA:
                    swaps_[k]->setPricingEngine(ext::make_shared<IsdaCdsEngine>(
                        adjustedTermStructures_[k], recoveryRates_[k], discountCurve_, false,
                        IsdaCdsEngine::Taylor, IsdaCdsEngine::HalfDayBias,
                        IsdaCdsEngine::Piecewise));
                    break;
                case CreditDefaultSwap::Midpoint:
                    swaps_[k]->setPricingEngine(ext::make_shared<MidPointCdsEngine>(
                        adjustedTermStructures_[k], recoveryRates_[k], discountCurve_));
                    break;
                default:
                    QL_FAIL("CDS pricing model must be IsdaCdsEngine or MidPointCdsEngine: " << model_);
            }
        }
    }

    Real SpreadCdsIndexHelper::impliedQuote() const {
        for (Size k = 0; k < swaps_.size(); ++k) {
            swaps_[k]->recalculate();
        }
        Real couponLegNPV = 0.0;
        Real accrualRebateNPV = 0.0;
        Real defaultLegNPV = 0.0;
        for (Size k = 0; k < swaps_.size(); ++k) {
            couponLegNPV += weights_[k] * swaps_[k]->couponLegNPV();
            accrualRebateNPV += weights_[k] * swaps_[k]->accrualRebateNPV();
            defaultLegNPV += weights_[k] * swaps_[k]->defaultLegNPV();
        }
        const Real riskyAnnity = -(couponLegNPV + accrualRebateNPV) / runningSpread_;
        // save for inspection
        couponLegNPV_ = couponLegNPV;
        accrualRebateNPV_ = accrualRebateNPV;
        defaultLegNPV_ = defaultLegNPV;
        riskyAnnuity_ = riskyAnnity;
        //
        QL_REQUIRE(riskyAnnity > 0.0, "riskyAnnity > 0.0 required");
        return defaultLegNPV / riskyAnnity;
    }

}
