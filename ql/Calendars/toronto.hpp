
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it under the
 terms of the QuantLib license.  You should have received a copy of the
 license along with this program; if not, please email quantlib-dev@lists.sf.net
 The license is also available online at http://quantlib.org/html/license.html

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file toronto.hpp
    \brief Toronto calendar
*/

#ifndef quantlib_toronto_calendar_h
#define quantlib_toronto_calendar_h

#include <ql/calendar.hpp>

namespace QuantLib {

    //! %Toronto calendar
    /*! Holidays:
      <ul>
      <li>Saturdays</li>
      <li>Sundays</li>
      <li>New Year's Day, January 1st (possibly moved to Monday)</li>
      <li>Good Friday</li>
      <li>Easter Monday</li>
      <li>Victoria Day, The Monday on or preceding 24 May</li>
      <li>Canada Day, July 1st (possibly moved to Monday)</li>
      <li>Provincial Holiday, first Monday of August</li>
      <li>Labour Day, first Monday of September</li>
      <li>Thanksgiving Day, second Monday of October</li>
      <li>Remembrance Day, November 11th</li>
      <li>Christmas, December 25th (possibly moved to Monday or Tuesday)</li>
      <li>Boxing Day, December 26th (possibly moved to Monday or Tuesday)</li>
      </ul>
    */
    class Toronto : public Calendar {
      private:
        class Impl : public Calendar::WesternImpl {
          public:
            std::string name() const { return "Toronto"; }
            bool isBusinessDay(const Date&) const;
        };
      public:
        Toronto()
        : Calendar(Handle<Calendar::Impl>(new Toronto::Impl)) {}
    };

}


#endif
