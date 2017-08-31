/* -*- mode: c++ -*- */
//# Try.h
//# Copyright (C) 2017
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU Library General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
//# License for more details.
//#
//# You should have received a copy of the GNU Library General Public License
//# along with this library; if not, write to the Free Software Foundation,
//# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
#ifndef TRY_H_
#define TRY_H_

#include <exception>
#include <functional>
#include <iostream>

namespace casa {

namespace vi {

class NoSuchElementException :
		public std::runtime_error {
public:
	NoSuchElementException()
		: std::runtime_error("no such element")
		{}
};

template <typename A>
class Try {

public:

	// Constructors

	Try()
		: is_success(true)
		, value() {}

	Try(const A& a)
		: is_success(true)
		, value(a) {}

	Try(A&& a)
		: is_success(true)
		, value(std::move(a)) {}

	Try(const std::exception_ptr& e)
		: is_success(false)
		, exception(e) {}

	Try(std::exception_ptr&& e)
		: is_success(false)
		, exception(std::move(e)) {}

	// Static methods

	/* from()
	 *
	 * Type F should be callable () -> A
	 */
	template <typename F>
	static Try<A>
	from(F&& f) {
		try {
			return Try<A>(std::forward<F>(f)());
		} catch (const std::exception &) {
			return Try<A>(std::current_exception());
		}
	}

	/* lift
	 *
	 * Type F should be callable const A& -> B
	 */
	template <typename B, typename F>
#if __cplusplus >= 201402L
	auto
#else
	static std::function<Try<B>(const Try<A>&)>
#endif
	lift(F f) {
		return [f](const Try<A>& ta) {
			return ta.map<B>(f);
		};
	}

	// Instance methods

	/* equality operator
	 *
	 * Value equality for successes, failures are never equal
	 */
	bool operator==(const Try<A> &ta) const {
		return (
			is_success == ta.is_success
			&& ((is_success && value == ta.value)
			    || (!is_success && exception == ta.exception)));
	}

	bool operator!=(const Try<A> &ta) const {
		return !(*this == ta);
	}

	/* filter()
	 *
	 * Type F should be callable const A& -> bool
	 */
	template <typename F>
	Try<A>
	filter(F&& f) const {
		if (is_success) {
			try {
				if (std::forward<F>(f)(value)) return *this;
				return Try<A>(
					std::make_exception_ptr(NoSuchElementException()));

			} catch (std::exception &) {
				return Try<A>(std::current_exception());
			}
		} else {
			return *this;
		}
	}

	/* flatMap()
	 *
	 * Type F should be callable const A& -> Try<B>
	 */
	template <typename B, typename F>
	Try<B>
	flatMap(F&& f) const {
		return map<Try<B> >(std::forward<F>(f)).flatten<B>();
	};

	/* flatten()
	 *
	 * May only flatten instances of Try<Try<...> >, removes one level of Try
	 */
	template <typename B,
	          class = typename std::enable_if<std::is_base_of<Try<B>, A>::value, B> >
	Try<B>
	flatten() const {
		if (is_success) return value;
		else return Try<B>(exception);
	}

	/* fold()
	 *
	 * Type Err should be callable const std::exception_ptr & -> B
	 * Type Val should be callable const A& -> B
	 */
	template <typename B, typename Err, typename Val>
	B
	fold(Err&& err, Val&& val) const {
		if (is_success) return std::forward<Val>(val)(value);
		else return std::forward<Err>(err)(exception);
	};

	/* foreach()
	 *
	 * Type F should be callable const A& -> anything
	 */
	template <typename F>
	void
	foreach(F&& f) const {
		if (is_success) std::forward<F>(f)(value);
	}

	/* get()
	 *
	 * Get value from Try instance, throws exception if failure
	 */
	A
	get() const {
		if (is_success) return value;
		else std::rethrow_exception(exception);
	}

	/* getOrElse()
	 *
	 * Type F should be callable () -> B
	 */
	template <typename B,
	          typename F,
	          class = typename std::enable_if<std::is_base_of<B, A>::value, B> >
	B
	getOrElse(F&& f) const {
		if (is_success) return value;
		else return std::forward<F>(f)();
	};

	bool
	isFailure() const {
		return !is_success;
	}

	bool
	isSuccess() const {
		return is_success;
	}

	/* map()
	 *
	 * Type F should be callable const A& -> B
	 */
	template <typename B, typename F>
	Try<B>
	map(F&& f) const {
		if (is_success) {
			try {
				return Try<B>(std::forward<F>(f)(value));
			} catch (std::exception &) {
				return Try<B>(std::current_exception());
			}
		} else {
			return Try<B>(exception);
		}
	};

	/* orElse
	 *
	 * Type F should be callable () -> Try<B>
	 */
	template <typename B,
	          typename F,
	          class = typename std::enable_if<std::is_base_of<B, A>::value, B> >
	Try<B>
	orElse(F&& f) const {
		if (is_success) return Try<B>(value);
		else return Try<Try<B> >::from(std::forward<F>(f)).flatten<B>();
	};

	/* transform()
	 *
	 * Type Err should be callable const std::exception_ptr & -> Try<B>
	 * Type Val should be callable const A& -> Try<B>
	 */
	template <typename B, typename Err, typename Val>
	Try<B>
	transform(Err&& err, Val&& val) const {
		return fold<Try<B> >(
			std::forward<Err>(err),
			std::forward<Val>(val)).
			flatten<B>();
	};

private:
	bool is_success;
	A value;
	std::exception_ptr exception;
};

} // end namespace vi

} // end namespace casa

#endif /* TRY_H_ */
