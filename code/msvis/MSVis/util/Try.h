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
#include <stdexcept>
#include <memory>

namespace casa {

namespace vi {

class NoSuchElementException :
		public std::runtime_error {
public:
	NoSuchElementException()
		: std::runtime_error("no such element")
		{}
};

// TryBase exists to enable tests in Try template methods of whether
// Try::value_type is a Try type
struct TryBase {};

template <typename A>
class Try
	: protected TryBase {

public:

	typedef A value_type;

	// Constructors

	Try()
		: m_isSuccess(false)
		, m_exception(std::make_exception_ptr(NoSuchElementException())) {}

	Try(const A& a)
		: m_isSuccess(true) {
		m_value.reset(new A(a));
	}

	Try(A&& a)
		: m_isSuccess(true) {
		m_value.reset(new A(std::move(a)));
	}

	Try(const std::exception_ptr& e)
		: m_isSuccess(false)
		, m_exception(e) {}

	Try(std::exception_ptr&& e)
		: m_isSuccess(false)
		, m_exception(std::move(e)) {}

	Try(const Try<A>& ta) {
		*this = ta;
	}

	Try(Try<A>&& ta) {
		*this = ta;
	}

	// Static methods

	/* from()
	 *
	 * Type F should be callable () -> A
	 */
	template <typename F,
	          class = std::enable_if<
		          std::is_convertible<std::result_of<F()>,A>::value,A> >
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
	 * Type F should be callable const A& -> B (for some B)
	 */
	template <typename F,
	          typename B = typename std::result_of<F(const A&)>::type>
#if __cplusplus >= 201402L
	static auto
#else
	static std::function<Try<B>(const Try<A>&)>
#endif
	lift(F f) {
		return [f](const Try<A>& ta) {
			return ta.map(f);
		};
	}

	// Instance methods

	/* copy assignment
	 */
	Try<A>& operator=(const Try<A>& ta) {
		m_isSuccess = ta.m_isSuccess;
		if (ta.m_isSuccess) m_value.reset(new A(*ta.m_value));
		else m_exception = ta.m_exception;
		return *this;
	}

	/* move assignment
	 */
	Try<A>& operator=(Try<A>&& ta) {
		m_isSuccess = ta.m_isSuccess;
		m_value = ta.m_value;
		m_exception = std::move(ta.m_exception);
		return *this;
	}

	/* equality operator
	 */
	bool operator==(const Try<A> &ta) const {
		return (
			m_isSuccess == ta.m_isSuccess
			&& ((m_isSuccess && *m_value == *ta.m_value)
			    || (!m_isSuccess && m_exception == ta.m_exception)));
	}

	bool operator!=(const Try<A> &ta) const {
		return !(*this == ta);
	}

	/* andThen()
	 *
	 * Type F should be callable () -> Try<B> (for some B)
	 */
	template <typename F,
	          typename TB = typename std::result_of<F()>::type,
	          class = std::enable_if<std::is_base_of<TryBase,TB>::value> >
	TB
	andThen(F&& f) const {
		return flatMap([&f](const A&){ return f(); });
	}

	/* filter()
	 *
	 * Type F should be callable const A& -> bool
	 */
	template <typename F,
	          class = std::enable_if<
		          std::is_convertible<std::result_of<F(const A&)>,bool>::value> >
	Try<A>
	filter(F&& f) const {
		if (m_isSuccess) {
			try {
				if (std::forward<F>(f)(*m_value)) return *this;
				return Try<A>();
			} catch (std::exception &) {
				return Try<A>(std::current_exception());
			}
		} else {
			return *this;
		}
	}

	/* flatMap()
	 *
	 * Type F should be callable const A& -> Try<B> (for some B)
	 */
	template <typename F,
	          typename TB = typename std::result_of<F(const A&)>::type,
	          class = std::enable_if<std::is_base_of<TryBase,TB>::value> >
	TB
	flatMap(F&& f) const {
		return map(std::forward<F>(f)).flatten();
	};

	/* flatten()
	 *
	 * May only flatten instances of Try<Try<...> >, removes one level of Try
	 */
	template <class = std::enable_if<std::is_base_of<TryBase,A>::value> >
	A
	flatten() const {
		if (m_isSuccess) return *m_value;
		else return Try<typename value_type::value_type>(m_exception);
	}

	/* fold()
	 *
	 * Type Err should be callable const std::exception_ptr & -> B (for some B)
	 *
	 * Type Val should be callable const A& -> B
	 */
	template <typename Err,
	          typename Val,
	          typename B = typename std::result_of<Err(const std::exception_ptr&)>::type,
	          typename C = typename std::result_of<Val(const A&)>::type,
	          class = std::enable_if<std::is_same<B,C>::value> >
	B
	fold(Err&& err, Val&& val) const {
		if (m_isSuccess) return std::forward<Val>(val)(*m_value);
		else return std::forward<Err>(err)(m_exception);
	};

	/* foreach()
	 *
	 * Type F should be callable const A& -> anything
	 */
	template <typename F>
	void
	foreach(F&& f) const {
		if (m_isSuccess) std::forward<F>(f)(*m_value);
	}

	/* get()
	 *
	 * Get value from Try instance, throws exception if failure
	 */
	A
	get() const {
		if (m_isSuccess) return *m_value;
		else std::rethrow_exception(m_exception);
	}

	/* getOrElse()
	 *
	 * Type F should be callable () -> B (for some B)
	 */
	template <typename F,
	          typename B = typename std::result_of<F()>::type,
	          class = typename std::enable_if<std::is_base_of<B, A>::value> >
	B
	getOrElse(F&& f) const {
		if (m_isSuccess) return *m_value;
		else return std::forward<F>(f)();
	};

	/* getOrElse_()
	 *
	 * Non-lazy version of getOrElse
	 */
	template <typename B,
	          class = typename std::enable_if<std::is_base_of<B, A>::value> >
	B
	getOrElse_(const B& b) const {
		return getOrElse([&b](){ return b; });
	};

	bool
	isSuccess() const {
		return m_isSuccess;
	}

	/* map()
	 *
	 * Type F should be callable const A& -> B (for some B)
	 */
	template <typename F,
	          typename B = typename std::result_of<F(const A&)>::type>
	Try<B>
	map(F&& f) const {
		if (m_isSuccess) {
			try {
				return Try<B>(std::forward<F>(f)(*m_value));
			} catch (std::exception &) {
				return Try<B>(std::current_exception());
			}
		} else {
			return Try<B>(m_exception);
		}
	};

	/* orElse
	 *
	 * Type F should be callable () -> Try<B> (for some B)
	 */
	template <typename F,
	          typename TB = typename std::result_of<F()>::type,
	          class = std::enable_if<std::is_base_of<TryBase,TB>::value>,
	          class = std::enable_if<
		          std::is_base_of<typename TB::value_type,A>::value> >
	TB
	orElse(F&& f) const {
		if (m_isSuccess)
			return Try<typename TB::value_type>(*m_value);
		else
			return std::forward<F>(f)();
	};

	/* orElse_
	 *
	 * Non-lazy version of orElse
	 */
	template <typename TB,
	          class = std::enable_if<std::is_base_of<TryBase,TB>::value>,
	          class = std::enable_if<
		          std::is_base_of<typename TB::value_type,A>::value> >
	TB
	orElse_(const TB& tb) const {
		return orElse([&tb](){ return tb; });
	};

	/* recoverWith()
	 *
	 * Type F should be callable const std::exception_ptr& -> A
	 */
	template <typename F,
	          class = std::enable_if<
		          std::is_convertible<
			          std::result_of<F(const std::exception&)>,A>::value> >
	Try<A>
	recoverWith(F&& f) const {
		if (m_isSuccess) {
			return *this;
		} else {
			auto e = m_exception;
			return Try<A>::from([f, e](){ return f(e); });
		}
	}

	/* transform()
	 *
	 * Type Err should be callable const std::exception_ptr & -> Try<B> (for
	 * some B)
	 *
	 * Type Val should be callable const A& -> Try<B>
	 */
	template <typename Err, typename Val,
	          typename TB = typename std::result_of<Err(const std::exception_ptr&)>::type,
	          typename TC = typename std::result_of<Val(const A&)>::type,
	          class = std::enable_if<std::is_base_of<TryBase,TB>::value>,
	          class = std::enable_if<std::is_base_of<TryBase,TC>::value>,
	          class = std::enable_if<std::is_same<TB,TC>::value > >
	TB
	transform(Err&& err, Val&& val) const {
		Try<TB> folded = fold(
			std::forward<Err>(err),
			std::forward<Val>(val));
		return folded.flatten();
	};

	/* operator|()
	 *
	 * Type F should be callable const A& -> B (for some B)
	 */
	template <typename F,
	          typename B = typename std::result_of<F(const A&)>::type>
	Try<B>
	operator|(F&&f) const {
		return map(std::forward<F>(f));
	}

	/* operator>>=()
	 *
	 * Type F should be callable const A& -> Try<B> (for some B)
	 */
	template <typename F,
	          typename TB = typename std::result_of<F(const A&)>::type,
	          class = std::enable_if<std::is_base_of<TryBase,TB>::value> >
	TB
	operator>>=(F&& f) const {
		return flatMap(std::forward<F>(f));
	}

	/* operator>>()
	 *
	 * Type F should be callable () -> Try<B> (for some B)
	 */
	template <typename F,
	          typename TB = typename std::result_of<F()>::type,
	          class = std::enable_if<std::is_base_of<TryBase,TB>::value> >
	TB
	operator>>(F&& f) const {
		return andThen(std::forward<F>(f));
	}

private:
	bool m_isSuccess;

	std::unique_ptr<A> m_value;

	std::exception_ptr m_exception;
};

template <typename F,
          typename A = typename std::result_of<F()>::type>
Try<A>
try_(F&& f) {
	return Try<A>::from(f);
}

} // end namespace vi

} // end namespace casa

#endif /* TRY_H_ */
