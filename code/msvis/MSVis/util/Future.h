/* -*- mode: c++ -*- */
//# Future.h
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
#ifndef FUTURE_H_
#define FUTURE_H_

#include <future>
#include <msvis/MSVis/util/Try.h>

namespace casa {

namespace vi {

// FutureBase exists to enable tests in Future template methods of whether
// Future::value_type is a Future type
struct FutureBase {};

template <typename A>
class Future {
public:

	typedef A value_type;

	// Constructors

	Future() {
		auto p = std::promise<Try<A> >();
		p.set_value(Try<A>());
		m_fta = p.get_future();
	}

	Future(const A& a) {
		auto p = std::promise<Try<A> >();
		p.set_value(Try<A>(a));
		m_fta = p.get_future();
	}

	Future(A&& a) {
		auto p = std::promise<Try<A> >();
		p.set_value(Try<A>(std::move(a)));
		m_fta = p.get_future();
	}

	Future(const std::exception_ptr& e) {
		auto p = std::promise<Try<A> >();
		p.set_value(Try<A>(e));
		m_fta = p.get_future();
	}

	Future(std::exception_ptr&& e) {
		auto p = std::promise<Try<A> >();
		p.set_value(Try<A>(std::move(e)));
		m_fta = p.get_future();
	}

	Future(const std::shared_future<Try<A> >& f)
		: m_fta(f) {}

	Future(std::shared_future<Try<A> >&& f)
		: m_fta(std::move(f)) {}

	// NB: the following method would be erroneous, since calling
	// std::future::share() is not const on std::future
	//
	// Future(const std::future<Try<A> >& f)
	//  : m_fta(f.share()) {};

	Future(std::future<Try<A> >&& f)
		: m_fta(f.share()) {}

	// Static methods

	/* from()
	 *
	 * Type F should be callable () -> A
	 */
	template <typename F,
	          class = std::enable_if<
		          std::is_convertible<std::result_of<F()>,A>::value,A> >
	static Future<A>
	from(F f) {
		return Future<A>(
			std::async(std::launch::async, [f](){ return try_(f); }));
	}

	/* lift()
	 *
	 * Type F should be callable const A& -> B (for some B)
	 */
	template <typename F,
	          typename B = typename std::result_of<F(const A&)>::type>
#if __cplusplus >= 201402L
	static auto
#else
	static std::function<Future<B>(const Future<A>&)>
#endif
	lift(F f) {
		return [f](const Future<A>& fa) {
			return fa.map(f);
		};
	}

	// Instance methods

	/* flatMap()
	 *
	 * Type F should be callable const &A -> Future<B> (for some B)
	 */
	template <typename F,
	          typename FB = typename std::result_of<F(const A&)>::type,
	          class = std::enable_if<std::is_base_of<FutureBase,FB>::value> >
	FB
	flatMap(F&& fn) const {
		return map(std::forward<F>(fn)).flatten();
	}

	/* flatten()
	 *
	 * May only flatten instances of Future<Future<...> >, removes one level of
	 * Future
	 */
	template <class = std::enable_if<std::is_base_of<FutureBase,A>::value> >
	A
	flatten() const {
		auto fta = m_fta;
		return Future<typename A::value_type>(
			std::async(
				std::launch::async,
				[fta](){
					return fta.get().recoverWith(
						[](const std::exception_ptr& e) {
							return Future<typename A::value_type>(e);
						}).get().get();
				}));
	}

	/* followedBy()
	 *
	 * Type F should be callable () -> Future<B> (for some B)
	 */
	template <typename F,
	          typename FB = typename std::result_of<F()>::type,
	          class = std::enable_if<std::is_base_of<FutureBase,FB>::value> >
	FB
	followedBy(F fn) const {
		return flatMap([fn](const A&){ return fn(); });
	}

	/* forEffect()
	 *
	 * Type F should be callable () -> Future<anything>
	 */
	template <typename F,
	          typename FB = typename std::result_of<F()>::type,
	          class = std::enable_if<std::is_base_of<FutureBase,FB>::value > >
	Future<A>
	forEffect(F fn) const {
		return flatMap(
			[fn](const A& a){
				return fn().map([a](const typename FB::value_type &){
						return a;
					});
			});
	}

	/* get()
	 */
	Try<A>
	get() const {
		return m_fta.get();
	}

	/* iterateUntil()
	 *
	 * Type P should be callable const A& -> bool
	 *
	 * Type F should be callable const A& -> Future<A>
	 */
	template <typename P,
	          typename F,
	          class = std::enable_if<
		          std::is_convertible<std::result_of<P(const A&)>,bool>::value>,
	          class = std::enable_if<
		          std::is_same<Future<A>,std::result_of<F(const A&)> >::value> >
	Future<A>
	iterateUntil(P p, F&& f) const {
		return iterateWhile(
			[p](const A& a){ return !p(a); },
			std::forward<F>(f));
	}

	/* iterateWhile()
	 *
	 * Type P should be callable const A& -> bool
	 *
	 * Type F should be callable const A& -> Future<A>
	 */
	template <typename P,
	          typename F,
	          class = std::enable_if<
		          std::is_convertible<std::result_of<P(const A&)>,bool>::value>,
	          class = std::enable_if<
		          std::is_same<Future<A>,std::result_of<F(const A&)> >::value> >
	Future<A>
	iterateWhile(P p, F f) const {
		Future<A> fa = *this;
		return
			Future<A>(
				std::async(
					std::launch::async,
					[fa, p, f]() {
						auto result = fa;
						while (result.get().map(p).getOrElse_(false))
							result = result.flatMap(f);
						return result.get();
					}));
	}

	/* map()
	 *
	 * Type F should be callable const A& -> B (for some B)
	 */
	template <typename F,
	          typename B = typename std::result_of<F(const A&)>::type>
	Future<B>
	map(F fn) const {
		auto fta = m_fta;
		return Future<B>(
			std::async(
				std::launch::async,
				[fta, fn](){ return fta.get().map(fn); }));
	}

	/* recoverWith()
	 *
	 * Type F should be callable const std::exception_ptr & -> A
	 */
	template <typename F,
	          class = std::enable_if<
		          std::is_convertible<
			          std::result_of<F(std::exception_ptr &)>,A>::value> >
	Future<A>
	recoverWith(F fn) const {
		auto fta = m_fta;
		return Future<A>(
			std::async(
				std::launch::async,
				[fta, fn]() {
					return fta.get().recoverWith(fn);
				}));
	}

	/* wait()
	 */
	void
	wait() const {
		m_fta.wait();
	}

	/* waitFor()
	 */
	template <class Rep, class Period>
	std::future_status
	waitFor(const std::chrono::duration<Rep,Period>& rel_time) const {
		return m_fta.wait_for(rel_time);
	}

	/* waitUntil()
	 */
	template <class Clock, class Duration>
	std::future_status
	waitUntil(const std::chrono::time_point<Clock,Duration>& abs_time) const {
		return m_fta.wait_until(abs_time);
	}

	/* operator>>=()
	 *
	 * Type F should be callable const A& -> Future<B> (for some B)
	 */
	template <typename F,
	          typename FB = typename std::result_of<F(const A&)>::type,
	          class = std::enable_if<std::is_base_of<FutureBase,FB>::value> >
	FB
	operator>>=(F&& fn) const {
		return flatMap(std::forward<F>(fn));
	}

	/* operator>>()
	 *
	 * Type F should be callable () -> Future<B> (for some B)
	 */
	template <typename F,
	          typename FB = typename std::result_of<F()>::type,
	          class = std::enable_if<std::is_base_of<FutureBase,FB>::value> >
	FB
	operator>>(F&& fn) const {
		return followedBy(std::forward<F>(fn));
	}


	/* operator|()
	 *
	 * Type F should be callable const A& -> B (for some B)
	 */
	template <typename F,
	          typename B = typename std::result_of<F(const A&)>::type>
	Future<B>
	operator|(F&& fn) const {
		return map(std::forward<F>(fn));
	}

private:

	std::shared_future<Try<A> > m_fta;

};

template <typename F,
          typename A = typename std::result_of<F()>::type>
Future<A>
future(F&& f) {
	return Future<A>::from(std::forward<F>(f));
}


} // end namespace vi

} // end namespace casa


#endif /* FUTURE_H_ */
