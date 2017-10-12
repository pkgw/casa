//# Future_Test.h
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
#include <msvis/MSVis/util/Future.h>
#include <gtest/gtest.h>
#include <chrono>

using namespace casa::vi;

int
main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

class FutureTest
	: public ::testing::Test {
public:
	struct BadIntException : public std::runtime_error {
		BadIntException()
			: std::runtime_error("BadIntException")
			{}
	};

	struct BadOuterIntException : public std::runtime_error {
		BadOuterIntException()
			: std::runtime_error("BadOuterIntException")
			{}
	};

	static const std::chrono::milliseconds initWait;

	static int
	plus2(const int& i) {
		return i + 2;
	}

	static bool
	isEven(const int& i) {
		return i % 2 == 0;
	}
};

const std::chrono::milliseconds FutureTest::initWait =
	std::chrono::milliseconds(500);

TEST_F(FutureTest, Constructors) {
	// construct with a value
	int four = 4;
	Future<int> a(four);
	// the future should be completed immediately, the following is an
	// approximation of "immediately"
	std::chrono::nanoseconds zero(0);
	EXPECT_EQ(a.waitFor(zero), std::future_status::ready);
	EXPECT_EQ(a.get().get(), four);

	// construct with an exception
	Future<int> c = Future<int>(std::make_exception_ptr(BadIntException()));
	EXPECT_THROW(c.get().get(), BadIntException);

	// construct with a deferred value
	std::promise<int> p;
	Future<int> d = future([&p](){ return p.get_future().get(); });
	EXPECT_EQ(d.waitFor(initWait), std::future_status::timeout);
	p.set_value(four);
	EXPECT_EQ(d.get().get(), four);

	// construct with a deferred exception
	std::promise<int> q;
	Future<int> e = future([&q](){ return q.get_future().get(); });
	q.set_exception(std::make_exception_ptr(BadIntException()));
	EXPECT_THROW(e.get().get(), BadIntException);
}

TEST_F(FutureTest, Map) {
	// map a value through two non-failing transformations
	std::promise<Try<int> > p1;
	Future<int> f1(p1.get_future());
	Future<int> f1m = f1.map(plus2);
	Future<bool> f1mm = f1m.map(isEven);
	// nothing should be completed yet
	EXPECT_EQ(f1m.waitFor(initWait), std::future_status::timeout);
	EXPECT_EQ(f1mm.waitFor(initWait), std::future_status::timeout);
	int init = 2;
	Try<int> initVal(init);
	p1.set_value(initVal);
	// get values in "reverse" order to additionally check multiple accesses to
	// a single result
	EXPECT_EQ(f1mm.get().get(), (init + 2) % 2 == 0);
	EXPECT_EQ(f1m.get().get(), init + 2);
	EXPECT_EQ(f1.get(), init);

	// map a value through one non-failing, one failing, and one non-failing
	// transformation
	std::promise<Try<int> > p2;
	Future<int> f2(p2.get_future());
	Future<int> f2m = f2.map([](const int& i){ return i + 2; });
	Future<int> f2mm =
		f2m.map([](const int& ){ throw BadIntException(); return 0; });
	Future<bool> f2mmm = f2mm.map(isEven);
	p2.set_value(initVal);
	EXPECT_THROW(f2mmm.get().get(), BadIntException);
	EXPECT_THROW(f2mm.get().get(), BadIntException);
	EXPECT_EQ(f2m.get().get(), init + 2);
	EXPECT_EQ(f2.get().get(), init);

	// the same as above, through some temporary values
	std::promise<Try<int> > p3;
	Future<bool> f3 =
		Future<int>(p3.get_future())
		.map(plus2)
		.map([](const int& ){ throw BadIntException(); return 0; })
		.map(isEven);
	p3.set_value(initVal);
	EXPECT_THROW(f3.get().get(), BadIntException);

	// check that a mapped function/callable object is being copied into the
	// asynchronous scope
	std::promise<Try<int> > p4;
	Future<int> f4(p4.get_future());
	Future<int> f4m;
	{
		auto i2 = [](const int& i){ return i + 2; };
		f4m = f4.map(i2);
	}
	p4.set_value(initVal);
	EXPECT_EQ(f4m.get().get(), init + 2);
}

TEST_F(FutureTest, MapOperator) {
	std::promise<Try<int> > p1;
	Future<int> f1(p1.get_future());
	Future<int> f1m = f1 | plus2;
	Future<bool> f1mm = f1m | isEven;
	int initVal = 11;
	p1.set_value(Try<int>(initVal));
	EXPECT_EQ(f1m.get().get(), plus2(initVal));
	EXPECT_EQ(f1mm.get().get(), isEven(plus2(initVal)));

	/* same as above, but with mapping chained into a single statement */
	std::promise<Try<int> > p2;
	Future<int> f2(p2.get_future());
	Future<bool> f2mm = (f2 | plus2) | isEven;
	p2.set_value(Try<int>(initVal));
	EXPECT_EQ(f2mm.get().get(), isEven(plus2(initVal)));
}

TEST_F(FutureTest, Flatten) {
	// flatten success of success
	std::promise<Try<Future<int> > > p1Outer;
	Future<Future<int> > f1Outer(p1Outer.get_future());
	Future<int> f1 = f1Outer.flatten();
	std::promise<Try<int> > p1Inner;
	Future<int> f1Inner(p1Inner.get_future());
	p1Outer.set_value(f1Inner);
	EXPECT_EQ(f1.waitFor(initWait), std::future_status::timeout);
	EXPECT_EQ(f1Inner.waitFor(initWait), std::future_status::timeout);
	int init = 13;
	Try<int> initVal(init);
	p1Inner.set_value(initVal);
	EXPECT_EQ(f1.get(), initVal);

	// flatten success of failure
	std::promise<Try<Future<int> > > p2Outer;
	Future<Future<int> > f2Outer(p2Outer.get_future());
	std::promise<Try<int> > p2Inner;
	Future<int> f2Inner(p2Inner.get_future());
	p2Outer.set_value(f2Inner);
	Try<int> initFail(std::make_exception_ptr(BadIntException()));
	p2Inner.set_value(initFail);
	EXPECT_THROW(f2Outer.flatten().get().get(), BadIntException);

	// flatten failure of success
	std::promise<Try<Future<int> > > p3Outer;
	Future<Future<int> > f3Outer(p3Outer.get_future());
	std::promise<Try<int> > p3Inner;
	Future<int> f3Inner(p3Inner.get_future());
	p3Outer.set_value(
		f3Inner | [](const int &){ throw BadOuterIntException(); return 0; });
	p3Inner.set_value(initVal);
	EXPECT_THROW(f3Outer.flatten().get().get(), BadOuterIntException);

	// flatten failure of failure
	std::promise<Try<Future<int> > > p4Outer;
	Future<Future<int> > f4Outer(p4Outer.get_future());
	std::promise<Try<int> > p4Inner;
	Future<int> f4Inner(p4Inner.get_future());
	p4Outer.set_value(
		Try<Future<int> >(std::make_exception_ptr(BadOuterIntException())));
	p4Inner.set_value(initFail);
	EXPECT_THROW(f4Outer.flatten().get().get(), BadOuterIntException);
}

TEST_F(FutureTest, FlatMap) {
	// simple flatMap, w.o. failure
	std::promise<Try<int> > p;
	Future<int> a(p.get_future());
	std::promise<Try<int> > q;
	Future<int> b(q.get_future());
	Future<int> c =
		a.flatMap(
			[&b](const int &i) {
				return b | [&i](const int &j){ return i + j; };
			});
	Try<int> jVal(2);
	q.set_value(jVal);
	EXPECT_EQ(c.waitFor(initWait), std::future_status::timeout);
	Try<int> iVal(1);
	p.set_value(iVal);
	EXPECT_EQ(c.get().get(), iVal.get() + jVal.get());

	// flatMap to failure
	Future<int> d =
		a.flatMap(
			[](const int&){
				return Future<int>(std::make_exception_ptr(BadIntException()));
			});
	EXPECT_THROW(d.get().get(), BadIntException);

	// flatMap from failure
	std::promise<Try<int> > r;
	Future<int> e(r.get_future());
	Future<int> f =
		e.flatMap(
			[](const int &i){
				return Future<int>(i + 2);
			});
	r.set_value(Try<int>(std::make_exception_ptr(BadOuterIntException())));
	EXPECT_THROW(f.get().get(), BadOuterIntException);

	// flatMap from failure to failure
	Future<int> g =
		e.flatMap(
			[](const int &){
				return Future<int>(std::make_exception_ptr(BadIntException()));
			});
	EXPECT_THROW(g.get().get(), BadOuterIntException);
}

TEST_F(FutureTest, FlatMapOperator) {
	std::promise<Try<int> > p1;
	Future<int> f1(p1.get_future());
	Future<int> f1m = f1 >>= [](const int &i){ return Future<int>(plus2(i)); };
	Future<bool> f1mm =
		f1m >>= [](const int &i){ return Future<bool>(isEven(i)); };
	int initVal = 14;
	p1.set_value(Try<int>(initVal));
	EXPECT_EQ(f1m.get().get(), plus2(initVal));
	EXPECT_EQ(f1mm.get().get(), isEven(plus2(initVal)));

	/* same as above, but with mapping chained into a single statement */
	std::promise<Try<int> > p2;
	Future<int> f2(p2.get_future());
	Future<bool> f2mm
		((f2
		  >>= [](const int &i){ return Future<int>(plus2(i)); })
		 >>= [](const int &i){ return Future<bool>(isEven(i)); });
	p2.set_value(Try<int>(initVal));
	EXPECT_EQ(f2mm.get().get(), isEven(plus2(initVal)));
}

TEST_F(FutureTest, FollowedBy) {
	std::promise<Try<int> > p1;
	Future<int> f1(p1.get_future());
	int init = 15;
	Future<int> g1 = f1.followedBy([&init](){ return Future<int>(init + 1); });
	EXPECT_EQ(g1.waitFor(initWait), std::future_status::timeout);
	p1.set_value(Try<int>(init));
	EXPECT_EQ(g1.get().get(), init + 1);
}

TEST_F(FutureTest, FollowedByOperator) {
	std::promise<Try<int> > p1;
	Future<int> f1(p1.get_future());
	int init = 15;
	Future<int> g1 = f1 >> [&init](){ return Future<int>(init + 1); };
	p1.set_value(Try<int>(init));
	EXPECT_EQ(g1.get().get(), init + 1);
}

TEST_F(FutureTest, ForEffect) {
	std::promise<Try<int> > p1;
	Future<int> f1(p1.get_future());
	int init = 15;
	int init1 = init;
	Future<int> g1 = f1.forEffect(
		[&init, &init1](){ init1 = init + 1;  return Future<int>(init1); });
	EXPECT_EQ(g1.waitFor(initWait), std::future_status::timeout);
	p1.set_value(Try<int>(init));
	EXPECT_EQ(g1.get().get(), init);
	EXPECT_EQ(init1, init + 1);

	// with failing effect
	std::promise<Try<int> > p2;
	Future<int> f2(p2.get_future());
	Future<int> g2 = f2.forEffect(
		[](){ throw BadIntException(); return Future<int>(0); });
	p2.set_value(Try<int>(init));
	EXPECT_THROW(g2.get().get(), BadIntException);
}

TEST_F(FutureTest, IterateUntil) {
	std::promise<Try<int> > p1;
	Future<int> f1(p1.get_future());
	Future<int> g1 =
		f1.iterateUntil(
			[](const int& i) { return i < 0; },
			[](const int& i) { return Future<int>(i - 1); });
	p1.set_value(Try<int>(12));
	EXPECT_EQ(g1.get().get(), -1);
}

TEST_F(FutureTest, IterateWhile) {
	std::promise<Try<int> > p1;
	Future<int> f1(p1.get_future());
	Future<int> g1 =
		f1.iterateWhile(
			[](const int& i) { return i > 0; },
			[](const int& i) { return Future<int>(i - 1); });
	p1.set_value(Try<int>(12));
	EXPECT_EQ(g1.get().get(), 0);

	// with an error along the way
	std::promise<Try<int> > p2;
	int init = 12;
	Future<int> f2(p2.get_future());
	Future<int> g2 =
		f2.iterateWhile(
			[](const int& i) { return i > 0; },
			[&init](const int& i) {
				if (i > init / 2)
					return Future<int>(i - 1);
				else
					return Future<int>(
						std::make_exception_ptr(BadIntException()));
			});
	p2.set_value(Try<int>(init));
	EXPECT_THROW(g2.get().get(), BadIntException);
}

TEST_F(FutureTest, Lift) {
	auto liftedPlus2 = Future<int>::lift(plus2);
	std::promise<Try<int> > p1;
	Future<int> f1(p1.get_future());
	Future<int> g1 = liftedPlus2(f1);
	int initVal = 10;
	p1.set_value(Try<int>(initVal));
	EXPECT_EQ(g1.get().get(), plus2(initVal));

	// lift a function that gets deleted, it should be copied into the lifted
	// function
	std::promise<Try<int> > p2;
	Future<int> f2(p2.get_future());
	std::function<Future<bool>(const Future<int>&)> lf;
	{
		auto e = [](const int& i){ return i % 2 == 0; };
		lf = Future<int>::lift(e);
	}
	p2.set_value(Try<int>(22));
	EXPECT_TRUE(lf(f2).get().get());
}

TEST_F(FutureTest, RecoverWith) {
	std::promise<Try<int> > p1;
	int init = 13;
	int alt = 31;
	Future<int> f1(p1.get_future());
	Future<int> f1m =
		f1 | [](const int &) { throw BadIntException(); return 0; };
	Future<int> f1mr =
		f1m.recoverWith(
			[&alt](const std::exception_ptr &){
				return alt;
			});
	EXPECT_EQ(f1mr.waitFor(initWait), std::future_status::timeout);
	p1.set_value(Try<int>(init));
	EXPECT_EQ(f1mr.get().get(), alt);
	EXPECT_THROW(f1m.get().get(), BadIntException);

	// example with selective recovery
	std::promise<Try<int> > p2;
	Future<int> f2(p2.get_future());
	Future<int> f2m =
		f2 | [](const int &) { throw BadIntException(); return 0; };
	Future<int> f2mr =
		f2m.recoverWith(
			[&alt](const std::exception_ptr &e){
				try {
					std::rethrow_exception(e);
					return -1;
				} catch (const BadOuterIntException &) {
					// this won't be executed, and exception will propagate
					return alt;
				}
			});
	p2.set_value(Try<int>(init));
	EXPECT_THROW(f2m.get().get(), BadIntException);
	EXPECT_THROW(f2mr.get().get(), BadIntException);
}
