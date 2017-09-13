//# Try_Test.h
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
#include <msvis/MSVis/util/Try.h>
#include <gtest/gtest.h>

using namespace casa::vi;

int
main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

class TryTest
	: public ::testing::Test {
public:
	struct BadIntException : public std::runtime_error {
		BadIntException()
			: std::runtime_error("BadIntException")
			{}
	};

	static int
	goodInt() {
		return 42;
	}

	static int
	badInt() {
		throw BadIntException();
		return 0;
	}

	static Try<int>
	goodTryInt() {
		return try_(goodInt);
	}

	struct BadTryIntException : public std::runtime_error {
		BadTryIntException()
			: std::runtime_error("BadTryIntException")
			{}
	};

	static Try<int>
	badTryInt() {
		throw BadTryIntException();
		return try_(goodInt);
	}

	static bool
	isEven(const int &i) {
		return i % 2 == 0;
	}

	static Try<int>
	goodTryAdd2Int(const int &i) {
		return Try<int>(i + 2);
	}

	static Try<int>
	badTryAdd2Int(const int &) {
		return Try<int>(std::make_exception_ptr(BadTryIntException()));
	}

	static const std::string fail;
	static const std::string pass;

	static std::string err(const std::exception_ptr&) { return fail; }

	static std::string val(const int &) { return pass; };

	static Try<std::string> tryErr(const std::exception_ptr&) {
		return Try<std::string>(fail);
	}

	static Try<std::string> tryVal(const int &) {
		return Try<std::string>(pass);
	};

	struct BadTryXformException : public std::runtime_error {
		BadTryXformException()
			: std::runtime_error("BadTryXformException")
			{}
	};

	static Try<std::string> tryErrFail(const std::exception_ptr&) {
		throw BadTryXformException();
		return Try<std::string>(fail);
	}

	static Try<std::string> tryValFail(const int &) {
		throw BadTryXformException();
		return Try<std::string>(pass);
	};

	class Foo {
	public:
		Foo(int i) : m_i(i) {}
		int m_i;
		bool operator==(const Foo& f) const {
			return m_i == f.m_i;
		}
	};

	struct BadRecoveryException : public std::runtime_error {
		BadRecoveryException()
			: std::runtime_error("BadRecoveryException")
			{}
	};
};

const std::string TryTest::fail = "fail";
const std::string TryTest::pass = "pass";

TEST_F(TryTest, Constructors) {
	int four = 4;
	Try<int> a(four);
	EXPECT_TRUE(a.isSuccess());
	EXPECT_EQ(a.get(), four);

	Try<int> c = Try<int>(std::make_exception_ptr(std::runtime_error("bad")));
	EXPECT_FALSE(c.isSuccess());
}

TEST_F(TryTest, Eq) {
	Try<int> a(4);
	Try<int> b(4);
	EXPECT_TRUE(a == b);
	Try<int> c(5);
	EXPECT_TRUE(a != c);
	Try<int> d(std::make_exception_ptr(std::runtime_error("error")));
	EXPECT_FALSE(a == d);
	Try<int> e(d);
	EXPECT_TRUE(d == e);
	Try<int> f(std::make_exception_ptr(std::runtime_error("error")));
	EXPECT_FALSE(d == f);
}

TEST_F(TryTest, From) {
	Try<int> good = try_(goodInt);
	Try<int> bad = try_(badInt);
	EXPECT_EQ(good.get(), goodInt());
	EXPECT_THROW(bad.get(), BadIntException);
}

TEST_F(TryTest, Flatten) {
	// flatten success of success
	Try<int> good(6);
	Try<Try<int> > a(good);
	EXPECT_EQ(good, a.flatten());

	// flatten success of failure
	Try<int> bad = try_(badInt);
	Try<Try<int> > b(bad);
	EXPECT_THROW(b.flatten().get(), BadIntException);

	// flatten failure of success
	Try<Try<int> > c = try_([&](){ badTryInt(); return good; });
	EXPECT_THROW(c.flatten().get(), BadTryIntException);

	// flatten failure of failure
	Try<Try<int> > d = try_([&](){ badTryInt(); return bad; });
	EXPECT_THROW(d.flatten().get(), BadTryIntException);
}

TEST_F(TryTest, Filter) {
	// filter true
	Try<int> a(8);
	EXPECT_EQ(a.filter(isEven), a);

	// filter false
	Try<int> b(9);
	EXPECT_THROW(b.filter(isEven).get(), NoSuchElementException);
}

TEST_F(TryTest, Map) {
	// simple map, w.o. failure
	Try<int> a(8);
	Try<bool> t(true);
	EXPECT_EQ(a.map(isEven), t);

	// create false value
	Try<bool> f = t.map([](const bool &b){ return !b; });
	EXPECT_EQ(f, Try<bool>(false));

	// another simple map
	Try<int> b = a.map([](const int &i){ return i + 1; });
	EXPECT_EQ(b.map(isEven), f);

	// map failure
	Try<int> c = try_(badInt);
	EXPECT_FALSE(c.map(isEven).isSuccess());

	// map value to failure
	Try<int> d = a.map([](const int &){ badInt(); return 0; });
	EXPECT_THROW(d.map(isEven).get(), BadIntException);
}

TEST_F(TryTest, MapOperator) {
	// simple map, w.o. failure
	Try<int> a(8);
	Try<bool> t(true);
	EXPECT_EQ(a | isEven, t);

	// create false value
	Try<bool> f = t | [](const bool &b){ return !b; };
	EXPECT_EQ(f, Try<bool>(false));

	// another simple map
	Try<int> b = a.map([](const int &i){ return i + 1; });
	EXPECT_EQ(b | isEven, f);

	// map failure
	Try<int> c = try_(badInt);
	EXPECT_FALSE((c | isEven).isSuccess());

	// map value to failure
	Try<int> d = a | [](const int &){ badInt(); return 0; };
	EXPECT_THROW((d | isEven).get(), BadIntException);

	// map sequence
	auto a2e = ((a | [](const int &i){ return i + 2; }) | isEven);
	EXPECT_EQ(a2e, t);
}

TEST_F(TryTest, FlatMap) {
	// simple flatMap, w.o. failure
	Try<int> a(10);
	EXPECT_EQ(a.flatMap(goodTryAdd2Int),
	          a.map([](const int& i){ return i + 2; }));

	// flatMap to failure
	EXPECT_THROW(a.flatMap(badTryAdd2Int).get(), BadTryIntException);

	// flatMap from failure
	Try<int> b = try_(badInt);
	EXPECT_THROW(b.flatMap(goodTryAdd2Int).get(), BadIntException);

	// flatMap from failure to failure
	EXPECT_THROW(b.flatMap(badTryAdd2Int).get(), BadIntException);
}

TEST_F(TryTest, FlatMapOperator) {
	// simple flatMap, w.o. failure
	Try<int> a(10);
	EXPECT_EQ(a >>= goodTryAdd2Int,
	          a.map([](const int& i){ return i + 2; }));

	// flatMap to failure
	EXPECT_THROW((a >>= badTryAdd2Int).get(), BadTryIntException);

	// flatMap from failure
	Try<int> b = try_(badInt);
	EXPECT_THROW((b >>= goodTryAdd2Int).get(), BadIntException);

	// flatMap from failure to failure
	EXPECT_THROW((b >>= badTryAdd2Int).get(), BadIntException);

	// flatMap sequence
	auto a2e =
		((a
		  >>= goodTryAdd2Int)
		 >>= [](const int &i){ return Try<bool>(i % 2 == 0); });
	EXPECT_EQ(a2e, Try<bool>(true));
}

TEST_F(TryTest, AndThen) {
	// simple andThen, w.o. failure
	Try<int> a(10);
	EXPECT_EQ(a.andThen(goodTryInt), goodTryInt());

	// andThen to failure
	EXPECT_THROW(a.andThen(badTryInt).get(), BadTryIntException);

	// andThen from failure
	Try<int> b = try_(badInt);
	EXPECT_THROW(b.andThen(goodTryInt).get(), BadIntException);

	// andThen from failure to failure
	EXPECT_THROW(b.andThen(badTryInt).get(), BadIntException);
}

TEST_F(TryTest, AndThenOperator) {
	// simple andThen, w.o. failure
	Try<int> a(10);
	EXPECT_EQ(a >> goodTryInt, goodTryInt());

	// andThen to failure
	EXPECT_THROW((a >> badTryInt).get(), BadTryIntException);

	// andThen from failure
	Try<int> b = try_(badInt);
	EXPECT_THROW((b >> goodTryInt).get(), BadIntException);

	// andThen from failure to failure
	EXPECT_THROW((b >> badTryInt).get(), BadIntException);
}

TEST_F(TryTest, Fold) {
	// fold on success
	EXPECT_EQ(try_(goodInt).fold(err, val), pass);

	// fold on failure
	EXPECT_EQ(try_(badInt).fold(err, val), fail);
}

TEST_F(TryTest, Foreach) {
	// foreach on success
	std::string s = fail;
	Try<int> a(12);
	a.foreach([&s](const int &i) { s = std::to_string(i); });
	EXPECT_EQ(s, std::to_string(a.get()));

	// foreach on failure
	s = fail;
	try_(badInt).foreach([&s](const int &) { s = pass; });
	EXPECT_EQ(s, fail);
}

TEST_F(TryTest, Get) {
	// get on success
	int val = 14;
	Try<int> a(val);
	EXPECT_EQ(a.get(), val);

	// get on failure
	Try<int> b = try_(badInt);
	EXPECT_THROW(b.get(), BadIntException);

	// throw exception on each get of a failure
	EXPECT_THROW(b.get(), BadIntException);
}

TEST_F(TryTest, GetOrElse) {
	// getOrElse on success
	int val = 16;
	Try<int> a(val);
	EXPECT_EQ(a.getOrElse(goodInt), val);
	EXPECT_EQ(a.getOrElse_(goodInt()), val);
	EXPECT_EQ(a.getOrElse_(val), val);

	// getOrElse on failure
	Try<int> b = try_(badInt);
	EXPECT_EQ(b.getOrElse(goodInt), goodInt());
	EXPECT_EQ(b.getOrElse_(goodInt()), goodInt());
}

TEST_F(TryTest, OrElse) {
	// orElse on success
	Try<int> a(18);
	Try<int> b = a.orElse([&](){ return try_(goodInt); });
	Try<int> c = a.orElse_(try_(goodInt));
	EXPECT_EQ(b, a);
	EXPECT_EQ(c, a);

	// orElse on failure
	Try<int> d = try_(badInt);
	Try<int> e = d.orElse([&a](){ return a; });
	Try<int> f = d.orElse_(a);
	EXPECT_EQ(e, a);
	EXPECT_EQ(f, a);
}

TEST_F(TryTest, Transform) {
	// transform on success
	EXPECT_EQ(
		try_(goodInt).transform(tryErr, tryVal),
		Try<std::string>(pass));

	// transform on failure
	EXPECT_EQ(
		try_(badInt).transform(tryErr, tryVal),
		Try<std::string>(fail));

	// failed transform on success
	EXPECT_THROW(
		try_(goodInt).transform(tryErr, tryValFail).get(),
		BadTryXformException);

	// failed transform on failure
	EXPECT_THROW(
		try_(badInt).transform(tryErrFail, tryVal).get(),
		BadTryXformException);
}

TEST_F(TryTest, Lift) {
	auto liftedIsEven = Try<int>::lift(isEven);
	Try<int> a(20);
	EXPECT_TRUE(liftedIsEven(a).getOrElse_(false));
	auto liftedPlusOne = Try<int>::lift([](const int &i){ return i + 1; });
	EXPECT_FALSE(liftedIsEven(liftedPlusOne(a)).getOrElse_(true));
	// lift a function that gets deleted, it should be copied into the lifted
	// function
	std::function<Try<int>(const Try<int>&)> lf;
	{
		auto i2 = [](const int& i){ return i + 2; };
		lf = Try<int>::lift(i2);
	}
	EXPECT_EQ(lf(a), a.map([](const int &i){ return i + 2; }));
}

TEST_F(TryTest, NoDfltCtor) {
	Foo f(88);
	Foo g(88);
	EXPECT_EQ(Try<Foo>(f).get(), g);
}

TEST_F(TryTest, RecoverWith) {
	Try<int> a(20);
	Try<int> b(std::make_exception_ptr(BadIntException()));
	int altValue = 40;
	auto alt = [&altValue](const std::exception_ptr &) { return altValue; };
	EXPECT_EQ(a.recoverWith(alt), a);
	EXPECT_EQ(b.recoverWith(alt), Try<int>(altValue));
	EXPECT_THROW(
		b.recoverWith([](const std::exception_ptr &) {
				throw BadRecoveryException();
				return 0;
			}).get(),
		BadRecoveryException);
}
