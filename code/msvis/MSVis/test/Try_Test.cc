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

	struct BadTryIntException : public std::runtime_error {
		BadTryIntException()
			: std::runtime_error("BadTryIntException")
			{}
	};

	static Try<int>
	badTryInt() {
		throw BadTryIntException();
		return Try<int>::from(goodInt);
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
	Try<int> good = Try<int>::from(goodInt);
	Try<int> bad = Try<int>::from(badInt);
	EXPECT_EQ(good.get(), goodInt());
	EXPECT_THROW(bad.get(), BadIntException);
}

TEST_F(TryTest, Flatten) {
	// flatten success of success
	Try<int> good(6);
	Try<Try<int> > a(good);
	EXPECT_EQ(good, a.flatten<int>());

	// flatten success of failure
	Try<int> bad = Try<int>::from(badInt);
	Try<Try<int> > b(bad);
	EXPECT_THROW(b.flatten<int>().get(), BadIntException);

	// flatten failure of success
	Try<Try<int> > c =
		Try<Try<int> >::from([&](){ badTryInt(); return good; });
	EXPECT_THROW(c.flatten<int>().get(), BadTryIntException);

	// flatten failure of failure
	Try<Try<int> > d =
		Try<Try<int> >::from([&](){ badTryInt(); return bad; });
	EXPECT_THROW(d.flatten<int>().get(), BadTryIntException);
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
	EXPECT_EQ(a.map<bool>(isEven), t);

	// create false value
	Try<bool> f = t.map<bool>([](const bool &b){ return !b; });
	EXPECT_EQ(f, Try<bool>(false));

	// another simple map
	Try<int> b = a.map<int>([](const int &i){ return i + 1; });
	EXPECT_EQ(b.map<bool>(isEven), f);

	// map failure
	Try<int> c = Try<int>::from(badInt);
	EXPECT_FALSE(c.map<bool>(isEven).isSuccess());

	// map value to failure
	Try<int> d = a.map<int>([](const int &){ badInt(); return 0; });
	EXPECT_THROW(d.map<bool>(isEven).get(), BadIntException);
}

TEST_F(TryTest, FlatMap) {
	// simple flatMap, w.o. failure
	Try<int> a(10);
	EXPECT_EQ(a.flatMap<int>(goodTryAdd2Int),
	          a.map<int>([](const int& i){ return i + 2; }));

	// flatMap to failure
	EXPECT_THROW(a.flatMap<int>(badTryAdd2Int).get(),
	             BadTryIntException);

	// flatMap from failure
	Try<int> b = Try<int>::from(badInt);
	EXPECT_THROW(b.flatMap<int>(goodTryAdd2Int).get(),
	             BadIntException);

	// flatMap from failure to failure
	EXPECT_THROW(b.flatMap<int>(badTryAdd2Int).get(),
	             BadIntException);
}

TEST_F(TryTest, Fold) {
	// fold on success
	EXPECT_EQ(Try<int>::from(goodInt).fold<std::string>(err, val), pass);

	// fold on failure
	EXPECT_EQ(Try<int>::from(badInt).fold<std::string>(err, val), fail);
}

TEST_F(TryTest, Foreach) {
	// foreach on success
	std::string s = fail;
	Try<int> a(12);
	a.foreach([&s](const int &i) { s = std::to_string(i); });
	EXPECT_EQ(s, std::to_string(a.get()));

	// foreach on failure
	s = fail;
	Try<int>::from(badInt).foreach([&s](const int &) { s = pass; });
	EXPECT_EQ(s, fail);
}

TEST_F(TryTest, Get) {
	// get on success
	int val = 14;
	Try<int> a(val);
	EXPECT_EQ(a.get(), val);

	// get on failure
	Try<int> b = Try<int>::from(badInt);
	EXPECT_THROW(b.get(), BadIntException);

	// throw exception on each get of a failure
	EXPECT_THROW(b.get(), BadIntException);
}

TEST_F(TryTest, GetOrElse) {
	// getOrElse on success
	int val = 16;
	Try<int> a(val);
	EXPECT_EQ(a.getOrElse<int>(goodInt), val);

	// getOrElse on failure
	Try<int> b = Try<int>::from(badInt);
	EXPECT_EQ(b.getOrElse<int>(goodInt), goodInt());
}

TEST_F(TryTest, OrElse) {
	// orElse on success
	Try<int> a(18);
	Try<int> b = a.orElse<int>([&](){ return Try<int>::from(goodInt); });
	EXPECT_EQ(b, a);

	// orElse on failure
	Try<int> c = Try<int>::from(badInt);
	Try<int> d = c.orElse<int>([&a](){ return a; });
	EXPECT_EQ(d, a);
}

TEST_F(TryTest, Transform) {
	// transform on success
	EXPECT_EQ(
		Try<int>::from(goodInt).transform<std::string>(tryErr, tryVal),
		Try<std::string>(pass));

	// transform on failure
	EXPECT_EQ(
		Try<int>::from(badInt).transform<std::string>(tryErr, tryVal),
		Try<std::string>(fail));

	// failed transform on success
	EXPECT_THROW(
		Try<int>::from(goodInt).transform<std::string>(tryErr, tryValFail).get(),
		BadTryXformException);

	// failed transform on failure
	EXPECT_THROW(
		Try<int>::from(badInt).transform<std::string>(tryErrFail, tryVal).get(),
		BadTryXformException);
}

TEST_F(TryTest, Lift) {
	auto liftedIsEven = Try<int>::lift<bool>(isEven);
	Try<int> a(20);
	EXPECT_TRUE(liftedIsEven(a).get());
	EXPECT_FALSE(
		liftedIsEven(a.map<int>([](const int& i){ return i + 1; })).get());
}
