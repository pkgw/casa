#include <gtest/gtest.h>

#include <casacore/coordinates/Coordinates/CoordinateUtil.h>

#include <imageanalysis/Images/test/ComponentListImageTest.h>
#include <imageanalysis/Images/ComponentListImage.h>

using namespace std;

using namespace casacore;
using namespace casa;

namespace test {

TEST_F(ComponentListImageTest, constructorTest) {
    auto csys = CoordinateUtil::defaultCoords3D();
    IPosition shape(3, 5, 5, 5);
    ComponentList cl;
    ASSERT_THROW(ComponentListImage cli(cl, csys, shape), AipsError);
}

}

int main (int nArgs, char * args []) {
    ::testing::InitGoogleTest(& nArgs, args);
    cout << "ComponentListImage test " << endl;
    return RUN_ALL_TESTS();
}

