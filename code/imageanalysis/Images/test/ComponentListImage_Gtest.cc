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
    ComponentListImage cl0(cl, csys, shape);
    csys = CoordinateUtil::defaultCoords4D();
    shape.resize(4);
    shape = IPosition(4, 5, 5, 5, 5);
    ASSERT_THROW(ComponentListImage(cl, csys, shape), AipsError);
    shape[2] = 4;
    shape[0] = 0;
    ASSERT_THROW(ComponentListImage(cl, csys, shape), AipsError);
    shape[0] = 5;
    ComponentListImage cli(cl, csys, shape);
    ASSERT_FALSE(cli.isPersistent());
}

TEST_F(ComponentListImageTest, constructorTest2) {
    auto csys = CoordinateUtil::defaultCoords4D();
    IPosition shape(4, 4);
    ComponentList cl = oneGaussianCL();
    String imageName = "my.im";
    Unit unit("Jy/pixel");
    Record miscInfo;
    miscInfo.define("me", "you");
    ImageInfo ii;
    String objectName = "Pluto";
    ii.setObjectName(objectName);
    {
        ComponentListImage cli(cl, csys, shape, imageName);
        cli.setUnits(unit);
        cli.setMiscInfo(miscInfo);
        cli.setImageInfo(ii);
        ASSERT_TRUE(cli.isPersistent());
    }
    {
        ComponentListImage cli(imageName);
        ASSERT_TRUE(cli.units() == unit);
        ASSERT_TRUE(cli.miscInfo().asString("me") == "you");
        ASSERT_TRUE(cli.imageInfo().objectName() == objectName);
        ASSERT_TRUE(cli.shape() == shape);
        ASSERT_TRUE(cli.componentList().nelements() == 1);
        ASSERT_TRUE(cli.isPersistent());
    }
    Table::deleteTable(imageName);
}

TEST_F(ComponentListImageTest, copyConstructorTest) {
    auto csys = CoordinateUtil::defaultCoords4D();
    IPosition shape(4, 4);
    ComponentList cl;
    ComponentListImage cli(cl, csys, shape);
    ComponentListImage copy(cli);
    ASSERT_TRUE(cli.shape() == copy.shape());
    ASSERT_TRUE(cli.coordinates().near(copy.coordinates()));
    ASSERT_FALSE(cli.isPersistent());
    ASSERT_FALSE(copy.isPersistent());
}

TEST_F(ComponentListImageTest, cloneIITest) {
    auto csys = CoordinateUtil::defaultCoords4D();
    IPosition shape(4, 4);
    ComponentList cl;
    ComponentListImage cli(cl, csys, shape);
    unique_ptr<ImageInterface<Float>> clone(cli.cloneII());
    ASSERT_TRUE(cli.shape() == clone->shape());
    ASSERT_TRUE(cli.coordinates().near(clone->coordinates()));
    ASSERT_TRUE(cli.imageType() == clone->imageType());
}

TEST_F(ComponentListImageTest, nameTest) {
    auto csys = CoordinateUtil::defaultCoords4D();
    IPosition shape(4, 4);
    ComponentList cl;
    ComponentListImage cli(cl, csys, shape);
    ASSERT_TRUE(cli.name(True) == "Temporary ComponentListImage");
    ASSERT_TRUE(cli.name(False) == "Temporary ComponentListImage");
    String imageName = "nameTest.im";
    {
        ComponentListImage cli2(cl, csys, shape, imageName);
        ASSERT_TRUE(cli2.name(True) == imageName);
    }
    Table::deleteTable(imageName);
}

TEST_F(ComponentListImageTest, setUnitsTest) {
    auto csys = CoordinateUtil::defaultCoords4D();
    IPosition shape(4, 4);
    ComponentList cl;
    ComponentListImage cli(cl, csys, shape);
    ASSERT_THROW(cli.setUnits("Jy/beam"), AipsError);
    ASSERT_TRUE(cli.setUnits("Jy/pixel"));
}

TEST_F(ComponentListImageTest, setImageInfoTest) {
    auto csys = CoordinateUtil::defaultCoords4D();
    IPosition shape(4, 4);
    ComponentList cl;
    ComponentListImage cli(cl, csys, shape);
    auto ii = cli.imageInfo();
    ii.setRestoringBeam(Quantity(5, "arcmin"), Quantity(4, "arcmin"), Quantity(0, "deg"));
    ASSERT_THROW(cli.setImageInfo(ii), AipsError);
    ii.removeRestoringBeam();
    ASSERT_TRUE(cli.setImageInfo(ii));
}

}

int main (int nArgs, char * args []) {
    ::testing::InitGoogleTest(& nArgs, args);
    cout << "ComponentListImage test " << endl;
    return RUN_ALL_TESTS();
}

