#include <catch2/catch_all.hpp>
#include <BinnedEDGenerator.h>
#include <DistTools.h>
#include <Gaussian.h>
#include <cmath>

TEST_CASE("Recover a 1D gaussian"){
    Gaussian gaus(0, 10);
    AxisCollection axes;
    axes.AddAxis(BinAxis("test", -100, 100, 200));
    
    BinnedED inPdf("test", DistTools::ToHist(gaus, axes));
    
    BinnedEDGenerator edGen;
    edGen.SetRates(std::vector<double> (1, 1000000));
    edGen.SetPdfs(std::vector<BinnedED> (1, inPdf));

    BinnedED outPdf = edGen.ExpectedRatesED();
    REQUIRE(outPdf.Integral() == 1000000);
    outPdf.Normalise();

    const double outMean = outPdf.Means().at(0);
    const double inMean = inPdf.Means().at(0);

    REQUIRE_THAT(outMean, Catch::Matchers::WithinAbs(inMean, 0.01));
}
