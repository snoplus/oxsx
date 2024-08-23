#include <DataSetGenerator.h>
#include <OXSXDataSet.h>
#include <ROOTNtuple.h>
int main()
{
    // (filename, tree name)
    ROOTNtuple zeroNuMC("<filename1>", "treename");
    ROOTNtuple twoNuMC("<filename2>", "treename");
    ROOTNtuple b8NuMC("<filename3>", "treename");

    std::vector<bool> bootstraps;

    // DataSetGenerator
    DataSetGenerator dataGen;
    dataGen.AddDataSet(&zeroNuMC, 5, false); // second arg is the expected rate
    bootstraps.push_back(true);
    dataGen.AddDataSet(&twoNuMC, 50, false);
    bootstraps.push_back(true);
    dataGen.AddDataSet(&b8NuMC, 100, false);
    bootstraps.push_back(true);

    dataGen.SetBootstrap(bootstraps);

    // Generate a data set
    OXSXDataSet fakeData = dataGen.PoissonFluctuatedDataSet();
    OXSXDataSet fakeData2 = dataGen.ExpectedRatesDataSet();

    return 0;
}
