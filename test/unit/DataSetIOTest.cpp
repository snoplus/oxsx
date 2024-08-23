#include <catch2/catch_all.hpp>
#include <catch2/catch_approx.hpp>
#include <ROOTTree.h>
#include <IO.h>
#include <OXSXDataSet.h>
#include <TTree.h>
#include <TFile.h>
#include <math.h>
#include <iostream>

TEST_CASE("Writing a data set to disk in the HDF5 format, and reading back"){
    
    // Create fake data set
    OXSXDataSet origDataSet;
    std::vector<double> eventObs(4);

    for(unsigned i = 0; i < 21811; i++){
        for(size_t j = 0; j < eventObs.size(); j++)
            eventObs[j] = j;
        
        origDataSet.AddEntry(Event(eventObs));
    }
    
    // name the observables for the save
    std::vector<std::string> names;
    names.push_back("var0");
    names.push_back("var1");
    names.push_back("var2");    
    names.push_back("var3");
    origDataSet.SetObservableNames(names);
    
    IO::SaveDataSet(origDataSet, "data_set_io_root_test.h5");
    OXSXDataSet* loadedSet = IO::LoadDataSet("data_set_io_root_test.h5");
    size_t nEntries = origDataSet.GetNEntries();

    SECTION("Names copied correctly"){
        REQUIRE(origDataSet.GetObservableNames() == loadedSet->GetObservableNames());
    }

    SECTION("Same Data dimension"){
        REQUIRE(origDataSet.GetNEntries() == loadedSet->GetNEntries());
        REQUIRE(origDataSet.GetNObservables() == loadedSet->GetNObservables());
    }
	SECTION("Same number of events"){
	  REQUIRE(loadedSet->GetNEntries() == nEntries);
	}

    SECTION("Same first and last data"){
        REQUIRE(origDataSet.GetEntry(0).GetData() == loadedSet->GetEntry(0).GetData());
        REQUIRE(origDataSet.GetEntry(nEntries -1).GetData() == loadedSet->GetEntry(nEntries -1).GetData());
    }
    
    remove("data_set_io_root_test.h5");
}

TEST_CASE("Read TTree file in from disk") {
    // First - we must create the TTree file!
    const std::string treename = "T";
    TTree tree(treename.c_str(), "");
    double energy;
    unsigned int nhits;
    tree.Branch("energy", &energy, "energy/D");
    tree.Branch("nhits", &nhits, "nhits/i");
    const std::vector<std::string> observable_names = {"energy", "nhits"};

    for (size_t i = 0; i < 10000; i++) {
        energy = std::sin(i);
        nhits = (i % 13) + (i % 53);
        tree.Fill();
    }
    // Write tree to a file
    const std::string filename = "data_set_io_ttree_test.root";
    TFile outfile(filename.c_str(), "RECREATE");
    tree.Write();
    outfile.Close();

    // Now - let's try and read in that file, using the ROOTTree class!
    ROOTTree intree(filename, treename);

    SECTION("Names copied correctly") {
        REQUIRE(intree.GetObservableNames() == observable_names);
    }
    SECTION("Same data dimension") {
        REQUIRE(intree.GetNObservables() == observable_names.size());
        REQUIRE(intree.GetNEntries() == tree.GetEntries());
    }
    SECTION("All data matches, to within a floating point conversion") {
        tree.SetBranchAddress("energy", &energy);
        tree.SetBranchAddress("nhits", &nhits);
        for (size_t i = 0; i < intree.GetNEntries(); i++) {
            Event evt = intree.GetEntry(i);
            tree.GetEntry(i);
            REQUIRE(evt.GetDatum("energy") == Catch::Approx(energy));
            REQUIRE(evt.GetDatum("nhits") == Catch::Approx(static_cast<double>(nhits)));
        }
    }

    // Finally, delete the ROOT file
    remove(filename.c_str());
}
