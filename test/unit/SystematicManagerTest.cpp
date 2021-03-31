//The point of this is to test the interace no do a fit.
#include <catch.hpp>
#include <BinnedEDManager.h>
#include <DistTools.h>
#include <Gaussian.h>
#include <SystematicManager.h>
#include <Systematic.h>
#include <Convolution.h>
#include <Scale.h>
#include <iostream>
#include <ContainerTools.hpp>
#include <numeric>

class FakeSystematic : public Systematic{
public:
    FakeSystematic(const std::string& name_)  {name = name_; }
    ~FakeSystematic()  {}

    // BinnedED 
    // operator()(const BinnedED& pdf_) const{;}

    // Don't overload these methods as they are good in Systematic.
    // void SetResponse(const SparseMatrix& responseMatrix_){resp= responseMatrix_;}
    // const SparseMatrix& GetResponse() const{return resp;}
        
    // void   SetTransformationObs(const ObsSet&){;}

    // ObsSet GetTransformationObs() const{return ObsSet(3);} 
    //
    // void   SetDistributionObs(const ObsSet&){;}
    // ObsSet GetDistributionObs() const{return ObsSet(3);}

    void Construct() {;}

    void   SetParameter(const std::string& name_, double value) {;}
    double GetParameter(const std::string& name_) const { return 1;}

    void   SetParameters(const ParameterDict&) {;}
    ParameterDict GetParameters() const {ParameterDict temp;return temp;}
    size_t GetParameterCount() const  {;}

    std::set<std::string> GetParameterNames() const {std::set<std::string> temp;return temp;}
    void   RenameParameter(const std::string& old_, const std::string& new_) {;}

    std::string GetName() const {return name;}
    void SetName(const std::string& name_) {name = name_;}
private:
    std::string name;
};

TEST_CASE("SystematicManager"){

    SparseMatrix id(50,50);
    id.SetToIdentity();
    SparseMatrix zero(50,50);
    zero.SetZeros();

    FakeSystematic sys1("sys1");
    sys1.SetResponse(id);
    FakeSystematic sys2("sys2");
    sys2.SetResponse(zero);
    FakeSystematic sys3("sys3");
    sys3.SetResponse(id);
    FakeSystematic sys4("sys4");
    sys4.SetResponse(id);

    AxisCollection axes;
    axes.AddAxis(BinAxis("axis1", -10, 10 ,50));

    BinnedED pdf1("pdf1",axes);
    BinnedED pdf2("pdf2",axes);

    pdf1.SetBinContent(10,0.5);
    pdf1.SetBinContent(40,0.5);

    pdf2.SetBinContent(10,0.5);
    pdf2.SetBinContent(40,0.5);

    std::vector<BinnedED> pdfs;
    pdfs.push_back(pdf1);
    pdfs.push_back(pdf2);

    SystematicManager man;

    man.Add(&sys1);
    man.Add(&sys4);
    man.Add(&sys2,"groupName");
    man.Add(&sys3,"groupName");
    man.Construct();

    std::vector<std::string> appliedGroups;
    appliedGroups.push_back("groupName");
    man.AddDist(pdf2,appliedGroups); 

    std::vector<BinnedED> OrignalPdfs(pdfs);

    man.DistortEDs(OrignalPdfs,pdfs);

    SECTION("Check number of groups"){
        REQUIRE(  man.GetNGroups() == 2 );
    }

    SECTION("Set group names correctly"){
        std::set<std::string> name =  man.GetGroupNames();
        REQUIRE( name.find("groupName") != name.end() );
    }

    SECTION("Counting systematics"){
        int n =  man.GetNSystematics();
        REQUIRE( n == 4);
        n =  man.GetNSystematicsInGroup("");
        REQUIRE( n == 2);
        n =  man.GetNSystematicsInGroup("groupName");
        REQUIRE( n == 2);
    }

    SECTION("Applying identity does nothing"){
        REQUIRE( pdfs.at(0).GetBinContent(10) == 0.5 );
    }

    SECTION("Applying null matrix destroys everything"){
        std::vector<double> vec= pdfs.at(1).GetBinContents();
        REQUIRE( std::accumulate(vec.begin(), vec.end(), 0) == 0 );
    }

    SECTION("Picking up emtpy groups"){
        bool flag = 0;
        try {
            std::vector<std::string> anotherGroup;
            anotherGroup.push_back("anotherGroup");
            man.AddDist(pdf2,anotherGroup); 
            man.DistortEDs(OrignalPdfs, pdfs);
        }catch(const NotFoundError& e_){
            flag=1;
        }

        REQUIRE( flag == 1 );
    }

}
