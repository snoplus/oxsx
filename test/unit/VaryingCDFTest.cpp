#include <catch.hpp>
#include <iostream>

#include <VaryingCDF.h>
#include <Exceptions.h>
#include <ContainerTools.hpp>
#include <ParameterDict.h>
#include <Function.h>
#include <Gaussian.h>
#include <math.h>


using ContainerTools::ToString;
using ContainerTools::GetValues;
using ContainerTools::GetKeys;

class Ploy : public Function{
    public:
        // Constructory things
        Ploy(const std::string& name_,const double grad, const double offset){
            fName=name_;
            parameters["grad"]=grad;
            parameters["offset"]=offset;
        }

        Ploy(const Ploy& copy_){
            fName=copy_.fName;
            parameters=copy_.parameters;
        }


        Ploy& operator=(const Ploy& copy_){
            fName=copy_.fName;
            parameters=copy_.parameters;
            return *this;
        }

        // Probability
        double operator()(const std::vector<double>& vals_) const{
            return parameters.at("grad")*vals_[0]+parameters.at("offset");
        }

        int GetNDims() const{
            return 1;
        }

        Function* Clone() const{
            return static_cast<Function*> (new Ploy(*this));
        }

        void SetParameter(const std::string& name_, double value_){
            parameters[name_]=value_;
        }

        double GetParameter(const std::string& name_) const{
            return parameters.at(name_);
        }

        void SetParameters(const ParameterDict& paraDict_){
            for (ParameterDict::const_iterator function =paraDict_.begin(); function != paraDict_.end(); ++function) {
                std::set<std::string> holder=GetKeys(parameters);

                if(holder.find(function->first)!=holder.end())
                    SetParameter(function->first,function->second);
            }
        }

        ParameterDict GetParameters() const{
            return parameters;
        }

        size_t GetParameterCount() const{
            return 2;
        }

        std::set<std::string> GetParameterNames() const {
            std::set<std::string> names_;
            names_.insert("grad");
            names_.insert("offset");
            return names_;
        }

        void RenameParameter(const std::string& old_, const std::string& new_){
            parameters[new_]=parameters[old_];
            parameters.erase(old_);
        }

        std::string GetName() const{
            return fName;   
        }

        void SetName(const std::string& name_){
            fName= name_;
        }
    private:
        std::string fName;
        ParameterDict parameters;
};

TEST_CASE("Varying CDF 1 Parameter", "[VaryingCDF]"){
    Ploy ploy("ploy",1,0);
    Gaussian gaus(0, 8);
    VaryingCDF smearer("smearer");
    smearer.SetKernel(&gaus);
    smearer.SetDependance("stddevs_0",&ploy);

    SECTION("Check name"){
        REQUIRE(smearer.GetName()=="smearer");
    }
    SECTION("Check probability"){
        REQUIRE(smearer.Integral(std::vector<double>(1,0),std::vector<double>(1,2),std::vector<double>(1,1))==Approx(0.6827));
        REQUIRE(smearer.Integral(std::vector<double>(1,0),std::vector<double>(1,4),std::vector<double>(1,2))==Approx(0.6827));
    }
}

TEST_CASE("Varying CDF 2 parameter", "[VaryingCDF]"){
    Ploy ployMean("ployMean",0,0);
    Ploy ployStd("ployStd",1,0);
    Gaussian gaus(0, 8);
    VaryingCDF smearer("smearer");
    smearer.SetKernel(&gaus);
    smearer.SetDependance("means_0",&ployMean);
    smearer.SetDependance("stddevs_0",&ployStd);

    SECTION("Check probability"){
        REQUIRE(smearer.Integral(std::vector<double>(1,0),std::vector<double>(1,2),std::vector<double>(1,1))==Approx(0.6827));
        REQUIRE(smearer.Integral(std::vector<double>(1,0),std::vector<double>(1,4),std::vector<double>(1,2))==Approx(0.6827));
    }
}

TEST_CASE("Varying CDF 2 parameter TWO", "[VaryingCDF]"){
    Ploy ployMean("ployMean",2,0);
    Ploy ployStd("ployStd",1,0);
    Gaussian gaus(0, 8);
    VaryingCDF smearer("smearer");
    smearer.SetKernel(&gaus);
    smearer.SetDependance("means_0",&ployMean);
    smearer.SetDependance("stddevs_0",&ployStd);

    SECTION("Check probability"){
        // You are centered at 1, the mean shifts by 2 there the gaussian is about 3  with a stddev of 1. Therefore integrate between [2,4].
        REQUIRE(smearer.Integral(std::vector<double>(1,2),std::vector<double>(1,4),std::vector<double>(1,1))==Approx(0.6827));
    }
}

