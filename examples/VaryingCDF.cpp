#include <string>        
#include <vector>        
#include <math.h>	
#include <fstream>

#include <Rand.h>
#include <Exceptions.h>
#include <Convolution.h>
#include <Gaussian.h>
#include <ParameterDict.h>
#include <ContainerTools.hpp>
#include <SparseMatrix.h>
#include <DistTools.h>

#include <VaryingCDF.h>
#include <ContainerParameter.h>
#include <Formatter.hpp>
#include <Function.h>
#include <Gaussian.h>

#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

#include <set>
#include <map>
#include <vector>

using ContainerTools::ToString;
using ContainerTools::GetValues;
using ContainerTools::GetKeys;

//Quick function ploy(x) = a*sqrt(abs(x)) + b
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
            return parameters.at("grad")*sqrt(abs(vals_[0]))+parameters.at("offset");
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

int main(int argc, char *argv[])
{
    AxisCollection axes;
    axes.AddAxis(BinAxis("axis1", 10, 30 ,200));

    ObsSet  obsSet(0);
    ObsSet  obsSetToTransform(0);

    Convolution* conv_a = new Convolution("conv_a");
    VaryingCDF smear("smear");

    Gaussian* gaus_a = new Gaussian(0,10e4,"gaus_a"); 
    Ploy* ploy = new Ploy("line",0.5,0);

    gaus_a->RenameParameter("means_0","mean");
    gaus_a->RenameParameter("stddevs_0","std");
    //Set the kernal.
    smear.SetKernel(gaus_a);

    //Parameter std now runs like ploy.
    smear.SetDependance("std",ploy);

    conv_a->SetConditionalPDF(&smear);

    conv_a->SetAxes(axes);
    conv_a->SetTransformationObs(obsSetToTransform);
    conv_a->SetDistributionObs(obsSet);
    conv_a->Construct();

    SparseMatrix mat = conv_a->GetResponse();
    mat.Print("");

    Gaussian gaus(20, 3);
    BinnedED pdf("pdf", DistTools::ToHist(gaus, axes));
    pdf.SetObservables(0);

    pdf.Scale(40000);


    TH1D Pdf = DistTools::ToTH1D(pdf);
    TH1D PdfAfterSmear = DistTools::ToTH1D(conv_a->operator()(pdf));
    TLegend leg(0.6,0.6,0.9,0.9); 
    leg.AddEntry(&Pdf,"Before","lf"); 
    leg.AddEntry(&PdfAfterSmear,"After","lf"); 

    TCanvas c1;
    gStyle->SetOptStat(kFALSE);  
    Pdf.SetTitle("Before/After smear");
    Pdf.SetLineColor(kRed);
    Pdf.SetLineWidth(2);
    Pdf.Scale(1./pdf.Integral());
    Pdf.Draw();
    PdfAfterSmear.SetLineColor(kBlue);
    PdfAfterSmear.SetLineWidth(2);
    PdfAfterSmear.Draw("same");
    leg.Draw();
    c1.Print("Smeared.png");
                   


    return 0;
}
