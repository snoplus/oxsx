#include <Minuit.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinimize.h>
#include <Minuit2/MnSimplex.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/MnUserCovariance.h>
#include <FitResult.h>
#include <Exceptions.h>
#include <Formatter.hpp>
#include <math.h> //sqrt

using ContainerTools::GetValues;
using ContainerTools::HasSameKeys;
using ContainerTools::GetKeys;
using ContainerTools::ToString;


using ROOT::Minuit2::MnMigrad;
using ROOT::Minuit2::MnMinimize;
using ROOT::Minuit2::MnSimplex;
using ROOT::Minuit2::MnUserParameters;

// FIXME:: CHeck the length of the vectors as you set them match fNparams

Minuit::~Minuit(){
    delete fMinimiser;
}

void
Minuit::SetMethod(const std::string& method_){
    fMethod = method_;   
}

std::string
Minuit::GetMethod() const{
    return fMethod;
}

void
Minuit::SetInitialErrors(const ParameterDict& errs_){
    fInitialErrors = errs_;
}

void
Minuit::SetInitialValues(const ParameterDict& vals_){
    fInitialValues = vals_;
}

void
Minuit::SetUpperContourEdge(double val_){
    fMinuitFCN.SetUp(val_);
}

double
Minuit::GetUpperContourEdge() const{
    return fMinuitFCN.GetUp();
}

void
Minuit::SetMaxima(const ParameterDict& maxima_) {fMaxima = maxima_;}

ParameterDict 
Minuit::GetMaxima() const {return fMaxima;}

void
Minuit::SetMinima(const ParameterDict& minima_) {fMinima = minima_;}

ParameterDict 
Minuit::GetMinima() const {return fMinima;}

void 
Minuit::Initialise(TestStatistic * testStat_){
    delete fMinimiser;
    fMinimiser = NULL;

    // check everything makes sense
    if( !HasSameKeys(fMinima, fMaxima))
        throw LogicError(Formatter()
                         << "Minuit initialisation error "
                         << " minima/maxima parameters dont't match:\n"
                         << "Minima for :\n" << ToString(GetKeys(fMinima)) << "\n"
                         << "Maxima for :\n" << ToString(GetKeys(fMaxima)) << "\n"
                         );

    if( !HasSameKeys(fInitialErrors, fInitialValues))
        throw LogicError(Formatter()
                         << "Minuit initialisation error "
                         << "Initial value/error parameters dont't match:\n"
                         << "Initial Values for :\n" << ToString(GetKeys(fInitialValues)) << "\n"
                         << "Initial Errors for :\n" << ToString(GetKeys(fInitialErrors)) << "\n"
                         );

    for(ParameterDict::iterator it = fMinima.begin(); it != fMinima.end();
        ++it){
        if(fMinima[it->first] >= fMaxima[it->first])
            throw LogicError(Formatter()
                             << "Minuit initialisation error "
                             << "Invalid fit range for " << it->first
                             << " : " << it->second << " - " << fMaxima[it->first]
                             << "\n"
                             );
    }


    // take a copy of these here to make sure we hit the same order each time - minuit uses vector we use set/map...
    fParameterNames = GetKeys(fInitialErrors);

    // pass this ordering on to the called function so it knows the right order too
    fMinuitFCN = MinuitFCN(testStat_, fParameterNames);

    // max or min?
    if(fMaximising)
        fMinuitFCN.SetSignFlip(true);


    // Create parameters and set limits
    MnUserParameters params(GetValues(fInitialValues, fParameterNames), GetValues(fInitialErrors, fParameterNames));

    if(fMinima.size() && fMaxima.size()){
        int i = 0;
        for(std::set<std::string>::iterator it = fParameterNames.begin();
            it != fParameterNames.end(); ++it){
            params.SetLimits(i++, fMinima[*it], fMaxima[*it]);
        }
    }

    if("Migrad" == fMethod)
        fMinimiser = new MnMigrad(fMinuitFCN, params);
    
    else if ("Minimize" == fMethod)
        fMinimiser = new MnMinimize(fMinuitFCN, params);
    
    else if ("Simplex" == fMethod)
        fMinimiser = new MnSimplex(fMinuitFCN, params);
    else 
        throw ValueError(Formatter() << "Minuit::Unrecognised method - "
                         << fMethod
                         );
}

void
Minuit::Fix(const std::string& name_){
    fFixedParameters.insert(name_);
}

void
Minuit::Release(const std::string& name_){
    fFixedParameters.erase(name_);
}

void
Minuit::SetMaxCalls(unsigned max_) {
    fMaxCalls = max_;
}

unsigned
Minuit::GetMaxCalls() const {
    return fMaxCalls;
}

const FitResult&
Minuit::Optimise(TestStatistic* testStat_){
    testStat_ -> RegisterFitComponents();

    Initialise(testStat_);

    if(testStat_->GetParameterNames() != fParameterNames)
        throw LogicError(Formatter() << "Minuit config parameters don't match the test statistic!\n" 
                         << "TestStatistic: " << ToString(testStat_->GetParameterNames()) << "\n"
                         << "Minuit: " << ToString(fParameterNames)
                         );

    // fix the requested parameters
    // first work out which minuit index this refers to
    std::set<std::string>::iterator it = fFixedParameters.begin();
    for(; it != fFixedParameters.end(); ++it){
        size_t pos = std::distance(std::find(fParameterNames.begin(), fParameterNames.end(), *it), fParameterNames.begin());
        fMinimiser -> Fix(pos);
    }

    // defaults are same as ROOT defaults
    ROOT::Minuit2::FunctionMinimum fnMin  = fMinimiser -> operator()(fMaxCalls, fTolerance);

    fFitResult.SetBestFit(ContainerTools::CreateMap(fParameterNames, fMinimiser -> Params()));
    fFitResult.SetValid(fnMin.IsValid());

    if(fMaximising)
        fFitResult.SetExtremeVal(-fnMin.Fval());
    else
        fFitResult.SetExtremeVal(fnMin.Fval());

    DenseMatrix covarianceMatrix = CalcCovarianceMatrix(fnMin);

    fFitResult.SetCovarianceMatrix(covarianceMatrix);

    return fFitResult;
}

void 
Minuit::SetTolerance(double tol_) {
  fTolerance = tol_;
}

double
Minuit::GetTolerance() const {
  return fTolerance;
}

FitResult
Minuit::GetFitResult() const{
    return fFitResult;
    
}

DenseMatrix
Minuit::CalcCovarianceMatrix(ROOT::Minuit2::FunctionMinimum functionMin) {
    ROOT::Minuit2::MnUserCovariance minCovariance = functionMin.UserCovariance();
    const std::vector<double>& minCovarianceVector = minCovariance.ROOT::Minuit2::MnUserCovariance::Data();
    int noVectorEntries = minCovarianceVector.size();
    double noRows = (sqrt((8*noVectorEntries) + 1) - 1) / 2;
    DenseMatrix minCovarianceMatrix (noRows, noRows);
    minCovarianceMatrix.SetSymmetricMatrix(minCovarianceVector);
    return minCovarianceMatrix;
}
