/*
  Demonstrates the use of an Analytic pdf, using a 2D gaussian function
*/
#include <AnalyticED.h>
#include <Gaussian.h>
int main(){
    /* We can use any function that implements the IntegrableFunction 
       interface, here just use a gaussian
     */
    Gaussian gaus(std::vector<double>(2, 0), std::vector<double>(2, 1));
    
    AnalyticED gausPdf("gausPDF", &gaus); // takes ownership of a copy of gaussian

    // The dimensionality of the pdf now matches the underlying function
    std::cout << "Pdf is " << gausPdf.GetNDims() << std::endl;

    // now we can call it directly
    std::cout << " gausPdf(2, 0) = "
              << gausPdf.Probability(std::vector<double>(2, 0)) << std::endl;

    // or use a data representation to query an event
    // this pdf picks out the 0th and 2nd observables in an event
    std::vector<size_t> relevantIndices;
    relevantIndices.push_back(0);
    relevantIndices.push_back(2);
    
    gausPdf.SetObservables(ObsSet(relevantIndices));

    // fake event with ten observables in it, our pdf knows where to look
    Event event(std::vector<double>(10, 0));
    std::cout << "Probability = "
              << gausPdf.Probability(event) << std::endl;

    return 0;
}
