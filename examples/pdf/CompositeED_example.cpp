/*
  Pdfs containing independent components for example:
        P(E, r, t) = P(E) * P(r) * P(t)
  are dealt with using the CompositeED class. Just create the pdfs 
  independently and then multiply them together to get a composite pdf.

  If you've correctly set the data representation of all the constituents
  you can call compositeED.Probability(anEvent) and each of the consitent
  pdfs will work out what to do!

  In this example we combine a 2D gaussian analytic pdf, 
  with a 1D gaussian analytic pdf to produce  a 3D pdf. 
  ANY OTHER class implementing the PDF interface would work 
  in exactly the same way.
 */

#include <Gaussian.h>
#include <AnalyticED.h>
#include <CompositeED.h>

int main(){
    // Build the constituent parts independently
    Gaussian gaus2d(2); // normal parameters by default
    AnalyticED pdf2d("analytic2D", &gaus2d);
    // 2d pdf will look at observables 1 and 3
    // you can address these by names if they come from a data set
    std::vector<size_t> indices;
    indices.push_back(1);
    indices.push_back(3);

    pdf2d.SetObservables(ObsSet(indices));

    Gaussian gaus1d(1);
    AnalyticED pdf1d("analytic1D", &gaus1d);
    // 1d pdf will look at observable 4
    pdf1d.SetObservables(ObsSet(4));


    // Now we combine them into a single composite pdf
    CompositeED combinedPdf = pdf2d * pdf1d;
    

    std::cout << "From a "<< pdf2d.GetNDims() << "D pdf\n"
              << "and a " << pdf1d.GetNDims() << "D pdf\n"
              << "we have created a " << combinedPdf.GetNDims()
              << "D pdf" << std::endl;

    std::cout << "Like any other pdf it can be normalised " << std::endl;
    
    Event fakeEvent(std::vector<double>(10, 2));

    std::cout << "Calling probability on the composite pdf "
              << "gives you " 
              << combinedPdf.Probability(fakeEvent)
              << "\n individually the pdfs give you " 
              << pdf2d.Probability(fakeEvent) << " and " 
              << pdf1d.Probability(fakeEvent)
              << "\n their product is " 
              << pdf2d.Probability(fakeEvent) * pdf1d.Probability(fakeEvent)
              << "\n so we have combined two independent components!"
              << std::endl;
        
    
}
