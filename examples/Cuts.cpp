/*
Cuts are functions that take an event and return true or false, depending on whether the event passes the cut or not.

Cuts can be grouped together into a cut collection, which returns true if all cuts are satisfied.
The CutCollection object is used to FillPdfs or generate data sets, only taking events that pass 
*/
#include <BoxCut.h>
#include <BoolCut.h>
#include <LineCut.h>
#include <CutCollection.h>
#include <Event.h>
#include <iostream>

int main(){
    // cut if first observable is bigger than two i.e. a lower limit
    LineCut lineCut(0, 2, "lower");
    
    // cut if second observable is between 5 and 10
    BoxCut boxCut(1, 5, 10);

    // cut if the third value is 1 exactly
    BoolCut boolCut(2, 1);

    // Combine the three cuts
    CutCollection combinedCuts;
    combinedCuts.AddCut(lineCut);
    combinedCuts.AddCut(boxCut);
    combinedCuts.AddCut(boolCut);

    // Create a fake event to test
    std::vector<double> observations;
    observations.push_back(1);
    observations.push_back(6);
    observations.push_back(1);

    Event fakeEvent(observations);

    // test it
    std::cout << "Event passes line cut " 
              << lineCut.PassesCut(fakeEvent)
              << std::endl;

    std::cout << "Event passes box cut " 
              << boxCut.PassesCut(fakeEvent)
              << std::endl;

    std::cout << "Event passes bool cut " 
              << boolCut.PassesCut(fakeEvent)
              << std::endl;

    std::cout << "Event passes all cuts " 
              << combinedCuts.PassesCuts(fakeEvent)
              << std::endl;

    return 0;
    
}
