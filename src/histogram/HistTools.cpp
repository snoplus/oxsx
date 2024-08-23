#include <HistTools.h>
#include <Histogram.h>
#include <ContainerTools.hpp>

std::map<std::string, Histogram>
HistTools::MakeAllHists(const AxisCollection &axes_,
                        const std::set<std::pair<std::string, std::string>> &combinations_)
{

    typedef std::set<std::pair<std::string, std::string>>::iterator SetIt;
    typedef std::map<std::string, Histogram> HistMap;

    HistMap hists;
    for (SetIt it = combinations_.begin(); it != combinations_.end(); ++it)
    {
        AxisCollection axes;
        axes.AddAxis(axes_.GetAxis(it->first));
        axes.AddAxis(axes_.GetAxis(it->second));

        hists[it->first + "_" + it->second] = (Histogram(axes));
    }
    return hists;
}

void HistTools::FillAllHists(std::map<std::string, Histogram> &hists_, const std::map<std::string, double> &fillVals_)
{
    // note this assumes that the histograms have axes that match the keys of fillvals throws otherwise
    typedef std::map<std::string, Histogram>::iterator It;
    for (It it = hists_.begin(); it != hists_.end(); ++it)
        it->second.Fill(fillVals_);
}
