#include <ContainerTools.hpp>

std::vector<std::string>
ContainerTools::Split(const std::string &s, char delim)
{
    std::vector<std::string> elems;
    ContainerTools::Split(s, delim, std::back_inserter(elems));
    return elems;
}
