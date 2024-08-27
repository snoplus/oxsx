#include <ConfigLoader.hh>
#include <ContainerTools.hpp>
#include <Exceptions.h>
#include <sstream>

INI::Parser *ConfigLoader::fParser = NULL;

void ConfigLoader::Load(const std::string &section_, const std::string &fieldName_, std::string &loadVal_)
{

    CheckExists(section_, fieldName_);
    loadVal_ = fParser->top()(section_)[fieldName_];
}

void ConfigLoader::Open(const std::string &fileName_)
{
    std::ifstream ss;
    ss.open(fileName_.c_str());
    if (!ss)
    {
        throw IOError("Can't find config file " + fileName_);
    }
    delete fParser;
    fParser = new INI::Parser(ss);
}

void ConfigLoader::Close()
{
    delete fParser;
    fParser = NULL;
}

std::set<std::string>
ConfigLoader::ListSections()
{
    if (!fParser)
        return std::set<std::string>();

    return ContainerTools::GetKeys(fParser->top().sections);
}

void ConfigLoader::CheckExists(const std::string &section_, const std::string &field_)
{

    if (fParser->top().sections.find(section_) == fParser->top().sections.end())
        throw ConfigSectionMissing("Can't find section " + section_);

    if (fParser->top()(section_).values.find(field_) == fParser->top()(section_).values.end())
        throw ConfigFieldMissing("Can't find field " + field_ + " in section " + section_);
}
