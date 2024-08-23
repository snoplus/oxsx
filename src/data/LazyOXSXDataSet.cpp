#include <LazyOXSXDataSet.h>
#include <IO.h>
#include <algorithm>

LazyOXSXDataSet::~LazyOXSXDataSet()
{
  delete fDataSet;
  Deregister();
}

LazyOXSXDataSet::LazyOXSXDataSet(const LazyOXSXDataSet &other_)
{
  fFileName = other_.fFileName;
  fDataSet = NULL;
}

LazyOXSXDataSet
LazyOXSXDataSet::operator=(const LazyOXSXDataSet &other_)
{
  fFileName = other_.fFileName;
  fDataSet = NULL;
  return *this;
}

unsigned
LazyOXSXDataSet::GetNEntries() const
{
  if (!fDataSet)
    Load();
  return fDataSet->GetNEntries();
}

unsigned
LazyOXSXDataSet::GetNObservables() const
{
  if (!fDataSet)
    Load();
  return fDataSet->GetNObservables();
}

std::vector<std::string>
LazyOXSXDataSet::GetObservableNames() const
{
  if (!fDataSet)
    Load();
  return fDataSet->GetObservableNames();
}

void LazyOXSXDataSet::Load() const
{
  for (size_t i = 0; i < fFiles.size(); i++)
    fFiles[i]->Close();
  fDataSet = IO::LoadDataSet(fFileName);
}

void LazyOXSXDataSet::Close() const
{
  delete fDataSet;
  fDataSet = NULL;
}

void LazyOXSXDataSet::Register()
{
  fFiles.push_back(this);
}

void LazyOXSXDataSet::Deregister()
{
  fFiles.erase(std::remove(fFiles.begin(), fFiles.end(), this), fFiles.end());
}

Event LazyOXSXDataSet::GetEntry(size_t i_) const
{
  if (!fDataSet)
    Load();
  return fDataSet->GetEntry(i_);
}
std::vector<LazyOXSXDataSet *> LazyOXSXDataSet::fFiles;
