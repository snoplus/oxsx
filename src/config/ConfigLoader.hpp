#include <sstream>
#include <Exceptions.h>
#include <ContainerTools.hpp>

template<typename TargetType>
void
ConfigLoader::Load(const std::string& section_, const std::string& fieldName_,  TargetType& loadVal_, typename enable_if<is_number<TargetType>::value, int>::type){
  CheckExists(section_, fieldName_);

  loadVal_ = Converter<TargetType>()(fParser->top()(section_)[fieldName_]);
}


template<typename TargetType>
void
ConfigLoader::Load(const std::string& section_, const std::string& fieldName_, TargetType& loadVal_, typename enable_if<is_container<TargetType>::value, int>::type){
    CheckExists(section_, fieldName_);
    ConvertContainer(ContainerTools::Split(fParser->top()(section_)[fieldName_], ','), loadVal_, Converter<TargetType>());
}

template<typename TargetType>
TargetType
Converter<TargetType>::operator()(const std::string& s_) const{
  std::istringstream buffer(s_);
  TargetType val;
  buffer >> val;

  if (buffer.fail())
      throw ValueError("String conversion failed : value is invalid (this can mean the numerical type isn't big enough)");
	
  return val;
}

template<typename InContainer, typename OutContainer, typename ConverterSp>
void
ConvertContainer(const InContainer& incntr_, OutContainer& outcntr_, const ConverterSp& cnvtr_){
    for(typename InContainer::const_iterator it = incntr_.begin(); it != incntr_.end(); ++it){
      outcntr_.insert(outcntr_.end(), *it);
    }
}
