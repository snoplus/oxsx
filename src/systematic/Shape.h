/****************************************************************/
/* A systematic that applies an arbitrary shape bin-by-bin onto */
/* an event distribution.                                       */
/* This shape is user-provided, and scales each bin by a        */
/* certain factor.                                              */
/****************************************************************/
#ifndef __OXSX_SHAPE__
#define __OXSX_SHAPE__
#include <Systematic.h>
#include <functional>

/*
* Firstly, for my own sanity, create a typedef for a general function that
* takes in parameters + point in an event distribution's space, and returns
* a value: this is a "shape function".
*/
typedef std::function<double(const ParameterDict&,
                         const std::vector<double>&)> ShapeFunction;

class Shape : public Systematic{
public:
    Shape(const std::string& name_) : fShapeFunc(), fParamDict(), fName(name_) {}

    void SetShapeFunction(ShapeFunction& shape_func, std::vector<std::string>& param_names);

    void Construct();

    // FitComponent interface: handling systematic's parameters
    void   SetParameter(const std::string& name_, double value);
    double GetParameter(const std::string& name_) const;

    void   SetParameters(const ParameterDict&);
    ParameterDict GetParameters() const { return fParamDict; }
    size_t GetParameterCount() const { return fParamDict.size(); }

    std::set<std::string> GetParameterNames() const;
    void   RenameParameter(const std::string& old_, const std::string& new_);

    std::string GetName() const { return fName; }
    void SetName(const std::string& name_) { fName = name_; }
private:
    ShapeFunction fShapeFunc;
    ParameterDict fParamDict;
    std::string fName;
};

#endif
